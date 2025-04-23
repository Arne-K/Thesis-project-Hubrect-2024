import mmap
import lmdb
from typing import Optional, List, Tuple
import argparse
import re
import struct
import collections # For OrderedDict
# File and directory handling
import sys
from pathlib import Path
import shutil
import tempfile
import time
# logging and error handling
import logging
import traceback
# Multiprocessing
import multiprocessing as mp

# --- Import Python Utilities ---
# -------------------------------
relative_target_path = Path("utils") / "python_utils"
project_root = Path.cwd()
found_root = None
while True:
    if (project_root / relative_target_path).is_dir():
        found_root = project_root
        break # Found it!
    # Stop if we reach the filesystem root
    if project_root == project_root.parent:
        raise FileNotFoundError(
            f"Could not find the directory structure '{relative_target_path}'"
        )
    # Go one level up for the next iteration
    project_root = project_root.parent
# Add the found project root to sys.path if it's not already there
if found_root:
    path_str = str(found_root)
    if path_str not in sys.path:
        sys.path.append(path_str)

from utils.python_utils import (
    setup_logging,
    setup_interrupt_handling,
    define_multiprocess_chunks,
    read_barcode_file,
    compression_utility
    )

# =======================================
# Setup logging and interruption handling
# =======================================
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('R2_demultiplex')

_interrupt_event = None

def is_interrupted():
    """Check if interruption has been requested."""
    global _interrupt_event
    return _interrupt_event and _interrupt_event.is_set()

# =======================
# FASTQ parsing functions  
# =======================
def encode_fastq_header_id(
    header_line: bytes,
    colon_character: bytes = b':',
    space_character: bytes = b' '
) -> Optional[bytes]:
    """
    1. Extracts the ID line from the FASTQ header (before first space).
    2. Extracts the illumina chip sequencing coordinates from the ID.
    3. Encodes the values into a single 64-bit integer (packed as bytes).

    Assumes Illumina format like: @instrument:run_id:flowcell:lane:tile:x:y read:is_filtered:control_num:sample_num
    Or simply @instrument:run_id:flowcell:lane:tile:x:y
    """
    try:
        # Find the end of the ID part (first space or end of line)
        space_pos = header_line.find(space_character)
        if space_pos != -1:
            id_portion = header_line[:space_pos]
        else:
            # If no space, the whole line (minus potential newline) is the ID.
            if header_line.endswith(b'\n'):
                 id_portion = header_line[:-1]
                 if id_portion.endswith(b'\r'): # Handle Windows CRLF
                     id_portion = id_portion[:-1]
            else:
                id_portion = header_line

        # Split the ID portion by the last 4 colons
        parts = id_portion.rsplit(colon_character, 4)

        # Expecting 5 parts: [prefix, lane, tile, x_coord, y_coord]
        if len(parts) != 5:
             # logger.warning(f"Could not split header ID into 5 parts: {id_portion!r}") # Optional: log malformed ID
            return None

        # parts[0] is everything before lane
        # parts[1] is lane
        # parts[2] is tile
        # parts[3] is x_coord
        # parts[4] is y_coord

        lane = int(parts[1])
        tile = int(parts[2])
        x_coord = int(parts[3])
        y_coord = int(parts[4])

        # Total bits: 2 (lane) + 20 (tile) + 21 (x) + 21 (y) = 64 bits
        # Mask values to fit into allocated bits (use &)
        lane_masked = lane & 0x3           # 2 bits (max value 3)
        tile_masked = tile & 0xFFFFF       # 20 bits (max value 1,048,575)
        x_coord_masked = x_coord & 0x1FFFFF # 21 bits (max value 2,097,151)
        y_coord_masked = y_coord & 0x1FFFFF # 21 bits (max value 2,097,151)

        # Combine using bit shifts into a 64-bit integer
        # Order: lane | tile | x | y
        encoded_value = (
            (lane_masked << 62) |       # Shift lane by 20 + 21 + 21 = 62 bits
            (tile_masked << 42) |       # Shift tile by 21 + 21 = 42 bits
            (x_coord_masked << 21) |    # Shift x by 21 bits
            y_coord_masked              # Y occupies the lowest 21 bits
        )

        # Pack the 64-bit integer into 8 bytes (big-endian) for the LMDB key
        return struct.pack('>Q', encoded_value)

    except (ValueError, IndexError, TypeError):
        # Log parsing errors if needed, e.g.:
        # logger.warning(f"Failed to parse/encode header: {header_line!r}", exc_info=True)
        return None # Return None if parsing or encoding fails

# ===============================
# Functions for building the LMDB  
# ===============================
def extract_barcode_info(filename: str) -> Tuple[str, int]:
    """
    Extract barcode and barcode_id from a filename.
    
    The barcode is an 8-letter DNA sequence that is either:
    1. At the beginning of the filename followed by an underscore, or
    2. Flanked by two underscores
    
    The barcode_id is an integer that is either:
    1. Flanked by two underscores, or
    2. After an underscore and before the file extension (as in the examples)
    
    Args:
        filename (str): The filename to extract information from
    
    Returns:
        tuple: (barcode, barcode_id) where barcode is a string and barcode_id is an integer
    """
    # Define patterns for barcode
    barcode_patterns = [
        r'^([ACGT]{8})_',  # Barcode at the beginning followed by underscore
        r'_([ACGT]{8})_'   # Barcode flanked by two underscores
    ]
    
    # Define patterns for barcode_id
    barcode_id_patterns = [
        r'bcode_(\d+)'        # Barcode_id flanked by two underscores
    ]
    
    # Check for barcode
    barcode = None
    for pattern in barcode_patterns:
        match = re.search(pattern, filename)
        if match:
            barcode = match.group(1)
            break
    
    # Check for barcode_id
    barcode_id = None
    for pattern in barcode_id_patterns:
        match = re.search(pattern, filename)
        if match:
            barcode_id = int(match.group(1))
            break
    
    return barcode, barcode_id

def build_lmdb_subprocess(
        input_file_path_list: List[Path],
        lmdb_path: Path,
        worker_id: int,
        map_size: int = 20 * 1024 * 1024 * 1024,
        header_line_tag: bytes = b'@NS',
        max_readers: int = 20,
        batch_size: int = 200000 
) -> None:
    """
    Processes a list of FASTQ R1 files assigned to this worker to build the LMDB index.
    Optimized for speed using batch processing and efficient LMDB handling.
    """
    env = None
    total_headers_processed = 0
    total_invalid_headers = 0
    temp_decompressed_file_path: Optional[Path] = None # Track temp file for cleanup

    try:
        # --- Open LMDB environment one per worker ---
        env = lmdb.open(
            str(lmdb_path),
            map_size=map_size,
            writemap=True,    # Generally faster for writes, but less crash-resilient without frequent syncs
            sync=False,       # Disable synchronous commits for speed (rely on writemap and eventual close)
            map_async=True,   # Enable asynchronous flushing
            max_readers=max_readers,
            lock=True         # Assume exclusive write access per process or external locking
        )

        # List to hold (key_bytes, value_bytes) tuples for batch insertion
        batch_data = []

        for input_file_path in input_file_path_list:
            try:
                barcode, _ = extract_barcode_info(input_file_path.name)

                if barcode is None:
                    logger.warning(f"Could not extract barcode from {input_file_path.name}. Skipping file.")
                    continue

                # --- Encode barcode one per file ---
                try:
                    barcode_bytes = barcode.encode('utf-8')
                except Exception as e:
                    logger.error(f"Failed to encode barcode '{barcode}' from {input_file_path.name}: {e}")
                    continue # Skip file if barcode encoding fails

                headers_in_file = 0
                invalid_headers_in_file = 0

                # Use a single transaction for potentially multiple batches within a file for efficiency
                # Note: Long-running transactions can increase lock contention if other *readers* exist.
                # Committing per batch strikes a balance. Let's commit per batch_size.
                txn = env.begin(write=True)

                # --- Decompress the input file if needed ---
                # --------------------------------------------
                temp_decompressed_file_path: Optional[Path] = None # Track temp file for cleanup
                if input_file_path.suffix in ['.gz', '.zst']:
                    logger.info(f"Input file {input_file_path.name} is compressed. Decompressing...")
                    start_time = time.time()
                    with tempfile.NamedTemporaryFile(
                        mode='wb', # Need binary mode for compression utility
                        suffix=".fastq", # Keep extension for clarity
                        dir=input_file_path.parent, # Place temp file near output
                        delete=False # Prevent auto-deletion on close, we manage it
                    ) as tmp_out:
                        temp_decompressed_file_path = Path(tmp_out.name)

                    logger.info(f"Decompressing to temporary file: {temp_decompressed_file_path.name}")
                    decompression_success = compression_utility(
                        input_path=input_file_path,
                        output_path=temp_decompressed_file_path,
                        compress=False
                    )

                    if not decompression_success \
                        or not temp_decompressed_file_path.is_file() \
                            or temp_decompressed_file_path.stat().st_size == 0:
                        # Clean up partially created temp file if it exists
                        if temp_decompressed_file_path:
                                temp_decompressed_file_path.unlink(missing_ok=True)
                        raise RuntimeError(f"Failed to decompress {input_file_path} to temporary file.")

                    input_r1_to_process = Path(temp_decompressed_file_path) # Workers will use the decompressed temp file
                    elapsed_time = time.time() - start_time
                    logger.info(f"Decompression complete ({elapsed_time:.2f}s). Processing temporary file.")
                else:
                    logger.info("Input file is not compressed.")
                    input_r1_to_process = Path(input_file_path) # Process original directly

                try:
                    with open(input_r1_to_process, 'rb') as f:
                        # Use iter and next for potentially slightly faster line reading in some CPython versions
                        line_iterator = iter(f)
                        while True: # Process lines until StopIteration
                            try:
                                line = next(line_iterator)

                                # Check for interruption periodically
                                # Check less frequently than original to reduce overhead
                                if headers_in_file > 0 and headers_in_file % 5000 == 0:
                                    if is_interrupted():
                                        logger.warning(f"Processing interrupted during file {input_file_path.name}")
                                        # Commit any pending data before breaking
                                        if batch_data:
                                            with env.begin(write=True) as final_txn: # Use separate txn for safety
                                                final_txn.cursor().putmulti(batch_data, append=True)
                                            batch_data = []
                                        raise InterruptedError("Processing interrupted by signal") # Custom exception or just break

                                # Identify header lines (check start + potentially length)
                                # Using only startswith should be sufficient and fast
                                if line.startswith(header_line_tag):
                                    headers_in_file += 1
                                    total_headers_processed += 1

                                    # Directly encode the header (line already includes potential newline)
                                    encoded_id_bytes = encode_fastq_header_id(line) # Handles potential newline inside

                                    if encoded_id_bytes:
                                        # Add (key, value) tuple to batch list
                                        batch_data.append((encoded_id_bytes, barcode_bytes))
                                    else:
                                        invalid_headers_in_file += 1
                                        total_invalid_headers += 1
                                        # Optional: More aggressive error checking for malformed files
                                        if invalid_headers_in_file >= 5:
                                            logger.error(f"Too many invalid headers ({invalid_headers_in_file}) early in {input_file_path.name}. Stopping processing of this file.")
                                            # Discard current batch for this file
                                            batch_data = []
                                            break # Stop processing this file

                                    # --- Commit the batch if full ---
                                    if len(batch_data) >= batch_size:
                                        # Use cursor().putmulti() for optimized batch writing
                                        # append=True is efficient if keys are roughly sorted (often true for FASTQ)
                                        txn.cursor().putmulti(batch_data, append=False)
                                        batch_data = [] # Clear the batch list
                                        # Commit and start a new transaction for the next batch
                                        txn.commit()
                                        txn = env.begin(write=True)


                                    # Skip the next 3 lines (sequence, +, quality) for speed
                                    try:
                                        next(line_iterator) # Sequence
                                        next(line_iterator) # +
                                        next(line_iterator) # Quality
                                    except StopIteration:
                                        # End of file reached unexpectedly after a header
                                        break # Exit the while loop for this file

                            except StopIteration:
                                # End of file reached normally
                                break # Exit the while loop for this file
                            except Exception as line_err:
                                logger.error(f"Error processing a line in {input_file_path.name}: {line_err}", exc_info=True)
                                # Decide whether to continue or stop processing the file
                                # For robustness, maybe skip the problematic record and continue
                                continue


                    # --- Commit any remaining items in the batch for this file ---
                    if batch_data:
                        txn.cursor().putmulti(batch_data, append=False)
                        batch_data = [] # Clear list
                    txn.commit() # Commit the final batch for the file
                    txn = None # Ensure we don't accidentally use it after commit

                    logger.info(f"Finished {input_file_path.name} - {headers_in_file} headers processed (detected {invalid_headers_in_file} invalid).")

                except InterruptedError:
                    logger.warning(f"Interruption handled for {input_file_path.name}.")
                    # Re-raise or handle as needed for the overall process
                    raise # Or break outer loop if desired

                except Exception as file_proc_err:
                    logger.error(f"Error processing file {input_file_path.name}: {file_proc_err}")
                    logger.error(traceback.format_exc())
                    if txn: # Abort transaction if an error occurred during file processing
                        try:
                            txn.abort()
                            txn = None
                        except Exception as abort_err:
                            logger.error(f"Failed to abort transaction for Worker {worker_id} after error: {abort_err}")
                    # Continue to the next file
            finally:
                if temp_decompressed_file_path and temp_decompressed_file_path.is_file():
                    temp_decompressed_file_path.unlink()


    except Exception as e:
        logger.error(f"A critical error occurred in the main process: {e}")
        logger.error(traceback.format_exc())

    finally:
        # Ensure the environment is properly closed
        if env is not None:
            logger.info(f"Total headers processed by worker: {total_headers_processed}. Total invalid headers detected by worker: {total_invalid_headers}.")
            try:
                # Ensure any final transaction is handled (though it should be committed/aborted)
                if txn:
                    logger.warning(f"Transaction was still active at final close for Worker {worker_id}. Aborting.")
                    txn.abort()

                # Commit any remaining data in batch_data (e.g., if interruption happened outside file loop)
                if batch_data:
                     logger.warning(f"Committing final batch data ({len(batch_data)} items) during close.")
                     try:
                         with env.begin(write=True) as final_txn:
                             final_txn.cursor().putmulti(batch_data, append=False)
                     except Exception as final_commit_err:
                         logger.error(f"Error committing final batch during close: {final_commit_err}")

                env.sync() # Force sync before closing with writemap=True
                env.close()
                logger.info(f"LMDB environment closed for Worker {worker_id}.")
            except Exception as e:
                logger.error(f"Error closing LMDB environment for Worker {worker_id}: {e}")
                logger.error(traceback.format_exc())
        # Clean up the temporary decompressed file if one was created
        if temp_decompressed_file_path and temp_decompressed_file_path.exists():
            logger.info(f"Cleaning up temporary decompressed file: {temp_decompressed_file_path.name}")
            temp_decompressed_file_path.unlink(missing_ok=True)

# =========================================
# Functions for demultiplexing the R2 files  
# =========================================
def demultiplex_r2_subprocess(
        r2_file_path: Path, 
        lmdb_path: Path,
        all_workers_output_dir: Path,
        barcode_file: Path, 
        worker_id: int,
        chunk_start: int, 
        chunk_end: int,
        buffer_size: int = 1 * 1024 * 1024,
        max_open_handles: int = 40) -> None:
    """
    This function should be implemented to read the R2 file, look up IDs in LMDB,
    and write demultiplexed records to output files.
    1. Reads in a header line from the R2 FASTQ file and finds a matching barcode in the LMDB.
    2. Writes the R2 record into a buffer for the corresponding barcode.
    3. If the buffer is full, writes the records to the temporary worker output files corresponding to the barcodes.
    (The worker temp files are later combined into the final output files)

    Parameters
    ----------
    r2_file_path : Path
        Path to the R2 FASTQ file to demultiplex.
    lmdb_path : Path
        Path to the LMDB database containing the barcode mapping.
    all_workers_output_dir : Path
        Temporary output directory for all workers. 
        The function further creates a subdirectory for each worker.
    barcode_file : Path
        Path to the barcode file used for demultiplexing.
        This is used to create temporary output files for each barcode
        and byte value mappings to the output paths.
    worker_id : int
        ID of the worker process.
    chunk_start : int
        Start byte position of the chunk to process.
    chunk_end : int
        End byte position of the chunk to process.
    batch_size : int (default = 1MB)
        Write buffer threshold before flushing.
    max_open_handles : int (default = 40)
        Maximum number of open file handles for output files.
    """
    # Create a temporary directory for this worker's output. Clean up any existing temp directory if it exists
    worker_temp_output_dir = Path(all_workers_output_dir / f"temp_output_worker_{worker_id}")
    if worker_temp_output_dir.exists():
        try:
            shutil.rmtree(worker_temp_output_dir)
        except Exception as e:
            logger.error(f"Failed to clean up existing temporary output directory: {e}")
            return
    try:
        worker_temp_output_dir.mkdir(parents=True, exist_ok=True)
        logger.info(f"Created temporary output directory: {worker_temp_output_dir}")
    except OSError as e:
        logger.error(f"Failed to create temporary output directory: {e}")
        return
    
    # Create temporary output files for each barcode and a dictionary mapping of barcode byte values to file handles
    lru_handles_map = collections.OrderedDict() # LRU cache for file handles
    barcode_bytes_to_path_map = {} # Standard dict for file paths
    barcode_bytes_to_buffer_map = {} # Standard dict for write buffers

    try:
        # Read the list of expected barcodes and their IDs
        input_barcode_list = read_barcode_file(barcode_file)
        if not input_barcode_list:
             logger.error(f"No valid barcode information found in {barcode_file}. Stopping.")
             raise

        # Create the mapping from barcode bytes to the expected temporary output file path
        logger.info(f"Mapping barcode bytes to output paths for {len(input_barcode_list)} barcodes...")
        for barcode_str, barcode_id in input_barcode_list:
            # Encode the barcode string to bytes (using utf-8 is standard)
            barcode_key = barcode_str.encode('utf-8')
            temp_output_file_path = worker_temp_output_dir / f"worker_{worker_id}_{barcode_str}_{barcode_id}_R2.fastq"
            barcode_bytes_to_path_map[barcode_key] = temp_output_file_path
            temp_output_file_path.touch(exist_ok=True)
            if not temp_output_file_path.is_file():
                logger.error(f"Failed to create temporary output file: {temp_output_file_path}")
                raise OSError(f"Failed to create temporary output file: {temp_output_file_path}")
    except Exception as e:
        logger.error(f"Failed during initialization (reading barcode file or mapping paths): {e}", exc_info=True)
        return # Stop worker if initialization fails
    
    # --- Resource Management Variables ---
    mmapped_file = None
    env = None
    r2_file_handle = None # Keep track of the file handle to close it in finally

    try:
        # --- Open R2 File and LMDB ---
        logger.debug(f"Opening R2 file {r2_file_path} and LMDB {lmdb_path}")
        # Create memory map of the R2 input file
        r2_file_handle = open(r2_file_path, 'rb')
        mmapped_file = mmap.mmap(r2_file_handle.fileno(), 0, access=mmap.ACCESS_READ)
        # Open LMDB environment for reading
        env = lmdb.open(str(lmdb_path), readonly=True, lock=False, readahead=False, max_readers=126) # Adjust max_readers if needed

        # --- Find Effective Start Position of the first full record in the chunk ---
        mmapped_file.seek(chunk_start)
        effective_start = -1
        current_search_pos = chunk_start
        # Scan forward from chunk_start to find the beginning ('@') of a FASTQ header
        while current_search_pos < chunk_end:
            mmapped_file.seek(current_search_pos)
            line = mmapped_file.readline()
            if not line: break # Reached end of file before finding a header
            if line.startswith(b'@'):
                effective_start = current_search_pos # Found the start of a record
                break
            # If not a header, advance position to the start of the next line
            current_search_pos = mmapped_file.tell()

        # If no header is found within the chunk, there's nothing to process
        if effective_start == -1:
            logger.warning(f"No FASTQ header found in chunk [{chunk_start}, {chunk_end}).")
            return # Exit cleanly

        # Position the memory map at the start of the first valid record
        mmapped_file.seek(effective_start)
        logger.info(f"Processing chunk from effective start {effective_start} to {chunk_end}")

        # --- Main Processing Loop ---
        # Use a read transaction for LMDB lookups
        with env.begin(write=False, buffers=True) as txn: # buffers=True returns memoryviews
            # Continue reading records as long as the current position is within the chunk bounds
            while mmapped_file.tell() < chunk_end:
                current_record_start_pos = mmapped_file.tell()

                # Read the 4 lines of a FASTQ record
                header_line = mmapped_file.readline()
                # Check for EOF immediately after reading header
                if not header_line: break
                seq_line = mmapped_file.readline()
                plus_line = mmapped_file.readline()
                qual_line = mmapped_file.readline()

                # Check if a complete record was read
                if not seq_line or not plus_line or not qual_line:
                    # This indicates a potentially corrupted file or reaching EOF mid-record
                    logger.warning(f"Incomplete FASTQ record near pos {current_record_start_pos}. Header: {header_line[:50]!r}")
                    break # Stop processing this chunk if file seems corrupt

                # --- LMDB Lookup and Barcode Key Preparation ---
                # Encode the header to get the key for LMDB lookup
                encoded_id = encode_fastq_header_id(header_line)
                barcode_key = None # Initialize barcode key for this record

                if encoded_id:
                    # Look up the encoded ID in the LMDB transaction
                    barcode_bytes_view = txn.get(encoded_id) # Returns memoryview due to buffers=True
                    if barcode_bytes_view:
                        try:
                            # Convert the memoryview to bytes for use as a dictionary key
                            barcode_key = bytes(barcode_bytes_view)
                        except TypeError as e:
                             # Log if conversion fails, but continue processing other records
                             logger.warning(f"Could not convert LMDB value to bytes for header {header_line[:50]!r}. Error: {e}")
                             continue # Skip this record

                # If no valid barcode was found for this record, skip to the next one
                if barcode_key is None:
                    continue

                # --- Get/Manage Output File Handle using LRU Cache ---
                output_handle = lru_handles_map.get(barcode_key)

                if output_handle is not None:
                    # Cache Hit: Handle exists, mark it as most recently used
                    lru_handles_map.move_to_end(barcode_key)
                else:
                    # Cache Miss: Handle not in cache, need to open it
                    # Look up the pre-defined output path for this barcode key
                    output_path = barcode_bytes_to_path_map.get(barcode_key)
                    if output_path:
                        # Check if the cache is full and we need to close the least recently used handle
                        if len(lru_handles_map) >= max_open_handles:
                            # Remove the oldest item (key, handle tuple) from the OrderedDict
                            lru_barcode_key, lru_handle = lru_handles_map.popitem(last=False)
                            try:
                                # Close the evicted file handle
                                lru_handle.close()
                                logger.debug(f"Closed LRU handle for {lru_barcode_key!r}.")
                            except OSError as e:
                                # Log error but continue, as we primarily need to open the new file
                                logger.error(f"Error closing LRU handle for {lru_barcode_key!r}: {e}")

                        # Open the required output file in append-binary mode
                        try:
                            output_handle = open(output_path, "ab")
                            # Add the new handle to the LRU cache (automatically becomes most recent)
                            lru_handles_map[barcode_key] = output_handle
                            logger.debug(f"Opened handle for {barcode_key!r}. Cache size: {len(lru_handles_map)}")
                        except OSError as e:
                            # If opening fails, log error and skip this record as we can't write it
                            logger.error(f"Failed to open output file {output_path}: {e}")
                            continue
                    else:
                        # This indicates a barcode found in LMDB was not in the initial barcode file/map
                        # This might suggest an inconsistency or error in LMDB build or barcode file
                        logger.warning(f"Barcode key {barcode_key!r} from LMDB not found in initial path map. Skipping record.")
                        continue # Skip this record

                # --- Write Record to the appropriate Buffer ---
                # Get the buffer for this barcode, creating it if it's the first time
                output_buffer = barcode_bytes_to_buffer_map.get(barcode_key)
                if output_buffer is None:
                    output_buffer = bytearray()
                    barcode_bytes_to_buffer_map[barcode_key] = output_buffer

                # Append the 4 lines of the FASTQ record to the buffer
                output_buffer.extend(header_line)
                output_buffer.extend(seq_line)
                output_buffer.extend(plus_line)
                output_buffer.extend(qual_line)

                # --- Flush Buffer to File if it Exceeds Threshold ---
                if len(output_buffer) >= buffer_size:
                    try:
                        # Write the buffer contents to the correct open file handle
                        output_handle.write(output_buffer)
                        output_buffer.clear() # Reset the buffer after writing
                    except OSError as write_err:
                        # Log write errors, decide on recovery strategy (e.g., skip, retry, stop?)
                        logger.error(f"Failed writing buffer for barcode {barcode_key!r} to {output_path}: {write_err}")
                        # Clearing buffer might lose data, but prevents infinite loop if disk is full
                        output_buffer.clear()

                # --- End of while loop (finished processing chunk) ---

        # --- Final Flush: Write any remaining data in buffers ---
        logger.info(f"Flushing remaining output buffers ({len(barcode_bytes_to_buffer_map)} buffers)...")
        for barcode_key, output_buffer in barcode_bytes_to_buffer_map.items():
            if len(output_buffer) > 0:
                # Get or re-open the file handle if needed
                output_handle = lru_handles_map.get(barcode_key)
                if not output_handle:
                    # Handle not found in cache - need to reopen
                    output_path = barcode_bytes_to_path_map.get(barcode_key)
                    if output_path:
                        try:
                            # Reopen the file handle since it was evicted from LRU cache
                            output_handle = open(output_path, "ab")
                            # Add back to the LRU cache
                            lru_handles_map[barcode_key] = output_handle
                            logger.debug(f"Reopened handle for final flush: {barcode_key!r}")
                        except OSError as e:
                            logger.error(f"Failed to reopen output file for final flush {output_path}: {e}")
                            continue
                    else:
                        logger.error(f"CRITICAL - Cannot find output path for barcode {barcode_key!r}")
                        continue
                
                # Write the buffer now that we have an output handle
                try:
                    output_handle.write(output_buffer)
                    logger.debug(f"Final flush for {barcode_key!r} ({len(output_buffer)} bytes)")
                except OSError as write_err:
                    logger.error(f"Failed during final buffer write for barcode {barcode_key!r}: {write_err}")
        
        # Clear buffers after final flush attempt
        barcode_bytes_to_buffer_map.clear()


    except Exception as e:
        # Catch any unexpected errors during processing
        logger.error(f"An unexpected error occurred during R2 processing: {e}", exc_info=True)
        # Re-raise the exception to ensure the main process knows the worker failed
        raise

    finally:
        # --- Cleanup: Ensure all resources are closed ---
        logger.info(f"Cleaning up resources...")

        # Close LMDB environment if it was opened
        if env is not None:
            try:
                env.close()
                logger.debug(f"Closed LMDB environment for Worker {worker_id}.")
            except Exception as lmdb_close_err:
                 logger.error(f"Error closing LMDB env: {lmdb_close_err}")

        # Close memory-mapped file if it was created
        if mmapped_file is not None:
            try:
                mmapped_file.close()
                logger.debug(f"Closed memory mapped file for Worker {worker_id}.")
            except Exception as mmap_close_err:
                 logger.error(f"Error closing mmap file: {mmap_close_err}")

        # Close the original R2 file handle if it was opened
        if r2_file_handle is not None and not r2_file_handle.closed:
             try:
                 r2_file_handle.close()
                 logger.debug(f"Closed R2 file handle.")
             except Exception as fh_close_err:
                  logger.error(f"Error closing R2 file handle: {fh_close_err}")


        # Close all file handles remaining in the LRU cache
        logger.info(f"Closing {len(lru_handles_map)} remaining open file handles in LRU cache.")
        # Use items() for clarity, though values() is sufficient if only closing
        for barcode_key, handle in lru_handles_map.items():
            try:
                if not handle.closed:
                    handle.close()
            except Exception as close_err:
                # Log error but continue closing others
                logger.error(f"Error closing file handle for {barcode_key!r}: {close_err}")
        # Clear the map after attempting to close all handles
        lru_handles_map.clear()

# ===============================
# Main R2 demultiplexing function
# ===============================
def main():
    """Main function for the demultiplexer."""
    parser = argparse.ArgumentParser(description="FASTQ R2 Demultiplexer using LMDB")

    # Required arguments
    parser.add_argument("--r1-dmplex-dir", required=True, help="Directory containing R1 demultiplexed FASTQ files")
    parser.add_argument("--r2-file", required=True, help="Path to R2 FASTQ file to demultiplex")
    parser.add_argument("--barcode-file", required=True, help="Path to tab-separated file (.tsv) with barcode sequences in first column")
    parser.add_argument("--dmplex-output-dir", required=True, help="Directory for the demultiplexed R2 output FASTQ files")
    # Optional arguments
    parser.add_argument("--log-dir", help="Directory for log files (default: log to console)")
    parser.add_argument("--processes", type=int, default=4, # Default to number of CPU cores
                        help=f"Number of parallel processes (default: 4)")
    parser.add_argument("--batch-size", type=int, default=100000,
                        help=f"Records per LMDB write batch (default: 100,000)")
    parser.add_argument("--lmdb-map-size-gb", type=int, default=20,
                        help="Estimated max size (GB) for the LMDB database (default: 20)")

    args = parser.parse_args()

    # Convert string paths to Path objects
    r1_dmplex_dir = Path(args.r1_dmplex_dir)
    r2_file = Path(args.r2_file)
    barcode_file = Path(args.barcode_file)
    dmplex_output_dir = Path(args.dmplex_output_dir)
    log_dir = Path(args.log_dir) if args.log_dir else None
    processes_to_use = args.processes # Renamed for clarity
    lmdb_batch_size = args.batch_size
    lmdb_map_size = args.lmdb_map_size_gb * 1024 * 1024 * 1024 # Convert GB to bytes
    temp_decompressed_file_path: Optional[Path] = None # Track temp file for cleanup


    # Initialize variables
    lmdb_path = None

    try:
        #
        # Setup logging and signal handling (only in the main process initially)
        # -----------------------------------------------------------------
        r2_file_basename = r2_file.stem
        r2_file_basename = r2_file_basename.split(".")[0]
        setup_logging(f"R2_demultiplex_{r2_file_basename}", log_dir)
        setup_interrupt_handling() # Setup for main process

        logger.info("===================================================")
        logger.info("Starting R2 Demultiplexer")
        logger.info(f"Python version: {sys.version}")
        logger.info(f"System CPU count: {mp.cpu_count()}")
        logger.info(f"Using {processes_to_use} worker processes")
        logger.info(f"LMDB Batch Size: {lmdb_batch_size:,}")
        logger.info(f"LMDB Map Size: {args.lmdb_map_size_gb} GB")
        logger.info(f"R1 Input Dir: {r1_dmplex_dir}")
        logger.info(f"R2 Input File: {r2_file}")
        logger.info(f"Barcode File: {barcode_file}")
        logger.info(f"Demultiplexed R2 Output Dir: {dmplex_output_dir}")
        logger.info(f"Log Directory: {log_dir if log_dir else 'Console'}")
        logger.info("===================================================")

        # Validate inputs
        if not r1_dmplex_dir.is_dir():
            logger.error(f"R1 directory not found or not a directory: {r1_dmplex_dir}")
            return 1
        if not r2_file.is_file():
            logger.error(f"R2 file not found or not a file: {r2_file}")
            return 1
        if not barcode_file.is_file():
            logger.error(f"Barcode file not found or not a file: {barcode_file}")
            return 1

        # Create output directory if it doesn't exist
        dmplex_output_dir.mkdir(parents=True, exist_ok=True)

        #
        # Map the barcodes and their IDs to input R1 file paths
        # -----------------------------------------------------
        logger.info(f"Reading barcodes from: {barcode_file}")
        input_barcode_info_list = read_barcode_file(barcode_file)
        if not input_barcode_info_list:
             logger.error(f"No valid barcode information found in {barcode_file}. Exiting.")
             return 1
        logger.info(f"Read {len(input_barcode_info_list)} barcode entries.")

        barcode_info_to_r1_file = {}
        logger.info(f"Scanning for R1 files in: {r1_dmplex_dir}")
        found_r1_files = 0
        # Search for input R1 files
        for pattern in ["*.dmplx.fastq*", "*.dmplx.fq*"]:
            for r1_file_path in r1_dmplex_dir.rglob(pattern):
                # Basic check to skip files likely R2 
                if '_R2_' in r1_file_path.name:
                    continue
                # Extract barcode info (barcode DNA string and barcode_id integer) from filename
                r1_barcode_info = extract_barcode_info(r1_file_path.name)
                if r1_barcode_info:
                     # Check if the extracted info matches one from the barcode file
                    if r1_barcode_info in input_barcode_info_list:
                        if r1_barcode_info in barcode_info_to_r1_file:
                             logger.warning(f"Duplicate R1 file found for barcode {r1_barcode_info}: {r1_file_path.name} and {barcode_info_to_r1_file[r1_barcode_info].name}. Using the first one found.")
                        else:
                            barcode_info_to_r1_file[r1_barcode_info] = r1_file_path
                            found_r1_files += 1
                    # else: # Optional: Log files found but not matching barcode list
                    #    logger.debug(f"R1 file found but barcode {r1_barcode_info} not in list: {r1_file_path.name}")

        # Check if any matching R1 files were found
        if not barcode_info_to_r1_file:
            logger.error(f"No R1 files found in {r1_dmplex_dir} that match the barcodes/IDs in {barcode_file}")
            return 1
        else:
            logger.info(f"Found {len(barcode_info_to_r1_file)} matching R1 files for barcodes.")
            if len(barcode_info_to_r1_file) < len(input_barcode_info_list):
                 logger.warning(f"Could not find matching R1 files for {len(input_barcode_info_list) - len(barcode_info_to_r1_file)} barcodes listed in {barcode_file}")

        #
        # Divide the R1 file paths between the worker processes
        # -----------------------------------------------------
        r1_file_paths = list(barcode_info_to_r1_file.values())
        num_r1_files = len(r1_file_paths)

        if num_r1_files < processes_to_use:
            logger.warning(f"Number of R1 files ({num_r1_files}) is less than the number of processes ({processes_to_use}). Reducing processes to {num_r1_files}.")
            processes_to_use = num_r1_files

        worker_assignments = [[] for _ in range(processes_to_use)]
        for i, file_path in enumerate(r1_file_paths):
             worker_assignments[i % processes_to_use].append(file_path)

        logger.info(f"Distributing {num_r1_files} R1 files among {processes_to_use} LMDB worker processes.")
        for i, assignment in enumerate(worker_assignments):
             logger.debug(f"Worker {i} assigned {len(assignment)} files.")

        #
        # Build the LMDB that maps encoded IDs to barcodes
        # ------------------------------------------------
        lmdb_path = dmplex_output_dir / f"{r1_dmplex_dir.name}_readID_to_barcode.lmdb"
        logger.info(f"Preparing to build LMDB database at: {lmdb_path}")

        # Clean up existing LMDB if it exists
        if lmdb_path and lmdb_path.exists():
             logger.warning(f"Existing LMDB found at {lmdb_path}. Removing it before rebuilding.")
             try:
                 shutil.rmtree(lmdb_path)
             except OSError as e:
                 logger.error(f"Error removing existing LMDB directory {lmdb_path}: {e}")
                 return 1 # Stop if cleanup fails

        # Ensure the parent directory exists (LMDB open will create the final dir)
        lmdb_path.parent.mkdir(parents=True, exist_ok=True)

        processes_list = []
        logger.info(f"Starting {processes_to_use} worker processes to build LMDB...")

        # --- Start LMDB build worker processes ---
        for i in range(processes_to_use):
            worker_files = worker_assignments[i]
            if not worker_files: # Should not happen with the distribution logic above, but check anyway
                logger.info(f"Skipping LMDB-Worker-{i} as it has no assigned files.")
                continue

            # Define kwargs for the target function build_lmdb_subprocess
            process_kwargs = {
                'input_file_path_list': worker_files,
                'lmdb_path': lmdb_path,
                'worker_id': i,
                'map_size': lmdb_map_size,
                'header_line_tag': b'@NS',
                'batch_size': lmdb_batch_size
                # 'max_readers': 20 # add this if needed, otherwise default is used
            }

            p = mp.Process(
                target=build_lmdb_subprocess,
                kwargs=process_kwargs,
                name=f"LMDB-Worker-{i}"
            )
            processes_list.append(p)
            p.start()
            logger.info(f"Started {p.name} with {len(worker_files)} files.")

        # --- Wait for all LMDB builder processes to finish ---
        logger.info("Waiting for LMDB build processes to complete...")
        successful_completion = True
        interrupted_during_wait = False
        try:
            for p in processes_list:
                p.join() # Wait for the process to terminate
                if is_interrupted(): # Check if main process was interrupted while waiting
                    interrupted_during_wait = True
                    # Don't break immediately, allow logging of exit codes if possible
                if p.exitcode != 0:
                    logger.error(f"Process {p.name} exited with non-zero code: {p.exitcode}")
                    successful_completion = False # Mark failure if any worker fails
                else:
                     logger.info(f"Process {p.name} finished successfully.")

        except KeyboardInterrupt: # Catch interrupt specifically during the join loop
            logger.warning("Interruption detected while waiting for workers. Signaling workers to stop.")
            if _interrupt_event:
                 _interrupt_event.set() # Signal workers
            # Assume join timeout handles it.
            interrupted_during_wait = True
            successful_completion = False # Mark as unsuccessful due to interruption

        if interrupted_during_wait:
             logger.error(f"LMDB build process was interrupted. Removing potentially incomplete/corrupted LMDB at {lmdb_path}")
             shutil.rmtree(lmdb_path, ignore_errors=True)
             return 1 # Exit due to interruption

        if not successful_completion:
             logger.error("One or more LMDB worker processes failed. Aborting.")
             shutil.rmtree(lmdb_path, ignore_errors=True)
             return 1 # Exit with error code

        # --- Final check for LMDB existence ---
        if not lmdb_path.exists() or not (lmdb_path / "data.mdb").is_file():
             logger.error(f"LMDB database was not created correctly at {lmdb_path} despite workers reporting success.")
             return 1
        logger.info(f"LMDB database successfully created at {lmdb_path}")

        # --- Decompress the input R2 file if needed ---
        # -----------------------------------------------
        temp_decompressed_file_path: Optional[Path] = None # Track temp file for cleanup
        if r2_file.suffix in ['.gz', '.zst']:
            logger.info(f"Input file {r2_file.name} is compressed. Decompressing...")
            start_time = time.time()
            with tempfile.NamedTemporaryFile(
                mode='wb', # Need binary mode for compression utility
                suffix=".fastq", # Keep extension for clarity
                dir=r2_file.parent, # Place temp file near output
                delete=False # Prevent auto-deletion on close, we manage it
            ) as tmp_out:
                temp_decompressed_file_path = Path(tmp_out.name)

            logger.info(f"Decompressing to temporary file: {temp_decompressed_file_path.name}")
            decompression_success = compression_utility(
                input_path=r2_file,
                output_path=temp_decompressed_file_path,
                compress=False
            )

            if not decompression_success \
                or not temp_decompressed_file_path.is_file() \
                    or temp_decompressed_file_path.stat().st_size == 0:
                # Clean up partially created temp file if it exists
                if temp_decompressed_file_path:
                        temp_decompressed_file_path.unlink(missing_ok=True)
                raise RuntimeError(f"Failed to decompress {r2_file} to temporary file.")
            else:
                input_r2_to_process = Path(temp_decompressed_file_path) # Workers will use the decompressed temp file
                elapsed_time = time.time() - start_time
                logger.info(f"Decompression complete ({elapsed_time:.2f}s). Processing temporary file.")
        else:
            logger.info("Input file is not compressed.")
            input_r2_to_process = Path(r2_file) # Process original directly
        
        # Demultiplex the R2 files using the built LMDB database
        # ------------------------------------------------------
        logger.info("=" * 20 + " Starting R2 Demultiplexing Phase " + "=" * 20)
        worker_assignment_list = define_multiprocess_chunks(input_r2_to_process, num_processes=processes_to_use)
        logger.info(f"Distributing R2 file {input_r2_to_process.name} into {len(worker_assignment_list)} chunks for {processes_to_use} worker processes.")

        # --- Start R2 demultiplex worker processes ---
        processes_list = []    
        for i in range(processes_to_use):
            start_pos, end_pos = worker_assignment_list[i] # Tuple of start and end byte positions in the R2 file
            # Output directory for the temporary worker demultiplexed R2 output
            r2_dmplex_output_dir = Path(dmplex_output_dir) / f"R2_dmplex_temp_output"
            # Define args for the target function
            process_kwargs = {
                'r2_file_path': input_r2_to_process,                   # R2 file to demultiplex (Path object)
                'lmdb_path': lmdb_path,                 # Path object for the LMDB directory containing the barcode mapping
                'all_workers_output_dir': r2_dmplex_output_dir,      # Output directory for all workers
                'barcode_file': barcode_file,             # Path to the barcode file
                'worker_id': i,                         # Worker ID (integer)
                'chunk_start': start_pos,                 # Start byte position of the chunk to process
                'chunk_end': end_pos,                   # End byte position of the chunk to process
                # buffer_size              # FASTQ record writing buffer size for commits (in bytes; default = 1 MB)
                # max_open_handles: int = 40
            }

            p = mp.Process(
                target=demultiplex_r2_subprocess,
                kwargs=process_kwargs,
                name=f"R2-Demultiplex-Worker-{i}" # Assign a worker name
            )
            processes_list.append(p)
            p.start()

        # --- Wait for all R2 demultiplexing processes to finish ---
        logger.info("Waiting for R2 demultiplexing processes to complete...")
        successful_completion = True
        interrupted_during_wait = False
        try:
            for p in processes_list:
                p.join() # Wait for the process to terminate
                if is_interrupted(): # Check if main process was interrupted while waiting
                    interrupted_during_wait = True
                    # Don't break immediately, allow logging of exit codes if possible
                if p.exitcode != 0:
                    logger.error(f"Process {p.name} exited with non-zero code: {p.exitcode}")
                    successful_completion = False # Mark failure if any worker fails
                else:
                     logger.info(f"Process {p.name} finished successfully.")

        except KeyboardInterrupt: # Catch interrupt specifically during the join loop
            logger.warning("Interruption detected while waiting for workers. Signaling workers to stop.")
            if _interrupt_event:
                 _interrupt_event.set() # Signal workers
            # Assume join timeout handles it.
            interrupted_during_wait = True
            successful_completion = False # Mark as unsuccessful due to interruption
        
        if interrupted_during_wait:
            logger.error(f"R2 demultiplex process was interrupted. Cleaning up temporary files.")
            # Cleanup temporary files and terminate processes
            shutil.rmtree(r2_dmplex_output_dir, ignore_errors=True)
            for p in processes_list:
                if p.is_alive():
                    p.terminate()
            return 1
        
        # --- Clean up the temporary decompressed file if one was created ---
        if temp_decompressed_file_path and temp_decompressed_file_path.exists():
            logger.info(f"Cleaning up temporary decompressed file: {temp_decompressed_file_path.name}")
            temp_decompressed_file_path.unlink(missing_ok=True)

        # --- Delete output files that are empty ---
        r2_dmplex_output_file_paths = list(r2_dmplex_output_dir.glob("**/worker*.fastq"))
        for temp_file in r2_dmplex_output_file_paths:
            if temp_file.is_file() and temp_file.stat().st_size == 0:
                try:
                    temp_file.unlink()
                except OSError as unlink_err:
                    logger.warning(f"Could not delete empty temporary file {temp_file}: {unlink_err}")
        
        # --- Merge the temporary output of the R2 demultiplexing workers into a single file for each barcode ---
        temp_missing_barcode_list = []
        logger.info("Merging temporary output files for each barcode...")
        for barcode, barcode_id in input_barcode_info_list:
            # Search recursively within the temp directory structure
            temp_output_files = list(r2_dmplex_output_dir.glob(f"**/worker_*_{barcode}_{barcode_id}_*.fastq"))
            if not temp_output_files:
                temp_missing_barcode_list.append(barcode)
                continue
            r2_file_prefix = r2_file.stem.split(".")[0]
            merged_output_file_path = dmplex_output_dir / f"bcode_{barcode_id}_{barcode}_{r2_file_prefix}.dmplx.fastq"

            try:
                # Merge into the final destination file
                with open(merged_output_file_path, 'wb') as merged_file:
                    for temp_file in temp_output_files:
                        if temp_file.is_file() and temp_file.stat().st_size > 0:
                            with open(temp_file, 'rb') as f:
                                shutil.copyfileobj(f, merged_file)
                        else:
                            logger.warning(f"Skipping empty or missing temporary file: {temp_file}")
                # Checkpoint: Verify merged file is not empty before deleting source temp files
                merged_size = merged_output_file_path.stat().st_size
                if merged_size > 0:
                    logger.info(f"Successfully merged {len(temp_output_files)} temporary files ({merged_size} bytes) for barcode {barcode} ({barcode_id}) into {merged_output_file_path.name}")
                    # Deletion loop
                    for temp_file in temp_output_files:
                        try:
                            temp_file.unlink()
                        except OSError as unlink_err:
                            logger.warning(f"Could not delete temporary file {temp_file}: {unlink_err}")
                else:
                    # Merged file is empty - log warning and do NOT delete temp files
                    logger.warning(f"Merged file {merged_output_file_path.name} for barcode {barcode} is empty. Temporary files were NOT deleted.")
                    # Delete the empty merged file itself
                    try:
                        merged_output_file_path.unlink()
                    except OSError as empty_unlink_err:
                        logger.warning(f"Could not delete empty merged file {merged_output_file_path.name}: {empty_unlink_err}")
            
            except Exception as merge_err:
                logger.error(f"Error during merging for barcode {barcode}: {merge_err}")
                logger.error(f"Temporary files for barcode {barcode} were NOT deleted due to merge error.")
        
        # Log the barcodes that were found in the input file but for which no temporary output files were created
        # Likely due to no matching records in the R2 file
        if temp_missing_barcode_list:
            formatted_barcode_list = '\n'.join([', '.join(temp_missing_barcode_list[i:i+5])
                           for i in range(0, len(temp_missing_barcode_list), 5)])
            formatted_barcode_list = '\n' + formatted_barcode_list
            logger.warning(f"Missing temporary output files for {len(temp_missing_barcode_list)} input barcodes: {formatted_barcode_list}")
        
        # If we've reached this point without returning an error code, consider it a success
        return 0
    
    except KeyboardInterrupt:
        logger.warning("Interrupted by user in main process, shutting down gracefully.")
        # Set interrupt event to signal any potential running workers (though LMDB build should be done)
        if _interrupt_event:
            _interrupt_event.set()
        return 1 # Exit code for interruption

    except Exception as e:
        logger.error(f"An unexpected error occurred in the main function: {e}")
        logger.error(traceback.format_exc())
        return 1 # Exit code for general error
    
    # --- Cleanup any remaining resources or temporary files ---
    finally:
        # --- Cleanup the LMDB database if it exists ---
        if lmdb_path and lmdb_path.is_dir():
            try:
                shutil.rmtree(lmdb_path, ignore_errors=True)
                logger.info(f"Removed LMDB database at {lmdb_path}")
            except OSError as lmdb_cleanup_err:
                logger.error(f"Error cleaning up LMDB database: {lmdb_cleanup_err}")
        
        # --- Cleanup the decompressed files if they were created ---
        if temp_decompressed_file_path and temp_decompressed_file_path.exists():
            logger.info(f"Cleaning up temporary decompressed file: {temp_decompressed_file_path.name}")
            temp_decompressed_file_path.unlink(missing_ok=True)
        
        # --- Cleanup the R2 output directory if it does not contain any files ---
        if 'r2_dmplex_output_dir' in locals():
            if not any(item.is_file() for item in r2_dmplex_output_dir.rglob('*')):
                try:
                    shutil.rmtree(r2_dmplex_output_dir, ignore_errors=True)
                    logger.info(f"Removed temporary output directory: {r2_dmplex_output_dir}")
                except OSError as cleanup_err:
                    logger.error(f"Error cleaning up temporary output directory: {cleanup_err}")

if __name__ == "__main__":
    # Ensure multiprocessing context is set up early if needed (e.g., for 'spawn' method)
    # mp.set_start_method('spawn') # Uncomment if needed, e.g., on macOS/Windows or for specific libraries
    exit_code = main()
    if exit_code == 0:
        logger.info("R2 Demultiplexer finished successfully.")
    else:
        logger.error(f"R2 Demultiplexer exited with error code {exit_code}.")
    sys.exit(exit_code)