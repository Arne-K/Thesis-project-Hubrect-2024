#!/usr/bin/env python3

import argparse
import logging
import shutil
import sys
import time
import tempfile # Added for safe temporary file creation
from pathlib import Path
from typing import Set, Optional, List
import mmap
import multiprocessing as mp
from multiprocessing import Manager


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
    # logging
    setup_logging,
    setup_interrupt_handling, # Assuming this sets up signal handlers
    format_short_path,
    # Barcode
    read_barcode_file,
    # Multiprocessing
    define_multiprocess_chunks,
    # FASTQ parsing
    count_fq_headers, # This function caused the TypeError
    compression_utility
    )


# =======================================
# Setup logging and interruption handling
# =======================================
# Basic config as fallback if setup_logging isn't called or fails
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('fastq_match_edit')

# Placeholder for interrupt event - setup_interrupt_handling should populate this
# Example: _interrupt_event = mp.Event()
_interrupt_event = None

def is_interrupted():
    """Check if interruption has been requested."""
    # This function relies on setup_interrupt_handling correctly setting _interrupt_event
    global _interrupt_event
    return _interrupt_event and _interrupt_event.is_set()

# --- Constants ---
NEWLINE_BIN = b'\n'
DEFAULT_QUALITY_CHAR = 'I'  # Phred+33 score of 40 (high quality)
# Define buffer size once
BUFFER_SIZE = 16 * 1024 * 1024  # 16 MB in bytes

# --- Helper Functions ---
def _load_query_sequences(query_seq_str: Optional[str],
                          query_seq_file: Optional[Path]) -> Set[bytes]:
    """
    Loads query sequences from command-line string and/or a file.

    Args:
        query_seq_str: Comma-separated string of sequences.
        query_seq_file: Path to a tab-separated file (1st column is sequence).

    Returns:
        A set of unique query sequences as bytes.
    """
    query_sequences_bin = set()
    valid_chars = set('ACGTN') # Define valid chars once

    # Load from string
    if query_seq_str:
        sequences = [seq.strip().upper() for seq in query_seq_str.split(',') if seq.strip()]
        loaded_count = 0
        for seq in sequences:
            if all(c in valid_chars for c in seq):
                    query_sequences_bin.add(seq.encode('ascii'))
                    loaded_count += 1
            else:
                    logger.warning(f"Ignoring invalid sequence from --query-seq: {seq}")
        logger.info(f"Loaded {loaded_count} valid sequences from --query-seq argument "
                    f"({len(query_sequences_bin)} unique).")

    # Load from file
    if query_seq_file:
        if not query_seq_file.is_file():
            logger.error(f"Query sequence file not found: {query_seq_file}")
            raise FileNotFoundError(f"Query sequence file not found: {query_seq_file}")
        try:
            logger.info(f"Loading query sequences from file: {format_short_path(query_seq_file)}")
            # Assuming read_barcode_file returns list of (barcode_str, other_info)
            barcode_info_list = read_barcode_file(query_seq_file)
            # Extract the sequence strings from the list of tuples
            file_sequences = [barcode for barcode, _ in barcode_info_list]
            initial_count = len(query_sequences_bin)
            valid_file_seq_count = 0
            ignored_file_seq_count = 0
            for seq in file_sequences:
                seq_upper = seq.strip().upper()
                if seq_upper: # Ensure not empty
                    if all(c in valid_chars for c in seq_upper):
                        query_sequences_bin.add(seq_upper.encode('ascii'))
                        valid_file_seq_count += 1
                    else:
                        # Log only a few examples to avoid flooding logs
                        if ignored_file_seq_count < 5:
                            logger.warning(f"Ignoring invalid sequence from file {query_seq_file.name}: {seq}")
                        ignored_file_seq_count += 1
            if ignored_file_seq_count > 5:
                 logger.warning(f"Ignored a total of {ignored_file_seq_count} invalid sequences from file {query_seq_file.name}.")
            newly_added = len(query_sequences_bin) - initial_count
            logger.info(f"Loaded {valid_file_seq_count} valid sequences from file "
                        f"({newly_added} added, {len(query_sequences_bin)} total unique).")

        except Exception as e:
            logger.error(f"Failed to read or parse query sequence file {query_seq_file}: {e}")
            raise

    if not query_sequences_bin:
        logger.error("No valid query sequences were loaded. Cannot proceed.")
        raise ValueError("No valid query sequences provided or loaded.")

    return query_sequences_bin

def _process_fastq_subprocess_wrapper(
    input_fastq_path: Path, # Renamed for clarity
    temp_worker_dir: Path,
    query_sequences_bin: Set[bytes],
    mode: str,
    byte_range: tuple[int, int],
    replacement_seq_bin: Optional[bytes] = None,
    header_start_tag=None,
    results_queue=None,
    worker_id=None):
    """
    Wrapper function to call the subprocess, handle exceptions, and put results in the queue.
    """
    worker_name = f"Worker-{worker_id:02d}" # Consistent naming
    try:
        # Setup logging within the worker if needed (e.g., different format/level)
        # logger.info(f"{worker_name}: Starting processing.") # Optional: worker start log

        match_count, miss_count = _process_fastq_subprocess(
            input_fastq_path=input_fastq_path,
            temp_worker_dir=temp_worker_dir,
            query_sequences_bin=query_sequences_bin,
            mode=mode,
            byte_range=byte_range,
            replacement_seq_bin=replacement_seq_bin,
            header_start_tag=header_start_tag,
            worker_name=worker_name # Pass worker name for logging
        )
        # Put results in queue
        results_queue.put({
            'worker_id': worker_id,
            'match_count': match_count,
            'miss_count': miss_count,
            'error': None # Indicate success
        })
        # logger.info(f"{worker_name}: Processing finished.") # Optional: worker end log
    except Exception as e:
        error_message = f"{worker_name} error: {type(e).__name__}: {str(e)}"
        # Log full traceback for unexpected errors in worker is helpful
        logger.exception(f"{worker_name} encountered an unhandled exception:")
        # Put error result in queue to avoid blocking and signal failure
        results_queue.put({
            'worker_id': worker_id,
            'match_count': 0,
            'miss_count': 0,
            'error': error_message # Include error message
        })
        # Do not re-raise here; allow main process to handle based on queue results


def _process_fastq_subprocess(
    input_fastq_path: Path,
    temp_worker_dir: Path,
    query_sequences_bin: Set[bytes],
    mode: str,
    byte_range: tuple[int, int],
    replacement_seq_bin: Optional[bytes] = None,
    header_start_tag=None,
    worker_name="Worker"): # Added for logging context
    """
    Processes a chunk of the input FASTQ file using memory mapping.

    Args:
        input_fastq_path: Path to the (potentially decompressed) input FASTQ file.
        temp_worker_dir: Directory for this worker's temporary output files.
        query_sequences_bin: Set of query sequences (bytes).
        mode: Processing mode ('filter' or 'replace').
        byte_range: Tuple (start_byte, end_byte) defining the chunk to process.
        replacement_seq_bin: The replacement sequence (bytes) if mode is 'replace'.
        header_start_tag: Bytes indicating the start of a FASTQ header line.
        worker_name: Identifier for logging messages.

    Returns:
        Tuple (match_count, miss_count) for the processed chunk.
    """
    if not temp_worker_dir.is_dir():
        # This should ideally not happen if main() creates dirs correctly
        raise FileNotFoundError(f"{worker_name}: Temporary worker directory not found: {temp_worker_dir}")
    if header_start_tag is None:
        raise ValueError(f"{worker_name}: Header start tag is None")

    worker_match_count = 0
    worker_miss_count = 0
    worker_processed_records = 0
    log_interval_records = 3_000_000 # Log progress every N records

    chunk_start, chunk_end = byte_range
    logger.info(f"{worker_name}: Processing bytes {chunk_start:,} to {chunk_end:,}")

    # Prepare replacement quality string if needed (same length as replacement sequence)
    # Ensure replacement sequence and quality include the newline for writing
    replacement_qual_bin_with_nl = None
    replacement_seq_bin_with_nl = None
    if mode == 'replace' and replacement_seq_bin:
        qual_char_bin = DEFAULT_QUALITY_CHAR.encode('ascii')
        # Create quality string *without* newline first for length matching
        replacement_qual_bin = qual_char_bin * len(replacement_seq_bin)
        replacement_qual_bin_with_nl = replacement_qual_bin + NEWLINE_BIN
        replacement_seq_bin_with_nl = replacement_seq_bin + NEWLINE_BIN # Add newline for writing

    # --- Initialize worker output ---
    fq_base = input_fastq_path.stem # Use stem from the actual input being processed
    # Use unique names based on worker ID (passed via temp_worker_dir name convention)
    worker_id_str = temp_worker_dir.name.split('_')[-1] # Extract worker ID
    match_output_file = temp_worker_dir / f"{fq_base}_match_temp_{worker_id_str}.fastq"
    miss_output_file = temp_worker_dir / f"{fq_base}_miss_temp_{worker_id_str}.fastq"

    # Ensure files are clean before starting
    match_output_file.unlink(missing_ok=True)
    miss_output_file.unlink(missing_ok=True)

    # Open output files for writing with buffering
    # Use try-with-resources (with statement) for automatic closing
    try:
        with open(match_output_file, 'wb', buffering=BUFFER_SIZE) as match_out_handle, \
             open(miss_output_file, 'wb', buffering=BUFFER_SIZE) as miss_out_handle, \
             open(input_fastq_path, 'rb') as input_file_handle: # No buffering needed here due to mmap

            # Memory map the input FASTQ file
            # Using length=0 maps the entire file. Access is read-only.
            with mmap.mmap(input_file_handle.fileno(), length=0, access=mmap.ACCESS_READ) as mmapped_file:

                # --- Position the worker to its chunk start ---
                # Find the start of the first complete record at or after chunk_start
                current_pos = chunk_start
                if chunk_start > 0:
                    # Seek near the start and find the beginning of the next line
                    mmapped_file.seek(chunk_start - 1)
                    # Readline moves to the end of the current line
                    mmapped_file.readline()
                    current_pos = mmapped_file.tell() # Now at the start of a line

                # Scan forward from current_pos to find the start of a FASTQ header
                effective_start = -1
                while current_pos < chunk_end:
                    mmapped_file.seek(current_pos)
                    line = mmapped_file.readline()
                    if not line: break # Reached end of file
                    if line.startswith(header_start_tag):
                        effective_start = current_pos # Found the start of a record
                        break
                    # If not a header, advance position to the start of the next line
                    current_pos = mmapped_file.tell()

                # If no header is found within the chunk range, there's nothing to process
                if effective_start == -1 or effective_start >= chunk_end:
                    logger.warning(f"{worker_name}: No FASTQ header found to start processing in chunk range "
                                   f"[{chunk_start:,}, {chunk_end:,}).")
                    return 0, 0 # Exit cleanly

                # Position the memory map at the start of the first valid record for this chunk
                mmapped_file.seek(effective_start)

                # --- Main Processing Loop with Batch Writing ---
                # Process records as long as the *start* of the record is within the chunk bounds
                BATCH_SIZE_RECORDS = 100_000 # Process N records before writing
                match_records_batch = []
                miss_records_batch = []
                last_log_time = time.time()

                while True:
                    # Check for interruption *before* reading next record
                    if is_interrupted():
                        logger.warning(f"{worker_name}: Interruption detected. Stopping processing.")
                        break

                    current_record_start_pos = mmapped_file.tell()
                    # Stop if the start of the next record is outside the designated chunk end
                    if current_record_start_pos >= chunk_end:
                        break

                    # Read the 4 lines of a FASTQ record
                    header_line = mmapped_file.readline()
                    # Check for EOF immediately after reading header
                    if not header_line:
                        logger.info(f"{worker_name}: Reached EOF.")
                        break # End of file

                    seq_line = mmapped_file.readline()
                    plus_line = mmapped_file.readline()
                    qual_line = mmapped_file.readline()

                    # Check if a complete record was read
                    if not seq_line or not plus_line or not qual_line:
                        # This indicates a potentially corrupted file or reaching EOF mid-record
                        logger.warning(f"{worker_name}: Incomplete FASTQ record near byte offset "
                                       f"{current_record_start_pos:,}. Header: {header_line[:50]!r}. Stopping chunk.")
                        # Treat as end of processing for this chunk
                        break

                    # Ensure header starts correctly (sanity check)
                    if not header_line.startswith(header_start_tag):
                         logger.warning(f"{worker_name}: Expected header starting with {header_start_tag!r} "
                                        f"but found {header_line[:50]!r} at offset {current_record_start_pos:,}. Skipping record.")
                         continue # Skip this potentially malformed record

                    worker_processed_records += 1

                    # Check for match (use the sequence line excluding the newline)
                    seq_bin = seq_line.rstrip(NEWLINE_BIN) # More specific than rstrip()
                    if seq_bin in query_sequences_bin:
                        worker_match_count += 1
                        if mode == 'filter':
                            match_records_batch.append((header_line, seq_line, plus_line, qual_line))
                        elif mode == 'replace':
                            # Use the pre-calculated replacement seq/qual with newlines
                            match_records_batch.append((header_line, replacement_seq_bin_with_nl, plus_line, replacement_qual_bin_with_nl))
                    else:
                        # No match, collect for mismatch file
                        worker_miss_count += 1
                        miss_records_batch.append((header_line, seq_line, plus_line, qual_line))

                    # Write batches periodically to avoid memory buildup
                    if len(match_records_batch) >= BATCH_SIZE_RECORDS:
                        # *** FIX: Pass iterable of lines directly to writelines ***
                        items_to_write = [item for record in match_records_batch for item in record]
                        match_out_handle.writelines(items_to_write)
                        match_records_batch = [] # Clear the list

                    if len(miss_records_batch) >= BATCH_SIZE_RECORDS:
                        # *** FIX: Pass iterable of lines directly to writelines ***
                        items_to_write = [item for record in miss_records_batch for item in record]
                        miss_out_handle.writelines(items_to_write)
                        miss_records_batch = [] # Clear the list

                    # Log progress periodically based on record count
                    if worker_processed_records % log_interval_records == 0:
                        current_time = time.time()
                        rate = log_interval_records / (current_time - last_log_time + 1e-6) # Avoid division by zero
                        logger.info(f"{worker_name}: Processed {worker_processed_records:,} records... ({rate:,.0f} records/s)")
                        last_log_time = current_time

                # Write any remaining records in the last batch
                if match_records_batch:
                    # *** FIX: Pass iterable of lines directly to writelines ***
                    items_to_write = [item for record in match_records_batch for item in record]
                    match_out_handle.writelines(items_to_write)
                if miss_records_batch:
                    # *** FIX: Pass iterable of lines directly to writelines ***
                    items_to_write = [item for record in miss_records_batch for item in record]
                    miss_out_handle.writelines(items_to_write)

    # Handles are closed automatically by 'with' statement, even if errors occur
    except Exception as e:
        # Log the exception with traceback here before raising
        logger.exception(f"{worker_name}: Unexpected error during FASTQ processing loop.")
        raise # Re-raise to be caught by the wrapper

    finally:
        # Final log message for the worker
        logger.info(f"{worker_name}: Finished processing chunk. "
                    f"Total records processed in chunk: {worker_processed_records:,}. "
                    f"Matches: {worker_match_count:,}, Mismatches: {worker_miss_count:,}")

        # Clean up empty output files (optional, but good practice)
        try:
            if match_output_file.exists() and match_output_file.stat().st_size == 0:
                logger.info(f"{worker_name}: Removing empty match file: {match_output_file.name}")
                match_output_file.unlink()
            if miss_output_file.exists() and miss_output_file.stat().st_size == 0:
                logger.info(f"{worker_name}: Removing empty mismatch file: {miss_output_file.name}")
                miss_output_file.unlink()
        except OSError as e:
             logger.warning(f"{worker_name}: Error during empty file cleanup: {e}")
        # NOTE: Do NOT clean up the temp_worker_dir here. Main process handles it.

    # Return counts for the caller (wrapper function)
    return worker_match_count, worker_miss_count


# --- Main Execution ---
def main(input_fastq_orig: Path, # Keep original path separate
         output_prefix: str,
         operational_mode: str,
         query_seq_str: Optional[str],
         query_seq_file: Optional[Path],
         replacement_seq: Optional[str],
         log_dir: Optional[Path] = None,
         output_dir: Optional[Path] = None,
         keep_original: bool = False,
         processes: int = 4,
         header_start_tag_str: str = None): # Allow configuring header start via arg
    """
    Main function for the FASTQ Match/Edit Tool.
    Handles setup, multiprocessing, file management, concatenation, and validation.
    """
    global _interrupt_event # Allow modification if setup_interrupt_handling needs it

    # --- Validate input arguments ---
    if operational_mode not in ['filter', 'replace']:
        raise ValueError(f"Invalid operational mode: {operational_mode}")
    if not input_fastq_orig.is_file():
         raise FileNotFoundError(f"Input FASTQ file not found: {input_fastq_orig}")
    # Basic FASTQ validation (can be slow for large files, consider skipping or sampling)
    # if not validate_fastq(input_fastq_orig):
    #     logger.warning(f"Input file {input_fastq_orig} might not be a valid FASTQ.")
        # raise ValueError(f"Not a valid FASTQ file based on initial check: {input_fastq_orig}")

    # Validate header start tag
    if not header_start_tag_str or not isinstance(header_start_tag_str, str):
        raise ValueError("Header start tag must be a non-empty string.")
    header_start_tag_bin = header_start_tag_str.encode('ascii')

    # --- Initialize output paths ---
    # Determine output directory
    if output_dir:
        final_output_dir = output_dir
    else:
        final_output_dir = input_fastq_orig.parent
    # Ensure output directory exists
    try:
        final_output_dir.mkdir(parents=True, exist_ok=True)
        logger.info(f"Using output directory: {final_output_dir}")
    except OSError as e:
        logger.error(f"Error creating output directory {final_output_dir}: {e}")
        sys.exit(1)

    # Define final output filenames
    miss_filename = f"{output_prefix}.seq_mismatch.fastq"
    if operational_mode == 'filter':
        match_filename = f"{output_prefix}.seq_match.fastq"
    else: # operational_mode == 'replace'
        match_filename = f"{output_prefix}.seq_replace.fastq"

    output_match_path = final_output_dir / match_filename
    output_miss_path = final_output_dir / miss_filename
    if output_match_path.is_file():
        output_match_path.unlink()
    if output_miss_path.is_file():
        output_miss_path.unlink()

    # --- Setup Logging and Interrupt Handling ---
    log_file_prefix = f"{output_prefix}.fq_match_edit"
    try:
        # Pass the global event object if needed by setup_interrupt_handling
        # setup_interrupt_handling(_interrupt_event) # Example
        setup_logging(log_file_prefix, log_dir) # Configure logging (file or console)
        setup_interrupt_handling() # Setup signal handlers (ensure it works with multiprocessing)
    except Exception as e:
        logger.error(f"Failed to set up logging or interrupt handling: {e}")
        # Continue without advanced logging/interrupts if necessary, or exit
        # sys.exit(1)

    logger.info("===================================================")
    logger.info("Starting FASTQ Match/Edit Tool")
    logger.info(f"Python version: {sys.version.split()[0]}")
    logger.info(f"Operational mode: {operational_mode}")
    logger.info(f"Input File: {format_short_path(input_fastq_orig)}")
    logger.info(f"Output Prefix: {output_prefix}")
    logger.info(f"Output Directory: {final_output_dir}")
    logger.info(f"Log Directory: {log_dir if log_dir else 'Console Only'}")
    logger.info(f"Match/Replaced Output: {output_match_path.name}")
    logger.info(f"Mismatch Output: {output_miss_path.name}")
    logger.info(f"Keep Original Input: {keep_original}")
    logger.info(f"Processes: {processes}")
    logger.info(f"Header Start Tag: {header_start_tag_str!r}")

    # --- Load Query Sequences ---
    logger.info("Loading query sequences...")
    try:
        query_sequences_bin = _load_query_sequences(query_seq_str, query_seq_file)
        logger.info(f"Total unique valid query sequences loaded: {len(query_sequences_bin)}")
    except (FileNotFoundError, ValueError) as e:
         logger.error(f"Failed to load query sequences: {e}")
         return 1 # Exit with error code

    if operational_mode == 'replace':
        if not replacement_seq:
             logger.error("Replacement sequence (--replacement-seq) is required for 'replace' mode.")
             return 1
        replacement_seq_bytes = replacement_seq.encode('ascii')
        logger.info(f"Replacement sequence: {replacement_seq}")
    else:
        replacement_seq_bytes = None # Not needed for filter mode

    logger.info(f"System CPU count: {mp.cpu_count()}")
    logger.info("===================================================")

    # --- State Variables ---
    input_fastq_to_process = input_fastq_orig # Path to the file workers will read
    temp_decompressed_file_path: Optional[Path] = None # Track temp file for cleanup
    workers_dir: Optional[Path] = None # Path to the main temporary dir for workers
    success = False # Overall success flag
    total_match_count = 0
    total_miss_count = 0
    input_record_count = 0
    worker_results = [] # Store results from queue
    process_list: List[mp.Process] = [] # Define process_list here to be accessible in finally

    try:
        # --- Decompress Input File if Needed (using temp file) ---
        if input_fastq_orig.suffix in ['.gz', '.zst']:
            logger.info(f"Input file {input_fastq_orig.name} is compressed. Decompressing...")
            start_time = time.time()
            # Create a temporary file in the output directory (or system temp)
            with tempfile.NamedTemporaryFile(
                mode='wb', # Need binary mode for compression utility
                suffix=".fastq", # Keep extension for clarity
                dir=final_output_dir, # Place temp file near output
                delete=False # Prevent auto-deletion on close, we manage it
            ) as tmp_out:
                temp_decompressed_file_path = Path(tmp_out.name)

            logger.info(f"Decompressing to temporary file: {temp_decompressed_file_path.name}")
            decompression_success = compression_utility(
                input_path=input_fastq_orig,
                output_path=temp_decompressed_file_path,
                compress=False
            )

            if not decompression_success or not temp_decompressed_file_path.is_file() or temp_decompressed_file_path.stat().st_size == 0:
                # Clean up partially created temp file if it exists
                if temp_decompressed_file_path:
                     temp_decompressed_file_path.unlink(missing_ok=True)
                raise RuntimeError(f"Failed to decompress {input_fastq_orig} to temporary file.")

            input_fastq_to_process = temp_decompressed_file_path # Workers will use the decompressed temp file
            elapsed_time = time.time() - start_time
            logger.info(f"Decompression complete ({elapsed_time:.2f}s). Processing temporary file.")
        else:
             logger.info("Input file is not compressed.")
             input_fastq_to_process = input_fastq_orig # Process original directly


        # --- Count Records and Define Chunks ---
        logger.info(f"Counting records in: {input_fastq_to_process.name}")
        start_time = time.time()
        # *** Call count_fq_headers without the unsupported header_start_tag argument ***
        input_record_count = count_fq_headers(input_fastq_to_process)
        elapsed_time = time.time() - start_time
        if input_record_count == 0:
            logger.warning(f"Input file {input_fastq_to_process.name} contains 0 records. Nothing to process.")
            success = True # Technically successful, just no work done
            return 0 # Exit gracefully

        logger.info(f"Input file contains {input_record_count:,} records (counted in {elapsed_time:.2f}s).")

        # Define byte chunks for multiprocessing
        mp_chunks = define_multiprocess_chunks(file_path=input_fastq_to_process, num_processes=processes)
        logger.info(f"Divided input into {len(mp_chunks)} chunks for {processes} processes.")

        # --- Setup Multiprocessing ---
        results_queue = mp.Queue() # Queue to collect results dicts from workers
        # process_list defined above try block

        # Create a *main* temporary directory for all worker outputs
        workers_dir = final_output_dir / f"fq_match_edit_workers_temp_{output_prefix}"
        if workers_dir.exists():
            logger.warning(f"Removing existing temporary worker directory: {workers_dir}")
            shutil.rmtree(workers_dir)
        workers_dir.mkdir(parents=True, exist_ok=False)
        logger.info(f"Created main temporary worker directory: {workers_dir}")

        temp_output_dir_list: List[Path] = [] # Keep track of individual worker dirs

        # --- Start Worker Processes ---
        logger.info("Starting worker processes...")
        for i, chunk in enumerate(mp_chunks):
            # Create a dedicated temp directory for this worker's output
            tmp_dir_path = workers_dir / f"worker_{i:03d}"
            try:
                tmp_dir_path.mkdir(exist_ok=False)
                temp_output_dir_list.append(tmp_dir_path)
            except OSError as e:
                logger.error(f"Failed to create temporary directory for worker {i}: {e}")
                raise RuntimeError(f"Failed to create worker temp dir {tmp_dir_path}") from e

            # Create and start a process for each chunk
            p = mp.Process(
                target=_process_fastq_subprocess_wrapper,
                kwargs={
                    'input_fastq_path': input_fastq_to_process,
                    'temp_worker_dir': tmp_dir_path,
                    'query_sequences_bin': query_sequences_bin,
                    'mode': operational_mode,
                    'byte_range': chunk,
                    'replacement_seq_bin': replacement_seq_bytes,
                    'header_start_tag': header_start_tag_bin,
                    'results_queue': results_queue,
                    'worker_id': i
                },
                name=f"Worker-{i:02d}"
            )
            p.start()
            process_list.append(p)
            logger.debug(f"Started process {p.name} (PID {p.pid}) for chunk {i}.")

        # --- Wait for Workers and Collect Results ---
        logger.info("Waiting for worker processes to complete...")
        for p in process_list:
            p.join()
            logger.debug(f"Process {p.name} (PID {p.pid}) joined.")
        logger.info("All worker processes finished. Collecting results...")
        
        # --- Cleanup the temporary decompressed file if one was created ---
        if temp_decompressed_file_path and temp_decompressed_file_path.exists():
            logger.info(f"Cleaning up temporary decompressed file: {temp_decompressed_file_path.name}")
            temp_decompressed_file_path.unlink(missing_ok=True)
        
        num_workers_finished = 0
        worker_errors_found = False
        while num_workers_finished < len(process_list):
             try:
                 result = results_queue.get(timeout=10)
                 worker_results.append(result)
                 total_match_count += result['match_count']
                 total_miss_count += result['miss_count']
                 if result['error']:
                     logger.error(f"Error reported by Worker-{result['worker_id']}: {result['error']}")
                     worker_errors_found = True
                 num_workers_finished += 1
             except mp.queues.Empty:
                 logger.warning("Results queue is empty, but not all workers have reported. Waiting...")
                 alive_procs = [p.name for p in process_list if p.is_alive()]
                 if alive_procs:
                      logger.error(f"Timeout waiting for results. Processes still alive: {alive_procs}. Aborting.")
                      worker_errors_found = True
                      break
                 else:
                      logger.warning("Queue empty and all processes joined. Possible issue?")
                      if num_workers_finished < len(process_list):
                           worker_errors_found = True
                           logger.error(f"Missing results from {len(process_list) - num_workers_finished} workers.")
                      break

        if worker_errors_found:
            logger.error("Errors occurred in one or more worker processes. Processing failed.")
            success = False
        else:
            # --- Preliminary Count Check ---
            total_processed_by_workers = total_match_count + total_miss_count
            logger.info(f"Total records processed by workers: {total_processed_by_workers:,}")
            logger.info(f"Total matches reported: {total_match_count:,}")
            logger.info(f"Total mismatches reported: {total_miss_count:,}")

            if total_processed_by_workers != input_record_count:
                 logger.warning(f"Mismatch between initial record count ({input_record_count:,}) and "
                               f"total processed by workers ({total_processed_by_workers:,}). "
                               f"This might indicate issues with chunking, record parsing or the initial count method.")

            success = True


    except (FileNotFoundError, ValueError, RuntimeError, MemoryError) as e:
        logger.error(f"Error during setup or processing: {e}")
        success = False
    except KeyboardInterrupt:
         logger.warning("Keyboard interrupt detected. Attempting graceful shutdown...")
         if _interrupt_event:
             _interrupt_event.set()
         time.sleep(2)
         for p in process_list:
             if p.is_alive():
                 logger.warning(f"Terminating process {p.name} (PID {p.pid}).")
                 p.terminate()
                 p.join(timeout=5)
         success = False
         logger.error("Processing aborted due to keyboard interrupt.")
    except Exception as e:
        logger.exception(f"An unexpected error occurred in the main process: {e}")
        success = False

    finally:
        # --- Cleanup and Final Steps ---
        logger.debug("Entering final cleanup phase...")

        if success:
            # --- Concatenate Worker Outputs ---
            logger.info("Concatenating worker output files...")
            try:
                # Concatenate MATCH files
                if workers_dir and workers_dir.is_dir():
                    worker_match_files = sorted(list(workers_dir.glob("worker_*/*_match_temp_*.fastq")))
                    if worker_match_files:
                        logger.debug(f"Found {len(worker_match_files)} temporary match files to concatenate.")
                        output_match_path.unlink(missing_ok=True)
                        with open(output_match_path, 'wb', buffering=BUFFER_SIZE) as f_out:
                            for i, worker_file in enumerate(worker_match_files):
                                logger.debug(f"Concatenating match file {i+1}/{len(worker_match_files)}: {worker_file.name}")
                                with open(worker_file, 'rb', buffering=BUFFER_SIZE) as f_in:
                                    shutil.copyfileobj(f_in, f_out, length=BUFFER_SIZE)
                        logger.debug(f"Successfully created final match/replaced file: {output_match_path}")
                        if not output_match_path.is_file() or output_match_path.stat().st_size == 0:
                             logger.warning(f"Final match file {output_match_path} appears empty or was not created.")
                             if total_match_count > 0:
                                  logger.error("Concatenation failed for match file (expected matches).")
                                  success = False

                    elif total_match_count > 0:
                         logger.error(f"Expected {total_match_count} matches, but found no temporary match files in {workers_dir}.")
                         success = False
                    else:
                         logger.info("No temporary match files found (as expected, 0 matches reported).")
                         output_match_path.touch(exist_ok=True)


                    # Concatenate MISMATCH files (only if concatenation hasn't failed yet)
                    if success:
                        worker_miss_files = sorted(list(workers_dir.glob("worker_*/*_miss_temp_*.fastq")))
                        if worker_miss_files:
                            logger.debug(f"Found {len(worker_miss_files)} temporary mismatch files to concatenate.")
                            output_miss_path.unlink(missing_ok=True)
                            with open(output_miss_path, 'wb', buffering=BUFFER_SIZE) as f_out:
                                for i, worker_file in enumerate(worker_miss_files):
                                    logger.debug(f"Concatenating mismatch file {i+1}/{len(worker_miss_files)}: {worker_file.name}")
                                    with open(worker_file, 'rb', buffering=BUFFER_SIZE) as f_in:
                                        shutil.copyfileobj(f_in, f_out, length=BUFFER_SIZE)
                            logger.debug(f"Successfully created final mismatch file: {output_miss_path}")
                            if not output_miss_path.is_file() or output_miss_path.stat().st_size == 0:
                                logger.warning(f"Final mismatch file {output_miss_path} appears empty or was not created.")
                                if total_miss_count > 0:
                                    logger.error("Concatenation failed for mismatch file (expected mismatches).")
                                    success = False

                        elif total_miss_count > 0:
                            logger.error(f"Expected {total_miss_count} mismatches, but found no temporary mismatch files in {workers_dir}.")
                            success = False
                        else:
                            logger.debug("No temporary mismatch files found (as expected, 0 mismatches reported).")
                            output_miss_path.touch(exist_ok=True)
                else:
                     logger.error("Temporary worker directory not found during concatenation phase.")
                     success = False


            except Exception as e:
                logger.exception(f"Error during file concatenation: {e}")
                success = False

        # --- Final Validation (only if concatenation succeeded) ---
        if success:
            logger.info("Concatenation complete. Performing final validation...")
            try:
                output_match_count = 0
                output_miss_count = 0

                expected_match_count = total_match_count
                expected_miss_count = total_miss_count

                if output_match_path.is_file() and output_match_path.stat().st_size > 0:
                     logger.debug(f"Counting records in final match file: {output_match_path.name}")
                     # *** Call count_fq_headers without the unsupported header_start_tag argument ***
                     output_match_count = count_fq_headers(output_match_path)
                     logger.debug(f"Found {output_match_count:,} records in {output_match_path.name}")
                elif expected_match_count > 0:
                     logger.error(f"Final match file {output_match_path.name} not found or empty, but expected {expected_match_count:,} records based on worker reports.")
                     success = False

                if output_miss_path.is_file() and output_miss_path.stat().st_size > 0:
                     logger.debug(f"Counting records in final mismatch file: {output_miss_path.name}")
                     # *** Call count_fq_headers without the unsupported header_start_tag argument ***
                     output_miss_count = count_fq_headers(output_miss_path)
                     logger.debug(f"Found {output_miss_count:,} records in {output_miss_path.name}")
                elif expected_miss_count > 0:
                     logger.error(f"Final mismatch file {output_miss_path.name} not found or empty, but expected {expected_miss_count:,} records based on worker reports.")
                     success = False

                # Compare counts only if validation hasn't failed yet
                if success:
                    total_output_count = output_match_count + output_miss_count
                    logger.debug(f"Total records counted in final output files: {total_output_count:,}")
                    logger.debug(f"Total records processed by workers: {total_processed_by_workers:,}") # Use worker total for comparison
                    # Compare final count with worker total
                    if total_output_count == total_processed_by_workers:
                        logger.debug("SUCCESS: Final output record count matches total processed by workers.")
                        success = True # Confirm success
                        # Handle original file deletion *only* on full success
                        if not keep_original:
                            logger.debug(f"Attempting to delete original input file as requested: {input_fastq_orig}")
                            try:
                                input_fastq_orig.unlink()
                                logger.debug(f"Successfully deleted original input file: {format_short_path(input_fastq_orig)}")
                            except OSError as e:
                                logger.error(f"Failed to delete original input file {input_fastq_orig}: {e}")
                    else:
                        logger.error("CRITICAL ERROR: Final output record count mismatch!")
                        logger.error(f"Total Processed by Workers: {total_processed_by_workers:,}, Final Output Count: {total_output_count:,}")
                        logger.error("There might be issues with processing, concatenation or the final count method.")
                        success = False # Mark as failure

            except FileNotFoundError as e:
                 logger.error(f"File not found during final validation: {e}")
                 success = False
            except Exception as e:
                 logger.exception(f"Error during final validation: {e}")
                 success = False


        # --- General Cleanup ---
        # Clean up the main temporary worker directory if it exists
        if workers_dir and workers_dir.exists():
            logger.info(f"Cleaning up temporary worker directory: {workers_dir}")
            shutil.rmtree(workers_dir, ignore_errors=True)

        # Clean up the temporary decompressed file if one was created
        if temp_decompressed_file_path and temp_decompressed_file_path.exists():
            logger.info(f"Cleaning up temporary decompressed file: {temp_decompressed_file_path.name}")
            temp_decompressed_file_path.unlink(missing_ok=True)

        # Final status message
        if success:
            logger.info("FASTQ Match/Edit Tool finished successfully.")
        else:
            logger.error("FASTQ Match/Edit Tool finished with errors.")
            logger.warning("Output files might be incomplete or incorrect due to errors.")

        logging.shutdown()
        return 0 if success else 1


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Filter or replace sequences in a FASTQ file based on a query set using multiprocessing.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    # --- Input Arguments ---
    parser.add_argument(
        "--input-fastq",
        required=True,
        type=Path,
        help="Path to the input FASTQ file (.fastq, .fq, .fastq.gz, .fq.gz, .fastq.zst, .fq.zst)."
    )
    query_group = parser.add_argument_group('Query Sequences (at least one required)')
    query_group.add_argument(
        "--query-seq",
        type=str,
        help="Comma-separated string of query DNA sequences (e.g., 'AGTC,GGCC,TTA'). Case-insensitive."
    )
    query_group.add_argument(
        "--query-seq-file",
        type=Path,
        help="Path to a text file containing query sequences (one per line or tab-separated with sequence in first column). Assumes header if .tsv/.csv."
    )
    # --- Operational Mode Arguments ---
    mode_group = parser.add_mutually_exclusive_group(required=True)
    mode_group.add_argument(
        "--seq-filter",
        action='store_true',
        help="Mode: Filter records. Output matching and non-matching records to separate files."
    )
    mode_group.add_argument(
        "--seq-replace",
        action='store_true',
        help="Mode: Replace sequences. Replace matching sequences with --replacement-seq."
    )
    parser.add_argument(
        "--replacement-seq",
        type=str,
        help="The DNA sequence (ACGTN) to use for replacement (Required if --seq-replace is used). Case-insensitive."
    )
    # --- Output Arguments ---
    parser.add_argument(
        "--output-prefix",
        required=True,
        type=str,
        help="Prefix for output file names (e.g., 'processed_sample')."
    )
    parser.add_argument(
        "--output-dir",
        type=Path, # Use Path directly for type checking
        default=None,
        help="Directory to write output files and logs. Defaults to the input file's directory."
    )
    parser.add_argument(
        "--keep-original",
        action='store_true',
        help="Keep the original input FASTQ file after successful processing (default is to delete it if processing succeeds)."
    )
    parser.add_argument(
        "--log-dir",
        type=Path, # Use Path directly
        default=None, # Default handled in main
        help="Directory to write log file ({output_prefix}.fq_match_edit.log). Defaults to --output-dir if specified, otherwise console only."
    )
    # --- Performance Arguments ---
    parser.add_argument(
        "--processes",
        type=int,
        default=max(1, mp.cpu_count() // 2), # Default to half the CPUs (min 1)
        help="Number of worker processes to use."
    )
    parser.add_argument(
        "--header-start",
        type=str,
        default="@",
        help="Character(s) indicating the start of a FASTQ header line."
    )

    args = parser.parse_args()

    # --- Argument Validation ---
    if not args.query_seq and not args.query_seq_file:
        parser.error("At least one of --query-seq or --query-seq-file must be provided.")
    if args.seq_replace and not args.replacement_seq:
        parser.error("--replacement-seq is required when using --seq-replace.")

    # Validate replacement sequence characters if provided
    validated_replacement_seq = None
    if args.replacement_seq:
        validated_replacement_seq = args.replacement_seq.strip().upper()
        if not validated_replacement_seq:
             parser.error("--replacement-seq cannot be empty.")
        if not all(c in 'ACGTN' for c in validated_replacement_seq):
            parser.error(f"Invalid characters found in --replacement-seq: {args.replacement_seq}. Only A, C, G, T, N allowed.")
        if args.seq_filter:
             # Use logger if available, otherwise print
             (logger or logging).warning("--replacement-seq is ignored when using --seq-filter.")


    if args.seq_filter:
        mode_str = "filter"
    else: # args.seq_replace must be True
        mode_str = "replace"

    # Determine log directory based on args.log_dir and args.output_dir
    actual_log_dir = args.log_dir if args.log_dir else args.output_dir

    # --- Execute Main Function ---
    exit_code = 1 # Default to error
    try:
        # Resolve paths early
        input_path = args.input_fastq.resolve()
        query_file_path = args.query_seq_file.resolve() if args.query_seq_file else None
        output_dir_path = args.output_dir.resolve() if args.output_dir else None
        log_dir_path = actual_log_dir.resolve() if actual_log_dir else None

        exit_code = main(
            input_fastq_orig = input_path,
            output_prefix = args.output_prefix,
            operational_mode = mode_str,
            query_seq_str = args.query_seq,
            query_seq_file = query_file_path,
            replacement_seq = validated_replacement_seq,
            output_dir = output_dir_path,
            log_dir = log_dir_path,
            keep_original = args.keep_original,
            processes = args.processes,
            header_start_tag_str = args.header_start
        )
    except Exception as e:
         # Catch errors during argument processing or main() call setup
         (logger or logging).exception(f"Critical error during script execution setup: {e}")
         exit_code = 1

    sys.exit(exit_code)
