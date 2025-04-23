#!/usr/bin/env python3
import mmap
import argparse
import collections # For orderedDict and defaultDict
from typing import Optional
# File and directory handling
import sys
from pathlib import Path
import shutil
# logging and error handling
import logging
import traceback
import time
import tempfile
# Multiprocessing
import multiprocessing as mp
# Custom utilities
# Find the project root by searching upwards for the 'utils' directory
script_path = Path(__file__).resolve()
project_root = script_path.parent
while not (project_root / 'python_utils').is_dir():
    if project_root == project_root.parent: # Reached filesystem root
        raise FileNotFoundError("Could not find the 'python_utils' directory in parent paths.")
    project_root = project_root.parent
# Add the project root to the Python path
sys.path.append(str(project_root))
from python_utils import (
    setup_logging,
    setup_interrupt_handling,
    define_multiprocess_chunks,
    read_barcode_file,
    compression_utility,
    validate_fastq,
    format_short_path
    )

# =======================================
# Setup logging and interruption handling
# =======================================
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('R1_demultiplex')

_interrupt_event = None

def is_interrupted():
    """Check if interruption has been requested."""
    global _interrupt_event
    return _interrupt_event and _interrupt_event.is_set()


# ============================
# R1 Demultiplexing Subprocess
# ============================
def demultiplex_r1_subprocess(
        r1_file_path: Path, 
        all_workers_output_dir: Path,
        barcode_file: Path, 
        worker_id: int,
        chunk_start: int, 
        chunk_end: int,
        buffer_size: int = 1 * 1024 * 1024,
        max_open_handles = 40) -> None:
    """
    This function should be implemented to read the R1 file, look up IDs in LMDB,
    and write demultiplexed records to output files.
    1. Reads in a header line from the R1 FASTQ file and finds a matching barcode in the LMDB.
    2. Writes the R1 record into a buffer for the corresponding barcode.
    3. If the buffer is full, writes the records to the temporary worker output files corresponding to the barcodes.
    (The worker temp files are later combined into the final output files)

    Parameters
    ----------
    r1_file_path : Path
        Path to the R1 FASTQ file to demultiplex.
    all_workers_output_dir: Path
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
        Write buffer threshold before flushing to output files.
    max_open_handles : int (default = 40)
        Maximum concurrent number of open output file handles.
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
    r1_file_handle = None # Keep track of the file handle to close it in finally
    input_barcode_set = set(barcode_bytes_to_path_map.keys()) # Set of barcodes to filter for
    unassigned_read_count = 0 # Count of reads that do not match any barcode

    try:
        # --- Open R1 File and LMDB ---
        # Create memory map of the R1 input file
        r1_file_handle = open(r1_file_path, 'rb')
        mmapped_file = mmap.mmap(r1_file_handle.fileno(), 0, access=mmap.ACCESS_READ)

        # --- Find Effective Start Position of the first full record in the chunk ---
        mmapped_file.seek(chunk_start)
        effective_start = -1
        current_search_pos = chunk_start
        # Scan forward from chunk_start to find the beginning ('@') of a FASTQ header
        while current_search_pos < chunk_end:
            mmapped_file.seek(current_search_pos)
            line = mmapped_file.readline()
            if not line: break # Reached end of file before finding a header
            if line.startswith(b'@NS'):
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

            # If no valid barcode was found for this record, skip to the next one
            barcode_key = seq_line.strip()
            if not barcode_key in input_barcode_set:
                unassigned_read_count += 1
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
        logger.error(f"An unexpected error occurred during R1 processing: {e}", exc_info=True)
        # Re-raise the exception to ensure the main process knows the worker failed
        raise

    finally:
        # --- Cleanup: Ensure all resources are closed ---
        logger.info(f"Cleaning up resources...")

        # Close memory-mapped file if it was created
        if mmapped_file is not None:
            try:
                mmapped_file.close()
                logger.debug(f"Closed memory mapped file for Worker {worker_id}.")
            except Exception as mmap_close_err:
                 logger.error(f"Error closing mmap file: {mmap_close_err}")

        # Close the original file handle if it was opened
        if r1_file_handle is not None and not r1_file_handle.closed:
             try:
                 r1_file_handle.close()
                 logger.debug(f"Closed R1 file handle.")
             except Exception as fh_close_err:
                  logger.error(f"Error closing R1 file handle: {fh_close_err}")
       
        # Close all file handles remaining in the LRU cache
        logger.info(f"Closing {len(lru_handles_map)} remaining open file handles in LRU cache.")
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
# Main R1 demultiplexing function
# ===============================
def main():
    """Main function for the demultiplexer."""
    parser = argparse.ArgumentParser(description="FASTQ R1 Demultiplexer using memory-mapped files")

    # Required arguments
    parser.add_argument("--r1-file", required=True, help="Path to R1 FASTQ file to demultiplex")
    parser.add_argument("--barcode-file", required=True, help="Path to tab-separated file (.tsv) with barcode sequences in first column")
    parser.add_argument("--dmplex-output-dir", required=True, help="Directory for the demultiplexed R1 output FASTQ files")
    # Optional arguments
    parser.add_argument("--log-dir", help="Directory for log files (default: log to console)")
    parser.add_argument("--processes", type=int, default=4, # Default to number of CPU cores
                        help=f"Number of parallel processes (default: 4)")

    args = parser.parse_args()

    # Convert string paths to Path objects
    r1_file = Path(args.r1_file)
    barcode_file = Path(args.barcode_file)
    dmplex_output_dir = Path(args.dmplex_output_dir)
    log_dir = Path(args.log_dir) if args.log_dir else None
    processes_to_use = args.processes # Renamed for clarity
    temp_decompressed_file_path: Optional[Path] = None

    try:
        #
        # Setup logging and signal handling (only in the main process initially)
        # -----------------------------------------------------------------
        r1_file_basename = r1_file.stem.split('.')[0]
        setup_logging(r1_file_basename, log_dir)
        setup_interrupt_handling() # Setup for main process

        logger.info("===================================================")
        logger.info("Starting R1 Demultiplexer")
        logger.info(f"Python version: {sys.version}")
        logger.info(f"System CPU count: {mp.cpu_count()}")
        logger.info(f"Using {processes_to_use} worker processes")
        logger.info(f"R1 Input File: {format_short_path(r1_file)}")
        logger.info(f"Barcode File: {format_short_path(barcode_file)}")
        logger.info(f"Demultiplexed R1 Output Dir: {format_short_path(dmplex_output_dir)}")
        logger.info(f"Log Directory: {log_dir if log_dir else 'Console'}")
        logger.info("===================================================")

        # Validate inputs
        if not validate_fastq(r1_file):
            logger.error(f"R1 file not found or not a file: {r1_file}")
            return 1
        if not barcode_file.is_file():
            logger.error(f"Barcode file not found or not a file: {barcode_file}")
            return 1

        # Create output directory if it doesn't exist
        dmplex_output_dir.mkdir(parents=True, exist_ok=True)
        # Read in the barcodes from the input file
        input_barcode_info_list = read_barcode_file(barcode_file)
        if not input_barcode_info_list:
            logger.error(f"No valid barcodes found in {barcode_file}.")
            return 1
        
        # --- Decompress the input file if needed ---
        # -------------------------------------------
        temp_decompressed_file_path: Optional[Path] = None # Track temp file for cleanup
        if r1_file.suffix in ['.gz', '.zst']:
            logger.info(f"Input file {r1_file.name} is compressed. Decompressing...")
            start_time = time.time()
            with tempfile.NamedTemporaryFile(
                mode='wb', # Need binary mode for compression utility
                suffix=".fastq", # Keep extension for clarity
                dir=r1_file.parent, # Place temp file near output
                delete=False # Prevent auto-deletion on close, we manage it
            ) as tmp_out:
                temp_decompressed_file_path = Path(tmp_out.name)

            logger.info(f"Decompressing to temporary file: {temp_decompressed_file_path.name}")
            decompression_success = compression_utility(
                input_path=r1_file,
                output_path=temp_decompressed_file_path,
                compress=False
            )

            if not decompression_success \
                or not temp_decompressed_file_path.is_file() \
                    or temp_decompressed_file_path.stat().st_size == 0:
                # Clean up partially created temp file if it exists
                if temp_decompressed_file_path:
                        temp_decompressed_file_path.unlink(missing_ok=True)
                raise RuntimeError(f"Failed to decompress {r1_file} to temporary file.")

            input_fastq_to_process = Path(temp_decompressed_file_path) # Workers will use the decompressed temp file
            elapsed_time = time.time() - start_time
            logger.info(f"Decompression complete ({elapsed_time:.2f}s). Processing temporary file.")
        else:
            logger.info("Input file is not compressed.")
            input_fastq_to_process = Path(r1_file) # Process original directly

        #
        # Demultiplex the R1 files using memory-mapped I/O
        # ------------------------------------------------
        worker_assignment_list = define_multiprocess_chunks(input_fastq_to_process, num_processes=processes_to_use)
        logger.info(f"Distributing R1 file {input_fastq_to_process.name} into {len(worker_assignment_list)} chunks for {processes_to_use} worker processes.")

        # --- Start R1 demultiplex worker processes ---
        processes_list = []    
        for i in range(processes_to_use):
            start_pos, end_pos = worker_assignment_list[i] # Tuple of start and end byte positions in the file
            # Output directory for the temporary worker demultiplexed output
            r1_dmplex_output_dir = Path(dmplex_output_dir) / f"R1_dmplex_temp_output"

            process_kwargs = {
                'r1_file_path': input_fastq_to_process,         # R1 file to demultiplex (Path object) 
                'all_workers_output_dir': r1_dmplex_output_dir, # Output directory for all workers
                'barcode_file': barcode_file,                   # Path to the barcode file
                'worker_id': i,                                 # Worker ID (integer)
                'chunk_start': start_pos,                       # Start byte position of the chunk to process
                'chunk_end': end_pos,                           # End byte position of the chunk to process
                # buffer_size                                   # FASTQ record writing buffer size for commits (in bytes; default = 1 MB)
                # max_open_handles: int = 40
            }

            p = mp.Process(
                target=demultiplex_r1_subprocess,
                kwargs=process_kwargs,
                name=f"R1-Demultiplex-Worker-{i}" # Assign a worker name
            )
            processes_list.append(p)
            p.start()

        # --- Wait for all demultiplexing processes to finish ---
        logger.info("Waiting for R1 demultiplexing processes to complete...")
        interrupted_during_wait = False
        try:
            for p in processes_list:
                p.join() # Wait for the process to terminate
                if is_interrupted(): # Check if main process was interrupted while waiting
                    interrupted_during_wait = True
                    # Don't break immediately, allow logging of exit codes if possible
                if p.exitcode != 0:
                    logger.error(f"Process {p.name} exited with non-zero code: {p.exitcode}")
                else:
                     logger.info(f"Process {p.name} finished successfully.")

        except KeyboardInterrupt: # Catch interrupt specifically during the join loop
            logger.warning("Interruption detected while waiting for workers. Signaling workers to stop.")
            if _interrupt_event:
                 _interrupt_event.set() # Signal workers
            # Assume join timeout handles it.
            interrupted_during_wait = True
        
        if interrupted_during_wait:
            logger.error(f"Demultiplex process was interrupted. Cleaning up temporary files.")
            # Cleanup temporary files and terminate processes
            shutil.rmtree(r1_dmplex_output_dir, ignore_errors=True)
            for p in processes_list:
                if p.is_alive():
                    p.terminate()
            return 1
        
        # --- Delete output files that are empty ---
        r1_dmplex_output_file_paths = list(r1_dmplex_output_dir.glob("**/worker*.fastq"))
        for temp_file in r1_dmplex_output_file_paths:
            if temp_file.is_file() and temp_file.stat().st_size == 0:
                try:
                    temp_file.unlink()
                except OSError as unlink_err:
                    logger.warning(f"Could not delete empty temporary file {temp_file}: {unlink_err}")
        
        # --- Merge the temporary output of the R1 demultiplexing workers into a single file for each barcode ---
        temp_missing_barcode_list = []
        logger.info("Merging temporary output files for each barcode...")
        for barcode, barcode_id in input_barcode_info_list:
            # Search recursively within the temp directory structure
            temp_output_files = list(r1_dmplex_output_dir.glob(f"**/worker_*_{barcode}_{barcode_id}_*.fastq"))
            if not temp_output_files:
                temp_missing_barcode_list.append(barcode)
                continue
            r1_file_basename = r1_file.stem.split('.')[0]
            merged_output_file_path = dmplex_output_dir / f"bcode_{barcode_id}_{barcode}_{r1_file_basename}.dmplx.fastq"

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
                    logger.info(f"Successfully merged {len(temp_output_files)} temporary files ({merged_size} bytes) for barcode {barcode} into {merged_output_file_path.name}")
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
        # Likely due to no matching records
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
        # Cleanup the R1 output directory if it exists and does not contain any files
        if 'r1_dmplex_output_dir' in locals():
            if r1_dmplex_output_dir.exists() and not any(item.is_file() for item in r1_dmplex_output_dir.rglob('*')):
                try:
                    shutil.rmtree(r1_dmplex_output_dir, ignore_errors=True)
                    logger.info(f"Removed temporary output directory: {r1_dmplex_output_dir}")
                except OSError as cleanup_err:
                    logger.error(f"Error cleaning up temporary output directory: {cleanup_err}")
        
        # Clean up the temporary decompressed file if one was created
        if temp_decompressed_file_path and temp_decompressed_file_path.exists():
            logger.info(f"Cleaning up temporary decompressed file: {temp_decompressed_file_path.name}")
            temp_decompressed_file_path.unlink(missing_ok=True)

if __name__ == "__main__":
    # Ensure multiprocessing context is set up early if needed (e.g., for 'spawn' method)
    # mp.set_start_method('spawn') # Uncomment if needed, e.g., on macOS/Windows or for specific libraries
    exit_code = main()
    if exit_code == 0:
        logger.info("R1 Demultiplexer finished successfully.")
    else:
        logger.error(f"R1 Demultiplexer exited with error code {exit_code}.")
    sys.exit(exit_code)