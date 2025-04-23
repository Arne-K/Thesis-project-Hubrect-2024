import os
from pathlib import Path
import sys
from typing import List, Iterator, Tuple, Optional
import subprocess
import shutil
import errno
import tempfile
# Multiprocessing
import multiprocessing as mp
from multiprocessing import Manager
# For memory-mapping
import mmap

# -- import utilities --
from .logging_utils import format_short_path
from .file_management_utils import retain_decompressed_duplicate_paths, compression_utility
from .command_utils import execute_command

def extract_subseq_from_fastq(input_fastq: Path,
                              output_fastq: Path,
                              subseq_positions: Tuple[int, int]) -> bool:
    """
    The output is a new FASTQ file where the records contain a positional subsequence.
    This function can used to extract the barcode sequence from the R1 file (position 1-8).
    
    Parameters
    ----------
    input_fastq : Path
        Path to the input FASTQ file.
    output_dir : Path
        Path to the output directory where the file will be created.
        The output file will be named as '{input_file_stem}.bcode.fastq'.
    subseq_positions: Tuple
        Tuple containing the start and end positions of the subsequence to extract.
        The counting starts at 1.
    """
    # --- Check input and initialize output ---
    # -----------------------------------------
    # check directories
    if not isinstance(input_fastq, Path):
        print(f"input file has to be a Path object")
        return False
    if not isinstance(output_fastq, Path):
        print(f"output directory has to be a Path object")
        return False
    
    # check subseq positions input
    if not isinstance(subseq_positions, tuple):
        print(f"subseq_positions has to be a tuple")
        return False
    elif len(subseq_positions) != 2:
        print(f"subseq_positions has to be a tuple of length 2")
        return False
    elif not isinstance(subseq_positions[0], int) or not isinstance(subseq_positions[1], int):
        print(f"subseq_positions has to be a tuple of integers")
        return False
    elif subseq_positions[0] > subseq_positions[1]:
        print(f"subseq_positions[0] has to be smaller than subseq_positions[1]")
        return False
    
    # check fastq input
    is_valid_fastq = validate_fastq(input_fastq)
    if not is_valid_fastq:
        print(f"Error: Input file '{input_fastq}' is not a valid FASTQ file.")
        return False
    
    
    # initialize output
    if output_fastq.is_file():
        output_fastq.unlink()

    # --- Decompress the input file if needed ---
    if input_fastq.suffix.lower() in ['.gz', '.gzip', '.zst']:
        decompressed_fastq = input_fastq.with_suffix('')
        if compression_utility(input_fastq, decompressed_fastq, compress=False):
            input_fastq = decompressed_fastq
        else:
            return False
    
    # --- Extract the subsequence ---
    # -------------------------------
    subseq_start, subseq_end = subseq_positions
    seqkit_cmd = [
        "seqkit", "subseq",
        "-r", f"{subseq_start}:{subseq_end}"
        ]
    try:
        # Open the output file in binary write mode ("wb")
        with open(output_fastq, "wb") as f_out:
            # Open the input file in binary read mode ("rb")
            with open(input_fastq, "rb") as f_in:
                print(f"Extracting sequence positions {subseq_start}:{subseq_end} from {format_short_path(input_fastq)}")
                print(f"Writing to {format_short_path(output_fastq)}")
                # Start the seqkit process
                # - stdin=subprocess.PIPE: We provide input via pipe
                # - stdout=f_out: Redirect seqkit's standard output DIRECTLY to our output file handle
                # - stderr=None: Let seqkit's errors print to the console (or redirect to a file/PIPE)
                process = subprocess.Popen(
                    seqkit_cmd,
                    stdin=subprocess.PIPE,
                    stdout=f_out,
                    stderr=None
                )

                # Stream the input data from f_in to the process's stdin
                try:
                    shutil.copyfileobj(f_in, process.stdin)
                except BrokenPipeError:
                    # This might happen if seqkit finishes or errors before reading all input
                    print(f"Warning: Pipe to seqkit's stdin broke. Seqkit might have exited before reading all input.")
                except Exception as e:
                    print(f"Error during input streaming to seqkit: {e}")
                    # Try to stop the process if streaming failed
                    try:
                        process.terminate()
                    except ProcessLookupError:
                        pass # Process already finished
                    raise # Re-raise the exception that caused the streaming failure

                finally:
                    # Close the process's stdin to signal EOF (End Of File)
                    # This tells seqkit that no more input is coming.
                    if process.stdin:
                        try:
                            process.stdin.close()
                        except OSError as e:
                            # Ignore EBADF (Bad file descriptor) error, which means it's likely already closed
                            if e.errno != errno.EBADF:
                                print(f"Error closing seqkit stdin: {e}")
                                # Don't raise here usually, allow process.wait() to report exit code
                        except Exception as e:
                            print(f"Non-OS error closing process stdin: {e}")

            # After the input file is processed and stdin is closed,
            # wait for the seqkit process to complete its work and exit.
            # The output is written directly to f_out during this time.
            return_code = process.wait()

            # Check if seqkit reported an error exit code
            if return_code != 0:
                print(f"Error: seqkit command failed with exit code {return_code}")
                return False
            # Check if the output file is empty
            elif not output_fastq.stat().st_size > 0:
                print(f"Error: Output file '{output_fastq}' is empty.")
                return False 
            else:
                print(f"Sub-sequence extraction completed successfully.")
                return True
    except FileNotFoundError:
        print(f"Error: Input file '{input_fastq}' not found or output path invalid.")
        return False
    except Exception as e:
        print(f"An unexpected script error occurred: {e}")
        return False


def concatenate_fastq_files(
        input_dir: Path,
        target_file_pattern: str,
        sample_prefix: str):
    """
    Concatenate FASTQ files in a directory.
    
    Parameters
    ----------
    input_dir : Path
        Path to the directory containing the FASTQ files to be concatenated.
    target_file_pattern : str
        Pattern to match the target file names. 
        The pattern will be searched as a substring flanked by two underscores in the file name:
        '*_{target_file_pattern}_*.fastq*'
        '*_{target_file_pattern}_*.fq*'
        '*_{target_file_pattern}_*.fastq.gz*'
        '*_{target_file_pattern}_*.fq.gz*'
    sample_prefix : str
        Prefix of the sample name. The concatenated FASTQ file will be named as '{sample_prefix}_{target_file_pattern}.fastq.gz'.

    Returns
    -------
    bool
        True if the concatenation was successful, False otherwise.
    """
    # --- Check input and initialize output ---
    if not input_dir.is_dir():
        print(f"Error: Input directory '{input_dir}' is not a directory.")
        return False
    target_file_list = None
    seqkit_paths_file = None
    
    try:
        # --- Search the input directory for FASTQ files ---
        target_file_list = input_dir.glob(f"*{target_file_pattern}*")
        if not target_file_list:
            print(f"Error: No files found in '{input_dir}' with target file pattern '{target_file_pattern}'.")
            return False
        target_file_list = [path for path in target_file_list if validate_fastq(path)]
        if not target_file_list:
            print(f"Error: No valid FASTQ files found in '{input_dir}' with target file pattern '{target_file_pattern}'.")
            return False
        print(f"Found {len(target_file_list)} FASTQ files in '{input_dir}' with target file pattern '{target_file_pattern}'.")

        # --- Decompress the input files if needed ---
        # Make sure the file list contains only decompressed file duplicates, if possible
        original_target_file_list = target_file_list.copy()
        target_file_list = retain_decompressed_duplicate_paths(target_file_list)
        decompressed_target_file_list = []
        for file in target_file_list:
            if file.suffix in ['.gz','.zst']:
                decompressed_file = file.with_suffix('')
                decompression_success = compression_utility(
                    input_path=file,
                    output_path=decompressed_file,
                    compress=False
                )
                # Check for failure
                if not decompression_success \
                    or not decompressed_file.is_file() \
                        or not decompressed_file.stat().st_size > 0:
                    print(f"Failed to decompress {file} to {decompressed_file}")
                    return False
                # Update the decompressed file list
                decompressed_target_file_list.append(decompressed_file)
            else:
                decompressed_target_file_list.append(file)
        if len(decompressed_target_file_list) == len(target_file_list):
            target_file_list = decompressed_target_file_list
        else:
            print(f"Error: Failed to decompress the input files.")
            return False

        # --- Create the path list file for seqkit '--infile-list' parameter ---
        seqkit_paths_file = Path(input_dir) / "seqkit_paths.txt"
        if seqkit_paths_file.is_file():
            seqkit_paths_file.unlink()
        with open(seqkit_paths_file, "w") as f:
            f.writelines('\n'.join(str(path) for path in target_file_list))
        if not seqkit_paths_file.is_file() or not seqkit_paths_file.stat().st_size > 0:
            print(f"Error: Failed to create temporary input file '{seqkit_paths_file}'.")
            return False
        
        # --- Concatenate the files using seqkit ---
        concat_fastq_out = Path(input_dir) / f"{target_file_pattern}_{sample_prefix}.fastq"
        if concat_fastq_out.is_file():
            concat_fastq_out.unlink()
        seqkit_cmd = [
            "seqkit", "seq", 
            "--threads", "6", 
            "--infile-list", str(seqkit_paths_file), 
            "-o", str(concat_fastq_out)
            ]
        print(f"Running seqkit command: {' '.join(seqkit_cmd)}")
        exit_status = execute_command(seqkit_cmd)
        # Check for failure
        if exit_status != 0:
            if concat_fastq_out.is_file():
                concat_fastq_out.unlink()
            print(f"Error: seqkit command failed with exit code {exit_status}")
            return False
        if not concat_fastq_out.is_file() or not concat_fastq_out.stat().st_size > 0:
            if concat_fastq_out.is_file():
                concat_fastq_out.unlink()
            print(f"Error: Failed to concatenate files in '{input_dir}'.")
            return False
        print(f"Successfully concatenated files in '{input_dir.name}'")
        # Delete the parent files that are no longer needed
        files_to_delete = []
        for path in original_target_file_list:
            if path != concat_fastq_out:
                files_to_delete.append(path)
        for path in target_file_list:
            if path != concat_fastq_out:
                files_to_delete.append(path)
        for f in files_to_delete:
            if f.is_file():
                f.unlink()
        return True
    finally:
        # Clean up the temporary file list
        if seqkit_paths_file and seqkit_paths_file.is_file():
            seqkit_paths_file.unlink()

def validate_fastq(
        input_file: Path,
        valid_file_extensions: List[str] = [".fastq", ".fq", ".fastq.gz", ".fq.gz", ".fastq.zst", ".fq.zst"]):
    """
    Validate the input FASTQ file.
    
    Parameters
    ----------
    input_file : Path
        Path to the file to be validated.
    valid_file_extensions : List[str], optional
        List of valid file extensions for the input FASTQ file. Default is [".fastq", ".fq", ".fastq.gz", ".fq.gz", ".fastq.zst", ".fq.zst"].
    
    Returns
    -------
    bool
        True if the input file is a valid FASTQ file, False otherwise.
    """
    if not input_file.exists() or not input_file.is_file():
        print(f"Error: Input file '{input_file}' not found.", file=sys.stderr)
        return False
    elif not input_file.stat().st_size > 0:
        print(f"Error: Input file '{input_file}' is empty.", file=sys.stderr)
        return False
    for ext in valid_file_extensions:
        if input_file.name.endswith(ext):
            return True
    print(f"Error: Input file '{input_file}' does not have a valid extension. Valid extensions are: {valid_file_extensions}", file=sys.stderr)
    return False

# ===================
# Count FASTQ Headers
# ===================
def count_fq_headers_subprocess(fq_inpath:Path, 
                                process_name:str, 
                                byte_range:tuple, 
                                return_dict,
                                header_start_tag=b'@N'):
    start_pos_byte = byte_range[0]
    end_pos_byte = byte_range[1]
    
    print(f"Worker {process_name} checking bytes {start_pos_byte} to {end_pos_byte}")
    header_count = 0
    mmapped_file = None
    input_file_handle = None
    
    try:
        # Define buffer size for I/O operations (16 MB)
        BUFFER_SIZE = 16 * 1024 * 1024  # 16 MB in bytes
        
        # --- Open file as memory-mapped with large buffer ---
        input_file_handle = open(fq_inpath, 'rb', buffering=BUFFER_SIZE)
        mmapped_file = mmap.mmap(input_file_handle.fileno(), 0, access=mmap.ACCESS_READ)

        # --- Position the worker to its chunk start ---
        # ----------------------------------------------
        mmapped_file.seek(start_pos_byte)
        effective_start = -1
        current_search_pos = start_pos_byte
        # Scan forward from chunk_start to find the beginning ('@') of a FASTQ header
        while current_search_pos < end_pos_byte:
            mmapped_file.seek(current_search_pos)
            line = mmapped_file.readline()
            if not line: break # Reached end of file before finding a header
            if line.startswith(header_start_tag):
                effective_start = current_search_pos # Found the start of a record
                break
            # If not a header, advance position to the start of the next line
            current_search_pos = mmapped_file.tell()
        # If no header is found within the chunk, there's nothing to process
        if effective_start == -1:
            return_dict[process_name] = 0
            return # Exit cleanly
        # Position the memory map at the start of the first valid record
        mmapped_file.seek(effective_start)

        # --- Main Processing Loop ---
        # ----------------------------
        while mmapped_file.tell() < end_pos_byte:
            line = mmapped_file.readline()
            if line.startswith(header_start_tag):
                header_count += 1 
                mmapped_file.readline() # Sequence
                mmapped_file.readline() # +
                mmapped_file.readline() # Quality
    finally:
        # Close memory-mapped file if it was created
        if mmapped_file is not None:
            try:
                mmapped_file.close()
            except Exception as mmap_close_err:
                print(f"Error closing mmap file: {mmap_close_err}")
        # Close the original R2 file handle if it was opened
        if input_file_handle is not None and not input_file_handle.closed:
             try:
                 input_file_handle.close()
             except Exception as fh_close_err:
                  print(f"Error closing file handle: {fh_close_err}")
        
        # Store the result in the shared dictionary instead of updating a shared counter
        return_dict[process_name] = header_count


def count_fq_headers(input_fastq: Path, processes: int=4):
    
    # --- Decompress the file if needed ---
    temp_decompressed_file_path: Optional[Path] = None
    total_headers = 0
    try:
        # --- Decompress Input File if Needed (using temp file) ---
        if input_fastq.suffix in ['.gz', '.zst']:
            # Create a temporary file in the output directory (or system temp)
            with tempfile.NamedTemporaryFile(
                mode='wb', # Need binary mode for compression utility
                suffix=".fastq", # Keep extension for clarity
                dir=input_fastq.parent, # Place temp file near output
                delete=False # Prevent auto-deletion on close, we manage it
            ) as tmp_out:
                temp_decompressed_file_path = Path(tmp_out.name)

            decompression_success = compression_utility(
                input_path=input_fastq,
                output_path=temp_decompressed_file_path,
                compress=False
            )

            if not decompression_success or not temp_decompressed_file_path.is_file() or temp_decompressed_file_path.stat().st_size == 0:
                # Clean up partially created temp file if it exists
                if temp_decompressed_file_path:
                     temp_decompressed_file_path.unlink(missing_ok=True)
                raise RuntimeError(f"Failed to decompress {input_fastq.name} to temporary file.")

            input_fastq_to_process = temp_decompressed_file_path # Workers will use the decompressed temp file
        else:
             input_fastq_to_process = input_fastq # Process original directly

    
        # Define chunks for multiprocessing
        chunks = define_multiprocess_chunks(file_path=input_fastq_to_process, num_processes=processes)
        
        # Use a Manager.dict() instead of Manager.Value() for less contention
        with Manager() as manager:
            return_dict = manager.dict()
            processes_list = []
            
            for i, chunk in enumerate(chunks):
                process_name = f"temp_{i+1:03d}"
                # Create and start a process for each chunk
                p = mp.Process(
                    target=count_fq_headers_subprocess,
                    args=(
                        input_fastq_to_process,
                        process_name,
                        chunk,
                        return_dict,
                    )
                )
                p.start()
                processes_list.append(p)
            
            # Wait for all processes to complete
            for p in processes_list:
                p.join()
            if return_dict:
                total_headers = sum(return_dict.values())
            else:
                total_headers = 0
    finally:
            if temp_decompressed_file_path and temp_decompressed_file_path.exists():
                  temp_decompressed_file_path.unlink()
            return total_headers

# =========================
# Multiprocessing utilities
# =========================
def define_multiprocess_chunks(file_path: Path, num_processes: int):
    """
    Divides a FASTQ file into byte ranges for multiprocessing, ensuring
    each chunk starts at the beginning of a record (@ line) and ends
    at the end of a record (quality line).
    """
    if num_processes <= 0:
        raise ValueError("Number of processes must be positive.")

    try:
        file_size = os.path.getsize(file_path)
    except FileNotFoundError:
        raise FileNotFoundError(f"File not found: {file_path}")

    if file_size == 0:
        return [] # No chunks for an empty file

    print(f"File size: {file_size} bytes")

    # Assign each worker an approximate proportion of the file size
    # Handle potential division by zero if num_processes was not validated
    # Although validation is added above, this check makes it more robust
    if num_processes == 0: num_processes = 1 # Avoid division by zero, process sequentially
    chunk_size = file_size // num_processes
    # Ensure minimum chunk size if file is very small or num_processes is large
    if chunk_size == 0: chunk_size = 1

    chunk_list = []
    current_pos = 0

    for i in range(num_processes):
        # Find the start of the next complete record boundary
        # We search from current_pos
        chunk_start = find_chunk_bound(
            file_path=file_path,
            byte_pos=current_pos,
            find_start_boundary=True # We want the start byte of the '@' line
        )

        # If find_chunk_bound returns None (EOF reached before finding header), stop
        if chunk_start is None:
             # Add the last chunk if there was a valid start before EOF
            if chunk_list and chunk_list[-1][1] < file_size and current_pos < file_size:
                 # Find the true end of the previous record if possible
                last_record_end = find_chunk_bound(
                    file_path=file_path,
                    byte_pos=chunk_list[-1][0], # Start searching from last known record start
                    find_start_boundary=False # Find the end of this record
                )
                if last_record_end is not None and last_record_end <= file_size:
                     chunk_list[-1] = (chunk_list[-1][0], last_record_end)
                else: # Couldn't find proper end, just go to file end
                     chunk_list[-1] = (chunk_list[-1][0], file_size)

            print(f"EOF reached while searching for start of chunk {i+1} near byte {current_pos}. Stopping.")
            break # No more chunks can start

        if i == num_processes - 1:
            # Last chunk: starts at chunk_start and goes to the end of the file
            chunk_end = file_size
        else:
            # Calculate a target end position based on ideal chunk size
            target_end_pos = chunk_start + chunk_size
            # Ensure target is not beyond file size for search efficiency
            if target_end_pos >= file_size:
                 target_end_pos = file_size -1 # Seek somewhere before the end

            # Find the end boundary of a record near the target position
            # We search from target_end_pos
            chunk_end = find_chunk_bound(
                file_path=file_path,
                byte_pos=target_end_pos,
                find_start_boundary=False # We want the position AFTER the quality line
            )

            # If EOF was reached before finding a suitable end boundary near target,
            # just end this chunk at the end of the file. The loop will break next.
            if chunk_end is None or chunk_end <= chunk_start :
                 chunk_end = file_size


        # Only add chunk if it's valid (start is before end)
        if chunk_start < chunk_end:
             #print(f"Defined chunk {i+1}: Start={chunk_start}, End={chunk_end}")
             chunk_list.append((chunk_start, chunk_end))
             # Start position for the next chunk search IS the end of the current one
             current_pos = chunk_end
        elif chunk_start == file_size: # If the only place to start is EOF, we are done
             print(f"Chunk start search for chunk {i+1} resulted in EOF. Stopping.")
             break
        else:
             # This might happen if find_chunk_bound seeking end wraps around or finds
             # the same record start was on. End chunk here and stop.
              print(f"Chunk start ({chunk_start}) >= chunk end ({chunk_end}). Ending chunk {i+1} at EOF and stopping.")
              chunk_list.append((chunk_start, file_size))
              break


        # Prevent infinite loop if current_pos isn't advancing
        if current_pos <= chunk_start and i < num_processes -1:
             print(f"Warning: current_pos ({current_pos}) did not advance past chunk_start ({chunk_start}). Forcing advance.")
             # Try to find the end of the current record to advance past it
             actual_end = find_chunk_bound(file_path, chunk_start, find_start_boundary=False)
             if actual_end is not None and actual_end > chunk_start:
                  current_pos = actual_end
             else:
                  # Cannot determine end, break to avoid infinite loop
                  print("Error: Cannot advance position. Stopping chunk definition.")
                  # Optionally adjust the last defined chunk to EOF
                  if chunk_list: chunk_list[-1] = (chunk_list[-1][0], file_size)
                  break


        # Stop if we've processed the whole file
        if current_pos >= file_size:
             break
        
    # If the last defined chunk doesn't reach the end of the file, adjust it.
    # This handles cases where the loop exited early but didn't set the last chunk to file_size
    if chunk_list and chunk_list[-1][1] < file_size:
         last_start = chunk_list[-1][0]
         # Try to find the real end boundary if the last chunk wasn't the designated 'last chunk'
         final_end = find_chunk_bound(file_path, last_start, find_start_boundary=False)
         if final_end is None or final_end < last_start:
             final_end = file_size # Default to EOF if finding end fails

         chunk_list[-1] = (last_start, min(final_end, file_size)) # Ensure not past EOF
         
    return chunk_list


def find_chunk_bound(
    file_path: Path,
    byte_pos: int, # Byte position should be an integer offset because that's how we seek
    find_start_boundary: bool,
    header_start_char: bytes = b'@' # More general FASTQ check
) -> int | None: # Return int offset or None if EOF reached before finding
    """
    Finds the byte boundary for a FASTQ record chunk.

    Args:
        file_path: Path to the FASTQ file.
        byte_pos: The approximate byte position to start searching from.
        find_start_boundary: If True, find the start of the header ('@') line.
                             If False, find the end of the quality line (start of next record).
        header_start_char: The character indicating the start of a header line.

    Returns:
        The exact byte offset of the boundary, or None if EOF is reached
        before a suitable boundary is found.
    """
    try:
        # Define buffer size for I/O operations (16 MB)
        BUFFER_SIZE = 16 * 1024 * 1024  # 16 MB in bytes
        
        with open(file_path, 'rb', buffering=BUFFER_SIZE) as file:
            file.seek(byte_pos)

            # If searching for end boundary, read the partial line we landed in
            # to ensure we start search on the *next* line. Doesn't apply
            # when searching for start, as we might land on the '@'.
            if not find_start_boundary and byte_pos != 0:
                 file.readline()


            while True:
                # Record position *before* reading the potential header line
                start_of_line_pos = file.tell()
                line = file.readline()

                # Check for EOF
                if not line:
                    # If looking for start, EOF means no more records.
                    # If looking for end, EOF is a valid end boundary.
                    return None if find_start_boundary else start_of_line_pos

                # Check if this is a header line
                if line.startswith(header_start_char):
                    if find_start_boundary:
                        # Found the start boundary, return the position *of* the '@'
                        return start_of_line_pos
                    else:
                        # Found a header, now read the remaining 3 lines of THIS record
                        # to find the position AFTER the quality line.
                        try:
                            # Read Sequence line
                            if not file.readline(): return file.tell() # EOF is valid end
                            # Read '+' line
                            if not file.readline(): return file.tell() # EOF is valid end
                            # Read Quality line
                            if not file.readline(): return file.tell() # EOF is valid end

                            # After reading the quality line, the current position
                            # is the start of the NEXT record or EOF. This is our END boundary.
                            return file.tell()
                        except Exception as e:
                             # Handle potential errors during reads if needed
                             print(f"Error reading record lines after header at {start_of_line_pos}: {e}")
                             return file.tell() # Return current position on error


    except FileNotFoundError:
         print(f"Error: File not found in find_chunk_bound: {file_path}")
         raise # Re-raise the exception
    except Exception as e:
         print(f"An unexpected error occurred in find_chunk_bound: {e}")
         return None # Indicate failure
