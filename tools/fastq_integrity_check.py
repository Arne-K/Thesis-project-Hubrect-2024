#!/usr/bin/env python3
import re
from pathlib import Path
import argparse
# Multiprocessing
import multiprocessing as mp
from multiprocessing import Manager
from collections import Counter
import sys
import shutil
import mmap
import logging

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
        print(f"Added '{path_str}' to sys.path")


from utils.python_utils import (
    setup_logging,
    setup_interrupt_handling,
    # Multiprocessing
    define_multiprocess_chunks
    )

# =======================================
# Setup logging and interruption handling
# =======================================
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger('fastq_integrity_check')

_interrupt_event = None

def is_interrupted():
    """Check if interruption has been requested."""
    global _interrupt_event
    return _interrupt_event and _interrupt_event.is_set()

# 
# ========================== Functions ==========================
#
def record_fmt_check(
        record_lines: list[bytes],
        header_line_start_tag=b'@NS', 
        plus_character=b'+', 
        valid_seq_characters=re.compile(b"^[ACGTN]+$")):
    """
    check the correctness of a single 4-line FASTQ record
    Parameters
    ----------
    record_lines : list[bytes]
        The lines of the FASTQ record to check.
    header_line_start_tag : bytes
        The expected starting tag of the header line.
    plus_character : bytes
        The expected character for the plus line.
    valid_seq_characters : re.Pattern
        A compiled regex pattern for valid sequence characters. 
    Returns
    -------
    bool
        True if the record is valid, False otherwise.
    """
    
    # Check if we have the correct number of lines
    if len(record_lines) != 4:
        return False, f"A valid record has 4 lines, received {len(record_lines)}"
    
    # Check if the header line starts with the expected tag
    if not record_lines[0].startswith(header_line_start_tag):
        return False, f"Header line does not start with {header_line_start_tag.decode()}"
    # 3rd line has to only contain +
    if record_lines[2].strip() != plus_character:
        return False, "Incorrect plus line"
    
    # Sequence line (2nd) has to contain valid characters
    seq_line = record_lines[1].strip()
    if not valid_seq_characters.match(seq_line):
        return False, "Sequence contains invalid characters"
    
    # Sequence and quality lines have to be of equal length. 
    quality_line = record_lines[3].strip()
    if len(seq_line) != len(quality_line):
        return False, "Sequence and quality line lengths not equal"
    
    # If we get here, the record is valid
    return True, None

def integrity_subprocess(
        fq_inpath: Path, 
        malformed_outdir: Path, 
        tmp_dirname: str, 
        byte_range: tuple[int, int], 
        store_malformed: bool = True,
        header_start_tag=b'@NS',
        shared_counter=None, 
        shared_error_list=None):
    """
    A subprocess for integrity checking
    Parameters
    ----------
    fq_inpath : Path
        Path to the input FASTQ file to be checked.
    malformed_outdir : Path
        Directory to store the malformed records output file.
    tmp_dirname : str
        Temporary directory name for the subprocess.
    byte_range : tuple
        Byte range to process in the FASTQ file.
    store_malformed : bool
        Whether to store malformed records in the output directory.
    header_start_tag : bytes
        The starting tag of the header line in the FASTQ file.
    shared_counter : multiprocessing.Value
        Shared counter for tracking malformed records across processes.
    shared_error_list : multiprocessing.Manager().list
        Shared list for storing error messages across processes.
    """
    
    chunk_start = byte_range[0]
    chunk_end = byte_range[1]
    logger.info(f"Processing bytes {chunk_start} to {chunk_end}")
    
    malformed_record_count = 0
    malformed_outfile = None
    if store_malformed:
        fq_base = fq_inpath.stem
        tmp_dir = malformed_outdir / tmp_dirname
        tmp_dir.mkdir(exist_ok=True, parents=True)
        malformed_record_output_file = tmp_dir / f"{fq_base}_malformed_temp.fastq"
        if malformed_record_output_file.is_file():
            malformed_record_output_file.unlink()
        malformed_record_output_file.touch()
        malformed_outfile = open(malformed_record_output_file, 'wb')
    
    try:
        # Memory map the input FASTQ file for efficient reading
        input_file_handle = open(fq_inpath, 'rb')
        mmapped_file = mmap.mmap(input_file_handle.fileno(), 0, access=mmap.ACCESS_READ)
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
        
        # --- Main Processing Loop ---
        # Continue reading records as long as the current position is within the chunk bounds
        while mmapped_file.tell() < chunk_end:
            current_record_start_pos = mmapped_file.tell()
            # Read the 4 lines of a FASTQ record
            header_line = mmapped_file.readline()
            # Check for EOF immediately after reading header
            if not header_line: break
            elif not header_line.startswith(header_start_tag):
                # If the header line doesn't start with the expected tag, this is a potentially malformed record
                logger.warning(f"Malformed header line near pos {current_record_start_pos}: {header_line[:50]!r}")
            
            # Read the next 3 lines to complete the record
            seq_line = mmapped_file.readline()
            plus_line = mmapped_file.readline()
            qual_line = mmapped_file.readline()
            # Check if a complete record was read
            if not seq_line or not plus_line or not qual_line:
                # This indicates a potentially corrupted file or reaching EOF mid-record
                logger.warning(f"Incomplete FASTQ record near pos {current_record_start_pos}. Header: {header_line[:50]!r}")
                malformed_outfile.write(header_line + b'\n')
            else:
                record_lines = []
                record_lines = [header_line, seq_line, plus_line, qual_line]
                is_valid_record, error_message = record_fmt_check(record_lines=record_lines)
                        
                if not is_valid_record:
                    if malformed_outfile:
                        for line in record_lines:
                            malformed_outfile.write(line + b'\n')
                    # Add error to shared list if it exists
                    if shared_error_list is not None:
                        shared_error_list.append(error_message)
                    malformed_record_count += 1
                    # Update the shared counter if it exists
                    if shared_counter is not None:
                        shared_counter.value += 1               
    finally:
        if input_file_handle:
            input_file_handle.close()
        if mmapped_file:
            mmapped_file.close()
        logger.info(f"Finished processing. Malformed records found: {malformed_record_count}")
        if malformed_outfile:
            malformed_outfile.close()
            # Delete the temporary file if no malformed records were found
            if malformed_record_count == 0:
                malformed_record_output_file.unlink()
                # Clean up the temporary directory if empty
                if not any(tmp_dir.iterdir()):
                    shutil.rmtree(tmp_dir)
            else:
                logger.info(f"Malformed records written to {malformed_record_output_file}")

def main(
        input_fastq: Path, 
        malformed_output_dir: Path = Path("./integrity_check_malformed"), 
        store_malformed:bool = True,
        processes:int = 4,
        log_dir:Path = None):
    """
    Main function to check FASTQ file integrity using multiprocessing.
    Parameters
    ----------
    fq_inpath : Path
        Path to the input FASTQ file to be checked.
    malformed_output_dir : Path
        Directory to store the malformed records output file.
    store_malformed : bool
        Whether to store malformed records in the output directory.
    processes : int
        Number of processes to use for multiprocessing.
    log_dir : Path
        Directory to store log files. If None, logs will be printed to console.
    """
    #
    # Setup logging and signal handling (only in the main process initially)
    # -----------------------------------------------------------------
    setup_logging(f"integrity_check_{input_fastq.stem}", log_dir)
    setup_interrupt_handling() # Setup for main process

    logger.info("===================================================")
    logger.info("Starting FASTQ integrity check")
    logger.info(f"Python version: {sys.version}")
    logger.info(f"System CPU count: {mp.cpu_count()}")
    logger.info(f"Using {processes} worker processes")
    logger.info(f"Input File: {input_fastq}")
    logger.info(f"Log Directory: {log_dir if log_dir else 'Console'}")
    if store_malformed:
        logger.info(f"Malformed records will be stored in {malformed_output_dir}")
    logger.info("===================================================")

    # Check if the input file exists
    if not input_fastq.is_file():
        logger.error(f"Input file {input_fastq} does not exist.")
        raise FileNotFoundError(f"Input file {input_fastq} does not exist.")
    
    try:
        # Create the output directory if it doesn't exist
        malformed_output_dir = Path(malformed_output_dir)
        malformed_output_dir.mkdir(parents=True, exist_ok=True)
        
        # Define multiprocessing chunks
        mp_chunks = define_multiprocess_chunks(file_path=input_fastq, num_processes=processes)
        
        # Create a shared counter for tracking malformed records across all processes
        with Manager() as manager:
            total_malformed_counter = manager.Value('i', 0)
            shared_error_list = manager.list()  # Create a shared list for error messages
            
            # Create directories for each subprocess
            processes = []
            for i, chunk in enumerate(mp_chunks):
                tmp_dirname = f"temp_{i+1:03d}"
                tmp_dir_path = malformed_output_dir / tmp_dirname
                if tmp_dir_path.is_dir():
                    shutil.rmtree(tmp_dir_path)
                tmp_dir_path.mkdir(exist_ok=True)
                
                # Create and start a process for each chunk
                p = mp.Process(
                    target=integrity_subprocess,
                    kwargs={
                        'fq_inpath': Path(input_fastq),
                        'malformed_outdir': malformed_output_dir,
                        'tmp_dirname': tmp_dirname,
                        'byte_range': chunk,
                        'store_malformed': store_malformed,
                        'shared_counter': total_malformed_counter,
                        'shared_error_list': shared_error_list  # Pass the shared list
                        # 'header_start_tag': b'@NS' # This is the default value
                    }
                )
                p.start()
                processes.append(p)
            
            # Wait for all processes to complete
            for p in processes:
                p.join()
            
            print(">>> All integrity check processes completed successfully\n")
            
            # Display the total count of malformed records
            if total_malformed_counter.value > 0:
                print(f"\n===== SUMMARY =====")
                print(f"Total malformed records found: {total_malformed_counter.value}")
                error_counter = Counter(shared_error_list)  # Use the shared list of errors
                for key, value in error_counter.items():
                    print(f"{key}: {value}")
                print(f"==================\n")
            # Store the total malformed records count in a variable
            total_malformed_records = total_malformed_counter.value
        
        # get malformed file paths in tmp directories 
        malformed_temp_output_file_list = [file for file in malformed_output_dir.rglob(f"{input_fastq.stem}*malformed*temp*.fastq") 
                        if file.parent.name.startswith("temp_")]
        # Merge the malformed output files
        if malformed_temp_output_file_list:
            malformed_output_file = malformed_output_dir / f"{input_fastq.stem}_malformed.fastq"
            if malformed_output_file.is_file():
                malformed_output_file.unlink()
            malformed_output_file.touch()
            # open the merged file for writing
            with open(malformed_output_file, 'wb') as merged_out_f:
                for temp_file in malformed_temp_output_file_list:
                    with open(temp_file, 'rb') as temp_f:
                        shutil.copyfileobj(temp_f, merged_out_f)
            # If the merge was successful, clean up the individual temp files
            if not malformed_output_file.is_file() or not malformed_output_file.stat().st_size > 0:
                print(f"Failed to merge malformed records into {malformed_output_file}")
                if malformed_output_file.is_file():
                    malformed_output_file.unlink()
                raise RuntimeError(f"Failed to merge malformed records into {malformed_output_file}")
            else:
                # Consider it a successful merge
                print(f"Malformed records are found in {malformed_output_file}")
                for temp_file in malformed_temp_output_file_list:
                    temp_file.unlink()
    finally:
        # Clean up the temporary directories
        if malformed_output_dir.is_dir():
            for tmp_dir in malformed_output_dir.iterdir():
                if tmp_dir.is_dir() and tmp_dir.name.startswith("temp_"):
                    shutil.rmtree(tmp_dir)
            # Remove the main output directory if empty
            if not any(malformed_output_dir.iterdir()):
                malformed_output_dir.rmdir()
        if total_malformed_records:
            return total_malformed_records
        else:
            return None


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Check FASTQ file integrity using multiprocessing")
    parser.add_argument("--input-fastq", type=Path, required=True, help="Path to input FASTQ file")
    parser.add_argument("--malformed-output-dir", type=Path, default="./integrity_check_malformed", help="Output directory for malformed records")
    parser.add_argument("--store-malformed", action="store_true", help="Store malformed records in output directory")
    parser.add_argument("--processes", type=int, default=4, help="Number of processes to use")
    parser.add_argument("--log-dir", type=Path, default=None, help="Directory to store log files")
    
    args = parser.parse_args()
    main(
        input_fastq=Path(args.input_fastq),
        malformed_output_dir=Path(args.malformed_output_dir),
        store_malformed=args.store_malformed,
        processes=args.processes,
        log_dir=Path(args.log_dir) if args.log_dir else None
    )