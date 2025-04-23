from pathlib import Path
import subprocess
import shutil
import tempfile # Added for safe temporary file creation
import time

# --- Import Utilities ---
from .file_management_utils import compression_utility, format_short_path
from .command_utils import execute_command
from .fastq_processing_utils import validate_fastq

# ======================
# Functions to run fastp 
# ======================
def run_fastp(
        input_fastq: Path,
        output_fastq: Path,
        log_dir: Path,
        read_type: str) -> bool:
    """
    Runs fastp on the input FASTQ file.
    
    Parameters:
    -----------
    input_fastq: Path
        Path to the input FASTQ file.
    output_dir:
        Path to the output directory where the prefix-specific output directory will be created.
    read_type: str
        Either 'R1' or 'R2'.
    log_dir: Path
        Path to the log directory where the log files will be created.
    
    Output:
    -------
    QC output structure: <output_dir> / <prefix_dir> / <{read_type}_prefix_dir> / output files
    log output structure: <log_dir> / <prefix_dir> / <.json_dir> + <.html_dir>

    Returns:
    --------
    True for successful completion, False otherwise.
    """
    # --- Check and initiate input/output ---
    # ---------------------------------------
    # Check input fastq
    if not input_fastq.is_file():
        print(f"Input file not found: {input_fastq}")
        return False
    if not validate_fastq(input_fastq):
        print(f"Invalid input FASTQ file: {input_fastq}")
        return False
    # Check read type value
    if not read_type in ['R1', 'R2']:
        print(f"Invalid read type: {read_type}. Accepted values: 'R1' or 'R2'.")
        return False
    
    # Initialize variables
    was_successful = False
    temp_output_dir = None
    temp_decompressed_file_path = None
    
    try:
        # Temporary directory for fastp split output
        temp_output_dir = output_fastq.parent / "temp_output"
        if temp_output_dir.is_dir():
            shutil.rmtree(temp_output_dir)
        temp_output_dir.mkdir(parents=True, exist_ok=False)

        # Initialize output file paths
        if output_fastq.is_file():
            output_fastq.unlink()
        file_prefix = output_fastq.name.split(".")[0]
        # temporary split output
        temp_split_output = temp_output_dir / f"{file_prefix}.split.fastq"
        
        # log output
        html_log_dir = log_dir / "fastp_html_logs"
        json_log_dir = log_dir / "fastp_json_logs"
        html_log_dir.mkdir(parents=True, exist_ok=True)
        json_log_dir.mkdir(parents=True, exist_ok=True)
        output_html = html_log_dir / f"{file_prefix}.log.html"
        output_json = json_log_dir / f"{file_prefix}.log.json"
        if output_html.is_file():
            output_html.unlink()
        if output_json.is_file():
            output_json.unlink()

        # --- Decompress Input File if Needed (using temp file) ---
        if input_fastq.suffix in ['.gz', '.zst']:
            print(f"Input file {input_fastq.name} is compressed. Decompressing...")
            start_time = time.time()
            temp_decompressed_file_path = output_fastq.with_suffix("")
            print(f"Decompressing to temporary file: {temp_decompressed_file_path.name}")
            decompression_success = compression_utility(
                input_path=input_fastq,
                output_path=temp_decompressed_file_path,
                compress=False
            )
            if not decompression_success or not temp_decompressed_file_path.is_file() or temp_decompressed_file_path.stat().st_size == 0:
                # Clean up partially created temp file if it exists
                if temp_decompressed_file_path:
                     temp_decompressed_file_path.unlink(missing_ok=True)
                raise RuntimeError(f"Failed to decompress {input_fastq} to temporary file.")

            input_fastq = temp_decompressed_file_path # Workers will use the decompressed temp file
            elapsed_time = time.time() - start_time
            print(f"Decompression complete ({elapsed_time:.2f}s). Processing temporary file.")
        else:
             print("Input file is not compressed.")
        
        # --- Execute fastp ---
        # ---------------------
        print("Running fastp...")
        if read_type == "R1":
            fastp_cmd = ["fastp", 
                        "-i", str(input_fastq), 
                        "-o", str(temp_split_output), 
                        "--split", "4",
                        "-h", str(output_html), 
                        "-j", str(output_json),
                        "--length_required", "8", 
                        "--length_limit", "8",
                        "--qualified_quality_phred", "20", 
                        "--unqualified_percent_limit", "25", 
                        "--average_qual", "27", 
                        "--n_base_limit", "1",
                        "--thread", "4"
                        ]
        
        elif read_type == "R2":
            fastp_cmd = ["fastp", 
                        "-i", str(input_fastq), 
                        "-o", str(temp_split_output), 
                        "--split", "4",
                        "-h", str(output_html), 
                        "-j", str(output_json),
                        "--length_required", "20", 
                        "--cut_front", 
                        "--cut_tail", 
                        "--cut_window_size", "2", 
                        "--cut_mean_quality", "20",
                        "--qualified_quality_phred", "20", 
                        "--unqualified_percent_limit", "25", 
                        "--average_qual", "27", 
                        "--n_base_limit", "3",
                        "--thread", "4"
                        ]
        # Execute fastp
        exit_status = execute_command(fastp_cmd)
        if exit_status != 0:
            raise subprocess.CalledProcessError(exit_status, fastp_cmd)
        
        # --- Check output files ---
        if next(temp_output_dir.iterdir(), None) is None:
            was_successful = False
            print(f"Error: No split output files found in {temp_output_dir}")
            shutil.rmtree(temp_output_dir)
        else:
            for path in temp_output_dir.iterdir():
                if not validate_fastq(path):
                    path.unlink()
                    print(f"non-FASTQ file found in {format_short_path(temp_output_dir)}: {path.name}")
                    print(f"Removed {path.name}")
                elif not path.stat().st_size > 0:
                    path.unlink()
                    print(f"Removed empty FASTQ file in {format_short_path(temp_output_dir)}: {path.name}")
            dir_is_empty = next(temp_output_dir.iterdir(), None) is None
            if not dir_is_empty:
                was_successful = True
                        
    finally:
        # --- Remove the temporary decompressed file if one was created ---
        if temp_decompressed_file_path and temp_decompressed_file_path.exists():
            temp_decompressed_file_path.unlink(missing_ok=True)
            print(f"Removed temporary decompressed file: {temp_decompressed_file_path.name}")
        
        # --- Combine the split output files if there are any ---
        if temp_output_dir and temp_output_dir.is_dir():
            dir_is_empty = next(temp_output_dir.iterdir(), None) is None
            if dir_is_empty:
                print(f"No split output files found in {format_short_path(temp_output_dir)}")
                shutil.rmtree(temp_output_dir)
                print(f"Removed empty directory {format_short_path(temp_output_dir)}")
            elif was_successful:
                num_split_files = len(list(temp_output_dir.iterdir()))
                print(f"Combining ({num_split_files}) split output files in {format_short_path(temp_output_dir)}")
                with open(output_fastq, 'wb') as f_out:
                    for path in temp_output_dir.iterdir():
                        with open(path, 'rb') as f_in:
                            shutil.copyfileobj(f_in, f_out)
                if not output_fastq.stat().st_size > 0:
                    print(f"Error: Failed to combine split output files in {temp_output_dir}")
                    was_successful = False
                else:
                    print(f"Completed fastp run on {input_fastq.name}")
                    shutil.rmtree(temp_output_dir)
        return was_successful

