import sys
from pathlib import Path
from typing import List
import zstandard as zstd
import tarfile
import gzip

# --- Import Utilities ---
from .logging_utils import format_short_path

# =============================
# Functions for searching files
# =============================
def fetch_r_files(read_type: str,
                input_dir: Path) -> List[Path]:
    """
    Finds all FASTQ files of a given read type in a directory.
    
    Parameters:
    -----------
    read_type: str
        Either 'R1' or 'R2'.
        
    Returns:
    --------
    List[Path]
        A list of Path objects for the FASTQ files of the given read type.
    """
    if read_type not in ['R1', 'R2']:
        raise ValueError(f"Invalid read type: {read_type}")
    if not input_dir.is_dir():
        raise FileNotFoundError(f"Input directory not found: {input_dir}")
    
    all_fastq_files = []
    target_files = []
    for file in input_dir.rglob("*"):
        if file.is_file():
            all_fastq_files.append(file)
    if not all_fastq_files:
        raise FileNotFoundError(f"No valid FASTQ files found in {input_dir}")
    
    if read_type == "R1":
        opposite_read_type = "R2"
    elif read_type == "R2":
        opposite_read_type = "R1"
    
    # --- Positive search conditions ---
    # ----------------------------------
    for file in all_fastq_files:
        file_name = file.name # Get the filename string
        if (file_name.startswith(f"{read_type}_") or    # Condition 1: Starts with R1_
            f"_{read_type}_" in file_name or          # Condition 2: Contains _R1_
            f"_{read_type}." in file_name):           # Condition 3: Contains _R1.
            target_files.append(file)
    if not target_files:
        raise FileNotFoundError(f"No FASTQ files found in {input_dir} with read type {read_type}")
    
    # --- Negative search conditions ---
    # ----------------------------------
    for file in target_files:
        file_name = file.name # Get the filename string
        if (file_name.startswith(f"{opposite_read_type}_") or    
            f"_{opposite_read_type}_" in file_name or         
            f"_{opposite_read_type}." in file_name):          
            target_files.remove(file)
    if not target_files:
        raise FileNotFoundError(f"No FASTQ files found in {input_dir} with read type {read_type}")
    
    print(f"Found {len(target_files)} FASTQ files in {input_dir.name} with read type {read_type}")
    return target_files

def fetch_input_file_list(
        input_dir: Path,
        target_file_extensions: List[str]
        ) -> List[Path]:
    """
    Fetches a list of files from a directory based on file extensions.
    
    Parameters:
    -----------
    input_dir: Path
        The directory to fetch files from.
    target_file_extensions: List[str]
        A list of valid file extensions to filter the files by.
    """
    target_file_list = []
    all_file_list = list(input_dir.rglob("*"))
    all_file_list = [path for path in all_file_list if path.is_file()]
    if not all_file_list:
        print(f"No files found in {input_dir}")
        return None
    for ext in target_file_extensions:
        for file in all_file_list:
            if str(file).endswith(ext):
                target_file_list.append(file)
    # Raise an error if no input genome FASTA files were found
    if not target_file_list or len(target_file_list) == 0:
        print(f"No files were found in {input_dir} with extensions {target_file_extensions}")
        return None
    return target_file_list

def retain_decompressed_duplicate_paths(
        file_path_list: List[Path]) -> List[Path]:
    """
    Function to check if the input contains compressed-decompressed duplicates.
    If there is any, retain only the decompressed duplicate.
    
    Parameters:
    -----------
    file_path_list: List[Path]
        List of file paths to check for duplicates.
    """
    paths_to_remove = []
    checked_file_path_list = file_path_list.copy()
    for path in file_path_list:
        if path.suffix in ['.gz','.zst']:
            decompressed_path = path.with_suffix('')
            # If the decompressed version is found among the input file list, remove the compressed version
            if decompressed_path in file_path_list:
                paths_to_remove.append(path)
    if len(paths_to_remove) > 0:
        checked_file_path_list = [path for path in file_path_list 
                                if path not in paths_to_remove]
    return checked_file_path_list

# ==========================
# Compressed data management 
# ==========================
def manage_directory_content_compression(
        input_dir_path: Path,
        compression_threads: int = 4) -> bool:
    """
    This function makes sure that files in a directory are maintained in a compressed state, 
    particularly in .zst format because it is fast, reliable and customizable. The function is especially meant
    to be used in functions that work with large files to maintain them in compressed state and only decompress them if needed.
    This function does the following: 
        1. accepts a directory path
        2. checks for compressed-decompressed file pairs - of file type of any kind, as long as the format is .<file_ext>.<compression_ext>
        3. if there is no compressed pair, the file is compressed to .zst
        4. if the compressed pair is in .zst format, the decompressed version is removed.
        5. if the compressed pair is in .gz format, or if there is just the .gz file, the file is compressed to .zst
    ***other compression formats are not supported currently and are ignored***

    Parameters:
    -----------
    input_dir_path: Path
        The directory to manage the content of.
    compression_threads: int
        Threads to use for zst compression by the compression_utility() function.
    
    Returns:
    --------
    True for successful completion, False otherwise.
    """
    compression_success = None
    
    # Check if they are type Path
    if not isinstance(input_dir_path, Path):
        print(f"Error: Input directory '{input_dir_path}' is not a Path object.")
        return False
    # Check the directory path
    if not input_dir_path.is_dir():
        print(f"Error: Input directory '{input_dir_path}' is not a valid directory.")
        return False
    # Check if the directory is empty
    if next(input_dir_path.iterdir(), None) is None:
        print(f"Error: Input directory '{input_dir_path}' is empty.")
        return False
    
    print('\n--------------------------')
    print(f"Directory content manager: {format_short_path(input_dir_path)}")
    # Fetch all files from the directory
    all_file_list = list(input_dir_path.rglob("*"))
    all_file_list = [path for path in all_file_list if path.is_file()]
    if not all_file_list:
        print(f"Error: No files found in {input_dir_path}")
        return False
    
    # Step 1: First handle all .zst files - remove any uncompressed duplicates
    zst_files = [f for f in all_file_list if f.suffix == '.zst']
    for file in zst_files:
        decompressed_file = file.with_suffix('')
        if decompressed_file.is_file():
            decompressed_file.unlink()
            print(f"Removed {decompressed_file.name}")
            print(f"Retained compressed file: {file.name}")
            # Update the main list
            if decompressed_file in all_file_list:
                all_file_list.remove(decompressed_file)
    
    # Step 2: Handle all .gz files - convert to .zst
    gz_files = [f for f in all_file_list if f.suffix == '.gz' and f.is_file()]
    for file in gz_files:
        # Look for a decompressed version of the .gz file
        decompressed_file = file.with_suffix('')
        if decompressed_file.is_file():
            # Remove the decompressed version
            decompressed_file.unlink()
            print(f"Removed {decompressed_file.name}")
        
        # --- Convert the .gz file to .zst ---
        # First decompress the .gz
        print(f"Converting .gz file {file.name} to .zst")
        decompression_success = compression_utility(
            input_path=file,
            output_path=decompressed_file,
            compress=False
        )
        if not decompression_success or not decompressed_file.is_file() or not decompressed_file.stat().st_size > 0:
            print(f"Failed to decompress {file} to {decompressed_file}")
            if decompressed_file.is_file():
                decompressed_file.unlink()
            return False
        
        # Then compress the decompressed file to .zst
        compressed_file = file.with_suffix('.zst')
        compression_success = compression_utility(
            input_path=decompressed_file,
            output_path=compressed_file,
            compress=True,
            threads=compression_threads
        )
        if not compression_success \
            or not compressed_file.is_file() \
                or not compressed_file.stat().st_size > 0:
            print(f"Failed to compress {decompressed_file} to {compressed_file}")
            # remove any remnants of the failed compression
            if compressed_file.is_file():
                compressed_file.unlink()
            return False
        # After successful compression to .zst, remove the parent .gz file
        file.unlink()
        print(f"Removed {file.name}")
        # Remove the decompressed file
        decompressed_file.unlink()
        print(f"Removed {decompressed_file.name}")
    
    # Step 3: Compress any remaining uncompressed files
    # Refresh the file list to get the latest state of the directory
    remaining_files = list(input_dir_path.rglob("*"))
    remaining_files = [path for path in remaining_files if path.is_file() and path.suffix != '.zst']
    
    for file in remaining_files:
        # Skip any .zst files or files that don't exist anymore
        if not file.is_file() or file.suffix == '.zst':
            continue
            
        # Check if a compressed version already exists
        compressed_file = Path(f"{file}.zst")
        if compressed_file.is_file():
            # If the compressed version exists, remove the original
            file.unlink()
            print(f"Removed duplicate {file.name} (compressed version exists)")
        else:
            # Compress the file
            print(f"Compressing {file.name} to .zst")
            compression_success = compression_utility(
                input_path=file,
                output_path=compressed_file,
                compress=True,
                threads=compression_threads
            )
            if not compression_success or not compressed_file.is_file() or not compressed_file.stat().st_size > 0:
                print(f"Failed to compress {file} to {compressed_file}")
                if compressed_file.is_file():
                    compressed_file.unlink()
                return False
            # Remove the original after successful compression
            file.unlink()
            print(f"Removed {file.name}")
            print(f"Retained compressed file: {compressed_file.name}")
    
    # Step 4: Final verification - check if any uncompressed files remain
    final_check = list(input_dir_path.rglob("*"))
    final_check = [f for f in final_check if f.is_file() and f.suffix != '.zst']
    if final_check:
        print(f"Warning: Found {len(final_check)} uncompressed files after processing:")
        for f in final_check:
            print(f"  {f.name}")
    
    print("Directory content management completed")
    return True

# =============================
# Compression utility functions
# =============================
def _compress_file(
        input_path: Path, 
        output_path: Path, 
        level: int, 
        threads: int, 
        chunk_size: int):
    """
    Helper function to compress a single file.
    """
    print(f"Compressing file '{input_path.name}' to '{output_path.name}'")
    print(f"  Level: {level}, Threads: {threads if threads != -1 else 'Auto'}")
    cctx = zstd.ZstdCompressor(level=level, threads=threads)
    try:
        # Ensure output directory exists
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(input_path, 'rb') as infile, open(output_path, 'wb') as outfile:
            with cctx.stream_writer(outfile) as compressor:
                while True:
                    chunk = infile.read(chunk_size)
                    if not chunk:
                        break
                    compressor.write(chunk)
        print("File compression successful.")
        return True
    except FileNotFoundError:
        print(f"Error: Input file '{input_path}' not found.", file=sys.stderr)
    except Exception as e:
        print(f"An error occurred during file compression: {e}", file=sys.stderr)
    return False

def _decompress_file(
        input_path: Path, 
        output_path: Path, 
        chunk_size: int):
    """
    Helper function to decompress a single zstandard file.
    """
    print(f"Decompressing zstd file '{input_path.name}' to '{output_path.name}'...")
    dctx = zstd.ZstdDecompressor()
    try:
         # Ensure output directory exists
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(input_path, 'rb') as infile, open(output_path, 'wb') as outfile:
            with dctx.stream_reader(infile) as decompressor:
                while True:
                    chunk = decompressor.read(chunk_size)
                    if not chunk:
                        break
                    outfile.write(chunk)
        print("File decompression successful.")
        return True
    except FileNotFoundError:
        print(f"Error: Input file '{input_path}' not found.", file=sys.stderr)
    except zstd.ZstdError as e:
        print(f"Zstandard decompression error: {e}. Is '{input_path}' a valid Zstandard file?", file=sys.stderr)
    except Exception as e:
        print(f"An error occurred during file decompression: {e}", file=sys.stderr)
    return False

def _decompress_gzip_file(
        input_path: Path, 
        output_path: Path, 
        chunk_size: int):
    """
    Helper function to decompress a single gzip file.
    """
    print(f"Decompressing gzip file '{input_path.name}' to '{output_path.name}'")
    try:
        # Ensure output directory exists
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with gzip.open(input_path, 'rb') as infile, open(output_path, 'wb') as outfile:
            while True:
                chunk = infile.read(chunk_size)
                if not chunk:
                    break
                outfile.write(chunk)
        print("Gzip file decompression successful.")
        return True
    except FileNotFoundError:
        print(f"Error: Input file '{input_path}' not found.", file=sys.stderr)
    except gzip.BadGzipFile as e:
        print(f"Gzip decompression error: {e}. Is '{input_path}' a valid gzip file?", file=sys.stderr)
    except Exception as e:
        print(f"An error occurred during gzip file decompression: {e}", file=sys.stderr)
    return False

def _compress_directory(
        input_path: Path, 
        output_path: Path, 
        level: int, 
        threads: int):
    """
    Helper function to compress a directory into a .tar.zst file.
    """
    # Ensure output path has the correct extension
    if not output_path.name.endswith(".tar.zst"):
       output_path = output_path.with_name(output_path.name + ".tar.zst")

    print(f"Compressing directory '{input_path.name}' \n\tto archive '{output_path.name}'...")
    print(f"  Level: {level}, Threads: {threads if threads != -1 else 'Auto'}")
    cctx = zstd.ZstdCompressor(level=level, threads=threads)
    try:
         # Ensure output directory exists
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, 'wb') as outfile:
            # Stream data through the zstd compressor into the tar writer
            with cctx.stream_writer(outfile) as compressor_stream:
                # Use tarfile in stream write mode ('w|') with the compressor stream
                with tarfile.open(fileobj=compressor_stream, mode='w|') as tar:
                    # Add the directory content; arcname='.' keeps paths relative inside tar
                    tar.add(input_path, arcname=input_path.name) # Use dir name as root inside tar
                    # Or use arcname='.' if you want files directly in the root of the tar
                    # tar.add(input_path, arcname='.')
        print("Directory compression successful.")
        return True
    except FileNotFoundError:
         print(f"Error: Input directory '{input_path}' not found.", file=sys.stderr)
    except Exception as e:
        print(f"An error occurred during directory compression: {e}", file=sys.stderr)
    return False

def _decompress_archive(
        input_path: Path, 
        output_path: Path):
    """
    Helper function to decompress a .tar.zst archive into a directory.
    """
    print(f"Decompressing archive '{input_path.name}' to directory '{output_path.name}'...")
    # Use specified threads, default to 1 if 0 was passed, -1 for auto
    dctx = zstd.ZstdDecompressor()
    try:
        # Ensure output directory exists
        output_path.mkdir(parents=True, exist_ok=True)
        with open(input_path, 'rb') as infile:
            # Stream data through the zstd decompressor into the tar reader
            with dctx.stream_reader(infile) as decompressor_stream:
                 # Use tarfile in stream read mode ('r|') with the decompressor stream
                with tarfile.open(fileobj=decompressor_stream, mode='r|') as tar:
                    tar.extractall(path=output_path)
        print("Archive decompression successful.")
        return True
    except FileNotFoundError:
        print(f"Error: Input archive '{input_path}' not found.", file=sys.stderr)
    except (zstd.ZstdError, tarfile.TarError) as e:
        print(f"Decompression error: {e}. Is '{input_path}' a valid .tar.zst archive?", file=sys.stderr)
    except Exception as e:
        print(f"An error occurred during archive decompression: {e}", file=sys.stderr)
    return False


def compression_utility(
        input_path: Path,
        output_path: Path,
        compress: bool,
        level: int = 4,
        threads: int = 4, # Default to 4, use -1 for auto, or specify number
        chunk_size: int = 1024 * 1024 * 16 # Default 16MB chunk size
        ) -> bool:
    """
    Compress/decompress a file or directory using Zstandard or Gzip.
    Directories are handled using the '.tar.zst' format.
    For decompression, supports both .zst (Zstandard) and .gz (Gzip) formats.

    Parameters:
    -----------
        input_path: Path
            Path to the file or directory to be compressed/decompressed.
        output_path: Path
            Path where the compressed/decompressed file or directory will be saved.
            If compressing, the output path must have a '.tar.zst' (for archives) or '.zst' extension (for single files).
            If decompressing an archive, this should be the target directory path.
        compress: bool
            True to compress, False to decompress.
        level: int
            Compression level (1-22 or negative for speed).
            Only used for Zstandard compression not decompression.
        threads: int
            Number of threads to use. Use -1 for auto-detect CPU cores.
            These are only used for compression, decompression is always single-threaded.
            Not bound by Python's GIL.
        chunk_size: int
            Size of chunks in bytes to read/write for file operations, default is 16MB.

    Returns:
    --------
        bool: True if the operation was successful, False otherwise.
    """
    # Input validation
    # Check if they are type Path
    if not isinstance(input_path, Path):
        print(f"Error: Input path '{input_path}' is not a Path object.")
        return False
    if not isinstance(output_path, Path):
        print(f"Error: Output path '{output_path}' is not a Path object.")
        return False
    if not input_path.exists():
        print(f"Error: Input path '{input_path}' does not exist.")
        return False
    # Output path validation
    if compress and output_path.suffix.lower() not in ['.tar.zst','.zst']:
        print(f"Error: Output path '{output_path}' must have a '.tar.zst' (for archives) or '.zst' extension (for single files).", file=sys.stderr)
        return False
    
    # --- Dispatch based on input type and operation ---
    if compress:
        if input_path.is_dir():
            # Compress Directory
            return _compress_directory(input_path, output_path, level, threads)
        elif input_path.is_file():
            # Compress File
            return _compress_file(input_path, output_path, level, threads, chunk_size)
        else:
            print(f"Error: Input path '{input_path}' is not a file or directory.", file=sys.stderr)
            return False
    else: # Decompress
        if input_path.is_dir():
             # Cannot decompress a directory
             print(f"Error: Cannot decompress a directory ('{input_path}'). Provide a '.zst', '.tar.zst', or '.gz' file.", file=sys.stderr)
             return False
        elif input_path.is_file():
            # Check file extension to determine compression type
            if input_path.name.endswith(".tar.zst"):
                # Decompress Archive (directory)
                return _decompress_archive(input_path, output_path, threads)
            elif input_path.name.endswith(".zst"):
                # Decompress Single zstd File
                # Ensure output_path is treated as a file path
                if output_path.exists() and output_path.is_dir():
                    print(f"Error: Output path '{output_path}' is a directory, expected a file path for single file decompression.", file=sys.stderr)
                    return False
                return _decompress_file(input_path, output_path, chunk_size)
            elif input_path.name.endswith(".gz"):
                # Decompress Single gzip File
                # Ensure output_path is treated as a file path
                if output_path.exists() and output_path.is_dir():
                    print(f"Error: Output path '{output_path}' is a directory, expected a file path for single file decompression.", file=sys.stderr)
                    return False
                return _decompress_gzip_file(input_path, output_path, chunk_size)
            else:
                print(f"Error: Input file '{input_path}' does not have a '.zst', '.tar.zst', or '.gz' extension.", file=sys.stderr)
                return False
        else:
             print(f"Error: Input path '{input_path}' is not a file.", file=sys.stderr)
             return False

