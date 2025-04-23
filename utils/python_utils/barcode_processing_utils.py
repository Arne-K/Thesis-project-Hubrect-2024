from pathlib import Path
from typing import Tuple, List, Optional
import csv
import re


def read_barcode_file(file_path: Path) -> List[Tuple[str, int]]:
    """
    Reads a tab-separated file where barcodes are in the first column,
    and barcode ID (like spatial slice number for Tomo-seq) is in the second column.
    The first row is treated as a header and skipped.
    The function returns a list of tuples, where each tuple contains a barcode and its ID.
    
    Parameters:
    -----------
    file_path : Path
        Path to the tab-separated file containing barcodes and IDs.
        
    Returns:
    --------
    List[Tuple[str, str]]
        a list of tuples, where each tuple contains a barcode and its ID
    """
    barcode_info_list = []
    
    with open(file_path, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        
        # Check if file is empty
        first_row = next(reader, None)
        if not first_row:
            raise ValueError(f"File {file_path} is empty.")
        
        # Check if the first column contains valid DNA barcodes in rest of the rows
        for row in reader:
            first_col = row[0].strip()
            barcode_seq = first_col.upper()
            is_valid_barcode = all(c in 'ACGTN' for c in barcode_seq)
            if not is_valid_barcode:
                raise ValueError(f"Invalid barcode sequence '{barcode_seq}' found in file {file_path}.")
            else:
                barcode_id = row[1].strip() if len(row) > 1 and row[1] else None
                barcode_id = int(barcode_id) if barcode_id and barcode_id.isdigit() else None
                barcode_info_list.append((barcode_seq, barcode_id))
    return barcode_info_list

def extract_barcode_from_filename(
        filename: Path) -> Optional[str]:
    """
    Extract barcode from filename.
    """
    # First try to match barcode at the beginning followed by underscore
    matches = re.findall(r'^([ACGTN]+)_', filename)
    if matches:
        return matches[0]
    
    # If that doesn't work, try to match barcode between underscores
    matches = re.findall(r'_([ACGTN]+)_', filename)
    if matches:
        return matches[0]
        
    return None

def barcode_info_from_filename(filename: str) -> Tuple[str, int]:
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
        r'_([ACGT]{8})_',   # Barcode flanked by two underscores
        r'_([ACGT]{8})$'   # Barcode preceded by an underscore at the end of path
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