from pathlib import Path
import shutil
import re

# --- Import Utilities ---
from .file_management_utils import retain_decompressed_duplicate_paths, compression_utility, fetch_input_file_list
from .command_utils import execute_command

# ===================================
# Functions for creating genome index
# ===================================
def extract_decoys_from_fasta(
        fasta_file: Path, 
        output_file: Path):
    """
    Extract sequence identifiers from a FASTA file, 
    removing the '>' character and keeping only the first field (before any spaces). 
    Part of Salmon's index generation.
    
    Equivalent to:
    grep "^>" GRCm38.primary_assembly.genome.fa | cut -d " " -f 1 > decoys.txt
    sed -i.bak -e 's/>//g' decoys.txt
    
    Parameters:
    -----------
        fasta_file (Path): Path to the input FASTA file
        output_file (Path): Path to the output file to write the decoys
    """
    with open(fasta_file, 'rb') as infile, open(output_file, 'wb') as outfile:
        for line in infile:
            if line.startswith(b'>'):
                # Extract the first field (before any spaces) and remove '>'
                identifier = line.strip().split(b' ')[0][1:]
                outfile.write(identifier + b'\n')

# ---------------------------------------------
# --- Main genome index generation function ---
# ---------------------------------------------
def salmon_genome_index_for_species(
        species_name: str,
        input_genome_dir: Path,
        output_genome_index_dir: Path,
        genome_indexing_threads: int = 4):
    """
    Generates a Salmon index for a given species.
    
    Parameters:
    -----------
    species_name: str
        The name of the species. This is used to name the output directory and files.
        Whitespaces, hyphens, and periods are accepted but will be replaced with underscores.
    input_genome_dir: Path
        The directory containing the input files for the species. 
        The sequence files must be in FASTA format, ending with '.fa' or '.fasta', 
        and can be compressed with gzip (.gz) or zstd. 
    output_genome_index_dir: Path
        The directory where the species index directory will be generated.
        The output structure is: <output_genome_index_dir> / <species_dir> / <index files>
    input_compression_threads: int
        The number of threads to use for zstd compression of the input files.
    genome_indexing_threads: int
        The number of threads to use for the genome generation step.
    """
    # Initialize variables
    cdna_file_list = []
    genome_file_list = []
    temp_index_dir = None
    sp_index_dir = None
    decompressed_genome_file = None
    decompressed_cdna_file = None
    print("Generating genome index for species:", species_name)
    
    try:
        # Check input arguments 
        if not input_genome_dir.is_dir():
            raise FileNotFoundError(f"Input genome directory not found: {input_genome_dir}")

        # --- Fetch the input FASTA files ---
        # -----------------------------------
        # Locate all .fa and .fasta files in the input genome directory 
        fasta_file_extensions = ['.fa.gz', '.fa.zst', '.fasta.gz', '.fasta.zst', '.fa', '.fasta']
        fasta_file_list = fetch_input_file_list(
            input_dir = input_genome_dir,
            target_file_extensions = fasta_file_extensions
        )
        if not fasta_file_list or not len(fasta_file_list) > 0:
            raise ValueError(f"No input FASTA files found in {input_genome_dir}")
        
        # Retain only the decompressed duplicates if possible
        fasta_file_list = retain_decompressed_duplicate_paths(fasta_file_list)
        
        # Separate out the cDNA files. Expect one cDNA file.
        for file in fasta_file_list:
            if 'cdna' in file.name:
                cdna_file_list.append(file)
        if not cdna_file_list:
            raise ValueError(f"No cDNA files found in {input_genome_dir}")
        elif len(cdna_file_list) > 1:
            raise ValueError(f"Multiple cDNA files found in {input_genome_dir}: {cdna_file_list}")
        cdna_file = Path(cdna_file_list[0])
        
        # Separate out the genome files. Expect only one genome file.
        for file in fasta_file_list:
            if 'primary_assembly' in file.name:
                genome_file_list.append(file)
        if not genome_file_list:
            raise ValueError(f"No genome files found in {input_genome_dir}")
        elif len(genome_file_list) > 1:
            raise ValueError(f"Multiple genome files found in {input_genome_dir}: {genome_file_list}")
        genome_file = Path(genome_file_list[0])

        # Genome and cDNA files cannot be the same
        if genome_file == cdna_file:
            raise ValueError(f"The genome and cDNA files cannot be the same: {genome_file}")
        print("Using genome sequence file:", genome_file)
        print("Using cDNA sequence file:", cdna_file)
        
        # --- Decompress the input files if needed ---
        # --------------------------------------------
        # Decompress the genome if needed
        if genome_file.suffix in ['.gz','.zst']:
            decompressed_genome_file = genome_file.with_suffix('')
            print(f"Decompressing genome file to: {decompressed_genome_file}... ", end='')
            decompression_success = compression_utility(
                input_path=genome_file,
                output_path=decompressed_genome_file,
                compress=False,
            )
            if not decompression_success \
                or not decompressed_genome_file.is_file() \
                    or not decompressed_genome_file.stat().st_size > 0:
                if decompressed_genome_file.is_file():
                    decompressed_genome_file.unlink()
                raise RuntimeError(f"Failed to decompress {genome_file} to {decompressed_genome_file}")
            print("Done.")
            genome_file = decompressed_genome_file
        
        # Decompress the cDNA if needed
        if cdna_file.suffix in ['.gz','.zst']:
            decompressed_cdna_file = cdna_file.with_suffix('')
            print(f"Decompressing cDNA file to: {decompressed_cdna_file}... ", end='')
            decompression_success = compression_utility(
                input_path=cdna_file,
                output_path=decompressed_cdna_file,
                compress=False,
            )
            if not decompression_success \
                or not decompressed_cdna_file.is_file() \
                    or not decompressed_cdna_file.stat().st_size > 0:
                if decompressed_cdna_file.is_file():
                    decompressed_cdna_file.unlink()
                raise RuntimeError(f"Failed to decompress {cdna_file} to {decompressed_cdna_file}")
            print("Done.")
            cdna_file = decompressed_cdna_file
        
        
        # --- Generate metadata files ---
        # -------------------------------
        # Create the species genome index directory
        # In species_name, replace whitespaces, hyphens, and periods with underscores
        species_name = re.sub(r'\s', '_', species_name)
        species_name = re.sub(r'[-.]+', '_', species_name)
        # Directory for temporary files
        temp_index_dir = output_genome_index_dir / f"{species_name}_genome_index_temp"
        if temp_index_dir.is_dir():
            shutil.rmtree(temp_index_dir)
        temp_index_dir.mkdir(parents=True, exist_ok=False)
        # index output directory
        sp_index_dir = output_genome_index_dir / f"{species_name}_genome_index"
        if sp_index_dir.is_dir():
            shutil.rmtree(sp_index_dir)
        sp_index_dir.mkdir(parents=True, exist_ok=False)

        # Create the decoy file
        print("Creating decoy file... ", end='')
        decoy_file = Path(temp_index_dir) / 'decoys.txt'
        if decoy_file.is_file():
            decoy_file.unlink()
        extract_decoys_from_fasta(fasta_file=genome_file, output_file=decoy_file)
        if not decoy_file.is_file() or not decoy_file.stat().st_size > 0:
            raise RuntimeError(f"Failed to create decoy file {decoy_file}")
        print("Done.")

        # Concatenate the genome and cDNA files
        print("Concatenating genome and cDNA files... ", end='')
        concatenated_file = Path(temp_index_dir) / 'gentrome.fa'
        if concatenated_file.is_file():
            concatenated_file.unlink()
        with open(concatenated_file, 'wb') as outfile:
            # Genome targets (decoys) should come after the transcriptome
            with open(cdna_file, 'rb') as infile:
                outfile.write(infile.read())
            with open(genome_file, 'rb') as infile:
                outfile.write(infile.read())
        if not concatenated_file.is_file() or not concatenated_file.stat().st_size > 0:
            raise RuntimeError(f"Failed to concatenate {genome_file} and {cdna_file} to {concatenated_file}")
        print('Done.')
        
        # --- Generate genome index ---
        # -----------------------------
        print("\nRunning Salmon genome index generation")
        print("--------------------------------------")
        salmon_cmd = [
            'salmon',
            'index',
            '-t', str(concatenated_file),
            '-d', str(decoy_file),
            '-p', str(genome_indexing_threads),
            '-i', str(sp_index_dir),
            '-k', '19' # Must be an odd value
        ]
        # From Salmon documentation about -k parameter:
        # "We find that a k of 31 seems to work well for reads of 75bp or longer, 
        # but you might consider a smaller k if you plan to deal with shorter reads." 
        
        exit_status = execute_command(salmon_cmd)
        if exit_status != 0:
            raise RuntimeError(f"Genome generation failed with exit status {exit_status}")
        # Check if the genome index directory is empty
        elif next(sp_index_dir.iterdir(), None) is None:
            shutil.rmtree(sp_index_dir)
            raise RuntimeError(f"Genome generation resulted in empty output directory {sp_index_dir}")

        # If we have reached this point, the genome index generation was successful
        print("Genome index generation successful.")
    
    finally:
        # --- Remove temporary files ---
        if temp_index_dir and temp_index_dir.is_dir():
            shutil.rmtree(temp_index_dir)
        if sp_index_dir and next(sp_index_dir.iterdir(), None) is None:
            shutil.rmtree(sp_index_dir)
        if decompressed_cdna_file and decompressed_cdna_file.is_file():
            decompressed_cdna_file.unlink()
        if decompressed_genome_file and decompressed_genome_file.is_file():
            decompressed_genome_file.unlink()

# =======
# Mapping
# =======
def salmon_map_reads(
        input_reads_fastq: Path,
        input_index_dir: Path,
        output_path_prefix: str,
        mapping_threads: int = 4):
    """
    
    """
    input_fastq_to_process = None
    decompressed_reads_file = None
    try:
        # --- Decompress the reads input if needed ---
        # --------------------------------------------
        if input_reads_fastq.suffix in ['.gz','.zst']:
            decompressed_reads_file = input_reads_fastq.with_suffix('')
            print(f"Decompressing reads file to: {decompressed_reads_file}... ", end='')
            decompression_success = compression_utility(
                input_path=input_reads_fastq,
                output_path=decompressed_reads_file,
                compress=False,
            )
            if not decompression_success \
                or not decompressed_reads_file.is_file() \
                    or not decompressed_reads_file.stat().st_size > 0:
                if decompressed_reads_file.is_file():
                    decompressed_reads_file.unlink()
                raise RuntimeError(f"Failed to decompress {input_reads_fastq} to {decompressed_reads_file}")
            print("Done.")
            input_fastq_to_process = decompressed_reads_file
        else:
            input_fastq_to_process = input_reads_fastq
        
        # Salmon mapping command
        mapping_command = [
            'salmon', 'quant',
            '-i', str(input_index_dir),
            '-p', str(mapping_threads),
            '-l', 'A',
            '-r', str(input_fastq_to_process),
            '-o', str(output_path_prefix),
            '--posBias'
            ]
        exit_status = execute_command(mapping_command)
        if exit_status != 0:
            raise RuntimeError(f"Mapping failed with exit status {exit_status}")
    finally:
        if decompressed_reads_file and decompressed_reads_file.is_file():
            decompressed_reads_file.unlink()