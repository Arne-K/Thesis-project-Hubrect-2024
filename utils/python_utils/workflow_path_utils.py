from pathlib import Path

utils_dir = Path(__file__).resolve().parent.parent
if utils_dir.name != "utils":
    raise ValueError(f"Invalid utils directory: {utils_dir}")

# --- Main directory paths ---
workflow_dir = utils_dir.parent
datasets_dir = workflow_dir / "datasets"
#
# Input/Output Directory Paths
# ----------------------------
tools_dir = workflow_dir / "tools"
log_dir = workflow_dir / "logs"
raw_reads_dir = datasets_dir / "Tomoseq_reads_raw"
qc_reads_dir = datasets_dir / "Tomoseq_reads_qc"
#
# Support Data Paths
# ------------------
barcode_list_file = datasets_dir / "Barcodes_CelSeq1.tsv"
bcode_permut_dict_file = tools_dir / "BcCorrect_Dictionary.json"
#
# Program Paths
# -------------
unambig_bcode_permuts = tools_dir / "unambig_bcode_permuts.py"
r1_demultiplex_program = tools_dir / "fastq_demultiplexer" / "R1_demultiplexer.py"
r2_demultiplex_program = tools_dir / "fastq_demultiplexer" / "R2_demultiplexer.py"
fastq_integrity_check = tools_dir / "fastq_integrity_check.py"
fastq_match_edit = tools_dir / "fastq_match_edit_tool.py"

datasets_dir.mkdir(exist_ok=True)
log_dir.mkdir(exist_ok=True)
tools_dir.mkdir(exist_ok=True)
raw_reads_dir.mkdir(exist_ok=True)
qc_reads_dir.mkdir(exist_ok=True)