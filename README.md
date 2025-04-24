# Tomo-seq analysis for major research project (Hubrecht-2024)
This repository contains the code for creating count matrices from Tomo-seq data for major research project (Bioinformatics and Biocomplexity MSc, Hubrecht-2024).
## 1. From reads to counts
### Quality control and demultiplexing
1. Quality control of the raw reads with [fastp](https://github.com/OpenGene/fastp)
2. Barcode error-correction within Hamming distance = 1.
3. Demultiplexing the reads between 96 barcodes. The demultiplexing follows a two-step process: the error-corrected R1 files are first demultiplexed using the barcodes provided in the text file, and the R2 files are then demultiplexed based on the header line IDs in the demultiplexed R1 files. 
### Read Mapping and quantification 
- [Salmon](https://salmon.readthedocs.io/en/latest/): the reads are mapped to reference transcriptome, and transcript counts are summarized to gene-level counts with [tximport](https://bioconductor.org/packages/release/bioc/html/tximport.html).
- [STAR](https://github.com/alexdobin/STAR): the reads are mapped to reference genome with `--quantMode TranscriptomeSAM` and the gene-level quantification is done with [RSEM](https://github.com/deweylab/RSEM). 
## 2. Orthologous gene inference
Doing RNA-seq analysis between two different species means that we will compare the transcriptional expression between orthologs - genes that have originated from a common ancestor. It is easiest to compare single-copy orthologs (that have not undergone duplication after phylogenetic split), as neo-functionalization and hypo-functionalization make the comparison of multi-copy orthologs more complex. The only way to distinguish between single-copy and multi-copy orthologs is through phylogenetic reconstruction methods, and it cannot be done with sequence similarity alone (e. g. reciprocal BLAST). This project used [OrthoFinder](https://github.com/davidemms/OrthoFinder)
# Environment Setup
- Operating system: Linux or MacOS. Many tools are from Bioconda which does not support Windows. This codebase was created and run in WSL2 (Linux 5.15.167.4-microsoft-standard-WSL2) with Ubuntu 24.04.2 LTS.
- Python and R installed.
- At least 32 GB of RAM. The more the merrier. Multiple processor cores wouldn't hurt. 
- Conda (Miniconda3) with [Bioconda channel configuration](https://bioconda.github.io/) + Jupyter Notebook installed (or the Jupyter extension if running through VSCode). 
- The main code is run from Jupyter notebooks that use [Conda](https://www.anaconda.com/docs/getting-started/anaconda/install) environments as kernels. The notebooks can be run through VScode (with the Jupyter extension) or through the Jupyter web interface. The required notebook's conda environment (kernel) installations are displayed at the beginning of each notebook.
- When using Jupyter web interface, the Python conda environments have to be registered as Jupyter kernels from the command line within the active conda environment: `python -m ipykernel install --user --name=<env_name>`. For R, the kernel activation has to be done from within an R instance: `IRkernel::installspec(name = '<env_name>', displayname = '<display_kernel_name>')`
- When using the notebooks through VScode, the Python kernel registration is automatic after installing the ipykernel conda package, and it is not needed to run the ipykernel command manually. For R kernels, this is also true but requires the installation of `jupyter_client` conda package in addition to `IRkenel`.
## Tomo-seq Data
- The [Tomo-seq](https://pubmed.ncbi.nlm.nih.gov/27443932/) method is essentially bulk RNA-seq of thin cryosections, and the library amplification is done with barcoded primers to retain spatial information. This data was from infarcted animal hearts to perform differential gene expression analysis between two species: *Mus musculus* and *Acomys cahirinus*. The animals were sacrificed at time points 14- and 42 days after MI (n=2), and the hearts were cryosectioned in the direction of atrium --> apex into 96 × 80 μm slices. 
- After total RNA extraction, the mRNA was amplified with polyT primers with slice-specific 8-
bp barcodes. The 96 × 8-bp barcodes were designed to have a minimum Hamming distance of 2 between any
pair.
- The libraries were sequenced with NextSeq500, by paired-end sequencing which produced a pair of sequences per each
amplified fragment: a 26 bp (read 1) and 62 bp (read 2). Read 1 contained only the barcode
sequence and read 2 contained the sequence used for alignment, so practically from the
alignment standpoint these were single-end reads and Read 1 is only used for demultiplexing.
### Practical comments on the data
- the poly-T read amplification produced 3' end sequence fragments with an average length of ~50 bp, which means there is no gene length bias in the data and normalization methods like RPKM must not be used. There is, however, bias in the read start position on the fragments, so Salmon and RSEM must be run with settings that infer read start bias from the data. Reads map to forward strand (SF). 
- The minimum Hamming distance between any two barcodes is 2 (to be more precise, there are 28 barcodes that have a minimum HD of 3 to any other of the 96 barcodes, the rest have minimum HD=2). This is relevant to barcode error correction, where we can correct sequencing errors by looking for permutations of the original sequences. For an 8-base DNA sequence, the set size of permutations within HD=1 is 33. For two 8-base DNA sequences that are HD=2 apart, there will be 2 permutations that are identical between the sets. Thus, the barcode error correction has to account for collisions by excluding permutations that are not unique to the barcode. Among the 96 barcodes, there are 16 barcodes that have an unique permutation set size of 27 (the lowest), because they have HD=2 with 3 other barcodes: 33-(3*2)=27. The `unambig_bcode_permuts.py` program accounts for these collisions. From practical standpoint, the bias from the uneven set sizes is minimal, as most (>90%) of the barcodes do not need error correction after QC anyway.
  - 28 barcodes have the full set of 33
  - 24 barcodes have 31 unambiguous sequences
  - 28 barcodes have 29 unambiguous sequences
  - 16 barcodes have 27 unambiguous sequences
