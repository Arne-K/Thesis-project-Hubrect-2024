library(tximport)
library(dplyr)
library(fs)

salmon_summarize_gene_level <- function(gtf_df, salmon_results_dir) {
    
    # ==========================================================
    # Summarize gene counts from Salmon mapping results.

    # Parameters:
    # -----------
    # gtf_df: Tibble of an annotation GTF file loaded with rtracklayer
    # salmon_results_dir: Path to directory containing Salmon .sf files. The file names should contain the barcode.
    # output_file: Path to output file to write the summarized counts to. The file will be a tab-separated
    # table with columns 'gene_id' and one column for each barcode.
    # ==========================================================
    
    # Check if the input directory exists
    if (!dir.exists(salmon_results_dir)) {
        stop(paste("Salmon directory not found:", salmon_results_dir))
    }

    # Fetch all .sf files in the directory
    sf_files <- dir_ls(salmon_results_dir, glob = "*.sf", recurse = TRUE)
    if (length(sf_files) == 0) {
        stop(paste("No Salmon .sf files found in directory:", salmon_results_dir))
    } else {
        message(paste("Found", length(sf_files), "Salmon .sf files in directory:", salmon_results_dir))
    }

    # Create a vector of barcode IDs for a named file vector. The names will be the column names in the final output.
    sf_file_names <- basename(sf_files)
    sf_file_names <- str_replace(sf_file_names, "_Salmon_quant.sf", "")
    names(sf_files) <- sf_file_names
    
    # Gene ID to transcript ID mapping
    tx2gene<- gtf_df %>%
        filter(type == "transcript") %>%
        select(transcript_id, gene_id) %>%
        unique()
    
    # Summarize gene counts
    txi.salmon <- tximport(sf_files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE)

    return(txi.salmon)
}