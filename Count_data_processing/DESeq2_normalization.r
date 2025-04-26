library("DESeq2")

# Logic outline:
# Step 1: estimate size factors on expression-filtered count data. 
# Step 2: apply the size factors to un-filtered count data to obtain normalized counts.
# The point is that the median of ratios normalization factor estimation benefits from expression filtering,
# However, in the downstream analysis we might be interested in genes whose expression is zero in one species and non-zero in another. 

input_dir_path <- "/home/a.kukkonen_cbs-niob.local/PartB/deseq2_input_data/tomoseq_datasets/"
out_dir_path="/home/a.kukkonen_cbs-niob.local/PartB/deseq2_output/tomoseq_normalized_counts_2/"
if (!dir.exists(out_dir_path)) {
    stop("output directory not found")
}

# Load the count data and meta data for both species
ac_count_data <- readRDS(paste0(input_dir_path, "AC_Tomoseq_counts_1-1orth_unnorm.RData"))
mm_count_data <- readRDS(paste0(input_dir_path, "MM_Tomoseq_counts_1-1orth_unnorm.RData"))

ac_meta_data <- readRDS(paste0(input_dir_path, "AC_Tomoseq_meta_1-1orth_unnorm.RData"))
mm_meta_data <- readRDS(paste0(input_dir_path, "MM_Tomoseq_meta_1-1orth_unnorm.RData"))

# Ensure that the Condition column in the meta data is class factor
mm_meta_data$Condition <- as.factor(mm_meta_data$Condition)
mm_meta_data$Condition <- as.factor(mm_meta_data$Condition)

# function for checking dataset format correctness
check_input_format <- function(count_data, meta_data) {
    dataset_name <- deparse(substitute(count_data))
    
    count_check <- colnames(count_data)
    col_check <- rownames(meta_data)
    check <- identical(count_check, col_check) # colData rownames and countData colnames must be identical
    if (check!=TRUE) {
        stop(paste0(dataset_name, ": count data column names are not identical to meta data row names"))
    } else {
        cat(paste0(dataset_name, ": dataset format correctness confirmed\n"))
    }

    cat(paste0("Any missing values in", dataset_name, ": ", anyNA(count_data), "\n"))
    cat(paste0("Samples in", dataset_name, ": ", length(colnames(count_data)), "\n"))
}

# Function for expression filtering. Only for the size factor estimation step
expression_filter <- function(count_data) {
    dataset_name <- deparse(substitute(count_data))
    n_before <- nrow(count_data)
    cat(paste(dataset_name, "#genes total:", n_before,"\n"))
    keep <- rowSums(count_data >= 2) >= 90
    count_data_filt <- count_data[keep,]
    n_after <- nrow(count_data_filt)
    cat(paste(dataset_name, "#genes after expression filter:", n_after,"\n"))
    return(count_data_filt)
}

# Function to estimate normalization size factors
estimate_size_factors <- function(count_data, meta_data) {
    # Create DESEq dataset
    dds <- DESeqDataSetFromMatrix(countData = count_data, 
        colData = meta_data, design = ~ 1)
    # Apply the median-of-ratios normalization
    dds <- estimateSizeFactors(dds)
    # Obtain normalized counts
    size_factors <- sizeFactors(dds)

    return(size_factors)
}

# Function to normalize the counts with pre-calculated size factors
normalize_counts <- function(count_data, meta_data, size_factors) {
    # Create DESEq dataset
    dds <- DESeqDataSetFromMatrix(countData = count_data, 
        colData = meta_data, design = ~ 1)
    # Supply the pre-calculated size factors
    sizeFactors(dds) <- size_factors
    # Obtain normalized counts
    normalized_counts <- counts(dds, normalized=TRUE)

    return(normalized_counts)
}

# ---- Species normalized separately ---- #

# Check that the row names of the meta data match column names of the count data
check_input_format(mm_count_data, mm_meta_data)
check_input_format(ac_count_data, ac_meta_data)

# 1. Filter the count data for expressed genes (for size-factor calculation)
ac_count_data_filt <- expression_filter(ac_count_data)
mm_count_data_filt <- expression_filter(mm_count_data) 

# 2. Obtain the size factors from the filtered data 
mm_sizefactors <- estimate_size_factors(mm_count_data_filt, mm_meta_data)
ac_sizefactors <- estimate_size_factors(ac_count_data_filt, ac_meta_data)

# 3. Normalize the unfiltered counts with the calculated size factors
mm_normalized_counts <- normalize_counts(mm_count_data, mm_meta_data, mm_sizefactors)
ac_normalized_counts <- normalize_counts(ac_count_data, ac_meta_data, ac_sizefactors)

saveRDS(ac_normalized_counts, paste0(out_dir_path, "AC_Tomoseq_counts_1-1orth_normalized.RData"))
saveRDS(mm_normalized_counts, paste0(out_dir_path, "MM_Tomoseq_counts_1-1orth_normalized.RData"))

