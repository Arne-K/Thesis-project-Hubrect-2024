library(DESeq2)
library(ggplot2)
library(dplyr)

input_dir <- "/home/a.kukkonen_cbs-niob.local/PartB/deseq2_input_data/tomoseq_datasets/"
output_dir <- "/home/a.kukkonen_cbs-niob.local/PartB/deseq2_output/tomoseq_sham_MMvAC/"

# Load the combined unnormalized counts and meta data
count_data <- readRDS(paste0(input_dir, "MM-AC_Tomoseq_counts_1-1orth_unnorm.RData"))
meta_data <- readRDS(paste0(input_dir, "MM-AC_Tomoseq_meta_1-1orth_unnorm.RData"))

# Subset the MI14D condition
count_data <- count_data %>% select(matches("_SHAM."))
meta_data <- meta_data[colnames(count_data), ]

# Ensure species is of class factor
meta_data$Species <- as.factor(meta_data$Species)

###################
#### FUNCTIONS ####
###################

# Function for subsetting slices
slice_subset <- function(count_data, meta_data, slice_numbers) {
    result_list <- list()
    count_data_subset <- count_data %>% select(matches(paste0("\\.slice(", paste(slice_numbers, collapse="|"), ")$")))
    meta_data_subset <- meta_data[colnames(count_data_subset), ]
    result_list[["count_data_subset"]] <- count_data_subset
    result_list[["meta_data_subset"]] <- meta_data_subset
    return(result_list)
}

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

    cat(paste("Any missing values in ", dataset_name, ": ", anyNA(count_data), "\n"))
    cat(paste0("Samples in ", dataset_name, ": ", length(colnames(count_data)), "\n"))
}


# Function for expression filtering. Only for the size factor estimation step
expression_filter <- function(count_data, count_threshold, sample_n_threshold) {
    dataset_name <- deparse(substitute(count_data))
    n_before <- nrow(count_data)
    cat(paste(dataset_name, "#genes total:", n_before,"\n"))
    keep <- rowSums(count_data >= count_threshold) >= sample_n_threshold
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

check_input_format(count_data = count_data, meta_data = meta_data)

# Apply an expression filter for normalization step
counts_for_norm <- expression_filter(count_data = count_data, count_threshold = 2, sample_n_threshold = 80)
sample_size_factors <- estimate_size_factors(counts_for_norm, meta_data)

# Subset the slices for DE
slice_numbers <- c(5:20)
subset_for_de <- slice_subset(count_data = count_data, meta_data = meta_data, slice_numbers = slice_numbers)
counts_for_de <- subset_for_de[["count_data_subset"]]
meta_for_de <- subset_for_de[["meta_data_subset"]]
rm(subset_for_de)

check_input_format(count_data = counts_for_de, meta_data = meta_for_de)

# subset the size factors
sample_size_factors <- sample_size_factors[colnames(counts_for_de)]

# Apply expression filter for DE analysis step
counts_for_de <- expression_filter(counts_for_de, count_threshold = 10, sample_n_threshold = 25)

# save the unfiltered input counts 
saveRDS(counts_for_de, paste0(output_dir, "DESeq2_input_counts_unfilt_Tomoseq_sham_MMvAC.RData"))

# Run differential expression analysis with the size factors
dds <- DESeqDataSetFromMatrix(countData = counts_for_de, 
    colData = meta_for_de, design = ~ Species)
# Set reference levels for variable
dds$Species = relevel( dds$Species, "AC" )
# Display the levels for checking
dds$Species
# Supply the normalization size factors
sizeFactors(dds) <- sample_size_factors
# Perform the DE analysis
dds <- DESeq(dds)

###############
#### PLOTS ####
###############

# Plot dispersion estimates
png(file = paste0(output_dir, "DispEst_Tomoseq_sham_MMvAC.png"), width = 1200, height = 900)
plotDispEsts(dds)
dev.off()

# Variance-stabilizing transformation (VST)
#vsd <- vst(dds, blind = FALSE)
rld <- rlog(dds, blind = FALSE)

# Make pca plot to model distances between species
png(file = paste0(output_dir, "PCA_Tomoseq_sham_MMvAC_RLog.png"), width = 600, height = 350)
plotPCA(rld, intgroup = c("Species"))
dev.off()

#################
#### RESULTS ####
#################

# Display all the variable names
resultsNames(dds)
# Get results
res <- results(dds, name = "Species_MM_vs_AC")
write.table(res, paste0(output_dir, "DESeq2_output_Tomoseq_sham_MMvAC.tsv"), sep="\t", col.names=T, quote=F, append=F, row.names=T)




