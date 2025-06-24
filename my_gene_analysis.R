# 1. Install BiocManager (if you don't have it already):
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
  print("BiocManager installed. Please proceed with GEOquery and limma installation.")
}

# 2. Install GEOquery, limma, dplyr, and ggrepel:
#    UNCOMMENT THE NEXT LINES AND RUN THEM IN YOUR R CONSOLE:
# BiocManager::install("GEOquery")
# BiocManager::install("limma")     # Essential for generating p-values
# install.packages("dplyr")         # For data manipulation, especially the pipe operator %>%
# BiocManager::install("ggrepel")   # Needed for labeling points on the volcano plot

# 3. Load the necessary libraries:
#    These lines should run without errors after successful installation.
library(GEOquery)
library(ggplot2)     # For visualization
library(dplyr)       # For data manipulation, especially the pipe operator %>% and mutate/case_when
library(limma)       # For differential expression analysis (generating p-values)
library(ggrepel)     # For non-overlapping text labels on plots (for Volcano Plot)

# ==============================================================================
# END OF PACKAGE INSTALLATION/LOADING INSTRUCTIONS.
# You can now run the rest of the script.
# ==============================================================================


# Define the GEO accession ID for the dataset we want to analyze
# GSE33146: Gene expression data for EMT in human mammary epithelial cells (HMLE)
geo_accession_id <- "GSE33146"

# Fetch the GEO dataset with error handling
print(paste("Attempting to fetch GEO dataset:", geo_accession_id, "... This may take a moment."))
gse <- tryCatch({
  getGEO(geo_accession_id, GSEMatrix = TRUE, AnnotGPL = FALSE)
}, error = function(e) {
  message(paste("ERROR: Failed to fetch GEO dataset:", e$message,
                "\n  Possible issues:",
                "\n  1. Internet connection problem.",
                "\n  2. Incorrect GEO accession ID ('", geo_accession_id, "').",
                "\n  3. 'GEOquery' package not properly installed or loaded.", sep=""))
  return(NULL) # Return NULL if an error occurs during fetching
})

# Crucial check: If 'gse' is NULL, it means data fetching failed. Terminate script.
if (is.null(gse)) {
  stop("Script terminated: Data acquisition failed. Please review error messages above and ensure internet connection/package installation.")
}

# 'getGEO' can return a list of ExpressionSet objects if a GSE ID contains data from multiple platforms.
# For simplicity, we'll take the first ExpressionSet object, which usually contains the main data.
if (is.list(gse) && length(gse) > 1) {
  print(paste("INFO: Multiple platforms found for GSE ID. Using the first ExpressionSet (gse[[1]]). Total platforms:", length(gse)))
  gse <- gse[[1]] # Extract the first ExpressionSet
} else if (is.list(gse) && length(gse) == 1) {
  gse <- gse[[1]] # If it's a list with one element, unlist it
}

# Verify that 'gse' is indeed an ExpressionSet object before proceeding
if (!inherits(gse, "ExpressionSet")) {
  stop(paste("ERROR: Fetched data is not a valid ExpressionSet object. Its class is:", class(gse),
             "\n  This might indicate an issue with the GEO dataset structure or GEOquery's processing."))
}
print("SUCCESS: Dataset fetched and confirmed as a valid ExpressionSet object.")
print(paste("Class of 'gse' object:", class(gse)))


# Extract the main components:
expression_data <- exprs(gse) # Matrix of gene expression values
pheno_data <- pData(gse)     # Data frame containing metadata/phenotypic information

# Display the dimensions of the expression data
print(paste("Expression data dimensions (genes x samples):", dim(expression_data)[1], "genes x", dim(expression_data)[2], "samples"))

# Display the first few rows and columns of the expression data matrix
print("--- Raw Expression Data (first 5 genes, first 5 samples) ---")
print(head(expression_data[, 1:min(5, ncol(expression_data)), drop = FALSE]))


# Display all column names from the phenotypic data to identify the condition column.
print("--- All Phenotypic Data Column Names (look for 'treatment', 'condition', 'group', 'medium', etc.) ---")
print(colnames(pheno_data))

# --- IMPORTANT: Identify the condition column and its values for Control and Treated ---
# For GSE33146, the 'characteristics_ch1.1' column contains the 'growth medium' information.
print("--- Inspecting Phenotypic Data for Condition Column and Unique Values ---")
condition_column <- NULL
if ("characteristics_ch1.1" %in% colnames(pheno_data)) {
  condition_column <- "characteristics_ch1.1"
  print(paste("Identified 'characteristics_ch1.1' as potential condition column."))
  print("Unique values in 'characteristics_ch1.1':")
  print(unique(pheno_data[[condition_column]]))
} else {
  condition_candidates <- grep("treatment|condition|group|medium|cell line|genotype", colnames(pheno_data), ignore.case = TRUE, value = TRUE)
  if (length(condition_candidates) > 0) {
    condition_column <- condition_candidates[1]
    print(paste("INFO: 'characteristics_ch1.1' not found. Using first general potential match:", condition_column))
    print(paste("Unique values in", condition_column, ":"))
    print(unique(pheno_data[[condition_column]]))
  } else {
    stop(paste("ERROR: No obvious condition column found in phenotypic data.",
               "\n  Please manually inspect 'pheno_data' (e.g., by running 'View(pData(gse))' in RStudio) ",
               "to find the column that specifies your control/treated conditions. ",
               "Then, modify the 'condition_column' assignment in the script (e.g., 'condition_column <- \"Your_Actual_Column_Name\"')."))
  }
}

# --- Part 2: Data Preparation and Limma Analysis for P-values ---

# Define control and treated samples based on the identified condition column and their specific values.
# For GSE33146, 'growth medium: MEGM' corresponds to control, and 'growth medium: SCGM' to treated.
control_label <- "MEGM"
treated_label <- "SCGM"

control_samples_indices <- which(grepl(control_label, pheno_data[[condition_column]], ignore.case = TRUE))
treated_samples_indices <- which(grepl(treated_label, pheno_data[[condition_column]], ignore.case = TRUE))

if (length(control_samples_indices) == 0) {
  stop(paste("ERROR: No '", control_label, "' samples found in column '", condition_column, "'.",
             "\n  Please verify the exact control condition value in your 'pheno_data' and adjust 'grepl(\"", control_label, "\", ...)' accordingly.", sep=""))
}
if (length(treated_samples_indices) == 0) {
  stop(paste("ERROR: No '", treated_label, "' samples found in column '", condition_column, "'.",
             "\n  Please verify the exact treated condition value in your 'pheno_data' and adjust 'grepl(\"", treated_label, "\", ...)' accordingly.", sep=""))
}

print(paste("Number of control samples (", control_label, "):", length(control_samples_indices)))
print(paste("Number of treated samples (", treated_label, "):", length(treated_samples_indices)))

# Subset the expression data to get control and treated sample matrices
# Ensure sample order matches the design matrix for limma.
samples_to_analyze_indices <- c(control_samples_indices, treated_samples_indices)
expression_data_subset <- expression_data[, samples_to_analyze_indices, drop = FALSE]

# Create a factor for experimental groups, ensuring levels are correctly ordered for contrasts
groups <- factor(c(rep(control_label, length(control_samples_indices)),
                   rep(treated_label, length(treated_samples_indices))),
                 levels = c(control_label, treated_label))

# Create a design matrix for limma
design <- model.matrix(~0 + groups)
colnames(design) <- levels(groups) # Rename columns to group names
print("--- Limma Design Matrix (first 6 rows) ---")
print(head(design))

# Fit the linear model to the expression data
fit <- lmFit(expression_data_subset, design)

# Define the contrast to compare Treated vs Control
# The 'coef' argument in topTable must match the name of the contrast produced here.
# For makeContrasts, the column name is usually derived directly from the 'levels' paste0 string.
contrast_name <- paste0(treated_label, "-", control_label) # This will be "SCGM-MEGM"
contrast_matrix <- makeContrasts(Contrasts = contrast_name, levels = design) # Use 'Contrasts' as generic name for this structure
print("--- Limma Contrast Matrix ---")
print(contrast_matrix)

# Apply contrasts to the fitted model
fit2 <- contrasts.fit(fit, contrast_matrix)

# Apply empirical Bayes moderation to estimate gene-wise variances
fit2 <- eBayes(fit2)

# Extract differential expression results
# CRITICAL FIX: Use the actual contrast name (e.g., "SCGM-MEGM") for the 'coef' argument.
limma_results <- topTable(fit2, coef = contrast_name, number = Inf, adjust.method = "fdr")

# Prepare the final analysis results data frame
analysis_results <- limma_results %>%
  as_tibble(rownames = "Gene") %>% # Convert to tibble and keep rownames as 'Gene'
  rename(Log2_Fold_Change = logFC, p_value = P.Value, Adjusted_P_Value = adj.P.Val) %>%
  mutate(
    Avg_Expression = AveExpr # Average expression across all samples
  ) %>%
  select(Gene, Avg_Expression, Log2_Fold_Change, p_value, Adjusted_P_Value)


# Display the first few rows of the analysis results with real p-values
print("--- Analysis Results (first 10 genes with calculated values from Limma) ---")
print(head(analysis_results, 10))


# --- Part 3: Identifying Differentially Expressed Genes with Real P-values ---

# Define thresholds for significance
log2fc_threshold <- 1      # Fold change threshold (2-fold change)
p_value_threshold <- 0.05 # Adjusted P-value threshold (False Discovery Rate - FDR)

# Classify genes based on both Log2 Fold Change and Adjusted P-value for plotting
analysis_results <- analysis_results %>%
  mutate(
    Significance = case_when(
      Log2_Fold_Change > log2fc_threshold & Adjusted_P_Value < p_value_threshold ~ "Upregulated",
      Log2_Fold_Change < -log2fc_threshold & Adjusted_P_Value < p_value_threshold ~ "Downregulated",
      TRUE ~ "Not Significant"
    )
  )

# Identify significantly upregulated and downregulated genes
upregulated_genes <- analysis_results %>%
  filter(Log2_Fold_Change > log2fc_threshold, Adjusted_P_Value < p_value_threshold)

downregulated_genes <- analysis_results %>%
  filter(Log2_Fold_Change < -log2fc_threshold, Adjusted_P_Value < p_value_threshold)


print(paste("--- Genes with |Log2 Fold Change| >", log2fc_threshold, "AND Adjusted P-value <", p_value_threshold, "---"))

# Print information about upregulated genes
if (nrow(upregulated_genes) > 0) {
  print("Upregulated Genes (Top 10 by Log2 Fold Change):")
  upregulated_genes_ordered <- upregulated_genes %>% arrange(desc(Log2_Fold_Change))
  print(head(upregulated_genes_ordered[, c("Gene", "Log2_Fold_Change", "p_value", "Adjusted_P_Value")], 10))
} else {
  print("No genes found significantly upregulated with the current thresholds.")
}

# Print information about downregulated genes
if (nrow(downregulated_genes) > 0) {
  print("Downregulated Genes (Top 10 by Log2 Fold Change):")
  downregulated_genes_ordered <- downregulated_genes %>% arrange(Log2_Fold_Change) # Ascending for downregulated
  print(head(downregulated_genes_ordered[, c("Gene", "Log2_Fold_Change", "p_value", "Adjusted_P_Value")], 10))
} else {
  print("No genes found significantly downregulated with the current thresholds.")
}


# --- Part 4: Simple Visualization (Top N Bar Plot) ---

# Prepare data for plotting: Select the top N most changed genes (either up or down).
plot_data_bar <- analysis_results %>%
  arrange(desc(abs(Log2_Fold_Change))) # Order by absolute fold change
num_genes_to_plot <- min(50, nrow(plot_data_bar))
plot_data_subset_bar <- head(plot_data_bar, num_genes_to_plot)

if (nrow(plot_data_subset_bar) == 0) {
  print("WARNING: Not enough genes with significant changes to create a bar plot with the current filters. Consider adjusting 'log2fc_threshold' or 'num_genes_to_plot'.")
} else {
  p_bar <- ggplot(plot_data_subset_bar, aes(x = reorder(Gene, Log2_Fold_Change), y = Log2_Fold_Change)) +
    geom_bar(stat = "identity", aes(fill = Significance)) +
    geom_hline(yintercept = log2fc_threshold, linetype = "dashed", color = "red", linewidth = 0.8) +
    geom_hline(yintercept = -log2fc_threshold, linetype = "dashed", color = "red", linewidth = 0.8) +
    scale_fill_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey70")) +
    coord_flip() +
    labs(
      title = paste0("Top Genes: Log2 Fold Change (GEO ID: ", geo_accession_id, ")"),
      x = "Gene (Probe ID)",
      y = "Log2 Fold Change (Treated / Control)",
      fill = "Differential Expression"
    ) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 7),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "bottom"
    )
  
  print(p_bar)
  ggsave("real_gene_expression_barplot.png", plot = p_bar, width = 10, height = 8, dpi = 300)
  print("Bar plot generated and SAVED to 'real_gene_expression_barplot.png'.")
}


# --- Part 5: Advanced Visualization (Volcano Plot) ---

# Check if there's enough data for a meaningful volcano plot
if (nrow(analysis_results) < 2) {
  print("WARNING: Not enough data points to create a meaningful Volcano Plot.")
} else {
  # Create a Volcano Plot using Adjusted_P_Value for Y-axis
  p_volcano <- ggplot(analysis_results, aes(x = Log2_Fold_Change, y = -log10(Adjusted_P_Value))) +
    geom_point(aes(color = Significance), alpha = 0.6, size = 2) +
    scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey70")) +
    
    # Add lines for significance thresholds
    geom_hline(yintercept = -log10(p_value_threshold), linetype = "dashed", color = "black", linewidth = 0.8) +
    geom_vline(xintercept = c(-log2fc_threshold, log2fc_threshold), linetype = "dashed", color = "black", linewidth = 0.8) +
    
    # Label top significant genes for clarity
    geom_text_repel(
      data = analysis_results %>% filter(Significance != "Not Significant") %>% 
        arrange(desc(abs(Log2_Fold_Change))) %>% head(15),
      aes(label = Gene),
      box.padding = 0.5, point.padding = 0.5,
      segment.color = 'grey50', max.overlaps = 20
    ) +
    
    labs(
      title = paste0("Volcano Plot (GEO ID: ", geo_accession_id, ")"),
      x = "Log2 Fold Change (Treated vs. Control)",
      y = "-log10(Adjusted p-value)", # Updated label
      color = "Gene Significance"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    )
  
  print(p_volcano)
  ggsave("real_gene_expression_volcanoplot.png", plot = p_volcano, width = 10, height = 8, dpi = 300)
  print("Volcano plot generated and SAVED to 'real_gene_expression_volcanoplot.png'.")
}


print("--- Real Data Analysis Complete! Please check your R console for output and plots. ---")

# Print the working directory so you know where the files are saved.
print("Your plot files are saved in your R working directory, which is:")
print(getwd())