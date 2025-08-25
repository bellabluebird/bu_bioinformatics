# R/fct_data_processing.R
# data loading and processing functions for GSE64810

#' parse GEO Series Matrix File
#' 
#' @param file_path Path to the series matrix file
#' @return Data frame with parsed sample metadata
parse_geo_series_matrix <- function(file_path) {
  # read the entire file
  lines <- readLines(file_path)
  
  # helper function to extract values from a metadata row
  extract_values <- function(field_name) {
    row <- lines[grepl(paste0("^", field_name), lines)]
    if (length(row) == 0) return(character(0))
    values <- unlist(strsplit(row, "\t"))[-1]
    gsub('"', '', values)
  }
  
  # extract basic metadata
  sample_titles <- extract_values("!Sample_title")
  geo_accessions <- extract_values("!Sample_geo_accession")
  organisms <- extract_values("!Sample_organism_ch1")
  platforms <- extract_values("!Sample_platform_id")
  instruments <- extract_values("!Sample_instrument_model")
  
  # get characteristics data
  characteristics_lines <- lines[grepl("^!Sample_characteristics_ch1", lines)]
  characteristics <- list()
  
  for (line in characteristics_lines) {
    values <- unlist(strsplit(line, "\t"))[-1]
    values <- gsub('"', '', values)
    
    first_value <- values[values != "" & !is.na(values)][1]
    if (!is.na(first_value) && grepl(":", first_value)) {
      field_name <- gsub(":.*", "", first_value)
      field_values <- sapply(values, function(x) {
        if (x == "" || is.na(x)) return("")
        if (grepl(":", x)) return(gsub(".*: ", "", x))
        return(x)
      })
      characteristics[[field_name]] <- field_values
    }
  }
  
  # create metadata dataframe
  n_samples <- length(sample_titles)
  metadata_df <- data.frame(
    sample_id = sample_titles,
    geo_accession = if(length(geo_accessions) == n_samples) geo_accessions else rep(NA, n_samples),
    organism = if(length(organisms) == n_samples) organisms else rep(NA, n_samples),
    platform = if(length(platforms) == n_samples) platforms else rep(NA, n_samples),
    instrument = if(length(instruments) == n_samples) instruments else rep(NA, n_samples),
    stringsAsFactors = FALSE
  )
  
  # add characteristics with standardized column names
  for (field in names(characteristics)) {
    col_name <- switch(field,
                       "tissue" = "tissue",
                       "diagnosis" = "diagnosis",
                       "pmi" = "pmi",
                       "age of death" = "age_of_death",
                       "rin" = "rin",
                       "mrna-seq reads" = "mrna_seq_reads",
                       "age of onset" = "age_of_onset", 
                       "duration" = "duration",
                       "cag" = "cag_repeats",
                       "vonsattel grade" = "vonsattel_grade",
                       "h-v striatal score" = "hv_striatal_score",
                       "h-v cortical score" = "hv_cortical_score",
                       gsub(" ", "_", field))
    
    if (field %in% c("pmi", "age of death", "rin", "mrna-seq reads", "age of onset", 
                     "duration", "cag", "vonsattel grade", "h-v striatal score", "h-v cortical score")) {
      values <- characteristics[[field]]
      values[values == "" | values == "NA"] <- NA
      metadata_df[[col_name]] <- as.numeric(values)
    } else {
      metadata_df[[col_name]] <- characteristics[[field]]
    }
  }
  
  # create group variable and standardize
  if ("diagnosis" %in% names(metadata_df)) {
    metadata_df$group <- ifelse(grepl("Huntington", metadata_df$diagnosis, ignore.case = TRUE), "HD", "Control")
  }
  
  # convert reads to millions for better visualization
  if ("mrna_seq_reads" %in% names(metadata_df)) {
    metadata_df$reads_millions <- metadata_df$mrna_seq_reads / 1e6
  }
  
  return(metadata_df)
}

#' Validate Sample Metadata
#' 
#' @param metadata_df Parsed metadata dataframe
#' @return list with validation results
validate_sample_metadata <- function(metadata_df) {
  validation <- list(
    valid = TRUE,
    messages = character(0),
    warnings = character(0)
  )
  
  # check required columns
  required_cols <- c("sample_id", "diagnosis")
  missing_cols <- setdiff(required_cols, names(metadata_df))
  
  if (length(missing_cols) > 0) {
    validation$valid <- FALSE
    validation$messages <- c(validation$messages, 
                             paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
  }
  
  # check sample counts
  if ("group" %in% names(metadata_df)) {
    group_counts <- table(metadata_df$group, useNA = "ifany")
    validation$messages <- c(validation$messages,
                             paste("Sample counts:", paste(names(group_counts), group_counts, sep = ": ", collapse = ", ")))
  }
  
  # check for hd-specific data
  hd_cols <- c("age_of_onset", "cag_repeats", "duration", "vonsattel_grade")
  available_hd_cols <- intersect(hd_cols, names(metadata_df))
  
  if (length(available_hd_cols) > 0) {
    non_na_counts <- sapply(available_hd_cols, function(col) sum(!is.na(metadata_df[[col]])))
    validation$messages <- c(validation$messages,
                             paste("HD-specific data available for:", 
                                   paste(names(non_na_counts)[non_na_counts > 0], collapse = ", ")))
  }
  
  return(validation)
}

#' Get Sample Metadata Summary Statistics
#' 
#' @param metadata_df Parsed metadata dataframe
#' @return list with summary statistics
get_metadata_summary <- function(metadata_df) {
  library(dplyr)
  
  summary_stats <- list()
  
  # basic counts
  summary_stats$total_samples <- nrow(metadata_df)
  
  if ("group" %in% names(metadata_df)) {
    group_counts <- table(metadata_df$group, useNA = "ifany")
    summary_stats$hd_samples <- ifelse("HD" %in% names(group_counts), group_counts[["HD"]], 0)
    summary_stats$control_samples <- ifelse("Control" %in% names(group_counts), group_counts[["Control"]], 0)
  }
  
  # age statistics
  age_col <- intersect(c("age_of_death", "age_death"), names(metadata_df))[1]
  if (!is.na(age_col) && age_col %in% names(metadata_df)) {
    if ("group" %in% names(metadata_df)) {
      age_stats <- metadata_df %>%
        group_by(group) %>%
        summarise(
          mean_age = round(mean(.data[[age_col]], na.rm = TRUE), 1),
          median_age = round(median(.data[[age_col]], na.rm = TRUE), 1),
          sd_age = round(sd(.data[[age_col]], na.rm = TRUE), 1),
          .groups = 'drop'
        )
      summary_stats$age_stats <- age_stats
    }
  }
  
  # hd-specific statistics
  if ("group" %in% names(metadata_df)) {
    hd_data <- metadata_df %>% filter(group == "HD")
    if (nrow(hd_data) > 0) {
      hd_summary <- list()
      
      if ("age_of_onset" %in% names(hd_data)) {
        hd_summary$mean_onset <- round(mean(hd_data$age_of_onset, na.rm = TRUE), 1)
      }
      if ("duration" %in% names(hd_data)) {
        hd_summary$mean_duration <- round(mean(hd_data$duration, na.rm = TRUE), 1)
      }
      if ("cag_repeats" %in% names(hd_data)) {
        hd_summary$mean_cag <- round(mean(hd_data$cag_repeats, na.rm = TRUE), 1)
        hd_summary$cag_range <- paste(range(hd_data$cag_repeats, na.rm = TRUE), collapse = "-")
      }
      
      summary_stats$hd_specific <- hd_summary
    }
  }
  
  return(summary_stats)
}

#' Calculate summary statistics for numeric variables by group
#' 
#' @param data Data frame with samples
#' @param group_var Column name for grouping variable
#' @param numeric_vars Vector of numeric column names
#' @return Data frame with summary statistics
calculate_group_summaries <- function(data, group_var = "group", numeric_vars = NULL) {
  library(dplyr)
  
  if (is.null(numeric_vars)) {
    numeric_vars <- names(data)[sapply(data, is.numeric)]
  }
  
  data %>%
    group_by(across(all_of(group_var))) %>%
    summarise(
      across(all_of(numeric_vars), 
             list(mean = ~mean(.x, na.rm = TRUE),
                  median = ~median(.x, na.rm = TRUE),
                  sd = ~sd(.x, na.rm = TRUE),
                  n = ~sum(!is.na(.x))),
             .names = "{.col}_{.fn}"),
      .groups = "drop"
    )
}

#' Validate uploaded file format
#' 
#' @param file_path Path to uploaded file
#' @param expected_type Type of file expected ("series_matrix", "deseq_results", "counts")
#' @return List with validation results
validate_file_format <- function(file_path, expected_type = "series_matrix") {
  
  result <- list(
    valid = FALSE,
    message = "",
    warnings = character(0)
  )
  
  tryCatch({
    # read first few lines
    lines <- readLines(file_path, n = 20)
    
    if (expected_type == "series_matrix") {
      # check for geo series matrix format
      if (any(grepl("^!Series_title", lines))) {
        result$valid <- TRUE
        result$message <- "Valid GEO series matrix file detected"
      } else {
        result$message <- "File does not appear to be a GEO series matrix format"
      }
    } else if (expected_type == "deseq_results") {
      # check for deseq2 results format
      header <- lines[1]
      if (grepl("baseMean.*log2FoldChange.*pvalue", header)) {
        result$valid <- TRUE
        result$message <- "Valid DESeq2 results file detected"
      } else {
        result$message <- "File does not appear to be DESeq2 results format"
      }
    } else if (expected_type == "counts") {
      # check for count matrix format
      header <- lines[1]
      if (grepl("^\\w+\\t.*GSM", header)) {
        result$valid <- TRUE
        result$message <- "Valid count matrix file detected"
      } else {
        result$message <- "File does not appear to be a count matrix format"
      }
    }
    
  }, error = function(e) {
    result$message <- paste("Error reading file:", e$message)
  })
  
  return(result)
}

#' Load DESeq2 Results
#' 
#' @param file_path Path to DESeq2 results file
#' @return Data frame with processed DESeq2 results
load_deseq_results <- function(file_path) {
  library(readr)
  library(dplyr)
  
  cat("Loading DESeq2 results from:", file_path, "\n")
  
  # read gse64810 deseq2 results
  deseq_results <- read_tsv(file_path, show_col_types = FALSE) %>%
    # handle potential gene id column naming variations
    rename_with(~case_when(
      . %in% c("X1", "gene_id", "Gene", "Geneid", "symbol") ~ "gene_id",
      . %in% c("baseMean", "AveExpr") ~ "baseMean", 
      . %in% c("log2FoldChange", "logFC") ~ "log2FoldChange",
      . %in% c("padj", "adj.P.Val") ~ "padj",
      . %in% c("pvalue", "P.Value") ~ "pvalue",
      TRUE ~ .
    )) %>%
    # add derived columns for visualization
    mutate(
      neg_log10_pval = -log10(pvalue),
      abs_log2fc = abs(log2FoldChange),
      significance = case_when(
        padj < 0.001 & abs_log2fc > 2 ~ "High significance",
        padj < 0.01 & abs_log2fc > 1 ~ "Medium significance", 
        padj < 0.05 & abs_log2fc > 0.5 ~ "Low significance",
        TRUE ~ "Not significant"
      ),
      regulation = case_when(
        log2FoldChange > 0 & padj < 0.05 ~ "Up-regulated",
        log2FoldChange < 0 & padj < 0.05 ~ "Down-regulated", 
        TRUE ~ "No change"
      )
    ) %>%
    # remove genes with missing critical values
    filter(!is.na(pvalue), !is.na(padj), !is.na(log2FoldChange))
  
  cat("✓ Loaded", nrow(deseq_results), "genes successfully\n")
  
  return(deseq_results)
}

#' Load Normalized Counts
#' 
#' @param file_path Path to normalized counts file
#' @param max_genes Maximum number of genes to load (for performance)
#' @return Matrix with normalized counts
load_normalized_counts <- function(file_path, max_genes = 10000) {
  library(data.table)
  library(dplyr)
  
  cat("Loading normalized counts from:", file_path, "\n")
  
  # use data.table for efficient reading of large count matrices
  counts <- fread(file_path, data.table = FALSE)
  
  # handle gene id column naming (first column is typically gene names)
  gene_names <- counts[[1]]
  counts_matrix <- as.matrix(counts[, -1])
  rownames(counts_matrix) <- gene_names
  
  # remove genes with all zero counts
  keep_genes <- rowSums(counts_matrix > 0) > 0
  counts_matrix <- counts_matrix[keep_genes, ]
  
  # for prototyping, sample genes if too large
  if (nrow(counts_matrix) > max_genes) {
    # keep top variable genes
    gene_vars <- apply(counts_matrix, 1, var, na.rm = TRUE)
    keep_genes <- names(sort(gene_vars, decreasing = TRUE)[1:max_genes])
    counts_matrix <- counts_matrix[keep_genes, ]
    cat("Subsampled to top", max_genes, "most variable genes\n")
  }
  
  cat("✓ Loaded", nrow(counts_matrix), "genes x", ncol(counts_matrix), "samples\n")
  
  return(counts_matrix)
}

#' Prepare Heatmap Data
#' 
#' @param counts_matrix Matrix of normalized counts
#' @param deseq_results DESeq2 results data frame
#' @param selected_genes Vector of selected gene names
#' @param n_top_genes Number of top DE genes to use if no genes selected
#' @return Matrix prepared for heatmap visualization
prepare_heatmap_data <- function(counts_matrix, deseq_results, 
                                 selected_genes = NULL, n_top_genes = 50) {
  library(dplyr)
  
  # if specific genes selected, use those; otherwise use top de genes
  if (!is.null(selected_genes) && length(selected_genes) > 0) {
    heatmap_genes <- intersect(selected_genes, rownames(counts_matrix))
    cat("Using", length(heatmap_genes), "selected genes for heatmap\n")
  } else {
    # get top differentially expressed genes
    top_genes <- deseq_results %>%
      filter(!is.na(padj), padj < 0.05) %>%
      arrange(padj) %>%
      slice_head(n = n_top_genes) %>%
      pull(gene_id)
    
    heatmap_genes <- intersect(top_genes, rownames(counts_matrix))
    cat("Using top", length(heatmap_genes), "DE genes for heatmap\n")
  }
  
  if (length(heatmap_genes) == 0) {
    stop("No genes found for heatmap visualization")
  }
  
  # subset count matrix
  heatmap_matrix <- counts_matrix[heatmap_genes, , drop = FALSE]
  
  # log2 transform and scale for visualization
  heatmap_matrix <- log2(heatmap_matrix + 1)
  
  # z-score scaling by gene (row)
  heatmap_matrix <- t(scale(t(heatmap_matrix)))
  
  return(heatmap_matrix)
}

#' Create Diagnostic Plots
#' 
#' @param counts_matrix Matrix of normalized counts
#' @return List of diagnostic plots
create_diagnostic_plots <- function(counts_matrix) {
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  
  # convert to data frame for ggplot
  counts_df <- as.data.frame(counts_matrix)
  
  # 1. library size plot
  library_sizes <- colSums(counts_matrix)
  lib_size_plot <- data.frame(
    sample = names(library_sizes),
    library_size = library_sizes
  ) %>%
    ggplot(aes(x = reorder(sample, library_size), y = library_size)) +
    geom_col(fill = "steelblue", alpha = 0.7) +
    coord_flip() +
    labs(title = "Library Sizes per Sample", 
         x = "Sample", y = "Total Counts") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 8))
  
  # 2. count distribution plot
  count_dist_data <- counts_df %>%
    pivot_longer(everything(), names_to = "sample", values_to = "counts") %>%
    filter(counts > 0) %>%
    mutate(log_counts = log2(counts + 1))
  
  count_dist_plot <- count_dist_data %>%
    ggplot(aes(x = log_counts)) +
    geom_density(alpha = 0.6, fill = "lightblue") +
    facet_wrap(~sample, scales = "free", ncol = 4) +
    labs(title = "Count Distribution per Sample (log2)", 
         x = "log2(counts + 1)", y = "Density") +
    theme_minimal() +
    theme(axis.text = element_text(size = 6),
          strip.text = element_text(size = 7))
  
  # 3. sample correlation heatmap data
  sample_cors <- cor(counts_matrix, use = "complete.obs")
  
  # 4. mean vs variance plot
  gene_means <- rowMeans(counts_matrix)
  gene_vars <- apply(counts_matrix, 1, var)
  
  mean_var_plot <- data.frame(
    mean = gene_means,
    variance = gene_vars
  ) %>%
    filter(mean > 0, variance > 0) %>%
    ggplot(aes(x = log10(mean), y = log10(variance))) +
    geom_point(alpha = 0.5, size = 0.5) +
    geom_smooth(method = "lm", color = "red") +
    labs(title = "Mean-Variance Relationship", 
         x = "log10(Mean counts)", y = "log10(Variance)") +
    theme_minimal()
  
  return(list(
    library_size = lib_size_plot,
    count_distribution = count_dist_plot,
    sample_correlation = sample_cors,
    mean_variance = mean_var_plot
  ))
}

#' Perform PCA Analysis
#' 
#' @param counts_matrix Matrix of normalized counts
#' @param n_genes Number of most variable genes to use for PCA
#' @return List with PCA results
perform_pca_analysis <- function(counts_matrix, n_genes = 1000) {
  library(dplyr)
  
  # select most variable genes for pca
  gene_vars <- apply(counts_matrix, 1, var, na.rm = TRUE)
  top_var_genes <- names(sort(gene_vars, decreasing = TRUE)[1:min(n_genes, length(gene_vars))])
  
  # subset and log transform
  pca_matrix <- counts_matrix[top_var_genes, ]
  pca_matrix <- log2(pca_matrix + 1)
  
  # remove genes with zero variance
  gene_vars_subset <- apply(pca_matrix, 1, var, na.rm = TRUE)
  pca_matrix <- pca_matrix[gene_vars_subset > 0, ]
  
  # perform pca on transposed matrix (samples as rows)
  pca_result <- prcomp(t(pca_matrix), scale. = TRUE, center = TRUE)
  
  # calculate variance explained
  var_explained <- (pca_result$sdev^2 / sum(pca_result$sdev^2)) * 100
  
  # create pca data frame
  pca_df <- data.frame(
    sample = rownames(pca_result$x),
    PC1 = pca_result$x[, 1],
    PC2 = pca_result$x[, 2],
    PC3 = pca_result$x[, 3]
  )
  
  # extract loadings (top contributing genes)
  loadings <- pca_result$rotation[, 1:3]
  top_loadings_pc1 <- loadings[order(abs(loadings[, 1]), decreasing = TRUE)[1:10], 1]
  top_loadings_pc2 <- loadings[order(abs(loadings[, 2]), decreasing = TRUE)[1:10], 2]
  
  return(list(
    pca_data = pca_df,
    variance_explained = var_explained,
    pca_object = pca_result,
    top_genes_pc1 = top_loadings_pc1,
    top_genes_pc2 = top_loadings_pc2,
    n_genes_used = nrow(pca_matrix)
  ))
}

#' Create PCA Plots
#' 
#' @param pca_analysis PCA analysis results from perform_pca_analysis
#' @return List of PCA plots
create_pca_plots <- function(pca_analysis) {
  library(ggplot2)
  library(plotly)
  
  pca_data <- pca_analysis$pca_data
  var_exp <- pca_analysis$variance_explained
  
  # 2d pca plot
  pca_2d <- pca_data %>%
    ggplot(aes(x = PC1, y = PC2, text = paste("Sample:", sample))) +
    geom_point(size = 3, alpha = 0.7, color = "steelblue") +
    labs(title = "PCA Analysis of Gene Expression",
         x = paste0("PC1 (", round(var_exp[1], 1), "% variance)"),
         y = paste0("PC2 (", round(var_exp[2], 1), "% variance)")) +
    theme_minimal() +
    theme(aspect.ratio = 1)
  
  # variance explained plot (scree plot)
  var_plot <- data.frame(
    PC = paste0("PC", 1:min(10, length(var_exp))),
    Variance = var_exp[1:min(10, length(var_exp))]
  ) %>%
    ggplot(aes(x = factor(PC, levels = PC), y = Variance)) +
    geom_col(fill = "lightcoral", alpha = 0.7) +
    geom_text(aes(label = paste0(round(Variance, 1), "%")), 
              vjust = -0.5, size = 3) +
    labs(title = "Variance Explained by Principal Components",
         x = "Principal Component", y = "% Variance Explained") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(list(
    pca_2d = pca_2d,
    variance_plot = var_plot
  ))
}
