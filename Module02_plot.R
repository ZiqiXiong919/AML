############################################################
##
## Module02B_QC_and_SaveObjects.R
##
## Project:
## Lysine metabolism-related genes predict prognosis in AML
##
## Purpose:
## Add QC analysis after Module02 / Module03 preprocessing
##
## Input:
## result/TCGA_expression.csv
## result/TCGA_survival.csv
## result/TCGA_clinical.csv
## result/GSE37642_expression.csv
## result/GSE37642_pheno.csv
## result/GSE12417_expression.csv
## result/GSE12417_pheno.csv
##
## Output:
## figure/QC_*.pdf
## result/Module02_QC_Report.csv
## object/Module02_Preprocessed_Data.RData
##
############################################################

rm(list = ls())
options(stringsAsFactors = FALSE)

cat("=============================================\n")
cat("Module02B: QC and Save Preprocessed Objects\n")
cat("=============================================\n\n")

############################################################
## 1. Packages
############################################################

cran_pkgs <- c(
  "data.table",
  "dplyr",
  "ggplot2",
  "reshape2",
  "pheatmap"
)

for(pkg in cran_pkgs){
  if(!requireNamespace(pkg, quietly = TRUE)){
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

############################################################
## 2. Create folders
############################################################

dir.create("figure", showWarnings = FALSE, recursive = TRUE)
dir.create("result", showWarnings = FALSE, recursive = TRUE)
dir.create("object", showWarnings = FALSE, recursive = TRUE)

############################################################
## 3. Helper functions
############################################################

read_expr_csv <- function(file){
  
  dat <- data.table::fread(
    file,
    data.table = FALSE,
    check.names = FALSE
  )
  
  rownames(dat) <- dat[, 1]
  dat <- dat[, -1, drop = FALSE]
  
  dat <- as.matrix(dat)
  mode(dat) <- "numeric"
  
  return(dat)
}

qc_summary <- function(expr, dataset_name){
  
  data.frame(
    Dataset = dataset_name,
    Genes = nrow(expr),
    Samples = ncol(expr),
    Missing_values = sum(is.na(expr)),
    Missing_rate = sum(is.na(expr)) / (nrow(expr) * ncol(expr)),
    Min = min(expr, na.rm = TRUE),
    Q1 = quantile(expr, 0.25, na.rm = TRUE),
    Median = median(expr, na.rm = TRUE),
    Mean = mean(expr, na.rm = TRUE),
    Q3 = quantile(expr, 0.75, na.rm = TRUE),
    Max = max(expr, na.rm = TRUE)
  )
}

plot_boxplot <- function(expr, dataset_name, outfile){
  
  plot_df <- data.frame(
    Expression = as.vector(expr),
    Sample = rep(colnames(expr), each = nrow(expr))
  )
  
  ## If too many samples, boxplot all but suppress x labels
  p <- ggplot(plot_df, aes(x = Sample, y = Expression)) +
    geom_boxplot(outlier.size = 0.2) +
    theme_bw() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.title = element_text(hjust = 0.5)
    ) +
    labs(
      title = paste0(dataset_name, " expression boxplot"),
      x = "Samples",
      y = "Expression"
    )
  
  ggsave(outfile, p, width = 10, height = 5)
}

plot_density <- function(expr, dataset_name, outfile){
  
  plot_df <- data.frame(
    Expression = as.vector(expr),
    Sample = rep(colnames(expr), each = nrow(expr))
  )
  
  p <- ggplot(plot_df, aes(x = Expression, group = Sample)) +
    geom_density(alpha = 0.25) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5)
    ) +
    labs(
      title = paste0(dataset_name, " expression density"),
      x = "Expression",
      y = "Density"
    )
  
  ggsave(outfile, p, width = 8, height = 5)
}

plot_pca <- function(expr, dataset_name, outfile){
  
  ## Remove genes with zero variance
  gene_sd <- apply(expr, 1, sd, na.rm = TRUE)
  expr2 <- expr[gene_sd > 0, , drop = FALSE]
  
  ## Use top variable genes to make PCA more stable
  gene_var <- apply(expr2, 1, var, na.rm = TRUE)
  
  top_n <- min(5000, length(gene_var))
  
  top_genes <- names(sort(gene_var, decreasing = TRUE))[1:top_n]
  
  pca_input <- t(expr2[top_genes, , drop = FALSE])
  
  pca <- prcomp(
    pca_input,
    center = TRUE,
    scale. = TRUE
  )
  
  pca_df <- data.frame(
    Sample = rownames(pca$x),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2]
  )
  
  pc1_var <- round(summary(pca)$importance[2, 1] * 100, 2)
  pc2_var <- round(summary(pca)$importance[2, 2] * 100, 2)
  
  p <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
    geom_point(size = 2) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5)
    ) +
    labs(
      title = paste0(dataset_name, " PCA"),
      x = paste0("PC1 (", pc1_var, "%)"),
      y = paste0("PC2 (", pc2_var, "%)")
    )
  
  ggsave(outfile, p, width = 6, height = 5)
  
  return(pca_df)
}

missing_report <- function(expr, dataset_name){
  
  sample_missing <- data.frame(
    Dataset = dataset_name,
    Sample = colnames(expr),
    Missing_values = colSums(is.na(expr)),
    Missing_rate = colSums(is.na(expr)) / nrow(expr)
  )
  
  gene_missing <- data.frame(
    Dataset = dataset_name,
    Gene = rownames(expr),
    Missing_values = rowSums(is.na(expr)),
    Missing_rate = rowSums(is.na(expr)) / ncol(expr)
  )
  
  return(
    list(
      sample_missing = sample_missing,
      gene_missing = gene_missing
    )
  )
}

############################################################
## 4. Read input files
############################################################

cat("Reading expression matrices...\n")

tcga_expr <- read_expr_csv(
  "C:/Users/Xiong/Desktop/AML/result/TCGA_expression.csv"
)

gse37642_expr <- read_expr_csv(
  "C:/Users/Xiong/Desktop/AML/result/GSE37642_expression.csv"
)

gse12417_expr <- read_expr_csv(
  "C:/Users/Xiong/Desktop/AML/result/GSE12417_expression.csv"
)

cat("Reading clinical and phenotype files...\n")

tcga_survival <- data.table::fread(
  "C:/Users/Xiong/Desktop/AML/result/TCGA_survival.csv",
  data.table = FALSE
)

tcga_clinical <- data.table::fread(
  "C:/Users/Xiong/Desktop/AML/result/TCGA_clinical.csv",
  data.table = FALSE
)

gse37642_pheno <- data.table::fread(
  "C:/Users/Xiong/Desktop/AML/result/GSE37642_pheno.csv",
  data.table = FALSE,
  check.names = FALSE
)

gse12417_pheno <- data.table::fread(
  "C:/Users/Xiong/Desktop/AML/result/GSE12417_pheno.csv",
  data.table = FALSE,
  check.names = FALSE
)

cat("Input data loaded.\n\n")

############################################################
## 5. Basic QC summary
############################################################

qc_tcga <- qc_summary(tcga_expr, "TCGA")
qc_gse37642 <- qc_summary(gse37642_expr, "GSE37642")
qc_gse12417 <- qc_summary(gse12417_expr, "GSE12417")

qc_report <- dplyr::bind_rows(
  qc_tcga,
  qc_gse37642,
  qc_gse12417
)

write.csv(
  qc_report,
  "C:/Users/Xiong/Desktop/AML/result/Module02_QC_Report.csv",
  row.names = FALSE
)

cat("QC summary saved:\n")
cat("result/Module02_QC_Report.csv\n\n")

print(qc_report)

############################################################
## 6. Missing value QC
############################################################

miss_tcga <- missing_report(tcga_expr, "TCGA")
miss_gse37642 <- missing_report(gse37642_expr, "GSE37642")
miss_gse12417 <- missing_report(gse12417_expr, "GSE12417")

sample_missing_report <- dplyr::bind_rows(
  miss_tcga$sample_missing,
  miss_gse37642$sample_missing,
  miss_gse12417$sample_missing
)

gene_missing_report <- dplyr::bind_rows(
  miss_tcga$gene_missing,
  miss_gse37642$gene_missing,
  miss_gse12417$gene_missing
)

write.csv(
  sample_missing_report,
  "C:/Users/Xiong/Desktop/AML/result/Module02_SampleMissing_Report.csv",
  row.names = FALSE
)

write.csv(
  gene_missing_report,
  "C:/Users/Xiong/Desktop/AML/result/Module02_GeneMissing_Report.csv",
  row.names = FALSE
)

cat("Missing value reports saved.\n\n")

############################################################
## 7. Boxplot QC
############################################################

cat("Drawing boxplot QC...\n")

plot_boxplot(
  tcga_expr,
  "TCGA",
  "C:/Users/Xiong/Desktop/AML/figure/QC_TCGA_Boxplot.pdf"
)

plot_boxplot(
  gse37642_expr,
  "GSE37642",
  "C:/Users/Xiong/Desktop/AML/figure/QC_GSE37642_Boxplot.pdf"
)

plot_boxplot(
  gse12417_expr,
  "GSE12417",
  "C:/Users/Xiong/Desktop/AML/figure/QC_GSE12417_Boxplot.pdf"
)

cat("Boxplot QC finished.\n\n")

############################################################
## 8. Density QC
############################################################

cat("Drawing density QC...\n")

plot_density(
  tcga_expr,
  "TCGA",
  "C:/Users/Xiong/Desktop/AML/figure/QC_TCGA_Density.pdf"
)

plot_density(
  gse37642_expr,
  "GSE37642",
  "C:/Users/Xiong/Desktop/AML/figure/QC_GSE37642_Density.pdf"
)

plot_density(
  gse12417_expr,
  "GSE12417",
  "C:/Users/Xiong/Desktop/AML/figure/QC_GSE12417_Density.pdf"
)

cat("Density QC finished.\n\n")

############################################################
## 9. PCA QC
############################################################

cat("Drawing PCA QC...\n")

pca_tcga <- plot_pca(
  tcga_expr,
  "TCGA",
  "C:/Users/Xiong/Desktop/AML/figure/QC_TCGA_PCA.pdf"
)

pca_gse37642 <- plot_pca(
  gse37642_expr,
  "GSE37642",
  "C:/Users/Xiong/Desktop/AML/figure/QC_GSE37642_PCA.pdf"
)

pca_gse12417 <- plot_pca(
  gse12417_expr,
  "GSE12417",
  "C:/Users/Xiong/Desktop/AML/figure/QC_GSE12417_PCA.pdf"
)

write.csv(
  pca_tcga,
  "C:/Users/Xiong/Desktop/AML/result/QC_TCGA_PCA_coordinates.csv",
  row.names = FALSE
)

write.csv(
  pca_gse37642,
  "C:/Users/Xiong/Desktop/AML/result/QC_GSE37642_PCA_coordinates.csv",
  row.names = FALSE
)

write.csv(
  pca_gse12417,
  "C:/Users/Xiong/Desktop/AML/result/QC_GSE12417_PCA_coordinates.csv",
  row.names = FALSE
)

cat("PCA QC finished.\n\n")

############################################################
## 10. Correlation heatmap QC
############################################################

cat("Drawing sample correlation heatmap...\n")

plot_correlation_heatmap <- function(expr, dataset_name, outfile){
  
  gene_sd <- apply(expr, 1, sd, na.rm = TRUE)
  expr2 <- expr[gene_sd > 0, , drop = FALSE]
  
  gene_var <- apply(expr2, 1, var, na.rm = TRUE)
  top_n <- min(2000, length(gene_var))
  top_genes <- names(sort(gene_var, decreasing = TRUE))[1:top_n]
  
  cor_mat <- cor(
    expr2[top_genes, , drop = FALSE],
    method = "pearson",
    use = "pairwise.complete.obs"
  )
  
  pdf(outfile, width = 8, height = 8)
  pheatmap::pheatmap(
    cor_mat,
    show_rownames = FALSE,
    show_colnames = FALSE,
    main = paste0(dataset_name, " sample correlation")
  )
  dev.off()
}

plot_correlation_heatmap(
  tcga_expr,
  "TCGA",
  "C:/Users/Xiong/Desktop/AML/figure/QC_TCGA_CorrelationHeatmap.pdf"
)

plot_correlation_heatmap(
  gse37642_expr,
  "GSE37642",
  "C:/Users/Xiong/Desktop/AML/figure/QC_GSE37642_CorrelationHeatmap.pdf"
)

plot_correlation_heatmap(
  gse12417_expr,
  "GSE12417",
  "C:/Users/Xiong/Desktop/AML/figure/QC_GSE12417_CorrelationHeatmap.pdf"
)

cat("Correlation heatmap QC finished.\n\n")

############################################################
## 11. Save multiple RData objects
############################################################

cat("Saving RData objects...\n")

save(
  tcga_expr,
  tcga_survival,
  tcga_clinical,
  file = "C:/Users/Xiong/Desktop/AML/object/Module02_TCGA_Preprocessed.RData"
)

save(
  gse37642_expr,
  gse37642_pheno,
  file = "C:/Users/Xiong/Desktop/AML/object/Module02_GSE37642_Preprocessed.RData"
)

save(
  gse12417_expr,
  gse12417_pheno,
  file = "C:/Users/Xiong/Desktop/AML/object/Module02_GSE12417_Preprocessed.RData"
)

save(
  tcga_expr,
  tcga_survival,
  tcga_clinical,
  gse37642_expr,
  gse37642_pheno,
  gse12417_expr,
  gse12417_pheno,
  qc_report,
  sample_missing_report,
  gene_missing_report,
  pca_tcga,
  pca_gse37642,
  pca_gse12417,
  file = "C:/Users/Xiong/Desktop/AML/object/Module02_All_Preprocessed_QC.RData"
)

cat("RData objects saved.\n\n")

############################################################
## 12. Generate text QC report
############################################################

report_file <- "C:/Users/Xiong/Desktop/AML/result/Module02_QC_Report.txt"

sink(report_file)

cat("Module02 QC Report\n")
cat("==================\n\n")

cat("Project: Lysine metabolism-related genes predict prognosis in AML\n\n")

cat("Input expression matrices:\n")
cat("TCGA: result/TCGA_expression.csv\n")
cat("GSE37642: result/GSE37642_expression.csv\n")
cat("GSE12417: result/GSE12417_expression.csv\n\n")

cat("Basic QC summary:\n")
print(qc_report)

cat("\nMissing value summary by dataset:\n")
cat("\nSample missing value maximum:\n")
print(
  sample_missing_report %>%
    dplyr::group_by(Dataset) %>%
    dplyr::summarise(
      Max_missing_values = max(Missing_values),
      Max_missing_rate = max(Missing_rate)
    )
)

cat("\nGene missing value maximum:\n")
print(
  gene_missing_report %>%
    dplyr::group_by(Dataset) %>%
    dplyr::summarise(
      Max_missing_values = max(Missing_values),
      Max_missing_rate = max(Missing_rate)
    )
)

cat("\nGenerated QC figures:\n")
cat("figure/QC_TCGA_Boxplot.pdf\n")
cat("figure/QC_TCGA_Density.pdf\n")
cat("figure/QC_TCGA_PCA.pdf\n")
cat("figure/QC_TCGA_CorrelationHeatmap.pdf\n\n")

cat("figure/QC_GSE37642_Boxplot.pdf\n")
cat("figure/QC_GSE37642_Density.pdf\n")
cat("figure/QC_GSE37642_PCA.pdf\n")
cat("figure/QC_GSE37642_CorrelationHeatmap.pdf\n\n")

cat("figure/QC_GSE12417_Boxplot.pdf\n")
cat("figure/QC_GSE12417_Density.pdf\n")
cat("figure/QC_GSE12417_PCA.pdf\n")
cat("figure/QC_GSE12417_CorrelationHeatmap.pdf\n\n")

cat("Generated RData objects:\n")
cat("object/Module02_TCGA_Preprocessed.RData\n")
cat("object/Module02_GSE37642_Preprocessed.RData\n")
cat("object/Module02_GSE12417_Preprocessed.RData\n")
cat("object/Module02_All_Preprocessed_QC.RData\n\n")

cat("Session information:\n")
print(sessionInfo())

sink()

cat("Text QC report saved:\n")
cat(report_file, "\n\n")

############################################################
## 13. Finish
############################################################

cat("=============================================\n")
cat("Module02B QC finished successfully.\n")
cat("=============================================\n")

