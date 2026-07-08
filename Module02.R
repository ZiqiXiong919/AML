## Module02
## TCGA & GEO preprocessing
##
rm(list=ls())
options(stringsAsFactors = FALSE)

library(data.table)
library(dplyr)
library(GEOquery)
library(org.Hs.eg.db)
library(AnnotationDbi)

############################################################
## Create folders
############################################################

dir.create("result",showWarnings = FALSE)
dir.create("figure",showWarnings = FALSE)

############################################################
##==============================
## Part 1 TCGA-LAML
##==============================
############################################################

cat("Reading TCGA expression...\n")

expr_raw <- fread(
  "C:/Users/Xiong/Desktop/AML/data/TCGA/TCGA-LAML.star_counts.tsv",
  data.table=FALSE
)

cat("Dimension:\n")
dim(expr_raw)

head(colnames(expr_raw))

############################################################
## Remove version number
############################################################

expr_raw$Ensembl_ID <- gsub("\\..*","",expr_raw$Ensembl_ID)

############################################################
## Ensembl -> Gene Symbol
############################################################

geneAnno <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys=expr_raw$Ensembl_ID,
  keytype="ENSEMBL",
  columns="SYMBOL"
)

geneAnno <- geneAnno[!duplicated(geneAnno$ENSEMBL),]

expr_raw <- merge(
  expr_raw,
  geneAnno,
  by.x="Ensembl_ID",
  by.y="ENSEMBL"
)

############################################################
## Remove NA
############################################################

expr_raw <- expr_raw[
  !is.na(expr_raw$SYMBOL),
]

############################################################
## Keep expression matrix
############################################################

expr <- expr_raw[, !(names(expr_raw) %in%
                       c("Ensembl_ID"))]

############################################################
## Average duplicated genes
############################################################

expr <- expr %>%
  group_by(SYMBOL) %>%
  summarise(across(everything(),mean))

expr <- as.data.frame(expr)

rownames(expr) <- expr$SYMBOL

expr$SYMBOL <- NULL

############################################################
## log2 transform
############################################################

expr <- log2(expr+1)

cat("TCGA expression dimension:\n")

dim(expr)

############################################################
## Save
############################################################

write.csv(
  expr,
  "C:/Users/Xiong/Desktop/AML/result/TCGA_expression.csv",
  quote=FALSE
)

############################################################
##==============================
## Part 2 Survival
##==============================
############################################################

survival <- fread(
  "C:/Users/Xiong/Desktop/AML/data/TCGA/TCGA-LAML.survival.tsv",
  data.table=FALSE
)

head(survival)

############################################################
## Find useful columns
############################################################

colnames(survival)

############################################################
## Example
############################################################

survival2 <- survival %>%
  dplyr::select(
    sample,
    OS.time,
    OS
  )

write.csv(
  survival2,
  "C:/Users/Xiong/Desktop/AML/result/TCGA_survival.csv",
  row.names=FALSE
)

############################################################
##==============================
## Part3 Clinical
##==============================
############################################################

clinical <- fread(
  "C:/Users/Xiong/Desktop/AML/data/TCGA/TCGA-LAML.clinical.tsv",
  data.table=FALSE
)

write.csv(
  clinical,
  "C:/Users/Xiong/Desktop/AML/result/TCGA_clinical.csv",
  row.names=FALSE
)

############################################################
##==============================
## Part4 GSE37642
##==============================
############################################################
############################################################
## GSE37642 preprocessing - corrected version
############################################################

library(GEOquery)
library(limma)
library(dplyr)

gse37642 <- getGEO(
  filename = "C:/Users/Xiong/Desktop/AML/data/GEO/GSE37642/GSE37642-GPL570_series_matrix.txt.gz"
)

expr37642 <- exprs(gse37642)
pdata37642 <- pData(gse37642)
fdata37642 <- fData(gse37642)

cat("Raw expression dimension:\n")
print(dim(expr37642))

cat("Feature annotation columns:\n")
print(colnames(fdata37642))

############################################################
## Automatically detect gene symbol column
############################################################

possible_gene_cols <- c(
  "Gene Symbol",
  "Gene symbol",
  "GENE_SYMBOL",
  "Gene.symbol",
  "gene_symbol",
  "SYMBOL",
  "Symbol"
)

gene_col <- intersect(possible_gene_cols, colnames(fdata37642))

if(length(gene_col) == 0){
  stop("No Gene Symbol column found. Please run colnames(fdata37642) and check the correct column name.")
}

gene_col <- gene_col[1]

cat("Using gene symbol column: ", gene_col, "\n")

gene37642 <- fdata37642[[gene_col]]

############################################################
## Clean gene symbols
############################################################

gene37642 <- as.character(gene37642)

# 有些探针对应多个基因，例如 "A /// B"，只保留第一个
gene37642 <- sapply(strsplit(gene37642, " /// "), `[`, 1)

# 去掉前后空格
gene37642 <- trimws(gene37642)

############################################################
## Remove empty genes
############################################################

keep <- !is.na(gene37642) & gene37642 != ""

expr37642 <- expr37642[keep, ]
gene37642 <- gene37642[keep]

cat("After removing empty gene symbols:\n")
print(dim(expr37642))

############################################################
## Probe to Gene
############################################################

rownames(expr37642) <- gene37642

expr37642 <- limma::avereps(expr37642)

cat("After averaging duplicated genes:\n")
print(dim(expr37642))

############################################################
## Check log2 status
############################################################

cat("Expression summary:\n")
print(summary(as.vector(expr37642)))

if(max(expr37642, na.rm = TRUE) > 100){
  expr37642 <- log2(expr37642 + 1)
  cat("Log2 transformation performed.\n")
}else{
  cat("Expression data already appear to be log2-transformed.\n")
}

############################################################
## Save
############################################################

write.csv(
  expr37642,
  "C:/Users/Xiong/Desktop/AML/result/GSE37642_expression.csv",
  quote = FALSE
)

write.csv(
  pdata37642,
  "C:/Users/Xiong/Desktop/AML/result/GSE37642_pheno.csv",
  row.names = FALSE
)

cat("GSE37642 preprocessing finished.\n")
cat("Final GSE37642 expression dimension:\n")
print(dim(expr37642))

############################################################
##==============================
## Part5 GSE12417
##==============================
############################################################

############################################################
## Part5 GSE12417 preprocessing - corrected version
############################################################

library(GEOquery)
library(limma)
library(dplyr)

cat("Reading GSE12417...\n")

gse12417 <- getGEO(
  filename = "C:/Users/Xiong/Desktop/AML/data/GEO/GSE12417/GSE12417-GPL570_series_matrix.txt.gz"
)

expr12417 <- exprs(gse12417)
pdata12417 <- pData(gse12417)
fdata12417 <- fData(gse12417)

cat("Raw expression dimension:\n")
print(dim(expr12417))

cat("Feature annotation columns:\n")
print(colnames(fdata12417))

############################################################
## Automatically detect Gene Symbol column
############################################################

possible_gene_cols <- c(
  "Gene Symbol",
  "Gene symbol",
  "GENE_SYMBOL",
  "Gene.symbol",
  "gene_symbol",
  "SYMBOL",
  "Symbol"
)

gene_col <- intersect(possible_gene_cols, colnames(fdata12417))

if(length(gene_col) == 0){
  stop("No Gene Symbol column found. Please run colnames(fdata12417) and check the correct column name.")
}

gene_col <- gene_col[1]

cat("Using gene symbol column: ", gene_col, "\n")

gene12417 <- fdata12417[[gene_col]]

############################################################
## Clean Gene Symbols
############################################################

gene12417 <- as.character(gene12417)

# Some probes map to multiple genes, such as "A /// B".
# Keep the first symbol to avoid ambiguous mapping.
gene12417 <- sapply(strsplit(gene12417, " /// "), `[`, 1)

# Remove spaces
gene12417 <- trimws(gene12417)

############################################################
## Remove empty / NA genes
############################################################

keep <- !is.na(gene12417) & gene12417 != ""

expr12417 <- expr12417[keep, ]
gene12417 <- gene12417[keep]

cat("After removing empty gene symbols:\n")
print(dim(expr12417))

############################################################
## Probe to Gene
############################################################

rownames(expr12417) <- gene12417

expr12417 <- limma::avereps(expr12417)

cat("After averaging duplicated genes:\n")
print(dim(expr12417))

############################################################
## Check log2 status
############################################################

cat("Expression summary:\n")
print(summary(as.vector(expr12417)))

if(max(expr12417, na.rm = TRUE) > 100){
  expr12417 <- log2(expr12417 + 1)
  cat("Log2 transformation performed.\n")
}else{
  cat("Expression data already appear to be log2-transformed.\n")
}

############################################################
## Save
############################################################

write.csv(
  expr12417,
  "C:/Users/Xiong/Desktop/AML/result/GSE12417_expression.csv",
  quote = FALSE
)

write.csv(
  pdata12417,
  "C:/Users/Xiong/Desktop/AML/result/GSE12417_pheno.csv",
  row.names = FALSE
)

cat("GSE12417 preprocessing finished.\n")
cat("Final GSE12417 expression dimension:\n")
print(dim(expr12417))
############################################################

cat("\n")
cat("=============================\n")
cat("Module02 Finished\n")
cat("=============================\n")