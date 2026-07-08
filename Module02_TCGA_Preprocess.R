Part01
Section1
Load Data
Section2
Gene Annotation
Section3
Gene QC

Part02
Section4
Clinical
Section5
Survival
Section6
Merge

Part03
Section7
QC
Section8
Save

############################################################
##
## Module02
## TCGA & GEO preprocessing
##
############################################################

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
    "data/TCGA/TCGA-LAML.star_counts.tsv",
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
    "result/TCGA_expression.csv",
    quote=FALSE
)

############################################################
##==============================
## Part 2 Survival
##==============================
############################################################

survival <- fread(
    "data/TCGA/TCGA-LAML.survival.tsv",
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

    select(sample=timepoint_submitter_id,
           OS.time,
           OS)

write.csv(
    survival2,
    "result/TCGA_survival.csv",
    row.names=FALSE
)

############################################################
##==============================
## Part3 Clinical
##==============================
############################################################

clinical <- fread(
    "data/TCGA/TCGA-LAML.clinical.tsv",
    data.table=FALSE
)

write.csv(
    clinical,
    "result/TCGA_clinical.csv",
    row.names=FALSE
)

############################################################
##==============================
## Part4 GSE37642
##==============================
############################################################

cat("Reading GSE37642...\n")

gse37642 <- getGEO(
    filename="data/GEO/GSE37642-GPL570_series_matrix.txt"
)

expr37642 <- exprs(gse37642)

pdata37642 <- pData(gse37642)

fdata37642 <- fData(gse37642)

############################################################
## probe -> gene
############################################################

gene37642 <- fdata37642$Gene.symbol

expr37642 <- as.data.frame(expr37642)

expr37642$Gene <- gene37642

expr37642 <- expr37642[
    expr37642$Gene!="",
]

expr37642 <- expr37642 %>%

    group_by(Gene) %>%

    summarise(across(everything(),mean))

expr37642 <- as.data.frame(expr37642)

rownames(expr37642) <- expr37642$Gene

expr37642$Gene <- NULL

write.csv(
    expr37642,
    "result/GSE37642_expression.csv"
)

write.csv(
    pdata37642,
    "result/GSE37642_pheno.csv"
)

############################################################
##==============================
## Part5 GSE12417
##==============================
############################################################

cat("Reading GSE12417...\n")

gse12417 <- getGEO(
    filename="data/GEO/GSE12417-GPL570_series_matrix.txt"
)

expr12417 <- exprs(gse12417)

pdata12417 <- pData(gse12417)

fdata12417 <- fData(gse12417)

############################################################
## probe -> gene
############################################################

gene12417 <- fdata12417$Gene.symbol

expr12417 <- as.data.frame(expr12417)

expr12417$Gene <- gene12417

expr12417 <- expr12417[
    expr12417$Gene!="",
]

expr12417 <- expr12417 %>%

    group_by(Gene) %>%

    summarise(across(everything(),mean))

expr12417 <- as.data.frame(expr12417)

rownames(expr12417) <- expr12417$Gene

expr12417$Gene <- NULL

write.csv(
    expr12417,
    "result/GSE12417_expression.csv"
)

write.csv(
    pdata12417,
    "result/GSE12417_pheno.csv"
)

############################################################

cat("\n")
cat("=============================\n")
cat("Module02 Finished\n")
cat("=============================\n")
