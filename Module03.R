############################################################
##
## Module03
## Data Integration
##
############################################################

rm(list=ls())

library(data.table)
library(dplyr)

dir.create("result",showWarnings=FALSE)

############################################################
## TCGA Expression
############################################################

expr <- read.csv(
  "C:/Users/Xiong/Desktop/AML/result/TCGA_expression.csv",
  row.names=1,
  check.names=FALSE
)

############################################################
## TCGA Survival
############################################################

survival <- read.csv(
  "C:/Users/Xiong/Desktop/AML/result/TCGA_survival.csv"
)

############################################################
## Clinical
############################################################

clinical <- read.csv(
  "C:/Users/Xiong/Desktop/AML/result/TCGA_clinical.csv"
)

############################################################
## sample ID
############################################################

colnames(expr) <- substr(
  colnames(expr),
  1,
  16
)

survival$sample <- substr(
  survival$sample,
  1,
  16
)

############################################################
## common samples
############################################################

common <- intersect(
  colnames(expr),
  survival$sample
)

length(common)

############################################################
## expression
############################################################

expr <- expr[,common]

############################################################
## survival
############################################################

survival <- survival[
  match(common,survival$sample),
]

############################################################
## check
############################################################

all(colnames(expr)==survival$sample)

############################################################
## remove NA
############################################################

keep <- !is.na(survival$OS.time)

expr <- expr[,keep]

survival <- survival[keep,]

############################################################
## save
############################################################

write.csv(
  expr,
  "C:/Users/Xiong/Desktop/AML/result/TCGA_train_expression.csv",
  quote=FALSE
)

write.csv(
  survival,
  "C:/Users/Xiong/Desktop/AML/result/TCGA_train_survival.csv",
  row.names=FALSE
)

############################################################
## GSE37642
############################################################

expr37642 <- read.csv(
  "C:/Users/Xiong/Desktop/AML/result/GSE37642_expression.csv",
  row.names=1,
  check.names=FALSE
)

pdata37642 <- read.csv(
  "C:/Users/Xiong/Desktop/AML/result/GSE37642_pheno.csv"
)

write.csv(
  expr37642,
  "C:/Users/Xiong/Desktop/AML/result/GSE37642_validation_expression.csv",
  quote=FALSE
)

write.csv(
  pdata37642,
  "C:/Users/Xiong/Desktop/AML/result/GSE37642_validation_pheno.csv",
  row.names=FALSE
)

############################################################
## GSE12417
############################################################

expr12417 <- read.csv(
  "C:/Users/Xiong/Desktop/AML/result/GSE12417_expression.csv",
  row.names=1,
  check.names=FALSE
)

pdata12417 <- read.csv(
  "C:/Users/Xiong/Desktop/AML/result/GSE12417_pheno.csv"
)

write.csv(
  expr12417,
  "C:/Users/Xiong/Desktop/AML/result/GSE12417_validation_expression.csv",
  quote=FALSE
)

write.csv(
  pdata12417,
  "C:/Users/Xiong/Desktop/AML/result/GSE12417_validation_pheno.csv",
  row.names=FALSE
)

############################################################
## Summary
############################################################

cat("\n==============================\n")

cat("TCGA samples:",ncol(expr),"\n")

cat("TCGA genes:",nrow(expr),"\n")

cat("GSE37642 samples:",ncol(expr37642),"\n")

cat("GSE12417 samples:",ncol(expr12417),"\n")

cat("==============================\n")

cat("Module03 Finished.\n")
