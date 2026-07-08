###############################################################
##
## Module01_DataDownload.R
##
## Project:
## Lysine metabolism-related genes predict prognosis in AML
##
## Description:
## Project initialization
## Check packages
## Create project directories
## Check TCGA / GEO files
##
## Author:
##
## Version: 1.0
##
###############################################################

rm(list = ls())

options(stringsAsFactors = FALSE)

cat("=============================================\n")
cat(" AML Lysine Project : Module 01\n")
cat(" Data Download & Project Initialization\n")
cat("=============================================\n\n")

###############################################################
## Section 1
## Load Packages
###############################################################

cran_packages <- c(

    "data.table",
    "tidyverse",
    "stringr",
    "dplyr",
    "tibble"

)

bioc_packages <- c(

    "GEOquery",
    "limma",
    "edgeR",
    "AnnotationDbi",
    "org.Hs.eg.db",
    "SummarizedExperiment",
    "Biobase"

)

###############################################################
## Install CRAN Packages
###############################################################

for(pkg in cran_packages){

    if(!require(pkg,
                character.only = TRUE)){

        install.packages(pkg)

        library(pkg,
                character.only = TRUE)

    }

}

###############################################################
## Install Bioconductor Packages
###############################################################

if(!requireNamespace("BiocManager",
                     quietly = TRUE))

    install.packages("BiocManager")

for(pkg in bioc_packages){

    if(!require(pkg,
                character.only = TRUE)){

        BiocManager::install(pkg,
                             ask = FALSE)

        library(pkg,
                character.only = TRUE)

    }

}

###############################################################
## Session Information
###############################################################

cat("\nLoaded Packages:\n")

sessionInfo()

###############################################################
## Section 2
## Set Working Directory
###############################################################

## Please modify this path

project_dir <- "D:/AML_Lysine_Project"

setwd(project_dir)

cat("\nWorking directory:\n")

print(getwd())


###############################################################
## Section 3
## Create Project Folder
###############################################################

dir_list <- c(

"data",

"data/TCGA",

"data/GEO",

"data/GeneSet",

"object",

"result",

"figure",

"script",

"log"

)

for(i in dir_list){

    if(!dir.exists(i)){

        dir.create(i,
                  recursive = TRUE)

        cat(i," created.\n")

    }

}

cat("\nDirectory checking finished.\n")

###############################################################
## Section 4
## TCGA Data Check
###############################################################

cat("\nChecking TCGA files...\n")

tcga_expression <-
"data/TCGA/TCGA-LAML.star_counts.tsv"

tcga_clinical <-
"data/TCGA/TCGA-LAML.clinical.tsv"

tcga_survival <-
"data/TCGA/TCGA-LAML.survival.tsv"

if(file.exists(tcga_expression)){

    cat("Expression file found.\n")

}else{

    stop("Cannot find TCGA expression file.")

}

if(file.exists(tcga_clinical)){

    cat("Clinical file found.\n")

}else{

    stop("Cannot find clinical file.")

}

if(file.exists(tcga_survival)){

    cat("Survival file found.\n")

}else{

    stop("Cannot find survival file.")

}

###############################################################
## Section 5
## GEO Data Check
###############################################################

geo37642 <-
"data/GEO/GSE37642/GSE37642-GPL570_series_matrix.txt.gz"

geo12417 <-
"data/GEO/GSE12417/GSE12417-GPL570_series_matrix.txt.gz"

if(file.exists(geo37642)){

    cat("GSE37642 found.\n")

}else{

    stop("Cannot find GSE37642.")

}

if(file.exists(geo12417)){

    cat("GSE12417 found.\n")

}else{

    stop("Cannot find GSE12417.")

}
