############################################################
##
## Module04_LysineGeneSet.R
##
## Project:
## Lysine metabolism-related genes predict prognosis in AML
##
## Purpose:
## 1. Build a comprehensive lysine metabolism-related gene set
## 2. Sources:
##    - MSigDB KEGG / GO / Reactome / GOMF
##    - KEGG official pathway hsa00310
## 3. Generate core, extended and full lysine gene lists
## 4. Extract lysine-related expression matrices from TCGA and GEO
##
## Input:
## result/TCGA_train_expression.csv
## result/TCGA_train_survival.csv
## result/GSE37642_validation_expression.csv
## result/GSE12417_validation_expression.csv
##
## Output:
## result/lysine_gene_list.csv
## result/lysine_gene_list_core.csv
## result/lysine_gene_list_extended.csv
## result/lysine_gene_list_full.csv
## result/Supplementary_Table_S1_LysineGeneSource.csv
## result/Supplementary_Table_S1_LysineGeneSource_WithCategory.csv
## result/Module04_GeneSet_Category.csv
## result/Module04_GeneSet_SourceSummary.csv
## result/TCGA_lysine_expression.csv
## result/GSE37642_lysine_expression.csv
## result/GSE12417_lysine_expression.csv
## result/Module04_Unmatched_LysineGenes.csv
## result/Module04_summary.csv
## result/Module04_Report.txt
## object/Module04_LysineGeneSet.RData
##
############################################################

rm(list = ls())
options(stringsAsFactors = FALSE)

cat("============================================\n")
cat("Module04: Lysine metabolism gene set\n")
cat("============================================\n\n")

############################################################
## 1. Package installation and loading
############################################################

if(!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}

cran_pkgs <- c(
  "data.table",
  "dplyr",
  "stringr",
  "msigdbr",
  "ggplot2"
)

bioc_pkgs <- c(
  "KEGGREST",
  "org.Hs.eg.db",
  "AnnotationDbi"
)

for(pkg in cran_pkgs){
  if(!requireNamespace(pkg, quietly = TRUE)){
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

for(pkg in bioc_pkgs){
  if(!requireNamespace(pkg, quietly = TRUE)){
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
  library(pkg, character.only = TRUE)
}

############################################################
## 2. Create folders
############################################################

dir.create("result", showWarnings = FALSE, recursive = TRUE)
dir.create("object", showWarnings = FALSE, recursive = TRUE)
dir.create("figure", showWarnings = FALSE, recursive = TRUE)
dir.create("log", showWarnings = FALSE, recursive = TRUE)

############################################################
## 3. Define input files
############################################################

tcga_expr_file <- "C:/Users/Xiong/Desktop/AML/result/TCGA_train_expression.csv"
tcga_surv_file <- "C:/Users/Xiong/Desktop/AML/result/TCGA_train_survival.csv"

gse37642_expr_file <- "C:/Users/Xiong/Desktop/AML/result/GSE37642_validation_expression.csv"
gse12417_expr_file <- "C:/Users/Xiong/Desktop/AML/result/GSE12417_validation_expression.csv"

input_files <- c(
  tcga_expr_file,
  tcga_surv_file,
  gse37642_expr_file,
  gse12417_expr_file
)

missing_files <- input_files[!file.exists(input_files)]

if(length(missing_files) > 0){
  stop(
    paste0(
      "The following input files are missing:\n",
      paste(missing_files, collapse = "\n")
    )
  )
}

cat("All Module03 input files found.\n\n")

############################################################
## 4. Helper function: read expression CSV
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

############################################################
## 5. Load expression matrices
############################################################

cat("Reading TCGA expression matrix...\n")

tcga_expr <- read_expr_csv(tcga_expr_file)

cat("TCGA expression dimension:\n")
print(dim(tcga_expr))

cat("Reading TCGA survival data...\n")

tcga_survival <- data.table::fread(
  tcga_surv_file,
  data.table = FALSE
)

cat("TCGA survival dimension:\n")
print(dim(tcga_survival))

cat("Reading GSE37642 expression matrix...\n")

gse37642_expr <- read_expr_csv(gse37642_expr_file)

cat("GSE37642 expression dimension:\n")
print(dim(gse37642_expr))

cat("Reading GSE12417 expression matrix...\n")

gse12417_expr <- read_expr_csv(gse12417_expr_file)

cat("GSE12417 expression dimension:\n")
print(dim(gse12417_expr))

############################################################
## 6. Get MSigDB human gene sets
############################################################

cat("\nDownloading MSigDB human gene sets using msigdbr...\n")

msig <- msigdbr::msigdbr(
  species = "Homo sapiens"
)

cat("MSigDB version:\n")
print(unique(msig$db_version))

cat("\nMSigDB collections:\n")
print(unique(msig$gs_collection))

############################################################
## 7. Search all lysine-related gene sets in MSigDB
############################################################

cat("\nSearching all lysine-related gene sets in MSigDB...\n")

lysine_set_table <- msig %>%
  dplyr::filter(
    grepl("LYSINE", gs_name, ignore.case = TRUE)
  ) %>%
  dplyr::select(
    gs_name,
    gs_collection,
    gs_subcollection
  ) %>%
  dplyr::distinct() %>%
  dplyr::arrange(gs_collection, gs_subcollection, gs_name)

cat("\nAll lysine-related gene sets found in current MSigDB version:\n")
print(lysine_set_table)

write.csv(
  lysine_set_table,
  "C:/Users/Xiong/Desktop/AML/result/Available_Lysine_GeneSets_in_MSigDB.csv",
  row.names = FALSE
)

############################################################
## 8. Define core and extended source labels
############################################################

## Core lysine metabolism gene sets:
## Directly related to lysine metabolism, degradation, or catabolism.
core_source_labels <- c(
  "KEGG_LYSINE_DEGRADATION",
  "GOBP_LYSINE_METABOLIC_PROCESS",
  "GOBP_L_LYSINE_CATABOLIC_PROCESS"
)

## Extended lysine-related gene sets:
## Related to lysine modification, acetylation, methylation,
## deacetylation, deacylation, ADP-ribosylation, and histone lysine regulation.
extended_source_labels <- c(
  "GOBP_PEPTIDYL_LYSINE_ACETYLATION",
  "GOBP_PEPTIDYL_LYSINE_METHYLATION",
  "GOBP_PEPTIDYL_LYSINE_MODIFICATION",
  "REACTOME_PKMTS_METHYLATE_HISTONE_LYSINES",
  "GOMF_NADPLUS_PROTEIN_LYSINE_ADP_RIBOSYLTRANSFERASE_ACTIVITY",
  "GOMF_NAD_DEPENDENT_PROTEIN_LYSINE_DEACYLASE_ACTIVITY",
  "GOMF_PROTEIN_LYSINE_6_OXIDASE_ACTIVITY",
  "GOMF_PROTEIN_LYSINE_DEACETYLASE_ACTIVITY"
)

core_source_labels_available <- intersect(
  core_source_labels,
  unique(msig$gs_name)
)

extended_source_labels_available <- intersect(
  extended_source_labels,
  unique(msig$gs_name)
)

cat("\nAvailable core source labels:\n")
print(core_source_labels_available)

cat("\nAvailable extended source labels:\n")
print(extended_source_labels_available)

selected_msig_sets <- unique(
  c(
    core_source_labels_available,
    extended_source_labels_available
  )
)

cat("\nFinal selected MSigDB gene sets:\n")
print(selected_msig_sets)

write.csv(
  data.frame(GeneSet = selected_msig_sets),
  "C:/Users/Xiong/Desktop/AML/result/Selected_Lysine_GeneSets_for_Modeling.csv",
  row.names = FALSE
)

############################################################
## 9. Extract MSigDB lysine genes with source information
############################################################

msig_lysine_source <- msig %>%
  dplyr::filter(
    gs_name %in% selected_msig_sets
  ) %>%
  dplyr::select(
    Gene = gene_symbol,
    Source = gs_name,
    Collection = gs_collection,
    Subcollection = gs_subcollection
  ) %>%
  dplyr::distinct()

cat("\nNumber of MSigDB lysine-related gene-source records:\n")
print(nrow(msig_lysine_source))

cat("\nNumber of unique MSigDB lysine-related genes:\n")
print(length(unique(msig_lysine_source$Gene)))

############################################################
## 10. Get KEGG official pathway hsa00310
############################################################

cat("\nDownloading KEGG official pathway hsa00310...\n")

kegg_official <- character(0)
kegg_ok <- TRUE

tryCatch({
  
  pathway <- KEGGREST::keggGet("hsa00310")
  
  gene_info <- pathway[[1]]$GENE
  
  if(!is.null(gene_info) && length(gene_info) > 0){
    
    ## KEGG GENE vector format:
    ## odd positions = Entrez IDs
    ## even positions = "GENE_SYMBOL; description"
    kegg_official <- gene_info[seq(2, length(gene_info), 2)]
    
    kegg_official <- gsub(";.*", "", kegg_official)
    kegg_official <- trimws(kegg_official)
    kegg_official <- unique(kegg_official)
    
  }
  
}, error = function(e){
  
  kegg_ok <<- FALSE
  message("KEGGREST failed to retrieve hsa00310.")
  message(e$message)
  
})

if(kegg_ok){
  cat("Number of KEGG official hsa00310 genes:\n")
  print(length(kegg_official))
} else {
  cat("KEGG official genes not retrieved. Continuing with MSigDB only.\n")
}

kegg_official_source <- data.frame(
  Gene = kegg_official,
  Source = "KEGG_OFFICIAL_hsa00310",
  Collection = "KEGG",
  Subcollection = "hsa00310",
  stringsAsFactors = FALSE
)

############################################################
## 11. Merge MSigDB + KEGG official sources
############################################################

source_table_raw <- dplyr::bind_rows(
  msig_lysine_source,
  kegg_official_source
) %>%
  dplyr::filter(!is.na(Gene), Gene != "") %>%
  dplyr::distinct()

cat("\nMerged source table records before QC:\n")
print(nrow(source_table_raw))

lysine_gene_raw <- unique(source_table_raw$Gene)

cat("\nRaw merged lysine gene number before QC:\n")
print(length(lysine_gene_raw))

############################################################
## 12. Gene Symbol QC using org.Hs.eg.db
############################################################

cat("\nChecking official gene symbols using org.Hs.eg.db...\n")

symbol_check <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = lysine_gene_raw,
  keytype = "SYMBOL",
  column = "SYMBOL",
  multiVals = "first"
)

valid_genes <- names(symbol_check)[!is.na(symbol_check)]

invalid_genes <- setdiff(lysine_gene_raw, valid_genes)

cat("\nValid gene symbols:\n")
print(length(valid_genes))

if(length(invalid_genes) > 0){
  cat("\nInvalid / unmapped gene symbols removed:\n")
  print(invalid_genes)
}

source_table <- source_table_raw %>%
  dplyr::filter(Gene %in% valid_genes) %>%
  dplyr::distinct()

############################################################
## 13. Build core, extended and full gene lists
############################################################

core_source_labels_final <- c(
  core_source_labels_available,
  "KEGG_OFFICIAL_hsa00310"
)

extended_source_labels_final <- extended_source_labels_available

lysine_gene_core <- source_table %>%
  dplyr::filter(Source %in% core_source_labels_final) %>%
  dplyr::pull(Gene) %>%
  unique() %>%
  sort()

lysine_gene_extended <- source_table %>%
  dplyr::filter(Source %in% extended_source_labels_final) %>%
  dplyr::pull(Gene) %>%
  unique() %>%
  sort()

lysine_gene_full <- unique(
  c(
    lysine_gene_core,
    lysine_gene_extended
  )
) %>%
  sort()

## Default downstream gene list:
## Use full gene set for the main analysis.
## Keep core gene set for sensitivity analysis.
lysine_gene <- lysine_gene_full

cat("\nCore lysine metabolism gene number:\n")
print(length(lysine_gene_core))

cat("\nExtended lysine-related gene number:\n")
print(length(lysine_gene_extended))

cat("\nFull lysine metabolism-related gene number:\n")
print(length(lysine_gene_full))

cat("\nDefault downstream lysine gene number:\n")
print(length(lysine_gene))

############################################################
## 14. Save gene lists and source tables
############################################################

write.csv(
  data.frame(Gene = lysine_gene_core),
  "C:/Users/Xiong/Desktop/AML/result/lysine_gene_list_core.csv",
  row.names = FALSE
)

write.csv(
  data.frame(Gene = lysine_gene_extended),
  "C:/Users/Xiong/Desktop/AML/result/lysine_gene_list_extended.csv",
  row.names = FALSE
)

write.csv(
  data.frame(Gene = lysine_gene_full),
  "C:/Users/Xiong/Desktop/AML/result/lysine_gene_list_full.csv",
  row.names = FALSE
)

write.csv(
  data.frame(Gene = lysine_gene),
  "C:/Users/Xiong/Desktop/AML/result/lysine_gene_list.csv",
  row.names = FALSE
)

gene_set_category <- data.frame(
  Source = c(
    core_source_labels_final,
    extended_source_labels_final
  ),
  Category = c(
    rep("Core_lysine_metabolism", length(core_source_labels_final)),
    rep("Extended_lysine_related_modification", length(extended_source_labels_final))
  ),
  stringsAsFactors = FALSE
)

source_table_with_category <- source_table %>%
  dplyr::left_join(
    gene_set_category,
    by = "Source"
  )

write.csv(
  source_table,
  "C:/Users/Xiong/Desktop/AML/result/Supplementary_Table_S1_LysineGeneSource.csv",
  row.names = FALSE
)

write.csv(
  source_table_with_category,
  "C:/Users/Xiong/Desktop/AML/result/Supplementary_Table_S1_LysineGeneSource_WithCategory.csv",
  row.names = FALSE
)

write.csv(
  gene_set_category,
  "C:/Users/Xiong/Desktop/AML/result/Module04_GeneSet_Category.csv",
  row.names = FALSE
)

source_summary <- source_table_with_category %>%
  dplyr::group_by(Source, Collection, Subcollection, Category) %>%
  dplyr::summarise(
    GeneNumber = dplyr::n_distinct(Gene),
    .groups = "drop"
  ) %>%
  dplyr::arrange(Category, Source)

write.csv(
  source_summary,
  "C:/Users/Xiong/Desktop/AML/result/Module04_GeneSet_SourceSummary.csv",
  row.names = FALSE
)

cat("\nGene list and source files saved.\n")

############################################################
## 15. Extract lysine expression matrices
############################################################

cat("\nExtracting lysine expression matrix from TCGA...\n")

tcga_lysine_genes <- intersect(
  rownames(tcga_expr),
  lysine_gene
)

cat("Matched lysine genes in TCGA:\n")
print(length(tcga_lysine_genes))

tcga_lysine_expr <- tcga_expr[
  tcga_lysine_genes,
  ,
  drop = FALSE
]

write.csv(
  tcga_lysine_expr,
  "C:/Users/Xiong/Desktop/AML/result/TCGA_lysine_expression.csv",
  quote = FALSE
)

cat("TCGA lysine expression dimension:\n")
print(dim(tcga_lysine_expr))

cat("\nExtracting lysine expression matrix from GSE37642...\n")

gse37642_lysine_genes <- intersect(
  rownames(gse37642_expr),
  lysine_gene
)

cat("Matched lysine genes in GSE37642:\n")
print(length(gse37642_lysine_genes))

gse37642_lysine_expr <- gse37642_expr[
  gse37642_lysine_genes,
  ,
  drop = FALSE
]

write.csv(
  gse37642_lysine_expr,
  "C:/Users/Xiong/Desktop/AML/result/GSE37642_lysine_expression.csv",
  quote = FALSE
)

cat("GSE37642 lysine expression dimension:\n")
print(dim(gse37642_lysine_expr))

cat("\nExtracting lysine expression matrix from GSE12417...\n")

gse12417_lysine_genes <- intersect(
  rownames(gse12417_expr),
  lysine_gene
)

cat("Matched lysine genes in GSE12417:\n")
print(length(gse12417_lysine_genes))

gse12417_lysine_expr <- gse12417_expr[
  gse12417_lysine_genes,
  ,
  drop = FALSE
]

write.csv(
  gse12417_lysine_expr,
  "C:/Users/Xiong/Desktop/AML/result/GSE12417_lysine_expression.csv",
  quote = FALSE
)

cat("GSE12417 lysine expression dimension:\n")
print(dim(gse12417_lysine_expr))

############################################################
## 16. Save unmatched gene information
############################################################

unmatched_tcga <- setdiff(lysine_gene, rownames(tcga_expr))
unmatched_gse37642 <- setdiff(lysine_gene, rownames(gse37642_expr))
unmatched_gse12417 <- setdiff(lysine_gene, rownames(gse12417_expr))

all_unmatched <- unique(c(
  unmatched_tcga,
  unmatched_gse37642,
  unmatched_gse12417
))

unmatched_table <- data.frame(
  Gene = all_unmatched,
  Missing_in_TCGA = all_unmatched %in% unmatched_tcga,
  Missing_in_GSE37642 = all_unmatched %in% unmatched_gse37642,
  Missing_in_GSE12417 = all_unmatched %in% unmatched_gse12417
)

write.csv(
  unmatched_table,
  "C:/Users/Xiong/Desktop/AML/result/Module04_Unmatched_LysineGenes.csv",
  row.names = FALSE
)

############################################################
## 17. Save RData objects
############################################################

save(
  lysine_gene,
  lysine_gene_core,
  lysine_gene_extended,
  lysine_gene_full,
  source_table,
  source_table_with_category,
  source_summary,
  gene_set_category,
  tcga_lysine_expr,
  gse37642_lysine_expr,
  gse12417_lysine_expr,
  tcga_survival,
  file = "C:/Users/Xiong/Desktop/AML/object/Module04_LysineGeneSet.RData"
)

cat("\nRData object saved:\n")
cat("object/Module04_LysineGeneSet.RData\n")

############################################################
## 18. Generate summary report
############################################################

summary_df <- data.frame(
  Item = c(
    "Core lysine metabolism genes",
    "Extended lysine-related genes",
    "Full lysine metabolism-related genes",
    "Default downstream lysine genes",
    "TCGA matched lysine genes",
    "GSE37642 matched lysine genes",
    "GSE12417 matched lysine genes",
    "TCGA samples",
    "GSE37642 samples",
    "GSE12417 samples"
  ),
  Number = c(
    length(lysine_gene_core),
    length(lysine_gene_extended),
    length(lysine_gene_full),
    length(lysine_gene),
    length(tcga_lysine_genes),
    length(gse37642_lysine_genes),
    length(gse12417_lysine_genes),
    ncol(tcga_expr),
    ncol(gse37642_expr),
    ncol(gse12417_expr)
  )
)

write.csv(
  summary_df,
  "C:/Users/Xiong/Desktop/AML/result/Module04_summary.csv",
  row.names = FALSE
)

cat("\nModule04 summary:\n")
print(summary_df)

############################################################
## 19. Text report
############################################################

report_file <- "C:/Users/Xiong/Desktop/AML/result/Module04_Report.txt"

sink(report_file)

cat("Module04 Report: Lysine metabolism gene set\n")
cat("===========================================\n\n")

cat("Project:\n")
cat("Lysine metabolism-related genes predict prognosis in AML\n\n")

cat("MSigDB version:\n")
print(unique(msig$db_version))

cat("\nAll lysine-related gene sets found in MSigDB:\n")
print(lysine_set_table)

cat("\nSelected core MSigDB gene sets:\n")
print(core_source_labels_available)

cat("\nSelected extended MSigDB gene sets:\n")
print(extended_source_labels_available)

cat("\nKEGG official pathway:\n")
cat("hsa00310: Lysine degradation\n")
cat("Retrieved genes:", length(kegg_official), "\n\n")

cat("Core lysine metabolism gene number:\n")
cat(length(lysine_gene_core), "\n\n")

cat("Extended lysine-related gene number:\n")
cat(length(lysine_gene_extended), "\n\n")

cat("Full lysine metabolism-related gene number:\n")
cat(length(lysine_gene_full), "\n\n")

cat("Default downstream gene number:\n")
cat(length(lysine_gene), "\n\n")

cat("Matched genes and samples:\n")
print(summary_df)

cat("\nSource summary:\n")
print(source_summary)

cat("\nGenerated files:\n")
cat("result/Available_Lysine_GeneSets_in_MSigDB.csv\n")
cat("result/Selected_Lysine_GeneSets_for_Modeling.csv\n")
cat("result/lysine_gene_list.csv\n")
cat("result/lysine_gene_list_core.csv\n")
cat("result/lysine_gene_list_extended.csv\n")
cat("result/lysine_gene_list_full.csv\n")
cat("result/Supplementary_Table_S1_LysineGeneSource.csv\n")
cat("result/Supplementary_Table_S1_LysineGeneSource_WithCategory.csv\n")
cat("result/Module04_GeneSet_Category.csv\n")
cat("result/Module04_GeneSet_SourceSummary.csv\n")
cat("result/TCGA_lysine_expression.csv\n")
cat("result/GSE37642_lysine_expression.csv\n")
cat("result/GSE12417_lysine_expression.csv\n")
cat("result/Module04_Unmatched_LysineGenes.csv\n")
cat("result/Module04_summary.csv\n")
cat("object/Module04_LysineGeneSet.RData\n\n")

cat("Session information:\n")
print(sessionInfo())

sink()

cat("\nText report saved:\n")
cat(report_file, "\n")

############################################################
## 20. Finish
############################################################

cat("\n============================================\n")
cat("Module04 finished successfully.\n")
cat("============================================\n")

