############################################################
##
## Module05_UnivariateCox.R
##
## Project:
## Lysine-associated gene signature predicts prognosis in AML
##
## Purpose:
## 1. Perform univariate Cox regression for lysine-associated genes
## 2. Identify prognosis-related lysine genes
## 3. Generate Cox result tables
## 4. Generate forest plot and volcano plot
##
## Input:
## result/TCGA_lysine_expression.csv
## result/TCGA_train_survival.csv
##
## Output:
## result/Module05_UnivariateCox_AllGenes.csv
## result/Module05_UnivariateCox_SignificantGenes.csv
## result/Module05_UnivariateCox_TopGenes.csv
## figure/Module05_UnivariateCox_ForestPlot.pdf
## figure/Module05_UnivariateCox_Volcano.pdf
## object/Module05_UnivariateCox.RData
##
############################################################

rm(list = ls())
options(stringsAsFactors = FALSE)

cat("============================================\n")
cat("Module05: Univariate Cox regression\n")
cat("============================================\n\n")

############################################################
## 1. Package installation and loading
############################################################

cran_pkgs <- c(
  "data.table",
  "dplyr",
  "survival",
  "ggplot2",
  "forestplot"
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

dir.create("result", showWarnings = FALSE, recursive = TRUE)
dir.create("figure", showWarnings = FALSE, recursive = TRUE)
dir.create("object", showWarnings = FALSE, recursive = TRUE)

############################################################
## 3. Define input files
############################################################

expr_file <- "C:/Users/Xiong/Desktop/AML/result/TCGA_lysine_expression.csv"
surv_file <- "C:/Users/Xiong/Desktop/AML/result/TCGA_train_survival.csv"

if(!file.exists(expr_file)){
  stop("Cannot find input file: result/TCGA_lysine_expression.csv")
}

if(!file.exists(surv_file)){
  stop("Cannot find input file: result/TCGA_train_survival.csv")
}

############################################################
## 4. Read expression matrix
############################################################

cat("Reading lysine expression matrix...\n")

expr <- data.table::fread(
  expr_file,
  data.table = FALSE,
  check.names = FALSE
)

rownames(expr) <- expr[, 1]
expr <- expr[, -1, drop = FALSE]

expr <- as.matrix(expr)
mode(expr) <- "numeric"

cat("Expression matrix dimension:\n")
print(dim(expr))

############################################################
## 5. Read survival data
############################################################

cat("Reading survival data...\n")

survival_data <- data.table::fread(
  surv_file,
  data.table = FALSE
)

cat("Survival data dimension:\n")
print(dim(survival_data))

cat("Survival columns:\n")
print(colnames(survival_data))

############################################################
## 6. Check survival columns
############################################################

required_cols <- c("sample", "OS.time", "OS")

missing_cols <- setdiff(required_cols, colnames(survival_data))

if(length(missing_cols) > 0){
  stop(
    paste0(
      "Missing required survival columns: ",
      paste(missing_cols, collapse = ", ")
    )
  )
}

############################################################
## 7. Match expression and survival samples
############################################################

cat("Matching expression samples with survival samples...\n")

colnames(expr) <- substr(colnames(expr), 1, 16)
survival_data$sample <- substr(survival_data$sample, 1, 16)

common_samples <- intersect(
  colnames(expr),
  survival_data$sample
)

cat("Number of matched samples:\n")
print(length(common_samples))

if(length(common_samples) < 30){
  stop("Too few matched samples. Please check sample IDs.")
}

expr <- expr[, common_samples, drop = FALSE]

survival_data <- survival_data[
  match(common_samples, survival_data$sample),
]

if(!all(colnames(expr) == survival_data$sample)){
  stop("Sample order between expression and survival data is not matched.")
}

cat("Sample matching finished.\n")

############################################################
## 8. Clean survival data
############################################################

survival_data$OS.time <- as.numeric(survival_data$OS.time)
survival_data$OS <- as.numeric(survival_data$OS)

keep <- !is.na(survival_data$OS.time) &
  !is.na(survival_data$OS) &
  survival_data$OS.time > 0

expr <- expr[, keep, drop = FALSE]
survival_data <- survival_data[keep, ]

cat("Final samples for Cox analysis:\n")
print(ncol(expr))

cat("Survival event table:\n")
print(table(survival_data$OS))

############################################################
## 9. Remove genes with zero variance
############################################################

gene_sd <- apply(expr, 1, sd, na.rm = TRUE)

expr <- expr[gene_sd > 0, , drop = FALSE]

cat("Genes retained after removing zero-variance genes:\n")
print(nrow(expr))

############################################################
## 10. Univariate Cox regression
############################################################

cat("Running univariate Cox regression...\n")

cox_result <- data.frame()

for(gene in rownames(expr)){
  
  gene_exp <- as.numeric(expr[gene, ])
  
  cox_data <- data.frame(
    OS.time = survival_data$OS.time,
    OS = survival_data$OS,
    expression = gene_exp
  )
  
  fit <- tryCatch(
    {
      survival::coxph(
        survival::Surv(OS.time, OS) ~ expression,
        data = cox_data
      )
    },
    error = function(e){
      return(NULL)
    }
  )
  
  if(is.null(fit)){
    next
  }
  
  fit_sum <- summary(fit)
  
  tmp <- data.frame(
    Gene = gene,
    HR = fit_sum$conf.int[1, "exp(coef)"],
    Lower95CI = fit_sum$conf.int[1, "lower .95"],
    Upper95CI = fit_sum$conf.int[1, "upper .95"],
    Coef = fit_sum$coefficients[1, "coef"],
    Z = fit_sum$coefficients[1, "z"],
    Pvalue = fit_sum$coefficients[1, "Pr(>|z|)"],
    stringsAsFactors = FALSE
  )
  
  cox_result <- rbind(cox_result, tmp)
}

cat("Univariate Cox finished.\n")
cat("Number of genes analyzed:\n")
print(nrow(cox_result))

############################################################
## 11. Multiple testing correction
############################################################

cox_result$FDR <- p.adjust(
  cox_result$Pvalue,
  method = "BH"
)

cox_result <- cox_result %>%
  dplyr::arrange(Pvalue)

############################################################
## 12. Define significant genes
############################################################

## 主分析建议先用 Pvalue < 0.05。
## FDR < 0.05 对样本量132可能过严，容易筛不出足够LASSO候选基因。
sig_cutoff <- 0.05

cox_sig <- cox_result %>%
  dplyr::filter(Pvalue < sig_cutoff)

cat("Significant genes with P < 0.05:\n")
print(nrow(cox_sig))

if(nrow(cox_sig) < 5){
  cat("\nWarning: fewer than 5 significant genes were found.\n")
  cat("You may consider using P < 0.1 for candidate selection, but report this clearly.\n")
}

############################################################
## 13. Add risk/protective label
############################################################

cox_result$RiskType <- ifelse(
  cox_result$HR > 1,
  "Risk",
  "Protective"
)

cox_sig$RiskType <- ifelse(
  cox_sig$HR > 1,
  "Risk",
  "Protective"
)

############################################################
## 14. Save Cox results
############################################################

write.csv(
  cox_result,
  "C:/Users/Xiong/Desktop/AML/result/Module05_UnivariateCox_AllGenes.csv",
  row.names = FALSE
)

write.csv(
  cox_sig,
  "C:/Users/Xiong/Desktop/AML/result/Module05_UnivariateCox_SignificantGenes.csv",
  row.names = FALSE
)

top_n <- min(30, nrow(cox_result))

cox_top <- cox_result[1:top_n, ]

write.csv(
  cox_top,
  "C:/Users/Xiong/Desktop/AML/result/Module05_UnivariateCox_TopGenes.csv",
  row.names = FALSE
)

cat("Cox result files saved.\n")

############################################################
## 15. Forest plot for top significant genes
############################################################

cat("Drawing forest plot...\n")

forest_data <- cox_result %>%
  dplyr::filter(Pvalue < 0.05) %>%
  dplyr::arrange(Pvalue)

if(nrow(forest_data) > 20){
  forest_data <- forest_data[1:20, ]
}

if(nrow(forest_data) >= 2){
  
  forest_data$Gene <- factor(
    forest_data$Gene,
    levels = rev(forest_data$Gene)
  )
  
  p_forest <- ggplot(
    forest_data,
    aes(
      x = Gene,
      y = HR,
      ymin = Lower95CI,
      ymax = Upper95CI
    )
  ) +
    geom_pointrange(size = 0.6) +
    geom_hline(
      yintercept = 1,
      linetype = "dashed"
    ) +
    coord_flip() +
    theme_bw() +
    labs(
      title = "Univariate Cox analysis of lysine-associated genes",
      x = "",
      y = "Hazard ratio (95% CI)"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.y = element_text(size = 9)
    )
  
  ggsave(
    "C:/Users/Xiong/Desktop/AML/figure/Module05_UnivariateCox_ForestPlot.pdf",
    p_forest,
    width = 7,
    height = max(5, nrow(forest_data) * 0.3)
  )
  
}else{
  
  cat("Not enough significant genes for forest plot.\n")
  
}

############################################################
## 16. Volcano-style Cox plot
############################################################

cat("Drawing Cox volcano plot...\n")

cox_result$logHR <- log2(cox_result$HR)

cox_result$minusLog10P <- -log10(cox_result$Pvalue)

cox_result$Significance <- "Not significant"

cox_result$Significance[
  cox_result$Pvalue < 0.05 & cox_result$HR > 1
] <- "Risk"

cox_result$Significance[
  cox_result$Pvalue < 0.05 & cox_result$HR < 1
] <- "Protective"

p_volcano <- ggplot(
  cox_result,
  aes(
    x = logHR,
    y = minusLog10P,
    color = Significance
  )
) +
  geom_point(alpha = 0.8, size = 2) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed"
  ) +
  theme_bw() +
  labs(
    title = "Univariate Cox volcano plot",
    x = "log2(HR)",
    y = "-log10(P value)"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

ggsave(
  "C:/Users/Xiong/Desktop/AML/figure/Module05_UnivariateCox_Volcano.pdf",
  p_volcano,
  width = 7,
  height = 6
)

############################################################
## 17. Extract candidate expression matrix for LASSO
############################################################

if(nrow(cox_sig) > 0){
  
  candidate_genes <- cox_sig$Gene
  
  candidate_expr <- expr[
    rownames(expr) %in% candidate_genes,
    ,
    drop = FALSE
  ]
  
  write.csv(
    candidate_expr,
    "C:/Users/Xiong/Desktop/AML/result/Module05_CandidateGene_Expression.csv",
    quote = FALSE
  )
  
}else{
  
  candidate_genes <- character(0)
  candidate_expr <- matrix(
    nrow = 0,
    ncol = ncol(expr)
  )
  
  colnames(candidate_expr) <- colnames(expr)
  
}

############################################################
## 18. Save RData object
############################################################

save(
  expr,
  survival_data,
  cox_result,
  cox_sig,
  cox_top,
  candidate_genes,
  candidate_expr,
  file = "C:/Users/Xiong/Desktop/AML/object/Module05_UnivariateCox.RData"
)

############################################################
## 19. Generate text report
############################################################

report_file <- "C:/Users/Xiong/Desktop/AML/result/Module05_Report.txt"

sink(report_file)

cat("Module05 Report: Univariate Cox regression\n")
cat("==========================================\n\n")

cat("Project:\n")
cat("Lysine-associated gene signature predicts prognosis in AML\n\n")

cat("Input files:\n")
cat(expr_file, "\n")
cat(surv_file, "\n\n")

cat("Expression matrix dimension:\n")
print(dim(expr))

cat("\nSurvival data dimension:\n")
print(dim(survival_data))

cat("\nSurvival event table:\n")
print(table(survival_data$OS))

cat("\nGenes analyzed by Cox regression:\n")
cat(nrow(cox_result), "\n\n")

cat("Significant genes with P < 0.05:\n")
cat(nrow(cox_sig), "\n\n")

cat("Top Cox genes:\n")
print(head(cox_result, 20))

cat("\nRisk genes among significant genes:\n")
print(sum(cox_sig$HR > 1))

cat("\nProtective genes among significant genes:\n")
print(sum(cox_sig$HR < 1))

cat("\nGenerated files:\n")
cat("result/Module05_UnivariateCox_AllGenes.csv\n")
cat("result/Module05_UnivariateCox_SignificantGenes.csv\n")
cat("result/Module05_UnivariateCox_TopGenes.csv\n")
cat("result/Module05_CandidateGene_Expression.csv\n")
cat("figure/Module05_UnivariateCox_ForestPlot.pdf\n")
cat("figure/Module05_UnivariateCox_Volcano.pdf\n")
cat("object/Module05_UnivariateCox.RData\n")

cat("\nSession information:\n")
print(sessionInfo())

sink()

cat("\nText report saved:\n")
cat(report_file, "\n")

############################################################
## 20. Finish
############################################################

cat("\n============================================\n")
cat("Module05 finished successfully.\n")
cat("============================================\n")