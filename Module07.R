############################################################
##
## Module07_MultivariateCox_Stepwise.R
##
## Project:
## Lysine-associated gene signature predicts prognosis in AML
##
## Purpose:
## 1. Use LASSO-selected genes from Module06
## 2. Perform multivariate Cox regression
## 3. Apply stepwise Cox regression based on AIC
## 4. Generate final prognostic gene signature
## 5. Output final coefficients and risk score formula
##
## Input:
## result/Module06_LASSO_Expression_lambda_min.csv
## result/TCGA_train_survival.csv
##
## Output:
## result/Module07_MultivariateCox_FullModel.csv
## result/Module07_StepwiseCox_FinalModel.csv
## result/Module07_FinalModel_Genes.csv
## result/Module07_RiskScore_Formula.txt
## result/Module07_FinalModel_Expression.csv
## figure/Module07_FinalModel_ForestPlot.pdf
## figure/Module07_FinalModel_CoefficientBarplot.pdf
## object/Module07_MultivariateCox_Stepwise.RData
##
############################################################

rm(list = ls())
options(stringsAsFactors = FALSE)

cat("============================================\n")
cat("Module07: Multivariate Cox + Stepwise Cox\n")
cat("============================================\n\n")

############################################################
## 1. Packages
############################################################

cran_pkgs <- c(
  "data.table",
  "dplyr",
  "survival",
  "MASS",
  "ggplot2",
  "stringr"
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
## 3. Input files
############################################################

expr_file <- "C:/Users/Xiong/Desktop/AML/result/Module06_LASSO_Expression_lambda_min.csv"
surv_file <- "C:/Users/Xiong/Desktop/AML/result/TCGA_train_survival.csv"

if(!file.exists(expr_file)){
  stop("Cannot find result/Module06_LASSO_Expression_lambda_min.csv")
}

if(!file.exists(surv_file)){
  stop("Cannot find result/TCGA_train_survival.csv")
}

############################################################
## 4. Read LASSO-selected expression matrix
############################################################

cat("Reading LASSO-selected expression matrix...\n")

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

required_cols <- c("sample", "OS.time", "OS")

missing_cols <- setdiff(required_cols, colnames(survival_data))

if(length(missing_cols) > 0){
  stop(
    paste0(
      "Missing survival columns: ",
      paste(missing_cols, collapse = ", ")
    )
  )
}

survival_data$sample <- substr(survival_data$sample, 1, 16)
survival_data$OS.time <- as.numeric(survival_data$OS.time)
survival_data$OS <- as.numeric(survival_data$OS)

############################################################
## 6. Match samples
############################################################

cat("Matching samples...\n")

colnames(expr) <- substr(colnames(expr), 1, 16)

common_samples <- intersect(
  colnames(expr),
  survival_data$sample
)

cat("Matched samples:\n")
print(length(common_samples))

if(length(common_samples) < 30){
  stop("Too few matched samples for multivariate Cox.")
}

expr <- expr[, common_samples, drop = FALSE]

survival_data <- survival_data[
  match(common_samples, survival_data$sample),
]

if(!all(colnames(expr) == survival_data$sample)){
  stop("Expression and survival samples are not in the same order.")
}

############################################################
## 7. Remove missing survival data
############################################################

keep <- !is.na(survival_data$OS.time) &
  !is.na(survival_data$OS) &
  survival_data$OS.time > 0

expr <- expr[, keep, drop = FALSE]
survival_data <- survival_data[keep, ]

cat("Final samples for multivariate Cox:\n")
print(ncol(expr))

cat("Event table:\n")
print(table(survival_data$OS))

############################################################
## 8. Prepare Cox dataframe
############################################################

cox_df <- data.frame(
  sample = survival_data$sample,
  OS.time = survival_data$OS.time,
  OS = survival_data$OS,
  t(expr),
  check.names = FALSE
)

gene_names <- rownames(expr)

cat("Genes entered into multivariate Cox:\n")
print(gene_names)

############################################################
## 9. Rename genes to syntactically valid variable names
############################################################

safe_gene_names <- make.names(gene_names)

name_map <- data.frame(
  OriginalGene = gene_names,
  SafeGene = safe_gene_names,
  stringsAsFactors = FALSE
)

colnames(cox_df)[4:ncol(cox_df)] <- safe_gene_names

write.csv(
  name_map,
  "C:/Users/Xiong/Desktop/AML/result/Module07_GeneName_Map.csv",
  row.names = FALSE
)

############################################################
## 10. Full multivariate Cox model
############################################################

cat("Running full multivariate Cox model...\n")

formula_full <- as.formula(
  paste0(
    "Surv(OS.time, OS) ~ ",
    paste(safe_gene_names, collapse = " + ")
  )
)

full_cox <- survival::coxph(
  formula_full,
  data = cox_df
)

full_summary <- summary(full_cox)

full_result <- data.frame(
  Gene = rownames(full_summary$coefficients),
  Coef = full_summary$coefficients[, "coef"],
  HR = full_summary$conf.int[, "exp(coef)"],
  Lower95CI = full_summary$conf.int[, "lower .95"],
  Upper95CI = full_summary$conf.int[, "upper .95"],
  Pvalue = full_summary$coefficients[, "Pr(>|z|)"],
  stringsAsFactors = FALSE
)

full_result <- full_result %>%
  dplyr::left_join(
    name_map,
    by = c("Gene" = "SafeGene")
  ) %>%
  dplyr::mutate(
    Gene = OriginalGene
  ) %>%
  dplyr::select(
    Gene,
    Coef,
    HR,
    Lower95CI,
    Upper95CI,
    Pvalue
  )

write.csv(
  full_result,
  "C:/Users/Xiong/Desktop/AML/result/Module07_MultivariateCox_FullModel.csv",
  row.names = FALSE
)

cat("Full multivariate Cox finished.\n")

############################################################
## 11. Stepwise Cox regression by AIC
############################################################

cat("Running stepwise Cox regression based on AIC...\n")

set.seed(12345)

step_cox <- MASS::stepAIC(
  full_cox,
  direction = "both",
  trace = TRUE
)

step_summary <- summary(step_cox)

############################################################
## 12. Extract final stepwise model result
############################################################

step_result <- data.frame(
  Gene = rownames(step_summary$coefficients),
  Coef = step_summary$coefficients[, "coef"],
  HR = step_summary$conf.int[, "exp(coef)"],
  Lower95CI = step_summary$conf.int[, "lower .95"],
  Upper95CI = step_summary$conf.int[, "upper .95"],
  Pvalue = step_summary$coefficients[, "Pr(>|z|)"],
  stringsAsFactors = FALSE
)

step_result <- step_result %>%
  dplyr::left_join(
    name_map,
    by = c("Gene" = "SafeGene")
  ) %>%
  dplyr::mutate(
    Gene = OriginalGene
  ) %>%
  dplyr::select(
    Gene,
    Coef,
    HR,
    Lower95CI,
    Upper95CI,
    Pvalue
  ) %>%
  dplyr::arrange(Pvalue)

step_result$RiskType <- ifelse(
  step_result$Coef > 0,
  "Risk",
  "Protective"
)

final_genes <- step_result$Gene

cat("Final genes selected by stepwise Cox:\n")
print(final_genes)

cat("Number of final genes:\n")
print(length(final_genes))

############################################################
## 13. Save final model tables
############################################################

write.csv(
  step_result,
  "C:/Users/Xiong/Desktop/AML/result/Module07_StepwiseCox_FinalModel.csv",
  row.names = FALSE
)

write.csv(
  data.frame(Gene = final_genes),
  "C:/Users/Xiong/Desktop/AML/result/Module07_FinalModel_Genes.csv",
  row.names = FALSE
)

############################################################
## 14. Save final model expression matrix
############################################################

final_expr <- expr[
  rownames(expr) %in% final_genes,
  ,
  drop = FALSE
]

write.csv(
  final_expr,
  "C:/Users/Xiong/Desktop/AML/result/Module07_FinalModel_Expression.csv",
  quote = FALSE
)

############################################################
## 15. Generate risk score formula
############################################################

formula_terms <- paste0(
  "(",
  round(step_result$Coef, 6),
  " × ",
  step_result$Gene,
  ")"
)

risk_formula <- paste(
  formula_terms,
  collapse = " + "
)

risk_formula_text <- paste0(
  "Risk score = ",
  risk_formula
)

cat("\nRisk score formula:\n")
cat(risk_formula_text, "\n")

writeLines(
  risk_formula_text,
  "C:/Users/Xiong/Desktop/AML/result/Module07_RiskScore_Formula.txt"
)

############################################################
## 16. Calculate preliminary risk score in TCGA
############################################################

coef_vec <- step_result$Coef
names(coef_vec) <- step_result$Gene

final_expr_ordered <- final_expr[
  names(coef_vec),
  ,
  drop = FALSE
]

risk_score <- as.numeric(
  t(final_expr_ordered) %*% coef_vec
)

risk_df <- data.frame(
  sample = colnames(final_expr_ordered),
  RiskScore = risk_score,
  OS.time = survival_data$OS.time,
  OS = survival_data$OS,
  stringsAsFactors = FALSE
)

write.csv(
  risk_df,
  "C:/Users/Xiong/Desktop/AML/result/Module07_TCGA_Preliminary_RiskScore.csv",
  row.names = FALSE
)

############################################################
## 17. Forest plot for final model
############################################################

cat("Drawing final model forest plot...\n")

forest_df <- step_result %>%
  dplyr::mutate(
    Gene = factor(Gene, levels = rev(Gene))
  )

p_forest <- ggplot(
  forest_df,
  aes(
    x = Gene,
    y = HR,
    ymin = Lower95CI,
    ymax = Upper95CI,
    color = RiskType
  )
) +
  geom_pointrange(size = 0.7) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    linewidth = 0.4
  ) +
  coord_flip() +
  scale_y_log10() +
  scale_color_manual(
    values = c(
      "Risk" = "#D73027",
      "Protective" = "#4575B4"
    )
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Final multivariate Cox model",
    x = "",
    y = "Hazard ratio (log scale)"
  )

ggsave(
  "C:/Users/Xiong/Desktop/AML/figure/Module07_FinalModel_ForestPlot.pdf",
  p_forest,
  width = 7,
  height = max(4, nrow(forest_df) * 0.45)
)

ggsave(
  "C:/Users/Xiong/Desktop/AML/figure/Module07_FinalModel_ForestPlot.png",
  p_forest,
  width = 7,
  height = max(4, nrow(forest_df) * 0.45),
  dpi = 300
)

############################################################
## 18. Coefficient barplot
############################################################

coef_plot_df <- step_result %>%
  dplyr::mutate(
    Gene = factor(Gene, levels = Gene[order(Coef)])
  )

p_coef <- ggplot(
  coef_plot_df,
  aes(
    x = Gene,
    y = Coef,
    fill = RiskType
  )
) +
  geom_bar(
    stat = "identity",
    width = 0.7
  ) +
  coord_flip() +
  scale_fill_manual(
    values = c(
      "Risk" = "#D73027",
      "Protective" = "#4575B4"
    )
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_blank(),
    panel.grid = element_blank()
  ) +
  labs(
    title = "Final model coefficients",
    x = "",
    y = "Coefficient"
  )

ggsave(
  "C:/Users/Xiong/Desktop/AML/figure/Module07_FinalModel_CoefficientBarplot.pdf",
  p_coef,
  width = 7,
  height = max(4, nrow(coef_plot_df) * 0.45)
)

ggsave(
  "C:/Users/Xiong/Desktop/AML/figure/Module07_FinalModel_CoefficientBarplot.png",
  p_coef,
  width = 7,
  height = max(4, nrow(coef_plot_df) * 0.45),
  dpi = 300
)

############################################################
## 19. Save RData object
############################################################

save(
  expr,
  survival_data,
  cox_df,
  full_cox,
  step_cox,
  full_result,
  step_result,
  final_genes,
  final_expr,
  coef_vec,
  risk_df,
  risk_formula_text,
  file = "C:/Users/Xiong/Desktop/AML/object/Module07_MultivariateCox_Stepwise.RData"
)

############################################################
## 20. Generate report
############################################################

report_file <- "C:/Users/Xiong/Desktop/AML/result/Module07_Report.txt"

sink(report_file)

cat("Module07 Report: Multivariate Cox and stepwise Cox regression\n")
cat("============================================================\n\n")

cat("Project:\n")
cat("Lysine-associated gene signature predicts prognosis in AML\n\n")

cat("Input files:\n")
cat(expr_file, "\n")
cat(surv_file, "\n\n")

cat("Input LASSO-selected expression dimension:\n")
print(dim(expr))

cat("\nSurvival data dimension:\n")
print(dim(survival_data))

cat("\nEvent table:\n")
print(table(survival_data$OS))

cat("\nGenes entered into full multivariate Cox model:\n")
print(gene_names)

cat("\nFull multivariate Cox result:\n")
print(full_result)

cat("\nFinal genes selected by stepwise Cox:\n")
print(final_genes)

cat("\nFinal model result:\n")
print(step_result)

cat("\nRisk score formula:\n")
cat(risk_formula_text, "\n")

cat("\nGenerated files:\n")
cat("result/Module07_MultivariateCox_FullModel.csv\n")
cat("result/Module07_StepwiseCox_FinalModel.csv\n")
cat("result/Module07_FinalModel_Genes.csv\n")
cat("result/Module07_RiskScore_Formula.txt\n")
cat("result/Module07_FinalModel_Expression.csv\n")
cat("result/Module07_TCGA_Preliminary_RiskScore.csv\n")
cat("figure/Module07_FinalModel_ForestPlot.pdf\n")
cat("figure/Module07_FinalModel_CoefficientBarplot.pdf\n")
cat("object/Module07_MultivariateCox_Stepwise.RData\n")

cat("\nSession information:\n")
print(sessionInfo())

sink()

cat("\nText report saved:\n")
cat(report_file, "\n")

############################################################
## 21. Finish
############################################################

cat("\n============================================\n")
cat("Module07 finished successfully.\n")
cat("============================================\n")

