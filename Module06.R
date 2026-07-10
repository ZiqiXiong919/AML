############################################################
##
## Module06_LASSO.R
##
## Project:
## Lysine-associated gene signature predicts prognosis in AML
##
## Purpose:
## 1. Perform LASSO Cox regression using candidate genes from Module05
## 2. Select prognostic model genes
## 3. Output LASSO coefficients
## 4. Generate LASSO coefficient path and cross-validation plots
##
## Input:
## result/Module05_CandidateGene_Expression.csv
## result/TCGA_train_survival.csv
##
## Output:
## result/Module06_LASSO_Coefficients_lambda_min.csv
## result/Module06_LASSO_Coefficients_lambda_1se.csv
## result/Module06_LASSO_SelectedGenes_lambda_min.csv
## result/Module06_LASSO_SelectedGenes_lambda_1se.csv
## figure/Module06_LASSO_CoefficientPath.pdf
## figure/Module06_LASSO_CV.pdf
## object/Module06_LASSO.RData
##
############################################################

rm(list = ls())
options(stringsAsFactors = FALSE)

cat("============================================\n")
cat("Module06: LASSO Cox regression\n")
cat("============================================\n\n")

############################################################
## 1. Packages
############################################################

cran_pkgs <- c(
  "data.table",
  "dplyr",
  "survival",
  "glmnet",
  "ggplot2"
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

expr_file <- "C:/Users/Xiong/Desktop/AML/result/Module05_CandidateGene_Expression.csv"
surv_file <- "C:/Users/Xiong/Desktop/AML/result/TCGA_train_survival.csv"

if(!file.exists(expr_file)){
  stop("Cannot find result/Module05_CandidateGene_Expression.csv")
}

if(!file.exists(surv_file)){
  stop("Cannot find result/TCGA_train_survival.csv")
}

############################################################
## 4. Read candidate gene expression matrix
############################################################

cat("Reading candidate gene expression matrix...\n")

expr <- data.table::fread(
  expr_file,
  data.table = FALSE,
  check.names = FALSE
)

rownames(expr) <- expr[, 1]
expr <- expr[, -1, drop = FALSE]

expr <- as.matrix(expr)
mode(expr) <- "numeric"

cat("Candidate expression dimension:\n")
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
  stop("Too few matched samples for LASSO Cox.")
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

cat("Final samples for LASSO:\n")
print(ncol(expr))

cat("Event table:\n")
print(table(survival_data$OS))

############################################################
## 8. Remove zero-variance genes
############################################################

gene_sd <- apply(expr, 1, sd, na.rm = TRUE)

expr <- expr[gene_sd > 0, , drop = FALSE]

cat("Genes retained after removing zero-variance genes:\n")
print(nrow(expr))

if(nrow(expr) < 2){
  stop("Too few genes for LASSO.")
}

############################################################
## 9. Prepare LASSO input
############################################################

## glmnet requires:
## x = samples × genes
## y = survival object

x <- t(expr)

y <- survival::Surv(
  time = survival_data$OS.time,
  event = survival_data$OS
)

############################################################
## 10. LASSO Cox regression
############################################################

cat("Running LASSO Cox regression...\n")

set.seed(12345)

cvfit <- glmnet::cv.glmnet(
  x = x,
  y = y,
  family = "cox",
  alpha = 1,
  nfolds = 10,
  type.measure = "deviance"
)

cat("LASSO finished.\n")

cat("lambda.min:\n")
print(cvfit$lambda.min)

cat("lambda.1se:\n")
print(cvfit$lambda.1se)

############################################################
## 11. Fit final glmnet model
############################################################

fit <- glmnet::glmnet(
  x = x,
  y = y,
  family = "cox",
  alpha = 1
)

############################################################
## 12. Extract coefficients: lambda.min
############################################################

coef_min <- coef(
  cvfit,
  s = "lambda.min"
)

coef_min_mat <- as.matrix(coef_min)

coef_min_df <- data.frame(
  Gene = rownames(coef_min_mat),
  Coefficient = as.numeric(coef_min_mat[, 1]),
  stringsAsFactors = FALSE
)

coef_min_df <- coef_min_df %>%
  dplyr::filter(Coefficient != 0) %>%
  dplyr::arrange(desc(abs(Coefficient)))

selected_genes_min <- coef_min_df$Gene

cat("Selected genes by lambda.min:\n")
print(selected_genes_min)

cat("Number of genes selected by lambda.min:\n")
print(length(selected_genes_min))

############################################################
## 13. Extract coefficients: lambda.1se
############################################################

coef_1se <- coef(
  cvfit,
  s = "lambda.1se"
)

coef_1se_mat <- as.matrix(coef_1se)

coef_1se_df <- data.frame(
  Gene = rownames(coef_1se_mat),
  Coefficient = as.numeric(coef_1se_mat[, 1]),
  stringsAsFactors = FALSE
)

coef_1se_df <- coef_1se_df %>%
  dplyr::filter(Coefficient != 0) %>%
  dplyr::arrange(desc(abs(Coefficient)))

selected_genes_1se <- coef_1se_df$Gene

cat("Selected genes by lambda.1se:\n")
print(selected_genes_1se)

cat("Number of genes selected by lambda.1se:\n")
print(length(selected_genes_1se))

############################################################
## 14. Save coefficient tables
############################################################

write.csv(
  coef_min_df,
  "C:/Users/Xiong/Desktop/AML/result/Module06_LASSO_Coefficients_lambda_min.csv",
  row.names = FALSE
)

write.csv(
  coef_1se_df,
  "C:/Users/Xiong/Desktop/AML/result/Module06_LASSO_Coefficients_lambda_1se.csv",
  row.names = FALSE
)

write.csv(
  data.frame(Gene = selected_genes_min),
  "C:/Users/Xiong/Desktop/AML/result/Module06_LASSO_SelectedGenes_lambda_min.csv",
  row.names = FALSE
)

write.csv(
  data.frame(Gene = selected_genes_1se),
  "C:/Users/Xiong/Desktop/AML/result/Module06_LASSO_SelectedGenes_lambda_1se.csv",
  row.names = FALSE
)

############################################################
## 15. Save selected expression matrices
############################################################

if(length(selected_genes_min) > 0){
  
  lasso_expr_min <- expr[
    rownames(expr) %in% selected_genes_min,
    ,
    drop = FALSE
  ]
  
  write.csv(
    lasso_expr_min,
    "C:/Users/Xiong/Desktop/AML/result/Module06_LASSO_Expression_lambda_min.csv",
    quote = FALSE
  )
  
}else{
  
  lasso_expr_min <- NULL
  
}

if(length(selected_genes_1se) > 0){
  
  lasso_expr_1se <- expr[
    rownames(expr) %in% selected_genes_1se,
    ,
    drop = FALSE
  ]
  
  write.csv(
    lasso_expr_1se,
    "C:/Users/Xiong/Desktop/AML/result/Module06_LASSO_Expression_lambda_1se.csv",
    quote = FALSE
  )
  
}else{
  
  lasso_expr_1se <- NULL
  
}

############################################################
## 16. LASSO cross-validation plot
############################################################

cat("Drawing LASSO CV plot...\n")

pdf(
  "C:/Users/Xiong/Desktop/AML/figure/Module06_LASSO_CV.pdf",
  width = 7,
  height = 6
)

plot(cvfit)

abline(
  v = log(cvfit$lambda.min),
  lty = 2,
  col = "red"
)

abline(
  v = log(cvfit$lambda.1se),
  lty = 2,
  col = "blue"
)

legend(
  "topright",
  legend = c("lambda.min", "lambda.1se"),
  col = c("red", "blue"),
  lty = 2,
  bty = "n"
)

dev.off()

png(
  "C:/Users/Xiong/Desktop/AML/figure/Module06_LASSO_CV.png",
  width = 1800,
  height = 1500,
  res = 300
)

plot(cvfit)

abline(
  v = log(cvfit$lambda.min),
  lty = 2,
  col = "red"
)

abline(
  v = log(cvfit$lambda.1se),
  lty = 2,
  col = "blue"
)

legend(
  "topright",
  legend = c("lambda.min", "lambda.1se"),
  col = c("red", "blue"),
  lty = 2,
  bty = "n"
)

dev.off()

############################################################
## 17. LASSO coefficient path plot
############################################################

cat("Drawing LASSO coefficient path plot...\n")

pdf(
  "C:/Users/Xiong/Desktop/AML/figure/Module06_LASSO_CoefficientPath.pdf",
  width = 7,
  height = 6
)

plot(
  fit,
  xvar = "lambda",
  label = FALSE
)

abline(
  v = log(cvfit$lambda.min),
  lty = 2,
  col = "red"
)

abline(
  v = log(cvfit$lambda.1se),
  lty = 2,
  col = "blue"
)

legend(
  "topright",
  legend = c("lambda.min", "lambda.1se"),
  col = c("red", "blue"),
  lty = 2,
  bty = "n"
)

dev.off()

png(
  "C:/Users/Xiong/Desktop/AML/figure/Module06_LASSO_CoefficientPath.png",
  width = 1800,
  height = 1500,
  res = 300
)

plot(
  fit,
  xvar = "lambda",
  label = FALSE
)

abline(
  v = log(cvfit$lambda.min),
  lty = 2,
  col = "red"
)

abline(
  v = log(cvfit$lambda.1se),
  lty = 2,
  col = "blue"
)

legend(
  "topright",
  legend = c("lambda.min", "lambda.1se"),
  col = c("red", "blue"),
  lty = 2,
  bty = "n"
)

dev.off()

############################################################
## 18. Coefficient barplot for lambda.min
############################################################

if(nrow(coef_min_df) > 0){
  
  coef_plot_df <- coef_min_df %>%
    dplyr::mutate(
      Direction = ifelse(Coefficient > 0, "Risk", "Protective"),
      Gene = factor(Gene, levels = Gene[order(Coefficient)])
    )
  
  p_coef <- ggplot(
    coef_plot_df,
    aes(
      x = Gene,
      y = Coefficient,
      fill = Direction
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
      title = "LASSO coefficients at lambda.min",
      x = "",
      y = "Coefficient"
    )
  
  ggsave(
    "C:/Users/Xiong/Desktop/AML/figure/Module06_LASSO_CoefficientBar_lambda_min.pdf",
    p_coef,
    width = 7,
    height = max(5, nrow(coef_plot_df) * 0.35)
  )
  
  ggsave(
    "C:/Users/Xiong/Desktop/AML/figure/Module06_LASSO_CoefficientBar_lambda_min.png",
    p_coef,
    width = 7,
    height = max(5, nrow(coef_plot_df) * 0.35),
    dpi = 300
  )
  
}

############################################################
## 19. Save RData object
############################################################

save(
  expr,
  survival_data,
  x,
  y,
  cvfit,
  fit,
  coef_min_df,
  coef_1se_df,
  selected_genes_min,
  selected_genes_1se,
  lasso_expr_min,
  lasso_expr_1se,
  file = "C:/Users/Xiong/Desktop/AML/object/Module06_LASSO.RData"
)

############################################################
## 20. Generate report
############################################################

report_file <- "C:/Users/Xiong/Desktop/AML/result/Module06_Report.txt"

sink(report_file)

cat("Module06 Report: LASSO Cox regression\n")
cat("=====================================\n\n")

cat("Project:\n")
cat("Lysine-associated gene signature predicts prognosis in AML\n\n")

cat("Input files:\n")
cat(expr_file, "\n")
cat(surv_file, "\n\n")

cat("Input candidate expression dimension:\n")
print(dim(expr))

cat("\nSurvival data dimension:\n")
print(dim(survival_data))

cat("\nEvent table:\n")
print(table(survival_data$OS))

cat("\nlambda.min:\n")
print(cvfit$lambda.min)

cat("\nlambda.1se:\n")
print(cvfit$lambda.1se)

cat("\nGenes selected by lambda.min:\n")
print(selected_genes_min)

cat("\nNumber of genes selected by lambda.min:\n")
print(length(selected_genes_min))

cat("\nGenes selected by lambda.1se:\n")
print(selected_genes_1se)

cat("\nNumber of genes selected by lambda.1se:\n")
print(length(selected_genes_1se))

cat("\nCoefficient table at lambda.min:\n")
print(coef_min_df)

cat("\nCoefficient table at lambda.1se:\n")
print(coef_1se_df)

cat("\nGenerated files:\n")
cat("result/Module06_LASSO_Coefficients_lambda_min.csv\n")
cat("result/Module06_LASSO_Coefficients_lambda_1se.csv\n")
cat("result/Module06_LASSO_SelectedGenes_lambda_min.csv\n")
cat("result/Module06_LASSO_SelectedGenes_lambda_1se.csv\n")
cat("result/Module06_LASSO_Expression_lambda_min.csv\n")
cat("result/Module06_LASSO_Expression_lambda_1se.csv\n")
cat("figure/Module06_LASSO_CV.pdf\n")
cat("figure/Module06_LASSO_CoefficientPath.pdf\n")
cat("figure/Module06_LASSO_CoefficientBar_lambda_min.pdf\n")
cat("object/Module06_LASSO.RData\n")

cat("\nSession information:\n")
print(sessionInfo())

sink()

cat("\nText report saved:\n")
cat(report_file, "\n")

############################################################
## 21. Finish
############################################################

cat("\n============================================\n")
cat("Module06 finished successfully.\n")
cat("============================================\n")

