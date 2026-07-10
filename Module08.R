############################################################
##
## Module08_RiskScore_TCGA.R
##
## Project:
## Lysine-associated gene signature predicts prognosis in AML
##
## Purpose:
## 1. Calculate TCGA risk score using final genes from Module07
## 2. Divide patients into high- and low-risk groups
## 3. Perform Kaplan-Meier survival analysis
## 4. Perform time-dependent ROC analysis
## 5. Generate risk plot, survival status plot, and heatmap
##
## Input:
## result/Module07_StepwiseCox_FinalModel.csv
## result/TCGA_train_expression.csv
## result/TCGA_train_survival.csv
##
## Output:
## result/Module08_TCGA_RiskScore.csv
## result/Module08_TCGA_RiskGroup_Summary.csv
## figure/Module08_TCGA_KM.pdf
## figure/Module08_TCGA_TimeROC.pdf
## figure/Module08_TCGA_RiskScore_Distribution.pdf
## figure/Module08_TCGA_SurvivalStatus.pdf
## figure/Module08_TCGA_ModelGene_Heatmap.pdf
## object/Module08_RiskScore_TCGA.RData
##
############################################################

rm(list = ls())
options(stringsAsFactors = FALSE)

cat("============================================\n")
cat("Module08: TCGA Risk Score Construction\n")
cat("============================================\n\n")

############################################################
## 1. Packages
############################################################

cran_pkgs <- c(
  "data.table",
  "dplyr",
  "survival",
  "survminer",
  "timeROC",
  "ggplot2",
  "pheatmap",
  "RColorBrewer"
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

model_file <- "C:/Users/Xiong/Desktop/AML/result/Module07_StepwiseCox_FinalModel.csv"
expr_file <- "C:/Users/Xiong/Desktop/AML/result/TCGA_train_expression.csv"
surv_file <- "C:/Users/Xiong/Desktop/AML/result/TCGA_train_survival.csv"

if(!file.exists(model_file)){
  stop("Cannot find result/Module07_StepwiseCox_FinalModel.csv")
}

if(!file.exists(expr_file)){
  stop("Cannot find result/TCGA_train_expression.csv")
}

if(!file.exists(surv_file)){
  stop("Cannot find result/TCGA_train_survival.csv")
}

############################################################
## 4. Read final model
############################################################

cat("Reading final model from Module07...\n")

model <- data.table::fread(
  model_file,
  data.table = FALSE
)

cat("Final model genes:\n")
print(model$Gene)

required_model_cols <- c("Gene", "Coef")

missing_model_cols <- setdiff(required_model_cols, colnames(model))

if(length(missing_model_cols) > 0){
  stop(
    paste0(
      "Missing columns in model file: ",
      paste(missing_model_cols, collapse = ", ")
    )
  )
}

model_genes <- model$Gene
coef_vec <- model$Coef
names(coef_vec) <- model_genes

############################################################
## 5. Read TCGA expression matrix
############################################################

cat("Reading TCGA expression matrix...\n")

expr <- data.table::fread(
  expr_file,
  data.table = FALSE,
  check.names = FALSE
)

rownames(expr) <- expr[, 1]
expr <- expr[, -1, drop = FALSE]

expr <- as.matrix(expr)
mode(expr) <- "numeric"

cat("TCGA expression dimension:\n")
print(dim(expr))

############################################################
## 6. Read survival data
############################################################

cat("Reading survival data...\n")

survival_data <- data.table::fread(
  surv_file,
  data.table = FALSE
)

required_surv_cols <- c("sample", "OS.time", "OS")

missing_surv_cols <- setdiff(required_surv_cols, colnames(survival_data))

if(length(missing_surv_cols) > 0){
  stop(
    paste0(
      "Missing survival columns: ",
      paste(missing_surv_cols, collapse = ", ")
    )
  )
}

survival_data$sample <- substr(survival_data$sample, 1, 16)
survival_data$OS.time <- as.numeric(survival_data$OS.time)
survival_data$OS <- as.numeric(survival_data$OS)

############################################################
## 7. Check model genes in expression matrix
############################################################

missing_genes <- setdiff(model_genes, rownames(expr))

if(length(missing_genes) > 0){
  stop(
    paste0(
      "These final model genes are missing in TCGA expression matrix: ",
      paste(missing_genes, collapse = ", ")
    )
  )
}

model_expr <- expr[
  model_genes,
  ,
  drop = FALSE
]

############################################################
## 8. Match samples
############################################################

cat("Matching expression samples with survival data...\n")

colnames(model_expr) <- substr(colnames(model_expr), 1, 16)

common_samples <- intersect(
  colnames(model_expr),
  survival_data$sample
)

cat("Matched samples:\n")
print(length(common_samples))

if(length(common_samples) < 30){
  stop("Too few matched samples for risk score analysis.")
}

model_expr <- model_expr[, common_samples, drop = FALSE]

survival_data <- survival_data[
  match(common_samples, survival_data$sample),
]

if(!all(colnames(model_expr) == survival_data$sample)){
  stop("Expression samples and survival samples are not matched.")
}

############################################################
## 9. Remove missing survival data
############################################################

keep <- !is.na(survival_data$OS.time) &
  !is.na(survival_data$OS) &
  survival_data$OS.time > 0

model_expr <- model_expr[, keep, drop = FALSE]
survival_data <- survival_data[keep, ]

cat("Final TCGA samples for risk model:\n")
print(ncol(model_expr))

cat("Event table:\n")
print(table(survival_data$OS))

############################################################
## 10. Calculate risk score
############################################################

cat("Calculating risk score...\n")

model_expr_ordered <- model_expr[
  names(coef_vec),
  ,
  drop = FALSE
]

risk_score <- as.numeric(
  t(model_expr_ordered) %*% coef_vec
)

risk_df <- data.frame(
  sample = colnames(model_expr_ordered),
  RiskScore = risk_score,
  OS.time = survival_data$OS.time,
  OS = survival_data$OS,
  stringsAsFactors = FALSE
)

############################################################
## 11. Define high- and low-risk groups
############################################################

cutoff <- median(
  risk_df$RiskScore,
  na.rm = TRUE
)

risk_df$RiskGroup <- ifelse(
  risk_df$RiskScore >= cutoff,
  "High",
  "Low"
)

risk_df$RiskGroup <- factor(
  risk_df$RiskGroup,
  levels = c("Low", "High")
)

cat("Risk score cutoff used:\n")
print(cutoff)

cat("Risk group table:\n")
print(table(risk_df$RiskGroup))

write.csv(
  risk_df,
  "C:/Users/Xiong/Desktop/AML/result/Module08_TCGA_RiskScore.csv",
  row.names = FALSE
)

risk_summary <- data.frame(
  Cohort = "TCGA",
  Samples = nrow(risk_df),
  Cutoff = cutoff,
  LowRisk = sum(risk_df$RiskGroup == "Low"),
  HighRisk = sum(risk_df$RiskGroup == "High"),
  Events = sum(risk_df$OS == 1),
  Censored = sum(risk_df$OS == 0)
)

write.csv(
  risk_summary,
  "C:/Users/Xiong/Desktop/AML/result/Module08_TCGA_RiskGroup_Summary.csv",
  row.names = FALSE
)

############################################################
## 12. Kaplan-Meier survival analysis
############################################################

cat("Performing Kaplan-Meier analysis...\n")

fit_km <- survival::survfit(
  survival::Surv(OS.time, OS) ~ RiskGroup,
  data = risk_df
)

cox_group <- survival::coxph(
  survival::Surv(OS.time, OS) ~ RiskGroup,
  data = risk_df
)

cox_group_summary <- summary(cox_group)

group_hr <- cox_group_summary$conf.int[1, "exp(coef)"]
group_lower <- cox_group_summary$conf.int[1, "lower .95"]
group_upper <- cox_group_summary$conf.int[1, "upper .95"]
group_p <- cox_group_summary$coefficients[1, "Pr(>|z|)"]

km_plot <- survminer::ggsurvplot(
  fit_km,
  data = risk_df,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  palette = c("#4575B4", "#D73027"),
  legend.title = "Risk group",
  legend.labs = c("Low risk", "High risk"),
  xlab = "Time (days)",
  ylab = "Overall survival probability",
  title = "TCGA-LAML overall survival by risk group",
  risk.table.height = 0.25,
  ggtheme = theme_bw()
)

pdf(
  "C:/Users/Xiong/Desktop/AML/figure/Module08_TCGA_KM.pdf",
  width = 7,
  height = 7
)

print(km_plot)

dev.off()

png(
  "C:/Users/Xiong/Desktop/AML/figure/Module08_TCGA_KM.png",
  width = 2100,
  height = 2100,
  res = 300
)

print(km_plot)

dev.off()

############################################################
## 13. Time-dependent ROC
############################################################

cat("Performing time-dependent ROC analysis...\n")

## TCGA-LAML survival time is usually in days.
## 1-year = 365 days, 3-year = 1095 days, 5-year = 1825 days.

roc_times <- c(
  365,
  1095,
  1825
)

roc_obj <- timeROC::timeROC(
  T = risk_df$OS.time,
  delta = risk_df$OS,
  marker = risk_df$RiskScore,
  cause = 1,
  times = roc_times,
  iid = TRUE
)

auc_df <- data.frame(
  Time = c("1-year", "3-year", "5-year"),
  Days = roc_times,
  AUC = as.numeric(roc_obj$AUC)
)

write.csv(
  auc_df,
  "C:/Users/Xiong/Desktop/AML/result/Module08_TCGA_TimeROC_AUC.csv",
  row.names = FALSE
)

pdf(
  "C:/Users/Xiong/Desktop/AML/figure/Module08_TCGA_TimeROC.pdf",
  width = 7,
  height = 6
)

plot(
  roc_obj,
  time = 365,
  col = "#D73027",
  title = FALSE
)

plot(
  roc_obj,
  time = 1095,
  add = TRUE,
  col = "#4575B4"
)

plot(
  roc_obj,
  time = 1825,
  add = TRUE,
  col = "#1A9850"
)

legend(
  "bottomright",
  legend = paste0(
    auc_df$Time,
    " AUC = ",
    sprintf("%.3f", auc_df$AUC)
  ),
  col = c("#D73027", "#4575B4", "#1A9850"),
  lwd = 2,
  bty = "n"
)

title(
  main = "Time-dependent ROC curve in TCGA-LAML"
)

dev.off()

png(
  "C:/Users/Xiong/Desktop/AML/figure/Module08_TCGA_TimeROC.png",
  width = 2100,
  height = 1800,
  res = 300
)

plot(
  roc_obj,
  time = 365,
  col = "#D73027",
  title = FALSE
)

plot(
  roc_obj,
  time = 1095,
  add = TRUE,
  col = "#4575B4"
)

plot(
  roc_obj,
  time = 1825,
  add = TRUE,
  col = "#1A9850"
)

legend(
  "bottomright",
  legend = paste0(
    auc_df$Time,
    " AUC = ",
    sprintf("%.3f", auc_df$AUC)
  ),
  col = c("#D73027", "#4575B4", "#1A9850"),
  lwd = 2,
  bty = "n"
)

title(
  main = "Time-dependent ROC curve in TCGA-LAML"
)

dev.off()

############################################################
## 14. Risk score distribution plot
############################################################

cat("Drawing risk score distribution plot...\n")

risk_order <- risk_df %>%
  dplyr::arrange(RiskScore)

risk_order$Index <- 1:nrow(risk_order)

p_risk <- ggplot(
  risk_order,
  aes(
    x = Index,
    y = RiskScore,
    color = RiskGroup
  )
) +
  geom_point(size = 2) +
  geom_hline(
    yintercept = cutoff,
    linetype = "dashed"
  ) +
  scale_color_manual(
    values = c(
      "Low" = "#4575B4",
      "High" = "#D73027"
    )
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_blank()
  ) +
  labs(
    title = "Risk score distribution in TCGA-LAML",
    x = "Patients ranked by risk score",
    y = "Risk score"
  )

ggsave(
  "C:/Users/Xiong/Desktop/AML/figure/Module08_TCGA_RiskScore_Distribution.pdf",
  p_risk,
  width = 7,
  height = 5
)

ggsave(
  "C:/Users/Xiong/Desktop/AML/figure/Module08_TCGA_RiskScore_Distribution.png",
  p_risk,
  width = 7,
  height = 5,
  dpi = 300
)

############################################################
## 15. Survival status plot
############################################################

cat("Drawing survival status plot...\n")

risk_order$Status <- ifelse(
  risk_order$OS == 1,
  "Dead",
  "Alive/Censored"
)

p_status <- ggplot(
  risk_order,
  aes(
    x = Index,
    y = OS.time,
    color = Status
  )
) +
  geom_point(size = 2) +
  scale_color_manual(
    values = c(
      "Dead" = "#D73027",
      "Alive/Censored" = "#4575B4"
    )
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_blank()
  ) +
  labs(
    title = "Survival status distribution in TCGA-LAML",
    x = "Patients ranked by risk score",
    y = "Survival time (days)"
  )

ggsave(
  "C:/Users/Xiong/Desktop/AML/figure/Module08_TCGA_SurvivalStatus.pdf",
  p_status,
  width = 7,
  height = 5
)

ggsave(
  "C:/Users/Xiong/Desktop/AML/figure/Module08_TCGA_SurvivalStatus.png",
  p_status,
  width = 7,
  height = 5,
  dpi = 300
)

############################################################
## 16. Heatmap of final model gene expression
############################################################

cat("Drawing final model gene heatmap...\n")

heat_expr <- model_expr_ordered[
  ,
  risk_order$sample,
  drop = FALSE
]

## Z-score by gene
heat_expr_z <- t(
  scale(
    t(heat_expr)
  )
)

annotation_col <- data.frame(
  RiskGroup = risk_order$RiskGroup
)

rownames(annotation_col) <- risk_order$sample

ann_colors <- list(
  RiskGroup = c(
    Low = "#4575B4",
    High = "#D73027"
  )
)

pdf(
  "C:/Users/Xiong/Desktop/AML/figure/Module08_TCGA_ModelGene_Heatmap.pdf",
  width = 9,
  height = 5
)

pheatmap::pheatmap(
  heat_expr_z,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_colnames = FALSE,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  color = colorRampPalette(
    rev(
      RColorBrewer::brewer.pal(
        n = 11,
        name = "RdBu"
      )
    )
  )(100),
  main = "Final model gene expression in TCGA-LAML"
)

dev.off()

png(
  "C:/Users/Xiong/Desktop/AML/figure/Module08_TCGA_ModelGene_Heatmap.png",
  width = 2700,
  height = 1500,
  res = 300
)

pheatmap::pheatmap(
  heat_expr_z,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_colnames = FALSE,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  color = colorRampPalette(
    rev(
      RColorBrewer::brewer.pal(
        n = 11,
        name = "RdBu"
      )
    )
  )(100),
  main = "Final model gene expression in TCGA-LAML"
)

dev.off()

############################################################
## 17. Save model performance summary
############################################################

performance_summary <- data.frame(
  Cohort = "TCGA",
  Samples = nrow(risk_df),
  Cutoff_method = "Median",
  Cutoff = cutoff,
  Low_risk = sum(risk_df$RiskGroup == "Low"),
  High_risk = sum(risk_df$RiskGroup == "High"),
  Event_dead = sum(risk_df$OS == 1),
  Event_censored = sum(risk_df$OS == 0),
  HR_high_vs_low = group_hr,
  HR_lower95CI = group_lower,
  HR_upper95CI = group_upper,
  Cox_Pvalue = group_p,
  AUC_1year = auc_df$AUC[auc_df$Time == "1-year"],
  AUC_3year = auc_df$AUC[auc_df$Time == "3-year"],
  AUC_5year = auc_df$AUC[auc_df$Time == "5-year"]
)

write.csv(
  performance_summary,
  "C:/Users/Xiong/Desktop/AML/result/Module08_TCGA_ModelPerformance.csv",
  row.names = FALSE
)

############################################################
## 18. Save RData object
############################################################

save(
  model,
  model_genes,
  coef_vec,
  expr,
  model_expr,
  model_expr_ordered,
  survival_data,
  risk_df,
  risk_order,
  cutoff,
  fit_km,
  cox_group,
  roc_obj,
  auc_df,
  performance_summary,
  file = "C:/Users/Xiong/Desktop/AML/object/Module08_RiskScore_TCGA.RData"
)

############################################################
## 19. Generate text report
############################################################

report_file <- "C:/Users/Xiong/Desktop/AML/result/Module08_Report.txt"

sink(report_file)

cat("Module08 Report: TCGA risk score construction\n")
cat("============================================\n\n")

cat("Project:\n")
cat("Lysine-associated gene signature predicts prognosis in AML\n\n")

cat("Input files:\n")
cat(model_file, "\n")
cat(expr_file, "\n")
cat(surv_file, "\n\n")

cat("Final model genes:\n")
print(model_genes)

cat("\nModel coefficients:\n")
print(model[, c("Gene", "Coef", "HR", "Pvalue", "RiskType")])

cat("\nRisk score cutoff method:\n")
cat("Median risk score\n\n")

cat("Risk score cutoff:\n")
print(cutoff)

cat("\nRisk group table:\n")
print(table(risk_df$RiskGroup))

cat("\nSurvival event table:\n")
print(table(risk_df$OS))

cat("\nHigh vs Low risk Cox result:\n")
cat("HR =", group_hr, "\n")
cat("95% CI =", group_lower, "-", group_upper, "\n")
cat("P value =", group_p, "\n\n")

cat("Time-dependent ROC AUC:\n")
print(auc_df)

cat("\nPerformance summary:\n")
print(performance_summary)

cat("\nGenerated files:\n")
cat("result/Module08_TCGA_RiskScore.csv\n")
cat("result/Module08_TCGA_RiskGroup_Summary.csv\n")
cat("result/Module08_TCGA_TimeROC_AUC.csv\n")
cat("result/Module08_TCGA_ModelPerformance.csv\n")
cat("figure/Module08_TCGA_KM.pdf\n")
cat("figure/Module08_TCGA_TimeROC.pdf\n")
cat("figure/Module08_TCGA_RiskScore_Distribution.pdf\n")
cat("figure/Module08_TCGA_SurvivalStatus.pdf\n")
cat("figure/Module08_TCGA_ModelGene_Heatmap.pdf\n")
cat("object/Module08_RiskScore_TCGA.RData\n")

cat("\nSession information:\n")
print(sessionInfo())

sink()

cat("\nText report saved:\n")
cat(report_file, "\n")

############################################################
## 20. Finish
############################################################

cat("\n============================================\n")
cat("Module08 finished successfully.\n")
cat("============================================\n")

