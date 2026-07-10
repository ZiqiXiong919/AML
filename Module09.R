############################################################
##
## Module09_GEO_ExternalValidation.R
##
## Project:
## Lysine-associated gene signature predicts prognosis in AML
##
## Purpose:
## 1. Validate TCGA-derived risk score model in external GEO cohorts
## 2. Cohorts:
##    - GSE37642
##    - GSE12417
## 3. Calculate risk score using final Module07 coefficients
## 4. Perform KM survival analysis
## 5. Perform time-dependent ROC analysis
## 6. Generate risk score plot, survival status plot, heatmap
##
## Important update:
## Corrected GSE12417 OS.time extraction.
## Avoids extracting contact_zip.postal_code = 81377 as survival time.
##
############################################################

rm(list = ls())
options(stringsAsFactors = FALSE)

cat("============================================\n")
cat("Module09: GEO external validation\n")
cat("============================================\n\n")

############################################################
## 0. Project directory
############################################################

project_dir <- "C:/Users/Xiong/Desktop/AML"

result_dir <- file.path(project_dir, "result")
figure_dir <- file.path(project_dir, "figure")
object_dir <- file.path(project_dir, "object")

dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(object_dir, showWarnings = FALSE, recursive = TRUE)

############################################################
## 1. Packages
############################################################

cran_pkgs <- c(
  "data.table",
  "dplyr",
  "stringr",
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
## 2. Input files
############################################################

model_file <- file.path(
  result_dir,
  "Module07_StepwiseCox_FinalModel.csv"
)

tcga_risk_file <- file.path(
  result_dir,
  "Module08_TCGA_RiskScore.csv"
)

gse37642_expr_file <- file.path(
  result_dir,
  "GSE37642_lysine_expression.csv"
)

gse12417_expr_file <- file.path(
  result_dir,
  "GSE12417_lysine_expression.csv"
)

gse37642_pheno_file <- file.path(
  result_dir,
  "GSE37642_validation_pheno.csv"
)

gse12417_pheno_file <- file.path(
  result_dir,
  "GSE12417_validation_pheno.csv"
)

input_files <- c(
  model_file,
  tcga_risk_file,
  gse37642_expr_file,
  gse12417_expr_file,
  gse37642_pheno_file,
  gse12417_pheno_file
)

missing_files <- input_files[!file.exists(input_files)]

if(length(missing_files) > 0){
  stop(
    paste0(
      "Missing input files:\n",
      paste(missing_files, collapse = "\n")
    )
  )
}

cat("All input files found.\n\n")

############################################################
## 3. Helper function: read expression CSV
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
## 4. Survival extraction for GSE37642
############################################################

clean_gse37642_survival <- function(pheno){
  
  cat("\nCleaning GSE37642 survival information...\n")
  
  sample <- if("geo_accession" %in% colnames(pheno)){
    pheno$geo_accession
  }else{
    rownames(pheno)
  }
  
  all_text <- apply(
    pheno,
    1,
    function(x){
      paste(x, collapse = " ; ")
    }
  )
  
  ##########################################################
  ## Extract OS time
  ##########################################################
  
  os_time_raw <- stringr::str_extract(
    all_text,
    "(?i)(OS|overall survival|survival)[^0-9]{0,50}[0-9]+\\.?[0-9]*"
  )
  
  os_time <- stringr::str_extract(
    os_time_raw,
    "[0-9]+\\.?[0-9]*"
  )
  
  os_time <- as.numeric(os_time)
  
  ##########################################################
  ## Extract OS status
  ##########################################################
  
  os_status <- rep(NA, length(all_text))
  
  os_status[
    grepl("dead|deceased|death|event|relapse", all_text, ignore.case = TRUE)
  ] <- 1
  
  os_status[
    grepl("alive|living|censored|censor", all_text, ignore.case = TRUE)
  ] <- 0
  
  status_raw <- stringr::str_extract(
    all_text,
    "(?i)(status|vital|event)[^0-9]{0,30}[01]"
  )
  
  status_num <- stringr::str_extract(
    status_raw,
    "[01]$"
  )
  
  status_num <- as.numeric(status_num)
  
  os_status[!is.na(status_num)] <- status_num[!is.na(status_num)]
  
  ##########################################################
  ## Convert time if clearly years/months
  ##########################################################
  
  if(all(os_time < 100, na.rm = TRUE)){
    if(any(grepl("year", all_text, ignore.case = TRUE))){
      os_time <- os_time * 365
    }else if(any(grepl("month", all_text, ignore.case = TRUE))){
      os_time <- os_time * 30.44
    }
  }
  
  surv <- data.frame(
    sample = sample,
    OS.time = os_time,
    OS = os_status,
    stringsAsFactors = FALSE
  )
  
  surv <- surv[
    !is.na(surv$OS.time) &
      !is.na(surv$OS) &
      surv$OS.time > 0,
  ]
  
  cat("GSE37642 survival extracted:\n")
  print(dim(surv))
  print(table(surv$OS))
  print(head(surv))
  
  write.csv(
    surv,
    file.path(result_dir, "Module09_GSE37642_ExtractedSurvival.csv"),
    row.names = FALSE
  )
  
  if(nrow(surv) == 0){
    stop("GSE37642 survival extraction failed.")
  }
  
  if(length(unique(surv$OS)) < 2){
    stop("GSE37642 survival extraction failed: only one survival status detected.")
  }
  
  return(surv)
}

############################################################
## 5. Corrected survival extraction for GSE12417
############################################################

clean_gse12417_survival <- function(pheno){
  
  cat("\nCleaning GSE12417 survival information...\n")
  
  ##########################################################
  ## Sample IDs
  ##########################################################
  
  sample <- if("geo_accession" %in% colnames(pheno)){
    pheno$geo_accession
  }else{
    rownames(pheno)
  }
  
  ##########################################################
  ## Select real clinical OS columns
  ##
  ## Correct columns usually look like:
  ## AML.normal.karyotype...OS...1280.days..status..0.alive.1.dead..ch1
  ##
  ## Required:
  ## - contains OS
  ## - contains days
  ## - contains status
  ##
  ## Exclude contact/postal/zip columns.
  ##########################################################
  
  clinical_cols <- colnames(pheno)[
    grepl("OS", colnames(pheno), ignore.case = TRUE) &
      grepl("days", colnames(pheno), ignore.case = TRUE) &
      grepl("status", colnames(pheno), ignore.case = TRUE)
  ]
  
  clinical_cols <- clinical_cols[
    !grepl(
      "contact|postal|zip|address|phone|fax|email|submission|update|public",
      clinical_cols,
      ignore.case = TRUE
    )
  ]
  
  cat("\nCandidate true clinical OS columns:\n")
  print(clinical_cols)
  
  if(length(clinical_cols) == 0){
    stop(
      "No true GSE12417 clinical OS columns found. Please inspect colnames(gse12417_pheno)."
    )
  }
  
  write.csv(
    data.frame(Clinical_OS_Column = clinical_cols),
    file.path(result_dir, "Module09_GSE12417_OS_ClinicalColumns.csv"),
    row.names = FALSE
  )
  
  ##########################################################
  ## Extract row-specific OS.time and OS
  ##
  ## For each sample:
  ## - one clinical column usually has cell value 0 or 1
  ## - column name contains that sample's OS time
  ## - cell value contains that sample's status
  ##########################################################
  
  os_time <- rep(NA_real_, nrow(pheno))
  os_status <- rep(NA_real_, nrow(pheno))
  source_column <- rep(NA_character_, nrow(pheno))
  source_value <- rep(NA_character_, nrow(pheno))
  
  for(i in seq_len(nrow(pheno))){
    
    row_values <- as.character(
      pheno[i, clinical_cols, drop = TRUE]
    )
    
    row_values <- trimws(row_values)
    
    active_idx <- which(
      !is.na(row_values) &
        row_values != "" &
        row_values != "NA"
    )
    
    if(length(active_idx) == 0){
      next
    }
    
    this_col <- clinical_cols[active_idx[1]]
    this_value <- row_values[active_idx[1]]
    
    source_column[i] <- this_col
    source_value[i] <- this_value
    
    ########################################################
    ## Extract OS.time from column name
    ## Example:
    ## OS...1280.days
    ########################################################
    
    time_raw <- stringr::str_extract(
      this_col,
      "(?i)OS[^0-9]{0,40}[0-9]+\\.?[0-9]*[^;]{0,30}days"
    )
    
    time_num <- stringr::str_extract(
      time_raw,
      "[0-9]+\\.?[0-9]*"
    )
    
    os_time[i] <- as.numeric(time_num)
    
    ########################################################
    ## Extract OS status from cell value
    ## 0 = alive
    ## 1 = dead
    ########################################################
    
    status_num <- suppressWarnings(
      as.numeric(this_value)
    )
    
    if(!is.na(status_num) && status_num %in% c(0, 1)){
      os_status[i] <- status_num
    }
    
  }
  
  ##########################################################
  ## Build survival table
  ##########################################################
  
  surv_all <- data.frame(
    sample = sample,
    OS.time = os_time,
    OS = os_status,
    SourceColumn = source_column,
    SourceValue = source_value,
    stringsAsFactors = FALSE
  )
  
  write.csv(
    surv_all,
    file.path(result_dir, "Module09_GSE12417_ExtractedSurvival_AllRows.csv"),
    row.names = FALSE
  )
  
  surv <- surv_all[
    !is.na(surv_all$OS.time) &
      !is.na(surv_all$OS) &
      surv_all$OS.time > 0,
  ]
  
  surv_simple <- surv[, c("sample", "OS.time", "OS")]
  
  write.csv(
    surv_simple,
    file.path(result_dir, "Module09_GSE12417_ExtractedSurvival.csv"),
    row.names = FALSE
  )
  
  ##########################################################
  ## QC output
  ##########################################################
  
  cat("\nGSE12417 survival extracted:\n")
  print(dim(surv_simple))
  
  cat("\nOS status table:\n")
  print(table(surv_simple$OS))
  
  cat("\nOS.time summary:\n")
  print(summary(surv_simple$OS.time))
  
  cat("\nNumber of unique OS.time values:\n")
  print(length(unique(surv_simple$OS.time)))
  
  cat("\nHead of extracted survival:\n")
  print(head(surv_simple))
  
  ##########################################################
  ## Safety checks
  ##########################################################
  
  if(nrow(surv_simple) == 0){
    stop("GSE12417 survival extraction failed: no valid survival records.")
  }
  
  if(length(unique(surv_simple$OS)) < 2){
    stop("GSE12417 survival extraction failed: only one survival status detected.")
  }
  
  if(length(unique(surv_simple$OS.time)) < 5){
    stop("GSE12417 survival extraction failed: OS.time values look incorrect.")
  }
  
  if(max(surv_simple$OS.time, na.rm = TRUE) > 30000){
    warning(
      "Some OS.time values are extremely large. Please inspect Module09_GSE12417_ExtractedSurvival_AllRows.csv"
    )
  }
  
  return(surv_simple)
}

############################################################
## 6. Generic validation function
############################################################

validate_geo <- function(
    cohort_name,
    expr,
    pheno,
    model,
    tcga_cutoff,
    cutoff_mode = c("TCGA", "Median")
){
  
  cutoff_mode <- match.arg(cutoff_mode)
  
  cat("\n--------------------------------------------\n")
  cat("Validating cohort: ", cohort_name, " | cutoff: ", cutoff_mode, "\n", sep = "")
  cat("--------------------------------------------\n")
  
  model_genes <- model$Gene
  coef_vec <- model$Coef
  names(coef_vec) <- model_genes
  
  ##########################################################
  ## Gene matching
  ##########################################################
  
  missing_genes <- setdiff(model_genes, rownames(expr))
  
  if(length(missing_genes) > 0){
    cat("Warning: missing model genes in ", cohort_name, ":\n", sep = "")
    print(missing_genes)
  }
  
  common_genes <- intersect(model_genes, rownames(expr))
  
  cat("Matched model genes:\n")
  print(length(common_genes))
  
  if(length(common_genes) < 3){
    stop(paste0("Too few model genes matched in ", cohort_name))
  }
  
  coef_vec <- coef_vec[common_genes]
  
  model_expr <- expr[
    common_genes,
    ,
    drop = FALSE
  ]
  
  model_expr <- model_expr[
    names(coef_vec),
    ,
    drop = FALSE
  ]
  
  ##########################################################
  ## Survival extraction
  ##########################################################
  
  if(cohort_name == "GSE12417"){
    
    surv <- clean_gse12417_survival(
      pheno = pheno
    )
    
  }else if(cohort_name == "GSE37642"){
    
    surv <- clean_gse37642_survival(
      pheno = pheno
    )
    
  }else{
    
    stop("Unknown cohort survival extraction rule.")
    
  }
  
  ##########################################################
  ## Sample matching
  ##########################################################
  
  colnames(model_expr) <- as.character(colnames(model_expr))
  surv$sample <- as.character(surv$sample)
  
  common_samples <- intersect(
    colnames(model_expr),
    surv$sample
  )
  
  cat("Matched samples:\n")
  print(length(common_samples))
  
  if(length(common_samples) < 20){
    stop(paste0("Too few matched samples in ", cohort_name))
  }
  
  model_expr <- model_expr[
    ,
    common_samples,
    drop = FALSE
  ]
  
  surv <- surv[
    match(common_samples, surv$sample),
  ]
  
  if(!all(colnames(model_expr) == surv$sample)){
    stop(paste0("Sample order mismatch in ", cohort_name))
  }
  
  surv$OS.time <- as.numeric(surv$OS.time)
  surv$OS <- as.numeric(surv$OS)
  
  keep <- !is.na(surv$OS.time) &
    !is.na(surv$OS) &
    surv$OS.time > 0
  
  model_expr <- model_expr[, keep, drop = FALSE]
  surv <- surv[keep, ]
  
  cat("Samples after survival filtering:\n")
  print(ncol(model_expr))
  
  cat("Event table:\n")
  print(table(surv$OS))
  
  if(ncol(model_expr) < 20){
    stop(paste0("Too few samples after survival filtering in ", cohort_name))
  }
  
  if(length(unique(surv$OS)) < 2){
    stop(paste0("Only one survival status present in ", cohort_name))
  }
  
  ##########################################################
  ## Calculate risk score
  ##########################################################
  
  risk_score <- as.numeric(
    t(model_expr) %*% coef_vec
  )
  
  risk_df <- data.frame(
    sample = colnames(model_expr),
    RiskScore = risk_score,
    OS.time = surv$OS.time,
    OS = surv$OS,
    stringsAsFactors = FALSE
  )
  
  ##########################################################
  ## Define cutoff
  ##########################################################
  
  if(cutoff_mode == "TCGA"){
    cutoff <- tcga_cutoff
  }else{
    cutoff <- median(risk_df$RiskScore, na.rm = TRUE)
  }
  
  risk_df$RiskGroup <- ifelse(
    risk_df$RiskScore >= cutoff,
    "High",
    "Low"
  )
  
  risk_df$RiskGroup <- factor(
    risk_df$RiskGroup,
    levels = c("Low", "High")
  )
  
  if(length(unique(risk_df$RiskGroup)) < 2){
    
    cat("Warning: one risk group is empty using ", cutoff_mode, " cutoff.\n", sep = "")
    cat("Using GEO median cutoff as fallback.\n")
    
    cutoff <- median(risk_df$RiskScore, na.rm = TRUE)
    
    risk_df$RiskGroup <- ifelse(
      risk_df$RiskScore >= cutoff,
      "High",
      "Low"
    )
    
    risk_df$RiskGroup <- factor(
      risk_df$RiskGroup,
      levels = c("Low", "High")
    )
    
    cutoff_mode <- paste0(cutoff_mode, "_fallbackMedian")
    
  }
  
  cat("Cutoff used:\n")
  print(cutoff)
  
  cat("Risk group table:\n")
  print(table(risk_df$RiskGroup))
  
  out_prefix <- paste0(
    "Module09_",
    cohort_name,
    "_",
    cutoff_mode
  )
  
  write.csv(
    risk_df,
    file.path(result_dir, paste0(out_prefix, "_RiskScore.csv")),
    row.names = FALSE
  )
  
  ##########################################################
  ## KM survival analysis
  ##########################################################
  
  fit_km <- survival::survfit(
    survival::Surv(OS.time, OS) ~ RiskGroup,
    data = risk_df
  )
  
  cox_group <- survival::coxph(
    survival::Surv(OS.time, OS) ~ RiskGroup,
    data = risk_df
  )
  
  cox_sum <- summary(cox_group)
  
  group_hr <- cox_sum$conf.int[1, "exp(coef)"]
  group_lower <- cox_sum$conf.int[1, "lower .95"]
  group_upper <- cox_sum$conf.int[1, "upper .95"]
  group_p <- cox_sum$coefficients[1, "Pr(>|z|)"]
  
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
    title = paste0(cohort_name, " survival by risk group"),
    risk.table.height = 0.25,
    ggtheme = theme_bw()
  )
  
  pdf(
    file.path(figure_dir, paste0(out_prefix, "_KM.pdf")),
    width = 7,
    height = 7
  )
  print(km_plot)
  dev.off()
  
  png(
    file.path(figure_dir, paste0(out_prefix, "_KM.png")),
    width = 2100,
    height = 2100,
    res = 300
  )
  print(km_plot)
  dev.off()
  
  ##########################################################
  ## Time-dependent ROC
  ##########################################################
  
  roc_times <- c(365, 1095, 1825)
  
  max_time <- max(risk_df$OS.time, na.rm = TRUE)
  
  valid_times <- roc_times[roc_times < max_time]
  
  if(length(valid_times) >= 1){
    
    roc_obj <- timeROC::timeROC(
      T = risk_df$OS.time,
      delta = risk_df$OS,
      marker = risk_df$RiskScore,
      cause = 1,
      times = valid_times,
      iid = TRUE
    )
    
    auc_df <- data.frame(
      Cohort = cohort_name,
      CutoffMode = cutoff_mode,
      Time = paste0(valid_times / 365, "-year"),
      Days = valid_times,
      AUC = as.numeric(roc_obj$AUC)
    )
    
    write.csv(
      auc_df,
      file.path(result_dir, paste0(out_prefix, "_TimeROC_AUC.csv")),
      row.names = FALSE
    )
    
    pdf(
      file.path(figure_dir, paste0(out_prefix, "_TimeROC.pdf")),
      width = 7,
      height = 6
    )
    
    plot(
      roc_obj,
      time = valid_times[1],
      col = "#D73027",
      title = FALSE
    )
    
    if(length(valid_times) >= 2){
      plot(
        roc_obj,
        time = valid_times[2],
        add = TRUE,
        col = "#4575B4"
      )
    }
    
    if(length(valid_times) >= 3){
      plot(
        roc_obj,
        time = valid_times[3],
        add = TRUE,
        col = "#1A9850"
      )
    }
    
    legend_cols <- c("#D73027", "#4575B4", "#1A9850")[seq_along(valid_times)]
    
    legend(
      "bottomright",
      legend = paste0(
        auc_df$Time,
        " AUC = ",
        sprintf("%.3f", auc_df$AUC)
      ),
      col = legend_cols,
      lwd = 2,
      bty = "n"
    )
    
    title(
      main = paste0(cohort_name, " time-dependent ROC")
    )
    
    dev.off()
    
  }else{
    
    roc_obj <- NULL
    
    auc_df <- data.frame(
      Cohort = cohort_name,
      CutoffMode = cutoff_mode,
      Time = NA,
      Days = NA,
      AUC = NA
    )
    
  }
  
  ##########################################################
  ## Risk score distribution
  ##########################################################
  
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
      title = paste0(cohort_name, " risk score distribution"),
      x = "Patients ranked by risk score",
      y = "Risk score"
    )
  
  ggsave(
    file.path(figure_dir, paste0(out_prefix, "_RiskScore_Distribution.pdf")),
    p_risk,
    width = 7,
    height = 5
  )
  
  ggsave(
    file.path(figure_dir, paste0(out_prefix, "_RiskScore_Distribution.png")),
    p_risk,
    width = 7,
    height = 5,
    dpi = 300
  )
  
  ##########################################################
  ## Survival status plot
  ##########################################################
  
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
      title = paste0(cohort_name, " survival status distribution"),
      x = "Patients ranked by risk score",
      y = "Survival time (days)"
    )
  
  ggsave(
    file.path(figure_dir, paste0(out_prefix, "_SurvivalStatus.pdf")),
    p_status,
    width = 7,
    height = 5
  )
  
  ggsave(
    file.path(figure_dir, paste0(out_prefix, "_SurvivalStatus.png")),
    p_status,
    width = 7,
    height = 5,
    dpi = 300
  )
  
  ##########################################################
  ## Heatmap
  ##########################################################
  
  heat_expr <- model_expr[
    ,
    risk_order$sample,
    drop = FALSE
  ]
  
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
    file.path(figure_dir, paste0(out_prefix, "_ModelGene_Heatmap.pdf")),
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
    main = paste0(cohort_name, " model gene expression")
  )
  
  dev.off()
  
  ##########################################################
  ## Performance summary
  ##########################################################
  
  performance <- data.frame(
    Cohort = cohort_name,
    CutoffMode = cutoff_mode,
    Samples = nrow(risk_df),
    Cutoff = cutoff,
    LowRisk = sum(risk_df$RiskGroup == "Low"),
    HighRisk = sum(risk_df$RiskGroup == "High"),
    Events = sum(risk_df$OS == 1),
    Censored = sum(risk_df$OS == 0),
    HR_high_vs_low = group_hr,
    HR_lower95CI = group_lower,
    HR_upper95CI = group_upper,
    Cox_Pvalue = group_p,
    AUC_1year = ifelse(any(auc_df$Days == 365), auc_df$AUC[auc_df$Days == 365], NA),
    AUC_3year = ifelse(any(auc_df$Days == 1095), auc_df$AUC[auc_df$Days == 1095], NA),
    AUC_5year = ifelse(any(auc_df$Days == 1825), auc_df$AUC[auc_df$Days == 1825], NA)
  )
  
  return(
    list(
      risk_df = risk_df,
      model_expr = model_expr,
      surv = surv,
      fit_km = fit_km,
      cox_group = cox_group,
      roc_obj = roc_obj,
      auc_df = auc_df,
      performance = performance
    )
  )
  
}

############################################################
## 7. Read model and TCGA cutoff
############################################################

cat("Reading final model...\n")

model <- data.table::fread(
  model_file,
  data.table = FALSE
)

required_model_cols <- c("Gene", "Coef")

missing_model_cols <- setdiff(
  required_model_cols,
  colnames(model)
)

if(length(missing_model_cols) > 0){
  stop(
    paste0(
      "Missing columns in model file: ",
      paste(missing_model_cols, collapse = ", ")
    )
  )
}

cat("Final model genes:\n")
print(model$Gene)

cat("Reading TCGA risk score...\n")

tcga_risk <- data.table::fread(
  tcga_risk_file,
  data.table = FALSE
)

tcga_cutoff <- median(
  tcga_risk$RiskScore,
  na.rm = TRUE
)

cat("TCGA median cutoff:\n")
print(tcga_cutoff)

############################################################
## 8. Read GEO data
############################################################

cat("Reading GSE37642 expression and phenotype...\n")

gse37642_expr <- read_expr_csv(
  gse37642_expr_file
)

gse37642_pheno <- data.table::fread(
  gse37642_pheno_file,
  data.table = FALSE,
  check.names = FALSE
)

cat("GSE37642 expression dimension:\n")
print(dim(gse37642_expr))

cat("GSE37642 phenotype dimension:\n")
print(dim(gse37642_pheno))

cat("Reading GSE12417 expression and phenotype...\n")

gse12417_expr <- read_expr_csv(
  gse12417_expr_file
)

gse12417_pheno <- data.table::fread(
  gse12417_pheno_file,
  data.table = FALSE,
  check.names = FALSE
)

cat("GSE12417 expression dimension:\n")
print(dim(gse12417_expr))

cat("GSE12417 phenotype dimension:\n")
print(dim(gse12417_pheno))

############################################################
## 9. Validate each GEO cohort
############################################################

gse37642_tcga <- validate_geo(
  cohort_name = "GSE37642",
  expr = gse37642_expr,
  pheno = gse37642_pheno,
  model = model,
  tcga_cutoff = tcga_cutoff,
  cutoff_mode = "TCGA"
)

gse37642_median <- validate_geo(
  cohort_name = "GSE37642",
  expr = gse37642_expr,
  pheno = gse37642_pheno,
  model = model,
  tcga_cutoff = tcga_cutoff,
  cutoff_mode = "Median"
)

gse12417_tcga <- validate_geo(
  cohort_name = "GSE12417",
  expr = gse12417_expr,
  pheno = gse12417_pheno,
  model = model,
  tcga_cutoff = tcga_cutoff,
  cutoff_mode = "TCGA"
)

gse12417_median <- validate_geo(
  cohort_name = "GSE12417",
  expr = gse12417_expr,
  pheno = gse12417_pheno,
  model = model,
  tcga_cutoff = tcga_cutoff,
  cutoff_mode = "Median"
)

############################################################
## 10. Save combined performance table
############################################################

performance_all <- dplyr::bind_rows(
  gse37642_tcga$performance,
  gse37642_median$performance,
  gse12417_tcga$performance,
  gse12417_median$performance
)

write.csv(
  performance_all,
  file.path(result_dir, "Module09_ExternalValidation_Performance.csv"),
  row.names = FALSE
)

cat("\nExternal validation performance:\n")
print(performance_all)

############################################################
## 11. Save RData object
############################################################

save(
  model,
  tcga_cutoff,
  gse37642_expr,
  gse12417_expr,
  gse37642_pheno,
  gse12417_pheno,
  gse37642_tcga,
  gse37642_median,
  gse12417_tcga,
  gse12417_median,
  performance_all,
  file = file.path(object_dir, "Module09_GEO_ExternalValidation.RData")
)

############################################################
## 12. Generate report
############################################################

report_file <- file.path(
  result_dir,
  "Module09_Report.txt"
)

sink(report_file)

cat("Module09 Report: GEO external validation\n")
cat("========================================\n\n")

cat("Project:\n")
cat("Lysine-associated gene signature predicts prognosis in AML\n\n")

cat("Project directory:\n")
cat(project_dir, "\n\n")

cat("Input files:\n")
cat(model_file, "\n")
cat(tcga_risk_file, "\n")
cat(gse37642_expr_file, "\n")
cat(gse12417_expr_file, "\n")
cat(gse37642_pheno_file, "\n")
cat(gse12417_pheno_file, "\n\n")

cat("Final model genes:\n")
print(model$Gene)

cat("\nTCGA median cutoff:\n")
print(tcga_cutoff)

cat("\nExternal validation performance:\n")
print(performance_all)

cat("\nGenerated main files:\n")
cat(file.path(result_dir, "Module09_ExternalValidation_Performance.csv"), "\n")
cat(file.path(object_dir, "Module09_GEO_ExternalValidation.RData"), "\n\n")

cat("Diagnostic files:\n")
cat(file.path(result_dir, "Module09_GSE12417_OS_ClinicalColumns.csv"), "\n")
cat(file.path(result_dir, "Module09_GSE12417_ExtractedSurvival_AllRows.csv"), "\n")
cat(file.path(result_dir, "Module09_GSE12417_ExtractedSurvival.csv"), "\n\n")

cat("Session information:\n")
print(sessionInfo())

sink()

cat("\nText report saved:\n")
cat(report_file, "\n")

############################################################
## 13. Finish
############################################################

cat("\n============================================\n")
cat("Module09 finished successfully.\n")
cat("============================================\n")

