############################################################
##
## Module05B_Cox_Visualization.R
##
## Project:
## Lysine-associated gene signature predicts prognosis in AML
##
## Purpose:
## Additional visualization for Module05 univariate Cox results
##
## Input:
## result/Module05_UnivariateCox_AllGenes.csv
## result/Module05_UnivariateCox_SignificantGenes.csv
##
## Output:
## figure/Module05B_Cox_Volcano_Labelled.pdf
## figure/Module05B_Cox_Forest_Top20.pdf
## figure/Module05B_Cox_Lollipop_Top30.pdf
## figure/Module05B_Cox_RiskProtective_Barplot.pdf
##
############################################################

rm(list = ls())
options(stringsAsFactors = FALSE)

cat("============================================\n")
cat("Module05B: Cox visualization with gene labels\n")
cat("============================================\n\n")

############################################################
## 1. Packages
############################################################

cran_pkgs <- c(
  "data.table",
  "dplyr",
  "ggplot2",
  "ggrepel",
  "stringr"
)

for(pkg in cran_pkgs){
  if(!requireNamespace(pkg, quietly = TRUE)){
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

dir.create("figure", showWarnings = FALSE, recursive = TRUE)
dir.create("result", showWarnings = FALSE, recursive = TRUE)

############################################################
## 2. Read Module05 results
############################################################

cox_all_file <- "C:/Users/Xiong/Desktop/AML/result/Module05_UnivariateCox_AllGenes.csv"
cox_sig_file <- "C:/Users/Xiong/Desktop/AML/result/Module05_UnivariateCox_SignificantGenes.csv"

if(!file.exists(cox_all_file)){
  stop("Cannot find result/Module05_UnivariateCox_AllGenes.csv")
}

if(!file.exists(cox_sig_file)){
  stop("Cannot find result/Module05_UnivariateCox_SignificantGenes.csv")
}

cox_all <- data.table::fread(
  cox_all_file,
  data.table = FALSE
)

cox_sig <- data.table::fread(
  cox_sig_file,
  data.table = FALSE
)

cat("All Cox genes:\n")
print(nrow(cox_all))

cat("Significant Cox genes:\n")
print(nrow(cox_sig))

############################################################
## 3. Prepare variables
############################################################

cox_all <- cox_all %>%
  dplyr::mutate(
    log2HR = log2(HR),
    minusLog10P = -log10(Pvalue),
    CoxType = dplyr::case_when(
      Pvalue < 0.05 & HR > 1 ~ "Risk",
      Pvalue < 0.05 & HR < 1 ~ "Protective",
      TRUE ~ "Not significant"
    )
  )

cox_sig <- cox_sig %>%
  dplyr::mutate(
    log2HR = log2(HR),
    minusLog10P = -log10(Pvalue),
    CoxType = dplyr::case_when(
      HR > 1 ~ "Risk",
      HR < 1 ~ "Protective",
      TRUE ~ "Not significant"
    )
  )

############################################################
## 4. Select genes to label
############################################################

## 自动标注：
## 1. P值最小的前10个基因
## 2. 风险基因中P值最小的前5个
## 3. 保护基因中P值最小的前5个

top_overall <- cox_all %>%
  dplyr::arrange(Pvalue) %>%
  dplyr::slice_head(n = 10) %>%
  dplyr::pull(Gene)

top_risk <- cox_all %>%
  dplyr::filter(Pvalue < 0.05, HR > 1) %>%
  dplyr::arrange(Pvalue) %>%
  dplyr::slice_head(n = 5) %>%
  dplyr::pull(Gene)

top_protective <- cox_all %>%
  dplyr::filter(Pvalue < 0.05, HR < 1) %>%
  dplyr::arrange(Pvalue) %>%
  dplyr::slice_head(n = 5) %>%
  dplyr::pull(Gene)

## 你也可以手动添加你特别想标注的基因
manual_genes <- c(
  "PARP3",
  "KMT2E",
  "DESI1",
  "SMC5",
  "PARP10",
  "LOXL4",
  "SIRT6",
  "HDAC11",
  "SIRT2"
)

label_genes <- unique(
  c(
    top_overall,
    top_risk,
    top_protective,
    manual_genes
  )
)

label_df <- cox_all %>%
  dplyr::filter(Gene %in% label_genes)

cat("Genes labelled in plots:\n")
print(label_df$Gene)

write.csv(
  label_df,
  "C:/Users/Xiong/Desktop/AML/result/Module05B_Labelled_Genes.csv",
  row.names = FALSE
)

############################################################
## 5. Labelled Cox volcano plot
############################################################

p_volcano <- ggplot(
  cox_all,
  aes(
    x = log2HR,
    y = minusLog10P
  )
) +
  geom_point(
    aes(color = CoxType),
    alpha = 0.85,
    size = 2.2
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.4
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    linewidth = 0.4
  ) +
  ggrepel::geom_text_repel(
    data = label_df,
    aes(label = Gene),
    size = 3.5,
    max.overlaps = Inf,
    box.padding = 0.5,
    point.padding = 0.3,
    segment.size = 0.3
  ) +
  scale_color_manual(
    values = c(
      "Risk" = "#D73027",
      "Protective" = "#4575B4",
      "Not significant" = "grey70"
    )
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.title = element_blank(),
    panel.grid = element_blank()
  ) +
  labs(
    title = "Univariate Cox analysis of lysine-associated genes",
    x = "log2(Hazard ratio)",
    y = "-log10(P value)"
  )

ggsave(
  "C:/Users/Xiong/Desktop/AML/figure/Module05B_Cox_Volcano_Labelled.pdf",
  p_volcano,
  width = 8,
  height = 6
)

ggsave(
  "C:/Users/Xiong/Desktop/AML/figure/Module05B_Cox_Volcano_Labelled.png",
  p_volcano,
  width = 8,
  height = 6,
  dpi = 300
)

cat("Labelled Cox volcano plot saved.\n")

############################################################
## 6. Forest plot for top 20 significant genes
############################################################

forest_df <- cox_sig %>%
  dplyr::arrange(Pvalue) %>%
  dplyr::slice_head(n = 20)

forest_df <- forest_df %>%
  dplyr::mutate(
    Gene = factor(Gene, levels = rev(Gene)),
    HR_label = paste0(
      sprintf("%.2f", HR),
      " (",
      sprintf("%.2f", Lower95CI),
      "-",
      sprintf("%.2f", Upper95CI),
      ")"
    )
  )

p_forest <- ggplot(
  forest_df,
  aes(
    x = Gene,
    y = HR,
    ymin = Lower95CI,
    ymax = Upper95CI
  )
) +
  geom_pointrange(
    aes(color = RiskType),
    size = 0.6
  ) +
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
    title = "Top prognostic lysine-associated genes",
    x = "",
    y = "Hazard ratio (log scale)"
  )

ggsave(
  "C:/Users/Xiong/Desktop/AML/figure/Module05B_Cox_Forest_Top20.pdf",
  p_forest,
  width = 7,
  height = 6
)

ggsave(
  "C:/Users/Xiong/Desktop/AML/figure/Module05B_Cox_Forest_Top20.png",
  p_forest,
  width = 7,
  height = 6,
  dpi = 300
)

cat("Forest plot saved.\n")

############################################################
## 7. Lollipop plot for top 30 genes
############################################################

lollipop_df <- cox_all %>%
  dplyr::arrange(Pvalue) %>%
  dplyr::slice_head(n = 30) %>%
  dplyr::mutate(
    Gene = factor(Gene, levels = rev(Gene)),
    Direction = ifelse(HR > 1, "Risk", "Protective")
  )

p_lollipop <- ggplot(
  lollipop_df,
  aes(
    x = Gene,
    y = minusLog10P
  )
) +
  geom_segment(
    aes(
      x = Gene,
      xend = Gene,
      y = 0,
      yend = minusLog10P
    ),
    linewidth = 0.5
  ) +
  geom_point(
    aes(color = Direction, size = abs(log2HR)),
    alpha = 0.9
  ) +
  coord_flip() +
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
    title = "Top 30 lysine-associated prognostic genes",
    x = "",
    y = "-log10(P value)",
    size = "|log2(HR)|"
  )

ggsave(
  "C:/Users/Xiong/Desktop/AML/figure/Module05B_Cox_Lollipop_Top30.pdf",
  p_lollipop,
  width = 7,
  height = 8
)

ggsave(
  "C:/Users/Xiong/Desktop/AML/figure/Module05B_Cox_Lollipop_Top30.png",
  p_lollipop,
  width = 7,
  height = 8,
  dpi = 300
)

cat("Lollipop plot saved.\n")

############################################################
## 8. Risk / protective gene count bar plot
############################################################

bar_df <- cox_sig %>%
  dplyr::count(RiskType)

p_bar <- ggplot(
  bar_df,
  aes(
    x = RiskType,
    y = n,
    fill = RiskType
  )
) +
  geom_bar(
    stat = "identity",
    width = 0.65
  ) +
  geom_text(
    aes(label = n),
    vjust = -0.5,
    size = 5
  ) +
  scale_fill_manual(
    values = c(
      "Risk" = "#D73027",
      "Protective" = "#4575B4"
    )
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none",
    panel.grid = element_blank()
  ) +
  labs(
    title = "Risk and protective genes from univariate Cox analysis",
    x = "",
    y = "Number of genes"
  )

ggsave(
  "C:/Users/Xiong/Desktop/AML/figure/Module05B_Cox_RiskProtective_Barplot.pdf",
  p_bar,
  width = 5,
  height = 5
)

ggsave(
  "C:/Users/Xiong/Desktop/AML/figure/Module05B_Cox_RiskProtective_Barplot.png",
  p_bar,
  width = 5,
  height = 5,
  dpi = 300
)

cat("Risk/protective bar plot saved.\n")

############################################################
## 9. Save visualization object
############################################################

save(
  cox_all,
  cox_sig,
  label_df,
  forest_df,
  lollipop_df,
  bar_df,
  file = "object/Module05B_Cox_Visualization.RData"
)

############################################################
## 10. Finish
############################################################

cat("\n============================================\n")
cat("Module05B Cox visualization finished.\n")
cat("Generated figures:\n")
cat("figure/Module05B_Cox_Volcano_Labelled.pdf\n")
cat("figure/Module05B_Cox_Forest_Top20.pdf\n")
cat("figure/Module05B_Cox_Lollipop_Top30.pdf\n")
cat("figure/Module05B_Cox_RiskProtective_Barplot.pdf\n")
cat("============================================\n")