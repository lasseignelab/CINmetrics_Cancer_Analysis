### Functions for CINmetrics Cancer Analysis - Clinical Variable: Sex
# Orignal Code: Sasha Thalluri
# Code Adaptations: Tabea M. Soelter

## get_cinmetrics
# A function which calculates 6 CIN metrics for all cancer projects of interest. Sasha's original code did not include information on gender/sex, so this is adapted to include that information.
get_cinmetrics <- function(all_cancer_df) {
  # Calculating CINmetrics
  modified.tai.cancer <- taiModified(all_cancer_df)
  cinmetrics.cancer <- CINmetrics(all_cancer_df)
  cinmetrics.cancer <- inner_join(cinmetrics.cancer,
                                  modified.tai.cancer,
                                  by = "sample_id")
  # Adding Patient ID, Project, and Sample Type Columns to Dataframe
  cinmetrics.cancer <- cinmetrics.cancer %>% separate(sample_id,
                                                      "patient_id",
                                                      12,
                                                      remove = FALSE)
  all.cinmetrics.df <- inner_join(cinmetrics.cancer,
                                  all_cancer_df[,c(7,10,11,18)],
                                  by = c("sample_id" = "Sample")) %>%
    unique(.) %>%
    select("sample_id",
           "patient_id",
           "project",
           "gender",
           "sample_type",
           "tai",
           "modified_tai",
           "cna",
           "base_segments",
           "break_points",
           "fga")
  # We store the output as a csv
  write.csv(all.cinmetrics.df,
            file = here("data", "processed", "all_cinmetrics_sex_df.csv"),
            row.names = FALSE)
  return(all.cinmetrics.df)
}

## get_utest_data
# A function which takes a df with all cin metrics calculations and performs a Mann-Whitney U-test on each cancer accounting for sex. Orignal code by Sasha was adapted by Tabea to only include sex and not other clinical variables.
get_utest_data <- function(all_cinmetrics_sex_df, cancer_project, metrics) {
  # We create an empty dataframe to store the U-test data
  sex_utest_df <- data.frame(.y. = character(),
                             group1 = character(),
                             group2 = character(),
                             n1 = numeric(),
                             n2 = numeric(),
                             p = numeric(),
                             p.signif = character(),
                             p.adj = numeric(),
                             p.adj.signif = character())

  for(i in 1:length(cancer_project)) {
    sex.df <- subset(all_cinmetrics_sex_df, project == cancer_project[i])

    # We filter the dataframe to only have the sexes and tumor samples we want included in the analysis
    sex.df <- sex.df %>%
      filter(sex %in% c("male", "female") &
               sample_type %in% c("Metastatic",
                                  "Primary Blood Derived Cancer",
                                  "Primary Tumor",
                                  "Recurrent Tumor",
                                  "Additional - New Primary",
                                  "Primary Blood Derived Cancer - Peripheral Blood"))

    for(j in metrics){
      # We performed the Mann-Whitney U-test (otherwise known as the Wilcoxon test)
      cinmetrics <- sex.df[, j]
      sex <- sex.df[,"sex"]
      cinmetrics.sex.df <- cbind(cinmetrics, sex) %>% as.data.frame(.)
      colnames(cinmetrics.sex.df)[1] <- "cinmetrics"
      cinmetrics.sex.df$cinmetrics <- cinmetrics.sex.df$cinmetrics %>% as.numeric()
      utest.results <- cinmetrics.sex.df %>% wilcox_test(cinmetrics ~ sex, paired = FALSE, p.adjust.method = "bonferroni")
      utest.results$project <- cancer_project[i]
      utest.results$cin <- j
      sex_utest_df <- rbind(sex_utest_df, utest.results)
    }
  }
  # We will store this output in a CSV file
  sex_utest_df <- sex_utest_df[ , -1]
  write.csv(sex_utest_df, file = here("data", "processed", "sex_utest_df.csv"), row.names = FALSE)
  return(sex_utest_df)
}

## get_utest_graphs
# A function which plots a barplot for all cancers for each CIN metric's adjusted p-value. Orignal code by Sasha was adapted by Tabea for sex bias plotting purposes.
get_utest_graphs <- function(cancer_project){
  # Filter columns and -Log10 transformation of p values
  filtered.sex.df<- sex_utest_df %>% select(cin, project, group1, group2)
  # Creating a function to input into sapply
  neglog10 <- function(x){
    neglog10 <- -log10(x)
    return(neglog10)
  }
  log.sex.utest <- sex_utest_df %>% select(p) %>% sapply(., neglog10)
  log.sex.utest <- cbind(filtered.sex.df, log.sex.utest)
  log.sex.utest$project <- str_remove(log.sex.utest$project, "TCGA-")

  # Base Segments Plot
  base.segments.utest <- filter(log.sex.utest, cin == "base_segments")
  base.segments.barplot <- ggbarplot(base.segments.utest,
                                     x = "project", y = "p",
                                     rotate = TRUE,
                                     xlab = "Cancer",
                                     ylab = "-Log10 p-value (female v. male)",
                                     width = 0.4,
                                     sort.val = "asc",
                                     fill = "black",
                                     title = "Base Segments",
                                     font.main = c(11, "bold"),
                                     font.x = c(10, "bold"),
                                     font.y = c(10, "bold"),
                                     font.xtickslab = c(10, "bold"),
                                     font.ytickslab = c(6, "bold"),) + geom_hline(yintercept = -log10(0.05/(length(base.segments.utest$p))), linetype = "dashed", color = "red")
  # Breakpoints Plot
  break.points.utest <- filter(log.sex.utest, cin == "break_points")
  break.points.barplot <- ggbarplot(break.points.utest,
                                    x = "project", y = "p",
                                    rotate = TRUE,
                                    xlab = "Cancer",
                                    ylab = "-Log10 p-value (female v. male)",
                                    width = 0.4,
                                    sort.val = "asc",
                                    fill = "black",
                                    title = "Breakpoints",
                                    font.main = c(11, "bold"),
                                    font.x = c(10, "bold"),
                                    font.y = c(10, "bold"),
                                    font.xtickslab = c(10, "bold"),
                                    font.ytickslab = c(6, "bold"),) + geom_hline(yintercept = -log10(0.05/(length(break.points.utest$p))), linetype = "dashed", color = "red")
  # TAI bar plot
  tai.utest <- filter(log.sex.utest, cin == "tai")
  tai.barplot <- ggbarplot(tai.utest,
                           x = "project", y = "p",
                           rotate = TRUE,
                           xlab = "Cancer",
                           ylab = "-Log10 p-value (female v. male)",
                           width = 0.4,
                           sort.val = "asc",
                           fill = "black",
                           title = "TAI",
                           font.main = c(11, "bold"),
                           font.x = c(10, "bold"),
                           font.y = c(10, "bold"),
                           font.xtickslab = c(10, "bold"),
                           font.ytickslab = c(6, "bold"),) + geom_hline(yintercept = -log10(0.05/(length(tai.utest$p))), linetype = "dashed", color = "red")
  # FGA bar plot
  fga.utest <- filter(log.sex.utest, cin == "fga")
  fga.barplot <- ggbarplot(fga.utest,
                           x = "project", y = "p",
                           rotate = TRUE,
                           xlab = "Cancer",
                           ylab = "-Log10 p-value (female v. male)",
                           width = 0.4,
                           sort.val = "asc",
                           fill = "black",
                           title = "FGA",
                           font.main = c(11, "bold"),
                           font.x = c(10, "bold"),
                           font.y = c(10, "bold"),
                           font.xtickslab = c(10, "bold"),
                           font.ytickslab = c(6, "bold"),) + geom_hline(yintercept = -log10(0.05/(length(fga.utest$p))), linetype = "dashed", color = "red")
  # CNA bar plot
  cna.utest <- filter(log.sex.utest, cin == "cna")
  cna.barplot <- ggbarplot(cna.utest,
                           x = "project", y = "p",
                           rotate = TRUE,
                           xlab = "Cancer",
                           ylab = "-Log10 p-value (female v. male)",
                           width = 0.4,
                           sort.val = "asc",
                           fill = "black",
                           title = "CNA",
                           font.main = c(11, "bold"),
                           font.x = c(10, "bold"),
                           font.y = c(10, "bold"),
                           font.xtickslab = c(10, "bold"),
                           font.ytickslab = c(6, "bold"),) + geom_hline(yintercept = -log10(0.05/(length(cna.utest$p))), linetype = "dashed", color = "red")

  # Modified TAI bar plot
  modified.tai.utest <- filter(log.sex.utest, cin == "modified_tai")
  modified.tai.barplot <- ggbarplot(modified.tai.utest,
                                    x = "project", y = "p",
                                    rotate = TRUE,
                                    xlab = "Cancer",
                                    ylab = "-Log10 p-value (female v. male)",
                                    width = 0.4,
                                    sort.val = "asc",
                                    fill = "black",
                                    title = "Modified TAI",
                                    font.main = c(11, "bold"),
                                    font.x = c(10, "bold"),
                                    font.y = c(10, "bold"),
                                    font.xtickslab = c(10, "bold"),
                                    font.ytickslab = c(6, "bold"),) + geom_hline(yintercept = -log10(0.05/(length(modified.tai.utest$p))), linetype = "dashed", color = "red")

  # Tumor Stage Bar Plot
  sex.barplot <- ggarrange(base.segments.barplot + rremove("xlab"),
                           fga.barplot + rremove("xlab"),
                           break.points.barplot + rremove("xlab"),
                           cna.barplot + rremove("xlab"),
                           tai.barplot + rremove("xlab"), modified.tai.barplot + rremove("xlab"),
                           ncol = 2, nrow = 3)
  sex.barplot <- annotate_figure(sex.barplot,
                                 bottom = text_grob("-Log10 p-value (female v. male)"))
  return(sex.barplot)

}

## get_utest_combined_graphs
# A function which plots a barplot of the top 10 most significant cancers with sex bias along with a raincloud plot for the top most cancer to show its sex distribution. Code adapted by Tabea from Sasha.
get_utest_combined_graphs <- function(cancer_project, sex_utest_df, all_cinmetrics_sex_df){
  # Filter columns and -Log10 transformation of p values
  filtered.sex.df <- sex_utest_df %>% select(cin, project, group1, group2)
  # Creating a function to input into sapply
  neglog10 <- function(x){
    neglog10 <- -log10(x)
    return(neglog10)
  }
  log.sex.utest <- sex_utest_df %>% select(p) %>% sapply(., neglog10)
  log.sex.utest <- cbind(filtered.sex.df, log.sex.utest)
  log.sex.utest$project <- str_remove(log.sex.utest$project, "TCGA-")

  # Base Segments Plot
  base.segments.utest <- filter(log.sex.utest, cin == "base_segments") %>% top_n(10)
  base.segments.barplot <- ggbarplot(base.segments.utest,
                                     x = "project", y = "p",
                                     rotate = TRUE,
                                     xlab = "Cancer",
                                     ylab = "-Log10 p-value (female v. male)",
                                     width = 0.4,
                                     sort.val = "asc",
                                     fill = "black",
                                     title = "Base Segments",
                                     font.main = c(11, "bold"),
                                     font.x = c(10, "bold"),
                                     font.y = c(10, "bold"),
                                     font.xtickslab = c(10, "bold"),
                                     font.ytickslab = c(6, "bold"),) +
    geom_hline(yintercept = -log10(0.05/(length(base.segments.utest$p))), linetype = "dashed", color = "red")

  # log10 transformation of Base Segments Values
  all_cinmetrics_sex_df$log10base_segments <- log10(all_cinmetrics_sex_df$base_segments)
  # subset for only TCGA-HNSC samples
  hnsc <- filter(all_cinmetrics_sex_df, project == "TCGA-HNSC" & sample_type == c("Metastatic", "Primary Tumor"))
  # Filter for HNSC and base segments
  hnsc.sex_utest_df <- filter(sex_utest_df, project == "TCGA-HNSC" & cin == "base_segments")
  # filter out any values below 7 for plotting purposes from hsnc df
  hnsc <- filter(hnsc, log10base_segments > "7")
  # raincloud plot
  hnsc_base_segments_raincloud <- ggplot(hnsc, aes(x = sex, y = log10base_segments, color = sex)) +
    ggtitle("Sex Distribution in TCGA-HNSC") +
    ylab("Base Segments") +
    xlab("TCGA-HNSC") +
    theme_cowplot() +
    scale_shape_identity() +
    theme(legend.position = "right",
          plot.title = element_text(size = 11),
          axis.title = element_text(size = 11, face = "bold"),
          axis.text = element_text(size = 11, face = "bold"),
          axis.text.x = element_text(angle = 0,
                                     hjust = 0.5,
                                     vjust = 0,
                                     face = "bold")) +
    scale_fill_manual(values = c("#FFA319FF", "#800000FF")) +
    scale_color_manual(values = c("#FFA319FF", "#800000FF")) +
    gghalves::geom_half_point(side = "l", range_scale = .4, alpha = .5, pch = c(21,25), fill = c("#FFA319FF", "#800000FF"), size = 2, transformation = position_jitter(width = 0.1), show.legend = FALSE) +
    guides(colour = guide_legend(override.aes = list(shape = c(21, 25)))) +
    geom_flat_violin(aes(sex, log10base_segments, fill = sex), # add the x and y parameters again
                     adjust = 2,
                     alpha = 0.6,
                     trim = TRUE,
                     scale = "width",
                     show.legend = FALSE) +
    geom_boxplot(fill = c("#FFA319FF", "#800000FF"),
                 notch = FALSE,
                 width = 0.06,
                 varwidth = FALSE,
                 outlier.shape = NA,
                 alpha = 0.3,
                 colour = "black",
                 show.legend = FALSE)

  # Breakpoints Plot
  break.points.utest <- filter(log.sex.utest, cin == "break_points") %>% top_n(10)
  break.points.barplot <- ggbarplot(break.points.utest,
                                    x = "project", y = "p",
                                    rotate = TRUE,
                                    xlab = "Cancer",
                                    ylab = "-Log10 p-value (female v. male)",
                                    width = 0.4,
                                    sort.val = "asc",
                                    fill = "black",
                                    title = "Breakpoints",
                                    font.main = c(11, "bold"),
                                    font.x = c(10, "bold"),
                                    font.y = c(10, "bold"),
                                    font.xtickslab = c(10, "bold"),
                                    font.ytickslab = c(6, "bold"),) + geom_hline(yintercept = -log10(0.05/(length(break.points.utest$p))), linetype = "dashed", color = "red")
  # Filter for HNSC and break points
  hnsc.sex_utest_df <- filter(sex_utest_df, project == "TCGA-HNSC" & cin == "break_points")
  # Raincloud plot
  hnsc_break_points_raincloud <- ggplot(hnsc, aes(x = sex, y = break_points, color = sex)) +
    ggtitle("Sex Distribution in TCGA-HNSC") +
    ylab("Break Points") +
    xlab("TCGA-HNSC") +
    theme_cowplot() +
    scale_shape_identity() +
    theme(legend.position = "right",
          plot.title = element_text(size = 11),
          axis.title = element_text(size = 11, face = "bold"),
          axis.text = element_text(size = 11, face = "bold"),
          axis.text.x = element_text(angle = 0,
                                     hjust = 0.5,
                                     vjust = 0,
                                     face = "bold")) +
    scale_fill_manual(values = c("#FFA319FF", "#800000FF")) +
    scale_color_manual(values = c("#FFA319FF", "#800000FF")) +
    gghalves::geom_half_point(side = "l", range_scale = .4, alpha = .5, pch = c(21,25), fill = c("#FFA319FF", "#800000FF"), size = 2, transformation = position_jitter(width = 0.1), show.legend = FALSE) +
    guides(colour = guide_legend(override.aes = list(shape = c(21, 25)))) +
    geom_flat_violin(aes(sex, break_points, fill = sex),
                     adjust = 2,
                     alpha = 0.6,
                     trim = TRUE,
                     scale = "width",
                     show.legend = FALSE) +
    geom_boxplot(fill = c("#FFA319FF", "#800000FF"),
                 notch = FALSE,
                 width = 0.06,
                 varwidth = FALSE,
                 outlier.shape = NA,
                 alpha = 0.3,
                 colour = "black",
                 show.legend = FALSE)
  # TAI bar plot
  tai.utest <- filter(log.sex.utest, cin == "tai") %>% top_n(10)
  tai.barplot <- ggbarplot(tai.utest,
                           x = "project", y = "p",
                           rotate = TRUE,
                           xlab = "Cancer",
                           ylab = "-Log10 p-value (female v. male)",
                           width = 0.4,
                           sort.val = "asc",
                           fill = "black",
                           title = "TAI",
                           font.main = c(11, "bold"),
                           font.x = c(10, "bold"),
                           font.y = c(10, "bold"),
                           font.xtickslab = c(10, "bold"),
                           font.ytickslab = c(6, "bold"),) + geom_hline(yintercept = -log10(0.05/(length(tai.utest$p))), linetype = "dashed", color = "red")
  # Filter for COAD and tai
  coad <- filter(all_cinmetrics_sex_df, project == "TCGA-COAD" & sample_type == c("Metastatic", "Primary Tumor", "Recurrent Tumor"))
  # filter out any values above 1 for plotting purposes from coad df
  coad <- filter(coad, tai < "1")

  coad.sex_utest_df <- filter(sex_utest_df, project == "TCGA-COAD" & cin == "tai")

  coad_tai_raincloud <- ggplot(coad, aes(x = sex, y = tai, color = sex)) +
    ggtitle("Sex Distribution in TCGA-COAD") +
    ylab("TAI") +
    xlab("TCGA-COAD") +
    theme_cowplot() +
    scale_shape_identity() +
    theme(legend.position = "right",
          plot.title = element_text(size = 11),
          axis.title = element_text(size = 11, face = "bold"),
          axis.text = element_text(size = 11, face = "bold"),
          axis.text.x = element_text(angle = 0,
                                     hjust = 0.5,
                                     vjust = 0,
                                     face = "bold")) +
    scale_fill_manual(values = c("#FFA319FF", "#800000FF")) +
    scale_color_manual(values = c("#FFA319FF", "#800000FF")) +
    gghalves::geom_half_point(side = "l", range_scale = .4, alpha = .5, pch = c(21,25), fill = c("#FFA319FF", "#800000FF"), size = 2, transformation = position_jitter(width = 0.1), show.legend = FALSE) +
    guides(colour = guide_legend(override.aes = list(shape = c(21, 25)))) +
    geom_flat_violin(aes(sex, tai, fill = sex),
                     adjust = 2,
                     alpha = 0.6,
                     trim = TRUE,
                     scale = "width",
                     show.legend = FALSE) +
    geom_boxplot(fill = c("#FFA319FF", "#800000FF"),
                 notch = FALSE,
                 width = 0.06,
                 varwidth = FALSE,
                 outlier.shape = NA,
                 alpha = 0.3,
                 colour = "black",
                 show.legend = FALSE)

  # FGA bar plot
  fga.utest <- filter(log.sex.utest, cin == "fga") %>% top_n(10)
  fga.barplot <- ggbarplot(fga.utest,
                           x = "project", y = "p",
                           rotate = TRUE,
                           xlab = "Cancer",
                           ylab = "-Log10 p-value (female v. male)",
                           width = 0.4,
                           sort.val = "asc",
                           fill = "black",
                           title = "FGA",
                           font.main = c(11, "bold"),
                           font.x = c(10, "bold"),
                           font.y = c(10, "bold"),
                           font.xtickslab = c(10, "bold"),
                           font.ytickslab = c(6, "bold"),) + geom_hline(yintercept = -log10(0.05/(length(fga.utest$p))), linetype = "dashed", color = "red")

  # Filter for HNSC and FGA
  hnsc.sex_utest_df <- filter(sex_utest_df, project == "TCGA-HNSC" & cin == "fga")

  hnsc_fga_raincloud <- ggplot(hnsc, aes(x = sex, y = fga, color = sex)) +
    ggtitle("Sex Distribution in TCGA-HNSC") +
    ylab("FGA") +
    xlab("TCGA-HNSC") +
    theme_cowplot() +
    scale_shape_identity() +
    theme(legend.position = "right",
          plot.title = element_text(size = 11),
          axis.title = element_text(size = 11, face = "bold"),
          axis.text = element_text(size = 11, face = "bold"),
          axis.text.x = element_text(angle = 0,
                                     hjust = 0.5,
                                     vjust = 0,
                                     face = "bold")) +
    scale_fill_manual(values = c("#FFA319FF", "#800000FF")) +
    scale_color_manual(values = c("#FFA319FF", "#800000FF")) +
    gghalves::geom_half_point(side = "l", range_scale = .4, alpha = .5, pch = c(21,25), fill = c("#FFA319FF", "#800000FF"), size = 2, transformation = position_jitter(width = 0.1), show.legend = FALSE) +
    guides(colour = guide_legend(override.aes = list(shape = c(21, 25)))) +
    geom_flat_violin(aes(sex, fga, fill = sex),
                     adjust = 2,
                     alpha = 0.6,
                     trim = TRUE,
                     scale = "width",
                     show.legend = FALSE) +
    geom_boxplot(fill = c("#FFA319FF", "#800000FF"),
                 notch = FALSE,
                 width = 0.06,
                 varwidth = FALSE,
                 outlier.shape = NA,
                 alpha = 0.3,
                 colour = "black",
                 show.legend = FALSE)

  # CNA bar plot
  cna.utest <- filter(log.sex.utest, cin == "cna") %>% top_n(10)
  cna.barplot <- ggbarplot(cna.utest,
                           x = "project", y = "p",
                           rotate = TRUE,
                           xlab = "Cancer",
                           ylab = "-Log10 p-value (female v. male)",
                           width = 0.4,
                           sort.val = "asc",
                           fill = "black",
                           title = "CNA",
                           font.main = c(11, "bold"),
                           font.x = c(10, "bold"),
                           font.y = c(10, "bold"),
                           font.xtickslab = c(10, "bold"),
                           font.ytickslab = c(6, "bold"),) + geom_hline(yintercept = -log10(0.05/(length(cna.utest$p))), linetype = "dashed", color = "red")
  # Filter for HNSC and CNA
  hnsc.sex_utest_df <- filter(sex_utest_df, project == "TCGA-HNSC" & cin == "cna")
  # Raincloud plot
  hnsc_cna_raincloud <- ggplot(hnsc, aes(x = sex, y = cna, color = sex)) +
    ggtitle("Sex Distribution in TCGA-HNSC") +
    ylab("CNA") +
    xlab("TCGA-HNSC") +
    theme_cowplot() +
    scale_shape_identity() +
    theme(legend.position = "right",
          plot.title = element_text(size = 11),
          axis.title = element_text(size = 11, face = "bold"),
          axis.text = element_text(size = 11, face = "bold"),
          axis.text.x = element_text(angle = 0,
                                     hjust = 0.5,
                                     vjust = 0,
                                     face = "bold")) +
    scale_fill_manual(values = c("#FFA319FF", "#800000FF")) +
    scale_color_manual(values = c("#FFA319FF", "#800000FF")) +
    gghalves::geom_half_point(side = "l", range_scale = .4, alpha = .5, pch = c(21,25), fill = c("#FFA319FF", "#800000FF"), size = 2, transformation = position_jitter(width = 0.1), show.legend = FALSE) +
    guides(colour = guide_legend(override.aes = list(shape = c(21, 25)))) +
    geom_flat_violin(aes(sex, cna, fill = sex),
                     adjust = 2,
                     alpha = 0.6,
                     trim = TRUE,
                     scale = "width",
                     show.legend = FALSE) +
    geom_boxplot(fill = c("#FFA319FF", "#800000FF"),
                 notch = FALSE,
                 width = 0.06,
                 varwidth = FALSE,
                 outlier.shape = NA,
                 alpha = 0.3,
                 colour = "black",
                 show.legend = FALSE)

  # Modified TAI bar plot
  modified.tai.utest <- filter(log.sex.utest, cin == "modified_tai") %>% top_n(10)
  modified.tai.barplot <- ggbarplot(modified.tai.utest,
                                    x = "project", y = "p",
                                    rotate = TRUE,
                                    xlab = "Cancer",
                                    ylab = "-Log10 p-value (female v. male)",
                                    width = 0.4,
                                    sort.val = "asc",
                                    fill = "black",
                                    title = "Modified TAI",
                                    font.main = c(11, "bold"),
                                    font.x = c(10, "bold"),
                                    font.y = c(10, "bold"),
                                    font.xtickslab = c(10, "bold"),
                                    font.ytickslab = c(6, "bold"),) + geom_hline(yintercept = -log10(0.05/(length(modified.tai.utest$p))), linetype = "dashed", color = "red")

  # Filter for GBM and modified tai
  gbm <- filter(all_cinmetrics_sex_df, project == "TCGA-GBM" & sample_type == c("Primary Tumor", "Recurrent Tumor"))

  # filter out any values below 0 for plotting purposes from gbm df
  gbm <- filter(gbm, modified_tai > "0")

  gbm.sex_utest_df <- filter(sex_utest_df, project == "TCGA-GBM" & cin == "modified_tai")

  gbm_modified_tai_raincloud <- ggplot(gbm, aes(x = sex, y = modified_tai, color = sex)) +
    ggtitle("Sex Distribution in TCGA-GBM") +
    ylab("Modified TAI") +
    xlab("TCGA-GBM") +
    theme_cowplot() +
    scale_shape_identity() +
    theme(legend.position = "right",
          plot.title = element_text(size = 11),
          axis.title = element_text(size = 11, face = "bold"),
          axis.text = element_text(size = 11, face = "bold"),
          axis.text.x = element_text(angle = 0,
                                     hjust = 0.5,
                                     vjust = 0,
                                     face = "bold")) +
    scale_fill_manual(values = c("#FFA319FF", "#800000FF")) +
    scale_color_manual(values = c("#FFA319FF", "#800000FF")) +
    gghalves::geom_half_point(side = "l", range_scale = .4, alpha = .5, pch = c(21,25), fill = c("#FFA319FF", "#800000FF"), size = 2, transformation = position_jitter(width = 0.1), show.legend = FALSE) +
    guides(colour = guide_legend(override.aes = list(shape = c(21, 25)))) +
    geom_flat_violin(aes(sex, modified_tai, fill = sex),
                     adjust = 2,
                     alpha = 0.6,
                     trim = TRUE,
                     scale = "width",
                     show.legend = FALSE) +
    geom_boxplot(fill = c("#FFA319FF", "#800000FF"),
                 notch = FALSE,
                 width = 0.06,
                 varwidth = FALSE,
                 outlier.shape = NA,
                 alpha = 0.3,
                 colour = "black",
                 show.legend = FALSE)

  # Tumor Stage Bar Plot
  sex.barplot <- ggarrange(base.segments.barplot,
                           hnsc_base_segments_raincloud + rremove("xlab"),
                           fga.barplot,
                           hnsc_fga_raincloud + rremove("xlab"),
                           break.points.barplot,
                           hnsc_break_points_raincloud + rremove("xlab"),
                           cna.barplot,
                           hnsc_cna_raincloud + rremove("xlab"),
                           tai.barplot,
                           coad_tai_raincloud + rremove("xlab"),
                           modified.tai.barplot,
                           gbm_modified_tai_raincloud + rremove("xlab"),
                           ncol = 2, nrow = 6)
  #sex.barplot <- annotate_figure(sex.barplot,
  #bottom = text_grob("-Log10 p-value (female v. male)"))
  return(sex.barplot)

}
