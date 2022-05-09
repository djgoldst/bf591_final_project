## -- BF591 Final Project -- ##
## Counts Matrix Exploration Tab
## Author: Daniel Goldstein
## Email: djgoldst@bu.edu

library(tidyverse)
library(RColorBrewer)
library(pheatmap)
library(grid)

options(scipen=0)

filter_counts <- function(counts, perc_var = 0, num_nonzero = 0){
  #variance percentile refers to distribution of gene variances 
  counts_var <- apply(counts,1,var)
  var_percentile <- as.vector(counts_var) %>%
    quantile(probs = perc_var/100)
  counts[counts_var >= var_percentile,] %>%
    filter(rowSums(. > 0) >= num_nonzero)
}

counts_summary_table <- function(counts, perc_var = 0, num_nonzero = 0){
  filt_counts <- filter_counts(counts, perc_var, num_nonzero)
  tibble(`Total samples` = ncol(counts),
         `Total genes` = nrow(counts),
         `Genes passing filter (% passed)` = paste0(nrow(filt_counts)," (",
                                                   format((nrow(filt_counts)/nrow(counts))*100,digits=4),"%)"),
         `Genes failing filter (% failed)` = paste0(nrow(counts)-nrow(filt_counts)," (",
                                                   format(((nrow(counts)-nrow(filt_counts))/nrow(counts))*100,digits=4),"%)"))
}

diagnostic_plot <- function(counts, plot_type = c("variance","non_zeros"), perc_var = 0, num_nonzero = 0, log_scale = FALSE){
  if(plot_type == "variance"){
    medians = apply(counts, 1, median)
    counts_var = apply(counts, 1, var)
    var_percentile <- as.vector(counts_var) %>%
      quantile(probs = perc_var/100)
    plot_data <- tibble(median = medians, 
                        var = counts_var) %>%
      mutate(pass_filt = factor(var >= var_percentile, levels = c(TRUE,FALSE)),
             rank = rank(median))
    
    median_var_plot <- ggplot(plot_data, 
                              aes(x = rank,
                                  y = var)) +
      geom_point(aes(color = pass_filt)) +
      geom_smooth(method='gam', formula = y ~ s(x, bs = "cs"), color = "red") +
      scale_color_manual(values = brewer.pal(n = 12, "Paired")[2:1])+
      theme_classic() +
      theme(text = element_text(size = 9)) +
      labs(x = "Rank(Median)",
           y = "Variance",
           color = paste0("Variance Percentile >= ",perc_var,"%"))
    if(log_scale){
      median_var_plot + 
        scale_y_log10()
    }
  } else if(plot_type == "non_zeros"){
    plot_data <- counts %>%
      rownames_to_column("ensgid") %>%
      mutate(nonzero_count = rowSums(.[-1] != 0)) %>%
      pivot_longer(-c(ensgid, nonzero_count),
                   names_to = "sample",
                   values_to = "count") %>%
      group_by(nonzero_count, ensgid) %>%
      dplyr::summarize(median_nonzero = median(count[count!=0]), 
                       .groups = "keep") %>%
      ungroup() %>%
      mutate(rank_median = rank(median_nonzero),
             pass_filt = factor(nonzero_count >= num_nonzero, levels = c(TRUE,FALSE)))
    
    ggplot(plot_data, 
           aes(x = rank_median,
               y = nonzero_count)) +
      geom_point(aes(color = pass_filt)) +
      geom_smooth(method='gam', formula = y ~ s(x, bs = "cs"), color = "red") +
      scale_y_continuous(breaks = seq(0,70,10)) +
      scale_color_manual(values = brewer.pal(n = 12, "Paired")[4:3]) +
      theme_classic() +
      theme(text = element_text(size = 9)) +
      labs(x = "Rank(Median)",
           y = "Non-Zero Sample Count",
           color = paste0("Non-Zero Sample Count >= ",num_nonzero))
  }
}

filt_counts_heatmap <- function(filt_counts){
  setHook("grid.newpage", function() pushViewport(viewport(x=0.95,y=0.95,width=0.85,height=0.85,name="vp",
                                                           just=c("right","top"))))
  pheatmap(log10(filt_counts+1),
          col = brewer.pal(n=9,"Blues"),
          legend_breaks = c(0,1,2,3,4,5,6,6.48),
          legend_labels = c("0","1","2","3","4","5","6","log10(counts)"),
          drop_levels = TRUE,
          show_rownames = FALSE,
          fontsize = 6,
          angle_col = 45)
  setHook("grid.newpage", NULL, "replace")
  grid.text("Genes", x=-0.03, rot=90, gp=gpar(fontsize=8))
  grid.text("Samples", y = -0.03, gp=gpar(fontsize=8))
}

pca_biplot <- function(filt_counts, metadata, first_PC = "PC1", second_PC = "PC2"){
  pca_res <- prcomp(scale(t(filt_counts)), center=FALSE, scale=FALSE)
  variance <- pca_res$sdev**2
  var_tbl <- tibble(PC = factor(str_c("PC",1:69),levels = str_c("PC",1:69)),
                    VE = variance/sum(variance))
  
  pca_condition <- data.frame(pca_res$x) %>%
    rownames_to_column("Sample_ID") %>%
    right_join(metadata[,c("Sample_ID","Condition")],
              by = "Sample_ID") %>%
    arrange(Condition)

  ggplot(pca_condition) +
    geom_point(aes(x = !!sym(first_PC),
                   y = !!sym(second_PC),
                   color = Condition)) +
    scale_color_brewer(type = "qual", palette = "Set2") +
    theme_classic() +
    theme(text = element_text(size = 9)) +
    labs(x = paste0(first_PC," (",round(var_tbl[var_tbl$PC == first_PC, "VE"],4)*100,"% Variance Explained)"),
         y = paste0(second_PC," (",round(var_tbl[var_tbl$PC == second_PC, "VE"],4)*100,"% Variance Explained)"))
}

