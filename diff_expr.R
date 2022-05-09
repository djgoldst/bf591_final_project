## -- BF591 Final Project -- ##
## Differential Expression Tab
## Author: Daniel Goldstein
## Email: djgoldst@bu.edu

library(tidyverse)

options(scipen=0)

de_volcano_plot <- function(dataf, x_name, y_name, slider, color1, color2) {
  dataf %>%
    select(c(x_name,y_name)) %>%
    mutate(volc_status = ifelse(.[,y_name] < 10^slider, "TRUE", "FALSE")) %>%
    ggplot() +
    geom_point(aes(x = !!sym(x_name),
                   y = -log10(!!sym(y_name)),
                   color = factor(volc_status, levels = c("TRUE","FALSE")))) +
    scale_color_manual(values = c(color2,color1)) +
    theme_classic() +
    theme(legend.position = "bottom") +
    labs(x = x_name,
         y = paste0("-log10(",y_name,")"),
         color = paste0(y_name," < 10^",slider)) %>%
    return()
}

de_summary_table <- function(dataf, slider) {
  dataf %>%
    rename(Symbol = symbol) %>%
    arrange(pvalue) %>%
    filter(padj < 10^slider) %>%
    mutate(pvalue = formatC(.$pvalue, digits = 2, format = "e"),
           padj = formatC(.$padj, digits = 2, format = "e")) %>%
    return()
}

