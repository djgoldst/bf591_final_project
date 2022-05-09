## -- BF591 Final Project -- ##
## GSEA Tab
## Author: Daniel Goldstein
## Email: djgoldst@bu.edu

library(tidyverse)
# library(fgsea)

options(scipen=0)

# run_fgsea <- function(de_res, gmt, min_size = 15, max_size = 500) {
#   gene_ids <- de_res %>%
#     arrange(padj) %>%
#     pull(symbol)
# 
#   genelist <- de_res %>%
#     drop_na(symbol, log2FoldChange) %>%
#     distinct(symbol, log2FoldChange, .keep_all=TRUE) %>%
#     arrange(desc(log2FoldChange)) %>%
#     dplyr::select(symbol, log2FoldChange) %>%
#     deframe()
# 
#   gmt_path <- gmtPathways(gmt)
# 
#   fgsea_res <- fgsea(gmt_path, genelist, minSize=min_size, maxSize=max_size) %>% as_tibble()
# 
#   return(fgsea_res)
# }

fgsea_NES_barplot <- function(fgsea_results, pathway_name, num_paths){
  fgsea_NES <- arrange(fgsea_results, desc(NES))
  topgene_pathways <- head(fgsea_NES, num_paths) %>% .$pathway
  botgene_pathways <- tail(fgsea_NES, num_paths) %>% .$pathway
  
  fgsea_subset <- fgsea_NES %>%
    filter(pathway %in% c(topgene_pathways, botgene_pathways)) %>%
    mutate(pathway = factor(pathway),
           path_labs = str_replace_all(pathway,"_"," "),
           path_labs = fct_reorder(factor(path_labs), NES))
  
  ggplot(fgsea_subset, aes(x = path_labs,
                        y = NES,
                        fill = NES < 0)) +
    geom_bar(stat = "identity",
             show.legend = FALSE) +
    scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "darkblue")) +
    labs(x = "",
         y = "Normalized Enrichment Score (NES)",
         title = paste0("FGSEA Results for ",pathway_name," Genesets from MSigDB")) +
    theme_classic() +
    theme(text = element_text(size = 6)) +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 70)) +
    coord_flip() %>%
    return()
}

fgsea_summary_table <- function(fgsea_results, pathway_direction = c("positive","negative","all"), padj_filter){
  fgsea_summary <- fgsea_results %>%
    rename(Pathway = pathway) %>%
    arrange(padj) %>%
    filter(padj < 10^padj_filter) %>%
    mutate(pval = formatC(.$pval, digits = 2, format = "e"),
           padj = formatC(.$padj, digits = 2, format = "e")) %>%
    select(-leadingEdge)
  
  if(pathway_direction == "positive"){
    fgsea_summary %>%
      filter(NES > 0) %>% 
      return()
  } else if(pathway_direction == "negative"){
    fgsea_summary %>%
      filter(NES < 0) %>%
      return()
  } else if(pathway_direction == "all"){
    fgsea_summary %>%
      return()
  }
}

fgsea_NES_scatterplot <- function(fgsea_results, padj_filter){
  fgsea_res_padj <- fgsea_results %>%
    mutate(pass_filt = padj < 10^padj_filter)
  
  ggplot(fgsea_res_padj) +
    geom_point(aes(x = NES,
                   y = -log10(padj),
                   color = factor(pass_filt, levels = c(TRUE,FALSE)))) +
    scale_color_manual(values = c("#1827C3","#CCCCCC")) +
    theme_classic() +
    theme(text = element_text(size=9)) +
    labs(x = "Normalized Enrichment Score (NES)",
         y = "-log10(Adjusted P-Value)",
         color = paste0("Adjusted P-Value < 10^",padj_filter))
}

