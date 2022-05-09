## -- BF591 Final Project -- ##
## Sample Information Tab
## Author: Daniel Goldstein
## Email: djgoldst@bu.edu

library(tidyverse)
library(data.table)
library(Hmisc)

options(scipen=999, warn=-1)

label_sample_data <- function(data){
  var.labels <- c(Age_of_Onset="HD only",
                  Duration="HD only",
                  CAG="HD only",
                  Vonsattel_Grade="HD only",
                  H_V_Striatal_Score="HD only",
                  H_V_Cortical_Score="HD_only")
  label(data) <- as.list(var.labels[match(names(data), names(var.labels))])
  return(data)
}

subset_data <- function(data, cond = c("Control","HD")){
  if(cond == "Control"){
    filter(data, Condition == cond) %>%
      select(which(colSums(is.na(.)) != nrow(.))) %>%
      return()
  } else{
    filter(data, Condition == cond) %>%
      return()
  }
}

mean_or_distinct <- function(data){
  if(typeof(data) == "double" | typeof(data) == "integer"){
    data_mean <- format(mean(data,na.rm=TRUE),digits=4)
    data_sd <- format(sd(data, na.rm=TRUE),digits=4)
    paste0(data_mean," (+/- ",data_sd,")")
  } else if(typeof(data) == "character"){
    paste0(unique(data),collapse=",")
  }
}

sample_summary_table <- function(data, cols){
  tibble(`Column Name` = cols,
         Type = c("character",rep("numeric", length(cols)-1)),
         `Mean (+/- SD) or Distinct Values` = c(mean_or_distinct(data[,cols[1]]),apply(data[,cols[-1]], 2, mean_or_distinct)))
}

sample_summary_plot <- function(data, ctrl_data, hd_data, sum_var){
  if(sum_var %in% colnames(ctrl_data)){
    ggplot(data) +
      geom_violin(aes(x=Condition,
                      y=!!sym(sum_var),
                      fill=Condition)) +
      theme_classic() +
      scale_fill_brewer(type="qual",palette="Paired") +
      labs(y = gsub("_"," ",sum_var),
           title = "Metadata Summary Distribution") %>%
      return()
  } else{
    hd_data %>%
      ggplot() +
      geom_violin(aes(x=Condition,
                      y=!!sym(sum_var), 
                      fill = Condition)) +
      theme_classic() +
      scale_fill_brewer(type="qual",palette="Set2") +
      labs(y = gsub("_"," ",sum_var),
           title = "Metadata Summary Distribution (HD only)")
  }

}

