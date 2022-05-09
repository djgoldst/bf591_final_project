# bf591_final_project

Dataset: mRNA-Seq Expression profiling of human post-mortem BA9 brain tissue for Huntington’s Disease and neurologically normal individuals

Project Details: RShiny app with four major components
1. Sample Information Exploration (samples.R)
  
    Input - Sample Information Table (CSV)
    
    a. Sample Information Summary Table
    
    b. Sample Information Table
    
    c. Sample Information Distribution Plots 
  
2. Counts Matrix Exploration (counts.R)
  
    Input - Normalized Counts Matrix (CSV)
    
    a. Count Filtering Summary Table
    
    b. Diagnostic Scatterplots
       
      i. Median Count vs. Variance
       
      ii. Medican Count vs. Number of Non-Zero Samples
    
    c. Clustered Heatmap of Counts vs. Samples
    
    d. PCA Biplot of User-Selected PCs
  
3. Differential Expression Analysis (diff_expr.R)

    Input - Differential Expression Results (CSV)
  
    a. Volcano Plot with User-Select Adjusted P-Value Threshold
  
    b. Table of DEGs with Significant Adjusted P-Value
  
4. Gene Set Enrichment Analysis (gsea.R)

    Input - FGSEA Results (C2: Canonical Pathway) (CSV)
    
    a. Barplot of Top Pathways by Normalized Enrichment Score (NES)
  
    b. Downloadable Table of Significantly Enriched Pathways
  
    c. Scatterplot of NES vs. -log10(Adjusted P-Value)

