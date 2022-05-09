## -- BF591 Final Project -- ##
## RShiny App
## Author: Daniel Goldstein
## Email: djgoldst@bu.edu

library(shiny)
library(bslib)
library(colourpicker)

options(shiny.maxRequestSize=30*1024^2)

source("samples.R")
source("counts.R")
source("diff_expr.R")
source("gsea.R")

ui <- fluidPage(
  theme = bs_theme(version = 5, bootswatch = "zephyr"),
  titlePanel("BF591 Final Project"),
  h4("Author: Daniel Goldstein"),
  h6("Dataset: Post-mortem Huntington's Disease (HD) prefrontal cortex compared with neurologically healthy controls"),
  tabsetPanel(
    tabPanel( #Sample Information Tab
      title = "Sample",
      sidebarLayout(
        sidebarPanel( #Inputs
          fileInput(inputId = "sample_info_upload",
                    label = "Sample Information (CSV):",
                    accept = ".csv",
                    placeholder = "metadata.csv"),
          width = 3
        ),
        mainPanel( #Outputs
          tabsetPanel(
            tabPanel( #Sample Information Summary Table
              title = "Summary",
              tableOutput("sample_info_summary_table")
            ),
            tabPanel( #Sample Information Table
              title = "Table",
              dataTableOutput("sample_info_table")
            ),
            tabPanel( #Sample Summary Violin Plot
              title = "Plot",
              sidebarLayout(
                sidebarPanel(
                  radioButtons(inputId = "sample_info_summary_variable",
                               label = "Select a variable for summary plot:",
                               choices = c("Age_of_Death","mRNA_Seq_reads","Age_of_Onset","Duration",
                                           "CAG","Vonsattel_Grade","H_V_Striatal_Score","H_V_Cortical_Score"),
                               selected = "Age_of_Death"),
                  submitButton(text = "Submit",
                               width = "100%")
                ),
                mainPanel(
                  plotOutput("sample_info_plot")
                )
              )
            )
          )
        )
      )
    ),
    tabPanel( #Counts Exploration Tab
      title = "Counts",
      sidebarLayout(
        sidebarPanel( #Inputs
          fileInput(inputId = "counts_upload",
                    label = "Counts Matrix (CSV):",
                    accept = ".csv",
                    placeholder = "HD_norm_counts.csv"),
          sliderInput(inputId = "counts_percent_var",
                      label = "Filter genes above the following percentile variance:",
                      min = 0,
                      max = 100,
                      value = 50,
                      step = 1),
          sliderInput(inputId = "counts_nonzero",
                      label = "Filter genes with the following number of non-zero samples:",
                      min = 0,
                      max = 69,
                      value = 35,
                      step = 1),
          submitButton(text = "Submit",
                       width = "100%"),
          width = 3
        ),
        mainPanel( #Outputs
          tabsetPanel(
            tabPanel( #Summary Filtered Counts Table
              title = "Table",
              tableOutput("filtered_counts_summary_table")
            ),
            tabPanel( #Diagnostic Scatter Plots
              title = "Diagnostic Plots",
              tabsetPanel(
                tabPanel(
                  title = "Median vs. Variance",
                  plotOutput("median_var_plot")
                  ),
                tabPanel(
                  title = "Median vs. Non-Zero Samples",
                  plotOutput("median_nonzero_plot")
                  )
                ) 
            ),
            tabPanel( #Clustered Heatmap Plot
              title = "Heatmap",
              plotOutput("counts_heatmap") #clustered heatmap output name
            ),
            tabPanel( #PCA
              title = "PCA",
              sidebarLayout(
                sidebarPanel( #Allow user to select PCs
                  #TODO: Radiobutton or slider??
                  fileInput(inputId = "counts_metadata",
                            label = "Sample Metadata (CSV):",
                            accept = ".csv",
                            placeholder = "metadata.csv"),
                  sliderInput(inputId = "first_PC",
                              label = "Select the first principal component (x-axis):",
                              min = 1,
                              max = 69,
                              value = 1,
                              step = 1),
                  sliderInput(inputId = "second_PC",
                              label = "Select the second principal component (y-axis):",
                              min = 1,
                              max = 69,
                              value = 2,
                              step = 1),
                  submitButton(text = "Submit",
                               width = "100%")
                ),
                mainPanel( #PCA Biplot
                  plotOutput("counts_pca_biplot")
                )
              )
            )
          )
        )
      )
    ),
    tabPanel( #Differential Expression Tab
      title = "DE",
      sidebarLayout(
        sidebarPanel( #Inputs
          fileInput(inputId = "de_res_upload",
                    label = "Differential Expression Results (CSV):",
                    accept = ".csv",
                    placeholder = "HD_DESEq2_DE_res.csv"),
          radioButtons(inputId = "volc_plot_x",
                       label = "Choose the x-axis column:",
                       choices = c("baseMean","HD.mean","Control.mean","log2FoldChange",
                                   "lfcSE","stat","pvalue","padj"),
                       selected = "log2FoldChange"),
          radioButtons(inputId = "volc_plot_y",
                       label = "Choose the y-axis column:",
                       choices = c("baseMean","HD.mean","Control.mean","log2FoldChange",
                                   "lfcSE","stat","pvalue","padj"),
                       selected = "padj"),
          colourInput(inputId = "volc_plot_base_col",
                      label = "Base point color:",
                      value = "#000000"),
          colourInput(inputId = "volc_plot_high_col",
                     label = "Highlight point color:",
                     value = "#FF0900"),
          sliderInput(inputId = "volc_plot_padj_slider",
                      label = "Select the magnitude of the adjust p-value coloring threshold:",
                      min = -35,
                      max = 0,
                      value = -15,
                      step = 1),
          submitButton(text = "Submit",
                       width = "100%"),
          width = 3
        ),
        mainPanel( #Outputs
          tabsetPanel(
            tabPanel( #Volcano Plot
              title = "Plot",
              plotOutput("volcano_plot") 
            ),
            tabPanel( #DE Results Table
              title = "Table",
              dataTableOutput("de_summary") 
            )
          )
        )
      )
    ),
    tabPanel( #GSEA Tab
      title = "GSEA",
      sidebarLayout(
        sidebarPanel( #Inputs
          fileInput(inputId = "fgsea_res_upload",
                    label = "FGSEA Results (CSV):",
                    accept = ".csv",
                    placeholder = "fgsea_c2path_res.csv"),
          width = 3
        ),
        mainPanel( #Outputs
          tabsetPanel(
            tabPanel( #Barplot of fgsea NES for top pathways
              title = "Top Results",
              sidebarLayout(
                sidebarPanel( #Slider: Adjust number of top pathways
                  sliderInput(inputId = "fgsea_top_path_slider",
                              label = "Select number of top pathways:",
                              min = 0,
                              max = 25,
                              value = 10,
                              step = 1),
                  submitButton(text = "Submit",
                               width = "100%")
                ),
                mainPanel(
                  plotOutput("fgsea_NES_barplot")
                )
              )
            ),
            tabPanel( #Table of FGSEA Pathways
              title = "Table",
              sidebarLayout(
                sidebarPanel(
                  sliderInput(inputId = "fgsea_table_padj_slider",
                              label = "Select adjusted p-value filtering threshold:",
                              min = -35,
                              max = 0,
                              value = -4,
                              step = 1),
                  radioButtons(inputId = "fgsea_NES_direction",
                               label = "Select NES direction (all, positive, negative):",
                               choices = c("all","positive","negative"),
                               selected = "all"),
                  submitButton(text = "Submit",
                               width = "100%"),
                ),
                mainPanel(
                  downloadButton(outputId = "fgsea_table_download",
                                 label = "Download Table"),
                  dataTableOutput("fgsea_sum_table")
                )
              )
            ),
            tabPanel( #Scatterplot of NES vs. -log10(padj)
              title = "Plot",
              sidebarLayout(
                sidebarPanel( #Slider: adjusted p-value filter
                  sliderInput(inputId = "NES_padj_slider",
                              label = "Select adjusted p-value filtering threshold:",
                              min = -35,
                              max = 0,
                              value = -5,
                              step = 1),
                  submitButton(text = "Submit",
                               width = "100%")
                ),
                mainPanel(
                  plotOutput("fgsea_NES_scatter")
                )
              )
            )
          )
        )
      )
    )
  )
)




server <- function(input, output) {
  ## Sample Information Tab
  
  sample_data <- reactive({
    read.csv(file = input$sample_info_upload$datapath,row.names = 1) %>%
      setnames(old = colnames(.), new = gsub("\\.","_", colnames(.))) %>%
      mutate(PMI = as.integer(PMI),
             mRNA_Seq_reads = gsub(",","",mRNA_Seq_reads),
             mRNA_Seq_reads = as.integer(mRNA_Seq_reads),
             Condition = ifelse(str_detect(Sample_ID, "C"),
                                "Control",
                                "HD")) %>%
      relocate(Condition, .after="Sample_ID") %>%
      arrange(Sample_ID) %>%
      return()
  })
  
  # Sample Information Summary Table
  output$sample_info_summary_table <- renderTable({
    req(input$sample_info_upload$datapath)
    
    sample_summary_cols <- colnames(sample_data())[-c(1,3,5)]
    sample_summary_table(data = sample_data(),
                         cols = sample_summary_cols)
  })
  
  # Sample Information Table
  output$sample_info_table <- renderDataTable({
    req(input$sample_info_upload$datapath)
    
    #Labels not shown on RShiny
    label_sample_data(data = sample_data())
  })
  
  # Sample Summary Violin Plot
  output$sample_info_plot <- renderPlot({
    req(input$sample_info_upload$datapath)
    req(input$sample_info_summary_variable)
    
    control_data <- subset_data(data = sample_data(),
                                cond = "Control")
    HD_data <- subset_data(data = sample_data(),
                           cond = "HD")
    
    sample_summary_plot(data = sample_data(),
                        ctrl_data = control_data,
                        hd_data = HD_data,
                        sum_var = input$sample_info_summary_variable)
  })
  
  
  ## Counts Exploration Tab
  
  counts_data <- reactive({
    read.csv(input$counts_upload$datapath, row.names=1) %>%
      return()
  })
  
  counts_metadata <- reactive({
    read.csv(file = input$counts_metadata$datapath,row.names = 1) %>%
      setnames(old = colnames(.), new = gsub("\\.","_", colnames(.))) %>%
      mutate(PMI = as.integer(PMI),
             mRNA_Seq_reads = gsub(",","",mRNA_Seq_reads),
             mRNA_Seq_reads = as.integer(mRNA_Seq_reads),
             Condition = ifelse(str_detect(Sample_ID, "C"),
                                "Control",
                                "HD")) %>%
      relocate(Condition, .after="Sample_ID") %>%
      arrange(Sample_ID) %>%
      return()
  })
  
  # Filtered Counts Table
  output$filtered_counts_summary_table <- renderTable({
    req(input$counts_upload$datapath)
    
    counts_summary_table(counts = counts_data(),
                         perc_var = input$counts_percent_var,
                         num_nonzero = input$counts_nonzero)
  })
  
  #Diagnostic Scatter Plots
  output$median_var_plot <- renderPlot({
    req(input$counts_upload$datapath)
    
    diagnostic_plot(counts = counts_data(),
                    plot_type = "variance",
                    perc_var = input$counts_percent_var,
                    log_scale = TRUE)
  })
  
  output$median_nonzero_plot <- renderPlot({
    req(input$counts_upload$datapath)
    
    diagnostic_plot(counts = counts_data(),
                    plot_type = "non_zeros",
                    num_nonzero = input$counts_nonzero)
  })
  
  #Clustered Heatmap
  output$counts_heatmap <- renderPlot({
    req(input$counts_upload$datapath)
    
    filtered_counts <- filter_counts(counts = counts_data(),
                                 perc_var = input$counts_percent_var,
                                 num_nonzero = input$counts_nonzero)
    filt_counts_heatmap(filt_counts = filtered_counts)
    
  })
  
  #PCA Biplot
  output$counts_pca_biplot <- renderPlot({
    req(input$counts_upload$datapath)
    req(input$counts_metadata$datapath)
    
    filtered_counts <- filter_counts(counts = counts_data(),
                                     perc_var = input$counts_percent_var,
                                     num_nonzero = input$counts_nonzero)
    pca_biplot(filt_counts = filtered_counts,
               metadata = counts_metadata(),
               first_PC = paste0("PC",input$first_PC),
               second_PC = paste0("PC",input$second_PC)) 
  })
  
  
  ## Differential Expression Tab
  
  de_data <- reactive({
    read.csv(file = input$de_res_upload$datapath, row.names=1) %>%
      rownames_to_column("ENSGID") %>%
      return()
  })
  
  #DE Volcano Plot
  output$volcano_plot <- renderPlot({
    req(input$de_res_upload$datapath)
    
    de_volcano_plot(dataf = de_data(),
                    x_name = input$volc_plot_x,
                    y_name = input$volc_plot_y,
                    slider = input$volc_plot_padj_slider,
                    color1 = input$volc_plot_base_col,
                    color2 = input$volc_plot_high_col)
  })
  
  #DE Summary Table
  output$de_summary <- renderDataTable({
    req(input$de_res_upload$datapath)
    
    de_summary_table(dataf = de_data(),
                     slider = input$volc_plot_padj_slider)
  })
  
  
  ## GSEA Tab
  
  fgsea_data <- reactive({
    read.csv(input$fgsea_res_upload$datapath) %>%
      return()
  })
  
  #Top Pathways (NES barplot)
  output$fgsea_NES_barplot <- renderPlot({
    req(input$fgsea_res_upload$datapath)
    
    fgsea_NES_barplot(fgsea_results = fgsea_data(),
                      pathway_name = "C2: Canonical Pathway",
                      num_paths = input$fgsea_top_path_slider)
  })
  
  fgsea_sum_padj_slider <- reactive({
    input$fgsea_table_padj_slider
  })
  
  #FGSEA Summary Table
  output$fgsea_sum_table <- renderDataTable({
    req(input$fgsea_res_upload$datapath)
    
    fgsea_summary_table(fgsea_results = fgsea_data(),
                        pathway_direction = input$fgsea_NES_direction,
                        padj_filter = fgsea_sum_padj_slider())
  })
  
  #FGSEA Summary Table Download
  output$fgsea_table_download <- downloadHandler(
    filename = function() {"fgsea_c2pathway_res.csv"},
    content = function(file){
      write_csv(fgsea_summary_table(fgsea_results = fgsea_data(),
                                    pathway_direction = input$fgsea_NES_direction,
                                    padj_filter = input$fgsea_table_padj_slider), 
                file)
    }
  )
  
  fgsea_scatter_padj_slider <- reactive({
    input$NES_padj_slider
  })
  
  #FGSEA NES Scatterplot
  output$fgsea_NES_scatter <- renderPlot({
    req(input$fgsea_res_upload$datapath)
    
    fgsea_NES_scatterplot(fgsea_results = fgsea_data(),
                          padj_filter = fgsea_scatter_padj_slider())
  })
}


shinyApp(ui, server)
