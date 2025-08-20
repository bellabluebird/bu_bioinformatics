# expression analysis module for RNA-seq data visualization
# last updated: 08/20/25

# expression analysis module UI
expression_analysis_ui <- function(id) {
  ns <- NS(id)
  
  fluidPage(
    sidebarLayout(
      sidebarPanel(width = 3,
                   # file upload section - expects normalized counts matrix
                   wellPanel(
                     h4("Data Upload"),
                     fileInput(ns("counts_file"), "Choose Normalized Counts File:",
                               accept = c(".txt", ".tsv", ".csv")),
                     
                     # show dataset info once data is loaded
                     conditionalPanel(
                       condition = paste0("output['", ns("data_loaded"), "']"),
                       div(style = "background-color: #f8f9fa; padding: 10px; border-radius: 5px;",
                           h5("Dataset Info"),
                           textOutput(ns("data_info"))
                       )
                     )
                   ),
                   
                   # conditional panels for different analysis types
                   # each panel shows only when that specific tab is active
                   
                   # diagnostic plots options - for QC assessment
                   conditionalPanel(
                     condition = paste0("input['", ns("analysis_tabs"), "'] == 'diagnostics'"),
                     wellPanel(
                       h4("Diagnostic Options"),
                       p("Diagnostic plots help assess data quality and distribution."),
                       actionButton(ns("generate_diagnostics"), "Generate Diagnostics", 
                                    class = "btn-info", width = "100%")
                     )
                   ),
                   
                   # PCA options - for dimensionality reduction and sample clustering
                   conditionalPanel(
                     condition = paste0("input['", ns("analysis_tabs"), "'] == 'pca'"),
                     wellPanel(
                       h4("PCA Options"),
                       numericInput(ns("pca_genes"), "Genes for PCA:", 
                                    value = 1000, min = 500, max = 5000, step = 100),
                       p("PCA uses the most variable genes to identify patterns within samples."),
                       actionButton(ns("generate_pca"), "Run PCA Analysis", 
                                    class = "btn-success", width = "100%")
                     )
                   ),
                   
                   # heatmap options - for expression pattern visualization
                   conditionalPanel(
                     condition = paste0("input['", ns("analysis_tabs"), "'] == 'heatmap'"),
                     wellPanel(
                       h4("Heatmap Options"),
                       numericInput(ns("heatmap_genes"), "Number of genes:", 
                                    value = 50, min = 10, max = 200, step = 10),
                       checkboxInput(ns("cluster_rows"), "Cluster genes", TRUE),
                       checkboxInput(ns("cluster_cols"), "Cluster samples", TRUE),
                       checkboxInput(ns("scale_data"), "Z-score scaling", TRUE),
                       actionButton(ns("generate_heatmap"), "Generate Heatmap", 
                                    class = "btn-primary", width = "100%")
                     )
                   )
      ),
      
      mainPanel(width = 9,
                # show analysis tabs only when data is loaded
                conditionalPanel(
                  condition = paste0("output['", ns("data_loaded"), "']"),
                  
                  # main analysis tabs
                  tabsetPanel(
                    id = ns("analysis_tabs"),
                    
                    # diagnostic plots tab - quality control visualizations
                    tabPanel(
                      title = "Diagnostic Plots",
                      value = "diagnostics",
                      br(),
                      # show plots when diagnostics are ready
                      conditionalPanel(
                        condition = paste0("output['", ns("diagnostics_ready"), "']"),
                        fluidRow(
                          # library size distribution - shows sequencing depth variation
                          column(6, 
                                 h4("Library Sizes"),
                                 plotOutput(ns("library_size_plot"), height = "400px")
                          ),
                          # mean-variance relationship - shows count distribution properties
                          column(6, 
                                 h4("Mean-Variance Relationship"),
                                 plotOutput(ns("mean_var_plot"), height = "400px")
                          )
                        ),
                        br(),
                        fluidRow(
                          # sample correlation heatmap - shows sample similarity
                          column(12,
                                 h4("Sample Correlation Heatmap"),
                                 plotOutput(ns("correlation_heatmap"), height = "500px")
                          )
                        )
                      ),
                      # placeholder when no diagnostics generated yet
                      conditionalPanel(
                        condition = paste0("!output['", ns("diagnostics_ready"), "']"),
                        div(style = "text-align: center; margin-top: 50px;",
                            h4("Click 'Generate Diagnostics' to create quality control plots"),
                            p("These plots help assess library quality, count distributions, and sample relationships.")
                        )
                      )
                    ),
                    
                    # PCA analysis tab - principal component analysis for sample clustering
                    tabPanel(
                      title = "PCA Analysis", 
                      value = "pca",
                      br(),
                      # show PCA results when analysis is complete
                      conditionalPanel(
                        condition = paste0("output['", ns("pca_ready"), "']"),
                        fluidRow(
                          # interactive PCA plot - shows sample relationships in 2D space
                          column(6,
                                 h4("PCA Plot"),
                                 plotlyOutput(ns("pca_plot"), height = "500px")
                          ),
                          column(6,
                                 # variance explained plot - shows contribution of each PC
                                 h4("Variance Explained"),
                                 plotOutput(ns("variance_plot"), height = "300px"),
                                 br(),
                                 # top contributing genes - shows which genes drive PC1
                                 h5("Top Contributing Genes"),
                                 DT::dataTableOutput(ns("pca_loadings"))
                          )
                        )
                      ),
                      # placeholder when no PCA generated yet
                      conditionalPanel(
                        condition = paste0("!output['", ns("pca_ready"), "']"),
                        div(style = "text-align: center; margin-top: 50px;",
                            h4("Click 'Run PCA Analysis' to identify sample patterns"),
                            p("PCA reveals the main sources of variation across samples using the most variable genes.")
                        )
                      )
                    ),
                    
                    # heatmap tab - expression pattern visualization
                    tabPanel(
                      title = "Expression Heatmap",
                      value = "heatmap", 
                      br(),
                      # show heatmap when generated
                      conditionalPanel(
                        condition = paste0("output['", ns("heatmap_ready"), "']"),
                        h4("Expression Heatmap"),
                        # bordered container for the heatmap
                        div(style = "border: 1px solid #ddd; padding: 10px;",
                            plotOutput(ns("heatmap_plot"), height = "600px")),
                        br(),
                        # table showing genes displayed in heatmap
                        h4("Displayed Genes"),
                        DT::dataTableOutput(ns("heatmap_genes_table"))
                      ),
                      # placeholder when no heatmap generated yet
                      conditionalPanel(
                        condition = paste0("!output['", ns("heatmap_ready"), "']"),
                        div(style = "text-align: center; margin-top: 50px;",
                            h4("Click 'Generate Heatmap' to visualize expression patterns"),
                            p("Heatmaps show expression levels of the most variable genes across all samples.")
                        )
                      )
                    )
                  )
                ),
                
                # no data loaded state - landing page with instructions
                conditionalPanel(
                  condition = paste0("!output['", ns("data_loaded"), "']"),
                  div(
                    style = "text-align: center; margin-top: 100px;",
                    h3("Upload Normalized Counts File"),
                    p("Expected file: GSE64810_mlhd_DESeq2_norm_counts_adjust.txt"),
                    br(),
                    div(style = "color: #28a745;",
                        h5("This module provides:"),
                        tags$ul(style = "text-align: left; display: inline-block;",
                                tags$li("Quality control diagnostic plots"),
                                tags$li("PCA analysis for sample clustering"),
                                tags$li("Expression heatmaps with clustering")
                        )
                    )
                  )
                )
      )
    )
  )
}

# expression analysis module server logic
expression_analysis_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    
    # load and validate counts data from uploaded file
    counts_data <- reactive({
      req(input$counts_file)
      tryCatch({
        load_normalized_counts(input$counts_file$datapath)
      }, error = function(e) {
        showNotification(paste("Error loading file:", e$message), type = "error")
        return(NULL)
      })
    })
    
    # reactive flag to check if data is successfully loaded
    output$data_loaded <- reactive({
      !is.null(counts_data())
    })
    outputOptions(output, "data_loaded", suspendWhenHidden = FALSE)
    
    # display basic dataset information (genes x samples)
    output$data_info <- renderText({
      req(counts_data())
      paste0("Genes: ", format(nrow(counts_data()), big.mark = ","), "\n",
             "Samples: ", ncol(counts_data()))
    })
    
    # DIAGNOSTIC PLOTS SECTION
    # generates QC plots when button is clicked
    diagnostic_results <- eventReactive(input$generate_diagnostics, {
      req(counts_data())
      withProgress(message = "Generating diagnostic plots...", value = 0.5, {
        create_diagnostic_plots(counts_data())
      })
    })
    
    # reactive flag for diagnostic plot availability
    output$diagnostics_ready <- reactive({
      !is.null(diagnostic_results())
    })
    outputOptions(output, "diagnostics_ready", suspendWhenHidden = FALSE)
    
    # render library size distribution plot
    output$library_size_plot <- renderPlot({
      req(diagnostic_results())
      diagnostic_results()$library_size
    })
    
    # render mean-variance relationship plot
    output$mean_var_plot <- renderPlot({
      req(diagnostic_results())
      diagnostic_results()$mean_variance
    })
    
    # render sample correlation heatmap using ComplexHeatmap
    output$correlation_heatmap <- renderPlot({
      req(diagnostic_results())
      
      library(ComplexHeatmap)
      library(circlize)
      
      cor_matrix <- diagnostic_results()$sample_correlation
      
      # create correlation heatmap with custom color scheme
      Heatmap(
        cor_matrix,
        name = "Correlation",
        col = colorRamp2(c(0.7, 0.85, 1), c("white", "lightblue", "darkblue")),
        show_row_names = TRUE,
        show_column_names = TRUE,
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 8),
        column_title = "Sample-Sample Correlation Matrix",
        heatmap_legend_param = list(direction = "vertical")
      )
    })
    
    # PCA ANALYSIS SECTION
    # performs PCA analysis when button is clicked
    pca_results <- eventReactive(input$generate_pca, {
      req(counts_data())
      withProgress(message = "Performing PCA analysis...", value = 0.5, {
        perform_pca_analysis(counts_data(), n_genes = input$pca_genes)
      })
    })
    
    # reactive flag for PCA result availability
    output$pca_ready <- reactive({
      !is.null(pca_results())
    })
    outputOptions(output, "pca_ready", suspendWhenHidden = FALSE)
    
    # render interactive PCA plot using plotly
    output$pca_plot <- renderPlotly({
      req(pca_results())
      
      pca_plots <- create_pca_plots(pca_results())
      ggplotly(pca_plots$pca_2d, tooltip = "text")
    })
    
    # render variance explained plot - shows PC contribution
    output$variance_plot <- renderPlot({
      req(pca_results())
      
      pca_plots <- create_pca_plots(pca_results())
      pca_plots$variance_plot
    })
    
    # render table of top genes contributing to PC1
    output$pca_loadings <- DT::renderDataTable({
      req(pca_results())
      
      # create loadings dataframe with PC1 and PC2 contributions
      loadings_df <- data.frame(
        Gene = names(pca_results()$top_genes_pc1),
        PC1_Loading = round(pca_results()$top_genes_pc1, 3),
        PC2_Loading = round(pca_results()$top_genes_pc2[names(pca_results()$top_genes_pc1)], 3)
      )
      
      DT::datatable(
        loadings_df,
        options = list(pageLength = 5, dom = 't'),
        rownames = FALSE,
        caption = "Top 10 genes contributing to PC1"
      )
    })
    
    # HEATMAP SECTION  
    # prepares heatmap data when button is clicked
    heatmap_data <- eventReactive(input$generate_heatmap, {
      req(counts_data())
      
      withProgress(message = "Preparing heatmap...", value = 0.5, {
        # identify most variable genes for heatmap display
        gene_vars <- apply(counts_data(), 1, var, na.rm = TRUE)
        top_genes <- names(sort(gene_vars, decreasing = TRUE)[1:input$heatmap_genes])
        
        # subset to top variable genes and log transform
        matrix_subset <- counts_data()[top_genes, , drop = FALSE]
        matrix_subset <- log2(matrix_subset + 1)
        
        # apply z-score scaling if requested
        if (input$scale_data) {
          matrix_subset <- t(scale(t(matrix_subset)))
        }
        
        return(matrix_subset)
      })
    })
    
    # reactive flag for heatmap data availability
    output$heatmap_ready <- reactive({
      !is.null(heatmap_data())
    })
    outputOptions(output, "heatmap_ready", suspendWhenHidden = FALSE)
    
    # render expression heatmap using ComplexHeatmap
    output$heatmap_plot <- renderPlot({
      req(heatmap_data())
      
      library(ComplexHeatmap)
      library(circlize)
      
      # set color scheme based on scaling option
      if (input$scale_data) {
        col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
        legend_name <- "Z-score"
      } else {
        max_val <- max(heatmap_data(), na.rm = TRUE)
        col_fun <- colorRamp2(c(0, max_val/2, max_val), c("white", "pink", "red"))
        legend_name <- "log2(counts+1)"
      }
      
      # create heatmap with clustering options
      ht <- Heatmap(
        heatmap_data(),
        name = legend_name,
        col = col_fun,
        cluster_rows = input$cluster_rows,
        cluster_columns = input$cluster_cols,
        show_row_names = TRUE,
        show_column_names = TRUE,
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 9),
        column_title = paste("Expression Heatmap -", nrow(heatmap_data()), "most variable genes")
      )
      
      draw(ht)
    })
    
    # render table of genes displayed in heatmap with expression stats
    output$heatmap_genes_table <- DT::renderDataTable({
      req(heatmap_data())
      
      # create summary table for genes in heatmap
      genes_df <- data.frame(
        Gene = rownames(heatmap_data()),
        Mean_Expression = round(rowMeans(heatmap_data()), 2),
        Variance = round(apply(heatmap_data(), 1, var), 2)
      )
      
      DT::datatable(
        genes_df,
        options = list(pageLength = 10, scrollX = TRUE),
        rownames = FALSE,
        caption = paste("Genes displayed in heatmap:", nrow(genes_df))
      )
    })
  })
}
