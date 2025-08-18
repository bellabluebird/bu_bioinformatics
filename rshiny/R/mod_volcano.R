# R/mod_volcano.R
# volcano plot module for GSE64810 visualization

# volcano plot module UI
volcano_plot_ui <- function(id) {
  ns <- NS(id)
  
  fluidPage(    
    sidebarLayout(
      sidebarPanel(width = 3,
                   wellPanel(
                     h4("Data Upload"),
                     fileInput(ns("deseq_file"), "Choose DESeq2 Results File:",
                               accept = c(".txt", ".tsv", ".csv")),
                     
                     h4("Significance Thresholds"),
                     numericInput(ns("pval_threshold"), "P-value cutoff:",
                                  value = 0.05, min = 0.001, max = 0.1, step = 0.005),
                     numericInput(ns("fc_threshold"), "Fold change cutoff:",
                                  value = 1.5, min = 1, max = 5, step = 0.1),
                     
                     h4("Display Options"),
                     selectInput(ns("color_scheme"), "Color by:",
                                 choices = c("Regulation" = "regulation", 
                                             "Significance" = "significance")),
                     checkboxInput(ns("show_labels"), "Show gene labels", FALSE),
                     
                     hr(),
                     downloadButton(ns("download_genes"), "Download Selected", 
                                    class = "btn-primary")
                   )
      ),
      
      mainPanel(width = 9,
                conditionalPanel(
                  condition = paste0("output['", ns("data_loaded"), "']"),
                  plotlyOutput(ns("volcano_plot"), height = "600px"),
                  br(),
                  h4("Selected/Significant Genes"),
                  DT::dataTableOutput(ns("gene_table"))
                ),
                conditionalPanel(
                  condition = paste0("!output['", ns("data_loaded"), "']"),
                  div(
                    style = "text-align: center; margin-top: 100px;",
                    h3("Please upload your GSE64810 DESeq2 results file"),
                    p("Expected file: GSE64810_mlhd_DESeq2_diffexp_DESeq2_outlier_trimmed_adjust.txt")
                  )
                )
      )
    )
  )
}

# volcano plot module server
volcano_plot_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    
    # reactive data loading (named to avoid conflict with base::data)
    deseq_data <- reactive({
      req(input$deseq_file)
      
      tryCatch({
        load_deseq_results(input$deseq_file$datapath)
      }, error = function(e) {
        showNotification(paste("Error loading file:", e$message), type = "error")
        return(NULL)
      })
    })
    
    # check if data is loaded
    output$data_loaded <- reactive({
      !is.null(deseq_data())
    })
    outputOptions(output, "data_loaded", suspendWhenHidden = FALSE)
    
    # reactive data with applied thresholds
    filtered_data <- reactive({
      req(deseq_data())
      
      deseq_data() %>%
        mutate(
          is_significant = padj < input$pval_threshold & 
            abs(log2FoldChange) > log2(input$fc_threshold),
          regulation_filtered = case_when(
            log2FoldChange > log2(input$fc_threshold) & 
              padj < input$pval_threshold ~ "Up-regulated",
            log2FoldChange < -log2(input$fc_threshold) & 
              padj < input$pval_threshold ~ "Down-regulated",
            TRUE ~ "No change"
          )
        )
    })
    
    # interactive volcano plot with plotly
    output$volcano_plot <- renderPlotly({
      req(filtered_data())
      
      plot_data <- filtered_data()
      
      p <- ggplot(plot_data, aes(x = log2FoldChange, y = neg_log10_pval)) +
        geom_point(aes(color = get(input$color_scheme),
                       text = paste("Gene:", gene_id,
                                    "<br>Log2FC:", round(log2FoldChange, 3),
                                    "<br>P-value:", format(pvalue, scientific = TRUE),
                                    "<br>Adj. P-value:", format(padj, scientific = TRUE))),
                   alpha = 0.7, size = 1.2) +
        
        # add threshold lines
        geom_hline(yintercept = -log10(input$pval_threshold), 
                   linetype = "dashed", color = "red", alpha = 0.7) +
        geom_vline(xintercept = c(-log2(input$fc_threshold), 
                                  log2(input$fc_threshold)), 
                   linetype = "dashed", color = "red", alpha = 0.7) +
        
        scale_color_manual(values = c(
          "Up-regulated" = "#e74c3c", "Down-regulated" = "#3498db", 
          "No change" = "#95a5a6", "High significance" = "#e74c3c",
          "Medium significance" = "#f39c12", "Low significance" = "#f1c40f",
          "Not significant" = "#95a5a6"
        )) +
        
        labs(x = "Log2 Fold Change", y = "-Log10 P-value", 
             color = tools::toTitleCase(gsub("_", " ", input$color_scheme))) +
        theme_minimal() +
        theme(legend.position = "right")
      
      ggplotly(p, tooltip = "text") %>%
        layout(dragmode = "select") %>%
        config(displayModeBar = TRUE, displaylogo = FALSE)
    })
    
    # gene table for selected/significant genes
    output$gene_table <- DT::renderDataTable({
      req(filtered_data())
      
      # check for selected points from plotly
      selected_points <- event_data("plotly_selected")
      
      if (!is.null(selected_points)) {
        display_data <- filtered_data()[selected_points$pointNumber + 1, ]
      } else {
        # show only significant genes by default
        display_data <- filtered_data() %>%
          filter(is_significant) %>%
          arrange(padj)
      }
      
      display_data %>%
        select(gene_id, log2FoldChange, pvalue, padj, regulation_filtered) %>%
        DT::datatable(
          options = list(pageLength = 15, scrollX = TRUE),
          rownames = FALSE
        ) %>%
        DT::formatRound(columns = c("log2FoldChange", "pvalue", "padj"), digits = 4)
    })
    
    # download handler
    output$download_genes <- downloadHandler(
      filename = function() paste0("GSE64810_selected_genes_", Sys.Date(), ".csv"),
      content = function(file) {
        selected_points <- event_data("plotly_selected")
        
        if (!is.null(selected_points)) {
          download_data <- filtered_data()[selected_points$pointNumber + 1, ]
        } else {
          download_data <- filtered_data() %>% filter(is_significant)
        }
        
        write.csv(download_data, file, row.names = FALSE)
      }
    )
  })
}
