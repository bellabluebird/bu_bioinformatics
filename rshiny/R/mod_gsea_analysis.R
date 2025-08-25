# R/mod_gsea_analysis.R

#' gsea Analysis UI
#' 
#' @param id module namespace ID
gsea_analysis_ui <- function(id) {
  ns <- NS(id)
  
  fluidPage(
    # file upload section
    fluidRow(
      column(6,
             wellPanel(
               h4("Upload Gene Data", icon("upload")),
               radioButtons(
                 ns("input_type"),
                 "Choose input type:",
                 choices = list(
                   "DESeq2 Results File (.csv/.txt)" = "deseq",
                   "Pre-ranked Genes File (.rnk/.txt)" = "ranked"
                 ),
                 selected = "deseq"
               ),
               
               # DESeq2 file input
               conditionalPanel(
                 condition = paste0("input['", ns("input_type"), "'] == 'deseq'"),
                 fileInput(
                   ns("deseq_file"), 
                   "Choose DESeq2 Results File",
                   accept = c(".csv", ".txt", ".tsv"),
                   placeholder = "deseq_results.csv"
                 ),
                 helpText("Expected columns: symbol, log2FoldChange, pvalue, padj, stat")
               ),
               
               # pre-ranked file input
               conditionalPanel(
                 condition = paste0("input['", ns("input_type"), "'] == 'ranked'"),
                 fileInput(
                   ns("ranked_genes_file"), 
                   "Choose Ranked Genes File",
                   accept = c(".txt", ".tsv", ".csv", ".rnk"),
                   placeholder = "ranked_genes.rnk"
                 ),
                 helpText("Two columns: gene symbol and ranking score")
               ),
               
               div(id = ns("gene_upload_status"))
             )
      ),
      column(6,
             wellPanel(
               h4("Upload Gene Sets File", icon("upload")),
               fileInput(
                 ns("gmt_file"), 
                 "Choose Gene Sets File (.gmt)",
                 accept = c(".gmt"),
                 placeholder = "gene_sets.gmt"
               ),
               div(id = ns("gmt_upload_status"))
             )
      )
    ),
    
    # DESeq2 processing parameters
    conditionalPanel(
      condition = paste0("input['", ns("input_type"), "'] == 'deseq' && output['", ns("deseq_uploaded"), "']"),
      fluidRow(
        column(12,
               wellPanel(
                 h4("⚙️ DESeq2 Processing Parameters", icon("cogs")),
                 fluidRow(
                   column(3,
                          selectInput(ns("ranking_method"), "Ranking Method",
                                      choices = list(
                                        "Wald Statistic (Recommended)" = "wald",
                                        "Log2 Fold Change" = "log2fc",
                                        "Signed -log10(padj)" = "signed_padj"
                                      ),
                                      selected = "wald")
                   ),
                   column(3,
                          numericInput(ns("padj_threshold"), "P-adj Threshold", 
                                       value = 0.05, min = 0.001, max = 1, step = 0.01)
                   ),
                   column(3,
                          checkboxInput(ns("filter_by_padj"), "Filter by P-adj", value = TRUE)
                   ),
                   column(3,
                          actionButton(ns("process_deseq"), "Process DESeq2 Data", 
                                       class = "btn-warning", style = "margin-top: 25px;")
                   )
                 ),
                 conditionalPanel(
                   condition = paste0("output['", ns("deseq_processed"), "']"),
                   hr(),
                   fluidRow(
                     column(12,
                            div(style = "background-color: #d4edda; padding: 10px; border-radius: 5px;",
                                h5("DESeq2 Data Processed Successfully"),
                                textOutput(ns("processing_summary"))
                            )
                     )
                   )
                 )
               )
        )
      )
    ),
    
    # standard GSEA parameters
    conditionalPanel(
      condition = paste0("output['", ns("files_ready"), "']"),
      fluidRow(
        column(12,
               wellPanel(
                 h4("⚙️ GSEA Parameters", icon("cogs")),
                 fluidRow(
                   column(3,
                          numericInput(ns("min_size"), "Min Gene Set Size", value = 15, min = 5)
                   ),
                   column(3,
                          numericInput(ns("max_size"), "Max Gene Set Size", value = 500, min = 50)
                   ),
                   column(3,
                          numericInput(ns("nperm"), "Permutations", value = 10000, min = 1000)
                   ),
                   column(3,
                          br(),
                          actionButton(ns("run_gsea"), "Run GSEA", class = "btn-primary btn-lg")
                   )
                 )
               )
        )
      )
    ),
    
    # data preview section
    conditionalPanel(
      condition = paste0("output['", ns("show_preview"), "']"),
      fluidRow(
        column(12,
               wellPanel(
                 h4("👀 Data Preview", icon("eye")),
                 tabsetPanel(
                   tabPanel("Input Data", DT::dataTableOutput(ns("input_preview"))),
                   tabPanel("Ranked Genes", DT::dataTableOutput(ns("ranked_preview")))
                 )
               )
        )
      )
    ),
    
    # main content (shown after analysis)
    conditionalPanel(
      condition = paste0("output['", ns("analysis_complete"), "']"),
      
      # summary cards
      fluidRow(
        column(3,
               div(class = "card bg-primary text-white mb-3",
                   div(class = "card-body",
                       h5(class = "card-title", "Total Pathways"),
                       h3(class = "card-text", textOutput(ns("total_pathways")))
                   )
               )
        ),
        column(3,
               div(class = "card bg-success text-white mb-3",
                   div(class = "card-body",
                       h5(class = "card-title", "Significant"),
                       h3(class = "card-text", textOutput(ns("sig_pathways")))
                   )
               )
        ),
        column(3,
               div(class = "card bg-danger text-white mb-3",
                   div(class = "card-body",
                       h5(class = "card-title", "Upregulated"),
                       h3(class = "card-text", textOutput(ns("up_pathways")))
                   )
               )
        ),
        column(3,
               div(class = "card bg-info text-white mb-3",
                   div(class = "card-body",
                       h5(class = "card-title", "Downregulated"),
                       h3(class = "card-text", textOutput(ns("down_pathways")))
                   )
               )
        )
      ),
      
      # main visualization tabs
      tabsetPanel(
        
        # results overview tab
        tabPanel(
          "Results Overview",
          br(),
          fluidRow(
            column(6,
                   wellPanel(
                     h4("Volcano Plot", icon("chart-line")),
                     plotlyOutput(ns("volcano_plot"), height = "400px")
                   )
            ),
            column(6,
                   wellPanel(
                     h4("Top Pathways", icon("chart-bar")),
                     plotlyOutput(ns("top_pathways_plot"), height = "400px")
                   )
            )
          ),
          fluidRow(
            column(12,
                   wellPanel(
                     h4("Pathway Summary Statistics", icon("table")),
                     DT::dataTableOutput(ns("summary_table"))
                   )
            )
          )
        ),
        
        # enrichment plots tab
        tabPanel(
          "Enrichment Plots",
          br(),
          fluidRow(
            column(4,
                   wellPanel(
                     h4("Select Pathway", icon("list")),
                     selectizeInput(
                       ns("selected_pathway"),
                       "Choose pathway for enrichment plot:",
                       choices = NULL,
                       options = list(placeholder = "Select a pathway...")
                     ),
                     br(),
                     div(id = ns("pathway_info"))
                   )
            ),
            column(8,
                   wellPanel(
                     h4("Enrichment Plot", icon("line-chart")),
                     plotlyOutput(ns("enrichment_plot"), height = "450px")
                   )
            )
          ),
          fluidRow(
            column(12,
                   wellPanel(
                     h4("Leading Edge Genes", icon("dna")),
                     DT::dataTableOutput(ns("leading_edge_table"))
                   )
            )
          )
        ),
        
        # results table tab
        tabPanel(
          "Results Table",
          br(),
          fluidRow(
            column(12,
                   wellPanel(
                     h4("Complete GSEA Results", icon("table")),
                     div(style = "margin-top: 10px;",
                         DT::dataTableOutput(ns("results_table"))
                     )
                   )
            )
          )
        ),
        
        # download results tab
        tabPanel(
          "Download",
          br(),
          fluidRow(
            column(6,
                   wellPanel(
                     h4("Download Results", icon("download")),
                     br(),
                     downloadButton(ns("download_results"), "Download Results Table", 
                                    class = "btn-primary", style = "width: 100%; margin-bottom: 10px;"),
                     br(),
                     downloadButton(ns("download_leading_edge"), "Download Leading Edge Genes", 
                                    class = "btn-info", style = "width: 100%; margin-bottom: 10px;"),
                     br(),
                     downloadButton(ns("download_ranked_genes"), "Download Ranked Genes", 
                                    class = "btn-warning", style = "width: 100%;")
                   )
            ),
            column(6,
                   wellPanel(
                     h4("Analysis Summary", icon("info-circle")),
                     verbatimTextOutput(ns("analysis_summary"))
                   )
            )
          )
        )
      )
    ),
    
    # loading/status messages
    conditionalPanel(
      condition = paste0("!output['", ns("files_ready"), "']"),
      div(class = "text-center", style = "margin-top: 50px;",
          h4("Upload gene data and gene sets file to begin"),
          conditionalPanel(
            condition = paste0("input['", ns("input_type"), "'] == 'deseq'"),
            p("Upload your DESeq2 results file (CSV/TXT) with columns: symbol, log2FoldChange, pvalue, padj, stat"),
            p("The system will automatically process it into a ranked gene list for GSEA.")
          ),
          conditionalPanel(
            condition = paste0("input['", ns("input_type"), "'] == 'ranked'"),
            p("Upload a pre-ranked genes file with gene symbols and ranking scores."),
            p("Gene sets file should be in GMT format (e.g., from MSigDB).")
          )
      )
    ),
    
    conditionalPanel(
      condition = paste0("output['", ns("files_ready"), "'] && !output['", ns("analysis_complete"), "']"),
      div(class = "text-center", style = "margin-top: 50px;",
          h4("⚙️ Configure parameters and click 'Run GSEA' to start analysis"),
          p("Review the parameter settings above before running the analysis.")
      )
    )
  )
}

#' GSEA Analysis Server
#' 
#' @param id Module namespace ID
gsea_analysis_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    
    # load required libraries
    library(fgsea)
    library(ggplot2)
    library(plotly)
    library(DT)
    library(dplyr)
    
    # reactive values
    values <- reactiveValues(
      deseq_data = NULL,
      ranked_genes = NULL,
      pathways = NULL,
      gsea_results = NULL,
      deseq_uploaded = FALSE,
      deseq_processed = FALSE,
      files_ready = FALSE,
      analysis_complete = FALSE,
      processing_info = NULL
    )
    
    # helper function to process DESeq2 data into ranked genes
    process_deseq_to_ranked <- function(deseq_df, method = "wald", filter_padj = TRUE, padj_thresh = 0.05) {
      
      # validate required columns
      required_cols <- c("symbol")
      if (method == "wald") required_cols <- c(required_cols, "stat")
      if (method == "log2fc") required_cols <- c(required_cols, "log2FoldChange")
      if (method == "signed_padj") required_cols <- c(required_cols, "log2FoldChange", "padj")
      if (filter_padj) required_cols <- c(required_cols, "padj")
      
      missing_cols <- setdiff(required_cols, colnames(deseq_df))
      if (length(missing_cols) > 0) {
        stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
      }
      
      # clean data
      df_clean <- deseq_df %>%
        filter(!is.na(symbol) & symbol != "") %>%
        filter(!duplicated(symbol))
      
      # apply p-adj filtering if requested
      if (filter_padj) {
        original_count <- nrow(df_clean)
        df_clean <- df_clean %>%
          filter(!is.na(padj) & padj <= padj_thresh)
        filtered_count <- nrow(df_clean)
        filter_info <- paste0("Filtered from ", original_count, " to ", filtered_count, " genes (padj ≤ ", padj_thresh, ")")
      } else {
        filter_info <- paste0("No p-adj filtering applied. ", nrow(df_clean), " genes retained.")
      }
      
      # create ranking metric based on method
      if (method == "wald") {
        df_clean <- df_clean %>%
          filter(!is.na(stat)) %>%
          mutate(rank_metric = stat)
        method_description <- "Wald statistic"
        
      } else if (method == "log2fc") {
        df_clean <- df_clean %>%
          filter(!is.na(log2FoldChange)) %>%
          mutate(rank_metric = log2FoldChange)
        method_description <- "Log2 fold change"
        
      } else if (method == "signed_padj") {
        df_clean <- df_clean %>%
          filter(!is.na(padj) & !is.na(log2FoldChange)) %>%
          filter(padj > 0) %>%  # Avoid infinite values
          mutate(rank_metric = -log10(padj) * sign(log2FoldChange))
        method_description <- "Signed -log10(padj)"
      }
      
      # create named vector and sort
      ranked_genes <- setNames(df_clean$rank_metric, df_clean$symbol)
      ranked_genes <- sort(ranked_genes, decreasing = TRUE)
      
      # return results with metadata
      list(
        ranked_genes = ranked_genes,
        method = method_description,
        filter_info = filter_info,
        final_count = length(ranked_genes)
      )
    }
    
    # file upload observers
    observeEvent(input$deseq_file, {
      req(input$deseq_file, input$input_type == "deseq")
      
      tryCatch({
        file_ext <- tools::file_ext(input$deseq_file$name)
        
        if (file_ext %in% c("csv")) {
          deseq_df <- read.csv(input$deseq_file$datapath, stringsAsFactors = FALSE)
        } else {
          deseq_df <- read.table(input$deseq_file$datapath, header = TRUE, 
                                 sep = "\t", stringsAsFactors = FALSE)
        }
        
        values$deseq_data <- deseq_df
        values$deseq_uploaded <- TRUE
        values$deseq_processed <- FALSE
        values$files_ready <- FALSE
        
        showNotification(
          paste("Loaded DESeq2 data with", nrow(deseq_df), "genes"),
          type = "message",
          duration = 3
        )
        
      }, error = function(e) {
        showNotification(
          paste("Error loading DESeq2 file:", e$message),
          type = "error",
          duration = 5
        )
      })
    })
    
    observeEvent(input$ranked_genes_file, {
      req(input$ranked_genes_file, input$input_type == "ranked")
      
      tryCatch({
        file_ext <- tools::file_ext(input$ranked_genes_file$name)
        
        if (file_ext %in% c("csv")) {
          ranks_df <- read.csv(input$ranked_genes_file$datapath, stringsAsFactors = FALSE)
        } else {
          ranks_df <- read.table(input$ranked_genes_file$datapath, header = FALSE, 
                                 sep = "\t", stringsAsFactors = FALSE)
        }
        
        # create named vector of ranks
        ranks <- setNames(ranks_df[[2]], ranks_df[[1]])
        ranks <- sort(ranks, decreasing = TRUE)
        
        values$ranked_genes <- ranks
        values$deseq_processed <- TRUE 
        
        showNotification(
          paste("Loaded", length(ranks), "ranked genes"),
          type = "message",
          duration = 3
        )
        
        check_files_ready()
        
      }, error = function(e) {
        showNotification(
          paste("Error loading ranked genes file:", e$message),
          type = "error",
          duration = 5
        )
      })
    })
    
    observeEvent(input$gmt_file, {
      req(input$gmt_file)
      
      tryCatch({
        pathways <- gmtPathways(input$gmt_file$datapath)
        values$pathways <- pathways
        
        showNotification(
          paste("Loaded", length(pathways), "gene sets"),
          type = "message",
          duration = 3
        )
        
        check_files_ready()
        
      }, error = function(e) {
        showNotification(
          paste("Error loading GMT file:", e$message),
          type = "error",
          duration = 5
        )
      })
    })
    
    # process DESeq2 data
    observeEvent(input$process_deseq, {
      req(values$deseq_data)
      
      withProgress(message = 'Processing DESeq2 data...', {
        
        tryCatch({
          incProgress(0.3, detail = "Validating data...")
          
          processed_result <- process_deseq_to_ranked(
            values$deseq_data,
            method = input$ranking_method,
            filter_padj = input$filter_by_padj,
            padj_thresh = input$padj_threshold
          )
          
          incProgress(0.8, detail = "Creating ranked gene list...")
          
          values$ranked_genes <- processed_result$ranked_genes
          values$processing_info <- processed_result
          values$deseq_processed <- TRUE
          
          incProgress(1.0, detail = "Processing complete!")
          
          showNotification(
            paste("DESeq2 data processed!", processed_result$final_count, "genes ranked"),
            type = "message",
            duration = 5
          )
          
          check_files_ready()
          
        }, error = function(e) {
          showNotification(
            paste("Error processing DESeq2 data:", e$message),
            type = "error",
            duration = 5
          )
        })
      })
    })
    
    # check if files are ready for GSEA
    check_files_ready <- function() {
      values$files_ready <- values$deseq_processed && !is.null(values$pathways) && !is.null(values$ranked_genes)
    }
    
    # GSEA analysis
    observeEvent(input$run_gsea, {
      req(values$ranked_genes, values$pathways)
      
      withProgress(message = 'Running GSEA analysis...', {
        
        tryCatch({
          incProgress(0.2, detail = "Preparing data...")
          
          incProgress(0.3, detail = "Running enrichment analysis...")
          gsea_results <- fgsea(
            pathways = values$pathways,
            stats = values$ranked_genes,
            minSize = input$min_size,
            maxSize = input$max_size,
            nperm = input$nperm
          )
          
          incProgress(0.8, detail = "Processing results...")
          
          # sort results and add additional columns
          gsea_results <- gsea_results[order(gsea_results$NES, decreasing = TRUE)]
          gsea_results$significant <- ifelse(gsea_results$padj < 0.05, "Yes", "No")
          gsea_results$direction <- ifelse(gsea_results$NES > 0, "Upregulated", "Downregulated")
          
          values$gsea_results <- gsea_results
          values$analysis_complete <- TRUE
          
          # update pathway choices for enrichment plots
          sig_pathways <- gsea_results[gsea_results$padj < 0.05, ]$pathway
          all_pathways <- gsea_results$pathway
          pathway_choices <- c(sig_pathways, all_pathways[!all_pathways %in% sig_pathways])
          
          updateSelectizeInput(session, "selected_pathway", choices = pathway_choices)
          
          incProgress(1.0, detail = "Analysis complete!")
          
          showNotification(
            paste("GSEA analysis complete!", nrow(gsea_results), "pathways analyzed"),
            type = "message",
            duration = 5
          )
          
        }, error = function(e) {
          showNotification(
            paste("Error in GSEA analysis:", e$message),
            type = "error",
            duration = 5
          )
        })
      })
    })
    
    # outtput reactive indicators
    output$deseq_uploaded <- reactive({
      values$deseq_uploaded
    })
    outputOptions(output, "deseq_uploaded", suspendWhenHidden = FALSE)
    
    output$deseq_processed <- reactive({
      values$deseq_processed
    })
    outputOptions(output, "deseq_processed", suspendWhenHidden = FALSE)
    
    output$files_ready <- reactive({
      values$files_ready
    })
    outputOptions(output, "files_ready", suspendWhenHidden = FALSE)
    
    output$analysis_complete <- reactive({
      values$analysis_complete
    })
    outputOptions(output, "analysis_complete", suspendWhenHidden = FALSE)
    
    output$show_preview <- reactive({
      !is.null(values$deseq_data) || !is.null(values$ranked_genes)
    })
    outputOptions(output, "show_preview", suspendWhenHidden = FALSE)
    
    # processing summary
    output$processing_summary <- renderText({
      req(values$processing_info)
      paste(
        paste("Ranking method:", values$processing_info$method),
        values$processing_info$filter_info,
        paste("Final ranked genes:", values$processing_info$final_count),
        sep = " | "
      )
    })
    
    # data previews
    output$input_preview <- DT::renderDataTable({
      req(values$deseq_data)
      DT::datatable(
        head(values$deseq_data, 100),
        options = list(pageLength = 10, scrollX = TRUE),
        rownames = FALSE
      )
    })
    
    output$ranked_preview <- DT::renderDataTable({
      req(values$ranked_genes)
      
      ranked_df <- data.frame(
        Gene = names(values$ranked_genes),
        Score = values$ranked_genes,
        Rank = 1:length(values$ranked_genes)
      )
      
      DT::datatable(
        head(ranked_df, 100),
        options = list(pageLength = 10, scrollX = TRUE),
        rownames = FALSE
      ) %>%
        DT::formatRound(columns = c("Score"), digits = 4)
    })
    
    # reactive expressions
    gsea_results <- reactive({
      req(values$gsea_results)
      values$gsea_results
    })
    
    # summary statistics
    output$total_pathways <- renderText({
      req(gsea_results())
      nrow(gsea_results())
    })
    
    output$sig_pathways <- renderText({
      req(gsea_results())
      sum(gsea_results()$padj < 0.05, na.rm = TRUE)
    })
    
    output$up_pathways <- renderText({
      req(gsea_results())
      sum(gsea_results()$NES > 0 & gsea_results()$padj < 0.05, na.rm = TRUE)
    })
    
    output$down_pathways <- renderText({
      req(gsea_results())
      sum(gsea_results()$NES < 0 & gsea_results()$padj < 0.05, na.rm = TRUE)
    })
    
    # volcano plot
    output$volcano_plot <- renderPlotly({
      req(gsea_results())
      
      p <- ggplot(gsea_results(), aes(x = NES, y = -log10(padj))) +
        geom_point(aes(color = significant, text = pathway), alpha = 0.6) +
        scale_color_manual(values = c("No" = "grey", "Yes" = "red")) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
        geom_vline(xintercept = 0, linetype = "solid", color = "black") +
        labs(
          title = "GSEA Volcano Plot",
          x = "Normalized Enrichment Score (NES)",
          y = "-log10(Adjusted P-value)",
          color = "Significant"
        ) +
        theme_minimal()
      
      ggplotly(p, tooltip = c("text", "x", "y")) %>%
        layout(hovermode = "closest")
    })
    
    # top pathways plot
    output$top_pathways_plot <- renderPlotly({
      req(gsea_results())
      
      top_pathways <- gsea_results()[gsea_results()$padj < 0.05, ]
      
      if (nrow(top_pathways) == 0) {
        return(plotly_empty() %>% 
                 add_annotations(text = "No significant pathways", showarrow = FALSE))
      }
      
      top_pathways <- head(top_pathways, 20)
      
      # truncate long pathway names
      top_pathways$pathway_short <- ifelse(nchar(top_pathways$pathway) > 50,
                                           paste0(substr(top_pathways$pathway, 1, 47), "..."),
                                           top_pathways$pathway)
      
      p <- ggplot(top_pathways, aes(x = reorder(pathway_short, NES), y = NES)) +
        geom_bar(stat = "identity", aes(fill = direction, text = pathway)) +
        scale_fill_manual(values = c("Upregulated" = "#e74c3c", "Downregulated" = "#3498db")) +
        coord_flip() +
        labs(
          title = "Top Enriched Pathways",
          x = "Pathway",
          y = "Normalized Enrichment Score",
          fill = "Direction"
        ) +
        theme_minimal() +
        theme(axis.text.y = element_text(size = 8))
      
      ggplotly(p, tooltip = c("text", "y")) %>%
        layout(hovermode = "closest")
    })
    
    # enrichment plot for selected pathway
    output$enrichment_plot <- renderPlotly({
      req(input$selected_pathway, values$pathways, values$ranked_genes)
      
      pathway_name <- input$selected_pathway
      
      tryCatch({
        p <- plotEnrichment(values$pathways[[pathway_name]], values$ranked_genes) +
          labs(
            title = paste("Enrichment Plot:", pathway_name),
            subtitle = paste("Pathway genes:", length(values$pathways[[pathway_name]]))
          ) +
          theme_minimal()
        
        ggplotly(p)
        
      }, error = function(e) {
        plotly_empty() %>% 
          add_annotations(text = "Error generating enrichment plot", showarrow = FALSE)
      })
    })
    
    # summary table
    output$summary_table <- DT::renderDataTable({
      req(gsea_results())
      
      total_pathways <- nrow(gsea_results())
      sig_pathways <- sum(gsea_results()$padj < 0.05, na.rm = TRUE)
      up_pathways <- sum(gsea_results()$NES > 0 & gsea_results()$padj < 0.05, na.rm = TRUE)
      down_pathways <- sum(gsea_results()$NES < 0 & gsea_results()$padj < 0.05, na.rm = TRUE)
      mean_nes <- ifelse(sig_pathways > 0, 
                         round(mean(gsea_results()[gsea_results()$padj < 0.05, ]$NES, na.rm = TRUE), 3),
                         0)
      median_size <- round(median(gsea_results()$size, na.rm = TRUE), 0)
      
      summary_stats <- data.frame(
        Metric = c("Total pathways tested", "Significantly enriched (padj < 0.05)", 
                   "Upregulated pathways", "Downregulated pathways",
                   "Mean NES (significant)", "Median pathway size"),
        Value = c(
          total_pathways,
          sig_pathways,
          up_pathways,
          down_pathways,
          mean_nes,
          median_size
        )
      )
      
      DT::datatable(summary_stats, options = list(dom = 't', pageLength = 10), rownames = FALSE)
    })
    
    # leading edge table
    output$leading_edge_table <- DT::renderDataTable({
      req(input$selected_pathway, gsea_results())
      
      pathway_result <- gsea_results()[gsea_results()$pathway == input$selected_pathway, ]
      
      if (nrow(pathway_result) == 0) {
        return(data.frame(Message = "No data available"))
      }
      
      leading_edge <- pathway_result$leadingEdge[[1]]
      if (length(leading_edge) == 0) {
        return(data.frame(Message = "No leading edge genes"))
      }
      
      # get gene scores
      gene_scores <- values$ranked_genes[leading_edge]
      
      leading_edge_df <- data.frame(
        Gene = leading_edge,
        Score = gene_scores,
        Rank = match(leading_edge, names(values$ranked_genes))
      ) %>%
        arrange(desc(Score))
      
      DT::datatable(
        leading_edge_df,
        options = list(pageLength = 15, scrollX = TRUE),
        rownames = FALSE
      ) %>%
        DT::formatRound(columns = c("Score"), digits = 3)
    })
    
    # results table
    output$results_table <- DT::renderDataTable({
      req(gsea_results())
      
      # prepare table data (exclude leadingEdge column for display)
      display_results <- gsea_results() %>%
        select(-leadingEdge) %>%
        arrange(padj)
      
      DT::datatable(
        display_results,
        options = list(
          pageLength = 25,
          scrollX = TRUE,
          dom = 'Bfrtip',
          buttons = c('copy', 'csv', 'excel')
        ),
        extensions = 'Buttons',
        rownames = FALSE
      ) %>%
        DT::formatRound(columns = c("pval", "padj", "ES", "NES"), digits = 4) %>%
        DT::formatStyle(
          "padj",
          backgroundColor = DT::styleInterval(0.05, c("lightcoral", "white"))
        )
    })
    
    # analysis summary
    output$analysis_summary <- renderText({
      req(gsea_results())
      
      total <- nrow(gsea_results())
      sig <- sum(gsea_results()$padj < 0.05, na.rm = TRUE)
      up <- sum(gsea_results()$NES > 0 & gsea_results()$padj < 0.05, na.rm = TRUE)
      down <- sum(gsea_results()$NES < 0 & gsea_results()$padj < 0.05, na.rm = TRUE)
      
      # add DESeq2 processing info if available
      processing_text <- ""
      if (!is.null(values$processing_info)) {
        processing_text <- paste(
          "",
          "DESeq2 Processing:",
          paste("- Ranking method:", values$processing_info$method),
          paste("- Filter info:", values$processing_info$filter_info),
          paste("- Genes used for GSEA:", values$processing_info$final_count),
          "",
          sep = "\n"
        )
      }
      
      paste(
        "GSEA Analysis Summary:",
        "===================",
        paste("Total pathways analyzed:", total),
        paste("Significantly enriched pathways:", sig, paste0("(", round(sig/total*100, 1), "%)")),
        paste("Upregulated pathways:", up),
        paste("Downregulated pathways:", down),
        processing_text,
        "GSEA Parameters:",
        paste("- Minimum gene set size:", input$min_size),
        paste("- Maximum gene set size:", input$max_size),
        paste("- Number of permutations:", input$nperm),
        paste("- Significance threshold: 0.05 (adjusted p-value)"),
        sep = "\n"
      )
    })
    
    # download handlers
    output$download_results <- downloadHandler(
      filename = function() {
        paste0("gsea_results_", Sys.Date(), ".txt")
      },
      content = function(file) {
        write.table(gsea_results(), file, sep = "\t", row.names = FALSE, quote = FALSE)
      }
    )
    
    output$download_leading_edge <- downloadHandler(
      filename = function() {
        paste0("leading_edge_genes_", Sys.Date(), ".txt")
      },
      content = function(file) {
        req(input$selected_pathway)
        
        pathway_result <- gsea_results()[gsea_results()$pathway == input$selected_pathway, ]
        leading_edge <- pathway_result$leadingEdge[[1]]
        gene_scores <- values$ranked_genes[leading_edge]
        
        leading_edge_df <- data.frame(
          Pathway = input$selected_pathway,
          Gene = leading_edge,
          Score = gene_scores,
          Rank = match(leading_edge, names(values$ranked_genes))
        )
        
        write.table(leading_edge_df, file, sep = "\t", row.names = FALSE, quote = FALSE)
      }
    )
    
    # download ranked genes
    output$download_ranked_genes <- downloadHandler(
      filename = function() {
        method_suffix <- switch(input$ranking_method,
                                "wald" = "wald_stat",
                                "log2fc" = "log2fc", 
                                "signed_padj" = "signed_padj")
        paste0("ranked_genes_", method_suffix, "_", Sys.Date(), ".rnk")
      },
      content = function(file) {
        req(values$ranked_genes)
        
        ranked_df <- data.frame(
          Gene = names(values$ranked_genes),
          Score = values$ranked_genes
        )
        
        write.table(ranked_df, file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
      }
    )
    
  })
}
