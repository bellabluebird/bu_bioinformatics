# R/mod_sample_metadata.R 

#' Sample Metadata UI
#' 
#' @param id Module namespace ID
sample_metadata_ui <- function(id) {
  ns <- NS(id)
  
  fluidPage(
    # file upload section
    fluidRow(
      column(12,
             wellPanel(
               h4("Upload Sample Metadata", icon("upload")),
               fileInput(
                 ns("metadata_file"), 
                 "Choose GSE64810 Series Matrix File (.txt) or Processed CSV",
                 accept = c(".txt", ".tsv", ".csv"),
                 placeholder = "GSE64810_series_matrix.txt or huntington_metadata.csv"
               ),
               div(id = ns("upload_status"))
             )
      )
    ),
    
    # main content (only shown after file upload)
    conditionalPanel(
      condition = paste0("output['", ns("file_uploaded"), "']"),
      
      # summary cards
      fluidRow(
        column(3,
               div(class = "card bg-primary text-white mb-3",
                   div(class = "card-body",
                       h5(class = "card-title", "Total Samples"),
                       h3(class = "card-text", textOutput(ns("total_samples")))
                   )
               )
        ),
        column(3,
               div(class = "card bg-success text-white mb-3",
                   div(class = "card-body",
                       h5(class = "card-title", "HD Patients"),
                       h3(class = "card-text", textOutput(ns("hd_samples")))
                   )
               )
        ),
        column(3,
               div(class = "card bg-info text-white mb-3",
                   div(class = "card-body",
                       h5(class = "card-title", "Controls"),
                       h3(class = "card-text", textOutput(ns("control_samples")))
                   )
               )
        ),
        column(3,
               div(class = "card bg-warning text-white mb-3",
                   div(class = "card-body",
                       h5(class = "card-title", "Variables"),
                       h3(class = "card-text", textOutput(ns("total_variables")))
                   )
               )
        )
      ),
      
      # main visualization tabs
      tabsetPanel(
        
        # sample overview tab
        tabPanel(
          "Sample Overview",
          br(),
          fluidRow(
            column(6,
                   wellPanel(
                     h4("Age Distribution", icon("chart-line")),
                     plotlyOutput(ns("age_plot"), height = "350px")
                   )
            ),
            column(6,
                   wellPanel(
                     h4("RNA Quality (RIN)", icon("dna")),
                     plotlyOutput(ns("rin_plot"), height = "350px")
                   )
            )
          ),
          fluidRow(
            column(6,
                   wellPanel(
                     h4("Post-Mortem Interval", icon("clock")),
                     plotlyOutput(ns("pmi_plot"), height = "350px")
                   )
            ),
            column(6,
                   wellPanel(
                     h4("Sequencing Depth", icon("microscope")),
                     plotlyOutput(ns("reads_plot"), height = "350px")
                   )
            )
          )
        ),
        
        # hd-specific tab
        tabPanel(
          "HD Analysis",
          br(),
          conditionalPanel(
            condition = paste0("output['", ns("has_hd_data"), "']"),
            fluidRow(
              column(6,
                     wellPanel(
                       h4("Disease Progression Overview", icon("project-diagram")),
                       plotlyOutput(ns("progression_plot"), height = "400px")
                     )
              ),
              column(6,
                     wellPanel(
                       h4("Age of Onset vs CAG Repeats", icon("scatter-chart")),
                       plotlyOutput(ns("onset_cag_plot"), height = "400px")
                     )
              )
            ),
            fluidRow(
              column(4,
                     wellPanel(
                       h4("Disease Duration", icon("hourglass-half")),
                       plotlyOutput(ns("duration_plot"), height = "300px")
                     )
              ),
              column(4,
                     wellPanel(
                       h4("CAG Repeat Length", icon("repeat")),
                       plotlyOutput(ns("cag_plot"), height = "300px")
                     )
              ),
              column(4,
                     wellPanel(
                       h4("Vonsattel Grade", icon("chart-pie")),
                       plotlyOutput(ns("vonsattel_plot"), height = "300px")
                     )
              )
            ),
            fluidRow(
              column(12,
                     wellPanel(
                       h4("HD Severity Scores", icon("thermometer-half")),
                       plotlyOutput(ns("severity_plot"), height = "350px")
                     )
              )
            )
          ),
          conditionalPanel(
            condition = paste0("!output['", ns("has_hd_data"), "']"),
            div(class = "alert alert-info",
                h4("ℹ️ No HD-specific data available"),
                p("Upload a file containing HD patient data to see disease-specific analyses.")
            )
          )
        ),
        
        # stats analysis tab
        tabPanel(
          "Statistical Analysis",
          br(),
          fluidRow(
            column(6,
                   wellPanel(
                     h4("Group Comparisons", icon("balance-scale")),
                     verbatimTextOutput(ns("statistical_tests"))
                   )
            ),
            column(6,
                   wellPanel(
                     h4("HD Correlations", icon("link")),
                     verbatimTextOutput(ns("hd_correlations"))
                   )
            )
          ),
          fluidRow(
            column(12,
                   wellPanel(
                     h4("Sample Characteristic Correlations", icon("project-diagram")),
                     div(style = "text-align: center;",
                         plotlyOutput(ns("correlation_plot"), height = "500px")
                     )
                   )
            )
          )
        ),
        
        # data table tab
        tabPanel(
          "Data Table",
          br(),
          fluidRow(
            column(12,
                   wellPanel(
                     h4("Sample Metadata Table", icon("table")),
                     div(style = "margin-top: 10px;",
                         DT::dataTableOutput(ns("metadata_table"))
                     )
                   )
            )
          )
        )
      )
    ),
    
    # loading message
    conditionalPanel(
      condition = paste0("!output['", ns("file_uploaded"), "']"),
      div(class = "text-center", style = "margin-top: 50px;",
          h4("Please upload the GSE64810 series matrix file or processed CSV to begin"),
          p("The file should contain sample metadata from the GEO dataset.")
      )
    )
  )
}

#' Sample Metadata Server
#' 
#' @param id Module namespace ID
sample_metadata_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    
    # reactive values
    values <- reactiveValues(
      metadata_df = NULL,
      file_uploaded = FALSE
    )
    
    # file upload observer
    observeEvent(input$metadata_file, {
      req(input$metadata_file)
      
      tryCatch({
        file_ext <- tools::file_ext(input$metadata_file$name)
        
        if (file_ext == "csv") {
          values$metadata_df <- read.csv(input$metadata_file$datapath, 
                                         stringsAsFactors = FALSE, 
                                         check.names = FALSE,
                                         na.strings = c("", "NA", "N/A"))
        } else {
          values$metadata_df <- parse_geo_series_matrix(input$metadata_file$datapath)
        }
        
        # basic data processing
        df <- values$metadata_df
        if ("diagnosis" %in% names(df)) {
          df$group <- ifelse(grepl("Huntington", df$diagnosis, ignore.case = TRUE), "HD", "Control")
        }
        if ("mrna_seq_reads" %in% names(df)) {
          df$reads_millions <- df$mrna_seq_reads / 1e6
        }
        
        values$metadata_df <- df
        values$file_uploaded <- TRUE
        
        showNotification(
          paste("✅ Loaded metadata for", nrow(values$metadata_df), "samples"),
          type = "message",
          duration = 3
        )
        
      }, error = function(e) {
        showNotification(
          paste("❌ Error loading file:", e$message),
          type = "error",
          duration = 5
        )
        values$file_uploaded <- FALSE
      })
    })
    
    # reactive expressions
    metadata_clean <- reactive({
      req(values$metadata_df)
      values$metadata_df
    })
    
    # file uploaded status
    output$file_uploaded <- reactive({
      values$file_uploaded
    })
    outputOptions(output, "file_uploaded", suspendWhenHidden = FALSE)
    
    # check if hd data is available
    output$has_hd_data <- reactive({
      req(metadata_clean())
      any(!is.na(metadata_clean()$age_of_onset)) || any(!is.na(metadata_clean()$cag_repeats))
    })
    outputOptions(output, "has_hd_data", suspendWhenHidden = FALSE)
    
    # summary statistics
    output$total_samples <- renderText({
      req(metadata_clean())
      nrow(metadata_clean())
    })
    
    output$hd_samples <- renderText({
      req(metadata_clean())
      if ("group" %in% names(metadata_clean())) {
        sum(metadata_clean()$group == "HD", na.rm = TRUE)
      } else {
        "N/A"
      }
    })
    
    output$control_samples <- renderText({
      req(metadata_clean())
      if ("group" %in% names(metadata_clean())) {
        sum(metadata_clean()$group == "Control", na.rm = TRUE)
      } else {
        "N/A"
      }
    })
    
    output$total_variables <- renderText({
      req(metadata_clean())
      ncol(metadata_clean())
    })
    
    # age distribution plot
    output$age_plot <- renderPlotly({
      req(metadata_clean())
      df <- metadata_clean()
      
      age_col <- intersect(c("age_of_death", "age_death"), names(df))[1]
      if (is.na(age_col) || !age_col %in% names(df) || !"group" %in% names(df)) {
        return(plotly_empty() %>% add_annotations(text = "Age data not available", showarrow = FALSE))
      }
      
      p <- ggplot(df, aes_string(x = age_col, fill = "group")) +
        geom_histogram(alpha = 0.7, bins = 20, position = "identity") +
        scale_fill_manual(values = c("HD" = "#e74c3c", "Control" = "#3498db")) +
        labs(
          title = "Age at Death Distribution",
          x = "Age at Death (years)",
          y = "Count",
          fill = "Group"
        ) +
        theme_minimal() +
        theme(legend.position = "bottom")
      
      ggplotly(p, tooltip = c("x", "y", "fill"))
    })
    
    # RIN plot
    output$rin_plot <- renderPlotly({
      req(metadata_clean())
      df <- metadata_clean()
      
      if (!"rin" %in% names(df) || !"group" %in% names(df)) {
        return(plotly_empty() %>% add_annotations(text = "RIN data not available", showarrow = FALSE))
      }
      
      p <- ggplot(df, aes(x = group, y = rin, fill = group)) +
        geom_boxplot(alpha = 0.7) +
        geom_jitter(width = 0.2, alpha = 0.6) +
        scale_fill_manual(values = c("HD" = "#e74c3c", "Control" = "#3498db")) +
        labs(
          title = "RNA Integrity Number (RIN)",
          x = "Group",
          y = "RIN Score",
          fill = "Group"
        ) +
        theme_minimal() +
        theme(legend.position = "none")
      
      ggplotly(p, tooltip = c("y", "fill"))
    })
    
    # PMI plot
    output$pmi_plot <- renderPlotly({
      req(metadata_clean())
      df <- metadata_clean()
      
      if (!"pmi" %in% names(df) || !"group" %in% names(df)) {
        return(plotly_empty() %>% add_annotations(text = "PMI data not available", showarrow = FALSE))
      }
      
      df <- df %>% filter(!is.na(pmi))
      
      p <- ggplot(df, aes(x = group, y = pmi, fill = group)) +
        geom_boxplot(alpha = 0.7) +
        geom_jitter(width = 0.2, alpha = 0.6) +
        scale_fill_manual(values = c("HD" = "#e74c3c", "Control" = "#3498db")) +
        labs(
          title = "Post-Mortem Interval",
          x = "Group",
          y = "PMI (hours)",
          fill = "Group"
        ) +
        theme_minimal() +
        theme(legend.position = "none")
      
      ggplotly(p, tooltip = c("y", "fill"))
    })
    
    # reads plot
    output$reads_plot <- renderPlotly({
      req(metadata_clean())
      df <- metadata_clean()
      
      reads_col <- intersect(c("reads_millions", "mrna_seq_reads"), names(df))[1]
      if (is.na(reads_col) || !"group" %in% names(df)) {
        return(plotly_empty() %>% add_annotations(text = "Sequencing data not available", showarrow = FALSE))
      }
      
      if (reads_col == "mrna_seq_reads") {
        df$reads_millions <- df$mrna_seq_reads / 1e6
        reads_col <- "reads_millions"
      }
      
      p <- ggplot(df, aes_string(x = "group", y = reads_col, fill = "group")) +
        geom_boxplot(alpha = 0.7) +
        geom_jitter(width = 0.2, alpha = 0.6) +
        scale_fill_manual(values = c("HD" = "#e74c3c", "Control" = "#3498db")) +
        labs(
          title = "Sequencing Depth",
          x = "Group",
          y = "Reads (Millions)",
          fill = "Group"
        ) +
        theme_minimal() +
        theme(legend.position = "none")
      
      ggplotly(p, tooltip = c("y", "fill"))
    })
    
    # HD progression plot
    output$progression_plot <- renderPlotly({
      req(metadata_clean())
      df <- metadata_clean() %>% 
        filter(!is.na(age_of_onset) & !is.na(duration) & !is.na(cag_repeats))
      
      if (nrow(df) == 0) {
        return(plotly_empty() %>% add_annotations(text = "HD progression data not available", showarrow = FALSE))
      }
      
      p <- ggplot(df, aes(x = age_of_onset, y = duration)) +
        geom_point(aes(size = cag_repeats, color = vonsattel_grade), alpha = 0.8) +
        scale_color_viridis_c(name = "Vonsattel\nGrade") +
        scale_size_continuous(name = "CAG\nRepeats", range = c(3, 8)) +
        labs(
          title = "HD Disease Progression",
          x = "Age of Onset (years)",
          y = "Disease Duration (years)"
        ) +
        theme_minimal()
      
      ggplotly(p, tooltip = c("x", "y", "size", "colour"))
    })
    
    # duration plot
    output$duration_plot <- renderPlotly({
      req(metadata_clean())
      df <- metadata_clean() %>% filter(!is.na(duration))
      
      if (nrow(df) == 0) {
        return(plotly_empty() %>% add_annotations(text = "Duration data not available", showarrow = FALSE))
      }
      
      p <- ggplot(df, aes(x = duration)) +
        geom_histogram(fill = "#e74c3c", alpha = 0.7, bins = 10) +
        labs(
          title = "Disease Duration",
          x = "Duration (years)",
          y = "Count"
        ) +
        theme_minimal()
      
      ggplotly(p, tooltip = c("x", "y"))
    })
    
    # CAG plot
    output$cag_plot <- renderPlotly({
      req(metadata_clean())
      df <- metadata_clean() %>% filter(!is.na(cag_repeats))
      
      if (nrow(df) == 0) {
        return(plotly_empty() %>% add_annotations(text = "CAG data not available", showarrow = FALSE))
      }
      
      p <- ggplot(df, aes(x = cag_repeats)) +
        geom_histogram(fill = "#9b59b6", alpha = 0.7, bins = 10) +
        labs(
          title = "CAG Repeats",
          x = "CAG Repeats",
          y = "Count"
        ) +
        theme_minimal()
      
      ggplotly(p, tooltip = c("x", "y"))
    })
    
    # Vonsattel plot
    output$vonsattel_plot <- renderPlotly({
      req(metadata_clean())
      df <- metadata_clean() %>% filter(!is.na(vonsattel_grade))
      
      if (nrow(df) == 0) {
        return(plotly_empty() %>% add_annotations(text = "Vonsattel data not available", showarrow = FALSE))
      }
      
      grade_counts <- df %>% count(vonsattel_grade)
      
      p <- ggplot(grade_counts, aes(x = factor(vonsattel_grade), y = n)) +
        geom_col(fill = "#f39c12", alpha = 0.8) +
        labs(
          title = "Vonsattel Grade",
          x = "Grade",
          y = "Count"
        ) +
        theme_minimal()
      
      ggplotly(p, tooltip = c("x", "y"))
    })
    
    # onset vs CAG plot
    output$onset_cag_plot <- renderPlotly({
      req(metadata_clean())
      df <- metadata_clean() %>% filter(!is.na(age_of_onset) & !is.na(cag_repeats))
      
      if (nrow(df) == 0) {
        return(plotly_empty() %>% add_annotations(text = "Onset/CAG data not available", showarrow = FALSE))
      }
      
      p <- ggplot(df, aes(x = cag_repeats, y = age_of_onset)) +
        geom_point(alpha = 0.7, size = 3, color = "#e74c3c") +
        geom_smooth(method = "lm", se = TRUE, color = "#3498db") +
        labs(
          title = "Age of Onset vs CAG Repeats",
          x = "CAG Repeats",
          y = "Age of Onset (years)"
        ) +
        theme_minimal()
      
      ggplotly(p, tooltip = c("x", "y"))
    })
    
    # severity scores plot
    output$severity_plot <- renderPlotly({
      req(metadata_clean())
      df <- metadata_clean() %>%
        filter(!is.na(hv_striatal_score) | !is.na(hv_cortical_score))
      
      if (nrow(df) == 0) {
        return(plotly_empty() %>% add_annotations(text = "Severity score data not available", showarrow = FALSE))
      }
      
      df_long <- df %>%
        select(sample_id, hv_striatal_score, hv_cortical_score) %>%
        pivot_longer(cols = c(hv_striatal_score, hv_cortical_score),
                     names_to = "score_type", values_to = "score") %>%
        filter(!is.na(score)) %>%
        mutate(score_type = case_when(
          score_type == "hv_striatal_score" ~ "Striatal",
          score_type == "hv_cortical_score" ~ "Cortical"
        ))
      
      p <- ggplot(df_long, aes(x = score_type, y = score, fill = score_type)) +
        geom_boxplot(alpha = 0.7) +
        geom_jitter(width = 0.2, alpha = 0.6) +
        scale_fill_brewer(palette = "Set2") +
        labs(
          title = "HD Severity Scores",
          x = "Brain Region",
          y = "Pathology Score",
          fill = "Region"
        ) +
        theme_minimal() +
        theme(legend.position = "none")
      
      ggplotly(p, tooltip = c("x", "y"))
    })
    
    # statistical tests
    output$statistical_tests <- renderText({
      req(metadata_clean())
      df <- metadata_clean()
      
      if (!"group" %in% names(df)) return("Group variable not available")
      
      results <- c()
      
      # age comparison
      age_col <- intersect(c("age_of_death", "age_death"), names(df))[1]
      if (!is.na(age_col) && age_col %in% names(df)) {
        age_test <- t.test(df[[age_col]] ~ df$group)
        results <- c(results, paste("Age at death t-test p-value:", format.pval(age_test$p.value)))
      }
      
      # RIN comparison
      if ("rin" %in% names(df)) {
        rin_test <- t.test(rin ~ group, data = df)
        results <- c(results, paste("RIN t-test p-value:", format.pval(rin_test$p.value)))
      }
      
      # PMI comparison
      if ("pmi" %in% names(df)) {
        pmi_test <- t.test(pmi ~ group, data = df)
        results <- c(results, paste("PMI t-test p-value:", format.pval(pmi_test$p.value)))
      }
      
      if (length(results) == 0) {
        "No statistical comparisons available"
      } else {
        paste(results, collapse = "\n")
      }
    })
    
    # HD correlations
    output$hd_correlations <- renderText({
      req(metadata_clean())
      df <- metadata_clean()
      
      results <- c()
      
      # CAG vs onset correlation
      if ("cag_repeats" %in% names(df) && "age_of_onset" %in% names(df)) {
        hd_data <- df %>% filter(!is.na(cag_repeats) & !is.na(age_of_onset))
        if (nrow(hd_data) > 2) {
          cor_test <- cor.test(hd_data$cag_repeats, hd_data$age_of_onset)
          results <- c(results, paste("CAG vs Onset correlation:", round(cor_test$estimate, 3)))
          results <- c(results, paste("p-value:", format.pval(cor_test$p.value)))
        }
      }
      
      # duration vs onset correlation
      if ("duration" %in% names(df) && "age_of_onset" %in% names(df)) {
        hd_data <- df %>% filter(!is.na(duration) & !is.na(age_of_onset))
        if (nrow(hd_data) > 2) {
          cor_test <- cor.test(hd_data$duration, hd_data$age_of_onset)
          results <- c(results, paste("Duration vs Onset correlation:", round(cor_test$estimate, 3)))
          results <- c(results, paste("p-value:", format.pval(cor_test$p.value)))
        }
      }
      
      if (length(results) == 0) {
        "No HD-specific correlations available"
      } else {
        paste(results, collapse = "\n")
      }
    })
    
    # correlation plot
    output$correlation_plot <- renderPlotly({
      req(metadata_clean())
      df <- metadata_clean()
      
      # select numeric columns for correlation
      numeric_cols <- df %>%
        select_if(is.numeric) %>%
        select(-matches("reads$")) %>%  # remove raw reads if reads_millions exists
        na.omit()
      
      if (ncol(numeric_cols) < 2) {
        return(plotly_empty() %>% add_annotations(text = "Insufficient numeric data for correlation", showarrow = FALSE))
      }
      
      cor_matrix <- cor(numeric_cols, use = "complete.obs")
      
      # convert to long format for plotting
      cor_long <- expand.grid(Var1 = rownames(cor_matrix), Var2 = colnames(cor_matrix)) %>%
        mutate(correlation = as.vector(cor_matrix))
      
      p <- ggplot(cor_long, aes(x = Var1, y = Var2, fill = correlation)) +
        geom_tile() +
        scale_fill_gradient2(low = "#3498db", mid = "white", high = "#e74c3c",
                             midpoint = 0, limits = c(-1, 1)) +
        labs(
          title = "Sample Characteristic Correlations",
          x = "", y = "", fill = "Correlation"
        ) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      
      ggplotly(p, tooltip = c("x", "y", "fill"))
    })
    
    # data table
    output$metadata_table <- DT::renderDataTable({
      req(metadata_clean())
      
      DT::datatable(
        metadata_clean(),
        options = list(
          pageLength = 25,
          scrollX = TRUE,
          dom = 'Bfrtip',
          buttons = c('copy', 'csv', 'excel')
        ),
        extensions = 'Buttons',
        rownames = FALSE
      ) %>%
        DT::formatRound(columns = intersect(names(metadata_clean()), 
                                            c("pmi", "rin", "reads_millions", "age_of_death", 
                                              "age_of_onset", "duration", "cag_repeats", 
                                              "vonsattel_grade", "hv_striatal_score", "hv_cortical_score")), 
                        digits = 2)
    })
    
  })
}
