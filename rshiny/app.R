# app.R - Multi-Module Genomics Explorer - Designed for GSE64810

# increase file upload size limit
options(shiny.maxRequestSize = 500*1024^2)

# load required libraries
library(shiny)
library(plotly)
library(ggplot2)
library(DT)
library(dplyr)
library(readr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(tools)

# source module functions
source("R/fct_data_processing.R")
source("R/mod_volcano.R")
source("R/mod_expression_analysis.R")
source("R/mod_sample_metadata.R")
source("R/mod_gsea_analysis.R")

# ui with tabs
ui <- fluidPage(
  titlePanel("RNASeq Genomics Explorer"),
  
  # css for styling
  tags$head(
    tags$style(HTML("
      .nav-tabs > li > a {
        background-color: #f8f9fa;
        border-color: #dee2e6;
      }
      .nav-tabs > li.active > a {
        background-color: #007bff;
        color: white;
        border-color: #007bff;
      }
      .tab-content {
        border: 1px solid #dee2e6;
        border-top: none;
        padding: 20px;
        background-color: #ffffff;
      }
      h4 {
        margin-top: 30px;
        margin-bottom: 15px;
      }
    "))
  ),
  
  # navigtaion tabs
  tabsetPanel(
    id = "main_tabs",
    
    # Sample Metadata Tab
    tabPanel(
      "Sample Metadata",
      sample_metadata_ui("sample_metadata")
    ),
    
    # Volcano Plot Tab
    tabPanel(
      title = "Volcano Plot",
      value = "volcano_tab",
      volcano_plot_ui("volcano")
    ),
    
    # expression analysis tab (diagnostic + pca + umap)
    tabPanel(
      title = "Expression Analysis",
      value = "expression_tab",
      expression_analysis_ui("expression")
    ),
    
    # gsea tab
    tabPanel(
      "GSEA Analysis",
      gsea_analysis_ui("gsea_module")
    ),
    
    # about tab
    tabPanel(
      title = "About",
      value = "about_tab",
      fluidRow(
        column(8, offset = 2,
               wellPanel(
                 h3("GSE64810 Huntington's Genomics Explorer"),
                 p("This interactive dashboard visualizes RNA-seq data from GSE64810, comparing gene expression 
              between Huntington's disease patients and healthy controls in prefrontal cortex brain tissue."),
                 
                 h4("Dataset Information:"),
                 tags$ul(
                   tags$li("Study: Huntington's disease vs. control brain samples"),
                   tags$li("Tissue: Prefrontal cortex (Brodmann area 9)"),
                   tags$li("Samples: 69 total (20 HD patients + 49 controls)"),
                   tags$li("Platform: Illumina RNA-sequencing"),
                   tags$li("Genes: ~28,000 in final dataset")
                 ),
                 h4("Required Files:"),
                 tags$ul(
                   tags$li(code("GSE64810_series_matrix.txt"), " - Sample metadata and clinical variables"),
                   tags$li(code("GSE64810_mlhd_DESeq2_diffexp_DESeq2_outlier_trimmed_adjust.txt"), " - Differential expression results"),
                   tags$li(code("GSE64810_mlhd_DESeq2_norm_counts_adjust.txt"), " - Normalized expression counts"),
                   tags$li(code("pathway_gene_sets.gmt"), " - Gene set collections (e.g., from MSigDB for GSEA)")
                 ),
                 
                 h4("Analysis Tabs:"),
                 div(style = "background-color: #f8f9fa; padding: 15px; border-radius: 5px;",
                     fluidRow(
                       column(2,
                              h6("Metadata Analysis"),
                              tags$ul(style = "font-size: 12px;",
                                      tags$li("Clinical characteristics"),
                                      tags$li("HD progression metrics"),
                                      tags$li("Quality control metrics")
                              )
                       ),
                       column(2,
                              h6("Quality Diagnostics"),
                              tags$ul(style = "font-size: 12px;",
                                      tags$li("Library size distribution"),
                                      tags$li("Mean-variance relationships"),
                                      tags$li("Sample correlation matrix")
                              )
                       ),
                       column(2,
                              h6("Differential Expression"),
                              tags$ul(style = "font-size: 12px;",
                                      tags$li("Interactive PCA plots"),
                                      tags$li("Variance explained analysis"),
                                      tags$li("Gene contribution rankings")
                              )
                       ),
                       column(2,
                              h6("Heatmap Analysis"),
                              tags$ul(style = "font-size: 12px;",
                                      tags$li("Hierarchical clustering"),
                                      tags$li("Z-score normalization"),
                                      tags$li("Variable gene selection")
                              )
                       ),
                       column(2,
                              h6("GSEA Analysis"),
                              tags$ul(style = "font-size: 12px;",
                                      tags$li("Pathway enrichment"),
                                      tags$li("Multiple ranking methods"),
                                      tags$li("Interactive visualizations")
                              )
                       ),
                       column(2,
                              h6("Data Export"),
                              tags$ul(style = "font-size: 12px;",
                                      tags$li("Results tables"),
                                      tags$li("Pathway gene lists"),
                                      tags$li("Ranked gene files")
                              )
                       )
                     )
                 ),
                 
                 hr(),
                 p(em("Built with RShiny, ComplexHeatmap, fgsea, and other modern genomics visualization tools. Want a 
                      dashboard for your data? Reach out here!"),
                   style = "text-align: center; color: #6c757d;")
               )
        )
      )
    )
  )
)

# server
server <- function(input, output, session) {
  
  # BP LATER - add modules for sharing data between modules
  
  # calling modules
  sample_metadata_server("sample_metadata")
  volcano_plot_server("volcano")
  expression_analysis_server("expression")
  gsea_analysis_server("gsea_module")
}

# run the app!
shinyApp(ui = ui, server = server)
