#
# Shiny App for Multi-Omics Biomarker Discovery Pipeline
#

# 1. --- SETUP: LOAD LIBRARIES AND HELPERS ---
# -------------------------------------------------

suppressPackageStartupMessages({
  library(shiny)
  library(shinyjs)
})

load_app_packages <- function(pkgs) {
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop(
      "Missing required packages for the app: ",
      paste(missing_pkgs, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(lapply(pkgs, function(p) {
    suppressPackageStartupMessages(library(p, character.only = TRUE))
  }))
}

# All required packages for the app
required_pkgs <- c(
  "WGCNA", "dplyr", "readr", "clusterProfiler", "org.Hs.eg.db", "enrichplot",
  "ReactomePA", "ggplot2", "openxlsx", "fgsea", "msigdbr", "tidyr",
  "tibble", "pheatmap", "stringr", "caret", "glmnet", "pROC", "DT", "sva", "imputeLCMD", "ggrepel", "impute"
)

message("Loading required packages...")
load_app_packages(required_pkgs)
message("All packages loaded.")

# Source the helper file containing the analysis functions
# (Functions expected: run_module_scoring, run_ml_analysis, run_gsea_analysis, plot_gsea_results)
if (file.exists("helpers.R")) { 
  source("helpers.R")
} else {
  warning("helpers.R not found in working directory. Some features may not work.")
}

# Increase maximum file upload size to 500MB (adjust as needed)
options(shiny.maxRequestSize = 500 * 1024^2)

# 2. --- UI DEFINITION ---
# --------------------------

ui <- tagList(
  tags$head(
    tags$style(HTML("
      .navbar { margin-bottom: 16px; }
      .well { border-radius: 8px; }
      .tab-pane { padding-top: 12px; }
      .pipeline-footer { margin: 20px 0 12px; text-align: center; font-size: 12px; color: #6b7280; }
    "))
  ),
  useShinyjs(),
  navbarPage(
             title = "ADNI-ADRC Proteomics Pipeline",
             # --- Tab 1: Data Preprocessing ---
             tabPanel("1. Preprocessing",
                      sidebarLayout(
                        sidebarPanel(
                          h4("Input Data"),
                          p("Upload raw proteomics and phenotype data (Excel/CSV)."),
                          fileInput("raw_proto_file", "Proteomics Data (Expression)", accept = c(".xlsx", ".csv")),
                          fileInput("raw_pheno_file", "Phenotype Data", accept = c(".xlsx", ".csv")),
                          selectInput("impute_method", "Imputation Method", choices = c("Hybrid", "KNN", "Minimum", "None"), selected = "Hybrid"),
                          actionButton("run_preprocess_btn", "Run Preprocessing", class = "btn-primary"),
                          hr(),
                          wellPanel(
                            h4("About This Step"),
                            p("This step handles missing value imputation (using QRILC or MinProb) and filters out low-quality samples/features."),
                            p("It ensures the expression data is numeric and aligned with the phenotype metadata.")
                          )
                        ),
                        mainPanel(
                          h3("Preprocessing Diagnostics"),
                          plotOutput("missingness_plot", height = "400px"),
                          downloadButton("dl_missing_plot", "Download Plot"),
                          hr(),
                          h4("Imputation Effect"),
                          plotOutput("impute_dist_plot", height = "400px"),
                          downloadButton("dl_impute_plot", "Download Plot"),
                          hr(),
                          h4("Processed Data Preview"),
                          DT::dataTableOutput("preview_processed_data")
                        )
                      )
             ),
             
             # --- Tab 2: Harmonization ---
             tabPanel("2. Harmonization",
                      sidebarLayout(
                        sidebarPanel(
                          h4("Batch Correction"),
                          p("Use ComBat to remove batch effects (e.g., Cohort)."),
                          selectInput("batch_col", "Select Batch Column", choices = NULL),
                          actionButton("run_combat_btn", "Run ComBat", class = "btn-primary"),
                          hr(),
                          downloadButton("dl_harmonized_csv", "Save Harmonized Data (.csv)"),
                          downloadButton("dl_harmonized_rdata", "Save Harmonized Data (02_output.RData)"),
                          hr(),
                          wellPanel(
                            h4("About This Step"),
                            p("Harmonization merges datasets from different sources (e.g., ADNI and ADRC) by removing systematic technical variations (batch effects) while preserving biological signal.")
                          )
                        ),
                        mainPanel(
                          h3("Harmonization Results"),
                          fluidRow(
                            column(6, h4("Before Correction"), plotOutput("pca_before_combat")),
                            column(6, h4("After Correction"), plotOutput("pca_after_combat"))
                          ),
                          downloadButton("dl_combat_plots", "Download Comparison Plots")
                        )
                      )
             ),
             
             # --- Tab 3: PCA ---
             tabPanel("3. PCA",
                      sidebarLayout(
                        sidebarPanel(
                          h4("Trait PCA"),
                          p("Decompose phenotype variables into Principal Components (e.g., WMH traits)."),
                          selectInput("pca_vars", "Select Traits for PCA", choices = NULL, multiple = TRUE),
                          selectInput("pca_color_col", "Color Samples by (Optional)", choices = NULL),
                          actionButton("run_pca_btn", "Run PCA", class = "btn-primary"),
                          hr(),
                          wellPanel(
                            h4("About This Step"),
                            p("This step performs PCA on selected clinical traits (e.g., white matter hyperintensity volumes) to reduce dimensionality and create composite scores (PC1, PC2, etc.) for downstream analysis.")
                          )
                        ),
                        mainPanel(
                          h3("PCA Outputs"),
                          tabsetPanel(
                            tabPanel("Scree Plot", plotOutput("pca_scree_plot", height = "400px")),
                            tabPanel("Biplot", plotOutput("pca_biplot", height = "600px")),
                            tabPanel("Loadings", plotOutput("pca_loadings_plot", height = "600px"))
                          ),
                          hr(),
                          downloadButton("dl_pca_scores", "Download PCA Scores (CSV)"),
                          downloadButton("dl_pca_plots", "Download Plots (PDF)")
                        )
                      )
             ),
             
             # --- Tab 4: Setup & Load Data (WGCNA) ---
             tabPanel("4. WGCNA Setup",
                      sidebarLayout(
                        sidebarPanel(
                          h4("Load WGCNA Data"),
                          p("Select the core output file from the WGCNA step (e.g., WGCNA2_core_outputs.RData)."),
                          fileInput("wgcna_file_input", "Choose .RData File", accept = ".RData"),
                          fileInput("pca_file_input", "Optional: Choose PCA Scores CSV", accept = ".csv"),
                          actionButton("load_data_btn", "Load Data", class = "btn-primary"),
                          hr(),
                          wellPanel(
                              h4("About This Step"),
                              p("This initial step is for loading the necessary data to run the pipeline. The primary input is an `.RData` file generated from a prior WGCNA run (for example from `03_WGCNA_running.Rmd`). This file should contain essential objects like the expression matrix (`expr`), trait data (`traits`), and pre-calculated module assignments (`moduleColors`)."),
                              p("Loading a pre-computed file saves significant time by skipping the computationally intensive network construction steps.")
                          )
                        ),
                        mainPanel(
                          h3("Data Loading Status"),
                          verbatimTextOutput("data_status_output")
                        )
                      )
             ),
             
             # --- Tab 5: WGCNA (Re-run or Parameter Check) ---
             tabPanel("5. WGCNA",
                      sidebarLayout(
                        sidebarPanel(
                          h4("1. Soft Threshold Analysis"),
                          p("Analyze network topology to select an appropriate soft-thresholding power."),
                          actionButton("run_sft_btn", "Analyze Soft Threshold"),
                          hr(),
                          h4("2. Network Construction"),
                          p("Run the main WGCNA pipeline with the chosen parameters."),
                          numericInput("wgcna_softPower", "Soft Threshold Power", value = 6, min = 1, max = 20),
                          numericInput("wgcna_minClusterSize", "Min. Cluster Size", value = 30, min = 10),
                          sliderInput("wgcna_cutHeight", "Module Merge Cut Height", value = 0.30, min = 0, max = 1),
                          selectInput("wgcna_traits", "Select Traits for Heatmap", choices = NULL, multiple = TRUE),
                          actionButton("run_wgcna_btn", "Run WGCNA Step", class = "btn-primary"),
                          hr(),
                          downloadButton("dl_wgcna_results", "Save WGCNA Results (.RData)"),
                          hr(),
                          wellPanel(
                              h4("About This Step"),
                              p("Weighted Gene Co-expression Network Analysis (WGCNA) is a systems biology method for describing correlation patterns among genes or proteins across multiple samples. It is used to find clusters (modules) of highly correlated molecules."),
                              h5("Key Concepts:"),
                              tags$ul(
                                  tags$li(strong("Soft Thresholding:"), " A power is chosen to amplify strong correlations and penalize weak ones, ensuring the network has 'scale-free' properties, which is common in biological networks."),
                                  tags$li(strong("Topological Overlap Matrix (TOM):"), " A robust measure of network interconnectedness that considers both direct and indirect connections between molecules."),
                                  tags$li(strong("Module Eigengene (ME):"), " The first principal component of a module, representing the overall expression profile of all molecules within that module.")
                              )
                          )
                        ),
                        mainPanel(
                          h3("WGCNA Outputs"),
                          div(style = "height: 800px; overflow-y: scroll; margin-top: 20px;",
                              
                              h4("A. Sample Clustering (QC)"),
                              plotOutput("wgcna_sample_tree", height = "400px"),
                              downloadButton("download_sample_tree", "Download Plot"),
                              br(), hr(),
                              
                              h4("B. Scale-Free Topology"),
                              plotOutput("wgcna_sft_fit_plot", height = "400px"),
                              downloadButton("download_sft_fit", "Download Plot"),
                              br(), hr(),
                              
                              h4("B. Mean Connectivity"),
                              plotOutput("wgcna_sft_mean_plot", height = "400px"),
                              downloadButton("download_sft_mean", "Download Plot"),
                              br(), hr(),
                              
                              h4("C. Gene Dendrogram and Module Colors"),
                              plotOutput("wgcna_dendro_plot", height = "600px"),
                              downloadButton("download_dendro", "Download Dendrogram"),
                              br(), hr(),
                              
                              h4("D. Eigengene Clustering"),
                              plotOutput("wgcna_eigengene_plot", height = "400px"),
                              downloadButton("download_eigengene", "Download Eigengene Plot"),
                              br(), hr(),
                              h4("E. Module-Trait Relationships"),
                              plotOutput("wgcna_trait_heatmap", height = "800px"),
                              downloadButton("download_heatmap", "Download Heatmap")
                          )
                        )
                      ),
                      hr(),
                      wellPanel(
                          h5("Citations:"),
                          tags$ul(
                              tags$li("Langfelder, P., & Horvath, S. (2008). WGCNA: an R package for weighted correlation network analysis. ", tags$i("BMC Bioinformatics, 9"), ", 559. ", tags$a(href="https://doi.org/10.1186/1471-2105-9-559", target="_blank", "DOI: 10.1186/1471-2105-9-559")),
                              tags$li("Zhang, B., & Horvath, S. (2005). A general framework for weighted gene co-expression network analysis. ", tags$i("Statistical Applications in Genetics and Molecular Biology, 4"), "(1). ", tags$a(href="https://doi.org/10.2202/1544-6115.1128", target="_blank", "DOI: 10.2202/1544-6115.1128"))
                          )
                      )
             ),
             
             # --- Tab 6: Module Scoring ---
             tabPanel("6. Module Scoring",
                      sidebarLayout(
                        sidebarPanel(
                          h4("Scoring Parameters"),
                          p("These parameters follow the canonical module-scoring notebook `05_Module_Scoring_Selection.Rmd`."),
                          sliderInput("scoring_r_thresh", "Correlation Threshold (r)", value = 0.12, min = 0, max = 1),
                          sliderInput("scoring_q_thresh", "Significance Threshold (q-value)", value = 0.05, min = 0.01, max = 0.2),
                          actionButton("run_scoring_btn", "Run Scoring", class = "btn-primary"),
                          hr(),
                          wellPanel(
                              h4("About This Step"),
                              p("This step quantifies the relevance of each module to different groups of clinical traits. A 'hit' is counted if a module's eigengene has a statistically significant correlation (below the chosen q-value threshold) with a trait, and the correlation strength is above the chosen 'r' threshold."),
                              p("This allows for a high-level biological categorization of modules, such as identifying modules primarily associated with 'vascular', 'CAA', or 'AD' traits.")
                          )
                        ),
                        mainPanel(
                          h3("Module Scoring Results"),
                          DT::dataTableOutput("scoring_table_output")
                        )
                      ),
                      hr(),
                      wellPanel(
                          h5("Citation:"),
                          p("This scoring method is a custom implementation based on the principles of module-trait relationship analysis described in the core WGCNA papers.")
                      )
             ),
             
             # --- Tab 7: Machine Learning ---
             tabPanel("7. Machine Learning",
                      sidebarLayout(
                        sidebarPanel(
                          h4("ML Parameters"),
                          p("These parameters follow the canonical machine-learning notebook `06_ADNI_ADRC_ML.Rmd`."),
                          selectInput("ml_target_col", "Target Variable", choices = c("stroke", "AD_status", "caa")),
                          sliderInput("ml_split_ratio", "Train/Test Split Ratio", value = 0.7, min = 0.5, max = 0.9),
                          sliderInput("ml_alpha", "Elastic Net Alpha", value = 0.5, min = 0, max = 1),
                          actionButton("run_ml_btn", "Run ML Model", class = "btn-primary"),
                          hr(),
                          wellPanel(
                              h4("About This Step"),
                              p("This step uses the module eigengenes (MEs) as features to predict a clinical outcome (e.g., 'stroke'). Using MEs instead of thousands of individual proteins drastically reduces the number of features, mitigates multiple testing issues, and incorporates biological structure into the model."),
                              p("An Elastic Net (a combination of LASSO and Ridge regression) is used for its ability to handle correlated predictors and perform automatic feature selection, identifying the most predictive modules.")
                          )
                        ),
                        mainPanel(
                          div(style = "height: 800px; overflow-y: scroll; margin-top: 20px;",
                              h3("Model Performance"),
                              verbatimTextOutput("ml_cm_output"),
                              plotOutput("ml_roc_plot"),
                              h3("Top Predictors"),
                              DT::dataTableOutput("ml_coeffs_output")
                          )
                        )
                      ),
                      hr(),
                      wellPanel(
                          h5("Citations:"),
                          tags$ul(
                              tags$li("Zou, H., & Hastie, T. (2005). Regularization and variable selection via the elastic net. ", tags$i("Journal of the Royal Statistical Society: Series B, 67"), "(2), 301-320. ", tags$a(href="https://doi.org/10.1111/j.1467-9868.2005.00503.x", target="_blank", "DOI: 10.1111/j.1467-9868.2005.00503.x"))
                          )
                      )
             ),
             
             # --- Tab 8: Functional Analysis (GSEA) ---
             tabPanel("8. GSEA",
                      sidebarLayout(
                        sidebarPanel(
                          h4("GSEA Parameters"),
                          p("These parameters follow the canonical GSEA notebook `07a_Functional_Analysis_GSEA.Rmd`."),
                          selectInput("gsea_module", "Select Module to Analyze", choices = NULL),
                          numericInput("gsea_minSize", "Min. Gene Set Size", value = 15, min = 5),
                          numericInput("gsea_maxSize", "Max. Gene Set Size", value = 500, min = 100),
                          sliderInput("gsea_padj_cutoff", "Significance Cutoff (padj)", value = 0.05, min = 0.01, max = 0.25),
                          # Optional: allow selecting analyte info instead of hardcoded path
                          fileInput("analyte_info_file", "Analyte Info (TSV)", accept = c(".tsv", ".txt")),
                          actionButton("run_gsea_btn", "Run GSEA", class = "btn-primary"),
                          hr(),
                          wellPanel(
                              h4("About This Step"),
                              p("Gene Set Enrichment Analysis (GSEA) is used to identify biological pathways or functions that are significantly enriched within a module. Unlike simple Over-Representation Analysis (ORA), GSEA considers the entire ranked list of proteins."),
                              p("In this pipeline, proteins are ranked by their Module Membership (kME) value for a selected module. GSEA then determines if known biological pathways (from databases like GO, Reactome, or Hallmark) are enriched at the top of this ranked list.")
                          )
                        ),
                        mainPanel(
                          div(style = "height: 800px; overflow-y: scroll; margin-top: 20px;",
                              h3("GSEA Results"),
                              plotOutput("gsea_dotplot_output"),
                              DT::dataTableOutput("gsea_table_output")
                          )
                        )
                      ),
                      hr(),
                      wellPanel(
                          h5("Citations:"),
                          tags$ul(
                              tags$li("Subramanian, A., et al. (2005). Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles. ", tags$i("PNAS, 102"), "(43), 15545-15550. ", tags$a(href="https://doi.org/10.1073/pnas.0506580102", target="_blank", "DOI: 10.1073/pnas.0506580102")),
                              tags$li("Sergushichev, A. A. (2016). An algorithm for fast preranked gene set enrichment analysis using cumulative statistic calculation. ", tags$i("bioRxiv"), ". ", tags$a(href="https://doi.org/10.1101/060012", target="_blank", "DOI: 10.1101/060012"))
                          )
                      )
             ),
             
             # --- Tab 9: Functional Analysis (ORA) ---
             tabPanel("9. ORA",
                      sidebarLayout(
                        sidebarPanel(
                          h4("ORA Parameters"),
                          p("These parameters follow the canonical ORA notebook `07b_Funtional_Analysis_ORA.Rmd`."),
                          selectInput("ora_module", "Select Module to Analyze", choices = NULL),
                          sliderInput("ora_padj_cutoff", "Significance Cutoff (padj)", value = 0.05, min = 0.01, max = 0.25),
                          actionButton("run_ora_btn", "Run ORA", class = "btn-primary"),
                          hr(),
                          wellPanel(
                              h4("About This Step"),
                              p("Over-Representation Analysis (ORA) determines whether known biological functions or pathways are over-represented in a given list of genes/proteins. Unlike GSEA, ORA requires a discrete list of interesting molecules (e.g., all proteins within a module) and compares it against a background universe."),
                              p("This analysis tests for enrichment in Gene Ontology (GO), KEGG, and Reactome pathway databases.")
                          )
                        ),
                        mainPanel(
                          div(style = "height: 800px; overflow-y: scroll; margin-top: 20px;",
                              h3("ORA Results"),
                              plotOutput("ora_dotplot_output"),
                              DT::dataTableOutput("ora_table_output")
                          )
                        )
                      ),
                      hr(),
                      wellPanel(
                          h5("Citations:"),
                          tags$ul(
                              tags$li("Ashburner, M., et al. (2000). Gene Ontology: tool for the unification of biology. ", tags$i("Nature Genetics, 25"), "(1), 25-29. ", tags$a(href="https://doi.org/10.1038/75556", target="_blank", "DOI: 10.1038/75556")),
                              tags$li("Yu, G., et al. (2012). clusterProfiler: an R package for comparing biological themes among gene clusters. ", tags$i("OMICS: A Journal of Integrative Biology, 16"), "(5), 284-287. ", tags$a(href="https://doi.org/10.1089/omi.2011.0118", target="_blank", "DOI: 10.1089/omi.2011.0118"))
                          )
                      )
             )
  ),
  div(class = "pipeline-footer", "ADNI-ADRC proteomics analysis app")
)


# 3. --- SERVER LOGIC ---
# --------------------------

server <- function(input, output, session) {
  
  # Reactive values to store data and results across steps
  pipeline_data <- reactiveValues(
    loaded_data = NULL,
    raw_proto = NULL,
    raw_pheno = NULL,
    processed_data = NULL,
    sft_results = NULL,
    wgcna_filtered_data = NULL, # Store QC'd data
    pca_results = NULL,
    pca_samples = NULL,
    wgcna_results = NULL,
    trait_cor_results = NULL,
    ora_results = NULL,
    status = "App started. All packages loaded. Please load data."
  )
  
  # Helper to check if expected objects exist in loaded RData
  has_objects <- function(obj_list, needed) {
    all(needed %in% names(obj_list))
  }
  
  # --- Server Logic for Tab 1: Preprocessing ---
  
  observeEvent(input$run_preprocess_btn, {
    req(input$raw_proto_file, input$raw_pheno_file)
    
    tryCatch({
      # Read files
      ext_proto <- tools::file_ext(input$raw_proto_file$name)
      if(ext_proto == "xlsx") {
        proto <- openxlsx::read.xlsx(input$raw_proto_file$datapath)
      } else {
        proto <- read.csv(input$raw_proto_file$datapath)
      }
      
      ext_pheno <- tools::file_ext(input$raw_pheno_file$name)
      if(ext_pheno == "xlsx") {
        pheno <- openxlsx::read.xlsx(input$raw_pheno_file$datapath)
      } else {
        pheno <- read.csv(input$raw_pheno_file$datapath)
      }
      
      # Run preprocessing
      res <- run_preprocessing(proto, pheno, input$impute_method)
      pipeline_data$processed_data <- res
      
      # Clear downstream results to prevent mismatch
      pipeline_data$harmonized_expr <- NULL
      pipeline_data$pca_results <- NULL
      
      # Update UI choices for next steps
      batch_choices <- colnames(res$pheno)
      default_batch <- if("Cohort" %in% batch_choices) "Cohort" else batch_choices[1]
      updateSelectInput(session, "batch_col", choices = batch_choices, selected = default_batch)
      
      showNotification("Preprocessing complete.", type = "message")
      
    }, error = function(e) {
      showNotification(paste("Preprocessing failed:", conditionMessage(e)), type = "error")
    })
  })
  
  # --- Reset Harmonization when Batch Column Changes ---
  observeEvent(input$batch_col, {
    pipeline_data$harmonized_expr <- NULL
    pipeline_data$pca_results <- NULL
  })
  
  output$missingness_plot <- renderPlot({
    req(input$raw_proto_file) # Just to trigger on upload, ideally use reactive
    # For demo, we plot missingness of the *uploaded* data before full processing if possible,
    # or just plot the processed data's stats. Let's plot the processed data for now.
    req(pipeline_data$processed_data)
    plot_missingness(pipeline_data$processed_data$expr)
  })
  
  output$dl_missing_plot <- downloadHandler(
    filename = "missingness_plot.pdf",
    content = function(file) {
      req(pipeline_data$processed_data)
      pdf(file, width = 8, height = 6)
      print(plot_missingness(pipeline_data$processed_data$expr))
      dev.off()
    }
  )
  
  output$impute_dist_plot <- renderPlot({
    req(pipeline_data$processed_data)
    # Check if we have raw data to compare
    if (!is.null(pipeline_data$processed_data$expr_raw)) {
      plot_imputation_distribution(pipeline_data$processed_data$expr_raw, pipeline_data$processed_data$expr)
    }
  })
  
  output$dl_impute_plot <- downloadHandler(
    filename = "imputation_distribution.pdf",
    content = function(file) {
      req(pipeline_data$processed_data$expr_raw)
      pdf(file, width = 8, height = 6)
      print(plot_imputation_distribution(pipeline_data$processed_data$expr_raw, pipeline_data$processed_data$expr))
      dev.off()
    }
  )
  
  output$preview_processed_data <- DT::renderDataTable({
    req(pipeline_data$processed_data)
    DT::datatable(head(pipeline_data$processed_data$expr, 20), options = list(scrollX = TRUE))
  })
  
  # --- Server Logic for Tab 2: Harmonization ---
  
  observeEvent(input$run_combat_btn, {
    req(pipeline_data$processed_data, input$batch_col)
    
    tryCatch({
      harmonized <- run_harmonization(
        expr_data = pipeline_data$processed_data$expr,
        pheno_data = pipeline_data$processed_data$pheno,
        batch_col = input$batch_col
      )
      pipeline_data$harmonized_expr <- harmonized
      showNotification("Harmonization (ComBat) complete.", type = "message")
    }, error = function(e) {
      showNotification(paste("Harmonization failed:", conditionMessage(e)), type = "error")
    })
  })
  
  output$pca_before_combat <- renderPlot({
    req(pipeline_data$processed_data, input$batch_col)
    plot_pca_batch(pipeline_data$processed_data$expr, pipeline_data$processed_data$pheno, input$batch_col, "Before ComBat")
  })
  
  output$pca_after_combat <- renderPlot({
    req(pipeline_data$harmonized_expr, pipeline_data$processed_data, input$batch_col)
    # Re-construct a data frame for plotting if needed, or just pass matrix
    plot_pca_batch(as.data.frame(pipeline_data$harmonized_expr), pipeline_data$processed_data$pheno, input$batch_col, "After ComBat")
  })
  
  output$dl_combat_plots <- downloadHandler(
    filename = "harmonization_plots.pdf",
    content = function(file) {
      req(pipeline_data$harmonized_expr)
      pdf(file, width = 12, height = 6)
      print(plot_pca_batch(pipeline_data$processed_data$expr, pipeline_data$processed_data$pheno, input$batch_col, "Before"))
      print(plot_pca_batch(as.data.frame(pipeline_data$harmonized_expr), pipeline_data$processed_data$pheno, input$batch_col, "After"))
      dev.off()
    }
  )
  
  output$dl_harmonized_csv <- downloadHandler(
    filename = function() {
      paste0("Harmonized_Expression_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(pipeline_data$harmonized_expr)
      write.csv(as.data.frame(pipeline_data$harmonized_expr), file, row.names = TRUE)
    }
  )
  
  output$dl_harmonized_rdata <- downloadHandler(
    filename = function() { "02_output.RData" },
    content = function(file) {
      req(pipeline_data$harmonized_expr, pipeline_data$processed_data)
      # Save with variable names expected by 03_WGCNA_running.Rmd
      Proto_combat_z <- as.data.frame(pipeline_data$harmonized_expr)
      # Ensure pheno matches harmonized data (in case samples were dropped/reordered)
      pheno_combat <- pipeline_data$processed_data$pheno[rownames(Proto_combat_z), , drop=FALSE]
      save(Proto_combat_z, pheno_combat, file = file)
    }
  )
  
  # --- Server Logic for Tab 3: PCA ---
  
  # Helper to get current phenotype data
  current_pheno <- reactive({
    if (!is.null(pipeline_data$processed_data$pheno)) {
      return(pipeline_data$processed_data$pheno)
    } else if (!is.null(pipeline_data$loaded_data$traits)) {
      return(pipeline_data$loaded_data$traits)
    }
    return(NULL)
  })
  
  # Update PCA inputs when data changes
  observeEvent(current_pheno(), {
    pheno <- current_pheno()
    req(pheno)
    nums <- names(pheno)[sapply(pheno, is.numeric)]
    
    # Default selection: try to match WMH vars if present, else first few
    wmh_defaults <- c("wmh_juxta", "wmh_frontal", "wmh_pv", "wmh_parietal", "wmh_posterior")
    selected <- intersect(wmh_defaults, nums)
    if(length(selected) == 0) selected <- nums[1:min(5, length(nums))]
    
    updateSelectInput(session, "pca_vars", choices = nums, selected = selected)
    updateSelectInput(session, "pca_color_col", choices = c("None", names(pheno)))
    
    # Also update WGCNA trait selector if we have data but haven't run WGCNA yet
    # This allows users to see what traits are available before running
    updateSelectInput(session, "wgcna_traits", choices = nums, selected = selected)
  })
  
  observeEvent(input$run_pca_btn, {
    req(current_pheno(), input$pca_vars)
    pheno_sub <- na.omit(current_pheno()[, input$pca_vars, drop = FALSE])
    pca_res <- run_pca_analysis(pheno_sub)
    pipeline_data$pca_results <- pca_res
    pipeline_data$pca_samples <- rownames(pheno_sub)
  })
  
  output$pca_scree_plot <- renderPlot({
    req(pipeline_data$pca_results)
    plot_pca_scree(pipeline_data$pca_results)
  })
  
  output$pca_biplot <- renderPlot({
    req(pipeline_data$pca_results)
    plot_pca_biplot(pipeline_data$pca_results)
  })
  
  output$pca_loadings_plot <- renderPlot({
    req(pipeline_data$pca_results)
    grid::grid.draw(plot_pca_loadings_heatmap(pipeline_data$pca_results))
  })
  
  output$dl_pca_plots <- downloadHandler(
    filename = "pca_plots.pdf",
    content = function(file) {
      req(pipeline_data$pca_results)
      pdf(file, width = 10, height = 8)
      
      # Scree
      print(plot_pca_scree(pipeline_data$pca_results))
      
      # Biplot
      print(plot_pca_biplot(pipeline_data$pca_results))
      
      # Loadings
      grid::grid.newpage()
      grid::grid.draw(plot_pca_loadings_heatmap(pipeline_data$pca_results))
      
      dev.off()
    }
  )

  output$dl_pca_scores <- downloadHandler(
    filename = function() { "PCA_scores.csv" },
    content = function(file) {
      req(pipeline_data$pca_results)
      scores <- as.data.frame(pipeline_data$pca_results$x[, 1:min(5, ncol(pipeline_data$pca_results$x))])
      # Add SampleID column
      scores$SampleID <- pipeline_data$pca_samples
      scores <- scores %>% dplyr::select(SampleID, everything())
      write.csv(scores, file, row.names = FALSE)
    }
  )
  
  # --- Server Logic for Tab 4: WGCNA Setup ---
  
  observeEvent(input$load_data_btn, {
    req(input$wgcna_file_input)
    
    # Create a new environment to load the RData into
    load_env <- new.env()
    tryCatch({
      load(input$wgcna_file_input$datapath, envir = load_env)
      pipeline_data$loaded_data <- as.list(load_env)
      
      # Handle pheno_combat -> traits alias if traits missing (common in 02_output.RData)
      if(is.null(pipeline_data$loaded_data$traits) && !is.null(pipeline_data$loaded_data$pheno_combat)) {
          pipeline_data$loaded_data$traits <- pipeline_data$loaded_data$pheno_combat
      }
      
      # Merge PCA scores if provided
      if (!is.null(input$pca_file_input)) {
        pca_df <- read.csv(input$pca_file_input$datapath, stringsAsFactors = FALSE)
        
        # Ensure SampleID column exists
        if(!"SampleID" %in% colnames(pca_df)) colnames(pca_df)[1] <- "SampleID"
        
        # Deduplicate
        pca_df <- pca_df[!duplicated(pca_df$SampleID), ]
        
        if(!is.null(pipeline_data$loaded_data$traits)) {
          traits <- pipeline_data$loaded_data$traits
          matched <- match(rownames(traits), pca_df$SampleID)
          
          # Merge PC columns
          for(col in intersect(colnames(pca_df), c("PC1", "PC2", "PC3"))) {
            traits[[col]] <- pca_df[[col]][matched]
          }
          pipeline_data$loaded_data$traits <- traits
          showNotification("Merged PCA scores into phenotype data.", type = "message")
        }
      }
      
      loaded_vars <- names(pipeline_data$loaded_data)
      pipeline_data$status <- paste(
        "Successfully loaded:", input$wgcna_file_input$name,
        "\nContains objects:", paste(loaded_vars, collapse = ", ")
      )
      
      # Reset results from previous runs to clear plots
      pipeline_data$sft_results <- NULL
      pipeline_data$wgcna_results <- NULL
      pipeline_data$trait_cor_results <- NULL
      
      # Update GSEA module selector if moduleColors exist
      if ("moduleColors" %in% loaded_vars) {
        modules <- setdiff(unique(pipeline_data$loaded_data$moduleColors), "grey")
        updateSelectInput(session, "gsea_module", choices = sort(modules), selected = modules[1])
        updateSelectInput(session, "ora_module", choices = sort(modules), selected = modules[1])
      } else {
        updateSelectInput(session, "gsea_module", choices = character(0))
        updateSelectInput(session, "ora_module", choices = character(0))
      }
      
      # Setup WGCNA heatmap if results exist
      if ("moduleTraitCor" %in% loaded_vars && "moduleTraitP" %in% loaded_vars) {
        pipeline_data$trait_cor_results <- list(
          cor = pipeline_data$loaded_data$moduleTraitCor,
          pval = pipeline_data$loaded_data$moduleTraitP
        )
        all_traits <- colnames(pipeline_data$loaded_data$moduleTraitCor)
        want_traits <- c("Sex","Age_at_draw","PC1","PC2","PC3","SBP","d_cmb","fsrp","htnscore","lacunar_stroke","dm","bmi","stroke","caa_prob","caa","AD_status","cog_status","E4status")
        updateSelectInput(session, "wgcna_traits", choices = all_traits, selected = intersect(want_traits, all_traits))
      }
      
      showNotification("Data loaded successfully!", type = "message")
    }, error = function(e) {
      pipeline_data$status <- paste("ERROR loading file:", conditionMessage(e))
      showNotification(pipeline_data$status, type = "error")
    })
  })
  
  output$data_status_output <- renderText({
    pipeline_data$status
  })
  
  # --- (Optional) WGCNA tab placeholders ---
  output$wgcna_sample_tree <- renderPlot({
    req(pipeline_data$wgcna_filtered_data)
    plot_sample_tree(pipeline_data$wgcna_filtered_data$expr)
  })

  output$wgcna_sft_fit_plot <- renderPlot({
    req(pipeline_data$sft_results)
    plot_scale_free_topology(pipeline_data$sft_results)
  })
  
  output$wgcna_sft_mean_plot <- renderPlot({
    req(pipeline_data$sft_results)
    plot_mean_connectivity(pipeline_data$sft_results)
  })
  
  output$wgcna_dendro_plot <- renderPlot({
    req(pipeline_data$wgcna_results)
    plot_dendro_and_colors(
      pipeline_data$wgcna_results$geneTree,
      pipeline_data$wgcna_results$moduleColors
    )
  })
  
  output$wgcna_eigengene_plot <- renderPlot({
    req(pipeline_data$wgcna_results)
    plot_eigengene_network(pipeline_data$wgcna_results$MEs)
  })
  
  output$wgcna_trait_heatmap <- renderPlot({
    req(pipeline_data$trait_cor_results)
    
    cor_mat <- pipeline_data$trait_cor_results$cor
    pval_mat <- pipeline_data$trait_cor_results$pval
    
    if (!is.null(input$wgcna_traits) && length(input$wgcna_traits) > 0) {
      keep_cols <- intersect(input$wgcna_traits, colnames(cor_mat))
      if (length(keep_cols) > 0) {
        cor_mat <- cor_mat[, keep_cols, drop = FALSE]
        pval_mat <- pval_mat[, keep_cols, drop = FALSE]
      }
    }
    plot_module_trait_heatmap(cor_mat, pval_mat)
  })
  
  # --- WGCNA Download Handlers ---
  output$download_sample_tree <- downloadHandler(
    filename = function() { "wgcna_sample_clustering.pdf" },
    content = function(file) {
      req(pipeline_data$wgcna_filtered_data)
      pdf(file, width = 10, height = 7)
      plot_sample_tree(pipeline_data$wgcna_filtered_data$expr)
      dev.off()
    }
  )

  output$download_sft_fit <- downloadHandler(
    filename = function() { "wgcna_scale_free_fit.pdf" },
    content = function(file) {
      req(pipeline_data$sft_results)
      pdf(file, width = 8, height = 6)
      plot_scale_free_topology(pipeline_data$sft_results)
      dev.off()
    }
  )
  
  output$download_sft_mean <- downloadHandler(
    filename = function() { "wgcna_mean_connectivity.pdf" },
    content = function(file) {
      req(pipeline_data$sft_results)
      pdf(file, width = 8, height = 6)
      plot_mean_connectivity(pipeline_data$sft_results)
      dev.off()
    }
  )
  
  output$download_dendro <- downloadHandler(
    filename = function() { "wgcna_dendrogram.pdf" },
    content = function(file) {
      req(pipeline_data$wgcna_results)
      pdf(file, width = 10, height = 7)
      plot_dendro_and_colors(
        pipeline_data$wgcna_results$geneTree,
        pipeline_data$wgcna_results$moduleColors
      )
      dev.off()
    }
  )
  
  output$download_eigengene <- downloadHandler(
    filename = function() { "wgcna_eigengene_clustering.pdf" },
    content = function(file) {
      req(pipeline_data$wgcna_results)
      pdf(file, width = 8, height = 6)
      plot_eigengene_network(pipeline_data$wgcna_results$MEs)
      dev.off()
    }
  )
  
  output$download_heatmap <- downloadHandler(
    filename = function() { "wgcna_heatmap.pdf" },
    content = function(file) {
      req(pipeline_data$trait_cor_results)
      pdf(file, width = 12, height = 10)
      
      cor_mat <- pipeline_data$trait_cor_results$cor
      pval_mat <- pipeline_data$trait_cor_results$pval
      
      if (!is.null(input$wgcna_traits) && length(input$wgcna_traits) > 0) {
        keep_cols <- intersect(input$wgcna_traits, colnames(cor_mat))
        if (length(keep_cols) > 0) {
          cor_mat <- cor_mat[, keep_cols, drop = FALSE]
          pval_mat <- pval_mat[, keep_cols, drop = FALSE]
        }
      }
      plot_module_trait_heatmap(cor_mat, pval_mat)
      dev.off()
    }
  )
  
  # --- Server Logic for Tab 5: WGCNA ---
  observeEvent(input$run_sft_btn, {
    # Determine input data (prefer pipeline data, fallback to loaded data)
    expr_to_use <- if (!is.null(pipeline_data$harmonized_expr)) {
      as.data.frame(pipeline_data$harmonized_expr)
    } else if (!is.null(pipeline_data$loaded_data$expr)) {
      pipeline_data$loaded_data$expr
    } else {
      NULL
    }
    
    if (is.null(expr_to_use)) {
      showNotification("No expression data found. Run Harmonization (Tab 2) or load WGCNA data (Tab 4).", type = "error")
      return()
    }
    
    # Get traits for alignment/QC
    traits_to_use <- if (!is.null(pipeline_data$processed_data$pheno)) {
      pipeline_data$processed_data$pheno
    } else {
      pipeline_data$loaded_data$traits
    }
    
    withProgress(message = 'Analyzing Soft Threshold', value = 0, {
      update_prog <- function(detail) {
        incProgress(0.5, detail = detail)
      }
      
      sft_res <- tryCatch({
        # Run QC first
        qc_res <- perform_wgcna_qc(expr_to_use, traits_to_use)
        
        # Store filtered data
        pipeline_data$wgcna_filtered_data <- qc_res
        
        analyze_soft_threshold(
          expr = qc_res$expr,
          updateProgress = update_prog
        )
      }, error = function(e) {
        showNotification(paste("Soft Threshold analysis failed:", conditionMessage(e)), type = "error")
        return(NULL)
      })
      
      pipeline_data$sft_results <- sft_res
      if (!is.null(sft_res)) showNotification("Soft Threshold analysis complete.", type = "message")
    })
  })
  
  observeEvent(input$run_wgcna_btn, {
    # Determine input data
    expr_to_use <- if (!is.null(pipeline_data$harmonized_expr)) {
      as.data.frame(pipeline_data$harmonized_expr)
    } else if (!is.null(pipeline_data$loaded_data$expr)) {
      pipeline_data$loaded_data$expr
    } else {
      NULL
    }
    
    if (is.null(expr_to_use)) {
      showNotification("No expression data found. Run Harmonization (Tab 2) or load WGCNA data (Tab 4).", type = "error")
      return()
    }
    
    withProgress(message = 'Running WGCNA', value = 0, {
      
      # Define update function to pass to helper
      update_prog <- function(detail) {
        incProgress(0.2, detail = detail)
      }
      
      tryCatch({
        # Get traits for alignment/QC
        traits_to_use <- if (!is.null(pipeline_data$processed_data$pheno)) {
          pipeline_data$processed_data$pheno
        } else {
          pipeline_data$loaded_data$traits
        }
        
        # Run QC (or use cached if available from SFT step, but safer to re-run to ensure consistency)
        qc_res <- perform_wgcna_qc(expr_to_use, traits_to_use)
        pipeline_data$wgcna_filtered_data <- qc_res
        expr_to_use <- qc_res$expr
        traits_to_use <- qc_res$traits
        
        res <- run_wgcna(
          expr = expr_to_use,
          softPower = input$wgcna_softPower,
          minClusterSize = input$wgcna_minClusterSize,
          cutHeight = input$wgcna_cutHeight,
          updateProgress = update_prog
        )
        
        # Update pipeline data
        pipeline_data$wgcna_results <- res
        
        # Merge PCA results if available (e.g. from Tab 3)
        if (!is.null(pipeline_data$pca_results) && !is.null(traits_to_use)) {
          pca_scores <- as.data.frame(pipeline_data$pca_results$x)
          # Keep top 5 PCs
          pca_scores <- pca_scores[, 1:min(5, ncol(pca_scores)), drop=FALSE]
          
          # Match samples using merge to handle potential rowname mismatches or subsets
          # First ensure SampleID column exists or use rownames
          if(!"SampleID" %in% colnames(traits_to_use)) traits_to_use$SampleID <- rownames(traits_to_use)
          pca_scores$SampleID <- rownames(pca_scores)
          
          # Remove existing PC cols from traits to avoid duplicates
          pc_cols <- setdiff(colnames(pca_scores), "SampleID")
          traits_to_use <- traits_to_use[, !colnames(traits_to_use) %in% pc_cols, drop=FALSE]
          
          # Merge
          traits_to_use <- merge(traits_to_use, pca_scores, by="SampleID", all.x=TRUE)
          rownames(traits_to_use) <- traits_to_use$SampleID
          traits_to_use$SampleID <- NULL
        }
        
        if (!is.null(traits_to_use)) {
          # Calculate correlations for ALL numeric traits available
          # Filtering happens at the visualization stage
          
          incProgress(0.1, detail = "Calculating Module-Trait Correlations")
          cor_res <- calculate_module_trait_relationships(res$MEs, traits_to_use)
          if (!is.null(cor_res)) {
            # Store results for heatmap
            pipeline_data$trait_cor_results <- cor_res
            # Store traits used for saving later
            pipeline_data$wgcna_results$traits <- traits_to_use
            
            # Update UI selection with defaults
            all_traits <- colnames(cor_res$cor)
            want_traits <- c("Sex","Age_at_draw","PC1","PC2","PC3","SBP","d_cmb","fsrp","htnscore","lacunar_stroke","dm","bmi","stroke","caa_prob","caa","AD_status","cog_status","E4status")
            selected <- intersect(want_traits, all_traits)
            if(length(selected) == 0) selected <- all_traits[1:min(5, length(all_traits))]
            
            updateSelectInput(session, "wgcna_traits", choices = all_traits, selected = selected)
          }
        }
        
        # Update GSEA module selector
        modules <- setdiff(unique(res$moduleColors), "grey")
        updateSelectInput(session, "gsea_module", choices = sort(modules), selected = modules[1])
        updateSelectInput(session, "ora_module", choices = sort(modules), selected = modules[1])

        showNotification("WGCNA analysis completed successfully.", type = "message")
        
      }, error = function(e) {
        showNotification(paste("WGCNA failed:", conditionMessage(e)), type = "error")
      })
    })
  })
  
  output$dl_wgcna_results <- downloadHandler(
    filename = function() {
      paste0("WGCNA_Results_", Sys.Date(), ".RData")
    },
    content = function(file) {
      req(pipeline_data$wgcna_results)
      
      # Prepare objects to match the expected input for downstream steps
      expr <- if(!is.null(pipeline_data$harmonized_expr)) as.data.frame(pipeline_data$harmonized_expr) else pipeline_data$loaded_data$expr
      MEs <- pipeline_data$wgcna_results$MEs
      moduleColors <- pipeline_data$wgcna_results$moduleColors
      geneTree <- pipeline_data$wgcna_results$geneTree
      
      # Traits
      traits <- pipeline_data$wgcna_results$traits
      traits_interest <- traits # For simplicity, save all as interest if not subsetted
      
      # Correlations
      moduleTraitCor <- pipeline_data$trait_cor_results$cor
      moduleTraitP <- pipeline_data$trait_cor_results$pval
      
      save(expr, MEs, moduleColors, geneTree, traits, traits_interest, moduleTraitCor, moduleTraitP, file = file)
    }
  )
  
  # --- Server Logic for Tab 6: Module Scoring ---
  observeEvent(input$run_scoring_btn, {
    # Ensure data is loaded
    if (is.null(pipeline_data$loaded_data)) {
      showNotification("Please load WGCNA data in Tab 4 first.", type = "error")
      return()
    }
    
    needed <- c("moduleTraitCor", "moduleTraitP")
    if (!has_objects(pipeline_data$loaded_data, needed)) {
      showNotification(paste("Missing objects for scoring:", paste(setdiff(needed, names(pipeline_data$loaded_data)), collapse = ", ")),
                       type = "error")
      return()
    }
    
    if (!exists("run_module_scoring")) {
      showNotification("Function 'run_module_scoring' not found (helpers.R).", type = "error")
      return()
    }
    
    id <- showNotification("Running module scoring...", duration = NULL, closeButton = FALSE)
    on.exit(removeNotification(id), add = TRUE)
    
    scoring_results <- tryCatch({
      run_module_scoring(
        moduleTraitCor = pipeline_data$loaded_data$moduleTraitCor,
        moduleTraitP   = pipeline_data$loaded_data$moduleTraitP,
        r_thresh = input$scoring_r_thresh,
        q_thresh = input$scoring_q_thresh
      )
    }, error = function(e) {
      showNotification(paste("Scoring failed:", conditionMessage(e)), type = "error")
      return(NULL)
    })
    
    if (!is.null(scoring_results)) {
      output$scoring_table_output <- DT::renderDataTable({
        DT::datatable(scoring_results, options = list(pageLength = 10, scrollX = TRUE))
      })
    }
  })
  
  
  # --- Server Logic for Tab 7: Machine Learning ---
  observeEvent(input$run_ml_btn, {
    if (is.null(pipeline_data$loaded_data)) {
      showNotification("Please load WGCNA data in Tab 4 first.", type = "error")
      return()
    }
    
    needed <- c("MEs", "traits")
    if (!has_objects(pipeline_data$loaded_data, needed)) {
      showNotification(paste("Missing objects for ML:", paste(setdiff(needed, names(pipeline_data$loaded_data)), collapse = ", ")),
                       type = "error")
      return()
    }
    
    if (!exists("run_ml_analysis")) {
      showNotification("Function 'run_ml_analysis' not found (helpers.R).", type = "error")
      return()
    }
    
    id <- showNotification("Training and evaluating ML model...", duration = NULL, closeButton = FALSE)
    on.exit(removeNotification(id), add = TRUE)
    
    ml_results <- tryCatch({
      run_ml_analysis(
        MEs        = pipeline_data$loaded_data$MEs,
        traits     = pipeline_data$loaded_data$traits,
        target_col = input$ml_target_col,
        split_ratio = input$ml_split_ratio,
        alpha       = input$ml_alpha
      )
    }, error = function(e) {
      showNotification(paste("ML failed:", conditionMessage(e)), type = "error")
      return(NULL)
    })
    
    if (!is.null(ml_results)) {
      output$ml_cm_output <- renderPrint({
        print(ml_results$confusion_matrix)
      })
      
      output$ml_roc_plot <- renderPlot({
        plot(ml_results$roc_object, main = paste0("ROC Curve (AUC = ", round(ml_results$auc_value, 3), ")"))
      })
      
      output$ml_coeffs_output <- DT::renderDataTable({
        DT::datatable(ml_results$coefficients, options = list(pageLength = 5))
      })
    }
  })
  
  
  # --- Server Logic for Tab 8: GSEA ---
  observeEvent(input$run_gsea_btn, {
    if (is.null(pipeline_data$loaded_data)) {
      showNotification("Please load WGCNA data in Tab 4 first.", type = "error")
      return()
    }
    
    if (!exists("run_gsea_analysis")) {
      showNotification("Function 'run_gsea_analysis' not found (helpers.R).", type = "error")
      return()
    }
    if (!exists("plot_gsea_results")) {
      showNotification("Function 'plot_gsea_results' not found (helpers.R).", type = "error")
      return()
    }
    
    # Require user to upload the analyte info file for this session
    if (is.null(input$analyte_info_file)) {
      showNotification("Please upload the Analyte Info TSV file to run GSEA.", type = "error")
      return()
    }
    analyte_file_path <- input$analyte_info_file$datapath
    
    if (is.null(input$gsea_module) || !nzchar(input$gsea_module)) {
      showNotification("Please select a module for GSEA.", type = "error")
      return()
    }
    
    id <- showNotification(paste("Running GSEA for module:", input$gsea_module), duration = NULL, closeButton = FALSE)
    on.exit(removeNotification(id), add = TRUE)
    
    gsea_results <- tryCatch({
      run_gsea_analysis(
        wgcna_objects    = pipeline_data$loaded_data,
        analyte_info_path = analyte_file_path,
        target_module     = input$gsea_module,
        min_size = input$gsea_minSize,
        max_size = input$gsea_maxSize
      )
    }, error = function(e) {
      showNotification(paste("GSEA failed:", conditionMessage(e)), type = "error")
      return(NULL)
    })
    
    if (is.null(gsea_results)) return()
    
    # Filter and render results
    sig_gsea <- tryCatch({
      dplyr::filter(gsea_results, padj < input$gsea_padj_cutoff)
    }, error = function(e) {
      showNotification("Unexpected GSEA results format.", type = "error")
      return(NULL)
    })
    
    if (!is.null(sig_gsea)) {
      output$gsea_table_output <- DT::renderDataTable({
        DT::datatable(sig_gsea, options = list(pageLength = 10, scrollX = TRUE))
      })
      
      output$gsea_dotplot_output <- renderPlot({
        plot_gsea_results(sig_gsea, input$gsea_module)
      })
    }
  })
  
  
  # --- Server Logic for Tab 9: ORA ---
  observeEvent(input$run_ora_btn, {
    if (is.null(pipeline_data$loaded_data)) {
      showNotification("Please load WGCNA data in Tab 4 first.", type = "error")
      return()
    }
    
    if (!exists("run_ora_analysis")) {
      showNotification("Function 'run_ora_analysis' not found (helpers.R).", type = "error")
      return()
    }
    
    # Require user to upload the analyte info file for this session
    if (is.null(input$analyte_info_file)) {
      showNotification("Please upload the Analyte Info TSV file on the GSEA tab to run ORA.", type = "error")
      return()
    }
    analyte_file_path <- input$analyte_info_file$datapath
    
    if (is.null(input$ora_module) || !nzchar(input$ora_module)) {
      showNotification("Please select a module for ORA.", type = "error")
      return()
    }
    
    id <- showNotification(paste("Running ORA for module:", input$ora_module), duration = NULL, closeButton = FALSE)
    on.exit(removeNotification(id), add = TRUE)
    
    ora_results <- tryCatch({
      run_ora_analysis(
        wgcna_objects    = pipeline_data$loaded_data,
        analyte_info_path = analyte_file_path,
        target_module     = input$ora_module
      )
    }, error = function(e) {
      showNotification(paste("ORA failed:", conditionMessage(e)), type = "error")
      return(NULL)
    })
    
    if (is.null(ora_results)) return()
    
    # Filter and render results
    sig_ora <- dplyr::filter(ora_results, p.adjust < input$ora_padj_cutoff)
    
    output$ora_table_output <- DT::renderDataTable({
      DT::datatable(sig_ora, options = list(pageLength = 10, scrollX = TRUE))
    })
    
    output$ora_dotplot_output <- renderPlot({
      plot_ora_results(sig_ora, input$ora_module)
    })
  })
}

# 4. --- RUN THE APP ---
# -----------------------
shinyApp(ui = ui, server = server)
