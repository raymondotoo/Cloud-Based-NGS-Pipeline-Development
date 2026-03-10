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
    message("Installing missing packages: ", paste(missing_pkgs, collapse = ", "))
    utils::install.packages(missing_pkgs, repos = "https://cloud.r-project.org")
  }

  invisible(lapply(pkgs, function(p) {
    suppressPackageStartupMessages(library(p, character.only = TRUE))
  }))
}

# All required packages for the app
required_pkgs <- c(
  "WGCNA", "dplyr", "readr", "clusterProfiler", "org.Hs.eg.db", "enrichplot",
  "ReactomePA", "ggplot2", "openxlsx", "fgsea", "msigdbr", "tidyr",
  "tibble", "pheatmap", "stringr", "caret", "glmnet", "pROC", "DT", "sva",
  "imputeLCMD", "ggrepel", "impute", "randomForest"
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
      .analysis-card { background: #f8fafc; border: 1px solid #dbe4ee; border-radius: 12px; padding: 14px 16px; margin-bottom: 14px; }
      .analysis-card h4 { margin-top: 0; font-weight: 700; }
      .section-note { color: #475569; font-size: 13px; }
      .pipeline-footer { margin: 20px 0 12px; text-align: center; font-size: 12px; color: #6b7280; }
    "))
  ),
  useShinyjs(),
  navbarPage(
             title = "ADNI-ADRC Proteomics Pipeline",
             # --- Tab 1: Inputs & Preprocessing ---
             tabPanel("1. Inputs",
                      sidebarLayout(
                        sidebarPanel(
                          h4("Required Inputs"),
                          p("Upload the source files once here. Downstream steps derive PCA, ComBat, and WGCNA objects from this run."),
                          fileInput("raw_proto_file", "Proteomics Data (Expression)", accept = c(".xlsx", ".csv")),
                          fileInput("raw_pheno_file", "Phenotype Data", accept = c(".xlsx", ".csv")),
                          fileInput("analyte_info_init", "Analyte Annotation (TSV)", accept = c(".tsv", ".txt")),
                          fileInput("background_init", "Background Gene List (Optional CSV)", accept = c(".csv")),
                          selectInput("impute_method", "Imputation Method", choices = c("Hybrid", "KNN", "Minimum", "None"), selected = "Hybrid"),
                          actionButton("run_preprocess_btn", "Run Preprocessing", class = "btn-primary"),
                          hr(),
                          wellPanel(
                            h4("About This Step"),
                            p("This step reproduces the original workflow: ID cleanup, duplicate removal, SOMA panel filtering when available, WMH PCA derivation from phenotype, imputation of proteomics missingness, and alignment of phenotype with expression samples."),
                            p("The canonical WMH PCA uses wmh_juxta, wmh_frontal, wmh_pv, wmh_parietal, and wmh_posterior, with PC1 sign flipped before merging into phenotype.")
                          )
                        ),
                        mainPanel(
                          h3("Input Status"),
                          verbatimTextOutput("input_status_output"),
                          hr(),
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
                            column(6, h4("Before Correction"), plotOutput("pca_before_combat"), downloadButton("dl_combat_before_plot", "Download Before Plot")),
                            column(6, h4("After Correction"), plotOutput("pca_after_combat"), downloadButton("dl_combat_after_plot", "Download After Plot"))
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
             
             # --- Tab 4: Pipeline Status ---
             tabPanel("4. Pipeline Status",
                      sidebarLayout(
                        sidebarPanel(
                          wellPanel(
                              h4("About This Step"),
                              p("This tab reports the in-app pipeline state. No prior WGCNA or PCA output upload is required."),
                              p("Run the tabs in order: preprocessing, harmonization, optional PCA review, WGCNA, then downstream scoring, ML, GSEA, and ORA.")
                          )
                        ),
                        mainPanel(
                          h3("Workflow Overview"),
                          tags$img(src = "workflow.svg", style = "max-width: 100%; border: 1px solid #e2e8f0; border-radius: 12px;"),
                          hr(),
                          h3("Pipeline Status"),
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
                          tabsetPanel(
                            tabPanel("Quality Control",
                                     div(class = "analysis-card",
                                         h4("Sample Clustering"),
                                         p(class = "section-note", "Hierarchical clustering of samples after WGCNA input QC."),
                                         plotOutput("wgcna_sample_tree", height = "420px"),
                                         downloadButton("download_sample_tree", "Download Plot")
                                     )),
                            tabPanel("Soft Threshold",
                                     div(class = "analysis-card",
                                         h4("Scale-Free Topology"),
                                         p(class = "section-note", "Used to choose a soft-threshold power consistent with scale-free network behavior."),
                                         plotOutput("wgcna_sft_fit_plot", height = "420px"),
                                         downloadButton("download_sft_fit", "Download Plot")
                                     ),
                                     div(class = "analysis-card",
                                         h4("Mean Connectivity"),
                                         p(class = "section-note", "Connectivity profile across candidate powers."),
                                         plotOutput("wgcna_sft_mean_plot", height = "420px"),
                                         downloadButton("download_sft_mean", "Download Plot")
                                     )),
                            tabPanel("Modules",
                                     div(class = "analysis-card",
                                         h4("Gene Dendrogram and Module Colors"),
                                         p(class = "section-note", "Dynamic tree cut modules and merged color assignments."),
                                         plotOutput("wgcna_dendro_plot", height = "620px"),
                                         downloadButton("download_dendro", "Download Dendrogram")
                                     ),
                                     div(class = "analysis-card",
                                         h4("Eigengene Clustering"),
                                         p(class = "section-note", "Relationships among module eigengenes after module construction."),
                                         plotOutput("wgcna_eigengene_plot", height = "420px"),
                                         downloadButton("download_eigengene", "Download Eigengene Plot")
                                     )),
                            tabPanel("Trait Heatmap",
                                     div(class = "analysis-card",
                                         h4("Module-Trait Relationships"),
                                         p(class = "section-note", "Correlations between module eigengenes and selected phenotype traits."),
                                         plotOutput("wgcna_trait_heatmap", height = "820px"),
                                         downloadButton("download_heatmap", "Download Heatmap")
                                     ))
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
                          selectInput(
                            "ml_model",
                            "Model",
                            choices = c(
                              "Elastic Net (glmnet)" = "elastic_net",
                              "Logistic Regression" = "logistic",
                              "Random Forest" = "random_forest"
                            ),
                            selected = "elastic_net"
                          ),
                          selectInput("ml_target_col", "Target Variable", choices = c("stroke", "AD_status", "caa")),
                          sliderInput("ml_split_ratio", "Train/Test Split Ratio", value = 0.7, min = 0.5, max = 0.9),
                          conditionalPanel(
                            condition = "input.ml_model == 'elastic_net'",
                            sliderInput("ml_alpha", "Elastic Net Alpha", value = 0.5, min = 0, max = 1)
                          ),
                          actionButton("run_ml_btn", "Run ML Model", class = "btn-primary"),
                          hr(),
                          wellPanel(
                              h4("About This Step"),
                              p("This step uses the module eigengenes (MEs) as features to predict a clinical outcome (e.g., 'stroke'). Using MEs instead of thousands of individual proteins drastically reduces the number of features, mitigates multiple testing issues, and incorporates biological structure into the model."),
                              p("Multiple model options are available. Elastic Net remains the default because it balances interpretability with predictive performance on correlated module eigengenes."),
                              h4("Computation Notes"),
                              tags$ul(
                                tags$li("Elastic Net: 5-fold CV on standardized eigengenes using glmnet; coefficients are reported at lambda.min."),
                                tags$li("Logistic Regression: GLM with binomial link on standardized eigengenes; coefficients are reported directly."),
                                tags$li("Random Forest: default randomForest settings on standardized eigengenes; variable importance is reported.")
                              ),
                              p("All models use a train/test split and are evaluated with ROC AUC and a 0.5 decision threshold for the confusion matrix.")
                          )
                        ),
                        mainPanel(
                          h3("Model Performance"),
                          tabsetPanel(
                            tabPanel("Summary",
                                     div(class = "analysis-card",
                                         h4("Metrics"),
                                         tableOutput("ml_metrics_output")
                                     ),
                                     div(class = "analysis-card",
                                         h4("Confusion Matrix"),
                                         verbatimTextOutput("ml_cm_output")
                                     )),
                            tabPanel("ROC Curve",
                                     div(class = "analysis-card",
                                         h4("ROC Curve"),
                                         plotOutput("ml_roc_plot", height = "420px"),
                                         downloadButton("dl_ml_roc_plot", "Download ROC Plot")
                                     )),
                            tabPanel("Top Features",
                                     div(class = "analysis-card",
                                         h4("Feature Weights or Importance"),
                                         DT::dataTableOutput("ml_coeffs_output")
                                     ))
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
                          selectInput("gsea_source", "Database View", choices = c("All", "Hallmark", "Reactome", "GO_BP"), selected = "All"),
                          numericInput("gsea_minSize", "Min. Gene Set Size", value = 15, min = 5),
                          numericInput("gsea_maxSize", "Max. Gene Set Size", value = 500, min = 100),
                          sliderInput("gsea_padj_cutoff", "Significance Cutoff (padj)", value = 0.05, min = 0.01, max = 0.25),
                          selectInput("gsea_pathway", "Pathway for Enrichment Curve", choices = NULL),
                          actionButton("run_gsea_btn", "Run GSEA", class = "btn-primary"),
                          hr(),
                          wellPanel(
                              h4("About This Step"),
                              p("Gene Set Enrichment Analysis (GSEA) is used to identify biological pathways or functions that are significantly enriched within a module. Unlike simple Over-Representation Analysis (ORA), GSEA considers the entire ranked list of proteins."),
                              p("In this pipeline, proteins are ranked by their Module Membership (kME) value for a selected module. GSEA then determines if known biological pathways (from databases like GO, Reactome, or Hallmark) are enriched at the top of this ranked list."),
                              p("Analyte annotation comes from the upfront input tab.")
                          )
                        ),
                        mainPanel(
                          div(style = "height: 800px; overflow-y: scroll; margin-top: 20px;",
                              h3("GSEA Results"),
                              div(class = "analysis-card",
                                  h4("Pathway Overview"),
                                  p(class = "section-note", "A ranked summary of significant pathway enrichments across Hallmark, Reactome, and GO biological process gene sets."),
                                  plotOutput("gsea_dotplot_output", height = "520px"),
                                  downloadButton("dl_gsea_overview_plot", "Download Overview Plot")
                              ),
                              tabsetPanel(
                                tabPanel("Network",
                                         div(class = "analysis-card",
                                             h4("Leading-Edge Network"),
                                             p(class = "section-note", "Connections between the strongest enriched pathways and the proteins driving them."),
                                             plotOutput("gsea_network_output", height = "620px"),
                                             downloadButton("dl_gsea_network_plot", "Download Network Plot")
                                         )),
                                tabPanel("Enrichment Curve",
                                         div(class = "analysis-card",
                                             h4("Running Score Plot"),
                                             p(class = "section-note", "The FGSEA running enrichment curve for a selected pathway."),
                                             plotOutput("gsea_curve_output", height = "500px"),
                                             downloadButton("dl_gsea_curve_plot", "Download Enrichment Plot")
                                         )),
                                tabPanel("Results Table",
                                         div(class = "analysis-card",
                                             h4("Detailed Results"),
                                             DT::dataTableOutput("gsea_table_output")
                                         ))
                              )
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
                          selectInput("ora_source", "Database View", choices = c("All", "GO_BP", "KEGG", "Reactome"), selected = "All"),
                          sliderInput("ora_padj_cutoff", "Significance Cutoff (padj)", value = 0.05, min = 0.01, max = 0.25),
                          actionButton("run_ora_btn", "Run ORA", class = "btn-primary"),
                          hr(),
                          wellPanel(
                              h4("About This Step"),
                              p("Over-Representation Analysis (ORA) determines whether known biological functions or pathways are over-represented in a given list of genes/proteins. Unlike GSEA, ORA requires a discrete list of interesting molecules (e.g., all proteins within a module) and compares it against a background universe."),
                              p("This analysis tests for enrichment in Gene Ontology (GO), KEGG, and Reactome pathway databases."),
                              p("Analyte annotation comes from the upfront input tab.")
                          )
                        ),
                        mainPanel(
                          div(style = "height: 800px; overflow-y: scroll; margin-top: 20px;",
                              h3("ORA Results"),
                              div(class = "analysis-card",
                                  h4("Term Overview"),
                                  p(class = "section-note", "Prioritized GO, KEGG, and Reactome terms for the selected module."),
                                  plotOutput("ora_dotplot_output", height = "520px"),
                                  downloadButton("dl_ora_overview_plot", "Download Overview Plot")
                              ),
                              tabsetPanel(
                                tabPanel("Concept Network",
                                         div(class = "analysis-card",
                                             h4("Functional Network"),
                                             p(class = "section-note", "A cnet-style view linking enriched terms to the genes in the module."),
                                             plotOutput("ora_cnet_output", height = "620px"),
                                             downloadButton("dl_ora_cnet_plot", "Download Network Plot")
                                         )),
                                tabPanel("Hub Genes",
                                         div(class = "analysis-card",
                                             h4("Pathway Hub Genes"),
                                             p(class = "section-note", "Proteins with strong module membership supporting significant GO terms."),
                                             DT::dataTableOutput("ora_hub_table_output")
                                         )),
                                tabPanel("Results Table",
                                         div(class = "analysis-card",
                                             h4("Detailed Results"),
                                             DT::dataTableOutput("ora_table_output")
                                         ))
                              )
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
             ),
             tabPanel("10. Supplementary",
                      sidebarLayout(
                        sidebarPanel(
                          h4("Supplementary Analyses"),
                          p("These outputs reproduce the core elements of `08_Supplementary_analysis_v3.Rmd`."),
                          selectInput("supp_modules", "Key Modules", choices = NULL, multiple = TRUE),
                          numericInput("supp_permutations", "Preservation Permutations", value = 100, min = 25, step = 25),
                          actionButton("run_supp_btn", "Run Supplementary Analysis", class = "btn-primary"),
                          hr(),
                          wellPanel(
                            h4("Included Outputs"),
                            p("Module preservation across ADNI and ADRC, AD-adjusted sensitivity models, module-membership hub proteins, and cohort-stratified module-trait heatmaps.")
                          )
                        ),
                        mainPanel(
                          div(style = "height: 850px; overflow-y: scroll; margin-top: 20px;",
                              h3("Supplementary Results"),
                              tabsetPanel(
                                tabPanel("Preservation",
                                         div(class = "analysis-card",
                                             h4("Module Preservation"),
                                             plotOutput("supp_preservation_plot", height = "500px"),
                                             downloadButton("dl_supp_preservation_plot", "Download Preservation Plot"),
                                             DT::dataTableOutput("supp_preservation_table")
                                         )),
                                tabPanel("Sensitivity",
                                         div(class = "analysis-card",
                                             h4("AD-Adjusted Sensitivity"),
                                             plotOutput("supp_forest_plot", height = "500px"),
                                             downloadButton("dl_supp_forest_plot", "Download Sensitivity Plot"),
                                             DT::dataTableOutput("supp_sensitivity_table")
                                         )),
                                tabPanel("Hub Proteins",
                                         div(class = "analysis-card",
                                             h4("Module Membership Hubs"),
                                             plotOutput("supp_hub_plot", height = "560px"),
                                             downloadButton("dl_supp_hub_plot", "Download Hub Plot"),
                                             DT::dataTableOutput("supp_kme_table")
                                         )),
                                tabPanel("Cohort Heatmaps",
                                         fluidRow(
                                           column(12,
                                                  div(class = "analysis-card",
                                                      h4("ADNI"),
                                                      plotOutput("supp_heatmap_adni", height = "520px"),
                                                      downloadButton("dl_supp_heatmap_adni", "Download ADNI Heatmap")
                                                  )),
                                           column(12,
                                                  div(class = "analysis-card",
                                                      h4("ADRC"),
                                                      plotOutput("supp_heatmap_adrc", height = "520px"),
                                                      downloadButton("dl_supp_heatmap_adrc", "Download ADRC Heatmap")
                                                  ))
                                         ))
                              )
                          )
                        )
                      ))
  ),
  div(class = "pipeline-footer", "ADNI-ADRC proteomics analysis app")
)


# 3. --- SERVER LOGIC ---
# --------------------------

server <- function(input, output, session) {
  save_current_plot <- function(file, expr, width = 10, height = 7) {
    grDevices::pdf(file, width = width, height = height, useDingbats = FALSE)
    on.exit(grDevices::dev.off(), add = TRUE)
    force(expr)
  }
  
  # Reactive values to store data and results across steps
  pipeline_data <- reactiveValues(
    loaded_data = NULL,
    raw_proto = NULL,
    raw_pheno = NULL,
    analyte_info = NULL,
    background = NULL,
    processed_data = NULL,
    harmonized = NULL,
    sft_results = NULL,
    wgcna_filtered_data = NULL, # Store QC'd data
    pca_results = NULL,
    pca_samples = NULL,
    wgcna_results = NULL,
    trait_cor_results = NULL,
    ml_results = NULL,
    gsea_results = NULL,
    ora_results = NULL,
    supp_results = NULL,
    status = "App started. Upload phenotype, proteomics, and annotation inputs on Tab 1."
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
      
      pipeline_data$raw_proto <- proto
      pipeline_data$raw_pheno <- pheno
      pipeline_data$analyte_info <- if (!is.null(input$analyte_info_init)) input$analyte_info_init$datapath else NULL
      pipeline_data$background <- if (!is.null(input$background_init)) input$background_init$datapath else NULL

      # Run preprocessing
      res <- run_preprocessing(proto, pheno, input$impute_method)
      pipeline_data$processed_data <- res
      pipeline_data$pca_results <- if (!is.null(res$wmh_pca)) res$wmh_pca$pca else NULL
      pipeline_data$pca_samples <- if (!is.null(res$wmh_pca)) res$wmh_pca$scores$SampleID else NULL
      
      # Clear downstream results to prevent mismatch
      pipeline_data$harmonized <- NULL
      pipeline_data$harmonized_expr <- NULL
      pipeline_data$loaded_data <- NULL
      pipeline_data$wgcna_results <- NULL
      pipeline_data$trait_cor_results <- NULL
      pipeline_data$sft_results <- NULL
      pipeline_data$ml_results <- NULL
      pipeline_data$gsea_results <- NULL
      pipeline_data$ora_results <- NULL
      pipeline_data$supp_results <- NULL
      
      # Update UI choices for next steps
      batch_choices <- colnames(res$pheno)
      default_batch <- if("Cohort" %in% batch_choices) "Cohort" else batch_choices[1]
      updateSelectInput(session, "batch_col", choices = batch_choices, selected = default_batch)

      pheno_for_pca <- if (!is.null(res$pheno_raw)) res$pheno_raw else res$pheno
      nums <- names(pheno_for_pca)[sapply(pheno_for_pca, is.numeric)]
      default_pca_vars <- intersect(canonical_wmh_vars(), nums)
      if (length(default_pca_vars) == 0) default_pca_vars <- nums[1:min(5, length(nums))]
      updateSelectInput(session, "pca_vars", choices = nums, selected = default_pca_vars)
      updateSelectInput(session, "pca_color_col", choices = c("None", names(pheno_for_pca)))
      
      pipeline_data$status <- paste(
        "Inputs loaded.",
        "\nProteomics rows:", nrow(res$expr),
        "\nProteomics features:", ncol(res$expr),
        "\nPhenotype rows:", nrow(res$pheno),
        "\nWMH PCA derived:", !is.null(res$wmh_pca),
        "\nAnalyte annotation uploaded:", !is.null(pipeline_data$analyte_info)
      )
      showNotification("Preprocessing complete.", type = "message")
      
    }, error = function(e) {
      showNotification(paste("Preprocessing failed:", conditionMessage(e)), type = "error")
    })
  })
  
  # --- Reset Harmonization when Batch Column Changes ---
  observeEvent(input$batch_col, {
      pipeline_data$harmonized <- NULL
      pipeline_data$harmonized_expr <- NULL
      pipeline_data$loaded_data <- NULL
      pipeline_data$ml_results <- NULL
      pipeline_data$gsea_results <- NULL
      pipeline_data$ora_results <- NULL
      pipeline_data$supp_results <- NULL
  })

  output$input_status_output <- renderText({
    parts <- c(
      paste("Proteomics uploaded:", !is.null(input$raw_proto_file)),
      paste("Phenotype uploaded:", !is.null(input$raw_pheno_file)),
      paste("Analyte annotation uploaded:", !is.null(input$analyte_info_init)),
      paste("Background uploaded:", !is.null(input$background_init)),
      paste("Preprocessing complete:", !is.null(pipeline_data$processed_data))
    )
    paste(parts, collapse = "\n")
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
      pipeline_data$harmonized <- harmonized
      pipeline_data$harmonized_expr <- as.data.frame(harmonized$expr)
      pipeline_data$status <- paste(
        "Preprocessing complete.",
        "\nHarmonization complete.",
        "\nHarmonized samples:", nrow(harmonized$expr),
        "\nHarmonized features:", ncol(harmonized$expr)
      )
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
    plot_pca_batch(as.data.frame(pipeline_data$harmonized_expr), pipeline_data$harmonized$pheno, input$batch_col, "After ComBat")
  })
  
  output$dl_combat_plots <- downloadHandler(
    filename = "harmonization_plots.pdf",
    content = function(file) {
      req(pipeline_data$harmonized_expr)
      pdf(file, width = 12, height = 6)
      print(plot_pca_batch(pipeline_data$processed_data$expr, pipeline_data$processed_data$pheno, input$batch_col, "Before"))
      print(plot_pca_batch(as.data.frame(pipeline_data$harmonized_expr), pipeline_data$harmonized$pheno, input$batch_col, "After"))
      dev.off()
    }
  )

  output$dl_combat_before_plot <- downloadHandler(
    filename = "combat_before_pca.pdf",
    content = function(file) {
      req(pipeline_data$processed_data, input$batch_col)
      save_current_plot(file, print(plot_pca_batch(
        pipeline_data$processed_data$expr,
        pipeline_data$processed_data$pheno,
        input$batch_col,
        "Before ComBat"
      )), width = 8, height = 6)
    }
  )

  output$dl_combat_after_plot <- downloadHandler(
    filename = "combat_after_pca.pdf",
    content = function(file) {
      req(pipeline_data$harmonized_expr, input$batch_col)
      save_current_plot(file, print(plot_pca_batch(
        as.data.frame(pipeline_data$harmonized_expr),
        pipeline_data$harmonized$pheno,
        input$batch_col,
        "After ComBat"
      )), width = 8, height = 6)
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
      pheno_combat <- pipeline_data$harmonized$pheno[rownames(Proto_combat_z), , drop=FALSE]
      save(Proto_combat_z, pheno_combat, file = file)
    }
  )
  
  # --- Server Logic for Tab 3: PCA ---
  
  # Helper to get current phenotype data
  current_pheno <- reactive({
    if (!is.null(pipeline_data$processed_data$pheno_raw)) {
      return(pipeline_data$processed_data$pheno_raw)
    } else if (!is.null(pipeline_data$processed_data$pheno)) {
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
  
  output$data_status_output <- renderText({
    paste(
      pipeline_data$status,
      paste("Upfront analyte annotation available:", !is.null(pipeline_data$analyte_info)),
      paste("Preprocessed data available:", !is.null(pipeline_data$processed_data)),
      paste("Harmonized data available:", !is.null(pipeline_data$harmonized_expr)),
      paste("WGCNA objects available:", !is.null(pipeline_data$loaded_data)),
      sep = "\n"
    )
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
    expr_to_use <- if (!is.null(pipeline_data$harmonized_expr)) as.data.frame(pipeline_data$harmonized_expr) else NULL
    
    if (is.null(expr_to_use)) {
      showNotification("No harmonized expression data found. Run Harmonization first.", type = "error")
      return()
    }
    
    # Get traits for alignment/QC
    traits_to_use <- if (!is.null(pipeline_data$harmonized$pheno)) pipeline_data$harmonized$pheno else pipeline_data$processed_data$pheno
    
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
    expr_to_use <- if (!is.null(pipeline_data$harmonized_expr)) as.data.frame(pipeline_data$harmonized_expr) else NULL
    
    if (is.null(expr_to_use)) {
      showNotification("No harmonized expression data found. Run Harmonization first.", type = "error")
      return()
    }
    
    withProgress(message = 'Running WGCNA', value = 0, {
      
      # Define update function to pass to helper
      update_prog <- function(detail) {
        incProgress(0.2, detail = detail)
      }
      
      tryCatch({
        # Get traits for alignment/QC
        traits_to_use <- if (!is.null(pipeline_data$harmonized$pheno)) pipeline_data$harmonized$pheno else pipeline_data$processed_data$pheno
        
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

        pipeline_data$loaded_data <- list(
          expr = expr_to_use,
          MEs = res$MEs,
          moduleColors = res$moduleColors,
          geneTree = res$geneTree,
          traits = traits_to_use,
          moduleTraitCor = if (!is.null(pipeline_data$trait_cor_results)) pipeline_data$trait_cor_results$cor else NULL,
          moduleTraitP = if (!is.null(pipeline_data$trait_cor_results)) pipeline_data$trait_cor_results$pval else NULL
        )
        pipeline_data$gsea_results <- NULL
        pipeline_data$ora_results <- NULL
        pipeline_data$supp_results <- NULL
        pipeline_data$ml_results <- NULL
        
        # Update GSEA module selector
        modules <- setdiff(unique(res$moduleColors), "grey")
        updateSelectInput(session, "gsea_module", choices = sort(modules), selected = modules[1])
        updateSelectInput(session, "ora_module", choices = sort(modules), selected = modules[1])
        default_supp_modules <- intersect(c("red", "lightyellow", "black", "magenta"), modules)
        if (length(default_supp_modules) == 0) default_supp_modules <- head(sort(modules), 4)
        updateSelectInput(session, "supp_modules", choices = sort(modules), selected = default_supp_modules)
        pipeline_data$status <- paste(
          "Preprocessing complete.",
          "\nHarmonization complete.",
          "\nWGCNA complete.",
          "\nModules detected:", length(unique(res$moduleColors))
        )
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
      expr <- as.data.frame(pipeline_data$harmonized_expr)
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
      showNotification("Run WGCNA in this app before scoring modules.", type = "error")
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
      showNotification("Run WGCNA in this app before training the ML model.", type = "error")
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
        alpha       = if (!is.null(input$ml_alpha)) input$ml_alpha else 0.5,
        model_type  = input$ml_model
      )
    }, error = function(e) {
      showNotification(paste("ML failed:", conditionMessage(e)), type = "error")
      return(NULL)
    })
    
    if (!is.null(ml_results)) {
      pipeline_data$ml_results <- ml_results
      output$ml_metrics_output <- renderTable({
        ml_results$metrics
      })
      
      output$ml_cm_output <- renderPrint({
        print(ml_results$confusion_matrix)
      })
      
      output$ml_roc_plot <- renderPlot({
        plot(ml_results$roc_object, main = paste0(ml_results$model_label, " ROC (AUC = ", round(ml_results$auc_value, 3), ")"))
      })
      
      output$ml_coeffs_output <- DT::renderDataTable({
        DT::datatable(ml_results$coefficients, options = list(pageLength = 5))
      })
    }
  })

  output$dl_ml_roc_plot <- downloadHandler(
    filename = function() { paste0("ml_roc_", input$ml_target_col, "_", input$ml_model, ".pdf") },
    content = function(file) {
      req(pipeline_data$ml_results)
      save_current_plot(file, {
        plot(
          pipeline_data$ml_results$roc_object,
          main = paste0(pipeline_data$ml_results$model_label, " ROC (AUC = ", round(pipeline_data$ml_results$auc_value, 3), ")")
        )
      }, width = 8, height = 6)
    }
  )
  
  
  # --- Server Logic for Tab 8: GSEA ---
  observeEvent(input$run_gsea_btn, {
    if (is.null(pipeline_data$loaded_data)) {
      showNotification("Run WGCNA in this app before GSEA.", type = "error")
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
    
    analyte_file_path <- pipeline_data$analyte_info
    if (is.null(analyte_file_path) && file.exists(file.path("data", "CSF_SOMAscan7k_analyte_information.tsv"))) {
      analyte_file_path <- file.path("data", "CSF_SOMAscan7k_analyte_information.tsv")
    }
    if (is.null(analyte_file_path)) {
      showNotification("Upload analyte annotation on Tab 1 before running GSEA.", type = "error")
      return()
    }
    
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
    pipeline_data$gsea_results <- gsea_results
    sig_gsea <- get_sig_gsea_results(gsea_results, input$gsea_padj_cutoff, input$gsea_source)
    pathway_choices <- unique(sig_gsea$pathway)
    selected_pathway <- if (length(pathway_choices) > 0) pathway_choices[[1]] else character(0)
    updateSelectInput(session, "gsea_pathway", choices = pathway_choices, selected = selected_pathway)
  })

  observe({
    req(pipeline_data$gsea_results)
    sig_gsea <- get_sig_gsea_results(pipeline_data$gsea_results, input$gsea_padj_cutoff, input$gsea_source)
    pathway_choices <- unique(sig_gsea$pathway)
    current_choice <- isolate(input$gsea_pathway)
    if (!length(pathway_choices)) {
      updateSelectInput(session, "gsea_pathway", choices = character(0), selected = character(0))
    } else {
      selected <- if (!is.null(current_choice) && current_choice %in% pathway_choices) current_choice else pathway_choices[[1]]
      updateSelectInput(session, "gsea_pathway", choices = pathway_choices, selected = selected)
    }
  })

  output$gsea_table_output <- DT::renderDataTable({
    req(pipeline_data$gsea_results)
    sig_gsea <- get_sig_gsea_results(pipeline_data$gsea_results, input$gsea_padj_cutoff, input$gsea_source)
    DT::datatable(sig_gsea, options = list(pageLength = 10, scrollX = TRUE))
  })

  output$gsea_dotplot_output <- renderPlot({
    req(pipeline_data$gsea_results)
    plot_gsea_results(pipeline_data$gsea_results, input$gsea_padj_cutoff, input$gsea_source)
  })

  output$gsea_network_output <- renderPlot({
    req(pipeline_data$gsea_results)
    source_name <- if (is.null(input$gsea_source) || input$gsea_source == "All") "Hallmark" else input$gsea_source
    plot_gsea_network(pipeline_data$gsea_results, source_name, input$gsea_padj_cutoff)
  })

  output$gsea_curve_output <- renderPlot({
    req(pipeline_data$gsea_results)
    source_name <- if (is.null(input$gsea_source) || input$gsea_source == "All") "Hallmark" else input$gsea_source
    plot_gsea_enrichment_curve(
      pipeline_data$gsea_results,
      source_name = source_name,
      pathway_name = input$gsea_pathway,
      padj_cutoff = input$gsea_padj_cutoff
    )
  })

  output$dl_gsea_overview_plot <- downloadHandler(
    filename = function() { paste0("gsea_overview_", input$gsea_module, ".pdf") },
    content = function(file) {
      req(pipeline_data$gsea_results)
      save_current_plot(file, print(plot_gsea_results(
        pipeline_data$gsea_results,
        input$gsea_padj_cutoff,
        input$gsea_source
      )), width = 12, height = 8)
    }
  )

  output$dl_gsea_network_plot <- downloadHandler(
    filename = function() { paste0("gsea_network_", input$gsea_module, ".pdf") },
    content = function(file) {
      req(pipeline_data$gsea_results)
      source_name <- if (is.null(input$gsea_source) || input$gsea_source == "All") "Hallmark" else input$gsea_source
      save_current_plot(file, print(plot_gsea_network(
        pipeline_data$gsea_results,
        source_name,
        input$gsea_padj_cutoff
      )), width = 10, height = 8)
    }
  )

  output$dl_gsea_curve_plot <- downloadHandler(
    filename = function() { paste0("gsea_curve_", input$gsea_module, ".pdf") },
    content = function(file) {
      req(pipeline_data$gsea_results)
      source_name <- if (is.null(input$gsea_source) || input$gsea_source == "All") "Hallmark" else input$gsea_source
      save_current_plot(file, print(plot_gsea_enrichment_curve(
        pipeline_data$gsea_results,
        source_name = source_name,
        pathway_name = input$gsea_pathway,
        padj_cutoff = input$gsea_padj_cutoff
      )), width = 9, height = 6)
    }
  )
  
  
  # --- Server Logic for Tab 9: ORA ---
  observeEvent(input$run_ora_btn, {
    if (is.null(pipeline_data$loaded_data)) {
      showNotification("Run WGCNA in this app before ORA.", type = "error")
      return()
    }
    
    if (!exists("run_ora_analysis")) {
      showNotification("Function 'run_ora_analysis' not found (helpers.R).", type = "error")
      return()
    }
    
    analyte_file_path <- pipeline_data$analyte_info
    if (is.null(analyte_file_path) && file.exists(file.path("data", "CSF_SOMAscan7k_analyte_information.tsv"))) {
      analyte_file_path <- file.path("data", "CSF_SOMAscan7k_analyte_information.tsv")
    }
    if (is.null(analyte_file_path)) {
      showNotification("Upload analyte annotation on Tab 1 before running ORA.", type = "error")
      return()
    }
    
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
    pipeline_data$ora_results <- ora_results
  })

  output$ora_table_output <- DT::renderDataTable({
    req(pipeline_data$ora_results)
    sig_ora <- get_sig_ora_results(pipeline_data$ora_results, input$ora_padj_cutoff, input$ora_source)
    DT::datatable(sig_ora, options = list(pageLength = 10, scrollX = TRUE))
  })

  output$ora_dotplot_output <- renderPlot({
    req(pipeline_data$ora_results)
    plot_ora_results(pipeline_data$ora_results, input$ora_padj_cutoff, input$ora_source)
  })

  output$ora_cnet_output <- renderPlot({
    req(pipeline_data$ora_results)
    source_name <- if (is.null(input$ora_source) || input$ora_source == "All") "GO_BP" else input$ora_source
    plot_ora_cnet(pipeline_data$ora_results, source_name, input$ora_padj_cutoff)
  })

  output$ora_hub_table_output <- DT::renderDataTable({
    req(pipeline_data$ora_results)
    hub_df <- extract_ora_hub_genes(
      pipeline_data$ora_results,
      source_name = "GO_BP",
      padj_cutoff = input$ora_padj_cutoff
    )
    DT::datatable(hub_df, options = list(pageLength = 10, scrollX = TRUE))
  })

  output$dl_ora_overview_plot <- downloadHandler(
    filename = function() { paste0("ora_overview_", input$ora_module, ".pdf") },
    content = function(file) {
      req(pipeline_data$ora_results)
      save_current_plot(file, print(plot_ora_results(
        pipeline_data$ora_results,
        input$ora_padj_cutoff,
        input$ora_source
      )), width = 12, height = 8)
    }
  )

  output$dl_ora_cnet_plot <- downloadHandler(
    filename = function() { paste0("ora_network_", input$ora_module, ".pdf") },
    content = function(file) {
      req(pipeline_data$ora_results)
      source_name <- if (is.null(input$ora_source) || input$ora_source == "All") "GO_BP" else input$ora_source
      save_current_plot(file, print(plot_ora_cnet(
        pipeline_data$ora_results,
        source_name,
        input$ora_padj_cutoff
      )), width = 10, height = 8)
    }
  )

  # --- Server Logic for Tab 10: Supplementary ---
  observeEvent(input$run_supp_btn, {
    if (is.null(pipeline_data$loaded_data)) {
      showNotification("Run WGCNA in this app before supplementary analysis.", type = "error")
      return()
    }

    id <- showNotification("Running supplementary analyses...", duration = NULL, closeButton = FALSE)
    on.exit(removeNotification(id), add = TRUE)

    supp_results <- tryCatch({
      run_supplementary_analysis(
        wgcna_objects = pipeline_data$loaded_data,
        selected_modules = input$supp_modules,
        n_permutations = input$supp_permutations
      )
    }, error = function(e) {
      showNotification(paste("Supplementary analysis failed:", conditionMessage(e)), type = "error")
      return(NULL)
    })

    if (is.null(supp_results)) return()
    pipeline_data$supp_results <- supp_results
    showNotification("Supplementary analysis complete.", type = "message")
  })

  output$supp_preservation_plot <- renderPlot({
    req(pipeline_data$supp_results)
    print(plot_supp_preservation(pipeline_data$supp_results))
  })

  output$supp_preservation_table <- DT::renderDataTable({
    req(pipeline_data$supp_results)
    df <- pipeline_data$supp_results$preservation_summary
    DT::datatable(if (is.null(df)) data.frame() else df, options = list(pageLength = 10, scrollX = TRUE))
  })

  output$supp_forest_plot <- renderPlot({
    req(pipeline_data$supp_results)
    print(plot_supp_forest(pipeline_data$supp_results))
  })

  output$supp_sensitivity_table <- DT::renderDataTable({
    req(pipeline_data$supp_results)
    df <- pipeline_data$supp_results$comparison_table
    DT::datatable(if (is.null(df)) data.frame() else df, options = list(pageLength = 10, scrollX = TRUE))
  })

  output$supp_hub_plot <- renderPlot({
    req(pipeline_data$supp_results)
    print(plot_supp_top_hubs(pipeline_data$supp_results))
  })

  output$supp_kme_table <- DT::renderDataTable({
    req(pipeline_data$supp_results)
    top_kme <- pipeline_data$supp_results$protein_kme
    DT::datatable(utils::head(top_kme, 100), options = list(pageLength = 10, scrollX = TRUE))
  })

  output$supp_heatmap_adni <- renderPlot({
    req(pipeline_data$supp_results)
    plot_supp_cohort_heatmap(pipeline_data$supp_results, "ADNI")
  })

  output$supp_heatmap_adrc <- renderPlot({
    req(pipeline_data$supp_results)
    plot_supp_cohort_heatmap(pipeline_data$supp_results, "ADRC")
  })

  output$dl_supp_preservation_plot <- downloadHandler(
    filename = "supplementary_preservation.pdf",
    content = function(file) {
      req(pipeline_data$supp_results)
      save_current_plot(file, print(plot_supp_preservation(pipeline_data$supp_results)), width = 8, height = 6)
    }
  )

  output$dl_supp_forest_plot <- downloadHandler(
    filename = "supplementary_sensitivity.pdf",
    content = function(file) {
      req(pipeline_data$supp_results)
      save_current_plot(file, print(plot_supp_forest(pipeline_data$supp_results)), width = 9, height = 6)
    }
  )

  output$dl_supp_hub_plot <- downloadHandler(
    filename = "supplementary_hub_proteins.pdf",
    content = function(file) {
      req(pipeline_data$supp_results)
      save_current_plot(file, print(plot_supp_top_hubs(pipeline_data$supp_results)), width = 10, height = 8)
    }
  )

  output$dl_supp_heatmap_adni <- downloadHandler(
    filename = "supplementary_heatmap_adni.pdf",
    content = function(file) {
      req(pipeline_data$supp_results)
      save_current_plot(file, plot_supp_cohort_heatmap(pipeline_data$supp_results, "ADNI"), width = 12, height = 8)
    }
  )

  output$dl_supp_heatmap_adrc <- downloadHandler(
    filename = "supplementary_heatmap_adrc.pdf",
    content = function(file) {
      req(pipeline_data$supp_results)
      save_current_plot(file, plot_supp_cohort_heatmap(pipeline_data$supp_results, "ADRC"), width = 12, height = 8)
    }
  )
}

# 4. --- RUN THE APP ---
# -----------------------
shinyApp(ui = ui, server = server)
