#
# Helper functions for the Multi-Omics Shiny App
# This file contains the core analysis logic adapted from the R Markdown scripts.
#


# --- Helper for WGCNA (aligned to 03_WGCNA_running.Rmd) ---

run_wgcna <- function(expr, softPower, minClusterSize, cutHeight, updateProgress = NULL) {
  
  # 1. Adjacency
  if (is.function(updateProgress)) updateProgress("Calculating Adjacency")
  adjacency <- WGCNA::adjacency(expr, power = softPower, type = "signed hybrid")
  
  # 2. TOM
  if (is.function(updateProgress)) updateProgress("Calculating TOM (this may take a while)")
  TOM <- WGCNA::TOMsimilarity(adjacency, TOMType = "signed")
  dissTOM <- 1 - TOM
  
  # 3. Clustering
  if (is.function(updateProgress)) updateProgress("Clustering Genes")
  geneTree <- hclust(as.dist(dissTOM), method = "average")
  
  # 4. Dynamic Modules
  if (is.function(updateProgress)) updateProgress("Detecting Dynamic Modules")
  dynamicMods <- dynamicTreeCut::cutreeDynamic(
    dendro = geneTree, 
    distM = dissTOM,
    deepSplit = 3, 
    pamRespectsDendro = FALSE,
    minClusterSize = minClusterSize
  )
  dynamicColors <- WGCNA::labels2colors(dynamicMods)
  
  # 5. Merging
  if (is.function(updateProgress)) updateProgress("Merging Close Modules")
  merge <- WGCNA::mergeCloseModules(expr, dynamicColors, cutHeight = cutHeight)
  
  return(list(
    geneTree = geneTree,
    dynamicColors = dynamicColors,
    moduleColors = merge$colors,
    MEs = merge$newMEs
  ))
}

# --- Helper for WGCNA QC (aligned to 03_WGCNA_running.Rmd) ---

perform_wgcna_qc <- function(expr, traits = NULL) {
  # 1. Orientation Check & Alignment
  if (!is.null(traits)) {
    # Check if samples are columns (match trait rownames to expr colnames)
    if (!is.null(rownames(traits)) && !is.null(colnames(expr)) && 
        all(rownames(traits) %in% colnames(expr)) && 
        !all(rownames(traits) %in% rownames(expr))) {
      expr <- t(expr)
    }
    
    # Intersect samples
    common <- intersect(rownames(expr), rownames(traits))
    if(length(common) == 0) stop("No common samples between expression and traits.")
    expr <- expr[common, , drop=FALSE]
    traits <- traits[common, , drop=FALSE]
  }
  
  # 2. GoodSamplesGenes Check
  gsg <- WGCNA::goodSamplesGenes(expr, verbose = 0)
  if (!gsg$allOK) {
    expr <- expr[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
    if (!is.null(traits)) {
      traits <- traits[gsg$goodSamples, , drop = FALSE]
    }
  }
  
  return(list(expr = expr, traits = traits, gsg = gsg))
}

# --- Helper for WGCNA Soft Threshold Analysis ---

analyze_soft_threshold <- function(expr, updateProgress = NULL) {
  if (is.function(updateProgress)) updateProgress("Calculating soft-threshold powers...")
  powers <- c(1:20)
  sft <- WGCNA::pickSoftThreshold(expr, powerVector = powers, verbose = 3, networkType = "signed hybrid")
  return(sft)
}

plot_scale_free_topology <- function(sft) {
  if (is.null(sft)) {
    plot.new()
    text(0.5, 0.5, "Run Soft Threshold Analysis first.")
    return()
  }
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=sft$fitIndices[,1],cex=0.9,col="red");
  abline(h=0.85,col="red")
}

plot_mean_connectivity <- function(sft) {
  if (is.null(sft)) {
    plot.new()
    text(0.5, 0.5, "Run Soft Threshold Analysis first.")
    return()
  }
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=sft$fitIndices[,1], cex=0.9,col="red")
}

plot_sample_tree <- function(expr) {
  sampleTree <- hclust(dist(expr), method = "average")
  plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
}

# --- Helper for WGCNA Plotting and Metrics ---

calculate_module_trait_relationships <- function(MEs, traits) {
  # Ensure traits is a data frame
  if (!is.data.frame(traits)) return(NULL)
  
  # Select numeric columns only
  numeric_traits <- traits %>% dplyr::select(where(is.numeric))
  
  # Handle case where no numeric traits
  if (ncol(numeric_traits) == 0) return(NULL)
  
  # Align samples
  common_samples <- intersect(rownames(MEs), rownames(numeric_traits))
  if (length(common_samples) == 0) return(NULL)
  
  MEs_sub <- MEs[common_samples, , drop = FALSE]
  traits_sub <- numeric_traits[common_samples, , drop = FALSE]
  
  moduleTraitCor <- WGCNA::cor(MEs_sub, traits_sub, use = "pairwise.complete.obs")
  moduleTraitP <- WGCNA::corPvalueStudent(moduleTraitCor, length(common_samples))
  
  return(list(cor = moduleTraitCor, pval = moduleTraitP))
}

plot_dendro_and_colors <- function(geneTree, moduleColors) {
  if (is.null(geneTree) || is.null(moduleColors)) {
    plot.new()
    text(0.5, 0.5, "Gene dendrogram or module colors missing.\n(Run WGCNA or load data containing 'geneTree')")
    return()
  }
  WGCNA::plotDendroAndColors(
    geneTree, 
    moduleColors, 
    "Module Colors", 
    dendroLabels = FALSE, 
    hang = 0.03, 
    addGuide = TRUE, 
    guideHang = 0.05,
    main = "Gene Dendrogram and Module Colors"
  )
}

plot_eigengene_network <- function(MEs) {
  if (is.null(MEs)) {
    plot.new()
    text(0.5, 0.5, "Module Eigengenes missing.")
    return()
  }
  # Calculate dissimilarity of module eigengenes
  MEDiss = 1 - cor(MEs, use = "p")
  METree = hclust(as.dist(MEDiss), method = "average")
  
  plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
  abline(h = 0.25, col = "red")
}

plot_module_trait_heatmap <- function(moduleTraitCor, moduleTraitP) {
  if (is.null(moduleTraitCor) || is.null(moduleTraitP)) {
    plot.new()
    text(0.5, 0.5, "Module-Trait correlations missing.\n(Load data with traits or run WGCNA)")
    return()
  }
  
  textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                      signif(moduleTraitP, 1), ")", sep = "")
  dim(textMatrix) <- dim(moduleTraitCor)
  
  par(mar = c(6, 8.5, 3, 3))
  
  WGCNA::labeledHeatmap(
    Matrix = moduleTraitCor,
    xLabels = colnames(moduleTraitCor),
    yLabels = rownames(moduleTraitCor),
    ySymbols = rownames(moduleTraitCor),
    colorLabels = FALSE,
    colors = WGCNA::blueWhiteRed(50),
    textMatrix = textMatrix,
    setStdMargins = FALSE,
    cex.text = 0.5,
    zlim = c(-1, 1),
    main = paste("Module-Trait Relationships")
  )
}

# --- Helper for Module Scoring (aligned to 05_Module_Scoring_Selection.Rmd) ---

run_module_scoring <- function(moduleTraitCor, moduleTraitP, r_thresh, q_thresh) {
  
  # Define trait groups (could also be passed as arguments)
  vascular_traits <- c("PC1", "fsrp", "htnscore", "SBP", "d_cmb", "lacunar_stroke")
  caa_traits <- c("caa_prob", "caa")
  ad_traits <- c("AD_status", "E4status", "cog_status")
  
  # Ensure trait groups only contain columns present in the data
  vascular_traits <- intersect(vascular_traits, colnames(moduleTraitCor))
  caa_traits <- intersect(caa_traits, colnames(moduleTraitCor))
  ad_traits <- intersect(ad_traits, colnames(moduleTraitCor))

  # Calculate FDR-Corrected P-values
  qval_mat <- moduleTraitP
  for (trait in colnames(moduleTraitP)) {
    qval_mat[, trait] <- p.adjust(moduleTraitP[, trait], method = "BH")
  }
  
  # Discrete scoring function
  score_discrete <- function(trait_set, direction = "positive") {
    sapply(rownames(moduleTraitCor), function(mod) {
      r_vec <- as.numeric(moduleTraitCor[mod, trait_set])
      q_vec <- as.numeric(qval_mat[mod, trait_set])
      if (direction == "positive") {
        sum(r_vec >= r_thresh & q_vec < q_thresh, na.rm = TRUE)
      } else {
        sum(r_vec <= -r_thresh & q_vec < q_thresh, na.rm = TRUE)
      }
    })
  }
  
  # Perform scoring
  module_scores <- data.frame(
    module = rownames(moduleTraitCor),
    vascular_hits = score_discrete(vascular_traits, direction = "positive"),
    caa_hits = score_discrete(caa_traits, direction = "positive"),
    ad_hits_pos = score_discrete(ad_traits, direction = "positive"),
    ad_hits_neg = score_discrete(ad_traits, direction = "negative")
  )
  
  # Assign labels
  module_scores$label <- "none"
  module_scores$label[module_scores$vascular_hits >= 2] <- "vascular"
  module_scores$label[module_scores$caa_hits >= 1] <- "CAA"
  module_scores$label[module_scores$ad_hits_neg >= 1] <- "protective_healthy"
  
  return(module_scores)
}


# --- Helper for Machine Learning (aligned to 06_ADNI_ADRC_ML.Rmd) ---

run_ml_analysis <- function(MEs, traits, target_col, split_ratio, alpha) {
  
  # Align samples
  common_samples <- intersect(rownames(MEs), rownames(traits))
  MEs <- MEs[common_samples, , drop = FALSE]
  traits <- traits[common_samples, , drop = FALSE]
  
  # Prepare data frame for modeling
  ml_df <- data.frame(
    MEs,
    outcome = factor(traits[[target_col]])
  )
  ml_df <- ml_df[!is.na(ml_df$outcome), ]
  
  if (length(unique(ml_df$outcome)) < 2) {
    stop("Outcome variable has fewer than 2 classes.")
  }
  
  # Train/Test Split
  set.seed(42)
  train_idx <- caret::createDataPartition(ml_df$outcome, p = split_ratio, list = FALSE)
  train_df <- ml_df[train_idx, ]
  test_df  <- ml_df[-train_idx, ]
  
  # Feature Scaling
  preproc <- caret::preProcess(train_df %>% dplyr::select(-outcome), method = c("center", "scale"))
  X_train <- predict(preproc, train_df %>% dplyr::select(-outcome))
  X_test  <- predict(preproc, test_df %>% dplyr::select(-outcome))
  y_train <- train_df$outcome
  y_test  <- test_df$outcome
  
  # Elastic Net Model
  y01 <- as.numeric(y_train == levels(y_train)[2])
  cv_fit <- glmnet::cv.glmnet(
    x = as.matrix(X_train),
    y = y01,
    family = "binomial",
    alpha = alpha,
    nfolds = 5,
    standardize = FALSE
  )
  
  best_lambda <- cv_fit$lambda.min
  
  # Evaluation
  test_probs <- as.numeric(predict(cv_fit, newx = as.matrix(X_test), s = best_lambda, type = "response"))
  levels_y <- levels(y_test)
  test_pred_class <- ifelse(test_probs > 0.5, levels_y[2], levels_y[1])
  test_pred_factor <- factor(test_pred_class, levels = levels_y)
  
  cm <- caret::confusionMatrix(test_pred_factor, y_test, positive = levels_y[2])
  roc_obj <- pROC::roc(response = y_test, predictor = as.numeric(test_probs), levels = levels_y, quiet = TRUE)
  
  # Coefficients
  coef_mat <- coef(cv_fit, s = best_lambda)
  coef_df <- data.frame(
    Module = rownames(coef_mat),
    Coefficient = as.numeric(coef_mat)
  ) %>%
    dplyr::filter(Module != "(Intercept)", Coefficient != 0) %>%
    dplyr::arrange(desc(abs(Coefficient)))
  
  return(list(
    confusion_matrix = cm,
    roc_object = roc_obj,
    auc_value = pROC::auc(roc_obj),
    coefficients = coef_df
  ))
}

# --- Helper for GSEA (aligned to 07a_Functional_Analysis_GSEA.Rmd) ---

run_gsea_analysis <- function(wgcna_objects, analyte_info_path, target_module, min_size, max_size) {
    # Unpack objects
    expr <- wgcna_objects[["expr"]]
    MEs <- wgcna_objects[["MEs"]]
    moduleColors <- wgcna_objects[["moduleColors"]]
    if(is.null(names(moduleColors))) names(moduleColors) <- colnames(expr)

    # Compute kME (module membership)
    datKME <- as.data.frame(cor(expr, MEs, use = "p"))
    colnames(datKME) <- paste0("kME_", sub("^ME", "", colnames(MEs)))
    rownames(datKME) <- colnames(expr)

    # Annotation Mapping (SOMAmer -> Entrez)
    analyte_info <- resolve_analyte_info(analyte_info_path)
    
    feature_key <- analyte_info %>%
      dplyr::filter(Analytes %in% colnames(expr)) %>%
      dplyr::select(Analytes, EntrezGeneSymbol) %>%
      dplyr::distinct() %>%
      dplyr::rename(SOMAmer = Analytes, SYMBOL = EntrezGeneSymbol)

    symbol_map <- clusterProfiler::bitr(
      feature_key$SYMBOL,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db
    ) %>% dplyr::distinct(SYMBOL, .keep_all = TRUE)

    feature_map <- dplyr::inner_join(feature_key, symbol_map, by = "SYMBOL")

    # Load Gene Sets
    hm <- msigdbr::msigdbr(species = "Homo sapiens", category = "H")
    react <- msigdbr::msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")
    go <- msigdbr::msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")

    gene_sets <- list(
      Hallmark = split(hm$entrez_gene, as.character(hm$gs_name)),
      Reactome = split(react$entrez_gene, as.character(react$gs_name)),
      GO_BP = split(go$entrez_gene, as.character(go$gs_name))
    )

    # GSEA Function to run for a specific database
    run_fgsea_for_db <- function(pathway_sets) {
        kme_col <- paste0("kME_", target_module)
        if (!kme_col %in% colnames(datKME)) return(NULL)

        ranks_df <- data.frame(
            SOMAmer = rownames(datKME),
            kME = datKME[[kme_col]]
        ) %>%
        dplyr::inner_join(feature_map, by = "SOMAmer") %>%
        dplyr::filter(is.finite(kME)) %>%
        dplyr::distinct(ENTREZID, .keep_all = TRUE)

        if (nrow(ranks_df) < 15) return(NULL)

        gene_ranks <- setNames(ranks_df$kME, ranks_df$ENTREZID)
        gene_ranks <- sort(gene_ranks, decreasing = TRUE)

        fgsea::fgseaMultilevel(
            pathways = pathway_sets,
            stats = gene_ranks,
            minSize = min_size,
            maxSize = max_size,
            nproc = 1
        )
    }

    # Run GSEA for each database
    h_res <- run_fgsea_for_db(gene_sets$Hallmark)
    r_res <- run_fgsea_for_db(gene_sets$Reactome)
    g_res <- run_fgsea_for_db(gene_sets$GO_BP)

    # Combine and return results
    combined_results <- dplyr::bind_rows(
        if (!is.null(h_res)) tibble::as_tibble(h_res) %>% dplyr::mutate(Source = "Hallmark"),
        if (!is.null(r_res)) tibble::as_tibble(r_res) %>% dplyr::mutate(Source = "Reactome"),
        if (!is.null(g_res)) tibble::as_tibble(g_res) %>% dplyr::mutate(Source = "GO_BP")
    )

    return(combined_results)
}

plot_gsea_results <- function(sig_gsea_df, module_name) {
    if (nrow(sig_gsea_df) == 0) return(NULL)
    
    top_res <- sig_gsea_df %>% dplyr::slice_max(order_by = abs(NES), n = 15)
    
    ggplot2::ggplot(top_res, ggplot2::aes(x = NES, y = reorder(pathway, NES), size = -log10(padj), color = NES)) +
        ggplot2::geom_point() +
        ggplot2::scale_color_gradient2(low = "blue", mid = "white", high = "red") +
        ggplot2::theme_bw() +
        ggplot2::labs(title = paste("GSEA Results – Module", module_name), y = "Pathway", x = "NES")
}


# --- Helper for ORA (aligned to 07b_Funtional_Analysis_ORA.Rmd) ---

run_ora_analysis <- function(wgcna_objects, analyte_info_path, target_module) {
    # Unpack objects
    expr <- wgcna_objects[["expr"]]
    moduleColors <- wgcna_objects[["moduleColors"]]
    if(is.null(names(moduleColors))) names(moduleColors) <- colnames(expr)

    # --- ID Mapping (similar to GSEA helper) ---
    analyte_info <- resolve_analyte_info(analyte_info_path)
    
    feature_key <- analyte_info %>%
      dplyr::filter(Analytes %in% colnames(expr)) %>%
      dplyr::select(Analytes, EntrezGeneSymbol) %>%
      dplyr::distinct() %>%
      dplyr::rename(SOMAmer = Analytes, SYMBOL = EntrezGeneSymbol)

    symbol_map <- clusterProfiler::bitr(
      feature_key$SYMBOL,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db
    ) %>% dplyr::distinct(SYMBOL, .keep_all = TRUE)

    feature_map <- dplyr::inner_join(feature_key, symbol_map, by = "SYMBOL")
    
    # --- Get Gene Lists ---
    universe_entrez <- unique(feature_map$ENTREZID)
    
    module_somamers <- colnames(expr)[moduleColors == target_module]
    module_entrez <- feature_map %>%
        dplyr::filter(SOMAmer %in% module_somamers) %>%
        dplyr::pull(ENTREZID) %>%
        unique()
        
    if (length(module_entrez) < 10) {
        stop("Not enough mapped genes in module to run ORA (less than 10).")
    }

    # --- Run ORA for each database (qvalueCutoff=1 to get all results before filtering in app) ---
    go_res <- clusterProfiler::enrichGO(
        gene = module_entrez, universe = universe_entrez, OrgDb = org.Hs.eg.db,
        keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 1, readable = TRUE
    )
    
    kegg_res <- clusterProfiler::enrichKEGG(
        gene = module_entrez, universe = universe_entrez, organism = "hsa",
        pAdjustMethod = "BH", qvalueCutoff = 1
    )
    
    reactome_res <- ReactomePA::enrichPathway(
        gene = module_entrez, universe = universe_entrez, organism = "human",
        pAdjustMethod = "BH", qvalueCutoff = 1, readable = TRUE
    )
    
    # --- Combine Results ---
    safe_as_tibble <- function(res, source_name) {
        if (!is.null(res) && nrow(as.data.frame(res)) > 0) {
            tibble::as_tibble(res) %>% dplyr::mutate(Source = source_name)
        } else { NULL }
    }
    
    combined <- dplyr::bind_rows(
        safe_as_tibble(go_res, "GO_BP"),
        safe_as_tibble(kegg_res, "KEGG"),
        safe_as_tibble(reactome_res, "Reactome")
    )
    
    return(combined)
}

plot_ora_results <- function(sig_ora_df, module_name) {
    if (is.null(sig_ora_df) || nrow(sig_ora_df) == 0) {
        return(ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::ggtitle(paste("No significant ORA results for module", module_name)))
    }
    
    top_res <- sig_ora_df %>% dplyr::group_by(Source) %>% dplyr::slice_min(order_by = p.adjust, n = 10, with_ties = FALSE) %>% dplyr::ungroup()
    top_res$GeneRatioNumeric <- sapply(top_res$GeneRatio, function(x) eval(parse(text=x)))

    ggplot2::ggplot(top_res, ggplot2::aes(x = GeneRatioNumeric, y = reorder(Description, GeneRatioNumeric))) +
        ggplot2::geom_point(ggplot2::aes(size = Count, color = p.adjust)) +
        ggplot2::facet_wrap(~Source, scales = "free_y") +
        ggplot2::scale_color_viridis_c(guide = ggplot2::guide_colorbar(reverse = TRUE)) +
        ggplot2::theme_bw(base_size = 11) +
        ggplot2::labs(title = paste("ORA Results – Module", module_name), x = "Gene Ratio", y = "Pathway / Term", color = "Adj. P-value") +
        ggplot2::theme(axis.text.y = ggplot2::element_text(size = 9))
}

# --- Helper for Data Preprocessing (aligned to 00_Data_preprocessing.Rmd) ---

canonical_wmh_vars <- function() {
  c("wmh_juxta", "wmh_frontal", "wmh_pv", "wmh_parietal", "wmh_posterior")
}

compute_wmh_pca <- function(pheno_df) {
  wmh_vars <- canonical_wmh_vars()
  if (!all(wmh_vars %in% colnames(pheno_df))) {
    return(NULL)
  }

  pca_input <- pheno_df[, wmh_vars, drop = FALSE]
  keep <- stats::complete.cases(pca_input)
  if (sum(keep) < 3) {
    return(NULL)
  }

  pca_res <- stats::prcomp(scale(pca_input[keep, , drop = FALSE]), center = TRUE, scale. = TRUE)
  pca_scores <- as.data.frame(pca_res$x[, 1:3, drop = FALSE])
  pca_scores$PC1 <- -pca_scores$PC1
  pca_scores$SampleID <- rownames(pca_input)[keep]
  pca_scores <- pca_scores[, c("SampleID", "PC1", "PC2", "PC3")]

  list(
    pca = pca_res,
    scores = pca_scores,
    vars = wmh_vars
  )
}

build_canonical_pheno <- function(pheno_df) {
  expected_cols <- c(
    "Cohort", "Sex", "Age_at_draw", "PC1", "PC2", "PC3", "SBP", "DM2",
    "HTNscore", "bmi", "fsrp", "d_cmb", "lacunar_stroke", "stroke",
    "caa_diag", "caa", "E4", "closest_cdr", "cdr_sob", "AD_status",
    "cog_status"
  )

  keep <- intersect(expected_cols, colnames(pheno_df))
  pheno_out <- pheno_df[, keep, drop = FALSE]

  if ("Sex" %in% colnames(pheno_out)) {
    pheno_out$Sex <- ifelse(as.character(pheno_out$Sex) == "Male", 1, 0)
  }
  if ("Cohort" %in% colnames(pheno_out)) {
    pheno_out$Cohort_bin <- ifelse(as.character(pheno_out$Cohort) == "ADNI", 1, 0)
  }
  if ("caa_diag" %in% colnames(pheno_out)) {
    pheno_out$caa_prob <- ifelse(as.character(pheno_out$caa_diag) == "probable", 1, 0)
  }
  if ("closest_cdr" %in% colnames(pheno_out)) {
    pheno_out$closest_cdr <- ifelse(as.character(pheno_out$closest_cdr) == "0" | pheno_out$closest_cdr == 0, 0, 1)
  }
  if ("lacunar_stroke" %in% colnames(pheno_out)) {
    pheno_out$lacunar_stroke <- ifelse(as.character(pheno_out$lacunar_stroke) == "0" | pheno_out$lacunar_stroke == 0, 0, 1)
  }
  if ("E4" %in% colnames(pheno_out)) {
    pheno_out$E4status <- ifelse(as.character(pheno_out$E4) == "0" | pheno_out$E4 == 0, 0, 1)
  }

  pheno_out
}

run_preprocessing <- function(proto_df, pheno_df, impute_method = "Hybrid") {
  # 0. Standardize ID columns
  get_id_col <- function(df) {
    if("SampleID" %in% colnames(df)) return("SampleID")
    if("ID" %in% colnames(df)) return("ID")
    return(colnames(df)[1]) # Fallback to first column
  }
  
  proto_id <- get_id_col(proto_df)
  pheno_id <- get_id_col(pheno_df)
  
  # Rename to SampleID for consistency
  colnames(proto_df)[colnames(proto_df) == proto_id] <- "SampleID"
  colnames(pheno_df)[colnames(pheno_df) == pheno_id] <- "SampleID"
  
  # Clean IDs and set rownames
  proto_df$SampleID <- trimws(as.character(proto_df$SampleID))
  pheno_df$SampleID <- trimws(as.character(pheno_df$SampleID))
  
  # Remove duplicates
  proto_df <- proto_df[!duplicated(proto_df$SampleID), ]
  pheno_df <- pheno_df[!duplicated(pheno_df$SampleID), ]
  
  rownames(proto_df) <- proto_df$SampleID
  rownames(pheno_df) <- pheno_df$SampleID

  # 0b. Canonical WMH PCA from phenotype data, matching 02b + 00 merge/flip logic
  wmh_pca <- compute_wmh_pca(pheno_df)
  if (!is.null(wmh_pca)) {
    matched <- match(pheno_df$SampleID, wmh_pca$scores$SampleID)
    for (pc in c("PC1", "PC2", "PC3")) {
      pheno_df[[pc]] <- wmh_pca$scores[[pc]][matched]
    }
  }

  # 0c. Restrict to SOMA panel if available, matching the original workflow
  if ("Panel.SOMA" %in% colnames(proto_df)) {
    proto_df <- proto_df[proto_df$Panel.SOMA == 1, , drop = FALSE]
  }
  if ("Panel.SOMA" %in% colnames(pheno_df)) {
    pheno_df <- pheno_df[pheno_df$Panel.SOMA == 1, , drop = FALSE]
  }
  
  # 1. Prepare Expression Matrix (Logic from 02a_ADNI_ADRC_harmonized.Rmd)
  
  # Remove non-numeric columns specified in script + SampleID
  cols_to_remove <- c("SampleID", "drawDate", "Panel.OLINK")
  cols_to_keep <- setdiff(colnames(proto_df), cols_to_remove)
  expr_data_for_matrix <- proto_df[, cols_to_keep, drop = FALSE]
  
  # Convert to strictly numeric matrix
  # Use suppressWarnings to avoid noise from potential non-numeric coercion
  expr_matrix <- suppressWarnings(as.matrix(sapply(expr_data_for_matrix, as.numeric)))
  rownames(expr_matrix) <- proto_df$SampleID 
  
  # Transpose (Proteins as Rows, Samples as Columns)
  expr_matrix_transposed <- t(expr_matrix)
  
  # Handle NAs: Remove proteins that are all NA
  all_na_rows <- which(rowSums(is.na(expr_matrix_transposed)) == ncol(expr_matrix_transposed))
  if(length(all_na_rows) > 0) {
      expr_matrix_clean <- expr_matrix_transposed[-all_na_rows, ]
      message(paste("Removed", length(all_na_rows), "proteins that were all NA."))
  } else {
      expr_matrix_clean <- expr_matrix_transposed
  }
  
  # Prepare data for imputation (Samples x Proteins)
  # We keep a copy of the raw data (filtered but with NAs) for plotting
  expr_raw_filtered <- t(expr_matrix_clean)
  
  # 2. Imputation
  imputed_data <- expr_raw_filtered # Default if no imputation needed
  
  if (any(is.na(expr_raw_filtered))) {
    if (impute_method == "Hybrid") {
      if (requireNamespace("imputeLCMD", quietly = TRUE)) {
        tryCatch({
          # Hybrid method from 00_Data_preprocessing.Rmd
          # Uses MLE for MAR and QRILC for MNAR
          imputed_data <- imputeLCMD::impute.MAR.MNAR(
            expr_raw_filtered,
            model.selector = 1,
            method.MAR = "MLE",
            method.MNAR = "QRILC"
          )
        }, error = function(e) {
          warning("Hybrid imputation failed: ", e$message, ". Falling back to Minimum.")
          impute_method <<- "Minimum"
        })
      } else {
        warning("Package 'imputeLCMD' not found. Falling back to Minimum.")
        impute_method <<- "Minimum"
      }
    }
    
    if (impute_method == "KNN") {
      if (requireNamespace("impute", quietly = TRUE)) {
        tryCatch({
          # impute.knn expects Features x Samples
          imputed_res <- impute::impute.knn(t(expr_raw_filtered))
          imputed_data <- t(imputed_res$data)
        }, error = function(e) {
          warning("KNN imputation failed. Falling back to Minimum.")
          impute_method <<- "Minimum"
        })
      }
    }
    
    if (impute_method == "Minimum") {
      # Replace NA with column minimum (since we are in Samples x Proteins orientation)
      imputed_data <- apply(expr_raw_filtered, 2, function(x) {
        x[is.na(x)] <- min(x, na.rm = TRUE)
        return(x)
      })
    }
  }
  
  # Transpose back to Samples x Proteins for the app pipeline
  processed_proto <- as.data.frame(imputed_data)
  
  # 3. Align Phenotype
  common_ids <- intersect(rownames(processed_proto), rownames(pheno_df))
  if (length(common_ids) == 0) stop("No common sample IDs between proteomics and phenotype data.")
  
  processed_proto <- processed_proto[common_ids, ]
  expr_raw_filtered <- expr_raw_filtered[common_ids, ]
  pheno_aligned <- pheno_df[common_ids, , drop = FALSE]
  pheno_model <- build_canonical_pheno(pheno_aligned)
  
  return(list(
    expr = processed_proto,
    expr_raw = expr_raw_filtered,
    pheno = pheno_model,
    pheno_raw = pheno_aligned,
    wmh_pca = wmh_pca
  ))
}

plot_missingness <- function(df) {
  is_num <- sapply(df, is.numeric)
  numeric_data <- df[, is_num, drop = FALSE]
  missing_counts <- colSums(is.na(numeric_data))
  
  df_plot <- data.frame(Protein = names(missing_counts), MissingCount = missing_counts)
  
  ggplot2::ggplot(df_plot, ggplot2::aes(x = MissingCount)) +
    ggplot2::geom_histogram(binwidth = 1, fill = "steelblue", color = "white") +
    ggplot2::theme_bw() +
    ggplot2::labs(title = "Distribution of Missing Values per Protein", x = "Number of Missing Samples", y = "Count of Proteins")
}

plot_imputation_distribution <- function(expr_raw, expr_imputed) {
  # Flatten matrices to vectors
  raw_vals <- as.vector(as.matrix(expr_raw))
  raw_vals <- raw_vals[!is.na(raw_vals)] # Remove NAs for plotting
  
  imp_vals <- as.vector(as.matrix(expr_imputed))
  
  df_plot <- data.frame(
    Value = c(raw_vals, imp_vals),
    Type = c(rep("Before Imputation", length(raw_vals)), rep("After Imputation", length(imp_vals)))
  )
  
  ggplot2::ggplot(df_plot, ggplot2::aes(x = Value, fill = Type)) +
    ggplot2::geom_density(alpha = 0.5) +
    ggplot2::theme_bw() +
    ggplot2::labs(title = "Data Distribution Before and After Imputation", x = "Intensity", y = "Density") +
    ggplot2::scale_fill_manual(values = c("Before Imputation" = "grey", "After Imputation" = "red"))
}

# --- Helper for Harmonization (ComBat) ---

run_harmonization <- function(expr_data, pheno_data, batch_col) {
  if (!requireNamespace("sva", quietly = TRUE)) stop("sva package required for ComBat")
  
  # Critical Alignment (ensure samples match exactly)
  common_samples <- intersect(rownames(expr_data), rownames(pheno_data))
  if (length(common_samples) == 0) stop("No common samples found for harmonization.")
  
  expr_data <- expr_data[common_samples, , drop = FALSE]
  pheno_data <- pheno_data[common_samples, , drop = FALSE]
  
  # Ensure exact order match
  match_idx <- match(rownames(expr_data), rownames(pheno_data))
  pheno_data <- pheno_data[match_idx, , drop = FALSE]
  
  # Ensure expr_data is numeric matrix (samples x features)
  # ComBat expects features x samples
  is_num <- sapply(expr_data, is.numeric)
  dat <- t(as.matrix(expr_data[, is_num]))
  
  # Check for NAs (ComBat doesn't like them)
  if (any(is.na(dat))) stop("Expression data contains NAs. Run preprocessing first.")
  
  # Batch vector
  if (!batch_col %in% colnames(pheno_data)) stop(paste("Batch column", batch_col, "not found in phenotype data"))
  batch <- pheno_data[[batch_col]]
  
  if(length(unique(batch)) < 2) stop("Batch variable must have at least 2 levels.")
  if(any(is.na(batch))) stop("Batch variable contains NAs.")
  
  # Run ComBat
  mod <- model.matrix(~1, data = pheno_data)
  pca_before <- stats::prcomp(expr_data[, is_num, drop = FALSE], scale. = TRUE)
  combat_edata <- sva::ComBat(dat = dat, batch = batch, mod = mod, par.prior = TRUE, prior.plots = FALSE)
  harmonized_expr <- t(combat_edata)
  pca_after <- stats::prcomp(harmonized_expr, scale. = TRUE)
  
  list(
    expr = harmonized_expr,
    pheno = pheno_data,
    pca_before = pca_before,
    pca_after = pca_after
  )
}

plot_pca_batch <- function(expr_data, pheno_data, batch_col, title_prefix = "PCA") {
  # Align pheno_data to expr_data
  common <- intersect(rownames(expr_data), rownames(pheno_data))
  if(length(common) == 0) return(ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::ggtitle("No common samples"))
  
  expr_data <- expr_data[common, , drop=FALSE]
  pheno_data <- pheno_data[common, , drop=FALSE]
  
  # Ensure exact order match (redundant but safe)
  match_idx <- match(rownames(expr_data), rownames(pheno_data))
  pheno_data <- pheno_data[match_idx, , drop=FALSE]
  
  if (!batch_col %in% colnames(pheno_data)) return(ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::ggtitle("Batch column not found"))
  
  # Simple PCA
  is_num <- sapply(expr_data, is.numeric)
  # Remove zero variance columns
  expr_data_num <- expr_data[, is_num, drop=FALSE]
  expr_data_num <- expr_data_num[, apply(expr_data_num, 2, var) > 0, drop=FALSE]
  
  if(ncol(expr_data_num) < 2) return(ggplot2::ggplot() + ggplot2::theme_void() + ggplot2::ggtitle("Not enough numeric features"))

  pca <- prcomp(expr_data_num, scale. = TRUE)
  
  df_plot <- data.frame(
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    Batch = as.factor(pheno_data[[batch_col]])
  )
  
  ggplot2::ggplot(df_plot, ggplot2::aes(x = PC1, y = PC2, color = Batch)) +
    ggplot2::geom_point(alpha = 0.7, size = 2) +
    ggplot2::theme_bw() +
    ggplot2::labs(title = paste(title_prefix, "- Colored by", batch_col), color = batch_col)
}

# --- Helper for PCA (Full Analysis) ---

run_pca_analysis <- function(expr_data) {
  is_num <- sapply(expr_data, is.numeric)
  pca <- prcomp(expr_data[, is_num], scale. = TRUE)
  return(pca)
}

resolve_analyte_info <- function(analyte_info_input) {
  if (is.data.frame(analyte_info_input)) {
    return(analyte_info_input)
  }
  read.delim(analyte_info_input, sep = "\t", header = TRUE, check.names = FALSE)
}

plot_pca_scree <- function(pca_obj) {
  var_explained <- pca_obj$sdev^2 / sum(pca_obj$sdev^2)
  df_plot <- data.frame(PC = 1:length(var_explained), Variance = var_explained)
  
  # Top 10 PCs
  df_plot <- head(df_plot, 10)
  
  ggplot2::ggplot(df_plot, ggplot2::aes(x = factor(PC), y = Variance)) +
    ggplot2::geom_bar(stat = "identity", fill = "steelblue") +
    ggplot2::geom_line(ggplot2::aes(group = 1), color = "red") +
    ggplot2::geom_point(color = "red") +
    ggplot2::theme_bw() +
    ggplot2::labs(title = "Scree Plot (Top 10 PCs)", x = "Principal Component", y = "Proportion of Variance Explained")
}

plot_pca_scatter_custom <- function(pca_obj, pheno_data, color_col) {
  df_plot <- data.frame(
    PC1 = pca_obj$x[, 1],
    PC2 = pca_obj$x[, 2]
  )
  if (!is.null(color_col) && color_col %in% colnames(pheno_data)) {
    df_plot$Color <- as.factor(pheno_data[[color_col]])
    p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = PC1, y = PC2, color = Color)) + ggplot2::labs(color = color_col)
  } else {
    p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = PC1, y = PC2))
  }
  p + ggplot2::geom_point(alpha = 0.7, size = 2) + ggplot2::theme_bw() + ggplot2::labs(title = "PCA: PC1 vs PC2")
}

plot_pca_biplot <- function(pca_obj, title_prefix = "PCA Biplot") {
  # Extract scores
  pca_scores <- as.data.frame(pca_obj$x)
  
  # Extract loadings
  loadings_raw <- as.data.frame(pca_obj$rotation)
  loadings_raw$Trait <- rownames(loadings_raw)
  
  # Scaling factor for arrows (simple heuristic to match score scale)
  r <- min(
    (max(pca_scores$PC1) - min(pca_scores$PC1)) / (max(loadings_raw$PC1) - min(loadings_raw$PC1)),
    (max(pca_scores$PC2) - min(pca_scores$PC2)) / (max(loadings_raw$PC2) - min(loadings_raw$PC2))
  )
  arrow_scale <- r * 0.8
  
  loadings_plot <- loadings_raw
  loadings_plot[, c("PC1", "PC2")] <- loadings_plot[, c("PC1", "PC2")] * arrow_scale
  
  # Variance explained
  var_expl <- summary(pca_obj)$importance[2, 1:2] * 100
  
  p <- ggplot2::ggplot(pca_scores, ggplot2::aes(x = PC1, y = PC2)) +
    ggplot2::geom_point(alpha = 0.5, color = "grey50") +
    ggplot2::geom_segment(
      data = loadings_plot,
      ggplot2::aes(x = 0, y = 0, xend = PC1, yend = PC2),
      arrow = ggplot2::arrow(length = ggplot2::unit(0.2, "cm")),
      color = "red"
    ) +
    ggrepel::geom_text_repel(
      data = loadings_plot,
      ggplot2::aes(x = PC1, y = PC2, label = Trait),
      color = "red",
      size = 4
    ) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = title_prefix,
      x = paste0("PC1 (", round(var_expl[1], 1), "%)"),
      y = paste0("PC2 (", round(var_expl[2], 1), "%)")
    )
  return(p)
}

plot_pca_loadings_heatmap <- function(pca_obj) {
  # Extract loadings for first 3 PCs
  loadings <- pca_obj$rotation[, 1:min(3, ncol(pca_obj$rotation)), drop = FALSE]
  
  # Plot using pheatmap
  pheatmap::pheatmap(
    loadings,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    display_numbers = TRUE,
    number_format = "%.2f",
    color = colorRampPalette(c("blue", "white", "red"))(100),
    main = "PCA Loadings (PC1-PC3)",
    silent = TRUE
  )$gtable
}
