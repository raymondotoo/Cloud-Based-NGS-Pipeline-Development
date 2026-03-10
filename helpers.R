#
# Helper functions for the Multi-Omics Shiny App
# This file contains the core analysis logic adapted from the R Markdown scripts.
#

app_theme <- function(base_size = 12) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_line(color = "#e2e8f0"),
      panel.grid.major.y = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold"),
      plot.subtitle = ggplot2::element_text(color = "#475569"),
      strip.background = ggplot2::element_rect(fill = "#e2e8f0", color = NA),
      legend.title = ggplot2::element_text(face = "bold")
    )
}


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
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  par(mar = c(5, 5, 3, 1))
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=sft$fitIndices[,1],cex=0.9,col="#9f1239");
  abline(h=0.85,col="#9f1239", lty = 2)
}

plot_mean_connectivity <- function(sft) {
  if (is.null(sft)) {
    plot.new()
    text(0.5, 0.5, "Run Soft Threshold Analysis first.")
    return()
  }
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  par(mar = c(5, 5, 3, 1))
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=sft$fitIndices[,1], cex=0.9,col="#0f766e")
}

plot_sample_tree <- function(expr) {
  sampleTree <- hclust(dist(expr), method = "average")
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  par(mar = c(5, 4, 3, 1))
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

# --- Shared helpers for Functional Analysis ---

compute_datkme <- function(expr, MEs) {
  datKME <- as.data.frame(stats::cor(expr, MEs, use = "pairwise.complete.obs"))
  colnames(datKME) <- paste0("kME_", sub("^ME", "", colnames(MEs)))
  rownames(datKME) <- colnames(expr)
  datKME
}

build_feature_map <- function(expr, analyte_info_input) {
  analyte_info <- resolve_analyte_info(analyte_info_input)

  feature_key <- analyte_info %>%
    dplyr::filter(Analytes %in% colnames(expr)) %>%
    dplyr::select(Analytes, EntrezGeneSymbol) %>%
    dplyr::filter(!is.na(EntrezGeneSymbol), EntrezGeneSymbol != "") %>%
    dplyr::distinct() %>%
    dplyr::rename(SOMAmer = Analytes, SYMBOL = EntrezGeneSymbol)

  symbol_map <- clusterProfiler::bitr(
    feature_key$SYMBOL,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  ) %>%
    dplyr::distinct(SYMBOL, .keep_all = TRUE)

  dplyr::inner_join(feature_key, symbol_map, by = "SYMBOL")
}

clean_pathway_label <- function(x, width = 42) {
  stringr::str_wrap(
    gsub("_", " ", gsub("^(GO_BP_|HALLMARK_|REACTOME_)", "", x)),
    width = width
  )
}

empty_plot <- function(title, subtitle = NULL) {
  ggplot2::ggplot() +
    ggplot2::theme_void() +
    ggplot2::labs(title = title, subtitle = subtitle)
}

# --- Helper for GSEA (aligned to 07a_Functional_Analysis_GSEA.Rmd) ---

run_gsea_analysis <- function(wgcna_objects, analyte_info_path, target_module, min_size, max_size) {
  expr <- wgcna_objects[["expr"]]
  MEs <- wgcna_objects[["MEs"]]
  moduleColors <- wgcna_objects[["moduleColors"]]
  if (is.null(names(moduleColors))) names(moduleColors) <- colnames(expr)

  datKME <- compute_datkme(expr, MEs)
  feature_map <- build_feature_map(expr, analyte_info_path)
  kme_col <- paste0("kME_", target_module)
  if (!kme_col %in% colnames(datKME)) {
    stop("Selected module is not available in the kME matrix.")
  }

  hm <- msigdbr::msigdbr(species = "Homo sapiens", category = "H")
  react <- msigdbr::msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")
  go <- msigdbr::msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")

  gene_sets <- list(
    Hallmark = split(hm$entrez_gene, as.character(hm$gs_name)),
    Reactome = split(react$entrez_gene, as.character(react$gs_name)),
    GO_BP = split(go$entrez_gene, as.character(go$gs_name))
  )

  module_genes <- names(moduleColors[moduleColors == target_module])
  ranks_df <- data.frame(
    SOMAmer = rownames(datKME),
    kME = datKME[[kme_col]],
    stringsAsFactors = FALSE
  ) %>%
    dplyr::filter(SOMAmer %in% module_genes) %>%
    dplyr::inner_join(feature_map, by = "SOMAmer") %>%
    dplyr::filter(is.finite(kME), !is.na(ENTREZID)) %>%
    dplyr::distinct(ENTREZID, .keep_all = TRUE)

  if (nrow(ranks_df) < 15) {
    stop("Too few mapped genes in the selected module to run GSEA.")
  }

  gene_ranks <- setNames(ranks_df$kME, ranks_df$ENTREZID)
  gene_ranks <- sort(gene_ranks, decreasing = TRUE)

  run_fgsea_for_db <- function(pathway_sets, source_name) {
    res <- fgsea::fgseaMultilevel(
      pathways = pathway_sets,
      stats = gene_ranks,
      minSize = min_size,
      maxSize = max_size,
      scoreType = "pos",
      nproc = 1
    )

    tibble::as_tibble(res) %>%
      dplyr::mutate(
        Source = source_name,
        pathway_clean = clean_pathway_label(pathway)
      ) %>%
      dplyr::arrange(padj)
  }

  results_by_source <- list(
    Hallmark = run_fgsea_for_db(gene_sets$Hallmark, "Hallmark"),
    Reactome = run_fgsea_for_db(gene_sets$Reactome, "Reactome"),
    GO_BP = run_fgsea_for_db(gene_sets$GO_BP, "GO_BP")
  )

  combined_results <- dplyr::bind_rows(results_by_source)

  list(
    module = target_module,
    combined = combined_results,
    results_by_source = results_by_source,
    feature_map = feature_map,
    datKME = datKME,
    gene_ranks = gene_ranks,
    gene_sets = gene_sets,
    ranks_df = ranks_df
  )
}

get_sig_gsea_results <- function(gsea_bundle, padj_cutoff = 0.05, source_name = NULL) {
  if (is.null(gsea_bundle) || is.null(gsea_bundle$combined)) {
    return(tibble::tibble())
  }

  res <- gsea_bundle$combined %>%
    dplyr::filter(!is.na(padj), padj <= padj_cutoff)

  if (!is.null(source_name) && nzchar(source_name) && source_name != "All") {
    res <- res %>% dplyr::filter(Source == source_name)
  }

  res %>% dplyr::arrange(Source, padj)
}

plot_gsea_results <- function(gsea_bundle, padj_cutoff = 0.05, source_name = "All", top_n = 18) {
  sig_gsea_df <- get_sig_gsea_results(gsea_bundle, padj_cutoff, source_name)
  if (nrow(sig_gsea_df) == 0) {
    return(empty_plot(
      title = paste("No significant GSEA pathways for module", gsea_bundle$module),
      subtitle = "Try a less stringent adjusted p-value cutoff."
    ))
  }

  top_res <- sig_gsea_df %>%
    dplyr::group_by(Source) %>%
    dplyr::slice_max(order_by = abs(NES), n = top_n, with_ties = FALSE) %>%
    dplyr::ungroup()

  top_res$pathway_clean <- reorder(top_res$pathway_clean, top_res$NES)

  ggplot2::ggplot(
    top_res,
    ggplot2::aes(
      x = NES,
      y = pathway_clean,
      fill = NES
    )
  ) +
    ggplot2::geom_col(width = 0.8) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("padj %.3g", padj)),
      hjust = ifelse(top_res$NES >= 0, -0.1, 1.1),
      size = 3.2,
      color = "#334155"
    ) +
    ggplot2::facet_wrap(~Source, scales = "free_y") +
    ggplot2::scale_fill_gradient2(low = "#155e75", mid = "#e2e8f0", high = "#9f1239") +
    app_theme() +
    ggplot2::labs(
      title = paste("FGSEA pathway overview:", gsea_bundle$module, "module"),
      x = "Normalized enrichment score",
      y = NULL,
      fill = "NES"
    ) +
    ggplot2::coord_cartesian(clip = "off")
}

plot_gsea_network <- function(gsea_bundle, source_name, padj_cutoff = 0.05, n_pathways = 5) {
  sig_res <- get_sig_gsea_results(gsea_bundle, padj_cutoff, source_name) %>%
    dplyr::slice_max(order_by = abs(NES), n = n_pathways, with_ties = FALSE)

  if (nrow(sig_res) == 0) {
    return(empty_plot(
      title = paste("No network plot available for", source_name),
      subtitle = "No significant leading-edge pathways matched the current cutoff."
    ))
  }

  id_to_symbol <- gsea_bundle$feature_map$SYMBOL
  names(id_to_symbol) <- gsea_bundle$feature_map$ENTREZID

  gene_list_symbols <- lapply(sig_res$leadingEdge, function(x) {
    symbols <- unname(id_to_symbol[as.character(x)])
    unique(symbols[!is.na(symbols)])
  })
  names(gene_list_symbols) <- clean_pathway_label(sig_res$pathway, width = 20)

  enrichplot::cnetplot(
    gene_list_symbols,
    showCategory = min(n_pathways, length(gene_list_symbols)),
    layout = "circle",
    cex_label_gene = 0.7,
    cex_label_category = 0.9,
    colorEdge = TRUE,
    color_edge = "category",
    shadowtext = "gene"
  ) +
    ggplot2::labs(
      title = paste("Leading-edge network:", gsea_bundle$module, "module"),
      subtitle = source_name
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5),
      plot.margin = ggplot2::margin(20, 20, 20, 20)
    )
}

plot_gsea_enrichment_curve <- function(gsea_bundle, source_name, pathway_name = NULL, padj_cutoff = 0.05) {
  sig_res <- get_sig_gsea_results(gsea_bundle, padj_cutoff, source_name)
  if (nrow(sig_res) == 0) {
    return(empty_plot(
      title = paste("No enrichment curve available for", source_name),
      subtitle = "No significant pathways matched the current cutoff."
    ))
  }

  if (is.null(pathway_name) || !pathway_name %in% sig_res$pathway) {
    pathway_name <- sig_res$pathway[[1]]
  }

  current_sets <- gsea_bundle$gene_sets[[source_name]]
  pathway_stats <- sig_res %>% dplyr::filter(pathway == pathway_name) %>% dplyr::slice(1)

  fgsea::plotEnrichment(current_sets[[pathway_name]], gsea_bundle$gene_ranks) +
    ggplot2::labs(
      title = paste("Enrichment curve:", gsea_bundle$module, "module"),
      subtitle = clean_pathway_label(pathway_name, width = 60),
      caption = paste0(
        "NES = ", round(pathway_stats$NES, 2),
        " | padj = ", format.pval(pathway_stats$padj, digits = 3)
      )
    ) +
    app_theme() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5),
      plot.caption = ggplot2::element_text(hjust = 0.5, color = "#7f1d1d")
    )
}

# --- Helper for ORA (aligned to 07b_Funtional_Analysis_ORA.Rmd) ---

run_ora_analysis <- function(wgcna_objects, analyte_info_path, target_module) {
  expr <- wgcna_objects[["expr"]]
  MEs <- wgcna_objects[["MEs"]]
  moduleColors <- wgcna_objects[["moduleColors"]]
  if (is.null(names(moduleColors))) names(moduleColors) <- colnames(expr)

  feature_map <- build_feature_map(expr, analyte_info_path)
  universe_entrez <- unique(feature_map$ENTREZID)

  module_somamers <- colnames(expr)[moduleColors == target_module]
  module_entrez <- feature_map %>%
    dplyr::filter(SOMAmer %in% module_somamers) %>%
    dplyr::pull(ENTREZID) %>%
    unique()

  if (length(module_entrez) < 10) {
    stop("Not enough mapped genes in module to run ORA (less than 10).")
  }

  datKME <- compute_datkme(expr, MEs)

  go_res <- clusterProfiler::enrichGO(
    gene = module_entrez,
    universe = universe_entrez,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    qvalueCutoff = 1,
    readable = TRUE
  )

  kegg_res <- clusterProfiler::enrichKEGG(
    gene = module_entrez,
    universe = universe_entrez,
    organism = "hsa",
    pAdjustMethod = "BH",
    qvalueCutoff = 1
  )

  reactome_res <- ReactomePA::enrichPathway(
    gene = module_entrez,
    universe = universe_entrez,
    organism = "human",
    pAdjustMethod = "BH",
    qvalueCutoff = 1,
    readable = TRUE
  )

  safe_as_tibble <- function(res, source_name) {
    if (is.null(res) || nrow(as.data.frame(res)) == 0) {
      return(NULL)
    }

    tibble::as_tibble(res) %>%
      dplyr::mutate(
        Source = source_name,
        pathway_clean = stringr::str_wrap(Description, width = 50)
      )
  }

  combined <- dplyr::bind_rows(
    safe_as_tibble(go_res, "GO_BP"),
    safe_as_tibble(kegg_res, "KEGG"),
    safe_as_tibble(reactome_res, "Reactome")
  )

  list(
    module = target_module,
    combined = combined,
    enrich_objects = list(GO_BP = go_res, KEGG = kegg_res, Reactome = reactome_res),
    feature_map = feature_map,
    datKME = datKME,
    module_somamers = module_somamers
  )
}

get_sig_ora_results <- function(ora_bundle, padj_cutoff = 0.05, source_name = NULL) {
  if (is.null(ora_bundle) || is.null(ora_bundle$combined)) {
    return(tibble::tibble())
  }

  res <- ora_bundle$combined %>%
    dplyr::filter(!is.na(p.adjust), p.adjust <= padj_cutoff)

  if (!is.null(source_name) && nzchar(source_name) && source_name != "All") {
    res <- res %>% dplyr::filter(Source == source_name)
  }

  res %>% dplyr::arrange(Source, p.adjust)
}

plot_ora_results <- function(ora_bundle, padj_cutoff = 0.05, source_name = "All", top_n = 10) {
  sig_ora_df <- get_sig_ora_results(ora_bundle, padj_cutoff, source_name)
  if (nrow(sig_ora_df) == 0) {
    return(empty_plot(
      title = paste("No significant ORA terms for module", ora_bundle$module),
      subtitle = "Try a less stringent adjusted p-value cutoff."
    ))
  }

  top_res <- sig_ora_df %>%
    dplyr::group_by(Source) %>%
    dplyr::slice_min(order_by = p.adjust, n = top_n, with_ties = FALSE) %>%
    dplyr::ungroup()

  top_res$GeneRatioNumeric <- vapply(top_res$GeneRatio, function(x) eval(parse(text = x)), numeric(1))

  ggplot2::ggplot(
    top_res,
    ggplot2::aes(
      x = GeneRatioNumeric,
      y = reorder(pathway_clean, GeneRatioNumeric),
      fill = -log10(p.adjust)
    )
  ) +
    ggplot2::geom_col(width = 0.8) +
    ggplot2::geom_text(
      ggplot2::aes(label = paste0("n=", Count)),
      hjust = -0.1,
      size = 3.2,
      color = "#334155"
    ) +
    ggplot2::facet_wrap(~Source, scales = "free_y") +
    ggplot2::scale_fill_gradient(low = "#fbbf24", high = "#9f1239") +
    app_theme() +
    ggplot2::labs(
      title = paste("ORA term overview:", ora_bundle$module, "module"),
      x = "Gene ratio",
      y = NULL,
      fill = "-log10(adj. p)"
    ) +
    ggplot2::coord_cartesian(clip = "off")
}

plot_ora_cnet <- function(ora_bundle, source_name, padj_cutoff = 0.05, top_n = 5) {
  enrich_obj <- ora_bundle$enrich_objects[[source_name]]
  if (is.null(enrich_obj) || nrow(as.data.frame(enrich_obj)) == 0) {
    return(empty_plot(title = paste("No cnet plot available for", source_name)))
  }

  df <- as.data.frame(enrich_obj)
  df <- df[df$p.adjust <= padj_cutoff, , drop = FALSE]
  if (nrow(df) == 0) {
    return(empty_plot(
      title = paste("No cnet plot available for", source_name),
      subtitle = "No significant enriched terms matched the current cutoff."
    ))
  }

  enrich_obj@result <- enrich_obj@result[enrich_obj@result$p.adjust <= padj_cutoff, , drop = FALSE]

  enrichplot::cnetplot(
    enrich_obj,
    showCategory = min(top_n, nrow(enrich_obj@result)),
    circular = FALSE,
    colorEdge = TRUE
  ) +
    ggplot2::labs(
      title = paste("ORA concept network:", ora_bundle$module, "module"),
      subtitle = source_name
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(hjust = 0.5)
    )
}

extract_ora_hub_genes <- function(ora_bundle, source_name = "GO_BP", padj_cutoff = 0.05, kME_threshold = 0.7) {
  ora_res <- ora_bundle$enrich_objects[[source_name]]
  if (is.null(ora_res) || nrow(as.data.frame(ora_res)) == 0) {
    return(tibble::tibble())
  }

  df <- as.data.frame(ora_res) %>%
    dplyr::filter(p.adjust <= padj_cutoff)

  if (nrow(df) == 0) {
    return(tibble::tibble())
  }

  kme_col <- paste0("kME_", ora_bundle$module)
  module_genes <- ora_bundle$feature_map %>%
    dplyr::filter(SOMAmer %in% ora_bundle$module_somamers) %>%
    dplyr::left_join(
      data.frame(
        SOMAmer = rownames(ora_bundle$datKME),
        kME = ora_bundle$datKME[, kme_col],
        stringsAsFactors = FALSE
      ),
      by = "SOMAmer"
    ) %>%
    dplyr::filter(kME >= kME_threshold)

  hub_rows <- lapply(seq_len(nrow(df)), function(i) {
    pathway_genes <- unlist(strsplit(df$geneID[i], "/"))
    hubs <- module_genes %>%
      dplyr::filter(SYMBOL %in% pathway_genes) %>%
      dplyr::arrange(dplyr::desc(kME))

    if (nrow(hubs) == 0) {
      return(NULL)
    }

    hubs %>%
      dplyr::mutate(
        Pathway = df$Description[i],
        Source = source_name
      ) %>%
      dplyr::select(Source, Pathway, SOMAmer, SYMBOL, ENTREZID, kME)
  })

  dplyr::bind_rows(hub_rows)
}

# --- Helper for Supplementary Analysis (aligned to 08_Supplementary_analysis_v3.Rmd) ---

infer_cohort_labels <- function(sample_ids, traits = NULL) {
  if (!is.null(traits) && "Cohort" %in% colnames(traits)) {
    cohort <- as.character(traits$Cohort)
    names(cohort) <- rownames(traits)
    cohort <- cohort[sample_ids]
    cohort[cohort %in% c("MAP", "ADRC")] <- "ADRC"
    return(cohort)
  }

  ifelse(
    grepl("^ADNI", sample_ids), "ADNI",
    ifelse(grepl("^MAP|^ADRC", sample_ids), "ADRC", NA_character_)
  )
}

prepare_preservation_set <- function(expr_subset, module_colors) {
  if (is.null(expr_subset) || nrow(expr_subset) < 3 || ncol(expr_subset) < 10) {
    return(NULL)
  }

  gsg <- WGCNA::goodSamplesGenes(expr_subset, verbose = 0)
  expr_clean <- expr_subset
  colors_clean <- module_colors

  if (!gsg$allOK) {
    expr_clean <- expr_subset[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
    colors_clean <- module_colors[gsg$goodGenes]
  }

  keep_var <- apply(expr_clean, 2, stats::var, na.rm = TRUE) > 0
  keep_var[is.na(keep_var)] <- FALSE
  expr_clean <- expr_clean[, keep_var, drop = FALSE]
  colors_clean <- colors_clean[keep_var]

  if (nrow(expr_clean) < 3 || ncol(expr_clean) < 10) {
    return(NULL)
  }

  list(expr = expr_clean, colors = as.character(colors_clean))
}

run_supplementary_analysis <- function(wgcna_objects, selected_modules = NULL, n_permutations = 100) {
  expr <- wgcna_objects[["expr"]]
  MEs <- wgcna_objects[["MEs"]]
  moduleColors <- wgcna_objects[["moduleColors"]]
  traits <- wgcna_objects[["traits"]]
  if (is.null(names(moduleColors))) names(moduleColors) <- colnames(expr)

  if (is.null(selected_modules) || length(selected_modules) == 0) {
    selected_modules <- intersect(c("red", "lightyellow", "black", "magenta"), unique(moduleColors))
  }
  selected_modules <- unique(selected_modules)

  datKME <- WGCNA::signedKME(expr, MEs)
  final_protein_kME <- data.frame(
    ProteinID = colnames(expr),
    ModuleColor = moduleColors,
    datKME,
    check.names = FALSE
  )

  sample_ids <- rownames(expr)
  cohort <- infer_cohort_labels(sample_ids, traits)
  valid_cohort <- !is.na(cohort)
  preservation_summary <- NULL
  preservation_raw <- NULL

  if (sum(cohort == "ADNI", na.rm = TRUE) > 10 && sum(cohort == "ADRC", na.rm = TRUE) > 10) {
    adni_set <- prepare_preservation_set(
      expr[cohort == "ADNI" & valid_cohort, , drop = FALSE],
      moduleColors
    )
    adrc_set <- prepare_preservation_set(
      expr[cohort == "ADRC" & valid_cohort, , drop = FALSE],
      moduleColors
    )

    if (!is.null(adni_set) && !is.null(adrc_set)) {
      common_genes <- intersect(colnames(adni_set$expr), colnames(adrc_set$expr))

      if (length(common_genes) >= 10) {
        adni_expr <- adni_set$expr[, common_genes, drop = FALSE]
        adrc_expr <- adrc_set$expr[, common_genes, drop = FALSE]
        color_idx <- match(common_genes, names(moduleColors))
        common_colors <- as.character(moduleColors[color_idx])

        multiExpr <- list(
          ADNI = list(data = adni_expr),
          ADRC = list(data = adrc_expr)
        )
        multiColor <- list(
          ADNI = common_colors,
          ADRC = common_colors
        )

        pres <- WGCNA::modulePreservation(
          multiExpr,
          multiColor = multiColor,
          referenceNetworks = 1,
          nPermutations = n_permutations,
          networkType = "signed",
          randomSeed = 123,
          verbose = 0
        )

        adrc_slot <- grep("ADRC", names(pres$preservation$Z$ref.ADNI), value = TRUE)[1]
        if (!is.na(adrc_slot) && length(selected_modules) > 0) {
          keep_modules <- intersect(selected_modules, rownames(pres$preservation$Z$ref.ADNI[[adrc_slot]]))
          if (length(keep_modules) > 0) {
            z_stats <- pres$preservation$Z$ref.ADNI[[adrc_slot]][keep_modules, , drop = FALSE]
            obs_stats <- pres$preservation$observed$ref.ADNI[[adrc_slot]][keep_modules, , drop = FALSE]
            preservation_summary <- data.frame(
              Module = rownames(z_stats),
              Size = z_stats$moduleSize,
              Z_ADRC = z_stats$Zsummary.pres,
              medianRank_ADRC = obs_stats$medianRank.pres,
              stringsAsFactors = FALSE
            ) %>%
              dplyr::mutate(
                Preservation_Status = dplyr::case_when(
                  Z_ADRC > 10 ~ "Strongly Preserved",
                  Z_ADRC > 2 ~ "Moderately Preserved",
                  TRUE ~ "Weak/No Preservation"
                )
              ) %>%
              dplyr::arrange(medianRank_ADRC)

            preservation_raw <- pres
          }
        }
      }
    }
  }

  analysis_map <- list(
    list(m = "MEred", t = "fsrp"),
    list(m = "MEred", t = "htnscore"),
    list(m = "MElightyellow", t = "PC1"),
    list(m = "MEblack", t = "caa"),
    list(m = "MEmagenta", t = "bmi")
  )

  run_adjusted_model <- function(module_name, trait_name) {
    req_cols <- c(module_name, trait_name, "Age_at_draw", "Sex", "AD_status")
    if (!all(req_cols %in% c(colnames(MEs), colnames(traits)))) {
      return(NULL)
    }

    df_reg <- data.frame(
      ME = MEs[[module_name]],
      Trait = traits[[trait_name]],
      Age = traits$Age_at_draw,
      Sex = traits$Sex,
      AD = traits$AD_status
    ) %>%
      stats::na.omit()

    if (nrow(df_reg) < 10) {
      return(NULL)
    }

    fit <- stats::lm(ME ~ Trait + Age + Sex + AD, data = df_reg)
    res <- coef(summary(fit))["Trait", ]
    data.frame(
      Module = module_name,
      Trait = trait_name,
      Beta = res[1],
      SE = res[2],
      P = res[4],
      stringsAsFactors = FALSE
    )
  }

  adj_results <- dplyr::bind_rows(lapply(analysis_map, function(x) run_adjusted_model(x$m, x$t)))
  comparison_table <- NULL
  if (nrow(adj_results) > 0) {
    comparison_table <- adj_results %>%
      dplyr::mutate(
        Original_Cor = vapply(
          seq_len(n()),
          function(i) stats::cor(MEs[[Module[i]]], traits[[Trait[i]]], use = "pairwise.complete.obs"),
          numeric(1)
        )
      )
  }

  top_hubs <- final_protein_kME %>%
    dplyr::filter(ModuleColor %in% selected_modules) %>%
    tidyr::pivot_longer(
      cols = dplyr::starts_with("kME"),
      names_to = "kME_Type",
      values_to = "kME_Value"
    ) %>%
    dplyr::filter(gsub("^kME", "", kME_Type) == ModuleColor) %>%
    dplyr::group_by(ModuleColor) %>%
    dplyr::slice_max(order_by = kME_Value, n = 10, with_ties = FALSE) %>%
    dplyr::ungroup()

  key_cor <- list(ADNI = NULL, ADRC = NULL)
  key_p <- list(ADNI = NULL, ADRC = NULL)

  if (!is.null(traits) && sum(cohort == "ADNI", na.rm = TRUE) > 5 && sum(cohort == "ADRC", na.rm = TRUE) > 5) {
    numeric_traits <- traits %>% dplyr::select(where(is.numeric))
    common_modules <- intersect(paste0("ME", selected_modules), colnames(MEs))
    if (ncol(numeric_traits) > 0 && length(common_modules) > 0) {
      cor_adni <- WGCNA::cor(
        MEs[cohort == "ADNI" & !is.na(cohort), common_modules, drop = FALSE],
        numeric_traits[cohort == "ADNI" & !is.na(cohort), , drop = FALSE],
        use = "pairwise.complete.obs"
      )
      p_adni <- WGCNA::corPvalueStudent(cor_adni, sum(cohort == "ADNI", na.rm = TRUE))

      cor_adrc <- WGCNA::cor(
        MEs[cohort == "ADRC" & !is.na(cohort), common_modules, drop = FALSE],
        numeric_traits[cohort == "ADRC" & !is.na(cohort), , drop = FALSE],
        use = "pairwise.complete.obs"
      )
      p_adrc <- WGCNA::corPvalueStudent(cor_adrc, sum(cohort == "ADRC", na.rm = TRUE))

      key_cor <- list(ADNI = cor_adni, ADRC = cor_adrc)
      key_p <- list(ADNI = p_adni, ADRC = p_adrc)
    }
  }

  list(
    selected_modules = selected_modules,
    preservation_summary = preservation_summary,
    preservation_raw = preservation_raw,
    adjusted_results = adj_results,
    comparison_table = comparison_table,
    protein_kme = final_protein_kME,
    top_hubs = top_hubs,
    cohort_cor = key_cor,
    cohort_p = key_p
  )
}

plot_supp_preservation <- function(supp_bundle) {
  df <- supp_bundle$preservation_summary
  if (is.null(df) || nrow(df) == 0) {
    return(empty_plot(
      title = "Module preservation unavailable",
      subtitle = "Both ADNI and ADRC cohorts with sufficient samples are required."
    ))
  }

  plot_data <- df %>%
    dplyr::mutate(is_key = ifelse(Module %in% supp_bundle$selected_modules, "Key", "Other"))

  ggplot2::ggplot(plot_data, ggplot2::aes(x = Size, y = Z_ADRC)) +
    ggplot2::geom_hline(yintercept = 2, linetype = "dashed", color = "#0f766e") +
    ggplot2::geom_hline(yintercept = 10, linetype = "dashed", color = "#9f1239") +
    ggplot2::geom_point(ggplot2::aes(color = is_key), size = 3) +
    ggrepel::geom_text_repel(ggplot2::aes(label = Module), size = 3.5, max.overlaps = 15) +
    ggplot2::scale_color_manual(values = c(Key = "#9f1239", Other = "#475569")) +
    app_theme() +
    ggplot2::labs(
      title = "Module preservation: ADNI to ADRC",
      subtitle = "Zsummary versus module size",
      x = "Module size",
      y = "Zsummary preservation score",
      color = NULL
    )
}

plot_supp_forest <- function(supp_bundle) {
  df <- supp_bundle$comparison_table
  if (is.null(df) || nrow(df) == 0) {
    return(empty_plot(
      title = "Sensitivity analysis unavailable",
      subtitle = "Required traits for the adjusted regression models were not all present."
    ))
  }

  forest_data <- df %>%
    dplyr::select(Module, Trait, Unadj_Beta = Original_Cor, Adj_Beta = Beta, SE) %>%
    tidyr::pivot_longer(
      cols = c(Unadj_Beta, Adj_Beta),
      names_to = "Model",
      values_to = "Estimate"
    )

  ggplot2::ggplot(
    forest_data,
    ggplot2::aes(x = Estimate, y = interaction(Module, Trait), color = Model)
  ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dotted", color = "#94a3b8") +
    ggplot2::geom_errorbarh(
      ggplot2::aes(xmin = Estimate - (1.96 * SE), xmax = Estimate + (1.96 * SE)),
      position = ggplot2::position_dodge(width = 0.5),
      height = 0.2
    ) +
    ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.5), size = 3) +
    ggplot2::scale_color_manual(values = c(Unadj_Beta = "#0f766e", Adj_Beta = "#9f1239")) +
    app_theme() +
    ggplot2::labs(
      title = "Sensitivity analysis with AD adjustment",
      subtitle = "95% confidence intervals around adjusted effect estimates",
      x = "Effect size",
      y = "Module-trait pair",
      color = NULL
    )
}

plot_supp_top_hubs <- function(supp_bundle) {
  df <- supp_bundle$top_hubs
  if (is.null(df) || nrow(df) == 0) {
    return(empty_plot(title = "Hub protein plot unavailable"))
  }

  ggplot2::ggplot(df, ggplot2::aes(x = reorder(ProteinID, kME_Value), y = kME_Value)) +
    ggplot2::geom_segment(
      ggplot2::aes(xend = reorder(ProteinID, kME_Value), yend = 0),
      color = "#cbd5e1"
    ) +
    ggplot2::geom_point(ggplot2::aes(color = ModuleColor), size = 3.5) +
    ggplot2::coord_flip() +
    ggplot2::facet_wrap(~ModuleColor, scales = "free_y") +
    ggplot2::scale_color_identity() +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank()) +
    ggplot2::labs(
      title = "Top hub proteins across selected modules",
      subtitle = "Ranked by module membership (kME)",
      x = "Protein",
      y = "kME strength"
    )
}

plot_supp_cohort_heatmap <- function(supp_bundle, cohort_name) {
  cor_matrix <- supp_bundle$cohort_cor[[cohort_name]]
  p_matrix <- supp_bundle$cohort_p[[cohort_name]]
  if (is.null(cor_matrix) || is.null(p_matrix)) {
    plot.new()
    text(0.5, 0.5, paste("No cohort heatmap available for", cohort_name))
    return(invisible(NULL))
  }

  textMatrix <- paste(signif(cor_matrix, 2), "\n(", signif(p_matrix, 1), ")", sep = "")
  dim(textMatrix) <- dim(cor_matrix)
  par(mar = c(10, 12, 4, 4))
  WGCNA::labeledHeatmap(
    Matrix = cor_matrix,
    xLabels = colnames(cor_matrix),
    yLabels = rownames(cor_matrix),
    ySymbols = rownames(cor_matrix),
    colorLabels = FALSE,
    colors = WGCNA::blueWhiteRed(50),
    textMatrix = textMatrix,
    setStdMargins = FALSE,
    cex.text = 0.6,
    cex.lab.x = 0.8,
    cex.lab.y = 0.8,
    zlim = c(-1, 1),
    main = paste(cohort_name, "module-trait relationships")
  )
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
    ggplot2::geom_histogram(binwidth = 1, fill = "#0f766e", color = "white") +
    app_theme() +
    ggplot2::labs(
      title = "Distribution of missing values per protein",
      subtitle = "Quality-control view before downstream harmonization",
      x = "Number of missing samples",
      y = "Count of proteins"
    )
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
    app_theme() +
    ggplot2::labs(
      title = "Signal distribution before and after imputation",
      subtitle = "Comparing observed intensities with the imputed expression matrix",
      x = "Intensity",
      y = "Density"
    ) +
    ggplot2::scale_fill_manual(values = c("Before Imputation" = "#94a3b8", "After Imputation" = "#9f1239"))
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
  expr_numeric <- expr_data[, is_num, drop = FALSE]
  dat <- t(as.matrix(expr_numeric))
  
  # Check for NAs (ComBat doesn't like them)
  if (any(is.na(dat))) stop("Expression data contains NAs. Run preprocessing first.")
  
  # Batch vector
  if (!batch_col %in% colnames(pheno_data)) stop(paste("Batch column", batch_col, "not found in phenotype data"))
  batch <- pheno_data[[batch_col]]
  
  if(length(unique(batch)) < 2) stop("Batch variable must have at least 2 levels.")
  if(any(is.na(batch))) stop("Batch variable contains NAs.")
  
  # Run ComBat
  mod <- model.matrix(~1, data = pheno_data)
  pca_input_before <- expr_numeric[, apply(expr_numeric, 2, stats::var, na.rm = TRUE) > 0, drop = FALSE]
  if (ncol(pca_input_before) < 2) {
    stop("Not enough non-constant numeric proteomic features available for harmonization PCA diagnostics.")
  }
  pca_before <- stats::prcomp(pca_input_before, scale. = TRUE)
  combat_edata <- sva::ComBat(dat = dat, batch = batch, mod = mod, par.prior = TRUE, prior.plots = FALSE)
  harmonized_expr <- t(combat_edata)
  pca_input_after <- harmonized_expr[, apply(harmonized_expr, 2, stats::var, na.rm = TRUE) > 0, drop = FALSE]
  if (ncol(pca_input_after) < 2) {
    stop("Not enough non-constant numeric proteomic features available after ComBat for PCA diagnostics.")
  }
  pca_after <- stats::prcomp(pca_input_after, scale. = TRUE)
  
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
    app_theme() +
    ggplot2::labs(
      title = paste(title_prefix, "PCA"),
      subtitle = paste("Samples colored by", batch_col),
      color = batch_col
    )
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
    ggplot2::geom_col(fill = "#0f766e") +
    ggplot2::geom_line(ggplot2::aes(group = 1), color = "#9f1239") +
    ggplot2::geom_point(color = "#9f1239") +
    app_theme() +
    ggplot2::labs(
      title = "Scree plot",
      subtitle = "Variance explained by the first 10 principal components",
      x = "Principal component",
      y = "Proportion of variance explained"
    )
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
  p + ggplot2::geom_point(alpha = 0.7, size = 2) + app_theme() + ggplot2::labs(
    title = "PCA projection",
    subtitle = "Samples projected onto the first two components"
  )
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
    app_theme() +
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
