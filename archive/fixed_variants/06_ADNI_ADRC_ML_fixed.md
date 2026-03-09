---
title: "06_ADNI_ADRC_ML (fixed)"
author: "automated-patch"
date: "2026-03-04"
output: pdf_document
---

# Setup



# Libraries


``` r
pkgs <- c("caret", "glmnet", "pROC", "dplyr", "here")
# Check for missing packages and install them
missing_pkgs <- pkgs[!(pkgs %in% installed.packages()[,"Package"])]
if(length(missing_pkgs) > 0){
  message("Installing missing packages: ", paste(missing_pkgs, collapse=", "))
  install.packages(missing_pkgs, dependencies = TRUE, repos = "https://cloud.r-project.org")
}

# Now, load all packages
lapply(pkgs, library, character.only = TRUE)
```

```
## [[1]]
##  [1] "caret"     "lattice"   "ggplot2"   "here"      "stats"     "graphics" 
##  [7] "grDevices" "utils"     "datasets"  "methods"   "base"     
## 
## [[2]]
##  [1] "glmnet"    "Matrix"    "caret"     "lattice"   "ggplot2"   "here"     
##  [7] "stats"     "graphics"  "grDevices" "utils"     "datasets"  "methods"  
## [13] "base"     
## 
## [[3]]
##  [1] "pROC"      "glmnet"    "Matrix"    "caret"     "lattice"   "ggplot2"  
##  [7] "here"      "stats"     "graphics"  "grDevices" "utils"     "datasets" 
## [13] "methods"   "base"     
## 
## [[4]]
##  [1] "dplyr"     "pROC"      "glmnet"    "Matrix"    "caret"     "lattice"  
##  [7] "ggplot2"   "here"      "stats"     "graphics"  "grDevices" "utils"    
## [13] "datasets"  "methods"   "base"     
## 
## [[5]]
##  [1] "dplyr"     "pROC"      "glmnet"    "Matrix"    "caret"     "lattice"  
##  [7] "ggplot2"   "here"      "stats"     "graphics"  "grDevices" "utils"    
## [13] "datasets"  "methods"   "base"
```

# Load Data


``` r
# Load WGCNA outputs
# Expecting 03_output.RData in results/wgcna/ (from step 03)
wgcna_file <- file.path(base_dir, "results", "wgcna", "03_output.RData")

if(!file.exists(wgcna_file)){
  # Fallback to root results if not in wgcna subdir (legacy compatibility)
  wgcna_file_legacy <- file.path(base_dir, "results", "03_output.RData")
  if(file.exists(wgcna_file_legacy)){
    wgcna_file <- wgcna_file_legacy
  } else {
    stop("WGCNA output file not found. Expected at: ", wgcna_file)
  }
}
message("Loading: ", wgcna_file)
```

```
## Loading: /Users/informatics/Desktop/Multi-Omics-Biomarker-Discovery/Proteomics/results/wgcna/03_output.RData
```

``` r
load(wgcna_file)

# Check objects
req_objs <- c("MEs", "traits", "expr", "moduleColors")
missing_objs <- setdiff(req_objs, ls())
if(length(missing_objs) > 0) stop("Missing objects in loaded RData: ", paste(missing_objs, collapse=", "))
```

# Feature Matrix Preparation


``` r
# Align samples
common_samples <- intersect(rownames(MEs), rownames(traits))
if(length(common_samples) == 0) stop("No common samples between MEs and traits.")

MEs <- MEs[common_samples, , drop=FALSE]
traits <- traits[common_samples, , drop=FALSE]
expr <- expr[common_samples, , drop=FALSE]

# Define target
target_col <- "stroke" 
if(!target_col %in% colnames(traits)) stop("Target column '", target_col, "' not found in traits.")

# Combine features + outcome
ml_df <- data.frame(
  MEs,
  outcome = factor(traits[[target_col]])
)

# Remove samples with missing outcome
ml_df <- ml_df[!is.na(ml_df$outcome), ]
message("Samples for ML: ", nrow(ml_df))
```

```
## Samples for ML: 995
```

``` r
message("Outcome distribution:")
```

```
## Outcome distribution:
```

``` r
print(table(ml_df$outcome))
```

```
## 
##   0   1 
## 899  96
```

``` r
# Ensure at least 2 classes
if(length(unique(ml_df$outcome)) < 2) stop("Outcome has fewer than 2 classes.")
```

# Train/Test Split


``` r
set.seed(42)
train_idx <- createDataPartition(ml_df$outcome, p = 0.7, list = FALSE)
train_df <- ml_df[train_idx, ]
test_df  <- ml_df[-train_idx, ]
```

# Feature Scaling


``` r
# Preprocess (center/scale) based on training data
preproc <- preProcess(train_df %>% select(-outcome), method = c("center", "scale"))

X_train <- predict(preproc, train_df %>% select(-outcome))
X_test  <- predict(preproc, test_df %>% select(-outcome))

y_train <- train_df$outcome
y_test  <- test_df$outcome
```

# Elastic Net Model


``` r
# glmnet requires numeric 0/1 outcome for binomial when y is not a two-column matrix
y01 <- as.numeric(y_train == levels(y_train)[2])

# CV for lambda
cv_fit <- cv.glmnet(
  x = as.matrix(X_train),
  y = y01,
  family = "binomial",
  alpha = 0.5, # Elastic Net
  nfolds = 5,
  standardize = FALSE # already centered/scaled via caret preProcess
)

best_lambda <- cv_fit$lambda.min
message("Best Lambda: ", best_lambda)
```

```
## Best Lambda: 0.046206567333412
```

``` r
plot(cv_fit)
```

![plot of chunk model](figure/model-1.png)

# Evaluation


``` r
# Predict probabilities
test_probs <- as.numeric(predict(cv_fit, newx = as.matrix(X_test), s = best_lambda, type = "response"))

# Convert to class labels (threshold 0.5)
levels_y <- levels(y_test)
test_pred_class <- ifelse(test_probs > 0.5, levels_y[2], levels_y[1])
test_pred_factor <- factor(test_pred_class, levels = levels_y)

# Confusion Matrix
cm <- confusionMatrix(test_pred_factor, y_test, positive = levels_y[2])
print(cm)
```

```
## Confusion Matrix and Statistics
## 
##           Reference
## Prediction   0   1
##          0 269  28
##          1   0   0
##                                           
##                Accuracy : 0.9057          
##                  95% CI : (0.8666, 0.9364)
##     No Information Rate : 0.9057          
##     P-Value [Acc > NIR] : 0.5501          
##                                           
##                   Kappa : 0               
##                                           
##  Mcnemar's Test P-Value : 3.352e-07       
##                                           
##             Sensitivity : 0.00000         
##             Specificity : 1.00000         
##          Pos Pred Value :     NaN         
##          Neg Pred Value : 0.90572         
##              Prevalence : 0.09428         
##          Detection Rate : 0.00000         
##    Detection Prevalence : 0.00000         
##       Balanced Accuracy : 0.50000         
##                                           
##        'Positive' Class : 1               
## 
```

``` r
# ROC (ensure response is a factor with expected ordering)
roc_obj <- roc(response = y_test, predictor = test_probs, levels = levels_y, direction = "<", quiet = TRUE)
auc_val <- as.numeric(auc(roc_obj))
message("AUC: ", auc_val)
```

```
## AUC: 0.5
```

``` r
png(file.path(results_dir, "ROC_Curve.png"))
plot(roc_obj, main = paste0("ROC Curve (AUC = ", round(auc_val, 3), ")"))
dev.off()
```

```
## quartz_off_screen 
##                 2
```

# Interpretation


``` r
# Extract coefficients
coef_mat <- coef(cv_fit, s = best_lambda)

coef_df <- data.frame(
  Module = rownames(coef_mat),
  Coefficient = as.numeric(coef_mat),
  row.names = NULL,
  stringsAsFactors = FALSE
)

# Remove intercept and filter non-zero
coef_df <- coef_df %>%
  filter(Module != "(Intercept)", Coefficient != 0) %>%
  arrange(desc(abs(Coefficient)))

print(coef_df)
```

```
## [1] Module      Coefficient
## <0 rows> (or 0-length row.names)
```

``` r
write.csv(coef_df, file.path(results_dir, "Model_Coefficients.csv"), row.names = FALSE)

# Map top module to proteins (moduleColors contains per-protein module assignment; values are colors)
if (nrow(coef_df) > 0) {
  top_module <- coef_df$Module[1]

  if (!startsWith(top_module, "ME")) {
    message("Top feature is not a module eigengene (expected prefix 'ME'): ", top_module)
  } else {
    top_color <- sub("^ME", "", top_module)
    message("Top Module Color: ", top_color)

    if (top_color %in% unique(moduleColors)) {
      proteins_in_top <- colnames(expr)[moduleColors == top_color]
      message("Number of proteins in top module: ", length(proteins_in_top))
      print(utils::head(proteins_in_top))

      write.csv(
        data.frame(Protein = proteins_in_top),
        file.path(results_dir, paste0("Proteins_in_", top_color, ".csv")),
        row.names = FALSE
      )
    } else {
      message("Module color not found in moduleColors values: ", top_color)
    }
  }
} else {
  message("No non-zero coefficients found.")
}
```

```
## No non-zero coefficients found.
```
