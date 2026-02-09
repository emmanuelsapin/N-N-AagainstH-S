# =============================================================================
# GTNAall2.R -- GTNA pipeline: merge, filter, GMM, supervised & 4-D mixture
# =============================================================================
#
# Purpose:
#   Pipeline for predicting Half-Sib (H-S) vs Not Applicable (N/A) pairs from
#   genotype/pi-hat features. Combines: (1) data read & merge with maximums,
#   (2) exclusions and Pihat/Pihat2 filters, (3) univariate 2-component GMMs
#   per variable, (4) supervised logistic classifier (4-D features + calibration),
#   (5) multivariate 2-component Gaussian mixture on 4-D features → P(H-S),
#   (6) cutoff selection (sensitivity=1 + max specificity, or max accuracy),
#   (7) mixture/supervised plots and TSV outputs, (8) scatter/diagnostics,
#   (9) final output table with predictions.
#
# Usage:
#   Rscript GTNAall2.R <input_file_path>
#   Example: Rscript GTNAall2.R /path/to/pairs_with_type.tsv
#
# Dependencies: mixtools, data.table; optional: glmnet, mvtnorm
#
# Main outputs (paths are configurable below):
#   - GMM_plot_*.pdf          : univariate GMM histograms (Case/Control, All)
#   - GMM_plot_supervised.pdf : supervised classifier calibration, P(HS), ROC
#   - GMM_mixture_4D_params.rds: saved 4-D mixture parameters (reusable)
#   - GMM_plot_mixture.pdf     : P(HS) histograms, cutoff, calibration, ROC
#   - GMM_mixture_PHS_by_group.txt : P(HS) by group, FN/FP at chosen cutoff
#   - GMM_plot_PHS_vs_age_diff.pdf : P(H-S) vs age difference
#   - GTNA_predictions.txt    : full table with P_HS_supervised, P_HS_mixture, etc.
#
# =============================================================================

require("mixtools")
require("data.table")

# --- Optional: plot before/after improvement (if col39_pair_debug.tsv exists) ---
# Plot real before vs after values from compareimprovement output
improvement_data_path <- "col39_pair_debug.tsv"
if (file.exists(improvement_data_path)) {
    improvement_data <- tryCatch(
        read.table(improvement_data_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE),
        error = function(e) NULL
    )
    if (!is.null(improvement_data) && all(c("before", "after") %in% names(improvement_data))) {
        before_vals <- as.numeric(improvement_data$before)
        after_vals <- as.numeric(improvement_data$after)
        valid_pairs <- !is.na(before_vals) & !is.na(after_vals)
        if (any(valid_pairs)) {
            bx <- before_vals[valid_pairs]
            ay <- after_vals[valid_pairs]
            fit_temp <- lm(ay ~ bx)
            resid_abs <- abs(residuals(fit_temp))
            outlier_idx <- which.max(resid_abs)
            keep_plot <- seq_along(bx) != outlier_idx
            bx_plot <- bx[keep_plot]
            ay_plot <- ay[keep_plot]
            pdf("improvement.pdf", width = 8, height = 6)
            par(cex.lab = 1.35, cex.axis = 1.2, cex.main = 1.25)
            plot(
                x = bx_plot,
                y = ay_plot,
                xlab = expression("Niece/nephew" ~ italic(ACPA) ~ "excluding avuncular"),
                ylab = expression("Niece/nephew" ~ italic(ACPA) ~ "keeping avuncular"),
                main = "",
                pch = 19,
                cex = 0.6
            )
            grid(col = "lightgray", lty = "dotted")
            abline(0, 1, col = "gray", lty = 2, lwd = 1.5)
            fit <- lm(ay_plot ~ bx_plot)
            abline(fit, col = "red", lwd = 2)
            dev.off()
            cat("Improvement plot saved to improvement.pdf using", sum(valid_pairs), "pairs\n")
        } else {
            cat("No valid before/after pairs found in", improvement_data_path, "\n")
        }
    } else {
        cat("Could not read required columns from", improvement_data_path, "\n")
    }
} else {
    cat("Improvement data file not found:", improvement_data_path, "\n")
}

# ============================================================================
# 1. Read data and parse columns
# ============================================================================
# Input: path to tab-separated file (pairs with type, Pihat, sum_pihat_* columns).
# Get file path from first command-line argument.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
    stop("Error: Please provide the input file path as the first argument")
}
file_path <- args[1]
cat("Reading data from:", file_path, "\n")

# Read main input file (no header; columns assigned below).
data <- read.table(file_path, header = FALSE, stringsAsFactors = FALSE, fill = TRUE)
cat("Data read from file_path: ", nrow(data), " rows\n")

# Merge with maximums file if present (adds max_pihatagainstall* columns by UKBID1/UKBID2).
# Read and merge GTNA_predictions_maximums.txt with data if it exists
maximums_path <- "/pl/active/KellerLab/Emmanuel/NAvHS/GTNA_predictions_maximums.txt"
data_maximums <- NULL
if (file.exists(maximums_path)) {
    cat("\nReading maximums file:", maximums_path, "\n")
    data_maximums <- tryCatch(
        read.table(maximums_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE),
        error = function(e) {
            cat("Impossible de lire maximums_path:", maximums_path, "\n")
            return(NULL)
        }
    )
    
    # Remove identical lines in data_maximums
    if (!is.null(data_maximums) && nrow(data_maximums) > 0) {
        n_before <- nrow(data_maximums)
        data_maximums <- unique(data_maximums)
        n_after <- nrow(data_maximums)
        if (n_before > n_after) {
            cat("Removed ", n_before - n_after, " duplicate lines from data_maximums\n")
        }
    }
    
    if (!is.null(data_maximums) && nrow(data_maximums) > 0) {
        cat("Maximums file read: ", nrow(data_maximums), " lignes, ", ncol(data_maximums), " colonnes\n")
    } else {
        cat("Maximums file is empty or could not be read\n")
        data_maximums <- NULL
    }
} else {
    cat("Maximums file does not exist, skipping merge\n")
    data_maximums <- NULL
}

# Assign column names: base columns (IDs, ages, Pihat, genotype probs), then
# sum_pihat_* or pihat_* depending on format, and finally "type".
# Nommer les colonnes de base
base_cols <- c("numtrio", "UKBID1", "UKBID2", "ID1", "ID2", "age1", "age2",
               "Pihat", 
               "p1_00", "p1_01", "p1_02", "p1_10", "p1_11", "p1_12", "p1_20", "p1_21", "p1_22",
               "Pihat2",
               "p2_00", "p2_01", "p2_02", "p2_10", "p2_11", "p2_12", "p2_20", "p2_21", "p2_22",
               "nsnp_00", "nsnp_01", "nsnp_02", "nsnp_10", "nsnp_11", "nsnp_12", "nsnp_20", "nsnp_21", "nsnp_22")

# Nommer les colonnes restantes
remaining_cols <- ncol(data) - length(base_cols) - 1  # -1 pour "type"
col_names <- base_cols

if (!is.na(remaining_cols) && (remaining_cols == 40 || remaining_cols == 41)) {
    # Format avec sommes
    col_names <- c(col_names, "sum_pihat_max")
    for (j in 1:9) col_names <- c(col_names, paste0("sum_pihat_max_", j))
    col_names <- c(col_names, "sum_pihat_opposite")
    for (j in 1:9) col_names <- c(col_names, paste0("sum_pihat_opposite_", j))
    col_names <- c(col_names, "sum_pihat_remaining_max") 
    for (j in 1:9) col_names <- c(col_names, paste0("sum_pihat_remaining_max_", j))
    col_names <- c(col_names, "sum_pihat_remaining_min")
    for (j in 1:9) col_names <- c(col_names, paste0("sum_pihat_remaining_min_", j))
    if (remaining_cols == 41) col_names <- c(col_names, "V_extra")
} else if (!is.na(remaining_cols) && remaining_cols == 36) {
    # Format sans sommes
    for (i in 1:4) {
        for (j in 1:9) col_names <- c(col_names, paste0("pihat_", i, "_", j))
    }
} else {
    # Format inattendu
    if (!is.na(remaining_cols) && remaining_cols > 0) {
        for (i in 1:remaining_cols) col_names <- c(col_names, paste0("V", i))
    }
}

col_names <- c(col_names, "type")
if (length(col_names) != ncol(data)) {
    if (length(col_names) > ncol(data)) {
        col_names <- col_names[1:ncol(data)]
        col_names[ncol(data)] <- "type"
    } else {
        extra_cols <- paste0("V", (length(col_names) + 1):(ncol(data) - 1))
        col_names <- c(col_names, extra_cols, "type")
    }
}
colnames(data) <- col_names
cat("Colonnes nommées:", ncol(data), "\n")

# Genotype-scale correction for Pihat/Pihat2 (chromosome length ratio).
# Apply correction factor to Pihat and Pihat2 (multiply by 330005/327744)
correction_factor <- 330005 / 327744
if ("Pihat" %in% names(data)) {
    if (!is.numeric(data$Pihat)) data$Pihat <- as.numeric(as.character(data$Pihat))
    n_pihat_before <- sum(!is.na(data$Pihat))
    data$Pihat <- data$Pihat * correction_factor
    cat("Applied correction factor (", correction_factor, ") to Pihat: ", n_pihat_before, " non-NA values corrected\n", sep="")
}
if ("Pihat2" %in% names(data)) {
    if (!is.numeric(data$Pihat2)) data$Pihat2 <- as.numeric(as.character(data$Pihat2))
    n_pihat2_before <- sum(!is.na(data$Pihat2))
    data$Pihat2 <- data$Pihat2 * correction_factor
    cat("Applied correction factor (", correction_factor, ") to Pihat2: ", n_pihat2_before, " non-NA values corrected\n", sep="")
}

# Print type counts after reading data 
if ("type" %in% names(data)) {
    cat("\n=== Type counts AFTER reading data ===\n")
    type_counts_after_read <- table(data$type, useNA = "ifany")
    print(type_counts_after_read)
    cat("Total observations: ", nrow(data), "\n")
} else {
    cat("Warning: 'type' column not found after column naming\n")
}

# Merge data_maximums with data using UKBID1/UKBID2 = indiv1/indiv2
if (!is.null(data_maximums) && nrow(data_maximums) > 0) {
    cat("\nMerging maximums data with data...\n")
    
    # Merge on UKBID1/UKBID2 = indiv1/indiv2
    # Convert UKBID1 and UKBID2 to numeric for merging
    if ("UKBID1" %in% names(data) && "UKBID2" %in% names(data)) {
        data$UKBID1 <- as.numeric(as.character(data$UKBID1))
        data$UKBID2 <- as.numeric(as.character(data$UKBID2))
        data_maximums$indiv1 <- as.numeric(as.character(data_maximums$indiv1))
        data_maximums$indiv2 <- as.numeric(as.character(data_maximums$indiv2)) 
        
        # Perform merge
        data_merged <- merge(data, data_maximums, 
                            by.x = c("UKBID1", "UKBID2"), 
                            by.y = c("indiv1", "indiv2"), 
                            all.x = TRUE)
        data_merged1 <- merge(data, data_maximums, 
                            by.x = c("UKBID1", "UKBID2"), 
                            by.y = c("indiv2", "indiv1"), 
                            all.x = TRUE) 
        # Append data_merged1 to data_merged, but only keep unique rows
        if (!is.null(data_merged1) && nrow(data_merged1) > 0) {
            # Ensure both data frames have the same columns
            common_cols <- intersect(names(data_merged), names(data_merged1))
            data_merged_temp <- rbind(data_merged[, common_cols], data_merged1[, common_cols])
            cat("Appended data_merged1 to data_merged: ", nrow(data_merged_temp), " total lignes (before deduplication)\n")
            
            # Remove duplicates based on original data identifiers
            # Use UKBID1 and UKBID2 as primary keys (these should uniquely identify each row)
            n_before_dedup <- nrow(data_merged_temp)
            
            if ("UKBID1" %in% names(data_merged_temp) && "UKBID2" %in% names(data_merged_temp)) {
                # Sort to prioritize rows with maximums data (non-NA values in maximums columns)
                max_cols <- setdiff(names(data_maximums), c("indiv1", "indiv2"))
                max_cols_in_data <- intersect(max_cols, names(data_merged_temp))
                
                if (length(max_cols_in_data) > 0) {
                    # Create a priority score: rows with more non-NA maximums values come first
                    na_count <- rowSums(is.na(data_merged_temp[, max_cols_in_data, drop = FALSE]))
                    data_merged_temp$priority_score <- -na_count  # Negative so fewer NAs come first
                    data_merged_temp <- data_merged_temp[order(data_merged_temp$priority_score, decreasing = TRUE), ]
                    data_merged_temp$priority_score <- NULL
                } 
                
                # Remove duplicates based on UKBID1 and UKBID2, keeping first occurrence
                data_merged <- data_merged_temp[!duplicated(data_merged_temp[, c("UKBID1", "UKBID2")]), ]
            } else {
                # Fallback: use all original columns from data (before merge) as key
                original_cols <- names(data)
                key_cols <- intersect(original_cols, names(data_merged_temp))
                if (length(key_cols) > 0) {
                    data_merged <- data_merged_temp[!duplicated(data_merged_temp[, key_cols]), ]
                } else {
                    data_merged <- data_merged_temp
                }
            }
            
            n_after_dedup <- nrow(data_merged)
            if (n_before_dedup > n_after_dedup) {
                cat("Removed ", n_before_dedup - n_after_dedup, " duplicate rows after merge\n")
            } else {
                cat("No duplicates found to remove\n")
            }
        }
        
        cat("Merge completed: ", nrow(data_merged), " lignes (", nrow(data), " original, ", 
            sum(!is.na(data_merged$max_pihatagainstall1_idx)), " avec maximums)\n")
        
        # Replace data with merged data
        data <- data_merged
        cat("Data updated with maximums columns\n")
    } else {
        cat("Warning: UKBID1 or UKBID2 not found in data, cannot merge\n")
    }
}

# Print type counts before exclusions
if ("type" %in% names(data)) {
    cat("\n=== Type counts BEFORE exclusions ===\n")
    type_counts_before <- table(data$type, useNA = "ifany")
    print(type_counts_before)
    cat("Total observations: ", nrow(data), "\n")
}

# ============================================================================
# 2a. Exclusions based on type and maximums (type 11/7 rules)
# ============================================================================
# Exclude lines with type = 11 (and some type = 7) using max_pihatagainstall*
# criteria (thresholds seuil / seuil2). Apply exclusions to all data if maximums exist.
if ("max_pihatagainstall1_below_033_pihat1" %in% names(data)) {

    # Create exclusion condition for all rows in data
    exclude_condition_data <- rep(FALSE, nrow(data))
    
    if ("type" %in% names(data) && "max_pihatagainstall1_below_033_pihat1" %in% names(data)) {
        # Get values for easier handling
        type_vals <- data$type
        max1_1 <- as.numeric(data$max_pihatagainstall1_below_033_pihat1)
        max2_1 <- as.numeric(data$max_pihatagainstall2_below_033_pihat1)
        max1_2 <- as.numeric(data$max_pihatagainstall1_below_033_pihat2) 
        max2_2 <- as.numeric(data$max_pihatagainstall2_below_033_pihat2)
        
        # Handle NAs in comparisons - replace NA with FALSE for comparisons
        seuil=0.115
        seuil2=0.14
       
        exclude_condition_data <- 
            (is.na(max1_1) & (type_vals == 11 | type_vals == 7)) |
            (	 
                (type_vals == 11) &
                (!is.na(max1_1)) &
                (
                    (!is.na(max2_1) & max2_1 < seuil) |
                    (!is.na(max1_1) & max1_1 < seuil) | 
                    (!is.na(max2_2) & max2_2 < seuil) |
                    (!is.na(max1_2) & max1_2 < seuil) |
                    (!is.na(max2_1) & max2_1 > seuil2) |
                    (!is.na(max1_1) & max1_1 > seuil2) |
                    (!is.na(max2_2) & max2_2 > seuil2) |
                    (!is.na(max1_2) & max1_2 > seuil2)
                )
            ) |
            (	
                (type_vals == 7) &
                (!is.na(max1_1)) &
                (
                    (!is.na(max2_1) & max2_1 > 0.025) |
                    (!is.na(max1_2) & max1_2 > 0.025) 
                )
            )
        # Replace any remaining NAs with FALSE
        exclude_condition_data[is.na(exclude_condition_data)] <- FALSE
    }
    
    n_before_exclusion <- nrow(data)
    # Remove excluded rows from data
    data <- data[!exclude_condition_data, ]
    
    n_excluded <- n_before_exclusion - nrow(data)
    if (n_excluded > 0) {
        cat("\nExcluded ", n_excluded, " lines from data (type=11/7 exclusion criteria)\n")
        cat("Remaining lines in data: ", nrow(data), "\n")
    }
}

# Print type counts after exclusions
if ("type" %in% names(data)) {
    cat("\n=== Type counts AFTER exclusions ===\n")
    type_counts_after_excl <- table(data$type, useNA = "ifany")
    print(type_counts_after_excl)
    cat("Total observations: ", nrow(data), "\n")
}

# ============================================================================
# 2b. Filtering by Pihat and Pihat2
# ============================================================================
# Keep pairs with Pihat in (0.2, 0.325] and Pihat2 in (-0.005, 0.03].
cat("\n=== Filtering (Pihat / Pihat2) ===\n")
cat("Using data for filtering: ", nrow(data), " lignes\n")
n_before <- nrow(data)

# Filtrage Pihat > 0.3
if ("Pihat" %in% names(data)) {
    if (!is.numeric(data$Pihat)) data$Pihat <- as.numeric(as.character(data$Pihat))
    keep_rows <- is.na(data$Pihat) | data$Pihat <= 0.325
    data <- data[keep_rows, ]
    cat("Après filtrage Pihat <= 0.3: ", nrow(data), " lignes (", n_before - nrow(data), " exclues)\n")
    if ("type" %in% names(data)) {
        cat("Type counts after Pihat <= 0.3 filter:\n")
        print(table(data$type, useNA = "ifany"))
    }
}
n_before1 <- nrow(data)

# Filtrage Pihat > 0.2
if ("Pihat" %in% names(data)) {
    if (!is.numeric(data$Pihat)) data$Pihat <- as.numeric(as.character(data$Pihat))
    keep_rows <- is.na(data$Pihat) | data$Pihat > 0.20
    data <- data[keep_rows, ]
    cat("Après filtrage Pihat > 0.2: ", nrow(data), " lignes (", n_before1 - nrow(data), " exclues)\n")
    if ("type" %in% names(data)) {
        cat("Type counts after Pihat > 0.2 filter:\n")
        print(table(data$type, useNA = "ifany"))
    }
}
# Filtrage Pihat2 > 0.02
n_before_pihat2 <- nrow(data)
if ("Pihat2" %in% names(data)) {
    if (!is.numeric(data$Pihat2)) data$Pihat2 <- as.numeric(as.character(data$Pihat2))
    keep_rows <- is.na(data$Pihat2) | data$Pihat2 <= 0.03
    data <- data[keep_rows, ]
    cat("Après filtrage Pihat2 <= 0.02: ", nrow(data), " lignes (", n_before_pihat2 - nrow(data), " exclues)\n")
    if ("type" %in% names(data)) {
        cat("Type counts after Pihat2 <= 0.02 filter:\n")
        print(table(data$type, useNA = "ifany"))
    }
}
# Filtrage Pihat2 > -0.005
n_before_pihat2_min <- nrow(data)
if ("Pihat2" %in% names(data)) {
    if (!is.numeric(data$Pihat2)) data$Pihat2 <- as.numeric(as.character(data$Pihat2))
    keep_rows <- is.na(data$Pihat2) | data$Pihat2 > -0.005
    data <- data[keep_rows, ]
    cat("Après filtrage Pihat2 > -0.005: ", nrow(data), " lignes (", n_before_pihat2_min - nrow(data), " exclues)\n")
    if ("type" %in% names(data)) {
        cat("Type counts after Pihat2 > -0.005 filter:\n")
        print(table(data$type, useNA = "ifany"))
    }
}

# Print final type counts after all filtering
if ("type" %in% names(data)) {
    cat("\n=== Type counts AFTER all filtering ===\n")
    type_counts_final <- table(data$type, useNA = "ifany")
    print(type_counts_final)
    cat("Total observations: ", nrow(data), "\n")
}

# ============================================================================
# 3. Prepare data and define groups for GMM and downstream
# ============================================================================
# Build 7 derived variables from sum_pihat_max_5, sum_pihat_opposite_5, etc.;
# define groups: type 100 = Control (N/A), type 2/7/11 = Case (HS), else Undetermined.
cat("\n=== Data preparation and group definition ===\n")

# Extract 4-D sum columns and build composite variables.
# Prepare base data
sum_pihat_remaining_max_5_data <- data[["sum_pihat_remaining_max_5"]]
sum_pihat_opposite_5_data <- data[["sum_pihat_opposite_5"]]
sum_pihat_remaining_min_5_data <- data[["sum_pihat_remaining_min_5"]]
sum_pihat_max_5_data <- data[["sum_pihat_max_5"]]
type_data <- data$type

# Create 7 variables
var1_data <- sum_pihat_max_5_data
var2_data <- sum_pihat_opposite_5_data
var3_data <- sum_pihat_remaining_max_5_data
var4_data <- sum_pihat_remaining_min_5_data
var5_data <- sum_pihat_max_5_data + sum_pihat_opposite_5_data
var6_data <- sum_pihat_remaining_max_5_data + sum_pihat_remaining_min_5_data
var7_data <- sum_pihat_max_5_data + sum_pihat_opposite_5_data - sum_pihat_remaining_max_5_data - sum_pihat_remaining_min_5_data

# Filter valid data for each variable
idx_valid <- !is.na(var1_data) & is.finite(var1_data) & 
             !is.na(var2_data) & is.finite(var2_data) &
             !is.na(var3_data) & is.finite(var3_data) &
             !is.na(var4_data) & is.finite(var4_data) &
             !is.na(var5_data) & is.finite(var5_data) &
             !is.na(var6_data) & is.finite(var6_data) &
             !is.na(var7_data) & is.finite(var7_data) 

var1_clean <- var1_data[idx_valid]
var2_clean <- var2_data[idx_valid]
var3_clean <- var3_data[idx_valid]
var4_clean <- var4_data[idx_valid]
var5_clean <- var5_data[idx_valid]
var6_clean <- var6_data[idx_valid]
var7_clean <- var7_data[idx_valid]
type_valid <- type_data[idx_valid]
type_valid[is.na(type_valid)] <- 0

cat("Valid data for GMM: ", length(var1_clean), " observations\n")

# Define groups based on type BEFORE GMM:
# - type = 100 → Control (N/A)
# - type in (2, 7, 11) → Case (HS)  
# - others → Undetermined
group_labels <- rep("Undetermined", length(type_valid))
group_labels[type_valid == 100] <- "Control (N/A)"
group_labels[type_valid %in% c(2, 7, 11)] <- "Case (HS)"

# Colors for groups
group_colors <- c("Control (N/A)" = "lightskyblue1", "Case (HS)" = "#FF80FF", "Undetermined" = "lightgray")

cat("\nGroup distribution:\n")
print(table(group_labels))

# Define variable names and data
# Using expression() for axis labels with p^Het format
variable_list <- list(
    list(name = "π̂_max_Het", label = expression(pi[Max]^{hh}), data = var1_clean),
    list(name = "π̂_opposite_Het", label = expression(Opposite((hat(p)^Het))), data = var2_clean),
    list(name = "π̂_remaining_max_Het", label = expression(Remaining_max((hat(p)^Het))), data = var3_clean),
    list(name = "π̂_remaining_min_Het", label = expression(Remaining_min((hat(p)^Het))), data = var4_clean),
    list(name = "π̂_max_Het + π̂_opposite_Het", label = expression(Max((hat(p)^Het)) + Opposite((hat(p)^Het))), data = var5_clean),
    list(name = "π̂_remaining_max_Het + π̂_remaining_min_Het", label = expression(Remaining_max((hat(p)^Het)) + Remaining_min((hat(p)^Het))), data = var6_clean),
    list(name = "π̂_max_Het + π̂_opposite_Het - π̂_remaining_max_Het - π̂_remaining_min_Het", label = expression(Max((hat(p)^Het)) + Opposite((hat(p)^Het)) - Remaining_max((hat(p)^Het)) - Remaining_min((hat(p)^Het))), data = var7_clean)
)

# ============================================================================
# 4. Univariate 2-component GMM per variable
# ============================================================================
# For each variable: (GMM 1) fit on Case+Control only; (GMM 2) fit on all data.
# Used for histograms, cutoffs, and later TP/FN coloring in scatter plots.
calculate_gmm_for_variable <- function(var_data, var_name) {
    cat("\n=== Variable:", var_name, "===\n")

    # Restrict to labeled Case and Control for GMM 1.
    # Filter for Case and Control only
    idx_case_control <- group_labels %in% c("Control (N/A)", "Case (HS)")
    var_case_control <- var_data[idx_case_control]
    group_labels_case_control <- group_labels[idx_case_control]
    
    # GMM 1: Only Case and Control
    gmm_result_case_control <- NULL
    if (length(var_case_control) >= 2 && length(unique(var_case_control)) >= 2) {
    tryCatch({
            gmm_result_case_control <- mixtools::normalmixEM(var_case_control, k = 2, verb = FALSE, maxit = 1000)
            cat("GMM 1 (Case/Control only) successful:\n")
            cat("  Component 1: μ=", round(gmm_result_case_control$mu[1], 4), 
                ", σ=", round(gmm_result_case_control$sigma[1], 4), 
                ", λ=", round(gmm_result_case_control$lambda[1], 4), "\n")
            cat("  Component 2: μ=", round(gmm_result_case_control$mu[2], 4), 
                ", σ=", round(gmm_result_case_control$sigma[2], 4), 
                ", λ=", round(gmm_result_case_control$lambda[2], 4), "\n")
            # Print number of observations for each population
            cat("  Number of observations per population:\n")
            pop_counts <- table(group_labels_case_control)
            for (pop_name in names(pop_counts)) {
                cat("    ", pop_name, ": ", pop_counts[pop_name], "\n", sep = "")
            }
            # Assign observations to clusters and determine case/control detection
            prob_cluster1 <- gmm_result_case_control$lambda[1] * dnorm(var_case_control, 
                                                                       mean = gmm_result_case_control$mu[1], 
                                                                       sd = gmm_result_case_control$sigma[1])
            prob_cluster2 <- gmm_result_case_control$lambda[2] * dnorm(var_case_control, 
                                                                       mean = gmm_result_case_control$mu[2], 
                                                                       sd = gmm_result_case_control$sigma[2])
            total_prob <- prob_cluster1 + prob_cluster2
            prob_cluster1 <- prob_cluster1 / total_prob
            prob_cluster2 <- prob_cluster2 / total_prob
            cluster_labels <- ifelse(prob_cluster1 > prob_cluster2, 1, 2)
            # Determine which cluster corresponds to "Case (HS)"
            case_in_cluster1 <- sum(group_labels_case_control == "Case (HS)" & cluster_labels == 1)
            case_in_cluster2 <- sum(group_labels_case_control == "Case (HS)" & cluster_labels == 2)
            # The cluster with more Case observations is the "case" cluster
            if (case_in_cluster2 > case_in_cluster1) {
                detected_case <- sum(cluster_labels == 2)
                detected_control <- sum(cluster_labels == 1)
            } else {
                detected_case <- sum(cluster_labels == 1)
                detected_control <- sum(cluster_labels == 2)
            }
            cat("  Detected by GMM:\n")
            cat("    Case: ", detected_case, "\n", sep = "")
            cat("    Control: ", detected_control, "\n", sep = "")
            cat("    Total: ", length(var_case_control), "\n", sep = "")
            # Calculate cutoff between the two populations
            # The cutoff is where λ₁ * dnorm(x, μ₁, σ₁) = λ₂ * dnorm(x, μ₂, σ₂)
            # We'll find this numerically by searching between the two means
            mu1 <- gmm_result_case_control$mu[1]
            mu2 <- gmm_result_case_control$mu[2]
            sigma1 <- gmm_result_case_control$sigma[1]
            sigma2 <- gmm_result_case_control$sigma[2]
            lambda1 <- gmm_result_case_control$lambda[1]
            lambda2 <- gmm_result_case_control$lambda[2]
            # Search range: between the two means, with some padding
            search_min <- min(mu1, mu2) - 3 * max(sigma1, sigma2)
            search_max <- max(mu1, mu2) + 3 * max(sigma1, sigma2)
            # Function to find where the difference is zero
            find_cutoff <- function(x) {
                prob1 <- lambda1 * dnorm(x, mean = mu1, sd = sigma1)
                prob2 <- lambda2 * dnorm(x, mean = mu2, sd = sigma2)
                return(prob1 - prob2)
            }
            # Find the root using uniroot; if no sign change, use midpoint between means
            tryCatch({
                cutoff <- uniroot(find_cutoff, interval = c(search_min, search_max))$root
                cat("  Cutoff between populations: ", round(cutoff, 6), "\n", sep = "")
            }, error = function(e) {
                cutoff_fallback <- (mu1 + mu2) / 2
                cat("  Cutoff (fallback, midpoint of means): ", round(cutoff_fallback, 6), " [uniroot failed: ", e$message, "]\n", sep = "")
            })
    }, error = function(e) {
            cat("Error in GMM 1:", e$message, "\n")
            })
        } else {
        cat("Not enough data for GMM 1\n")
        }
        
    # GMM 2: All Data
    gmm_result_all <- NULL
    if (length(var_data) >= 2 && length(unique(var_data)) >= 2) {
        tryCatch({
            gmm_result_all <- mixtools::normalmixEM(var_data, k = 2, verb = FALSE, maxit = 1000)
            cat("GMM 2 (All Data) successful:\n")
            cat("  Component 1: μ=", round(gmm_result_all$mu[1], 4), 
                ", σ=", round(gmm_result_all$sigma[1], 4), 
                ", λ=", round(gmm_result_all$lambda[1], 4), "\n")
            cat("  Component 2: μ=", round(gmm_result_all$mu[2], 4), 
                ", σ=", round(gmm_result_all$sigma[2], 4), 
                ", λ=", round(gmm_result_all$lambda[2], 4), "\n")
            # Print number of observations for each population
            cat("  Number of observations per population:\n")
            # Get group_labels for all data (need to pass it to the function or use the global variable)
            # Since var_data is the full data, we need to use group_labels which corresponds to var_data
            # But group_labels is defined outside the function, so we need to filter it to match var_data
            # Actually, var_data is already filtered by idx_valid, so group_labels should match
            # But wait, var_data in this context is var_data[idx_valid], so group_labels should match
            # Let me check: var_data is passed as parameter, and group_labels is global
            # We need to get the corresponding group_labels for var_data
            # Since calculate_gmm_for_variable receives var_data which is already filtered, 
            # we need to pass group_labels as well or use the global one
            # Actually, looking at the code, var_data is var_info$data which is var1_clean, var2_clean, etc.
            # These are already filtered by idx_valid, and group_labels is also filtered by idx_valid
            # So group_labels should match var_data in length
            if (length(group_labels) == length(var_data)) {
                pop_counts_all <- table(group_labels)
                for (pop_name in names(pop_counts_all)) {
                    cat("    ", pop_name, ": ", pop_counts_all[pop_name], "\n", sep = "")
                }
                # Assign observations to clusters and determine case/control detection
                prob_cluster1_all <- gmm_result_all$lambda[1] * dnorm(var_data, 
                                                                      mean = gmm_result_all$mu[1], 
                                                                      sd = gmm_result_all$sigma[1])
                prob_cluster2_all <- gmm_result_all$lambda[2] * dnorm(var_data, 
                                                                      mean = gmm_result_all$mu[2], 
                                                                      sd = gmm_result_all$sigma[2])
                total_prob_all <- prob_cluster1_all + prob_cluster2_all
                prob_cluster1_all <- prob_cluster1_all / total_prob_all
                prob_cluster2_all <- prob_cluster2_all / total_prob_all
                cluster_labels_all <- ifelse(prob_cluster1_all > prob_cluster2_all, 1, 2)
                # Determine which cluster corresponds to "Case (HS)"
                case_in_cluster1_all <- sum(group_labels == "Case (HS)" & cluster_labels_all == 1)
                case_in_cluster2_all <- sum(group_labels == "Case (HS)" & cluster_labels_all == 2)
                # The cluster with more Case observations is the "case" cluster
                if (case_in_cluster2_all > case_in_cluster1_all) {
                    detected_case_all <- sum(cluster_labels_all == 2)
                    detected_control_all <- sum(cluster_labels_all == 1)
                } else {
                    detected_case_all <- sum(cluster_labels_all == 1)
                    detected_control_all <- sum(cluster_labels_all == 2)
                }
                cat("  Detected by GMM:\n")
                cat("    Case: ", detected_case_all, "\n", sep = "")
                cat("    Control: ", detected_control_all, "\n", sep = "")
                cat("    Total: ", length(var_data), "\n", sep = "")
                # Calculate cutoff between the two populations
                # The cutoff is where λ₁ * dnorm(x, μ₁, σ₁) = λ₂ * dnorm(x, μ₂, σ₂)
                # We'll find this numerically by searching between the two means
                mu1_all <- gmm_result_all$mu[1]
                mu2_all <- gmm_result_all$mu[2]
                sigma1_all <- gmm_result_all$sigma[1]
                sigma2_all <- gmm_result_all$sigma[2]
                lambda1_all <- gmm_result_all$lambda[1]
                lambda2_all <- gmm_result_all$lambda[2]
                # Search range: between the two means, with some padding
                search_min_all <- min(mu1_all, mu2_all) - 3 * max(sigma1_all, sigma2_all)
                search_max_all <- max(mu1_all, mu2_all) + 3 * max(sigma1_all, sigma2_all)
                # Function to find where the difference is zero
                find_cutoff_all <- function(x) {
                    prob1 <- lambda1_all * dnorm(x, mean = mu1_all, sd = sigma1_all)
                    prob2 <- lambda2_all * dnorm(x, mean = mu2_all, sd = sigma2_all)
                    return(prob1 - prob2)
                }
                # Find the root using uniroot; if no sign change, use midpoint between means
                tryCatch({
                    cutoff_all <- uniroot(find_cutoff_all, interval = c(search_min_all, search_max_all))$root
                    cat("  Cutoff between populations: ", round(cutoff_all, 6), "\n", sep = "")
                }, error = function(e) {
                    cutoff_fallback <- (mu1_all + mu2_all) / 2
                    cat("  Cutoff (fallback, midpoint of means): ", round(cutoff_fallback, 6), " [uniroot failed: ", e$message, "]\n", sep = "")
                })
            } else {
                cat("    Warning: group_labels length does not match var_data length\n")
            }
        }, error = function(e) {
            cat("Error in GMM 2:", e$message, "\n")
        })
        } else {
        cat("Not enough data for GMM 2\n")
    }
    
    return(list(
        gmm_case_control = gmm_result_case_control,
        gmm_all = gmm_result_all,
        var_case_control = var_case_control,
        group_labels_case_control = group_labels_case_control
    ))
}

# ============================================================================
# 5. GMM plot and accuracy (stacked histogram + density curves, TP/TN/FP/FN)
# ============================================================================
# Draws stacked histogram by group and GMM density curves; computes accuracy,
# sensitivity, specificity on labeled Case/Control; optionally breaks down by type.
create_gmm_plot <- function(gmm_result, data_values, group_labels_plot, title, var_name, type_data_plot = NULL) {
    if (is.null(gmm_result)) return(NULL)

    # Posterior: assign each observation to the component with higher probability.
    # Assign each observation to cluster with highest probability
    prob_cluster1 <- gmm_result$lambda[1] * dnorm(data_values, 
                                                   mean = gmm_result$mu[1], 
                                                   sd = gmm_result$sigma[1])
    prob_cluster2 <- gmm_result$lambda[2] * dnorm(data_values, 
                                                   mean = gmm_result$mu[2], 
                                                   sd = gmm_result$sigma[2])
    total_prob <- prob_cluster1 + prob_cluster2
    prob_cluster1 <- prob_cluster1 / total_prob
    prob_cluster2 <- prob_cluster2 / total_prob
    
    cluster_labels <- ifelse(prob_cluster1 > prob_cluster2, 1, 2)
    
    # Calculate accuracy only on Case/Control (exclude Undetermined)
    idx_labeled <- group_labels_plot %in% c("Control (N/A)", "Case (HS)")
        accuracy <- NA
        sensitivity <- NA
        specificity <- NA
    TP <- NA
    TN <- NA
    FP <- NA
    FN <- NA
    
    # Determine cluster alignment first (needed for all type calculations)
    cluster_aligned <- rep(NA, length(cluster_labels))
    if (sum(idx_labeled) > 0) {
        # Create binary labels: Case (HS) = 1, Control (N/A) = 0
        true_labels <- ifelse(group_labels_plot[idx_labeled] == "Case (HS)", 1, 0)
        pred_labels <- cluster_labels[idx_labeled]
        
        # Find which cluster best corresponds to "Case (HS)"
        case_in_cluster1 <- sum(group_labels_plot[idx_labeled] == "Case (HS)" & pred_labels == 1)
        case_in_cluster2 <- sum(group_labels_plot[idx_labeled] == "Case (HS)" & pred_labels == 2)
        
        # Align clusters: the cluster with the most "Case (HS)" becomes 1
        if (case_in_cluster2 > case_in_cluster1) {
            pred_labels_aligned <- ifelse(pred_labels == 2, 1, 0)
            # Apply alignment to all data
            cluster_aligned <- ifelse(cluster_labels == 2, 1, 0)
            } else {
            pred_labels_aligned <- ifelse(pred_labels == 1, 1, 0)
            # Apply alignment to all data
            cluster_aligned <- ifelse(cluster_labels == 1, 1, 0)
        }
        
        accuracy <- mean(true_labels == pred_labels_aligned)
        
        # Calculate sensitivity and specificity
        TP <- sum(true_labels == 1 & pred_labels_aligned == 1, na.rm = TRUE)
        TN <- sum(true_labels == 0 & pred_labels_aligned == 0, na.rm = TRUE)
        FP <- sum(true_labels == 0 & pred_labels_aligned == 1, na.rm = TRUE)
        FN <- sum(true_labels == 1 & pred_labels_aligned == 0, na.rm = TRUE)
        tot_pos <- TP + FN; tot_neg <- TN + FP
        sensitivity <- if (is.finite(tot_pos) && tot_pos > 0) TP / tot_pos else 0
        specificity <- if (is.finite(tot_neg) && tot_neg > 0) TN / tot_neg else 0
    }
    
    # Calculate TP, TN, FP, FN for all types if type_data_plot is provided
    if (!is.null(type_data_plot) && length(type_data_plot) == length(data_values)) {
        cat("\n=== TP, TN, FP, FN breakdown by type ===\n")
        
        # Get unique types
        unique_types <- sort(unique(type_data_plot[!is.na(type_data_plot)]))
        
        # For each type, calculate TP, TN, FP, FN
        # Case types (2, 7, 11): predicted as Case (1) = TP, predicted as Control (0) = FN
        # Control type (100): predicted as Control (0) = TN, predicted as Case (1) = FP
        # Other types: not included in TP/TN/FP/FN calculation
        
        for (type_val in unique_types) {
            idx_type <- type_data_plot == type_val & !is.na(cluster_aligned)
            
            if (sum(idx_type) > 0) {
                if (type_val %in% c(2, 7, 11)) {
                    # Case types: TP = predicted as Case (1), FN = predicted as Control (0)
                    tp_type <- sum(cluster_aligned[idx_type] == 1)
                    fn_type <- sum(cluster_aligned[idx_type] == 0)
                    tn_type <- 0
                    fp_type <- 0
                    cat("Type ", type_val, " (Case): TP=", tp_type, ", FN=", fn_type, 
                        ", TN=0, FP=0, Total=", sum(idx_type), "\n", sep="")
                } else if (type_val == 100) {
                    # Control type: TN = predicted as Control (0), FP = predicted as Case (1)
                    tn_type <- sum(cluster_aligned[idx_type] == 0)
                    fp_type <- sum(cluster_aligned[idx_type] == 1)
                    tp_type <- 0
                    fn_type <- 0
                    cat("Type ", type_val, " (Control): TN=", tn_type, ", FP=", fp_type, 
                        ", TP=0, FN=0, Total=", sum(idx_type), "\n", sep="")
                } else {
                    # Other types: not classified
                    cat("Type ", type_val, " (Other): Not classified, Total=", sum(idx_type), "\n", sep="")
                }
            }
        }
        
        # Print overall summary
        cat("\nOverall summary (Case vs Control only):\n")
        cat("  TP:", TP, "| TN:", TN, "| FP:", FP, "| FN:", FN, "\n")
    }
    
    # Single plot combining histogram and GMM
    par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1, cex.lab = 1.35, cex.axis = 1.2, cex.main = 1.25)
    
    # Colors for points
    point_colors_plot <- group_colors[group_labels_plot]
    
    # Prepare histogram
    hist_breaks <- seq(min(data_values), max(data_values), length.out = 31)
    group_names <- c("Control (N/A)", "Case (HS)", "Undetermined")
    n_bins <- length(hist_breaks) - 1
    stacked_counts <- matrix(0, nrow = n_bins, ncol = length(group_names))
    colnames(stacked_counts) <- group_names
    
    # Calculate counts for each group in each bin
    for (i in 1:n_bins) {
        bin_data <- data_values[data_values >= hist_breaks[i] & data_values < hist_breaks[i+1]]
        if (i == n_bins) {
            bin_data <- data_values[data_values >= hist_breaks[i] & data_values <= hist_breaks[i+1]]
        }
        if (length(bin_data) > 0) {
            bin_idx <- which(data_values %in% bin_data)
            for (grp in group_names) {
                stacked_counts[i, grp] <- sum(group_labels_plot[bin_idx] == grp)
            }
        }
    }
    
    max_total <- max(rowSums(stacked_counts), na.rm = TRUE)
    if (is.na(max_total) || max_total == 0) max_total <- 1
    
    # Calculate GMM density curves for scaling
    x_dens <- seq(min(data_values), max(data_values), length.out = 500)
    y_dens1 <- gmm_result$lambda[1] * dnorm(x_dens, mean = gmm_result$mu[1], sd = gmm_result$sigma[1])
    y_dens2 <- gmm_result$lambda[2] * dnorm(x_dens, mean = gmm_result$mu[2], sd = gmm_result$sigma[2])
    y_dens_total <- y_dens1 + y_dens2
    
    # Normalize GMM curves to match histogram scale
    max_dens <- max(y_dens_total, na.rm = TRUE)
    scale_factor <- max_total / max_dens * 0.8  # Scale to ~80% of histogram height
    
    # Determine which GMM component corresponds to Control (lightskyblue1)
    # Check which cluster has more Control observations
    idx_control <- group_labels_plot == "Control (N/A)"
    if (sum(idx_control) > 0) {
        control_in_cluster1 <- sum(idx_control & cluster_labels == 1)
        control_in_cluster2 <- sum(idx_control & cluster_labels == 2)
        
        # If cluster 2 has more Control, swap the components
        if (control_in_cluster2 > control_in_cluster1) {
            # Swap components
            y_dens_control <- y_dens2
            y_dens_case <- y_dens1
            mu_control <- gmm_result$mu[2]
            mu_case <- gmm_result$mu[1]
            col_control <- "blue"
            col_case <- "red"
        } else {
            # Keep original order
            y_dens_control <- y_dens1
            y_dens_case <- y_dens2
            mu_control <- gmm_result$mu[1]
            mu_case <- gmm_result$mu[2]
            col_control <- "blue"
            col_case <- "red"
        }
    } else {
        # No Control data, use default assignment
        y_dens_control <- y_dens1
        y_dens_case <- y_dens2
        mu_control <- gmm_result$mu[1]
        mu_case <- gmm_result$mu[2]
        col_control <- "blue"
        col_case <- "red"
    }
    
    # Draw stacked histogram (no title)
    plot(hist_breaks[1], 0, type = "n",
         xlim = range(hist_breaks),
         ylim = c(0, max_total * 1.15),
         xlab = var_name,
         ylab = "Frequency",
         main = "")
    
    # Add grid
    grid(col = "lightgray", lty = "dotted", lwd = 0.5)
    
    # Draw stacked histogram bars
    for (i in 1:n_bins) {
        bottom <- 0
        for (grp in group_names) {
            count <- stacked_counts[i, grp]
            if (count > 0) {
                rect(hist_breaks[i], bottom,
                     hist_breaks[i+1], bottom + count,
                     col = group_colors[grp],
                     border = "black")
                bottom <- bottom + count
            }
        }
    }
    
    # Add GMM density curves on top of histogram (blue for Control/lightskyblue1, red for Case)
    lines(x_dens, y_dens_control * scale_factor, col = col_control, lwd = 3, lty = 1)
    lines(x_dens, y_dens_case * scale_factor, col = col_case, lwd = 3, lty = 1)
    lines(x_dens, y_dens_total * scale_factor, col = "purple", lwd = 2, lty = 2)
    
    # Add GMM means as vertical lines
    abline(v = mu_control, col = col_control, lwd = 2, lty = 2)
    abline(v = mu_case, col = col_case, lwd = 2, lty = 2)
    
legend("top",
           legend = c(paste("N/A, n=", sum(group_labels_plot == "Control (N/A)")),
                     paste("HS, n=", sum(group_labels_plot == "Case (HS)")),
                     paste("Undetermined, n=", sum(group_labels_plot == "Undetermined")),
                     "GMM μ (population 1)", "GMM μ (population 2)", "GMM Distribution (population 1)", "GMM Distribution (population 2)", "GMM total"),
           col = c(group_colors, col_control, col_case, col_control, col_case, "purple"),
           pch = c(19, 19, 19, NA, NA, NA, NA, NA),
           lty = c(NA, NA, NA, 2, 2, 1, 1, 2),
           lwd = c(NA, NA, NA, 2, 2, 3, 3, 2),
           cex = 1.3)

    return(list(accuracy = accuracy, sensitivity = sensitivity, specificity = specificity,
                TP = TP, TN = TN, FP = FP, FN = FN, cluster_labels = cluster_labels))
}

# ============================================================================
# 6. Loop over 7 variables: fit GMMs and write one PDF per plot (14 PDFs total)
# ============================================================================
# GMM 1 = Case/Control only; GMM 2 = All data. Each gets a separate PDF file.
plot_counter <- 0

for (var_idx in 1:length(variable_list)) {
    var_info <- variable_list[[var_idx]]
    var_name <- var_info$name  # For printing
    var_label <- var_info$label  # For plotting (expression)
    var_data <- var_info$data
    
    cat("\n", rep("=", 70), "\n", sep = "")
    cat("Processing variable", var_idx, "of", length(variable_list), ":", var_name, "\n")
    cat(rep("=", 70), "\n", sep = "")
    
    # Calculate GMMs for this variable
    gmm_results <- calculate_gmm_for_variable(var_data, var_name)
    
    # Plot 1: GMM on Case/Control only 
    if (!is.null(gmm_results$gmm_case_control)) {
        plot_counter <- plot_counter + 1
        pdf_file <- paste0("/pl/active/KellerLab/Emmanuel/NAvHS/GMM_plot_", plot_counter, ".pdf")
        pdf(pdf_file, width = 12, height = 8)
        
        cat("\n=== Creating plot for GMM 1 (Case/Control only) ===\n")
        # For GMM 1, filter type data to match case/control only
        idx_case_control <- group_labels %in% c("Control (N/A)", "Case (HS)")
        type_data_case_control <- type_valid[idx_case_control]
        
        stats1 <- create_gmm_plot(gmm_results$gmm_case_control, 
                                   gmm_results$var_case_control, 
                                   gmm_results$group_labels_case_control,
                                   "",
                                   var_label,
                                   type_data_plot = type_data_case_control)
        
        dev.off()
        cat("Saved plot to:", pdf_file, "\n")
        
        if (!is.null(stats1) && !is.na(stats1$accuracy)) {
            cat("GMM 1 - Accuracy (Case vs Control):", round(stats1$accuracy, 4), "\n")
            cat("  Sensitivity:", round(stats1$sensitivity, 4), "\n")
            cat("  Specificity:", round(stats1$specificity, 4), "\n")
            cat("  TP:", stats1$TP, "| TN:", stats1$TN, "| FP:", stats1$FP, "| FN:", stats1$FN, "\n")
        }
    }
    
    # Plot 2: GMM on all data
    if (!is.null(gmm_results$gmm_all)) {
        plot_counter <- plot_counter + 1
        pdf_file <- paste0("/pl/active/KellerLab/Emmanuel/NAvHS/GMM_plot_", plot_counter, ".pdf")
        pdf(pdf_file, width = 12, height = 8)
        
        cat("\n=== Creating plot for GMM 2 (All Data) ===\n")
        stats2 <- create_gmm_plot(gmm_results$gmm_all, 
                                   var_data, 
                                   group_labels,
                                   "",
                                   var_label,
                                   type_data_plot = type_valid)
        
        dev.off()
        cat("Saved plot to:", pdf_file, "\n")
        
        if (!is.null(stats2) && !is.na(stats2$accuracy)) {
            cat("GMM 2 - Accuracy (Case vs Control):", round(stats2$accuracy, 4), "\n")
            cat("  Sensitivity:", round(stats2$sensitivity, 4), "\n")
            cat("  Specificity:", round(stats2$specificity, 4), "\n")
            cat("  TP:", stats2$TP, "| TN:", stats2$TN, "| FP:", stats2$FP, "| FN:", stats2$FN, "\n")
        }
    }
}

# ============================================================================
# 6b. Supervised classifier: 4-D features, ridge logistic, Platt calibration, CV
# ============================================================================
# Features: pi_hat_max, pi_hat_R_max, pi_hat_R_min, pi_hat_opp + gap, sum_top2,
# ratio_second_max, entropy, gini_conc. Labels: type 2/7/11 = HS (1), type 100 = N/A (0).
build_features_sup <- function(data, eps = 1e-6) {
  x1 <- as.numeric(data[["sum_pihat_max_5"]])
  x2 <- as.numeric(data[["sum_pihat_remaining_max_5"]])
  x3 <- as.numeric(data[["sum_pihat_remaining_min_5"]])
  x4 <- as.numeric(data[["sum_pihat_opposite_5"]])
  gap_max_second   <- x1 - x2
  sum_top2         <- x1 + x2
  ratio_second_max <- x2 / (x1 + eps)
  tot <- x1 + x2 + x3 + x4
  p1 <- x1 / (tot + eps); p2 <- x2 / (tot + eps); p3 <- x3 / (tot + eps); p4 <- x4 / (tot + eps)
  entropy <- -(p1*log(p1+eps) + p2*log(p2+eps) + p3*log(p3+eps) + p4*log(p4+eps))
  gini_conc <- 1 - (p1^2 + p2^2 + p3^2 + p4^2)
  data.frame(pi_hat_max = x1, pi_hat_R_max = x2, pi_hat_R_min = x3, pi_hat_opp = x4,
             gap_max_second = gap_max_second, sum_top2 = sum_top2, ratio_second_max = ratio_second_max,
             entropy = entropy, gini_conc = gini_conc)
}
# Map type to binary: HS (2,7,11) -> 1, N/A (100) -> 0, else NA.
get_binary_label_sup <- function(type) {
  y <- rep(NA_integer_, length(type))
  y[type %in% c(2, 7, 11)] <- 1L
  y[type == 100] <- 0L
  y
}
# Fit logistic (or glmnet ridge) on labeled rows; optional class weights.
fit_logistic_sup <- function(X, y, lambda_ridge = NULL, use_weights = TRUE) {
  idx <- !is.na(y)
  X_fit <- X[idx, , drop = FALSE]
  y_fit <- y[idx]
  # Class weights: inverse frequency so HS (minority) counts more
  wts <- NULL
  if (use_weights) {
    n1 <- sum(y_fit == 1); n0 <- sum(y_fit == 0)
    if (n1 > 0 && n0 > 0) wts <- ifelse(y_fit == 1, 1 / n1, 1 / n0)
  }
  if (is.null(lambda_ridge) || lambda_ridge == 0) {
    df <- data.frame(y = factor(y_fit, levels = c(0, 1)), X_fit)
    fit <- if (is.null(wts)) glm(y ~ ., data = df, family = binomial) else glm(y ~ ., data = df, family = binomial, weights = wts)
    return(list(model = fit, type = "glm", lambda = NULL))
  }
  if (requireNamespace("glmnet", quietly = TRUE)) {
    X_mat <- as.matrix(X_fit)
    fit <- if (is.null(wts)) glmnet::glmnet(X_mat, y_fit, family = "binomial", alpha = 0, lambda = lambda_ridge) else glmnet::glmnet(X_mat, y_fit, family = "binomial", alpha = 0, lambda = lambda_ridge, weights = wts)
    return(list(model = fit, type = "glmnet", lambda = lambda_ridge, X_scale = list(mean = colMeans(X_mat), sd = apply(X_mat, 2, sd))))
  }
  df <- data.frame(y = factor(y_fit, levels = c(0, 1)), X_fit)
  fit <- if (is.null(wts)) glm(y ~ ., data = df, family = binomial) else glm(y ~ ., data = df, family = binomial, weights = wts)
  list(model = fit, type = "glm", lambda = NULL)
}
predict_logistic_sup <- function(fit, X_new) {
  if (fit$type == "glm") return(as.numeric(predict(fit$model, newdata = X_new, type = "response")))
  if (fit$type == "glmnet") {
    X_mat <- as.matrix(X_new)
    for (j in seq_len(ncol(X_mat))) if (fit$X_scale$sd[j] > 0) X_mat[, j] <- (X_mat[, j] - fit$X_scale$mean[j]) / fit$X_scale$sd[j]
    return(as.numeric(predict(fit$model, newx = X_mat, type = "response")))
  }
  NA
}
calibrate_platt_sup <- function(pred, y) {
  idx <- !is.na(y) & is.finite(pred) & pred > 0 & pred < 1
  if (sum(idx) < 10) return(function(p) p)
  logit <- qlogis(pred[idx]); logit[is.infinite(logit)] <- NA
  ok <- is.finite(logit); if (sum(ok) < 10) return(function(p) p)
  cal <- glm(y[idx][ok] ~ logit[ok], family = binomial)
  function(p) { lp <- qlogis(p); lp[is.infinite(lp)] <- max(logit[ok], na.rm = TRUE); as.numeric(predict(cal, newdata = data.frame(logit = lp), type = "response")) }
}
calibrate_isotonic_sup <- function(pred, y) {
  idx <- !is.na(y) & is.finite(pred) & pred >= 0 & pred <= 1
  if (sum(idx) < 10) return(function(p) p)
  ord <- order(pred[idx]); x_iso <- pred[idx][ord]; y_iso <- y[idx][ord]
  iso <- stats::isoreg(x_iso, y_iso); x_vals <- iso$x; y_vals <- iso$yf
  function(p) { p <- pmax(0, pmin(1, p)); i <- findInterval(p, x_vals, left.open = FALSE); i <- pmin(pmax(i, 1), length(y_vals)); y_vals[i] }
}
evaluate_cv_sup <- function(X, y, K = 5, method = "glm", lambda_ridge = 0, calibrate = "platt", use_weights = TRUE) {
  idx <- !is.na(y); X <- X[idx, , drop = FALSE]; y <- y[idx]; n <- nrow(X)
  if (n < 2*K) return(NULL)
  set.seed(42); folds <- sample(rep(seq_len(K), length.out = n))
  pred_cv <- numeric(n); pred_cv[] <- NA_real_
  for (k in seq_len(K)) {
    i_tr <- folds != k; i_te <- folds == k
    fit <- fit_logistic_sup(X[i_tr, , drop = FALSE], y[i_tr], lambda_ridge = if (method == "glmnet") lambda_ridge else NULL, use_weights = use_weights)
    pred_cv[i_te] <- predict_logistic_sup(fit, X[i_te, , drop = FALSE])
  }
  if (calibrate == "platt") { cal_fun <- calibrate_platt_sup(pred_cv, y); pred_cv <- cal_fun(pred_cv) }
  if (calibrate == "isotonic") { cal_fun <- calibrate_isotonic_sup(pred_cv, y); pred_cv <- cal_fun(pred_cv) }
  pred_cv <- pmax(0, pmin(1, pred_cv)); pred_class <- as.integer(pred_cv >= 0.5)
  idx_ok <- !is.na(pred_class)
  TP <- sum(y == 1 & pred_class == 1, na.rm = TRUE); TN <- sum(y == 0 & pred_class == 0, na.rm = TRUE); FP <- sum(y == 0 & pred_class == 1, na.rm = TRUE); FN <- sum(y == 1 & pred_class == 0, na.rm = TRUE)
  tot_pos <- TP + FN; tot_neg <- TN + FP
  sensitivity <- if (is.finite(tot_pos) && tot_pos > 0) TP / tot_pos else NA_real_
  specificity <- if (is.finite(tot_neg) && tot_neg > 0) TN / tot_neg else NA_real_
  acc <- if (sum(idx_ok)) mean(pred_class[idx_ok] == y[idx_ok]) else NA_real_
  bins <- cut(pred_cv[idx_ok], breaks = seq(0, 1, length.out = 11), include.lowest = TRUE)
  ece <- if (sum(idx_ok)) sum(sapply(levels(bins), function(b) { i <- which(bins == b); if (length(i)) length(i)/sum(idx_ok) * abs(mean(pred_cv[idx_ok][i]) - mean(y[idx_ok][i])) else 0 })) else NA_real_
  list(pred_cv = pred_cv, y = y, accuracy = acc, sensitivity = sensitivity, specificity = specificity, TP = TP, TN = TN, FP = FP, FN = FN, ECE = ece)
}
# Run full pipeline: fit on labeled complete rows, predict on all complete rows,
# optionally calibrate (Platt/isotonic), run K-fold CV, optionally write TSV.
run_supervised_classifier_inline <- function(data, output_file = NULL, use_ridge = FALSE, lambda_ridge = 1e-4, calibrate = "platt", K_cv = 5, use_weights = TRUE) {
  if (!"type" %in% names(data)) stop("Column 'type' required")
  X <- build_features_sup(data); y <- get_binary_label_sup(data$type)
  idx_labeled <- !is.na(y); ok <- complete.cases(X) & idx_labeled
  X_fit <- X[ok, , drop = FALSE]; y_fit <- y[ok]
  if (nrow(X_fit) < 10) { warning("Too few labeled rows for supervised fit"); return(list(ok_all = rep(FALSE, nrow(data)), prob_all = rep(NA_real_, nrow(data)), cv = NULL)) }
  fit_final <- fit_logistic_sup(X_fit, y_fit, lambda_ridge = if (use_ridge) lambda_ridge else NULL, use_weights = use_weights)
  ok_all <- complete.cases(X); X_all <- X[ok_all, , drop = FALSE]
  prob_all <- predict_logistic_sup(fit_final, X_all)
  if (calibrate == "platt") { cal_fun <- calibrate_platt_sup(prob_all[ok_all & idx_labeled], y[ok_all & idx_labeled]); prob_all <- cal_fun(prob_all) }
  if (calibrate == "isotonic") { cal_fun <- calibrate_isotonic_sup(prob_all[ok_all & idx_labeled], y[ok_all & idx_labeled]); prob_all <- cal_fun(prob_all) }
  prob_all <- pmax(0, pmin(1, prob_all))
  if (!is.null(output_file) && nzchar(output_file)) { out <- data[ok_all, ]; out$P_HS_supervised <- prob_all; out$pred_class_supervised <- as.integer(prob_all >= 0.5); write.table(out, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE); cat("Wrote supervised predictions to:", output_file, "\n") }
  cv <- evaluate_cv_sup(X, y, K = K_cv, method = if (use_ridge) "glmnet" else "glm", lambda_ridge = lambda_ridge, calibrate = calibrate, use_weights = use_weights)
  if (!is.null(cv)) { cat("\n--- Supervised classifier CV ---\n"); cat("Accuracy:", round(cv$accuracy, 4), "Sens:", round(cv$sensitivity, 4), "Spec:", round(cv$specificity, 4), "| ECE:", round(cv$ECE, 4), "\n") }
  invisible(list(ok_all = ok_all, prob_all = prob_all, pred_class_supervised = as.integer(prob_all >= 0.5), cv = cv))
}

res_sup <- NULL
data_for_supervised <- data[idx_valid, ]
if (nrow(data_for_supervised) > 0 && all(c("sum_pihat_max_5", "sum_pihat_remaining_max_5", "sum_pihat_remaining_min_5", "sum_pihat_opposite_5") %in% names(data_for_supervised))) {
  res_sup <- run_supervised_classifier_inline(data_for_supervised, output_file = NULL, use_ridge = TRUE, lambda_ridge = 1e-4, calibrate = "platt", K_cv = 5)
}

# ============================================================================
# 6c. Multivariate 2-component Gaussian mixture (4-D) → P(H-S|x)
# ============================================================================
# Fit EM on ALL pairs with complete 4-D (unsupervised). Use labeled subset only
# to identify which component = HS (comp_HS) and to evaluate accuracy/sens/spec.
# Model (π̂_Max, π̂_R_Max, π̂_R_Min, π̂_Opp) with 2-component MVN mixture; use posterior P(HS|x).
#
# Formula for P(H-S) from the 4-D vector x = (x1,x2,x3,x4) (sum_pihat_max_5, sum_pihat_remaining_max_5, sum_pihat_remaining_min_5, sum_pihat_opposite_5):
#   d1(x) = dmvnorm(x, mean = mu1, sigma = Sigma1)   # density of component 1
#   d2(x) = dmvnorm(x, mean = mu2, sigma = Sigma2)   # density of component 2
#   post1(x) = lambda1 * d1(x) / (lambda1 * d1(x) + lambda2 * d2(x))
#   P(H-S|x) = post1(x) if component 1 = HS, else 1 - post1(x)
# (then clamped to [0,1].)
# Saved parameters (see GMM_mixture_4D_params.rds) allow recomputing P(H-S) from any 4-D vector.
# Example: p <- readRDS("/pl/active/KellerLab/Emmanuel/NAvHS/GMM_mixture_4D_params.rds")
#   x <- c(sum_pihat_max_5, sum_pihat_remaining_max_5, sum_pihat_remaining_min_5, sum_pihat_opposite_5)
#   d1 <- mvtnorm::dmvnorm(x, p$mu[[1]], p$sigma[[1]]); d2 <- mvtnorm::dmvnorm(x, p$mu[[2]], p$sigma[[2]])
#   post1 <- p$lambda[1]*d1 / (p$lambda[1]*d1 + p$lambda[2]*d2)
#   P_HS <- if (p$comp_HS == 1) post1 else (1 - post1); P_HS <- pmax(0, pmin(1, P_HS))
res_mix <- NULL
if (nrow(data_for_supervised) > 0 && all(c("sum_pihat_max_5", "sum_pihat_remaining_max_5", "sum_pihat_remaining_min_5", "sum_pihat_opposite_5") %in% names(data_for_supervised))) {
  # 4-D feature matrix and binary labels (HS=1, N/A=0).
  X_mv <- as.matrix(data_for_supervised[, c("sum_pihat_max_5", "sum_pihat_remaining_max_5", "sum_pihat_remaining_min_5", "sum_pihat_opposite_5")])
  y_mv <- get_binary_label_sup(data_for_supervised$type)
  idx_mv <- !is.na(y_mv) & complete.cases(X_mv)
  idx_complete_4d <- complete.cases(X_mv)
  n_complete_4d <- sum(idx_complete_4d)
  n_lab <- sum(idx_mv)
  n_hs <- sum(y_mv[idx_mv] == 1)
  n_control <- sum(y_mv[idx_mv] == 0)
  X_fit_mv <- X_mv[idx_complete_4d, , drop = FALSE]
  cat("\n--- Multivariate mixture (4-D) setup ---\n")
  cat("EM fit on ALL pairs with complete 4-D (unsupervised). n_complete_4d =", n_complete_4d, "| Labeled subset (for comp_HS): n_lab =", n_lab, "| HS:", n_hs, "| Control:", n_control, "\n")
  if (n_lab < 20 || n_hs < 5 || n_control < 5) {
    cat("Skipping multivariate mixture: need n_lab >= 20 and >= 5 per class (HS and Control).\n")
  } else if (!requireNamespace("mixtools", quietly = TRUE)) {
    cat("Skipping multivariate mixture: package mixtools not installed.\n")
  } else {
    tryCatch({
      fit_mv <- mixtools::mvnormalmixEM(X_fit_mv, k = 2, verb = FALSE, maxit = 1000)
      # Which component is HS? Use labeled subset posteriors only
      idx_lab_in_fit <- idx_mv[idx_complete_4d]
      post_lab <- fit_mv$posterior[idx_lab_in_fit, , drop = FALSE]
      y_fit_lab <- y_mv[idx_complete_4d][idx_lab_in_fit]
      hs_in_1 <- sum(post_lab[y_fit_lab == 1, 1]); hs_in_2 <- sum(post_lab[y_fit_lab == 1, 2])
      comp_HS <- if (hs_in_1 >= hs_in_2) 1 else 2
      # P(HS|x) for all rows with complete 4-D
      ok_all_mv <- complete.cases(X_mv)
      X_all_mv <- X_mv[ok_all_mv, , drop = FALSE]
      if (!requireNamespace("mvtnorm", quietly = TRUE)) {
        cat("Multivariate mixture: mvtnorm required for posterior; install.packages('mvtnorm'). Skipping.\n")
      } else {
        d1 <- mvtnorm::dmvnorm(X_all_mv, fit_mv$mu[[1]], fit_mv$sigma[[1]])
        d2 <- mvtnorm::dmvnorm(X_all_mv, fit_mv$mu[[2]], fit_mv$sigma[[2]])
        post1 <- fit_mv$lambda[1] * d1 / (fit_mv$lambda[1] * d1 + fit_mv$lambda[2] * d2)
        P_HS_mix <- if (comp_HS == 1) post1 else (1 - post1)
        P_HS_mix <- pmax(0, pmin(1, P_HS_mix))
        pred_mix <- as.integer(P_HS_mix >= 0.5)
        # Accuracy on labeled
        idx_lab_mv <- !is.na(y_mv) & ok_all_mv
        pred_lab <- pred_mix[ok_all_mv][idx_lab_mv[ok_all_mv]]
        y_lab <- y_mv[idx_lab_mv]
        if (length(pred_lab) == length(y_lab) && length(y_lab) > 0) {
          acc_mv <- mean(pred_lab == y_lab)
          sens_mv <- if (sum(y_lab == 1) > 0) mean(pred_lab[y_lab == 1] == 1) else NA
          spec_mv <- if (sum(y_lab == 0) > 0) mean(pred_lab[y_lab == 0] == 0) else NA
          cat("--- Multivariate 2-component mixture (4-D) ---\n")
          cat("Accuracy (labeled):", round(acc_mv, 4), "| Sens:", round(sens_mv, 4), "| Spec:", round(spec_mv, 4), "\n")
        }
        res_mix <- list(ok_all = ok_all_mv, prob_all = P_HS_mix, pred_class_mixture = pred_mix)
        # Save mixture parameters so P(H-S) can be recomputed from any 4-D vector
        mixture_params_file <- "/pl/active/KellerLab/Emmanuel/NAvHS/GMM_mixture_4D_params.rds"
        mixture_params <- list(
          lambda = fit_mv$lambda,
          mu = fit_mv$mu,
          sigma = fit_mv$sigma,
          comp_HS = comp_HS,
          var_names = c("sum_pihat_max_5", "sum_pihat_remaining_max_5", "sum_pihat_remaining_min_5", "sum_pihat_opposite_5")
        )
        tryCatch({
          saveRDS(mixture_params, mixture_params_file)
          cat("Mixture parameters saved to:", mixture_params_file, "(use to compute P(H-S) from 4-D vector)\n")
        }, error = function(e) cat("Could not save mixture params:", e$message, "\n"))
      }
    }, error = function(e) { cat("Multivariate mixture fit failed:", e$message, "\n") })
  }
}

# --- Supervised classifier: PDF with calibration, P(HS) by group, ROC ---
supervised_pdf <- "/pl/active/KellerLab/Emmanuel/NAvHS/GMM_plot_supervised.pdf"
if (!is.null(res_sup) && any(res_sup$ok_all)) {
  pdf(supervised_pdf, width = 10, height = 4)
  par(mfrow = c(1, 3), mar = c(4, 4, 3, 1), cex.lab = 1.35, cex.axis = 1.2, cex.main = 1.25)
  P_HS <- res_sup$prob_all
  y_sup <- get_binary_label_sup(data_for_supervised$type)
  idx_lab <- !is.na(y_sup)

  # 1. Calibration (reliability) diagram
  n_bins <- 10
  br <- seq(0, 1, length.out = n_bins + 1)
  bins <- cut(P_HS[idx_lab], breaks = br, include.lowest = TRUE)
  pred_bin <- sapply(levels(bins), function(b) { i <- which(bins == b); if (length(i)) mean(P_HS[idx_lab][i]) else NA })
  obs_bin  <- sapply(levels(bins), function(b) { i <- which(bins == b); if (length(i)) mean(y_sup[idx_lab][i]) else NA })
  plot(c(0, 1), c(0, 1), type = "n", xlab = "Mean predicted P(HS)", ylab = "Observed fraction HS", main = "Calibration (reliability)")
  abline(0, 1, col = "gray", lty = 2)
  points(pred_bin, obs_bin, pch = 19, col = "darkblue", cex = 1.2)
  grid(col = "lightgray", lty = "dotted")

  # 2. Histogram of P(HS) by group (HS vs N/A)
  idx_hs <- y_sup == 1; idx_na <- y_sup == 0
  h_na <- hist(P_HS[idx_na], breaks = seq(0, 1, 0.05), plot = FALSE)
  h_hs <- hist(P_HS[idx_hs], breaks = seq(0, 1, 0.05), plot = FALSE)
  ymax <- max(c(h_na$counts, h_hs$counts, 1), na.rm = TRUE) * 1.2
  hist(P_HS[idx_na], breaks = seq(0, 1, 0.05), col = "lightskyblue1", border = "gray", main = "P(HS) by group", xlab = "P(H-S)", xlim = c(0, 1), ylim = c(0, ymax))
  hist(P_HS[idx_hs], breaks = seq(0, 1, 0.05), col = "#FF80FF", border = "gray", add = TRUE)
  legend("topright", legend = c("N/A", "HS"), fill = c("lightskyblue1", "#FF80FF"), bty = "n", cex = 1.3)

  # 3. ROC curve
  th <- sort(unique(c(0, P_HS[idx_lab], 1)))
  sens <- spec <- numeric(length(th))
  for (i in seq_along(th)) {
    pred_pos <- P_HS >= th[i]
    TP <- sum(y_sup[idx_lab] == 1 & pred_pos[idx_lab], na.rm = TRUE); FN <- sum(y_sup[idx_lab] == 1 & !pred_pos[idx_lab], na.rm = TRUE)
    TN <- sum(y_sup[idx_lab] == 0 & !pred_pos[idx_lab], na.rm = TRUE); FP <- sum(y_sup[idx_lab] == 0 & pred_pos[idx_lab], na.rm = TRUE)
    tot_pos <- TP + FN; tot_neg <- TN + FP
    sens[i] <- if (is.finite(tot_pos) && tot_pos > 0) TP / tot_pos else 0
    spec[i] <- if (is.finite(tot_neg) && tot_neg > 0) TN / tot_neg else 0
  }
  plot(1 - spec, sens, type = "l", lwd = 2, col = "darkblue", xlab = "1 - Specificity", ylab = "Sensitivity", main = "ROC (supervised)")
  abline(0, 1, col = "gray", lty = 2)
  grid(col = "lightgray", lty = "dotted")
  dev.off()
  cat("Supervised classifier plots saved to:", supervised_pdf, "\n")
}

# ============================================================================
# 6d. Multivariate mixture: TSV (P(HS) by group), cutoff choice, PDFs
# ============================================================================
# Build table with P_HS_mixture, group (HS/N/A/Other), FN/FP at chosen cutoff.
# Cutoff: if any threshold gives sensitivity=1, pick that with max specificity;
# else use threshold that maximizes accuracy. Re-write TSV with cut_show.
mixture_pdf <- "/pl/active/KellerLab/Emmanuel/NAvHS/GMM_plot_mixture.pdf"
mixture_tsv <- "/pl/active/KellerLab/Emmanuel/NAvHS/GMM_mixture_PHS_by_group.txt"
if (!is.null(res_mix) && any(res_mix$ok_all)) {
  y_mix <- get_binary_label_sup(data_for_supervised$type)[res_mix$ok_all]
  idx_lab_mix <- !is.na(y_mix)
  P_HS_mix_all <- res_mix$prob_all
  type_mix <- data_for_supervised$type[res_mix$ok_all]
  group_mix <- ifelse(type_mix %in% c(2, 7, 11), "HS", ifelse(type_mix == 100, "N/A", "Other"))
  out_mix <- data.frame(
    type = type_mix,
    group = group_mix,
    P_HS_mixture = P_HS_mix_all,
    pred_class_mixture = res_mix$pred_class_mixture
  )
  if ("Pihat" %in% names(data_for_supervised)) out_mix$pihat <- as.numeric(data_for_supervised$Pihat[res_mix$ok_all])
  if ("numtrio" %in% names(data_for_supervised)) out_mix$numtrio <- data_for_supervised$numtrio[res_mix$ok_all]
  if ("UKBID1" %in% names(data_for_supervised)) out_mix$UKBID1 <- data_for_supervised$UKBID1[res_mix$ok_all]
  if ("UKBID2" %in% names(data_for_supervised)) out_mix$UKBID2 <- data_for_supervised$UKBID2[res_mix$ok_all]
  if ("ID1" %in% names(data_for_supervised)) out_mix$ID1 <- data_for_supervised$ID1[res_mix$ok_all]
  if ("ID2" %in% names(data_for_supervised)) out_mix$ID2 <- data_for_supervised$ID2[res_mix$ok_all]
  # FN/FP with temporary cutoff; will overwrite with cut_show after cutoff computation
  cut_tsv <- 0.999996
  pred_at_cut <- P_HS_mix_all >= cut_tsv
  is_FN <- (out_mix$group == "HS" & !pred_at_cut)
  is_FP <- (out_mix$group == "N/A" & pred_at_cut)
  out_mix$FN_FP <- ifelse(is_FN, "FN", ifelse(is_FP, "FP", ""))
  ord <- order(!is_FN, !is_FP)
  out_mix <- out_mix[ord, , drop = FALSE]
  write.table(out_mix, file = mixture_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("P(HS) by group written to:", mixture_tsv, "(", nrow(out_mix), "rows)\n")

  cut_show <- 0.999996
  # Compute best cutoff on labeled N/A vs HS: prefer sensitivity=1 + max specificity.
  if (any(idx_lab_mix)) {
    P_HS_lab <- P_HS_mix_all[idx_lab_mix]
    y_lab_mix <- y_mix[idx_lab_mix]
    p_na <- P_HS_lab[y_lab_mix == 0]
    p_hs <- P_HS_lab[y_lab_mix == 1]
    max_na <- max(p_na)
    min_hs <- min(p_hs)
    min_na <- min(p_na)
    max_hs <- max(p_hs)
    if (max_na < min_hs) {
      cutoff_mid <- (max_na + min_hs) / 2
      pred_mid <- P_HS_lab >= cutoff_mid
      TP_mid <- sum(y_lab_mix == 1 & pred_mid, na.rm = TRUE)
      TN_mid <- sum(y_lab_mix == 0 & !pred_mid, na.rm = TRUE)
      FP_mid <- sum(y_lab_mix == 0 & pred_mid, na.rm = TRUE)
      FN_mid <- sum(y_lab_mix == 1 & !pred_mid, na.rm = TRUE)
      cat("--- P(HS) perfect cutoff (labeled N/A vs HS) ---\n")
      cat("Yes. Any cutoff in (", round(max_na, 6), ", ", round(min_hs, 6), ") discriminates perfectly.\n")
      cat("Suggested cutoff (midpoint):", round(cutoff_mid, 6), "\n")
      cat("At midpoint cutoff: TP =", TP_mid, " TN =", TN_mid, " FP =", FP_mid, " FN =", FN_mid, "\n")
    } else {
      pred_05 <- P_HS_lab >= 0.5
      TP_05 <- sum(y_lab_mix == 1 & pred_05, na.rm = TRUE)
      TN_05 <- sum(y_lab_mix == 0 & !pred_05, na.rm = TRUE)
      FP_05 <- sum(y_lab_mix == 0 & pred_05, na.rm = TRUE)
      FN_05 <- sum(y_lab_mix == 1 & !pred_05, na.rm = TRUE)
      cat("--- P(HS) perfect cutoff (labeled N/A vs HS) ---\n")
      cat("No. Overlap: N/A P(HS) in [", round(min_na, 4), ", ", round(max_na, 4), "], HS P(HS) in [", round(min_hs, 4), ", ", round(max_hs, 4), "].\n")
      cat("At cutoff 0.5: TP =", TP_05, " TN =", TN_05, " FP =", FP_05, " FN =", FN_05, "\n")
    }
    # Scan all unique P(HS) values as thresholds; pick best accuracy or (if sens=1) max specificity.
    th_acc <- sort(unique(c(0, P_HS_lab, 1)))
    acc_vec <- sens_vec <- spec_vec <- numeric(length(th_acc))
    tp_vec <- tn_vec <- fp_vec <- fn_vec <- integer(length(th_acc))
    for (i in seq_along(th_acc)) {
      pred_1 <- P_HS_lab >= th_acc[i]
      TP <- sum(y_lab_mix == 1 & pred_1, na.rm = TRUE)
      TN <- sum(y_lab_mix == 0 & !pred_1, na.rm = TRUE)
      FP <- sum(y_lab_mix == 0 & pred_1, na.rm = TRUE)
      FN <- sum(y_lab_mix == 1 & !pred_1, na.rm = TRUE)
      tp_vec[i] <- TP; tn_vec[i] <- TN; fp_vec[i] <- FP; fn_vec[i] <- FN
      n_correct <- TP + TN
      acc_vec[i] <- n_correct / length(y_lab_mix)
      tot_pos <- TP + FN
      tot_neg <- TN + FP
      sens_vec[i] <- if (is.finite(tot_pos) && tot_pos > 0) TP / tot_pos else 0
      spec_vec[i] <- if (is.finite(tot_neg) && tot_neg > 0) TN / tot_neg else 0
    }
    best_i <- which.max(acc_vec)
    best_cutoff <- th_acc[best_i]
    idx_sens1 <- (sens_vec >= 1 - 1e-9)
    if (any(idx_sens1)) {
      spec_at_sens1 <- spec_vec
      spec_at_sens1[!idx_sens1] <- -1
      best_sens1_i <- which.max(spec_at_sens1)
      cut_show <- th_acc[best_sens1_i]
      cat("--- Cutoff with sensitivity 1 and maximum specificity ---\n")
      cat("Cutoff:", round(cut_show, 6), "| Sens: 1 | Spec:", round(spec_vec[best_sens1_i], 4), "| (chosen so Sens=1, Spec max)\n")
    } else {
      cut_show <- best_cutoff
    }
    cat("--- Best P(HS) cutoff for accuracy (labeled N/A vs HS) ---\n")
    cat("Cutoff:", round(best_cutoff, 6), "| Accuracy:", round(acc_vec[best_i], 4), "| Sens:", round(sens_vec[best_i], 4), "| Spec:", round(spec_vec[best_i], 4), "\n")
    cat("TP =", tp_vec[best_i], " TN =", tn_vec[best_i], " FP =", fp_vec[best_i], " FN =", fn_vec[best_i], "\n")
    cat("(Used cutoff for plots/TSV:", round(cut_show, 6), if (any(idx_sens1)) "[Sens=1, max Spec]" else "[max Accuracy]", ")\n")
    cat("\n--- Why some N/A have high P(HS), and why HS are almost always higher ---\n")
    cat("The 4-D mixture fits two Gaussians in (pi_Max, pi_R_Max, pi_R_Min, pi_Opp). Some N/A pairs\n")
    cat("fall close to the HS component (e.g. haplotype sharing by chance or structure), so they get\n")
    cat("high P(HS). HS pairs are concentrated in the HS component, so their P(HS) is almost always\n")
    cat("high; the overlap is why a high cutoff (e.g. ", round(best_cutoff, 4), ") is needed to maximize accuracy.\n")
    cat("\n--- Why use this cutoff? ---\n")
    cat("This cutoff was chosen to maximize accuracy on labeled N/A vs HS. See mixture PDF page\n")
    cat("\"Justification: accuracy is maximized at chosen cutoff\" for the curve.\n\n")
    # Table at selected cutoffs: 0.5, best, 0.99, 0.999
    cutoffs_out <- unique(c(0.5, best_cutoff, 0.99, 0.999))
    cutoffs_out <- sort(cutoffs_out)
    cat("--- P(HS) cutoffs: TP, TN, FP, FN, Acc, Sens, Spec ---\n")
    for (co in cutoffs_out) {
      i_co <- which(th_acc >= co)[1]
      if (is.na(i_co)) i_co <- length(th_acc)
      pred_co <- P_HS_lab >= co
      TP_co <- sum(y_lab_mix == 1 & pred_co, na.rm = TRUE)
      TN_co <- sum(y_lab_mix == 0 & !pred_co, na.rm = TRUE)
      FP_co <- sum(y_lab_mix == 0 & pred_co, na.rm = TRUE)
      FN_co <- sum(y_lab_mix == 1 & !pred_co, na.rm = TRUE)
      acc_co <- (TP_co + TN_co) / length(y_lab_mix)
      sens_co <- if (TP_co + FN_co > 0) TP_co / (TP_co + FN_co) else 0
      spec_co <- if (TN_co + FP_co > 0) TN_co / (TN_co + FP_co) else 0
      cat("Cutoff", round(co, 6), ": TP =", TP_co, " TN =", TN_co, " FP =", FP_co, " FN =", FN_co, "| Acc =", round(acc_co, 4), " Sens =", round(sens_co, 4), " Spec =", round(spec_co, 4), "\n")
    }
  }

  # Apply chosen cutoff to FN/FP and re-write TSV (rows ordered: FN, then FP, then rest).
  pred_at_cut <- out_mix$P_HS_mixture >= cut_show
  is_FN <- (out_mix$group == "HS" & !pred_at_cut)
  is_FP <- (out_mix$group == "N/A" & pred_at_cut)
  out_mix$FN_FP <- ifelse(is_FN, "FN", ifelse(is_FP, "FP", ""))
  ord <- order(!is_FN, !is_FP)
  out_mix <- out_mix[ord, , drop = FALSE]
  write.table(out_mix, file = mixture_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("P(HS) by group (cutoff ", round(cut_show, 6), ") written to:", mixture_tsv, "\n", sep = "")

  # Multi-page PDF: P(HS) histograms, by-group breakdown, calibration, ROC, justification.
  pdf(mixture_pdf, width = 10, height = 8)
  br <- seq(0, 1, 0.05)

  # Page 1: histogram of P(HS) for all samples with 4-D data (single color, cutoff line).
  par(mfrow = c(1, 1), mar = c(5, 5, 4, 2), cex.lab = 1.35, cex.axis = 1.2, cex.main = 1.25)
  hist(P_HS_mix_all, breaks = br, col = "lightgray", border = "darkgray", main = "P(HS) (4-D mixture, all samples)", xlab = "P(H-S)", xlim = c(0, 1), ylab = "Frequency")
  abline(v = cut_show, col = "red", lty = 2, lwd = 1.5)
  legend("topright", legend = paste("n =", length(P_HS_mix_all)), bty = "n", cex = 1.3)

  # Page 2: P(HS) histogram color-coded by GT class (HS, N/A, Other); last bin = P(H-S) >= cut_show, same width as others
  par(mfrow = c(1, 1), mar = c(5, 5, 4, 2), cex.lab = 1.35, cex.axis = 1.2, cex.main = 1.25)
  idx_na_gt <- group_mix == "N/A"
  idx_hs_gt <- group_mix == "HS"
  idx_other_gt <- group_mix == "Other"
  br_p2 <- c(seq(0, 0.95, 0.05), cut_show, 1)
  h_all_p2 <- hist(P_HS_mix_all, breaks = br_p2, plot = FALSE)
  h_na_p2 <- if (sum(idx_na_gt) > 0) hist(P_HS_mix_all[idx_na_gt], breaks = br_p2, plot = FALSE) else list(counts = rep(0, length(br_p2) - 1))
  h_hs_p2 <- if (sum(idx_hs_gt) > 0) hist(P_HS_mix_all[idx_hs_gt], breaks = br_p2, plot = FALSE) else list(counts = rep(0, length(br_p2) - 1))
  c_all <- c(h_all_p2$counts[1:18], h_all_p2$counts[19] + h_all_p2$counts[20], h_all_p2$counts[21])
  c_na <- c(h_na_p2$counts[1:18], h_na_p2$counts[19] + h_na_p2$counts[20], h_na_p2$counts[21])
  c_hs <- c(h_hs_p2$counts[1:18], h_hs_p2$counts[19] + h_hs_p2$counts[20], h_hs_p2$counts[21])
  n_display <- 20
  ymax_gt <- max(c(c_all, c_na, c_hs, 1), na.rm = TRUE) * 1.2
  n_pred_na <- sum(P_HS_mix_all < cut_show)
  n_pred_hs <- sum(P_HS_mix_all >= cut_show)
  plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, ymax_gt), xlab = "P(H-S)", ylab = "Frequency", main = "")
  x_left <- seq(0, 0.95, 0.05)
  x_right <- seq(0.05, 1, 0.05)
  for (i in 1:n_display) {
    bar_col <- if (x_right[i] <= cut_show) rgb(0.85, 0.9, 1) else rgb(1, 0.92, 0.96)
    if (c_all[i] > 0) rect(x_left[i], 0, x_right[i], c_all[i], col = bar_col, border = "gray")
  }
  for (i in 1:n_display) {
    if (c_na[i] > 0) rect(x_left[i], 0, x_right[i], c_na[i], col = rgb(0.53, 0.81, 0.98, 0.85), border = "gray")
  }
  for (i in 1:n_display) {
    if (c_hs[i] > 0) rect(x_left[i], 0, x_right[i], c_hs[i], col = rgb(1, 0.5, 1, 0.85), border = "gray")
  }
  abline(v = 0.95, col = "darkgreen", lty = 2, lwd = 1.5)
  text(0.95, ymax_gt * 0.5, paste0("last bin: pair with\n P(H-S) >= ", round(cut_show, 5), "\nPredicted as H-S"), col = "darkgreen", adj = c(0.5, 0.5), cex = 1, xpd = TRUE)
  pred_cut_leg <- P_HS_mix_all >= cut_show
  TN_leg <- sum(idx_na_gt & !pred_cut_leg)
  FP_leg <- sum(idx_na_gt & pred_cut_leg)
  TP_leg <- sum(idx_hs_gt & pred_cut_leg)
  FN_leg <- sum(idx_hs_gt & !pred_cut_leg)
  n_other_pred_na <- sum(idx_other_gt & !pred_cut_leg)
  n_other_pred_hs <- sum(idx_other_gt & pred_cut_leg)
  leg_text <- c(
    paste0("N/N-A (n=", sum(idx_na_gt), ")"),
    paste0("H-S (n=", sum(idx_hs_gt), ")"),
    "",
    paste0("N/N-A predicted N/N-A (n=", TN_leg, ")"),
    paste0("N/N-A predicted H-S (n=", FP_leg, ")"),
    paste0("H-S predicted H-S (n=", TP_leg, ")"),
    paste0("H-S predicted N/N-A (n=", FN_leg, ")"),
    "",
    paste0("Other Predicted N/N-A (n=", n_other_pred_na, ")"),
    paste0("Other Predicted H-S (n=", n_other_pred_hs, ")")
  )
  leg_fill <- c(rgb(0.53, 0.81, 0.98, 0.9), rgb(1, 0.5, 1, 0.9), NA, rgb(0.53, 0.81, 0.98, 0.9), rgb(0.53, 0.81, 0.98, 0.9), rgb(1, 0.5, 1, 0.9), rgb(1, 0.5, 1, 0.9), NA, rgb(0.85, 0.9, 1), rgb(1, 0.92, 0.96))
  leg_border <- rep(NA, 10)
  leg_lty <- c(NA, NA, 1, NA, NA, NA, NA, 1, NA, NA)
  leg_lwd <- c(NA, NA, 2, NA, NA, NA, NA, 2, NA, NA)
  leg_col <- c(rep(NA, 2), "gray", rep(NA, 4), "gray", rep(NA, 2))
  legend("top", legend = leg_text, fill = leg_fill, border = leg_border, lty = leg_lty, lwd = leg_lwd, col = leg_col, bty = "n", cex = 1.2, xpd = TRUE)

  # Page 2 by itself: save to separate PDF (same 20 equal-width bars, green line at 0.95)
  mixture_pdf_page2 <- "/pl/active/KellerLab/Emmanuel/NAvHS/GMM_plot_mixture_page2.pdf"
  pdf(mixture_pdf_page2, width = 10, height = 8)
  par(mfrow = c(1, 1), mar = c(5, 5, 4, 2), cex.lab = 1.35, cex.axis = 1.2, cex.main = 1.25)
  plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, ymax_gt), xlab = "P(H-S)", ylab = "Frequency", main = "")
  x_left_p2 <- seq(0, 0.95, 0.05)
  x_right_p2 <- seq(0.05, 1, 0.05)
  for (i in 1:n_display) {
    bar_col <- if (x_right_p2[i] <= cut_show) rgb(0.85, 0.9, 1) else rgb(1, 0.92, 0.96)
    if (c_all[i] > 0) rect(x_left_p2[i], 0, x_right_p2[i], c_all[i], col = bar_col, border = "gray")
  }
  for (i in 1:n_display) {
    if (c_na[i] > 0) rect(x_left_p2[i], 0, x_right_p2[i], c_na[i], col = rgb(0.53, 0.81, 0.98, 0.85), border = "gray")
  }
  for (i in 1:n_display) {
    if (c_hs[i] > 0) rect(x_left_p2[i], 0, x_right_p2[i], c_hs[i], col = rgb(1, 0.5, 1, 0.85), border = "gray")
  }
  abline(v = 0.95, col = "darkgreen", lty = 2, lwd = 1.5)
  text(0.95, ymax_gt * 0.5, paste0("last bin: pair with\n P(H-S) >= ", round(cut_show, 5), "\nPredicted as H-S"), col = "darkgreen", adj = c(0.5, 0.5), cex = 1, xpd = TRUE)
  legend("top", legend = leg_text, fill = leg_fill, border = leg_border, lty = leg_lty, lwd = leg_lwd, col = leg_col, bty = "n", cex = 1.2, xpd = TRUE)
  dev.off()
  cat("Page 2 (P(HS) by GT class) saved to:", mixture_pdf_page2, "\n")

  # Page 2b: P(HS) by GT class (N/A vs HS only) + TP, FP, TN, FN at cut_show
  par(mfrow = c(1, 1), mar = c(5, 5, 4, 2), cex.lab = 1.35, cex.axis = 1.2, cex.main = 1.25)
  pred_cut <- P_HS_mix_all >= cut_show
  TP <- sum(idx_hs_gt & pred_cut)
  FN <- sum(idx_hs_gt & !pred_cut)
  TN <- sum(idx_na_gt & !pred_cut)
  FP <- sum(idx_na_gt & pred_cut)
  ymax_lab <- max(c(if (sum(idx_na_gt) > 0) max(hist(P_HS_mix_all[idx_na_gt], breaks = br, plot = FALSE)$counts) else 0, if (sum(idx_hs_gt) > 0) max(hist(P_HS_mix_all[idx_hs_gt], breaks = br, plot = FALSE)$counts) else 0, 1), na.rm = TRUE) * 1.2
  plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, ymax_lab), xlab = "P(H-S)", ylab = "Frequency", main = "P(HS) by GT class (4-D mixture, N/A vs HS only)")
  if (sum(idx_na_gt) > 0) hist(P_HS_mix_all[idx_na_gt], breaks = br, col = rgb(0.53, 0.81, 0.98, 0.7), border = "gray", add = TRUE)
  if (sum(idx_hs_gt) > 0) hist(P_HS_mix_all[idx_hs_gt], breaks = br, col = rgb(1, 0.5, 1, 0.7), border = "gray", add = TRUE)
  abline(v = cut_show, col = "red", lty = 2, lwd = 1.5)
  legend("topright", legend = c(paste("N/A (n=", sum(idx_na_gt), ")"), paste("HS (n=", sum(idx_hs_gt), ")")), fill = c(rgb(0.53, 0.81, 0.98, 0.8), rgb(1, 0.5, 1, 0.8)), bty = "n", cex = 1.2)
  txt <- paste("Cutoff =", cut_show, "\nTP =", TP, "  FP =", FP, "\nTN =", TN, "  FN =", FN)
  text(0.02, 0.98 * ymax_lab, txt, adj = c(0, 1), font = 2, cex = 1.2, family = "mono")
  cat("--- P(HS) 4-D mixture (N/A vs HS only, cutoff ", cut_show, ") ---\n", sep = "")
  cat("TP =", TP, " FP =", FP, " TN =", TN, " FN =", FN, "\n")

  # Page 3: density of P(HS) by GT class
  par(mfrow = c(1, 1), mar = c(5, 5, 4, 2), cex.lab = 1.35, cex.axis = 1.2, cex.main = 1.25)
  den_na <- if (sum(idx_na_gt) > 1) density(P_HS_mix_all[idx_na_gt], from = 0, to = 1, n = 128) else NULL
  den_hs <- if (sum(idx_hs_gt) > 1) density(P_HS_mix_all[idx_hs_gt], from = 0, to = 1, n = 128) else NULL
  den_other <- if (sum(idx_other_gt) > 1) density(P_HS_mix_all[idx_other_gt], from = 0, to = 1, n = 128) else NULL
  ymax_den <- max(c(if (!is.null(den_na)) den_na$y, if (!is.null(den_hs)) den_hs$y, if (!is.null(den_other)) den_other$y, 0), na.rm = TRUE) * 1.1
  plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, ymax_den), xlab = "P(H-S)", ylab = "Density", main = "P(HS) density by GT class (4-D mixture)")
  if (!is.null(den_na)) lines(den_na, col = "dodgerblue", lwd = 2.5)
  if (!is.null(den_hs)) lines(den_hs, col = "magenta", lwd = 2.5)
  if (!is.null(den_other)) lines(den_other, col = "gray40", lwd = 2.5)
  abline(v = cut_show, col = "red", lty = 2, lwd = 1.5)
  legend("topright", legend = c(paste("N/A (n=", sum(idx_na_gt), ")"), paste("HS (n=", sum(idx_hs_gt), ")"), paste("Other (n=", sum(idx_other_gt), ")")), col = c("dodgerblue", "magenta", "gray40"), lwd = 2.5, bty = "n", cex = 1.2)

  # Page 4: boxplot of P(HS) by GT class
  par(mfrow = c(1, 1), mar = c(5, 5, 4, 2), cex.lab = 1.35, cex.axis = 1.2, cex.main = 1.25)
  group_fac <- factor(group_mix, levels = c("N/A", "HS", "Other"))
  boxplot(P_HS_mix_all ~ group_fac, col = c("lightskyblue1", "#FF80FF", "lightgray"), border = c("dodgerblue", "magenta", "gray40"), main = "P(HS) by GT class (4-D mixture)", xlab = "GT class", ylab = "P(H-S)", ylim = c(0, 1))
  abline(h = cut_show, col = "red", lty = 2, lwd = 1)
  grid(nx = NA, ny = NULL, col = "lightgray", lty = "dotted")

  # Page 4a: Scatter P(HS) vs P(N/A) = 1 - P(HS), color by GT (points on diagonal)
  par(mfrow = c(1, 1), mar = c(5, 5, 4, 2), cex.lab = 1.35, cex.axis = 1.2, cex.main = 1.25)
  P_NA_mix <- 1 - P_HS_mix_all
  plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1), xlab = "P(H-S)", ylab = "P(N/A) = 1 - P(HS)", main = "P(HS) vs P(N/A) by GT class (4-D mixture)")
  abline(0, -1, col = "gray70", lty = 2, lwd = 1.5)
  points(P_HS_mix_all[idx_na_gt], P_NA_mix[idx_na_gt], col = "dodgerblue", pch = 19, cex = 0.6)
  points(P_HS_mix_all[idx_hs_gt], P_NA_mix[idx_hs_gt], col = "magenta", pch = 19, cex = 0.6)
  points(P_HS_mix_all[idx_other_gt], P_NA_mix[idx_other_gt], col = "gray40", pch = 19, cex = 0.6)
  abline(v = cut_show, col = "red", lty = 2, lwd = 1.5)
  legend("topright", legend = c(paste("N/A (n=", sum(idx_na_gt), ")"), paste("HS (n=", sum(idx_hs_gt), ")"), paste("Other (n=", sum(idx_other_gt), ")")), col = c("dodgerblue", "magenta", "gray40"), pch = 19, bty = "n", cex = 1.2)
  grid(col = "lightgray", lty = "dotted")

  # Page 4a2: Zoom histogram P(HS) in [0.9999, 1], GT only (N/A vs HS)
  par(mfrow = c(1, 1), mar = c(5, 5, 4, 2), cex.lab = 1.35, cex.axis = 1.2, cex.main = 1.25)
  zoom_lo <- 0.9999
  zoom_hi <- 1
  idx_zoom <- (P_HS_mix_all >= zoom_lo & P_HS_mix_all <= zoom_hi) & (idx_na_gt | idx_hs_gt)
  n_zoom <- sum(idx_zoom)
  br_zoom <- seq(zoom_lo, zoom_hi, length.out = 21)
  if (n_zoom > 0) {
    idx_na_zoom <- idx_na_gt & (P_HS_mix_all >= zoom_lo & P_HS_mix_all <= zoom_hi)
    idx_hs_zoom <- idx_hs_gt & (P_HS_mix_all >= zoom_lo & P_HS_mix_all <= zoom_hi)
    h_na_z <- if (sum(idx_na_zoom) > 0) hist(P_HS_mix_all[idx_na_zoom], breaks = br_zoom, plot = FALSE) else list(counts = rep(0, length(br_zoom) - 1))
    h_hs_z <- if (sum(idx_hs_zoom) > 0) hist(P_HS_mix_all[idx_hs_zoom], breaks = br_zoom, plot = FALSE) else list(counts = rep(0, length(br_zoom) - 1))
    ymax_z <- max(c(h_na_z$counts, h_hs_z$counts, 1), na.rm = TRUE) * 1.2
    plot(0, 0, type = "n", xlim = c(zoom_lo, zoom_hi), ylim = c(0, ymax_z), xlab = "P(H-S)", ylab = "Frequency", main = "P(HS) zoom [0.9999, 1] by GT (N/A vs HS only)")
    if (sum(idx_na_zoom) > 0) hist(P_HS_mix_all[idx_na_zoom], breaks = br_zoom, col = rgb(0.53, 0.81, 0.98, 0.7), border = "gray", add = TRUE)
    if (sum(idx_hs_zoom) > 0) hist(P_HS_mix_all[idx_hs_zoom], breaks = br_zoom, col = rgb(1, 0.5, 1, 0.7), border = "gray", add = TRUE)
    abline(v = cut_show, col = "red", lty = 2, lwd = 1.5)
    legend("topright", legend = c(paste("N/A (n=", sum(idx_na_zoom), ")"), paste("HS (n=", sum(idx_hs_zoom), ")"), paste("Total:", n_zoom)), fill = c(rgb(0.53, 0.81, 0.98, 0.8), rgb(1, 0.5, 1, 0.8), NA), bty = "n", cex = 1.1)
  } else {
    plot(0, 0, type = "n", xlim = c(zoom_lo, zoom_hi), ylim = c(0, 1), xlab = "P(H-S)", ylab = "Frequency", main = "P(HS) zoom [0.9999, 1] by GT (N/A vs HS only)")
    text(0.5 * (zoom_lo + zoom_hi), 0.5, "No N/A or HS in this range", cex = 1.2)
  }
  grid(col = "lightgray", lty = "dotted")

  # Page 4b & 4c: Illustration at cut_show and justification (when we have labeled data)
  if (any(idx_lab_mix) && exists("best_cutoff") && exists("best_i")) {
    # Page 4b: Illustration at chosen cutoff (N/A vs HS only)
    par(mfrow = c(1, 1), mar = c(5, 5, 4, 2), cex.lab = 1.35, cex.axis = 1.2, cex.main = 1.25)
    ymax_best <- max(c(if (sum(idx_na_gt) > 0) max(hist(P_HS_mix_all[idx_na_gt], breaks = br, plot = FALSE)$counts) else 0, if (sum(idx_hs_gt) > 0) max(hist(P_HS_mix_all[idx_hs_gt], breaks = br, plot = FALSE)$counts) else 0, 1), na.rm = TRUE) * 1.2
    plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, ymax_best), xlab = "P(H-S)", ylab = "Frequency", main = paste0("P(HS) by GT (N/A vs HS): cutoff ", round(cut_show, 6)))
    if (sum(idx_na_gt) > 0) hist(P_HS_mix_all[idx_na_gt], breaks = br, col = rgb(0.53, 0.81, 0.98, 0.7), border = "gray", add = TRUE)
    if (sum(idx_hs_gt) > 0) hist(P_HS_mix_all[idx_hs_gt], breaks = br, col = rgb(1, 0.5, 1, 0.7), border = "gray", add = TRUE)
    abline(v = cut_show, col = "red", lty = 1, lwd = 2.5)
    legend("topright", legend = c(paste("N/A (n=", sum(idx_na_gt), ")"), paste("HS (n=", sum(idx_hs_gt), ")")), fill = c(rgb(0.53, 0.81, 0.98, 0.8), rgb(1, 0.5, 1, 0.8)), bty = "n", cex = 1.2)
    txt_best <- paste("Cutoff =", round(cut_show, 6), "\nTP =", TP, "  TN =", TN, "  FP =", FP, "  FN =", FN)
    text(0.02, 0.98 * ymax_best, txt_best, adj = c(0, 1), font = 2, cex = 1.1, family = "mono")
    # Page 4c: Justification — Accuracy (and Sens, Spec) vs cutoff
    par(mfrow = c(1, 1), mar = c(5, 5, 4, 2), cex.lab = 1.35, cex.axis = 1.2, cex.main = 1.25)
    plot(th_acc, acc_vec, type = "l", lwd = 2, col = "darkblue", xlab = "Cutoff (P(HS) >= cutoff -> predict HS)", ylab = "Accuracy", main = "Justification: accuracy is maximized at chosen cutoff", xlim = c(0, 1), ylim = c(0, 1))
    abline(v = cut_show, col = "red", lty = 2, lwd = 2)
    acc_at_cut <- (TP + TN) / (TP + TN + FP + FN)
    points(cut_show, acc_at_cut, pch = 19, col = "red", cex = 1.5)
    lines(th_acc, sens_vec, col = "magenta", lwd = 1.5, lty = 2)
    lines(th_acc, spec_vec, col = "dodgerblue", lwd = 1.5, lty = 2)
    legend("right", legend = c("Accuracy", "Sensitivity", "Specificity", paste("Cutoff =", cut_show)), col = c("darkblue", "magenta", "dodgerblue", "red"), lwd = c(2, 1.5, 1.5, 2), lty = c(1, 2, 2, 2), bty = "n", cex = 1.1)
    grid(col = "lightgray", lty = "dotted")
  }

  # Page 5: calibration, P(HS) by group (labeled), ROC (1x3)
  par(mfrow = c(1, 3), mar = c(4, 4, 3, 1), cex.lab = 1.35, cex.axis = 1.2, cex.main = 1.25)
  if (any(idx_lab_mix)) {
    P_HS_mix <- P_HS_mix_all[idx_lab_mix]
    y_lab_mix <- y_mix[idx_lab_mix]
    # 1. Calibration (reliability) for 4-D mixture
    n_bins_m <- 10
    br_m <- seq(0, 1, length.out = n_bins_m + 1)
    bins_m <- cut(P_HS_mix, breaks = br_m, include.lowest = TRUE)
    pred_bin_m <- sapply(levels(bins_m), function(b) { i <- which(bins_m == b); if (length(i)) mean(P_HS_mix[i]) else NA })
    obs_bin_m  <- sapply(levels(bins_m), function(b) { i <- which(bins_m == b); if (length(i)) mean(y_lab_mix[i]) else NA })
    plot(c(0, 1), c(0, 1), type = "n", xlab = "Mean predicted P(HS)", ylab = "Observed fraction HS", main = "Calibration (4-D mixture)")
    abline(0, 1, col = "gray", lty = 2)
    points(pred_bin_m, obs_bin_m, pch = 19, col = "darkgreen", cex = 1.2)
    grid(col = "lightgray", lty = "dotted")

    # 2. Histogram of P(HS) by group (4-D mixture)
    idx_hs_mix <- y_lab_mix == 1
    idx_na_mix <- y_lab_mix == 0
    h_na_mix <- hist(P_HS_mix[idx_na_mix], breaks = seq(0, 1, 0.05), plot = FALSE)
    h_hs_mix <- hist(P_HS_mix[idx_hs_mix], breaks = seq(0, 1, 0.05), plot = FALSE)
    ymax_mix <- max(c(h_na_mix$counts, h_hs_mix$counts, 1), na.rm = TRUE) * 1.2
    hist(P_HS_mix[idx_na_mix], breaks = seq(0, 1, 0.05), col = "lightskyblue1", border = "gray", main = "P(HS) by group (4-D mixture)", xlab = "P(H-S)", xlim = c(0, 1), ylim = c(0, ymax_mix))
    hist(P_HS_mix[idx_hs_mix], breaks = seq(0, 1, 0.05), col = "#FF80FF", border = "gray", add = TRUE)
    abline(v = cut_show, col = "red", lty = 2, lwd = 1.5)
    legend("topright", legend = c("N/A", "HS"), fill = c("lightskyblue1", "#FF80FF"), bty = "n", cex = 1.3)

    # 3. ROC curve (4-D mixture)
    th_m <- sort(unique(c(0, P_HS_mix, 1)))
    sens_m <- spec_m <- numeric(length(th_m))
    for (i in seq_along(th_m)) {
      pred_pos_m <- P_HS_mix >= th_m[i]
      TP <- sum(y_lab_mix == 1 & pred_pos_m, na.rm = TRUE)
      FN <- sum(y_lab_mix == 1 & !pred_pos_m, na.rm = TRUE)
      TN <- sum(y_lab_mix == 0 & !pred_pos_m, na.rm = TRUE)
      FP <- sum(y_lab_mix == 0 & pred_pos_m, na.rm = TRUE)
      tot_pos <- TP + FN; tot_neg <- TN + FP
      sens_m[i] <- if (is.finite(tot_pos) && tot_pos > 0) TP / tot_pos else 0
      spec_m[i] <- if (is.finite(tot_neg) && tot_neg > 0) TN / tot_neg else 0
    }
    plot(1 - spec_m, sens_m, type = "l", lwd = 2, col = "darkgreen", xlab = "1 - Specificity", ylab = "Sensitivity", main = "ROC (4-D mixture)")
    abline(0, 1, col = "gray", lty = 2)
    grid(col = "lightgray", lty = "dotted")
  } else {
    plot.new()
    text(0.5, 0.7, "Calibration (4-D mixture)", cex = 1.2)
    text(0.5, 0.5, "No labeled rows (Case/Control) with 4-D data.", cex = 1)
    hist(P_HS_mix_all, breaks = seq(0, 1, 0.05), col = "lightgray", border = "gray", main = "P(HS) (4-D mixture, all)", xlab = "P(H-S)", xlim = c(0, 1))
    plot.new()
    text(0.5, 0.5, "ROC requires labeled data.", cex = 1.2)
  }

  dev.off()
  cat("Multivariate mixture plots saved to:", mixture_pdf, "\n")

  # --- P(H-S) vs age difference (separate PDF) ---
  # Custom y-axis so high P(HS) (0.99–1) is spread out. Exclude age=0. Regression on "Others" only.
  trans_phs_axis <- function(p) {
    p <- pmax(0, pmin(1, p))
    coord <- numeric(length(p))
    for (i in seq_along(p)) {
      if (p[i] <= 0.9) { coord[i] <- 0.25 * (p[i] / 0.9)
      } else if (p[i] <= 0.99) { coord[i] <- 0.25 + 0.25 * (p[i] - 0.9) / 0.09
      } else if (p[i] <= 0.999) { coord[i] <- 0.5 + 0.25 * (p[i] - 0.99) / 0.009
      } else if (p[i] <= 0.9999) { coord[i] <- 0.75 + 0.125 * (p[i] - 0.999) / 0.0009
      } else if (p[i] <= 0.99999) { coord[i] <- 0.875 + 0.0625 * (p[i] - 0.9999) / 0.00009
      } else { coord[i] <- 0.9375 + 0.0625 * (p[i] - 0.99999) / 0.00001 }
    }
    coord
  }
  if (all(c("age1", "age2") %in% names(data_for_supervised))) {
    age1_num <- as.numeric(as.character(data_for_supervised$age1[res_mix$ok_all]))
    age2_num <- as.numeric(as.character(data_for_supervised$age2[res_mix$ok_all]))
    age_diff <- abs(age1_num - age2_num)
    idx_age_ok <- is.finite(age_diff) & !is.na(age_diff) & (age1_num != 0) & (age2_num != 0)
    pred_hs <- P_HS_mix_all >= cut_show
    col_na <- "dodgerblue"
    col_hs <- "magenta"
    col_other_na <- rgb(0.85, 0.9, 1)
    col_other_hs <- rgb(1, 0.75, 0.85)
    if (sum(idx_age_ok) > 0) {
      pdf_age <- "/pl/active/KellerLab/Emmanuel/NAvHS/GMM_plot_PHS_vs_age_diff.pdf"
      pdf(pdf_age, width = 10, height = 8)
      par(mar = c(5, 5, 4, 2), cex.lab = 1.35, cex.axis = 1.2, cex.main = 1.25)
      xr <- range(age_diff[idx_age_ok], na.rm = TRUE)
      plot(NA, NA, xlim = xr, ylim = c(0, 1), xlab = "Age difference", ylab = "P(H-S)", main = "", yaxt = "n")
      axis(2, at = c(0, 0.25, 0.5, 0.75, 0.875, 0.9375, 1), labels = c("0", "0.9", "0.99", "0.9999", "0.99999", "0.99999", "1"), las = 1)
      P_HS_trans <- trans_phs_axis(P_HS_mix_all)
      idx_na_age <- idx_na_gt & idx_age_ok
      idx_hs_age <- idx_hs_gt & idx_age_ok
      idx_other_na_age <- idx_other_gt & !pred_hs & idx_age_ok
      idx_other_hs_age <- idx_other_gt & pred_hs & idx_age_ok
      abline(h = trans_phs_axis(cut_show), col = "green", lty = 2, lwd = 1.5)
      if (sum(idx_other_na_age) > 0) points(age_diff[idx_other_na_age], P_HS_trans[idx_other_na_age], col = col_other_na, pch = 19, cex = 0.7)
      if (sum(idx_other_hs_age) > 0) points(age_diff[idx_other_hs_age], P_HS_trans[idx_other_hs_age], col = col_other_hs, pch = 19, cex = 0.7)
      if (sum(idx_na_age) > 0) points(age_diff[idx_na_age], P_HS_trans[idx_na_age], col = col_na, pch = 19, cex = 0.7)
      if (sum(idx_hs_age) > 0) points(age_diff[idx_hs_age], P_HS_trans[idx_hs_age], col = col_hs, pch = 19, cex = 0.7)
      idx_other_age <- idx_other_gt & idx_age_ok
      idx_gt_age <- (idx_na_gt | idx_hs_gt) & idx_age_ok
      has_fit <- sum(idx_other_age) >= 2
      if (has_fit) {
        fit_other <- lm(P_HS_mix_all[idx_other_age] ~ age_diff[idx_other_age])
        x_line <- seq(xr[1], xr[2], length.out = 100)
        pred_line <- pmax(0, pmin(1, predict(fit_other, newdata = data.frame(age_diff = x_line))))
        lines(x_line, trans_phs_axis(pred_line), col = "darkgreen", lwd = 2, lty = 1)
      }
      cor_other <- if (sum(idx_other_age) >= 2) cor(age_diff[idx_other_age], P_HS_mix_all[idx_other_age], use = "pairwise.complete.obs") else NA_real_
      cor_all   <- if (sum(idx_age_ok) >= 2) cor(age_diff[idx_age_ok], P_HS_mix_all[idx_age_ok], use = "pairwise.complete.obs") else NA_real_
      cor_gt    <- if (sum(idx_gt_age) >= 2) cor(age_diff[idx_gt_age], P_HS_mix_all[idx_gt_age], use = "pairwise.complete.obs") else NA_real_
      fit_all <- if (sum(idx_age_ok) >= 2) lm(P_HS_mix_all[idx_age_ok] ~ age_diff[idx_age_ok]) else NULL
      fit_gt  <- if (sum(idx_gt_age) >= 2) lm(P_HS_mix_all[idx_gt_age] ~ age_diff[idx_gt_age]) else NULL
      slope_other <- if (has_fit) coef(fit_other)[2] else NA_real_
      slope_all   <- if (!is.null(fit_all)) coef(fit_all)[2] else NA_real_
      slope_gt    <- if (!is.null(fit_gt)) coef(fit_gt)[2] else NA_real_
      r2_other <- if (has_fit) summary(fit_other)$r.squared else NA_real_
      r2_all   <- if (!is.null(fit_all)) summary(fit_all)$r.squared else NA_real_
      r2_gt    <- if (!is.null(fit_gt)) summary(fit_gt)$r.squared else NA_real_
      cat("--- P(H-S) vs age difference: correlation and regression ---\n")
      cat("Others (n=", sum(idx_other_age), "): r =", round(cor_other, 4), " slope =", round(slope_other, 6), " R2 =", round(r2_other, 4), "\n")
      cat("All (n=", sum(idx_age_ok), "): r =", round(cor_all, 4), " slope =", round(slope_all, 6), " R2 =", round(r2_all, 4), "\n")
      cat("GT only (n=", sum(idx_gt_age), "): r =", round(cor_gt, 4), " slope =", round(slope_gt, 6), " R2 =", round(r2_gt, 4), "\n")
      leg_lab <- c(paste0("Other Predicted N/N-A (n=", sum(idx_other_na_age), ")"), paste0("Other Predicted H-S (n=", sum(idx_other_hs_age), ")"), paste0("N/N-A (n=", sum(idx_na_age), ")"), paste0("H-S (n=", sum(idx_hs_age), ")"))
      leg_col <- c(col_other_na, col_other_hs, col_na, col_hs)
      leg_pch <- c(19, 19, 19, 19)
      leg_lty <- c(NA, NA, NA, NA)
      leg_lwd <- c(NA, NA, NA, NA)
      if (has_fit) {
        leg_lab <- c(leg_lab, "Regression (Others only)")
        leg_col <- c(leg_col, "darkgreen")
        leg_pch <- c(leg_pch, NA)
        leg_lty <- c(leg_lty, 1)
        leg_lwd <- c(leg_lwd, 2)
      }
      leg_lab <- c(leg_lab, "", paste0("Others: r = ", round(cor_other, 4), " slope = ", round(slope_other, 5), " R2 = ", round(r2_other, 4)), paste0("All: r = ", round(cor_all, 4), " slope = ", round(slope_all, 5), " R2 = ", round(r2_all, 4)), paste0("GT: r = ", round(cor_gt, 4), " slope = ", round(slope_gt, 5), " R2 = ", round(r2_gt, 4)))
      leg_col <- c(leg_col, "white", "black", "black", "black")
      leg_pch <- c(leg_pch, NA, NA, NA, NA)
      leg_lty <- c(leg_lty, NA, NA, NA, NA)
      leg_lwd <- c(leg_lwd, NA, NA, NA, NA)
      par(fig = c(0.02, 0.34, 0.18, 0.85), new = TRUE, mar = c(0.5, 0.5, 0.5, 0.5))
      plot.new()
      plot.window(xlim = c(0, 1), ylim = c(0, 1))
      legend(0.02, 0.98, legend = leg_lab, col = leg_col, pch = leg_pch, lty = leg_lty, lwd = leg_lwd, bty = "o", bg = "white", border = "gray40", xjust = 0, yjust = 1, cex = 0.82, y.intersp = 1.02)
      par(fig = c(0, 1, 0, 1))
      dev.off()
      cat("P(H-S) vs age difference (non-missing age only) saved to:", pdf_age, "\n")
    }
  }
} else if (nrow(data_for_supervised) > 0 && all(c("sum_pihat_max_5", "sum_pihat_remaining_max_5", "sum_pihat_remaining_min_5", "sum_pihat_opposite_5") %in% names(data_for_supervised))) {
  cat("Multivariate mixture PDF not created: fit failed or was skipped (see messages above).\n")
}

# ============================================================================
# Scatter: max(pihat against all) vs sum_pihat_remaining_max_5 (case types 2, 7, 11)
# ============================================================================
if ("type" %in% names(data) && 
    "max_pihatagainstall2_below_033_pihat1" %in% names(data) &&
    "max_pihatagainstall1_below_033_pihat2" %in% names(data) &&
    "sum_pihat_remaining_max_5" %in% names(data)) {
    
    # Use all data, but we'll filter to rows with non-NA maximums values
    data_for_plots <- data
    
    # Calculate max value for each row
    max_pihat_val <- pmax(
        as.numeric(data_for_plots$max_pihatagainstall2_below_033_pihat1),
        as.numeric(data_for_plots$max_pihatagainstall1_below_033_pihat2),
        na.rm = TRUE
    )
    # Set to NA if both are NA
    max_pihat_val[is.na(data_for_plots$max_pihatagainstall2_below_033_pihat1) & 
                 is.na(data_for_plots$max_pihatagainstall1_below_033_pihat2)] <- NA
    
    # Get sum_pihat_remaining_max_5 values
    sum_pihat_vals <- as.numeric(data_for_plots$sum_pihat_remaining_max_5)
    
    # Filter for different types (only rows with non-NA maximums values will be included)
    idx_type2 <- which(data_for_plots$type == 2 & !is.na(max_pihat_val) & !is.na(sum_pihat_vals))
    idx_type7 <- which(data_for_plots$type == 7 & !is.na(max_pihat_val) & !is.na(sum_pihat_vals))
    idx_type11 <- which(data_for_plots$type == 11 & !is.na(max_pihat_val) & !is.na(sum_pihat_vals))
    
    if (length(idx_type2) > 0 || length(idx_type7) > 0 || length(idx_type11) > 0) {
        par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1, cex.lab = 1.35, cex.axis = 1.2, cex.main = 1.25)
        
        # Determine plot limits
        all_x <- c(sum_pihat_vals[idx_type2], sum_pihat_vals[idx_type7], sum_pihat_vals[idx_type11])
        all_y <- c(max_pihat_val[idx_type2], max_pihat_val[idx_type7], max_pihat_val[idx_type11])
        xlim <- range(all_x, na.rm = TRUE)
        ylim <- range(all_y, na.rm = TRUE)
        
        # Create empty plot
        plot(0, 0, type = "n", 
             xlim = xlim, ylim = ylim,
             xlab = "sum_pihat_remaining_max_5",
             ylab = "max(max_pihatagainstall2_below_033_pihat1, max_pihatagainstall1_below_033_pihat2)",
             main = "Scatter plot: max pihat vs sum_pihat_remaining_max_5\n(type=2, 7, 11)")
        
        # Add grid with 0.02 spacing
        x_grid <- seq(floor(xlim[1] / 0.02) * 0.02, ceiling(xlim[2] / 0.02) * 0.02, by = 0.02)
        y_grid <- seq(floor(ylim[1] / 0.02) * 0.02, ceiling(ylim[2] / 0.02) * 0.02, by = 0.02)
        abline(v = x_grid, col = "lightgray", lty = "dotted", lwd = 0.5)
        abline(h = y_grid, col = "lightgray", lty = "dotted", lwd = 0.5)
        
        # Plot points for each type with different colors
        if (length(idx_type2) > 0) {
            points(sum_pihat_vals[idx_type2], max_pihat_val[idx_type2], 
                   col = "purple", pch = 19, cex = 0.8)
        }
        if (length(idx_type7) > 0) {
            points(sum_pihat_vals[idx_type7], max_pihat_val[idx_type7], 
                   col = "blue", pch = 19, cex = 0.8)
        }
        if (length(idx_type11) > 0) {
            points(sum_pihat_vals[idx_type11], max_pihat_val[idx_type11], 
                   col = "red", pch = 19, cex = 0.8)
        }
        
        # Add legend
        legend_items <- c()
        legend_cols <- c()
        if (length(idx_type2) > 0) {
            legend_items <- c(legend_items, paste("type=2 (n=", length(idx_type2), ")", sep=""))
            legend_cols <- c(legend_cols, "purple")
        }
        if (length(idx_type7) > 0) {
            legend_items <- c(legend_items, paste("type=7 (n=", length(idx_type7), ")", sep=""))
            legend_cols <- c(legend_cols, "blue")
        }
        if (length(idx_type11) > 0) {
            legend_items <- c(legend_items, paste("type=11 (n=", length(idx_type11), ")", sep=""))
            legend_cols <- c(legend_cols, "red")
        }
        if (length(legend_items) > 0) {
            legend("topright", 
                   legend = legend_items,
                   col = legend_cols, pch = 19, cex = 1.3)
        }
        
        cat("Scatter plot created: ", length(idx_type2), " type=2, ", 
            length(idx_type7), " type=7, ", length(idx_type11), " type=11 points\n")
    } else {
        cat("No valid data for scatter plot (no case types with non-NA maximums values)\n")
    }
} else {
    cat("Warning: Required columns not found in data for scatter plot\n")
}

# ============================================================================
# Pairwise scatter of 4 max_pihatagainstall* variables, colored by TP/FN (GMM-based)
# ============================================================================
# For types 2, 7, 11 separately: 6 pairs of variables, green=TP, red=FN (from univariate GMM on var1).
var_names <- c("max_pihatagainstall1_below_033_pihat1",
               "max_pihatagainstall1_below_033_pihat2",
               "max_pihatagainstall2_below_033_pihat1",
               "max_pihatagainstall2_below_033_pihat2")

# Check if all variables exist in data
all_vars_exist <- all(var_names %in% names(data))
has_sum_var <- "sum_pihat_remaining_max_5" %in% names(data)
has_type <- "type" %in% names(data)

if (all_vars_exist && has_sum_var && has_type) {
    # Use all data
    data_for_plots <- data
    # Extract variables as numeric
    var_data <- list()
    for (v in var_names) {
        var_data[[v]] <- as.numeric(data_for_plots[[v]])
    }
    sum_var <- as.numeric(data_for_plots$sum_pihat_remaining_max_5)
    type_var <- data_for_plots$type
    
    # Calculate TP/FN for coloring
    # Get sum_pihat_remaining_max_5 values for data_for_plots
    sum_pihat_vals_plot <- as.numeric(data_for_plots$sum_pihat_remaining_max_5)
    
    # Calculate GMM predictions for these observations
    # Use the GMM from the first variable (sum_pihat_remaining_max_5)
    tp_fn_status <- rep(NA, nrow(data_for_plots))
    if (length(variable_list) >= 1 && !is.null(variable_list[[1]]$data)) {
        var1_info <- variable_list[[1]]
        var1_data_all <- var1_info$data
        gmm_results_var1 <- calculate_gmm_for_variable(var1_data_all, var1_info$name)
        
        if (!is.null(gmm_results_var1$gmm_all)) {
            # Predict clusters for sum_pihat_remaining_max_5 values (only for non-NA values)
            valid_for_pred <- !is.na(sum_pihat_vals_plot) & is.finite(sum_pihat_vals_plot)
            predicted_cluster <- rep(NA, length(sum_pihat_vals_plot))
            
            if (sum(valid_for_pred) > 0) {
                prob_cluster1 <- gmm_results_var1$gmm_all$lambda[1] * dnorm(sum_pihat_vals_plot[valid_for_pred], 
                                                                            mean = gmm_results_var1$gmm_all$mu[1], 
                                                                            sd = gmm_results_var1$gmm_all$sigma[1])
                prob_cluster2 <- gmm_results_var1$gmm_all$lambda[2] * dnorm(sum_pihat_vals_plot[valid_for_pred], 
                                                                            mean = gmm_results_var1$gmm_all$mu[2], 
                                                                            sd = gmm_results_var1$gmm_all$sigma[2])
                total_prob <- prob_cluster1 + prob_cluster2
                prob_cluster1 <- prob_cluster1 / total_prob
                prob_cluster2 <- prob_cluster2 / total_prob
                
                predicted_cluster[valid_for_pred] <- ifelse(prob_cluster1 > prob_cluster2, 1, 2)
            }
            
            # Determine which cluster corresponds to "Case (HS)"
            # Use the same logic as in create_gmm_plot - find which cluster has more Case observations
            # For types 2, 7, 11, true_label = 1 (Case)
            true_labels_plot <- ifelse(type_var %in% c(2, 7, 11), 1, 
                                      ifelse(type_var == 100, 0, NA))
            
            # Find which cluster best corresponds to "Case (HS)" (only using valid predictions)
            valid_pred_idx <- !is.na(predicted_cluster) & !is.na(true_labels_plot)
            if (sum(valid_pred_idx) > 0) {
                case_in_cluster1 <- sum(true_labels_plot[valid_pred_idx] == 1 & predicted_cluster[valid_pred_idx] == 1)
                case_in_cluster2 <- sum(true_labels_plot[valid_pred_idx] == 1 & predicted_cluster[valid_pred_idx] == 2)
                
                # Align clusters: the cluster with the most "Case" becomes 1
                if (case_in_cluster2 > case_in_cluster1) {
                    pred_labels_aligned <- ifelse(predicted_cluster == 2, 1, 
                                                 ifelse(predicted_cluster == 1, 0, NA))
                } else {
                    pred_labels_aligned <- ifelse(predicted_cluster == 1, 1, 
                                                 ifelse(predicted_cluster == 2, 0, NA))
                }
                
                # Calculate TP/FN for case types (true_label = 1)
                # TP: predicted as Case (1) and actually Case (1)
                # FN: predicted as Control (0) but actually Case (1)
                for (k in 1:length(true_labels_plot)) {
                    if (!is.na(true_labels_plot[k]) && true_labels_plot[k] == 1) {
                        if (!is.na(pred_labels_aligned[k])) {
                            if (pred_labels_aligned[k] == 1) {
                                tp_fn_status[k] <- "TP"
                            } else {
                                tp_fn_status[k] <- "FN"
                            }
                        }
                    }
                }
            }
        }
    }
    
    # Generate all pairs (6 pairs for 4 variables)
    n_vars <- length(var_names)
    
    # Function to create scatter plots for a given type
    create_scatter_plots <- function(type_value, type_label) {
        # Filter by type
        type_filter <- type_var == type_value
        
        if (sum(type_filter) == 0) {
            cat("No observations with type=", type_value, " for scatter plots\n", sep="")
            return(0)
        }
        
        pair_count <- 0
        
        for (i in 1:(n_vars - 1)) {
            for (j in (i + 1):n_vars) {
                var1_name <- var_names[i]
                var2_name <- var_names[j]
                var1_data <- var_data[[var1_name]]
                var2_data <- var_data[[var2_name]]
                
                # Filter by type and remove NA values
                valid_idx <- type_filter & !is.na(var1_data) & !is.na(var2_data)
                
                if (sum(valid_idx) > 0) {
                    pair_count <- pair_count + 1
                    par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1, cex.lab = 1.35, cex.axis = 1.2, cex.main = 1.25)
                    
                    # Get colors based on TP/FN status
                    # TP = green, FN = red
                    colors <- rep("gray", sum(valid_idx))
                    colors[tp_fn_status[valid_idx] == "TP"] <- "green"
                    colors[tp_fn_status[valid_idx] == "FN"] <- "red"
                    
                    # Determine plot limits
                    xlim_plot <- range(var1_data[valid_idx], na.rm = TRUE)
                    ylim_plot <- range(var2_data[valid_idx], na.rm = TRUE)
                    
                    # Create empty plot to establish coordinate system
                    plot(0, 0, type = "n",
                         xlim = xlim_plot,
                         ylim = ylim_plot,
                         xlab = var1_name,
                         ylab = var2_name,
                         main = paste("Scatter plot (type=", type_label, "): ", var1_name, " vs ", var2_name, 
                                     "\n(green: TP, red: FN)", sep=""))
                    
                    # Add grid with 0.02 spacing
                    x_grid <- seq(floor(xlim_plot[1] / 0.02) * 0.02, ceiling(xlim_plot[2] / 0.02) * 0.02, by = 0.02)
                    y_grid <- seq(floor(ylim_plot[1] / 0.02) * 0.02, ceiling(ylim_plot[2] / 0.02) * 0.02, by = 0.02)
                    abline(v = x_grid, col = "lightgray", lty = "dotted", lwd = 0.5)
                    abline(h = y_grid, col = "lightgray", lty = "dotted", lwd = 0.5)
                    
                    # Add points on top of grid
                    points(var1_data[valid_idx], var2_data[valid_idx],
                           col = colors,
                           pch = 19,
                           cex = 0.8)
                    
                    # Add legend
                    n_tp <- sum(tp_fn_status[valid_idx] == "TP", na.rm = TRUE)
                    n_fn <- sum(tp_fn_status[valid_idx] == "FN", na.rm = TRUE)
                    n_other <- sum(is.na(tp_fn_status[valid_idx]) | (!tp_fn_status[valid_idx] %in% c("TP", "FN")))
                    legend_items <- c()
                    legend_cols <- c()
                    if (n_tp > 0) {
                        legend_items <- c(legend_items, paste("TP (n=", n_tp, ")", sep=""))
                        legend_cols <- c(legend_cols, "green")
                    }
                    if (n_fn > 0) {
                        legend_items <- c(legend_items, paste("FN (n=", n_fn, ")", sep=""))
                        legend_cols <- c(legend_cols, "red")
                    }
                    if (n_other > 0) {
                        legend_items <- c(legend_items, paste("Other (n=", n_other, ")", sep=""))
                        legend_cols <- c(legend_cols, "gray")
                    }
                    if (length(legend_items) > 0) {
                        legend("topright",
                               legend = legend_items,
                               col = legend_cols,
                               pch = 19,
                               cex = 1.3)
                    }
                    
                    cat("Scatter plot (type=", type_label, ") ", pair_count, ": ", var1_name, " vs ", var2_name, 
                        " (", sum(valid_idx), " points, TP=", n_tp, ", FN=", n_fn, ")\n", sep="")
                }
            }
        }
        
        return(pair_count)
    }
    
    # Create scatter plots for type=2
    cat("\nCreating scatter plots for type=2...\n")
    count_type2 <- create_scatter_plots(2, "2")
    
    # Create scatter plots for type=7
    cat("\nCreating scatter plots for type=7...\n")
    count_type7 <- create_scatter_plots(7, "7")
    
    # Create scatter plots for type=11
    cat("\nCreating scatter plots for type=11...\n")
    count_type11 <- create_scatter_plots(11, "11")
    
    cat("\nCreated ", count_type2, " scatter plots for type=2, ", 
        count_type7, " for type=7, and ", 
        count_type11, " scatter plots for type=11\n")
    } else {
        cat("Warning: Required columns not found in data for pair scatter plots\n")
        if (!all_vars_exist) {
            missing_vars <- var_names[!var_names %in% names(data_for_plots)]
            cat("Missing variables:", paste(missing_vars, collapse = ", "), "\n")
        }
        if (!has_sum_var) {
            cat("Missing variable: sum_pihat_remaining_max_5\n")
        }
        if (!has_type) {
            cat("Missing variable: type\n")
        }
    }

# ============================================================================
# Pihat vs Pihat2 scatter by type (2, 7, 11: TP/FN; type 100: TN/FP)
# ============================================================================
cat("\n=== Creating Pihat vs Pihat2 scatter plots ===\n")
if ("Pihat" %in% names(data) && "Pihat2" %in% names(data) && "type" %in% names(data)) {
    # Get Pihat and Pihat2 values
    pihat_vals <- as.numeric(data$Pihat)
    pihat2_vals <- as.numeric(data$Pihat2)
    type_vals_plot <- data$type
    
    # Calculate TP/FN for all data (same as in pairwise scatter plots)
    # Use the GMM from the first variable (sum_pihat_remaining_max_5)
    tp_fn_status_plot <- rep(NA, nrow(data))
    if (length(variable_list) >= 1 && !is.null(variable_list[[1]]$data)) {
        var1_info <- variable_list[[1]]
        var1_data_all <- var1_info$data
        gmm_results_var1 <- calculate_gmm_for_variable(var1_data_all, var1_info$name)
        
        if (!is.null(gmm_results_var1$gmm_all)) {
            # Get sum_pihat_remaining_max_5 for predictions
            sum_pihat_vals_plot_all <- as.numeric(data$sum_pihat_remaining_max_5)
            
            # Predict clusters for sum_pihat_remaining_max_5 values (only for non-NA values)
            valid_for_pred_all <- !is.na(sum_pihat_vals_plot_all) & is.finite(sum_pihat_vals_plot_all)
            predicted_cluster_all <- rep(NA, length(sum_pihat_vals_plot_all))
            
            if (sum(valid_for_pred_all) > 0) {
                prob_cluster1 <- gmm_results_var1$gmm_all$lambda[1] * dnorm(sum_pihat_vals_plot_all[valid_for_pred_all], 
                                                                            mean = gmm_results_var1$gmm_all$mu[1], 
                                                                            sd = gmm_results_var1$gmm_all$sigma[1])
                prob_cluster2 <- gmm_results_var1$gmm_all$lambda[2] * dnorm(sum_pihat_vals_plot_all[valid_for_pred_all], 
                                                                            mean = gmm_results_var1$gmm_all$mu[2], 
                                                                            sd = gmm_results_var1$gmm_all$sigma[2])
                total_prob <- prob_cluster1 + prob_cluster2
                prob_cluster1 <- prob_cluster1 / total_prob
                prob_cluster2 <- prob_cluster2 / total_prob
                
                predicted_cluster_all[valid_for_pred_all] <- ifelse(prob_cluster1 > prob_cluster2, 1, 2)
            }
            
            # Determine which cluster corresponds to "Case (HS)"
            true_labels_plot_all <- ifelse(type_vals_plot %in% c(2, 7, 11), 1, 
                                          ifelse(type_vals_plot == 100, 0, NA))
            
            # Find which cluster best corresponds to "Case (HS)" (only using valid predictions)
            valid_pred_idx_all <- !is.na(predicted_cluster_all) & !is.na(true_labels_plot_all)
            if (sum(valid_pred_idx_all) > 0) {
                case_in_cluster1 <- sum(true_labels_plot_all[valid_pred_idx_all] == 1 & predicted_cluster_all[valid_pred_idx_all] == 1)
                case_in_cluster2 <- sum(true_labels_plot_all[valid_pred_idx_all] == 1 & predicted_cluster_all[valid_pred_idx_all] == 2)
                
                # Align clusters: the cluster with the most "Case" becomes 1
                if (case_in_cluster2 > case_in_cluster1) {
                    pred_labels_aligned_all <- ifelse(predicted_cluster_all == 2, 1, 
                                                     ifelse(predicted_cluster_all == 1, 0, NA))
                } else {
                    pred_labels_aligned_all <- ifelse(predicted_cluster_all == 1, 1, 
                                                     ifelse(predicted_cluster_all == 2, 0, NA))
                }
                
                # Calculate TP/FN for case types (true_label = 1)
                for (k in 1:length(true_labels_plot_all)) {
                    if (!is.na(true_labels_plot_all[k]) && true_labels_plot_all[k] == 1) {
                        if (!is.na(pred_labels_aligned_all[k])) {
                            if (pred_labels_aligned_all[k] == 1) {
                                tp_fn_status_plot[k] <- "TP"
                            } else {
                                tp_fn_status_plot[k] <- "FN"
                            }
                        }
                    }
                }
                
                # Calculate TN/FP for type 100 (Control, true_label = 0)
                # TN = predicted as Control (0) and actually Control (0)
                # FP = predicted as Case (1) but actually Control (0)
                tn_fp_status_plot <- rep(NA, nrow(data))
                for (k in 1:length(true_labels_plot_all)) {
                    if (!is.na(true_labels_plot_all[k]) && true_labels_plot_all[k] == 0) {
                        if (!is.na(pred_labels_aligned_all[k])) {
                            if (pred_labels_aligned_all[k] == 0) {
                                tn_fp_status_plot[k] <- "TN"
                            } else {
                                tn_fp_status_plot[k] <- "FP"
                            }
                        }
                    }
                }
            }
        }
    }
    
    # Create scatter plots for each case type
    case_types <- c(2, 7, 11)
    for (type_val in case_types) {
        # Filter for this type and non-NA Pihat/Pihat2 values
        type_filter_plot <- type_vals_plot == type_val & !is.na(pihat_vals) & !is.na(pihat2_vals)
        
        if (sum(type_filter_plot) > 0) {
            par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1, cex.lab = 1.35, cex.axis = 1.2, cex.main = 1.25)
            
            # Get colors based on TP/FN status
            colors_plot <- rep("gray", sum(type_filter_plot))
            colors_plot[tp_fn_status_plot[type_filter_plot] == "TP"] <- "green"
            colors_plot[tp_fn_status_plot[type_filter_plot] == "FN"] <- "red"
            
            # Determine plot limits
            xlim_plot <- range(pihat_vals[type_filter_plot], na.rm = TRUE)
            ylim_plot <- range(pihat2_vals[type_filter_plot], na.rm = TRUE)
            
            # Create empty plot
            plot(0, 0, type = "n",
                 xlim = xlim_plot,
                 ylim = ylim_plot,
                 xlab = "Pihat",
                 ylab = "Pihat2",
                 main = paste("Scatter plot (type=", type_val, "): Pihat vs Pihat2\n(green: TP, red: FN)", sep=""))
            
            # Add grid with 0.02 spacing
            x_grid <- seq(floor(xlim_plot[1] / 0.02) * 0.02, ceiling(xlim_plot[2] / 0.02) * 0.02, by = 0.02)
            y_grid <- seq(floor(ylim_plot[1] / 0.02) * 0.02, ceiling(ylim_plot[2] / 0.02) * 0.02, by = 0.02)
            abline(v = x_grid, col = "lightgray", lty = "dotted", lwd = 0.5)
            abline(h = y_grid, col = "lightgray", lty = "dotted", lwd = 0.5)
            
            # Add points on top of grid
            points(pihat_vals[type_filter_plot], pihat2_vals[type_filter_plot],
                   col = colors_plot,
                   pch = 19,
                   cex = 0.8)
            
            # Add legend
            n_tp_plot <- sum(tp_fn_status_plot[type_filter_plot] == "TP", na.rm = TRUE)
            n_fn_plot <- sum(tp_fn_status_plot[type_filter_plot] == "FN", na.rm = TRUE)
            n_other_plot <- sum(is.na(tp_fn_status_plot[type_filter_plot]) | (!tp_fn_status_plot[type_filter_plot] %in% c("TP", "FN")))
            legend_items_plot <- c()
            legend_cols_plot <- c()
            if (n_tp_plot > 0) {
                legend_items_plot <- c(legend_items_plot, paste("TP (n=", n_tp_plot, ")", sep=""))
                legend_cols_plot <- c(legend_cols_plot, "green")
            }
            if (n_fn_plot > 0) {
                legend_items_plot <- c(legend_items_plot, paste("FN (n=", n_fn_plot, ")", sep=""))
                legend_cols_plot <- c(legend_cols_plot, "red")
            }
            if (n_other_plot > 0) {
                legend_items_plot <- c(legend_items_plot, paste("Other (n=", n_other_plot, ")", sep=""))
                legend_cols_plot <- c(legend_cols_plot, "gray")
            }
            if (length(legend_items_plot) > 0) {
                legend("topright",
                       legend = legend_items_plot,
                       col = legend_cols_plot,
                       pch = 19,
                       cex = 1.3)
            }
            
            cat("Pihat vs Pihat2 scatter plot (type=", type_val, "): ", sum(type_filter_plot), 
                " points (TP=", n_tp_plot, ", FN=", n_fn_plot, ")\n", sep="")
        } else {
            cat("No observations with type=", type_val, " and non-NA Pihat/Pihat2 values\n", sep="")
        }
    }
    
    # Create scatter plot for type 100 (Control), color-coded by detected class (TN vs FP)
    type_filter_plot_100 <- type_vals_plot == 100 & !is.na(pihat_vals) & !is.na(pihat2_vals)
    
    if (sum(type_filter_plot_100) > 0) {
        par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1, cex.lab = 1.35, cex.axis = 1.2, cex.main = 1.25)
        
        # Get colors based on TN/FP status
        # TN (predicted as Control) = blue, FP (predicted as Case) = orange/red
        colors_plot_100 <- rep("gray", sum(type_filter_plot_100))
        colors_plot_100[tn_fp_status_plot[type_filter_plot_100] == "TN"] <- "blue"
        colors_plot_100[tn_fp_status_plot[type_filter_plot_100] == "FP"] <- "orange"
        
        # Determine plot limits
        xlim_plot_100 <- range(pihat_vals[type_filter_plot_100], na.rm = TRUE)
        ylim_plot_100 <- range(pihat2_vals[type_filter_plot_100], na.rm = TRUE)
        
        # Create empty plot
        plot(0, 0, type = "n",
             xlim = xlim_plot_100,
             ylim = ylim_plot_100,
             xlab = "Pihat",
             ylab = "Pihat2",
             main = "Scatter plot (type=100, Control): Pihat vs Pihat2\n(blue: TN, orange: FP)")
        
        # Add grid with 0.02 spacing
        x_grid_100 <- seq(floor(xlim_plot_100[1] / 0.02) * 0.02, ceiling(xlim_plot_100[2] / 0.02) * 0.02, by = 0.02)
        y_grid_100 <- seq(floor(ylim_plot_100[1] / 0.02) * 0.02, ceiling(ylim_plot_100[2] / 0.02) * 0.02, by = 0.02)
        abline(v = x_grid_100, col = "lightgray", lty = "dotted", lwd = 0.5)
        abline(h = y_grid_100, col = "lightgray", lty = "dotted", lwd = 0.5)
        
        # Add points on top of grid
        points(pihat_vals[type_filter_plot_100], pihat2_vals[type_filter_plot_100],
               col = colors_plot_100,
               pch = 19,
               cex = 0.8)
        
        # Add legend
        n_tn_plot_100 <- sum(tn_fp_status_plot[type_filter_plot_100] == "TN", na.rm = TRUE)
        n_fp_plot_100 <- sum(tn_fp_status_plot[type_filter_plot_100] == "FP", na.rm = TRUE)
        n_other_plot_100 <- sum(is.na(tn_fp_status_plot[type_filter_plot_100]) | (!tn_fp_status_plot[type_filter_plot_100] %in% c("TN", "FP")))
        legend_items_plot_100 <- c()
        legend_cols_plot_100 <- c()
        if (n_tn_plot_100 > 0) {
            legend_items_plot_100 <- c(legend_items_plot_100, paste("TN (n=", n_tn_plot_100, ")", sep=""))
            legend_cols_plot_100 <- c(legend_cols_plot_100, "blue")
        }
        if (n_fp_plot_100 > 0) {
            legend_items_plot_100 <- c(legend_items_plot_100, paste("FP (n=", n_fp_plot_100, ")", sep=""))
            legend_cols_plot_100 <- c(legend_cols_plot_100, "orange")
        }
        if (n_other_plot_100 > 0) {
            legend_items_plot_100 <- c(legend_items_plot_100, paste("Other (n=", n_other_plot_100, ")", sep=""))
            legend_cols_plot_100 <- c(legend_cols_plot_100, "gray")
        }
        if (length(legend_items_plot_100) > 0) {
            legend("topright",
                   legend = legend_items_plot_100,
                   col = legend_cols_plot_100,
                   pch = 19,
                   cex = 1.3)
        }
        
        cat("Pihat vs Pihat2 scatter plot (type=100, Control): ", sum(type_filter_plot_100), 
            " points (TN=", n_tn_plot_100, ", FP=", n_fp_plot_100, ")\n", sep="")
    } else {
        cat("No observations with type=100 and non-NA Pihat/Pihat2 values\n")
    }
} else {
    cat("Warning: Required columns (Pihat, Pihat2, type) not found for Pihat vs Pihat2 scatter plots\n")
}

dev.off()
cat("\nAll plots saved in:", pdf_file, "\n")
cat("Total pages:", length(variable_list) * 2, "\n")

# ============================================================================
# 7. Final output table: valid rows + type_norm, predicted_cluster, P_HS_supervised, P_HS_mixture
# ============================================================================
# Writes GTNA_predictions.txt with all input columns plus predictions. Supervised
# and mixture columns are NA where features are incomplete.
cat("\n", rep("=", 70), "\n", sep = "")
cat("Generating output file for all valid rows\n")
cat(rep("=", 70), "\n", sep = "")

# idx_valid was set in section 3 (all 7 variables non-NA). Valid rows = data[idx_valid, ].
idx_valid_rows <- which(idx_valid)
cat("Number of valid rows after filtering and validation:", length(idx_valid_rows), "\n")

if (length(idx_valid_rows) > 0) {
    # Get the first two variables for prediction (var1 and var2)
    # var1_clean and var2_clean are already filtered by idx_valid
    var1_for_pred <- var1_clean
    var2_for_pred <- var2_clean
    
    # Get type data for these rows
    type_for_pred <- type_valid
    
    # Calculate type_norm: 2,7,11 → 2; 100 → 3; others → 0
    type_norm <- rep(0, length(type_for_pred))
    type_norm[type_for_pred %in% c(2, 7, 11)] <- 2
    type_norm[type_for_pred == 100] <- 3
    type_norm[is.na(type_for_pred)] <- 0
    
    # Calculate true_label: 2,3,4,5,7,11 → 1; 100 → 0; others → NA
    true_label <- rep(NA, length(type_for_pred))
    true_label[type_for_pred %in% c(2, 7, 11)] <- 1
    true_label[type_for_pred == 100] <- 0
    
    # Get GMM results for var1 (first variable) - use GMM 2 (all data)
    predicted_cluster <- rep(NA, length(type_for_pred))
    if (length(variable_list) >= 1) {
        var1_info <- variable_list[[1]]
        var1_data_all <- var1_info$data
        gmm_results_var1 <- calculate_gmm_for_variable(var1_data_all, var1_info$name)
        
        if (!is.null(gmm_results_var1$gmm_all)) {
            # Predict clusters for var1 values
            prob_cluster1 <- gmm_results_var1$gmm_all$lambda[1] * dnorm(var1_for_pred, 
                                                                        mean = gmm_results_var1$gmm_all$mu[1], 
                                                                        sd = gmm_results_var1$gmm_all$sigma[1])
            prob_cluster2 <- gmm_results_var1$gmm_all$lambda[2] * dnorm(var1_for_pred, 
                                                                        mean = gmm_results_var1$gmm_all$mu[2], 
                                                                        sd = gmm_results_var1$gmm_all$sigma[2])
            total_prob <- prob_cluster1 + prob_cluster2
            prob_cluster1 <- prob_cluster1 / total_prob
            prob_cluster2 <- prob_cluster2 / total_prob
            
            predicted_cluster <- ifelse(prob_cluster1 > prob_cluster2, 1, 2)
        }
    }
    
    # Create output data frame with valid rows and new columns
    # data[idx_valid, ] gives us the valid rows
    output_data <- data[idx_valid, ]
    
    # Add new columns
    output_data$type_norm <- type_norm
    output_data$predicted_cluster <- predicted_cluster
    output_data$true_label <- true_label
    output_data$variable1_value <- var1_for_pred
    output_data$variable2_value <- var2_for_pred

    # Merge supervised and mixture predictions (NA where 4-D features missing).
    output_data$P_HS_supervised <- NA_real_
    output_data$pred_class_supervised <- NA_integer_
    if (!is.null(res_sup) && any(res_sup$ok_all)) {
        output_data$P_HS_supervised[res_sup$ok_all] <- res_sup$prob_all
        output_data$pred_class_supervised[res_sup$ok_all] <- res_sup$pred_class_supervised
    }

    # Mixture columns: P_HS_mixture and pred_class_mixture (NA where 4-D incomplete).
    output_data$P_HS_mixture <- NA_real_
    output_data$pred_class_mixture <- NA_integer_
    if (!is.null(res_mix) && any(res_mix$ok_all)) {
        output_data$P_HS_mixture[res_mix$ok_all] <- res_mix$prob_all
        output_data$pred_class_mixture[res_mix$ok_all] <- res_mix$pred_class_mixture
    }

    # Write output file
    output_file <- "/pl/active/KellerLab/Emmanuel/NAvHS/GTNA_predictions.txt"
    write.table(output_data, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
    cat("Output file saved:", output_file, "\n")
    cat("Number of rows:", nrow(output_data), "\n")
    cat("Columns:", paste(colnames(output_data), collapse = ", "), "\n")
} else {
    cat("No valid rows found after filtering. Output file not created.\n")
}

cat("\nScript finished.\n")
