#' Helper function to filter out sparse proteins
filter_features <- function(DT_limma, control_samples, treatment_samples, alpha){
  n_treatment <- length(treatment_samples)
  n_control <- length(control_samples)
  if (n_treatment>3&n_control>3) {
    # Drop rows (protein groups) with > 50% missingness in samples
    DT_limma$missing_value= apply(DT_limma, 1, function(x) sum(x==0))
    DT_limma$missing_value_c= apply(DT_limma[,control_samples], 1, function(x) sum(x==0))
    DT_limma$missing_value_t= apply(DT_limma[,treatment_samples], 1, function(x) sum(x==0))
    DT_limma <- DT_limma[!(DT_limma$missing_value_t >(n_treatment * alpha)&DT_limma$missing_value_t < n_treatment),]
    DT_limma <- DT_limma[!(DT_limma$missing_value_c >(n_control * alpha)&DT_limma$missing_value_c <n_control),]
    DT_limma <- DT_limma[DT_limma$missing_value != (n_treatment+n_control),]
    DT_limma[,grep('missing_value',colnames(DT_limma))]=NULL
  } else{
    # Drop rows (protein groups) with missingness in samples
    DT_limma$missing_value= apply(DT_limma, 1, function(x) sum(x==0))
    DT_limma$missing_value_c= apply(DT_limma[,control_samples], 1, function(x) sum(x==0))
    DT_limma$missing_value_t= apply(DT_limma[,treatment_samples], 1, function(x) sum(x==0))
    DT_limma <- DT_limma[DT_limma$missing_value_t %in% c(0,n_treatment),]
    DT_limma <- DT_limma[DT_limma$missing_value_c %in% c(0,n_control),]
    DT_limma <- DT_limma[DT_limma$missing_value != (n_treatment+n_control),]
    DT_limma[,grep('missing_value',colnames(DT_limma))]=NULL
  }
  return(DT_limma)

}

#' Perform limma differential expression on a SummarizedExperiment
#'
#' This function takes a `SummarizedExperiment` object and uses group labels in
#' `colData(object)` to perform differential expression between a treatment group
#' and a control group. Covariates can be included in the model design.
#'
#' @param object A `SummarizedExperiment` object containing protein intensities and metadata.
#' @param condition String: column name in `colData(object)` that holds group labels.
#' @param treatment_group String: name of treatment group in `condition`.
#' @param control_group String: name of control group in `condition`.
#' @param covariates Optional character vector: column names in `colData(object)` to use as covariates.
#'
#' @return A data frame with metadata, intensities, log fold change, p-values, and adjusted p-values.
#' @export
#'
setGeneric("do_limma_binary", function(object,
                                condition,
                                treatment_group,
                                control_group,
                                covariates = NULL) {
  standardGeneric("do_limma_binary")
})

setMethod("do_limma_binary", "SummarizedExperiment",
          function(object, condition, treatment_group, control_group, covariates) {
            # -- prepare assay data --
            meta_cols <- names(rowData(object))
            if (ProtPipe::has_step(object, "log2_transform")) {
              data <- assay(object) %>% as.data.frame()
            } else {
              object <- ProtPipe::log2_transform(object)
              data <- assay(object) %>% as.data.frame()
            }
            if (anyNA(data)) stop("Missing values detected. Please impute before running limma.")

            DT <- cbind(rowData(object), data) %>% as.data.frame()
            meta <- colData(object) %>% as.data.frame()

            # -- check condition column --
            if (!(condition %in% colnames(meta))) {
              stop("`condition` must be a column name in colData(object).")
            }
            # -- ensure groups are different --
            if (treatment_group == control_group) {
              stop("Treatment group and control group must be different.")
            }

            # -- covariate check: must not include the grouping variable --
            if (!is.null(covariates) && condition %in% covariates) {
              stop("Covariates cannot include the grouping variable used for treatment/control.")
            }

            # -- extract sample groups --
            groups <- meta[[condition]]
            if (!(treatment_group %in% groups && control_group %in% groups)) {
              stop("Both treatment and control groups must exist in `condition`.")
            }

            treatment_samples <- rownames(meta[meta[[condition]] == treatment_group, , drop = FALSE])
            control_samples   <- rownames(meta[meta[[condition]] == control_group, , drop = FALSE])

            if (length(treatment_samples) < 2 || length(control_samples) < 2) {
              stop("Each group must contain at least 2 samples.")
            }

            DT_limma <- DT[, c(treatment_samples, control_samples)]

            # -- build design matrix --
            group_list <- factor(c(rep("treatment", length(treatment_samples)),
                                   rep("control",   length(control_samples))),
                                 levels = c("control", "treatment"))

            design_df <- data.frame(group = group_list)
            rownames(design_df) <- c(treatment_samples, control_samples)

            # add covariates if provided
            if (!is.null(covariates)) {
              for (cov in covariates) {
                if (!(cov %in% colnames(meta))) {
                  stop(paste("Covariate", cov, "not found in colData(object)."))
                }
                design_df[[cov]] <- meta[rownames(design_df), cov]
              }
            }

            limma_design <- model.matrix(~ 0 + ., data = design_df)
            colnames(limma_design) <- make.names(colnames(limma_design))

            # -- contrast matrix --
            cont.matrix <- limma::makeContrasts(grouptreatment - groupcontrol, levels = limma_design)

            # -- limma pipeline --
            fit <- limma::lmFit(DT_limma, limma_design)
            fit2 <- limma::contrasts.fit(fit, cont.matrix)
            fit2 <- limma::eBayes(fit2, trend = TRUE)

            result_limma <- limma::topTable(fit2, coef = 1, n = Inf)
            result_limma <- merge(DT[, c(meta_cols, treatment_samples, control_samples)],
                                  result_limma, by.x = 0, by.y = 0)
            result_limma$Row.names <- NULL

            # sort by adjusted p-value then effect size
            result_limma_sorted <- result_limma[order(result_limma$adj.P.Val, -abs(result_limma$logFC)), ]
            return(result_limma_sorted)
          })

#' Perform t-test differential expression on a SummarizedExperiment
#'
#' This function takes a `SummarizedExperiment` object and uses group labels in
#' `colData(object)` to perform differential expression between a treatment group
#' and a control group using per-protein t-tests. Optional covariates can be
#' regressed out before testing.
#'
#' @param object A `SummarizedExperiment` object containing protein intensities and metadata.
#' @param condition String: column name in `colData(object)` that holds group labels.
#' @param treatment_group String: name of treatment group in `condition`.
#' @param control_group String: name of control group in `condition`.
#' @param covariates Optional character vector: column names in `colData(object)` to regress out before testing.
#'
#' @return A data frame with metadata, intensities, log fold change, p-values, and adjusted p-values.
#' @export
#'
setGeneric("do_t_test_binary", function(object,
                                       condition,
                                       treatment_group,
                                       control_group,
                                       covariates = NULL) {
  standardGeneric("do_t_test_binary")
})

setMethod("do_t_test_binary", "SummarizedExperiment",
          function(object, condition, treatment_group, control_group, covariates) {
            meta_cols <- names(rowData(object))
            if (ProtPipe::has_step(object, "log2_transform")) {
              data <- assay(object) %>% as.data.frame()
            } else {
              object <- ProtPipe::log2_transform(object)
              data <- assay(object) %>% as.data.frame()
            }
            if (anyNA(data)) stop("Missing values detected. Please impute before running t-tests.")

            DT <- cbind(rowData(object), data) %>% as.data.frame()
            meta <- colData(object) %>% as.data.frame()

            if (!(condition %in% colnames(meta))) {
              stop("`condition` must be a column name in colData(object).")
            }
            if (treatment_group == control_group) {
              stop("Treatment group and control group must be different.")
            }
            if (!is.null(covariates) && condition %in% covariates) {
              stop("Covariates cannot include the grouping variable used for treatment/control.")
            }

            groups <- meta[[condition]]
            if (!(treatment_group %in% groups && control_group %in% groups)) {
              stop("Both treatment and control groups must exist in `condition`.")
            }

            treatment_samples <- rownames(meta[meta[[condition]] == treatment_group, , drop = FALSE])
            control_samples <- rownames(meta[meta[[condition]] == control_group, , drop = FALSE])

            if (length(treatment_samples) < 2 || length(control_samples) < 2) {
              stop("Each group must contain at least 2 samples.")
            }

            selected_samples <- c(treatment_samples, control_samples)
            feature_ids <- rownames(object)
            test_data <- as.matrix(DT[, selected_samples, drop = FALSE])
            if (is.null(feature_ids) || length(feature_ids) != nrow(test_data)) {
              feature_ids <- seq_len(nrow(test_data))
            }

            if (!is.null(covariates)) {
              for (cov in covariates) {
                if (!(cov %in% colnames(meta))) {
                  stop(paste("Covariate", cov, "not found in colData(object)."))
                }
              }

              covariate_df <- meta[selected_samples, covariates, drop = FALSE]
              covariate_design <- stats::model.matrix(~ ., data = covariate_df)

              residualized <- t(apply(test_data, 1, function(y) {
                fit <- stats::lm.fit(x = covariate_design, y = as.numeric(y))
                stats::residuals(fit) + stats::coef(fit)[1]
              }))
              colnames(residualized) <- selected_samples
              rownames(residualized) <- rownames(test_data)
              test_data <- residualized
            }

            test_results <- lapply(seq_len(nrow(test_data)), function(i) {
              treat_vals <- as.numeric(test_data[i, treatment_samples])
              ctrl_vals <- as.numeric(test_data[i, control_samples])

              tt <- stats::t.test(treat_vals, ctrl_vals)
              data.frame(
                Row.names = feature_ids[i],
                logFC = mean(treat_vals) - mean(ctrl_vals),
                AveExpr = mean(c(treat_vals, ctrl_vals)),
                t = unname(tt$statistic),
                P.Value = tt$p.value,
                B = NA_real_,
                stringsAsFactors = FALSE
              )
            })

            result_ttest <- do.call(rbind, test_results)
            result_ttest$adj.P.Val <- stats::p.adjust(result_ttest$P.Value, method = "BH")

            result_ttest <- merge(
              DT[, c(meta_cols, treatment_samples, control_samples), drop = FALSE],
              result_ttest,
              by.x = 0,
              by.y = "Row.names"
            )
            result_ttest$Row.names <- NULL

            result_ttest_sorted <- result_ttest[order(result_ttest$adj.P.Val, -abs(result_ttest$logFC)), ]
            return(result_ttest_sorted)
          })


#' Perform limma differential expression for a continuous outcome
#'
#' This function takes a SummarizedExperiment object, a numeric condition column,
#' and optional covariates, and performs limma-based differential expression.
#'
#' @param object A SummarizedExperiment object
#' @param condition Column name in colData(object) used as the continuous predictor
#'
#' @return A data frame with metadata, intensities, logFC, p-values
#' @export
#'
setGeneric("do_comparison_continuous",
           function(object, condition) standardGeneric("do_comparison_continuous"))
setMethod("do_comparison_continuous", "SummarizedExperiment",
          function(object, condition) {

            # -- Prepare assay data --
            meta_cols <- names(rowData(object))
            if (ProtPipe::has_step(object, "log2_transform")) {
              data <- assay(object) %>% as.data.frame()
            } else {
              object <- ProtPipe::log2_transform(object)
              data <- assay(object) %>% as.data.frame()
            }

            if (anyNA(data))
              stop("Missing values detected. Please impute before running correlation analysis.")

            DT <- cbind(rowData(object), data) %>% as.data.frame()
            meta <- colData(object) %>% as.data.frame()

            # -- Check condition column --
            if (!(condition %in% colnames(meta))) {
              stop("`condition` must be a column name in colData(object).")
            }

            # -- Condition must be numeric --
            if (!is.numeric(meta[[condition]])) {
              stop("For continuous outcome, `condition` must be numeric.")
            }

            # -- Run Spearman correlations for each protein --
            cor_results <- apply(data, 1, function(x) {
              test <- suppressWarnings(cor.test(x, meta[[condition]], method = "spearman"))
              c(rho = unname(test$estimate), pval = test$p.value)
            })

            cor_df <- as.data.frame(t(cor_results)) %>%
              dplyr::rename(P.Value = pval)
            cor_df$protein <- rownames(cor_df)

            # -- Add adjusted p-values --
            cor_df$adj.P.Val <- p.adjust(cor_df$P.Value, method = "BH")

            # -- Merge with metadata (if needed) --
            cor_df <- merge(DT, cor_df, by.x = "row.names", by.y = "protein", all.x = TRUE)
            cor_df$Row.names <- NULL
            # -- Sort by adjusted p-value (then by |rho|) --
            cor_df_sorted <- cor_df[order(cor_df$adj.P.Val, -abs(cor_df$rho)), ]

            return(cor_df_sorted)
          })










#' Plot a Volcano Plot for Differential Expression Results
#'
#' This function generates a volcano plot based on log fold changes (logFC) and
#' adjusted p-values from differential expression analysis. It highlights genes
#' that pass the specified thresholds for logFC and FDR (false discovery rate).
#' Optionally, it can label genes of interest or the top up/downregulated genes.
#'
#' @param DT.original A data frame containing the differential expression results.
#'        It should have columns `logFC` (log fold change) and `adj.P.Val` (adjusted p-value).
#' @param label_col The column name (as a string) used for labeling genes in the plot. Default is `NULL`.
#'        If `NULL`, the function will select the first column of `DT.original` for labeling.
#' @param lfc_threshold A numeric value representing the threshold for log fold change (default is 1).
#'        Genes with logFC greater than or equal to this value are labeled as "UP",
#'        and genes with logFC less than or equal to the negative of this value are labeled as "DOWN".
#' @param fdr_threshold A numeric value representing the false discovery rate (FDR) threshold (default is 0.01).
#'        Genes with an adjusted p-value greater than or equal to this threshold will be labeled as "Others".
#' @param labelgene A character vector of gene names to be labeled in the plot (default is `NULL`).
#'        If provided, only these genes will be labeled in the plot.
#' @param adj A boolean. Set to true to use adj.P.val for the y axis (default) and false to
#'        use P.value.
#' @return A `ggplot2` object representing the volcano plot.
#' @export
#'
plot_volcano <- function(DT.original, label_col = NULL, lfc_threshold=1, fdr_threshold=0.01, labelgene=NULL, adj=T) {

  if(is.null(label_col)){
    label_col = names(DT.original)[1]
  }

  options(ggrepel.max.overlaps = Inf)
  DT <- DT.original

  # 1. Determine which column to use based on the 'adj' parameter
  target_p_col <- if(adj) "adj.P.Val" else "P.Value"

  # Check if column exists to prevent crashing
  if(!target_p_col %in% names(DT)) stop(paste("Column", target_p_col, "not found in input data."))

  # Set initial group to 'Others' and update based on thresholds
  DT <- DT %>%
    dplyr::mutate(Group = 'Others',
                  Group = dplyr::if_else(logFC >= lfc_threshold, 'UP', Group),
                  Group = dplyr::if_else(logFC <= -lfc_threshold, 'DOWN', Group),
                  # 2. Use .data[[target_p_col]] to dynamically access the p-value column
                  Group = dplyr::if_else(.data[[target_p_col]] >= fdr_threshold, 'Others', Group),
                  labeltext = '')

  # If labelgene is provided, update labeltext accordingly
  if (!is.null(labelgene)) {
    DT <- DT %>%
      dplyr::mutate(labeltext = dplyr::if_else(!!rlang::sym(label_col) %in% labelgene, !!rlang::sym(label_col), labeltext))
  } else{
    # Select top 5 genes for UP and DOWN groups
    up_rows <- which(DT$Group == "UP")
    sorted_up_indices <- up_rows[order(DT$logFC[up_rows], decreasing = TRUE)]
    top_up_indices <- head(sorted_up_indices, 5)

    down_rows <- which(DT$Group == "DOWN")
    sorted_down_indices <- down_rows[order(DT$logFC[down_rows], decreasing = FALSE)]
    top_down_indices <- head(sorted_down_indices, 5)

    top_indices <- c(top_up_indices, top_down_indices)
    DT$labeltext[top_indices] <- DT[top_indices, label_col]
  }

  # plot
  # 3. Updated y mapping to use the dynamic column
  g <- ggplot2::ggplot(DT, ggplot2::aes(x = logFC, y = -log10(.data[[target_p_col]]))) +
    ggplot2::geom_point(ggplot2::aes(color = Group)) +
    ggplot2::scale_color_manual(breaks = c("DOWN", "Others", "UP"),
                                values = c("#67a9cf", "#969696", "#ef8a62")) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(legend.position = "bottom") +
    ggrepel::geom_label_repel(
      data = DT,
      ggplot2::aes(label = labeltext),
      size = 5,
      box.padding = grid::unit(0.35, "lines"),
      point.padding = grid::unit(0.3, "lines")) +
    ggplot2::geom_hline(yintercept = -log10(fdr_threshold), linetype = "dashed") +
    ggplot2::geom_vline(xintercept = lfc_threshold, linetype = "dashed") +
    ggplot2::geom_vline(xintercept = -lfc_threshold, linetype = "dashed") +
    ggplot2::theme_classic() +
    # 4. Optional: Update label so you know what you are looking at
    ggplot2::ylab(paste0("-log10(", target_p_col, ")"))+
    ggplot2::xlab(paste0("log2FC"))

  return(g)
}


#' Plot a Volcano Plot for Correlation Results
#'
#' This function generates a volcano plot based on the correlation coefficient rho and
#' adjusted p-values. It highlights genes
#' that pass the specified thresholds for rho and FDR (false discovery rate).
#' Optionally, it can label genes of interest or the top up/downregulated genes.
#'
#' @param DT.original A data frame containing the differential expression results.
#'        It should have columns `logFC` (log fold change) and `adj.P.Val` (adjusted p-value).
#' @param label_col The column name (as a string) used for labeling genes in the plot. Default is `NULL`.
#'        If `NULL`, the function will select the first column of `DT.original` for labeling.
#' @param rho_threshold A numeric value representing the threshold for rho (default is 0.35).
#'        Genes with logFC greater than or equal to this value are labeled as "UP",
#'        and genes with logFC less than or equal to the negative of this value are labeled as "DOWN".
#' @param fdr_threshold A numeric value representing the false discovery rate (FDR) threshold (default is 0.01).
#'        Genes with an adjusted p-value greater than or equal to this threshold will be labeled as "Others".
#' @param labelgene A character vector of gene names to be labeled in the plot (default is `NULL`).
#'        If provided, only these genes will be labeled in the plot.
#' @param adj A boolean. Set to true to use adj.P.val for the y axis (default) and false to
#'        use P.value.
#' @return A `ggplot2` object representing the volcano plot.
#' @export
#'
plot_correlation_volcano <- function(DT.original, label_col = NULL, rho_threshold = 0.35, fdr_threshold = 0.01, labelgene = NULL, adj = T) {

  # 1. Select the column to use based on the 'adj' parameter
  target_p_col <- if (adj) "adj.P.Val" else "P.Value"

  if(is.null(label_col)){
    label_col = names(DT.original)[1]
  }

  options(ggrepel.max.overlaps = Inf)
  DT <- DT.original

  # 2. Update Grouping logic to use the dynamic 'target_p_col'
  # We use .data[[target_p_col]] to access the column by string name within dplyr
  DT <- DT %>%
    dplyr::mutate(Group = 'Others',
                  Group = dplyr::if_else(rho >= rho_threshold, "Positive", Group),
                  Group = dplyr::if_else(rho <= -rho_threshold, "Negative", Group),
                  Group = dplyr::if_else(.data[[target_p_col]] >= fdr_threshold, 'Others', Group),
                  labeltext = '')

  # If labelgene is provided, update labeltext accordingly
  if (!is.null(labelgene)) {
    DT <- DT %>%
      dplyr::mutate(labeltext = dplyr::if_else(!!rlang::sym(label_col) %in% labelgene, !!rlang::sym(label_col), labeltext))
  } else {
    # Select top 5 genes for UP and DOWN groups
    up_rows <- which(DT$Group == "Positive")
    sorted_up_indices <- up_rows[order(DT$rho[up_rows], decreasing = TRUE)]
    top_up_indices <- head(sorted_up_indices, 5)

    down_rows <- which(DT$Group == "Negative")
    sorted_down_indices <- down_rows[order(DT$rho[down_rows], decreasing = FALSE)]
    top_down_indices <- head(sorted_down_indices, 5)

    top_indices <- c(top_up_indices, top_down_indices)
    DT$labeltext[top_indices] <- DT[top_indices, label_col]
  }

  # 3. Plot update:
  # - Y-axis uses the dynamic column
  # - Added labs(y=...) to change the label text
  g <- ggplot2::ggplot(DT, ggplot2::aes(x = rho, y = -log10(.data[[target_p_col]]))) +
    ggplot2::geom_point(ggplot2::aes(color = Group)) +
    ggplot2::scale_color_manual(breaks = c("Negative", "Others", "Positive"),
                                values = c("#67a9cf", "#969696", "#ef8a62")) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(legend.position = "bottom") +
    ggrepel::geom_label_repel(
      data = DT,
      ggplot2::aes(label = labeltext),
      size = 5,
      box.padding = grid::unit(0.35, "lines"),
      point.padding = grid::unit(0.3, "lines")) +
    ggplot2::geom_hline(yintercept = -log10(fdr_threshold), linetype = "dashed") +
    ggplot2::geom_vline(xintercept = rho_threshold, linetype = "dashed") +
    ggplot2::geom_vline(xintercept = -rho_threshold, linetype = "dashed") +
    ggplot2::labs(y = paste0("-log10(", target_p_col, ")")) + # Dynamic Label
    ggplot2::theme_classic()
  return(g)
}


#' Add Entrez Gene IDs to a Data Frame
#'
#' @description Maps gene symbols from a specified column to Entrez IDs using
#' clusterProfiler::bitr and joins them back to the original data frame. It also
#' cleans gene symbols that may contain multiple entries separated by semicolons.
#'
#' @param DE A data frame, such as one containing differential expression results.
#' @param org An organism annotation database object (e.g., org.Hs.eg.db).
#' @param gene_col A character string specifying the name of the column in DE
#'   that contains the gene symbols.
#'
#' @return A data frame with an appended ENTREZID column containing the mapped
#'   IDs. Returns NULL if no genes can be successfully mapped.
#'
#' @export
#'
#' @importFrom clusterProfiler bitr
#' @importFrom dplyr inner_join
#'
add_entrez <- function(DE, org = org.Hs.eg.db::org.Hs.eg.db, gene_col = "Genes") {
  DE[[gene_col]] <- sapply(strsplit(DE[[gene_col]], ";"), `[`, 1)

  # Remove NA or empty strings before mapping
  genes_to_map <- DE[[gene_col]]
  genes_to_map <- genes_to_map[!is.na(genes_to_map) & genes_to_map != ""]

  if (length(genes_to_map) == 0) {
    warning("No valid gene symbols found to map.")
    return(NULL)
  }

  # Use tryCatch to prevent app crash
  mapping <- tryCatch({
    clusterProfiler::bitr(
      genes_to_map,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org
    )
  }, error = function(e) {
    warning("Gene mapping failed: ", conditionMessage(e))
    return(NULL)
  })

  if (is.null(mapping) || nrow(mapping) == 0) {
    warning("No Entrez IDs mapped.")
    return(NULL)
  }

  dplyr::inner_join(
    DE,
    mapping,
    by = setNames("SYMBOL", gene_col)
  )
}

#' Read a Custom Ontology File
#'
#' @description
#' Reads a user-supplied ontology into the `TERM2GENE` and optional `TERM2NAME`
#' structures required by `clusterProfiler::enricher()` and
#' `clusterProfiler::GSEA()`.
#'
#' Supported formats are:
#' - GMT files, where each line is `term`, `name`, then one or more gene IDs
#' - Delimited text files (`csv` or `tsv`) with either:
#'   - 2 columns: `term`, `gene`
#'   - 3 columns: `term`, `name`, `gene`
#'
#' Gene identifiers in uploaded ontologies must match the identifiers used for
#' enrichment. In the ProtPipe Shiny app, this currently means Entrez gene IDs.
#'
#' @param file Path to a `.gmt`, `.csv`, or `.tsv` ontology file.
#'
#' @return A named list with `term2gene` and `term2name`.
#' @export
read_ontology <- function(file) {
  ext <- tolower(tools::file_ext(file))

  if (ext == "gmt") {
    lines <- readLines(file, warn = FALSE)
    lines <- lines[nzchar(lines)]

    parsed <- lapply(lines, function(line) strsplit(line, "\t", fixed = TRUE)[[1]])
    parsed <- parsed[lengths(parsed) >= 3]

    if (length(parsed) == 0) {
      stop("The ontology file did not contain any valid GMT records.")
    }

    term2gene <- do.call(rbind, lapply(parsed, function(x) {
      data.frame(term = x[[1]], gene = x[-c(1, 2)], stringsAsFactors = FALSE)
    }))
    term2name <- unique(do.call(rbind, lapply(parsed, function(x) {
      data.frame(term = x[[1]], name = x[[2]], stringsAsFactors = FALSE)
    })))
  } else if (ext %in% c("csv", "tsv")) {
    sep <- if (ext == "csv") "," else "\t"
    dat <- data.table::fread(file, sep = sep, data.table = FALSE)

    if (ncol(dat) < 2) {
      stop("Delimited ontology files must contain at least 2 columns.")
    }

    if (ncol(dat) == 2) {
      term2gene <- dat[, 1:2, drop = FALSE]
      colnames(term2gene) <- c("term", "gene")
      term2name <- NULL
    } else {
      term2gene <- dat[, c(1, 3), drop = FALSE]
      colnames(term2gene) <- c("term", "gene")
      term2name <- unique(dat[, 1:2, drop = FALSE])
      colnames(term2name) <- c("term", "name")
    }
  } else {
    stop("Unsupported ontology format. Use .gmt, .csv, or .tsv.")
  }

  term2gene <- unique(term2gene[!is.na(term2gene$term) & !is.na(term2gene$gene), , drop = FALSE])
  if (nrow(term2gene) == 0) {
    stop("The ontology file did not contain any valid term-to-gene mappings.")
  }

  if (!is.null(term2name)) {
    term2name <- unique(term2name[!is.na(term2name$term) & !is.na(term2name$name), , drop = FALSE])
    if (nrow(term2name) == 0) {
      term2name <- NULL
    }
  }

  list(term2gene = term2gene, term2name = term2name)
}

#' Run Over-Representation Analysis for Custom Gene Sets
#'
#' @description
#' A thin wrapper around `clusterProfiler::enricher()` for user-supplied gene
#' sets represented as `TERM2GENE` and optional `TERM2NAME` tables.
#'
#' @param gene_id A character vector of significant gene IDs.
#' @param all_gene_vector A character vector containing the background universe.
#' @param term2gene A two-column data frame with term IDs in the first column and
#'   gene IDs in the second.
#' @param term2name An optional two-column data frame with term IDs and readable
#'   term names.
#' @param enrich_pvalue The p-value and q-value cutoff. Defaults to `1`.
#'
#' @return An `enrichResult` object, or `NULL` if no terms are enriched.
#' @export
enrich_terms <- function(gene_id, all_gene_vector, term2gene, term2name = NULL, enrich_pvalue = 1) {
  enrichment <- tryCatch({
    clusterProfiler::enricher(
      gene = gene_id,
      universe = unique(all_gene_vector),
      TERM2GENE = term2gene,
      TERM2NAME = term2name,
      pAdjustMethod = "BH",
      pvalueCutoff = enrich_pvalue,
      qvalueCutoff = enrich_pvalue
    )
  }, error = function(e) {
    warning("Custom enrichment failed: ", conditionMessage(e))
    return(NULL)
  })

  if (!is.null(enrichment) && !is.null(enrichment@result) && nrow(enrichment@result) > 0) {
    return(enrichment)
  }

  NULL
}

#' Run GSEA for Custom Gene Sets
#'
#' @description
#' A thin wrapper around `clusterProfiler::GSEA()` for user-supplied gene sets
#' represented as `TERM2GENE` and optional `TERM2NAME` tables.
#'
#' @param gene_list A named numeric vector of ranked genes.
#' @param term2gene A two-column data frame with term IDs in the first column and
#'   gene IDs in the second.
#' @param term2name An optional two-column data frame with term IDs and readable
#'   term names.
#' @param enrich_pvalue The p-value cutoff. Defaults to `1`.
#'
#' @return A `gseaResult` object, or `NULL` if no terms are enriched.
#' @export
gse_terms <- function(gene_list, term2gene, term2name = NULL, enrich_pvalue = 1) {
  enrichment <- tryCatch({
    clusterProfiler::GSEA(
      geneList = gene_list,
      TERM2GENE = term2gene,
      TERM2NAME = term2name,
      pvalueCutoff = enrich_pvalue,
      verbose = FALSE
    )
  }, error = function(e) {
    warning("Custom GSEA failed: ", conditionMessage(e))
    return(NULL)
  })

  if (!is.null(enrichment) && !is.null(enrichment@result) && nrow(enrichment@result) > 0) {
    return(enrichment)
  }

  NULL
}

has_valid_gsea_ranking <- function(gene_list) {
  stats <- as.numeric(gene_list)
  stats <- stats[!is.na(stats)]

  if (length(stats) < 2) {
    return(FALSE)
  }

  if (length(unique(stats)) < 2) {
    return(FALSE)
  }

  if (!any(stats > 0) || !any(stats < 0)) {
    return(FALSE)
  }

  TRUE
}

has_enrichment_signal <- function(stats) {
  stats <- as.numeric(stats)
  stats <- stats[!is.na(stats)]

  if (length(stats) == 0) {
    return(FALSE)
  }

  length(unique(stats)) > 1 && any(stats != 0)
}


#' Perform Gene Ontology (GO) Over-Representation Analysis (ORA)
#'
#' @description A wrapper around clusterProfiler::enrichGO to perform ORA to
#' find enriched GO terms (Biological Process, Molecular Function, Cellular Component).
#'
#' @param gene_id A character vector of significant Entrez gene IDs to test.
#' @param all_gene_vector A character vector of all background (universe) Entrez
#'   gene IDs against which to test.
#' @param enrich_pvalue The p-value and q-value cutoff for enrichment. Defaults to 1.
#' @param org An organism annotation database (e.g., org.Hs.eg.db).
#' @param ont GO ontology to test. One of `"BP"`, `"MF"`, `"CC"`, or `"ALL"`.
#'
#' @return An enrichResult object if the analysis is successful and finds
#'   enriched terms, otherwise NULL.
#'
#' @export
#'
#' @importFrom clusterProfiler enrichGO
#'
enrich_go <- function(gene_id, all_gene_vector, enrich_pvalue = 1, org = org.Hs.eg.db::org.Hs.eg.db, ont = "BP") {
  cat("Processing GO\n")

  GO <- tryCatch({
    clusterProfiler::enrichGO(
      gene          = gene_id,
      universe      = unique(all_gene_vector),
      OrgDb         = org,
      ont           = ont,
      pAdjustMethod = "BH",
      pvalueCutoff  = enrich_pvalue,
      qvalueCutoff  = enrich_pvalue,
      readable      = TRUE
    )
  }, error = function(e) {
    warning("GO enrichment failed: ", conditionMessage(e))
    return(NULL)
  })

  if (!is.null(GO) && !is.null(GO@result) && nrow(GO@result) > 0) {
    return(GO)
  } else {
    return(NULL)
  }
}


#' Perform KEGG Over-Representation Analysis (ORA)
#'
#' @description A wrapper around clusterProfiler::enrichKEGG to perform ORA
#' to find enriched KEGG pathways from a list of significant genes.
#'
#' @param gene_id A character vector of significant Entrez gene IDs to test for enrichment.
#' @param all_gene_vector A character vector of all background (universe) Entrez
#'   gene IDs against which to test.
#' @param enrich_pvalue The p-value and q-value cutoff for enrichment. Defaults to 1.
#' @param org An organism annotation database (e.g., org.Hs.eg.db) used to
#'   map Entrez IDs to gene symbols.
#' @param organism A character string for the KEGG organism code (e.g., 'hsa' for human).
#'
#' @return An enrichResult object if the analysis is successful and finds
#'   enriched pathways, otherwise NULL.
#'
#' @export
#'
#' @importFrom clusterProfiler enrichKEGG
#' @importFrom DOSE setReadable
#'
enrich_kegg <- function(gene_id, all_gene_vector, enrich_pvalue = 1, org = org.Hs.eg.db::org.Hs.eg.db, organism = 'hsa') {
  cat("Processing KEGG\n")

  KEGG <- tryCatch({
    clusterProfiler::enrichKEGG(
      gene         = gene_id,
      organism     = organism,
      universe     = all_gene_vector,
      pvalueCutoff = enrich_pvalue,
      qvalueCutoff = enrich_pvalue
    )
  }, error = function(e) {
    warning("KEGG enrichment failed: ", conditionMessage(e))
    return(NULL)
  })

  if (!is.null(KEGG)) {
    KEGG <- DOSE::setReadable(KEGG, OrgDb = org, keyType = "ENTREZID")
    if (!is.null(KEGG@result) && nrow(KEGG@result) > 0) {
      return(KEGG)
    }
  }

  return(NULL)
}

#' Perform Gene Ontology (GO) Gene Set Enrichment Analysis (GSEA)
#'
#' @description A wrapper around clusterProfiler::gseGO to perform GSEA for
#' Gene Ontology terms on a ranked list of genes.
#'
#' @param gene_list A named numeric vector of genes. Names should be Entrez IDs
#'   and values should be the ranking metric (e.g., log2 fold change).
#' @param enrich_pvalue The p-value cutoff for enrichment. Defaults to 1 (no cutoff).
#' @param org An organism annotation database (e.g., org.Hs.eg.db).
#' @param ont GO ontology to test. One of `"BP"`, `"MF"`, `"CC"`, or `"ALL"`.
#'
#' @return A gseResult object if the analysis is successful and finds
#'   enriched terms, otherwise NULL.
#'
#' @export
#'
#' @importFrom clusterProfiler gseGO
#' @importFrom DOSE setReadable
#'
gse_go <- function(gene_list, enrich_pvalue = 1, org = org.Hs.eg.db::org.Hs.eg.db, ont = "BP") {
  cat("Processing GSEA GO\n")

  GO <- tryCatch({
    clusterProfiler::gseGO(
      geneList     = gene_list,
      OrgDb        = org,
      ont          = ont,
      pAdjustMethod = "BH",
      pvalueCutoff  = enrich_pvalue,
      verbose       = FALSE
    )
  }, error = function(e) {
    warning("GSEA GO failed: ", conditionMessage(e))
    return(NULL)
  })

  if (!is.null(GO)) {
    GO <- DOSE::setReadable(GO, OrgDb = org, keyType = "ENTREZID")
    if (!is.null(GO@result) && nrow(GO@result) > 0) {
      return(GO)
    }
  }

  return(NULL)
}

#' Perform KEGG Gene Set Enrichment Analysis (GSEA)
#'
#' @description A wrapper around clusterProfiler::gseKEGG to perform GSEA on a
#' ranked list of genes and convert IDs to readable gene symbols.
#'
#' @param gene_list A named numeric vector of genes. Names should be Entrez IDs
#'   and values should be the ranking metric (e.g., log2 fold change).
#' @param enrich_pvalue The p-value cutoff for enrichment. Defaults to 1 (no cutoff).
#' @param org An organism annotation database (e.g., org.Hs.eg.db) used to
#'   map Entrez IDs to gene symbols.
#' @param organism A character string for the KEGG organism code (e.g., 'hsa' for human).
#'
#' @return A gseKEGGResult object if the analysis is successful and finds
#'   enriched pathways, otherwise NULL.
#'
#' @export
#'
#' @importFrom clusterProfiler gseKEGG
#' @importFrom DOSE setReadable
#'
gse_kegg <- function(gene_list, enrich_pvalue = 1, org = org.Hs.eg.db::org.Hs.eg.db, organism = 'hsa') {
  cat("Processing GSEA KEGG\n")

  KEGG <- tryCatch({
    clusterProfiler::gseKEGG(
      geneList     = gene_list,
      organism     = organism,
      pvalueCutoff = enrich_pvalue,
      verbose      = FALSE
    )
  }, error = function(e) {
    warning("GSEA KEGG failed: ", conditionMessage(e))
    return(NULL)
  })

  if (!is.null(KEGG)) {
    KEGG <- DOSE::setReadable(KEGG, OrgDb = org, keyType = "ENTREZID")
    if (!is.null(KEGG@result) && nrow(KEGG@result) > 0) {
      return(KEGG)
    }
  }

  return(NULL)
}

###pathway analysis
#' Perform Comprehensive GO and KEGG Pathway Enrichment Analysis
#'
#' @description
#' This function conducts a full suite of pathway enrichment analyses on a given
#' differential expression (DE) result set. It performs Over-Representation
#' Analysis (ORA) for significantly up- and down-regulated genes, as well as
#' Gene Set Enrichment Analysis (GSEA) on the complete ranked list of genes.
#' It analyzes both Gene Ontology (GO) terms and KEGG pathways.
#'
#' @param DE A data frame containing differential expression results. Must include
#'   columns for gene identifiers, log fold change, and adjusted p-values.
#' @param lfc_threshold A numeric value for the absolute log2 fold change
#'   threshold to define significant genes for ORA. Defaults to `1`.
#' @param fdr_threshold A numeric value for the adjusted p-value (FDR) threshold
#'   to define significant genes for ORA. Defaults to `0.01`.
#' @param enrich_pvalue A numeric value for the p-value cutoff used to determine
#'   significant enrichment for pathways/terms. Defaults to `0.05`.
#' @param go_org An annotation database object (e.g., `org.Hs.eg.db` for human)
#'   used for GO analysis and mapping gene IDs.
#' @param kegg_org A character string specifying the KEGG organism code (e.g.,
#'   `'hsa'` for Homo sapiens).
#' @param gene_col A character string indicating the name of the column in the `DE`
#'   data frame that contains the gene symbols/identifiers. Defaults to `"Genes"`.
#' @param source Pathway source. Use `"go"` for Gene Ontology or `"custom"` for
#'   user-supplied ontologies.
#' @param go_ont GO ontology to test when `source = "go"`. Defaults to `"BP"`.
#' @param term2gene Optional two-column data frame of term-to-gene mappings for
#'   custom ontologies. Gene IDs must match the IDs used for enrichment,
#'   currently Entrez IDs in the Shiny app workflow.
#' @param term2name Optional two-column data frame of term-to-name mappings for
#'   custom ontologies.
#' @param run_ora Logical; run over-representation analysis.
#' @param run_gsea Logical; run gene set enrichment analysis.
#' @param run_kegg Logical; also run KEGG analyses. Defaults to `FALSE`.
#'
#' @return
#' A list containing two named elements:
#' \describe{
#'   \item{`results`}{A list of data frames with the detailed enrichment statistics for each analysis type (e.g., `go_up`, `kegg_down`, `gse_go`).}
#'   \item{`plots`}{A list of `ggplot` objects for visualizing the enrichment results (e.g., dot plots, bar plots).}
#' }
#' The function returns `NULL` if the initial mapping of gene identifiers to
#' Entrez IDs fails.
#'
#' @export
#'
#' @importFrom enrichplot dotplot pairwise_termsim
#' @importFrom ggplot2 facet_grid ggtitle
#' @importFrom clusterProfiler filter
#'
#' @examples
#' # Create a sample DE results dataframe
#' de_results <- data.frame(
#'   Genes = c("TP53", "EGFR", "BRCA1", "TNF", "MMP9", "IL6", "VEGFA", "JUN"),
#'   logFC = c(2.5, -1.8, 2.1, 1.6, -2.2, 3.0, -1.7, 2.2),
#'   adj.P.Val = c(0.001, 0.005, 0.002, 0.009, 0.003, 0.0001, 0.008, 0.001)
#' )
#'
#' # This example requires an internet connection and the org.Hs.eg.db package.
#' \dontrun{
#' if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
#'   enrichment_output <- enrich_pathways(
#'     DE = de_results,
#'     lfc_threshold = 1.5,
#'     fdr_threshold = 0.01,
#'     go_org = org.Hs.eg.db,
#'     kegg_org = 'hsa'
#'   )
#'
#'   # Check if the analysis produced results
#'   if (!is.null(enrichment_output)) {
#'     # View the head of the GO results for upregulated genes
#'     print(head(enrichment_output$results$go_up))
#'
#'     # Display one of the generated plots
#'     if (!is.null(enrichment_output$plots$go_up_dotplot)) {
#'       print(enrichment_output$plots$go_up_dotplot)
#'     }
#'   }
#' }
#' }
enrich_pathways = function(DE, lfc_threshold=1, fdr_threshold=0.01, enrich_pvalue=0.05,
                           go_org = org.Hs.eg.db, kegg_org = 'hsa', gene_col = "Genes", adj = TRUE,
                           source = c("go", "custom"), go_ont = "BP", term2gene = NULL,
                           term2name = NULL, run_ora = TRUE, run_gsea = TRUE, run_kegg = FALSE){
  source <- match.arg(source)
  datas <- list()
  plots <- list()

  stat_col <- if ("rho" %in% names(DE)) "rho" else "logFC"
  if (!stat_col %in% names(DE) || !has_enrichment_signal(DE[[stat_col]])) {
    datas$message <- data.frame(
      message = "Pathway enrichment was skipped because the comparison had no effect-size signal.",
      stringsAsFactors = FALSE
    )
    plots$ora_up_dotplot <- plots$ora_down_dotplot <- plots$gsea_dotplot <- NULL
    return(list(results = datas, plots = plots))
  }
  DT <- ProtPipe::add_entrez(DE, org = go_org, gene_col = gene_col)

  if("rho" %in% names(DT)){
    DT$logFC <- DT$rho
  }

  # Check if any genes were successfully mapped
  if (is.null(DT)) {
    warning("No genes were mapped to Entrez IDs. Skipping enrichment.")
    return(NULL)
  }

  #initialize all plots and dataframes
  datas$ora_up <- datas$ora_down <- datas$gsea <- NULL
  plots$ora_up_dotplot <- plots$ora_down_dotplot <- plots$gsea_dotplot <- NULL


  ## up and down regulated genes
  if(adj){
    up_genes=DT[which(DT$logFC>=lfc_threshold&DT$adj.P.Val<=fdr_threshold),]
    down_genes=DT[which(DT$logFC<=(-lfc_threshold)&DT$adj.P.Val<=fdr_threshold),]
  }else{
    up_genes=DT[which(DT$logFC>=lfc_threshold&DT$P.Value<=fdr_threshold),]
    down_genes=DT[which(DT$logFC<=(-lfc_threshold)&DT$P.Value<=fdr_threshold),]
  }


  # over representation enrichment for upregulated genes ###############################
  if (run_ora && nrow(up_genes) > 0){
    ora_up <- if (source == "go") {
      ProtPipe::enrich_go(up_genes$ENTREZID, DT$ENTREZID, org = go_org, enrich_pvalue = enrich_pvalue, ont = go_ont)
    } else {
      ProtPipe::enrich_terms(up_genes$ENTREZID, DT$ENTREZID, term2gene = term2gene, term2name = term2name, enrich_pvalue = enrich_pvalue)
    }

    if(!is.null(ora_up)){
      datas$ora_up <- ora_up@result
      plots$ora_up_dotplot <- enrichplot::dotplot(ora_up, showCategory = 10)
    }

    if (run_kegg) {
      kegg_up <- ProtPipe::enrich_kegg(up_genes$ENTREZID, DT$ENTREZID, org = go_org, organism = kegg_org, enrich_pvalue = enrich_pvalue)
      if(!is.null(kegg_up)){
        datas$kegg_up <- kegg_up@result
        plots$kegg_up_dotplot <- enrichplot::dotplot(kegg_up, showCategory = 10)
      }
    }
  }
  # over representation enrichment for downregulated genes ###############################
  if (run_ora && nrow(down_genes) > 0){
    ora_down <- if (source == "go") {
      ProtPipe::enrich_go(down_genes$ENTREZID, DT$ENTREZID, org = go_org, enrich_pvalue = enrich_pvalue, ont = go_ont)
    } else {
      ProtPipe::enrich_terms(down_genes$ENTREZID, DT$ENTREZID, term2gene = term2gene, term2name = term2name, enrich_pvalue = enrich_pvalue)
    }

    if(!is.null(ora_down)){
      datas$ora_down <- ora_down@result
      plots$ora_down_dotplot <- enrichplot::dotplot(ora_down, showCategory = 10)
    }

    if (run_kegg) {
      kegg_down <- ProtPipe::enrich_kegg(down_genes$ENTREZID, DT$ENTREZID, org = go_org, organism = kegg_org, enrich_pvalue = enrich_pvalue)
      if(!is.null(kegg_down)){
        datas$kegg_down <- kegg_down@result
        plots$kegg_down_dotplot <- enrichplot::dotplot(kegg_down, showCategory = 10)
      }
    }
  }

  #Gene Set Enrichment ##############
  df_ordered <- DT[order(DT$logFC, decreasing = TRUE), ]
  ordered_genes <- df_ordered$logFC
  names(ordered_genes) <- df_ordered$ENTREZID
  ordered_genes_unique <- ordered_genes[!duplicated(names(ordered_genes))]

  if (run_gsea) {
    if (!has_valid_gsea_ranking(ordered_genes_unique)) {
      datas$gsea_message <- data.frame(
        message = "GSEA was skipped because the ranked statistics had no directional signal.",
        stringsAsFactors = FALSE
      )
    } else {
      gse_result <- if (source == "go") {
        ProtPipe::gse_go(ordered_genes_unique, org = go_org, enrich_pvalue = enrich_pvalue, ont = go_ont)
      } else {
        ProtPipe::gse_terms(ordered_genes_unique, term2gene = term2gene, term2name = term2name, enrich_pvalue = enrich_pvalue)
      }

      if(!is.null(gse_result)){
        datas$gsea <- gse_result@result
        plots$gsea_dotplot <- enrichplot::dotplot(gse_result, showCategory=10, split=".sign") + facet_grid(.~.sign)
      }

      if (run_kegg) {
        gse_kegg <- ProtPipe::gse_kegg(ordered_genes_unique, org = go_org, organism = kegg_org, enrich_pvalue = enrich_pvalue)
        if(!is.null(gse_kegg)){
          datas$gse_kegg <- gse_kegg@result
          plots$gse_kegg_dotplot <- enrichplot::dotplot(gse_kegg, showCategory=10, split=".sign") + facet_grid(.~.sign)
        }
      }
    }
  }
  return(list(results = datas, plots = plots))
}
