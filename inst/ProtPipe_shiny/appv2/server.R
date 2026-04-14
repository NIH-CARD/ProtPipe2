source("shared.R")
source("helpers.R")

server <- function(input, output, session) {
  page_state <- reactiveVal(page_ids[[1]])
  rv <- reactiveValues(
    data = NULL,
    condition = NULL,
    number_samples = NULL,
    type = "Upload a proteomics dataset to begin."
  )

  get_page_index <- function(page_id) {
    match(page_id, page_ids)
  }

  observe({
    updateTextInput(session, "current_page", value = page_state())
  })

  observeEvent(input$back_page, {
    current_index <- get_page_index(page_state())
    if (!is.na(current_index) && current_index > 1) {
      page_state(page_ids[[current_index - 1]])
    }
  })

  observeEvent(input$next_page, {
    current_index <- get_page_index(page_state())
    if (!is.na(current_index) && current_index < length(page_ids)) {
      page_state(page_ids[[current_index + 1]])
    }
  })

  lapply(page_ids, function(id) {
    observeEvent(input[[paste0("tab_", id)]], {
      page_state(id)
    }, ignoreInit = TRUE)
  })

  output$top_tabs <- renderUI({
    tryCatch({
      active_id <- page_state()
      div(
        class = "appv2-tabs",
        lapply(seq_along(page_ids), function(i) {
          top_tab_button(
            id = page_ids[[i]],
            label = workflow_pages[[i]],
            active = identical(page_ids[[i]], active_id)
          )
        })
      )
    }, error = function(e) {
      validate(need(FALSE, paste("Navigation failed:", e$message)))
    })
  })

  # Input ---------------------------------------------------------------------

  intensity_file <- reactive({
    req(rv$data)
    rv$data
  })

  data_type <- reactive({
    rv$type
  })

  observeEvent(input$intensity_matrix, {
    file_info <- input$intensity_matrix

    if (is.null(file_info)) {
      rv$data <- NULL
      rv$condition <- NULL
      rv$number_samples <- NULL
      rv$type <- "Upload a proteomics dataset to begin."
      return()
    }

    tryCatch({
      ext <- tolower(tools::file_ext(file_info$datapath))

      validate(
        need(ext %in% c("csv", "tsv", "xlsx", "xls", "adat"), "Invalid file format.")
      )

      if (ext == "adat") {
        rv$type <- "SomaScan"
        out <- SomaDataIO::read_adat(file_info$datapath) |>
          ProtPipe::soma_all_output()
        rv$data <- out$data
        rv$condition <- out$condition
        rv$number_samples <- out$number_samples
      } else if (detect_olink_npx(file_info$datapath)) {
        rv$type <- "Olink"
        out <- OlinkAnalyze::read_NPX(file_info$datapath) |>
          ProtPipe::olink_all_output()
        rv$data <- out$data
        rv$condition <- out$condition
        rv$number_samples <- out$number_samples
      } else {
        rv$type <- "Standard Matrix"
        if (ext %in% c("csv", "tsv")) {
          rv$data <- data.table::fread(file_info$datapath, data.table = FALSE)
        } else {
          rv$data <- readxl::read_excel(file_info$datapath)
        }
        rv$condition <- NULL
        rv$number_samples <- NULL
      }
    }, error = function(e) {
      rv$data <- NULL
      rv$condition <- NULL
      rv$number_samples <- NULL
      rv$type <- paste("File Processing Error:", e$message)
      showNotification(rv$type, type = "error", duration = 10)
    })
  })

  metadata_df <- reactive({
    file_info <- input$sample_metadata
    if (is.null(file_info)) {
      return(NULL)
    }

    ext <- tolower(tools::file_ext(file_info$datapath))
    validate(
      need(ext %in% c("csv", "tsv", "xlsx", "xls"), "Invalid metadata file format.")
    )

    if (ext %in% c("csv", "tsv")) {
      data.table::fread(file_info$datapath, data.table = FALSE)
    } else {
      readxl::read_excel(file_info$datapath)
    }
  })

  condition_file <- reactive({
    metadata_df()
  })

  output$column_range_ui_v2 <- renderUI({
    tryCatch({
      req(intensity_file())
      df <- ProtPipe::convert_numeric_cols(intensity_file())
      choices <- names(df)

      if (data_type() == "Standard Matrix") {
        intensity_cols <- ProtPipe::detect_intensity_cols(df)
        if (length(intensity_cols) > 0) {
          first <- intensity_cols[[1]]
          last <- intensity_cols[[length(intensity_cols)]]
        } else {
          first <- 1
          last <- length(choices)
        }
      } else {
        num_cols <- length(choices)
        sample_count <- if (is.null(rv$number_samples)) 0 else rv$number_samples
        first <- max(1, num_cols - sample_count + 1)
        last <- num_cols
      }

      tagList(
        selectInput("lower_col_v2", "Intensity columns start at:", choices = choices, selected = choices[[first]]),
        selectInput("upper_col_v2", "Intensity columns end at:", choices = choices, selected = choices[[last]])
      )
    }, error = function(e) {
      validate(need(FALSE, paste("Preparing intensity column selection failed:", e$message)))
    })
  })

  raw_prot_data_result <- reactive({
    tryCatch({
      req(intensity_file(), input$lower_col_v2, input$upper_col_v2)
      df <- intensity_file()
      cols <- names(df)
      lower_idx <- match(input$lower_col_v2, cols)
      upper_idx <- match(input$upper_col_v2, cols)

      validate(
        need(!is.na(lower_idx) && !is.na(upper_idx), "Selected columns were not found."),
        need(lower_idx <= upper_idx, "The first intensity column must come before the last intensity column.")
      )

      se_obj <- if (data_type() == "Standard Matrix") {
        ProtPipe::create_se(
          dat = intensity_file(),
          intensity_cols = c(lower_idx:upper_idx),
          sample_metadata = condition_file()
        )
      } else {
        if (is.null(condition_file())) {
          condition <- rv$condition
        } else if (is.null(rv$condition)) {
          condition <- condition_file()
        } else {
          condition <- dplyr::left_join(rv$condition, condition_file(), by = "SampleID")
        }
        ProtPipe::create_se(
          dat = intensity_file(),
          intensity_cols = c(lower_idx:upper_idx),
          sample_metadata = condition
        )
      }

      list(data = se_obj, error = NULL)
    }, error = function(e) {
      list(data = NULL, error = e$message)
    })
  })

  raw_prot_data <- reactive(raw_prot_data_result()$data)
  raw_prot_data_error <- reactive(raw_prot_data_result()$error)

  validate_raw_prot_data <- function() {
    validate(need(is.null(raw_prot_data_error()), raw_prot_data_error()))
  }

  output$file_summary_v2 <- renderText({
    tryCatch({
      se_obj <- raw_prot_data()
      error_msg <- raw_prot_data_error()
      if (is.null(se_obj)) {
        paste(
          "Data type:", rv$type,
          "\nNumber of samples: 0",
          "\nNumber of proteins: 0",
          if (!is.null(error_msg)) paste0("\n", error_msg) else ""
        )
      } else {
        paste(
          "Data type:", rv$type,
          "\nNumber of samples:", ncol(se_obj),
          "\nNumber of proteins:", nrow(se_obj)
        )
      }
    }, error = function(e) {
      validate(need(FALSE, paste("Creating data summary failed:", e$message)))
    })
  })

  output$input_preview_table <- renderTable({
    tryCatch({
      validate_raw_prot_data()
      se_obj <- raw_prot_data()
      req(se_obj)

      switch(
        input$input_main_view_v2,
        "Protein Intensity" = utils::head(as.data.frame(SummarizedExperiment::assay(se_obj)), 12),
        "Sample Metadata" = {
          col_df <- as.data.frame(SummarizedExperiment::colData(se_obj))
          if (ncol(col_df) == 0) data.frame() else utils::head(col_df, 12)
        },
        "Protein Metadata" = {
          row_df <- as.data.frame(SummarizedExperiment::rowData(se_obj))
          if (ncol(row_df) == 0) data.frame() else utils::head(row_df, 12)
        }
      )
    }, error = function(e) {
      validate(need(FALSE, paste("Rendering input preview failed:", e$message)))
    })
  }, striped = TRUE, bordered = TRUE, spacing = "s")

  output$download_example_proteomics_v2 <- downloadHandler(
    filename = function() {
      "example_proteomics_dataset.csv"
    },
    content = function(file) {
      file.copy(file.path("www", "iPSC.csv"), file, overwrite = TRUE)
    }
  )

  selected_main_view <- reactive({
    current_id <- page_state()
    if (identical(current_id, "input")) {
      if (is.null(input$input_main_view_v2)) "Assay" else input$input_main_view_v2
    } else {
      view_value <- input[[paste0(current_id, "_main_view_v2")]]
      if (is.null(view_value)) "Plot" else view_value
    }
  })

  # Pre-Processing ------------------------------------------------------------

  output$lod_filtering_v2 <- renderUI({
    tryCatch({
      req(intensity_file())
      req(data_type() == "Olink" || data_type() == "SomaScan")

      label <- if (data_type() == "Olink") {
        "Filter values based on LOD values (if present)"
      } else {
        "Filter values based on average buffer values (if present)"
      }

      checkboxInput("lod_filter_v2", label = label, value = FALSE)
    }, error = function(e) {
      validate(need(FALSE, paste("Preparing LOD filtering options failed:", e$message)))
    })
  })

  output$preprocessing_min_filter_ui_v2 <- renderUI({
    tryCatch({
      left_block <- div(
        checkboxInput("min_int_filter_v2", "Apply minimum intensity filter", value = FALSE),
        numericInput("min_int_filter_lod_v2", "Minimum intensity", value = 0)
      )

      if (data_type() %in% c("Olink", "SomaScan")) {
        div(
          class = "appv2-subsection-grid",
          left_block,
          div(uiOutput("lod_filtering_v2"))
        )
      } else {
        left_block
      }
    }, error = function(e) {
      validate(need(FALSE, paste("Preparing minimum intensity filtering options failed:", e$message)))
    })
  })

  output$batch_correct_column_v2 <- renderUI({
    tryCatch({
      validate_raw_prot_data()
      req(raw_prot_data())
      if (is.null(condition_file())) {
        tags$p(class = "appv2-subtitle", "Must upload sample condition file.")
      } else {
        choices <- names(SummarizedExperiment::colData(raw_prot_data()))
        selectInput("batch_correct_column_v2", "Select condition for correction", choices = choices)
      }
    }, error = function(e) {
      validate(need(FALSE, paste("Preparing batch correction options failed:", e$message)))
    })
  })

  output$batch_correct_section_v2 <- renderUI({
    tryCatch({
      validate_raw_prot_data()
      req(raw_prot_data())
      if (is.null(condition_file())) {
        tags$p(class = "appv2-subtitle", "Must upload sample condition file.")
      } else {
        div(
          class = "appv2-subsection-grid",
          div(checkboxInput("batch_correct_v2", "Batch correct", value = FALSE)),
          div(uiOutput("batch_correct_column_v2"))
        )
      }
    }, error = function(e) {
      validate(need(FALSE, paste("Preparing batch correction section failed:", e$message)))
    })
  })

  prot_data_result <- reactive({
    tryCatch({
      validate_raw_prot_data()
      req(raw_prot_data())
      PD <- raw_prot_data()

      if (isTRUE(input$lod_filter_v2)) {
        lod_col <- NULL
        if (data_type() == "Olink") {
          lod_col <- "LOD"
        }
        if (data_type() == "SomaScan") {
          lod_col <- "Buffer"
        }
        if (!is.null(lod_col)) {
          PD <- ProtPipe::lod_filter(PD, lod_col)
        }
      }

      if (isTRUE(input$min_int_filter_v2)) {
        PD <- ProtPipe::apply_min_intenisty(PD, input$min_int_filter_lod_v2)
      }

      if (isTRUE(input$remove_outliers_v2)) {
        PD <- ProtPipe::filter_outlier_samples(PD, sds = input$outlier_sds_v2)
      }
      if (isTRUE(input$remove_sparse_proteins_v2)) {
        PD <- ProtPipe::filter_proteins_by_percent(PD, percent = input$sparse_protein_percent_v2)
      }

      if (isTRUE(input$log2_transform_v2)) {
        tryCatch({
          PD <- ProtPipe::log2_transform(PD)
        }, error = function(e) {
          message("Transformation failed: ", e$message)
        })
      }

      if (isTRUE(input$normalize_v2)) {
        tryCatch({
          if (input$normalize_method_v2 == "mean") {
            PD <- ProtPipe::mean_normalize(PD)
          } else if (input$normalize_method_v2 == "median") {
            PD <- ProtPipe::median_normalize(PD)
          }
        }, error = function(e) {
          message("Normalization failed: ", e$message)
        })
      }

      if (isTRUE(input$impute_v2)) {
        if (input$imputation_method_v2 == "fixed value") {
          PD <- ProtPipe::impute(PD, input$impute_fixed_value_v2)
        } else if (input$imputation_method_v2 == "minimum") {
          PD <- ProtPipe::impute_min(PD, input$impute_min_value_v2)
        } else if (input$imputation_method_v2 == "left-shifted distribution") {
          PD <- ProtPipe::impute_left_dist(PD, input$impute_left_dist_shift_v2, input$impute_left_dist_scale_v2)
        }
      }

      if (!is.null(input$batch_correct_column_v2) && isTRUE(input$batch_correct_v2)) {
        PD <- ProtPipe::batch_correct(PD, input$batch_correct_column_v2)
      }

      list(data = PD, error = NULL)
    }, error = function(e) {
      list(data = NULL, error = e$message)
    })
  })

  prot_data <- reactive(prot_data_result()$data)
  prot_data_error <- reactive(prot_data_result()$error)

  validate_prot_data <- function() {
    validate(need(is.null(prot_data_error()), prot_data_error()))
  }

  # Quality Control -----------------------------------------------------------

  output$quality_control_params_v2 <- renderUI({
    tryCatch({
      validate_raw_prot_data()
      req(raw_prot_data())

      qc_choices <- names(SummarizedExperiment::colData(raw_prot_data()))
      selected_view <- if (is.null(input$quality_control_main_view_v2)) "Protein Groups" else input$quality_control_main_view_v2

      switch(
        selected_view,
        "Protein Groups" = div(
          class = "appv2-params-stack",
          selectInput("qc_pg_condition_v2", "Group samples by", choices = qc_choices)
        ),
        "Protein Intensities" = div(
          class = "appv2-params-stack",
          tags$p(class = "appv2-subtitle", "No additional parameters for this view.")
        ),
        "Sample Correlations" = div(
          class = "appv2-params-stack",
          checkboxInput("qc_use_all_features_v2", "Use all proteins for correlation (slower)", value = FALSE),
          conditionalPanel(
            condition = "!input.qc_use_all_features_v2",
            numericInput("qc_num_features_v2", "Number of proteins for correlation (n)", value = 1000, min = 1, step = 100)
          ),
          tags$p(
            class = "appv2-subtitle",
            "The top n most variable proteins in the dataset will be used to calculate spearman correlation. This will decrease analysis time, but may result in lower correlation values compared to using the full set of proteins."
          )
        ),
        "Coefficient of Variation" = div(
          class = "appv2-params-stack",
          selectInput("qc_condition_v2", "Group samples by", choices = qc_choices),
          selectInput("cv_plot_type_v2", "Display style", choices = c("violin", "jitter"), selected = "violin")
        )
      )
    }, error = function(e) {
      validate(need(FALSE, paste("Preparing quality control parameters failed:", e$message)))
    })
  })

  qc_cvs_reactive <- reactive({
    req(raw_prot_data(), input$qc_condition_v2)
    list(
      data = ProtPipe::get_CVs(raw_prot_data(), condition = input$qc_condition_v2),
      plot = ProtPipe::plot_CVs(
        raw_prot_data(),
        condition = input$qc_condition_v2,
        plot_type = input$cv_plot_type_v2
      )
    )
  })

  qc_protein_groups_reactive <- reactive({
    req(raw_prot_data(), input$qc_pg_condition_v2)
    list(
      data = get_pg_counts(prot_data()),
      plot = ProtPipe::plot_pg_counts(raw_prot_data(), condition = input$qc_pg_condition_v2)
    )
  })

  qc_intensities_reactive <- reactive({
    req(raw_prot_data())
    list(
      data = as.data.frame(SummarizedExperiment::assay(raw_prot_data())),
      plot = ProtPipe::plot_pg_intensities(raw_prot_data())
    )
  })

  qc_correlations_reactive <- reactive({
    req(raw_prot_data())
    num_features <- if (isTRUE(input$qc_use_all_features_v2)) NULL else input$qc_num_features_v2
    if (!isTRUE(input$qc_use_all_features_v2)) {
      req(num_features)
    }
    list(
      data = ProtPipe::get_sample_correlation(raw_prot_data(), num_features = num_features),
      plot = ProtPipe::plot_correlation_heatmap(raw_prot_data(), num_features = num_features)
    )
  })

  qc_plot_dimensions_v2 <- reactive({
    list(
      width_in = if (is.null(input$qc_plot_width_v2)) 8 else input$qc_plot_width_v2,
      height_in = if (is.null(input$qc_plot_height_v2)) 6 else input$qc_plot_height_v2
    )
  })

  output$quality_control_plot_ui_v2 <- renderUI({
    dims <- qc_plot_dimensions_v2()
    width_px <- round(dims$width_in * 96)
    height_px <- round(dims$height_in * 96)

    div(
      style = paste0("width:", width_px, "px;"),
      plotOutput("quality_control_plot_v2", height = paste0(height_px, "px"))
    )
  })

  output$quality_control_plot_v2 <- renderPlot({
    tryCatch({
      validate_raw_prot_data()
      req(raw_prot_data())
      selected_view <- if (is.null(input$quality_control_main_view_v2)) "Protein Groups" else input$quality_control_main_view_v2

      p <- switch(
        selected_view,
        "Protein Groups" = qc_protein_groups_reactive()$plot,
        "Protein Intensities" = qc_intensities_reactive()$plot,
        "Sample Correlations" = qc_correlations_reactive()$plot,
        "Coefficient of Variation" = qc_cvs_reactive()$plot
      )

      print(p)
    }, error = function(e) {
      validate(need(FALSE, paste("Rendering quality control plot failed:", e$message)))
    })
  })

  output$download_qc_plot_v2 <- downloadHandler(
    filename = function() {
      selected_view <- if (is.null(input$quality_control_main_view_v2)) "Protein Groups" else input$quality_control_main_view_v2
      paste0(gsub(" ", "_", tolower(selected_view)), ".pdf")
    },
    content = function(file) {
      selected_view <- if (is.null(input$quality_control_main_view_v2)) "Protein Groups" else input$quality_control_main_view_v2
      dims <- qc_plot_dimensions_v2()
      p <- switch(
        selected_view,
        "Protein Groups" = qc_protein_groups_reactive()$plot,
        "Protein Intensities" = qc_intensities_reactive()$plot,
        "Sample Correlations" = qc_correlations_reactive()$plot,
        "Coefficient of Variation" = qc_cvs_reactive()$plot
      )
      ggsave(file, plot = p, device = "pdf", width = dims$width_in, height = dims$height_in, units = "in")
    }
  )

  output$download_qc_table_v2 <- downloadHandler(
    filename = function() {
      selected_view <- if (is.null(input$quality_control_main_view_v2)) "Protein Groups" else input$quality_control_main_view_v2
      paste0(gsub(" ", "_", tolower(selected_view)), ".tsv")
    },
    content = function(file) {
      selected_view <- if (is.null(input$quality_control_main_view_v2)) "Protein Groups" else input$quality_control_main_view_v2
      dat <- switch(
        selected_view,
        "Protein Groups" = qc_protein_groups_reactive()$data,
        "Protein Intensities" = qc_intensities_reactive()$data,
        "Sample Correlations" = qc_correlations_reactive()$data,
        "Coefficient of Variation" = qc_cvs_reactive()$data
      )
      write.table(dat, file = file, sep = "\t", quote = FALSE, row.names = FALSE)
    }
  )

  # Pre-Processing ------------------------------------------------------------

  output$download_preprocessing_data_v2 <- downloadHandler(
    filename = function() {
      "processed_data.tsv"
    },
    content = function(file) {
      data <- base::cbind(SummarizedExperiment::rowData(prot_data()), SummarizedExperiment::assay(prot_data()))
      write.table(data, file = file, sep = "\t", quote = FALSE, row.names = FALSE)
    }
  )

  output$download_preprocessing_report_v2 <- downloadHandler(
    filename = function() {
      "processing_report.md"
    },
    content = function(file) {
      req(prot_data())
      ProtPipe::generate_preprocessing_report(
        object = prot_data(),
        output_file = file
      )
    }
  )

  # Clustering / Dimensionality Reduction ------------------------------------

  output$clustering_params_v2 <- renderUI({
    tryCatch({
      validate_prot_data()
      req(prot_data())
      selected_view <- if (is.null(input$clustering_main_view_v2)) "Hierarchical Clustering" else input$clustering_main_view_v2
      choices <- names(SummarizedExperiment::colData(prot_data()))

      switch(
        selected_view,
        "Hierarchical Clustering" = div(
          class = "appv2-params-stack",
          tags$p(class = "appv2-subtitle", "No additional parameters for this view.")
        ),
        "PCA" = div(
          class = "appv2-params-stack",
          selectInput("cluster_condition_v2", "Group samples by", choices = choices)
        ),
        "UMAP" = div(
          class = "appv2-params-stack",
          selectInput("cluster_condition_umap_v2", "Group samples by", choices = choices),
          sliderInput("neighbors_v2", "Number of neighbors", min = 2, max = max(2, ncol(prot_data())), value = min(15, max(2, floor((ncol(prot_data()) + 2) / 2))), step = 1)
        )
      )
    }, error = function(e) {
      validate(need(FALSE, paste("Preparing clustering parameters failed:", e$message)))
    })
  })

  clustering_hcluster_reactive <- reactive({
    req(prot_data())
    ProtPipe::plot_hierarchical_cluster(prot_data())
  })

  clustering_pca_data_reactive <- reactive({
    req(prot_data(), input$cluster_condition_v2)
    get_PCs(prot_data(), condition = input$cluster_condition_v2)
  })

  clustering_pca_plot_reactive <- reactive({
    req(prot_data(), input$cluster_condition_v2)
    ProtPipe::plot_PCs(prot_data(), condition = input$cluster_condition_v2)
  })

  clustering_umap_results <- reactive({
    req(prot_data(), input$neighbors_v2, input$cluster_condition_umap_v2)
    set.seed(123)
    list(
      table = ProtPipe::get_umap(
        prot_data(),
        neighbors = input$neighbors_v2,
        condition = input$cluster_condition_umap_v2
      ),
      plot = ProtPipe::plot_umap(
        prot_data(),
        neighbors = input$neighbors_v2,
        condition = input$cluster_condition_umap_v2
      )
    )
  })

  clustering_plot_dimensions_v2 <- reactive({
    list(
      width_in = if (is.null(input$clustering_plot_width_v2)) 8 else input$clustering_plot_width_v2,
      height_in = if (is.null(input$clustering_plot_height_v2)) 6 else input$clustering_plot_height_v2
    )
  })

  output$clustering_plot_ui_v2 <- renderUI({
    dims <- clustering_plot_dimensions_v2()
    width_px <- round(dims$width_in * 96)
    height_px <- round(dims$height_in * 96)

    div(
      style = paste0("width:", width_px, "px;"),
      plotOutput("clustering_plot_v2", height = paste0(height_px, "px"))
    )
  })

  output$clustering_plot_v2 <- renderPlot({
    tryCatch({
      validate_prot_data()
      req(prot_data())
      selected_view <- if (is.null(input$clustering_main_view_v2)) "Hierarchical Clustering" else input$clustering_main_view_v2
      p <- switch(
        selected_view,
        "Hierarchical Clustering" = clustering_hcluster_reactive(),
        "PCA" = clustering_pca_plot_reactive(),
        "UMAP" = clustering_umap_results()$plot
      )
      print(p)
    }, error = function(e) {
      validate(need(FALSE, paste("Rendering clustering plot failed:", e$message)))
    })
  })

  output$download_clustering_plot_v2 <- downloadHandler(
    filename = function() {
      selected_view <- if (is.null(input$clustering_main_view_v2)) "Hierarchical Clustering" else input$clustering_main_view_v2
      paste0(gsub(" ", "_", tolower(selected_view)), ".pdf")
    },
    content = function(file) {
      dims <- clustering_plot_dimensions_v2()
      selected_view <- if (is.null(input$clustering_main_view_v2)) "Hierarchical Clustering" else input$clustering_main_view_v2
      p <- switch(
        selected_view,
        "Hierarchical Clustering" = clustering_hcluster_reactive(),
        "PCA" = clustering_pca_plot_reactive(),
        "UMAP" = clustering_umap_results()$plot
      )
      ggsave(file, plot = p, device = "pdf", width = dims$width_in, height = dims$height_in, units = "in")
    }
  )

  output$download_clustering_table_v2 <- downloadHandler(
    filename = function() {
      selected_view <- if (is.null(input$clustering_main_view_v2)) "Hierarchical Clustering" else input$clustering_main_view_v2
      paste0(gsub(" ", "_", tolower(selected_view)), ".tsv")
    },
    content = function(file) {
      selected_view <- if (is.null(input$clustering_main_view_v2)) "Hierarchical Clustering" else input$clustering_main_view_v2
      dat <- switch(
        selected_view,
        "Hierarchical Clustering" = as.data.frame(SummarizedExperiment::assay(prot_data())),
        "PCA" = clustering_pca_data_reactive()$components,
        "UMAP" = clustering_umap_results()$table
      )
      write.table(dat, file = file, sep = "\t", quote = FALSE, row.names = FALSE)
    }
  )

  output$workflow_preview_table <- renderTable({
    tryCatch({
      validate_prot_data()
      se_obj <- prot_data()
      req(se_obj)

      switch(
        selected_main_view(),
        "Assay" = utils::head(as.data.frame(SummarizedExperiment::assay(se_obj)), 12),
        "Sample Metadata" = {
          col_df <- as.data.frame(SummarizedExperiment::colData(se_obj))
          if (ncol(col_df) == 0) data.frame() else utils::head(col_df, 12)
        },
        "Protein Metadata" = {
          row_df <- as.data.frame(SummarizedExperiment::rowData(se_obj))
          if (ncol(row_df) == 0) data.frame() else utils::head(row_df, 12)
        }
      )
    }, error = function(e) {
      validate(need(FALSE, paste("Rendering preview table failed:", e$message)))
    })
  }, striped = TRUE, bordered = TRUE, spacing = "s")

  output$workflow_main_body <- renderUI({
    tryCatch({
      validate_prot_data()
      if (identical(selected_main_view(), "Plot")) {
        current_index <- get_page_index(page_state())
        return(
          div(
            class = "appv2-panel appv2-plot-card",
            h2(workflow_pages[[current_index]]),
            tags$p("Plot placeholder")
          )
        )
      }

      div(
        class = "appv2-panel appv2-plot",
        div(class = "appv2-table-wrap", tableOutput("workflow_preview_table"))
      )
    }, error = function(e) {
      validate(need(FALSE, paste("Preparing main view failed:", e$message)))
    })
  })

  # Differential Intensity ----------------------------------------------------

  differential_selected_view_v2 <- reactive({
    if (is.null(input$differential_main_view_v2)) "Volcano Plot" else input$differential_main_view_v2
  })

  differential_plot_dimensions_v2 <- reactive({
    list(
      width_in = if (is.null(input$differential_plot_width_v2)) 6 else input$differential_plot_width_v2,
      height_in = if (is.null(input$differential_plot_height_v2)) 4 else input$differential_plot_height_v2
    )
  })

  observe({
    if (!is.null(prot_data_error()) || is.null(prot_data())) {
      return()
    }
    col_choices <- names(SummarizedExperiment::colData(prot_data()))
    row_choices <- names(SummarizedExperiment::rowData(prot_data()))
    updateSelectInput(session, "de_condition_v2", choices = col_choices)
    updateSelectInput(session, "label_col_v2", choices = row_choices)
    updateSelectInput(session, "gene_col_v2", choices = row_choices)
  })

  output$de_groups_v2 <- renderUI({
    tryCatch({
      validate_prot_data()
      req(prot_data(), input$de_condition_v2, input$de_mode_v2 == "binary")
      groups <- unique(SummarizedExperiment::colData(prot_data())[[input$de_condition_v2]])
      tagList(
        selectInput("control_condition_v2", "Control group", choices = groups),
        selectInput("treatment_condition_v2", "Treatment group", choices = groups)
      )
    }, error = function(e) {
      validate(need(FALSE, paste("Preparing group selectors failed:", e$message)))
    })
  })

  output$de_covariates_v2 <- renderUI({
    tryCatch({
      validate_prot_data()
      req(prot_data(), input$de_mode_v2 == "binary")
      choices <- names(SummarizedExperiment::colData(prot_data()))
      selectInput("de_covariates_v2", "Covariates", choices = choices, multiple = TRUE, selected = NULL)
    }, error = function(e) {
      validate(need(FALSE, paste("Preparing covariate selector failed:", e$message)))
    })
  })

  output$logfc_v2 <- renderUI({
    tryCatch({
      if (input$de_mode_v2 == "continuous") {
        numericInput("logfc_v2", "Spearman rho cutoff", value = 0.35)
      } else {
        numericInput("logfc_v2", "Log2 fold-change cutoff", value = 1)
      }
    }, error = function(e) {
      validate(need(FALSE, paste("Preparing effect size input failed:", e$message)))
    })
  })

  output$pvalue_v2 <- renderUI({
    tryCatch({
      label <- if (isTRUE(input$use_adj_pval_v2)) "Adjusted P-value cutoff" else "P-value cutoff"
      numericInput("pvalue_v2", label = label, value = 0.01, min = 0)
    }, error = function(e) {
      validate(need(FALSE, paste("Preparing p-value input failed:", e$message)))
    })
  })

  gene_labels_v2 <- reactive({
    file_info <- input$gene_labels_v2
    if (is.null(file_info)) {
      return(NULL)
    }
    dat <- data.table::fread(file_info$datapath, data.table = FALSE)
    if (!"Genes" %in% names(dat)) {
      validate(need(FALSE, "The uploaded gene label file must contain a column called Genes."))
    }
    dat$Genes
  })

  dea_result_v2 <- reactiveVal(NULL)

  dea_v2 <- reactive({
    req(prot_data(), input$de_condition_v2)
    result <- tryCatch({
      if (input$de_mode_v2 == "continuous") {
        ProtPipe::do_comparison_continuous(prot_data(), condition = input$de_condition_v2)
      } else {
        ProtPipe::do_limma_binary(
          prot_data(),
          condition = input$de_condition_v2,
          control_group = input$control_condition_v2,
          treatment_group = input$treatment_condition_v2,
          covariates = input$de_covariates_v2
        )
      }
    }, error = function(e) {
      dea_result_v2(NULL)
      validate(need(FALSE, paste("Calculating differential intensity failed:", e$message)))
    })
    dea_result_v2(result)
    result
  })

  volcano_plot_v2 <- reactive({
    req(dea_v2(), input$label_col_v2, input$pvalue_v2, input$logfc_v2)
    tryCatch({
      if (input$de_mode_v2 == "continuous") {
        ProtPipe::plot_correlation_volcano(
          dea_v2(),
          label_col = input$label_col_v2,
          labelgene = gene_labels_v2(),
          fdr_threshold = input$pvalue_v2,
          rho_threshold = input$logfc_v2,
          adj = input$use_adj_pval_v2
        )
      } else {
        ProtPipe::plot_volcano(
          dea_v2(),
          label_col = input$label_col_v2,
          labelgene = gene_labels_v2(),
          fdr_threshold = input$pvalue_v2,
          lfc_threshold = input$logfc_v2,
          adj = input$use_adj_pval_v2
        )
      }
    }, error = function(e) {
      validate(need(FALSE, paste("Plotting volcano failed:", e$message)))
    })
  })

  selected_org_v2 <- reactive({
    req(input$organism_v2)
    organism_map[[input$organism_v2]]
  })

  selected_ontology_v2 <- reactive({
    req(input$pathway_source_v2)
    if (identical(input$pathway_source_v2, "GO")) {
      return(NULL)
    }
    req(input$ontology_file_v2)
    ProtPipe::read_ontology(input$ontology_file_v2$datapath)
  })

  enrichment_result_v2 <- reactiveVal(NULL)
  enrichment_message_v2 <- reactiveVal(NULL)

  observeEvent(input$run_enrichment_v2, {
    showNotification("Running enrichment analysis, please wait...", duration = 8)
    tryCatch({
      req(dea_result_v2(), input$gene_col_v2, input$enrich_pval_v2)
      ontology <- selected_ontology_v2()
      result <- ProtPipe::enrich_pathways(
        dea_result_v2(),
        lfc_threshold = input$logfc_v2,
        fdr_threshold = input$pvalue_v2,
        enrich_pvalue = input$enrich_pval_v2,
        go_org = selected_org_v2()$OrgDb,
        kegg_org = selected_org_v2()$kegg,
        gene_col = input$gene_col_v2,
        adj = input$use_adj_pval_v2,
        source = if (identical(input$pathway_source_v2, "GO")) "go" else "custom",
        go_ont = input$go_ontology_v2,
        term2gene = if (is.null(ontology)) NULL else ontology$term2gene,
        term2name = if (is.null(ontology)) NULL else ontology$term2name,
        run_ora = TRUE,
        run_gsea = TRUE
      )
      enrichment_result_v2(result)
      if (is.null(result)) {
        enrichment_message_v2("No genes were mapped to Entrez IDs.")
      } else if ("message" %in% names(result$results)) {
        enrichment_message_v2(result$results$message$message[[1]])
      } else {
        enrichment_message_v2(NULL)
      }
    }, error = function(e) {
      enrichment_result_v2(NULL)
      enrichment_message_v2(paste("Pathway enrichment failed:", e$message))
    })
  })

  output$differential_main_body_v2 <- renderUI({
    tryCatch({
      if (identical(differential_selected_view_v2(), "Volcano Plot")) {
        dims <- differential_plot_dimensions_v2()
        width_px <- round(dims$width_in * 96)
        height_px <- round(dims$height_in * 96)
        return(
          div(
            class = "appv2-panel appv2-plot",
            div(
              style = paste0("width:", width_px, "px;"),
              plotOutput("volcano_plot_v2", height = paste0(height_px, "px"))
            ),
            div(
              class = "appv2-download-row",
              downloadButton("download_volcano_v2", "Download PDF"),
              downloadButton("download_de_table_v2", "Download TSV")
            )
          )
        )
      }

      div(
        class = "appv2-panel appv2-plot",
        actionButton("run_enrichment_v2", "Run pathway enrichment", class = "appv2-actionbtn"),
        verbatimTextOutput("enrichment_message_v2", placeholder = TRUE),
        tags$p(class = "appv2-subsection-title", "GO Up"),
        uiOutput("ora_up_enrich_ui_v2"),
        div(
          class = "appv2-download-row",
          downloadButton("download_ora_up_plot_v2", "GO Up PDF"),
          downloadButton("download_ora_up_table_v2", "GO Up TSV")
        ),
        tags$p(class = "appv2-subsection-title", "GO Down"),
        uiOutput("ora_down_enrich_ui_v2"),
        div(
          class = "appv2-download-row",
          downloadButton("download_ora_down_plot_v2", "GO Down PDF"),
          downloadButton("download_ora_down_table_v2", "GO Down TSV")
        ),
        tags$p(class = "appv2-subsection-title", "GSEA"),
        uiOutput("gsea_enrich_ui_v2"),
        div(
          class = "appv2-download-row",
          downloadButton("download_gsea_plot_v2", "GSEA PDF"),
          downloadButton("download_gsea_table_v2", "GSEA TSV")
        ),
      )
    }, error = function(e) {
      validate(need(FALSE, paste("Preparing differential intensity view failed:", e$message)))
    })
  })

  output$volcano_plot_v2 <- renderPlot({
    tryCatch({
      print(volcano_plot_v2())
    }, error = function(e) {
      validate(need(FALSE, paste("Rendering volcano plot failed:", e$message)))
    })
  })

  output$enrichment_message_v2 <- renderText({
    tryCatch({
      enrichment_message_v2()
    }, error = function(e) {
      validate(need(FALSE, paste("Preparing enrichment message failed:", e$message)))
    })
  })

  output$ora_up_enrich_ui_v2 <- renderUI({
    dims <- differential_plot_dimensions_v2()
    div(
      style = paste0("width:", round(dims$width_in * 96), "px;"),
      plotOutput("ora_up_enrich_v2", height = paste0(round(dims$height_in * 96), "px"))
    )
  })

  output$ora_down_enrich_ui_v2 <- renderUI({
    dims <- differential_plot_dimensions_v2()
    div(
      style = paste0("width:", round(dims$width_in * 96), "px;"),
      plotOutput("ora_down_enrich_v2", height = paste0(round(dims$height_in * 96), "px"))
    )
  })

  output$gsea_enrich_ui_v2 <- renderUI({
    dims <- differential_plot_dimensions_v2()
    div(
      style = paste0("width:", round(dims$width_in * 96), "px;"),
      plotOutput("gsea_enrich_v2", height = paste0(round(dims$height_in * 96), "px"))
    )
  })

  output$ora_up_enrich_v2 <- renderPlot({
    tryCatch({
      req(enrichment_result_v2())
      validate(need(is.null(enrichment_message_v2()), enrichment_message_v2()))
      validate(need(!is.null(enrichment_result_v2()$plots$ora_up_dotplot), "No GO up terms were enriched."))
      enrichment_result_v2()$plots$ora_up_dotplot
    }, error = function(e) {
      validate(need(FALSE, paste("Rendering GO up plot failed:", e$message)))
    })
  })

  output$ora_down_enrich_v2 <- renderPlot({
    tryCatch({
      req(enrichment_result_v2())
      validate(need(is.null(enrichment_message_v2()), enrichment_message_v2()))
      validate(need(!is.null(enrichment_result_v2()$plots$ora_down_dotplot), "No GO down terms were enriched."))
      enrichment_result_v2()$plots$ora_down_dotplot
    }, error = function(e) {
      validate(need(FALSE, paste("Rendering GO down plot failed:", e$message)))
    })
  })

  output$gsea_enrich_v2 <- renderPlot({
    tryCatch({
      req(enrichment_result_v2())
      validate(need(is.null(enrichment_message_v2()), enrichment_message_v2()))
      validate(
        need(
          !("gsea_message" %in% names(enrichment_result_v2()$results)),
          enrichment_result_v2()$results$gsea_message$message[[1]]
        )
      )
      validate(need(!is.null(enrichment_result_v2()$plots$gsea_dotplot), "No GSEA terms were enriched."))
      enrichment_result_v2()$plots$gsea_dotplot
    }, error = function(e) {
      validate(need(FALSE, paste("Rendering GSEA plot failed:", e$message)))
    })
  })

  output$download_volcano_v2 <- downloadHandler(
    filename = function() { "volcano.pdf" },
    content = function(file) {
      dims <- differential_plot_dimensions_v2()
      ggsave(file, plot = volcano_plot_v2(), device = "pdf", width = dims$width_in, height = dims$height_in, units = "in")
    }
  )

  output$download_de_table_v2 <- downloadHandler(
    filename = function() { "differential_expression_results.tsv" },
    content = function(file) {
      write.table(dea_v2(), file = file, sep = "\t", row.names = FALSE, quote = FALSE)
    }
  )

  output$download_enrichment_plots_v2 <- downloadHandler(
    filename = function() { "pathway_enrichment_plots.zip" },
    content = function(file) {
      req(enrichment_result_v2())
      zip_dir <- tempfile("enrichment_plots_")
      dir.create(zip_dir, recursive = TRUE)
      plots <- enrichment_result_v2()$plots

      for (name in names(plots)) {
        p <- plots[[name]]
        if (!is.null(p)) {
          ggplot2::ggsave(
            filename = file.path(zip_dir, paste0(name, ".pdf")),
            plot = p,
            device = "pdf",
            width = 6,
            height = 4,
            units = "in"
          )
        }
      }

      old_wd <- getwd()
      on.exit(setwd(old_wd), add = TRUE)
      setwd(zip_dir)
      utils::zip(zipfile = file, files = list.files())
    }
  )

  output$download_enrichment_tables_v2 <- downloadHandler(
    filename = function() { "pathway_enrichment_tables.zip" },
    content = function(file) {
      req(enrichment_result_v2())
      zip_dir <- tempfile("enrichment_tables_")
      dir.create(zip_dir, recursive = TRUE)
      results <- enrichment_result_v2()$results

      for (name in names(results)) {
        dat <- results[[name]]
        if (is.data.frame(dat) && nrow(dat) > 0) {
          utils::write.table(
            dat,
            file = file.path(zip_dir, paste0(name, ".tsv")),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE
          )
        }
      }

      old_wd <- getwd()
      on.exit(setwd(old_wd), add = TRUE)
      setwd(zip_dir)
      utils::zip(zipfile = file, files = list.files())
    }
  )

  output$download_ora_up_plot_v2 <- downloadHandler(
    filename = function() { "ora_up_dotplot.pdf" },
    content = function(file) {
      req(enrichment_result_v2())
      dims <- differential_plot_dimensions_v2()
      ggplot2::ggsave(file, plot = enrichment_result_v2()$plots$ora_up_dotplot, device = "pdf", width = dims$width_in, height = dims$height_in, units = "in")
    }
  )

  output$download_ora_up_table_v2 <- downloadHandler(
    filename = function() { "ora_up.tsv" },
    content = function(file) {
      req(enrichment_result_v2())
      dat <- enrichment_result_v2()$results$ora_up
      utils::write.table(dat, file = file, sep = "\t", quote = FALSE, row.names = FALSE)
    }
  )

  output$download_ora_down_plot_v2 <- downloadHandler(
    filename = function() { "ora_down_dotplot.pdf" },
    content = function(file) {
      req(enrichment_result_v2())
      dims <- differential_plot_dimensions_v2()
      ggplot2::ggsave(file, plot = enrichment_result_v2()$plots$ora_down_dotplot, device = "pdf", width = dims$width_in, height = dims$height_in, units = "in")
    }
  )

  output$download_ora_down_table_v2 <- downloadHandler(
    filename = function() { "ora_down.tsv" },
    content = function(file) {
      req(enrichment_result_v2())
      dat <- enrichment_result_v2()$results$ora_down
      utils::write.table(dat, file = file, sep = "\t", quote = FALSE, row.names = FALSE)
    }
  )

  output$download_gsea_plot_v2 <- downloadHandler(
    filename = function() { "gsea_dotplot.pdf" },
    content = function(file) {
      req(enrichment_result_v2())
      dims <- differential_plot_dimensions_v2()
      ggplot2::ggsave(file, plot = enrichment_result_v2()$plots$gsea_dotplot, device = "pdf", width = dims$width_in, height = dims$height_in, units = "in")
    }
  )

  output$download_gsea_table_v2 <- downloadHandler(
    filename = function() { "gsea.tsv" },
    content = function(file) {
      req(enrichment_result_v2())
      dat <- enrichment_result_v2()$results$gsea
      utils::write.table(dat, file = file, sep = "\t", quote = FALSE, row.names = FALSE)
    }
  )

  # Abundance Profiling -------------------------------------------------------

  abundance_plot_dimensions_v2 <- reactive({
    list(
      width_in = if (is.null(input$abundance_plot_width_v2)) 7 else input$abundance_plot_width_v2,
      height_in = if (is.null(input$abundance_plot_height_v2)) 4 else input$abundance_plot_height_v2
    )
  })

  abundance_selected_view_v2 <- reactive({
    if (is.null(input$abundance_main_view_v2)) "Barplot" else input$abundance_main_view_v2
  })

  output$abundance_params_v2 <- renderUI({
    tryCatch({
      validate_prot_data()
      req(prot_data())
      row_choices <- names(SummarizedExperiment::rowData(prot_data()))
      col_choices <- names(SummarizedExperiment::colData(prot_data()))

      switch(
        abundance_selected_view_v2(),
        "Barplot" = div(
          class = "appv2-params-stack",
          selectInput("pv_prot_meta_v2", "Column used for protein labels", choices = row_choices),
          uiOutput("pv_protein_v2"),
          selectInput("pv_condition_v2", "Group samples by", choices = c("No grouping", col_choices)),
          uiOutput("barchart_selected_groups_v2")
        ),
        "Heatmap" = div(
          class = "appv2-params-stack",
          selectInput("protein_label_v2", "Column used for protein labels", choices = row_choices),
          fileInput("heatmap_labels_v2", "Upload heatmap labels"),
          selectInput("heatmap_condition_v2", "Group samples by", choices = c("No grouping", col_choices)),
          checkboxInput("cluster_cols_heatmap_v2", "Cluster columns", value = TRUE),
          checkboxInput("cluster_rows_heatmap_v2", "Cluster rows", value = TRUE)
        )
      )
    }, error = function(e) {
      validate(need(FALSE, paste("Preparing abundance profiling parameters failed:", e$message)))
    })
  })

  output$pv_protein_v2 <- renderUI({
    tryCatch({
      validate_prot_data()
      req(prot_data(), input$pv_prot_meta_v2)
      choices <- SummarizedExperiment::rowData(prot_data())[[input$pv_prot_meta_v2]]
      selectInput("pv_protein_v2", "Select a protein", choices = choices)
    }, error = function(e) {
      validate(need(FALSE, paste("Preparing protein selector failed:", e$message)))
    })
  })

  pv_selected_condition_v2 <- reactive({
    if (is.null(input$pv_condition_v2) || input$pv_condition_v2 == "No grouping") NULL else input$pv_condition_v2
  })

  output$barchart_selected_groups_v2 <- renderUI({
    tryCatch({
      validate_prot_data()
      req(prot_data(), pv_selected_condition_v2())
      choices <- SummarizedExperiment::colData(prot_data())[[pv_selected_condition_v2()]]
      selectInput("barchart_selected_groups_v2", "Groups to display", choices = choices, multiple = TRUE, selected = NULL)
    }, error = function(e) {
      validate(need(FALSE, paste("Preparing group selector failed:", e$message)))
    })
  })

  heatmap_labels_v2 <- reactive({
    file_info <- input$heatmap_labels_v2
    if (is.null(file_info)) {
      return(NULL)
    }
    dat <- data.table::fread(file_info$datapath, data.table = FALSE)
    if (!"Gene" %in% names(dat)) {
      validate(need(FALSE, "The uploaded file must contain a column called Gene."))
    }
    dat$Gene
  })

  heatmap_condition_v2 <- reactive({
    if (is.null(input$heatmap_condition_v2) || input$heatmap_condition_v2 == "No grouping") NULL else input$heatmap_condition_v2
  })

  abundance_barchart_reactive <- reactive({
    req(prot_data(), input$pv_protein_v2, input$pv_prot_meta_v2)
    ProtPipe::compare_protein(
      prot_data(),
      prot = input$pv_protein_v2,
      prot_meta_col = input$pv_prot_meta_v2,
      condition = pv_selected_condition_v2(),
      selected_groups = input$barchart_selected_groups_v2
    )
  })

  abundance_heatmap_reactive <- reactive({
    req(prot_data(), input$protein_label_v2)
    ProtPipe::plot_proteomics_heatmap(
      object = prot_data(),
      protmeta_col = input$protein_label_v2,
      genes = heatmap_labels_v2(),
      condition = heatmap_condition_v2(),
      cluster_cols = input$cluster_cols_heatmap_v2,
      cluster_rows = input$cluster_rows_heatmap_v2
    )
  })

  output$abundance_plot_ui_v2 <- renderUI({
    dims <- abundance_plot_dimensions_v2()
    width_px <- round(dims$width_in * 96)
    height_px <- round(dims$height_in * 96)

    div(
      style = paste0("width:", width_px, "px;"),
      plotOutput("abundance_plot_v2", height = paste0(height_px, "px"))
    )
  })

  output$abundance_plot_v2 <- renderPlot({
    tryCatch({
      validate_prot_data()
      req(prot_data())
      p <- switch(
        abundance_selected_view_v2(),
        "Barplot" = abundance_barchart_reactive(),
        "Heatmap" = abundance_heatmap_reactive()
      )
      print(p)
    }, error = function(e) {
      validate(need(FALSE, paste("Rendering abundance profiling plot failed:", e$message)))
    })
  })

  output$download_abundance_plot_v2 <- downloadHandler(
    filename = function() {
      paste0(gsub(" ", "_", tolower(abundance_selected_view_v2())), ".pdf")
    },
    content = function(file) {
      dims <- abundance_plot_dimensions_v2()
      p <- switch(
        abundance_selected_view_v2(),
        "Barplot" = abundance_barchart_reactive(),
        "Heatmap" = abundance_heatmap_reactive()
      )
      ggsave(file, plot = p, device = "pdf", width = dims$width_in, height = dims$height_in, units = "in")
    }
  )

  output$download_abundance_table_v2 <- downloadHandler(
    filename = function() {
      paste0(gsub(" ", "_", tolower(abundance_selected_view_v2())), ".tsv")
    },
    content = function(file) {
      dat <- switch(
        abundance_selected_view_v2(),
        "Barplot" = as.data.frame(SummarizedExperiment::assay(prot_data())),
        "Heatmap" = as.data.frame(SummarizedExperiment::assay(prot_data()))
      )
      write.table(dat, file = file, sep = "\t", quote = FALSE, row.names = FALSE)
    }
  )

  # Help ----------------------------------------------------------------------

  # Placeholder section for help-specific server logic.

}
