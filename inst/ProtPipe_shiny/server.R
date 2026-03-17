options(shiny.maxRequestSize=5000 * 1024^2)
server <- function(input, output, session) {

  # Update hidden input$select based on which button was clicked
  observeEvent(input$view_0, { updateTextInput(session, "select", value = "0") })
  observeEvent(input$view_1, { updateTextInput(session, "select", value = "1") })
  observeEvent(input$view_2, { updateTextInput(session, "select", value = "2") })
  observeEvent(input$view_3, { updateTextInput(session, "select", value = "3") })
  observeEvent(input$view_4, { updateTextInput(session, "select", value = "4") })
  observeEvent(input$view_5, { updateTextInput(session, "select", value = "5") })
  observeEvent(input$view_6, { updateTextInput(session, "select", value = "6") })

  # Create a temp working directory for this zip session
  zip_workspace <- file.path(tempdir(), "zip_workspace")
  dir.create(zip_workspace, showWarnings = FALSE, recursive = TRUE)
  plot_dirs <- c("quality_control", "clustering", "differential_expression", "pathway_enrichment", "protein_view")

  # Subfolders inside workspace
  subfolders <- file.path(zip_workspace, plot_dirs)

  # Create the folders (with optional .keep files)
  observe({
    unlink(zip_workspace, recursive = TRUE)  # clean old
    dir.create(zip_workspace, showWarnings = FALSE, recursive = TRUE)
    for (sub in subfolders) {
      dir.create(sub, recursive = TRUE)
      file.create(file.path(sub, ".keep"))
    }

    # Create the zip with correct folder structure
    relative_dirs <- plot_dirs  # these are folder names relative to root
    # zip::zip(
    #   zipfile = file.path(zip_workspace, "output.zip"),
    #   files = relative_dirs,
    #   root = zip_workspace,
    #   mode = "cherry-pick"
    # )
    # zip::zip(
    #   zipfile = file.path(zip_workspace, "pathways.zip"),
    #   files = c("pathway_enrichment"),
    #   root = zip_workspace,
    #   mode = "cherry-pick")
  })
  # Download handler
  output$downloadZip <- downloadHandler(
    filename = function() {
      "output.zip"
    },
    content = function(file) {
      relative_dirs <- plot_dirs  # these are folder names relative to root
      zip::zip(
        zipfile = file.path(zip_workspace, "output.zip"),
        files = relative_dirs,
        root = zip_workspace,
        mode = "cherry-pick"
      )
      file.copy(
        from = file.path(zip_workspace, "output.zip"),
        to = file,
        overwrite = TRUE
      )
    },
    contentType = "application/zip"
  )

  #read file uploads
  intensity <- fileUploadServer("intensity")
  sample_condition <- fileUploadServer("sample_condition")
  gene_labels_file <- fileUploadServer("gene_labels")
  heatmap_labels <- fileUploadServer("heatmap_labels")
  ontology_file <- fileUploadServer("ontology_file")

  #### Reactive functions ############################################################################################

  # intensity() is your reactive fileInput, e.g.,
  # intensity <- reactive(input$yourFileInput)

  # 1. Use reactiveValues to store state that can be changed.
  # This is the correct way to manage variables like data frames and their types.
  rv <- reactiveValues(
    data = NULL,
    condition = NULL,
    number_samples = NULL,
    type = "upload file first"
  )

  # 2. Use an observeEvent to perform the ACTION of reading the file
  #    and detecting its type. This runs whenever the file input or
  #    the example checkbox changes.
  observeEvent(list(intensity(), input$use_example), {

    # Logic for loading the example data
    if (input$use_example) {
      rv$type <- "Standard Matrix"
      rv$data <- data.table::fread("www/iPSC.csv", data.table = FALSE)
      return() # Stop execution here
    }

    # If not using an example, require a file upload
    req(intensity())
    file_info <- intensity()

    # Use a tryCatch block for robust error handling during file read
    tryCatch({
      ext <- tolower(tools::file_ext(file_info$datapath))

      # Validate file extension
      validate(
        need(ext %in% c("csv", "tsv", "xlsx", "xls", "adat"), "Invalid file format.")
      )

      # --- Simplified and Corrected Logic ---
      if (ext == "adat") {
        rv$type <- "SomaScan"
        out <- SomaDataIO::read_adat(file_info$datapath) %>%
          ProtPipe::soma_all_output()
        rv$data <- out$data
        rv$condition <- out$condition
        rv$number_samples <- out$number_samples
      } else if (detect_olink_npx(file_info$datapath)) {
        rv$type <- "Olink"
        out <- OlinkAnalyze::read_NPX(file_info$datapath) %>%
          ProtPipe::olink_all_output()
        rv$data <- out$data
        rv$condition <- out$condition
        rv$number_samples <- out$number_samples
      } else {
        # Fallback for standard CSV, TSV, or Excel files
        rv$type <- "Standard Matrix"
        if (ext %in% c("csv", "tsv")) {
          rv$data <- data.table::fread(file_info$datapath, data.table = FALSE)
        } else {
          rv$data <- readxl::read_excel(file_info$datapath)
        }
      }

    }, error = function(e) {
      # Show a user-friendly error if anything in the try block fails
      showNotification(paste("File Processing Error:", e$message), type = "error", duration = 10)
      # Reset reactive values on error
      rv$data <- NULL
      rv$type <- "error"
    })
  })



  # simple reactive expressions to safely access the results.
  intensity_file <- reactive({
    req(rv$data) # Require data to be non-NULL
    rv$data
  })

  data_type <- reactive({
    rv$type
  })

  intermediate_condition <- reactive({
    rv$condition
  })

  number_samples <- reactive({
    rv$number_samples
  })



  output$file_type_output <- renderText({
    # The output will display the string returned by the reactive
    paste("Detected File Type: ", data_type())
  })

  condition_file <- reactive({
    if (!is.null(sample_condition())) {
      # Extract the file extension and convert to lowercase for matching
      ext <- tolower(tools::file_ext(sample_condition()$datapath))
      validate(
        need(ext %in% c("csv", "tsv", "xlsx", "xls"), "Invalid file format. Please upload a .csv, .tsv, or .xlsx file.")
      )
      if (ext %in% c("csv", "tsv")) {
        data.table::fread(sample_condition()$datapath, data.table = FALSE)
      } else { # ext is "xlsx"
        readxl::read_excel(sample_condition()$datapath)
      }
    } else {
      return(NULL)
    }
  })
  # Dynamically generate dropdowns for column range selection
  output$column_range_ui <- renderUI({
    req(intensity_file())
    #req(data_type() == "Standard Matrix")
    df <- intensity_file() %>%
      ProtPipe::convert_numeric_cols()
    choices <- names(df)

    # get default intensity columns
    if(data_type() == "Standard Matrix"){
      intensity_cols <- ProtPipe::detect_intensity_cols(df)
      first <- intensity_cols[[1]]
      last <- intensity_cols[[length(intensity_cols)]]
    }else{
      num_cols <- length(choices)
      first <- num_cols-number_samples()+1
      last <- num_cols
    }

    tagList(
      selectInput("lower_col", "Intensity columns start at:", choices = choices, selected = choices[first]),
      selectInput("upper_col", "Intensity columns end at:", choices = choices, selected = choices[last])
    )
  })


  # Validate and report selection
  output$range_result <- renderPrint({
    req(data_type() == "Standard Matrix")
    df <- intensity_file() %>% ProtPipe::convert_numeric_cols()

    # Ensure both selections are made
    req(input$lower_col, input$upper_col)

    cols <- names(df)
    lower_idx <- match(input$lower_col, cols)
    upper_idx <- match(input$upper_col, cols)

    if (is.na(lower_idx) || is.na(upper_idx)) {
      return("❌ Column not found in file.")
    }

    # Ensure lower <= upper
    if (lower_idx > upper_idx) {
      return("❌ First column must come before or be the same as the last column.")
    }

    selected <- df[, lower_idx:upper_idx, drop = FALSE]

    if (!all(sapply(selected, is.numeric))) {
      return("❌ All selected columns must be numeric.")
    }

    paste("✅ Selected columns:", input$lower_col, "to", input$upper_col,
          "| Count:", ncol(selected))
  })

  output$download_ex <- downloadHandler(
    filename = function(){
      paste("neuron_differentiation.csv")
    },
    content = function(file){
      file.copy("www/iPSC.csv", file)
    }
  )

  raw_prot_data <- reactive({
    req(intensity_file())
    df <- intensity_file()

    # Ensure both selections are made
    #req(input$lower_col, input$upper_col)

    cols <- names(df)
    lower_idx <- match(input$lower_col, cols)
    upper_idx <- match(input$upper_col, cols)

    data_type <- data_type()

    if (data_type == "Standard Matrix") {
      PD <- ProtPipe::create_se(dat = intensity_file(), intensity_cols = c(lower_idx:upper_idx), sample_metadata = condition_file())
    } else {
      if(is.null(condition_file())){
        condition <-intermediate_condition()
      }else if(is.null(intermediate_condition())){
        condition <-condition_file()
      }else{
        condition <- dplyr::left_join(intermediate_condition(), condition_file(), by = "SampleID")
      }
      PD <- ProtPipe::create_se(dat = intensity_file(), intensity_cols = c(lower_idx:upper_idx), sample_metadata = condition)
    }
    return(PD)
  })


  prot_data <- reactive({
    req(raw_prot_data())
    PD <- raw_prot_data()

    #1 min intensity filtering
    if(isTRUE(input$lod_filter)){
      if(data_type() == "Olink"){
        lod_col <- "LOD"
      }
      if(data_type() == "SomaScan"){
        lod_col <- "Buffer"
      }
      PD <- ProtPipe::lod_filter(PD, lod_col)
    }

    if(isTRUE(input$min_int_filter)){
      PD <- ProtPipe::apply_min_intenisty(PD, input$min_int_filter_lod)
    }

    #2 outlier removal
    if(input$remove_outliers == TRUE){
      PD <- ProtPipe::filter_outlier_samples(PD, sds = input$outlier_sds)
    }
    if(input$remove_sparse_proteins == TRUE){
      PD <- ProtPipe::filter_proteins_by_percent(PD, percent = input$sparse_protein_percent)
    }

    #3 transformation
    if(input$log2_transform == TRUE){
      print(paste("Perform", input$normalize_method))
      tryCatch({
        PD <- ProtPipe::log2_transform(PD)
      }, error = function(e) {
        print("Transformation failed")
        print(e)
      })
    }

    #4 normalization
    if(input$normalize == TRUE){
      print(paste("Normalizing using", input$normalize_method))
      tryCatch({
        if (input$normalize_method == "mean") {
          PD <- ProtPipe::mean_normalize(PD)
        } else if (input$normalize_method == "median") {
          PD <- ProtPipe::median_normalize(PD)
        }
      }, error = function(e) {
        print("Normalization failed")
        print(e)
      })
    }

    #5 imputation
    if(input$impute == TRUE){
      if(input$imputation_method == "fixed value"){
        PD <- ProtPipe::impute(PD, input$impute_fixed_value)
      }else if(input$imputation_method == "minimum"){
        PD <- ProtPipe::impute_min(PD, input$impute_min_value)
      }else if(input$imputation_method == "left-shifted distribution"){
        PD <- ProtPipe::impute_left_dist(PD, input$impute_left_dist_shift, input$impute_left_dist_scale)
      }
    }

    #6 batch correction
    if(!is.null(input$batch_correct_column) && input$batch_correct == TRUE){
      PD <- ProtPipe::batch_correct(PD, input$batch_correct_column)
    }

    pdata <- base::cbind(rowData(PD), assay(PD)) %>% as.data.frame()
    add_zip_tabular(pdata, "processed_data.tsv", "quality_control", zip_workspace, "output.zip")

    return(PD)

  })

  output$download_data <- downloadHandler(
    filename = function(){
      paste("processed_data.tsv")
    },
    content = function(file){
      data <- base::cbind(rowData(prot_data()), assay(prot_data()))
      write.table(data, file = file, sep = "\t", quote = FALSE, row.names = FALSE)
    }
  )

  output$download_preprocessing_report <- downloadHandler(
    filename = function() {
      # This part is fine
      "processing_report.md"
    },

    content = function(file) {
      # 'file' is the special temporary path Shiny provides (e.g., "/tmp/RtmpsXYZ/file123.md")

      # 1. Get your reactive data object
      # It's good practice to add req() to ensure data exists before downloading.
      req(prot_data())

      # 2. Pass the special 'file' path from Shiny to your function's argument
      ProtPipe::generate_preprocessing_report(
        object = prot_data(),
        output_file = file  # <-- This is the fix
      )
    }
  )


  ### Pre-processing ############################################################################################

  output$lod_filtering <- renderUI({
    req(intensity_file())
    req(data_type() == "Olink" || data_type() == "SomaScan")
    if(data_type() == "Olink"){
      label <- "filter values based on LOD values (if present)"
    }
    if(data_type() == "SomaScan"){
      label <- "filter values based on average buffer values (if present)"
    }
    checkboxInput("lod_filter", label = label, value = FALSE)
  })

  output$imputation_parameters <- renderUI({
    req(intensity_file())
    if(input$imputation_method == "fixed value"){
      tagList(
        numericInput("impute_fixed_value", "value:", value = 0)
      )
    }else if(input$imputation_method == "minimum"){
      tagList(
        numericInput("impute_min_value", "scale minimum by:", value = 1)
      )
    }else if(input$imputation_method == "left-shifted distribution"){
      tagList(
        numericInput("impute_left_dist_shift", "shift mean of distribution by n standard deviations:", value = 1.8),
        numericInput("impute_left_dist_scale", "scale standard deviation of distribution by:", value = 0.3)
      )
    }
  })

  output$batch_correct_column <- renderUI({
    req(intensity_file())
    req(sample_condition())

    choices <- names(colData(raw_prot_data()))

    selectInput("batch_correct_column", "select condition for correction:", choices = choices)
  })


  #### QC ############################################################################################

  #select condition
  output$quality_control_condition <- renderUI({
    req(intensity_file())
    #req(sample_condition())

    choices <- names(colData(prot_data()))

    selectInput("qc_condition", "select condition to group by:", choices = choices)
  })

  #compute CVs
  cvs_reactive <- reactive({
    req(intensity_file())
    req(input$qc_condition)
    cvs <- tryCatch({
      # This is the "try" block. R will attempt to run this code.
      ProtPipe::get_CVs(prot_data(), condition = input$qc_condition)

    }, error = function(e) {
      # This is the "catch" block. It only runs if an error occurs.
      # We use validate() to display a user-friendly message in the plot area.
      validate(need(FALSE, paste("Calculating CVs failed:", e$message)))
    })
    p <- tryCatch({
      # This is the "try" block. R will attempt to run this code.
      ProtPipe::plot_CVs(prot_data(), condition = input$qc_condition, plot_type = input$cv_plot_type)

    }, error = function(e) {
      # This is the "catch" block. It only runs if an error occurs.
      # We use validate() to display a user-friendly message in the plot area.
      validate(need(FALSE, paste("Plotting CVs failed:", e$message)))
    })
    return(list(cvs = cvs, plot = p))
  })

  # CV plot
  output$cv_graph <- renderPlot({

    cv_data <- cvs_reactive()

    #save data
    add_zip_tabular(cv_data$cvs, "CVs.tsv", "quality_control", zip_workspace, "output.zip")
    add_zip_plot(cv_data$plot, "CV_plot.pdf", "quality_control", zip_workspace, "output.zip")

    #print
    print(cv_data$plot)
  })

  output$download_cv <- downloadHandler(
    filename = function(){
      paste("cv_plot.pdf")
    },
    content = function(file){
      p <- cvs_reactive()$plot
      ggsave(file, plot=p, device = "pdf")
    }
  )

  output$download_cv_tsv <- downloadHandler(
    filename = function(){
      paste("cv.tsv")
    },
    content = function(file){
      cvs <- cvs_reactive()$cvs
      write.table(cvs, file = file, sep = "\t", quote = FALSE, row.names = FALSE)
    }
  )

  #compute protein intensities
  intensities_reactive <- reactive({
    req(intensity_file())
    tryCatch({
      # This is the "try" block. R will attempt to run this code.
      ProtPipe::plot_pg_intensities(prot_data())

    }, error = function(e) {
      # This is the "catch" block. It only runs if an error occurs.
      # We use validate() to display a user-friendly message in the plot area.
      validate(need(FALSE, paste("Plotting intensity failed:", e$message)))
    })
  })

  # intensity graph
  output$intensity_graph <- renderPlot({
    #save plot
    p <-intensities_reactive()
    add_zip_plot(p, "intensities_plot.pdf", "quality_control", zip_workspace, "output.zip")

    #print
    print(p)
  })

  output$download_intensity <- downloadHandler(
    filename = function(){
      paste("intensities.pdf")
    },
    content = function(file){
      p <- intensities_reactive()
      ggsave(file, plot=p, device = "pdf")
    }
  )

  protein_group_counts_reactive <- reactive({
    req(intensity_file())
    pgcounts <- tryCatch({
      # This is the "try" block. R will attempt to run this code.
      get_pg_counts(prot_data())

    }, error = function(e) {
      # This is the "catch" block. It only runs if an error occurs.
      # We use validate() to display a user-friendly message in the plot area.
      validate(need(FALSE, paste("Calculating counts failed:", e$message)))
    })
    p <- tryCatch({
      # This is the "try" block. R will attempt to run this code.
      ProtPipe::plot_pg_counts(prot_data(), condition = input$qc_condition)
    }, error = function(e) {
      # This is the "catch" block. It only runs if an error occurs.
      # We use validate() to display a user-friendly message in the plot area.
      validate(need(FALSE, paste("Plotting protein counts failed:", e$message)))
    })
    return(list(counts = pgcounts, plot = p))
  })

  # protein group counts
  output$pgroup_graph <- renderPlot({
    req(intensity_file())
    pg_counts <- protein_group_counts_reactive()

    #save data
    add_zip_tabular(pg_counts$counts, "pg_counts.tsv", "quality_control", zip_workspace, "output.zip")
    add_zip_plot(pg_counts$plot, "pg_groups_plot.pdf", "quality_control", zip_workspace, "output.zip")
    print(pg_counts$plot)
  })

  output$download_pg <- downloadHandler(
    filename = function(){
      paste("protein_groups_nonzero_counts.pdf")
    },
    content = function(file){
      p <- protein_group_counts_reactive()$plot
      ggsave(file, plot=p, device = "pdf")
    }
  )

  output$download_pg_tsv <- downloadHandler(
    filename = function(){
      paste("protein_group_nonzero_counts.tsv")
    },
    content = function(file){
      pgcounts <- protein_group_counts_reactive()$counts
      write.table(pgcounts, file = file, sep = "\t", row.names = FALSE, quote = FALSE)
    }
  )

  #compute sample correlations
  correlation_heatmap_reactive <- reactive({
    req(intensity_file())
    dat.correlations <- tryCatch({
      # This is the "try" block. R will attempt to run this code.
      ProtPipe::get_sample_correlation(prot_data())
    }, error = function(e) {
      # This is the "catch" block. It only runs if an error occurs.
      # We use validate() to display a user-friendly message in the plot area.
      validate(need(FALSE, paste("Calculating sample correlation failed:", e$message)))
    })
    p <- tryCatch({
      # This is the "try" block. R will attempt to run this code.
      ProtPipe::plot_correlation_heatmap(prot_data())
    }, error = function(e) {
      # This is the "catch" block. It only runs if an error occurs.
      # We use validate() to display a user-friendly message in the plot area.
      validate(need(FALSE, paste("Plotting sample correlation failed:", e$message)))
    })
    return(list(correlations = dat.correlations, plot = p))
  })

  #correlation heatmap
  output$correlation_graph <- renderPlot({
    req(intensity_file())
    sample_correlations <- correlation_heatmap_reactive()

    #save data
    add_zip_tabular(sample_correlations$correlations, "sample_correlations.tsv", "quality_control", zip_workspace, "output.zip")
    add_zip_plot(sample_correlations$plot, "sample_correlation_heatmap.pdf", "quality_control", zip_workspace, "output.zip")

    print(sample_correlations$plot)
  })

  output$download_cor <- downloadHandler(
    filename = function(){
      paste("sample_correlation.pdf")
    },
    content = function(file){
      p<-correlation_heatmap_reactive()$plot
      ggsave(file, plot=p, device = "pdf")
    }
  )

  output$download_cor_tsv <- downloadHandler(
    filename = function(){
      paste("sample_correlation.tsv")
    },
    content = function(file){
      dat.correlations <- correlation_heatmap_reactive()$correlations
      write.table(dat.correlations , file = file, sep = "\t", quote = FALSE, row.names = FALSE)
    }
  )

  #### Clustering ############################################################################################

  #select condition
  folder_path <- file.path(zip_workspace, 'clustering') #used to add data to zip

  output$clustering_condition <- renderUI({
    req(intensity_file())
    #req(sample_condition())

    choices <- names(colData(prot_data()))

    selectInput("cluster_condition", "select condition to group by:", choices = choices)
  })

  hcluster_reactive <- reactive({
    req(intensity_file())
    tryCatch({
      # This is the "try" block. R will attempt to run this code.
      ProtPipe::plot_hierarchical_cluster(prot_data())

    }, error = function(e) {
      # This is the "catch" block. It only runs if an error occurs.
      # We use validate() to display a user-friendly message in the plot area.
      validate(need(FALSE, paste("Clustering failed:", e$message)))
    })
  })

  #hierarchical clustering
  output$hcluster <- renderPlot({
    req(intensity_file())  # Ensure file is uploaded
    p <- hcluster_reactive()

    #save data to zip
    add_zip_plot(p, "hierarchical_clustering.pdf", "clustering", zip_workspace, "output.zip")
    #print plot
    print(p)
  })

  output$download_hcluster <- downloadHandler(
    filename = function(){
      paste("hierarchical_clustering.pdf")
    },
    content = function(file){
      p <-  hcluster_reactive()
      ggsave(file, plot=p, device = "pdf")
    }
  )

  pca_reactive <- reactive({
    req(intensity_file())
    tryCatch({
      # This is the "try" block. R will attempt to run this code.
      ProtPipe::plot_PCs(prot_data(), condition = input$cluster_condition)

    }, error = function(e) {
      # This is the "catch" block. It only runs if an error occurs.
      # We use validate() to display a user-friendly message in the plot area.
      validate(need(FALSE, paste("Clustering failed:", e$message)))
    })
  })

  pca_data_reactive <- reactive({
    req(intensity_file())
    tryCatch({
      # This is the "try" block. R will attempt to run this code.
      get_PCs(prot_data(), condition = input$cluster_condition)

    }, error = function(e) {
      # This is the "catch" block. It only runs if an error occurs.
      # We use validate() to display a user-friendly message in the plot area.
      validate(need(FALSE, paste("Computing PCA failed:", e$message)))
    })
  })

  #PCA
  output$pca <- renderPlot({
    req(intensity_file())  # Ensure file is uploaded

    p <- pca_reactive()
    components <- pca_data_reactive()$components
    summary <- pca_data_reactive()$summary
    #save data to zip
    add_zip_plot(p, "PCA.pdf", "clustering", zip_workspace, "output.zip")
    add_zip_tabular(components, "pca_components.tsv", "clustering", zip_workspace, "output.zip")
    add_zip_tabular(summary, "pca_summary.tsv", "clustering", zip_workspace, "output.zip")
    #print plot
    print(p)

  })

  output$download_pca <- downloadHandler(
    filename = function(){
      paste("PCA.pdf")
    },
    content = function(file){
      p <- pca_reactive()
      ggsave(file, plot=p, device = "pdf")
    }
  )

  output$download_pca_tsv <- downloadHandler(
    filename = function(){
      paste('PCA.tsv')
    },
    content = function(file){
      req(intensity_file())  # Ensure file is uploaded
      data.table::fwrite(pca_data_reactive()$components, file, sep = "\t")
    }
  )
  output$download_pca_sum <- downloadHandler(
    filename = function(){
      paste('PCA_summary.tsv')
    },
    content = function(file){
      req(intensity_file())  # Ensure file is uploaded
      data.table::fwrite(pca_data_reactive()$summary, file, sep = "\t")
    }
  )
  output$neighbors_slider <- renderUI({
    req(intensity_file())
    max = ncol(prot_data())
    default = min(15, (max+2)/2)
    div(style = "display: flex; align-items: center; height: 100%;",
        sliderInput("neighbors", label = h3("Select number of neighbors for UMAP"),
                    min = 2, max = max, step = 1, value = default))
  })

  # 1. Create a reactive to handle the calculation and saving.
  #    This only runs when 'input$neighbors' or 'input$cluster_condition' changes,
  #    NOT when the user resizes the window.
  umap_results <- reactive({
    req(intensity_file())
    req(input$neighbors)

    # Set a seed so the Plot and the Table data are mathematically identical
    set.seed(123)

    tryCatch({
      # A. Generate the Data
      # We generate the data first to ensure consistency
      umap_table <- ProtPipe::get_umap(
        prot_data(),
        neighbors = input$neighbors,
        condition = input$cluster_condition
      )

      # B. Generate the Plot
      # (Assuming plot_umap re-runs calculations, we set the seed above to match.
      #  Ideally, plot_umap would accept 'umap_table' as input to be perfectly safe.)
      p <- ProtPipe::plot_umap(
        prot_data(),
        neighbors = input$neighbors,
        condition = input$cluster_condition
      )

      # C. Save to Zip (Only happens once per input change)
      add_zip_plot(p, "umap.pdf", "clustering", zip_workspace, "output.zip")
      add_zip_tabular(umap_table, "umap_summary.tsv", "clustering", zip_workspace, "output.zip")

      return(list(plot = p, table = umap_table))

    }, error = function(e) {
      # Return a specific error message string if it fails
      return(paste("Clustering failed:", e$message))
    })
  })

  # 2. The render function is now lightweight.
  #    It just shows the result of the reactive above.
  output$umap <- renderPlot({
    # Get the result from the reactive
    result <- umap_results()

    # Check if the result is a string (which means it's our error message)
    # If it is, trigger the validate/stop mechanism
    if (is.character(result)) {
      validate(need(FALSE, result))
    }

    print(result$plot)
  })

  output$download_umap <- downloadHandler(
    filename = function(){
      paste("umap.pdf")
    },
    content = function(file){
      req(intensity_file())  # Ensure file is uploaded
      ggsave(file, plot=umap_results()$plot, device = "pdf")
    }
  )

  output$download_umap_tsv <- downloadHandler(
    filename = function(){
      paste("umap.tsv")
    },
    content = function(file){
      req(intensity_file())  # Ensure file is uploaded
      data.table::fwrite(umap_results()$table, file, sep = "\t")
    }
  )

  #### Differential Intensity ############################################################################################

  #select condition
  output$de_condition <- renderUI({
    req(intensity_file())
    choices <- names(colData(prot_data()))
    selectInput("de_condition", "select the column used to compare groups:", choices = choices)
  })

  #select groups
  output$de_groups <- renderUI({
    req(intensity_file())
    req(input$de_condition)
    req(input$outcome_type == "binary")
    groups <- unique(colData(prot_data())[[input$de_condition]])

    tagList(
      selectInput("control_condition", "select the control groups:", choices = groups),
      selectInput("treatment_condition", "select the treatment groups:", choices = groups)
    )
  })

  output$de_covariates <- renderUI({
    req(intensity_file())
    req(input$outcome_type == "binary")
    choices <- names(colData(prot_data()))
    selectInput("de_covariates", "select the covariates:", choices = choices,multiple = TRUE,selected = NULL)
  })


  #select column to label the proteins
  output$label_col <- renderUI({
    req(intensity_file())
    choices <- names(rowData(prot_data()))
    selectInput("label_col", "select the column used to label proteins:", choices = choices)
  })
  output$logfc <- renderUI({
    if (input$outcome_type == "continuous"){
      x<-0.35
      label <- "Enter spearman rho cutoff"
    }else{
        x<-1
        label <- "Enter log2 fold-change cutoff"
    }
    numericInput("logfc", label = label, value = x)
  })

  output$pvalue <- renderUI({
    if (isTRUE(input$use_adj_pval)){
      label <- "Enter adjusted P-value cutoff"
    }else{
      label <- "Enter P-value cutoff"
    }
    numericInput("pvalue", label = label, value = 0.01)
  })


  #custom labels for volcano
  gene_labels <- reactive({
    req(intensity_file())  # Ensure file is uploaded
    if(is.null(gene_labels_file())){
      return(NULL)
    }
    dat <- data.table::fread(gene_labels_file()$datapath, data.table=FALSE)

    if (!'Genes' %in% names(dat)) {
      validate("The uploaded gene label file must contain a column called 'Genes'.")
      return(NULL)
    }

    return(dat$Genes)
  })

  dea_result <- reactiveVal(NULL)

  dea <- reactive({
    condition <- input$de_condition
    control_group <- input$control_condition
    treatment_group <- input$treatment_condition
    covariates <- input$de_covariates

    result <- tryCatch({
      if (input$outcome_type == "continuous") {
        ProtPipe::do_comparison_continuous(prot_data(), condition = condition)
      } else {
        ProtPipe::do_limma_binary(
          prot_data(),
          condition = condition,
          control_group = control_group,
          treatment_group = treatment_group,
          covariates = covariates
        )
      }
    }, error = function(e) {
      dea_result(NULL)
      validate(need(FALSE, paste("Calculating limma failed:", e$message)))
    })

    dea_result(result)
    result
  })

  output$dea_ready <- shiny::renderText({
    !is.null(dea_result())
  })
  outputOptions(output, "dea_ready", suspendWhenHidden = FALSE)

  # compute volcano plot
  volcano_reactive <- reactive({
    req(intensity_file())
    tryCatch({
      if (input$outcome_type == "continuous"){
        ProtPipe::plot_correlation_volcano(dea(), label_col = input$label_col,
                                           labelgene = gene_labels(),
                                           fdr_threshold = input$pvalue,
                                           rho_threshold = input$logfc,
                                           adj = input$use_adj_pval)
      }else{
        ProtPipe::plot_volcano(dea(), label_col = input$label_col,
                               labelgene = gene_labels(),
                               fdr_threshold = input$pvalue,
                               lfc_threshold = input$logfc,
                               adj = input$use_adj_pval)}
    }, error = function(e) {
      # This is the "catch" block. It only runs if an error occurs.
      # We use validate() to display a user-friendly message in the plot area.
      validate(need(FALSE, paste("Plotting vocano failed:", e$message)))
    })
  })


  #render volcano plot
  output$volcano <- renderPlot({
    req(intensity_file())
    p <- volcano_reactive()
    #save data to zip
    add_zip_plot(p, "volcano_plot.pdf", "differential_expression", zip_workspace, "output.zip")
    add_zip_tabular(dea(), "differential_expression.tsv", "differential_expression", zip_workspace, "output.zip")

    print(p)
  })

  #download volcano plot
  output$download_volcano <- downloadHandler(
    filename = function(){
      paste("volcano.pdf")
    },
    content = function(file){
      p <- volcano_reactive()
      ggsave(file, plot=p, device = "pdf")
    }
  )

  output$download_DE_tsv <- downloadHandler(
    filename = function(){
      paste("differential_expression_results.tsv")
    },
    content = function(file){
      req(intensity_file())  # Ensure file is uploaded
      write.table(dea(), file = file, sep = "\t", row.names = FALSE, quote = FALSE)
    }
  )

  #select condition
  output$gene_col <- renderUI({
    req(intensity_file())
    choices <- names(rowData(prot_data()))
    selectInput("gene_col", "select the column containing official gene symbols (e.g., TP53):", choices = choices)
  })

  #go and kegg db for different organisms
  selected_org <- reactive({
    req(input$organism)
    organism_map[[input$organism]]
  })

  selected_ontology <- reactive({
    req(input$pathway_source)

    if (identical(input$pathway_source, "GO")) {
      return(NULL)
    }

    req(ontology_file())
    ProtPipe::read_ontology(ontology_file()$datapath)
  })

  # Create a reactiveVal to store pathway enrichment results
  enrichment_result <- reactiveVal(NULL)
  enrichment_message <- reactiveVal(NULL)

  observeEvent(input$run_enrichment, {
    if (isTRUE(input$run_enrichment)) {
      # Disable the checkbox/button (if checkboxInput used as button, or use actionButton)
      shinyjs::disable("run_enrichment")  # requires shinyjs package and call to use it in UI

      # Optionally show a notification
      showNotification("Running enrichment analysis, please wait...", duration = NULL, id = "enrich_msg")

      tryCatch({
        req(dea_result())
        ontology <- selected_ontology()
        result <- ProtPipe::enrich_pathways(
          dea_result(),
          lfc_threshold = input$logfc,
          fdr_threshold = input$pvalue,
          enrich_pvalue = input$enrich_pval,
          go_org = selected_org()$OrgDb,
          kegg_org = selected_org()$kegg,
          gene_col = input$gene_col,
          adj = input$use_adj_pval,
          source = if (identical(input$pathway_source, "GO")) "go" else "custom",
          go_ont = input$go_ontology,
          term2gene = if (is.null(ontology)) NULL else ontology$term2gene,
          term2name = if (is.null(ontology)) NULL else ontology$term2name,
          run_ora = isTRUE(input$run_ora),
          run_gsea = isTRUE(input$run_gsea)
        )
        enrichment_result(result)
        if (is.null(result)) {
          enrichment_message("No genes were mapped to Entrez IDs. Check that the selected gene column contains official gene symbols.")
        } else if ("message" %in% names(result$results)) {
          enrichment_message(result$results$message$message[[1]])
        } else {
          enrichment_message(NULL)
        }
      }, error = function(e) {
        enrichment_result(NULL)
        enrichment_message(paste("Pathway enrichment failed:", e$message))
      }, finally = {
        removeNotification("enrich_msg")
        shinyjs::enable("run_enrichment")
      })
    }
  })

  #pathway enrichment plots
  output$ora_up_enrich <- renderPlot({
    req(intensity_file())
    validate(need(is.null(enrichment_message()), enrichment_message()))
    req(enrichment_result())
    enrichment_result()$plots$ora_up_dotplot
  })

  output$ora_down_enrich <- renderPlot({
    req(intensity_file())
    validate(need(is.null(enrichment_message()), enrichment_message()))
    req(enrichment_result())
    enrichment_result()$plots$ora_down_dotplot
  })

  output$gsea_enrich <- renderPlot({
    req(intensity_file())
    validate(need(is.null(enrichment_message()), enrichment_message()))
    req(enrichment_result())
    validate(
      need(
        !("gsea_message" %in% names(enrichment_result()$results)),
        enrichment_result()$results$gsea_message$message[[1]]
      )
    )
    enrichment_result()$plots$gsea_dotplot
  })

  #save all enrichment results to temp zip
  observe({
    req(intensity_file())
    req(enrichment_result())

    # Save data frames
    for (name in names(enrichment_result()$results)) {
      df <- enrichment_result()$results[[name]]
      if (!is.null(df) && nrow(df) > 0) {
        add_zip_tabular(df, paste0(name, ".tsv"), "pathway_enrichment", zip_workspace, "output.zip")
        gc()
      } else {
        message(paste("Skipping empty or null dataframe:", name))
      }
    }

    # Save plots
    for (name in names(enrichment_result()$plots)) {
      p <- enrichment_result()$plots[[name]]
      if (!is.null(p) && nrow(p$data) > 0) {
        add_zip_plot(p, paste0(name, ".pdf"), "pathway_enrichment", zip_workspace, "output.zip")
        gc()
      } else {
        message(paste("Skipping empty or null plot:", name))
      }
    }
  })

  # download zip file of all pathway enrichment plots and data
  output$download_enrichment <- downloadHandler(
    filename = function() {
      "pathwaysg.zip"
    },
    content = function(file) {
      relative_dirs <- plot_dirs  # these are folder names relative to root
      zip::zip(
        zipfile = file.path(zip_workspace, "pathways.zip"),
        files = c("pathway_enrichment"),
        root = zip_workspace,
        mode = "cherry-pick")
      file.copy(
        from = file.path(zip_workspace, "pathways.zip"),
        to = file,
        overwrite = TRUE
      )
    },
    contentType = "application/zip"
  )

  #### Heatmap ############################################################################################


  #select condition
  output$protein_label <- renderUI({
    req(intensity_file())
    choices <- names(rowData(prot_data()))
    selectInput("protein_label", "select the column used to label proteins:", choices = choices)
  })

  #heatmap subset
  prot_labels <- reactive({
    req(intensity_file())  # Ensure file is uploaded
    req(input$protein_label)  # Ensure file is uploaded
    if(is.null(heatmap_labels())){
      return(NULL)
    }
    dat <- data.table::fread(heatmap_labels()$datapath, data.table=FALSE)
    if('Gene' %in% names(dat)){
      print("The csv must contain a column called Gene containing the labels")
    }
    return(dat$Gene)
  })

  output$heatmap_condition <- renderUI({
    req(intensity_file())
    choices <- c("no goruping", names(colData(prot_data())))
    selectInput("heatmap_condition", "Select condition to group samples:", choices = choices)
  })

  heatmap_condition <- reactive({
    req(intensity_file())  # Ensure file is uploaded
    req(input$heatmap_condition)  # Ensure file is uploaded
    if(is.null(input$heatmap_condition) || input$heatmap_condition == "no goruping"){
      return(NULL)
    }
    return(input$heatmap_condition)
  })

  #compute plot reactively
  heatmap_reactive <- reactive({
    req(intensity_file())

    tryCatch({
      ProtPipe::plot_proteomics_heatmap(
        object = prot_data(),
        protmeta_col = input$protein_label,
        genes = prot_labels(),
        condition = heatmap_condition(),
        cluster_cols = input$cluster_cols_heatmap,
        cluster_rows = input$cluster_rows_heatmap
      )
    }, error = function(e) {
      # If this fails, it stops here and sends the message to any output using it
      validate(need(FALSE, paste("Plotting heatmap failed:", e$message)))
    })
  })

  # Render the Plot (Calls the reactive)
  output$h_map <- renderPlot({
    p <- heatmap_reactive() # Retrieves the cached plot

    # Save to zip (Side effect)
    add_zip_plot(p, "heatmap.pdf", "protein_view", zip_workspace, "output.zip")
    print(p)
  })

  #Download Handler (Calls the SAME reactive)
  output$download_hmap <- downloadHandler(
    filename = function(){
      paste("heatmap.pdf")
    },
    content = function(file){
      p <- heatmap_reactive() # Retrieves the cached plot
      ggsave(file, plot = p, device = "pdf")
    }
  )

  #### Protein Barchart ############################################################################################

  #select condition
  output$pv_prot_meta <- renderUI({
    req(intensity_file())
    choices <- names(rowData(prot_data()))
    selectInput("pv_prot_meta", "select the column used to label proteins:", choices = choices)
  })

  #select condition
  output$pv_protein <- renderUI({
    req(intensity_file())
    req(input$pv_prot_meta)
    choices <- rowData(prot_data())[[input$pv_prot_meta]]
    selectInput("pv_protein", "select a protein:", choices = choices)
  })

  #select condition
  output$pv_condition <- renderUI({
    req(intensity_file())
    choices <- c("No grouping", names(colData(prot_data())))
    selectInput("pv_condition", "select the column used to group samples:", choices = choices)
  })

  pv_selected_condition <-reactive({
    if(input$pv_condition == "No grouping"){
      NULL
    }else{
      input$pv_condition
    }
  })

  output$barchart_selected_groups <- renderUI({
    req(intensity_file())
    req(pv_selected_condition())
    choices <- colData(prot_data())[[pv_selected_condition()]]
    selectInput("barchart_selected_groups", "select groups to display:", choices = choices, multiple = TRUE,selected = NULL)
  })

  protein_barchart_reactive <- reactive({
    req(intensity_file())

    tryCatch({
      ProtPipe::compare_protein(prot_data(),
                                prot = input$pv_protein,
                                prot_meta_col = input$pv_prot_meta,
                                condition = pv_selected_condition(),
                                selected_groups = input$barchart_selected_groups
                                )
    }, error = function(e) {
      # If this fails, it stops here and sends the message to any output using it
      validate(need(FALSE, paste("Plotting protein barchart failed:", e$message)))
    })
  })

  #complete barchart
  output$protein_barchart <- renderPlot({
    p <- protein_barchart_reactive()
    #save to zip
    add_zip_plot(p, paste(input$pv_protein, "_levels.pdf"), "protein_view", zip_workspace, "output.zip")
    print(p)
  })

  output$download_protein_barchart <- downloadHandler(
    filename = function(){
      paste(input$pv_protein, "_levels.pdf")
    },
    content = function(file){
      p <- protein_barchart_reactive()
      ggsave(file, plot=p, device = "pdf")
    }
  )


}
