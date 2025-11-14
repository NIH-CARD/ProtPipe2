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
  plot_dirs <- c("quality_control", "clustering", "differential_expression", "pathway_enrichment")

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

  #### Reactive functions ############################################################################################

  # intensity() is your reactive fileInput, e.g.,
  # intensity <- reactive(input$yourFileInput)

  # 1. Use reactiveValues to store state that can be changed.
  # This is the correct way to manage variables like data frames and their types.
  rv <- reactiveValues(
    data = NULL,
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
        rv$data <- SomaDataIO::read_adat(file_info$datapath)
      } else if (detect_olink_npx(file_info$datapath)) {
        rv$type <- "Olink"
        # **CORRECTED**: Use an Olink reader or a generic one.
        # olinkanalyze::read_npx is the standard for Olink files.
        rv$data <- OlinkAnalyze::read_NPX(file_info$datapath)
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

  # 3. Create simple reactive expressions to safely access the results.
  #    Your original `intensity_file` now just returns the data.
  intensity_file <- reactive({
    req(rv$data) # Require data to be non-NULL
    rv$data
  })

  # A new reactive to access the type
  data_type <- reactive({
    rv$type
  })

  output$file_type_output <- renderText({
    # The output will display the string returned by the reactive
    data_type()
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
    req(data_type() == "Standard Matrix")
    df <- intensity_file() %>%
      ProtPipe::convert_numeric_cols()
    choices <- names(df)

    # get default intensity columns
    intensity_cols <- ProtPipe::detect_intensity_cols(df)
    #intensity_cols <<- intensity_cols
    first <- intensity_cols[[1]]
    last <- intensity_cols[[length(intensity_cols)]]

    tagList(
      selectInput("lower_col", "Intensity columns start at:", choices = choices, selected = choices[first]),
      selectInput("upper_col", "Intensity columns end at:", choices = choices, selected = choices[last])
    )
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
    } else if(data_type == "SomaScan"){
      PD <- ProtPipe::create_se_from_soma(adat = intensity_file(), condition = condition_file(), filter = T)
    } else if(data_type == "Olink"){
      PD <- ProtPipe::create_se_from_olink(npx = intensity_file(), condition = condition_file(), filter = T)
    }
    return(PD)
  })


  prot_data <- reactive({
    req(raw_prot_data())
    PD <- raw_prot_data()
    #1 outlier removal
    if(input$remove_outliers == TRUE){
      PD <- ProtPipe::filter_outlier_samples(PD, sds = input$outlier_sds)
    }
    if(input$remove_sparse_proteins == TRUE){
      PD <- ProtPipe::filter_proteins_by_percent(PD, percent = input$sparse_protein_percent)
    }

    #2 transformation
    if(input$log2_transform == TRUE){
      print(paste("Perform", input$normalize_method))
      tryCatch({
        PD <- ProtPipe::log2_transform(PD)
      }, error = function(e) {
        print("Transformation failed")
        print(e)
      })
    }

    #3 normalization
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

    #4 imputation
    if(input$impute == TRUE){
      if(input$imputation_method == "fixed value"){
        PD <- ProtPipe::impute(PD, input$impute_fixed_value)
      }else if(input$imputation_method == "minimum"){
        PD <- ProtPipe::impute_min(PD, input$impute_min_value)
      }else if(input$imputation_method == "left-shifted distribution"){
        PD <- ProtPipe::impute_left_dist(PD, input$impute_left_dist_shift, input$impute_left_dist_scale)
      }
    }

    #5 batch correction
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
      data <- base::cbind(prot_data()@prot_meta, prot_data()@data)
      write.table(data, file = file, sep = "\t", quote = FALSE, row.names = FALSE)
    }
  )





  ### Color Pallete ############################################################################################


  ### Batch Correction ############################################################################################

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

  # CV plot
  output$cv_graph <- renderPlot({
    req(intensity_file())
    #req(sample_condition())
    req(input$qc_condition)

    #save tabular data
    cvs <- tryCatch({
      # This is the "try" block. R will attempt to run this code.
      ProtPipe::get_CVs(prot_data(), condition = input$qc_condition)

    }, error = function(e) {
      # This is the "catch" block. It only runs if an error occurs.
      # We use validate() to display a user-friendly message in the plot area.
      validate(need(FALSE, paste("Calculating CVs failed:", e$message)))
    })

    tryCatch({
      # This is the "try" block. R will attempt to run this code.
      add_zip_tabular(cvs, "CVs.tsv", "quality_control", zip_workspace, "output.zip")

    }, error = function(e) {
      # This is the "catch" block. It only runs if an error occurs.
      # We use validate() to display a user-friendly message in the plot area.
      validate(need(FALSE, paste("Calculating CVs failed:", e$message)))
    })


    #save plot
    p <- tryCatch({
      # This is the "try" block. R will attempt to run this code.
      ProtPipe::plot_CVs(prot_data(), condition = input$qc_condition, plot_type = input$cv_plot_type)

    }, error = function(e) {
      # This is the "catch" block. It only runs if an error occurs.
      # We use validate() to display a user-friendly message in the plot area.
      validate(need(FALSE, paste("Plotting CVs failed:", e$message)))
    })

    add_zip_plot(p, "CV_plot.pdf", "quality_control", zip_workspace, "output.zip")

    #print
    print(p)
  })

  output$download_cv <- downloadHandler(
    filename = function(){
      paste("cv_plot.pdf")
    },
    content = function(file){
      p <- ProtPipe::plot_CVs(prot_data(), condition = input$qc_condition, plot_type = input$cv_plot_type)
      ggsave(file, plot=p, device = "pdf")
    }
  )

  output$download_cv_tsv <- downloadHandler(
    filename = function(){
      paste("cv.tsv")
    },
    content = function(file){
      cvs <- ProtPipe::get_CVs(prot_data(), condition = input$qc_condition)
      write.table(cvs, file = file, sep = "\t", quote = FALSE, row.names = FALSE)
    }
  )

  # intensity graph
  output$intensity_graph <- renderPlot({
    req(intensity_file())

    #save plot
    p <- tryCatch({
      # This is the "try" block. R will attempt to run this code.
      ProtPipe::plot_pg_intensities(prot_data())

    }, error = function(e) {
      # This is the "catch" block. It only runs if an error occurs.
      # We use validate() to display a user-friendly message in the plot area.
      validate(need(FALSE, paste("Plotting intensity failed:", e$message)))
    })

    add_zip_plot(p, "intensities_plot.pdf", "quality_control", zip_workspace, "output.zip")

    #print
    print(p)
  })

  output$download_intensity <- downloadHandler(
    filename = function(){
      paste("intensities.pdf")
    },
    content = function(file){
      p <- ProtPipe::plot_pg_intensities(prot_data())
      ggsave(file, plot=p, device = "pdf")
    }
  )

  # protein group counts
  output$pgroup_graph <- renderPlot({
    req(intensity_file())

    #save tabular data
    pgcounts <- tryCatch({
      # This is the "try" block. R will attempt to run this code.
      get_pg_counts(prot_data())

    }, error = function(e) {
      # This is the "catch" block. It only runs if an error occurs.
      # We use validate() to display a user-friendly message in the plot area.
      validate(need(FALSE, paste("Calculating counts failed:", e$message)))
    })

    add_zip_tabular(pgcounts, "pg_counts.tsv", "quality_control", zip_workspace, "output.zip")

    #save plot
    p <- tryCatch({
      # This is the "try" block. R will attempt to run this code.
      ProtPipe::plot_pg_counts(prot_data())
    }, error = function(e) {
      # This is the "catch" block. It only runs if an error occurs.
      # We use validate() to display a user-friendly message in the plot area.
      validate(need(FALSE, paste("Plotting counts failed:", e$message)))
    })

    add_zip_plot(p, "pg_groups_plot.pdf", "quality_control", zip_workspace, "output.zip")
    print(p)
  })

  output$download_pg <- downloadHandler(
    filename = function(){
      paste("protein_groups_nonzero_counts.pdf")
    },
    content = function(file){
      p <- ProtPipe::plot_pg_counts(prot_data())
      ggsave(file, plot=p, device = "pdf")
    }
  )

  output$download_pg_tsv <- downloadHandler(
    filename = function(){
      paste("protein_group_nonzero_counts.tsv")
    },
    content = function(file){
      pgcounts <- get_pg_counts(prot_data())
      write.table(pgcounts, file = file, sep = "\t", row.names = FALSE, quote = FALSE)
    }
  )

  #correlation heatmap
  output$correlation_graph <- renderPlot({
    req(intensity_file())

    #save tabular data
    dat.correlations <- tryCatch({
      # This is the "try" block. R will attempt to run this code.
      ProtPipe::get_sample_correlation(prot_data())
    }, error = function(e) {
      # This is the "catch" block. It only runs if an error occurs.
      # We use validate() to display a user-friendly message in the plot area.
      validate(need(FALSE, paste("Calculating sample correlation failed:", e$message)))
    })
    add_zip_tabular(dat.correlations, "sample_correlations.tsv", "quality_control", zip_workspace, "output.zip")

    #save plot
    p <- tryCatch({
      # This is the "try" block. R will attempt to run this code.
      ProtPipe::plot_correlation_heatmap(prot_data())
    }, error = function(e) {
      # This is the "catch" block. It only runs if an error occurs.
      # We use validate() to display a user-friendly message in the plot area.
      validate(need(FALSE, paste("Plotting sample correlation failed:", e$message)))
    })
    add_zip_plot(p, "sample_correlation_heatmap.pdf", "quality_control", zip_workspace, "output.zip")

    print(p)
  })

  output$download_cor <- downloadHandler(
    filename = function(){
      paste("sample_correlation.pdf")
    },
    content = function(file){
      p<-ProtPipe::plot_correlation_heatmap(prot_data())
      ggsave(file, plot=p, device = "pdf")
    }
  )

  output$download_cor_tsv <- downloadHandler(
    filename = function(){
      paste("sample_correlation.tsv")
    },
    content = function(file){
      dat.correlations <- ProtPipe::get_sample_correlation(prot_data())
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

  #hierarchical clustering
  output$hcluster <- renderPlot({
    req(intensity_file())  # Ensure file is uploaded

    p <- tryCatch({
      # This is the "try" block. R will attempt to run this code.
      ProtPipe::plot_hierarchical_cluster(prot_data())

    }, error = function(e) {
      # This is the "catch" block. It only runs if an error occurs.
      # We use validate() to display a user-friendly message in the plot area.
      validate(need(FALSE, paste("Clustering failed:", e$message)))
    })


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
      req(intensity_file())  # Ensure file is uploaded
      p <-  ProtPipe::plot_hierarchical_cluster(prot_data())
      ggsave(file, plot=p, device = "pdf")
    }
  )

  #PCA
  output$pca <- renderPlot({
    req(intensity_file())  # Ensure file is uploaded

    p <- tryCatch({
      # This is the "try" block. R will attempt to run this code.
      ProtPipe::plot_PCs(prot_data(), condition = input$cluster_condition)

    }, error = function(e) {
      # This is the "catch" block. It only runs if an error occurs.
      # We use validate() to display a user-friendly message in the plot area.
      validate(need(FALSE, paste("Clustering failed:", e$message)))
    })
    #save data to zip
    add_zip_plot(p, "PCA.pdf", "clustering", zip_workspace, "output.zip")
    add_zip_tabular(get_PCs(prot_data(), condition = input$cluster_condition)$components, "pca_components.tsv", "clustering", zip_workspace, "output.zip")
    add_zip_tabular(get_PCs(prot_data(), condition = input$cluster_condition)$summary, "pca_summary.tsv", "clustering", zip_workspace, "output.zip")
    #print plot
    print(p)

  })

  output$download_pca <- downloadHandler(
    filename = function(){
      paste("PCA.pdf")
    },
    content = function(file){
      p <- ProtPipe::plot_PCs(prot_data(), condition = input$cluster_condition)
      ggsave(file, plot=p, device = "pdf")
    }
  )

  output$download_pca_tsv <- downloadHandler(
    filename = function(){
      paste('PCA.tsv')
    },
    content = function(file){
      req(intensity_file())  # Ensure file is uploaded
      data.table::fwrite(get_PCs(prot_data(), condition = input$cluster_condition)$components, file, sep = "\t")
    }
  )
  output$download_pca_sum <- downloadHandler(
    filename = function(){
      paste('PCA_summary.tsv')
    },
    content = function(file){
      req(intensity_file())  # Ensure file is uploaded
      data.table::fwrite(get_PCs(prot_data(), condition = input$cluster_condition)$summary, file, sep = "\t")
    }
  )

  #UMAP
  output$umap <- renderPlot({
    req(intensity_file())  # Ensure file is uploaded'
    #req(sample_condition())

    p <- tryCatch({
      # This is the "try" block. R will attempt to run this code.
      ProtPipe::plot_umap(prot_data(), neighbors = input$neighbors, condition = input$cluster_condition)

    }, error = function(e) {
      # This is the "catch" block. It only runs if an error occurs.
      # We use validate() to display a user-friendly message in the plot area.
      validate(need(FALSE, paste("Clustering failed:", e$message)))
    })

    #save data to zip
    add_zip_plot(p, "umap.pdf", "clustering", zip_workspace, "output.zip")
    add_zip_tabular(get_umap(prot_data(), neighbors = input$neighbors, condition = input$cluster_condition), "umap_summary.tsv", "clustering", zip_workspace, "output.zip")

    #plot data
    print(p)
  })

  output$download_umap <- downloadHandler(
    filename = function(){
      paste("umap.pdf")
    },
    content = function(file){
      req(intensity_file())  # Ensure file is uploaded
      p <- ProtPipe::plot_umap(prot_data(), neighbors = input$neighbors, condition = input$cluster_condition)
      ggsave(file, plot=p, device = "pdf")
    }
  )

  output$download_umap_tsv <- downloadHandler(
    filename = function(){
      paste("umap.tsv")
    },
    content = function(file){
      req(intensity_file())  # Ensure file is uploaded
      data.table::fwrite(get_umap(prot_data(), neighbors = input$neighbors, condition = input$cluster_condition), file, sep = "\t")
    }
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

  #complete heatmap
  output$h_map <- renderPlot({
    req(intensity_file())  # Ensure file is uploaded
    p <- tryCatch({
      # This is the "try" block. R will attempt to run this code.
      ProtPipe::plot_proteomics_heatmap(prot_data(), protmeta_col = input$protein_label, genes = prot_labels())
    }, error = function(e) {
      # This is the "catch" block. It only runs if an error occurs.
      # We use validate() to display a user-friendly message in the plot area.
      validate(need(FALSE, paste("Plotting heatmap failed:", e$message)))
    })

    #save to zip
    add_zip_plot(p, "heatmap.pdf", "quality_control", zip_workspace, "output.zip")

    print(p)
    # grid::grid.newpage()
    # grid::grid.draw(p$gtable)
  })

  output$download_hmap <- downloadHandler(
    filename = function(){
      paste("heatmap.pdf")
    },
    content = function(file){
      req(intensity_file())  # Ensure file is uploaded
      p <- ProtPipe::plot_proteomics_heatmap(prot_data(), protmeta_col = input$protein_label, genes = prot_labels())
      ggsave(file, plot=p, device = "pdf")
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

  pv_condition <-reactive({
    if(input$pv_condition == "No grouping"){
      NULL
    }else{
      input$pv_condition
    }
  })

  #complete barchart
  output$protein_barchart <- renderPlot({
    req(intensity_file())  # Ensure file is uploaded
    p <- tryCatch({
      # This is the "try" block. R will attempt to run this code.
      ProtPipe::compare_protein(prot_data(), prot = input$pv_protein, prot_meta_col = input$pv_prot_meta, condition = pv_condition())
    }, error = function(e) {
      # This is the "catch" block. It only runs if an error occurs.
      # We use validate() to display a user-friendly message in the plot area.
      validate(need(FALSE, paste("Plotting barchart failed:", e$message)))
    })
    #save to zip
    add_zip_plot(p, paste(input$pv_protein, "_levels.pdf"), "quality_control", zip_workspace, "output.zip")

    print(p)
  })

  output$download_protein_barchart <- downloadHandler(
    filename = function(){
      paste(input$pv_protein, "_levels.pdf")
    },
    content = function(file){
      req(intensity_file())  # Ensure file is uploaded
      p <- ProtPipe::compare_protein(prot_data(), prot = input$pv_protein, prot_meta_col = input$pv_prot_meta, condition = pv_condition())
      ggsave(file, plot=p, device = "pdf")
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
    choices <- names(colData(prot_data()))
    selectInput("de_covariates", "select the covariaes:", choices = choices,multiple = TRUE,selected = NULL)
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


  #custom labels for volcano
  gene_labels <- reactive({
    req(intensity_file())  # Ensure file is uploaded
    if(is.null(gene_labels_file())){
      return(NULL)
    }
    dat <- data.table::fread(gene_labels_file()$datapath, data.table=FALSE)

    if (!'Gene' %in% names(dat)) {
      warning("The uploaded gene label file must contain a column called 'Gene'.")
      return(NULL)
    }

    return(dat$Gene)
  })

  dea <- reactive({
    condition <- input$de_condition
    control_group <- input$control_condition
    treatment_group <- input$treatment_condition
    covariates <- input$de_covariates

    dea <-tryCatch({
        # This is the "try" block. R will attempt to run this code.
      if (input$outcome_type == "continuous"){
        (ProtPipe::do_comparison_continuous(prot_data()
                                   ,condition = condition,
                                   covariates = covariates))
      }else{(ProtPipe::do_limma_binary(prot_data()
                                       ,condition = condition,
                                       control_group = control_group,
                                       treatment_group = treatment_group,
                                       covariates = covariates))}
      }, error = function(e) {
        # This is the "catch" block. It only runs if an error occurs.
        # We use validate() to display a user-friendly message in the plot area.
        validate(need(FALSE, paste("Calculating limma failed:", e$message)))
      })
  })

  #volcano plot
  output$volcano <- renderPlot({
    req(intensity_file())
    p <- tryCatch({
      if (input$outcome_type == "continuous"){
        ProtPipe::plot_correlation_volcano(dea(), label_col = input$label_col,
                                 labelgene = gene_labels(),
                                 fdr_threshold = input$pvalue,
                                 rho_threshold = input$logfc)
      }else{
        ProtPipe::plot_volcano(dea(), label_col = input$label_col,
                               labelgene = gene_labels(),
                               fdr_threshold = input$pvalue,
                               lfc_threshold = input$logfc)}
    }, error = function(e) {
      # This is the "catch" block. It only runs if an error occurs.
      # We use validate() to display a user-friendly message in the plot area.
      validate(need(FALSE, paste("Plotting vocano failed:", e$message)))
    })

    #save data to zip
    add_zip_plot(p, "volcano_plot.pdf", "differential_expression", zip_workspace, "output.zip")
    add_zip_tabular(dea(), "differential_expression.tsv", "differential_expression", zip_workspace, "output.zip")

    print(p)
  })

  output$download_volcano <- downloadHandler(
    filename = function(){
      paste("volcano.pdf")
    },
    content = function(file){
      req(intensity_file())  # Ensure file is uploaded
      p <- ProtPipe::plot_volcano(dea(), label_col = input$label_col, labelgene = gene_labels(), fdr_threshold = input$pvalue, lfc_threshold = input$logfc)
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

  # Create a reactiveVal to store pathway enrichment results
  enrichment_result <- reactiveVal(NULL)

  observeEvent(input$run_enrichment, {
    if (isTRUE(input$run_enrichment)) {
      # Disable the checkbox/button (if checkboxInput used as button, or use actionButton)
      shinyjs::disable("run_enrichment")  # requires shinyjs package and call to use it in UI

      # Optionally show a notification
      showNotification("Running enrichment analysis, please wait...", duration = NULL, id = "enrich_msg")

      # Run the long function (blocking)
      result <- ProtPipe::enrich_pathways(dea(), lfc_threshold=input$logfc, fdr_threshold=input$pvalue, enrich_pvalue=input$enrich_pval, go_org = selected_org()$OrgDb, kegg_org = selected_org()$kegg, gene_col = input$gene_col)
      enrichment_result(result)

      # Remove notification
      removeNotification("enrich_msg")

      # Re-enable button so user can rerun if needed
      shinyjs::enable("run_enrichment")
    }
  })

  #pathway enrichment plots
  output$go_up_enrich <- renderPlot({
    req(intensity_file())
    req(enrichment_result())
    enrichment_result()$plots$go_up_dotplot
  })

  output$kegg_up_enrich <- renderPlot({
    req(intensity_file())
    req(enrichment_result())
    enrichment_result()$plots$kegg_up_dotplot
  })
  output$go_down_enrich <- renderPlot({
    req(intensity_file())
    req(enrichment_result())
    enrichment_result()$plots$go_down_dotplot
  })

  output$kegg_down_enrich <- renderPlot({
    req(intensity_file())
    req(enrichment_result())
    enrichment_result()$plots$kegg_down_dotplot
  })
  output$go_gsea <- renderPlot({
    req(intensity_file())
    req(enrichment_result())
    enrichment_result()$plots$gse_go_dotplot
  })

  output$kegg_gsea <- renderPlot({
    req(intensity_file())
    req(enrichment_result())
    enrichment_result()$plots$gse_kegg_dotplot
  })

  #save all enrichment results to temp zip
  observe({
    req(intensity_file())
    req(enrichment_result())

    # Save data frames
    for (name in names(enrichment_result()$results)) {
      #df <- enrichment_result()$results[[name]]
      if (!is.null(df)) {
        add_zip_tabular(enrichment_result()$results[[name]], paste0(name, ".tsv"), "pathway_enrichment", zip_workspace, "output.zip")
        #add_zip_tabular(enrichment_result()$results[[name]], paste0(name, ".tsv"), "pathway_enrichment", zip_workspace, "pathways.zip")
        gc()
      }
    }

    # Save plots
    for (name in names(enrichment_result()$plots)) {
      #p <- enrichment_result()$plots[[name]]
      if (!is.null(p)) {
        add_zip_plot(enrichment_result()$plots[[name]], paste0(name, ".pdf"), "pathway_enrichment", zip_workspace, "output.zip")
        #add_zip_plot(enrichment_result()$plots[[name]], paste0(name, ".pdf"), "pathway_enrichment", zip_workspace, "pathways.zip")
        gc()
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


}
