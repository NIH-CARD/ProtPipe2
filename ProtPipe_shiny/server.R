
server <- function(input, output, session) {

  # Update hidden input$select based on which button was clicked
  observeEvent(input$view_0, { updateTextInput(session, "select", value = "0") })
  observeEvent(input$view_1, { updateTextInput(session, "select", value = "1") })
  observeEvent(input$view_2, { updateTextInput(session, "select", value = "2") })
  observeEvent(input$view_3, { updateTextInput(session, "select", value = "3") })
  observeEvent(input$view_4, { updateTextInput(session, "select", value = "4") })

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

  #if sample conditions are provided, the data is reformatted so that the intensity columns
  #match the condition rows

  intensity_file <- reactive({
    # This part is unchanged. It requires either the example checkbox
    # to be checked or a file to be uploaded.
    req(input$use_example || !is.null(intensity()))

    # This logic for loading the example data is also unchanged.
    if(input$use_example){
      return(data.table::fread("www/iPSC.csv", data.table=FALSE))
    } else {
      # --- START OF EDITED SECTION ---

      # Get information about the uploaded file
      file_info <- intensity() # Avoid calling the reactive multiple times

      # Extract the file extension and convert to lowercase for matching
      ext <- tolower(tools::file_ext(file_info$datapath))

      # Use validate() to stop execution and show a user-friendly message
      # if the file extension is not one of the allowed types.
      validate(
        need(ext %in% c("csv", "tsv", "xlsx", "adat"), "Invalid file format. Please upload a .csv, .tsv, or .xlsx file.")
      )

      # Use the correct function based on the file extension
      tryCatch({

        # First, validate the inputs to make sure the combination is valid
        # This provides specific feedback to the user.
        if (input$data_type == 1 && !(ext %in% c("csv", "tsv", "xlsx"))) {
          validate(need(FALSE, "For Mass Spec, please upload a .csv, .tsv, or .xlsx file."))
        } else if (input$data_type == 2 && ext != "adat") { # Assuming 2 is SomaScan
          validate(need(FALSE, "For SomaScan, please upload an .adat file."))
        } else if (input$data_type == 3 && !(ext %in% c("csv", "tsv", "xlsx"))) {
          validate(need(FALSE, "For Olink, please upload a .csv, .tsv, or .xlsx file."))
        }

        # If validation passes, proceed to read the file
        if (input$data_type == 1) { # Mass Spec
          if (ext %in% c("csv", "tsv")) {
            data.table::fread(file_info$datapath, data.table = FALSE)
          } else { # ext is "xlsx"
            readxl::read_excel(file_info$datapath)
          }
        } else if (input$data_type == 2) { # SomaScan
          SomaDataIO::read_adat(file_info$datapath)
        } else if (input$data_type == 3) { # Olink
          OlinkAnalyze::read_NPX(file_info$datapath)
        }

        # This 'error' function is the key to preventing crashes.
        # It catches any error from the reading functions (e.g., read_NPX)
        # and displays it as a safe validation message instead of crashing.
      }, error = function(e) {
        validate(need(FALSE, paste("File processing error:", e$message)))
      })
    }
  })
  condition_file <- reactive({
    if (!is.null(sample_condition())) {
      return(data.table::fread(sample_condition()$datapath, data.table=FALSE))
    } else {
      return(NULL)
    }
  })
  # Dynamically generate dropdowns for column range selection
  output$column_range_ui <- renderUI({
    req(intensity_file())
    df <- intensity_file()
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

  # Validate and report selection
  output$range_result <- renderPrint({
    df <- intensity_file()

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

  raw_prot_data <- reactive({
    req(intensity_file())
    df <- intensity_file()

    # Ensure both selections are made
    req(input$lower_col, input$upper_col)

    cols <- names(df)
    lower_idx <- match(input$lower_col, cols)
    upper_idx <- match(input$upper_col, cols)

    # if(any(grepl("NPX", intensity_file(), ignore.case = TRUE))){
    #   PD <- ProtPipe::create_protdata_from_olink(dat = intensity_file(), condition = condition_file())
    #   print("yooooooo")
    # }
    data_type <- input$data_type
    if(input$use_example){
      data_type <- 1
    }


    if (data_type == 1) {
      PD <- ProtPipe::create_protdata(dat = intensity_file(), intensity_cols = c(lower_idx:upper_idx), condition = condition_file())
    } else if(data_type == 2){
      PD <- ProtPipe::create_protdata_from_soma(adat = intensity_file(), condition = condition_file())
    } else if(data_type == 3){
      PD <- ProtPipe::create_protdata_from_olink(npx = intensity_file(), condition = condition_file())
    }
    return(PD)
  })
  prot_data <- reactive({
    req(raw_prot_data())
    PD <- raw_prot_data()
    #1 outlier removal
    if(input$remove_outliers == TRUE){
      PD <- ProtPipe::remove_outliers(PD, sds = input$outlier_sds)
    }

    #2 normalization
    if(input$normalize == TRUE){
      print(paste("Normalizing using", input$normalize_method))
      tryCatch({
        if (input$normalize_method == "mean") {
          PD <- ProtPipe::mean_normalize(PD)
        } else if (input$normalize_method == "median") {
          PD <- ProtPipe::median_normalize(PD)
        }
        PD1 <<- PD
      }, error = function(e) {
        print("Normalization failed")
        print(e)
      })
    }

    #3 imputation
    if(input$impute == TRUE){
      if(input$imputation_method == "zero"){
        PD <- ProtPipe::impute(PD, 0)
      }else if(input$imputation_method == "minimum"){
        PD <- ProtPipe::impute_min(PD, 1)
      }else if(input$imputation_method == "left-shifted distribution"){
        PD <- ProtPipe::impute_minimal(PD)
      }
    }

    #4 batch correction
    if(!is.null(input$batch_correct_column) && input$batch_correct == TRUE){
      PD <- ProtPipe::batch_correct(PD, input$batch_correct_column)
    }
    #PD2 <<- PD
    return(PD)

  })



  ### Color Pallete ############################################################################################


  ### Batch Correction ############################################################################################

  output$batch_correct_column <- renderUI({
    req(intensity_file())
    req(sample_condition())

    choices <- names(raw_prot_data()@condition)

    selectInput("batch_correct_column", "select condition for correction:", choices = choices)
  })


  #### QC ############################################################################################

  #select condition
  output$quality_control_condition <- renderUI({
    req(intensity_file())
    #req(sample_condition())

    choices <- names(prot_data()@condition)

    selectInput("qc_condition", "select condition to group by:", choices = choices)
  })

  # CV plot
  output$cv_graph <- renderPlot({
    req(intensity_file())
    #req(sample_condition())
    req(input$qc_condition)

    #save tabular data
    cvs <- ProtPipe::get_CVs(prot_data(), condition = input$qc_condition)
    add_zip_tabular(cvs, "CVs.tsv", "quality_control", zip_workspace, "output.zip")

    #save plot
    p <- ProtPipe::plot_CVs(prot_data(), condition = input$qc_condition, plot_type = input$cv_plot_type)
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
    p <- ProtPipe::plot_pg_intensities(prot_data())
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
    pgcounts <- get_pg_counts(prot_data())
    add_zip_tabular(pgcounts, "pg_counts.tsv", "quality_control", zip_workspace, "output.zip")

    #save plot
    p <- ProtPipe::plot_pg_counts(prot_data())
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
    dat.correlations <- ProtPipe::get_spearman(prot_data())
    add_zip_tabular(dat.correlations, "sample_correlations.tsv", "quality_control", zip_workspace, "output.zip")

    #save plot
    p <- ProtPipe::plot_correlation_heatmap(prot_data())
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
      dat.correlations <- ProtPipe::get_spearman(prot_data())
      write.table(dat.correlations , file = file, sep = "\t", quote = FALSE, row.names = FALSE)
    }
  )

  #### Clustering ############################################################################################

  #select condition
  folder_path <- file.path(zip_workspace, 'clustering') #used to add data to zip

  output$clustering_condition <- renderUI({
    req(intensity_file())
    #req(sample_condition())

    choices <- names(prot_data()@condition)

    selectInput("cluster_condition", "select condition to group by:", choices = choices)
  })

  #hierarchical clustering
  output$hcluster <- renderPlot({
    req(intensity_file())  # Ensure file is uploaded
    p <- ProtPipe::plot_hierarchical_cluster(prot_data())

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
    p <- ProtPipe::plot_PCs(prot_data(), condition = input$cluster_condition)
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
    p <- ProtPipe::plot_umap(prot_data(), neighbors = input$neighbors, condition = input$cluster_condition)

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
    choices <- names(prot_data()@prot_meta)
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
    p <- ProtPipe::plot_proteomics_heatmap(prot_data(), protmeta_col = input$protein_label, genes = prot_labels())

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
  #### Differential Intensity ############################################################################################

  #select condition
  output$de_condition <- renderUI({
    req(intensity_file())
    choices <- names(prot_data()@condition)
    selectInput("de_condition", "select the column used to compare groups:", choices = choices)
  })

  #select groups
  output$de_groups <- renderUI({
    req(intensity_file())
    req(input$de_condition)
    groups <- unique(prot_data()@condition[[input$de_condition]])

    tagList(
      selectInput("control_condition", "select the control groups:", choices = groups),
      selectInput("treatment_condition", "select the treatment groups:", choices = groups)
    )
  })

  #select column to label the proteins
  output$label_col <- renderUI({
    req(intensity_file())
    choices <- names(prot_data()@prot_meta)
    selectInput("label_col", "select the column used to label proteins:", choices = choices)
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
    df <- ProtPipe::log2_transform(prot_data())
    condition <- input$de_condition
    control_group <- input$control_condition
    treatment_group <- input$treatment_condition

    return(ProtPipe::do_limma_by_condition(df,condition = condition, control_group = control_group, treatment_group = treatment_group))
  })

  #volcano plot
  output$volcano <- renderPlot({
    req(intensity_file())
    p <- ProtPipe::plot_volcano(dea(), label_col = input$label_col, labelgene = gene_labels(), fdr_threshold = input$pvalue, lfc_threshold = input$logfc)

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
    choices <- names(prot_data()@prot_meta)
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
