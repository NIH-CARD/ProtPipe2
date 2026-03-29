
#options(shiny.maxRequestSize=5000 * 1024^2)
source("helpers.R")

workflow_header <- function(title) {
  card(
    card_header(
      div(
        style = "width: 100%;",
        div(
          style = "display: flex; align-items: center; justify-content: space-between; gap: 1rem; width: 100%;",
          div(
            actionButton("back_page", "Back", class = "btn-default")
          ),
          div(
            actionButton("next_page", "Next", class = "btn-primary")
          )
        ),
        div(
          style = "width: 100%; text-align: center; font-size: 1.75rem; font-weight: 700; line-height: 1.2; margin-top: 0.75rem;",
          title
        )
      )
    )
  )
}


ui <- page_sidebar(
  useShinyjs(),

  tags$head(
    # 2. Add the meta tag
    tags$meta(name="google-site-verification", content="_Y-vcB0obJKbpSV7gjLqwfiS-lPu-EYNaFJyFjsQUSk")
  ),

  # Title and subtitle
  title = tagList(
    h1("ProtPipe2", style = "margin-bottom: 0;"),
    tags$div("This website is free and open to all users. Shiny app made by Jacob Epstein", style = "font-size: 0.9em; color: #666; margin-top: 0.2em;")
  ),

  # SIDEBAR CONTENT GOES HERE
  sidebar = sidebar(
    open = "closed",
    h3("Select view"),

    # Hidden textInput used for conditionalPanel control
    tags$div(style = "display:none;", textInput("select", label = NULL, value = "0")),

    # Block-style buttons
    actionButton("view_0", "1. Input", class = "btn-block btn-primary mb-2"),
    actionButton("view_1", "2. Quality Control", class = "btn-block btn-primary mb-2"),
    actionButton("view_2", "3. Pre-Processing", class = "btn-block btn-primary mb-2"),
    actionButton("view_3", "4. Clustering / Dimensionality Reduction", class = "btn-block btn-primary mb-2"),
    actionButton("view_4", "5. Differential Intensity", class = "btn-block btn-primary mb-2"),
    actionButton("view_5", "6. Abundance Profiling", class = "btn-block btn-primary mb-2"),
    actionButton("view_6", "Help", class = "btn-block btn-primary mb-2"),

    hr(),
    downloadButton("downloadZip", "Download All Plots")
  ),

  # MAIN PANEL CONTENT
  main = tagList(
    conditionalPanel("input.select == '0'", h4("Input parameters content")),
    conditionalPanel("input.select == '1'", h4("Quality Control content")),
    conditionalPanel("input.select == '2'", h4("Pre Processing content")),
    conditionalPanel("input.select == '3'", h4("Clustering content")),
    conditionalPanel("input.select == '4'", h4("Differential Intensity content")),
    conditionalPanel("input.select == '5'", h4("Protein content")),
    conditionalPanel("input.select == '6'", h4("Help content"))
  ),

  ### Parameter input screen ############################################################################################
  conditionalPanel(condition = "input.select == 0",
                   fluidPage(
                       workflow_header("Input"),
                       card(
                         card_header(h4("Protein Intensity File")),
                         fluidRow(
                           column(width = 6, fileUploadUI("intensity", label = NULL)),
                           column(width = 6, verbatimTextOutput("file_type_output", T))),
                         fluidRow(
                           column(width = 6, checkboxInput("use_example", "Or use our iPSC to neuron differentiation example dataset", value = FALSE),),
                           column(width = 6, downloadButton("download_ex", "Download example dataset"))),
                         uiOutput("column_range_ui"),
                         verbatimTextOutput("range_result")),
                       card(
                         card_header(h4("Sample Condition File")),
                         p("make sure row names match the column names of the intensity file exactly"),
                         fileUploadUI("sample_condition", label = NULL)
                       )
                   )
  ),


  ### Quality control screen ############################################################################################
  conditionalPanel(condition = "input.select == 1",
                   fluidPage(
                     workflow_header("Quality Control"),
                     uiOutput("quality_control_condition"),
                     card(card_header("Coefficients of Variation (requires condition file)"),
                          selectInput("cv_plot_type", "Select format:", choices = c("violin", "jitter"), selected = "violin"),
                          plotOutput("cv_graph"),
                          downloadButton("download_cv", "Download Plot as PDF"),
                          downloadButton("download_cv_tsv", "Download data as tsv")
                     ),
                     card(card_header("Intensities"),
                          plotOutput("intensity_graph"),
                          downloadButton("download_intensity", "Download Plot as PDF")
                     ),
                     card(card_header("Non-zero Protein Group Counts"),
                          plotOutput("pgroup_graph"),
                          downloadButton("download_pg", "Download Plot as PDF"),
                          downloadButton("download_pg_tsv", "Download data as tsv")
                     ),
                     card(card_header("Correlation Heatmap"),
                          plotOutput("correlation_graph"),
                          downloadButton("download_cor", "Download Plot as PDF"),
                          downloadButton("download_cor_tsv", "Download data as tsv")
                     )
                   )
  ),
  ### Pre Processing Screen ############################################################################################
  conditionalPanel(condition = "input.select == 2",
                   fluidPage(
                     workflow_header("Pre-Processing"),
                          card(card_header("1. Minimum Intensity Filtering"),
                               fluidRow(
                                  column(width=6,
                                      checkboxInput("min_int_filter", label = "set minimum intensity", value = FALSE),
                                      numericInput("min_int_filter_lod", label="min: ", value = 0)),
                                  column(width=6,uiOutput("lod_filtering"))
                               )),
                          card(card_header("2. Outlier Removal"),
                            fluidRow(
                            column(width = 6,
                                   checkboxInput("remove_outliers", label = "remove outlier samples", value = FALSE),
                                   numericInput("outlier_sds", label = "Remove samples with protein groups outside n standard deviations from the mean", value = 3)),
                            column(width = 6,
                                   checkboxInput("remove_sparse_proteins", label = "remove outlier proteins", value = FALSE),
                                   numericInput("sparse_protein_percent", label = "Remove proteins present in less than n% of samples", value = 30)))),
                          card(card_header("3. Transformation"),
                               checkboxInput("log2_transform", label = "log2_transform", value = FALSE)),
                          card(card_header("4. Normalization"),
                               checkboxInput("normalize", label = "normalize", value = FALSE),
                               selectInput("normalize_method", label = "normalize_method", choices = c("mean", "median"), selected = "median")),
                          card(card_header("5. Imputation"),
                               fluidRow(
                                 column(width = 6,
                                        checkboxInput("impute", label = "impute", value = TRUE),
                                        selectInput("imputation_method", label = "imputation method", choices = c("fixed value", "minimum", "left-shifted distribution"), selected = "fixed value")),
                                 column(width = 6,
                                        conditionalPanel(
                                          condition = "input.imputation_method == 'fixed value'",
                                          numericInput("impute_fixed_value", "value:", value = 0)
                                        ),
                                        conditionalPanel(
                                          condition = "input.imputation_method == 'minimum'",
                                          numericInput("impute_min_value", "scale minimum by:", value = 1)
                                        ),
                                        conditionalPanel(
                                          condition = "input.imputation_method == 'left-shifted distribution'",
                                          numericInput("impute_left_dist_shift", "shift mean of distribution by n standard deviations:", value = 1.8),
                                          numericInput("impute_left_dist_scale", "scale standard deviation of distribution by:", value = 0.3)
                                        )
                                 ))),
                          card(card_header("6. Batch Correction"),
                               checkboxInput("batch_correct", label = "batch correct", value = FALSE),
                               uiOutput("batch_correct_column")),
                          card(
                            downloadButton("download_data", "Download pre-processed data"),
                            downloadButton("download_preprocessing_report", "Download pre-processing report"))
                    )),
  ### Clustering screen ############################################################################################
  conditionalPanel(condition = "input.select == 3",
                   fluidPage(
                     workflow_header("Clustering / Dimensionality Reduction"),
                     uiOutput("clustering_condition"),
                     card(card_header("heirarchial clustering"),
                          plotOutput("hcluster"),
                          downloadButton("download_hcluster", "Download Plot as PDF")
                     ),card(card_header("PCA (requires condition file)"),
                            plotOutput("pca"),
                            downloadButton("download_pca", "Download Plot as PDF"),
                            downloadButton("download_pca_tsv", "Download PCA as tsv"),
                            downloadButton("download_pca_sum", "Download PCA summary as tsv")
                     ),card(card_header("UMAP (requires condition file)"),
                            fluidRow(
                              column(width = 3, uiOutput("neighbors_slider")),
                              column(width = 9, plotOutput("umap"))
                            ),
                            downloadButton("download_umap", "Download Plot as PDF"),
                            downloadButton("download_umap_tsv", "Download UMAP as tsv")
                     )
                   )
  ),
  ### Differential Intensity ############################################################################################
  conditionalPanel(condition = "input.select == 4",
                   workflow_header("Differential Expression"),
                   card(card_header("Options"),
                        fluidPage(
                          radioButtons(
                            inputId = "outcome_type",
                            label = "Select outcome type:",
                            choices = c("Categorical" = "binary", "Continuous" = "continuous"),
                            selected = "binary"
                          ),
                          fluidRow(
                            column(width = 4,
                                   uiOutput("de_condition"),
                                   uiOutput("de_groups"),
                                   uiOutput("de_covariates")),
                            column(width = 4,
                                   # this is called logfc but can be spearman coef if outcome is continuous
                                   uiOutput("logfc"),
                                   checkboxInput(inputId = "use_adj_pval",
                                                 label = "Use Adjusted P-value (FDR)",
                                                 value = TRUE),
                                   uiOutput("pvalue")),
                                   #numericInput("pvalue", label = "Enter P-value cutoff", value = 0.01)),
                            column(width = 4,
                                   uiOutput("label_col"),
                                   p("optional: upload csv file with genes to label in volcano plot.
                                   Make sure column name is Genes. This must contain values present in the column
                                  selected above"),
                                   fileUploadUI("gene_labels", label = NULL))
                          )
                        )
                   ),
                   card(card_header("Volcano Plot"),
                        plotOutput("volcano"),
                        downloadButton("download_volcano", "Download Plot as PDF"),
                        downloadButton("download_DE_tsv", "Download differential expression results as tsv")
                   ),
                   conditionalPanel(
                     condition = "output.dea_ready",
                     card(card_header("Pathway Enrichment Options"),
                        fluidRow(
                          column(width=6,
                                 numericInput("enrich_pval", label = "Enter enrichment pvalue cutoff", value = 0.05)
                                 ,
                                 selectInput("organism", "Select Organism:", choices = names(organism_map), selected = "Human"),
                                 uiOutput("gene_col")
                          ),
                          column(width=6,
                                 actionButton("run_enrichment", "Run Pathway Analysis", class = "btn-primary"),
                                 div(style = "margin-top: 0.75rem;"),
                                 selectInput("pathway_source", "Ontology source:", choices = c("GO", "Custom ontology"), selected = "GO"),
                                 conditionalPanel(
                                   condition = "input.pathway_source == 'GO'",
                                   selectInput("go_ontology", "GO ontology:", choices = c("BP", "MF", "CC"), selected = "BP")
                                 ),
                                 conditionalPanel(
                                   condition = "input.pathway_source == 'Custom ontology'",
                                   p("Upload a GMT, TSV, or CSV ontology file using Entrez gene IDs."),
                                   p("Required formats: GMT with term, name, then genes; TSV/CSV with either term/gene or term/name/gene columns."),
                                   p("Gene IDs in the uploaded ontology must match the IDs used for enrichment. ProtPipe currently uses Entrez IDs for custom ontology analysis."),
                                   fileUploadUI("ontology_file", label = NULL)
                                 )
                          )
                        )
                     )
                   ),
                   card(card_header("Pathway Enrichment"),
                        fluidRow(
                          column(width=6 ,card(card_header("Upregulated Pathways"),plotOutput("ora_up_enrich"))),
                          column(width=6 ,card(card_header("Downregulated Pathways"),plotOutput("ora_down_enrich")))
                        ),
                        fluidRow(
                          column(width=12 ,card(card_header("Gene Set Enrichment"),plotOutput("gsea_enrich")))
                        ),
                        downloadButton("download_enrichment", "Download pathway enrichment results")
                     )
                   ),

  ### Protein View ############################################################################################
  conditionalPanel(condition = "input.select == 5",
                   workflow_header("Abundance Profiling"),
                   card(card_header("Heatmap"),
                     uiOutput("protein_label"),
                     uiOutput("heatmap_condition"),
                     checkboxInput(inputId = "cluster_cols_heatmap", label = "Cluster Columns", value = FALSE),
                     checkboxInput(inputId = "cluster_rows_heatmap", label = "Cluster Rows", value = FALSE),
                     p("optional: upload csv file with genes to include in heatmap subset.
                                   Make sure column name is Genes"),
                     fileUploadUI("heatmap_labels", label = NULL),
                     card(
                          plotOutput("h_map"),
                          downloadButton("download_hmap", "Download Plot as PDF")
                     )
                   ),

                   card(card_header("Protein Barchart"),
                     uiOutput("pv_prot_meta"),
                     uiOutput("pv_protein"),
                     uiOutput("pv_condition"),
                     uiOutput("barchart_selected_groups"),
                     card(
                          plotOutput("protein_barchart"),
                          downloadButton("download_protein_barchart", "Download Plot as PDF")
                     )
                   )
  ),

  ### Help ############################################################################################
  conditionalPanel(condition = "input.select == 6",
                   tags$iframe(
                     src = "help.html",
                     style = "width: 100%; height: 100vh; border: none; display: block;"
                   )
  )
)
