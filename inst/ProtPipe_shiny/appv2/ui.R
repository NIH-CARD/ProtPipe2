source("shared.R")

appv2_page_layout <- function(param_ui, main_ui, show_main = TRUE, param_card_class = NULL) {
  fluidRow(
    style = "height: 100%; margin-left: 0; margin-right: 0;",
    column(
      width = if (show_main) 4 else 12,
      class = "appv2-col",
      div(
        class = "appv2-panel appv2-sidebar",
        div(
          class = paste("appv2-card appv2-params-card", param_card_class),
          tags$p(class = "appv2-title", "Parameters"),
          div(class = "appv2-params-body", param_ui)
        )
      )
    ),
    if (show_main) {
      column(width = 8, class = "appv2-col", main_ui)
    }
  )
}

input_main_ui <- div(
  class = "appv2-panel appv2-plot",
  div(
    class = "appv2-view-shell",
    div(
      class = "appv2-panel appv2-top-card",
      tags$p(class = "appv2-title", "Data View"),
      div(
        class = "appv2-top-row",
        div(
          class = "appv2-top-control",
          selectInput("input_main_view_v2", "Select view", choices = c("Protein Intensity", "Sample Metadata", "Protein Metadata"))
        ),
        div(
          class = "appv2-top-copy",
          verbatimTextOutput("file_summary_v2", placeholder = TRUE)
        )
      )
    ),
    div(
      class = "appv2-main-body",
      div(
        class = "appv2-panel appv2-plot",
        div(class = "appv2-table-wrap", tableOutput("input_preview_table"))
      )
    )
  )
)

generic_main_ui <- function(page_id, title) {
  div(
    class = "appv2-panel appv2-plot",
    div(
      class = "appv2-view-shell",
      div(
        class = "appv2-panel appv2-top-card",
        tags$p(class = "appv2-title", title),
        div(
          class = "appv2-top-row",
          div(
            class = "appv2-top-control",
            selectInput(
              paste0(page_id, "_main_view_v2"),
              "Select view",
              choices = c("Plot", "Assay", "Sample Metadata", "Protein Metadata")
            )
          ),
          div(class = "appv2-top-copy")
        )
      ),
      div(class = "appv2-main-body", uiOutput("workflow_main_body"))
    )
  )
}

quality_control_main_ui <- div(
  class = "appv2-panel appv2-plot",
  div(
    class = "appv2-view-shell",
    div(
      class = "appv2-panel appv2-top-card",
      tags$p(class = "appv2-title", "Quality Control"),
      div(
        class = "appv2-top-row",
        div(
          class = "appv2-top-control",
          selectInput(
            "quality_control_main_view_v2",
            "Select view",
            choices = c(
              "Protein Groups",
              "Protein Intensities",
              "Sample Correlations",
              "Coefficient of Variation"
            )
          )
        ),
        div(
          class = "appv2-top-copy",
          div(
            style = "max-width: 260px;",
            sliderInput("qc_plot_width_v2", "Plot width", min = 4, max = 14, value = 6, step = 1),
            sliderInput("qc_plot_height_v2", "Plot height", min = 3, max = 10, value = 4, step = 1)
          )
        )
      )
    ),
    div(
      class = "appv2-main-body",
      div(
        class = "appv2-panel appv2-plot",
        uiOutput("quality_control_plot_ui_v2"),
        div(
          class = "appv2-download-row",
          downloadButton("download_qc_plot_v2", "Download PDF"),
          downloadButton("download_qc_table_v2", "Download TSV")
        )
      )
    )
  )
)

clustering_main_ui <- div(
  class = "appv2-panel appv2-plot",
  div(
    class = "appv2-view-shell",
    div(
      class = "appv2-panel appv2-top-card",
      tags$p(class = "appv2-title", "Clustering / Dimensionality Reduction"),
      div(
        class = "appv2-top-row",
        div(
          class = "appv2-top-control",
          selectInput(
            "clustering_main_view_v2",
            "Select view",
            choices = c("Hierarchical Clustering", "PCA", "UMAP")
          )
        ),
        div(
          class = "appv2-top-copy",
          div(
            style = "max-width: 260px;",
            sliderInput("clustering_plot_width_v2", "Plot width", min = 4, max = 14, value = 6, step = 1),
            sliderInput("clustering_plot_height_v2", "Plot height", min = 3, max = 10, value = 4, step = 1)
          )
        )
      )
    ),
    div(
      class = "appv2-main-body",
      div(
        class = "appv2-panel appv2-plot",
        uiOutput("clustering_plot_ui_v2"),
        div(
          class = "appv2-download-row",
          downloadButton("download_clustering_plot_v2", "Download PDF"),
          downloadButton("download_clustering_table_v2", "Download TSV")
        )
      )
    )
  )
)

abundance_profiling_main_ui <- div(
  class = "appv2-panel appv2-plot",
  div(
    class = "appv2-view-shell",
    div(
      class = "appv2-panel appv2-top-card",
      tags$p(class = "appv2-title", "Abundance Profiling"),
      div(
        class = "appv2-top-row",
        div(
          class = "appv2-top-control",
          selectInput(
            "abundance_main_view_v2",
            "Select view",
            choices = c("Barplot", "Heatmap")
          )
        ),
        div(
          class = "appv2-top-copy",
          div(
            style = "max-width: 260px;",
            sliderInput("abundance_plot_width_v2", "Plot width", min = 4, max = 14, value = 6, step = 1),
            sliderInput("abundance_plot_height_v2", "Plot height", min = 3, max = 10, value = 4, step = 1)
          )
        )
      )
    ),
    div(
      class = "appv2-main-body",
      div(
        class = "appv2-panel appv2-plot",
        uiOutput("abundance_plot_ui_v2"),
        div(
          class = "appv2-download-row",
          downloadButton("download_abundance_plot_v2", "Download PDF"),
          downloadButton("download_abundance_table_v2", "Download TSV")
        )
      )
    )
  )
)

differential_intensity_main_ui <- div(
  class = "appv2-panel appv2-plot",
  div(
    class = "appv2-view-shell",
    div(
      class = "appv2-panel appv2-top-card",
      tags$p(class = "appv2-title", "Differential Intensity"),
      div(
        class = "appv2-top-row",
        div(
          class = "appv2-top-control",
          selectInput(
            "differential_main_view_v2",
            "Select view",
            choices = c("Volcano Plot", "Pathway Enrichment")
          )
        ),
        div(
          class = "appv2-top-copy",
          div(
            style = "max-width: 260px;",
            sliderInput("differential_plot_width_v2", "Plot width", min = 4, max = 14, value = 6, step = 1),
            sliderInput("differential_plot_height_v2", "Plot height", min = 3, max = 10, value = 4, step = 1)
          )
        )
      )
    ),
    div(class = "appv2-main-body", uiOutput("differential_main_body_v2"))
  )
)

differential_params_ui <- div(
  class = "appv2-params-stack",
  radioButtons("de_mode_v2", "Outcome type", choices = c("Categorical" = "binary", "Continuous" = "continuous"), selected = "binary", inline = TRUE),
  selectInput("de_condition_v2", "Compare groups using", choices = character(0)),
  uiOutput("de_groups_v2"),
  uiOutput("de_covariates_v2"),
  selectInput("label_col_v2", "Protein label column", choices = character(0)),
  fileInput("gene_labels_v2", "Upload labeled genes"),
  uiOutput("logfc_v2"),
  checkboxInput("use_adj_pval_v2", "Use adjusted P-value (FDR)", value = TRUE),
  uiOutput("pvalue_v2"),
  conditionalPanel(
    condition = "input.differential_main_view_v2 == 'Pathway Enrichment'",
    div(
      class = "appv2-subsection",
      div(class = "appv2-subsection-title", "Pathway Enrichment"),
      selectInput("gene_col_v2", "Gene symbol column", choices = character(0)),
      selectInput("organism_v2", "Organism", choices = names(organism_map), selected = "Human"),
      selectInput("pathway_source_v2", "Pathway source", choices = c("GO", "Custom")),
      conditionalPanel(
        condition = "input.pathway_source_v2 == 'GO'",
        selectInput("go_ontology_v2", "GO ontology", choices = c("BP", "CC", "MF"), selected = "BP")
      ),
      conditionalPanel(
        condition = "input.pathway_source_v2 != 'GO'",
        fileInput("ontology_file_v2", "Upload ontology file")
      ),
      numericInput("enrich_pval_v2", "Enrichment p-value cutoff", value = 0.05, min = 0)
    )
  )
)

ui <- fluidPage(
  tags$head(
    tags$title("ProtPipe2 App V2"),
    tags$style(HTML("
      html, body {
        min-height: 100%;
        margin: 0;
        overflow-y: auto;
        background: linear-gradient(180deg, #f5f1e8 0%, #eef3f7 100%);
        font-family: 'Avenir Next', Avenir, 'Segoe UI', sans-serif;
        color: #1f2933;
      }
      .container-fluid {
        min-height: 100vh;
        padding: 16px 18px 18px 18px;
        max-width: none;
        width: 100%;
      }
      .appv2-shell {
        min-height: calc(100vh - 34px);
        display: flex;
        flex-direction: column;
        gap: 14px;
      }
      .appv2-header {
        display: flex;
        align-items: center;
        justify-content: space-between;
        gap: 18px;
        padding: 14px 18px;
        border-radius: 22px;
        background: rgba(255, 252, 247, 0.88);
        box-shadow: 0 16px 40px rgba(57, 72, 89, 0.08);
      }
      .appv2-brand h1 {
        margin: 0;
        font-size: 28px;
        font-weight: 800;
        letter-spacing: 0.02em;
      }
      .appv2-brand p {
        margin: 2px 0 0 0;
        color: #566373;
        font-size: 13px;
      }
      .appv2-tabs {
        display: grid;
        grid-template-columns: repeat(7, minmax(0, 1fr));
        gap: 10px;
        flex: 1;
      }
      .appv2-tab {
        width: 100%;
        min-height: 58px;
        border-radius: 16px;
        border: 1px solid #d7dfeb;
        background: #ffffff;
        color: #314154;
        font-weight: 700;
        font-size: 13px;
        line-height: 1.2;
        white-space: normal;
        padding: 10px 8px;
      }
      .appv2-tab.active-tab,
      .appv2-tab.active-tab:hover {
        background: #1f6f78;
        border-color: #1f6f78;
        color: #ffffff;
      }
      .appv2-main {
        flex: 1;
        min-height: 700px;
      }
      .appv2-panel {
        border-radius: 24px;
        background: rgba(255, 255, 255, 0.86);
        box-shadow: 0 18px 42px rgba(31, 41, 51, 0.09);
        border: 1px solid rgba(214, 220, 228, 0.9);
      }
      .appv2-col {
        min-height: 700px;
      }
      .appv2-sidebar {
        padding: 16px;
        min-height: 100%;
        overflow: auto;
      }
      .appv2-card {
        border-radius: 18px;
        background: #fbfcfd;
        border: 1px solid #e2e8f0;
        padding: 16px;
      }
      .appv2-card h3,
      .appv2-section-title {
        margin: 0 0 8px 0;
        font-size: 15px;
        font-weight: 800;
        text-transform: uppercase;
        letter-spacing: 0.05em;
        color: #4d6278;
      }
      .appv2-params-card .shiny-input-container {
        margin-bottom: 12px;
      }
      .appv2-params-card-tall {
        min-height: 800px;
      }
      .appv2-params-stack {
        display: flex;
        flex-direction: column;
        gap: 10px;
      }
      .appv2-subsection {
        border-top: 1px solid #dde4ec;
        padding-top: 12px;
      }
      .appv2-subsection:first-child {
        border-top: 0;
        padding-top: 0;
      }
      .appv2-subsection-title {
        margin: 0 0 10px 0;
        font-size: 13px;
        font-weight: 800;
        text-transform: uppercase;
        letter-spacing: 0.05em;
        color: #4d6278;
      }
      .appv2-subsection-grid {
        display: grid;
        grid-template-columns: 1fr 1fr;
        gap: 12px;
      }
      .appv2-view-shell {
        min-height: 100%;
        display: flex;
        flex-direction: column;
        gap: 14px;
      }
      .appv2-top-card,
      .appv2-plot {
        padding: 16px;
      }
      .appv2-top-row {
        display: flex;
        align-items: flex-start;
        justify-content: space-between;
        gap: 16px;
      }
      .appv2-top-copy {
        flex: 1 1 auto;
        min-width: 0;
      }
      .appv2-top-control {
        flex: 0 0 240px;
      }
      .appv2-top-control .shiny-input-container {
        margin-bottom: 0;
      }
      .appv2-main-body {
        flex: 1 1 auto;
        min-height: 0;
      }
      .appv2-plot {
        min-height: 0;
        overflow: auto;
      }
      .appv2-plot-card {
        min-height: 420px;
        display: flex;
        flex-direction: column;
        justify-content: center;
        align-items: center;
        text-align: center;
        border-radius: 18px;
        background:
          radial-gradient(circle at top left, rgba(31, 111, 120, 0.16), transparent 35%),
          linear-gradient(135deg, #f8fafc 0%, #eef6f7 100%);
        border: 1px dashed #90a4ae;
        padding: 24px;
      }
      .appv2-table-wrap {
        width: 100%;
        min-height: 500px;
        max-height: calc(100vh - 340px);
        overflow: auto;
        border-radius: 16px;
        border: 1px solid #d6dee6;
        background: #ffffff;
        padding: 8px 10px;
      }
      .appv2-table-wrap table {
        margin-bottom: 0;
        font-size: 13px;
      }
      .appv2-download-row {
        display: flex;
        gap: 12px;
        margin-top: 12px;
      }
      .appv2-actionbtn {
        min-height: 48px;
        min-width: 220px;
        font-size: 15px;
        font-weight: 700;
      }
      .appv2-title {
        margin: 0;
        font-size: 26px;
        font-weight: 800;
      }
      .appv2-subtitle {
        margin: 6px 0 0 0;
        color: #607082;
        font-size: 15px;
      }
      .appv2-footer {
        display: flex;
        justify-content: space-between;
        gap: 12px;
        padding-top: 4px;
      }
      .appv2-navbtn {
        min-width: 120px;
        min-height: 48px;
        border-radius: 14px;
        font-size: 15px;
        font-weight: 700;
      }
      @media (max-width: 1180px) {
        .appv2-header {
          flex-direction: column;
          align-items: stretch;
        }
        .appv2-tabs {
          grid-template-columns: repeat(4, minmax(0, 1fr));
        }
        .appv2-top-row {
          flex-wrap: wrap;
        }
        .appv2-top-control {
          width: 100%;
        }
      }
    "))
  ),
  div(style = "display:none;", textInput("current_page", label = NULL, value = page_ids[[1]])),
  div(
    class = "appv2-shell",
    div(
      class = "appv2-header",
      div(
        class = "appv2-brand",
        h1("ProtPipe2")
      ),
      uiOutput("top_tabs")
    ),
    div(
      class = "appv2-main",

      ## Input -----------------------------------------------------------------
      conditionalPanel(
        condition = "input.current_page == 'input'",
        appv2_page_layout(
          param_ui = div(
            class = "appv2-params-stack",
            fileInput("intensity_matrix", "Upload proteomics data"),
            fileInput("sample_metadata", "Upload sample metadata"),
            downloadButton("download_example_proteomics_v2", "Example proteomics dataset"),
            uiOutput("column_range_ui_v2")
          ),
          main_ui = input_main_ui
        )
      ),

      ## Quality Control -------------------------------------------------------
      conditionalPanel(
        condition = "input.current_page == 'quality_control'",
        appv2_page_layout(
          param_ui = uiOutput("quality_control_params_v2"),
          main_ui = quality_control_main_ui,
          param_card_class = "appv2-params-card-tall"
        )
      ),

      ## Pre-Processing --------------------------------------------------------
      conditionalPanel(
        condition = "input.current_page == 'pre_processing'",
        appv2_page_layout(
          param_ui = div(
            class = "appv2-params-stack",
            div(
              class = "appv2-subsection",
              div(class = "appv2-subsection-title", "1. Minimum Intensity Filtering"),
              uiOutput("preprocessing_min_filter_ui_v2")
            ),
            div(
              class = "appv2-subsection",
              div(class = "appv2-subsection-title", "2. Outlier Removal"),
              div(
                class = "appv2-subsection-grid",
                div(
                  checkboxInput("remove_outliers_v2", "Remove outlier samples", value = FALSE),
                  numericInput("outlier_sds_v2", "Remove samples with protein groups outside n standard deviations", value = 3)
                ),
                div(
                  checkboxInput("remove_sparse_proteins_v2", "Remove sparse proteins", value = FALSE),
                  numericInput("sparse_protein_percent_v2", "Remove proteins present in less than n% of samples", value = 30)
                )
              )
            ),
            div(
              class = "appv2-subsection",
              div(class = "appv2-subsection-title", "3. Transformation"),
              checkboxInput("log2_transform_v2", "Log2 transform", value = FALSE)
            ),
            div(
              class = "appv2-subsection",
              div(class = "appv2-subsection-title", "4. Normalization"),
              div(
                class = "appv2-subsection-grid",
                div(checkboxInput("normalize_v2", "Normalize", value = FALSE)),
                div(selectInput("normalize_method_v2", "Normalization method", choices = c("mean", "median"), selected = "median"))
              )
            ),
            div(
              class = "appv2-subsection",
              div(class = "appv2-subsection-title", "5. Imputation"),
              div(
                class = "appv2-subsection-grid",
                div(
                  checkboxInput("impute_v2", "Impute missing values", value = TRUE),
                  selectInput("imputation_method_v2", "Imputation method", choices = c("fixed value", "minimum", "left-shifted distribution"), selected = "fixed value")
                ),
                div(
                  conditionalPanel(
                    condition = "input.imputation_method_v2 == 'fixed value'",
                    numericInput("impute_fixed_value_v2", "Value", value = 0)
                  ),
                  conditionalPanel(
                    condition = "input.imputation_method_v2 == 'minimum'",
                    numericInput("impute_min_value_v2", "Scale minimum by", value = 1)
                  ),
                  conditionalPanel(
                    condition = "input.imputation_method_v2 == 'left-shifted distribution'",
                    numericInput("impute_left_dist_shift_v2", "Shift mean of distribution by n standard deviations", value = 1.8),
                    numericInput("impute_left_dist_scale_v2", "Scale standard deviation of distribution by", value = 0.3)
                  )
                )
              )
            ),
            div(
              class = "appv2-subsection",
              div(class = "appv2-subsection-title", "6. Batch Correction"),
              uiOutput("batch_correct_section_v2")
            ),
            div(
              class = "appv2-subsection",
              div(class = "appv2-subsection-title", "Downloads"),
              downloadButton("download_preprocessing_data_v2", "Download TSV"),
              downloadButton("download_preprocessing_report_v2", "Download Markdown")
            )
          ),
          main_ui = NULL,
          show_main = FALSE
        )
      ),

      ## Clustering ------------------------------------------------------------
      conditionalPanel(
        condition = "input.current_page == 'clustering'",
        appv2_page_layout(
          param_ui = uiOutput("clustering_params_v2"),
          main_ui = clustering_main_ui,
          param_card_class = "appv2-params-card-tall"
        )
      ),

      ## Differential Intensity -----------------------------------------------
      conditionalPanel(
        condition = "input.current_page == 'differential_intensity'",
        appv2_page_layout(
          param_ui = differential_params_ui,
          main_ui = differential_intensity_main_ui
        )
      ),

      ## Abundance Profiling ---------------------------------------------------
      conditionalPanel(
        condition = "input.current_page == 'abundance_profiling'",
        appv2_page_layout(
          param_ui = uiOutput("abundance_params_v2"),
          main_ui = abundance_profiling_main_ui,
          param_card_class = "appv2-params-card-tall"
        )
      ),

      ## Help ------------------------------------------------------------------
      conditionalPanel(
        condition = "input.current_page == 'help'",
        fluidRow(
          style = "height: 100%; margin-left: 0; margin-right: 0;",
          column(
            width = 12,
            class = "appv2-col",
            div(
              class = "appv2-panel appv2-plot",
              div(
                class = "appv2-view-shell",
                div(
                  class = "appv2-panel appv2-top-card",
                  tags$p(class = "appv2-title", "Help")
                ),
                div(
                  class = "appv2-main-body",
                  div(
                    class = "appv2-panel appv2-plot",
                    includeMarkdown("help.md")
                  )
                )
              )
            )
          )
        )
      )
    ),
    div(
      class = "appv2-footer",
      actionButton("back_page", "Back", class = "btn-default appv2-navbtn"),
      actionButton("next_page", "Next", class = "btn-primary appv2-navbtn")
    )
  )
)
