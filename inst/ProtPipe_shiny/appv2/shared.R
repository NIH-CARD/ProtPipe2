workflow_pages <- c(
  "Input",
  "Quality Control",
  "Pre-Processing",
  "Clustering / Dimensionality Reduction",
  "Differential Intensity",
  "Abundance Profiling",
  "Help"
)

page_ids <- c(
  "input",
  "quality_control",
  "pre_processing",
  "clustering",
  "differential_intensity",
  "abundance_profiling",
  "help"
)

page_descriptions <- c(
  "Upload the core analysis files and define how the dataset should enter the workflow.",
  "Inspect sample quality and data consistency before making analytical decisions.",
  "Configure filtering, transformation, normalization, imputation, and batch handling.",
  "View sample-level structure with clustering and dimensionality reduction outputs.",
  "Explore protein-level group differences and statistical significance.",
  "Inspect abundance trends for selected proteins across groups and samples.",
  "Find guidance for data formatting, workflow order, and interpretation."
)

top_tab_button <- function(id, label, active = FALSE) {
  actionButton(
    inputId = paste0("tab_", id),
    label = label,
    class = paste("appv2-tab", if (active) "active-tab" else "")
  )
}
