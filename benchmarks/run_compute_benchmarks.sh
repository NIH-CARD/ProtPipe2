#!/usr/bin/env bash
set -euo pipefail

ITERATIONS="${ITERATIONS:-3}"
OUT_DIR="benchmarks/results"
TIMESTAMP="$(date +%Y%m%d_%H%M%S)"
RAW_CSV="${OUT_DIR}/compute_benchmark_raw_${TIMESTAMP}.csv"
SUMMARY_CSV="${OUT_DIR}/compute_benchmark_summary_${TIMESTAMP}.csv"

mkdir -p "${OUT_DIR}"

benchmarks=(
  "se_creation"
  "qc_get_pg_counts"
  "qc_plot_pg_counts"
  "qc_plot_pg_intensities"
  "qc_get_cvs"
  "qc_plot_cvs"
  "qc_get_sample_correlation"
  "qc_plot_correlation_heatmap"
  "clustering_get_pcs"
  "clustering_plot_pcs"
  "clustering_plot_hierarchical"
  "clustering_get_umap"
  "clustering_plot_umap"
  "dea_day0_vs_day28"
  "dea_plot_volcano"
  "pathway_go_enrichment"
  "protein_barchart"
  "proteomics_heatmap"
  "total_workflow"
)

echo "benchmark,iteration,elapsed_sec,max_rss_bytes,max_rss_mib,max_rss_gib" > "${RAW_CSV}"

for benchmark in "${benchmarks[@]}"; do
  for ((i=1; i<=ITERATIONS; i++)); do
    log_file="$(mktemp)"
    if ! /usr/bin/time -l Rscript benchmarks/run_one_benchmark.R "${benchmark}" > "${log_file}" 2>&1; then
      cat "${log_file}" >&2
      rm -f "${log_file}"
      exit 1
    fi

    elapsed_sec="$(
      awk -F= '/^elapsed_sec=/ { print $2; exit }' "${log_file}"
    )"
    max_rss_bytes="$(
      awk '/maximum resident set size/ { print $1; exit }' "${log_file}"
    )"

    if [[ -z "${elapsed_sec}" || -z "${max_rss_bytes}" ]]; then
      cat "${log_file}" >&2
      rm -f "${log_file}"
      echo "Failed to parse timing output for ${benchmark}" >&2
      exit 1
    fi

    max_rss_mib="$(awk -v bytes="${max_rss_bytes}" 'BEGIN { printf "%.3f", bytes / (1024 * 1024) }')"
    max_rss_gib="$(awk -v bytes="${max_rss_bytes}" 'BEGIN { printf "%.3f", bytes / (1024 * 1024 * 1024) }')"

    echo "${benchmark},${i},${elapsed_sec},${max_rss_bytes},${max_rss_mib},${max_rss_gib}" >> "${RAW_CSV}"
    rm -f "${log_file}"
  done
done

Rscript -e '
df <- read.csv(commandArgs(trailingOnly = TRUE)[1], stringsAsFactors = FALSE)
split_df <- split(df, df$benchmark)
summary_df <- do.call(
  rbind,
  lapply(split_df, function(x) {
    data.frame(
      benchmark = x$benchmark[[1]],
      iterations = nrow(x),
      median_runtime_sec = stats::median(x$elapsed_sec),
      max_rss_bytes = max(x$max_rss_bytes),
      max_rss_mib = max(x$max_rss_mib),
      max_rss_gib = max(x$max_rss_gib),
      stringsAsFactors = FALSE
    )
  })
)
summary_df <- summary_df[order(summary_df$benchmark), ]
write.csv(summary_df, commandArgs(trailingOnly = TRUE)[2], row.names = FALSE)
' "${RAW_CSV}" "${SUMMARY_CSV}"

echo "Raw results: ${RAW_CSV}"
echo "Summary results: ${SUMMARY_CSV}"
