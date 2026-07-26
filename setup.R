#!/usr/bin/env Rscript
#
# setup.R -- install every package needed to develop, test, check, and run ProtPipe2.
#
# Usage (from the project root):
#
#   Rscript setup.R              # install whatever is missing
#   Rscript setup.R --force      # reinstall everything, even if already present
#   Rscript setup.R --check      # report status only, install nothing
#
# The dependency list is read from DESCRIPTION (Imports + Suggests) so this script
# stays in sync with the package. A short list of development-only tools
# (devtools, roxygen2, pkgdown, ...) is added on top.

args <- commandArgs(trailingOnly = TRUE)
force_reinstall <- "--force" %in% args
check_only <- "--check" %in% args

# ---------------------------------------------------------------------------
# Repositories and install options
# ---------------------------------------------------------------------------

repos <- getOption("repos")
if (is.null(repos[["CRAN"]]) || is.na(repos[["CRAN"]]) || repos[["CRAN"]] == "@CRAN@") {
  repos[["CRAN"]] <- "https://cloud.r-project.org"
}
options(
  repos = repos,
  Ncpus = max(1L, parallel::detectCores() - 1L),
  warn = 1
)

# On macOS and Windows, prefer precompiled binaries. Without this, R silently
# falls back to building from source whenever the source version is newer than
# the binary, which fails for packages needing a Fortran toolchain
# (RcppArmadillo, lme4, ...) unless gfortran is installed.
if (.Platform$OS.type == "windows" || Sys.info()[["sysname"]] == "Darwin") {
  options(
    pkgType = "binary",
    install.packages.check.source = "no",
    install.packages.compile.from.source = "never"
  )
}

message("R version:  ", getRversion())
message("Library:    ", .libPaths()[1])
message("CRAN repo:  ", getOption("repos")[["CRAN"]])

# ---------------------------------------------------------------------------
# BiocManager -- needed to resolve the Bioconductor dependencies
# ---------------------------------------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  message("\nInstalling BiocManager ...")
  install.packages("BiocManager")
}
message("Bioconductor: ", as.character(BiocManager::version()))

# ---------------------------------------------------------------------------
# Dependency list
# ---------------------------------------------------------------------------

desc_path <- file.path(getwd(), "DESCRIPTION")
if (!file.exists(desc_path)) {
  stop("DESCRIPTION not found. Run this script from the ProtPipe2 project root.", call. = FALSE)
}

parse_deps <- function(path, fields = c("Imports", "Suggests", "Depends")) {
  dcf <- read.dcf(path, fields = fields)
  entries <- unlist(strsplit(paste(na.omit(as.vector(dcf)), collapse = ","), ","))
  entries <- trimws(sub("\\(.*\\)", "", entries))     # drop version constraints
  entries <- entries[nzchar(entries)]
  unique(entries)
}

# Packages shipped with R itself -- never install these.
base_pkgs <- rownames(installed.packages(priority = "base"))

# Development / CI tooling that is not a package dependency.
dev_pkgs <- c(
  "devtools",     # load_all(), test(), check(), install()
  "roxygen2",     # document()
  "pkgdown",      # website build (see .github/workflows/pkgdown.yaml)
  "rcmdcheck",    # R CMD check driver used in CI
  "remotes",
  "testthat",
  "knitr",
  "rmarkdown"
)

# Used in R/heatmap.R and by the tidyverse-style code paths, but not currently
# declared in DESCRIPTION.
undeclared_pkgs <- c("tibble")

pkgs <- setdiff(
  unique(c(parse_deps(desc_path), dev_pkgs, undeclared_pkgs)),
  c(base_pkgs, "R")
)
pkgs <- sort(pkgs)

installed <- rownames(installed.packages())
missing <- setdiff(pkgs, installed)

message("\n", length(pkgs), " packages required, ", length(missing), " missing.")
if (length(missing)) {
  message("Missing: ", paste(missing, collapse = ", "))
}

if (check_only) {
  quit(status = if (length(missing)) 1L else 0L)
}

# ---------------------------------------------------------------------------
# Install
# ---------------------------------------------------------------------------

to_install <- if (force_reinstall) pkgs else missing

if (length(to_install) == 0) {
  message("\nNothing to install -- all dependencies are already present.")
} else {
  message("\nInstalling ", length(to_install), " packages ...\n")
  # BiocManager::install() resolves both CRAN and Bioconductor packages, and
  # pins Bioconductor to the release matching this R version.
  BiocManager::install(
    to_install,
    update = FALSE,
    ask = FALSE,
    checkBuilt = FALSE
  )
}

# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

installed <- rownames(installed.packages())
still_missing <- setdiff(pkgs, installed)

message("\n", strrep("-", 60))
if (length(still_missing) == 0) {
  message("All ", length(pkgs), " dependencies are installed.")
  message("\nNext steps:")
  message("  devtools::load_all()   # load the package for development")
  message("  devtools::test()       # run the test suite")
} else {
  message("Installed: ", length(pkgs) - length(still_missing), "/", length(pkgs))
  message("STILL MISSING (", length(still_missing), "): ",
          paste(still_missing, collapse = ", "))
  message("\nThese usually fail because of missing system libraries. ",
          "Re-run with --force after resolving, or install them individually ",
          "with BiocManager::install(\"<pkg>\").")
}
message(strrep("-", 60))

quit(status = if (length(still_missing)) 1L else 0L)
