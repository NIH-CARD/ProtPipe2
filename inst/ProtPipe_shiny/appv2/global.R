options(shiny.sanitize.errors = TRUE)
options(shiny.maxRequestSize = 5000 * 1024^2)

source("helpers.R")

library(SummarizedExperiment)
library(shiny)
library(ProtPipe)
library(clusterProfiler)
library(Hmisc)
library(markdown)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(org.Rn.eg.db)
library(org.Dm.eg.db)
library(org.Ce.eg.db)

organism_map <- list(
  "Human" = list(OrgDb = org.Hs.eg.db, kegg = "hsa"),
  "Mouse" = list(OrgDb = org.Mm.eg.db, kegg = "mmu"),
  "Rat" = list(OrgDb = org.Rn.eg.db, kegg = "rno"),
  "Fly" = list(OrgDb = org.Dm.eg.db, kegg = "dme"),
  "Nematode" = list(OrgDb = org.Ce.eg.db, kegg = "cel")
)
