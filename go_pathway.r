options(show.error.locations = TRUE)

p <- function (...) {
  cat(as.character(...), sep = "\n")
}

warn <- function(...) {
  output <- paste(
    "__LATCH_WARNING_START__",
    as.character(...),
    "__LATCH_WARNING_END__",
    sep = "\n"
  )
  cat(output)
}

p("Importing")
suppressMessages(suppressWarnings(library(vctrs)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(tibble)))
suppressMessages(suppressWarnings(library(readr)))
suppressMessages(suppressWarnings(library(readxl)))
suppressMessages(suppressWarnings(library(stringr)))

suppressMessages(suppressWarnings(library(purrr)))

suppressMessages(suppressWarnings(library(clusterProfiler)))
suppressMessages(suppressWarnings(library(DOSE)))
suppressMessages(suppressWarnings(library(ggridges)))
suppressMessages(suppressWarnings(library(enrichplot)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(msigdbr)))
suppressMessages(suppressWarnings(library(pathview)))

read_tabular = function (path) {
  tryCatch({
    read_excel(path)
  }, error = function (cond) {
    read_delim(path, trim_ws = TRUE)
  })
}

p("Loading gene annotations")
organism = "org.Hs.eg.db"
suppressMessages(suppressWarnings(library(organism, character.only=TRUE)))
msig_db <- msigdbr(species="Homo sapiens") %>% dplyr::select(gs_name, entrez_gene)

p("Loading contrast data")
args <- commandArgs(trailingOnly = TRUE)

res <- read_tabular(args[[1]])
num_pathways <- strtoi(args[[2]])

p("Preparing data")
dsD <- res %>%
  mutate(gene_name = .data[[vec_as_names("", repair = "unique")]]) %>%
  dplyr::select(gene_name, log2FoldChange) %>%
  na.omit %>%
  mutate(gene_name = mapIds(org.Hs.eg.db, keys=gene_name, keytype="ALIAS", column="ENTREZID")) %>%
  na.omit %>%
  distinct(gene_name, .keep_all=TRUE) %>%
  column_to_rownames("gene_name")
ds <- dsD$log2FoldChange
names(ds) <- rownames(dsD)
ds <- ds %>% sort(decreasing=TRUE)

p("Running MSig")
start <- Sys.time()
msig <- GSEA(ds, TERM2GENE=msig_db)
print(Sys.time() - start)

head(msig)
p("Running GO")
start <- Sys.time()
go <- gseGO(ds, ont="ALL", organism, keyType="ENTREZID")
print(Sys.time() - start)

head(go)

if (nrow(go) > 0) {
  p("  Plotting")
  dir.create("/root/res/Gene Ontology", showWarnings = FALSE, recursive = TRUE)
  png(file="/root/res/Gene Ontology/Dot Plot.png", width=960, height=900)
  print(dotplot(go, showCategory=20, split=".sign") + facet_grid(.~.sign))
  dev.off()
  
  png(file="/root/res/Gene Ontology/Ridge Plot.png", width=960, height=900)
  print(ridgeplot(go) + labs(x = "enrichment distribution"))
  dev.off()
} else {
  warn(paste(
    "No statistically significant enriched gene sets (with cutoff p ≤ 0.05) found after ",
    "running gene set enrichment analysis (GSEA) on Gene Ontology.",
    sep = ""
  ))
}

p("Running KEGG")
start <- Sys.time()
kks <- gseKEGG(ds, "hsa")
print(Sys.time() - start)


if (nrow(kks) > 0) {
  p("  Plotting")
  dir.create("/root/res/KEGG", showWarnings = FALSE, recursive = TRUE)
  png(file="/root/res/KEGG/Dot Plot.png", width=960, height=900)
  print(dotplot(kks, showCategory=20, split=".sign") + facet_grid(.~.sign))
  dev.off()
  
  png(file="/root/res/KEGG/Ridge Plot.png", width=960, height=900)
  print(ridgeplot(kks) + labs(x = "enrichment distribution"))
  dev.off()
  
  entrezIDsToGeneNames <- function (entrezIDs) {
    return(mapIds(org.Hs.eg.db, entrezIDs, "SYMBOL", "ENTREZID") %>% paste(collapse = " "))
  }
  
  pathways <- kks@result %>%
    slice_max(order_by = enrichmentScore, n = num_pathways) %>%
    mutate(coreEntrezIDs = core_enrichment) %>%
    mutate(entrezList = strsplit(core_enrichment, '/'), .keep = "unused") %>%
    mutate(coreEnrichedGenes = unlist(lapply(entrezList, entrezIDsToGeneNames)), .keep = "unused")

  write.csv(pathways, "/root/res/KEGG/table.csv", row.names = FALSE)

  dir.create("/root/tempres", showWarnings = FALSE, recursive = TRUE) 
  genesets_path <- "/root/tempres/genesets.txt"
  geneSets <- kks@geneSets
  write("PATHWAYIDS", genesets_path)
  lapply(names(geneSets), write, genesets_path, append=TRUE, ncolumns=100000)
  write("ENTREZIDS", genesets_path, append=TRUE)
  lapply(geneSets, write, genesets_path, append=TRUE, ncolumns=100000)
  write("NAMES", genesets_path, append=TRUE)
  geneNames <- lapply(geneSets, entrezIDsToGeneNames)
  lapply(geneNames, write, genesets_path, append=TRUE, ncolumns=100000)
  
  for (pathwayID in pathways$ID) {
    p(paste("  Running pathview on", pathwayID))
    pathview(gene.data=ds, pathway.id=pathwayID, species="Homo sapiens")
  }
} else {
  warn(paste(
    "No statistically significant enriched gene sets (with cutoff p ≤ 0.05) found after ",
    "running gene set enrichment analysis (GSEA) on KEGG.",
    sep = ""
  ))
}


p("  Plotting")
dir.create("/root/res/MSig", showWarnings = FALSE, recursive = TRUE)
png(file="/root/res/MSig/Dot Plot.png", width=960, height=900)
print(dotplot(msig, showCategory=20, split=".sign") + facet_grid(.~.sign))
dev.off()

png(file="/root/res/MSig/Ridge Plot.png", width=960, height=900)
print(ridgeplot(msig) + labs(x = "enrichment distribution"))
dev.off()
