options(show.error.locations = TRUE)

p <- function(...) {
  cat(as.character(...), sep = "\n")
}

latch_error <- function(...) {
  output <- paste(
    "__LATCH_ERROR_START__",
    as.character(...),
    "__LATCH_ERROR_END__",
    sep = "\n"
  )
  cat(output)
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

read_tabular <- function(path) {
  tryCatch(
    {
      read_excel(path)
    },
    error = function(cond) {
      read_delim(path, trim_ws = TRUE)
    }
  )
}

p("Loading contrast data")
args <- commandArgs(trailingOnly = TRUE)

res <- read_tabular(args[[1]])

num_pathways <- strtoi(args[[2]])
msig_species_id <- args[[3]]
go_species_id <- args[[4]]
kegg_species_id <- args[[5]]

p("Loading gene annotations")
suppressMessages(suppressWarnings(library("org.Hs.eg.db", character.only = TRUE)))
suppressMessages(suppressWarnings(library("org.Mm.eg.db", character.only = TRUE)))
suppressMessages(suppressWarnings(library("org.Sc.sgd.db", character.only = TRUE)))
suppressMessages(suppressWarnings(library("org.Rn.eg.db", character.only = TRUE)))

msig_db <- msigdbr(species = msig_species_id) %>% dplyr::select(gs_name, entrez_gene)

p("Preparing data")


get_names_package <- function(package_string) {
  if (package_string == "org.Hs.eg.db") {
    return(org.Hs.eg.db)
  } else if (package_string == "org.Mm.eg.db") {
    return(org.Mm.eg.db)
  } else if (package_string == "org.Sc.sgd.db") {
    return(org.Sc.sgd.db)
  } else if (package_string == "org.Rn.eg.db") {
    return(org.Rn.eg.db)
  } else {
    return(NULL)
  }
}

dsD <- tryCatch(
  {
    res %>%
      mutate(gene_name = .data[[vec_as_names("", repair = "unique")]]) %>%
      dplyr::select(gene_name, log2FoldChange) %>%
      na.omit() %>%
      mutate(gene_name = mapIds(get_names_package(go_species_id), keys = gene_name, keytype = "ALIAS", column = "ENTREZID", multiVals = "first")) %>%
      na.omit() %>%
      distinct(gene_name, .keep_all = TRUE) %>%
      column_to_rownames("gene_name")
  },
  error = function(cond) {
    return(
      res %>%
        mutate(gene_name = .data[[vec_as_names("", repair = "unique")]]) %>%
        dplyr::select(gene_name, log2FoldChange) %>%
        na.omit() %>%
        mutate(gene_name = gsub("\\.\\d+", "", gene_name)) %>%
        mutate(gene_name = mapIds(get_names_package(go_species_id), keys = gene_name, keytype = "ENSEMBL", column = "ENTREZID", multiVals = "first")) %>%
        na.omit() %>%
        distinct(gene_name, .keep_all = TRUE) %>%
        column_to_rownames("gene_name")
    )
  }
)

ds <- dsD$log2FoldChange
names(ds) <- rownames(dsD)
ds <- ds %>% sort(decreasing = TRUE)

p("Running GO GSE")
start <- Sys.time()
go <- gseGO(ds, ont = "ALL", go_species_id, keyType = "ENTREZID")
print(Sys.time() - start)

head(go)

if (nrow(go) > 0) {
  p("  Plotting")
  dir.create("./res/Gene Ontology", showWarnings = FALSE, recursive = TRUE)
  png(file = "./res/Gene Ontology/Dot Plot.png", width = 960, height = 900)
  print(dotplot(go, showCategory = 20, split = ".sign") + facet_grid(. ~ .sign))
  dev.off()

  png(file = "./res/Gene Ontology/Ridge Plot.png", width = 960, height = 900)
  print(ridgeplot(go) + labs(x = "enrichment distribution"))
  dev.off()
} else {
  warn(paste(
    "No statistically significant enriched gene sets (with cutoff p ≤ 0.05) found after ",
    "running gene set enrichment analysis (GSEA) on Gene Ontology.",
    sep = ""
  ))
}

p("Running GO enrichment")
start <- Sys.time()
enrich_go <- enrichGO(names(ds), go_species_id, ont = "ALL", keyType = "ENTREZID")
enrich_go <- pairwise_termsim(enrich_go)
print(Sys.time() - start)

head(enrich_go)

if (nrow(enrich_go) > 0) {
  p("  Plotting")
  png(file = "./res/Gene Ontology/Enrichment Map.png", width = 960, height = 900)
  print(emapplot(enrich_go))
  dev.off()
} else {
  warn(paste(
    "No statistically significant enrichment tests (with cutoff p ≤ 0.05) found after ",
    "running enrichment analysis on Gene Ontology.",
    sep = ""
  ))
}

p("Running KEGG GSE")
start <- Sys.time()
kks <- gseKEGG(ds, kegg_species_id)
print(Sys.time() - start)

if (nrow(kks) > 0) {
  p("  Plotting")
  dir.create("./res/KEGG", showWarnings = FALSE, recursive = TRUE)
  png(file = "./res/KEGG/Dot Plot.png", width = 960, height = 900)
  print(dotplot(kks, showCategory = 20, split = ".sign") + facet_grid(. ~ .sign))
  dev.off()

  png(file = "./res/KEGG/Ridge Plot.png", width = 960, height = 900)
  print(ridgeplot(kks) + labs(x = "enrichment distribution"))
  dev.off()

  entrezIDsToGeneNames <- function(entrezIDs) {
    return(mapIds(get_names_package(go_species_id), entrezIDs, "SYMBOL", "ENTREZID") %>% paste(collapse = " "))
  }

  pathways <- kks@result %>%
    slice_max(order_by = enrichmentScore, n = num_pathways) %>%
    mutate(coreEntrezIDs = core_enrichment) %>%
    mutate(entrezList = strsplit(core_enrichment, "/"), .keep = "unused") %>%
    mutate(coreEnrichedGenes = unlist(lapply(entrezList, entrezIDsToGeneNames)), .keep = "unused")

  write.csv(pathways, "./res/KEGG/table.csv", row.names = FALSE)

  dir.create("./tempres", showWarnings = FALSE, recursive = TRUE)
  genesets_path <- "./tempres/genesets.txt"
  geneSets <- kks@geneSets
  write("PATHWAYIDS", genesets_path)
  lapply(names(geneSets), write, genesets_path, append = TRUE, ncolumns = 100000)
  write("ENTREZIDS", genesets_path, append = TRUE)
  lapply(geneSets, write, genesets_path, append = TRUE, ncolumns = 100000)
  write("NAMES", genesets_path, append = TRUE)
  geneNames <- lapply(geneSets, entrezIDsToGeneNames)
  lapply(geneNames, write, genesets_path, append = TRUE, ncolumns = 100000)

  for (pathwayID in pathways$ID) {
    p(paste("  Running pathview on", pathwayID))
    pathview(gene.data = ds, pathway.id = pathwayID, species = msig_species_id)
  }
} else {
  warn(paste(
    "No statistically significant enriched gene sets (with cutoff p ≤ 0.05) found after ",
    "running gene set enrichment analysis (GSEA) on KEGG. ",
    "Only static plots will be available in the report",
    sep = ""
  ))
}

p("Running KEGG enrichment")
start <- Sys.time()
enrich_kegg <- enrichGO(names(ds), go_species_id, ont = "ALL", keyType = "ENTREZID")
enrich_kegg <- pairwise_termsim(enrich_kegg)
print(Sys.time() - start)

head(enrich_kegg)

if (nrow(enrich_kegg) > 0) {
  p("  Plotting")
  png(file = "./res/KEGG/Enrichment Map.png", width = 960, height = 900)
  print(emapplot(enrich_kegg))
  dev.off()
} else {
  warn(paste(
    "No statistically significant enrichment tests (with cutoff p ≤ 0.05) found after ",
    "running enrichment analysis on KEGG.",
    sep = ""
  ))
}

p("Running MSig")
start <- Sys.time()
msig <- GSEA(ds, TERM2GENE = msig_db)
print(Sys.time() - start)

head(msig)

p("  Plotting")
dir.create("./res/MSig", showWarnings = FALSE, recursive = TRUE)
png(file = "./res/MSig/Dot Plot.png", width = 960, height = 900)
print(dotplot(msig, showCategory = 20, split = ".sign") + facet_grid(. ~ .sign))
dev.off()

png(file = "./res/MSig/Ridge Plot.png", width = 960, height = 900)
print(ridgeplot(msig) + labs(x = "enrichment distribution"))
dev.off()

p("Running MSig enrichment")
start <- Sys.time()
enrich_msig <- enricher(names(ds), TERM2GENE = msig_db)
enrich_msig <- pairwise_termsim(enrich_msig)
print(Sys.time() - start)

head(enrich_msig)

p("  Plotting")
png(file = "./res/MSig/Enrichment Map.png", width = 960, height = 900)
print(emapplot(enrich_msig))
dev.off()
