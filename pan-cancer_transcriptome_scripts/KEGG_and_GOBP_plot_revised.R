# ==============================================================================
# GO Biological Process and KEGG enrichment analysis
#
# This script:
#   1. reads an Ensembl gene list;
#   2. converts Ensembl IDs to Entrez IDs;
#   3. performs GO Biological Process and KEGG enrichment analyses;
#   4. calculates enrichment folds;
#   5. exports result tables and publication-ready figures.
#
# Update the configuration section before running the script.
# ==============================================================================

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(dplyr)
  library(stringr)
})

# ------------------------------------------------------------------------------
# 1. Configuration
# ------------------------------------------------------------------------------

analysisName <- "magoh"

geneFile <- "/data1/DYY/bambu/all_ts_gene.csv"
resultDir <- path.expand("~/bambu_data_kegg_go")
plotDir <- "/data1/DYY/bambu/plot/enrichment_plot"

geneIdType <- "ENSEMBL"
goOntology <- "BP"

pValueCutoff <- 0.5
qValueCutoff <- 0.5
topTermCount <- 10

goPlotWidth <- 10
goPlotHeight <- 10
keggPlotWidth <- 10
keggPlotHeight <- 10

dir.create(resultDir, recursive = TRUE, showWarnings = FALSE)
dir.create(plotDir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 2. Helper functions
# ------------------------------------------------------------------------------

checkRequiredFile <- function(filePath) {
  if (!file.exists(filePath)) {
    stop("Required input file was not found: ", filePath)
  }
}

readGeneIds <- function(filePath) {
  checkRequiredFile(filePath)

  geneTable <- read.table(
    file = filePath,
    header = FALSE,
    sep = "",
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )

  if (ncol(geneTable) < 1L) {
    stop("The gene file does not contain any columns: ", filePath)
  }

  geneIds <- unique(trimws(as.character(geneTable[[1L]])))
  geneIds <- geneIds[!is.na(geneIds) & nzchar(geneIds)]

  if (length(geneIds) == 0L) {
    stop("No valid gene identifiers were found in: ", filePath)
  }

  geneIds
}

parseRatio <- function(ratioVector) {
  ratioParts <- str_split_fixed(as.character(ratioVector), "/", n = 2)

  numerator <- suppressWarnings(as.numeric(ratioParts[, 1L]))
  denominator <- suppressWarnings(as.numeric(ratioParts[, 2L]))

  ratioValues <- numerator / denominator
  ratioValues[!is.finite(ratioValues)] <- NA_real_

  ratioValues
}

addEnrichmentFold <- function(enrichmentTable) {
  requiredColumns <- c("GeneRatio", "BgRatio")

  missingColumns <- setdiff(requiredColumns, colnames(enrichmentTable))

  if (length(missingColumns) > 0L) {
    stop(
      "The enrichment result is missing required columns: ",
      paste(missingColumns, collapse = ", ")
    )
  }

  enrichmentTable %>%
    mutate(
      geneRatioValue = parseRatio(GeneRatio),
      backgroundRatioValue = parseRatio(BgRatio),
      enrichmentFold = geneRatioValue / backgroundRatioValue
    )
}

selectTopTerms <- function(enrichmentTable, termCount) {
  if (nrow(enrichmentTable) == 0L) {
    return(enrichmentTable)
  }

  enrichmentTable %>%
    filter(
      !is.na(enrichmentFold),
      is.finite(enrichmentFold),
      !is.na(p.adjust)
    ) %>%
    arrange(p.adjust, desc(enrichmentFold)) %>%
    slice_head(n = termCount)
}

saveResultTable <- function(resultTable, filePath) {
  write.table(
    resultTable,
    file = filePath,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

createKeggPlot <- function(keggTable) {
  plotTable <- selectTopTerms(keggTable, topTermCount)

  if (nrow(plotTable) == 0L) {
    warning("No KEGG terms are available for plotting.")
    return(NULL)
  }

  plotTable <- plotTable %>%
    mutate(
      Description = factor(
        Description,
        levels = rev(Description)
      ),
      negativeLog10AdjustedP = -log10(
        pmax(p.adjust, .Machine$double.xmin)
      )
    )

  ggplot(
    plotTable,
    aes(
      x = enrichmentFold,
      y = Description,
      colour = negativeLog10AdjustedP,
      size = Count
    )
  ) +
    geom_point() +
    scale_size(range = c(2, 10)) +
    scale_colour_gradient(
      low = "#4399ab",
      high = "#e5eca2"
    ) +
    labs(
      x = "Enrichment fold",
      y = "KEGG pathway",
      colour = expression(-log[10]("adjusted P")),
      size = "Gene count"
    ) +
    theme_bw(base_size = 14) +
    theme(
      axis.title.y = element_text(margin = margin(r = 20)),
      axis.title.x = element_text(margin = margin(t = 12)),
      axis.text.x = element_text(
        face = "bold",
        colour = "black"
      ),
      axis.text.y = element_text(colour = "black"),
      axis.line = element_line(
        linewidth = 0.5,
        colour = "black"
      ),
      legend.key = element_blank()
    )
}

createGoBpPlot <- function(goTable) {
  plotTable <- selectTopTerms(goTable, topTermCount)

  if (nrow(plotTable) == 0L) {
    warning("No GO Biological Process terms are available for plotting.")
    return(NULL)
  }

  plotTable <- plotTable %>%
    mutate(
      Description = factor(
        Description,
        levels = rev(Description)
      ),
      negativeLog10AdjustedP = -log10(
        pmax(p.adjust, .Machine$double.xmin)
      )
    )

  ggplot(
    plotTable,
    aes(
      x = Description,
      y = enrichmentFold,
      fill = negativeLog10AdjustedP
    )
  ) +
    geom_col() +
    coord_flip() +
    scale_fill_gradient(
      low = "#4399ab",
      high = "#e5eca2"
    ) +
    scale_x_discrete(
      labels = function(termLabels) {
        str_wrap(termLabels, width = 45)
      }
    ) +
    labs(
      x = "GO Biological Process term",
      y = "Enrichment fold",
      fill = expression(-log[10]("adjusted P"))
    ) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "right",
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      axis.title.y = element_text(size = 14)
    )
}

# ------------------------------------------------------------------------------
# 3. Read and map gene identifiers
# ------------------------------------------------------------------------------

geneIds <- readGeneIds(geneFile)

geneMapping <- bitr(
  geneIds,
  fromType = geneIdType,
  toType = c("SYMBOL", "ENTREZID"),
  OrgDb = org.Hs.eg.db
)

if (nrow(geneMapping) == 0L) {
  stop("None of the supplied gene identifiers could be mapped to Entrez IDs.")
}

geneMapping <- geneMapping %>%
  distinct(ENTREZID, .keep_all = TRUE)

entrezIds <- unique(geneMapping$ENTREZID)

mappingFile <- file.path(
  resultDir,
  paste0(analysisName, "_gene_id_mapping.tsv")
)

saveResultTable(geneMapping, mappingFile)

message(
  "Mapped ",
  length(entrezIds),
  " unique Entrez IDs from ",
  length(geneIds),
  " input genes."
)

# ------------------------------------------------------------------------------
# 4. GO Biological Process enrichment
# ------------------------------------------------------------------------------

goResult <- enrichGO(
  gene = entrezIds,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = goOntology,
  pAdjustMethod = "BH",
  pvalueCutoff = pValueCutoff,
  qvalueCutoff = qValueCutoff,
  readable = TRUE
)

goTable <- as.data.frame(goResult)

if (nrow(goTable) > 0L) {
  goTable <- addEnrichmentFold(goTable)

  goResultFile <- file.path(
    resultDir,
    paste0(analysisName, "_go_bp.tsv")
  )

  saveResultTable(goTable, goResultFile)

  goPlot <- createGoBpPlot(goTable)

  if (!is.null(goPlot)) {
    goPlotFile <- file.path(
      plotDir,
      paste0(analysisName, "_go_bp_plot.pdf")
    )

    ggsave(
      filename = goPlotFile,
      plot = goPlot,
      width = goPlotWidth,
      height = goPlotHeight,
      units = "in"
    )
  }
} else {
  warning("GO Biological Process enrichment returned no terms.")
}

# ------------------------------------------------------------------------------
# 5. KEGG enrichment
# ------------------------------------------------------------------------------

keggResult <- enrichKEGG(
  gene = entrezIds,
  organism = "hsa",
  keyType = "kegg",
  pAdjustMethod = "BH",
  pvalueCutoff = pValueCutoff,
  qvalueCutoff = qValueCutoff
)

keggTable <- as.data.frame(keggResult)

if (nrow(keggTable) > 0L) {
  keggResult <- setReadable(
    keggResult,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID"
  )

  keggTable <- as.data.frame(keggResult)
  keggTable <- addEnrichmentFold(keggTable)

  keggResultFile <- file.path(
    resultDir,
    paste0(analysisName, "_kegg.tsv")
  )

  saveResultTable(keggTable, keggResultFile)

  keggPlot <- createKeggPlot(keggTable)

  if (!is.null(keggPlot)) {
    keggPlotFile <- file.path(
      plotDir,
      paste0(analysisName, "_kegg_dotplot.pdf")
    )

    ggsave(
      filename = keggPlotFile,
      plot = keggPlot,
      width = keggPlotWidth,
      height = keggPlotHeight,
      units = "in"
    )
  }
} else {
  warning("KEGG enrichment returned no pathways.")
}

message("GO and KEGG enrichment analysis completed.")
