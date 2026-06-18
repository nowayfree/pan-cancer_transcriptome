#!/usr/bin/env Rscript

# ==============================================================================
# Transcript structure visualization with transPlotR
#
# Usage:
#   Rscript transplotr_revised.R <annotation.gtf> <gene_name> \
#     <transcript_id_1> [transcript_id_2 ...]
#
# Example:
#   Rscript transplotr_revised.R gencode.v42.annotation.gtf TP53 \
#     ENST00000269305 ENST00000420246
#
# The resulting PDF is written to the configured output directory.
# ==============================================================================

suppressPackageStartupMessages({
  library(transPlotR)
  library(rtracklayer)
  library(ggplot2)
})

# ------------------------------------------------------------------------------
# 1. Configuration
# ------------------------------------------------------------------------------

outputDir <- "/data1/DYY/bambu/plot/gene_structure"
figureHeight <- 8
figureWidth <- 12
arrowType <- "open"

# ------------------------------------------------------------------------------
# 2. Helper functions
# ------------------------------------------------------------------------------

printUsage <- function() {
  message(
    paste(
      "Usage:",
      "Rscript transplotr_revised.R <annotation.gtf> <gene_name>",
      "<transcript_id_1> [transcript_id_2 ...]"
    )
  )
}

checkRequiredFile <- function(filePath) {
  if (!file.exists(filePath)) {
    stop("Required GTF file was not found: ", filePath)
  }
}

removeEnsemblVersion <- function(identifierVector) {
  sub("\\.[0-9]+$", "", as.character(identifierVector))
}

findMatchingColumn <- function(columnNames, candidates) {
  matchedColumn <- candidates[candidates %in% columnNames]

  if (length(matchedColumn) == 0L) {
    return(NULL)
  }

  matchedColumn[[1L]]
}

prepareGtfTable <- function(gtfTable) {
  geneColumn <- findMatchingColumn(
    colnames(gtfTable),
    c("gene_name", "gene_id")
  )

  transcriptIdColumn <- findMatchingColumn(
    colnames(gtfTable),
    c("transcript_id", "transcript_name")
  )

  if (is.null(geneColumn)) {
    stop(
      "The GTF annotation contains neither 'gene_name' nor 'gene_id'."
    )
  }

  if (is.null(transcriptIdColumn)) {
    stop(
      paste(
        "The GTF annotation contains neither",
        "'transcript_id' nor 'transcript_name'."
      )
    )
  }

  if (!"gene_name" %in% colnames(gtfTable)) {
    gtfTable$gene_name <- as.character(gtfTable[[geneColumn]])
  }

  if (!"transcript_id" %in% colnames(gtfTable)) {
    gtfTable$transcript_id <- as.character(
      gtfTable[[transcriptIdColumn]]
    )
  }

  if (!"transcript_name" %in% colnames(gtfTable)) {
    gtfTable$transcript_name <- as.character(
      gtfTable$transcript_id
    )
  }

  gtfTable$gene_name <- as.character(gtfTable$gene_name)
  gtfTable$transcript_id <- as.character(gtfTable$transcript_id)
  gtfTable$transcript_name <- as.character(
    gtfTable$transcript_name
  )

  gtfTable$transcriptIdWithoutVersion <- removeEnsemblVersion(
    gtfTable$transcript_id
  )

  gtfTable
}

# ------------------------------------------------------------------------------
# 3. Parse command-line arguments
# ------------------------------------------------------------------------------

commandArguments <- commandArgs(trailingOnly = TRUE)

if (length(commandArguments) < 3L) {
  printUsage()
  stop(
    "At least one transcript ID must be supplied.",
    call. = FALSE
  )
}

gtfFile <- commandArguments[[1L]]
targetGene <- commandArguments[[2L]]
requestedTranscriptIds <- unique(commandArguments[3:length(commandArguments)])

requestedTranscriptIds <- requestedTranscriptIds[
  !is.na(requestedTranscriptIds) &
    nzchar(trimws(requestedTranscriptIds))
]

if (length(requestedTranscriptIds) == 0L) {
  stop("No valid transcript IDs were supplied.")
}

checkRequiredFile(gtfFile)
dir.create(outputDir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 4. Import and validate the GTF annotation
# ------------------------------------------------------------------------------

gtfRanges <- import(gtfFile)
gtfTable <- as.data.frame(gtfRanges)
gtfTable <- prepareGtfTable(gtfTable)

availableGenes <- unique(gtfTable$gene_name)

if (!targetGene %in% availableGenes) {
  stop(
    "Target gene was not found in the GTF annotation: ",
    targetGene
  )
}

targetGeneTable <- gtfTable[
  gtfTable$gene_name == targetGene,
  ,
  drop = FALSE
]

requestedTranscriptIdsWithoutVersion <- removeEnsemblVersion(
  requestedTranscriptIds
)

availableTranscriptIdsWithoutVersion <- unique(
  targetGeneTable$transcriptIdWithoutVersion
)

missingTranscriptIds <- requestedTranscriptIds[
  !requestedTranscriptIdsWithoutVersion %in%
    availableTranscriptIdsWithoutVersion
]

if (length(missingTranscriptIds) > 0L) {
  stop(
    "The following transcripts were not found for ",
    targetGene,
    ": ",
    paste(missingTranscriptIds, collapse = ", ")
  )
}

matchedTranscriptIds <- unique(
  targetGeneTable$transcript_id[
    targetGeneTable$transcriptIdWithoutVersion %in%
      requestedTranscriptIdsWithoutVersion
  ]
)

# Remove the temporary matching column before passing the table to transPlotR.
gtfTable$transcriptIdWithoutVersion <- NULL

# ------------------------------------------------------------------------------
# 5. Draw and save the transcript structure
# ------------------------------------------------------------------------------

transcriptPlot <- trancriptVis(
  gtfFile = gtfTable,
  gene = targetGene,
  arrowType = arrowType,
  myTranscript = matchedTranscriptIds
)

safeGeneName <- gsub(
  pattern = "[^A-Za-z0-9_.-]",
  replacement = "_",
  x = targetGene
)

outputFile <- file.path(
  outputDir,
  paste0(safeGeneName, ".pdf")
)

ggsave(
  filename = outputFile,
  plot = transcriptPlot,
  width = figureWidth,
  height = figureHeight,
  units = "in"
)

message("Transcript structure figure saved to: ", outputFile)
