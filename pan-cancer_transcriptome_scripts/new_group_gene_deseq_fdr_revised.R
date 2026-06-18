#!/usr/bin/env Rscript

# ==============================================================================
# Tissue-specific DESeq2, Levene-test, and FDR analysis
# ============================================================================== 

suppressPackageStartupMessages({
  library(DESeq2)
  library(car)
})

# -----------------------------
# 1. File configuration
# -----------------------------
baseDir <- "/data1/DYY/bambu"
countInputDir <- file.path(baseDir, "288_count_gene_by_tissue")
metadataInputDir <- file.path(baseDir, "group_by_tissue")
deseqOutputDir <- file.path(baseDir, "288_count_matrix_deseq2")
leveneOutputDir <- file.path(baseDir, "288_LEVENE")
fdrOutputDir <- file.path(baseDir, "288_FDR")

invisible(lapply(
  c(deseqOutputDir, leveneOutputDir, fdrOutputDir),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
))

# -----------------------------
# 2. Helper functions
# -----------------------------
readCountMatrix <- function(inputFile) {
  inputData <- read.delim(
    inputFile,
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  if (ncol(inputData) < 2L) {
    stop("Count file must contain a feature column and at least one sample column: ", inputFile)
  }

  featureIds <- inputData[[1L]]
  countMatrix <- as.matrix(inputData[-1L])
  rownames(countMatrix) <- featureIds
  storage.mode(countMatrix) <- "numeric"

  if (anyNA(countMatrix)) {
    stop("Missing or non-numeric count values were detected in: ", inputFile)
  }

  if (any(countMatrix < 0)) {
    stop("Negative count values were detected in: ", inputFile)
  }

  if (any(countMatrix %% 1 != 0)) {
    warning("Non-integer counts were rounded in: ", inputFile)
    countMatrix <- round(countMatrix)
  }

  countMatrix
}

alignMetadata <- function(metadataFile, sampleNames) {
  sampleMetadata <- read.delim(
    metadataFile,
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  if (!"Type" %in% colnames(sampleMetadata)) {
    stop("Metadata file is missing the Type column: ", metadataFile)
  }

  sampleIdColumns <- c("Sample.ID", "sample", "ID", "DATA_ID")
  matchedColumns <- intersect(sampleIdColumns, colnames(sampleMetadata))

  if (length(matchedColumns) > 0L) {
    sampleIdColumn <- matchedColumns[1L]
    sampleIds <- as.character(sampleMetadata[[sampleIdColumn]])

    if (anyDuplicated(sampleIds)) {
      stop("Duplicated sample identifiers were detected in: ", metadataFile)
    }

    missingSamples <- setdiff(sampleNames, sampleIds)
    if (length(missingSamples) > 0L) {
      stop(
        "Metadata is missing sample(s) required by the count matrix: ",
        paste(missingSamples, collapse = ", ")
      )
    }

    sampleMetadata <- sampleMetadata[
      match(sampleNames, sampleIds),
      ,
      drop = FALSE
    ]
  } else {
    if (nrow(sampleMetadata) != length(sampleNames)) {
      stop(
        "Metadata has no recognized sample-ID column and its row count does not ",
        "match the number of count-matrix columns: ", metadataFile
      )
    }
  }

  rownames(sampleMetadata) <- sampleNames
  sampleMetadata$Type <- factor(sampleMetadata$Type)

  if (nlevels(sampleMetadata$Type) < 2L) {
    stop("Type must contain at least two groups in: ", metadataFile)
  }

  sampleMetadata
}

calculateLevenePValue <- function(expressionValues, sampleGroups) {
  testResult <- tryCatch(
    car::leveneTest(expressionValues, sampleGroups),
    error = function(errorCondition) NULL
  )

  if (is.null(testResult)) {
    return(NA_real_)
  }

  as.numeric(testResult[["Pr(>F)"]][1L])
}

analyzeDataset <- function(datasetId) {
  countFile <- file.path(countInputDir, paste0(datasetId, ".csv"))
  metadataFile <- file.path(metadataInputDir, paste0(datasetId, ".csv"))

  if (!file.exists(countFile)) {
    stop("Count file not found: ", countFile)
  }

  if (!file.exists(metadataFile)) {
    stop("Metadata file not found: ", metadataFile)
  }

  message("Processing dataset: ", datasetId)

  countMatrix <- readCountMatrix(countFile)
  sampleMetadata <- alignMetadata(metadataFile, colnames(countMatrix))

  deseqData <- DESeqDataSetFromMatrix(
    countData = countMatrix,
    colData = sampleMetadata,
    design = ~ Type
  )

  deseqData <- DESeq(deseqData)
  deseqResults <- results(deseqData, tidy = TRUE, alpha = 0.05)

  deseqResultFile <- file.path(
    deseqOutputDir,
    paste0(datasetId, "_result.csv")
  )

  write.table(
    deseqResults,
    file = deseqResultFile,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )

  normalizedCounts <- counts(deseqData, normalized = TRUE)
  normalizedCountFile <- file.path(
    deseqOutputDir,
    paste0(datasetId, "_norm.counts.fn")
  )

  write.table(
    normalizedCounts,
    file = normalizedCountFile,
    sep = "\t",
    row.names = TRUE,
    col.names = NA,
    quote = FALSE
  )

  featureStandardDeviation <- apply(normalizedCounts, 1L, sd)
  variableFeatures <- normalizedCounts[
    !is.na(featureStandardDeviation) & featureStandardDeviation > 1,
    ,
    drop = FALSE
  ]

  levenePValues <- apply(
    variableFeatures,
    1L,
    calculateLevenePValue,
    sampleGroups = sampleMetadata$Type
  )

  leveneResults <- data.frame(
    featureId = rownames(variableFeatures),
    variableFeatures,
    pValue = levenePValues,
    check.names = FALSE
  )

  leveneResultFile <- file.path(leveneOutputDir, paste0(datasetId, ".csv"))

  write.table(
    leveneResults,
    file = leveneResultFile,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )

  leveneResults$fdr <- p.adjust(leveneResults$pValue, method = "BH")
  leveneResults <- leveneResults[order(leveneResults$pValue), , drop = FALSE]
  significantResults <- leveneResults[
    !is.na(leveneResults$fdr) & leveneResults$fdr < 0.05,
    ,
    drop = FALSE
  ]

  fdrResultFile <- file.path(fdrOutputDir, paste0(datasetId, ".csv"))

  write.table(
    significantResults,
    file = fdrResultFile,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )

  invisible(list(
    datasetId = datasetId,
    deseqResultFile = deseqResultFile,
    normalizedCountFile = normalizedCountFile,
    leveneResultFile = leveneResultFile,
    fdrResultFile = fdrResultFile
  ))
}

# -----------------------------
# 3. Identify and analyze datasets
# -----------------------------
countFiles <- list.files(
  countInputDir,
  pattern = "\\.csv$",
  full.names = FALSE
)

if (length(countFiles) == 0L) {
  stop("No CSV files were found in: ", countInputDir)
}

datasetIds <- unique(substr(countFiles, 1L, 3L))
analysisOutputs <- lapply(datasetIds, analyzeDataset)

message("Completed ", length(analysisOutputs), " dataset(s).")
