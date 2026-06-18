#!/usr/bin/env Rscript

# ==============================================================================
# Differential expression and variance analysis for the complete transcript set
# ============================================================================== 

suppressPackageStartupMessages({
  library(DESeq2)
  library(car)
})

# -----------------------------
# 1. File configuration
# -----------------------------
countFile <- "/data1/DYY/bambu/288count_transcript.csv"
metadataFile <- "/data1/DYY/bambu/288_group.csv"
deseqOutputDir <- "/data1/DYY/bambu/288_count_matrix_deseq2"
leveneOutputDir <- "/data1/DYY/bambu/288_LEVENE"
fdrOutputDir <- "/data1/DYY/bambu/288_FDR"

resultFile <- file.path(deseqOutputDir, "all_trans_result.csv")
normalizedCountFile <- file.path(deseqOutputDir, "all_trans_norm.counts.fn")
leveneResultFile <- file.path(leveneOutputDir, "all_trans.csv")
fdrResultFile <- file.path(fdrOutputDir, "all_trans.csv")

invisible(lapply(
  c(deseqOutputDir, leveneOutputDir, fdrOutputDir),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
))

# -----------------------------
# 2. Helper functions
# -----------------------------
checkRequiredColumns <- function(inputData, requiredColumns, objectName) {
  missingColumns <- setdiff(requiredColumns, colnames(inputData))

  if (length(missingColumns) > 0L) {
    stop(
      objectName,
      " is missing required column(s): ",
      paste(missingColumns, collapse = ", ")
    )
  }
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

# -----------------------------
# 3. Read and validate inputs
# -----------------------------
if (!file.exists(countFile)) {
  stop("Count file not found: ", countFile)
}

if (!file.exists(metadataFile)) {
  stop("Metadata file not found: ", metadataFile)
}

countData <- read.delim(
  countFile,
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

sampleMetadata <- read.delim(
  metadataFile,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

checkRequiredColumns(
  sampleMetadata,
  c("Sample.ID", "Type"),
  "sampleMetadata"
)

countMatrix <- as.matrix(countData)
storage.mode(countMatrix) <- "numeric"

if (anyNA(countMatrix)) {
  stop("The count matrix contains missing or non-numeric values.")
}

if (any(countMatrix < 0)) {
  stop("The count matrix contains negative values.")
}

if (any(countMatrix %% 1 != 0)) {
  warning("Non-integer values were detected and rounded before DESeq2 analysis.")
  countMatrix <- round(countMatrix)
}

if (anyDuplicated(sampleMetadata$Sample.ID)) {
  stop("Duplicated Sample.ID values were detected in the metadata.")
}

missingSamples <- setdiff(colnames(countMatrix), sampleMetadata$Sample.ID)
if (length(missingSamples) > 0L) {
  stop(
    "The following count-matrix samples are absent from the metadata: ",
    paste(missingSamples, collapse = ", ")
  )
}

sampleMetadata <- sampleMetadata[
  match(colnames(countMatrix), sampleMetadata$Sample.ID),
  ,
  drop = FALSE
]
rownames(sampleMetadata) <- sampleMetadata$Sample.ID

if (!identical(colnames(countMatrix), rownames(sampleMetadata))) {
  stop("Count-matrix columns and metadata rows could not be aligned.")
}

sampleMetadata$sampleGroup <- factor(substr(sampleMetadata$Sample.ID, 1L, 3L))
sampleMetadata$Type <- factor(sampleMetadata$Type)

if (nlevels(sampleMetadata$Type) < 2L) {
  stop("The Type column must contain at least two groups.")
}

# -----------------------------
# 4. DESeq2 analysis
# -----------------------------
deseqData <- DESeqDataSetFromMatrix(
  countData = countMatrix,
  colData = sampleMetadata,
  design = ~ sampleGroup + Type
)

deseqData <- DESeq(deseqData)
deseqResults <- results(deseqData, tidy = TRUE, alpha = 0.05)

write.table(
  deseqResults,
  file = resultFile,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

normalizedCounts <- counts(deseqData, normalized = TRUE)

write.table(
  normalizedCounts,
  file = normalizedCountFile,
  sep = "\t",
  row.names = TRUE,
  col.names = NA,
  quote = FALSE
)

# -----------------------------
# 5. Levene's test
# -----------------------------
sampleGroups <- sampleMetadata$Type

levenePValues <- apply(
  normalizedCounts,
  1L,
  calculateLevenePValue,
  sampleGroups = sampleGroups
)

leveneResults <- data.frame(
  featureId = rownames(normalizedCounts),
  normalizedCounts,
  pValue = levenePValues,
  check.names = FALSE
)

write.table(
  leveneResults,
  file = leveneResultFile,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

# -----------------------------
# 6. Multiple-testing correction
# -----------------------------
leveneResults$fdr <- p.adjust(leveneResults$pValue, method = "BH")
leveneResults <- leveneResults[order(leveneResults$pValue), , drop = FALSE]
significantResults <- leveneResults[
  !is.na(leveneResults$fdr) & leveneResults$fdr < 0.05,
  ,
  drop = FALSE
]

write.table(
  significantResults,
  file = fdrResultFile,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

message("Analysis completed.")
message("DESeq2 results: ", resultFile)
message("Normalized counts: ", normalizedCountFile)
message("Levene results: ", leveneResultFile)
message("FDR-filtered results: ", fdrResultFile)
