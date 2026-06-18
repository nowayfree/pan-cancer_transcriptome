# ==============================================================================
# Supplementary Figure 3: selected GO Biological Process terms
#
# This script visualizes selected downregulated GO terms across tissues.
# Point colour represents statistical significance, and point size represents
# the gene ratio associated with each enriched term.
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(stringr)
})

# ------------------------------------------------------------------------------
# 1. Configuration
# ------------------------------------------------------------------------------

inputFile <- "/home/bwu4/go_kegg/Supplementary_Table_2.csv"
outputFile <- "/data1/DYY/plot/supplementary_figure_3.pdf"

selectedTerms <- c(
  "anion transmembrane transport",
  "glucose metabolic process",
  "sodium ion homeostasis"
)

figureWidth <- 10
figureHeight <- 10

# ------------------------------------------------------------------------------
# 2. Helper functions
# ------------------------------------------------------------------------------

checkRequiredFile <- function(filePath) {
  if (!file.exists(filePath)) {
    stop("Required input file was not found: ", filePath)
  }
}

checkRequiredColumns <- function(dataTable, requiredColumns, tableName) {
  missingColumns <- setdiff(requiredColumns, colnames(dataTable))

  if (length(missingColumns) > 0L) {
    stop(
      tableName,
      " is missing required columns: ",
      paste(missingColumns, collapse = ", ")
    )
  }
}

parseRatio <- function(ratioVector) {
  ratioParts <- str_split_fixed(
    as.character(ratioVector),
    pattern = "/",
    n = 2
  )

  numerator <- suppressWarnings(as.numeric(ratioParts[, 1L]))
  denominator <- suppressWarnings(as.numeric(ratioParts[, 2L]))

  ratioValues <- numerator / denominator
  ratioValues[!is.finite(ratioValues)] <- NA_real_

  ratioValues
}

# ------------------------------------------------------------------------------
# 3. Read and validate enrichment results
# ------------------------------------------------------------------------------

checkRequiredFile(inputFile)

enrichmentTable <- read.table(
  file = inputFile,
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  quote = "",
  comment.char = ""
)

requiredColumns <- c(
  "tissue",
  "Description",
  "down_gene_ratio",
  "down_p.adjust"
)

checkRequiredColumns(
  dataTable = enrichmentTable,
  requiredColumns = requiredColumns,
  tableName = "Supplementary enrichment table"
)

# ------------------------------------------------------------------------------
# 4. Prepare plotting data
# ------------------------------------------------------------------------------

plotData <- enrichmentTable %>%
  select(
    tissue,
    Description,
    down_gene_ratio,
    down_p.adjust
  ) %>%
  mutate(
    geneRatio = parseRatio(down_gene_ratio),
    adjustedPValue = suppressWarnings(
      as.numeric(down_p.adjust)
    )
  ) %>%
  filter(
    Description %in% selectedTerms,
    !is.na(tissue),
    !is.na(Description),
    !is.na(geneRatio),
    !is.na(adjustedPValue),
    geneRatio > 0,
    adjustedPValue >= 0,
    adjustedPValue <= 1
  ) %>%
  mutate(
    negativeLog10AdjustedP = -log10(
      pmax(adjustedPValue, .Machine$double.xmin)
    ),
    Description = factor(
      Description,
      levels = rev(selectedTerms)
    )
  ) %>%
  arrange(Description, tissue)

if (nrow(plotData) == 0L) {
  stop(
    "No rows remained after filtering for the selected GO terms."
  )
}

missingTerms <- setdiff(
  selectedTerms,
  unique(as.character(plotData$Description))
)

if (length(missingTerms) > 0L) {
  warning(
    "The following selected GO terms were not found: ",
    paste(missingTerms, collapse = ", ")
  )
}

# ------------------------------------------------------------------------------
# 5. Draw and save Supplementary Figure 3
# ------------------------------------------------------------------------------

supplementaryFigure3 <- ggplot(
  plotData,
  aes(
    x = tissue,
    y = Description,
    colour = negativeLog10AdjustedP,
    size = geneRatio
  )
) +
  geom_point() +
  scale_size(
    range = c(2, 8),
    name = "Gene ratio"
  ) +
  scale_colour_gradient(
    low = "blue",
    high = "red",
    name = expression(-log[10]("adjusted P"))
  ) +
  labs(
    x = NULL,
    y = "GO Biological Process term"
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.title.y = element_text(
      margin = margin(r = 30)
    ),
    axis.text.x = element_text(
      face = "bold",
      colour = "black",
      angle = 0,
      vjust = 1
    ),
    axis.text.y = element_text(
      colour = "black"
    ),
    axis.line = element_line(
      linewidth = 0.5,
      colour = "black"
    ),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.key = element_blank()
  )

outputDirectory <- dirname(outputFile)
dir.create(
  outputDirectory,
  recursive = TRUE,
  showWarnings = FALSE
)

ggsave(
  filename = outputFile,
  plot = supplementaryFigure3,
  width = figureWidth,
  height = figureHeight,
  units = "in"
)

message(
  "Supplementary Figure 3 was saved to: ",
  outputFile
)
