#!/usr/bin/env Rscript

# ==============================================================================
# Bubble plot for pathway effect sizes and significance values
# ============================================================================== 

suppressPackageStartupMessages({
  library(ggplot2)
  library(reshape2)
})

# -----------------------------
# 1. File configuration
# -----------------------------
effectSizeFile <- "/data1/DYY/bambu/3f_data_generatio.csv"
pValueFile <- "/data1/DYY/bambu/3f_data_p.csv"
outputFile <- "/data1/DYY/bambu/plot/figure3/3f.pdf"

significanceCutoff <- 0.05
separatorPositions <- c(20.5, 39.5, 58.5, 78.5, 95.5, 114.5, 134.5)

outputDir <- dirname(outputFile)
dir.create(outputDir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 2. Read and validate inputs
# -----------------------------
if (!file.exists(effectSizeFile)) {
  stop("Effect-size file not found: ", effectSizeFile)
}

if (!file.exists(pValueFile)) {
  stop("P-value file not found: ", pValueFile)
}

effectSizeData <- read.delim(
  effectSizeFile,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

pValueData <- read.delim(
  pValueFile,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

requiredColumn <- "Description"
if (!requiredColumn %in% colnames(effectSizeData)) {
  stop("Effect-size data must contain a Description column.")
}

if (!requiredColumn %in% colnames(pValueData)) {
  stop("P-value data must contain a Description column.")
}

if (anyDuplicated(effectSizeData$Description)) {
  stop("Duplicated Description values were detected in the effect-size data.")
}

if (anyDuplicated(pValueData$Description)) {
  stop("Duplicated Description values were detected in the P-value data.")
}

# -----------------------------
# 3. Reshape and merge data
# -----------------------------
effectSizeLong <- reshape2::melt(
  effectSizeData,
  id.vars = "Description",
  variable.name = "group",
  value.name = "effectSize"
)

pValueLong <- reshape2::melt(
  pValueData,
  id.vars = "Description",
  variable.name = "group",
  value.name = "pValue"
)

effectSizeLong$effectSize <- as.numeric(effectSizeLong$effectSize)
pValueLong$pValue <- as.numeric(pValueLong$pValue)

if (anyNA(effectSizeLong$effectSize)) {
  stop("Non-numeric or missing effect-size values were detected.")
}

if (anyNA(pValueLong$pValue)) {
  stop("Non-numeric or missing P values were detected.")
}

if (any(pValueLong$pValue < 0 | pValueLong$pValue > 1)) {
  stop("P values must fall between 0 and 1.")
}

minimumPositivePValue <- .Machine$double.xmin
pValueLong$pValue <- pmax(pValueLong$pValue, minimumPositivePValue)
pValueLong$negativeLog10P <- -log10(pValueLong$pValue)

plotData <- merge(
  effectSizeLong,
  pValueLong[c("Description", "group", "negativeLog10P")],
  by = c("Description", "group"),
  all = FALSE,
  sort = FALSE
)

expectedRows <- nrow(effectSizeLong)
if (nrow(plotData) != expectedRows) {
  stop(
    "Effect-size and P-value tables do not contain the same Description-group combinations."
  )
}

plotData$Description <- factor(
  plotData$Description,
  levels = rev(unique(effectSizeData$Description))
)
plotData$group <- factor(
  plotData$group,
  levels = colnames(effectSizeData)[colnames(effectSizeData) != "Description"]
)

negativeLog10Cutoff <- -log10(significanceCutoff)
maximumSignificance <- max(plotData$negativeLog10P, na.rm = TRUE)

# -----------------------------
# 4. Generate bubble plot
# -----------------------------
figure5f <- ggplot() +
  geom_point(
    data = subset(plotData, negativeLog10P < negativeLog10Cutoff),
    mapping = aes(x = group, y = Description, size = effectSize),
    color = "lightgray"
  ) +
  geom_point(
    data = subset(plotData, negativeLog10P >= negativeLog10Cutoff),
    mapping = aes(
      x = group,
      y = Description,
      size = effectSize,
      color = negativeLog10P
    )
  ) +
  scale_color_gradientn(
    colors = c("darkgreen", "orange"),
    limits = c(negativeLog10Cutoff, maximumSignificance),
    oob = scales::squish,
    name = expression(-log[10](italic(P)))
  ) +
  scale_size_continuous(
    range = c(1, 10),
    name = "Effect size"
  ) +
  geom_vline(
    xintercept = separatorPositions,
    linewidth = 0.5
  ) +
  labs(
    title = "Normal",
    x = NULL,
    y = NULL
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(
      fill = NA,
      color = "black",
      linewidth = 1,
      linetype = "solid"
    ),
    axis.text.x = element_text(
      angle = 90,
      hjust = 0.5,
      vjust = 0.5
    )
  )

# -----------------------------
# 5. Save figure
# -----------------------------
ggsave(
  filename = outputFile,
  plot = figure5f,
  device = "pdf",
  width = 10,
  height = 5,
  units = "in",
  dpi = 300
)

message("Figure saved to: ", outputFile)
