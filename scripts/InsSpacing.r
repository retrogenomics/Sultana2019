#!/usr/bin/env Rscript

###################################################
# Plot distribution of spacing between insertions #
###################################################

library("ggpubr")
library("reshape2")
library("scales")

# import dataframe from file (InsSpacing.sh output)
mydata <- read.delim("InsSpacing.tab",
  header = TRUE,
  sep = "\t",
  fill = TRUE,
  comment.char = "#",
  check.names = FALSE,
)
mydata2 <- melt(mydata,
  measure.vars = c("l1neo", "random", "mrc"),
  variable.name = "datasets",
  value.name = "distance",
  na.rm = TRUE
)
fs <- 6
pdf.options(
  width = 2.5,
  height = 2.5,
  family = "Helvetica",
  useDingbats = FALSE
)
pdf(file = "InsSpacing.pdf")
p <- ggecdf(mydata2,
  x = "distance",
  color = "datasets",
  palette = "npg",
  size = 0.5
) +
scale_x_log10(
  breaks = scales::trans_breaks("log10", function(x) 10 ^ x),
  labels = scales::trans_format("log10", scales::math_format(10 ^ .x))
) +
annotation_logticks(sides = "b", size = 0.176389) +
font("legend.title", size = 6) +
font("legend.text", size = 6)
ggpar(p,
  xlab = "Insertion-to-insertion spacing",
  ylab = "Frequency",
  xlim = c(1e2, 1e7),
  ylim = c(0, 1),
  font.x = c(fs, "bold"),
  font.y = c(fs, "bold"),
  font.legend = 6,
  font.tickslab = 6,
  legend = c(0.2, 0.8)
) +
rremove("axis") +
rremove("legend.title") +
theme(
  axis.ticks = element_line(size = 0.25),
  panel.border = element_rect(size = 0.5, fill = NA)
)
dev.off()
