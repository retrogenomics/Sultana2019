#!/usr/bin/env Rscript

######################################################################################
# Plot correlation between number of insertions per chr and wgs coverage of each chr #
######################################################################################

library("ggpubr")

# import dataframe from .bed file
mydata <- read.delim("hg19.l1neo.soni.loc.helas3.bed",
  header = TRUE,
  sep = "\t",
  fill = TRUE,
  comment.char = "#",
  check.names = FALSE,
)

# calculate count and percentage of each strand
colnames(mydata) <- c("chr", "start", "stop", "name", "score", "strand")
orientation <- data.frame(
  strand = c("plus", "minus"),
  count = c(sum(mydata$strand == "+"), sum(mydata$strand == "-"))
)
orientation$perc <- round(prop.table(orientation$count) * 100, digit = 1)

# prepare labels of piechart
labs <- paste0(orientation$strand, "\n(", orientation$perc, "%)")

# test if proportion between strands deviates from randomness
binom.test(
  x = orientation$count[orientation$strand == "plus"],
  n = sum (orientation$count),
  p = 0.5,
  alternative = "two.sided"
)

# generate piechart
pdf.options(
  width = 1.5,
  height = 1.5,
  family = "Helvetica",
  useDingbats = FALSE
)
pdf(file = "StrandPieChart.pdf")
ggpie(orientation,
  "count",
  label = labs,
  fill = "strand",
  palette = "npg",
  color = "white",
  lab.pos = "in",
  lab.font = c(2, "bold", "white"),
  size = 1
) +
rremove("legend")
dev.off()
