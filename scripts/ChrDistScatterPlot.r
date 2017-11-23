#!/usr/bin/env Rscript

######################################################################################
# Plot correlation between number of insertions per chr and wgs coverage of each chr #
######################################################################################

library("ggpubr")

# import dataframe from ChrDistTableMaker.sh output file
mydata <- read.delim("chromosomal_distribution.tab",
  header = TRUE,
  sep = "\t",
  quote = "\"",
  dec = ".",
  fill = TRUE,
  comment.char = "",
  check.names = FALSE,
  row.names = 1
)

# create a table to run and store binomial test result and multitesting correction for each chromosome
a <- cbind(
  hits = mydata$loc,
  total = sum(mydata$loc),
  prob = (mydata$wgs / sum(mydata$wgs))
)
a <- as.data.frame(a, row.names = row.names(mydata))
test <- function(z) {
  binom.test(x = z[1], n = z[2], p = z[3])$p.value
}
a$pval <- cbind(apply(a, 1, test))
a <- cbind(a, padj = p.adjust(a$pval, method = "fdr"))

# add columns to original dataframe with:
#   - significance of binomial test
#   - name of chr
# this will be used for coloring and labelling in the scatterplot
mydata$sig <- ifelse(a$padj <= 0.05, TRUE, FALSE)
mydata$name <- row.names(mydata)

# generate scatterplot
pdf.options(
  width = 2.5,
  height = 2.5,
  family = "Helvetica",
  useDingbats = FALSE
)
pdf(file = "ChrDistScatterPlot.pdf")
fs <- 6
p <- ggscatter(mydata,
  x = "wgs",
  y = "loc",
  color = "sig",
  alpha = 0.6,
  palette = c("black", "red"),
  size = 1,
  add = "reg.line",
  add.params = list (color = "black", size = 0.5, linetype = "dashed"),
  fullrange = TRUE,
  label = "name",
  label.select = subset(mydata, sig == TRUE)$name,
  font.label = list (color = "black", size = fs),
  repel = TRUE,
  cor.coef = TRUE,
  cor.coef.size = 2,
  cor.coeff.args = list(
    method = "pearson",
    label.x.npc = "left",
    label.y.npc = "top",
    label.sep = "\n"
  )
)
ggpar(p,
  xlab = "HeLaS3 WGS coverage",
  ylab = "# target loci",
  font.x = c(fs, "bold"),
  font.y = c(fs, "bold"),
  font.tickslab = c(fs)
) +
scale_x_continuous(
  labels = scales::comma
) +
rremove("legend") +
rremove("axis") +
theme(
  axis.ticks = element_line(size = 0.25),
  panel.border = element_rect(size = 0.5, fill = NA)
)
dev.off()
