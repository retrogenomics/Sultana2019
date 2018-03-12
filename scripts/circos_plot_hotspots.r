#!/usr/bin/env Rscript

# circos plotting of insertion density with red bars corresponding to hotspots (Poisson, FDR<=0.05)

# load libraries
library(circlize)

# load and reformat data
file_list <- list.files(path = "./", pattern = "poisson.txt$")
bin_size <- sub("Mb_bin_poisson.txt","",sub("^nb_ins_per_","",file_list))
ordered_file_list <- file_list[order(bin_size)]


# prepare axis and colors
color <- c("black", "red")
spacing <- c(rep(1,22),5) # unequal spacing will give space for the y-axis scale bar

pdf("Hotspots.pdf")


# initialize plot
circos.par(start.degree = 90, gap.degree = spacing)
circos.initializeWithIdeogram(chromosome.index = c(paste0("chr", 1:22),"chrX"), plotType = "ideogram", ideogram.height = 0.03)

# add labels of chromosomes
for(chr in get.all.sector.index()) {
    xlim = get.cell.meta.data("cell.xlim", sector.index = chr)
    ylim = get.cell.meta.data("cell.ylim", sector.index = chr)
    circos.text(mean(xlim), ylim[2], gsub("chr", "", chr), niceFacing = TRUE, cex = 0.8, font = 2, adj = c(0.5, -0.5), sector.index = chr)
}

# generate track for each bin size
for (file in ordered_file_list) {
	
	
	# load and prepare data tables
	myfile <- read.delim(file, header = TRUE, sep = "\t", quote = "\"", dec = ".", fill = TRUE, comment.char = "", check.names = FALSE)
	file_ID <- sub("Mb_bin_poisson.txt","",sub("^nb_ins_per_","",file))
	colnames(myfile)[1] <- "chr"
	mydata <- cbind(myfile[,1:3],myfile[,7:9])
	
	# separate bins above and below FDR of 5%
	hotspot0 <- subset(mydata, mydata$poisson_fdr > 0.05)
	hotspot1 <- subset(mydata, mydata$poisson_fdr <= 0.05)
	if(dim(hotspot1)[1]>=1){
	  hotspots <- list(hotspot0, hotspot1)
	}else{
	   hotspots <- list(hotspot0)
	}
	
	# prepare y-axis ticks
	ticks <- seq(from=0, to=max(mydata$nb_ins), by=(max(mydata$nb_ins)%/%3)) # define ticks for y-axis and dotted lines
	track_legend <- 2.5*(max(mydata$nb_ins)%/%3)

	# plot bars and lines
	circos.genomicTrackPlotRegion(hotspots, panel.fun = function(region, value, ...) {
		i=getI(...)
		circos.genomicLines(region, value, numeric.column = 1, type="h", col = color[i], lwd = 1)
		cell.xlim = get.cell.meta.data("cell.xlim")
		for (i in ticks){
			circos.lines(cell.xlim, c(i, i), lty = 3, cex = 0.5, col = "#888888")
		}
	})
	
	# draw y-axis and labels
	first_sector = get.all.sector.index()[1]
	last_sector = get.all.sector.index()[length(get.all.sector.index())]
	circos.yaxis(labels.cex=0.5, side = c("left"), labels.niceFacing = TRUE, sector.index = first_sector, at = ticks)
	xlim = get.cell.meta.data("cell.xlim", sector.index = first_sector)
	ylim = get.cell.meta.data("cell.ylim", sector.index = first_sector)
	# circos.text(xlim[1], mean(ylim), file_ID, facing = "clockwise", niceFacing = FALSE, cex = 0.5, adj = c(1, degree(4)), sector.index = first_sector)
	circos.text(xlim[1], track_legend, paste(file_ID,"Mb", sep = " "), facing = "bending.inside", niceFacing = FALSE, cex = 0.5, font = 2, adj = c(-0.5, 0.25), sector.index = first_sector)
}

# end circos plot
circos.clear()


dev.off()

