#!/usr/bin/env Rscript

usage <- "
Usage:
    plot_sibs.R (treefile name)
" 

args <- commandArgs(TRUE)
if (length(args) != 1) {
    cat(usage)
    q('no')
}

library(rgdal)
library(sp)
library(sf)
library(rgeos)
library(raster)

outbase <- gsub(".trees$", "", args[1])
outfile <- paste0(outbase, ".bear_sibs.pdf")

load('layers/up_maps.RData')
sibs <- read.table(paste0(outbase, ".sib_locs.tsv"))
samples <- read.table(paste0(outbase, ".sample_locs.tsv"))

slim_coords <- function (xy) {
    xy[,1] <-  xy[,1] * (extent(up)[2] - extent(up)[1])/530 + extent(up)[1]
    xy[,2] <-  xy[,2] * (extent(up)[4] - extent(up)[3])/265 + extent(up)[3]
    return(xy)
}
sibs[,1:2] <- slim_coords(sibs[,1:2])
sibs[,3:4] <- slim_coords(sibs[,3:4])
samples <- slim_coords(samples)

pdf(file=outfile, width=6, height=3, pointsize=10)
par(mar=c(0,0,0,0)+.1)
    plot(up, col="#dedbc8")
    lines(roads, col='grey')
    plot(hydro, col="#79BCDF", add=TRUE, border=NA)
    plot(up, add=TRUE)
    points(samples, pch=20, cex=0.25)
    segments(x0=sibs[,1], x1=sibs[,3],
             y0=sibs[,2], y1=sibs[,4],
             col='red')
dev.off()
