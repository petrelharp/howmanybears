library(rgdal)
library(sp)
library(sf)
library(rgeos)
library(raster)

hydro <- readOGR(dsn="./Hydrography_Polygons.shp")
hydro <- crop(hydro,
              extent(c(165182, 791789.8, 5e5, 8e5)))

the_clip <- as(extent(c(4.8e5, 8e5, 5e5, 5.88e5)), "SpatialPolygons")
proj4string(the_clip) <- CRS(proj4string(hydro))
hydro <- gDifference(hydro, the_clip, byid=TRUE)

roads <- gDifference(crop(readOGR(dsn="./State_Owned_Roads_v17a.shp"),
                          extent(c(165182, 791789.8, 5e5, 8e5))),
                     the_clip, byid=TRUE)
counties <- gDifference(crop(readOGR(dsn="./Counties_v17a.shp"),
                             extent(c(165182, 791789.8, 5e5, 8e5))),
                        the_clip, byid=TRUE)

geo <- crop(readOGR(dsn="./Quaternary_Geology_Map.shp"),
                        extent(c(165182, 791789.8, 5e5, 8e5)))
geo <- geo[!as.vector(gIntersects(tmp, the_clip, byid=TRUE)),]

dup <- disaggregate(gUnaryUnion(counties))
up <- dup[which.max(sapply(1:length(dup), function (k) {x = dup[k]; x@bbox[3] - x@bbox[1]}))]

save(up, roads, hydro, file="./up_maps.RData")

png(file="./the_up.png", width=6*144, height=4*144, pointsize=10, res=144)
plot(up, col="#dedbc8")
lines(roads, col='grey')
plot(hydro, col="#79BCDF", add=TRUE, border=NA)
plot(up, add=TRUE)
dev.off()

# the UP has an aspect ratio of almost exactly 2:1
r <- raster(ncol=200, nrow=100)
extent(r) <- extent(up)

upraster <- rasterize(up, r)
georaster <- upraster * 0.1
# make up some values for habitat quality
bear_goodness <- rexp(length(unique(geo$CODE)))
for (k in 1:length(unique(geo$CODE))) {
    usethese <- (geo$CODE == k)
    if (sum(usethese) > 0) {
        georaster <- georaster + bear_goodness[k] * rasterize(geo[usethese,], r, getCover=TRUE)
    }
}

# output values for SLiM
x <- values(georaster)
x[is.na(x)] <- 0
dim(x) <- dim(upraster)[2:1]
x <- x[, ncol(x):1]
cat(x, file='up_raster.txt')

png(file="./bear_goodness.png", width=6*144, height=4*144, pointsize=10, res=144)
plot(georaster, xaxt='n', yaxt='n')
lines(roads, col='black')
plot(hydro, col="#79BCDF", add=TRUE, border=NA)
dev.off()
