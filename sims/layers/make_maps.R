library(rgdal)
library(sp)
library(sf)
library(rgeos)
library(mapshaper)

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

dup <- disaggregate(gUnaryUnion(counties))
up <- dup[which.max(sapply(1:length(dup), function (k) {x = dup[k]; x@bbox[3] - x@bbox[1]}))]

save(up, roads, hydro, file="./up_maps.RData")

png(file="./the_up.png", width=6*144, height=4*144, pointsize=10, res=144)
plot(up, col="#dedbc8")
lines(roads, col='grey')
plot(hydro, col="#79BCDF", add=TRUE, border=NA)
plot(up, add=TRUE)
dev.off()


