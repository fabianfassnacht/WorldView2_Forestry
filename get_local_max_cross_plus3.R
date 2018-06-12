require(raster)
require(rgdal)

# 0 = soil
# 1 = shadow
# 2 = coniferous
# 3 = broadleaved

# load rasterfile with classification map
classmap <- stack("/home/nonameornot/KIT/WV2_Paper/Tree_density_workflow/pan/local_max_win_9.tif")
# load shape with cluster positions
setwd("/media/nonameornot/intenso_1500/38_Daniel_Mangold_Biomasse_WV2/09_Density/06_bl_cf_class/shape")
shp <- readOGR(".", "crosses_polygons_plus3")
#plot(classmap)
#plot(shp, add=T, col="red")
#shp <- fishn
# extract raster values of all pixels in the cluster polygons

classes <- extract(classmap, shp, fun=sum)
locmax_7 <- classes

setwd("/home/nonameornot/KIT/WV2_Paper/Tree_density_workflow/pan/")

save(locmax_7, file="local_max_win7.RData")

hist(classes)
