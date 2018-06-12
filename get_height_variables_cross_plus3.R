require(raster)
require(rgdal)

# load rasterfile with height values
classmap <- stack("/home/nonameornot/KIT/WV2_Paper/Tree_density_workflow/ndsm/wv2_ndsm_toWVJan.tif")
# load shape with cluster positions
setwd("/media/nonameornot/intenso_1500/38_Daniel_Mangold_Biomasse_WV2/09_Density/06_bl_cf_class/shape")
shp <- readOGR(".", "crosses_polygons_plus3")

plot(classmap)
par(new=T)
plot(shp, add=T)

# create empty matrix to store results
results <- matrix(NA, nrow=length(shp), ncol=4)
# rename column names
colnames(results) <- c("sum", "median", "mean", "stdev") 

# extract raster values of all pixels in the cluster polygons
classes <- extract(classmap, shp)

# iter <- shp@data$plot
# reslist <- list()
# 
# for (i in 1:length(shp)) {
#   
#   sub_shp <- shp[shp@data$plot == iter[i],]
#   classes <- extract(classmap, sub_shp)
#   reslist[[i]] <- classes
#   print(i)
# }


# count pixels of all considered classes (soil, shadow, coniferous, broadleaved) in each cluster-polygon
# calculate percentages of each class

for(i in 1:length(shp)) {
  
  # get first polygon
  data <- classes[[i]]
  
  # calculate percentages of each class
  results[i,1] <- sum(data)
  results[i,2] <- median(data)
  results[i,3] <- mean(data)
  results[i,4] <- sd(data)
}

setwd("/home/nonameornot/KIT/WV2_Paper/Tree_density_workflow/ndsm/")

heights <- results
save(heights, file="height_pred.RData")

# attach results to shapefile
shp$sum_ht <- results[,1]
shp$med_ht <- results[,2]
shp$mean_ht <- results[,3]
shp$sd_ht <- results[,4]

# write out a new shapefile with attached results
writeOGR(shp, ".", "Kreise_gemessen_ndsm", driver="ESRI Shapefile")
