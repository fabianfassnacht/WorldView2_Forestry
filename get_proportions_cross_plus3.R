require(raster)
require(rgdal)

# 0 = soil
# 1 = shadow
# 2 = coniferous
# 3 = broadleaved

# load rasterfile with classification map
classmap <- stack("/home/nonameornot/KIT/WV2_Paper/Tree_density_workflow/result_class/SVM_classification_wv2tex_4cl.tif")
# load shape with cluster positions
setwd("/media/nonameornot/intenso_1500/38_Daniel_Mangold_Biomasse_WV2/09_Density/06_bl_cf_class/shape")
shp <- readOGR(".", "crosses_polygons_plus3")
#plot(classmap)
#plot(shp, add=T, col="red")
#shp <- fishn

# create empty matrix to store results
results <- matrix(NA, nrow=length(shp), ncol=4)
# rename column names
colnames(results) <- c("broadleaved", "coniferous", "shadow", "soil") 

# extract raster values of all pixels in the cluster polygons
classes <- extract(classmap, shp)

# count pixels of all considered classes (soil, shadow, coniferous, broadleaved) in each cluster-polygon
# calculate percentages of each class

for(i in 1:length(shp)) {
  
  # get first polygon
  data <- classes[[i]]
  
  # calculate percentages of each class
  results[i,1] <- length(data[data==1])/length(data)
  results[i,2] <- length(data[data==2])/length(data)
  results[i,3] <- length(data[data==3])/length(data)
  results[i,4] <- length(data[data==4])/length(data)
}

# attach results to shapefile
shp$percsoil <- results[,1]
shp$percshad <- results[,2]
shp$perccf <- results[,3]
shp$percbl <- results[,4]

setwd("/home/nonameornot/KIT/WV2_Paper/Tree_density_workflow/result_class/")

# write out a new shapefile with attached results
writeOGR(shp, ".", "Vektorgitter_class_metrics", driver="ESRI Shapefile")
save(results, file="bl_cf_fr_new.RData")

head(shp@data)

