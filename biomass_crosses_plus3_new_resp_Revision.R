require(randomForest)
require(autopls)
require(rgdal)


setwd("E:/KIT_Forschung/WV2_Paper/Tree_density_workflow")
#load("bl_cf_fr_crosses_orig.RData")
load("bl_cf_fr_new.RData")#
#load("local_max_win11.RData")
load("local_max_win9.RData")
#load("local_max_win7.RData")
load("localmax_topred.RData")
fracs <- results
load("bl_cf_fr_topred.RData")
#load("bl_cf_fr_new_clouded.RData")

setwd("E:/KIT_Forschung/WV2_Paper/Tree_density_workflow/shape_predict")
shp <- readOGR(".", "dens_pred")

setwd("E:/KIT_Forschung/WV2_Paper/Biomass_workflow")
#load("resp_b.RData")
load("height_pred.RData")
load("height_topred.RData")
heights_to_pred <- results

load("biom_ha_dup_med_qu.RData")
#load("biom_ha_dup_mean_qu.RData")
#bla <- biomass_fin
#load("biom_ha.RData")


resp_new <- data.frame(p = biomass_fin[,1], b=biomass_fin[,2])
resp_new <- resp_new[with(resp_new, order(p)), ]
resp_new2 <- resp_new[,2]

plot(resp_new2, resp_b)

# remove plots covered by clouds
# and outliers (2111, 2083)
cl_plots <- c(1113, 1114, 1143, 1144, 1145, 1146, 1148, 1149, 1163, 2016, 2063, 2091)#, 2111, 2083)
tab <- read.table("BM2.txt", header=T, sep=";")
pt <- tab[,1]%in%cl_plots # Auswahl der Plots in biom_ref
pt3 <- seq(1,80,1)[!pt]

# plot(tab[,1], resp_new[,1])

colnames(locmax_9) <- "locmax"

preds <- cbind(fracs, heights, locmax_9)
preds2 <- preds[pt3,]
#resp_b2 <- resp_b[pt3]
resp_b3 <- resp_new2[pt3]

cor(resp_b3, preds2)

# run random forest model
set.seed(5456)

rfm <- randomForest(preds2, resp_b3, mtry=3, ntree=500)
rfm

############################################
## export data for scatterplot / Teja
wo_outl <- matrix(nrow=68, ncol=2)
wo_outl[,1] <- rfm$predicted
wo_outl[,2] <- resp_b3
colnames(wo_outl) <- c("predicted", "response")
save(wo_outl, file="wo_outl_biomass.RData")
############################################

# create plot predicted vs. modelled

pdf(file="biomass_scatter_cl_rem_2outl.pdf", width=7.5, height=6)

  plot(resp_b3, rfm$predicted, ylab="predicted ABG biomass [t/ha]", xlab="reference ABG biomass [t/ha]", cex=0.9,cex.axis=1.3, cex.lab=1.4)
  #text(resp_b3, rfm$predicted, resp_new[pt3,1])
  cort <- cor(resp_b3, rfm$predicted)^2
  rmse <- sqrt(sum((resp_b3-rfm$predicted)^2)/length(resp_b3))
  te <- paste("r2 =", round(cort, digits=2))
  tex <- paste("rmse =", round(rmse, digits=3))
  mtext(te, side=3, adj=0.06, padj=3, cex=1.5)
  mtext(tex, side=3, adj=0.068, padj=5, cex=1.5)
  
  abline(0,1, col="black")
  m <- lm(resp_b3~ rfm$predicted)
  abline(m$coefficients[1],m$coefficients[2], col="black", lty=2)

dev.off()


# predict to full data
#setwd("/home/nonameornot/KIT/WV2_Paper/Biomass_workflow")

locmax_pred[locmax_pred%in%NA] <- 66

topred_total <- cbind(topred, heights_to_pred, locmax_pred)
colnames(topred_total) <- colnames(preds2)
head(topred_total)

predicted <- predict(rfm, topred_total)

shp@data$pred_biom <- predicted

writeOGR(shp, ".", "pred_biom_locmax_cl_rem_2outl", driver="ESRI Shapefile")




### check height predictors low correlation ###

pdf(file="mean_h_biom_cl_rem_2outl.pdf", height=5, width=5)
plot(resp_b3, preds2[,6], ylab="mean height", xlab="biomass [t/ha]")
dev.off()


### export data to create scatterplot

biom_scatter <- cbind(resp_b3, rfm$predicted)
colnames(biom_scatter) <- c("resp", "pred")

save(biom_scatter, file="scatter_biom.RData")
