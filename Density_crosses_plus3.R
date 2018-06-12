require(randomForest)
require(autopls)
require(hier.part)

setwd("/media/nonameornot/intenso_1500/WV2_Paper/Tree_density_workflow")
load("response.RData")
#load("bl_cf_fr_crosses_orig.RData")
load("bl_cf_fr_new.RData")#
load("bl_cf_fr_topred.RData")
load("local_max_win11.RData")
load("local_max_win9.RData")
load("local_max_win7.RData")
load("localmax_topred.RData")

#loc_max <- cbind(locmax_11, locmax_9, locmax_7)
#(cor(resp, results))^2
#(cor(resp, locmax_7))^2

table(is.na(topred))


locmax_pred[locmax_pred%in%NA] <- 66
colnames(locmax_pred) <- "locmax"
topred2 <- cbind(topred, locmax_pred)

preds <- cbind(results, locmax_9)
colnames(preds) <- colnames(topred2)
  


setwd("/home/nonameornot/KIT/WV2_Paper/Tree_density_workflow/shape_predict")
shp <- readOGR(".", "dens_pred")

# create rf model
set.seed(56)
rfm <- randomForest(preds, resp, mtry=2, ntree=500)
rfm

# create scatterplot
pdf(file="density_scatter.pdf", width=7.5, height=6)

  plot(resp, rfm$predicted, ylab="predicted trees/ha", xlab="reference trees/ha", cex=0.9,cex.axis=1.3, cex.lab=1.4)
  cort <- cor(resp, rfm$predicted)^2
  rmse <- sqrt(sum((resp-rfm$predicted)^2)/length(resp))
  te <- paste("r2 =", round(cort, digits=2))
  tex <- paste("rmse =", round(rmse, digits=3))
  mtext(te, side=3, adj=0.06, padj=3, cex=1.5)
  mtext(tex, side=3, adj=0.068, padj=5, cex=1.5)
  
  abline(0,1, col="black")
  m <- lm(resp~ rfm$predicted)
  abline(9.4569, 0.9433, col="black", lty=2)

dev.off()

# predict to full dataset

head(preds)
head(topred2)

predicted <- predict(rfm, topred2)
shp@data$pred_den_loc <- predicted
writeOGR(shp, ".", "predicted_dens_loc_max", driver="ESRI Shapefile")
