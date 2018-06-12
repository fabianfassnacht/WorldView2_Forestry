### load packages

pkgs<-c("rgdal","caret","raster","foreign", "kernlab", "landsat", "e1071", "Hmisc", "vegan", "doSNOW")
lapply(pkgs,require, character.only=T)
library(mlbench)
library(caret)
library(randomForest)

##########
########## run SVM classification
########## 



##########
########## load images to R
##########
#setwd("I:/KIT_Forschung/33_WV2_Jannika/1_Workflow/3_SVM_classification/Input")
#hym <- stack("WV2_MS_8channel_11_texture_measures.tif")
# setwd("I:/KIT_Forschung/33_WV2_Jannika/1_Workflow/3_SVM_classification/Input/MNF_MS_only")
# wv2_ms <-stack("WV2_MNF_8ch")
# 
# setwd("I:/KIT_Forschung/33_WV2_Jannika/1_Workflow/3_SVM_classification/Input/MNF_tex_only")
# wv2_tex <-stack("mnf_text_11_11")
# 
# setwd("I:/KIT_Forschung/33_WV2_Jannika/1_Workflow/3_SVM_classification/Input/MNF_MS_tex_stacked")
# wv2_ms_tex_mnf_stack <- stack("mnf_ms_tex_stack_2m")
# 
# setwd("I:/KIT_Forschung/33_WV2_Jannika/1_Workflow/3_SVM_classification/Input/MNF")
# wv2_mnf_tot <- stack("mnf")

setwd("/media/nonameornot/intenso_1500/33_WV2_Jannika/1_Workflow/3_SVM_classification/Input")
wv2 <- stack("WV2_MS_8channel_11_texture_measures.tif")
# 
# setwd("I:/KIT_Forschung/33_WV2_Jannika/1_Workflow/3_SVM_classification/Input/Original_img_ms")
# wv2_ms_only <- stack("13JUL22105730-M2AS-12EUSI-2318-01-recol.tif")


##########
########## load training samples (point-Shapefile)
##########

setwd("/home/nonameornot/KIT/WV2_Paper/Tree_species_workflow/tr_data")
vec<-readOGR(".","poly_comb1")
trainclass <-vec@data$Art
trainclass1 <- factor(c(as.character(trainclass),as.character(trainclass_sh)))

vec_red <- readOGR(".","poly_comb1wo_linde_hainb")
trainclass_red <-vec_red@data$Art
trainclass1_red <- factor(c(as.character(trainclass_red),as.character(trainclass_sh)))

vec_five <- readOGR(".", "poly_comb1_earl_stud")
trainclass_five <-vec_five@data$Art
trainclass1_five <- factor(c(as.character(trainclass_five),as.character(trainclass_sh)))

vec_sh <- readOGR(".", "shadow_soil")
trainclass_sh <-vec_sh@data$class
table(trainclass_sh)
##########
########## run SVM classification
########## 




# create subset of images (if needed)
# satImage <- wv2_ms[[1:6]]
# satImage <- wv2_tex[[1:6]]
# satImage <- wv2_wo_mnf
# satImage <- wv2_ms_tex_mnf_stack
 satImage <- wv2


#define function to get only high values
#quant <- function(x,...){quantile(x, c(.90))}
#gtmean <- function(x,...){mean(x[x>mean(x)])}


# Extract pixel values from the image for each sample point
trainval <- extract(satImage, vec, fun=mean)
trainval_red <- extract(satImage, vec_red, fun=mean)
trainval_five <- extract(satImage, vec_five, fun=mean)
trainval_sh <- extract(satImage, vec_sh, fun=mean)
#trainval1 <- extract(satImage, vec, fun=quant)
#trainval1 <- extract(satImage, vec, fun=gtmean)

save(trainval, file="trainval.RData")
save(trainval_red, file="trainval_red.RData")
save(trainval_five, file="trainval_five.RData")
save(trainval_sh, file="trainval_sh.RData")

# Get pixel DNs from the image for each sample point (in case of polygons)
#trainval1 <- extract(satImage,vec, fun=mean)

trainval1 <- rbind(trainval, trainval_sh)
trainval1_red <- rbind(trainval_red, trainval_sh)
trainval1_five <- rbind(trainval_five, trainval_sh)

# extract class definitions from shapefile - class definitions have to be stored in
# column names "class" or the vec$class has to be replaced with the corresponding column name
set.seed(1173)

#m1 <- randomForest(trainval1, trainclass, mtry=6, ntrees=500)
#m1

#####
####
# SVM classification (with or without reduced dataset)
####
#####


##############################
##############################
########### Scenario A #######
##############################
##############################

# set output directory
setwd("/home/nonameornot/KIT/WV2_Paper/Tree_species_workflow")

# define output filenames

out2 <- "boxplots_table_wv2tex.txt"
out3 <- "boxplots_hym_wv2tex.pdf"
out4 <- "summed_conf_wv2tex.txt"
outImage <-'SVM_classification_wv2tex.tif'


# prepare dataframe for classification
tl<-ncol(trainval1)-1
data.df<-cbind(trainclass1,trainval1)
d<-ncol(data.df)
N1<-length(data.df[,1])

# set the number of iterations for the iterative validation
N<-100

#############################################
# Apply SVM with kernel radial basis function for all classes 

# set.seed to allow for reproducible results
set.seed(1173)

# define tuning range of grid search
gammat = seq(.1, .9, by = .1)
costt=seq(1,128, by = 8)

# run parameter tuning
tune1 <- tune.svm(trainval1, as.factor(trainclass1), gamma = gammat, cost=costt)
plot(tune1)

# extract best parameters from parameter tuning
gamma <- tune1$best.parameters$gamma
cost <- tune1$best.parameters$cost

#gamma = 0.1
# cost = 41

# train model with all data (to use maximum information for the classification of the full image)
model <- svm(trainval1, as.factor(trainclass1), gamma = gamma, cost = cost, probability = TRUE)
# train model with 5-fold cross-validation to get first impression on accuracies
model2 <- svm(trainval1, as.factor(trainclass1), gamma = gamma, cost = cost, probability = TRUE, cross=5)
# display results of 5-fold cross-validation
summary(model2)

#apply model to image
svmPred <- predict(satImage, model, filename=outImage, na.rm=TRUE, progress='text', format='GTiff', datatype='INT1U',overwrite=TRUE)



##########Now store the results of accuracy across bootstrap samples using the best tuning parameters ##############
set.seed(1173)

############### Create empty slots of Kappa,OA,PA and UA ###################

results <- matrix(0, nrow=N, ncol=(((nlevels(trainclass1)*2)+2)))

CL = length(table(trainclass1))
CL2 <- CL^2
cmat <- matrix(0, ncol = CL2, nrow = N)

################ Run the process to store the value in the empty slots #################
require(kernlab)
for (i in 1:N){
  
  idx=sample(1:N1,N1,replace=TRUE)
  #split the data in test and train
  training.df<-data.df[idx,]
  testing.df<-data.df[-idx,]
  predTest<-testing.df[,c(2:d)]
  predClass<-testing.df[,c(1)]
  svmFit1 <-svm(training.df[,c(2:d)],as.factor(training.df[,c(1)]), gamma = gamma, cost = cost, probability = TRUE)
  pred<-predict(svmFit1,predTest)
  
  if (nlevels(pred) == nlevels(as.factor(predClass))) { 
    
    con.mat<-confusionMatrix(pred,predClass)
    
    for (i2 in 1:(nlevels(trainclass1))) {

        results[i,i2]<-con.mat$byClass[i2,1] # producer's accuracies
        results[i,(i2+nlevels(trainclass1))]<-con.mat$byClass[i2,3] # user's accuracies

      }
    
    results[i, ((nlevels(trainclass1)*2)+1)] <- con.mat$overall[1] # overall accuracy
    results[i, ((nlevels(trainclass1)*2)+2)] <- con.mat$overall[2] # kappa
    
    cm <- table(predClass, pred)
    cmat[i,1:CL2] <- cm[1:CL2] }
}



write.table(results, out2, sep="\t")
write.table(cmat, out4, sep="\t")

names1 <- paste0("PA_cl", 1:nlevels(trainclass1))
names2 <- paste0("UA_cl", 1:nlevels(trainclass1))
names <- matrix("",nrow=1,ncol=(nlevels(trainclass1)*2+2))
names[1:nlevels(trainclass1)] <- names1
names[(nlevels(trainclass1)+1):(nlevels(trainclass1)*2)] <- names2
names[(nlevels(trainclass1)*2)+1] <- "OA"
names[(nlevels(trainclass1)*2)+2] <- "kappa"

trainclass1

####################################### Boxplot of SVM results ##########################################
pdf(file = out3, width = 12, height = 6)
boxplot(results, notch=TRUE, ylim=c(0,1), las=2, names=names)
dev.off()

mean(results[,25])
sd(results[,25])
mean(results[,26])
sd(results[,26])






##############################
##############################
########### Scenario B #######
##############################
##############################

# set output directory
setwd("/home/nonameornot/KIT/WV2_Paper/Tree_species_workflow/results")

# define output filenames

out2 <- "boxplots_table_wv2tex_red.txt"
out3 <- "boxplots_hym_wv2tex_red.pdf"
out4 <- "summed_conf_wv2tex_red.txt"
outImage <-'SVM_classification_wv2tex_red.tif'


# prepare dataframe for classification
tl<-ncol(trainval1_red)-1
data.df<-cbind(trainclass1_red,trainval1_red)
d<-ncol(data.df)
N1<-length(data.df[,1])

# set the number of iterations for the iterative validation
N<-100

#############################################
# Apply SVM with kernel radial basis function for all classes 

# set.seed to allow for reproducible results
set.seed(1173)

# define tuning range of grid search
gammat = seq(.1, .9, by = .1)
costt=seq(1,128, by = 8)

# run parameter tuning
tune1_red <- tune.svm(trainval1_red, as.factor(trainclass1_red), gamma = gammat, cost=costt)
plot(tune1_red)

# extract best parameters from parameter tuning
gamma <- tune1_red$best.parameters$gamma
cost <- tune1_red$best.parameters$cost

#gamma = 0.2
# cost = 9

# train model with all data (to use maximum information for the classification of the full image)
model_red <- svm(trainval1_red, as.factor(trainclass1_red), gamma = gamma, cost = cost, probability = TRUE)
# train model with 5-fold cross-validation to get first impression on accuracies
model2_red <- svm(trainval1_red, as.factor(trainclass1_red), gamma = gamma, cost = cost, probability = TRUE, cross=5)
# display results of 5-fold cross-validation
summary(model2_red)

#apply model to image
svmPred_red <- predict(satImage, model_red, filename=outImage, na.rm=TRUE, progress='text', format='GTiff', datatype='INT1U',overwrite=TRUE)



##########Now store the results of accuracy across bootstrap samples using the best tuning parameters ##############
set.seed(1173)

############### Create empty slots of Kappa,OA,PA and UA ###################

results <- matrix(0, nrow=N, ncol=(((nlevels(trainclass1_red)*2)+2)))

CL = length(table(trainclass1_red))
CL2 <- CL^2
cmat <- matrix(0, ncol = CL2, nrow = N)

################ Run the process to store the value in the empty slots #################
require(kernlab)
for (i in 1:N){
  
  idx=sample(1:N1,N1,replace=TRUE)
  #split the data in test and train
  training.df<-data.df[idx,]
  testing.df<-data.df[-idx,]
  predTest<-testing.df[,c(2:d)]
  predClass<-testing.df[,c(1)]
  svmFit1 <-svm(training.df[,c(2:d)],as.factor(training.df[,c(1)]), gamma = gamma, cost = cost, probability = TRUE)
  pred<-predict(svmFit1,predTest)
  
  if (nlevels(pred) == nlevels(as.factor(predClass))) { 
    
    con.mat<-confusionMatrix(pred,predClass)
    
    for (i2 in 1:(nlevels(trainclass1_red))) {
      
      results[i,i2]<-con.mat$byClass[i2,1] # producer's accuracies
      results[i,(i2+nlevels(trainclass1_red))]<-con.mat$byClass[i2,3] # user's accuracies
      
    }
    
    results[i, ((nlevels(trainclass1_red)*2)+1)] <- con.mat$overall[1] # overall accuracy
    results[i, ((nlevels(trainclass1_red)*2)+2)] <- con.mat$overall[2] # kappa
    
    cm <- table(predClass, pred)
    cmat[i,1:CL2] <- cm[1:CL2] }
}



write.table(results, out2, sep="\t")
write.table(cmat, out4, sep="\t")

names1 <- paste0("PA_cl", 1:nlevels(trainclass1_red))
names2 <- paste0("UA_cl", 1:nlevels(trainclass1_red))
names <- matrix("",nrow=1,ncol=(nlevels(trainclass1_red)*2+2))
names[1:nlevels(trainclass1_red)] <- names1
names[(nlevels(trainclass1_red)+1):(nlevels(trainclass1_red)*2)] <- names2
names[(nlevels(trainclass1_red)*2)+1] <- "OA"
names[(nlevels(trainclass1_red)*2)+2] <- "kappa"

####################################### Boxplot of SVM results ##########################################
pdf(file = out3, width = 12, height = 6)
boxplot(results, notch=TRUE, ylim=c(0,1), las=2, names=names)
dev.off()

mean(results[,21])
sd(results[,21])
mean(results[,22])
sd(results[,22])






##############################
##############################
########### Scenario C #######
##############################
##############################

# set output directory
setwd("/home/nonameornot/KIT/WV2_Paper/Tree_species_workflow/results")

# define output filenames

out2 <- "boxplots_table_wv2tex_five.txt"
out3 <- "boxplots_hym_wv2tex_five.pdf"
out4 <- "summed_conf_wv2tex_five.txt"
outImage <-'SVM_classification_wv2tex_five.tif'


# prepare dataframe for classification
tl<-ncol(trainval1_five)-1
data.df<-cbind(trainclass1_five,trainval1_five)
d<-ncol(data.df)
N1<-length(data.df[,1])

# set the number of iterations for the iterative validation
N<-100

#############################################
# Apply SVM with kernel radial basis function for all classes 

# set.seed to allow for reproducible results
set.seed(1173)

# define tuning range of grid search
gammat = seq(.1, .9, by = .1)
costt=seq(1,128, by = 8)

# run parameter tuning
tune1_five <- tune.svm(trainval1_five, as.factor(trainclass1_five), gamma = gammat, cost=costt)
plot(tune1_five)

# extract best parameters from parameter tuning
gamma <- tune1_five$best.parameters$gamma
cost <- tune1_five$best.parameters$cost

#gamma = 0.1
# cost = 9

# train model with all data (to use maximum information for the classification of the full image)
model_five <- svm(trainval1_five, as.factor(trainclass1_five), gamma = gamma, cost = cost, probability = TRUE)
# train model with 5-fold cross-validation to get first impression on accuracies
model2_five <- svm(trainval1_five, as.factor(trainclass1_five), gamma = gamma, cost = cost, probability = TRUE, cross=5)
# display results of 5-fold cross-validation
summary(model2_five)

#apply model to image
svmPred_five <- predict(satImage, model_five, filename=outImage, na.rm=TRUE, progress='text', format='GTiff', datatype='INT1U',overwrite=TRUE)



##########Now store the results of accuracy across bootstrap samples using the best tuning parameters ##############
set.seed(1173)

############### Create empty slots of Kappa,OA,PA and UA ###################

results <- matrix(0, nrow=N, ncol=(((nlevels(trainclass1_five)*2)+2)))

CL = length(table(trainclass1_five))
CL2 <- CL^2
cmat <- matrix(0, ncol = CL2, nrow = N)

################ Run the process to store the value in the empty slots #################
require(kernlab)
for (i in 1:N){
  
  idx=sample(1:N1,N1,replace=TRUE)
  #split the data in test and train
  training.df<-data.df[idx,]
  testing.df<-data.df[-idx,]
  predTest<-testing.df[,c(2:d)]
  predClass<-testing.df[,c(1)]
  svmFit1 <-svm(training.df[,c(2:d)],as.factor(training.df[,c(1)]), gamma = gamma, cost = cost, probability = TRUE)
  pred<-predict(svmFit1,predTest)
  
  if (nlevels(pred) == nlevels(as.factor(predClass))) { 
    
    con.mat<-confusionMatrix(pred,predClass)
    
    for (i2 in 1:(nlevels(trainclass1_five))) {
      
      results[i,i2]<-con.mat$byClass[i2,1] # producer's accuracies
      results[i,(i2+nlevels(trainclass1_five))]<-con.mat$byClass[i2,3] # user's accuracies
      
    }
    
    results[i, ((nlevels(trainclass1_five)*2)+1)] <- con.mat$overall[1] # overall accuracy
    results[i, ((nlevels(trainclass1_five)*2)+2)] <- con.mat$overall[2] # kappa
    
    cm <- table(predClass, pred)
    cmat[i,1:CL2] <- cm[1:CL2] }
}



write.table(results, out2, sep="\t")
write.table(cmat, out4, sep="\t")

names1 <- paste0("PA_cl", 1:nlevels(trainclass1_five))
names2 <- paste0("UA_cl", 1:nlevels(trainclass1_five))
names <- matrix("",nrow=1,ncol=(nlevels(trainclass1_five)*2+2))
names[1:nlevels(trainclass1_five)] <- names1
names[(nlevels(trainclass1_five)+1):(nlevels(trainclass1_five)*2)] <- names2
names[(nlevels(trainclass1_five)*2)+1] <- "OA"
names[(nlevels(trainclass1_five)*2)+2] <- "kappa"

####################################### Boxplot of SVM results ##########################################
pdf(file = out3, width = 12, height = 6)
boxplot(results, notch=TRUE, ylim=c(0,1), las=2, names=names)
dev.off()

mean(results[,15])
sd(results[,15])
mean(results[,16])
sd(results[,16])


