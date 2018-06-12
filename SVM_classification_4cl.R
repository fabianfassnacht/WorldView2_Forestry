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

setwd("/home/nonameornot/KIT/WV2_Paper/WV2_scene_tex")
wv2 <- stack("WV2_MS_8channel_11_texture_measures.tif")
# 
# setwd("I:/KIT_Forschung/33_WV2_Jannika/1_Workflow/3_SVM_classification/Input/Original_img_ms")
# wv2_ms_only <- stack("13JUL22105730-M2AS-12EUSI-2318-01-recol.tif")


##########
########## load training samples (point-Shapefile)
##########

setwd("/home/nonameornot/KIT/WV2_Paper/Tree_density_workflow/shape")
vec<-readOGR(".","samples")
trainclass <-vec@data$class
table(trainclass)
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
trainval1 <- extract(satImage, vec, fun=mean)
trainclass1 <- vec@data$class

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
########### Scenario D #######
##############################
##############################

# set output directory
setwd("/home/nonameornot/KIT/WV2_Paper/Tree_density_workflow")

# define output filenames

out2 <- "boxplots_table_wv2tex_4cl.txt"
out3 <- "boxplots_hym_wv2tex_4cl.pdf"
out4 <- "summed_conf_wv2tex_4cl.txt"
outImage <-'SVM_classification_wv2tex_4cl.tif'


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

trainval1 <- trainval

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

mean(results[,9])
sd(results[,9])
mean(results[,10])
sd(results[,10])





