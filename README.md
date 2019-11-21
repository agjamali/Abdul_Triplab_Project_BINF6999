# CAMH_Krembil_Center_for_Neuroinformatics_Project
Using Machine Learning and Statistical Model fitting to predict electrophysiological parameters of Mus Musculus neuron cell types using gene expression data

## LOAD LIBRARIES ##

install.packages("glmnet")
install.packages("bootstrap")
install.packages("DAAG")
install.packages("bestglm")
install.packages("tidyverse")
install.packages("homologene")
library(glmnet)
library(caret)
library(tidyverse)
library(pROC)
library(bootstrap)
library(ggplot2)
library(tidyr)
library(DAAG)
library(reshape)
library(biomaRt)
library(org.Mm.eg.db)
library(annotate)
library(homologene)

## Installing Bioconductor for conversion of Entrez IDs to Gene names
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

## Loading correct packages for conversion of Entrez IDs to Gene names
BiocManager::install(c("biomaRt","org.Mm.eg.db", "annotate"))

browseVignettes("biomaRt")

## DATA MANIPULATION ##

# Set working directory to the one that includes our data set
setwd("C:/Users/abdul_jamali/Downloads")

# Import data as a csv file from directory
neuroEphysData <- read.csv("TableS6.csv")

# View data to ensure validity
# head(neuroEphysData)

# Make the first row of data into rownames
neuroEphysData2 <- neuroEphysData[,-1]
rownames(neuroEphysData2) <- neuroEphysData[,1]

# Transpose and convert data into a data frame
neuroEphysData2 <- t(neuroEphysData2)
neuroEphysData2 <- as.data.frame(neuroEphysData2)

## Converting all the gene entrez IDs to gene names
entrezIDs <- colnames(neuroEphysData2[,29:45796])

## Ensure the data base used for conversion is the Mus Musculus data set: "org.Mm.eg"
GeneNameDF <- getSYMBOL(entrezIDs, data='org.Mm.eg')

colnames(neuroEphysData2[,29:45796]) <- GeneNameDF

neuroEphysData3 <- neuroEphysData2[,29:45796]
colnames(neuroEphysData3) <- GeneNameDF

# Apply the converted Gene names to data frame format
neuroEphysData2[,29:45796] <- neuroEphysData3

# Combine new Gene name data frame with just the Electrophysiological responses
neuroEphysData4 <- data.frame(neuroEphysData2[,1:28],neuroEphysData3)

## Specific Gene and Electrophysiological Response data for simple visualization of Gene expression and response plots ##

cellInfo = rownames(neuroEphysData4) %>% as.data.frame()
colnames(cellInfo) = 'sample_name'

cellInfo = cellInfo %>% separate('sample_name', c('cre', 'layer', 'class'), sep = '__')


gene_expression_plot <- function(geneName,EphysFeatureName){
  
  GeneExpressionData <- neuroEphysData5[,geneName]
  EphysFeatureData <- neuroEphysData5[,EphysFeatureName]
  DF <- data.frame(GeneExpressionData,EphysFeatureData)
  DF <- data.frame(DF, cellInfo)
  
  print(ggplot(DF, aes(x=GeneExpressionData, y=EphysFeatureData)) +
    geom_point(aes(shape = class, color = class), size = 4) +
    scale_shape_manual(values = c(18, 16)) +
    scale_color_manual(values = c("#E69F00", "#009E73")) + 
    theme(axis.text=element_text(size=14), axis.text.x=element_text(angle=45), axis.title=element_text(size=14,face="bold"), legend.title=element_text(size=14), legend.text=element_text(size=14)) +
    ggtitle('Gene Expression Plot of Gm40627 vs f_i_curve_slope') +
    xlab('(log2 CPM+1)') +
    ylab('geneName'))
  
}


gene_expression_plot(geneName='Gm40627',EphysFeatureName='apamp')


AHP_Amplitude <- as.matrix(neuroEphysData4[5,])
AP_Amplitude <- as.matrix(neuroEphysData2[4,])
f_i_curve_slope <- as.matrix(neuroEphysData4[,13])
apthr <- as.matrix(neuroEphysData4[,3])
rmp <- as.matrix(neuroEphysData4[,1])

colnames(neuroEphysData4[,1:28]) <- 

Ggrasp_entrez_id = '67298'
Gprasp1 <- as.matrix(neuroEphysData2[Ggrasp_entrez_id,])

Camk2g_entrez_id = '12325'
Camk2g <- as.matrix(neuroEphysData2[Camk2g_entrez_id,])

Xxylt1_entrez_id = '268880'
Xxylt1 <- as.matrix(neuroEphysData2[Xxylt1_entrez_id,])

Clic4_gene_symbol = 'Clic4'
Clic4 <- as.matrix(neuroEphysData4[,Clic4_gene_symbol])

Gm40627_gene_symbol = 'Gm40627'
Gm40627 <- as.matrix(neuroEphysData4[,Gm40627_gene_symbol])
Gm40627 <- as.data.frame(Gm40627)

Ttc21a_gene_symbol = 'Ttc21a'
Ttc21a <- as.matrix(neuroEphysData4[,Ttc21a_gene_symbol])
Ttc21a <- as.data.frame(Ttc21a)

Lrrn1_gene_symbol = 'Lrrn1'
Lrrn1 <- as.matrix(neuroEphysData4[,Lrrn1_gene_symbol])
Lrrn1 <- as.data.frame(Lrrn1)

sample_meta <- sample_meta[-1,]
sample_meta[,4] = f_i_curve_slope
sample_meta[,5] = Gm40627
sample_meta[,6]= Ttc21a
sample_meta[,7] = apthr
sample_meta[,8] = rmp
sample_meta[,9] = Lrrn1

ggplot(sample_meta, aes(x=V5, y=V4)) +
  geom_point(aes(shape = class, color = class), size = 4) +
  scale_shape_manual(values = c(18, 16)) +
  scale_color_manual(values = c("#E69F00", "#009E73")) + 
  theme(axis.text=element_text(size=14), axis.text.x=element_text(angle=45), axis.title=element_text(size=14,face="bold"), legend.title=element_text(size=14), legend.text=element_text(size=14)) +
  ggtitle('Gene Expression Plot of Gm40627 vs f_i_curve_slope') +
  xlab('Gm40627 (log2 CPM+1)') +
  ylab('f_i_curve_slope')

ggplot(sample_meta, aes(x=V1, y=V7)) +
  geom_point(aes(shape = class, color = class), size = 4) +
  scale_shape_manual(values = c(18, 16)) +
  scale_color_manual(values = c("#E69F00", "#009E73")) + 
  theme(axis.text=element_text(size=14), axis.text.x=element_text(angle=45), axis.title=element_text(size=14,face="bold"), legend.title=element_text(size=14), legend.text=element_text(size=14)) +
  ggtitle('Gene Expression Plot of Ttc21a vs apthr') +
  xlab('Ttc21a (log2 CPM+1)') +
  ylab('apthr')


ggplot(sample_meta, aes(x=V1.1, y=V9)) +
  geom_point(aes(shape = class, color = class), size = 4) +
  scale_shape_manual(values = c(18, 16)) +
  scale_color_manual(values = c("#E69F00", "#009E73")) + 
  theme(axis.text=element_text(size=14), axis.text.x=element_text(angle=45), axis.title=element_text(size=14,face="bold"), legend.title=element_text(size=14), legend.text=element_text(size=14)) +
  ggtitle('Gene Expression Plot of Lrrn1 vs rmp') +
  xlab('Lrrn1 (log2 CPM+1)') +
  ylab('rmp')


Clic4DF <- as.data.frame(AHP_Amplitude)

ElectroPhysio <- neuroEphysData2[,1:28]
ElectroPhysio <- t(ElectroPhysio)

GeneExpression <- neuroEphysData2[,29:45796]
GeneExpression <- t(GeneExpression)

## Begin fitting the multivariate regression LASSO model

model <- lm(formula=GeneExpression ~ ElectroPhysio, data=neuroEphysData2)

sampleNeuroEphysData2 <- neuroEphysData2[1:40,1:40]
sampleNeuroEphysData2$cellType <- rownames(sampleNeuroEphysData2)
sampleNeuroEphysData2 = sampleNeuroEphysData2 %>% separate('cellType', c('cre', 'layer', 'class'), sep = '__')


#Splitting the data into 34 cell types for the training set and 14 cell types of the validation set
ind <- sample(1:dim(neuroEphysData2)[1], 36, replace = FALSE) 

train.neuroEphys <- neuroEphysData2[ind, ]
test.neuroEphys <- neuroEphysData2[-ind, ]

trainingY <- as.matrix(train.neuroEphys[,4])
testingY <- as.matrix(test.neuroEphys[,4])

trainingX <- as.matrix(train.neuroEphys[,29:45796])
testingX <- as.matrix(test.neuroEphys[,29:45796])

## Change Response to ONE of the electrophysiological parameters (ie. apamp)
Lasso.neuroEphys <- cv.glmnet(trainingX, trainingY, alpha = 1, type.measure="mse", family="gaussian")



set.seed(1235)
nrfolds <- 4
nEphysFeatures <- 16

folds <- rep_len(1:nrfolds, nrow(neuroEphysData2))
folds <- sample(folds, nrow(neuroEphysData2))
APamp <- list()

save_df = list()
save_df2 = list()
save_df3 = list()
r2 = list()
shuffledr2 = list()
save_df2_shuffled = list()
save_values = list()
pred_df=list()
shuffledTestingY = list()
shuffledPred_df = list()
EphysFeature = list()
rsqShuffled = list()
rsq=list()
rss=list()
tss=list()
rssShuffled=list()
tssShuffled=list()
rsqShuffled=list()
rsq4=list()
rsq4shuffled=list()

## Leave One Out Cross-Validation across all Electrophysiological features from Gene expression data
## For 4-fold CV, change nrfolds to 4
for(x in 1:nEphysFeatures){
  for(i in 4){
    ind = which(folds != i)
    
    train.neuroEphys <- neuroEphysData4[ind, ]
    test.neuroEphys <- neuroEphysData4[-ind, ]
    
    trainingY <- as.matrix(train.neuroEphys[,1])
    testingY <- as.matrix(test.neuroEphys[,1])
    #shuffledY <- as.matrix(train.neuroEphys[sample(nrow(trainingY)),x])
    
    trainingX <- as.matrix(train.neuroEphys[,29:45796])
    testingX <- as.matrix(test.neuroEphys[,29:45796])
    
    Lasso.neuroEphys <- cv.glmnet(trainingX, trainingY, alpha = 1, type.measure="mse", family="gaussian")
    #ShuffledLasso.neuroEphys <- cv.glmnet(trainingX, shuffledY, alpha = 1, type.measure="mse", family="gaussian")
    
    fit <- glmnet(trainingX,
                  trainingY,
                  alpha = 1,
                  lambda = Lasso.neuroEphys$lambda.min,
                  family="gaussian")
    
    #shuffledFit <- glmnet(trainingX,
    #                      shuffledY,
    #                      alpha = 1,
    #                      lambda = ShuffledLasso.neuroEphys$lambda.min,
    #                      family="gaussian")
    
    prds.test <- predict(fit,newx = as.matrix(testingX), type = "response", s=Lasso.neuroEphys$lambda.min)
    #Shuffledprds.test <- predict(shuffledFit,newx = as.matrix(testingX), type = "response", s=ShuffledLasso.neuroEphys$lambda.min)
    ##(cbind(prds.test, testingY))
    
    save_df[[i]] <- data.frame(prds.test, testingY)
    
    pred_df[i] <- prds.test
    
    save_df2[i] <- testingY
    
    #shuffledPred_df[i] <- Shuffledprds.test
    
    #shuffled_values_df <- data.frame(Shuffledprds.test, testingY)
    
    rss[i] <- sum((prds.test - testingY) ^ 2)  ## residual sum of squares
    tss[i] <- sum((testingY - mean(testingY)) ^ 2)  ## total sum of squares
    rssVal <- as.numeric(rss[i])
    tssVal <- as.numeric(tss[i])
    rsq[i] <- 1 - (rssVal/tssVal)
    
    #rssShuffled[i] <- sum((Shuffledprds.test - testingY) ^ 2)  ## residual sum of squares
    #tssShuffled[i] <- sum((testingY - mean(testingY)) ^ 2)  ## total sum of squares
    #rssValShuffled <- as.numeric(rssShuffled[i])
    #tssValShuffled <- as.numeric(tssShuffled[i])
    #rsqShuffled[i] <- 1 - (rssValShuffled/tssValShuffled)
  }
  EphysFeature[x] <- colnames(train.neuroEphys)[x]
  
  rsq4[[x]] <- data.frame(rsq, EphysFeature[x])
  rsq4shuffled[[x]] <- data.frame(rsqShuffled, EphysFeature[x])
  
  #LOOCV Code:
  #pred_df <- unlist(pred_df)
  
  #save_df2 <- as.numeric(unlist(t(save_df2)))
  #r2[x] <- summary(lm(pred_df ~ save_df2))$r.squared
  #r2 <- summary(lm(pred_df ~ save_df2))$r.squared
  
  # Generates an r2 value of 0.888292 for AHPamp
  # Generates an r2 value of 0.654871 for AHPw
  #rss[x] <- sum((pred_df - save_df2) ^ 2)  ## residual sum of squares
  #tss[x] <- sum((save_df2 - mean(save_df2)) ^ 2)  ## total sum of squares
  #rssVal <- as.numeric(rss[x])
  #tssVal <- as.numeric(tss[x])
  #rsq[x] <- 1 - (rssVal/tssVal)
  # Generates an r2 value of 0.872248 for AHPamp
  # Generates an r2 value of 0.647678 for AHPw
  #shuffledPred_df <- unlist(shuffledPred_df)
  
  #shuffledr2[x] <- summary(lm(shuffledPred_df ~ save_df2))$r.squared
  
  #rssShuffled[x] <- sum((shuffledPred_df - save_df2) ^ 2)  ## residual sum of squares
  #tssShuffled[x] <- sum((save_df2 - mean(save_df2)) ^ 2)  ## total sum of squares
  #rssValShuffled <- as.numeric(rssShuffled[x])
  #tssValShuffled <- as.numeric(tssShuffled[x])
  #rsqShuffled[x] <- 1 - (rssValShuffled/tssValShuffled)
  
}

rsq4DF <- rsq4
names(rsq4DF[[1]]) <- c('rsquared', 'rsquared', 'rsquared', 'rsquared', 'EphysFeatures')

for (i in 1:nEphysFeatures){
  names(rsq4DF[[i]]) <- names(rsq4DF[[1]])
}

rsq4DF <- do.call("rbind",rsq4DF)
rsq4DF$Data <- 'AIBS'

rsq4shuffledDF <- rsq4shuffled
names(rsq4shuffledDF[[1]]) <- c('rsquared', 'rsquared', 'rsquared', 'rsquared', 'EphysFeatures')

for (i in 1:nEphysFeatures){
  names(rsq4shuffledDF[[i]]) <- names(rsq4shuffledDF[[1]])
}

rsq4shuffledDF <- do.call("rbind",rsq4shuffledDF)
rsq4shuffledDF$Data <- 'Shuffled'

FinalFinal4DF <- rbind(rsq4DF, rsq4shuffledDF)

FinalFinal5DF <- cbind(FinalFinal4DF[5:6], stack(FinalFinal4DF[1:4]))

ggplot(FinalFinal5DF, aes(EphysFeatures, values, fill=Data)) +
  geom_boxplot() +
  scale_x_discrete() +
  scale_y_continuous(limits=c(-1, 1)) +
  theme(axis.text=element_text(size=14), axis.text.x=element_text(angle=45), axis.title=element_text(size=14,face="bold"), legend.title=element_text(size=14), legend.text=element_text(size=14)) +
  ggtitle('Predicting EphysFeatures using Gene Expression Data via 4-fold CV - Boxplot') +
  ylab('R^2 (4-fold CV)') +
  xlab('EphysFeature')

# Coefficient analysis on 4-fold CV model (Go back and run the 4-fold CV for specific Ephys Feature "x" and pick a random iteration "i" before executing the following code)
coef(Lasso.neuroEphys, s=Lasso.neuroEphys$lambda.min) %>%
  broom::tidy() %>%
  filter(row != "(Intercept)") %>%
  ggplot(aes(value, reorder(row, value), color = value > 0)) +
  geom_point(show.legend = FALSE, size=4) +
  theme(axis.text=element_text(size=14), axis.text.x=element_text(angle=45), axis.title=element_text(size=14,face="bold"), legend.title=element_text(size=14), legend.text=element_text(size=14)) +
  ggtitle("Coefficient Aanalysis in LASSO of EphysFeature: rmp prediction for iteration: 4 of the 4-fold CV") +
  xlab("Coefficient") +
  ylab("Gene Name") 

# Predicted vs Observed Plots for LASSO results
ggplot(save_df[[i]], aes(x=testingY, y=X1)) +
  geom_point(size = 4) +
  theme(axis.text=element_text(size=14), axis.text.x=element_text(angle=45), axis.title=element_text(size=14,face="bold"), legend.title=element_text(size=14), legend.text=element_text(size=14)) +
  geom_smooth(method = "lm") +
  ggtitle('Predicitng rmp with 4-fold CV iteration: 1') +
  xlab('Observed rmp') +
  ylab('Predicted rmp')

# Gene Expression plots from the Ceofficient analysis


# LOOCV data manipulation code to generate visualization of r^2 across all folds

'
df2 <- data.frame(rsq,EphysFeature)
rsquaredDF <- as.data.frame(t(df2[,1:nEphysFeatures]))
rsquaredDF <- rsquaredDF[,1:2]
colnames(rsquaredDF)[1] <- "rsquared"
rsquaredDF <- rsquaredDF$rsquared

FinalDataFrame <- data.frame(rsquaredDF, EphysFeature)
FinalDataFrame$Data <- "AIBS"
colnames(FinalDataFrame)[1] <- "rsquared"
FinalDataFrame <- FinalDataFrame[,1:2]
FinalDataFrame[,2] <- unlist(EphysFeature[])
colnames(FinalDataFrame)[2] <- "EphysFeature"

save_shuffled_df <- data.frame(rsqShuffled, EphysFeature)
rsquaredDF_shuffled <- as.data.frame(t(save_shuffled_df[,1:nEphysFeatures]))
rsquaredDF_shuffled <- rsquaredDF_shuffled[,1:2]
colnames(rsquaredDF_shuffled)[1] <- "rsquared"
rsquaredDF_shuffled <- rsquaredDF_shuffled$rsquared

FinalDataFrameShuffled <- data.frame(rsquaredDF_shuffled, EphysFeature)
FinalDataFrameShuffled$Data <- "Shuffled"
colnames(FinalDataFrameShuffled)[1] <- "rsquared"
FinalDataFrameShuffled <- FinalDataFrameShuffled[,1:2]
FinalDataFrameShuffled[,2] <- unlist(EphysFeature[])
colnames(FinalDataFrameShuffled)[2] <- "EphysFeature"

FinalFinalDF <- rbind(FinalDataFrame, FinalDataFrameShuffled)
'

ggplot(FinalDataFrame, aes(EphysFeature, rsquared)) +
  geom_boxplot() +
  scale_x_discrete(labels=abbreviate) +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  ggtitle('Predicting EphysFeatures using Gene expression Boxplot via LOOCV on Normal Data') +
  ylab('R^2 (LOOCV)') +
  xlab('EphysFeature')

ggplot(FinalDataFrameShuffled, aes(EphysFeature, rsquared)) +
  geom_boxplot() +
  scale_x_discrete(labels=abbreviate) +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  ggtitle('Predicting EphysFeatures using Gene expression Boxplot via LOOCV on Shuffled Data') +
  ylab('Shuffled R^2 (LOOCV)') +
  xlab('EphysFeature')

mse = function(x,y) { mean((x-y)^2)}

set.seed(1235)
## Declaring variables for cell class 
nrfolds <- 2
nEphysFeatures <- 16
r2 = list()
save_df = list()
save_df2 = list()
save_df3 = list()
r2 = list()
shuffledr2 = list()
save_values = list()
pred_df=list()
class=list()
FinalFinalDF=list()
FinalFinalDF2=list()
EphysFeature=list()
hyperParameter=list()
MSEValue=list()
ngenes=list()
rss1=list()
tss1=list()
rsq1=list()
zScoreMean=list()
zSD = list()
zScore = list()


## Checking Inhibitory vs Excitatory cell class
for(x in 1:nEphysFeatures){
  for(i in 1:2){
    if(i == 1){
      train.neuroEphys <- neuroEphysData4[grep("*inh",rownames(neuroEphysData4)),]
      test.neuroEphys <- neuroEphysData4[grep("*exc",rownames(neuroEphysData4)),]
      class[i] = "train:inh/test:exc"
    }
    else
    {
      train.neuroEphys <- neuroEphysData4[grep("*exc",rownames(neuroEphysData4)),]
      test.neuroEphys <- neuroEphysData4[grep("*inh",rownames(neuroEphysData4)),]
      class[i] = "train:exc/test:inh"
    }
    
    
    trainingY <- as.matrix(train.neuroEphys[,13])
    testingY <- as.matrix(test.neuroEphys[,13])
    
    #Zsdtrain <- sd(trainingY)*sqrt((length(trainingY)-1)/(length(trainingY)))
    #zScoretrain = (trainingY - mean(trainingY))/ Zsd
    
    trainingX <- as.matrix(train.neuroEphys[,29:45796])
    testingX <- as.matrix(test.neuroEphys[,29:45796])
    
    Lasso.neuroEphys <- cv.glmnet(trainingX, trainingY, alpha = 1, type.measure="mse", family="gaussian")
    
    fit <- glmnet(trainingX,
                  trainingY,
                  alpha = 1,
                  lambda = Lasso.neuroEphys$lambda.min,
                  family="gaussian")
    
    MSEValue[[i]] = Lasso.neuroEphys$cvm[Lasso.neuroEphys$lambda == Lasso.neuroEphys$lambda.1se]
  
    hyperParameter[[i]] = Lasso.neuroEphys$lambda.min
    #plot(Lasso.neuroEphys)
    
    ngenes[[i]]= as.matrix(coef(Lasso.neuroEphys, Lasso.neuroEphys$lambda.min))
    #rownames(ngenes[[i]][ngenes[[i]]!=0])
    
    prds.test <- predict(fit,newx = as.matrix(testingX), type = "response", s=Lasso.neuroEphys$lambda.min)
    
    save_df[[i]] <- data.frame(prds.test, testingY)
    
    pred_df[i] <- prds.test
    
    save_df2[i] <- testingY
    
    #r2[i] <- summary(lm(prds.test ~ testingY))$r.squared
    rss1[i] <- sum((prds.test - testingY) ^ 2)  ## residual sum of squares
    tss1[i] <- sum((testingY - mean(testingY)) ^ 2)  ## total sum of squares
    rssVal <- as.numeric(rss1[i])
    tssVal <- as.numeric(tss1[i])
    rsq1[i] <- 1 - (rssVal/tssVal)
    
    #Zsd <- sd(prds.test)*sqrt((length(prds.test)-1)/(length(prds.test)))
    #zScorePreds = (prds.test - mean(prds.test))/ Zsd
    

    #ZsdObs <- sd(testingY)*sqrt((length(testingY)-1)/(length(testingY)))
    #zScoreObs = (testingY - mean(testingY))/ ZsdObs
    
    #zScore[[i]] <- data.frame(zScorePreds, zScoreObs)
    
    
    FinalFinalDF[[i]] <- data.frame(rsq1[i], class[i])
    colnames(FinalFinalDF[[i]]) <- c("rsquared", "class split")
    
    
    
  }
  EphysFeature[x] <- colnames(train.neuroEphys)[x]
  FinalFinalDF2[[x]] <- as.data.frame(FinalFinalDF)
}

for (x in 1:nEphysFeatures){
  FinalFinalDF2[[x]] <- as.data.frame(FinalFinalDF2)
  FinalFinalDF2[[x]]$EphysFeatures <- EphysFeature[x]
}

ggplot(save_df[[1]], aes(x=X1, y=testingY)) +
  geom_point(size = 3) +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_smooth(method = "lm") +
  geom_text(x = 25, y = 300, label = r2[1], parse = TRUE) +
  ggtitle('Predicitng f_i_curve_slope with class split: exc/inh') +
  xlab('Observed AHPamp (mV)') +
  ylab('Predicted AHPamp (mV)')

ggplot(save_df[[2]], aes(x=X1, y=testingY)) +
  geom_point(size = 3) +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_smooth(method = "lm") +
  ggtitle('Predicitng AHPamp with class split: exc/inh') +
  xlab('Observed AHPamp (mV)') +
  ylab('Predicted AHPamp (mV)')

ggplot(zScore[[1]], aes(x=X1, y=zScoreObs)) +
  geom_point(size = 3) +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_smooth(method = "lm") +
  ggtitle('Predicitng Z-score AHPamp with class split: inh/exc') +
  xlab('Observed Z-score AHPamp (mV)') +
  ylab('Predicted Z-score AHPamp (mV)')

ggplot(zScore[[2]], aes(x=X1, y=zScoreObs)) +
  geom_point(size = 3) +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  geom_smooth(method = "lm") +
  ggtitle('Predicitng Z-score AHPamp with class split: exc/inh') +
  xlab('Observed Z-score AHPamp (mV)') +
  ylab('Predicted Z-score AHPamp (mV)')

coef(Lasso.neuroEphys, s=Lasso.neuroEphys$lambda.min) %>%
  broom::tidy() %>%
  filter(row != "(Intercept)") %>%
  ggplot(aes(value, reorder(row, value), color = value > 0)) +
  geom_point(show.legend = FALSE, size=4) +
  theme(axis.text=element_text(size=14), axis.text.x=element_text(angle=45), axis.title=element_text(size=14,face="bold"), legend.title=element_text(size=14), legend.text=element_text(size=14)) +
  ggtitle("Coefficient Aanalysis in LASSO of EphysFeature: f_i_curve_slope prediction for train:exc/test:inh of the Cell Class split") +
  xlab("Coefficient") +
  ylab("Gene Name") 

ggplot(save_df[[2]], aes(x=testingY, y=X1)) +
  geom_point(size = 4) +
  theme(axis.text=element_text(size=14), axis.text.x=element_text(angle=45), axis.title=element_text(size=14,face="bold"), legend.title=element_text(size=14), legend.text=element_text(size=14)) +
  geom_smooth(method = "lm") +
  ggtitle('Predicitng f_i_curve_slope with class split: exc/inh') +
  xlab('Observed f_i_curve_slope') +
  ylab('Predicted f_i_curve_slope')

##pred_df <- unlist(pred_df)
##save_df2 <- as.numeric(unlist(t(save_df2)))

view(FinalFinalDF2[[2]])

##shuffledPred_df <- unlist(shuffledPred_df)

##shuffledr2[x] <- summary(lm(shuffledPred_df ~ save_df2))$r.squared


FinalFinalDF3 <- FinalFinalDF2


for (i in 1:nEphysFeatures){
  names(FinalFinalDF3[[i]]) <- names(FinalFinalDF3[[]])
}

FinalFinalDF3 <- do.call("rbind",FinalFinalDF3)

FinalFinalDF3$EphysFeatures <- unlist(EphysFeature)


ggplot(FinalFinalDF3, aes(EphysFeatures, rsquared)) +
  geom_boxplot() +
  ggtitle('Predicting EphysFeature Correlation Coefficients using Gene expression Boxplot for train:inh/test:exc') +
  scale_y_continuous() +
  theme(axis.text=element_text(size=14), axis.text.x=element_text(angle=45), axis.title=element_text(size=14,face="bold"), legend.title=element_text(size=14), legend.text=element_text(size=14)) +
  ylab('R^2 (train:inh/test:exc)') +
  xlab('EphysFeature')

ggplot(FinalFinalDF3, aes(EphysFeatures, rsquared.1)) +
  geom_boxplot() +
  ggtitle('Predicting EphysFeature Correlation Coefficients using Gene expression Boxplot for train:exc/test:inh') +
  scale_y_continuous() +
  theme(axis.text=element_text(size=14), axis.text.x=element_text(angle=45), axis.title=element_text(size=14,face="bold"), legend.title=element_text(size=14), legend.text=element_text(size=14)) +
  ylab('R^2 (train:exc/test:inh)') +
  xlab('EphysFeature')

## Fid the most important variables contributing to the LASSO prediction model
coef(Lasso.neuroEphys,s=Lasso.neuroEphys$lambda.min)

## Plotting the most important variables coefficients of the LASSO model
coef(Lasso.neuroEphys, s=Lasso.neuroEphys$lambda.min) %>%
  broom::tidy() %>%
  filter(row != "(Intercept)") %>%
  ggplot(aes(value, reorder(row, value), color = value > 0)) +
  geom_point(show.legend = FALSE) +
  ggtitle("Variable Importance in LASSO") +
  xlab("Coefficient") +
  ylab(NULL)

coef(Lasso.neuroEphys, s=Lasso.neuroEphys$lambda.min) %>%
  broom::tidy() %>%
  filter(row != "(Intercept)") %>%
  ggplot(aes(value, reorder(row, value), color = value > 0)) +
  geom_point(show.legend = FALSE) +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=14,face="bold")) +
  ggtitle("Variable Importance in LASSO of class split for AHPamp inh/exc") +
  xlab("Coefficient") +
  ylab("Gene Name")


# Splitting data into Inhibitory and Excitatory splits
neuroEphysData4INH <- neuroEphysData4[grep("*inh",rownames(neuroEphysData4)),]
neuroEphysData4EXC <- neuroEphysData4[grep("*exc",rownames(neuroEphysData4)),]


set.seed(1235)
## Declaring variables for cell class
nrfolds <- nrow(neuroEphysData4INH)
nEphysFeatures <- 16

folds <- rep_len(1:nrfolds, nrow(neuroEphysData4INH))
folds <- sample(folds, nrow(neuroEphysData4INH))

r2 = list()
save_df = list()
save_df2 = list()
save_df3 = list()
r2 = list()
shuffledr2 = list()
save_values = list()
pred_df=list()
class=list()
FinalFinalDF=list()
FinalFinalDF2=list()
EphysFeature=list()
hyperParameter=list()
MSEValue=list()
ngenes=list()
rss1=list()
tss1=list()
rsq1=list()
rsq5=list()

# Inhibitory LOOCV
for(x in 1:nEphysFeatures){
  for(i in 1:nrfolds){
    ind = which(folds != i)
    
    train.neuroEphys <- neuroEphysData4INH[ind, ]
    test.neuroEphys <- neuroEphysData4INH[-ind, ]
    
    trainingY <- as.matrix(train.neuroEphys[,x])
    testingY <- as.matrix(test.neuroEphys[,x])
    #shuffledY <- as.matrix(train.neuroEphys[sample(nrow(trainingY)),x])
    
    trainingX <- as.matrix(train.neuroEphys[,29:45796])
    testingX <- as.matrix(test.neuroEphys[,29:45796])
    
    Lasso.neuroEphys <- cv.glmnet(trainingX, trainingY, alpha = 1, type.measure="mse", family="gaussian")
    #ShuffledLasso.neuroEphys <- cv.glmnet(trainingX, shuffledY, alpha = 1, type.measure="mse", family="gaussian")
    
    fit <- glmnet(trainingX,
                  trainingY,
                  alpha = 1,
                  lambda = Lasso.neuroEphys$lambda.min,
                  family="gaussian")
    
    #shuffledFit <- glmnet(trainingX,
    #                      shuffledY,
    #                      alpha = 1,
    #                      lambda = ShuffledLasso.neuroEphys$lambda.min,
    #                      family="gaussian")
    
    prds.test <- predict(fit,newx = as.matrix(testingX), type = "response", s=Lasso.neuroEphys$lambda.min)
    #Shuffledprds.test <- predict(shuffledFit,newx = as.matrix(testingX), type = "response", s=ShuffledLasso.neuroEphys$lambda.min)
    ##(cbind(prds.test, testingY))
    
    save_df[[i]] <- data.frame(prds.test, testingY)
    
    pred_df[i] <- prds.test
    
    save_df2[i] <- testingY
    
    #shuffledPred_df[i] <- Shuffledprds.test
    
    #shuffled_values_df <- data.frame(Shuffledprds.test, testingY)
    
    #rss[i] <- sum((prds.test - testingY) ^ 2)  ## residual sum of squares
    #tss[i] <- sum((testingY - mean(testingY)) ^ 2)  ## total sum of squares
    #rssVal <- as.numeric(rss[i])
    #tssVal <- as.numeric(tss[i])
    #rsq[i] <- 1 - (rssVal/tssVal)
    
    #rssShuffled[i] <- sum((Shuffledprds.test - testingY) ^ 2)  ## residual sum of squares
    #tssShuffled[i] <- sum((testingY - mean(testingY)) ^ 2)  ## total sum of squares
    #rssValShuffled <- as.numeric(rssShuffled[i])
    #tssValShuffled <- as.numeric(tssShuffled[i])
    #rsqShuffled[i] <- 1 - (rssValShuffled/tssValShuffled)
  }
  EphysFeature[x] <- colnames(train.neuroEphys)[x]
  
  #rsq4[[x]] <- data.frame(rsq, EphysFeature[x])
  #rsq4shuffled[[x]] <- data.frame(rsqShuffled, EphysFeature[x])
  
  #LOOCV Code:
  pred_df <- unlist(pred_df)
  
  save_df2 <- as.numeric(unlist(t(save_df2)))
  r2[x] <- summary(lm(pred_df ~ save_df2))$r.squared
  r2 <- summary(lm(pred_df ~ save_df2))$r.squared
  
  # Generates an r2 value of 0.888292 for AHPamp
  # Generates an r2 value of 0.654871 for AHPw
  rss[x] <- sum((pred_df - save_df2) ^ 2)  ## residual sum of squares
  tss[x] <- sum((save_df2 - mean(save_df2)) ^ 2)  ## total sum of squares
  rssVal <- as.numeric(rss[x])
  tssVal <- as.numeric(tss[x])
  rsq[x] <- 1 - (rssVal/tssVal)
  # Generates an r2 value of 0.872248 for AHPamp
  # Generates an r2 value of 0.647678 for AHPw
  #shuffledPred_df <- unlist(shuffledPred_df)
  
  #shuffledr2[x] <- summary(lm(shuffledPred_df ~ save_df2))$r.squared
  
  #rssShuffled[x] <- sum((shuffledPred_df - save_df2) ^ 2)  ## residual sum of squares
  #tssShuffled[x] <- sum((save_df2 - mean(save_df2)) ^ 2)  ## total sum of squares
  #rssValShuffled <- as.numeric(rssShuffled[x])
  #tssValShuffled <- as.numeric(tssShuffled[x])
  #rsqShuffled[x] <- 1 - (rssValShuffled/tssValShuffled)
  
}

rsq5DF <- unlist(rsq)
EphysFeature <- unlist(EphysFeature)
rsq5DF <- data.frame(rsq5DF, EphysFeature)

ggplot(rsq5DF, aes(EphysFeature, rsq5DF)) +
  geom_boxplot() +
  ggtitle('Predicting EphysFeature Correlation Coefficients using Gene expression Boxplot for only Inhibitory Cells via LOOCV') +
  scale_y_continuous() +
  theme(axis.text=element_text(size=14), axis.text.x=element_text(angle=45), axis.title=element_text(size=14,face="bold"), legend.title=element_text(size=14), legend.text=element_text(size=14)) +
  ylab('R^2 (LOOCV)') +
  xlab('EphysFeature')


set.seed(1235)
## Declaring variables for cell class
nrfolds <- nrow(neuroEphysData4EXC)
nEphysFeatures <- 16

folds <- rep_len(1:nrfolds, nrow(neuroEphysData4EXC))
folds <- sample(folds, nrow(neuroEphysData4EXC))

r2 = list()
save_df = list()
save_df2 = list()
save_df3 = list()
r2 = list()
shuffledr2 = list()
save_values = list()
pred_df=list()
class=list()
FinalFinalDF=list()
FinalFinalDF2=list()
EphysFeature=list()
hyperParameter=list()
MSEValue=list()
ngenes=list()
rss1=list()
tss1=list()
rsq1=list()
rsq5=list()

# Excitatory LOOCV
for(x in 1:nEphysFeatures){
  for(i in 1:nrfolds){
    ind = which(folds != i)
    
    train.neuroEphys <- neuroEphysData4EXC[ind, ]
    test.neuroEphys <- neuroEphysData4EXC[-ind, ]
    
    trainingY <- as.matrix(train.neuroEphys[,x])
    testingY <- as.matrix(test.neuroEphys[,x])
    #shuffledY <- as.matrix(train.neuroEphys[sample(nrow(trainingY)),x])
    
    trainingX <- as.matrix(train.neuroEphys[,29:45796])
    testingX <- as.matrix(test.neuroEphys[,29:45796])
    
    Lasso.neuroEphys <- cv.glmnet(trainingX, trainingY, alpha = 1, type.measure="mse", family="gaussian")
    #ShuffledLasso.neuroEphys <- cv.glmnet(trainingX, shuffledY, alpha = 1, type.measure="mse", family="gaussian")
    
    fit <- glmnet(trainingX,
                  trainingY,
                  alpha = 1,
                  lambda = Lasso.neuroEphys$lambda.min,
                  family="gaussian")
    
    #shuffledFit <- glmnet(trainingX,
    #                      shuffledY,
    #                      alpha = 1,
    #                      lambda = ShuffledLasso.neuroEphys$lambda.min,
    #                      family="gaussian")
    
    prds.test <- predict(fit,newx = as.matrix(testingX), type = "response", s=Lasso.neuroEphys$lambda.min)
    #Shuffledprds.test <- predict(shuffledFit,newx = as.matrix(testingX), type = "response", s=ShuffledLasso.neuroEphys$lambda.min)
    ##(cbind(prds.test, testingY))
    
    save_df[[i]] <- data.frame(prds.test, testingY)
    
    pred_df[i] <- prds.test
    
    save_df2[i] <- testingY
    
    #shuffledPred_df[i] <- Shuffledprds.test
    
    #shuffled_values_df <- data.frame(Shuffledprds.test, testingY)
    
    #rss[i] <- sum((prds.test - testingY) ^ 2)  ## residual sum of squares
    #tss[i] <- sum((testingY - mean(testingY)) ^ 2)  ## total sum of squares
    #rssVal <- as.numeric(rss[i])
    #tssVal <- as.numeric(tss[i])
    #rsq[i] <- 1 - (rssVal/tssVal)
    
    #rssShuffled[i] <- sum((Shuffledprds.test - testingY) ^ 2)  ## residual sum of squares
    #tssShuffled[i] <- sum((testingY - mean(testingY)) ^ 2)  ## total sum of squares
    #rssValShuffled <- as.numeric(rssShuffled[i])
    #tssValShuffled <- as.numeric(tssShuffled[i])
    #rsqShuffled[i] <- 1 - (rssValShuffled/tssValShuffled)
  }
  EphysFeature[x] <- colnames(train.neuroEphys)[x]
  
  #rsq4[[x]] <- data.frame(rsq, EphysFeature[x])
  #rsq4shuffled[[x]] <- data.frame(rsqShuffled, EphysFeature[x])
  
  #LOOCV Code:
  pred_df <- unlist(pred_df)
  
  save_df2 <- as.numeric(unlist(t(save_df2)))
  r2[x] <- summary(lm(pred_df ~ save_df2))$r.squared
  r2 <- summary(lm(pred_df ~ save_df2))$r.squared
  
  # Generates an r2 value of 0.888292 for AHPamp
  # Generates an r2 value of 0.654871 for AHPw
  rss[x] <- sum((pred_df - save_df2) ^ 2)  ## residual sum of squares
  tss[x] <- sum((save_df2 - mean(save_df2)) ^ 2)  ## total sum of squares
  rssVal <- as.numeric(rss[x])
  tssVal <- as.numeric(tss[x])
  rsq[x] <- 1 - (rssVal/tssVal)
  # Generates an r2 value of 0.872248 for AHPamp
  # Generates an r2 value of 0.647678 for AHPw
  #shuffledPred_df <- unlist(shuffledPred_df)
  
  #shuffledr2[x] <- summary(lm(shuffledPred_df ~ save_df2))$r.squared
  
  #rssShuffled[x] <- sum((shuffledPred_df - save_df2) ^ 2)  ## residual sum of squares
  #tssShuffled[x] <- sum((save_df2 - mean(save_df2)) ^ 2)  ## total sum of squares
  #rssValShuffled <- as.numeric(rssShuffled[x])
  #tssValShuffled <- as.numeric(tssShuffled[x])
  #rsqShuffled[x] <- 1 - (rssValShuffled/tssValShuffled)
  
}

rsq6DF <- unlist(rsq)
EphysFeature <- unlist(EphysFeature)
rsq6DF <- data.frame(rsq6DF, EphysFeature)

ggplot(rsq6DF, aes(EphysFeature, rsq6DF)) +
  geom_boxplot() +
  ggtitle('Predicting EphysFeature Correlation Coefficients using Gene expression Boxplot for only Excitatory Cells via LOOCV') +
  scale_y_continuous() +
  theme(axis.text=element_text(size=14), axis.text.x=element_text(angle=45), axis.title=element_text(size=14,face="bold"), legend.title=element_text(size=14), legend.text=element_text(size=14)) +
  ylab('R^2 (LOOCV)') +
  xlab('EphysFeature')


### Rough Notes for code
'''MSE_LOOCV <- cv.glm(train.neuroEphys, model1)
MSE_LOOCV$delta[1]

MSE_LOOCV = NULL

for (i in 1:5){
  model = glm(trainingX[,4] ~ poly(trainingY[,29:50], i), data=train.neuroEphys)
  MSE_LOOCV[i] <- cv.glm(train.neuroEphys, model1)$delta[1]
}



cv.lm(data = neuroEphysData2, form.lm = formula(ElectroPhysio ~ GeneExpression),
      m = 3, dots = FALSE, seed = 29, plotit = c("Observed","Residual"),
      main="Small symbols show cross-validation predicted values",
      legend.pos="topleft", printit = TRUE)

plot(Lasso.neuroEphys)

## PLOTTING ##

names(neuroEphysData2) %>% strsplit(., '__')

Inh_Exc_Classification <- grep("inh", names(neuroEphysData2), value=TRUE)

sample_name_df = names(neuroEphysData2) %>% as.data.frame()
colnames(sample_name_df) = 'sample_name'

sample_meta = sample_name_df %>% separate('sample_name', c('cre', 'layer', 'class'), sep = '__')

Gprasp1_transposed <- t(Gprasp1)
Camk2g_transposed <- t(Camk2g)
Xxylt1_transposed <- t(Xxylt1)
AHP_Amplitude_transposed <- t(AHP_Amplitude)

Ggprasp1.df <- as.data.frame(Gprasp1_transposed)
Camk2g.df <- as.data.frame(Camk2g_transposed)
Xxylt1.df <- as.data.frame(Xxylt1_transposed)
AHP_Amplitude.df <- as.data.frame(AHP_Amplitude_transposed)

sample_meta['Ggprasp1 (log2 CPM+1)']= Ggprasp1.df[Ggrasp_entrez_id]
sample_meta['Camk2g_(log2_CPM+1)']= Camk2g.df[Camk2g_entrez_id]
sample_meta['Xxylt1_(log2_CPM+1)']= Xxylt1.df[Xxylt1_entrez_id]
sample_meta['AHP Amplitude (mV)']= AHP_Amplitude.df['ahpamp']

sample_meta['Ggprasp1_(log2_CPM+1)']=sample_meta['Ggprasp1 (log2 CPM+1)']
sample_meta['AHP_Amplitude_(mV)']= sample_meta['AHP Amplitude (mV)']


ggplot(sample_meta, aes(x=`Ggprasp1_(log2_CPM+1)`, y=`AHP_Amplitude_(mV)`)) +
  geom_point(aes(shape = class, color = class), size = 3) +
  scale_shape_manual(values = c(18, 16)) +
  scale_color_manual(values = c("#00AFBB", "#E35334")) + 
  geom_smooth(method = "lm") +
  ggtitle('Class-driven') +
  xlab('Gprasp1 (log2 CPM+1)') +
  ylab('AHP Amplitude (mV)')


ggplot(sample_meta, aes(x=`Camk2g_(log2_CPM+1)`, y=`AHP_Amplitude_(mV)`)) +
  geom_point(aes(shape = class, color = class), size = 3) +
  scale_shape_manual(values = c(18, 16)) +
  scale_color_manual(values = c("#00AFBB", "#E35334")) + 
  geom_smooth(method = "lm") +
  ggtitle('Non-class-driven') +
  xlab('Camk2g (log2 CPM+1)') +
  ylab('AHP Amplitude (mV)')


ggplot(sample_meta, aes(x=`Xxylt1_(log2_CPM+1)`, y=`AHP_Amplitude_(mV)`)) +
  geom_point(aes(shape = class, color = class), size = 3) +
  scale_shape_manual(values = c(18, 16)) +
  scale_color_manual(values = c("#00AFBB", "#E35334")) +
  geom_smooth(method = "lm") +
  ggtitle('Non-class-driven; sig. in both models') +
  xlab('Xxylt1 (log2 CPM+1)') +
  ylab('AHP Amplitude (mV)')

'''
