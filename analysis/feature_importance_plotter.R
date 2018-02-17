# DSSL imports
source("dssFunctionLibrary.R")
source("dssPerformanceEvaluation.R")
source("mimicUsefulFunction.R")
source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/my.colorFct.R")
source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/dendroCol.R")

# imports from CRAN
require(gridExtra)
require(ROCR)
require(caret)
require(chron)
require(magrittr)
require(chron)
require(proxy)
require(ggplot2)
require(lubridate)
require(dplyr)
library(tidyverse)
require(hash)
require(outliers)
require(RSQLite)
require(sqldf)
require(data.table)
require(gtools)

############################################################################################################
# Determine which features DSSL weights the most across active test modalities when learning the mPDS 
# Included here is code to visualize feature importance
# Note that the data in this file (sorted_weight.csv) is out of date
# Refer to "Linear DSSL-v1-tap2-robust-10fcv-pt-dssv2-mPDS-short.ipynb" for up-to-date results
############################################################################################################

# read sorted weights on features from learned model
sorted_weights <- read.csv("sorted_weight.csv", header = FALSE, sep = ",", quote = "\"", dec = ".", fill = TRUE,
                           comment.char = "", na.strings = c("NA", "#DIV/0!", ""))
colnames(sorted_weights) <- c('test','weight')
sorted_weights$test <- as.character(sorted_weights$test)

# normalize weights to sum to percentages (to simplify relative comparison)
sorted_weights$weight <- sorted_weights$weight / 100
sorted_weights <- sorted_weights[which(sorted_weights$test != 'date'),]

# extract top 10 most informative features in DSSL
best10_features <- sorted_weights[1:10,1]
# extract 10 least informative features in DSSL
worst10_features <- sorted_weights[(nrow(sorted_weights)-10):nrow(sorted_weights),1]
# extract 10 features of "medium importance" (to round out the visualization)
mid10_features <- sorted_weights[20:30,1]

# extract weights associated with 10 most informative features
best10_weights <- sorted_weights[1:10,2]
# extract weights associated with 10 least informative features
worst10_weights <- sorted_weights[(nrow(sorted_weights)-10):nrow(sorted_weights),2]
# extract weights associated with "medium importance" features
mid10_weights <- sorted_weights[20:30,2]

# plotting code
x_features <- c(best10_features,mid10_features,worst10_features)
y_weights <- c(best10_weights,mid10_weights,worst10_weights)
x_labels <- 1:length(y_weights)
plot(x_labels,y_weights,xaxt='n',xlab='',ylab='Normalized Feature Weights - DSSL',
     pch=19, col=ifelse(x_labels < 10, "darkgreen", ifelse(x_labels > (length(x_labels)-10),
                                                           "red","gold")))
# plotting code
x_features <- c(as.character(best10_features),as.character(mid10_features),
                as.character(worst10_features))
axis(1, at=1:length(y_weights), labels=x_features, las=2, cex.axis=0.7, xlab='features')
title('HopkinsPD Feature Importance in Computing PDSS')
par(oma=c(0,3,0,0))

# plotting code
x_labels <- length(y_weights):1
plot(y_weights,x_labels,yaxt='n',xlab='Normalized Feature Weights - DSSL', ylab="",
     pch=19, col=ifelse(x_labels < 10, "red", ifelse(x_labels > (length(x_labels)-10),
                                                     "darkgreen","gold")))
x_features <- c(as.character(best10_features),as.character(mid10_features),
                as.character(worst10_features))
axis(2, at=length(y_weights):1, labels=x_features, las=2, cex.axis=0.6, ylab='features')
title('HopkinsPD Feature Importance in Computing PDSS')
legend("topleft", inset = .02, title = "Feature Importance",
       legend=c('10 worst','10 moderate','10 best'),
       fill = c('red','gold','green'), horiz=F)

# extract features and associated weights on tapping active test
tap.wgts <- wgtsDF[which(wgtsDF$colnames.wgtCopy. %in% colnames(roch_tap_1)),'weights']
# extract features and associated weights on postural instability active test
bal.wgts <- wgtsDF[which(wgtsDF$colnames.wgtCopy. %in% gsub(" ","",paste('balance_',colnames(roch_balance_1)))),'weights']
# extract features and associated weights on gait active test
gait.wgts <- wgtsDF[which(wgtsDF$colnames.wgtCopy. %in% gsub(" ","",paste('gait_',colnames(roch_gait_1)))),'weights']
# extract features and associated weights on pitch active test
aud.wgts <- wgtsDF[which(wgtsDF$colnames.wgtCopy. %in% colnames(roch_audio_1)),'weights']
# extract features and associated weights on reaction time active test
rct.wgts <- wgtsDF[which(wgtsDF$colnames.wgtCopy. %in% colnames(roch_react_1)),'weights']

# normalize weights to assess what "percentage of total weights" learned by DSSL was assigned to each modality
tot_sum <- sum(tap.wgts) + sum(bal.wgts) + sum(gait.wgts) + sum(aud.wgts) + sum(rct.wgts)
tap.con <- sum(tap.wgts) / tot_sum
bal.con <- sum(bal.wgts) / tot_sum
gait.con <- sum(gait.wgts) / tot_sum
aud.con <- sum(aud.wgts) / tot_sum
rct.con <- sum(rct.wgts) / tot_sum

# plotting code for legend
legend("bottomright", inset=.02, title = "Test Contribution",
       legend=c(paste('finger dexterity:',round(tap.con,3)),paste('postural instability:',round(bal.con,3)),
                paste('gait:',round(gait.con,3)),paste('voice',round(aud.con,3)),
                paste('reaction time',round(rct.con,3))),
       fill = c('grey','grey','grey','grey','grey'), horiz=FALSE)