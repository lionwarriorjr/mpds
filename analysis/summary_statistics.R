source("dssFunctionLibrary.R")
source("dssPerformanceEvaluation.R")
source("mimicUsefulFunction.R")
source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/my.colorFct.R")
source("http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/dendroCol.R")
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

###############################################################################
# Summary Statistics
###############################################################################

######################################################################
# Analysis of mPDS in individuals with PD vs. without PD
######################################################################
updrs_new <- read.csv("udall_assessments_1231.csv", header = TRUE, sep = ",", quote = "\"", dec = ".", fill = TRUE,
                      comment.char = "", na.strings = c("NA", "#DIV/0!", "")) # read UPDRS data for cross-sectional analysis
updrs_new$pd_yn <- as.character(updrs_new$pd_yn)
updrs_new_id <- sapply(updrs_new$rec_id,function(x){
  div_cnt <- 0; tmp_x <- x
  while(tmp_x > 0) {
    div_cnt = div_cnt + 1
    tmp_x = as.integer(tmp_x / 10)
  }
  padded <- paste(rep("0",times=(3-div_cnt)),collapse='')
  temp <- gsub(" ","",paste("udall",padded,x))
})
updrs_new$rec_id <- updrs_new_id
newScore <- read.csv('udall_scores_1231.csv')
newScore$id <- as.character(newScore$id)
yes_id <- updrs_new[which(updrs_new$pd_yn == 'Yes'),'rec_id']
no_id <- updrs_new[which(updrs_new$pd_yn == 'No'),'rec_id']
yes_scores <- newScore[which(newScore$id %in% yes_id),'scores_scaled']
no_scores <- newScore[which(newScore$id %in% no_id),'scores_scaled']

# return mean mPDS score in patients with PD
mean(yes_scores)
# return mean mPDS score in patients without PD
mean(no_scores)
#return standard deviation in mPDS for patients with PD
sd(yes_scores)
# return standard deviation in mPDS for patients without PD
sd(no_scores)

#############################################################################################################
# Analysis of intraday change in mPDS as measure of mPDS's sensitivity to daily symptom fluctuations
#############################################################################################################

### clinic cohort ###
tap2 <- read.csv('udall_scores_1231.csv')
tap2 <- data.frame(tap2)
tap2 <- tap2[order(tap2$tsp),]
tap2$id <- as.character(tap2$id)
tap2 <- tap2[which(tap2$id %in% yes_id),]

unique_ids <- as.character(unique(tap2$id))
for(udall.id in unique_ids){
  write.csv(tap2[which(tap2$id == udall.id),],gsub(" ","",paste(udall.id,'.csv')),row.names=F)
}

intraday <- intraday.ids <- NULL
count_assessments <- NULL
suppressWarnings(
  lapply(1:length(unique_ids),function(udall.index){
    tryCatch({
      udall.id <- unique_ids[udall.index]
      plot_dt <- read.csv(gsub(" ","",paste(udall.id,'.csv')))
      j <- 'scores_scaled'
      if(nrow(plot_dt) %% 2 == 1) plot_dt <- plot_dt[1:(nrow(plot_dt)-1),]
      plot_dt$tsp <- as.character(plot_dt$tsp)
      event_date <- strsplit(plot_dt$tsp," ")
      event_time <- unlist(lapply(event_date,function(x){x[[2]]}))
      event_date <- unlist(lapply(event_date,function(x){x[[1]]}))
      plot_dt$event_date <- as.factor(event_date)
      plot_dt$event_time <- as.factor(event_time)
      intraday.diffs <- sqldf('select event_date, count(*) count, max(scores_scaled)-min(scores_scaled) intraday_change 
                              from plot_dt group by event_date having count(*) > 1')
      intraday.ids <<- c(intraday.ids,udall.id)
      intraday <<- c(intraday,intraday.diffs$intraday_change)
      count_assessments <<- c(count_assessments, sum(intraday.diffs$count))
    } , error=function(e){
      cat("Too few data points\n")
    })
  })
)

# maximum number of mPDS assessments in clinic cohort
max(count_assessments)
# minimum number of mPDS assessments in clinic cohort 
min(count_assessments)
# mean intraday change in mPDS (as measurement of daily severity fluctuations) in clinic cohort
mean(intraday) #13.87762
# standard deviation in change in mPDS (as measurement of daily severity fluctuations) in clinic cohort
sd(intraday) #10.33258

############################################################################################################
# Run Wilcoxon signed-rank test to assess whether mPDS detects a decrease in severity
# in response to dopaminergic therapy in the clinic cohort
############################################################################################################

# read score differences associated with adminstration of dopaminergic therapy in clinic cohort
diff <- read_tsv("differences.txt")
diff <- data.frame(diff)
tot_diffs <- diff$score_diff

# run Wilcoxon signed-rank test to assess whether mPDS is generally lower after therapy
shapiro.test(tot_diffs)
ds <- wilcox.test(tot_diffs, alternative="greater")
# test statistic reported in the main text
ds$statistic
# associated p-value reported in the main text
ds$p.value