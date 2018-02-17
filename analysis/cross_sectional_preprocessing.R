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

############################################################################################################
# Cross-Sectional Validation of mPDS Against Traditional Rating Scales
# Read in learned mPDS data and perform preprocessing to get data in format for analysis
# The actual cross-sectional analysis is performed in the mPDS_cross_sectional_analysis-documented.ipynb
# notebook using the data generated in the two .csv files at the end of this script
############################################################################################################

# Read file containing data associated with traditional rating scales
updrs_new <- read.csv("udall_assessments_1231.csv", header = TRUE, sep = ",", quote = "\"", dec = ".", fill = TRUE,
                      comment.char = "", na.strings = c("NA", "#DIV/0!", "")) # read UPDRS data for cross-sectional analysis
updrs_new$pd_yn <- as.character(updrs_new$pd_yn)

colnames(updrs_new)[colnames(updrs_new) == 'rec_id'] <- 'Record.ID'
colnames(updrs_new)[colnames(updrs_new) == 'event'] <- 'Event.Name'
colnames(updrs_new)[colnames(updrs_new) == 'timed_up_and_go'] <- 'Timed.Up.and.Go.time..seconds.'
colnames(updrs_new)[colnames(updrs_new) == 'updrs_part3'] <- 'UPDRS.Part.3..Motor.Exam..clinician.rated..Score.'
colnames(updrs_new)[colnames(updrs_new) == 'updrs_total'] <- 'UPDRS.Total.Score'
colnames(updrs_new)[colnames(updrs_new) == 'hy'] <- 'Estimated.Hoehn.and.Yahr.stage'
updrs_new <- updrs_new[which(updrs_new$pd_yn == 'Yes'),]

# read file containing additional data associated with the traditional rating scales
all <- read_csv("AllData_01.02.18.csv")
all <- as.data.frame(all)
all <- all[which(all$`Do you have Parkinson disease?` == 'Yes' | is.na(all$`Do you have Parkinson disease?`)),]
t6 <- which(!is.na(all$`Timed Up and Go (seconds) [6-MONTH ON-STATE ONLY]`))
all[t6,'Timed Up and Go time (seconds)'] <- all$`Timed Up and Go (seconds) [6-MONTH ON-STATE ONLY]`[t6]
p3 <- which(!is.na(all$`UPDRS Part 3, Motor Exam 'ON' State Score`))
all[p3,'UPDRS Part 3, Motor Exam (clinician-rated) Score'] <- all$`UPDRS Part 3, Motor Exam 'ON' State Score`[p3]
hy <- which(!is.na(all$`Estimated Hoehn and Yahr stage 'ON' State Score`))
all[hy,'Estimated Hoehn and Yahr stage'] <- all$`Estimated Hoehn and Yahr stage 'ON' State Score`[hy]
all$`Timed Up and Go (seconds) [6-MONTH ON-STATE ONLY]` <- NULL
all$`Estimated Hoehn and Yahr stage 'ON' State Score` <- NULL
all$`UPDRS Part 3, Motor Exam 'ON' State Score` <- NULL

# read file containing learned mPDS scores
tap2 <- read.csv('udall_scores_1231.csv', header = TRUE, sep = ",", quote = "\"", dec = ".", fill = TRUE,
                 comment.char = "", na.strings = c("NA", "#DIV/0!", ""))
colnames(tap2)[colnames(tap2) == 'scores_scaled'] <- 'pdss'

updated_id <-
  sapply(updrs_new$Record.ID,function(x){
    div_cnt <- 0; tmp_x <- x
    while(tmp_x > 0) {
      div_cnt = div_cnt + 1
      tmp_x = as.integer(tmp_x / 10)
    }
    padded <- paste(rep("0",times=(3-div_cnt)),collapse='')
    temp <- gsub(" ","",paste("udall",padded,x))
  })

updrs_new$Record.ID <- updated_id
updrs_new <- updrs_new[which(updrs_new$Record.ID %in% tap2$id),] # filter only updrs ids on which we can compare

# groupby most afflicted side
sides <- tap2 %>% group_by(id,test) %>% summarise(pdss=mean(pdss))
sides.sub <- sides %>% group_by(id) %>% summarise(pdss=max(pdss))
sides.sub$test <- NA
sub.ids <- as.character(sides.sub$id)
sub.pdss <- as.numeric(sides.sub$pdss)
sides.test <- as.character(sides$test)

sub.test <- c()
lapply(1:nrow(sides.sub),function(x){
  test <- sides.test[which(sides$id == sub.ids[x] & sides$pdss == sub.pdss[x])]
  sub.test <<- c(sub.test,test)
})

sides.sub$test <- sub.test
sides.sub <- data.frame(sides.sub)

# match updrs_new to tap2 on timestamp and remove NAs
set.seed(2536)
timestamps <- read.csv("TIMESTAMP_12.7.17.csv", header = TRUE, sep = ",", quote = "\"", dec = ".", fill = TRUE,
                       comment.char = "", na.strings = c("NA", "#DIV/0!", "")) # timestamps to match PDSS to UPDRS
colnames(timestamps)[colnames(timestamps) == 'Smartphone.Activity_time'] <- 'tsp'
timestamps$tsp <- as.character(timestamps$tsp)
timestamps$tsp <- gsub(" ","",paste(timestamps$tsp,':00'))
timestamps$tsp <- chron(times=timestamps$tsp)
timestamps$Motor.Exam_Time <- NULL
timestamps <- timestamps[complete.cases(timestamps),]
updated_id <- sapply(timestamps$record,function(x){
  div_cnt <- 0; tmp_x <- x
  while(tmp_x > 0) {
    div_cnt = div_cnt + 1
    tmp_x = as.integer(tmp_x / 10)
  }
  padded <- paste(rep("0",times=(3-div_cnt)),collapse='')
  temp <- gsub(" ","",paste("udall",padded,x))
})

timestamps$record <- updated_id
tsps <- as.character(tap2$tsp)
dates <- NULL; times <- NULL

tsps <- sapply(tsps,function(x){
  res <- strsplit(x, " ", fixed=F)[[1]]
  dates <<- c(dates,res[[1]])
  times <<- c(times,res[[2]])
})

timestamps$event_date <- as.character(timestamps$event_date)
timestamps$event_date <- as.Date(timestamps$event_date, "%m/%d/%y")
times <- chron(times=times)
pdss <- tap2$pdss; tap2$pdss <- NULL; tap2$tsp <- NULL
tap2$event_date <- dates; tap2$event_time <- times
tap2$pdss <- pdss
tap2$event_date <- as.Date(tap2$event_date, "%Y-%m-%d")

tap2.filt <- tap2[which(tap2$event_date %in% timestamps$event_date),]
updrs_new$event_date <- NA; updrs_new$event_time <- NA
updrs_new$Event.Name <- as.character(updrs_new$Event.Name)
timestamps$descrip <- as.character(timestamps$descrip)
updrs.count <- nrow(updrs_new)
all_id <- sapply(all$record_id,function(x){
  div_cnt <- 0; tmp_x <- x
  while(tmp_x > 0) {
    div_cnt = div_cnt + 1
    tmp_x = as.integer(tmp_x / 10)
  }
  padded <- paste(rep("0",times=(3-div_cnt)),collapse='')
  temp <- gsub(" ","",paste("udall",padded,x))
})

all$record_id <- all_id
all <- all[which(all[['Event Name']] %in% c('Baseline (enrollment)','Month 3','Month 6')),]
colnames(all)[colnames(all) == "Estimated Hoehn and Yahr stage"] <- 'Estimated.Hoehn.and.Yahr.stage'
colnames(all)[colnames(all) == "Timed Up and Go time (seconds)"] <- 'Timed.Up.and.Go.time..seconds.'
colnames(all)[colnames(all) == "UPDRS Total Score"] <- 'UPDRS.Total.Score'
colnames(all)[colnames(all) == "UPDRS Part 3, Motor Exam (clinician-rated) Score"] <- "UPDRS.Part.3..Motor.Exam..clinician.rated..Score."
colnames(all)[colnames(all) == 'record_id'] <- "Record.ID"
colnames(all)[colnames(all) == "Event Name"] <- "Event.Name"
res_all <- smartbind(all,updrs_new)
res_all <- res_all[,colnames(res_all) %in% colnames(updrs_new)]
rownames(res_all) <- NULL
updrs_new <- res_all

for(ind in 1:updrs.count){
  curr.desc <- updrs_new[ind,'Event.Name']
  curr.record <- updrs_new[ind,'Record.ID']
  if(curr.record %in% timestamps$record){
    record.sub <- timestamps[which(timestamps$record == curr.record),]
    if(curr.desc %in% record.sub$descrip){
      updrs_new[ind,'event_date'] <- as.Date(record.sub[which(record.sub$descrip == curr.desc),'event_date'])
      updrs_new[ind,'event_time'] <- as.character(record.sub[which(record.sub$descrip == curr.desc),'tsp'])
    }
  }
}

updrs_new <- updrs_new[-which(is.na(updrs_new$event_date) | is.na(updrs_new$event_time)),]
updrs_new$event_time <- chron(times=updrs_new$event_time)
updrs_new$event_date <- as.Date(updrs_new$event_date)
colnames(updrs_new)[colnames(updrs_new) == 'UPDRS.Part.3..Motor.Exam..clinician.rated..Score.'] <- 'updrs_part3'
colnames(updrs_new)[colnames(updrs_new) == 'Timed.Up.and.Go.time..seconds.'] <- 'timed_up_and_go'
colnames(updrs_new)[colnames(updrs_new) == 'UPDRS.Total.Score'] <- 'updrs_total_score'
colnames(updrs_new)[colnames(updrs_new) == 'Estimated.Hoehn.and.Yahr.stage'] <- 'hoehn_yahr'
tap2.filt <- tap2.filt[which(tap2.filt$event_date %in% updrs_new$event_date),]
tap2.filt$event_time <- as.character(tap2.filt$event_time)
updrs_new$event_time <- as.character(updrs_new$event_time)
updrs_new$event_date <- as.character(updrs_new$event_date)
write.csv(tap2.filt,'__tap2_filt_openSource.csv',row.names=F)
write.csv(updrs_new,'__updrs_new_openSource.csv',row.names=F)

# Pass these two data files to the "mPDS_cross_sectional_analysis-documented.ipynb" notebook
# to run the cross-sectional analysis
