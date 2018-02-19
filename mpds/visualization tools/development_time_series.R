# DSSL imports
source("dssFunctionLibrary.R")
source("dssPerformanceEvaluation.R")
source("mimicUsefulFunction.R")

# imports from CRAN
require(gridExtra)
require(ROCR)
require(ISLR)
require(caret)
require(magrittr)
require(chron)
require(randomcoloR)
require(proxy)
require(tidyverse)
library(ggplot2)
require(plyr)

############################################################################################################
# Generate plots of efficacy of mPDS in tracking intraday fluctuation in the development set 
# In the development set, we have access to the actual on-state and off-state timestamps associated
# with each smartphone assessment as outlined in the main text
# This code is used to generate intraday fluctuation time series reported on the development set
# in the accompanying paper
############################################################################################################

# read in score matrix
#     id                 tsp scores_react scores_voice scores_gait  scores_scaled
#1  0022 2014-07-11 05:58:34  0.056703681    0.5807599   2.7295750       83.49790
#2  0022 2014-07-11 07:10:05  0.066239872   -0.8723297   0.2115280       68.49464
#3  0035 2014-07-12 05:38:39 -0.008428049    0.8452914   2.6596697       62.04420
#4  0022 2014-07-12 06:00:48 -0.043301346    0.3299085   0.5734744       64.10186
#5  0039 2014-07-12 07:10:18 -0.097928252   -2.1624151   2.5494258       46.20457
#6  0035 2014-07-12 07:18:18  0.011130728    0.7006828   2.7620053       61.39492

# this matrix simply contains the relevant ids, timestamps of collected measurements,
# and the learned scores by DSSL. Here we stored the learned scores for
# a subset of the smartphone modalities alongside the final, scaled mPDS (scores_scaled)
dt <- read_csv("roch_scores_1231.csv")
dt <- data.frame(dt)
dt <- dt[order(dt$tsp),]
patients <- dt$id %>% unique
patient_ids <- patients
tests = c( "scores_react", "scores_voice", "scores_balance",
           "scores_tap", "scores_gait", "scores", "scores_scaled")
mean_inc <- c()

# pick a score to visualize (here we visualize scores_scaled, i.e. mPDS)
j = tests[7]

# output generate plots to PDF
filename = gsub(" ","",paste("trajectories_PDSS_scaled.pdf"), fixed=T); print(filename) 
pdss_inc <- NULL; abs_diff <- NULL; tot_diffs <- NULL; durations <- NULL
pdf(file=filename, height=2.155, width=6.974)

for(i in 1:length(patient_ids)){
  
  # get mPDS scores for patient in on/off states
  patients_id <- patient_ids[i]
  plot_dt <- filter(dt, id == patients_id) %>% arrange(tsp)
  plot_dt <- plot_dt %>% mutate(color_pair = rep(c("A", "B"),nrow(plot_dt) / 2),
                                line_pair = rep(seq(1, nrow(plot_dt) / 2), each = 2))
  col.list <- rep(c("red","blue"), length(plot_dt$color_pair)/2)
  
  # mark even assessments as taken in off state and odd assessments as taken in on state
  plot_dt$color_pair[which(plot_dt$color_pair == "A")] = "Before Medication (OFF)"
  plot_dt$color_pair[which(plot_dt$color_pair == "B")] = "After Medication (ON)"
  plot_dt$color_pair <- as.factor(plot_dt$color_pair)
  time_diff <- plot_dt$tsp[length(plot_dt$tsp)] - plot_dt$tsp[1]
  
  # set off-state score to be plotted in red, on-state score in blue
  cbPalette <- c("blue", "red")
  # cast timestamps to consistent date format
  plot_dt$tsp <- as.Date(plot_dt$tsp)
  
  # generate mPDS scatterplot where lines connect off-state score to associated on-state score
  diffs <- ggplot(data = plot_dt) + 
    geom_line(aes_string(x = "tsp", y = j,
                         group = "line_pair")) +
    geom_point(aes_string(x = "tsp", y = j, color="color_pair"), size = 2) +
    xlab("Disease Timeline") + ylab("mPDS Score") +
    scale_color_discrete(labels = c("Before", "After"),
                         name = "Medication State")
  
  # get coefficients of best-fit line
  test_name <- paste('scores', substr(j, 8, nchar(j)), sep = "_")
  params <- c()
  if(j == "scores_react"){
    params <- coef(lm(scores_react ~ tsp, data = plot_dt))
  } else if(j == "scores_voice") {
    params <- coef(lm(scores_voice ~ tsp, data = plot_dt))
  } else if(j == "scores_balance") {
    params <- coef(lm(scores_balance ~ tsp, data = plot_dt))
  } else if(j == "scores_tap") {
    params <- coef(lm(scores_tap ~ tsp, data = plot_dt))
  } else {
    params <- coef(lm(scores_scaled ~ tsp, data = plot_dt))
  }
  pdss_inc <- append(pdss_inc, params[2])
  
  # add best-fit line to the plot
  diffs <- diffs + geom_abline(intercept = params[1], slope = params[2], col="darkgreen", lwd=1.25, linetype = "dashed") + 
    theme_bw() + theme(panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), 
                       axis.line = element_line(colour = "black")) +
    theme(axis.title = element_text(face="bold", size=8)) +
    theme(axis.text = element_text(size=8)) +
    ggtitle(paste("mPDS Severity Trajectory")) +
    theme(plot.title = element_text(color="black", face="bold", size=10)) + ylim(0,100) +
    scale_colour_manual(name="Status",values=cbPalette)
  
  # update plot legend
  diffs <- diffs + theme(legend.title = element_text(colour="black", size=6, face="bold")) +
    theme(legend.position=c(0.78,0.85),
          legend.direction="horizontal",
          legend.key.width=unit(1, "lines"), 
          legend.key.height=unit(1, "lines"),
          plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "lines"),
          legend.text=element_text(size=6))
  
  grid.arrange(diffs, nrow=1)
  plot_dt$color_pair <- as.character(plot_dt$color_pair)
  before <- plot_dt[which(plot_dt$color_pair == "Before Medication (OFF)"),]
  after <- plot_dt[which(plot_dt$color_pair == "After Medication (ON)"),]
  
  # can access off-state - on-state changes in mPDS with diff_scores
  diff_scores <- before$scores_scaled - after$scores_scaled
}
dev.off()