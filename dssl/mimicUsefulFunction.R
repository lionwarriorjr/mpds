#
# Library of Additional DSSL Functions (not pertaining to model training)
#
# If you use this code, please cite:
# Dyagilev, K. and Saria, S., 2016. Learning (predictive) risk scores in the
# presence of censoring due to interventions. Machine Learning, 102(3),
# pp.323-348.
# Dyagilev, K. and Saria, S., 2015. Learning a severity score for sepsis:
# A novel approach based on clinical comparisons. In AMIA Annual Symposium
# Proceedings (Vol. 2015, p. 1890). American Medical Informatics Association.


# Input argument:
# * (id, ht) - patient id and hospital time.
# * event - binary indicator if the current entry is the event of interest
# Return argument:
# * a list with fields timeToNextEvent, timeToFirstEvent, timeSinceLastEvent, timeSinceFirstEvent;
#
# For debugging
# id <- icustayID; ht <- hospTime; event <- cmiInd;
library(boot)
library(zoo)
dss.Auxiliary.CalculateTimeSinceAndToEvent <- function(id, ht, event){

  # get the number of entries
  numOfEntries = length(event);

  # initialization
  timeToNextEvent = rep(NA, numOfEntries);
  timeToFirstEvent = rep(NA, numOfEntries);
  timeSinceLastEvent = rep(NA, numOfEntries);
  timeSinceFirstEvent = rep(NA, numOfEntries);

  # parse forward to calculate time since event vectors
  lastId = NA; lastEventHT = NA; firstEventHT = NA;
  for(cnt in 1:numOfEntries){

    # get current values
    currId = id[cnt]; currHT = ht[cnt];

    # update if there was a change of user
    if(is.na(lastId) | currId != lastId){
      lastId = currId; lastEventHT = NA; firstEventHT = NA;
    }

    # if event was already encountered
    if(!is.na(lastEventHT) | !is.na(firstEventHT)){
      timeSinceLastEvent[cnt] = currHT - lastEventHT;
      timeSinceFirstEvent[cnt] = currHT - firstEventHT;
    }

    # if even is happening right now
    if(event[cnt] == 1){
      lastEventHT = currHT;
      if(is.na(firstEventHT)){
        firstEventHT = currHT;
      }
    }
  }

  # parse backward and calculate time to even vectors
  # cbind(id, ht, event, timeToNextEvent)[7700:8150,]
  lastId = NA; lastEventHT = NA;
  for(cnt in numOfEntries:1){
    #for(cnt in 8150:7700){

    # get current values
    currId = id[cnt]; currHT = ht[cnt];

    # update if there was a change of user
    if(is.na(lastId) | currId != lastId){
      lastId = currId; lastEventHT = NA;
    }

    # if event was already encountered
    if(!is.na(lastEventHT)){
      timeToNextEvent[cnt] = lastEventHT - currHT;
    }

    # if even is happening right now
    if(event[cnt] == 1){
      lastEventHT = currHT;
    }
  }
  cond = !is.na(timeToNextEvent) & is.na(timeSinceFirstEvent)
  timeToFirstEvent[cond] = timeToNextEvent[cond]

  res_list = list(timeToNextEvent = timeToNextEvent,
                  timeToFirstEvent = timeToFirstEvent,
                  timeSinceLastEvent = timeSinceLastEvent,
                  timeSinceFirstEvent = timeSinceFirstEvent)

  res_list
}


# counts the number of NA entries in the object
na.cnt <- function(obj){
  sum(is.na(obj))
}

# allows to treat correctly the case of vec containing a single value
mysample <- function(vec, num, replace = FALSE){
  if(length(vec) == 1){
    res_vec <- rep(vec, times = num)
  }else{
    res_vec <- sample(vec, num, replace = replace)
  }
  res_vec
}


bootstrapCIforAccuracy <- function(isCorrect){
  #x <- c(rep(1,1000), rep(0,1000))
  #b <- boot(x, function(u,i) mean(u[i]), R = 2000)
  b <- boot(isCorrect, function(u,i) mean(u[i]), R = 4000)
  ci <- boot.ci(b, type = c("basic"))
  ci
}


calculateCenteringAndRescalingDataForMatrix <- function(dataMat){

  # calculate the mean, the range, and the scale for each column
  mean = apply(dataMat, 2, FUN = mean)
  mx_feat = apply(dataMat, 2, max)
  mn_feat = apply(dataMat, 2, min)
  range_feat = mx_feat - mn_feat; range_feat[range_feat == 0] = 1;
  scale = 1/range_feat

  res_list <- list(mean = mean, scale = scale)
  res_list
}


centerAndRescaleMatrix <- function(dataMat, centeringCoeff){

  # get the sizes of the matrix
  numOfEntries = nrow(dataMat);
  numOfFeat = ncol(dataMat);

  # allocate auxiliary matrix and perform a basic arithmetic operations
  meanMat = matrix(data = centeringCoeff$mean, nrow = numOfEntries, ncol = numOfFeat, byrow = TRUE)
  scaleMat = matrix(data = centeringCoeff$scale, nrow = numOfEntries, ncol = numOfFeat, byrow = TRUE)
  dataMat_centered = (dataMat - meanMat)*scaleMat

  # return the result
  dataMat_centered
}


# swap entries of two vectors based on 0-1 indicator "ifswap"
vector_swap <- function(vec1, vec2, ifswap){
  vec1_swaped <- vec1*(1-ifswap) + vec2*ifswap
  vec2_swaped <- vec2*(1-ifswap) + vec1*ifswap
  res_list <- list(vec1 = vec1_swaped, vec2 = vec2_swaped)
  res_list
}


# This function might not generate the exact number of pairs
createRandomPairsOrderedByRank <- function(rankVec, numOfPairs){

  # get the number of entrie
  numOfEntries = length(rankVec)

  # potential on and off indices
  sampleOn = sample(numOfEntries, 2*numOfPairs, replace = T)
  sampleOff = sample(numOfEntries, 2*numOfPairs, replace = T)

  # remove rows with identical indices
  cond = (sampleOn == sampleOff) | (rankVec[sampleOn] == rankVec[sampleOff])
  sampleOn <- sampleOn[!cond]; sampleOff <- sampleOff[!cond];

  # swap the entries when needed
  ifswap = ifelse(rankVec[sampleOn] < rankVec[sampleOff], 1, 0)

  # swapping
  swapping_res <- vector_swap(sampleOn, sampleOff, ifswap)

  # get rid of repeating pairs
  orderingPairs <- data.frame(onIdx = swapping_res$vec1, offIdx = swapping_res$vec2)
  cond_duplicated = duplicated(orderingPairs)
  orderingPairs <- orderingPairs[!cond_duplicated,]

  # need to downsample into the required size
  currentNumberOfPairs = nrow(orderingPairs);
  if(currentNumberOfPairs>numOfPairs){
    orderingPairs <- orderingPairs[1:numOfPairs,];
  }

  orderingPairs
}


# Katie's function for calculating nice color palette
gg_color_hue<- function(n) {
  hues=seq(15,375,length=n+1)
  hcl(h=hues,l=65,c=100)[1:n]
}


# translate into binary vectors
aux.DecimalVecToBinMatrix <- function(decVec, numOfClasses){

  # checking the input
  uniqClasses = sort(unique(decVec));
  if(sum(uniqClasses != 0:(numOfClasses-1))){
    cat("aux.DecimalVecToBinMatrix: The list of unique classes is not concordant with", numOfClasses, "classes\n")
    stop()
  }

  # number of bits
  numOfBits <- ceil(log2(numOfClasses))

  # allocate matrix
  resMat <- matrix(NA, ncol = numOfBits, nrow = length(decVec))

  tmpVec = decVec
  for(cntB in numOfBits:1){
    currBit = mod(tmpVec, 2);
    resMat[,cntB] = currBit;
    tmpVec <- (tmpVec - currBit)/2
  }


  resMatColNames = c();
  for(cntB in 1:numOfBits){
    resMatColNames[cntB] <- sprintf("ClassBit%d", numOfBits - cntB + 1)
  }
  colnames(resMat) <- resMatColNames

  # return the value
  resMat
}


energyOfVector.total <- function(vec){
  sum(abs(vec))
}

energyOfVector.analyzefDifferenceOfTwoVectors <- function(vec1, vec2){
  diff_vec <- vec1 - vec2;
  totalDiffEnergy = energyOfVector.total(diff_vec)
  fractionOfDiffEnergyPerComponent <- 2*abs(diff_vec)/(abs(vec1) + abs(vec2))
  fractionOfTotalDiffEnergy <- 2*totalDiffEnergy/(energyOfVector.total(vec1) + energyOfVector.total(vec2))
  res_list = list(totalDiffEnergy = totalDiffEnergy,
                  fractionOfDiffEnergyPerComponent = fractionOfDiffEnergyPerComponent,
                  fractionOfTotalDiffEnergy = fractionOfTotalDiffEnergy)

  res_list
}


df.bindOrUpdate <- function(df, vec, vec_name){
  if(vec_name %in% names(df)){
    df[,vec_name]  <- vec
  }else{
    df <- cbind(df, vec);
    df_names <- names(df);
    df_names[length(df_names)] <- vec_name;
    names(df) <- df_names;
  }
  df
}


union.list <- function(set_list){
  curr_res <- set_list[[1]];
  if(length(set_list) > 1){
    for(cnt in 2:length(set_list)){
      curr_res <- union(curr_res, set_list[[cnt]])
    }
  }
  curr_res
}


union.all <- function(...){
  input_args <- list(...)
  union.list(input_args)
}


# this function proves very useful in various contexts
vector.getFirstValue <- function(vec){
  vec[1]
}


# Function mimic.SampleAndHoldPropagation.WithExpiration(id, ht, vec, valid_time)
# Input arguments
# * (id, ht) - icustay_id and hospitalization time. We assume that (id, ht) is ordered by id and then by ht.
# * vec - values to be interpolated, is either NA or numeric value
# * valid_time - how long does the measurement remains valid (in minutes)
#
# Returned value is a list with three fields:
# * locf - last one carry forward as it is
# * timeSinceLastSample - time since the last value was observed, NA if there was no value yet
# * locf_withExpir - locf with expiration
mimic.SampleAndHoldPropagation.WithExpiration <- function(id, ht, vec, valid_time, markerForFirstEntraceOfID = "FEOfID"){

  # create a working copy of the vector
  copy_vec = vec
  last_sample = ht;
  last_sample[is.na(copy_vec)] = NA

  # find indexes of all first entries of ids that are equal NA
  # replace them with temporary place holder
  first_idEntry_idx = which(diff(id) != 0)  + 1;
  is_firstIDEntry_na = is.na(copy_vec[first_idEntry_idx])
  copy_vec[first_idEntry_idx[is_firstIDEntry_na]] = markerForFirstEntraceOfID
  last_sample[first_idEntry_idx[is_firstIDEntry_na]] = markerForFirstEntraceOfID

  # now replace all NA with last one carry forward, except for the first values that are NA
  copy_vec = na.locf(copy_vec, na.rm = FALSE)
  copy_vec[copy_vec == markerForFirstEntraceOfID] = NA
  copy_vec <- as.numeric(copy_vec)
  last_sample = na.locf(last_sample, na.rm = FALSE)
  last_sample[last_sample == markerForFirstEntraceOfID] = NA;
  last_sample <- as.numeric(last_sample)

  # find out how much time passed since last sample
  time_since_last_sample = ht - last_sample;
  vec_with_expiration = copy_vec;
  cond_expiration = !is.na(time_since_last_sample) & (time_since_last_sample > valid_time)
  vec_with_expiration[cond_expiration] = NA

  res_list = list(locf = copy_vec, timeSinceLastSample = time_since_last_sample,
                  locf_withExpir = vec_with_expiration)
  res_list
}


# Calculate instantaneuous conditions for SIRS.
# Function mimic.CalculateInstantaneuousSIRSConditions
# Input arguments:
# * temperature, hr, resp_rate, paco2, wbc - numerical or NA
# Output argument: a list with the following fields, all TRUE/FALSE/NA
# * temperature - standard temperature condition
# * hr - standard hr condition
# * resp_rate - standard resp_rate condition
# * paco2 - standard paco2 condition
# * wbc - standard wbc condition
# * resp_rate_or_paco2 - if one of the conditions hold (the other one can be NA)
# * sirs_raw - instantaneuous SIRS signal, i.e., whether two out of four conditions hold (the rest can be NA)
# * none_holds - all data is non-NA but none of the conditions hold
mimic.CalculateInstantaneuousSIRSConditions <- function(temperature, hr, resp_rate, paco2, wbc){

  # calculate each condition separately (makes use of the fact that NA|TRUE == TRUE )
  temp_cond = temperature < 96.8  | temperature > 100.4;
  hr_cond = hr > 90
  rr_cond = resp_rate > 20;
  paco2_cond = paco2 < 32;
  wbc_cond = wbc < 4 | wbc > 12;
  rr_or_paco2_cond = rr_cond | paco2_cond

  # see if any of two of them hold simultaneously
  # i.e., this is where we have enough information to pronounce two of the conditions as currently holding
  tmp_cond = cbind(temp_cond, hr_cond, rr_or_paco2_cond, wbc_cond)
  #tmp_cond[is.na(temp_cond)] = 0;
  sirs_raw_cond = (rowSums(tmp_cond, na.rm = TRUE) >=2)

  # find entries where we have all the information but none of the conditions holds
  none_holds_cond = (temp_cond == FALSE) & (hr_cond == FALSE) & (rr_cond == FALSE) &
    (paco2_cond == FALSE | is.na(paco2_cond)) & (wbc_cond == FALSE | is.na(wbc_cond))

  # resulting list
  res_list = list(temperature = temp_cond, hr = hr_cond, resp_rate = rr_cond,
                  paco2 = paco2_cond, wbc = wbc_cond, resp_rate_or_paco2 = rr_or_paco2_cond,
                  sirs_raw = sirs_raw_cond, none_holds = none_holds_cond)

  res_list
}


# Function mimic.RemoveShortOnIntervals
# Input arguments:
# * id, ht - icustay_id, hosp_time. The data is assumed to be ordered by (id,ht)
# * vec - vector of 0/1/NA values to be processed
# * timeTh - required time length of an undisturbed window
# Output arguments:
# * a vector of processed value
mimic.RemoveShortOnIntervals <- function(id, ht, vec, timeTh){

  # calculate a vector indicating whether the current entry is the last
  # one for the current user
  is_last_user_entry = ifelse(c(diff(id) != 0, TRUE),1,0)
  sig_raw = vec; sig_raw[is.na(sig_raw)] = 0;

  numOfEntries <- length(ht)
  out_sig <- rep(0, numOfEntries)
  prevVal = -1; prevT = -timeTh - 1 + min(timeTh);
  risingEdgeT = -timeTh - 1; risingEdgeIdx = -1;
  for(cntR in 1:numOfEntries){
    currT = ht[cntR];
    currVal = sig_raw[cntR];
    is_last = is_last_user_entry[cntR];

    # rising edge
    if(currVal == 1 & prevVal <= 0){
      risingEdgeT = ht[cntR]; risingEdgeIdx = cntR;
    }

    # falling edge
    if(currVal == 0 & prevVal == 1){
      # if there were at least 5 hours
      # then declare the whole interval as sirs
      if(prevT - risingEdgeT >= timeTh){
        for(cnt in risingEdgeIdx:(cntR-1)){
          out_sig[cnt] = 1;
        }
      }
    }

    # this is the last one for the user
    if(currVal == 1 & prevVal == 1 & is_last == 1){
      if(currT - risingEdgeT >= timeTh){
        for(cnt in risingEdgeIdx:(cntR)){
          out_sig[cnt] = 1;
        }
      }
    }

    # make sure to update all the runners
    if(is_last == 1){
      prevVal = -1; prevT = -timeTh - 1;
    }else{
      prevVal= currVal;
      prevT = currT;
    }

  }

  out_sig
}


# Function mimic.UniteAdjacentOnIntervals unites "on" intervals
# that are separated by less than unionTh minutes
#
# Input arguments:
# * id, ht - icustay_id, hosp_time. The data is assumed to be ordered by (id,ht)
# * vec - vector of 0/1/NA values to be processed
# * unionTh - maximal separation for intervals to be uniteds
#
# Output arguments:
# * a vector of processed value
mimic.UniteAdjacentOnIntervals <- function(id, ht, vec, unionTh){

  # compatibility transformations
  is_last_user_entry = ifelse(c(diff(id) != 0, TRUE),1,0)
  sig_raw = vec; sig_raw[is.na(sig_raw)] = 0;
  ht[ht<0] = 0

  # initialization
  numOfEntries <- length(ht)
  prevVal = -1; prevT = -1;
  fallingEdgeT = -unionTh - 1; fallingEdgeIdx = -1;

  # run over all rows
  for(cntR in 1:numOfEntries){

    # get current time, value and whether this is the last entry for a user
    currT = ht[cntR];
    currVal = sig_raw[cntR];
    is_last = is_last_user_entry[cntR];

    # rising edge
    if(prevVal == 0 & currVal == 1){
      if(currT - fallingEdgeT <= unionTh){
        for(cnt in fallingEdgeIdx:cntR){
          sig_raw[cnt] = 1;
        }
      }
    }

    # falling edge
    if(prevVal == 1 & currVal == 0){
      fallingEdgeT = prevT; fallingEdgeIdx = cntR - 1;
    }

    # update running variables
    if(is_last){
      prevVal = -1; prevT = -1;
      fallingEdgeT = -unionTh - 1; fallingEdgeIdx = -1;
    }else{
      prevVal = currVal; prevT = currT;
    }

  }

  sig_raw
}


mimic.TranslateLogicalVectorToZeroOne <- function(vec){
  res_vec <- vec;
  res_vec[res_vec == TRUE] = 1; res_vec[res_vec == FALSE] = 0;
  res_vec
}


mimic.CoarseGradingSurroundedBySameStatus <-
  function(id, ht, trusted_status_raw, statusOfInterest, timeTh){

    status_noyes = rep(1, length(trusted_status_raw));
    status_noyes[trusted_status_raw == statusOfInterest] = 0;
    timeToData <- dss.Auxiliary.CalculateTimeSinceAndToEvent(id, ht, status_noyes);

    status_safe = trusted_status_raw == statusOfInterest & (is.na(timeToData$timeToNextEvent) | timeToData$timeToNextEvent > timeTh) &
      (is.na(timeToData$timeSinceLastEvent) | timeToData$timeSinceLastEvent > timeTh)

    res_list <- list(timeToData = timeToData, status_safe = status_safe)
    res_list
}


mimic.getSampleAndHoldSepsisStageMarkers <- function(dataMat){

  # extract all the data as usual
  severe_sepsis_sh = dataMat$severe_sepsis;
  septic_shock_sh = dataMat$septic_shock;
  sirs_inpt_sh = dataMat$sirs_intp;

  # status signals
  none_sh = (sirs_inpt_sh == 0) & (severe_sepsis_sh == 0) & (septic_shock_sh == 0);
  h_status_sh <- mimic.highest_status(sirs_inpt_sh, severe_sepsis_sh, septic_shock_sh);
  status_sh <- factor(h_status_sh, labels = c("none", "sirs", "severe", "shock"), ordered = TRUE)

  # sirs related signals
  sirs_raw_sh <- dataMat$sirs_raw
  sirs_temp_sh <- dataMat$sirs_temperature_oor
  sirs_hr_oor_sh <- dataMat$sirs_hr_oor
  sirs_resp_oor_sh <- dataMat$sirs_resp_oor
  sirs_wbc_oor_sh <- dataMat$sirs_wbc_oor

  # indicators used definition of infection suspicion
  angus_infection <- dataMat$angus_infection
  sepsis_note <- dataMat$sepsis_note
  time_to_ss = dataMat$minutes_to_shock_onset;

  # push all the signals into a single data frame
  res <- data.frame(severe_sepsis_sh = severe_sepsis_sh, septic_shock_sh = septic_shock_sh, sirs_inpt_sh = sirs_inpt_sh,
                    none_sh = none_sh,  status_sh = status_sh, sirs_raw_sh = sirs_raw_sh, sirs_temp_sh = sirs_temp_sh,
                    sirs_hr_oor_sh = sirs_hr_oor_sh, sirs_resp_oor_sh = sirs_resp_oor_sh, sirs_wbc_oor_sh = sirs_wbc_oor_sh,
                    sepsis_note = sepsis_note, angus_infection = angus_infection, time_to_ss = time_to_ss);
  res
}


# function written by Peter Schulam
mimic.highest_status <- function(sirs, severe, shock) {
  f <- ifelse(
    shock == 1,
    "shock",
    ifelse(
      severe == 1,
      "severe",
      ifelse(
        sirs == 1,
        "sirs",
        "none"
      )
    )
  )

  ordered(f, c("none", "sirs", "severe", "shock"))
}


mimic.dropDFColumnsByListOfNames <- function (df, killList){
  df_new <- df[,!(names(df) %in% killList)];
  df_new
}


# this function calculates defaults for different features
# We assume that the input matrix is unpropagated, thus
# the population median is a good starting point.
# For some of the features, the default value is manually corrected
# to their "normal/healthy" value.
mimic.calculateFeatureDefaults <- function(dataMat){

  # get population level median
  # can be solved more elegantly
  medianNA <- function(x){median(x, na.rm = TRUE)}
  populationMedian = apply(dataMat, 2, FUN = medianNA)

  # some of the values need to be defaulted to their "normal" state
  # rather than to some population median. E.g., aids when not mentioned is usally 0.
  correctionsToPopMed = list();
  correctionsToPopMed["aids"] <- 0;  # no aids by default
  correctionsToPopMed["met_carcinoma"] <- 0; # no carcinoma by default
  correctionsToPopMed["organ_insuff"] <- 0; # no organ_insufficiency by default
  correctionsToPopMed["hem_malig"] <- 0;
  correctionsToPopMed["immuno_comp"] <- 0;
  correctionsToPopMed["neurologic_sofa"] <- 0;
  correctionsToPopMed["cardio_sofa"] <- 0;
  correctionsToPopMed["resp_sofa"] <- 0;
  correctionsToPopMed["hepatic_sofa"] <- 0;
  correctionsToPopMed["vent"] <- 0;
  correctionsToPopMed["pacemkr"] <- 0;
  correctionsToPopMed["all_input_24hr"] <- 0;
  correctionsToPopMed["gcs"] <- 15;
  correctionsToPopMed["resp_rate"] <- 18;
  correctionsToPopMed["any_input"] <- 0;
  correctionsToPopMed["bilirubin"] <- 1.9;
  correctionsToPopMed["worst_sofa"] <- 0;

  corrFeatNames <- names(correctionsToPopMed);
  for(cnt in seq(along = correctionsToPopMed)){
    featName <- corrFeatNames[[cnt]];
    correctVal <- correctionsToPopMed[[cnt]];
    
    # if the feature exists in the current population median table
    if(featName %in% names(populationMedian)){
      populationMedian[[featName]] <- correctVal;
    }

  }

  if(sum(is.na(populationMedian))==1){
    stop("One or more of the values in populationMedian are NA")
  }

  populationMedian
}


mimic.aux.linear_approx <- function(x,y){
  
  # get the number of elements
  lngth = length(x);

  # get the indices of all non NA elements
  nonNA_idx = which(!is.na(y))
  numNonNA = length(nonNA_idx);

  # sanity checks
  if(numNonNA <= 1){
    error("linear_approx: less than one non-Na point in the vector")
  }
  if(is.na(y[1]) | is.na(y[lngth])){
    error("linear_approx: either the first element or the last one are NA")
  }

  # prepare current and next value
  prevIdx = nonNA_idx[1:(numNonNA - 1)];
  nextIdx = nonNA_idx[2:numNonNA];
  xPrev = x[prevIdx]; yPrev = y[prevIdx]
  xNext = x[nextIdx]; yNext = y[nextIdx]
  intCnt = 0;
  for(cnt in seq(along = x)){
    if(is.na(y[cnt])){
      yn = yNext[intCnt]; yp = yPrev[intCnt];
      xn = xNext[intCnt]; xp = xPrev[intCnt];
      xc = x[cnt];
      yc = (yn*(xc-xp) + yp*(xn-xc))/(xn-xp);
      y[cnt] = yc
    }else{
      intCnt = intCnt + 1;
    }
  }
  res <- y
  res
}


# function written by Kirill Dyagilev
mimic.singlePatientLinearInterpolation <- function(x, y, defVal){

  # get the number of elements
  lngth = length(x);

  # get the indices of all NA elements
  nonNA_idx = which(!is.na(y))

  # return the vector as it is
  if(lngth == length(nonNA_idx)){
    return(y)
  }

  # if all NA then exit, return a list of defaults
  if(length(nonNA_idx) == 0){
    res <- rep(defVal, lngth);
    return(res);
  }

  # make sure the last element is non NA
  if(is.na(tail(y, n=1))){
    last_nonNA_val = y[nonNA_idx[length(nonNA_idx)]]
    y[lngth] = last_nonNA_val
  }

  # make sure the first element is non NA
  if(is.na(head(y, n=1))){
    #first_nonNA_val = y[nonNA_idx[1]]
    y[1] = defVal
  }

  # use the approximation function
  cond = !is.na(y);
  na_idx = which(!cond);
  not_na_idx = which(cond);
  y_padded <- mimic.aux.linear_approx(x, y);
  y_padded
}


mimic.singlePatientSampleAndHoldInterpolation <- function(vec, defVal){

  if(is.na(vec[1])){
    vec[1] = defVal;
  }
  vec_padded <- na.locf(vec)
  vec_padded
}


mimic.sampleUniquePairs.withReplacement <- function(vec1, vec2, numOfPairs){

  tmp1 <- sample(vec1, 4*numOfPairs, replace = T)
  tmp2 <- sample(vec2, 4*numOfPairs, replace = T)


  tmp_df <- data.frame(vec1 = tmp1, vec2 = tmp2)
  tmp_df <- unique(tmp_df)

  if(nrow(tmp_df) > numOfPairs){
    tmp_df <- tmp_df[sample(nrow(tmp_df), numOfPairs) ,]
  }
  tmp_df
}
