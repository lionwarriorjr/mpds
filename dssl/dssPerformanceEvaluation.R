################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################
# Move later to a separate file
#
################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################
################################################################################################################################



# combine data on per-patient basis
# ids = test_info[,1]; time_to_e = test_info[,2]; cens = test_info[,3];
# res_df <- perUserIndCalc(ids, time_to_e, cens)
perUserIndicatorCalculation <- function(ids, time_to_e, cens){

  # checking input
  if(length(ids) != length(time_to_e) | length(time_to_e) != length(cens)){
    stop("perUserIndCalc: input vectors are of different length")
  }

  # max - will give NA iff all the values are NA
  tmp_fn <- function(v){v[1]}
  isPos = !is.na(tapply(time_to_e, ids, FUN = tmp_fn))
  isCens = !is.na(tapply(cens, ids, FUN = tmp_fn))

  uniq_ids = as.numeric(names(isCens))

  # for each id, create a label -1, 0, 1
  # -1 for negative censored,
  # 0 for negative uncensored
  # +1 for positive
  lbl = ifelse(isPos, 1, 0);
  lbl[ isCens & !isPos] = -1;


  tmp_df <- data.frame(id = uniq_ids, isPos = isPos, isCens = isCens, lbl = lbl)
  tmp_df
}

#
# This function finds the maximal obtained score value for each user
# Input arguments:
# - two vectors (ids, scores) with scores for different patients at different time steps.
# Returned value:
# data frame with two columns, uniqId - the list of all unique ids in the data,
# maxScore - a vector of maximal scores for each unique user
#
calcPerUserMaxScore <- function(ids, scores){

  # find the maximal obtained score value for each user
  maxScores <- tapply(scores, ids, FUN = max);
  uniq_usr <- as.numeric(names(maxScores))

  #
  tmp_df  = data.frame(uniqId = uniq_usr, maxScore = maxScores)
  tmp_df
}


#
# Calculate time-to-detection for each user in a given set
# Input arguments:
# * (ids, scores, time_to_event, time_to_organ_failure) - icustayId, the value of the score, and the corresponding time-to-event,
#   the value of time_to_event is always strictly positive prior to the first event and NA for negative patients
#   and after the first event. time_to_organ_failure - positive value indicates minutes to first organ failure (should less or equal to time to
#   event), it is negative after the first organ failure, and NA for patients without any organ failure.
# * positiveIds - list of ids of all relevant positive patients
# * (truePositiveTh, falseNegativeTh, detectionTh) - three lists that defines measurement regime, i.e., the value of
#   score detection threshold. Our algorithm assumes that the detection point is the earliest entrie of the patient
#   in which the value of the score exceeds the detection threshold. The value of the detection threshold can be
#   given as a value (detectionTh), or as desired false negative (falseNegativeRate) or true positive (truePositiveRate) rates.
# Return value:
# *
#
#
#
# For debugging:
# * ids = icustayID; scores = dssRanks; time_to_event = minToShock; positiveIds =  userInfo_test$id[userInfo_test$isPos == 1];
#   negativeIds = userInfo_test$id[userInfo_test$isPos == 0 & userInfo_test$isCens == 0];
#   truePositiveTh = c(0.85); falsePositiveTh = c(0.35); detectionTh = c(-3.03)
calcTimeToDetection <- function(ids, scores, time_to_event, time_to_organ_failure, positiveIds, negativeIds,
                                truePositiveTh = c(), falsePositiveTh = c(), detectionTh = c()){

  # calculate the maximal score for all users
  maxScores_df <- calcPerUserMaxScore(ids = ids, scores = scores)


  # get relevant entries
  rel_ids = c(positiveIds, negativeIds)
  rel_lbl = c(rep(1, length(positiveIds)), rep(0, length(negativeIds)))
  map_idx = match(rel_ids, maxScores_df$uniqId);
  rel_max_scores = maxScores_df$maxScore[map_idx]
  roc_obj <- roc(rel_lbl, rel_max_scores)
  ci_obj = ci(roc_obj)

  # augmenting roc object
  num_pos = length(roc_obj$cases)
  num_neg = length(roc_obj$controls)
  tpr = roc_obj$sensitivities;
  fpr = 1 - roc_obj$specificities;
  num_found = num_pos*tpr + num_neg*fpr

  recall = tpr;
  precision = num_pos*tpr/num_found
  f_measure = (2*precision*recall)/(precision + recall)
  roc_obj$precision = precision
  roc_obj$recall = recall
  roc_obj$f_measure = f_measure

  neg_id_ids <- maxScores_df$uniqId[maxScores_df$uniqId %in%  negativeIds]
  neg_id_max_scores <- maxScores_df$maxScore[maxScores_df$uniqId %in%  negativeIds]
  pos_id_ids <- maxScores_df$uniqId[maxScores_df$uniqId %in%  positiveIds]
  pos_id_max_scores <- maxScores_df$maxScore[maxScores_df$uniqId %in%  positiveIds]

  #
  # calulate the list of all thresholds and their names
  #
  th_list = c()
  th_names = c()
  th_type = c()
  cnt = 0;

  # calculate thresholds on the score that achieve the desired truePositiveRate
  for(cnt_tp_rt in seq_along(truePositiveTh)){
    tp_rt = truePositiveTh[cnt_tp_rt]
    cnt = cnt + 1
    th_list[cnt] = quantile(pos_id_max_scores, c(1-tp_rt))
    th_names[cnt] = sprintf("TruePosTh_%d", cnt_tp_rt)
    th_type[cnt] = "TruePositive"
  }

  # calculate thresholds on the score that achieve the desired falsePositiveRate
  for(cnt_fp_rt in seq_along(falsePositiveTh)){
    fp_rt = falsePositiveTh[cnt_fp_rt]
    cnt = cnt + 1
    th_list[cnt] = quantile(neg_id_max_scores, c(fp_rt))
    th_names[cnt] = sprintf("FalsePosTh_%d", cnt_fp_rt)
    th_type[cnt] = "FalsePositive"
  }

  # copy the manually defined thresholds
  for(cnt_det_th in seq_along(detectionTh)){
    det_th = detectionTh[cnt_det_th]
    cnt = cnt + 1
    th_list[cnt] = det_th
    th_names[cnt] = sprintf("DetTh_%d", cnt_det_th)
    th_type[cnt] = "DetectionThreshold"
  }


  # calculate detection time for all patients
  time_to_event[is.na(time_to_event)] = 0;
  tmp_df = data.frame(ids = ids, t2e = time_to_event, score = as.vector(scores), orgFail = time_to_organ_failure);
  cond = tmp_df$ids %in% positiveIds;
  tmp_df <- tmp_df[cond, ]

  detTime_df_list = list()
  for(cnt_th in seq_along(th_list)){
    # current detection threshold
    th = th_list[cnt_th]
    # 0-1 indicator whether the current score is above the detection threshold
    mult_coeff = ifelse(tmp_df$score >=  th,1,0)

    # make a temporary copy of the data frame
    curr_df <- tmp_df

    # calculate the time to event (recall that these are positive patients)
    # it is calculated as the per-user maximal time_to_event entrie for which the score is above detection threshold
    # it equals zero if the septic shock was not detected prior to its onset
    curr_df$t2e <- curr_df$t2e*mult_coeff
    time_between_detection_and_septic_shock <- tapply(curr_df$t2e, curr_df$ids, FUN = max)
    cond_detPriorSS = time_between_detection_and_septic_shock >0

    # we repeat a similar procedure for the time to organ failure
    # first_org_fail - positive for detection prior to organ failure, otherwise negative.
    # cannot be NA since we consider only positive patients.
    time2OrganFail = curr_df$orgFail*mult_coeff;
    first_org_fail <- tapply(time2OrganFail, curr_df$ids, FUN = max)
    # cond_detPriorToOF <- first_org_fail>= 0 & cond_detPriorSS

    detTime_df_list[[cnt_th]] = data.frame(uniqId = as.numeric(names(time_between_detection_and_septic_shock)[cond_detPriorSS]),
                                           detectionTime = time_between_detection_and_septic_shock[cond_detPriorSS],
                                           detToOrgFailTime = first_org_fail[cond_detPriorSS])
  }

  #
  # calculate some summary statistics for detection times
  #
  tmp_val = c(truePositiveTh, falsePositiveTh, detectionTh);
  fpRate = c(); tpRate = c(); detTime_median = c();
  detTime_IQR_l = c(); detTime_IQR_u = c();
  fractionOfDetectedPriorToOrgFail = c();
  detTimePriorToOrgFail_median = c();
  detTimePriorToOrgFail_IQR_l = c();
  detTimePriorToOrgFail_IQR_u = c();
  for(cnt_th in seq_along(th_list)){
    # current detection threshold
    th = th_list[cnt_th]

    # actual true and false positive rates for the current detection threshold
    tpRate[cnt_th] = mean(pos_id_max_scores >= th)
    fpRate[cnt_th] = mean(neg_id_max_scores >= th)

    # calculate summary statistics for detection time
    detTime_median[cnt_th] = median(detTime_df_list[[cnt_th]]$detectionTime);
    detTime_IQR_l[cnt_th] = quantile(detTime_df_list[[cnt_th]]$detectionTime,1/4);
    detTime_IQR_u[cnt_th] = quantile(detTime_df_list[[cnt_th]]$detectionTime,3/4);

    # calculate summary statistics for time between the detection and the first organ failure
    fractionOfDetectedPriorToOrgFail[cnt_th] = mean(detTime_df_list[[cnt_th]]$detToOrgFailTime>0)
    tmp_vec = detTime_df_list[[cnt_th]]$detToOrgFailTime; tmp_vec <- tmp_vec[tmp_vec>0];
    detTimePriorToOrgFail_median[cnt_th] = median(tmp_vec)
    detTimePriorToOrgFail_IQR_l[cnt_th] = quantile(tmp_vec,1/4)
    detTimePriorToOrgFail_IQR_u[cnt_th] = quantile(tmp_vec,3/4)
  }


  sumStats <- data.frame(type = th_type, value = tmp_val, th = th_list, falsePositiveRate = fpRate, truePositiveRate = tpRate,
                         detTime_median = detTime_median, detTime_IQR_l = detTime_IQR_l, detTime_IQR_u = detTime_IQR_u,
                         fractionOfDetectedPriorToOrgFail = fractionOfDetectedPriorToOrgFail,
                         detTimePriorToOrgFail_median = detTimePriorToOrgFail_median,
                         detTimePriorToOrgFail_IQR_l = detTimePriorToOrgFail_IQR_l, detTimePriorToOrgFail_IQR_u = detTimePriorToOrgFail_IQR_u)

  res <- list(detailedTime_df = detTime_df_list, summaryStatistics = sumStats, maxScores = maxScores_df,
              medianAUC = ci_obj[2], AUC_ci95_l = ci_obj[1],  AUC_ci95_u = ci_obj[3], roc_obj = roc_obj)
  res
}
