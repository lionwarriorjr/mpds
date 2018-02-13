#
# this is the library of DSS functions
#
library(pracma)
library(glmnet)
library(dplyr)
library(pROC)
library(permute)
library(pracma)
library(ggplot2)
library(grid)
library(rpart)
source("mimicUsefulFunction.R")

# for (package in c('<package1>', '<package2>')) {
#   if (!(require(package, character.only=T, quietly=T))) {
#     install.packages(package)
#     library(package, character.only=T)
#   }
# }


# A = ii1; B = ii2; N = numOfPoints;
subsampleSetOfPossiblePairs<- function(A,B, N){
  if(as.double(as.numeric(length(A))*as.numeric(length(B)))<=(10*N)){
    tmp <- expand.grid(A,B);
  }else{
    tmp <- data.frame(
      Var1 = mysample(A, 10*N, replace = TRUE),
      Var2 = mysample(B, 10*N, replace = TRUE)
    )
    tmp <- unique(tmp)
  }
  sampleNum = min(c(nrow(tmp), N))
  tmp <- tmp[sample(1:nrow(tmp), sampleNum),]
  res_df = data.frame(onIdx = tmp$Var1, offIdx = tmp$Var2)
  res_df
}


#tmp_pairs_df <- generateBalancedPairs(patD, numberOfIntraPairsPerPatientPerPair);
# inp_df = patD; numOfPoints = numberOfIntraPairsPerPatientPerPair; min_ht = 2*60

generateBalancedPairs <- function(inp_df, numOfPoints, min_ht = 2*60){

  # Extract the vectors
  orig_idx = inp_df$idx;
  stage = inp_df$stage;
  ht <- inp_df$hosp_time;

  resulting_idx = list()

  # create all cross-stage combinations
  uniq_stages = unique(stage);
  uniq_stages = uniq_stages[uniq_stages>=0];
  stage_pairs = expand.grid(uniq_stages,uniq_stages)
  cond = (stage_pairs$Var1 > stage_pairs$Var2);
  stage_pairs = stage_pairs[cond,]
  names(stage_pairs) <- c("high", "low")
  numberOfPairs = sum(cond)




  # generate pairs for each combination of stages
  for (cntP in 1:numberOfPairs){
    # get all the indexes
    ii1 = which(stage == stage_pairs$high[cntP]);
    ii2 = which(stage == stage_pairs$low[cntP]);

    if(length(ii2)>0 & length(ii1)>0){
      tmp_df <- subsampleSetOfPossiblePairs(ii1, ii2,numOfPoints);
      cond = abs(ht[tmp_df$onIdx] - ht[tmp_df$offIdx])>=min_ht;
      tmp_df = tmp_df[cond,]

      #
      resulting_idx[[length(resulting_idx)+1]] = tmp_df;
    }
  }

  # resulting pairs
  tmp_pairs = rbind_all(resulting_idx);

  # need to make sure that there are no repetitions
  if(nrow(tmp_pairs)>1){
    tmp_pairs = unique(tmp_pairs)
  }


  res_pairs <- data.frame(onIdx = orig_idx[tmp_pairs$onIdx],
                          offIdx = orig_idx[tmp_pairs$offIdx])

  res_pairs
}


#dssGenerateClinicalComparisons(icustayID, severity, hospTime, ids_list[[cnt]],
#                               numberOfIntraPairsPerPatientPerPair)
# patientIDVec = icustayID; levelsVec = severity;
# timeVec = hospTime;  consideredIDs = ids_list[[1]];
# numberOfIntraPairsPerPatientPerPair = numberOfIntraPairsPerPatientPerPair; numberOfInterPairsPerPair = NA
#  dssGenerateClinicalComparisons(conc_id_sorted, conc_sv_sorted_new, conc_ht_sorted, ids_list[[cnt]],
#                                            numberOfIntraPairsPerPatientPerPair)
# patientIDVec = conc_id_sorted; levelsVec = conc_sv_sorted_new; timeVec = conc_ht_sorted;
# consideredIDs = ids_list[[1]]; numberOfIntraPairsPerPatientPerPair = numberOfIntraPairsPerPatientPerPair;
# numberOfInterPairsPerPair = NA

# Generates both Intra- and Inter-Patient Clinical Comparisosns
dssGenerateClinicalComparisons <- function(
  patientIDVec, levelsVec,
  timeVec, consideredIDs,
  numberOfIntraPairsPerPatientPerPair,
  numberOfInterPairsPerPair = NA
  )
{

  # this function generates clincal comparisons
  numOfEntries = length(patientIDVec);
  inp_data_df = data.frame(idx = 1:numOfEntries, stage = levelsVec,
                           hosp_time = timeVec)

  # leave only relevant input data
  cond = patientIDVec %in% consideredIDs
  inp_data_df = inp_data_df[cond,]
  rel_ids = patientIDVec[cond]


  # get all intra-patient clinical comparisons
  all_pairs_list = list();
  inp_df_rgd <- split(inp_data_df, as.factor(rel_ids));
  for(patD in inp_df_rgd){
    tmp_pairs_df <- generateBalancedPairs(patD,
                                          numberOfIntraPairsPerPatientPerPair
                                          );

    # if weird pairs are found
    if(sum(levelsVec[tmp_pairs_df$onIdx] <= levelsVec[tmp_pairs_df$offIdx])){
      cat("Problematic id:", conc_id_sorted[patD$idx[1]], "\n")
    }

    all_pairs_list[[length(all_pairs_list)+1]] = tmp_pairs_df
  }
  all_pairs_df = rbind_all(all_pairs_list)


  # now get the inter-patient clinical comparisons
  if(is.na(numberOfInterPairsPerPair)){
    numOfStages = length(unique(levelsVec));
    numOfStagePairs = numOfStages*(numOfStages-1)/2
    numberOfInterPairsPerPair = round(nrow(all_pairs_df)/numOfStagePairs)
  }

  # generate pairs only between entries that were already used
  # it will save a lot of training time for GBRTs
  uniq_idx = unique(c(all_pairs_df$onIdx, all_pairs_df$offIdx))
  cond = inp_data_df$idx %in% uniq_idx
  inp_data_df = inp_data_df[cond, ]
  tmp_pairs_df <- generateBalancedPairs(inp_data_df, numberOfInterPairsPerPair);
  all_pairs_df <- unique(rbind(all_pairs_df, tmp_pairs_df))

  all_pairs_df

}



dssGenerateSmoothnessPairs <- function(
  patientIDVec, clinicalComparisons,
  numberOfPairs = NA)
{
  # generate everything, then subsample is needed
  all_indexes = sort(unique(c(clinicalComparisons$offIdx, clinicalComparisons$onIdx)));
  numOfSmoothPairs = length(all_indexes);
  smooth_pairs_df = data.frame(offIdx = all_indexes, onIdx = all_indexes+1);

  # check if the last entry is out of range
  badLastEntry = 0;
  numOfEntries = length(patientIDVec)
  if(smooth_pairs_df$onIdx[numOfSmoothPairs] > numOfEntries){
    smooth_pairs_df$onIdx[numOfSmoothPairs] = smooth_pairs_df$offIdx[numOfSmoothPairs];
    badLastEntry = 1;
  }

  cond = patientIDVec[smooth_pairs_df$offIdx] ==  patientIDVec[smooth_pairs_df$onIdx]
  if(badLastEntry){
    cond[numOfSmoothPairs] = FALSE;
  }
  smooth_pairs_df = smooth_pairs_df[cond,];

  if(!is.na(numberOfPairs)){
    if(numberOfPairs < nrow(smooth_pairs_df)){
      smpl_idx <- sort(sample(nrow(smooth_pairs_df), numberOfPairs))
      smooth_pairs_df <- smooth_pairs_df[smpl_idx,]
    }
  }

  smooth_pairs_df
}


########################################################################################
# SVM training
# res = trainSmoothDSS_quadSmoothness_diff(diffMat_ord, diffMat_smooth, diffT, timeVec, muO, orderingPairs,
#        smoothnessPairs, h_v, LambdaO_list[cntP], LambdaS_list[cntP]);
#
trainSmoothDSS_quadSmoothness_diff<-
  function(diffMatO, diffMatS, diffTVecS, hospTimeS, muO, pairsO, pairsS, h_v, Creg1, Creg2){

  # get dimensions
  numOfPairsO = dim(pairsO)[1];
  numOfPairsS = dim(pairsS)[1];
  dimOfFeatures = dim(diffMatO)[2];


  # rescaling Creg1 and Creg2
  normalizerC = Creg1
  Creg0 = 1/normalizerC;
  Creg1_n = (Creg1/numOfPairsO)/normalizerC
  Creg2_n = (Creg2/numOfPairsS)/normalizerC

  # loss function that allows vectors, currently it is Huber loss
  loss_fun_vec <- function(margin, indZ2, indZ3){

    # initialization
    res<- rep(0, length(margin))

    # In zone 1 it should be 0

    # if \mu-f>h then the loss grows linearly
    res[indZ3] <- margin[indZ3];

    # if |\mu-f|<= h then the loss is quadratic
    res[indZ2] <-(margin[indZ2] + h_v)^2/(4*h_v);

    res
  }

  # calculate matrices of differences and time deltas multiplied by Lipschitz constant
  #diffMatO <- calculateDiffMatrix(dataMatO, pairsO)
  #diffMatS <- calculateDiffMatrix(dataMatS, pairsS)
  #diffTVecS <- Lip*calculateDiffVec(hospTimeS, pairsS)
  #print(min(diffTVecS))
  invTVecS <- 1/diffTVecS

  w_init = rep(1,dimOfFeatures)
  # x = w_init
  # this function calculates target value for parameter x
  # as an attribute to its value it attaches gradient and hessian
  fgh <- function(x, ord){


    # get difference times parameter vector
    rankO = diffMatO%*%x;
    marginO = muO - rankO;
    indHuberZ1_O = marginO < -h_v;
    indHuberZ3_O = marginO > h_v;
    indHuberZ2_O = marginO <= h_v & marginO >= -h_v;
    penaltyO <- loss_fun_vec(marginO, indHuberZ2_O, indHuberZ3_O)



    rankS = diffMatS%*%x;
    #print(rankS)
    #print(invTVecS)
    rankS_overT = rankS*invTVecS;
    penaltyS = rankS_overT^2

    # calculate the number of elements in zones 2 and 3
    numZ2_O = sum(indHuberZ2_O)
    numZ3_O = sum(indHuberZ3_O)
    #print(penaltyS)
    # value calculation
    res <- Creg0*(1/2)*crossprod(x,x) + Creg1_n*sum(penaltyO) + Creg2_n*sum(penaltyS)
    res <- res[1]
    #print(res)
    #
    # calculation of the gradient
    #

    # regularization term
    gW = x;

    # ordering term
    auxO = (marginO + h_v)/(2*h_v);
    gO = 0;
    if(numZ3_O == 1){
      gO = gO - diffMatO[which(indHuberZ3_O),]
    }
    if(numZ3_O > 1){
      gO = gO - colSums(diffMatO[indHuberZ3_O,])
    }

    if(numZ2_O == 1){
      # Helping R to deal with vectors as matrices of size 1
      gO = gO - diffMatO[which(indHuberZ2_O),]*auxO[which(indHuberZ2_O)]
    }
    if(numZ2_O > 1){
      # Helping R to deal with vectors as matrices of size 1
      gO = gO - colSums(diffMatO[indHuberZ2_O,]*auxO[indHuberZ2_O])
    }

    # smoothness term
    gS = 2*colSums(diffMatS*as.vector(rankS_overT))

    grad_vec <- Creg0*gW + Creg1_n*gO + Creg2_n*gS;

    if(ord == 2){
      #
      # Calculating the Hessian - not currently implemented
      #
      featDim = length(x)
      HW = diag(1,featDim)


      z2Idx_O  = which(indHuberZ2_O)
      z2Idx_SP  = which(indHuberZ2_SP)
      z2Idx_SN  = which(indHuberZ2_SN)



      HO <- matrix(data = 0, ncol = featDim, nrow = featDim);
      for (cntP in z2Idx_O){
        curr_d = diffMatO[cntP, ];
        HO <- HO + curr_d%o%curr_d;
      }
      HO <- HO / (2*h_v)

      HSP <- matrix(data = 0, ncol = featDim, nrow = featDim);
      for (cntP in z2Idx_SP){
        curr_d = diffMatS[cntP, ];
        HSP <- HSP + curr_d%o%curr_d;
      }
      HSP <- HSP / (2*h_v)

      HSN <- matrix(data = 0, ncol = featDim, nrow = featDim);
      for (cntP in z2Idx_SN){
        curr_d = diffMatS[cntP, ];
        HSN <- HSN + curr_d%o%curr_d;
      }
      HSN <- HSN / (2*h_v)


      hss <- Creg0*HW + Creg1_n*HO + Creg2_n*(HSP + HSN)
    }
    if(ord >= 1){
      attr(res,"gradient") <- grad_vec
    }
    if(ord >= 2){
      attr(res,"hessian") <- hss
    }
    #

    res
  }



  # Rprof("out.out")
  #w <- nlm(f, w_init, iterlim = 40, print.level = 0, check.analyticals = TRUE)
  w <- nlm(fgh, w_init, ord = 1, print.level = 0, check.analyticals = FALSE,
           iterlim = 10000, stepmax = 10000, steptol = 1e-10, hessian = FALSE)

  debugMode = 0;
  if(debugMode){

    w0 <- nlm(fgh, w_init, ord = 0, print.level = 2, check.analyticals = FALSE, iterlim = 1, stepmax = 10000, steptol = 1e-8, hessian = FALSE)
    w_init <- w1$estimate
    w1 <- nlm(fgh, w_init, ord = 1, print.level = 2, check.analyticals = FALSE, iterlim = 1000, stepmax = 10000, steptol = 1e-8, hessian = FALSE)
    sum()
    w0
  }
  #
  #
  #
  #Rprof(NULL)
  #summaryRprof("out.out")

  w

}



calculateDiffMatrix <- function(dataMat, trainPairs){

  # get dimensions
  numOfTrainPairs = dim(trainPairs)[1];
  dimOfFeatures = dim(dataMat)[2];

  # calculate differences for all of the pairs
  # we assume that trainPairs are a collection of (i,j),
  # where i is always more severe than j
  diffMat = matrix(data=NA, nrow = numOfTrainPairs, ncol = dimOfFeatures);
  dataMat = as.matrix(dataMat)
  for(cntR in 1:numOfTrainPairs){
    i = trainPairs$onIdx[cntR];
    j = trainPairs$offIdx[cntR];
    diffMat[cntR, ] = dataMat[i,] -dataMat[j,];
  }

  diffMat
}

calculateDiffVec <- function(dataVec, trainPairs){

  # calculate differences for all of the pairs
  # we assume that trainPairs are a collection of (i,j),
  # where i is always more severe than j
  diffVec <- dataVec[trainPairs$onIdx] - dataVec[trainPairs$offIdx]
  diffVec
}


dssTrain.Linear <- function(dataMat_interpolated, orderingPairs, smoothnessPairs, timeVec,
                            LambdaO_list, LambdaS_list, save_fn_template, doParallel = 1)
{
 w_list = list();

  # calculate matrix of differences
  diffMat_ord <- calculateDiffMatrix(dataMat_interpolated, orderingPairs)
  diffMat_smooth <- calculateDiffMatrix(dataMat_interpolated, smoothnessPairs)
  diffT <- calculateDiffVec(timeVec, smoothnessPairs)

  h_v = 0.1
  numberOfPairs = dim(orderingPairs)[1];
  muO = rep(1, times = numberOfPairs)


  #for(cntP in seq(along = LambdaO_list)){
  if(doParallel){
    w_list <- foreach(cntP = seq(along = LambdaO_list))%dopar%{
      res = trainSmoothDSS_quadSmoothness_diff(diffMatO = diffMat_ord,
                                               diffMatS = diffMat_smooth,
                                               diffTVecS = diffT,
                                               hospTimeS = timeVec,
                                               muO = muO,
                                               pairsO = orderingPairs,
                                               pairsS = smoothnessPairs,
                                               h_v = h_v,
                                               Creg1 = LambdaO_list[cntP],
                                               Creg2 = LambdaS_list[cntP]);

      w_list[[cntP]] = res
      model_obj = res
      if(!is.na(save_fn_template)){
        save(model_obj, file = sprintf(save_fn_template, cntP))
      }
      res
    }
  }else{

    for(cntP in seq_along(LambdaO_list)){
      res = trainSmoothDSS_quadSmoothness_diff(diffMat_ord, diffMat_smooth, diffT, timeVec, muO, orderingPairs,
                                               smoothnessPairs, h_v, LambdaO_list[cntP], LambdaS_list[cntP]);
      w_list[[cntP]] = res
      model_obj = res
      if(!is.na(save_fn_template)){
        save(model_obj, file = sprintf(save_fn_template, cntP))
      }
    }
  }
  # return the list
  w_list

}


#
# Functions that are responsible for NL-DSS training
#
# This function trains severity score regressor and returns its parameters "w"
# based on data in
# dumb examples
# dataMat <- matrix(data = 1, ncol = 2, nrow = 3)
# dataMat <- rbind(c(1,2), c(2,3), c(6,8), c(5,9))
# trainPairs <- list(list(1,1), list(1,2), list(2, 3), list(3,2), list(3,1), list(4,3))
# h_val = 0.5 - parameter of Huber loss
# Creg = 1 - parameter of regularization
# w = trainSeverityScoreRegressor(truncDataMat_intp, trainPair_df, 0.5, 40000)
# dataMat = truncDataMat_intp; trainPairs = trainPair_df; h_val = 0.5; Creg = 40000
# dataMat = truncDataMat_intp; dataMat = truncDataMat_intp; hospTimeS = hosp_time; pairsO = trainPair_df; pairsS= smooth_pairs_df; h_v = 0.5; Creg1 = 1; Creg2 = 1
# w <- trainSmoothDSS(truncDataMat_intp, truncDataMat_intp, hosp_time, muO, trainPair_df, smooth_pairs_df, h_v, 1, 1, 1000)
#
#
# dataMat <- truncDataMat_intp; dataMat <- truncDataMat_intp; hospTimeS <- hosp_time;
# pairsO <- trainPair_df; pairsS <- trainSmoothPair_df;  Creg1 <- CregO[cntP];
# Creg2 <- CregS[cntP]; Lip <- LipList[cntP]
# res = trainSmoothDSS(truncDataMat_intp, truncDataMat_intp, hosp_time, muO, trainPair_df,
#       trainSmoothPair_df, h_v, CregO[cntP], CregS[cntP], LipList[cntP]);
#
# model <- trainSmoothDSS_gbrt_quadSmooth(dataMat, hospTime, muO, devCC_subsample, subsample_smooth_df,
#                                         h_v, 5e10, maxDepth = 6, maxTree = 200, maxStep = 1e6,
#                                         r_init = r_i)
# hospTimeS = hospTime; pairsS = subsample_smooth_df; pairsO = devCC_subsample; r_init = r_i; Creg2 = 48; maxTree = 200; maxDepth = 5; maxStep = 1e6
#model <- trainSmoothDSS_gbrt_quadSmooth(conc_dm_sorted, conc_ht_sorted, muO, devCC_subsample, subsample_smooth_df,
#                                        h_v, 250, maxDepth = 6, maxTree = 10, maxStep = 1e6,
#                                        r_init = r_i)
# dataMat_full = conc_dm_sorted; hospTimeS_full = conc_ht_sorted; pairsO_full = devCC_df; pairsS_full <- devSmoothness_df_1; r_init_full = r_i
# model <-
#   dssTrain.NonLinear(dataMat_full, hospTimeS_full, muO, pairsO_full, pairsS_full, h_v, 48, maxDepth = 5, maxTree = 5, r_init_full = rep(0, nrow(dataMat_full)) )
dssTrain.NonLinear <-
  function(dataMat_full,
           hospTimeS_full,
           muO,
           pairsO_full,
           pairsS_full,
           h_v,
           Creg2,
           r_init_full,
           maxDepth = 5,
           maxTree = 100,
           maxStep = 1e6,
           relImprovementStopTh = 1e-8,
           minStep = 1e-8,
           linSearchDamp = 0.5,
           interSave = 0,
           fn_inter_tmpl = NULL){

    # get dimensions
    numOfPairsO = dim(pairsO_full)[1];
    numOfPairsS = dim(pairsS_full)[1];

    # rescaling Creg1 and Creg2
    # for GBRT approach, C1 always equals 1
    Creg1_n = 1/numOfPairsO
    Creg2_n = Creg2/numOfPairsS
    if(0){
      Creg1_n = 1/numOfPairsO/Creg2
      Creg2_n = 1/numOfPairsS
    }

    # create a local copy of the reduced data matrix that contains only
    # entries that appear in in one of the pairs
    uniqEntries = unique(c(pairsO_full$onIdx, pairsO_full$offIdx, pairsS_full$onIdx, pairsS_full$offIdx))
    dataMat <- dataMat_full[uniqEntries, ]
    hospTimeS <- hospTimeS_full[uniqEntries]
    r_init <- r_init_full[uniqEntries]
    pairsO = data.frame(onIdx  = match(pairsO_full$onIdx, uniqEntries),
                        offIdx  = match(pairsO_full$offIdx, uniqEntries))
    pairsS = data.frame(onIdx  = match(pairsS_full$onIdx, uniqEntries),
                        offIdx  = match(pairsS_full$offIdx, uniqEntries))

    # dimensions of the
    dimOfFeatures = ncol(dataMat);
    numOfPatterns = nrow(dataMat);


    # Huber loss function that allows vector input
    # Here, margin is $mu - f_val$
    # indZ2 == 1 iff \size{margin} <= h_v
    # indZ3 == 1 iff margin > h_v
    loss_fun_vec <- function(margin, indZ2, indZ3){

      # initialization
      res<- rep(0, length(margin))

      # In zone 1 it should be 0

      # if \mu-f>h then the loss grows linearly
      res[indZ3] <- margin[indZ3];

      # if |\mu-f|<= h then the loss is quadratic
      res[indZ2] <-(margin[indZ2] + h_v)^2/(4*h_v);

      res
    }

    # Derivative of Huber loss function that allows vector input
    # Here, margin is $mu - f_val$
    # indZ2 == 1 iff \size{margin} <= h_v
    # indZ3 == 1 iff margin > h_v
    loss_deriv_vec <- function(margin, indZ2, indZ3){

      # initialization
      res<- rep(0, length(margin))

      # In zone 1 it should be 0

      # if \mu-f>h then the derivative is +1
      res[indZ3] <- 1;

      # if |\mu-f|<= h then the derivative is linear
      res[indZ2] <-(margin[indZ2] + h_v)^2/(2*h_v);

      res
    }


    # calculate matrices of differences and time deltas
    diffTVecS <- abs(calculateDiffVec(hospTimeS, pairsS))
    invDiffTVecS = 1/diffTVecS;
    sqInvDiffTVecS = (invDiffTVecS)^2


    # x = w_init
    # this function returns a list  of two elements
    # the first element is "val" - the value of the cost function
    # the second element is "grad" - the gradient vector with respect
    # to ranking of each point
    # the gradient is only calculated if (isGrad == 1)
    # rnk = currRnk; isGrad = 1
    fg <- function(rnk, isGrad){

      # initialize the output
      res = list(val = NA, grad = NA)

      # get differences of rankings for all the points
      rankO = calculateDiffVec(rnk, pairsO)
      marginO = muO - rankO
      indHuberZ1_O = marginO < -h_v;
      indHuberZ3_O = marginO > h_v;
      indHuberZ2_O = marginO <= h_v & marginO >= -h_v;
      penaltyO <- loss_fun_vec(marginO, indHuberZ2_O, indHuberZ3_O)


      # get differences times parameter vector
      rankS = calculateDiffVec(rnk, pairsS);
      penaltyS = (rankS)^2*sqInvDiffTVecS

      tmp <- Creg1_n*sum(penaltyO) + Creg2_n*sum(penaltyS);
      res$val = tmp[1]

      # if gradient calculation is not needed
      if(isGrad != 1){
        return(res)
      }

      #
      # Calculate the gradient with respect to every entry in the training data set
      # it is important that we are not playing with rankings of other points
      #
      derivO <- loss_deriv_vec(marginO, indHuberZ2_O, indHuberZ3_O)
      derivS <- 2*rankS*sqInvDiffTVecS


      # Fast calculation of per-entry derivative
      gradO = rep(0, numOfPatterns)
      i_com <- tapply(derivO, pairsO$onIdx, FUN = sum)
      i <- as.numeric(names(i_com))
      j_com <- tapply(derivO, pairsO$offIdx, FUN = sum)
      j <- as.numeric(names(j_com))
      gradO[i] = gradO[i] - i_com
      gradO[j] = gradO[j] + j_com

      gradS = rep(0, numOfPatterns)
      i_com <- tapply(derivS, pairsS$onIdx, FUN = sum)
      i <- as.numeric(names(i_com))
      j_com <- tapply(derivS, pairsS$offIdx, FUN = sum)
      j <- as.numeric(names(j_com))
      gradS[i] = gradS[i] + i_com
      gradS[j] = gradS[j] - j_com

      res$grad = Creg1_n*gradO + Creg2_n*gradS;

      #
      #

      res
    }

    # perform linear search
    dssGBRTLinearSearch <- function(v, currRnk, grad, descVec, aMax, aMin, tau, c, linSearchDamp){

      # init return value
      resList = list(nextV = NA, nextR = NA, conv = 1, a = NA);

      #
      curr_a <- aMax;
      m = crossprod(descVec, grad);

      repeat{
        nextR <- currRnk + curr_a*descVec;
        nextV <- fg(rnk = nextR, isGrad = 0); nextV <- nextV$val
        # v + curr_a*m*c - nextV
        if(nextV <= v + curr_a*m*c | curr_a <=aMin){
          break;
        }
        # take a smaller step
        curr_a <- curr_a*tau;
      }

      # add damping, i.e., take a smaller step than the one that is prescriberd by line search
      nextR <- currRnk + linSearchDamp*curr_a*descVec;
      nextV <- fg(rnk = nextR, isGrad = 0); nextV <- nextV$val


      # prepare and return the value
      resList$nextR = nextR;
      resList$nextV = nextV;
      resList$conv = ifelse(curr_a>aMin, 1, 0)
      resList$a = curr_a*linSearchDamp;
      resList

    }

    cat("Inside training function\n")
    # find all entires that are part of the training set
    #     uniqueTrainEntries = unique(c(pairsO$onIdx, pairsO$offIdx, pairsS$onIdx, pairsS$offIdx))
    #     isTrain = rep(0, numOfPatterns);
    #     isTrain[uniqueTrainEntries] = 1;
    #     numOfTrainPoints = sum(isTrain)
    rnk_init = r_init;
    #
    a_list = list(); # list of weights
    t_list = list(); # list of trees
    currRnk <- rnk_init;
    grad <- rep(0, numOfPatterns);
    dataMatWithRnk = cbind(dataMat, desc = -grad)

    # control parameters for the tree

    if(0){
      rcontrol <- rpart.control(minsplit = 2, cp = 1e-8,
                                maxcompete = 0, maxsurrogate = 0, usesurrogate = 0, xval = 0,
                                surrogatestyle = 0, maxdepth = maxDepth)
    }else{
      rcontrol <- rpart.control(minsplit = 100, cp = 1e-8,
                                maxcompete = 0, maxsurrogate = 0, usesurrogate = 0, xval = 0,
                                surrogatestyle = 0, maxdepth = maxDepth)
    }

    numTree_max = maxTree;
    currV <- 0;
    # repeat until convergence or until maximal number of trees reached
    repeat{

      # get the value and the gradient for the current ranking
      t1 <- Sys.time()
      vg_list <- fg(currRnk, isGrad = 1)
      Sys.time() - t1
      cat("Gradient direction found in ", Sys.time() - t1, "\n")
      # update the target value for the regression tree
      nof <- ncol(dataMatWithRnk)
      dataMatWithRnk$desc <- -vg_list$grad;

      t2 <- Sys.time()
      # fit a tree to the direction of gradient descend, only on training points
      currTree <- rpart( desc ~ . ,
                         data = dataMatWithRnk,
                         method = "anova", control = rcontrol, x= FALSE, y = FALSE)
      # remove data and pointer to the environment from the object thus drastically reducing its size
      currTree$where = c();
      attr(currTree$terms,".Environment") <- NULL
      cat("Gradient tree fitted\n")

      # calculate the actual direction of descend
      # (discrepancies from desc because of imperfect fit by regt.ree)
      descVec <- rep(0, numOfPatterns);
      descVec <- predict(currTree, dataMat, type = "vector")
      cat("Gradient approximation predicted\n")
      cat("Time: ", Sys.time() - t2, "\n")

      # find the best line step
      v <- vg_list$val; grad <- vg_list$grad; aMax = maxStep; aMin = minStep; tau = 1/2; c = 1/2
      resList <- dssGBRTLinearSearch(v, currRnk, grad, descVec, aMax, aMin, tau, c, linSearchDamp)
      currV <- v;

      if(resList$conv){
        # a step was found
        nextV <- resList$nextV
        nextR <- resList$nextR

        # report results
        cat("#Tree= ", length(a_list)+1, "Current value: ", currV, ", Next value: ", nextV, ", max|grad| = ", max(abs(vg_list$grad)),
            "Step a =", resList$a, " Relative improvement: ", (abs(currV-nextV)/max(abs(currV), abs(nextV))), "\n")

        # check stopping condition
        stopping_cond = (abs(currV-nextV)/max(abs(currV), abs(nextV))) < relImprovementStopTh | length(t_list) ==  numTree_max;

        # in case of stop
        if(is.na(stopping_cond) | stopping_cond){
          cat("The simulation converged\n")
          tmp <- list(a = a_list, t = t_list, uniqEntries = uniqEntries, currRnk = currRnk);
          return(tmp);
        }

      }else{
        cat("No descending direction was found - halt\n")
        # we are stuck - need to stop
        tmp <- list(a = a_list, t = t_list, uniqEntries = uniqEntries, currRnk = currRnk);
        return(tmp);
      }

      # save the tree, and proceed to next iteration
      t_list[[length(t_list)+1]] = currTree
      a_list[[length(a_list)+1]] = resList$a;
      currRnk <- nextR;
      currV <- nextV;

      if(mod(length(t_list), 5000) == 0){
        lbl <- length(t_list)/5000
        model_intermediate = list(t = t_list, a = a_list, uniqEntries = uniqEntries, currRnk = currRnk)
        save(model_intermediate, file = sprintf(fn_inter_tmpl, lbl))
      }

    }

    model = list(t = t_list, a = a_list, uniqEntries = uniqEntries, currRnk = currRnk)
  }

#
# just for testing
#
# registerDoMC(5); t1 <- Sys.time(); aaa <- dss.NL.calculateDSS(dataMat_full, model, 0); Sys.time() -t1
dss.NL.calculateDSS <- function(dataMat, res_forest, return_all){
  a_list = unlist(res_forest$a)
  t_list = res_forest$t

  numOfTrees = length(a_list)
  privateVecSum <- function(x,y){ (x+y) }
  res_rnk <- foreach(cntT = 1:numOfTrees, .combine = privateVecSum) %dopar% {
    cat(cntT,"\n")
    tmp <- a_list[[cntT]]*predict(t_list[[cntT]], dataMat, type = "vector")
    tmp
  }


  # return the result
  res_rnk
}



dss.NL.calculateDSS.Serial <- function(dataMat, res_forest, return_all){
  a_list = unlist(res_forest$a)
  t_list = res_forest$t

  numOfTrees = length(a_list)
  res_rnk <- rep(0, nrow(dataMat))
  if(return_all){
    res_rnk_all <- matrix(0, nrow = nrow(dataMat), ncol = length(a_list))
  }else{
    res_rnk_all = NULL
  }

  for(cntT in 1:numOfTrees){
    t1 <- Sys.time();
    tmp <- a_list[[cntT]]*predict(t_list[[cntT]], dataMat, type = "vector")
    res_rnk <- res_rnk + tmp
    if(return_all){
      res_rnk_all[,cntT] <- res_rnk
    }

    dt <- Sys.time() - t1; units(dt) <- "mins"
    cat("Tree", cntT,"out of", numOfTrees, "processed in", dt, "mins\n")
  }


  # return the result
  res_list = list(res_rnk = res_rnk, res_rnk_all = res_rnk_all)
  res_list
}


# Save the dss vector for all trees; for a single tree and for the soa maximizer on the cc-df
dss.NL.calculateDSS.Serial.ThreeVectors <- function(dataMat, res_forest, validCC_df){
  a_list = unlist(res_forest$a)
  t_list = res_forest$t

  numOfTrees = length(a_list)
  res_rnk <- rep(0, nrow(dataMat))
  res_rnk_all <- matrix(0, nrow = nrow(dataMat), ncol = 3)
  soa_maximum = -1;
  soa_maximizer = -1;


  for(cntT in 1:numOfTrees){
    t1 <- Sys.time();
    tmp <- a_list[[cntT]]*predict(t_list[[cntT]], dataMat, type = "vector")
    res_rnk <- res_rnk + tmp


    if(cntT == 1){
      res_rnk_all[,1] <- res_rnk
    }


    # calculate the soa
    cc_diff = calculateDiffVec(res_rnk, validCC_df);
    curr_soa = mean(cc_diff>=0)
    if(curr_soa> soa_maximum){
      res_rnk_all[,2] <- res_rnk
      soa_maximizer = cntT
      soa_maximum = curr_soa
    }


    if(cntT == numOfTrees){
      res_rnk_all[,3] = res_rnk
    }

    dt <- Sys.time() - t1; units(dt) <- "mins"
    cat("Tree", cntT,"out of", numOfTrees, "processed in", dt, "mins\n")
  }


  # return the result
  res_list = list(res_rnk = res_rnk, res_rnk_all = res_rnk_all, soa_maximizer = soa_maximizer, soa_maximum = soa_maximum)
  res_list
}





#
# Featurizing DSS
#
extractDSSFeaturesForEachPatient <- function(patD){
  numOfEntries  = nrow(patD)
  ht <- patD$ht
  ht_sq <- ht^2
  dss <- patD$dss
  lineNum = 1:numOfEntries
  cumsum_ht = cumsum(ht)
  cumsum_ht_sq = cumsum(ht_sq)
  cumsum_dss = cumsum(dss)
  cumsum_dss_ht = cumsum(dss*ht)
  cumsum_dss_ht_sq = cumsum(dss*ht_sq)
  abs_diff_dss = c(0,abs(diff(dss)));
  diff_ht = c(1, abs(diff(ht)))
  abs_dss_rate = abs_diff_dss/diff_ht
  cumsum_dss_abs_rate = cumsum(abs_dss_rate)
  cumsum_dss_abs_rate_ht = cumsum(abs_dss_rate*ht)
  diff_dss = c(0,diff(dss));
  dss_rate = diff_dss/diff_ht
  cumsum_dss_rate = cumsum(dss_rate)
  cumsum_dss_rate_ht = cumsum(dss_rate*ht)
  cummax_dss = cummax(dss)
  cummin_dss = cummin(dss)

  dss_avg = cumsum_dss/lineNum;
  dss_ht_avg = cumsum_dss_ht/cumsum_ht
  dss_ht_sq_avg = cumsum_dss_ht_sq/cumsum_ht_sq
  dss_abs_rate_avg = cumsum_dss_abs_rate/lineNum
  dss_abs_rate_ht_avg = cumsum_dss_abs_rate_ht/cumsum_ht
  dss_rate_avg = cumsum_dss_rate/lineNum
  dss_rate_ht_avg = cumsum_dss_rate_ht/cumsum_ht
  dss_max = cummax_dss

  res_df = data.frame(id = patD$id, ht = ht, dss = dss,
                      dss_avg = dss_avg, dss_ht_avg = dss_ht_avg,
                      dss_ht_sq_avg = dss_ht_sq_avg, dss_max = dss_max,
                      dss_rate_avg = dss_rate_avg, dss_rate_ht_avg = dss_rate_ht_avg,
                      dss_abs_rate_avg = dss_abs_rate_avg, dss_abs_rate_ht_avg = dss_abs_rate_ht_avg
  )

  res_df[is.na(res_df)] = 0;
  res_df
}

#dssWithDerivedFeatures <- extractDSSDerivedFeatures(icustayID, hospTime, dssRanks)
# patientIDVec = icustayID; timeVec = hospTime; probPat <- 22
# patD <- derDF_rgd[["22"]]
extractDSSDerivedFeatures <- function(patientIDVec, timeVec, dssRanks)
{
  derivationDF = data.frame(id = patientIDVec, ht = timeVec, dss = dssRanks)


  derDF_rgd <- split(derivationDF, as.factor(patientIDVec))


  res_list = list();
  for(patD in derDF_rgd){
    res_df <- extractDSSFeaturesForEachPatient(patD)
    res_list[[length(res_list)+1]] = res_df
  }
  dss_derived_feat <- rbind_all(res_list)

  dss_derived_feat

}



#
# Trains a logistic regression mapping from DSS+DerivedFeatures to adverse event using the glmnet package
# the training with the goal of maximizing AUC
# ---- ADD DESCRIPTION ----
# dssAndFeatures_obj = dssWithDerivedFeatures;
# user_level_label_valid - in lbl, '0' for negative, '1' for positive, '-1' to skip (censored), in id - user ids
#
# res_list <- trainLR.fromDSSAndDerivedFeaturesToAdverseEvent(dssWithDerivedFeatures, entrie_level_label, dev_ids, valid_ids, userInd)
#
# dssAndFeatures_obj = dssWithDerivedFeatures; user_level_label = userInd; lambda_list = logspace(100, 1e-7,11)
#
# lr_training_data <-
#trainLR.fromDSSAndDerivedFeaturesToAdverseEvent(dssWithDerivedFeatures, entrie_level_label,
#                                              dev_ids, valid_ids, userClassification,
#                                                lambda = logspace(-1, -3, 11))
# dssAndFeatures_obj <- dssWithDerivedFeatures; user_level_label <- userClassification; lambda_list = logspace(-1, -3, 11)

dss.Aux.TrainLRFromFeaturesToAdverseEvent <- function(ids, patternsMat, entrie_level_label, dev_ids, valid_ids, user_level_label,
                                                      lambda_list = logspace(2, -7,10)){

  # get the development and the validation set
  dev_idx = which(ids %in% dev_ids)
  valid_idx = which(ids %in% valid_ids)


  # train LR with L1-regularization (alpha = 1)
  mdl <- glmnet(x = patternsMat[dev_idx, ], entrie_level_label[dev_idx], family = "binomial",
                alpha = 0, lambda = lambda_list, intercept = TRUE)

  #
  # preliminary calculations for the validation set
  #

  uniq_usr_valid = unique(ids[valid_idx]); # unique users in validation set, in the same order as it ragged arrays

  mtch_idx = match(uniq_usr_valid, user_level_label$id) # fetch -1,0,1 labels for each user
  user_valid_lbl_all = user_level_label$lbl[mtch_idx]

  # remove users that have -1 label
  del_cond = user_valid_lbl_all == -1;
  user_valid_lbl <- user_valid_lbl_all[!del_cond]

  # compute the prediction on the validation set for all
  dssRanks_valid_auc_vec <- c()
  for(cnt in seq(along.with = mdl$lambda)){
    curr_lambda = mdl$lambda[cnt]
    tmp_ranks <- predict(mdl, patternsMat[valid_idx, ], s = curr_lambda, type = "response")
    usr_mx_score_valid_all <- tapply(tmp_ranks, ids[valid_idx], FUN = max);
    usr_mx_score_valid <- usr_mx_score_valid_all[!del_cond]

    roc_obj <- roc(user_valid_lbl, usr_mx_score_valid)
    dssRanks_valid_auc_vec[cnt]<- roc_obj$auc[1]
  }

  mx_idx = which(dssRanks_valid_auc_vec == max(dssRanks_valid_auc_vec))[1]
  best_lambda = mdl$lambda[mx_idx]
  best_auc <- dssRanks_valid_auc_vec[mx_idx]

  res_list <- list(mdl = mdl, best_lambda = best_lambda,
                   best_auc = best_auc, AUCs_valid = dssRanks_valid_auc_vec, best_idx = mx_idx)
  res_list

}


trainLR.fromDSSAndDerivedFeaturesToAdverseEvent <- function(dssAndFeatures_obj, entrie_level_label, dev_ids, valid_ids, user_level_label,
                                                lambda_list = logspace(2, -7,10)){

  # drop "id" and "ht" from the pattern and turn the result into a matrix
  patternsMat = as.matrix(dssAndFeatures_obj[,-c(1,2)])
  ids <- dssAndFeatures_obj$id;

  res_list <- dss.Aux.TrainLRFromFeaturesToAdverseEvent(ids, patternsMat,entrie_level_label, dev_ids, valid_ids, user_level_label,lambda_list)
  res_list

}

#
# This function will calculate the behavior of DSS prior to septic shock
# ids = icustayID; time_to_event = minToShock; score = dssRanks; test_ids = test_ids
dssTrendPriorToSepticShock <- function(ids, time_to_event, score, test_ids,
                                       typicalScoreRange, typicalInPatientScoreRange,
                                       basicWindowSize = 6*60,
                                       fn_delta1_plot
                                       ){

  # temporaty object
  tmp_df = data.frame(ids = ids, t2e = time_to_event, dss = as.vector(score))

  # truncate the temporary object
  cond = !is.na(tmp_df$t2e); tmp_df <- tmp_df[cond,]
  cond = tmp_df$t2e <= 3*basicWindowSize; tmp_df <- tmp_df[cond,]
  cond = tmp_df$ids %in% test_ids; tmp_df <- tmp_df[cond,]

  # indicators of different periods
  ind_1 = ifelse( tmp_df$t2e <= basicWindowSize,1,0)
  ind_2 = ifelse( tmp_df$t2e > basicWindowSize & tmp_df$t2e <= 2*basicWindowSize,1,0)
  ind_3 = ifelse( tmp_df$t2e > 2*basicWindowSize & tmp_df$t2e <= 3*basicWindowSize,1,0)

  # calculate the mean value of DSS for each patient in each time interval
  num_dss_1 = tapply(ind_1, tmp_df$ids, FUN = sum)
  num_dss_2 = tapply(ind_2, tmp_df$ids, FUN = sum)
  num_dss_3 = tapply(ind_3, tmp_df$ids, FUN = sum)
  sum_dss_1 = tapply(tmp_df$dss*ind_1, tmp_df$ids, FUN = sum)
  sum_dss_2 = tapply(tmp_df$dss*ind_2, tmp_df$ids, FUN = sum)
  sum_dss_3 = tapply(tmp_df$dss*ind_3, tmp_df$ids, FUN = sum)
  cond = !((num_dss_1 == 0) | (num_dss_2 == 0) |(num_dss_3 == 0))
  mean_dss_1 = sum_dss_1[cond]/num_dss_1[cond]
  mean_dss_2 = sum_dss_2[cond]/num_dss_2[cond]
  mean_dss_3 = sum_dss_3[cond]/num_dss_3[cond]

  delta_dss_2from1 = mean_dss_1 - mean_dss_2
  mean_delta_2from1 = mean(delta_dss_2from1)
  fractionOfPositive = mean(delta_dss_2from1 > 0)
  ci_fop = bootstrapCIforAccuracy(delta_dss_2from1 > 0)
  mean_delta_2from1_norm = mean_delta_2from1/typicalScoreRange

  # summarize the results into a list
  uniqIds = as.numeric(names(mean_dss_1))
  res_df <- data.frame(uniqIds = uniqIds, mean_dss_0_6 = mean_dss_1,
                       mean_dss_6_12 = mean_dss_2, mean_dss_12_18 = mean_dss_3)
  numEl = nrow(res_df);
  std_up = std(res_df$mean_dss_0_6 - res_df$mean_dss_6_12)
  std_acc = std(res_df$mean_dss_0_6 + res_df$mean_dss_12_18 - 2*res_df$mean_dss_6_12)


  sum_list = list();
  sum_list$mean_0_6 = mean(res_df$mean_dss_0_6);
  sum_list$mean_6_12 = mean(res_df$mean_dss_6_12);
  sum_list$mean_12_18 = mean(res_df$mean_dss_12_18);
  sum_list$delta_6_12 = sum_list$mean_0_6 - sum_list$mean_6_12;
  sum_list$delta_12_18 = sum_list$mean_6_12 - sum_list$mean_12_18;
  sum_list$delta_6_18 = sum_list$mean_0_6 - sum_list$mean_12_18;

  sum_list$p_trend_up = t.test(x = mean_dss_1, y = mean_dss_2, paired = T)$p.value
  sum_list$p_accel = t.test(x = mean_dss_1 -mean_dss_2,
                            y = mean_dss_2 - mean_dss_3, paired = T)$p.value

  sum_list$mean_delta1_norm_crossPat = mean_delta_2from1_norm;
  sum_list$mean_delta2_norm_crossPat = sum_list$delta_12_18/typicalScoreRange;
  sum_list$mean_delta3_norm_crossPat = sum_list$delta_6_18/typicalScoreRange;
  sum_list$mean_delta1_norm_inPat = mean(mean_dss_1)/typicalInPatientScoreRange;
  sum_list$mean_delta2_norm_inPat = sum_list$delta_12_18/typicalInPatientScoreRange;
  sum_list$mean_delta3_norm_inPat = sum_list$delta_6_18/typicalInPatientScoreRange;
  sum_list$fracOfPosDelta1_mean = fractionOfPositive;
  sum_list$fracOfPosDelta1_lowerB = ci_fop$basic[4]
  sum_list$fracOfPosDelta1_upperB = ci_fop$basic[5]
  sum_list$numberOfSamples = length(mean_dss_1)

  # plot the density of Delta1
  plot_df = data.frame(x = delta_dss_2from1)
  col <- gg_color_hue(4); col <- col[3]
  obj <- ggplot(data = plot_df, aes(x = x)) + geom_density(alpha = 0.3, color = col , fill = col) + geom_vline(xintercept = 0)
  obj <- obj + xlab("Value of Delta 1") + theme(text = element_text(size=20), plot.margin=unit(c(0.2,0.2,-0.0,0.0), "cm"))
  obj <- obj + ylab("Density of Delta 1")
  #obj <- obj + ggtitle("L-DSS")

  # save the actual plot
  pdf(file = fn_delta1_plot, width = 6, height = 6)
  print(obj)
  dev.off()


  list(sum_list = sum_list, res_df = res_df)


}





####################################################################
###################################################################
# Function that caclculates the optimal "r" for the training set.
# The returned values are "trainData" and "r_opt"
#dataMat_full = conc_dm_sorted;  hospTimeS_full = conc_ht_sorted; pairsO_full = devCC_df; pairsS_full = devSmoothness_df_1
#Creg2 = LambdaS_list[cnt]; r_init_full = rep(0, nrow(dataMat_full))
dssTrain.CalcOptimalDSSForTrainData <-
  function(dataMat_full, hospTimeS_full, muO, pairsO_full, pairsS_full, h_v, Creg2,
           r_init_full){

    # get dimensions
    numOfPairsO = dim(pairsO_full)[1];
    numOfPairsS = dim(pairsS_full)[1];



    # create a local copy of the reduced data matrix that contains only
    # entries that appear in in one of the pairs
    uniqEntries = unique(c(pairsO_full$onIdx, pairsO_full$offIdx, pairsS_full$onIdx, pairsS_full$offIdx))
    dataMat <- dataMat_full[uniqEntries, ]
    hospTimeS <- hospTimeS_full[uniqEntries]
    r_init <- r_init_full[uniqEntries]
    pairsO = data.frame(onIdx  = match(pairsO_full$onIdx, uniqEntries),
                        offIdx  = match(pairsO_full$offIdx, uniqEntries))
    pairsS = data.frame(onIdx  = match(pairsS_full$onIdx, uniqEntries),
                        offIdx  = match(pairsS_full$offIdx, uniqEntries))

    # dimensions of the
    dimOfFeatures = ncol(dataMat);
    numOfPatterns = nrow(dataMat);




    # rescaling Creg1 and Creg2
    # for GBRT approach, C1 always equals 1
    Creg1_n = 1/numOfPairsO
    Creg2_n = Creg2/numOfPairsS



    # Huber loss function that allows vector input
    # Here, margin is $mu - f_val$
    # indZ2 == 1 iff \size{margin} <= h_v
    # indZ3 == 1 iff margin > h_v
    loss_fun_vec <- function(margin, indZ2, indZ3){

      # initialization
      res<- rep(0, length(margin))

      # In zone 1 it should be 0

      # if \mu-f>h then the loss grows linearly
      res[indZ3] <- margin[indZ3];

      # if |\mu-f|<= h then the loss is quadratic
      res[indZ2] <-(margin[indZ2] + h_v)^2/(4*h_v);

      res
    }

    # Derivative of Huber loss function that allows vector input
    # Here, margin is $mu - f_val$
    # indZ2 == 1 iff \size{margin} <= h_v
    # indZ3 == 1 iff margin > h_v
    loss_deriv_vec <- function(margin, indZ2, indZ3){

      # initialization
      res<- rep(0, length(margin))

      # In zone 1 it should be 0

      # if \mu-f>h then the derivative is +1
      res[indZ3] <- 1;

      # if |\mu-f|<= h then the derivative is linear
      res[indZ2] <-(margin[indZ2] + h_v)^2/(2*h_v);

      res
    }


    # calculate matrices of differences and time deltas
    diffTVecS <- abs(calculateDiffVec(hospTimeS, pairsS))
    invDiffTVecS = 1/diffTVecS;
    sqInvDiffTVecS = (invDiffTVecS)^2


    # x = w_init
    # this function returns a list  of two elements
    # the first element is "val" - the value of the cost function
    # the second element is "grad" - the gradient vector with respect
    # to ranking of each point
    # the gradient is only calculated if (isGrad == 1)
    # rnk = r_init; isGrad = 1
    fg <- function(rnk, isGrad = 1){

      #cat("Function fg called\n")

      # initialize the output
      res = list(val = NA, grad = NA)

      # get differences of rankings for all the points
      rankO = calculateDiffVec(rnk, pairsO)
      marginO = muO - rankO
      indHuberZ1_O = marginO < -h_v;
      indHuberZ3_O = marginO > h_v;
      indHuberZ2_O = marginO <= h_v & marginO >= -h_v;
      penaltyO <- loss_fun_vec(marginO, indHuberZ2_O, indHuberZ3_O)


      # get differences times parameter vector
      rankS = calculateDiffVec(rnk, pairsS);
      penaltyS = (rankS)^2*sqInvDiffTVecS

      tmp <- Creg1_n*sum(penaltyO) + Creg2_n*sum(penaltyS);
      res$val = tmp[1]

      # if gradient calculation is not needed
      if(isGrad != 1){
        return(res)
      }

      #
      # Calculate the gradient with respect to every entry in the training data set
      # it is important that we are not playing with rankings of other points
      #
      derivO <- loss_deriv_vec(marginO, indHuberZ2_O, indHuberZ3_O)
      derivS <- 2*rankS*sqInvDiffTVecS


      # Fast calculation of per-entry derivative
      gradO = rep(0, numOfPatterns)
      i_com <- tapply(derivO, pairsO$onIdx, FUN = sum)
      i <- as.numeric(names(i_com))
      j_com <- tapply(derivO, pairsO$offIdx, FUN = sum)
      j <- as.numeric(names(j_com))
      gradO[i] = gradO[i] - i_com
      gradO[j] = gradO[j] + j_com

      gradS = rep(0, numOfPatterns)
      i_com <- tapply(derivS, pairsS$onIdx, FUN = sum)
      i <- as.numeric(names(i_com))
      j_com <- tapply(derivS, pairsS$offIdx, FUN = sum)
      j <- as.numeric(names(j_com))
      gradS[i] = gradS[i] + i_com
      gradS[j] = gradS[j] - j_com

      res$grad = Creg1_n*gradO + Creg2_n*gradS;

      res
    }



    cat("Inside training function\n")
    rnk_init = r_init;


    # perform linear search
    dssGBRTLinearSearch <- function(v, currRnk, grad, descVec, aMax, aMin, tau, c, linSearchDamp){

      # init return value
      resList = list(nextV = NA, nextR = NA, conv = 1, a = NA);

      #
      curr_a <- aMax;
      m = crossprod(descVec, grad);

      repeat{
        nextR <- currRnk + curr_a*descVec;
        nextV <- fg(rnk = nextR, isGrad = 0); nextV <- nextV$val

        if(nextV <= v + curr_a*m*c | curr_a <=aMin){
          break;
        }
        # take a smaller step
        curr_a <- curr_a*tau;
      }

      # add damping, i.e., take a smaller step than the one that is prescriberd by line search
      nextR <- currRnk + linSearchDamp*curr_a*descVec;
      nextV <- fg(rnk = nextR, isGrad = 0); nextV <- nextV$val


      # prepare and return the value
      resList$nextR = nextR;
      resList$nextV = nextV;
      resList$conv = ifelse(curr_a>aMin, 1, 0)
      resList$a = curr_a*linSearchDamp;
      resList

    }

    cat("Inside training function\n")

    currRnk <- r_init;
    grad <- rep(0, numOfPatterns);

    currV <- 0;
    maxIter = 200;
    cntIter = 0;
    relImprovementStopTh = 1e-8

    # repeat until convergence or until maximal number of iterations reached
    repeat{

      cntIter = cntIter + 1

      # get the value and the gradient for the current ranking
      t1 <- Sys.time()
      vg_list <- fg(currRnk, isGrad = 1)
      cat("Gradient direction found in ", Sys.time() - t1, "\n")

      # find the best line step
      v <- vg_list$val; grad <- vg_list$grad; aMax = 10^6; aMin = 1e-6; tau = 1/2; c = 1/2; descVec <- -grad
      resList <- dssGBRTLinearSearch(v, currRnk, grad, descVec, aMax, aMin, tau, c, linSearchDamp = 1)
      currV <- v;

      if(resList$conv){
        # a step was found
        nextV <- resList$nextV
        nextR <- resList$nextR

        # report results
        cat("#Interation = ", cntIter, "Current value: ", currV, ", Next value: ", nextV, ", max|grad| = ", max(abs(vg_list$grad)),
            "Step a =", resList$a, " Relative improvement: ", (abs(currV-nextV)/max(abs(currV), abs(nextV))), "\n")

        # check stopping condition
        stopping_cond = (abs(currV-nextV)/max(abs(currV), abs(nextV))) < relImprovementStopTh | cntIter ==  maxIter;

        # in case of stop
        if(stopping_cond){
          cat("The simulation converged\n")
          opt_scoring <- list(trainData = dataMat, dss_opt = nextR)
          return(opt_scoring);
        }

      }else{
        cat("No descending direction was found - halt\n")
        # we are stuck - need to stop
        opt_scoring <- list(trainData = dataMat, dss_opt = currRnk)
        return(opt_scoring);
      }

      currRnk <- nextR;
      currV <- nextV;
    }







    #
    # Here
    #

    #r_opt <- nlm(fg, isGrad = 0, r_init, print.level = 2, check.analyticals = TRUE, iterlim = 10, stepmax = 10000, steptol = 1e-8, hessian = FALSE)

#     f <- function(x){cat("Function f called\n"); tmp <- fg(x, 0); tmp}
#     g <- function(x){tmp <- fg(x, 1); cat("Function g called; value = ", tmp[1], "\n");tmp2 <- attr(tmp, "gradient"); tmp2}
#
#     r_opt <- optim(r_init, f, g, method = c("CG"), hessian = FALSE, control = list(trace = 0, maxit = 2))

    opt_scoring <- list(trainData = dataMat, dss_opt = currRnk, pairsO = pairsO,
                        pairsS = pairsS, uniqEntries = uniqEntries)
    opt_scoring

  }



###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
### Implementing analysis of the response to treatment ############################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################


dss.Aux.CalculateTimeWeightedAverage <- function(val, ht){

  # check the validity of inputs
  if(length(val) != length(ht)){
    cat("dss.Aux.CalculateTimeWeightedAverage: arguments are of different length\n")
    stop()
  }

  # if only one value
  if(length(val) == 1){
    return(val[1])
  }

  # if several
  val_but_last <-val[1:(length(val)-1)]
  dv <- diff(val)/2; dht <- diff(ht)

  return(sum((val_but_last + dv)*dht)/sum(dht))


}

# # This function analyzes data for a single patient
# # For every "clean" treatment administration point (no treatments "wind" minutes before
# # and no treatments "wind" minutes after),
# # we calculate the average value of DSS prior and after this event point.
# dss.Auxiliary.CalculateDSSResponseToTreatment <- function(id, ht, dss, treat, wind){
#
#
#
#   #
#   res_list <- list(); numOfEntries = length(ht)
#   for(cntP in seq_along(treat)){
#     # get current point of interest
#     curr_idx = treat[cntP]
#
#     # it's id and ht
#     curr_ht = ht[curr_idx]; curr_id = id[curr_idx]
#
#     # find relevant points prior to the current index
#     iter_bck = curr_idx;
#     while(ifelse(iter_bck>1, (curr_ht - ht[iter_bck-1] <= wind) & (id[iter_bck-1] == curr_id), FALSE)){
#       iter_bck = iter_bck -1;
#     }
#
#     # find relevant points after the current index
#     iter_fwd = curr_idx;
#     while(ifelse(iter_fwd < numOfEntries,  (ht[iter_fwd+1] - curr_ht <= wind) & (id[iter_fwd + 1] == curr_id), FALSE)){
#       iter_fwd = iter_fwd +1;
#     }
#
#     #avgBefore <- dss.Aux.CalculateTimeWeightedAverage(dss[iter_bck:curr_idx], ht[iter_bck:curr_idx])
#     #avgAfter <- dss.Aux.CalculateTimeWeightedAverage(dss[(curr_idx+1):iter_fwd], ht[(curr_idx+1):iter_fwd])
#     avgBefore <- mean(dss[iter_bck:curr_idx]); avgAfter <-mean(dss[(curr_idx+1):iter_fwd])
#     if(iter_fwd > curr_idx & iter_bck < curr_idx){
#       tmp_df =
#         data.frame(numBefore = (curr_idx - iter_bck + 1), numAfter = (iter_fwd - curr_idx),
#                    avgBefore = avgBefore, avgAfter = avgAfter)
#       res_list[[length(res_list)+1]] = tmp_df
#     }
#
#   }
#   res_df  = rbind_all(res_list)
#   res_df
# }


# This function analyzes data for a single patient
# For every "clean" treatment administration point (no treatments "wind" minutes before
# and no treatments "wind" minutes after),
# we calculate the average value of DSS prior and after this event point.
dss.Auxiliary.CalculateDSSResponseToTreatment <- function(id, ht, dss, treat, wind){



  #
  res_list <- list(); numOfEntries = length(ht)
  for(cntP in seq_along(treat)){
    # get current point of interest
    curr_idx = treat[cntP]

    # it's id and ht
    curr_ht = ht[curr_idx]; curr_id = id[curr_idx]

    # find relevant points prior to the current index
#     iter_bck = curr_idx;
#     while(ifelse(iter_bck>1, (curr_ht - ht[iter_bck-1] <= wind) & (id[iter_bck-1] == curr_id), FALSE)){
#       iter_bck = iter_bck -1;
#     }

    # find relevant points after the current index
    iter_fwd = curr_idx;
    while(ifelse(iter_fwd < numOfEntries,  (ht[iter_fwd+1] - curr_ht <= wind) & (id[iter_fwd + 1] == curr_id), FALSE)){
      iter_fwd = iter_fwd +1;
    }

    #avgBefore <- dss.Aux.CalculateTimeWeightedAverage(dss[iter_bck:curr_idx], ht[iter_bck:curr_idx])
    #avgAfter <- dss.Aux.CalculateTimeWeightedAverage(dss[(curr_idx+1):iter_fwd], ht[(curr_idx+1):iter_fwd])
    #avgBefore <- mean(dss[iter_bck:curr_idx]); avgAfter <-mean(dss[(curr_idx+1):iter_fwd])
    #avgBefore <- dss[curr_idx]; avgAfter <- mean(dss[(curr_idx+1):iter_fwd])
    if(iter_fwd > curr_idx){
      avgBefore <- dss[(curr_idx)];
      avgAfter <- dss.Aux.CalculateTimeWeightedAverage(dss[(curr_idx):iter_fwd], ht[(curr_idx):iter_fwd])
      tmp_df =
        data.frame(numBefore = 1, numAfter = (iter_fwd - curr_idx + 1),
                   avgBefore = avgBefore, avgAfter = avgAfter)

    }else{
      tmp_df = data.frame(numBefore = NA, numAfter = NA, avgBefore = NA, avgAfter = NA)
    }
    res_list[[length(res_list)+1]] = tmp_df
  }
  res_df  = rbind_all(res_list)
  res_df
}



# This function analyzes data for a single patient
# For every "clean" treatment administration point (no treatments "wind" minutes before
# and no treatments "wind" minutes after),
# we calculate the average value of DSS prior and after this event point.
dss.Auxiliary.CalculateDSSResponseToTreatment.PrePost <- function(id, ht, dss, treat, wind){

  #
  res_list <- list(); numOfEntries = length(ht)
  for(cntP in seq_along(treat)){
    # get current point of interest
    curr_idx = treat[cntP]

    # it's id and ht
    curr_ht = ht[curr_idx]; curr_id = id[curr_idx]

    # find relevant points prior to the current index
    iter_bck = curr_idx;
    while(ifelse(iter_bck>1, (curr_ht - ht[iter_bck-1] <= wind) & (id[iter_bck-1] == curr_id), FALSE)){
      iter_bck = iter_bck -1;
    }

    # find relevant points after the current index
    iter_fwd = curr_idx;
    while(ifelse(iter_fwd < numOfEntries,  (ht[iter_fwd+1] - curr_ht <= wind) & (id[iter_fwd + 1] == curr_id), FALSE)){
      iter_fwd = iter_fwd +1;
    }

    #avgBefore <- dss.Aux.CalculateTimeWeightedAverage(dss[iter_bck:curr_idx], ht[iter_bck:curr_idx])
    #avgAfter <- dss.Aux.CalculateTimeWeightedAverage(dss[(curr_idx+1):iter_fwd], ht[(curr_idx+1):iter_fwd])
    #avgBefore <- mean(dss[iter_bck:curr_idx]); avgAfter <-mean(dss[(curr_idx+1):iter_fwd])
    #avgBefore <- dss[curr_idx]; avgAfter <- mean(dss[(curr_idx+1):iter_fwd])
    if(iter_fwd > curr_idx & iter_bck < curr_idx){
      avgBefore <- dss.Aux.CalculateTimeWeightedAverage(dss[iter_bck:curr_idx], ht[iter_bck:curr_idx])
      avgAfter <- dss.Aux.CalculateTimeWeightedAverage(dss[(curr_idx):iter_fwd], ht[(curr_idx):iter_fwd])
      avgNow <- dss[curr_idx]
      tmp_df =
        data.frame(numBefore = (curr_idx - iter_bck + 1), numAfter = (iter_fwd - curr_idx),
                   avgBefore = avgBefore, avgAfter = avgAfter, avgNow  = avgNow,
                   trendBefore = avgNow - avgBefore, trendAfter = avgAfter - avgNow)

    }else{
      tmp_df = data.frame(numBefore = NA, numAfter = NA, avgBefore = NA, avgAfter = NA,
                          avgNow  = NA, trendBefore = NA, trendAfter = NA)
    }
    res_list[[length(res_list)+1]] = tmp_df
  }
  res_df  = rbind_all(res_list)
  res_df
}



# Train the treatment prediction model used for control;
dss.ResponseToTreatment.Preliminary <- function(icustayID, dataMat, treatmentTimeData, dev_ids, valid_ids, test_ids, wind, sepsis_related_cond){

  if(class(dataMat) != "matrix"){
    cat("dss.ResponseToTreatment.Preliminary: Input argument dataMat is expected to be a matrix and not", class(dataMat),"\n")
    stop("dss.ResponseToTreatment.Preliminary: input arguments are of the wrong format")
  }

  # classify data entries into positive (treatment within an hour)
#   pos_examples = (!is.na(treatmentTimeData$timeSinceLastEvent) & treatmentTimeData$timeSinceLastEvent <= 60) |
#     (!is.na(treatmentTimeData$timeToNextEvent) & treatmentTimeData$timeToNextEvent <= 60) | treatmentTimeData$indicator

  pos_examples = (!is.na(treatmentTimeData$timeToNextEvent) & treatmentTimeData$timeToNextEvent <= 60) | treatmentTimeData$indicator == 1
  pos_ex_idx = which(pos_examples == 1 & sepsis_related_cond)
  neg_ex_idx_all = which(pos_examples == 0 & sepsis_related_cond)
  neg_ex_idx = sample(neg_ex_idx_all, 10*length(pos_ex_idx))

  # find indices of dev, valid, and test data poitns
  dev_idx = which(icustayID %in% dev_ids)
  valid_idx = which(icustayID %in% valid_ids)
  test_idx = which(icustayID %in% test_ids)

  # find positive and negative examples for development set
  dev_pos_idx = intersect(dev_idx, pos_ex_idx)
  dev_neg_idx = intersect(dev_idx, neg_ex_idx);
  dev_all_idx = union(dev_pos_idx, dev_neg_idx)
  dev_labels = ifelse(dev_all_idx %in% dev_pos_idx, 1, 0)

  # train the predictor on the development set for a variety of regularization values
  lambda_list = logspace(0, -9, 11)
  alpha_list = linspace(0,1,11)
  #alpha_list = c(0)
  treat_model_list = list();
  t1 <- Sys.time()
  for(cntA in seq_along(alpha_list)){

    treat_model <- glmnet(x = dataMat[dev_all_idx, ],
                          dev_labels, family = "binomial",
                          alpha = alpha_list[cntA], lambda = lambda_list,
                          intercept = TRUE)
    treat_model_list[[cntA]] <- treat_model

    # reporting
    dt <- Sys.time() - t1; units(dt) <- "mins"
    cat("Treatment model training, cntA = ", cntA,"out of", length(alpha_list), "calculated in",dt, "mins\n")
    t1 <- Sys.time()
  }

  # evaluate the performance by calculating the AUC on the validation set
  valid_pos_idx = intersect(valid_idx, pos_ex_idx)
  valid_neg_idx = intersect(valid_idx, neg_ex_idx);  # we are taking ALL negative examples, since there is no need for balancing
  valid_all_idx = union(valid_pos_idx, valid_neg_idx)
  valid_labels = ifelse(valid_all_idx %in% valid_pos_idx, 1, 0)
  AUC_valid <- matrix(NA, ncol = length(lambda_list), nrow = length(alpha_list))
  t1 <- Sys.time()
  for(cntA in seq_along(alpha_list)){
    for(cntL in seq(along = lambda_list)){
      currL = lambda_list[cntL];
      valid_pred <- predict(treat_model_list[[cntA]], dataMat[valid_all_idx, ], s = currL, type = "response")
      # roc_obj <- roc(valid_labels, as.vector(valid_pred))
      #AUC_valid[cntA, cntL] <- roc_obj$auc[1]
      AUC_valid[cntA, cntL] <- performance(prediction(valid_pred, valid_labels), "auc")@y.values[[1]]
      # reporting
      dt <- Sys.time() - t1; units(dt) <- "mins"
      cat("Estimating accuracy of trained treatment models. For cntA = ", cntA,"out of", length(alpha_list), "and CntL =", cntL, "out of", length(lambda_list), "evaluated in",dt, "mins\n")
      t1 <- Sys.time()
    }
  }
  best_AUC = max(AUC_valid)
  best_ind = which(AUC_valid == best_AUC, arr.ind = TRUE)[1,]
  best_Alpha = alpha_list[best_ind[1]]
  best_Lambda = lambda_list[best_ind[2]]
  bestModel = treat_model_list[[best_ind[1]]]
  #  bestModel$beta[,best_ind[2]]

  # will be returned
  treatmentModelData = list(best_AUC = best_AUC, best_Alpha = best_Alpha, best_Lambda = best_Lambda, bestModel = bestModel)


  # predict treatment probability for all entries.
  treat_prob <- predict(bestModel, dataMat, s =  best_Lambda, type = "response")
  #th_75 <- quantile(treat_prob, c(0.5))[1]
  #th_75 <- quantile(treat_prob[treatmentTimeData$indicator == 1], c(0.5))[1]; th_75
  # th <- th_75
  th <- quantile(treat_prob[sepsis_related_cond], c(0.90))[1]

  hor = wind;
  # find the baseline examples
  cond_high_prob = treat_prob >= th & sepsis_related_cond;
  cond_test = icustayID %in% test_ids;
  cond_no_CMI_after = (treatmentTimeData$timeToNextEvent >= 2*hor) | is.na(treatmentTimeData$timeToNextEvent)
  cond_no_CMI_before = (treatmentTimeData$timeSinceLastEvent >= 2*hor) | is.na(treatmentTimeData$timeSinceLastEvent)
  baseline_idxs_untilSS = which(cond_no_CMI_after & cond_no_CMI_before & cond_high_prob & cond_test & treatmentTimeData$indicator == 0)

  # positive examples
  #pos_idxs_untilSS = which(cond_no_CMI_after & cond_no_CMI_before & cond_high_prob & cond_test & treatmentTimeData$indicator == 1)
  # all_positive = which(sepsis_related_cond & cond_test & treatmentTimeData$indicator == 1)
  pos_idxs_untilSS = which(cond_high_prob & cond_test & treatmentTimeData$indicator == 1)

  # baseline examples
  baseline_idxs_untilSS_subsample <- sample(baseline_idxs_untilSS, length(pos_idxs_untilSS), replace = FALSE)

  #
  retObj <- list(treatmentModelData = treatmentModelData, treat_ex_idx = pos_idxs_untilSS, baseline_ex_idx = baseline_idxs_untilSS_subsample,
                 baseline_idx_all = baseline_idxs_untilSS)
  # treatmentResponseAnalysisObj <- retObj

  retObj

}


dss.PerformanceAnalysis.AnalyzeTreatmentResponse <- function(icustayID, ht, dssRanks,
                                                             treatmentResponseAnalysisObj, wind,
                                                             fn_plot_treatmentResponse)
{
  # get the actual tables
  treat_stat <- dss.Auxiliary.CalculateDSSResponseToTreatment(icustayID, ht, dssRanks,
                                                              treatmentResponseAnalysisObj$treat_ex_idx,
                                                              wind)

  # debugging
  if(0){

    dssRanks <- (match(sepsisIndLin$status, c("none", "sirs", "severe", "shock")) -1)[until_ss]
    treat_stat <- dss.Auxiliary.CalculateDSSResponseToTreatment(icustayID, ht, dssRanks,
                                                                all_positive,
                                                                6*60)
   cat("Before = ", mean(treat_stat$avgBefore), ", After = ", mean(treat_stat$avgAfter))
  }

  # debugging
  if(0){
    all_positive = which(bolus_indicator == 1)
    #dssRanks <- truncDataMat_intp$anbp_sys
    dssRanks <- (match(sepsisIndLin$status, c("none", "sirs", "severe", "shock")) -1)
    wind = 12*60
    treat_stat <- dss.Auxiliary.CalculateDSSResponseToTreatment(icustay_id, hosp_time, dssRanks,
                                                                all_positive,
                                                                wind)
    cat("Before = ", mean(treat_stat$avgBefore), ", After = ", mean(treat_stat$avgAfter))
  }

  min_sbms = min(c(max(c(length(treatmentResponseAnalysisObj$treat_ex_idx), 2000)), length(treatmentResponseAnalysisObj$baseline_idx_all)))
  indx_smbs = sample(length(treatmentResponseAnalysisObj$baseline_idx_all), min_sbms, replace = FALSE)
  baseline_cmi_stat <- dss.Auxiliary.CalculateDSSResponseToTreatment(icustayID, ht, dssRanks,
                                                                     treatmentResponseAnalysisObj$baseline_idx_all[indx_smbs],
                                                                     wind)

#   baseline_cmi_stat <- dss.Auxiliary.CalculateDSSResponseToTreatment(icustayID, ht, dssRanks,
#                                                                      treatmentResponseAnalysisObj$baseline_idx_all,
#                                                                      wind)

  #
  before_cmi = mean(treat_stat$avgBefore)
  after_cmi = mean(treat_stat$avgAfter)
  before_no = mean(baseline_cmi_stat$avgBefore)
  after_no = mean(baseline_cmi_stat$avgAfter)
  delta_cmi = treat_stat$avgAfter - treat_stat$avgBefore
  delta_no = baseline_cmi_stat$avgAfter-baseline_cmi_stat$avgBefore
  #min_l = min(c(length(delta_cmi), length(delta_no)))
  #delta_tmp = delta_no[1:min_l] - delta_cmi[1:min_l]
  #p1 = pnorm( mean(delta_tmp)/(std(delta_tmp)/sqrt(length(delta_tmp))),lower.tail=FALSE)
  p1 = t.test(x = delta_cmi, y = delta_no, alternative = "t")$p.value
  numOfSamples = length(delta_cmi)
  fractionBelowControlMean = mean(delta_cmi <= mean(delta_no))
  ci_fbc = bootstrapCIforAccuracy(delta_cmi <= mean(delta_no))


  # plot the delta's to a file
  cols = gg_color_hue(2)
  lbl = c(rep(1, length(delta_cmi)), rep(2, length(delta_no)))
  uniq_stages = c("Treated", "Control")
#   plot_df = data.frame(scores = c(delta_cmi, delta_no),
#                        stage = uniq_stages[lbl],
#                        cols = cols[lbl]);
  stg = uniq_stages[lbl]

  plot_df = data.frame(scores = c(delta_cmi, delta_no),
                     stage = factor(stg, level = uniq_stages, ordered = TRUE),
                     cols = cols[lbl]);
  plotTitle = "L-DSS"
  # generate the actual plot
  obj <- ggplot(plot_df, aes(x = scores, fill = stage))
  obj <- obj + geom_density(alpha = 0.1)
  obj <-obj +
    scale_fill_manual(values = cols, guide=guide_legend(title='Group'),
                      labels = uniq_stages) +
    scale_color_manual(values=cols, guide=guide_legend(title='Detection status'),
                       labels = uniq_stages)
  obj <- obj + xlab("Value of Delta DSS")
  obj <- obj + ylab("Density of Delta DSS")
  obj <- obj + ggtitle(plotTitle)

  pdf(file = fn_plot_treatmentResponse, width = 12, height = 8)
  print(obj)
  dev.off()


  res_list = list(treat_stat = treat_stat, baseline_cmi_stat = baseline_cmi_stat,
                  before_treat = before_cmi, after_treat = after_cmi,
                  before_baseline = before_no, after_baseline = after_no,
                  meanDeltaTreatment = mean(delta_cmi),
                  meanDeltaControl = mean(delta_no),
                  #delta_table = delta_tmp,
                  p_treat = p1,
                  numOfSamples = numOfSamples,
                  fractionBelowControlMean = fractionBelowControlMean,
                  fractionBelowControlMean_lb = ci_fbc$basic[4],
                  fractionBelowControlMean_ub = ci_fbc$basic[5],
                  delta_no = delta_no, delta_cmi = delta_cmi
                  )

  res_list
}




dss.PerformanceAnalysis.PlotScoreDensityPerStage <- function(
  scores,
  sevStagePerEntry, sevStageToLabelMapping, plotTitle,
  fn_plot_save,
  fn_liftCurve){

  # get the number of unique stages
  numOfStages = nrow(sevStageToLabelMapping)
  cols=gg_color_hue(numOfStages)

  # translate severity stages into row numbers in the mapping data frame
  sevStages_inRows = match(sevStagePerEntry, sevStageToLabelMapping$sevStage)

  # plot data frames
  plot_df = data.frame(DSS = scores,
                       stage = sevStageToLabelMapping$sevStageLabel[sevStages_inRows],
                       cols = cols[sevStages_inRows]);

  # generate the actual plot
  obj <- ggplot(plot_df, aes(x = DSS, fill = factor(stage, level = sevStageToLabelMapping$sevStageLabel, ordered = TRUE)))
  obj <- obj + geom_density(alpha = 0.1) + theme(text = element_text(size=20), legend.position = "right")
  obj <-obj +
    scale_fill_manual(values=cols, guide=guide_legend(title='Stage'),
                      labels = sevStageToLabelMapping$sevStageLabel) +
    scale_color_manual(values=cols,guide=guide_legend(title='Detection status'),
                       labels = sevStageToLabelMapping$sevStageLabel)
  obj <- obj + xlab("Value of DSS")
  obj <- obj + ylab("Density of DSS")
  obj <- obj# + ggtitle(plotTitle)

  # save the actual plot
  pdf(file = fn_plot_save, width = 10, height = 6)
  print(obj)
  dev.off()

  # get the boundaries
  # dss_boundaries = quantile(scores, c(0.02, 0.92))
  # dss_lb = dss_boundaries[1]
  # dss_up = dss_boundaries[2]

#   dss_lb = quantile(scores[sevStagePerEntry < 4], c(0.01))
#   dss_up = quantile(scores[sevStagePerEntry == 3], c(0.98))
#   #bin_boundaries = linspace(dss_lb, dss_up, 101)
#   bin_boundaries = c(min(scores) - 1, bin_boundaries);
#   bin_boundaries[length(bin_boundaries)+1] = max(scores) + 1;
#   bin_boundaries = linspace(min(scores), max(scores), 101)
#
#
#   histInfo_notShock <- hist(scores[sevStagePerEntry != 3], breaks = bin_boundaries)
#   cumProb_notShock = cumsum(histInfo_notShock$counts)
#   cumProb_notShock_comp = sum(histInfo_notShock$counts) - cumsum(histInfo_notShock$counts)
#
#   histInfo_Shock <- hist(scores[sevStagePerEntry == 3], bin_boundaries)
#   cumProb_Shock = cumsum(histInfo_Shock$counts)
#   cumProb_Shock_comp = sum(histInfo_Shock$counts) - cumsum(histInfo_Shock$counts)
#
#   #oddsRatio = (cumProb_Shock_comp/cumProb_Shock)/(cumProb_notShock_comp/cumProb_notShock)
#   oddsRatio = (cumProb_Shock_comp/cumProb_notShock_comp)/(sum(histInfo_Shock$counts)/sum(histInfo_notShock$counts))
#   qplot(x = bin_boundaries[-1], y = oddsRatio)
#
#   #cumProb_notShock = cumsum(histInfo_notShock$density)/sum(histInfo_notShock$density)
#   odds_notShock = (1-cumProb_notShock)/cumProb_notShock
#   qplot(x = bin_boundaries[-1], y = log(odds_notShock))
#
#   #cumProb_Shock = cumsum(histInfo_Shock$density)/sum(histInfo_Shock$density)
#   odds_Shock = (1-cumProb_Shock)/cumProb_Shock
#   qplot(x = bin_boundaries[-1], y = log(odds_Shock))
#
#   oddsRatio = odds_Shock/odds_notShock
#   qplot(x = bin_boundaries[2:(length(bin_boundaries)-1)], y = log(oddsRatio[1:(length(oddsRatio)-1)]))
#
#
  # approximation

#   mean_shock = mean(scores[sevStagePerEntry == 3]); std_shock = std(scores[sevStagePerEntry == 3])
#   mean_noShock = mean(scores[sevStagePerEntry < 3]);
#   std_noShock = std(scores[sevStagePerEntry < 3])

  dss_lb = quantile(scores[sevStagePerEntry < 3], c(0.01))
  dss_up = quantile(scores[sevStagePerEntry == 3], c(0.98))
  bin_boundaries = linspace(dss_lb, dss_up, 1000)
  bin_boundaries = c(min(scores), bin_boundaries, max(scores))

  histInfo_notShock <- hist(scores[sevStagePerEntry != 3], breaks = bin_boundaries)
  cumProb_notShock = cumsum(histInfo_notShock$counts)
  cumProb_notShock_comp = sum(histInfo_notShock$counts) - cumsum(histInfo_notShock$counts)

  histInfo_Shock <- hist(scores[sevStagePerEntry == 3], bin_boundaries)
  cumProb_Shock = cumsum(histInfo_Shock$counts)
  cumProb_Shock_comp = sum(histInfo_Shock$counts) - cumsum(histInfo_Shock$counts)

  lift = (cumProb_Shock_comp/(cumProb_notShock_comp + cumProb_Shock_comp))/mean(sevStagePerEntry == 3)
  plot_df = data.frame(x = bin_boundaries[2:(length(lift))], y = lift[1:(length(lift)-1)]);
  obj <- ggplot(data = plot_df) + geom_line(aes(x = x, y = y))
  obj <- obj + xlab("Value of DSS")
  obj <- obj + ylab("Lift of fraction of Septic Shock points")
  obj <- obj + ggtitle("L-DSS")

  # save the actual plot
  pdf(file = fn_liftCurve, width = 12, height = 8)
  print(obj)
  dev.off()

#   print(obj)
#
#   x11(); qplot(x = bin_boundaries[2:(length(lift))],  y = lift[1:(length(lift)-1)])
#
#   odds_shock_app = (1-pnorm((bin_boundaries - mean_shock)/std_shock))/pnorm((bin_boundaries - mean_shock)/std_shock)
#   x11(); qplot(x = bin_boundaries, y = log(odds_shock_app))
#   odds_noshock_app = (1- pnorm((bin_boundaries - mean_noShock)/std_noShock))/pnorm((bin_boundaries - mean_noShock)/std_noShock)
#   x11(); qplot(x = bin_boundaries, y = log(odds_noshock_app))
#   odds_ratio_app = odds_shock_app/odds_noshock_app
#   x11(); qplot(x = bin_boundaries, y = log(odds_ratio_app))
#
#
#   log_likelihood_ratio = log(dnorm(x = bin_boundaries, mean = mean_shock, sd = std_shock)) - log(dnorm(x = bin_boundaries, mean = mean_noShock, sd = std_noShock))
#   x11(); qplot(x = bin_boundaries, y = log_likelihood_ratio)

}




dss.Aux.PValueOfPairedPositivityTestVsNumberOfSamples <- function(deltaVec, numOfSizePoints, numOfRepetitions = 100, file_name,
                                                                  xbreaks = NULL, ybreaks = NULL){

  # the total number of examples
  numOfSamples = length(deltaVec);
  valuesToTry = round(logseq(x1 = 2, x2 = numOfSamples, numOfSizePoints)); valuesToTry = unique(valuesToTry)
  resultsMat = matrix(NA, nrow = numOfRepetitions, ncol=length(valuesToTry))
  for(cntV in seq_along(valuesToTry)){
    currV = valuesToTry[cntV]
    for(cnt in 1:numOfRepetitions){
      currSample = sample(numOfSamples, currV, replace = TRUE)
      currData = deltaVec[currSample];

      if(length(unique(currData))>1){
        resultsMat[cnt, cntV] = t.test(x = deltaVec[currSample])$p.value
        if(mean(deltaVec[currSample])<0){
          resultsMat[cnt, cntV] = 1
        }
      }else{
        resultsMat[cnt, cntV] = 1
      }
    }
  }

  # extract mean and 95% CI
  meanPval = rep(0, length(valuesToTry))
  pVal_95ub = rep(0, length(valuesToTry))
  pVal_95lb = rep(0, length(valuesToTry))
  for(cntV in seq_along(valuesToTry)){
    meanPval[cntV] = median(resultsMat[,cntV])
    pVal_95ub[cntV] = quantile(resultsMat[,cntV], 0.95)
    pVal_95lb[cntV] = quantile(resultsMat[,cntV], 0.05)
  }

  ymin_val = min(pVal_95lb)
  x_onepercent = valuesToTry[which(meanPval < 0.01)[1]]

  if(is.null(xbreaks)){
    xbreaks = logseq(1, 10^ceil(log10(max(valuesToTry))), ceil(log10(max(valuesToTry)))+1);
  }
  xbreaks = sort(unique(c(xbreaks,  x_onepercent)))

  if(is.null(ybreaks)){
    ybreaks = logseq(1, 10^(4*floor(log10(ymin_val)/4)), abs(floor(log10(ymin_val)/4))+1);
  }
  ybreaks = sort(unique(c(ybreaks, 0.01)))

  # plot the graph
  plot_df <- data.frame(x = valuesToTry, y = meanPval, lb = pVal_95lb, ub = pVal_95ub)
  obj <- ggplot(plot_df, aes(x = x, y= y))+  geom_point()+ geom_line()+
    geom_ribbon(data=plot_df,aes(ymin=lb,ymax=ub),alpha=0.1) +
    scale_y_log10(limits = c(ymin_val, 1), breaks = ybreaks) + scale_x_log10(breaks = xbreaks, limits = c(1, 2*max(xbreaks))) +
    geom_vline(xintercept = x_onepercent, linetype = "longdash") + geom_hline(yintercept = 0.01, linetype = "longdash") +
    xlab("#Samples") + ylab("p-value") + theme(text = element_text(size=20), plot.margin=unit(c(0.2,0.2,-0.0,0.0), "cm"))
  pdf(file_name, height = 4, width = 6)
  print(obj)
  dev.off()

  res_list = list(meanPval = meanPval, pVal_95ub = pVal_95ub, pVal_95lb = pVal_95lb)
}



dss.Aux.PlotPrecisionRecall <- function(curr_roc_obj, file_name, xmin, xmax){

  # data frame for precision/recall/f-measure
  num_of_th = length(curr_roc_obj$thresholds)
  x = rep(curr_roc_obj$thresholds, 3)
  y = c(curr_roc_obj$precision, curr_roc_obj$recall, curr_roc_obj$f_measure)
  lbl = c(rep("Precision", num_of_th), rep("Recall", num_of_th), rep("F-Measure", num_of_th))
  precRecallFMeas_df = data.frame(x = x, y = y, label = lbl)
  obj <- ggplot(precRecallFMeas_df, aes(x = x, y = y, group = label, colour = label)) + geom_line()
  obj <- obj + theme(text = element_text(size=20), legend.position = "top")
  obj <- obj + ylim(0, 1) + xlab("Threshold") + ylab("") + xlim(xmin, xmax);

  pdf(file_name, height = 4, width = 6)
  print(obj)
  dev.off()

  1
}



dss.Aux.PlotSensitivityOneMinusSpecificity <- function(curr_roc_obj, file_name, xmin, xmax){

  # data frame for precision/recall/f-measure
  num_of_th = length(curr_roc_obj$thresholds)
  x = rep(curr_roc_obj$thresholds, 2)
  y = c(curr_roc_obj$sensitivities, 1 - curr_roc_obj$specificities)
  lbl = c(rep("Sensitivity", num_of_th), rep("1-Specificity", num_of_th))
  precRecallFMeas_df = data.frame(x = x, y = y, label = lbl)
  obj <- ggplot(precRecallFMeas_df, aes(x = x, y = y, group = label, colour = label)) + geom_line()
  obj <- obj + theme(text = element_text(size=20), legend.position = "top")
  obj <- obj + ylim(0, 1) + xlab("Threshold") + ylab("") + xlim(xmin, xmax);

  pdf(file_name, height = 4, width = 6)
  print(obj)
  dev.off()

  1
}


dss.Aux.PlotSensitivitySpecificity <- function(curr_roc_obj, file_name, xmin, xmax){

  # data frame for precision/recall/f-measure
  num_of_th = length(curr_roc_obj$thresholds)
  x = rep(curr_roc_obj$thresholds, 2)
  y = c(curr_roc_obj$sensitivities, curr_roc_obj$specificities)
  lbl = c(rep("Sensitivity", num_of_th), rep("Specificity", num_of_th))
  precRecallFMeas_df = data.frame(x = x, y = y, label = lbl)
  obj <- ggplot(precRecallFMeas_df, aes(x = x, y = y, group = label, colour = label)) + geom_line()
  obj <- obj + theme(text = element_text(size=20), legend.position = "top")
  obj <- obj + ylim(0, 1) + xlab("Threshold") + ylab("") + xlim(xmin, xmax);

  pdf(file_name, height = 4, width = 6)
  print(obj)
  dev.off()

  1
}


# this function slices GBRT model to smaller models with at most numOfTreesPerSlice in each.
# This way we can calculate their DSS faster
dss.Aux.SliceGBRT <- function(model_list, numOfTreesPerSlice){

  cntSlice = 0;
  cntAll = 0;
  sliced_models_list = list()
  bm_model = c()
  bm_slice = c();
  for(cntM in seq_along(model_list)){
    currModel = model_list[[cntM]]
    currNumOfT = length(currModel$t)
    currNumOfSlices = ceil(currNumOfT/numOfTreesPerSlice)
    for(cntSlice in 1:currNumOfSlices){
      slice_st = (cntSlice - 1)*numOfTreesPerSlice + 1
      slice_end = min(c(cntSlice*numOfTreesPerSlice, currNumOfT))
      slice_idxs = slice_st:slice_end
      #cat("Start = ", slice_st, ", End = ", slice_end, "\n")

      cntAll = cntAll + 1
      sliced_models_list[[cntAll]] <- list(a = currModel$a[slice_idxs], t = currModel$t[slice_idxs])
      bm_model[cntAll] <- cntM;
      bm_slice[cntAll] <- cntSlice;

    }
  }
  backwardMapping <- list(model = bm_model, slice = bm_slice)
  res_list = list(sliced_models_list = sliced_models_list, backwardMapping = backwardMapping)
  res_list
}


#
#
dss.Aux.CollectSlices <- function(backwardMapping, fn_modelSlices_template){

  # load files
  currModel = -1;
  curr_dss <- NULL;
  model_number_list = c()
  dss_list = list();
  for(cnt in seq_along(backwardMapping$model)){
    # load the data
    load(sprintf(fn_modelSlices_template, backwardMapping$model[cnt], backwardMapping$slice[cnt]))

    # first model ever - save results
    if(currModel == -1){
      currModel = backwardMapping$model[cnt];
      curr_dss <- dssRanks;
    }else{

      # not the first slice
      if(currModel != backwardMapping$model[cnt]){
        # new model
        dss_list[[length(dss_list)+1]] = curr_dss
        model_number_list[length(model_number_list)+1] = currModel;
        currModel = backwardMapping$model[cnt]
        curr_dss <- dssRanks;
      }else{
        curr_dss = curr_dss + dssRanks;
      }
    }
  }

  # remember to store the dss from the last model.
  dss_list[[length(dss_list)+1]] = curr_dss
  model_number_list[length(model_number_list)+1] = currModel;
  res_list = list(dss_list = dss_list, model_number = model_number_list)

  res_list
}


