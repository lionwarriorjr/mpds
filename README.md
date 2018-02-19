# mpds
Source Code for Learning the Mobile Parkinson Disease Score (mPDS)

## feature_extraction
This folder contains the Matlab feature library used to extract the relevant smartphone sensor features including but not limited to measurements from the accelerometer, in addition to those related to audio and screen touch. See paper "High Frequency Remote Monitoring of Parkinson's Disease via Smartphone: Platform Overview and Medication Response Detection" (at https://arxiv.org/abs/1601.00960) for more detail. 

## mPDS
This folder contains the Jupyter Notebook used to learn the mPDS in our accompanying paper using a smartphone based dataset derived from the feature extraction library above and DSSL, a ranked based machine learning framework described in "Learning (Predictive) Risk Scores in the Presence of Censoring due to Interventions." See the associated paper at https://arxiv.org/abs/1507.07295 for more detail. Included as well is a pair of useful visualization tools (in subfolder visualization tools) used to generate the learned severity trajectories as well as more easily visualize the features DSSL weights as most influential in the learning procedure.

## dssl
This folder contains the R source code for DSSL (Disease Severity Score Learning).

