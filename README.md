# mpds
Mobile Parkinson Disease Score

## feature_extraction
This folder contains the Matlab feature library to extract Parkinson disease related features from raw data samples collected from HopkinsPD Android Application.

## data_preparation
This folder contains the Jupyter Notebook to load feature files into a Postgres database called HopkinsPD. The tables in HopkinsPD database will be served for mPDS training.

## dssl
This folder contains the R source code of DSSL (Disease Severity Score Learning) and the Jupyter Notebook script to train the mPDS model using DSSL.

## analysis
This folder contains the R scripts used to generate results for each major component of the paper.
Specifically, one script is used for preprocessing of learned mPDS scores for the cross-sectional
analysis, another IPython notebook carries out the cross-sectional analysis, another generates
visualizations of the mPDS trajectories in the development cohort, while another computes
the relevant summary statistics reported in the main text.

## data
This folder contains the relevant data files referenced in the analysis.
