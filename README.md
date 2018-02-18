# mpds
Mobile Parkinson Disease Score

## feature_extraction
This folder contains the Matlab feature library to extract Parkinson disease related features from raw data samples collected from HopkinsPD Android Application.

## mPDS
This folder contains the Jupyter Notebook to load feature files into a Postgres database called HopkinsPD. The tables in HopkinsPD database will be served for mPDS training.

## dssl
This folder contains the R source code of DSSL (Disease Severity Score Learning) and the Jupyter Notebook script to train the mPDS model using DSSL.
