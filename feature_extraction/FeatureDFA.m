% Extract Feature - Detrended Fluctuation Analysis
%
% 8 July 2013
% -------------------------------------------------------------------------

function DFASignal = FeatureDFA(x)

DFASignal = fastdfa(x);
