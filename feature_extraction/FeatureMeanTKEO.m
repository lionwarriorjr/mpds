% Extract Feature - Mean TKEO
%
% 8 July 2013
% -------------------------------------------------------------------------

function MeanTKEOSignal = FeatureMeanTKEO(x)

m = length(x);
MeanTKEOSignal = mean((x(2:m-1,1).^2) - x(3:m,1).*x(1:m-2,1));

