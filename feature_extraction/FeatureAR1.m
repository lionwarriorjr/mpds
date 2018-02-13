% Extract Feature - AR1 coefficient
%
% 8 July 2013
% -------------------------------------------------------------------------

function AR_coeff = FeatureAR1(x)

out = x(2:length(x));
in  = x(1:length(x)-1);

AR_coeff = regress(out,in);
