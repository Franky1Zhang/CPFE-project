function Bvals = Bfunction(gamma_star,maturities)
%
% Function returns B(maturity) in the Vasicek model.
% P = exp(A(maturity) - B(maturity)*r)
%
% with the pseudo risk-neutral SDE parameterization 
%
% dr = gamma_star*(rbar_star - r)dt + sigma dX
%
% 'maturities' can be a vector.  If so, 'Bvals' is also a vector. 
%

Bvals = (1-exp(-gamma_star.*maturities))/gamma_star; 


return
