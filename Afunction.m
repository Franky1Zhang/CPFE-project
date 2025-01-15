function Avals = Afunction(rbar_star, gamma_star, sigma, maturities)
%
% Function returns A(maturity) in the Vasicek model.
% P = exp(A(maturity) - B(maturity)*r)
%
% with the pseudo risk-neutral SDE parameterization 
%
% dr = gamma_star (rbar-star - r)dt + \sigma dX
%
% 'maturities' can be a vector.  If so, 'Avals' is also a vector. 
%
sigma2 = sigma^2;

Bvals = (1-exp(-gamma_star.*maturities))/gamma_star; 

Term1=rbar_star-sigma2/(2*gamma_star^2);

Term2 = -sigma2*Bvals.^2/(4*gamma_star); 

Avals = (Bvals-maturities)*Term1+Term2; 

return
