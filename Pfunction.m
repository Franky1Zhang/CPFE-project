function [Pvals] = Pfunction(paravec, xdata, sigma, maturities)
%PFUNCTION is the function to calculate the price of Treasury Bond
%   Produce the bond prices given the input of parameters， initial interest
%   rate, sigma and maturities

rbar = paravec(1);
gamma = paravec(2);

Avals = Afunction(rbar, gamma, sigma, maturities);
Bvals = Bfunction(gamma, maturities);

Pvals = zeros(length(xdata), length(maturities));% T by d matrix

Pvals= 100 * exp(Avals- Bvals.*xdata);


end

