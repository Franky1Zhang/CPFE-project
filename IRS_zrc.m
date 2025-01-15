clc
clear
close all
dataFile = 'China Treasury Spot Yields.xlsx';

maturitiesInMonths = xlsread(dataFile, 'A1:M1');
maturitiesInYears = maturitiesInMonths / 12;
spotRates = 0.01 * xlsread(dataFile, 'A2:M206'); 

[numObs, numBonds] = size(spotRates);

%% a
instantrate = spotRates(:,1);

drt = instantrate(2:end) - instantrate(1:end-1);
dt = 1/252;
constant = ones(length(instantrate(1:end-1)),1);
[b,bint,r,rint,stats] = regress(drt,[constant,instantrate(1:end-1)]);

gamma = -b(2) / dt;
rbar = b(1) / gamma / dt;
sigma = sqrt(stats(4)) / sqrt(dt);
disp(gamma);
disp(rbar);
disp(sigma);

%% b 
% Initial Guess of Parameters
rbar_star = 0.04; 
gamma_star = 1.24;

paravec = [rbar_star, gamma_star];

xdata = instantrate;

myfun = @(parameters, xdata) Pfunction(parameters, xdata, sigma, maturitiesInYears);
options = optimset('MaxIter', 2000, 'MaxFunEvals', 2000);

bondPrices = zeros(numObs, numBonds);
for i = 1:numBonds
    bondPrices(:, i) = 100 ./ exp(spotRates(:, i) * maturitiesInYears(i));
end

[estimatedParas] = lsqcurvefit(myfun, paravec, xdata, bondPrices, [0 0], [10 4], options);
rbar_star = estimatedParas(1);
gamma_star = estimatedParas(2);
disp(gamma_star);
disp(rbar_star);

%% c
r0 = 1.8430 / 100;
fixrate = 0.03582;
FV = 1e6;
paravec = [rbar_star, gamma_star, sigma];

IRS_price = IRS_price_zrc(r0, fixrate, FV, paravec);
disp(IRS_price);

%% d
rng('default');

timeStepPerPeriod = 252;
simulationHorizon = 5;
totalNumSteps = timeStepPerPeriod * simulationHorizon;
numSimulation = 10000;
deltat = 1 / timeStepPerPeriod;

% 半年支付一次，持有IRS直到半年后
holdHorizon = 0.5; 
leftHorizon = simulationHorizon - holdHorizon;
holdNumSteps = holdHorizon * timeStepPerPeriod;
leftNumSteps = leftHorizon * timeStepPerPeriod;

holdSimNum = 1000;
leftSimNum = 1000;

% 存储持有期结束时的浮动利率
EndRate = zeros(holdSimNum, 1);

rpath_h = zeros(1, holdNumSteps); % 持有期的利率路径
rpath_l = zeros(1, leftNumSteps); % 剩余期的利率路径

% 存储每次模拟的IRS价格
allEndPrices = zeros(holdSimNum, leftSimNum);

% 存储剩余期的折现因子
discountFactor_l = zeros(1, leftNumSteps);

% 半年支付相关变量
rt_l = zeros(leftSimNum, leftHorizon * 2); % 剩余期每半年支付的浮动利率
dt_l = zeros(leftSimNum, leftHorizon * 2); % 剩余期每半年支付的折现因子
cashFlow_l = zeros(leftSimNum, leftHorizon * 2); % 剩余期每半年支付的现金流

% ----------------------------------------------------------------------
% 计算持有期结束时可能的浮动利率并存储到 EndRate 数组中
% ----------------------------------------------------------------------
for i = 1:holdSimNum
    randVars_h = randn(1, holdNumSteps);
    oldr = r0;
    for j = 1:holdNumSteps
        dr = gamma * (rbar - oldr) * deltat + sigma * sqrt(deltat) * randVars_h(j);
        rpath_h(j) = oldr + dr;
        oldr = rpath_h(j);        
    end
    EndRate(i) = rpath_h(holdNumSteps);
end

% ----------------------------------------------------------------------
% 计算持有期结束后IRS的所有可能价格
% ----------------------------------------------------------------------

for i = 1:holdSimNum
    for j = 1:leftSimNum
        randVars_l = randn(1, leftNumSteps);
        oldr = EndRate(i);
        olddiscountFactor = 1;
        for k = 1:leftNumSteps
            dr = gamma_star * (rbar_star - oldr) * deltat + sigma * sqrt(deltat) * randVars_l(k);
            rpath_l(k) = oldr + dr;
            discountFactor_l(k) = exp(-oldr * deltat) * olddiscountFactor;
            olddiscountFactor = discountFactor_l(k);
            oldr = rpath_l(k);
        end
        rt_l(j, :) = rpath_l((leftNumSteps / (leftHorizon * 2)) * [1 : leftHorizon * 2]);
        dt_l(j, :) = discountFactor_l((leftNumSteps / (leftHorizon * 2)) * [1 : leftHorizon * 2]);
        cashFlow_l(j, :) = FV * 0.5 * (0.03582 - 0.0001 * 0.8750 - rt_l(j, :)); % 半年现金流
        allEndPrices(i, j) = cashFlow_l(j, :) * dt_l(j, :)';
    end
end

endPrices = reshape(allEndPrices, holdSimNum * leftSimNum, 1);

VaR_1 = mean(endPrices) - prctile(endPrices, 1);
VaR_5 = mean(endPrices) - prctile(endPrices, 5);
disp(VaR_1);
disp(VaR_5);
