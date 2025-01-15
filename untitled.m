% 参数设定
clc; clear;

% Read data from 'China Treasury Spot Yields.xlsx'

dataFile = 'China Treasury Spot Yields.xlsx';
maturitiesInMonth = xlsread(dataFile, 'A1:M1');
maturitiesInYear = maturitiesInMonth / 12;

spotRates = 0.01 .* xlsread(dataFile, 'A2:M206');
[numObs numBonds] = size(spotRates);
spotRates = spotRates(:,1:12);

IRS_maturityInYear = 5;
IRS_maturityInMonth = IRS_maturityInYear * 12;

%% 1 - ZCB Prices

% Calculate zero coupon bond prices of each maturity in time series
tempMaturity = repmat(maturitiesInYear, numObs, 1);
priceMatrix = 100 * exp( - spotRates.* tempMaturity);
meanPrice = mean(priceMatrix);% calculate average of each columns
stdPrice = std(priceMatrix);

fprintf(1,'\n\n');
fprintf(1,'\ Mean and Std-Dev of Prices Derived from Spot Rates\n');
fprintf(1,'      Maturities         Mean         Std\n');
fprintf(1,'----------------------------------------------------\n');

n = min([length(maturitiesInMonth), length(meanPrice), length(stdPrice)]);
for i = 1:n
    fprintf(1, '%6d  %15.2f  %10.2f\n', maturitiesInMonth(i), meanPrice(i), stdPrice(i));
end


%% 2 - Risk Natural Parameter Estimates

% Parameter estimates for Vasicek Model
TreasSpotRate_1Month = spotRates(:,1);
FirstDifference = TreasSpotRate_1Month(2:end) - TreasSpotRate_1Month(1:end-1);
constant = ones(length(FirstDifference),1);

[b bint r rint stats] = regress(FirstDifference, [constant, TreasSpotRate_1Month(1:end-1)]);

% This leads to the following risk natrual parameters 
dt=1/252; 
alpha=b(1);
beta=b(2);
gamma=-beta/dt;
residSD=sqrt(stats(4));
sigma=residSD/sqrt(dt);
rbar=alpha/(gamma*dt);

%% 3 - Risk Neutral Parameter Estimates


rbar_star = 0.02;
gamma_star = 1.24;

paravec = [rbar_star gamma_star];

xdata = spotRates(:,1);

myfun = @(parameters, xdata) Pfunction(parameters,xdata,sigma,maturitiesInYear);

options = optimset('MaxIter',2000,'MaxFunEvals',2000);

[estimatedParas] = lsqcurvefit(myfun,paravec,xdata,priceMatrix,[0 0 0],[10 4 4],options);

rbar_star = estimatedParas(1);
gamma_star = estimatedParas(2);

fprintf('rbar_star_hat equals to %5.5f\n', estimatedParas(1));
fprintf('gamma_star_hat equals to %5.5f\n', estimatedParas(2));

%% 4 - Pricing by Monte Carlo Simulations
% For IRS, the Monte Carlo simulation methodology calls for the
% contemporaneous simulations of cash flows and discounts.

r0 = 1.8430 / 100; % instantaneous rate of floating rate path
FaceVal = 1e6;
fixedRate = 0.03582;

IRSPrice = irsPricing(r0,fixedRate,FaceVal,estimatedParas(1),estimatedParas(2),sigma);

fprintf('Price of the Interest Rate Swap equals to %5.5f\n', IRSPrice);

%% 5 - VaR in 0.5 year (半年付一次)

rng('default');

timeStepPerPeriod = 252; % 每年252个时间步
simulationHorizon = 5; % 模拟5年
totalNumSteps = timeStepPerPeriod * simulationHorizon; % 总步数
numSimulation = 10000;
deltat = 1 / timeStepPerPeriod;

holdHorizon = 0.5; % 持有期0.5年
leftHorizon = simulationHorizon - holdHorizon; % 剩余投资期
holdNumSteps = holdHorizon * timeStepPerPeriod; % 持有期时间步
leftNumSteps = leftHorizon * timeStepPerPeriod; % 剩余期时间步

holdSimNum = 1000; % 持有期模拟数量
leftSimNum = 1000; % 剩余期模拟数量

% 半年支付一次
paymentFrequency = 2; % 每年2次支付（半年付一次）

% 存储持有期结束时的浮动利率
EndRate = zeros(holdSimNum,1);

rpath_h = zeros(1,holdNumSteps); % 持有期利率路径
rpath_l = zeros(1,leftNumSteps); % 剩余期利率路径

% 存储每次模拟的价格
allEndPrices = zeros(holdSimNum,leftSimNum);

% 存储剩余期的折现因子
discountFactor_l = zeros(1,leftNumSteps);

% 存储每次支付的浮动利率
rt_l = zeros(leftSimNum,leftHorizon * paymentFrequency);
% 存储每次支付的折现因子
dt_l = zeros(leftSimNum,leftHorizon * paymentFrequency);
% 存储每次支付的现金流
cashFlow_l = zeros(leftSimNum,leftHorizon * paymentFrequency);

% ----------------------------------------------------------------------
% 计算持有期结束时可能的浮动利率并存储到 EndRate 数组中
% ----------------------------------------------------------------------
for i = 1 : holdSimNum
    randVars_h = randn(1,holdNumSteps);
    for j = 1 : holdNumSteps
        if j == 1
            dr = gamma * (rbar - r0) * deltat + sigma * sqrt(deltat) * randVars_h(j);
            rpath_h(j) = r0 + dr;
        else
            dr = gamma * (rbar - rpath_h(j-1)) * deltat + sigma * sqrt(deltat) * randVars_h(j);
            rpath_h(j) = rpath_h(j-1) + dr;
        end
    end
    EndRate(i) = rpath_h(holdNumSteps);
end

% ----------------------------------------------------------------------
% 计算持有期结束后 IRS 的所有可能价格
% ----------------------------------------------------------------------

for i = 1 : holdSimNum
    for j = 1 : leftSimNum
        randVars_l = randn(1,leftNumSteps);
        for k = 1 : leftNumSteps
            if k == 1
                dr = gamma_star * (rbar_star - EndRate(i)) * deltat + sigma * sqrt(deltat) * randVars_l(k);
                rpath_l(k) = EndRate(i) + dr;
                discountFactor_l(k) = exp(-EndRate(i) * deltat);
            else
                dr = gamma_star * (rbar_star - rpath_l(k-1)) * deltat + sigma * sqrt(deltat) * randVars_l(k);
                rpath_l(k) = rpath_l(k-1) + dr;
                discountFactor_l(k) = exp(-rpath_l(k-1) * deltat) * discountFactor_l(k-1);
            end
        end
        rt_l(j,:) = rpath_l((leftNumSteps / (leftHorizon * paymentFrequency)) * [1 : leftHorizon * paymentFrequency]);
        dt_l(j,:) = discountFactor_l((leftNumSteps / (leftHorizon * paymentFrequency)) * [1 : leftHorizon * paymentFrequency]);
        cashFlow_l(j,:) = FaceVal * 0.5 * (0.03582 - 0.0001 * 0.8750 - rt_l(j,:)); % 半年支付一次，现金流调整为0.5
        allEndPrices(i,j) = cashFlow_l(j,:) * dt_l(j,:)';
    end
end

endPrices = reshape(allEndPrices, holdSimNum * leftSimNum, 1);

sortedEndP = sort(endPrices);
% 1% 最坏情况下的投资组合价值
onePercent = sortedEndP(length(endPrices) * 0.01); disp(onePercent);
% 5% 最坏情况下的投资组合价值
fivePercent = sortedEndP(length(endPrices) * 0.05); disp(fivePercent);

% 1% VaR 和 5% VaR
onePercentVaR = mean(endPrices) - onePercent; disp(onePercentVaR);
fivePercentVaR = mean(endPrices) - fivePercent; disp(fivePercent);

