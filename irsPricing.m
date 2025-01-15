function IRSPrice = irsPricing(r0, fixedRate, FaceVal, rbar_star, gamma_star, sigma)
%IRSPrice is to calculate the price of an interest rate swap
%   r0 is the initial rate of floating rate base, which is Treasury Spot
%   Rate of 1 month at the beginning. Fixed rate is the fixed rate you need
%   to pay in an interest rate swap.

rng('default');

IRS_maturityInYear = 5; % IRS到期时间5年

timeStepPerPeriod = 252; % 每年252个时间步
simulationHorizon = IRS_maturityInYear; % 模拟总时间5年
totalNumSteps = timeStepPerPeriod * simulationHorizon; % 总步数
numSimulation = 10000; % 模拟次数
deltat = 1 / timeStepPerPeriod; % 时间步长

% 存储每次模拟的浮动利率
rpath = zeros(1, totalNumSteps);
% 存储每半年支付的浮动利率
rt = zeros(numSimulation, IRS_maturityInYear * 2);

% 存储每次模拟的折现因子
discountFactor = zeros(1, totalNumSteps);
% 存储每半年支付的折现因子
dt = zeros(numSimulation, IRS_maturityInYear * 2);

% 存储每半年支付的现金流
cashFlow = zeros(numSimulation, IRS_maturityInYear * 2);

% 存储每次模拟的IRS价格
allPrices = zeros(numSimulation, 1);

% Monte Carlo模拟
for i = 1:numSimulation
    randVars = randn(totalNumSteps, 1);
    for j = 1:totalNumSteps
        if j == 1
            dr = gamma_star * (rbar_star - r0) * deltat + sigma * sqrt(deltat) * randVars(j);
            rpath(j) = r0 + dr;
            discountFactor(j) = exp(-r0 * deltat);
        else
            dr = gamma_star * (rbar_star - rpath(j-1)) * deltat + sigma * sqrt(deltat) * randVars(j);
            rpath(j) = rpath(j-1) + dr;
            discountFactor(j) = discountFactor(j-1) * exp(-rpath(j-1) * deltat);
        end
    end
    % 提取半年支付点的浮动利率和折现因子
    rt(i, :) = rpath((totalNumSteps / (IRS_maturityInYear * 2)) * [1:IRS_maturityInYear * 2]);
    dt(i, :) = discountFactor((totalNumSteps / (IRS_maturityInYear * 2)) * [1:IRS_maturityInYear * 2]);
    % 半年现金流计算
    cashFlow(i, :) = FaceVal * 0.5 * (fixedRate - 0.0001 * 0.8750 - rt(i, :));
    % 计算模拟价格
    allPrices(i) = cashFlow(i, :) * dt(i, :)';
end

% IRS价格取所有模拟价格的平均值
IRSPrice = mean(allPrices);

return
