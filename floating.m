function [value,discountFactors] = floating(r0,FaceValue,paravec,simulationTime)
rbar_star = paravec(1);
gamma_star = paravec(2);
sigma = paravec(3);

rng('default');


timeStepPerPeriod = 252; % 每年252个时间步

totalNumSteps = timeStepPerPeriod * simulationTime; % 总时间步
numSimulation = 10000; % 模拟次数
deltat = 1 / timeStepPerPeriod; % 时间步长

% 初始化变量
rpath = zeros(1, totalNumSteps); % 利率路径
rt = zeros(numSimulation, simulationTime * 2); % 每半年支付的浮动利率
dt = zeros(numSimulation, simulationTime * 2); % 每半年支付的折现因子

discountFactor = zeros(1, totalNumSteps); % 折现因子
cashFlow = zeros(numSimulation, simulationTime * 2); % 每半年支付的现金流

allPrices = zeros(1, numSimulation); % 每次模拟的价格

% Monte Carlo 模拟
for i = 1:numSimulation
    RandomVars = randn(1, totalNumSteps);
    oldr = r0;
    olddiscountFactor = 1;
    for j = 1:totalNumSteps
        % 利率路径更新公式
        rpath(j) = oldr + gamma_star * (rbar_star - oldr) * deltat + sigma * sqrt(deltat) * RandomVars(j);       
        % 折现因子更新公式
        discountFactor(j) = olddiscountFactor * exp(-oldr * deltat);
        olddiscountFactor = discountFactor(j);
        oldr = rpath(j);
    end
    % 提取半年支付点的浮动利率和折现因子
    rt(i, :) = rpath((timeStepPerPeriod / 2) * [1:simulationTime * 2]);
    dt(i, :) = discountFactor((timeStepPerPeriod / 2) * [1:simulationTime* 2]);
    % 半年支付的现金流计算
    cashFlow(i, :) = FaceValue * 0.5 *rt(i, :);
    % 模拟的IRS价格
    allPrices(i) = cashFlow(i, :) * dt(i, :)';
end

% 计算平均价格
value = mean(allPrices);
discountFactors=mean(dt,1);
end
