function value = IRS_price_zrc(r0, fixrate, FV, paravec)
rbar_star = paravec(1);
gamma_star = paravec(2);
sigma = paravec(3);

rng('default');
IRS_maturityInYear = 5; % 到期时间5年

timeStepPerPeriod = 252; % 每年252个时间步
simulationHorizon = 5; % 总模拟时间为5年
totalNumSteps = timeStepPerPeriod * simulationHorizon; % 总时间步
numSimulation = 10000; % 模拟次数
deltat = 1 / timeStepPerPeriod; % 时间步长

% 初始化变量
rpath = zeros(1, totalNumSteps); % 利率路径
rt = zeros(numSimulation, IRS_maturityInYear * 2); % 每半年支付的浮动利率
dt = zeros(numSimulation, IRS_maturityInYear * 2); % 每半年支付的折现因子

discountFactor = zeros(1, totalNumSteps); % 折现因子
cashFlow = zeros(numSimulation, IRS_maturityInYear * 2); % 每半年支付的现金流

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
    rt(i, :) = rpath((timeStepPerPeriod / 2) * [1:IRS_maturityInYear * 2]);
    dt(i, :) = discountFactor((timeStepPerPeriod / 2) * [1:IRS_maturityInYear * 2]);
    % 半年支付的现金流计算
    cashFlow(i, :) = FV * 0.5 * (fixrate - 0.0001 * 0.8750 - rt(i, :));
    % 模拟的IRS价格
    allPrices(i) = cashFlow(i, :) * dt(i, :)';
end

% 计算平均价格
value = mean(allPrices);

end
