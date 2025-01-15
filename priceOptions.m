%%这是PS3作业bonus的延伸思考，当没有解析解的时候，如果要给这个期权定价应该如何做？有点类似于VaR
%%使用的是CIR模型


rng('default');
BondMaturity = 5;
OptionMaturity = 3;
LeftMaturity = 2;
numSim = 1000;
timeStepsPerPeriod = 252;

deltat =1/timeStepsPerPeriod;
sqrt_deltat = sqrt(deltat);
timeStepsPerSemi = timeStepsPerPeriod/2;
totalOptionSteps = timeStepsPerPeriod *OptionMaturity;
totalLeftSteps = timeStepsPerPeriod *LeftMaturity;
FaceVal = 100;

initialVal = 0.04;
sigma= 1.2373;
gamma_star= 1.24;
r_star = 0.04;
APR = 5/100;
EndRate=zeros(1,numSim);
discountFactor = zeros(1,totalLeftSteps);
rate_Left = zeros(1,totalLeftSteps);
discountFactor_semi = zeros(1,2*LeftMaturity);
rate_semi = zeros(1,2*LeftMaturity);
cashflow_semi = zeros(1,2*LeftMaturity);
OptionPrice = zeros(1,numSim);
discount = zeros(1,numSim);
EndPrices = zeros(numSim,numSim);
EndPrice = zeros(1,numSim);
%%先对于期权未到期时候进行模拟，只需记录时间末点的折现因子和利率，折现因子用来最后算期权的折现价格，利率是下一步继续模拟的起点
for i = 1:numSim
  discountfactor = 1;
  oldr = initialVal;
  randomVars = randn(1,totalOptionSteps);
  for thisStep = 1:totalOptionSteps
      discountfactor = discountfactor*exp(-deltat*oldr);
      thisShock =randomVars(thisStep);
      driftTerm = gamma_star*(r_star-oldr);
      diffusionTerm = sqrt(sigma*max(oldr,0));
      newr = oldr+driftTerm*deltat+diffusionTerm*thisShock*sqrt_deltat;
      oldr = newr;
  end
  EndRate(i) = oldr;
  discount(i) = discountfactor;
end
for i = 1:4
    if i == 4
        cashflow_semi(i) = FaceVal* APR/2+FaceVal;
    else
        cashflow_semi(i) = FaceVal* APR/2;
    end
end 
cashflow = repmat(cashflow_semi,numSim,1);
for i = 1:numSim
    for j = 1:numSim
        OlddiscountFactor = 1;
        oldr = EndRate(i);
        randomVars_2 = randn(1,totalLeftSteps);
        for thisStep = 1:totalLeftSteps
            NewdiscountFactor = OlddiscountFactor*exp(-deltat*oldr);
            OlddiscountFactor = NewdiscountFactor;
            discountFactor (thisStep) = NewdiscountFactor;
            thisShock =randomVars_2(thisStep);
            driftTerm = gamma_star*(r_star-oldr);
            diffusionTerm = sqrt(sigma*max(oldr,0));
            newr = oldr+driftTerm*deltat+diffusionTerm*thisShock*sqrt_deltat;
            oldr = newr;
            rate_Left(thisStep) = newr;
        end
    discountFactor_semi (j,:) = discountFactor(timeStepsPerSemi *[1:2*LeftMaturity]);
    rate_semi(j,:) = rate_Left(timeStepsPerSemi *[1:2*LeftMaturity]);
    EndPrices(i,j) = cashflow(j,:)*discountFactor_semi(j,:)';
    end
    EndPrice(i) = mean(EndPrices(i,:));
end
%%现在我们得到了option到期时，coupon的价格，可以根据这个价格来决定是否行权
for i = 1:numSim
    if EndPrice(i)>FaceVal
        OptionPayoff = max(0,EndPrice(i) - FaceVal);
    else
        OptionPayoff = 0;
    end
    OptionPrice(i) = OptionPayoff*discount(i);
end

FinalOptionPrice = mean(OptionPrice);
disp(FinalOptionPrice);
