clear
close all
%%
%Five factor model compared to Fama-French's three factor model.
%Data is from China's stock market.

%Read the data from .xlsx or .CSV to MATLAB.
stockData=readtable("StockData.CSV",'VariableNamingRule','preserve');
ROE=readtable("ROE.xlsx",'VariableNamingRule','preserve');
ATGR=readtable("ATGR.xlsx",'VariableNamingRule','preserve');
marketReturn=readtable("MarketReturn.CSV",'VariableNamingRule','preserve');
rf=readtable("rf.xlsx",'VariableNamingRule','preserve');
stockData.date.Format="yyyy/MM/dd";
ROE.date.Format="yyyy/MM/dd";
ATGR.date.Format="yyyy/MM/dd";
marketReturn.date.Format="yyyy/MM/dd";
rf.date.Format="yyyy/MM/dd";



%Data processing, consolidation, and cleaning.
stockData.beme=1./stockData.beme;
stockData1=stockData(~(isnan(stockData.return)),1:end);
rf.rf=LinearFillNaN(rf.rf);
rf.rf=rf.rf/100;
marketReturn.rm=marketReturn.rm/100;

%Join all the tables into the table "data".
data=outerjoin(stockData1,ROE,'Keys',{'stockid','date'},'MergeKeys',true,'Type','left');
data=outerjoin(data,ATGR,'Keys',{'stockid','date'},'MergeKeys',true,'Type','left');
data=outerjoin(data,marketReturn,'Keys',{'date'},'MergeKeys',true,'Type','left');
data=outerjoin(data,rf,'Keys',{'date'},'MergeKeys',true,'Type','left');

writetable(data,"test.xlsx")
%-------------------------------------------------------------
%上面这段代码要跑很长时间，可以直接data=readtable("test.xlsx")
%(这个表格刚开始有点空，往下拉就是填满的，数据问题不影响）
%data的格式如下
%|stockid | date | return | me | beme | roe | tagr | rm | rf |
%__________________________________________________________
%stockid:股票代码。date:日期。return:月收益率。me:市值（元）。beme:BM值，账市比。roe:净资产收益率，用于计算RMW因子。
%tagr:总资产增长率，用于计算CMA因子
%
%最后分别进行三因子和五因子回归，比较各个因子的系数和R^2。








%%
% Design a swap according to spot yields, Vasicek and Monte Carlo
% Simulation. Simulate the floating rate and then decide the fixed rate.

%
clc; clear;

% Read data from 'China Treasury Spot Yields.xlsx'
% Data reading and data cleaning

dataFile = 'China Treasury Spot Yields.xlsx';
maturitiesInMonth = readmatrix(dataFile, Range='B1:O1');
maturitiesInYear = maturitiesInMonth / 12;

spotRates = readmatrix(dataFile, Range='B2:O206')/100;
spotRates=rmmissing(spotRates);

[numObs,numBonds] = size(spotRates);

%% 1 - Risk Natural Parameter Estimates

TreasSpotRate_1Month = spotRates(:,1);
FirstDifference = TreasSpotRate_1Month(2:end) - TreasSpotRate_1Month(1:end-1);
constant = ones(length(FirstDifference),1);

[b,bint,r,rint,stats] = regress(FirstDifference, [constant, TreasSpotRate_1Month(1:end-1)]);

% This leads to the following risk natrual parameters 
dt=1/12; 
alpha=b(1);
beta=b(2);
gamma_1=-beta/dt;
residSD=sqrt(stats(4));
sigma_1=residSD/sqrt(dt);
rbar_1=alpha/(gamma_1*dt);


fprintf('-----------------------------------------------------------------\n')
fprintf('Parameter estimates for Vasicek Model in risk natural environment\n');
fprintf('rbar equals to %5.5f\n', rbar_1);
fprintf('gamma equals to %5.5f\n', gamma_1);
fprintf('sigma equals to %5.5f\n', sigma_1);



%% 2 - Risk Neutral Parameter Estimates


rbar_star = 0.02;
gamma_star = 1.24;
sigma=0.1;
priceMatrix=zeros(numObs,numBonds);
for i=1:numBonds
    priceMatrix(:,i)=exp( - spotRates(:,i).* maturitiesInYear(i));
end

paravec = [rbar_star,gamma_star,sigma];

xdata = spotRates(:,1);

myfun = @(parameters, xdata) Pfunction(parameters,xdata,maturitiesInYear);

options = optimset('MaxIter',2000,'MaxFunEvals',2000);

[estimatedParas] = lsqcurvefit(myfun,paravec,xdata,priceMatrix,[0 0 0],[10 4 4],options);

rbar_star = estimatedParas(1);
gamma_star = estimatedParas(2);

fprintf('-----------------------------------------------------------------\n')
fprintf('Parameter estimates for Vasicek Model in risk neutral environment\n');
fprintf('rbar_star_hat equals to %5.5f\n', estimatedParas(1));
fprintf('gamma_star_hat equals to %5.5f\n', estimatedParas(2));
fprintf('sigma_hat equals to %5.5f\n', estimatedParas(3));

%% 3 - Design swaps by Monte Carlo simulation
% Assume that we are now designing a swap for two firms. 
% As we know the initial rate is 1.843%(today's interest rate), and then
% use Monte Carlo simulation to simulate floating rate path. The NPV of floating side should be same for fixed side. 
% So we can get the fixed rate of the swap

r0 = 1.8430 / 100; % instantaneous rate of floating rate path
FaceValue = 1e6;


times=(1:10)';
gamma=gamma_1;
rbar=rbar_1;
sigma=sigma_1;
paravec=[rbar,gamma,sigma];
FixedRate=zeros(length(times),1);
for i=1:length(times)
    t=times(i);
    [NPV,discountFactors]=floating(r0,FaceValue,paravec,t);
    FixedRate(i)=NPV/sum(0.5*FaceValue*discountFactors);
end

result=table(times,FixedRate,'VariableNames',{'swap years','fixed rate'});
fprintf('-----------------------------------------------------------------\n')
disp(result)






