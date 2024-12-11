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
%Pricing swaps according to spot yields, Vasicek and CIR model and
%Monte Carlo Simulation.

%数据是China Treasury Spot Yields.xlsx，纵轴是时间，横轴是期限
%分别用Vasicek和CIR模型，在风险中性下用lsqcurvefit函数得到参数
%今天即期利率为1.881%，为不同到期期限的利率互换定价。
%最终输出一个表格，第一列是不同的到期期限（1个月、3个月、6个月、1年、2年、5年），第二列是利率互换对应的固定利率









