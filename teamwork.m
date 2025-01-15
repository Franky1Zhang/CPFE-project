clear
close all


%% Data reading and cleaning
data=readtable("test.xlsx","VariableNamingRule","preserve","ReadVariableNames",true);
data.return=data.return/100;
data=rmmissing(data);
data.date=12*year(data.date)+month(data.date);


% Creat Lagged me
data1=data(:,{'date','stockid','me'});
data1.wt=data1.me;
data1.me=[];
data1.date=data1.date+1;

data=outerjoin(data,data1,'Keys',{'date','stockid'},'MergeKeys',true,'Type','left');
data=rmmissing(data);

%% Creat groups of portfolios and 5 factors

% Creat labels
[G,date]=findgroups(data.date);

sizemedn =splitapply(@median,data.me,G);
breaks=table(date,sizemedn);

prctile_30=@(input)prctile(input,30);
prctile_70=@(input)prctile(input,70);

breaks.roe30 = splitapply(prctile_30,data.roe,G);
breaks.roe70 = splitapply(prctile_70,data.roe,G);

breaks.bm30=splitapply(prctile_30,data.beme,G);
breaks.bm70=splitapply(prctile_70,data.beme,G);

breaks.me30=splitapply(prctile_30,data.me,G);
breaks.me70=splitapply(prctile_70,data.me,G);

breaks.tagr30=splitapply(prctile_30,data.tagr,G);
breaks.tagr70=splitapply(prctile_70,data.tagr,G);

data_breakes=innerjoin(data,breaks,'Keys',{'date'});

% Assign size portfolios
szport=rowfun(@sz_bucket,data_breakes(:,{'me','sizemedn'}),'OutputFormat','cell');
data_breakes.szport=cell2mat(szport);

% Assign book-to-market portfolio
bmport=rowfun(@bm_bucket,data_breakes(:,{'beme','bm30','bm70'}),'OutputFormat','cell');
data_breakes.bmport=cell2mat(bmport);

% Assign RMW portfolio
roeport=rowfun(@roe_bucket,data_breakes(:,{'roe','roe30','roe70'}),'OutputFormat','cell');
data_breakes.roeport=cell2mat(roeport);

% Assign CMA portfolio
tagrport=rowfun(@tagr_bucket,data_breakes(:,{'tagr','tagr30','tagr70'}),'OutputFormat','cell');
data_breakes.tagrport=cell2mat(tagrport);



port=data_breakes(:,{'stockid','date','return','bmport','szport','roeport','tagrport','wt'});


% Create 54 groups of portfolios

[G,date,bmport,szport,roeport,tagrport]=findgroups(port.date,port.bmport,port.szport,port.roeport,port.tagrport);

vwret=splitapply(@wavg,port(:,{'return','wt'}),G);

% Create labels for the 54 portfolios

portname=strcat(szport,bmport,roeport,tagrport);
portname=cellstr(portname);
vwret_table=table(vwret,date,szport,bmport,roeport,tagrport,portname);

% Reshape the dataset into a a table for portfolio return organized with
% date in ascending order. This faciliates the creation of factors

ff_factors=unstack(vwret_table(:,{'vwret','date','portname'}),'vwret','portname');

SMB=['S','B'];
HML=['H','M','L'];
RMW=['R','M','W'];
CMA=['C','M','A'];
ff_factors = fillmissing(ff_factors, 'constant', 0);

% Creat MKT factor
data.mkt=data.rm-data.rf;
[G,date]=findgroups(data.date);
ff_factors.MKT=splitapply(@mean,data(:,{'mkt'}),G);

% Creat SMB factor
S={};
B={};
for i=1:3
    for j=1:3
        for k=1:3
            S=[S,strcat('S',HML(i),RMW(j),CMA(k))];
            B=[B,strcat('B',HML(i),RMW(j),CMA(k))];
        end
    end
end
T1=ff_factors(:,S);
T2=ff_factors(:,B);
ff_factors.SMB = mean(T1{:,:}, 2)-mean(T2{:,:}, 2);

% Creat HML factor
H={};
L={};
for i=1:2
    for j=1:3
        for k=1:3
            H=[H,strcat(SMB(i),'H',RMW(j),CMA(k))];
            L=[L,strcat(SMB(i),'L',RMW(j),CMA(k))];
        end
    end
end
T1=ff_factors(:,H);
T2=ff_factors(:,L);
ff_factors.HML = mean(T1{:,:}, 2)-mean(T2{:,:}, 2);

% Creat RMW factor
R={};
W={};
for i=1:2
    for j=1:3
        for k=1:3
            R=[R,strcat(SMB(i),HML(j),'R',CMA(k))];
            W=[W,strcat(SMB(i),HML(j),'W',CMA(k))];
        end
    end
end
T1=ff_factors(:,R);
T2=ff_factors(:,W);
ff_factors.RMW = mean(T1{:,:}, 2)-mean(T2{:,:}, 2);

% Creat CMA factor
C={};
A={};
for i=1:2
    for j=1:3
        for k=1:3
            C=[C,strcat(SMB(i),HML(j),RMW(k),'C')];
            A=[A,strcat(SMB(i),HML(j),RMW(k),'A')];
        end
    end
end
T1=ff_factors(:,C);
T2=ff_factors(:,A);
ff_factors.CMA = mean(T1{:,:}, 2)-mean(T2{:,:}, 2);

%% Fama-French 3 factors and 5 factors regression
numGroups=54;
ff5factor = table2array(ff_factors(:, { 'MKT', 'SMB', 'HML','RMW','CMA'}));
groupsReturn=table2array(ff_factors(:,2:numGroups+1));
constant = ones(size(ff5factor,1),1); 

% FF3 regression
ff3_gamma1=zeros(numGroups,1);
ff3_gamma2=zeros(numGroups,1);
ff3_gamma3=zeros(numGroups,1);
ff3_alpha = zeros(numGroups,1);
ff3_r_square = zeros(numGroups, 1);
for i = 1:numGroups
  [b,bint,r,rint,stats] = ...
      regress(groupsReturn(:,i), [constant ff5factor(:,1:3)]); 

   ff3_gamma1(i)=b(2);
   ff3_gamma2(i)=b(3);
   ff3_gamma3(i)=b(4);
   ff3_r_square(i) = stats(1); 
end
result1=table((1:numGroups)',ff3_gamma1,ff3_gamma2,ff3_gamma3,ff3_r_square,'VariableNames',{'Group','MKT','SMB','HML','R_square'});
disp(result1)

fprintf(1, ...
        '--------------------------------------------\n'); 
fprintf(1, ...
        'FF3 average R_square is %4.4f\n',mean(ff3_r_square)); 


% FF5 regression
ff5_gamma1=zeros(numGroups,1);
ff5_gamma2=zeros(numGroups,1);
ff5_gamma3=zeros(numGroups,1);
ff5_gamma4=zeros(numGroups,1);
ff5_gamma5=zeros(numGroups,1);
ff5_alpha = zeros(numGroups,1);
ff5_r_square = zeros(numGroups, 1);
for i = 1:numGroups
  [b,bint,r,rint,stats] = ...
      regress(groupsReturn(:,i), [constant ff5factor(:,1:5)]); 

   ff5_gamma1(i)=b(2);
   ff5_gamma2(i)=b(3);
   ff5_gamma3(i)=b(4);
   ff5_gamma4(i)=b(5);
   ff5_gamma5(i)=b(6);
   ff5_r_square(i) = stats(1); 
end

result2=table((1:numGroups)',ff5_gamma1,ff5_gamma2,ff5_gamma3,ff5_gamma4,ff5_gamma5,ff5_r_square,'VariableNames',{'Group','MKT','SMB','HML','RMW','CMA','R_square'});
disp(result2)

fprintf(1, ...
        '--------------------------------------------\n'); 
fprintf(1, ...
        'FF5 average R_square is %4.4f\n',mean(ff5_r_square)); 





