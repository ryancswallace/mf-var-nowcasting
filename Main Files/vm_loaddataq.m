% DATASET: 10 Variables          
% UNRATE,HOURS,CPI,INDPRO,PCE,FEDFUND,TBOND,SP500: 1964:M1 - 2011:M3 
YM = csvread('YM.csv',1,1);
YM(:,1)   = YM(:,1)/100;
YM(:,6:7) = YM(:,6:7)/100;
YM(:,2:5) = log(YM(:,2:5));
YM(:,8)   = log(YM(:,8)); 
% USE BALANCED PANEL
YM        = YM(1:end-3,:);
% GDP,INVFIX: 1964:M1 - 2010:M12
YQ = csvread('YQ.csv',1,1);         
YQ = log(YQ);
YMM     = zeros(size(YM,1)/3,size(YM,2));
for i=1:size(YM,1)/3
    YMM(i,:)=mean(YM(3*(i-1)+1:3*i,:));
end

YY = [YMM YQ];
%==========================================================================
%                           DATA SPECIFICATION
%==========================================================================

nlags   = 2;              % number of lags   
T0      = nlags+1;        % size of pre-sample 
nex     = 1;              % number of exogenous vars; 1 means intercept only 
nv      = size(YY,2);     % number of variables 
nobs    = size(YY,1)-T0;  % number of observations 

spec    = [nlags T0 nex nv nobs];