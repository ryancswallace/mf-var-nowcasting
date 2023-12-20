[YM0, YMX]  = xlsread('YM_May2020.xls');
[YQ0, YQX]  = xlsread('YQ_May2020.xls');

mdate       = YMX(2:end,1);
select      = [1 0 0 0 0 1 1 0];
YM0(:,select==1) = YM0(:,select==1)./100;
YM0(:,select==0) = log(YM0(:,select==0));
YQ0         = log(YQ0(:,:));
YM          = YM0;
YQ          = kron(YQ0,ones(3,1));

% allowing for possible release date mismatches within QF variables
% Tstar is the length of all available data (at least includes a non-NaN variable)
% T is the length of the balanced panel set
Tstar       = size(YM,1); 
Tdata       = find((sum(isnan(YQ),2))==0, 1, 'last' );    
Torigin     = Tdata;

% collect data
YDATA       = nan(Tstar,11);                      % 11 variables
YDATA(:,1:8)= YM; 
YDATA(1:length(YQ),9:end) = YQ;

YDATA_S = YDATA;
YM_S    = YM;
YQ_S    = YQ;

% set dimension, mark missing data
Nm          = size(YM,2);
Nq          = size(YQ,2);
index_NY    = isnan(YDATA(Tdata+1:Tstar,:))';

% drop a few covid samples: Can only do after dselect>=7
T           = find(datenum(YMX(2:end,1))==datenum('3/1/2020'));
% Tdrop1      = find(datenum(YMX(2:end,1))==datenum('4/1/2020'));
% Tdrop2      = find(datenum(YMX(2:end,1))==datenum('7/1/2020'));

index_NY    = [zeros(Nm+Nq,Tdata-T) index_NY];
%index_NY(:,Tdrop1-T:Tdrop2-T) = 1;
%YDATA(Tdrop1:Tdrop2,:)    = nan;
%YM(Tdrop1:Tdrop2,:)       = nan;
%YQ(Tdrop1:Tdrop2,:)       = nan;


% specification
nlags_  = 6;                  % number of lags   
T0      = 2*nlags_;           % size of pre-sample
nex     = 1;                  % number of exogenous vars;1(=intercept only) 
p       = nlags_;
nlags   = p;
nv      = Nm+Nq;
kq      = Nq*p;
nobs    = T-T0;
Tnew    = Tstar-T;
Tnobs   = Tstar-T0;
   