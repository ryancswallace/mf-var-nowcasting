if dselect==1
    [YM0, YMX]  = xlsread('YM_Jan2020.xls');
    [YQ0, YQX]  = xlsread('YQ_Jan2020.xls');
elseif dselect==3
    [YM0, YMX]  = xlsread('YM_March2020.xls');
    [YQ0, YQX]  = xlsread('YQ_March2020.xls');
elseif dselect==4
    [YM0, YMX]  = xlsread('YM_April2020.xls');
    [YQ0, YQX]  = xlsread('YQ_April2020.xls');
elseif dselect==5
    [YM0, YMX]  = xlsread('YM_May2020.xls');
    [YQ0, YQX]  = xlsread('YQ_May2020.xls');
elseif dselect==6
    [YM0, YMX]  = xlsread('YM_June2020.xls');
    [YQ0, YQX]  = xlsread('YQ_June2020.xls');
elseif dselect==7
    [YM0, YMX]  = xlsread('YM_July2020.xls');
    [YQ0, YQX]  = xlsread('YQ_July2020.xls');
elseif dselect==8
    [YM0, YMX]  = xlsread('YM_Aug2020.xls');
    [YQ0, YQX]  = xlsread('YQ_Aug2020.xls');
elseif dselect==9
    [YM0, YMX]  = xlsread('YM_Sept2020.xls');
    [YQ0, YQX]  = xlsread('YQ_Sept2020.xls');
elseif dselect==10
    [YM0, YMX]  = xlsread('YM_Oct2020.xls');
    [YQ0, YQX]  = xlsread('YQ_Oct2020.xls');
elseif dselect==11
    [YM0, YMX]  = xlsread('YM_Nov2020.xls');
    [YQ0, YQX]  = xlsread('YQ_Nov2020.xls');
elseif dselect==12
    [YM0, YMX]  = xlsread('YM_Dec2020.xls');
    [YQ0, YQX]  = xlsread('YQ_Dec2020.xls');
elseif dselect==13
    [YM0, YMX]  = xlsread('YM_Jan2021.xls');
    [YQ0, YQX]  = xlsread('YQ_Jan2021.xls');
elseif dselect==14
    [YM0, YMX]  = xlsread('YM_Feb2021.xls');
    [YQ0, YQX]  = xlsread('YQ_Feb2021.xls');
elseif dselect==15
    [YM0, YMX]  = xlsread('YM_Mar2021.xls');
    [YQ0, YQX]  = xlsread('YQ_Mar2021.xls');
elseif dselect==16
    [YM0, YMX]  = xlsread('YM_Apr2021.xls');
    [YQ0, YQX]  = xlsread('YQ_Apr2021.xls');
elseif dselect==17
    [YM0, YMX]  = xlsread('YM_May2021.xls');
    [YQ0, YQX]  = xlsread('YQ_May2021.xls');
elseif dselect==18
    [YM0, YMX]  = xlsread('YM_June2021.xls');
    [YQ0, YQX]  = xlsread('YQ_June2021.xls');
elseif dselect==19
    [YM0, YMX]  = xlsread('YM_July2021.xls');
    [YQ0, YQX]  = xlsread('YQ_July2021.xls');
elseif dselect==20
    [YM0, YMX]  = xlsread('YM_Aug2021.xls');
    [YQ0, YQX]  = xlsread('YQ_Aug2021.xls');
end

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

% set dimension, mark missing data
Nm          = size(YM,2);
Nq          = size(YQ,2);
index_NY    = isnan(YDATA(Tdata+1:Tstar,:))';

% drop a few covid samples: Can only do after dselect>=7
if dselect>=3
    T           = find(datenum(YMX(2:end,1))==datenum('12/1/2019'));
    index_NY    = [zeros(Nm+Nq,Tdata-T) index_NY];
else
    T = Tdata;
end

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
   



