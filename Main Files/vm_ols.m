%==========================================================================
%                                OLS                                     
%  Notes :  Generates variance of the residual term and actual observations
%  Output:  YYact, XXact, sbar
%==========================================================================

%==========================================================================
%                  Data Specification and Setting
%==========================================================================

nlags_  = spec(1);      % number of lags   
T0      = spec(2);      % size of pre-sample 
nex_    = spec(3);      % number of exogenous vars; 1 means intercept only 
nv      = spec(4);      % number of variables 
nobs    = spec(5);      % number of observations 

k       = nv*nlags_+1;

%==========================================================================
%                           Observation
%==========================================================================

% Actual observations

YYact = YY(T0+1:T0+nobs,:);
XXact = zeros(nobs,nv*nlags_);
i = 1;

while (i <= nlags_)
    XXact(:,(i-1)*nv+1:i*nv) = YY(T0-(i-1):T0+nobs-i,:);
    i = i+1;
end

% last column of XXact = constant
XXact = [XXact ones(nobs,1)];

% Obtain standard deviation from a 6-lag univariate autoregression
sbar  = zeros(nv,1);

for i=1:nv
    ydata  = YYact(:,i);
    xdata  = XXact(:,i);
    bols   = inv(xdata'*xdata)*(xdata'*ydata);
    sbar(i)= std(ydata-xdata*bols);
end





