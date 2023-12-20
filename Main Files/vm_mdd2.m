%==========================================================================                                                          
%                 COMPUTE MARGINAL DATA DENSITY                                                                        
%==========================================================================

function [mdd,YYact,YYdum,XXact,XXdum]=vm_mdd2(hyp,YY,spec,efficient)

%==========================================================================
%                  Data Specification and Setting
%==========================================================================

nlags_  = spec(1);      % number of lags   
T0      = spec(2);      % size of pre-sample 
nex_    = spec(3);      % number of exogenous vars; 1 means intercept only 
nv      = spec(4);      % number of variables 
nobs    = spec(5);      % number of observations 

%========================================================================== 
%                         Dummy Observations                                    
%==========================================================================

% Obtain mean and standard deviation from expandend pre-sample data
 
YY0     =   YY(1:T0+24,:);  
ybar    =   mean(YY0)';      
sbar    =   std(YY0)'; 
premom  =   [ybar sbar];

% Generate matrices with dummy observations

[YYdum, XXdum] = varprior(nv,nlags_,nex_,hyp,premom);

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

% dummy:  YYdum XXdum 
% actual: YYact XXact 

YY=[YYdum' YYact']';
XX=[XXdum' XXact']';

n_total = size(YY,1);
n_dummy = n_total-nobs;
nv     = size(YY,2);
k      = size(XX,2); 

if efficient==0

    %==========================================================================
    %         Compute the log marginal data density for the VAR model
    %==========================================================================
    
    [Up, Sp, Vp] = svd(XXdum'*XXdum);
    sp = diag(Sp); kp = sum(sp>1e-12);
    inv_sv_XXdum = (Up(:,1:kp)*diag(1./sp(1:kp))*Vp(:,1:kp)')';
    det_XXdum    = prod(sp(1:kp));
    if det_XXdum==0
        det_XXdum = 1e-300;
    end
    
    [Up, Sp, Vp] = svd(XX'*XX);
    sp = diag(Sp); kp = sum(sp>1e-12);
    inv_sv_XX    = (Up(:,1:kp)*diag(1./sp(1:kp))*Vp(:,1:kp)')';
    det_XX       = prod(sp(1:kp));
    
    
    
    S0     = (YYdum'*YYdum)-((YYdum'*XXdum)*(inv_sv_XXdum))*XXdum'*YYdum;
    
    S1     = (YY'*YY)-((YY'*XX)*(inv_sv_XX))*XX'*YY;
    
    % compute constants for integrals
    
    i=1;
    gam0=0;
    gam1=0;
    
    while i <= nv;
        %gam0=gam0+log(gamma(0.5*(n_dummy-k+1-i)));
        gam0=gam0+gammaln(0.5*(n_dummy-k+1-i));
        %gam1=gam1+log(gamma(0.5*(n_total-k+1-i)));
        gam1=gam1+gammaln(0.5*(n_total-k+1-i));
        i=i+1;
    end;
    
    % dummy observation
    
    lnpY0 = -nv*(n_dummy-k)*0.5*log(pi)-(nv/2)*log(abs(det_XXdum))...
        -(n_dummy-k)*0.5*log(abs(det(S0)))+nv*(nv-1)*0.25*log(pi)+gam0;
    
    % dummy and actual observations
    
    lnpY1 = -nv*(n_total-k)*0.5*log(pi)-(nv/2)*log(abs(det_XX))...
        -(n_total-k)*0.5*log(abs(det(S1)))+nv*(nv-1)*0.25*log(pi)+gam1;
    
    lnpYY = lnpY1-lnpY0;
    
    % marginal data density
    mdd = lnpYY;

else
    
    mdd = -1e+10;
    
end
