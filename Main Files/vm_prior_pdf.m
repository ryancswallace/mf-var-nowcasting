%==========================================================================                                                          
%                PRIOR DISTRIBUTION FROM PRESAMPLE                                                                       
%==========================================================================

function [MN_logpdf, IW_logpdf]=vm_prior_pdf(hyp,YY,spec,PHI,SIG)

%==========================================================================
%                  Data Specification and Setting
%==========================================================================

nlags_  = spec(1);      % number of lags   
T0      = spec(2);      % size of pre-sample 
nex_    = spec(3);      % number of exogenous vars; 1 means intercept only 
nv      = spec(4);      % number of variables 
Nm      = spec(5);      % number of monthly variables
%========================================================================== 
%                         Dummy Observations                                    
%==========================================================================

% Obtain mean and standard deviation from expandend pre-sample data
 
YY0     =   YY(1:T0,:);  
ybar    =   mean(YY0)';      
sbar    =   std(YY0)'; 
premom  =   [ybar sbar];

% Generate matrices with dummy observations

[YYdum, XXdum] = varprior(nv,nlags_,nex_,hyp,premom);
n = size(YYdum,2);
inv_x     = inv(XXdum'*XXdum);
Phi_tilde = (inv_x)*XXdum'*YYdum;
Sigma     = (YYdum-XXdum*Phi_tilde)'*(YYdum-XXdum*Phi_tilde);

MN_pdf = mvnpdf(reshape(PHI,n*(n*nlags_+1),1),reshape(Phi_tilde,n*(n*nlags_+1),1),kron(SIG,inv_x));
MN_logpdf = log(MN_pdf);    

% inverse wishart pdf
[IW_pdf, IW_logpdf ] = invwishpdf(SIG,Sigma,length(YYdum)-nv*nlags_-1); 