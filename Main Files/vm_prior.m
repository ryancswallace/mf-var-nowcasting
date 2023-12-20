%=========================================================================
%         Implementation of Equation-by-Equation Minnesota Prior  
%=========================================================================

clear Bprior 

h1  = hyp(1);
h2  = hyp(2);
h3  = hyp(3);
h4  = hyp(4);
h5  = hyp(5);
h6  = hyp(6);
h7  = hyp(7);
h8  = hyp(8);
h9  = hyp(9);

% discrete time representation of continuous time random walk
alpha  = sqrt(3)-2;

% prior mean
Bprior = zeros(k,1);

for ii=1:nlags_
    Bprior(rr+nv*(ii-1))=(1-alpha)*alpha^(ii-1);
end


% prior variance
Vprior_p = zeros(k-1,1);

for l=1:nlags_
    for j=1:nv        
        if j==rr
            Vprior_p(j+(l-1)*nv)=(h4*h1)/l;
        else
            Vprior_p(j+(l-1)*nv)=(h4*h2*SIGMA_i)/(l*SIGMA(j));
        end        
    end
end

Vprior = diag([Vprior_p;(h4*h3*SIGMA_i)]);

% Prior mean and variance for the FEDFUND is zero
if rr==10           
        Bprior      = zeros(k,1);
        Bprior(end) = log(0.2);
        Vprior      = zeros(k,k);
end

% variance of changes in the parameter vector
Sprior          = h9*Vprior;

% restriction on sum of coefficient
sig_i  = sqrt(SIGMA_i);
sumvec = zeros(k,1);

for kk=1:nv
    sumvec = zeros(k,1);
    for hh=1:nlags_
        sumvec(kk+(hh-1)*nv)=(h5*sqrt(SIGMA(kk)))/sig_i;
    end
    Vprior = Vprior - ((Vprior*sumvec)*(Vprior*sumvec)')/(1+sumvec'*(Vprior*sumvec));
end


%=========================================================================
%                        Dummy Observations
%     
%   Purpose: 1) introduce appropriate off-diagonal components
%            2) get price-level neutrality 
%=========================================================================


if additional1==1
% forecast over dummy observation using kalman filter

% Initialization 
At = Bprior;
Pt = Vprior;

% model's own forecast: kalman filter finds 0 error and sets posterior to
% prior

for kk=1:size(Xdm,1)

XX = h6*Xdm(kk,:);
YY = h6*Xdm(kk,:)*Bprior;

At1 = At;
Pt1 = Pt;

% Forecasting   
alphahat = h8 * At1+(1-h8)*Bprior;
Phat     = h8 * Pt1 * h8' + h9*Vprior;
Phat     = 0.5*(Phat+Phat');
      
yhat = XX*alphahat;
nut  = YY - yhat';   
Ft = XX*Phat*XX' + SIGMA_i;  
Ft = 0.5*(Ft+Ft');

% Updating 
At = alphahat + (Phat*XX')*inv(Ft)*nut';
Pt = Phat - (Phat*XX')*inv(Ft)*(Phat*XX')';

end

Bprior = At;
Vprior = Pt;

end

if additional2==1
% additional weight on dummy observations of nominal series to get
% neutrality
nomdata  = h7*[0 0 0 0 0 0 0 1 1 1 1 1 0 1 0 0];
nomdatax = [repmat(nomdata,1,nlags_) 0];
Ydm = nomdata(rr);                     % 1-by-1
Xdm = nomdatax;                        % 1-by-k

% Initialization 
At = Bprior;
Pt = Vprior;

for kk=1:nlags_

XX = Xdm;
YY = Ydm;

At1 = At;
Pt1 = Pt;

% Forecasting   
alphahat = h8 * At1+(1-h8)*Bprior;
Phat     = h8 * Pt1 * h8' + h9*Vprior;
Phat     = 0.5*(Phat+Phat');
      
yhat = XX*alphahat;
nut  = YY - yhat';   
Ft = XX*Phat*XX' + SIGMA_i;  
Ft = 0.5*(Ft+Ft');

% Updating 
At = alphahat + (Phat*XX')*inv(Ft)*nut';
Pt = Phat - (Phat*XX')*inv(Ft)*(Phat*XX')';

end

Bprior = At;
Vprior = Pt;

end












