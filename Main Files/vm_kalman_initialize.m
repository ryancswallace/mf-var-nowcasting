function [At Pt At_mat Pt_mat]=vm_kalman_initialize(Phi,sigma,YM,YQ,init_spec)

Tobs_ = init_spec(1);
kn    = init_spec(2);
p     = init_spec(3);
nv    = init_spec(4);
Nm    = init_spec(5);
Nq    = init_spec(6);
T0    = init_spec(7);

At_mat=zeros(Tobs_,kn);
Pt_mat=zeros(Tobs_,kn^2);
loglh = 0;
counter =0;

% STARTING VALUES
% Define Companion Form Matrix PHIF
PHIF         = zeros(nv*p,nv*p);          
IF           = eye(nv);

for i=1:p-1
    PHIF(i*nv+1:(i+1)*nv,(i-1)*nv+1:i*nv) = IF;
end

PHIF(1:nv,:) = Phi(1:end-1,:)';

% Define Constant Term CONF
CONF = [Phi(end,:)';zeros(nv*(p-1),1)];

% Define Covariance Term SIGF
SIGF = zeros(nv*p,nv*p);
SIGF(1:nv,1:nv) = sigma;

% Initialization 
At = zeros(kn,1);   
Pt = zeros(kn,kn);

for kk=1:5
    Pt = PHIF * Pt * PHIF' + SIGF;
end

% Measurement Equation
Z1         = zeros(Nm,kn);
Z1(:,1:Nm) = eye(Nm);

Z2         = zeros(Nq,kn);
for bb=1:Nq
    for ll=1:3
        Z2(bb,ll*Nm+(ll-1)*Nq+bb)=1/3;
    end
end

Z          = [Z1;Z2];


% Filter loop 
for t = 1:Tobs_
   
   if ((t+T0)/3)-floor((t+T0)/3)==0
       
       At1 = At;
       Pt1 = Pt;

       % Forecasting   
       alphahat = PHIF * At1 + CONF;
       Phat     = PHIF * Pt1 * PHIF' + SIGF;
       Phat     = 0.5*(Phat+Phat');
      
       yhat = Z*alphahat;
       nut  = [YM(t,:) YQ(t,:)]' - yhat;
       
       Ft = Z*Phat*Z';  
       Ft = 0.5*(Ft+Ft');  
 
       % Updating 
       At = alphahat + (Phat*Z')*inv(Ft)*nut;
       Pt = Phat - (Phat*Z')*inv(Ft)*(Phat*Z')';
       At_mat(t,:)  = At';
       Pt_mat(t,:)  = reshape(Pt,1,kn^2);
   
   else
       
       At1 = At;
       Pt1 = Pt;

       % Forecasting   
       alphahat = PHIF * At1 + CONF;
       Phat     = PHIF * Pt1 * PHIF' + SIGF;
       Phat     = 0.5*(Phat+Phat');
      
       yhat = Z1*alphahat;
       nut  = YM(t,:)' - yhat;
       
       Ft = Z1*Phat*Z1';  
       Ft = 0.5*(Ft+Ft');
      
  
       % Updating 
       At = alphahat + (Phat*Z1')*inv(Ft)*nut;
       Pt = Phat - (Phat*Z1')*inv(Ft)*(Phat*Z1')';
       At_mat(t,:)  = At';
       Pt_mat(t,:)  = reshape(Pt,1,kn^2);
       
   end
   
end