function [At_final,Pt_final] = vm_initialize(GAMMAs,GAMMAz,GAMMAc,GAMMAu,...
    LAMBDAs,LAMBDAz,LAMBDAc,LAMBDAu,LAMBDAs_t,LAMBDAz_t,LAMBDAc_t,LAMBDAu_t,...
    sig_qq,sig_mm,sig_qm,sig_mq,Zm,YDATA,init_mean,init_var,spec)

% Specification
p       = spec(1);      % number of lags   
T0      = spec(2);      % size of pre-sample 
nex_    = spec(3);      % number of exogenous vars; 1 means intercept only 
nv      = spec(4);      % number of variables 
Nm      = spec(5);      % number of monthly variables

% Initialization         
At   = init_mean;  
Pt   = init_var;
    
    
% Kalman Filter loop 
for t = p+2:T0            
       
   if ((t)/3)-floor((t)/3)==0
       
       At1 = At;
       Pt1 = Pt;

       % Forecasting   
       alphahat = GAMMAs * At1 + GAMMAz * Zm(t-p-1,:)' + GAMMAc;
       Phat     = GAMMAs * Pt1 * GAMMAs' + GAMMAu * sig_qq * GAMMAu';
       Phat     = 0.5*(Phat+Phat');
       
       yhat     = LAMBDAs * alphahat + LAMBDAz * Zm(t-p-1,:)' + LAMBDAc;
       nut      = [YDATA(t,:)'] - yhat;
       
       Ft       = LAMBDAs * Phat * LAMBDAs' + LAMBDAu * sig_mm * LAMBDAu'...
                  + LAMBDAs*GAMMAu*sig_qm*LAMBDAu'...
                  + LAMBDAu*sig_mq*GAMMAu'*LAMBDAs';
       Ft       = 0.5*(Ft+Ft');
       Xit      = LAMBDAs * Phat + LAMBDAu * sig_mq * GAMMAu';       
       
       At       = alphahat + Xit' * inv(Ft) * nut;
       Pt       = Phat     - Xit' * inv(Ft) * Xit;             
          
   else
       
       At1 = At;
       Pt1 = Pt;
       
       % Forecasting   
       alphahat = GAMMAs * At1 + GAMMAz * Zm(t-p-1,:)' + GAMMAc;
       Phat     = GAMMAs * Pt1 * GAMMAs' + GAMMAu * sig_qq * GAMMAu';
       Phat     = 0.5*(Phat+Phat');
       
       yhat     = LAMBDAs_t * alphahat + LAMBDAz_t * Zm(t-p-1,:)' + LAMBDAc_t;
       nut      = YDATA(t,1:8)' - yhat;
       
       Ft       = LAMBDAs_t * Phat * LAMBDAs_t' + LAMBDAu_t * sig_mm * LAMBDAu_t'...
                  + LAMBDAs_t*GAMMAu*sig_qm*LAMBDAu_t'...
                  + LAMBDAu_t*sig_mq*GAMMAu'*LAMBDAs_t';
       Ft       = 0.5*(Ft+Ft');
       Xit      = LAMBDAs_t * Phat + LAMBDAu_t * sig_mq * GAMMAu';
       
       At       = alphahat + Xit' * inv(Ft) * nut;
       Pt       = Phat     - Xit' * inv(Ft) * Xit; 
       
   end
   
end

At_final = At;
Pt_final = Pt;

end
