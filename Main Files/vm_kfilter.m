%=========================================================================
%                         Kalman Filter
% 
%      Measurement eq:    Yi(t)= X(t)*bi(t)+Ui(t)     
%      Transition eq:     bi(t)= hyp(8)*bi(t-1)+(1-hyp(8))*Bprior + mui(t)
%                                                 
%                        [ Ui(t)    ~   n( [0 , [SIGMA(i)     0                     
%                         mui(t) ]          0]       0    hyp(9)*Vprior] )
%=========================================================================

% hyperparameters

h1  = hyp(1);
h2  = hyp(2);
h3  = hyp(3);
h4  = hyp(4);
h5  = hyp(5);
h6  = hyp(6);
h7  = hyp(7);
h8  = hyp(8);
h9  = hyp(9);

%=========================================================================
% kalman filter
%=========================================================================
k  = nv*nlags_+1;

Tobs_=Tobs;

At_mat=zeros(Tobs_,k);
Pt_mat=zeros(k^2, Tobs_);
loglh = 0;
Res_=zeros(Tobs_,1);

% Initialization 
At = Bprior;
Pt = Vprior;

% Filter loop 
for t = 1:Tobs_

   At1 = At;
   Pt1 = Pt;

   % Forecasting   
   alphahat = h8 * At1+(1-h8)*Bprior;
   Phat     = h8 * Pt1 * h8' + Sprior;
   Phat     = 0.5*(Phat+Phat');
      
   yhat = X(t,:)*alphahat;
   nut  = Y(t) - yhat';
   Res_(t) = nut;
   Ft = X(t,:)*Phat*X(t,:)' + SIGMA_i;  
   Ft = 0.5*(Ft+Ft');
   
   
   % Log Likelihood
   loglh = loglh - 0.5*size(YY,2)*log(2*pi) - 0.5*log(det(Ft))...
                - 0.5*nut*inv(Ft)*nut';
   

   %** Updating **/
   At = alphahat + (Phat*X(t,:)')*inv(Ft)*nut';
   Pt = Phat - (Phat*X(t,:)')*inv(Ft)*(Phat*X(t,:)')';
   At_mat(t,:)  = At';
   Pt_mat(:,t)  = reshape(Pt,k^2,1);
end


loglh = real(loglh);

    
    
    
    
    
    
    
    
    
    
    
    
    
    