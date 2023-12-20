%======================================================================
%
%        MF BVAR : MIXED FREQ BAYESIAN VAR FORECASTING                     
%
%        Date : September 2023
%
%        Note : Given estimated MF-VAR coefficients, this code generates  
%               forecasts     
%
%        Code Written By:  Frank Schorfheide  schorf@ssc.upenn.edu
%                          Dongho Song        dongho.song@jhu.edu
%======================================================================

%======================================================================
%                          HOUSEKEEPING
%======================================================================
close all
clc;clear

ss = path;
path('Data',path);
path('Main Files',path);
path('Moment',path);
path('Gauss Files',path);
path('Output',path);


H        = 63;         % forecast horizon
nsim     = 20000;      % number of draws from Posterior Density
nburn    = 0.2*nsim;   % number of draws to discard
savetime = 1;          % if 1 don't compute mdd, 0 otherwise
model    = 1;          % 1 full sample estimation 2 drop subsample (bypass) estimation 3 Lenza-Primiceri (weigh) estimation

for dselect  = 7:7    % number refers to a specific real-time vintage available in "dselect" month

%======================================================================
%                          HYPERPARAMETERS
%======================================================================
load hyp.txt;           % selected hyperparameters from BVAR_PRIOR
lambda1=hyp(1);lambda2=hyp(2);lambda3=hyp(3);lambda4=hyp(4);


%======================================================================
%                    BAYESIAN ESTIMATION OF VAR
%======================================================================
% load data and set specification
vm_loaddata
Yq   = YQ(T0+1:nobs+T0,:);
Ym   = YM(T0+1:nobs+T0,:);

% prespecified federal fund rate/100
setff = 0.0005;

% read off parameter estimates
if model==1
    load(['Output\estimates_', num2str(dselect), '.mat'])
    nsim = size(Phip,1); SSpredsim = ones(nsim,H+1);
elseif model==2
    load(['Output\estimates_bypass_', num2str(dselect), '.mat'])
    nsim = size(Phip,1); SSpredsim = ones(nsim,H+1);
elseif model==3
    load(['Output\estimates_weigh_', num2str(dselect), '.mat'])
end


%======================================================================
%                     BLOCK 3: FORECASTING
%======================================================================
disp('                                                                ');
disp('                                                                ');
disp('                                                                ');
disp('                    MIXED FREQUENCY VAR: FORECASTING            ');
disp('                                                                ');
disp('                                                                ');
disp('                                                                ');

% store nowcasts(Tnew)+forecasts(H) in monthly frequency
YMforecastsim  = zeros(nsim,Tnew+H,Nm+Nq);
FMforecastsim  = zeros(nsim,Tnew+H,Nm+Nq);

gYMforecastsim = zeros(nsim,Tnew+H,Nm+Nq);
gFMforecastsim = zeros(nsim,Tnew+H,Nm+Nq);

% store nowcasts(Tnew)+forecasts(H) in quarterly frequency
HQ = floor((Tnew+H)/3);
YQforecastsim  = zeros(nsim,HQ,Nm+Nq);
FQforecastsim  = zeros(nsim,HQ,Nm+Nq);

gYQforecastsim = zeros(nsim,HQ,Nm+Nq);
gFQforecastsim = zeros(nsim,HQ,Nm+Nq);

% at the end of balanced panel for growth rate computation
YMactT = [Ym(end,:) Yq(end,:)];
YQactT = mean([Ym(end-2:end,:) Yq(end-2:end,:)]);
    
for jj=1:nsim
    
    
    YYact    = squeeze(YYactsim(jj,end,:));
    XXact    = squeeze(XXactsim(jj,end,:));    
    post_phi = squeeze(Phip(jj,:,:));
    post_sig = squeeze(Sigmap(jj,:,:));
    post_S   = SSpredsim(jj,:);
    
    % change elements in post_phi associated with federal funds rates
    Fpost_phi      = post_phi;
    Fpost_phi(:,6) = zeros(length(post_phi),1);
    Fpost_phi(6,:) = zeros(1,size(post_phi,2));
    Fpost_phi(6,6) = 1;
    
    YYpred              = zeros(H+1,nv);     % forecasts from VAR
    YYpred(1,:)         = YYact;
    XXpred              = zeros(H+1,(nv)*nlags+1);
    XXpred(:,end)       = ones(H+1,1);
    XXpred(1,:)         = XXact;
    
    FXXpred             = XXpred;
    FYYpred             = YYpred;
    
    % given posterior draw, draw #{H+1} random sequence
    error_pred = zeros(H+1,nv);     
    for h=1:H+1
        error_pred(h,:) = mvnrnd(zeros(nv,1), post_S(h)^2*post_sig);         
    end

    % given posterior draw, iterate forward to construct forecasts
    for h=2:H+1
        XXpred(h,nv+1:end-1)  = XXpred(h-1,1:end-nv-1);
        XXpred(h,1:nv)        = YYpred(h-1,:);
        YYpred(h,:)           = XXpred(h,:)*post_phi+error_pred(h,:);
        
        FXXpred(h,nv+1:end-1) = FXXpred(h-1,1:end-nv-1);
        FXXpred(h,1:nv)       = FYYpred(h-1,:);
        FYYpred(h,:)          = FXXpred(h,:)*Fpost_phi+error_pred(h,:);  
        if h==2
           FYYpred(h,:)       = YYpred(h,:);
           FYYpred(h,6)       = setff;
           FXXpred(h,:)       = XXpred(h,:); 
        end
    end
     
    if mod(jj,1000)==0
    disp('                                                               '); 
    disp(['                         FORECAST HORIZON:   ', num2str(H)]    );
    disp(['                         DRAW NUMBER:   ', num2str(jj)]        );
    disp(['                         REMAINING DRAWS:   ', num2str(nsim-jj)]);
    disp('                                                               ');
    end    
    
    
    % store in monthly frequency
    if Tnew==1
    YMall = [YMactT;[squeeze(YYactsim(jj,end,:))';YYpred(2:end,:)]];
    FMall = [YMactT;[squeeze(YYactsim(jj,end,:))';FYYpred(2:end,:)]];
    else
    YMall = [YMactT;[squeeze(YYactsim(jj,end-Tnew+1:end,:));YYpred(2:end,:)]];
    FMall = [YMactT;[squeeze(YYactsim(jj,end-Tnew+1:end,:));FYYpred(2:end,:)]];
    end

    % levels
    YMforecastsim(jj,:,:)  = YMall(2:end,:);
    FMforecastsim(jj,:,:)  = FMall(2:end,:);
    % growth rates
    gYMforecastsim(jj,:,:) = diff(YMall);
    gFMforecastsim(jj,:,:) = diff(FMall);

    % store in quarterly frequency: compute the averages
    YQall = zeros(HQ+1,size(YMall,2));
    FQall = zeros(HQ+1,size(YMall,2));
    YQall(1,:) = YQactT;
    FQall(1,:) = YQactT;    
    for hq=1:HQ
        YQall(hq+1,:) = mean(YMall(3*(hq-1)+2:3*hq+1,:));
        FQall(hq+1,:) = mean(FMall(3*(hq-1)+2:3*hq+1,:));
    end

    % levels
    YQforecastsim(jj,:,:)  = YQall(2:end,:);
    FQforecastsim(jj,:,:)  = FQall(2:end,:);
    
%     % growth rates: ln(X1)-ln(X0)
%     gYQforecastsim(jj,:,:) = diff(YQall);
%     gFQforecastsim(jj,:,:) = diff(FQall);   

    % growth rates: (X1-X0)/X0
    gYQforecastsim(jj,:,:) = vm_growthrates(YQall,select);
    gFQforecastsim(jj,:,:) = vm_growthrates(FQall,select);    

end

if model==1
    save(['Output\forecasts_', num2str(dselect), '.mat'])
elseif model==2
    save(['Output\forecasts_bypass_', num2str(dselect), '.mat'])
elseif model==3
    save(['Output\forecasts_weigh_', num2str(dselect), '.mat'])
end

end


path=ss;


