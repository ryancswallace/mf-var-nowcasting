%======================================================================
%
%        MF BVAR : MIXED FREQ BAYESIAN VAR INFERENCE                      
%
%        Date : September 2023
%
%        Note : Given selected hyperparameters, this code performs  
%               Bayesian Inference of MF-VARs     
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
fixpara  = 0;          % fixpara should be always 0 when running MF_BVAR_ESTIMATION.m
smodel   = 1;          % 1 full sample estimation 2 drop subsample (bypass) estimation 3 Lenza-Primiceri (weigh) estimation


for dselect  = 7:7     % number refers to a specific real-time vintage available in "dselect" month 

%======================================================================
%                          HYPERPARAMETERS
%======================================================================
load hyp.txt;           % selected hyperparameters from BVAR_PRIOR
hyp = [0.09 4.3 2.7 4.3]';
lambda1=hyp(1);lambda2=hyp(2);lambda3=hyp(3);lambda4=hyp(4);


%======================================================================
%                    BAYESIAN ESTIMATION OF VAR
%======================================================================
if smodel==1
    % load data and set specification
    vm_loaddata % "xlsread" will run properly in windows not in mac
    
    % estimation
    SUB_MF_BVAR_ESTIMATION_FULL
    
    save(['Output\estimates_', num2str(dselect), '.mat'],'Phip','Sigmap','YYactsim','XXactsim')
elseif smodel==2
    % load data and set specification
    vm_loaddata_bypass % "xlsread" will run properly in windows not in mac
    
    % estimation
    SUB_MF_BVAR_ESTIMATION_BYPASS
    
    save(['Output\estimates_bypass_', num2str(dselect), '.mat'],'Phip','Sigmap','YYactsim','XXactsim')
elseif smodel==3
    % load data and set specification
    vm_loaddata_primiceri % "xlsread" will run properly in windows not in mac
    
    % initialization
    SSold = ones(size(YDATA,1),1);
    sparaold = [1.5 5 3 0.8];
    % rho prior
    pmean  = 0.8; pstdd = 0.2; rhoprior.p = (1-pmean)*pmean^2/pstdd^2 - pmean; rhoprior.q = rhoprior.p*(1/pmean -1);
    lubound = [[1 1000];
        [0    1]];
    H_S     = 1e-4*diag([10000 10]);
    % estimation
    SUB_MF_BVAR_ESTIMATION_WEIGH
    
    save(['Output\estimates_weigh_', num2str(dselect), '.mat'],'Phip','Sigmap','SSsim','SSpredsim','Sparasim','YYactsim','XXactsim')
end

end
path=ss;


