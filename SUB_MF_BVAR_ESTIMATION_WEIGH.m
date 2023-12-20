disp('                                                                ');
disp('                    MIXED FREQUENCY VAR: ESTIMATION             ');
disp('                                                                ');
disp(['                           VINTAGE NUMBER:   ', num2str(dselect)] );
disp('                                                               ');
disp('                          PARAMETER ESTIMATION                  ');
disp('                                                                ');
%======================================================================
%                     BLOCK 1: PARAMETER ESTIMATION 
%======================================================================

% Matrices for collecting draws from Posterior Density


lstate    = zeros(nsim,Nq,Tnobs);
YYactsim  = zeros(nsim,Tnew,nv);
XXactsim  = zeros(nsim,Tnew,nv*p+1);
likesim   = zeros(nsim,1);
zmddsim   = zeros(nsim,1);    
priorsim  = zeros(nsim,2);
At_mat    = zeros(Tnobs,Nq*(p+1));
Pt_mat    = zeros(Tnobs,(Nq*(p+1))^2);
Atildemat = zeros(nsim,Nq*(p+1));
Ptildemat = zeros(nsim,Nq*(p+1),Nq*(p+1));

loglh   = 0;

spec_init   = [nlags_ T0 nex nv Nm];
if fixpara==0
    rej       = zeros(nsim,1);
    SSsim     = zeros(nsim,Tnobs);
    SSpredsim = zeros(nsim,H+1);
    Sparasim  = zeros(nsim,4);
    Sigmap    = zeros(nsim,nv,nv);
    Phip      = zeros(nsim,nv*p+1,nv);
    Cons      = zeros(nsim,nv);
    parasim   = zeros(nsim,(nv*p+1)*nv+nv*(nv+1)/2);
    % prior variance: initialization
    [Phi,sigma] = vm_prior_init(hyp,YDATA,spec_init);
else
    Phi   = squeeze(Phip(1,:,:));
    sigma = squeeze(Sigmap(1,:,:));
end

init_mean = [];
for it = 1:p+1
    init_mean = [YQ(it,:)';init_mean];
end
    

%% Define Coefficient Matrices
% Define phi(qm), phi(qq), phi(qc)
phi_qm = zeros(Nm*p,Nq);
for i=1:p
    phi_qm(Nm*(i-1)+1:Nm*i,:)=Phi((i-1)*(Nm+Nq)+1:(i-1)*(Nm+Nq)+Nm,Nm+1:end);
end
phi_qq = zeros(Nq*p,Nq);
for i=1:p
    phi_qq(Nq*(i-1)+1:Nq*i,:)=Phi((i-1)*(Nm+Nq)+Nm+1:i*(Nm+Nq),Nm+1:end);
end
phi_qc = Phi(end,Nm+1:end);

% Define phi(mm), phi(mq), phi(mc)
phi_mm = zeros(Nm*p,Nm);
for i=1:p
    phi_mm(Nm*(i-1)+1:Nm*i,:)=Phi((i-1)*(Nm+Nq)+1:(i-1)*(Nm+Nq)+Nm,1:Nm);
end
phi_mq = zeros(Nq*p,Nm);
for i=1:p
    phi_mq(Nq*(i-1)+1:Nq*i,:)=Phi((i-1)*(Nm+Nq)+Nm+1:i*(Nm+Nq),1:Nm);
end
phi_mc = Phi(end,1:Nm);

% Define Covariance Term sig_mm, sig_mq, sig_qm, sig_qq
sig_mm  = sigma(1:Nm,1:Nm);
sig_mq  = 0.5*(sigma(1:Nm,Nm+1:end)+sigma(Nm+1:end,1:Nm)');
sig_qm  = 0.5*(sigma(Nm+1:end,1:Nm)+sigma(1:Nm,Nm+1:end)');
sig_qq  = sigma(Nm+1:end,Nm+1:end);

% Define Transition Equation Matrices
GAMMAs  = [[phi_qq' zeros(Nq,Nq)];[eye(p*Nq,p*Nq) zeros(p*Nq,Nq)]];
GAMMAz  = [phi_qm'; zeros(p*Nq,p*Nm)];
GAMMAc  = [phi_qc'; zeros(p*Nq,1)];
GAMMAu  = [eye(Nq); zeros(p*Nq,Nq)];

% Define Measurement Equation Matrices
LAMBDAs = [[zeros(Nm,Nq) phi_mq'];(1/3)*[eye(Nq) eye(Nq) eye(Nq) zeros(Nq,Nq*(p-2))]];
LAMBDAz = [phi_mm'; zeros(Nq,p*Nm)];
LAMBDAc = [phi_mc'; zeros(Nq,1)];
LAMBDAu = [eye(Nm); zeros(Nq,Nm)];

% Define W matrix
Wmatrix   = [eye(Nm) zeros(Nm,Nq)];
LAMBDAs_t = Wmatrix * LAMBDAs;
LAMBDAz_t = Wmatrix * LAMBDAz;
LAMBDAc_t = Wmatrix * LAMBDAc;
LAMBDAu_t = Wmatrix * LAMBDAu;

%% Initialization   
init_var  = zeros(Nq*(p+1),Nq*(p+1));
for kk=1:p+1
    init_var = GAMMAs * init_var * GAMMAs' + GAMMAu * sig_qq * GAMMAu';
end

Zm_init   = zeros(T0-p-1,Nm*p);
i=1;
while i<=p
    Zm_init(:,(i-1)*Nm+1:i*Nm) = YM(p+1-(i-1):T0-i,:);
    i=i+1;
end

[init_At,init_Pt] = vm_initialize(GAMMAs,GAMMAz,GAMMAc,GAMMAu,...
LAMBDAs,LAMBDAz,LAMBDAc,LAMBDAu,LAMBDAs_t,LAMBDAz_t,LAMBDAc_t,LAMBDAu_t,...
sig_qq,sig_mm,sig_qm,sig_mq,Zm_init,YDATA,init_mean,init_var,spec_init);

%% Set Dataset
% Lagged Monthly Observations
Zm   = zeros(nobs,Nm*p);
i=1;
while i<=p
    Zm(:,(i-1)*Nm+1:i*Nm) = YM(T0-(i-1):T0+nobs-i,:);
    i=i+1;
end

% Observations in Monthly Freq
Yq   = YQ(T0+1:nobs+T0,:);
Ym   = YM(T0+1:nobs+T0,:);


%% Estimation
% Block Algorithm
for j=1:nsim
    
    if mod(j,1000)==0
    disp('                                                               '); 
    disp('                                                               ');
    disp(['                          DRAW NUMBER:   ', num2str(j)]        );
    disp(['                          REMAINING DRAWS:   ',num2str(nsim-j)]);
    disp('                                                               ');
    end 
    
    % Initialization         
    At   = init_At;  
    Pt   = init_Pt;
    loglh = 0;
    
% Kalman Filter loop 
for t = 1:nobs            % note that t=T0+t originally
   
   if ((t+T0)/3)-floor((t+T0)/3)==0
       
       At1 = At;
       Pt1 = Pt;

       % Forecasting   
       alphahat = GAMMAs * At1 + GAMMAz * Zm(t,:)' + GAMMAc;
       Phat     = GAMMAs * Pt1 * GAMMAs' + GAMMAu * sig_qq * GAMMAu';
       Phat     = 0.5*(Phat+Phat');
       
       yhat     = LAMBDAs * alphahat + LAMBDAz * Zm(t,:)' + LAMBDAc;
       nut      = [Ym(t,:)'; Yq(t,:)'] - yhat;
       
       Ft       = LAMBDAs * Phat * LAMBDAs' + LAMBDAu * sig_mm * LAMBDAu'...
                  + LAMBDAs*GAMMAu*sig_qm*LAMBDAu'...
                  + LAMBDAu*sig_mq*GAMMAu'*LAMBDAs';
       Ft       = 0.5*(Ft+Ft');
       
       loglh    = loglh-0.5*length([Ym(t,:)'; Yq(t,:)'])*log(2*pi)-0.5*log(det(Ft))-0.5*nut'*inv(Ft)*nut;    
       
       Xit      = LAMBDAs * Phat + LAMBDAu * sig_mq * GAMMAu';       
       
       At       = alphahat + Xit' * inv(Ft) * nut;
       Pt       = Phat     - Xit' * inv(Ft) * Xit;             
       
       At_mat(t,:)  = At';
       Pt_mat(t,:)  = reshape(Pt,1,(Nq*(p+1))^2);
   
   else
       
       At1 = At;
       Pt1 = Pt;
       
       % Forecasting   
       alphahat = GAMMAs * At1 + GAMMAz * Zm(t,:)' + GAMMAc;
       Phat     = GAMMAs * Pt1 * GAMMAs' + GAMMAu * sig_qq * GAMMAu';
       Phat     = 0.5*(Phat+Phat');
       
       yhat     = LAMBDAs_t * alphahat + LAMBDAz_t * Zm(t,:)' + LAMBDAc_t;
       nut      = Ym(t,:)' - yhat;
       
       Ft       = LAMBDAs_t * Phat * LAMBDAs_t' + LAMBDAu_t * sig_mm * LAMBDAu_t'...
                  + LAMBDAs_t*GAMMAu*sig_qm*LAMBDAu_t'...
                  + LAMBDAu_t*sig_mq*GAMMAu'*LAMBDAs_t';
       Ft       = 0.5*(Ft+Ft');
       
       loglh    = loglh-0.5*length(Ym(t,:)')*log(2*pi)-0.5*log(det(Ft))-0.5*nut'*inv(Ft)*nut; 
      
       Xit      = LAMBDAs_t * Phat + LAMBDAu_t * sig_mq * GAMMAu';
       
       At       = alphahat + Xit' * inv(Ft) * nut;
       Pt       = Phat     - Xit' * inv(Ft) * Xit; 

       At_mat(t,:)  = At';
       Pt_mat(t,:)  = reshape(Pt,1,(Nq*(p+1))^2);
       
   end
   
end

Atildemat(j,:) = At_mat(nobs,:);
Pt_last  = reshape(Pt_mat(nobs,:),Nq*(p+1),Nq*(p+1));
Ptildemat(j,:,:) = Pt_last;


%% unbalanced dataset

kn        = nv*(p+1);

% Measurement Equation
Z1         = zeros(Nm,kn);
Z1(:,1:Nm) = eye(Nm);

Z2         = zeros(Nq,kn);
for bb=1:Nq
    for ll=1:3
        Z2(bb,ll*Nm+(ll-1)*Nq+bb)=1/3;
    end
end
ZZ         = [Z1;Z2];

BAt = [];
for rr=1:p+1
    BAt=[BAt;[Ym(end-rr+1,:) squeeze(Atildemat(j,(rr-1)*Nq+1:rr*Nq))]'];
end

BPt = zeros(kn,kn);
for rr=1:p+1
    for vv=1:p+1
    BPt(rr*Nm+(rr-1)*Nq+1:rr*(Nm+Nq),vv*Nm+(vv-1)*Nq+1:vv*(Nm+Nq))=...
    squeeze(Ptildemat(j,(rr-1)*Nq+1:rr*Nq,(vv-1)*Nq+1:vv*Nq));
    end
end

BAt_mat=zeros(Tnobs,kn);
BPt_mat=zeros(Tnobs,kn^2);

BAt_mat(nobs,:) = BAt;
BPt_mat(nobs,:) = reshape(BPt,1,kn^2);

% Define Companion Form Matrix PHIF
PHIF         = zeros(kn,kn);          
IF           = eye(nv);
for i=1:p
    PHIF(i*nv+1:(i+1)*nv,(i-1)*nv+1:i*nv) = IF;
end
PHIF(1:nv,1:nv*p) = Phi(1:end-1,:)';

% Define Constant Term CONF
CONF = [Phi(end,:)';zeros(nv*p,1)];

% Define Covariance Term SIGF
SIGF = zeros(kn,kn);
SIGF(1:nv,1:nv) = sigma;

% Filter loop 
for t = nobs+1:Tnobs
       
       % New indicator
       kkk = t-nobs;
       
       % Define New Data (ND) and New Z matrix (NZ)
       ND  = [delif(YDATA(nobs+T0+kkk,:)',index_NY(:,kkk))];
       NZ  = delif(ZZ,index_NY(:,kkk));       
      
       BAt1 = BAt;
       BPt1 = BPt;

       % Forecasting   
       Balphahat = PHIF * BAt1 + CONF;
       BPhat     = PHIF * BPt1 * PHIF' + SSold(t)^2*SIGF; %% SS(t) is the Lenza-Primiceri weight
       BPhat     = 0.5*(BPhat+BPhat');
      
       if isempty(ND)==0
           Byhat = NZ*Balphahat;
           Bnut  = ND - Byhat;

           BFt = NZ*BPhat*NZ';
           BFt = 0.5*(BFt+BFt');

           loglh    = loglh-0.5*length(ND)*log(2*pi)-0.5*log(det(BFt))-0.5*Bnut'*inv(BFt)*Bnut;

           % Updating
           BAt = Balphahat + (BPhat*NZ')*inv(BFt)*Bnut;
           BPt = BPhat - (BPhat*NZ')*inv(BFt)*(BPhat*NZ')';
       else
           BAt = Balphahat;
           BPt = BPhat;
       end

       BAt_mat(t,:)  = BAt';
       BPt_mat(t,:)  = reshape(BPt,1,kn^2);
   
end

likesim(j,:) = loglh;

AT_draw = zeros(Tnew+1,kn);

% singular value decomposition
[u s v] = svd(reshape(BPt_mat(Tnobs,:),kn,kn));
Pchol = u*sqrt(s);
AT_draw(end,:) = BAt_mat(Tnobs,:)+(Pchol*randn(kn,1))';


% Kalman Smoother 
for i = 1:Tnew
    
    BAtt  = BAt_mat(Tnobs-i,:)';
    BPtt  = reshape(BPt_mat(Tnobs-i,:),kn,kn);

    BPhat = PHIF * BPtt * PHIF' + SSold(Tnobs-i)^2*SIGF; %% SS(Tnobs-i) is the Lenza-Primiceri weight
    BPhat = 0.5*(BPhat+BPhat');
    
    [Up, Sp, Vp] = svd(BPhat);
    sp = diag(Sp); kp = sum(sp>1e-15);
    inv_BPhat = (Up(:,1:kp)*diag(1./sp(1:kp))*Vp(:,1:kp)')';
    
    Bnut  = AT_draw(end-i+1,:)'-PHIF*BAtt -CONF;
      
    Amean = BAtt + (BPtt*PHIF')*inv_BPhat*Bnut;
    Pmean = BPtt - (BPtt*PHIF')*inv_BPhat*(BPtt*PHIF')';   
  
    % singular value decomposition
    [um sm vm] = svd(Pmean);
    Pmchol = um*sqrt(sm);
    AT_draw(end-i,:) = (Amean+Pmchol*randn(kn,1))';  
    
end


%% balanced dataset

At_draw = zeros(nobs,Nq*(p+1));

for kk=1:p+1
    At_draw(nobs,(kk-1)*Nq+1:kk*Nq) = ...
        AT_draw(1,kk*Nm+(kk-1)*Nq+1:kk*(Nm+Nq));
end

% Kalman Smoother 
for i = 1:nobs-1
    
    Att  = At_mat(nobs-i,:)';
    Ptt  = reshape(Pt_mat(nobs-i,:),Nq*(p+1),Nq*(p+1));

    Phat = GAMMAs * Ptt * GAMMAs' + GAMMAu * sig_qq * GAMMAu';
    Phat = 0.5*(Phat+Phat');
    
    [Up, Sp, Vp] = svd(Phat);
    sp = diag(Sp); kp = sum(sp>1e-12);
    inv_Phat = (Up(:,1:kp)*diag(1./sp(1:kp))*Vp(:,1:kp)')';
    
    nut  = At_draw(nobs-i+1,:)'-GAMMAs * Att - GAMMAz * Zm(nobs-i,:)'- GAMMAc;
      
    Amean = Att + (Ptt*GAMMAs')*inv_Phat*nut;
    Pmean = Ptt - (Ptt*GAMMAs')*inv_Phat*(Ptt*GAMMAs')';   
  
    % singular value decomposition
    [um sm vm] = svd(Pmean);
    Pmchol = um*sqrt(sm);
    At_draw(nobs-i,:) = (Amean+Pmchol*randn(Nq*(p+1),1))';  
end

%======================================================================
%                     BLOCK 2: MINNESOTA PRIOR
%======================================================================

% update Dataset YY
YY = [[Ym At_draw(:,1:Nq)];AT_draw(2:end,1:(Nm+Nq))];

% save latent states
for hh=1:Nq
    lstate(j,hh,1:nobs)=At_draw(:,hh);
    lstate(j,hh,nobs+1:end)=AT_draw(2:end,Nm+hh);
end 

nobs_   = size(YY,1)-T0;
spec    = [nlags_ T0 nex nv nobs_];

% dummy observations and actual observations
[mdd,YYact,YYdum,XXact,XXdum]=vm_mdd2(hyp,YY,spec,savetime);
zmddsim(j,:)    = mdd;
YYactsim(j,:,:) = YYact(end-Tnew+1:end,:);
XXactsim(j,:,:) = XXact(end-Tnew+1:end,:);
    
% draws from posterior distribution
[Tdummy,n] = size(YYdum);
[Tobs,n]   = size(YYact);
X          = [XXact./SSold(end-size(XXact,1)+1:end); XXdum];
Y          = [YYact./SSold(end-size(XXact,1)+1:end); YYdum];
p          = spec(1);                 % Number of lags in the VAR
T          = Tobs+Tdummy;             % Number of observations
F          = zeros(n*p,n*p);          % Matrix for Companion Form
I          = eye(n);

for i=1:p-1
    F(i*n+1:(i+1)*n,(i-1)*n+1:i*n) = I;
end



if fixpara==0
    inv_x     = inv(X'*X);
    Phi_tilde = (inv_x)*X'*Y;
    Sigma     = (Y-X*Phi_tilde)'*(Y-X*Phi_tilde);
    
    % Draws from the density Sigma | Y
    sigma   = iwishrnd(Sigma,T-n*p-1);
    
    % Draws from the density vec(Phi) |Sigma(j), Y
    phi_new = mvnrnd(reshape(Phi_tilde,n*(n*p+1),1),kron(sigma,inv_x));
    
    % Rearrange vec(Phi) into Phi
    Phi     = reshape(phi_new,n*p+1,n);
    
    % evaluate prior
    [MN_logpdf, IW_logpdf] = vm_prior_pdf(hyp,YY,spec,Phi,sigma);
    Sigmap(j,:,:) = sigma;
    Phip(j,:,:)   = Phi;
    Cons(j,:)     = Phi(end,:);
    priorsim(j,:) = [MN_logpdf IW_logpdf];
    parasim(j,1:(nv*p+1)*nv) = reshape(Phi,1,(nv*p+1)*nv);
    clear interm
    interm = [];
    for ii=1:nv; 
        interm = [interm sigma(ii,ii:end)];    
    end
    parasim(j,(nv*p+1)*nv+1:end) = interm;
else
    Phi   = squeeze(Phip(j,:,:)); 
    sigma = squeeze(Sigmap(j,:,:));
end

% Define phi(qm), phi(qq), phi(qc)
phi_qm = zeros(Nm*p,Nq);
for i=1:p
    phi_qm(Nm*(i-1)+1:Nm*i,:)=Phi((i-1)*(Nm+Nq)+1:(i-1)*(Nm+Nq)+Nm,Nm+1:end);
end
phi_qq = zeros(Nq*p,Nq);
for i=1:p
    phi_qq(Nq*(i-1)+1:Nq*i,:)=Phi((i-1)*(Nm+Nq)+Nm+1:i*(Nm+Nq),Nm+1:end);
end
phi_qc = Phi(end,Nm+1:end);

% Define phi(mm), phi(mq), phi(mc)
phi_mm = zeros(Nm*p,Nm);
for i=1:p
    phi_mm(Nm*(i-1)+1:Nm*i,:)=Phi((i-1)*(Nm+Nq)+1:(i-1)*(Nm+Nq)+Nm,1:Nm);
end
phi_mq = zeros(Nq*p,Nm);
for i=1:p
    phi_mq(Nq*(i-1)+1:Nq*i,:)=Phi((i-1)*(Nm+Nq)+Nm+1:i*(Nm+Nq),1:Nm);
end
phi_mc = Phi(end,1:Nm);

% Define Covariance Term sig_mm, sig_mq, sig_qm, sig_qq
sig_mm  = sigma(1:Nm,1:Nm);
sig_mq  = 0.5*(sigma(1:Nm,Nm+1:end)+sigma(Nm+1:end,1:Nm)');
sig_qm  = 0.5*(sigma(Nm+1:end,1:Nm)+sigma(1:Nm,Nm+1:end)');
sig_qq  = sigma(Nm+1:end,Nm+1:end);

% Define Transition Equation Matrices
GAMMAs  = [[phi_qq' zeros(Nq,Nq)];[eye(p*Nq,p*Nq) zeros(p*Nq,Nq)]];
GAMMAz  = [phi_qm'; zeros(p*Nq,p*Nm)];
GAMMAc  = [phi_qc'; zeros(p*Nq,1)];
GAMMAu  = [eye(Nq); zeros(p*Nq,Nq)];

% Define Measurement Equation Matrices
LAMBDAs = [[zeros(Nm,Nq) phi_mq'];(1/3)*[eye(Nq) eye(Nq) eye(Nq) zeros(Nq,Nq*(p-2))]];
LAMBDAz = [phi_mm'; zeros(Nq,p*Nm)];
LAMBDAc = [phi_mc'; zeros(Nq,1)];
LAMBDAu = [eye(Nm); zeros(Nq,Nm)];

LAMBDAs_t = Wmatrix * LAMBDAs;
LAMBDAz_t = Wmatrix * LAMBDAz;
LAMBDAc_t = Wmatrix * LAMBDAc;
LAMBDAu_t = Wmatrix * LAMBDAu;
  
% draw S parameters by Lenza and Primiceri 
if fixpara==0
   [lhold,alpha_S,beta_S,~,SSpredold] = drawSpara(sparaold,AT_draw,CONF,PHIF,SIGF,rhoprior,nv,Tnobs,Tnew,H);
    s0 = max(sqrt(gamrnd(alpha_S, beta_S(1)^(-1))^(-1)),1);
    s1 = max(sqrt(gamrnd(alpha_S, beta_S(2)^(-1))^(-1)),1);
    % propose   
    genacc=0;
    while genacc<1
        sparanew = mvnrnd(sparaold(3:4), H_S,1)';
        genacc  = all(sparanew < lubound(:,2) & sparanew > lubound(:,1));  
    end
    sparaold = [s0;s1;reshape(sparaold(3:4),[],1)]; 
    sparanew = [s0;s1;reshape(sparanew,[],1)];      
    [lhnew,~,~,SSnew,SSprednew] = drawSpara(sparanew,AT_draw,CONF,PHIF,SIGF,rhoprior,nv,Tnobs,Tnew,H);
 
    r = min([1 exp(lhnew-lhold)]);            
    if rand > r 
       % reject proposed jump             
       rej(j) = 1;
    else
       % accept proposed jump  
       sparaold  = sparanew;
       SSold     = SSnew;
       SSpredold = SSprednew;
    end
    Sparasim(j,:) = sparaold;
    SSsim(j,:)    = SSold;
    SSpredsim(j,:)= SSpredold;
end    
    
end