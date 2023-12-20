function [mdd, store1, store2, index] = vm_mdd_final(lstate,zmddsim,nsim,nobs,densfac,sD)
rng(sD);

Nq    = size(lstate,2);
lstate0  = lstate(:,:,1:nobs);
Nv       = size(lstate0,2);
NN       = nobs/3;

clear lstateA 
lstateA = [];

for nn=1:Nv
pick1    = repmat([1 0 0],1,nobs/3);
pick2    = repmat([0 1 0],1,nobs/3);
pick3    = repmat([0 0 1],1,nobs/3);

lstate1  = lstate0(:,nn,pick1==1);
lstate2  = lstate0(:,nn,pick2==1);
lstate3  = lstate0(:,nn,pick3==1);

lstate1  = reshape(lstate1,nsim,NN);
lstate2  = reshape(lstate2,nsim,NN);
lstate3  = reshape(lstate3,nsim,NN);

lstateA1 = lstate3-lstate2;
lstateA2 = lstate2-lstate1;

temp     = [lstateA1 lstateA2];

lstateA  = [lstateA temp];
end

parasim  = lstateA(0.5*nsim+1:end,:);
zmddsim  = zmddsim(0.5*nsim+1:end,:);


%==========================================================================
% Computation of the marginal likelihood 
%==========================================================================

%==========================================================================
% default parameters
%==========================================================================

% for data density (modified harmonic mean)
%densfac = 4000;
hmax    = 20;

%==========================================================================
% means and standard deviations
%==========================================================================

% initialize output
ndraws     = 0;
drawmean   = 0;
drawsqmean = 0;
sumdrawsq  = 0;
	 
postsim = zmddsim;
     
% number of simulations in each block
nsimul = size(parasim,1);
npara  = size(parasim,2);
%npara  = npara*3/2;

% collect simulations
drawblock = parasim;
drawdim   = size(drawblock,2);

% compute sums of x(i) and x(i)^2 to be used to calculate means and s.d.
drawmean   = drawmean + sum(drawblock,1);
drawsqmean = drawsqmean + sum(drawblock.^2,1);
sumdrawsq  = sumdrawsq + parasim'*parasim;
ndraws     = ndraws + nsimul;

% compute means and standard deviations
drawmean   = drawmean/ndraws;
drawsqmean = drawsqmean/ndraws;
drawstdd   = sqrt(drawsqmean - drawmean.^2);
%drawsig    = diag(drawstdd.^2);

drawsig    = sumdrawsq/ndraws - drawmean'*drawmean;

% [up sp vp] = svd(drawsig);
% drawsiginv = zeros(size(drawsig));
% drawsigdet = zeros(size(drawsig));
% 
% for j=1:size(drawsig,1)
%     if sp(j,j) > 1e-12
%         drawsiginv(j,j) = 1./sp(j,j);
%         drawsigdet(j,j) = sp(j,j);
%     else
%         drawsigdet(j,j) = 1;
%     end
% end

[Up, Sp, Vp] = svd(drawsig);
sp = diag(Sp); kp = sum(sp>1e-7);
drawsiginv = (Up(:,1:kp)*diag(1./sp(1:kp))*Vp(:,1:kp)')';
drawsigdet = prod(sp(1:kp));  
if drawsigdet==0
    drawsigdet = 1e-300;
end

drawsiglndet = sum(log(diag(drawsigdet)));

% drawsiglndet = sum(log(diag(drawsig)));
% drawsiginv   = inv(drawsig);

%==========================================================================
% marginal data density:  modified harmonic mean by Geweke (1999)
%==========================================================================

p = (.1:.1:.9)';
pcrit = chi2inv(p,ones(length(p),1)*npara);

ndraws  = 0;
suminvlike = zeros(length(p),1);
laginvlike = zeros(hmax,length(p));
gaminvlike = zeros(hmax,length(p));
	    
% number of simulations in each block
[nsimul,npara] = size(parasim);				
%npara  = npara*3/2;
paradev  = parasim-repmat(drawmean,nsimul,1);
quadpara = sum((paradev*drawsiginv).*paradev,2);

store1 = zeros(nsimul,length(p));
store2 = zeros(nsimul,length(p));

index = zeros(nsimul,length(p));
invlikevec = zeros(nsimul,length(p));

densfac = - 0.5*npara*log(2*pi) - 0.5*drawsiglndet - 0.5*quadpara(1) - log(p) - postsim(1);
densfac = -mean(densfac);
% simulation loop
for j=1:nsimul
	lnfpara = - 0.5*npara*log(2*pi) - 0.5*drawsiglndet - 0.5*quadpara(j) - log(p);
	indpara = (quadpara(j)<pcrit);
    store1(j,:) = (lnfpara );
    store2(j,:) = (postsim(j));
    index(j,:) = indpara;
	invlike = exp(lnfpara - postsim(j) + densfac).*indpara;
    invlikevec(j,:) = (lnfpara - postsim(j)).*indpara;
	laginvlike = [invlike' ; laginvlike(1:hmax-1,:)];
	gaminvlike = gaminvlike + laginvlike.*repmat(invlike',hmax,1);
	suminvlike = suminvlike + invlike;
end	% simulation loop

ndraws = ndraws + nsimul;
meaninvlike = suminvlike/ndraws;
%gaminvlike  = gaminvlike/ndraws - repmat((meaninvlike.^2)',hmax,1);

% % standard error
% suminvlikeerror = gaminvlike(1,:);
% for k=2:hmax
%    suminvlikeerror = suminvlikeerror + 2*gaminvlike(k,:)*(1-(k-1)/hmax);
% end
% 
% suminvlikeerror = 100*sqrt(suminvlikeerror/ndraws)./meaninvlike' ;

% marginalized data density
mdd = densfac-log(meaninvlike);
%mdd = [p mdd];

mdd = mean(mdd);