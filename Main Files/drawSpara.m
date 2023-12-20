function [lh,alpha,beta,SS,SSpred] = drawSpara(spara,AT_draw,CONF,PHIF,SIGF,rhoprior,nv,Tnobs,Tnew,H)

% Estimating s parameters
alpha = nv/2;
beta  = zeros(2,1);

lh    = 0;
SSall = ones(Tnobs+H,1);

for tt=Tnobs-Tnew+3:Tnobs+H
    if tt==Tnobs-Tnew+3
        SSall(tt,:) = spara(1);
    elseif tt==Tnobs-Tnew+4
        SSall(tt,:) = spara(2);
    elseif tt==Tnobs-Tnew+5
        SSall(tt,:) = spara(3);
    else
        SSall(tt,:) = 1+(spara(3)-1)*spara(4)^(tt-(Tnobs-Tnew+5));
    end  
end
SS     = SSall(1:Tnobs);
SSpred = SSall(Tnobs:end);

for tt=Tnobs-Tnew+3:Tnobs
    if tt<Tnobs-Tnew+5
       nut  = AT_draw(tt-(Tnobs-Tnew)+1,:)'-CONF-PHIF*AT_draw(tt-(Tnobs-Tnew),:)';
       beta(tt-(Tnobs-Tnew)-2) = 0.5*trace(inv(SIGF(1:nv,1:nv))*(nut(1:nv)*nut(1:nv)'));
    end
       nut  = AT_draw(tt-(Tnobs-Tnew)+1,:)'-CONF-PHIF*AT_draw(tt-(Tnobs-Tnew),:)';
       lh   = lh - log(SS(tt)^(nv/2)) - 0.5*(1/(SS(tt)^2))*trace(inv(SIGF(1:nv,1:nv))*(nut(1:nv)*nut(1:nv)'));
end

lh = lh - log(spara(3)^2) - lpdfbeta(spara(4), rhoprior.p, rhoprior.q);


