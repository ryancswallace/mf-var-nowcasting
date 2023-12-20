function [ pdf, logpdf ] = invwishpdf(IW,S,d) 

[k,k2] = size(IW) ;

logexpterm = -.5*trace(S/IW) ;
logdetIWterm = log(det(IW))*(d/2) ;
logdetSterm = log(det(S))*((d-k-1)/2) ;
logtwoterm = log(2)*((d-k-1)*k/2) ;
logpiterm = log(pi)*((k-1)*k/4) ;

klst = 1:k ;
dkk2 = (d-k-klst)/2 ;
% ind  = dkk2>0;
% dkk2 = dkk2(ind==1);
dkk2(end) = 1e-20;


gamln = gammaln(dkk2) ;
sumgamln = sum(gamln) ;

logpdf = logexpterm + logdetSterm - ...
         (logdetIWterm + logtwoterm + logpiterm + sumgamln  ) ;

pdf = exp(logpdf) ;
