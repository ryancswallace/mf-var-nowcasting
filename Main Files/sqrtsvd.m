function [sqrtPPP, ind] = sqrtsvd(PPP)

tol = 1e-15;

[a,b,c] = svd(PPP);

ind = sum(diag(b)>tol);

aa = a(:,1:ind);
bb = b(1:ind,1:ind);
% cc = c(:,1:ind);
% 
% sum(sum(abs(PPP - aa*bb*cc')))

sqrtPPP = aa*sqrt(bb);