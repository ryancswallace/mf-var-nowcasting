function ldens = lpdfbeta(x, a, b)
%  log beta pdf
% ldens = lngam(a+b) - lngam(a) - lngam(b) + (a-1)*log(x) + (b-1)*log(1-x);
ldens = gammaln(a+b) - gammaln(a) - gammaln(b) + (a-1)*log(x) + (b-1)*log(1-x);