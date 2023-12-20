function [y95,y80,y50,y20,y05] = moment(IRF)

[nsim,H,n] = size(IRF);
y05 = zeros(H,n);
y20 = zeros(H,n);
y50 = squeeze(median(IRF));
y80 = zeros(H,n);
y95 = zeros(H,n);

for i=1:n
  d= sort(squeeze(IRF(:,:,i))); 
  y05(:,i) = d(round(0.05*nsim),:);
  y20(:,i) = d(round(0.2*nsim),:);
  y80(:,i) = d(round(0.8*nsim),:);
  y95(:,i) = d(round(0.95*nsim),:);
end
end
