clc;clear  

[YM0, YMX]  = xlsread('YM_Jan2022.xls');
[YQ0, YQX]  = xlsread('YQ_Jan2022.xls');

qtm = 1964.25:.25:2022;

select      = [1 0 0 0 0 1 1 0];
YM0(:,select==1) = YM0(:,select==1)./100;
YM0(:,select==0) = log(YM0(:,select==0));
YQ0         = log(YQ0(:,:));

nM = floor(size(YM0,1)/3);

YM1 = nan(nM,8);
for t=1:nM
YM1(t,:) = nanmean(YM0(3*(t-1)+1:3*t,:));
end

QDATA = nan(nM,11);
QDATA(:,1:8) = YM1;
QDATA(1:size(YQ0,1),9:end) = YQ0;

save('recentdata.mat','qtm','QDATA')


