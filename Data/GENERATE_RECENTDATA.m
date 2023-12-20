clc;clear;
load recentdata

cgr_data = 100*cumsum(exp(QDATA(end-8:end,:))./exp(QDATA(end-9:end-1,:))-1); 
cgr_data = cgr_data-cgr_data(1,:);

gr_data = 400*(exp(QDATA(end-8:end,:))./exp(QDATA(end-9:end-1,:))-1); 
gr_data(1,:) = 0;

lv_data = QDATA(end-8:end,:);
lv_data(:,[1 6 7]) = 100*lv_data(:,[1 6 7]);

save('forecast_realdata.mat','cgr_data','gr_data','lv_data');