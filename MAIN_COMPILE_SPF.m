clc;clear
path('Data',path);
path('SPF',path);

%% real-time data vintage

time         = 2019.75:.25:2021.75;
qtm          = 1964.25:.25:2021;
% 2020:Q1 
[YM1, YMX1]  = xlsread('YM_Jan2020.xls');
[YQ1, YQX1]  = xlsread('YQ_Jan2020.xls');
QDATA1       = construct_recentdata(YM1,YQ1);
d2019q4_2020q1 = QDATA1(find(qtm==2020),:);
% 2020:Q2
[YM2, YMX2]  = xlsread('YM_April2020.xls');
[YQ2, YQX2]  = xlsread('YQ_April2020.xls');
QDATA2       = construct_recentdata(YM2,YQ2);
d2019q4_2020q2 = QDATA2(find(qtm==2020),:);
% 2020:Q3
[YM3, YMX3]  = xlsread('YM_July2020.xls');
[YQ3, YQX3]  = xlsread('YQ_July2020.xls');
QDATA3       = construct_recentdata(YM3,YQ3);
d2019q4_2020q3 = QDATA3(find(qtm==2020),:);
% 2020:Q4
[YM4, YMX4]  = xlsread('YM_Oct2020.xls');
[YQ4, YQX4]  = xlsread('YQ_Oct2020.xls');
QDATA4       = construct_recentdata(YM4,YQ4);
d2019q4_2020q4 = QDATA4(find(qtm==2020),:);
% 2021:Q1
[YM5, YMX5]  = xlsread('YM_Jan2021.xls');
[YQ5, YQX5]  = xlsread('YQ_Jan2021.xls');
QDATA5       = construct_recentdata(YM5,YQ5);
d2019q4_2021q1 = QDATA5(find(qtm==2020),:);
% 2021:Q2
[YM6, YMX6]  = xlsread('YM_Apr2021.xls');
[YQ6, YQX6]  = xlsread('YQ_Apr2021.xls');
QDATA6       = construct_recentdata(YM6,YQ6);
d2019q4_2021q2 = QDATA6(find(qtm==2020),:);
% 2021:Q3
[YM7, YMX7]  = xlsread('YM_July2021.xls');
[YQ7, YQX7]  = xlsread('YQ_July2021.xls');
QDATA7       = construct_recentdata(YM7,YQ7);
d2019q4_2021q3 = QDATA7(find(qtm==2020),:);
% 2021:Q4
[YM8, YMX8]  = xlsread('YM_Oct2021.xls');
[YQ8, YQX8]  = xlsread('YQ_Oct2021.xls');
QDATA7       = construct_recentdata(YM8,YQ8);
d2019q4_2021q4 = QDATA7(find(qtm==2020),:);

%% spf median forecasts
[spf_rgdp,etc_rgdp]     = xlsread('Median_RGDP_Level.xlsx');
[spf_cpi, etc_cpi]      = xlsread('Median_CPI_Level.xlsx');
[spf_indpro,etc_indpro] = xlsread('Median_INDPROD_Level.xlsx');
[spf_tbond,etc_tbond]   = xlsread('Median_TBOND_Level.xlsx');
[spf_unrate,etc_unrate] = xlsread('Median_UNEMP_Level.xlsx');
[spf_rpce,etc_rpce]     = xlsread('Median_RCONSUM_Level.xlsx');
[spf_govs,etc_govs]     = xlsread('Median_RFEDGOV_Level.xlsx');
[spf_nriv,etc_nriv]     = xlsread('Median_RNRESIN_Level.xlsx');
[spf_rinv,etc_rinv]     = xlsread('Median_RRESINV_Level.xlsx');

% Q1 forecasts based on Q4 information, Feb released to public (Jan questionnaire after 1st NIPA release)
% Q2 forecasts based on Q1 information, May (Apr)
% Q3 forecasts based on Q2 information, Aug (Jul)
% Q4 forecasts based on Q3 information, Nov (Oc)

% from 2020:q1 onward
rgdp_spf   = log(spf_rgdp((spf_rgdp(:,1)>2019)==1,find(ismember(etc_rgdp(1,:),'RGDP2')==1):find(ismember(etc_rgdp(1,:),'RGDP6')==1)));
cpi_spf    = spf_cpi((spf_cpi(:,1)>2019)==1,find(ismember(etc_cpi(1,:),'CPI2')==1):find(ismember(etc_cpi(1,:),'CPI6')==1));
indpro_spf = log(spf_indpro((spf_indpro(:,1)>2019)==1,find(ismember(etc_indpro(1,:),'INDPROD2')==1):find(ismember(etc_indpro(1,:),'INDPROD6')==1)));
tbond_spf  = spf_tbond((spf_tbond(:,1)>2019)==1,find(ismember(etc_tbond(1,:),'TBOND2')==1):find(ismember(etc_tbond(1,:),'TBOND6')==1))/100;
unrate_spf = spf_unrate((spf_unrate(:,1)>2019)==1,find(ismember(etc_unrate(1,:),'UNEMP2')==1):find(ismember(etc_unrate(1,:),'UNEMP6')==1))/100;
rpce_spf   = log(spf_rpce((spf_rpce(:,1)>2019)==1,find(ismember(etc_rpce(1,:),'RCONSUM2')==1):find(ismember(etc_rpce(1,:),'RCONSUM6')==1)));
govs_spf   = log(spf_govs((spf_govs(:,1)>2019)==1,find(ismember(etc_govs(1,:),'RFEDGOV2')==1):find(ismember(etc_govs(1,:),'RFEDGOV6')==1)));

nriv_spf   = spf_nriv((spf_nriv(:,1)>2019)==1,find(ismember(etc_nriv(1,:),'RNRESIN2')==1):find(ismember(etc_nriv(1,:),'RNRESIN6')==1));
rinv_spf   = spf_rinv((spf_rinv(:,1)>2019)==1,find(ismember(etc_rinv(1,:),'RRESINV2')==1):find(ismember(etc_rinv(1,:),'RRESINV6')==1));
inv_spf    = log(nriv_spf+rinv_spf);

vname = {'Unemployment Rate','Hours Worked','Consumer Price Index ',...
     'Industrial Production Index ','Personal Consumption Expenditure ',...
    'Federal Funds Rate','Treasury Bond Yield','S\&P 500 ',...
    'Gross Domestic Product ','Fixed Investment','Government Spending'};

% subtract 2019:Q4 values from real-time data, but don't seem to match with SPF level
% jj = 9;
% rgdp_2019q4   = [d2019q4_2020q1(jj);d2019q4_2020q2(jj);d2019q4_2020q3(jj);d2019q4_2020q4(jj);d2019q4_2021q1(jj);d2019q4_2021q2(jj);d2019q4_2021q3(jj);d2019q4_2021q4(jj)];
% jj = 1;
% unrate_2019q4 = [d2019q4_2020q1(jj);d2019q4_2020q2(jj);d2019q4_2020q3(jj);d2019q4_2020q4(jj);d2019q4_2021q1(jj);d2019q4_2021q2(jj);d2019q4_2021q3(jj);d2019q4_2021q4(jj)];
% jj = 3;
% cpi_2019q4    = [d2019q4_2020q1(jj);d2019q4_2020q2(jj);d2019q4_2020q3(jj);d2019q4_2020q4(jj);d2019q4_2021q1(jj);d2019q4_2021q2(jj);d2019q4_2021q3(jj);d2019q4_2021q4(jj)];
% jj = 4;
% indpro_2019q4 = [d2019q4_2020q1(jj);d2019q4_2020q2(jj);d2019q4_2020q3(jj);d2019q4_2020q4(jj);d2019q4_2021q1(jj);d2019q4_2021q2(jj);d2019q4_2021q3(jj);d2019q4_2021q4(jj)];
% jj = 7;
% tbond_2019q4  = [d2019q4_2020q1(jj);d2019q4_2020q2(jj);d2019q4_2020q3(jj);d2019q4_2020q4(jj);d2019q4_2021q1(jj);d2019q4_2021q2(jj);d2019q4_2021q3(jj);d2019q4_2021q4(jj)];
% jj = 5;
% rpce_2019q4   = [d2019q4_2020q1(jj);d2019q4_2020q2(jj);d2019q4_2020q3(jj);d2019q4_2020q4(jj);d2019q4_2021q1(jj);d2019q4_2021q2(jj);d2019q4_2021q3(jj);d2019q4_2021q4(jj)];
% jj = 11;
% govs_2019q4   = [d2019q4_2020q1(jj);d2019q4_2020q2(jj);d2019q4_2020q3(jj);d2019q4_2020q4(jj);d2019q4_2021q1(jj);d2019q4_2021q2(jj);d2019q4_2021q3(jj);d2019q4_2021q4(jj)];

rgdp_2019q4     = log(spf_rgdp((spf_rgdp(:,1)==2020 & spf_rgdp(:,2)==1)==1,find(ismember(etc_rgdp(1,:),'RGDP1')==1)))*ones(8,1);
indpro_2019q4   = log(spf_indpro((spf_indpro(:,1)==2020 & spf_indpro(:,2)==1)==1,find(ismember(etc_indpro(1,:),'INDPROD1')==1)))*ones(8,1);
rpce_2019q4     = log(spf_rpce((spf_rpce(:,1)==2020 & spf_rpce(:,2)==1)==1,find(ismember(etc_rpce(1,:),'RCONSUM1')==1)))*ones(8,1);
govs_2019q4     = log(spf_govs((spf_govs(:,1)==2020 & spf_govs(:,2)==1)==1,find(ismember(etc_govs(1,:),'RFEDGOV1')==1)))*ones(8,1);
unrate_2019q4   = spf_unrate((spf_unrate(:,1)==2020 & spf_unrate(:,2)==1)==1,find(ismember(etc_unrate(1,:),'UNEMP1')==1))*ones(8,1);
cpi_2019q4      = spf_cpi((spf_cpi(:,1)==2020 & spf_cpi(:,2)==1)==1,find(ismember(etc_cpi(1,:),'CPI1')==1))*ones(8,1);
tbond_2019q4    = spf_tbond((spf_tbond(:,1)==2020 & spf_tbond(:,2)==1)==1,find(ismember(etc_tbond(1,:),'TBOND1')==1))*ones(8,1);

nriv_2019q4     = spf_nriv((spf_nriv(:,1)==2020 & spf_nriv(:,2)==1)==1,find(ismember(etc_nriv(1,:),'RNRESIN1')==1))*ones(8,1);
rinv_2019q4     = spf_rinv((spf_rinv(:,1)==2020 & spf_rinv(:,2)==1)==1,find(ismember(etc_rinv(1,:),'RRESINV1')==1))*ones(8,1);
inv_2019q4      = log(nriv_2019q4+rinv_2019q4);

%
rgdp_spf20q2 = nan(length(time),1);
rgdp_spf20q2(find(time==2020.25):find(time==2021.25)) = 100*(rgdp_spf(2,:)-rgdp_2019q4(2));
rgdp_spf20q3 = nan(length(time),1);
rgdp_spf20q3(find(time==2020.5):find(time==2021.5))   = 100*(rgdp_spf(3,:)-rgdp_2019q4(3));
rgdp_spf20q4 = nan(length(time),1);
rgdp_spf20q4(find(time==2020.75):find(time==2021.75)) = 100*(rgdp_spf(4,:)-rgdp_2019q4(4));
rgdp_spf21q1 = nan(length(time),1);
rgdp_spf21q1(find(time==2021.00):find(time==2021.75)) = 100*(rgdp_spf(5,1:end-1)-rgdp_2019q4(5));
rgdp_spf21q2 = nan(length(time),1);
rgdp_spf21q2(find(time==2021.25):find(time==2021.75)) = 100*(rgdp_spf(6,1:end-2)-rgdp_2019q4(6));
rgdp_spf21q3 = nan(length(time),1);
rgdp_spf21q3(find(time==2021.50):find(time==2021.75)) = 100*(rgdp_spf(7,1:end-3)-rgdp_2019q4(7));

%
ip_spf20q2 = nan(length(time),1);
ip_spf20q2(find(time==2020.25):find(time==2021.25)) = 100*(indpro_spf(2,:)-indpro_2019q4(2));
ip_spf20q3 = nan(length(time),1);
ip_spf20q3(find(time==2020.5):find(time==2021.5))   = 100*(indpro_spf(3,:)-indpro_2019q4(3));
ip_spf20q4 = nan(length(time),1);
ip_spf20q4(find(time==2020.75):find(time==2021.75)) = 100*(indpro_spf(4,:)-indpro_2019q4(4));
ip_spf21q1 = nan(length(time),1);
ip_spf21q1(find(time==2021.00):find(time==2021.75)) = 100*(indpro_spf(5,1:end-1)-indpro_2019q4(5));
ip_spf21q2 = nan(length(time),1);
ip_spf21q2(find(time==2021.25):find(time==2021.75)) = 100*(indpro_spf(6,1:end-2)-indpro_2019q4(6));
ip_spf21q3 = nan(length(time),1);
ip_spf21q3(find(time==2021.50):find(time==2021.75)) = 100*(indpro_spf(7,1:end-3)-indpro_2019q4(7));

%
cpi_spf20q2 = nan(length(time),1);
cpi_spf20q2(find(time==2020.25):find(time==2021.25)) = cpi_spf(2,:);
cpi_spf20q3 = nan(length(time),1);
cpi_spf20q3(find(time==2020.5):find(time==2021.5))   = cpi_spf(3,:);
cpi_spf20q4 = nan(length(time),1);
cpi_spf20q4(find(time==2020.75):find(time==2021.75)) = cpi_spf(4,:);
cpi_spf21q1 = nan(length(time),1);
cpi_spf21q1(find(time==2021.00):find(time==2021.75)) = cpi_spf(5,1:end-1);
cpi_spf21q2 = nan(length(time),1);
cpi_spf21q2(find(time==2021.25):find(time==2021.75)) = cpi_spf(6,1:end-2);
cpi_spf21q3 = nan(length(time),1);
cpi_spf21q3(find(time==2021.50):find(time==2021.75)) = cpi_spf(7,1:end-3);

%
unrate_spf20q2 = nan(length(time),1);
unrate_spf20q2(find(time==2020.25):find(time==2021.25)) = 100*unrate_spf(2,:);
unrate_spf20q3 = nan(length(time),1);
unrate_spf20q3(find(time==2020.5):find(time==2021.5))   = 100*unrate_spf(3,:);
unrate_spf20q4 = nan(length(time),1);
unrate_spf20q4(find(time==2020.75):find(time==2021.75)) = 100*unrate_spf(4,:);
unrate_spf21q1 = nan(length(time),1);
unrate_spf21q1(find(time==2021.00):find(time==2021.75)) = 100*unrate_spf(5,1:end-1);
unrate_spf21q2 = nan(length(time),1);
unrate_spf21q2(find(time==2021.25):find(time==2021.75)) = 100*unrate_spf(6,1:end-2);
unrate_spf21q3 = nan(length(time),1);
unrate_spf21q3(find(time==2021.50):find(time==2021.75)) = 100*unrate_spf(7,1:end-3);

%
tbond_spf20q2 = nan(length(time),1);
tbond_spf20q2(find(time==2020.25):find(time==2021.25)) = 100*tbond_spf(2,:);
tbond_spf20q3 = nan(length(time),1);
tbond_spf20q3(find(time==2020.5):find(time==2021.5))   = 100*tbond_spf(3,:);
tbond_spf20q4 = nan(length(time),1);
tbond_spf20q4(find(time==2020.75):find(time==2021.75)) = 100*tbond_spf(4,:);
tbond_spf21q1 = nan(length(time),1);
tbond_spf21q1(find(time==2021.00):find(time==2021.75)) = 100*tbond_spf(5,1:end-1);
tbond_spf21q2 = nan(length(time),1);
tbond_spf21q2(find(time==2021.25):find(time==2021.75)) = 100*tbond_spf(6,1:end-2);
tbond_spf21q3 = nan(length(time),1);
tbond_spf21q3(find(time==2021.50):find(time==2021.75)) = 100*tbond_spf(7,1:end-3);

%
rpce_spf20q2 = nan(length(time),1);
rpce_spf20q2(find(time==2020.25):find(time==2021.25)) = 100*(rpce_spf(2,:)-rpce_2019q4(2));
rpce_spf20q3 = nan(length(time),1);
rpce_spf20q3(find(time==2020.5):find(time==2021.5))   = 100*(rpce_spf(3,:)-rpce_2019q4(3));
rpce_spf20q4 = nan(length(time),1);
rpce_spf20q4(find(time==2020.75):find(time==2021.75)) = 100*(rpce_spf(4,:)-rpce_2019q4(4));
rpce_spf21q1 = nan(length(time),1);
rpce_spf21q1(find(time==2021.00):find(time==2021.75)) = 100*(rpce_spf(5,1:end-1)-rpce_2019q4(5));
rpce_spf21q2 = nan(length(time),1);
rpce_spf21q2(find(time==2021.25):find(time==2021.75)) = 100*(rpce_spf(6,1:end-2)-rpce_2019q4(6));
rpce_spf21q3 = nan(length(time),1);
rpce_spf21q3(find(time==2021.50):find(time==2021.75)) = 100*(rpce_spf(7,1:end-3)-rpce_2019q4(7));

%
govs_spf20q2 = nan(length(time),1);
govs_spf20q2(find(time==2020.25):find(time==2021.25)) = 100*(govs_spf(2,:)-govs_2019q4(2));
govs_spf20q3 = nan(length(time),1);
govs_spf20q3(find(time==2020.5):find(time==2021.5))   = 100*(govs_spf(3,:)-govs_2019q4(3));
govs_spf20q4 = nan(length(time),1);
govs_spf20q4(find(time==2020.75):find(time==2021.75)) = 100*(govs_spf(4,:)-govs_2019q4(4));
govs_spf21q1 = nan(length(time),1);
govs_spf21q1(find(time==2021.00):find(time==2021.75)) = 100*(govs_spf(5,1:end-1)-govs_2019q4(5));
govs_spf21q2 = nan(length(time),1);
govs_spf21q2(find(time==2021.25):find(time==2021.75)) = 100*(govs_spf(6,1:end-2)-govs_2019q4(6));
govs_spf21q3 = nan(length(time),1);
govs_spf21q3(find(time==2021.50):find(time==2021.75)) = 100*(govs_spf(7,1:end-3)-govs_2019q4(7));

%
inv_spf20q2 = nan(length(time),1);
inv_spf20q2(find(time==2020.25):find(time==2021.25)) = 100*(inv_spf(2,:)-inv_2019q4(2));
inv_spf20q3 = nan(length(time),1);
inv_spf20q3(find(time==2020.5):find(time==2021.5))   = 100*(inv_spf(3,:)-inv_2019q4(3));
inv_spf20q4 = nan(length(time),1);
inv_spf20q4(find(time==2020.75):find(time==2021.75)) = 100*(inv_spf(4,:)-inv_2019q4(4));
inv_spf21q1 = nan(length(time),1);
inv_spf21q1(find(time==2021.00):find(time==2021.75)) = 100*(inv_spf(5,1:end-1)-inv_2019q4(5));
inv_spf21q2 = nan(length(time),1);
inv_spf21q2(find(time==2021.25):find(time==2021.75)) = 100*(inv_spf(6,1:end-2)-inv_2019q4(6));
inv_spf21q3 = nan(length(time),1);
inv_spf21q3(find(time==2021.50):find(time==2021.75)) = 100*(inv_spf(7,1:end-3)-inv_2019q4(7));

%% save files: cgr - rgdp & indpro; gr - cpi; lv - unrate & tbond


rgdp_spf_all = [repmat(rgdp_spf20q2,1,3)... % 4, 5, 6
                repmat(rgdp_spf20q3,1,3)... % 7, 8, 9
                repmat(rgdp_spf20q4,1,3)... % 10,11,12
                repmat(rgdp_spf21q1,1,3)... % 13,14,15
                repmat(rgdp_spf21q2,1,3)... % 16,17,18
                repmat(rgdp_spf21q3,1,3)];  % 19,20,21

indpro_spf_all = [repmat(ip_spf20q2,1,3)... % 4, 5, 6
                  repmat(ip_spf20q3,1,3)... % 7, 8, 9
                  repmat(ip_spf20q4,1,3)... % 10,11,12
                  repmat(ip_spf21q1,1,3)... % 13,14,15
                  repmat(ip_spf21q2,1,3)... % 16,17,18
                  repmat(ip_spf21q3,1,3)];  % 19,20,21

cpi_spf_all  = [repmat(cpi_spf20q2,1,3)... % 4, 5, 6
                repmat(cpi_spf20q3,1,3)... % 7, 8, 9
                repmat(cpi_spf20q4,1,3)... % 10,11,12
                repmat(cpi_spf21q1,1,3)... % 13,14,15
                repmat(cpi_spf21q2,1,3)... % 16,17,18
                repmat(cpi_spf21q3,1,3)];  % 19,20,21

unr_spf_all  = [repmat(unrate_spf20q2,1,3)... % 4, 5, 6
                repmat(unrate_spf20q3,1,3)... % 7, 8, 9
                repmat(unrate_spf20q4,1,3)... % 10,11,12
                repmat(unrate_spf21q1,1,3)... % 13,14,15
                repmat(unrate_spf21q2,1,3)... % 16,17,18
                repmat(unrate_spf21q3,1,3)];  % 19,20,21

tb_spf_all   = [repmat(tbond_spf20q2,1,3)... % 4, 5, 6
                repmat(tbond_spf20q3,1,3)... % 7, 8, 9
                repmat(tbond_spf20q4,1,3)... % 10,11,12
                repmat(tbond_spf21q1,1,3)... % 13,14,15
                repmat(tbond_spf21q2,1,3)... % 16,17,18
                repmat(tbond_spf21q3,1,3)];  % 19,20,21

con_spf_all  = [repmat(rpce_spf20q2,1,3)... % 4, 5, 6
                repmat(rpce_spf20q3,1,3)... % 7, 8, 9
                repmat(rpce_spf20q4,1,3)... % 10,11,12
                repmat(rpce_spf21q1,1,3)... % 13,14,15
                repmat(rpce_spf21q2,1,3)... % 16,17,18
                repmat(rpce_spf21q3,1,3)];  % 19,20,21
 
gov_spf_all  = [repmat(govs_spf20q2,1,3)... % 4, 5, 6
                repmat(govs_spf20q3,1,3)... % 7, 8, 9
                repmat(govs_spf20q4,1,3)... % 10,11,12
                repmat(govs_spf21q1,1,3)... % 13,14,15
                repmat(govs_spf21q2,1,3)... % 16,17,18
                repmat(govs_spf21q3,1,3)];  % 19,20,21
            
inv_spf_all  = [repmat(inv_spf20q2,1,3)... % 4, 5, 6
                repmat(inv_spf20q3,1,3)... % 7, 8, 9
                repmat(inv_spf20q4,1,3)... % 10,11,12
                repmat(inv_spf21q1,1,3)... % 13,14,15
                repmat(inv_spf21q2,1,3)... % 16,17,18
                repmat(inv_spf21q3,1,3)];  % 19,20,21
            
cgr_spf = []; gr_spf = []; lv_spf = [];

%             1    2    3   4  5   6   7   8    9    10    11     
% variables: UNR,HOURS,CPI,IP,PCE,FFR,TB,SP500,GDP,FIXINV,GOV

for vintage=4:20
    cgr_spf{vintage} = nan(length(time),11);
    cgr_spf{vintage}(:,4)  = indpro_spf_all(:,vintage-3);
    cgr_spf{vintage}(:,9)  = rgdp_spf_all(:,vintage-3);
    cgr_spf{vintage}(:,5)  = con_spf_all(:,vintage-3);
    cgr_spf{vintage}(:,10) = inv_spf_all(:,vintage-3);
    cgr_spf{vintage}(:,11) = gov_spf_all(:,vintage-3);

    gr_spf{vintage} = nan(length(time),11);
    gr_spf{vintage}(:,3) = cpi_spf_all(:,vintage-3);

    lv_spf{vintage} = nan(length(time),11);
    lv_spf{vintage}(:,1) = unr_spf_all(:,vintage-3);
    lv_spf{vintage}(:,7) = tb_spf_all(:,vintage-3);
end

save('Output\SPFforecast.mat','cgr_spf','gr_spf','lv_spf')

%%
function QDATA = construct_recentdata(YM0,YQ0)

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

end
