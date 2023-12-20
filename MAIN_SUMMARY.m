clc;clear;
close all

ss = path;
path('Data',path);
path('Main Files',path);
path('Moment',path);
path('Gauss Files',path);
path('Output',path);
path('Figures',path);



% load forecasts stored in "Output" folder and save graphics in "Figures"
% this is an illustration:
for dselect = 7:7


model = 1; % 1 full sample estimation 2 drop subsample (bypass) estimation 3 Lenza-Primiceri (weigh) estimation

if model==1
    load(['Output\forecasts_', num2str(dselect), '.mat'])
elseif model==2
    load(['Output\forecasts_bypass_', num2str(dselect), '.mat'])
elseif model==3
    load(['Output\forecasts_weigh_', num2str(dselect), '.mat'])
end


vname = {'Unemployment Rate','Hours Worked','Consumer Price Index ',...
     'Industrial Production Index ','Personal Consumption Expenditure ',...
    'Federal Funds Rate','Treasury Bond Yield','S\&P 500 ',...
    'Gross Domestic Product ','Fixed Investment','Government Spending'};

% scale-back
m100  = [1 0 0 0 0 1 1 0 0 0 0];
cpi4  = [ 0 0 1 0 0 0 0 0 0 0 0];


yname = {'Log value','Annualized percent'};

load recentdata

%% data
YD = [Ym Yq];
QD = zeros(size(YD,1)/3,size(YD,2));
for ii=1:size(YD,1)/3
   QD(ii,:) = mean(YD(3*(ii-1)+1:3*ii,:)); 
end
gYD = diff(YD);
gQD = diff(QD);

% log difference as growth rates
%QDATAgr = diff(QDATA);
% (Y1-Y0)/Y0 computation for growth rates
QDATAgr = vm_growthrates(QDATA,select);
QDATAgr = [nan*QDATAgr(1,:);100*QDATAgr];


% set time
origin  = datevec(mdate(Torigin));
bench_l = [1 4 7 10];
bench_u = [3 6 9 12];
bench   = [1 2 3 4]/4;
addQ    = bench((bench_l<=origin(2) & origin(2)<=bench_u));
qtime   = 2020:0.25:2022;

qloc    = find(qtime==(origin(1)+addQ));
q1      = length(qtime(1:qloc));
q2      = length(qtime(qloc+1:end));

% most recent data
recentLV = nan(length(qtime),size(YD,2));
recentGR = recentLV;

recentLV(1:size(QDATA(find(qtm==qtime(1)):end,:),1),:) = QDATA(find(qtm==qtime(1)):end,:);
recentGR(1:size(QDATA(find(qtm==qtime(1)):end,:),1),:) = QDATAgr(find(qtm==qtime(1)):end,:);

recentLV(:,m100==1) = 100*recentLV(:,m100==1);

 
%% forecasts: levels

YMforecastsim(:,:,m100==1) = 100*YMforecastsim(:,:,m100==1);
FMforecastsim(:,:,m100==1) = 100*FMforecastsim(:,:,m100==1);
YQforecastsim(:,:,m100==1) = 100*YQforecastsim(:,:,m100==1);
FQforecastsim(:,:,m100==1) = 100*FQforecastsim(:,:,m100==1);
QD(:,m100==1)              = 100*QD(:,m100==1);
% monthly forecast
[lym95,lym80,lym50,lym20,lym05] = moment(YMforecastsim);
% monthly forecast fixed FF
[lfm95,lfm80,lfm50,lfm20,lfm05] = moment(FMforecastsim);

% quarterly forecast
[lyq95,lyq80,lyq50,lyq20,lyq05] = moment(YQforecastsim);
% quarterly forecast fixed FF
[lfq95,lfq80,lfq50,lfq20,lfq05] = moment(FQforecastsim);

LYQ05 = vm_forecast_presentation(QD,lyq05,q1,q2);
LYQ20 = vm_forecast_presentation(QD,lyq20,q1,q2);
LYQ50 = vm_forecast_presentation(QD,lyq50,q1,q2);
LYQ80 = vm_forecast_presentation(QD,lyq80,q1,q2);
LYQ95 = vm_forecast_presentation(QD,lyq95,q1,q2);

nmlz  = repmat(QD(end-q1+find(qtime==2020),:),size(LYQ50,1),1);

LFQ05 = vm_forecast_presentation(QD,lfq05,q1,q2);
LFQ20 = vm_forecast_presentation(QD,lfq20,q1,q2);
LFQ50 = vm_forecast_presentation(QD,lfq50,q1,q2);
LFQ80 = vm_forecast_presentation(QD,lfq80,q1,q2);
LFQ95 = vm_forecast_presentation(QD,lfq95,q1,q2);


%% forecasts: growth rates
% monthly forecast
[gym95,gym80,gym50,gym20,gym05] = moment(gYMforecastsim);
% monthly forecast fixed FF
[gfm95,gfm80,gfm50,gfm20,gfm05] = moment(gFMforecastsim);
% quarterly forecast
[gyq95,gyq80,gyq50,gyq20,gyq05] = moment(gYQforecastsim);
% quarterly forecast fixed FF
[gfq95,gfq80,gfq50,gfq20,gfq05] = moment(gFQforecastsim);

GYQ05 = 100*vm_forecast_presentation(QDATAgr(1:Torigin/3,:)/100,gyq05,q1,q2);
GYQ20 = 100*vm_forecast_presentation(QDATAgr(1:Torigin/3,:)/100,gyq20,q1,q2);
GYQ50 = 100*vm_forecast_presentation(QDATAgr(1:Torigin/3,:)/100,gyq50,q1,q2);
GYQ80 = 100*vm_forecast_presentation(QDATAgr(1:Torigin/3,:)/100,gyq80,q1,q2);
GYQ95 = 100*vm_forecast_presentation(QDATAgr(1:Torigin/3,:)/100,gyq95,q1,q2);

GYQ05(1,:) = 0;
GYQ20(1,:) = 0;
GYQ50(1,:) = 0;
GYQ80(1,:) = 0;
GYQ95(1,:) = 0;


GYQ05(:,cpi4==1) = 4*GYQ05(:,cpi4==1);
GYQ20(:,cpi4==1) = 4*GYQ20(:,cpi4==1);
GYQ50(:,cpi4==1) = 4*GYQ50(:,cpi4==1);
GYQ80(:,cpi4==1) = 4*GYQ80(:,cpi4==1);
GYQ95(:,cpi4==1) = 4*GYQ95(:,cpi4==1);

GFQ05 = 100*vm_forecast_presentation(QDATAgr(1:Torigin/3,:)/100,gfq05,q1,q2);
GFQ20 = 100*vm_forecast_presentation(QDATAgr(1:Torigin/3,:)/100,gfq20,q1,q2);
GFQ50 = 100*vm_forecast_presentation(QDATAgr(1:Torigin/3,:)/100,gfq50,q1,q2);
GFQ80 = 100*vm_forecast_presentation(QDATAgr(1:Torigin/3,:)/100,gfq80,q1,q2);
GFQ95 = 100*vm_forecast_presentation(QDATAgr(1:Torigin/3,:)/100,gfq95,q1,q2);

GFQ05(1,:) = 0;
GFQ20(1,:) = 0;
GFQ50(1,:) = 0;
GFQ80(1,:) = 0;
GFQ95(1,:) = 0;

GFQ05(:,cpi4==1) = 4*GFQ05(:,cpi4==1);
GFQ20(:,cpi4==1) = 4*GFQ20(:,cpi4==1);
GFQ50(:,cpi4==1) = 4*GFQ50(:,cpi4==1);
GFQ80(:,cpi4==1) = 4*GFQ80(:,cpi4==1);
GFQ95(:,cpi4==1) = 4*GFQ95(:,cpi4==1);

recentGR(1,:) = 0;
recentGR(:,cpi4==1) = 4*recentGR(:,cpi4==1);

scrsz = get(0,'ScreenSize');
set(0,'defaulttextinterpreter','latex')
set(0,'defaultaxesfontname','times')  
set(0,'DefaultAxesFontSize',20)
set(0,'DefaultTextFontSize',20);

for ii=1:nv
    fig = figure('Position',[20,20,900,600],'Name','','Color','w',...
      'Position',[1 scrsz(4)/13 scrsz(3)/4.5 scrsz(4)/3]);

    for jj=1:length(qtime)
        hold on
        line([qtime(jj) qtime(jj)],[-500 500],'color',[0.9 0.9 0.9],'linestyle','-')
    end
    hold on
    shadedplot(qtime,LYQ95(:,ii)',LYQ05(:,ii)',[0.9 0.9 0.9],[0.7 0.7 0.7]);
%     hold on
%     shadedplot(qtime,LYQ80(:,ii)',LYQ20(:,ii)',[0.7 0.7 0.7],[0.7 0.7 0.7]);
    hold on
%    k1=plot(qtime,LFQ50(:,ii),'Color',[0.3 0.3 0.9],'marker','s','markerfacecolor',[.3 .3 .9],'markeredgecolor',[0.3 0.3 0.9],'markersize',8);
    k2=plot(qtime,recentLV(:,ii),'Color',[0.9 0.3 0.3],'marker','o','markerfacecolor',[0.9 0.3 0.3],'markeredgecolor',[0.9 0.3 0.3],'markersize',10,'linewidth',2);
    plot(qtime,LYQ50(:,ii),'Color',[0.1 0.1 0.1],'linewidth',2);
    if m100(ii)==1
       ylabel('Percent') 
    end
    xticks([qtime(1) qtime(2) qtime(3) qtime(4) qtime(5) qtime(6) qtime(7) qtime(8) qtime(9)])
    xticklabels({'19:Q4','','20:Q2','','20:Q4','','21:Q2','','21:Q4'})
    axis([qtime(1) qtime(end) min(min(min(LYQ05(:,ii)),nanmin(LFQ50(:,ii))),nanmin(recentLV(:,ii)))...
                              max(max(max(LYQ95(:,ii)),nanmax(LFQ50(:,ii))),nanmax(recentLV(:,ii)))])
    
    print('-depsc2',['Figures\lv_', num2str(dselect),'_', num2str(ii), '.eps']);

end

for ii=1:nv
    fig = figure('Position',[20,20,900,600],'Name','','Color','w',...
      'Position',[1 scrsz(4)/13 scrsz(3)/4.5 scrsz(4)/3]);

    for jj=1:length(qtime)
        hold on
        line([qtime(jj) qtime(jj)],[-500 500],'color',[0.9 0.9 0.9],'linestyle','-')
    end
    hold on
    shadedplot(qtime,GYQ95(:,ii)',GYQ05(:,ii)',[0.9 0.9 0.9],[0.7 0.7 0.7]);
%     hold on
%     shadedplot(qtime,GYQ80(:,ii)',GYQ20(:,ii)',[0.7 0.7 0.7],[0.7 0.7 0.7]);
    hold on
%    k1=plot(qtime,GFQ50(:,ii),'Color',[0.3 0.3 0.9],'marker','s','markerfacecolor',[.3 .3 .9],'markeredgecolor',[0.3 0.3 0.9],'markersize',8);
    k2=plot(qtime,recentGR(:,ii),'Color',[0.9 0.3 0.3],'marker','o','markerfacecolor',[0.9 0.3 0.3],'markeredgecolor',[0.9 0.3 0.3],'markersize',10,'linewidth',2);
    plot(qtime,GYQ50(:,ii),'Color',[0.1 0.1 0.1],'linewidth',2);
    xticks([qtime(1) qtime(2) qtime(3) qtime(4) qtime(5) qtime(6) qtime(7) qtime(8) qtime(9)])
    xticklabels({'19:Q4','','20:Q2','','20:Q4','','21:Q2','','21:Q4'})
    axis([qtime(1) qtime(end) min(min(min(GYQ05(:,ii)),nanmin(GFQ50(:,ii))),nanmin(recentGR(:,ii)))...
                              max(max(max(GYQ95(:,ii)),nanmax(GFQ50(:,ii))),nanmax(recentGR(:,ii)))])
    ylabel('Percent')
    print('-depsc2',['Figures\gr_', num2str(dselect),'_', num2str(ii), '.eps']);
    
end


for ii=1:nv
    fig = figure('Position',[20,20,900,600],'Name','','Color','w',...
      'Position',[1 scrsz(4)/13 scrsz(3)/4.5 scrsz(4)/3]);

    for jj=1:length(qtime)
        hold on
        line([qtime(jj) qtime(jj)],[-500 500],'color',[0.9 0.9 0.9],'linestyle','-')
    end
    hold on
    shadedplot(qtime,cumsum(GYQ95(:,ii))',cumsum(GYQ05(:,ii))',[0.9 0.9 0.9],[0.7 0.7 0.7]);
%     hold on
%     shadedplot(qtime,cumsum(GYQ80(:,ii))',cumsum(GYQ20(:,ii))',[0.7 0.7 0.7],[0.7 0.7 0.7]);
    hold on
%    k1=plot(qtime,nancumsum(GFQ50(:,ii)),'Color',[0.3 0.3 0.9],'marker','s','markerfacecolor',[.3 .3 .9],'markeredgecolor',[0.3 0.3 0.9],'markersize',8);
    k2=plot(qtime,nancumsum(recentGR(:,ii)),'Color',[0.9 0.3 0.3],'marker','o','markerfacecolor',[0.9 0.3 0.3],'markeredgecolor',[0.9 0.3 0.3],'markersize',10,'linewidth',2);
    plot(qtime,cumsum(GYQ50(:,ii)),'Color',[0.1 0.1 0.1],'linewidth',2);
    xticks([qtime(1) qtime(2) qtime(3) qtime(4) qtime(5) qtime(6) qtime(7) qtime(8) qtime(9)])
    xticklabels({'19:Q4','','20:Q2','','20:Q4','','21:Q2','','21:Q4'})
    axis([qtime(1) qtime(end) min(min(min(cumsum(GYQ05(:,ii))),nanmin(cumsum(GFQ50(:,ii)))),nanmin(nancumsum(recentGR(:,ii))))...
                              max(max(max(cumsum(GYQ95(:,ii))),nanmax(cumsum(GFQ50(:,ii)))),nanmax(nancumsum(recentGR(:,ii))))])
    ylabel('Percent')
    print('-depsc2',['Figures\cgr_', num2str(dselect),'_', num2str(ii), '.eps']);
end

close all

% store forecasts in spreadsheet
year    = 2020:1:2021; quarter = 1:4;
time    = [[2019 4];[kron(year,ones(1,4))' repmat(quarter,1,length(year))']];

% row cell array (for column labels)
col_header = {'Year','Quarter','Unemployment Rate','Hours Worked','Consumer Price Index ',...
     'Industrial Production Index ','Personal Consumption Expenditure ',...
    'Federal Funds Rate','Treasury Bond Yield','S\&P 500 ',...
    'Gross Domestic Product ','Fixed Investment','Government Spending',...
    'Unemployment Rate','Hours Worked','Consumer Price Index ',...
     'Industrial Production Index ','Personal Consumption Expenditure ',...
    'Federal Funds Rate','Treasury Bond Yield','S\&P 500 ',...
    'Gross Domestic Product ','Fixed Investment','Government Spending'};    


% store level forecasts
xlswrite(['Output\LV_forecast_', num2str(dselect),'_', num2str(model),'.xls'],col_header(1:13),'Actual','A1');     
xlswrite(['Output\LV_forecast_', num2str(dselect),'_', num2str(model),'.xls'],[time recentLV],'Actual','A2');  

xlswrite(['Output\LV_forecast_', num2str(dselect),'_', num2str(model),'.xls'],col_header,'Median','A1');     
xlswrite(['Output\LV_forecast_', num2str(dselect),'_', num2str(model),'.xls'],[time LYQ50 LFQ50],'Median','A2');   

xlswrite(['Output\LV_forecast_', num2str(dselect),'_', num2str(model),'.xls'],col_header,'5%','A1');     
xlswrite(['Output\LV_forecast_', num2str(dselect),'_', num2str(model),'.xls'],[time LYQ05 LFQ05],'5%','A2');   

xlswrite(['Output\LV_forecast_', num2str(dselect),'_', num2str(model),'.xls'],col_header,'20%','A1');     
xlswrite(['Output\LV_forecast_', num2str(dselect),'_', num2str(model),'.xls'],[time LYQ20 LFQ20],'20%','A2');   

xlswrite(['Output\LV_forecast_', num2str(dselect),'_', num2str(model),'.xls'],col_header,'80%','A1');     
xlswrite(['Output\LV_forecast_', num2str(dselect),'_', num2str(model),'.xls'],[time LYQ80 LFQ80],'80%','A2');   

xlswrite(['Output\LV_forecast_', num2str(dselect),'_', num2str(model),'.xls'],col_header,'95%','A1');     
xlswrite(['Output\LV_forecast_', num2str(dselect),'_', num2str(model),'.xls'],[time LYQ95 LFQ95],'95%','A2');   


% store growth rate forecasts
xlswrite(['Output\GR_forecast_', num2str(dselect),'_', num2str(model),'.xls'],col_header(1:13),'Actual','A1');     
xlswrite(['Output\GR_forecast_', num2str(dselect),'_', num2str(model),'.xls'],[time recentGR],'Actual','A2');  

xlswrite(['Output\GR_forecast_', num2str(dselect),'_', num2str(model),'.xls'],col_header,'Median','A1');     
xlswrite(['Output\GR_forecast_', num2str(dselect),'_', num2str(model),'.xls'],[time GYQ50 GFQ50],'Median','A2'); 

xlswrite(['Output\GR_forecast_', num2str(dselect),'_', num2str(model),'.xls'],col_header,'5%','A1');     
xlswrite(['Output\GR_forecast_', num2str(dselect),'_', num2str(model),'.xls'],[time GYQ05 GFQ05],'5%','A2'); 

xlswrite(['Output\GR_forecast_', num2str(dselect),'_', num2str(model),'.xls'],col_header,'20%','A1');     
xlswrite(['Output\GR_forecast_', num2str(dselect),'_', num2str(model),'.xls'],[time GYQ20 GFQ20],'20%','A2'); 

xlswrite(['Output\GR_forecast_', num2str(dselect),'_', num2str(model),'.xls'],col_header,'80%','A1');     
xlswrite(['Output\GR_forecast_', num2str(dselect),'_', num2str(model),'.xls'],[time GYQ80 GFQ80],'80%','A2'); 

xlswrite(['Output\GR_forecast_', num2str(dselect),'_', num2str(model),'.xls'],col_header,'95%','A1');     
xlswrite(['Output\GR_forecast_', num2str(dselect),'_', num2str(model),'.xls'],[time GYQ95 GFQ95],'95%','A2'); 

% store cumulative growth rate forecasts
xlswrite(['Output\CGR_forecast_', num2str(dselect),'_', num2str(model),'.xls'],col_header(1:13),'Actual','A1');     
xlswrite(['Output\CGR_forecast_', num2str(dselect),'_', num2str(model),'.xls'],[time cumsum(recentGR)],'Actual','A2');  

xlswrite(['Output\CGR_forecast_', num2str(dselect),'_', num2str(model),'.xls'],col_header,'Median','A1');     
xlswrite(['Output\CGR_forecast_', num2str(dselect),'_', num2str(model),'.xls'],[time cumsum(GYQ50) cumsum(GFQ50)],'Median','A2'); 

xlswrite(['Output\CGR_forecast_', num2str(dselect),'_', num2str(model),'.xls'],col_header,'5%','A1');     
xlswrite(['Output\CGR_forecast_', num2str(dselect),'_', num2str(model),'.xls'],[time cumsum(GYQ05) cumsum(GFQ05)],'5%','A2'); 

xlswrite(['Output\CGR_forecast_', num2str(dselect),'_', num2str(model),'.xls'],col_header,'20%','A1');     
xlswrite(['Output\CGR_forecast_', num2str(dselect),'_', num2str(model),'.xls'],[time cumsum(GYQ20) cumsum(GFQ20)],'20%','A2'); 

xlswrite(['Output\CGR_forecast_', num2str(dselect),'_', num2str(model),'.xls'],col_header,'80%','A1');     
xlswrite(['Output\CGR_forecast_', num2str(dselect),'_', num2str(model),'.xls'],[time cumsum(GYQ80) cumsum(GFQ80)],'80%','A2'); 

xlswrite(['Output\CGR_forecast_', num2str(dselect),'_', num2str(model),'.xls'],col_header,'95%','A1');     
xlswrite(['Output\CGR_forecast_', num2str(dselect),'_', num2str(model),'.xls'],[time cumsum(GYQ95) cumsum(GFQ95)],'95%','A2'); 


end

path=ss;




