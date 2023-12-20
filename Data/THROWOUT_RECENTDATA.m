clc;clear  

[YM0, YMX]  = xlsread('YM_Jan2022.xls');
[YQ0, YQX]  = xlsread('YQ_Jan2022.xls');

qtm = 1964.25:.25:2021.75;
mtm = 1964+1/12:1/12:2022;

select      = [1 0 0 0 0 1 1 0];
YM0(:,select==1) = YM0(:,select==1)./100;
YM0(:,select==0) = log(YM0(:,select==0));
YQ0         = log(YQ0(:,:));

YM0 = YM0(1:end-1,:);
gYM = 1200*diff(YM0);
gYQ = 400*diff(YQ0);


load recession.txt
xx = sum(recession,2)/2;
ww = recession(:,2) - recession(:,1);

folder_save = 'C:\Users\dsong14\Dropbox\MF-VAR Resuscitation\Documents\Figures';

vname = {'Unemployment Rate','Hours Worked','Consumer Price Index ',...
     'Industrial Production Index ','Personal Consumption Expenditure ',...
    'Federal Funds Rate','Treasury Bond Yield','S\&P 500 ',...
    'Gross Domestic Product ','Fixed Investment','Government Spending'};

scrsz = get(0,'ScreenSize');
set(0,'defaulttextinterpreter','latex')
set(0,'defaultaxesfontname','times')  
set(0,'DefaultAxesFontSize',20)
set(0,'DefaultTextFontSize',20);

for ii=1:1
    meann = mean(gYM(:,ii));
    stdev = std(gYM(:,ii));
    
    fig = figure('Position',[20,20,900,600],'Name','','Color','w',...
      'Position',[1 scrsz(4)/13 scrsz(3)/4.5 scrsz(4)/3]);
    hold on
    for kk = 1:length(xx), bar(xx(kk),15000,ww(kk),'FaceColor',[0.8 0.8 0.8],'EdgeColor','w'), end
    for kk = 1:length(xx), bar(xx(kk),-15000,ww(kk),'FaceColor',[0.8 0.8 0.8],'EdgeColor','w'), end
    hold on
    plot(mtm(2:end),ones(length(mtm)-1,1)*(meann-4*stdev),'r--');
    plot(mtm(2:end),ones(length(mtm)-1,1)*(meann+4*stdev),'r--');
    plot(mtm(2:end),gYM(:,ii),'Color',[0.3 0.3 0.9],'linewidth',2)
    axis([mtm(2) mtm(end) min(min(gYM(:,ii)),min(meann-5*stdev)) max(max(gYM(:,ii)),max(meann+5*stdev))])
    ylabel('\%')
    print('-depsc2',[fullfile(folder_save, 'mdatagr_'), num2str(ii), '.eps']);
    
    ind1 = gYM(:,ii)<meann-4*stdev; 
    ind2 = gYM(:,ii)>meann+4*stdev; 
    indd = [0;(ind1+ind2)];
    YM0(indd==1,ii) = nan; 
    
    fig = figure('Position',[20,20,900,600],'Name','','Color','w',...
      'Position',[1 scrsz(4)/13 scrsz(3)/4.5 scrsz(4)/3]);
    hold on
    for kk = 1:length(xx), bar(xx(kk),15000,ww(kk),'FaceColor',[0.8 0.8 0.8],'EdgeColor','w'), end
    for kk = 1:length(xx), bar(xx(kk),-15000,ww(kk),'FaceColor',[0.8 0.8 0.8],'EdgeColor','w'), end
    hold on
    plot(mtm,YM0(:,ii),'Color',[0.3 0.3 0.9],'linewidth',2)
    axis([mtm(1) mtm(end) min(YM0(:,ii)) max(YM0(:,ii))])
    ylabel('\%')
    print('-depsc2',[fullfile(folder_save, 'mdatadrop_'), num2str(ii), '.eps']);   
end






