clc;clear  

[YM0, YMX]  = xlsread('YM_Jan2022.xls');
[YQ0, YQX]  = xlsread('YQ_Jan2022.xls');


yname = {'Log value','Percent'};

mdate       = YMX(2:end,1); 
qdate       = YQX(2:end,1); 

mdate0      = 1964+1/12:1/12:2022+1/12;
qdate0      = 1964+.25:.25:2022;


select      = [1 0 0 0 0 1 1 0];
%YM0(:,select==1) = YM0(:,select==1)./100;
YM0(:,select==0) = log(YM0(:,select==0));
YQ0         = log(YQ0(:,:));
YM          = YM0;
YQ          = kron(YQ0,ones(3,1));

% allowing for possible release date mismatches within QF variables
% Tstar is the length of all available data (at least includes a non-NaN variable)
% T is the length of the balanced panel set
Tstar       = size(YM,1); 
T           = find((sum(isnan(YQ),2))==0, 1, 'last' );    
Torigin     = T;

% collect data
YDATA       = nan(Tstar,11);                      % 11 variables
YDATA(:,1:8)= YM; 
YDATA(1:length(YQ),9:end) = YQ;

gDATA = [nan(1,11);1200*diff(YDATA)];
gQ0   = [nan(1,3);400*diff(YQ0)];


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

for ii=1:8
    fig = figure('Position',[20,20,900,600],'Name','','Color','w',...
      'Position',[1 scrsz(4)/13 scrsz(3)/4.5 scrsz(4)/3]);
    hold on
    for kk = 1:length(xx), bar(xx(kk),15000,ww(kk),'FaceColor',[0.8 0.8 0.8],'EdgeColor','w'), end
    for kk = 1:length(xx), bar(xx(kk),-15000,ww(kk),'FaceColor',[0.8 0.8 0.8],'EdgeColor','w'), end
    hold on
    plot(mdate0,YDATA(:,ii),'Color',[0.3 0.3 0.9],'linewidth',2)
    axis([mdate0(1) mdate0(end) min(YDATA(:,ii)) max(YDATA(:,ii))])
    ylabel(yname(1+select(ii)))
    print('-depsc2',[fullfile(folder_save, 'mdata_'), num2str(ii), '.eps']);
end

for ii=9:11
    jj = ii-8;
    fig = figure('Position',[20,20,900,600],'Name','','Color','w',...
      'Position',[1 scrsz(4)/13 scrsz(3)/4.5 scrsz(4)/3]);
    hold on
    for kk = 1:length(xx), bar(xx(kk),15000,ww(kk),'FaceColor',[0.8 0.8 0.8],'EdgeColor','w'), end
    for kk = 1:length(xx), bar(xx(kk),-15000,ww(kk),'FaceColor',[0.8 0.8 0.8],'EdgeColor','w'), end
    hold on
    plot(qdate0,YQ0(:,jj),'Color',[0.3 0.3 0.9],'linewidth',2)
    axis([qdate0(1) mdate0(end) min(YQ0(:,jj)) max(YQ0(:,jj))])
    ylabel('Log value')
    print('-depsc2',[fullfile(folder_save, 'qdata_'), num2str(ii), '.eps']);
end


for ii=1:8
    fig = figure('Position',[20,20,900,600],'Name','','Color','w',...
      'Position',[1 scrsz(4)/13 scrsz(3)/4.5 scrsz(4)/3]);
    hold on
    for kk = 1:length(xx), bar(xx(kk),15000,ww(kk),'FaceColor',[0.8 0.8 0.8],'EdgeColor','w'), end
    for kk = 1:length(xx), bar(xx(kk),-15000,ww(kk),'FaceColor',[0.8 0.8 0.8],'EdgeColor','w'), end
    hold on
    plot(mdate0,gDATA(:,ii),'Color',[0.3 0.3 0.9],'linewidth',1.5)
    axis([mdate0(1) mdate0(end) min(gDATA(:,ii)) max(gDATA(:,ii))])
    ylabel('Percent')
    print('-depsc2',[fullfile(folder_save, 'mdatagr_'), num2str(ii), '.eps']);
end

for ii=9:11
    jj = ii-8;
    fig = figure('Position',[20,20,900,600],'Name','','Color','w',...
      'Position',[1 scrsz(4)/13 scrsz(3)/4.5 scrsz(4)/3]);
    hold on
    for kk = 1:length(xx), bar(xx(kk),15000,ww(kk),'FaceColor',[0.8 0.8 0.8],'EdgeColor','w'), end
    for kk = 1:length(xx), bar(xx(kk),-15000,ww(kk),'FaceColor',[0.8 0.8 0.8],'EdgeColor','w'), end
    hold on
    plot(qdate0,gQ0(:,jj),'Color',[0.3 0.3 0.8],'linewidth',1.5)
    axis([qdate0(1) mdate0(end) min(gQ0(:,jj)) max(gQ0(:,jj))])
    ylabel('Percent')
    print('-depsc2',[fullfile(folder_save, 'qdatagr_'), num2str(ii), '.eps']);
end






for ii=1:8
    fig = figure('Position',[20,20,900,600],'Name','','Color','w',...
      'Position',[1 scrsz(4)/13 scrsz(3)/4.5 scrsz(4)/3]);
    hold on
    for kk = 1:length(xx), bar(xx(kk),15000,ww(kk),'FaceColor',[0.8 0.8 0.8],'EdgeColor','w'), end
    for kk = 1:length(xx), bar(xx(kk),-15000,ww(kk),'FaceColor',[0.8 0.8 0.8],'EdgeColor','w'), end
    hold on
    plot(mdate0,YDATA(:,ii),'Color',[0.3 0.3 0.9],'linewidth',2)
    axis([mdate0(432) mdate0(end) min(YDATA(:,ii)) max(YDATA(:,ii))])
    ylabel(yname(1+select(ii)))
    print('-depsc2',[fullfile(folder_save, 'mdata_post2000_'), num2str(ii), '.eps']);
end

for ii=9:11
    jj = ii-8;
    fig = figure('Position',[20,20,900,600],'Name','','Color','w',...
      'Position',[1 scrsz(4)/13 scrsz(3)/4.5 scrsz(4)/3]);
    hold on
    for kk = 1:length(xx), bar(xx(kk),15000,ww(kk),'FaceColor',[0.8 0.8 0.8],'EdgeColor','w'), end
    for kk = 1:length(xx), bar(xx(kk),-15000,ww(kk),'FaceColor',[0.8 0.8 0.8],'EdgeColor','w'), end
    hold on
    plot(qdate0,YQ0(:,jj),'Color',[0.3 0.3 0.9],'linewidth',2)
    axis([qdate0(144) mdate0(end) min(YQ0(:,jj)) max(YQ0(:,jj))])
    ylabel('Log value')
    print('-depsc2',[fullfile(folder_save, 'qdata_post2000_'), num2str(ii), '.eps']);
end


for ii=1:8
    fig = figure('Position',[20,20,900,600],'Name','','Color','w',...
      'Position',[1 scrsz(4)/13 scrsz(3)/4.5 scrsz(4)/3]);
    hold on
    for kk = 1:length(xx), bar(xx(kk),15000,ww(kk),'FaceColor',[0.8 0.8 0.8],'EdgeColor','w'), end
    for kk = 1:length(xx), bar(xx(kk),-15000,ww(kk),'FaceColor',[0.8 0.8 0.8],'EdgeColor','w'), end
    hold on
    plot(mdate0,gDATA(:,ii),'Color',[0.3 0.3 0.9],'linewidth',1.5)
    axis([mdate0(432) mdate0(end) min(gDATA(:,ii)) max(gDATA(:,ii))])
    ylabel('Percent')
    print('-depsc2',[fullfile(folder_save, 'mdatagr_post2000_'), num2str(ii), '.eps']);
end

for ii=9:11
    jj = ii-8;
    fig = figure('Position',[20,20,900,600],'Name','','Color','w',...
      'Position',[1 scrsz(4)/13 scrsz(3)/4.5 scrsz(4)/3]);
    hold on
    for kk = 1:length(xx), bar(xx(kk),15000,ww(kk),'FaceColor',[0.8 0.8 0.8],'EdgeColor','w'), end
    for kk = 1:length(xx), bar(xx(kk),-15000,ww(kk),'FaceColor',[0.8 0.8 0.8],'EdgeColor','w'), end
    hold on
    plot(qdate0,gQ0(:,jj),'Color',[0.3 0.3 0.9],'linewidth',1.5)
    axis([qdate0(144) mdate0(end) min(gQ0(:,jj)) max(gQ0(:,jj))])
    ylabel('Percent')
    print('-depsc2',[fullfile(folder_save, 'qdatagr_post2000_'), num2str(ii), '.eps']);
end





