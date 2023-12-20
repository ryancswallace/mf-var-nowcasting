%% level series
for ii=1:nv

    % 1: full estimation for each vintage
    fig = figure('Position',[20,20,900,600],'Name','','Color','w',...
      'Position',[1 scrsz(4)/13 scrsz(3)/4.5 scrsz(4)/3]);
    for jj=1:length(time)
        hold on
        line([time(jj) time(jj)],[-5000 5000],'color',[0.9 0.9 0.9],'linestyle','-')
    end
    hold on
    shadedplot(time,lv95vec1{vintage}(:,ii)',lv05vec1{vintage}(:,ii)',[0.9 0.9 0.9],[0.6 0.6 0.6]);
    hold on
    k1=plot(time(1:end),lv_data(1:end,ii),'Color',[0.9 0.3 0.3],'marker','o','markerfacecolor',[0.9 0.3 0.3],'markeredgecolor',[0.9 0.3 0.3],'markersize',10,'linewidth',2);
    hold on
    plot(time,lv50vec1{vintage}(:,ii),'Color',[0.1 0.1 0.1],'linewidth',2);
    plot(time,lv_spf{vintage}(:,ii),'s','markersize',10,'MarkerEdgeColor',[0.1 0.6 0.1],'MarkerFaceColor',[0.1 0.6 0.1]);
    xticks([time(1) time(2) time(3) time(4) time(5) time(6) time(7) time(8) time(9)])
    xticklabels({'19:Q4','','20:Q2','','20:Q4','','21:Q2','','21:Q4'})
    axis([time(1) time(end) min([min(lv05vec3{vintage}(:,ii)) min(lv05vec1{vintage}(:,ii))  min(lv05vec2{vintage}(:,ii)) min(lv05vec2{latest_vintage}(1:end-2,ii))  min(lv_spf{vintage}(:,ii))])...
                            min(max([max(lv95vec3{vintage}(:,ii)) max(lv95vec1{vintage}(:,ii))  max(lv95vec2{vintage}(:,ii)) max(lv95vec2{latest_vintage}(1:end-2,ii))  max(lv_spf{vintage}(:,ii))]),60)])
    if ii==1 || ii==6 || ii==7; ylabel('\%'); end
    print('-depsc2',[fullfile(folder_save, 'lv1_'), num2str(vintage),'_', num2str(ii), '.eps']);  

    % 2: skip covid sample in the estimation: March, April, May, June 2020
    fig = figure('Position',[20,20,900,600],'Name','','Color','w',...
      'Position',[1 scrsz(4)/13 scrsz(3)/4.5 scrsz(4)/3]);
    for jj=1:length(time)
        hold on
        line([time(jj) time(jj)],[-5000 5000],'color',[0.9 0.9 0.9],'linestyle','-')
    end
    hold on
    shadedplot(time,lv95vec2{vintage}(:,ii)',lv05vec2{vintage}(:,ii)',[0.9 0.9 0.9],[0.6 0.6 0.6]);
    hold on
    k1=plot(time(1:end),lv_data(1:end,ii),'Color',[0.9 0.3 0.3],'marker','o','markerfacecolor',[0.9 0.3 0.3],'markeredgecolor',[0.9 0.3 0.3],'markersize',10,'linewidth',2);
    hold on
    plot(time,lv50vec2{vintage}(:,ii),'Color',[0.1 0.1 0.1],'linewidth',2);
    plot(time,lv_spf{vintage}(:,ii),'s','markersize',10,'MarkerEdgeColor',[0.1 0.6 0.1],'MarkerFaceColor',[0.1 0.6 0.1]);
    xticks([time(1) time(2) time(3) time(4) time(5) time(6) time(7) time(8) time(9)])
    xticklabels({'19:Q4','','20:Q2','','20:Q4','','21:Q2','','21:Q4'})
    axis([time(1) time(end) min([min(lv05vec3{vintage}(:,ii)) min(lv05vec1{vintage}(:,ii))  min(lv05vec2{vintage}(:,ii)) min(lv05vec2{latest_vintage}(1:end-2,ii))  min(lv_spf{vintage}(:,ii))])...
                            min(max([max(lv95vec3{vintage}(:,ii)) max(lv95vec1{vintage}(:,ii))  max(lv95vec2{vintage}(:,ii)) max(lv95vec2{latest_vintage}(1:end-2,ii))  max(lv_spf{vintage}(:,ii))]),60)])
    if ii==1 || ii==6 || ii==7; ylabel('\%'); end
    print('-depsc2',[fullfile(folder_save, 'lv2_'), num2str(vintage),'_', num2str(ii), '.eps']);   
    
    % 3: Lenza and Primiceri
    fig = figure('Position',[20,20,900,600],'Name','','Color','w',...
      'Position',[1 scrsz(4)/13 scrsz(3)/4.5 scrsz(4)/3]);
    for jj=1:length(time)
        hold on
        line([time(jj) time(jj)],[-5000 5000],'color',[0.9 0.9 0.9],'linestyle','-')
    end
    hold on
    shadedplot(time,lv95vec3{vintage}(:,ii)',lv05vec3{vintage}(:,ii)',[0.9 0.9 0.9],[0.6 0.6 0.6]);
    hold on
    k1=plot(time(1:end),lv_data(1:end,ii),'Color',[0.9 0.3 0.3],'marker','o','markerfacecolor',[0.9 0.3 0.3],'markeredgecolor',[0.9 0.3 0.3],'markersize',10,'linewidth',2);
    hold on
    plot(time,lv50vec3{vintage}(:,ii),'Color',[0.1 0.1 0.1],'linewidth',2);
    plot(time,lv_spf{vintage}(:,ii),'s','markersize',10,'MarkerEdgeColor',[0.1 0.6 0.1],'MarkerFaceColor',[0.1 0.6 0.1]);
    xticks([time(1) time(2) time(3) time(4) time(5) time(6) time(7) time(8) time(9)])
    xticklabels({'19:Q4','','20:Q2','','20:Q4','','21:Q2','','21:Q4'})
    axis([time(1) time(end) min([min(lv05vec3{vintage}(:,ii)) min(lv05vec1{vintage}(:,ii))  min(lv05vec2{vintage}(:,ii)) min(lv05vec2{latest_vintage}(1:end-2,ii))  min(lv_spf{vintage}(:,ii))])...
                            min(max([max(lv95vec3{vintage}(:,ii)) max(lv95vec1{vintage}(:,ii))  max(lv95vec2{vintage}(:,ii)) max(lv95vec2{latest_vintage}(1:end-2,ii))  max(lv_spf{vintage}(:,ii))]),60)])
    if ii==1 || ii==6 || ii==7; ylabel('\%'); end
    print('-depsc2',[fullfile(folder_save, 'lv3_'), num2str(vintage),'_', num2str(ii), '.eps']);      
end

%% growth rates 
for ii=1:nv

    % 1: full estimation for each vintage
    fig = figure('Position',[20,20,900,600],'Name','','Color','w',...
      'Position',[1 scrsz(4)/13 scrsz(3)/4.5 scrsz(4)/3]);
    for jj=1:length(time)
        hold on
        line([time(jj) time(jj)],[-5000 5000],'color',[0.9 0.9 0.9],'linestyle','-')
    end
    hold on
    shadedplot(time,gr95vec1{vintage}(:,ii)',gr05vec1{vintage}(:,ii)',[0.9 0.9 0.9],[0.6 0.6 0.6]);
    hold on
    k1=plot(time(1:end),gr_data(1:end,ii),'Color',[0.9 0.3 0.3],'marker','o','markerfacecolor',[0.9 0.3 0.3],'markeredgecolor',[0.9 0.3 0.3],'markersize',10,'linewidth',2);
    hold on
    plot(time,gr50vec1{vintage}(:,ii),'Color',[0.1 0.1 0.1],'linewidth',2);
    plot(time,gr_spf{vintage}(:,ii),'s','markersize',10,'MarkerEdgeColor',[0.1 0.6 0.1],'MarkerFaceColor',[0.1 0.6 0.1]);
    xticks([time(1) time(2) time(3) time(4) time(5) time(6) time(7) time(8) time(9)])
    xticklabels({'19:Q4','','20:Q2','','20:Q4','','21:Q2','','21:Q4'})
    axis([time(1) time(end) max(min([min(gr05vec3{vintage}(:,ii)) min(gr05vec1{vintage}(:,ii))  min(gr05vec2{vintage}(:,ii)) min(gr05vec2{latest_vintage}(1:end-2,ii)) min(gr_spf{vintage}(:,ii))]),-100)...
                            max([max(gr95vec3{vintage}(:,ii)) max(gr95vec1{vintage}(:,ii))  max(gr95vec2{vintage}(:,ii)) max(gr95vec2{latest_vintage}(1:end-2,ii)) max(gr_spf{vintage}(:,ii))]) ]);
    ylabel('\%');
    print('-depsc2',[fullfile(folder_save, 'gr1_'), num2str(vintage),'_', num2str(ii), '.eps']);  

    % 2: skip covid sample in the estimation: March, April, May, June 2020
    fig = figure('Position',[20,20,900,600],'Name','','Color','w',...
      'Position',[1 scrsz(4)/13 scrsz(3)/4.5 scrsz(4)/3]);
    for jj=1:length(time)
        hold on
        line([time(jj) time(jj)],[-5000 5000],'color',[0.9 0.9 0.9],'linestyle','-')
    end
    hold on
    shadedplot(time,gr95vec2{vintage}(:,ii)',gr05vec2{vintage}(:,ii)',[0.9 0.9 0.9],[0.6 0.6 0.6]);
    hold on
    k1=plot(time(1:end),gr_data(1:end,ii),'Color',[0.9 0.3 0.3],'marker','o','markerfacecolor',[0.9 0.3 0.3],'markeredgecolor',[0.9 0.3 0.3],'markersize',10,'linewidth',2);
    hold on
    plot(time,gr50vec2{vintage}(:,ii),'Color',[0.1 0.1 0.1],'linewidth',2);
    plot(time,gr_spf{vintage}(:,ii),'s','markersize',10,'MarkerEdgeColor',[0.1 0.6 0.1],'MarkerFaceColor',[0.1 0.6 0.1]);
    xticks([time(1) time(2) time(3) time(4) time(5) time(6) time(7) time(8) time(9)])
    xticklabels({'19:Q4','','20:Q2','','20:Q4','','21:Q2','','21:Q4'})
    axis([time(1) time(end) max(min([min(gr05vec3{vintage}(:,ii)) min(gr05vec1{vintage}(:,ii))  min(gr05vec2{vintage}(:,ii)) min(gr05vec2{latest_vintage}(1:end-2,ii)) min(gr_spf{vintage}(:,ii))]),-100)...
                            max([max(gr95vec3{vintage}(:,ii)) max(gr95vec1{vintage}(:,ii))  max(gr95vec2{vintage}(:,ii)) max(gr95vec2{latest_vintage}(1:end-2,ii)) max(gr_spf{vintage}(:,ii))]) ]);
    ylabel('\%');
    print('-depsc2',[fullfile(folder_save, 'gr2_'), num2str(vintage),'_', num2str(ii), '.eps']);    

    % 3: Lenza and Primiceri
    fig = figure('Position',[20,20,900,600],'Name','','Color','w',...
      'Position',[1 scrsz(4)/13 scrsz(3)/4.5 scrsz(4)/3]);
    for jj=1:length(time)
        hold on
        line([time(jj) time(jj)],[-5000 5000],'color',[0.9 0.9 0.9],'linestyle','-')
    end
    hold on
    shadedplot(time,gr95vec3{vintage}(:,ii)',gr05vec3{vintage}(:,ii)',[0.9 0.9 0.9],[0.6 0.6 0.6]);
    hold on
    k1=plot(time(1:end),gr_data(1:end,ii),'Color',[0.9 0.3 0.3],'marker','o','markerfacecolor',[0.9 0.3 0.3],'markeredgecolor',[0.9 0.3 0.3],'markersize',10,'linewidth',2);
    hold on
    plot(time,gr50vec3{vintage}(:,ii),'Color',[0.1 0.1 0.1],'linewidth',2);
    plot(time,gr_spf{vintage}(:,ii),'s','markersize',10,'MarkerEdgeColor',[0.1 0.6 0.1],'MarkerFaceColor',[0.1 0.6 0.1]);
    xticks([time(1) time(2) time(3) time(4) time(5) time(6) time(7) time(8) time(9)])
    xticklabels({'19:Q4','','20:Q2','','20:Q4','','21:Q2','','21:Q4'})
    axis([time(1) time(end) max(min([min(gr05vec3{vintage}(:,ii)) min(gr05vec1{vintage}(:,ii))  min(gr05vec2{vintage}(:,ii)) min(gr05vec2{latest_vintage}(1:end-2,ii)) min(gr_spf{vintage}(:,ii))]),-100)...
                            max([max(gr95vec3{vintage}(:,ii)) max(gr95vec1{vintage}(:,ii))  max(gr95vec2{vintage}(:,ii)) max(gr95vec2{latest_vintage}(1:end-2,ii)) max(gr_spf{vintage}(:,ii))]) ]);
    ylabel('\%');
    print('-depsc2',[fullfile(folder_save, 'gr3_'), num2str(vintage),'_', num2str(ii), '.eps']);  
end


%% cumulative growth rates 
for ii=1:nv


    % 1: full estimation for each vintage
    fig = figure('Position',[20,20,900,600],'Name','','Color','w',...
      'Position',[1 scrsz(4)/13 scrsz(3)/4.5 scrsz(4)/3]);
    for jj=1:length(time)
        hold on
        line([time(jj) time(jj)],[-5000 5000],'color',[0.9 0.9 0.9],'linestyle','-')
    end
    hold on
    shadedplot(time,cgr95vec1{vintage}(:,ii)',cgr05vec1{vintage}(:,ii)',[0.9 0.9 0.9],[0.6 0.6 0.6]);
    hold on
    k1=plot(time(1:end),cgr_data(1:end,ii),'Color',[0.9 0.3 0.3],'marker','o','markerfacecolor',[0.9 0.3 0.3],'markeredgecolor',[0.9 0.3 0.3],'markersize',10,'linewidth',2);
    hold on
    plot(time,cgr50vec1{vintage}(:,ii),'Color',[0.1 0.1 0.1],'linewidth',2);
    plot(time,cgr_spf{vintage}(:,ii),'s','markersize',10,'MarkerEdgeColor',[0.1 0.6 0.1],'MarkerFaceColor',[0.1 0.6 0.1]);
    xticks([time(1) time(2) time(3) time(4) time(5) time(6) time(7) time(8) time(9)])
    xticklabels({'19:Q4','','20:Q2','','20:Q4','','21:Q2','','21:Q4'})
    axis([time(1) time(end) max(min([min(cgr05vec3{vintage}(:,ii)) min(cgr05vec1{vintage}(:,ii))  min(cgr05vec2{vintage}(:,ii)) min(cgr05vec2{latest_vintage}(1:end-2,ii)) min(cgr_spf{vintage}(:,ii))]),-60)...
                            max([max(cgr95vec3{vintage}(:,ii)) max(cgr95vec1{vintage}(:,ii))  max(cgr95vec2{vintage}(:,ii)) max(cgr95vec2{latest_vintage}(1:end-2,ii)) max(cgr_spf{vintage}(:,ii))])])
    ylabel('\%');
    print('-depsc2',[fullfile(folder_save, 'cgr1_'), num2str(vintage),'_', num2str(ii), '.eps']);  

    % 2: skip covid sample in the estimation: March, April, May, June 2020
    fig = figure('Position',[20,20,900,600],'Name','','Color','w',...
      'Position',[1 scrsz(4)/13 scrsz(3)/4.5 scrsz(4)/3]);
    for jj=1:length(time)
        hold on
        line([time(jj) time(jj)],[-5000 5000],'color',[0.9 0.9 0.9],'linestyle','-')
    end
    hold on
    shadedplot(time,cgr95vec2{vintage}(:,ii)',cgr05vec2{vintage}(:,ii)',[0.9 0.9 0.9],[0.6 0.6 0.6]);
    hold on
    k1=plot(time(1:end),cgr_data(1:end,ii),'Color',[0.9 0.3 0.3],'marker','o','markerfacecolor',[0.9 0.3 0.3],'markeredgecolor',[0.9 0.3 0.3],'markersize',10,'linewidth',2);
    hold on
    plot(time,cgr50vec2{vintage}(:,ii),'Color',[0.1 0.1 0.1],'linewidth',2);
    plot(time,cgr_spf{vintage}(:,ii),'s','markersize',10,'MarkerEdgeColor',[0.1 0.6 0.1],'MarkerFaceColor',[0.1 0.6 0.1]);
    xticks([time(1) time(2) time(3) time(4) time(5) time(6) time(7) time(8) time(9)])
    xticklabels({'19:Q4','','20:Q2','','20:Q4','','21:Q2','','21:Q4'})
    axis([time(1) time(end) max(min([min(cgr05vec3{vintage}(:,ii)) min(cgr05vec1{vintage}(:,ii))  min(cgr05vec2{vintage}(:,ii)) min(cgr05vec2{latest_vintage}(1:end-2,ii)) min(cgr_spf{vintage}(:,ii))]),-60)...
                            max([max(cgr95vec3{vintage}(:,ii)) max(cgr95vec1{vintage}(:,ii))  max(cgr95vec2{vintage}(:,ii)) max(cgr95vec2{latest_vintage}(1:end-2,ii)) max(cgr_spf{vintage}(:,ii))])])
    ylabel('\%');
    print('-depsc2',[fullfile(folder_save, 'cgr2_'), num2str(vintage),'_', num2str(ii), '.eps']);   
    
    % 3: Lenza and Primiceri method
    fig = figure('Position',[20,20,900,600],'Name','','Color','w',...
      'Position',[1 scrsz(4)/13 scrsz(3)/4.5 scrsz(4)/3]);
    for jj=1:length(time)
        hold on
        line([time(jj) time(jj)],[-5000 5000],'color',[0.9 0.9 0.9],'linestyle','-')
    end
    hold on
    shadedplot(time,cgr95vec3{vintage}(:,ii)',cgr05vec3{vintage}(:,ii)',[0.9 0.9 0.9],[0.6 0.6 0.6]);
    hold on
    k1=plot(time(1:end),cgr_data(1:end,ii),'Color',[0.9 0.3 0.3],'marker','o','markerfacecolor',[0.9 0.3 0.3],'markeredgecolor',[0.9 0.3 0.3],'markersize',10,'linewidth',2);
    hold on
    plot(time,cgr50vec3{vintage}(:,ii),'Color',[0.1 0.1 0.1],'linewidth',2);
    plot(time,cgr_spf{vintage}(:,ii),'s','markersize',10,'MarkerEdgeColor',[0.1 0.6 0.1],'MarkerFaceColor',[0.1 0.6 0.1]);
    xticks([time(1) time(2) time(3) time(4) time(5) time(6) time(7) time(8) time(9)])
    xticklabels({'19:Q4','','20:Q2','','20:Q4','','21:Q2','','21:Q4'})
    axis([time(1) time(end) max(min([min(cgr05vec3{vintage}(:,ii)) min(cgr05vec1{vintage}(:,ii))  min(cgr05vec2{vintage}(:,ii)) min(cgr05vec2{latest_vintage}(1:end-2,ii)) min(cgr_spf{vintage}(:,ii))]),-60)...
                            max([max(cgr95vec3{vintage}(:,ii)) max(cgr95vec1{vintage}(:,ii))  max(cgr95vec2{vintage}(:,ii)) max(cgr95vec2{latest_vintage}(1:end-2,ii)) max(cgr_spf{vintage}(:,ii))])])
    ylabel('\%');
    print('-depsc2',[fullfile(folder_save, 'cgr3_'), num2str(vintage),'_', num2str(ii), '.eps']);   
end



