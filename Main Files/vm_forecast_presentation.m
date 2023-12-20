function output = vm_forecast_presentation(actual,forecast,t1,t2)

% if sum(isnan(actual(end,:)))==0
%     output = [actual(end-t1+1:end,:);forecast(1:t2,:)];
% else
%     output = [actual(end-t1+1:end-1,:);forecast(1:t2,:)];
% end


    output = [actual(end-t1+1:end,:);forecast(1:t2,:)];