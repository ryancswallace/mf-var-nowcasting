clc;clear;

path('Data',path);
path('SPF',path);
path('Main Files',path);
path('Moment',path);
path('Gauss Files',path);
path('Output',path);
path('Figures',path);

folder_save = 'C:\Users\dsong14\Dropbox\Projects\MFVAR_resuscitation\MFVAR-2023\Figures';

%% read MF-VAR forecasts
gr50vec0 = []; lv50vec0 = []; cgr50vec0 = [];
gr05vec0 = []; lv05vec0 = []; cgr05vec0 = [];
gr95vec0 = []; lv95vec0 = []; cgr95vec0 = [];

gr50vec1 = []; lv50vec1 = []; cgr50vec1 = [];
gr05vec1 = []; lv05vec1 = []; cgr05vec1 = [];
gr95vec1 = []; lv95vec1 = []; cgr95vec1 = [];

gr50vec2 = []; lv50vec2 = []; cgr50vec2 = [];
gr05vec2 = []; lv05vec2 = []; cgr05vec2 = [];
gr95vec2 = []; lv95vec2 = []; cgr95vec2 = [];

gr50vec3 = []; lv50vec3 = []; cgr50vec3 = [];
gr05vec3 = []; lv05vec3 = []; cgr05vec3 = [];
gr95vec3 = []; lv95vec3 = []; cgr95vec3 = [];


for dselect = 7:7
    % 1: full estimation for each vintage
    fixpara = 0;
    [gr50,gr05,gr95,lv50,lv05,lv95,cgr50,cgr05,cgr95]...
        = compile(['Output\GR_forecast_',  num2str(dselect),'_', num2str(1), '.xls'],...
        ['Output\LV_forecast_',  num2str(dselect),'_', num2str(1), '.xls'],...
        ['Output\CGR_forecast_', num2str(dselect),'_', num2str(1), '.xls'],3,4,7);

    gr50vec1{dselect}  = gr50;  gr05vec1{dselect}  = gr05;  gr95vec1{dselect}  = gr95;
    cgr50vec1{dselect} = cgr50; cgr05vec1{dselect} = cgr05; cgr95vec1{dselect} = cgr95;
    lv50vec1{dselect}  = lv50;  lv05vec1{dselect}  = lv05;  lv95vec1{dselect}  = lv95;

    % 2: skip covid sample in the estimation: March, April, May, June 2020
    [gr50,gr05,gr95,lv50,lv05,lv95,cgr50,cgr05,cgr95]...
        = compile(['Output\GR_forecast_',  num2str(dselect),'_', num2str(2), '.xls'],...
        ['Output\LV_forecast_',  num2str(dselect),'_', num2str(2), '.xls'],...
        ['Output\CGR_forecast_', num2str(dselect),'_', num2str(2), '.xls'],3,4,7);

    gr50vec2{dselect}  = gr50;  gr05vec2{dselect}  = gr05;  gr95vec2{dselect}  = gr95;
    cgr50vec2{dselect} = cgr50; cgr05vec2{dselect} = cgr05; cgr95vec2{dselect} = cgr95;
    lv50vec2{dselect}  = lv50;  lv05vec2{dselect}  = lv05;  lv95vec2{dselect}  = lv95;

    % 3: Lenza and Primiceri method
    [gr50,gr05,gr95,lv50,lv05,lv95,cgr50,cgr05,cgr95]...
        = compile(['Output\GR_forecast_',  num2str(dselect),'_', num2str(3), '.xls'],...
        ['Output\LV_forecast_',  num2str(dselect),'_', num2str(3), '.xls'],...
        ['Output\CGR_forecast_', num2str(dselect),'_', num2str(3), '.xls'],3,4,7);

    gr50vec3{dselect}  = gr50;  gr05vec3{dselect}  = gr05;  gr95vec3{dselect}  = gr95;
    cgr50vec3{dselect} = cgr50; cgr05vec3{dselect} = cgr05; cgr95vec3{dselect} = cgr95;
    lv50vec3{dselect}  = lv50;  lv05vec3{dselect}  = lv05;  lv95vec3{dselect}  = lv95;
end


%% load SPF forecasts
load('Output\SPFforecast.mat')

%% note
%             1    2    3   4  5   6   7   8    9    10    11     
% variables: UNR,HOURS,CPI,IP,PCE,FFR,TB,SP500,GDP,FIXINV,GOV

time           = 2019.75:.25:2021.75;
nv             = 11;
latest_vintage = 7; % August 2021 vintage which has data upto 2021:q2

% latest data Jan 2020 vintage
load forecast_realdata

scrsz = get(0,'ScreenSize');
set(0,'defaulttextinterpreter','latex')
set(0,'defaultaxesfontname','times')  
set(0,'DefaultAxesFontSize',20)
set(0,'DefaultTextFontSize',20);

for vintage = 7:latest_vintage
    SUB_COMPARE_FIGURE
    close all
end
















function [gr50,gr05,gr95,lv50,lv05,lv95,cgr50,cgr05,cgr95] = compile(inputGR,inputLV,inputCGR,ind50,ind05,ind95)
%             1    2    3   4  5   6   7   8    9    10    11     
% variables: UNR,HOURS,CPI,IP,PCE,FFR,TB,SP500,GDP,FIXINV,GOV

% what we report in the paper
% %              1 2 3 4 5 6 7  8  9 10 11 
% indLV  = [[0 0 1 0 0 0 0 1 1  0  0  0  0] zeros(1,11)];
% indGR  = [[0 0 0 0 1 0 0 0 0  0  0  0  0] zeros(1,11)];
% indCGR = [[0 0 0 1 0 1 1 0 0  1  1  1  1] zeros(1,11)];

ind = [0 0 ones(1,11) zeros(1,11)];

[gr50,etc] = xlsread(inputGR,ind50); gr50 = gr50(:,ind==1);
[gr05,~]   = xlsread(inputGR,ind05); gr05 = gr05(:,ind==1);
[gr95,~]   = xlsread(inputGR,ind95); gr95 = gr95(:,ind==1);

[lv50,~]   = xlsread(inputLV,ind50); lv50 = lv50(:,ind==1);
[lv05,~]   = xlsread(inputLV,ind05); lv05 = lv05(:,ind==1);
[lv95,~]   = xlsread(inputLV,ind95); lv95 = lv95(:,ind==1);

[cgr50,~]  = xlsread(inputCGR,ind50);cgr50 = cgr50(:,ind==1);
[cgr05,~]  = xlsread(inputCGR,ind05);cgr05 = cgr05(:,ind==1);
[cgr95,~]  = xlsread(inputCGR,ind95);cgr95 = cgr95(:,ind==1);

end

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