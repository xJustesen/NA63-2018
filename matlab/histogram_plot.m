clc; clear all; close all;
%% load data
datpath = '/home/christian/Dropbox/speciale/data/';
figpath = '/home/christian/Dropbox/speciale/figures/Align/';
filepath = strcat(datpath,'hits_coord_data_run44.txt');
formatSpec = '%f %f';
fileID = fopen(filepath,'r');
datblocks = 6;

hitdat = loaddat(fileID, formatSpec, datblocks);

%% Analyse data
binwidth = 25;

xhist = makehist(hitdat, binwidth, datblocks);
yhist = makehist(hitdat, binwidth, datblocks);

for i = 1:datblocks
    field = strcat('plane',num2str(i-1));  
    hitssorted = sort(hitdat.(field));
    x = linspace(hitssorted(1, 1), hitssorted(end, 1), length(xhist.(field)));
    y = linspace(hitssorted(1, 2), hitssorted(end, 2), length(yhist.(field)));
    
    disp(field)
%     fitx.(field) = fitgaussian(x.', xhist.(field));
%     fity.(field) = fitgaussian(y.', yhist.(field));
    
    figx = figure(1);
    figx.PaperUnits = 'centimeters';
    figx.PaperPosition = [0 0 22 26];
    figx.PaperSize = [22 26];
    spx = subplot(datblocks,1,i);
    hold on
    plot(x, xhist.(field),'k')
%     plot(x, fitx.(field).full(x),'linewidth',2,'color','blue','linestyle','-')
    xlabel('xpos'); ylabel('Total hits');
    xlim([-1.5, 1.5]*10^4);
    
    figy = figure(2);
    figy.PaperUnits = 'centimeters';
    figy.PaperPosition = [0 0 22 26];
    figy.PaperSize = [22 26];
    spy = subplot(datblocks,1,i);
    hold on
    plot(y, yhist.(field),'k')
%     plot(y, fity.(field).full(y),'linewidth',2,'color','blue','linestyle','-')
    xlabel('ypos'); ylabel('Total hits');
    xlim([-6, 11]*1000)
end

%% Function definitions
function fitparams = fitgaussian(x, y)
    initialfit = fit(x, y, 'gauss1');
    fitparams.full = initialfit;
    fitparams.center = initialfit.a1;
    fitparams.width = initialfit.c1;
end

function hitdat = loaddat(fileID, formatSpec, datblocks)
    for i = 1:datblocks
        field = strcat('plane',num2str(i-1));
        data = textscan(fileID,formatSpec,'HeaderLines',1);
        hitdat.(field) = [data{1},data{2}];
    end
end

function counts = makehist(hitdat, binwidth, datblocks)
    for i = 1:datblocks
        field = strcat('plane',num2str(i-1));
        counts.(field) = transpose(histcounts(hitdat.(field)(:,1),'binwidth',binwidth));
    end
end