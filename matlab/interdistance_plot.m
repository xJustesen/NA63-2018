clear all; clc; close all;
%% load data
datpath = '/home/christian/Dropbox/speciale/data/';
filepath = strcat(datpath,'interdistance_data_7.txt');
filepath2 = strcat(datpath,'interdistance_data_7.txt');
figpath = '/home/christian/Dropbox/speciale/figures/';
formatSpec = '%f';
fileID = fopen(filepath,'r');
fileID2 = fopen(filepath2,'r');
datblocks = 4;
dat = loaddat(fileID, formatSpec, datblocks);
dat2 = loaddat(fileID2, formatSpec, datblocks);

%% analyse data
for i = 1:datblocks
    field = strcat('plane',num2str(i-1));
    [counts, edges] = histcounts(dat.(field),'binwidth',5);
    bincenters = (edges(2:end) - edges(1:end-1))/2 + edges(1:end-1);
    [counts2, edges2] = histcounts(dat2.(field),'binwidth',5);
    bincenters2 = (edges2(2:end) - edges2(1:end-1))/2 + edges2(1:end-1);
    fig = figure(i);
    ax = gca;
    ax.XScale = 'log';
    title(strcat(num2str(i-1), ' + ', num2str(i),' -> ',num2str(i+1)))
    hold on
    plot(bincenters,counts,'-')
    xlabel('Distance (\mum)'); ylabel('Occurance')
    grid on
end
%% functions
function hitdat = loaddat(fileID, formatSpec, datblocks)
    for i = 1:datblocks
        field = strcat('plane',num2str(i-1));
        data = textscan(fileID,formatSpec,'HeaderLines',1);
        hitdat.(field) = data{1};
    end
end