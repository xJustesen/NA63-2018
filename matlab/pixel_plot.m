clc; clear all; close all

%% load data
datpath = '/home/christian/Dropbox/speciale/data/';
figpath = '/home/christian/Dropbox/speciale/figures/';

filepath = strcat(datpath,'pixeldata_run44.txt');
filepath2 = strcat(datpath,'pixeldata_run44_with_hot_pixels.txt');

formatSpec = '%f %f';
fileID = fopen(filepath,'r');
fileID2 = fopen(filepath2,'r');
datblocks = 6;
pixeldat = loaddat(fileID, formatSpec, datblocks);
pixeldat2 = loaddat(fileID2, formatSpec, datblocks);

%% plot data
for i = 1:datblocks
    field = strcat('plane',num2str(i-1));
    titlestr = strcat('Plane ',num2str(i-1),' hits');

    figure(i);
    hold on
    plot(1:length(pixeldat2.(field)),pixeldat2.(field));
    plot(1:length(pixeldat.(field)),pixeldat.(field));
    xlabel('Pixel');ylabel('Total hits');
    xlim([0, 576*1152])
    legend('Raw data','Hot pixels removed')
    ax = gca; ax.YScale = 'log';
    grid on
    box on

end

%% functions
function hitdat = loaddat(fileID, formatSpec, datblocks)
    for i = 1:datblocks
        field = strcat('plane',num2str(i-1));
        data = textscan(fileID,formatSpec,'HeaderLines',1);
        hitdat.(field) = data{1};
    end
end
