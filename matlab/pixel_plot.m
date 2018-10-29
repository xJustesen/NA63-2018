clc; clear all; close all

%% load data
datpath = '/home/christian/Dropbox/speciale/data/';
figpath = '/home/christian/Dropbox/speciale/figures/';

filepath = strcat(datpath,'pixeldata_run54.txt');
filepath2 = strcat(datpath,'pixeldata_run54_incl_hotpix.txt');

formatSpec = '%f %f';
fileID = fopen(filepath,'r');
fileID2 = fopen(filepath2,'r');
datblocks = 6;
pixeldat = loaddat(fileID, formatSpec, datblocks);
pixeldat2 = loaddat(fileID2, formatSpec, datblocks);

%% plot data
for i = 1:1
    field = strcat('plane',num2str(i-1));
    titlestr = strcat('Plane ',num2str(i-1),' hits');

    f= figure(i);
    hold on
    plot(1:length(pixeldat2.(field)),pixeldat2.(field),'linewidth',2.5);
    plot(1:length(pixeldat.(field)),pixeldat.(field),'linewidth',2.5);
    xlabel('Pixel','fontsize',22,'interpreter','latex');ylabel('Total hits','fontsize',22,'interpreter','latex');
    xlim([0, 576*1152])
    ylim([0, 1e6])
    legend({'Including hotpixels','Excluding hotpixels'},'fontsize',22,'interpreter','latex');
    ax = gca; ax.YScale = 'log';
    set(gca, 'FontSize', 14)
    grid on
    box on
    set(f, 'Units','centimeters','PaperUnits','centimeters', 'PaperSize',[18, 12],'PaperPosition',[0, 0, 18, 12],'Position',[0 0 18 12])
    print('../../figures/hotpix_run54.pdf', '-dpdf','-r600','-painters')


end

%% functions
function hitdat = loaddat(fileID, formatSpec, datblocks)
    for i = 1:datblocks
        field = strcat('plane',num2str(i-1));
        data = textscan(fileID,formatSpec,'HeaderLines',1);
        hitdat.(field) = data{1};
    end
end
