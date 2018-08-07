clear all; clc; close all
%% load data
datpath = '/home/christian/Dropbox/speciale/data/';
figpath = '/home/christian/Dropbox/speciale/figures/';
filename = 'hits_coord_data_sim_background1_80GeV_1.5mm';
filepath = strcat(datpath,filename,'.txt');
formatSpec = '%f %f';
fileID = fopen(filepath,'r');
datblocks = 6;

hitdat = loaddat(fileID, formatSpec, datblocks);

%% plot data
for i = 1:datblocks
    field = strcat('plane',num2str(i-1));
    disp(field);
    titlestr = strcat('Plane ',num2str(i-1),' hits');    
    
    figure(i)
    hold on
    title(titlestr);
    plot(hitdat.(field)(:,1),hitdat.(field)(:,2),'.','markersize',.001);
    xlabel('xpos (\mum)'); ylabel('ypos (\mum)');
    xlim([-10000, 10000]);
    box on
    grid on
    axis equal
   
    [xcounts, xbins] = hist(hitdat.(field)(:,1),1000);
    [ycounts, ybins] = hist(hitdat.(field)(:,2),1000);
%     
%     if i == 1
%         xbins = xbins.';
%         ybins = ybins.';
%         xcounts = xcounts.';
%         ycounts = ycounts.';
%         save('../cpp/xdat_80GeV_beam_params.txt','xbins','-ascii');
%         save('../cpp/ydat_80GeV_beam_params.txt','ybins','-ascii');
%         save('../cpp/xweight_80GeV_beam_params.txt','xcounts','-ascii');
%         save('../cpp/yweight_80GeV_beam_params.txt','ycounts','-ascii');
%     end
% %     % 

    fitx = fit(xbins.', xcounts.', 'gauss1')
    fity = fit(ybins.', ycounts.', 'gauss1')
  
    x = linspace(min(xbins), max(xbins), 1000);
    % Find the half max value.
    halfMax = (min(fitx(x)) + max(fitx(x))) / 2;
    % Find where the data first drops below half the max.
    xl = x(find(fitx(x) >= halfMax, 1, 'first'));
    xr = x(find(fitx(x) >= halfMax, 1, 'last'));
    fwhmx = xr - xl;
    
    y = linspace(min(ybins), max(ybins), 1000);
    halfMax = (min(fity(y)) + max(fity(y))) / 2;
    % Find where the data first drops below half the max.
    yl = y(find(fity(y) >= halfMax, 1, 'first'));
    yr = y(find(fity(y) >= halfMax, 1, 'last'));
    fwhmy = yr - yl;
    disp('')
    disp(['xmean = ',num2str(x(fitx(x) == max(fitx(x)))),' ; fwhmx = ', num2str(fwhmx), ' ; stddev = ', num2str(fwhmx/(2*sqrt(2*log(2))))]);
    disp(['ymean = ',num2str(y(fity(y) == max(fity(y)))),' ; fwhmy = ', num2str(fwhmy), ' ; stddev = ', num2str(fwhmy/(2*sqrt(2*log(2))))]);
    
    f = figure(i*10);
    hold on
    grid on
    box on
    title('Spatial distrobution ; 20GeV')
    xlabel('Position (\mu m)')
    ylabel('Counts')
    stairs(xbins, xcounts);
    stairs(ybins, ycounts);
%     plot(fitx);
%     plot(fity);
%     plot(xl,fitx(xl),'ro')
%     plot(xr,fitx(xr),'ro')
%     plot(yl,fity(yl),'ro')
%     plot(yr,fity(yr),'ro')
%     plot(x(fitx(x) == max(fitx(x))), max(fitx(x)),'bo')
%     plot(y(fity(y) == max(fity(y))), max(fity(y)),'bo')
    legend('x','y')
%     set(f, 'Units','centimeters','PaperUnits','centimeters', 'PaperSize',[18, 12],'PaperPosition',[0, 0, 18, 12],'Position',[0 0 18 12])
%     if i == 1
%         print('../../figures/spatial_distro_sim_amorphous1_80GeV', '-dpdf','-r600','-painters')
%     end
end

%% functions
function hitdat = loaddat(fileID, formatSpec, datblocks)
    for i = 1:datblocks
        field = strcat('plane',num2str(i-1));
        data = textscan(fileID,formatSpec,'HeaderLines',1);
        hitdat.(field) = [data{1},data{2}];
%         disp([min(data{1}), min(data{2}(data{1} == max(data{1})));
%               max(data{1}), max(data{2}(data{1} == min(data{1})));
%               max(data{1}(data{2} == max(data{2}))), max(data{2});
%               min(data{1}(data{2} == min(data{2}))), min(data{2})])
    end
end
