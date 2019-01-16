clear all; clc; close all
%% load data
filename = 'hits_coord_data_77';
filename1 = 'hits_coord_data_4'; % fake hits

hitdat   = loaddat(filename);
fakehitdat  = loaddat(filename1);
fakeevents = 130042
events = 201678

datpath = '/home/christian/Documents/cern2018/simdata/';

%% load data
input_pos   = load('/home/christian/Dropbox/Cern2018Experiment/kode/beam_position_4mm.txt');
input_ang   = load('/home/christian/Dropbox/Cern2018Experiment/kode/beam_direction_4mm.txt'); 

xpos = input_pos(:,3);
ypos = input_pos(:,4);
xang = input_ang(:,3);
yang = input_ang(:,4);

%% plot data
for i = 1:numel(fieldnames(hitdat))
    field = strcat('plane',num2str(i-1));
    disp(field);
    disp('  ')
    
    titlestr = strcat('Plane ',num2str(i-1),' hits');
       
    [xcounts, xbins] = hist(hitdat.(field)(:,1),1000);
    [ycounts, ybins] = hist(hitdat.(field)(:,2),1000);
    
    [fakexcounts, fakexbins] = hist(fakehitdat.(field)(:,1),1000);
    [fakeycounts, fakeybins] = hist(fakehitdat.(field)(:,2),1000);

    netcountsx = xcounts/events - fakexcounts/fakeevents;
    netcountsy = ycounts/events - fakeycounts/fakeevents;
    
    netcountsx = netcountsx(netcountsx > 0);
    netcountsy = netcountsy(netcountsy > 0);
    netxbins = xbins(netcountsx > 0);
    netybins = ybins(netcountsy > 0);

    f = figure;
    title(titlestr)
    plot(hitdat.(field)(:,1), hitdat.(field)(:,2));
    
    if i == 1
                
        f = figure;
        title('x')
        hold on
        plot(netxbins, netcountsx/max(netcountsx))
        plot(xpos, input_pos(:,1)/max(input_pos(:,1)))
        
        figure
        title('y')
        hold on
        plot(netybins, netcountsy/max(netcountsy))
        plot(ypos, input_pos(:,2)/max(input_pos(:,2)))
%         box on
%         grid on
%         legend({'x','y'},'interpreter','latex');
%         set(gca, 'FontSize', 24)
%         xlabel('position [$\mu$m]','fontsize',36,'interpreter','latex'); ylabel('Normalized counts','fontsize',36,'interpreter','latex');

    end
%     
    if i == 1
%         netcountsx = netcountsx.';
%         netcountsy = netcountsy.';
%         xbins = xbins.';
%         ybins = ybins.';
%         save('../beamParameters/xdat_alignment_beam_params.txt','xbins','-ascii');
%         save('../beamParameters/ydat_alignment_beam_params.txt','ybins','-ascii');
%         save('../beamParameters/xweight_alignment_beam_params.txt','netcountsx','-ascii');
%         save('../beamParameters/yweight_alignment_beam_params.txt','netcountsy','-ascii');
    end
end

%% functions
function hitdat = loaddat(filename)
    datpath = '/home/christian/Documents/cern2018/simdata/';
    filepath = strcat(datpath,filename,'.txt');
    formatSpec = '%f %f';
    fileID = fopen(filepath,'r');
    datblocks = 6;
    for i = 1:datblocks
        field = strcat('plane',num2str(i-1));
        data = textscan(fileID,formatSpec,'HeaderLines',1);
        hitdat.(field) = [data{1},data{2}];
    end
end
