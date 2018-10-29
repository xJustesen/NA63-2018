clear all; clc; close all
%% load data
filename1 = 'hits_coord_data_4'; % fake hits
filename2 = 'hits_coord_data_49';
filename3 = 'hits_coord_data_sim_alignment';

fakehitdat  = loaddat(filename1);
hitdat      = loaddat(filename2);
simhitdat      = loaddat(filename3);

fakeevents = 130042;
events = 928547;
simevents = 1e6;

%% plot data
for i = 1:numel(fieldnames(hitdat))
    field = strcat('plane',num2str(i-1));
    disp(field);
    disp('  ')
    disp(['Hits per event = ', num2str(length(hitdat.(field))/events)]);
    disp(['Sim hits per event = ', num2str(length(simhitdat.(field))/simevents)]);
    disp('  ')
    
    titlestr = strcat('Plane ',num2str(i-1),' hits');
       
    [xcounts, xbins] = hist(hitdat.(field)(:,1),100);
    [ycounts, ybins] = hist(hitdat.(field)(:,2),100);
    
    [fakexcounts, fakexbins] = hist(fakehitdat.(field)(:,1),100);
    [fakeycounts, fakeybins] = hist(fakehitdat.(field)(:,2),100);

    netcountsx = xcounts/events - fakexcounts/fakeevents;
    netcountsy = ycounts/events - fakeycounts/fakeevents;
    
    netcountsx = netcountsx(netcountsx > 0);
    netcountsy = netcountsy(netcountsy > 0);
    netxbins = xbins(netcountsx > 0);
    netybins = ybins(netcountsy > 0);
    
    if i == 1
                
        f = figure;
        subplot(1,2,1)
        hold on
        plot(xbins, xcounts/events,'linewidth',1.5)
        plot(ybins, ycounts/events,'linewidth',1.5)
        box on
        grid on
        ylim([0, 0.04]);
        xlim([-1.5e4, 1.5e4]);
        legend({'x','y'},'interpreter','latex');
        set(gca, 'FontSize', 18)
        xlabel('position [$\mu$m]','fontsize',22,'interpreter','latex'); ylabel('Normalized counts','fontsize',22,'interpreter','latex');
        title('a) incl fakes','fontsize',22,'interpreter','latex')
%         set(f, 'Units','centimeters','PaperUnits','centimeters', 'PaperSize',[18, 12],'PaperPosition',[0, 0, 18, 12],'Position',[0 0 18 12])
%         print('../../figures/beamspatial_80GeV_1_5mm_fakes', '-dpdf','-r600','-painters')

        subplot(1,2,2)
%         f = figure;
        hold on
        plot(netxbins, netcountsx,'-','linewidth',1.5)
        plot(netybins, netcountsy,'-','linewidth',1.5)
        box on
        grid on       
        ylim([0, 0.04]);
        xlim([-1.5e4, 1.5e4]);
        legend({'x','y'},'interpreter','latex');
        set(gca, 'FontSize', 18)
        xlabel('position [$\mu$m]','fontsize',22,'interpreter','latex'); ylabel('Normalized counts','fontsize',22,'interpreter','latex');
        title('b) excl fakes','fontsize',22,'interpreter','latex');
        
        set(f,'Units','centimeters','PaperUnits','centimeters', 'PaperSize',[36, 12],'PaperPosition',[0, 0, 36, 12],'Position',[0 0 36 12])
        print(f, '../../figures/beamspatial_80GeV_1_5mm', '-dpdf','-r600','-painters')

    end
%     
    if i == 1
        netcountsx = netcountsx.';
        netcountsy = netcountsy.';
        xbins = xbins.';
        ybins = ybins.';
%         save('../beamParameters/xdat_alignment_beam_params.txt','xbins','-ascii');
%         save('../beamParameters/ydat_alignment_beam_params.txt','ybins','-ascii');
%         save('../beamParameters/xweight_alignment_beam_params.txt','netcountsx','-ascii');
%         save('../beamParameters/yweight_alignment_beam_params.txt','netcountsy','-ascii');
    end
end

%% functions
function hitdat = loaddat(filename)
    datpath = '/home/christian/Dropbox/speciale/data/';
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
