clear all; clc; close all;
%% load data
datpath = '/home/christian/Dropbox/speciale/data/';
filepath = strcat(datpath,'interdistance_data_18.txt');
filepath1 = strcat(datpath,'interdistance_data_73_no_align.txt');
filepath2 = strcat(datpath,'interdistance_data_sim_alignment.txt');
figpath = '/home/christian/Dropbox/speciale/figures/';

formatSpec = '%f';
fileID = fopen(filepath,'r');
fileID1 = fopen(filepath1,'r');
fileID2 = fopen(filepath2,'r');

datblocks = 4;
dat = loaddat(fileID, formatSpec, datblocks);
dat1 = loaddat(fileID1, formatSpec, datblocks);
dat2 = loaddat(fileID2, formatSpec, datblocks);

colors =    [       0    0.4470    0.7410
                0.8500    0.3250    0.0980
                0.9290    0.6940    0.1250
                0.4940    0.1840    0.5560
                0.4660    0.6740    0.1880
                0.3010    0.7450    0.9330
                0.6350    0.0780    0.1840
            ];
    indx = 1;
%% analyse data
for i = 1:datblocks
    field = strcat('plane',num2str(i-1));
    [counts, edges] = histcounts(dat.(field),'binwidth',1);
    [counts1, edges1] = histcounts(dat1.(field),'binwidth',10);
    [counts2, edges2] = histcounts(dat2.(field),'binwidth',1);
        
    if 1
        k = 1;
        if i == 1
            filter = edges(2:end) < 350 & edges(2:end) > 150;
            filter2 = edges2(2:end) < 350 & edges2(2:end) > 150;
            d = edges(2:end);
            c = counts(2:end)/max(counts);
            k = 783.5052e-003;
            d2 = edges2(2:end);
            c2 = counts2/max(counts2);
            
            d2 = d2(filter2);
            c2 = k * c2(filter2);
            d = d(filter);
            c = c(filter);
        end
        
        fig = figure(1);
        ax = gca;
        box on
        hold on
        subplot(1,2,1)
        hold on
        grid on
        box on
        p(indx) = plot(edges1(2:end),counts1/max(counts1),'-','linewidth',1.5,'color',colors(i,:));
        xlabel('Distance $[\mu m]$','fontsize',22,'interpreter','latex'); ylabel('Normalized counts','fontsize',22,'interpreter','latex')
        xlim([0, 5e3])
        ylim([0, 1.1])
        title('a) without alignment','fontsize',22,'interpreter','latex');
        set(gca, 'fontsize',18)
        subplot(1,2,2)
        hold on
        plot(edges(2:end),counts/max(counts),'o','linewidth',1.5,'color',colors(i,:));
        plot(edges2(2:end),k * counts2/max(counts2),'-','linewidth',2.5,'color',colors(i,:));
        title('b) with alignment','fontsize',22,'interpreter','latex');ylabel('Normalized counts','fontsize',22,'interpreter','latex')
        xlabel('Distance $[\mu m]$','fontsize',22,'interpreter','latex');% ylabel('Normalized counts','fontsize',22,'interpreter','latex')
        xlim([0, 100])
        ylim([0, 1.1])
        set(gca, 'fontsize',18)
        grid on
        box on
        axes('position',[.71 .65 .15 .2])
        box on
        grid on
        hold on
        plot(d,c,'o','linewidth',1.5,'color',colors(1,:));
        plot(d2,c2,'-','linewidth',2.5,'color',colors(1,:));
        axis tight
        set(gca, 'FontSize', 14)

        set(fig, 'Units','centimeters','PaperUnits','centimeters', 'PaperSize',[36, 12],'PaperPosition',[0, 0, 36, 12],'Position',[0 0 36 12])
        indx = indx + 1;
    end
end
legend(p(:),{'M3','M4','M5','M6'},'fontsize',22,'interpreter','latex')

print(fig,'../../figures/alignment_comparison.pdf', '-dpdf','-r600','-painters')

%% functions
function hitdat = loaddat(fileID, formatSpec, datblocks)
    for i = 1:datblocks
        field = strcat('plane',num2str(i-1));
        data = textscan(fileID,formatSpec,'HeaderLines',1);
        hitdat.(field) = data{1};
    end
end