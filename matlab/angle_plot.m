clear all; clc; close all;
datpath = '/home/christian/Dropbox/speciale/data/';

colors =    [       0    0.4470    0.7410
                0.8500    0.3250    0.0980
                0.9290    0.6940    0.1250
                0.4940    0.1840    0.5560
                0.4660    0.6740    0.1880
                0.3010    0.7450    0.9330
                0.6350    0.0780    0.1840
            ];


ang_in_sim = [];
ang_out_sim = [];
ang_in_dat = [];
ang_out_dat = [];

sim = 1:5;
dat = [86, 87, 88, 89];

for i = sim
    filepath = strcat(datpath,['angles_sim_amorphous',num2str(i)','_80GeV_04012019.txt']);
   ang_mat = load(filepath);
   ang_in_sim = [ang_in_sim; ang_mat(:,1)];
   ang_out_sim = [ang_out_sim; ang_mat(:,2)];
end

for i = dat
   filepath = strcat(datpath,['angles_',num2str(i),'.txt']);
   ang_mat = load(filepath);
   ang_in_dat = [ang_in_dat; ang_mat(:,1)];
   ang_out_dat = [ang_out_dat; ang_mat(:,2)];
end

bins = linspace(-0.01, 0.01, 300);

h1 = histogram(ang_in_sim(ang_in_sim > bins(1) & ang_in_sim < bins(end)), bins);
counts_in_sim = h1.Values;
h2 = histogram(ang_out_sim(ang_out_sim > bins(1) & ang_out_sim < bins(end)), bins);
counts_out_sim = h2.Values;
h3 = histogram(ang_in_dat(ang_in_dat > bins(1) & ang_in_dat < bins(end)), bins);
counts_in_dat = h3.Values;
h4 = histogram(ang_out_dat(ang_out_dat > bins(1) & ang_out_dat < bins(end)), bins);
counts_out_dat = h4.Values;

shift = 0 * 3e-4;

f = figure;
subplot(2,1,1)
hold on;box on
plot(bins(2:end)+shift, counts_in_dat, 'o','color',colors(1,:),'linewidth',1.5)
plot(bins(2:end), max(counts_in_dat)/max(counts_in_sim) * counts_in_sim,'linewidth',2.5,'color',colors(1,:))
set(gca, 'FontSize', 14)
title('a) M3-M4 $\rightarrow$ MM','fontsize',22,'interpreter','latex')
legend({'Data','Simulation'},'fontsize',22,'interpreter','latex');
xlabel('ang (rad)','fontsize',22,'interpreter','latex'); ylabel('counts','fontsize',22,'interpreter','latex');
xlim([bins(1), bins(end)])
subplot(2,1,2)
hold on;box on
plot(bins(2:end)+shift, counts_out_dat, 's','color',colors(2,:),'linewidth',1.5)
plot(bins(2:end), max(counts_out_dat)/max(counts_out_sim) * counts_out_sim,'linewidth',2.5,'color',colors(2,:))
set(gca, 'FontSize', 14)
title('b) MM $\rightarrow$ M5-M6','fontsize',22,'interpreter','latex')
xlabel('ang (rad)','fontsize',22,'interpreter','latex'); ylabel('counts','fontsize',22,'interpreter','latex');
legend({'Data','Simulation'},'fontsize',22,'interpreter','latex');
xlim([bins(1), bins(end)])

% set(f, 'Units','centimeters','PaperUnits','centimeters', 'PaperSize',[18, 12],'PaperPosition',[0, 0, 18, 12],'Position',[0 0 18 12])
% print(f, '../../presentation/figures/angle_distro_40GeV_sim', '-dsvg','-r600','-painters')
