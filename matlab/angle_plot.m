clear all; clc; close all;
datpath = '/home/christian/Dropbox/speciale/data/';

ang_in_sim = [];
ang_out_sim = [];
ang_in_dat = [];
ang_out_dat = [];

sim = 1:5;
dat = [32, 33, 34, 39, 40, 41, 43];

for i = sim
    filepath = strcat(datpath,['angles_sim_amorphous',num2str(i)','_40GeV.txt']);
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

bins = linspace(-0.02, 0.02, 300);

h1 = histogram(ang_in_sim(ang_in_sim > -0.02 & ang_in_sim < 0.02), bins);
counts_in_sim = h1.Values;
h2 = histogram(ang_out_sim(ang_out_sim > -0.02 & ang_out_sim < 0.02), bins);
counts_out_sim = h2.Values;
h3 = histogram(ang_in_dat(ang_in_dat > -0.02 & ang_in_dat < 0.02), bins);
counts_in_dat = h3.Values;
h4 = histogram(ang_out_dat(ang_out_dat > -0.02 & ang_out_dat < 0.02), bins);
counts_out_dat = h4.Values;

f = figure;
hold on
box on
errorbar(bins(2:end), counts_in_dat, sqrt(counts_in_dat))
errorbar(bins(2:end), counts_out_dat, sqrt(counts_in_dat))
plot(bins(2:end), max(counts_in_dat)/max(counts_in_sim) * counts_in_sim,'linewidth',2)
plot(bins(2:end), max(counts_in_dat)/max(counts_in_sim) * counts_out_sim,'linewidth',2)
xlabel('ang (rad)'); ylabel('counts');
legend('Angle in dat','Angle out dat', 'Angle in sim','Angle out sim');
set(f, 'Units','centimeters','PaperUnits','centimeters', 'PaperSize',[18, 12],'PaperPosition',[0, 0, 18, 12],'Position',[0 0 18 12])
% print(f, '../../figures/angle_distro_40GeV_sim.pdf', '-dpdf','-r600','-painters')
