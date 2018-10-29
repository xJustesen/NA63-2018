clc; clear all; close all;
datpath = '/home/christian/Dropbox/speciale/data/';
energy = linspace(0, 80, 80);

filepath2 = strcat(datpath,'energy_sim_amorphous1_80GeV_1.5mm_monochrome_5.txt');
filepath3 = strcat(datpath,'energy_sim_amorphous1_80GeV_1.5mm_monochrome_10.txt');
filepath4 = strcat(datpath,'energy_sim_amorphous1_80GeV_1.5mm_monochrome_25.txt');
filepath5 = strcat(datpath,'energy_sim_amorphous1_80GeV_1.5mm_monochrome_40.txt');
filepath6 = strcat(datpath,'energy_sim_amorphous1_80GeV_1.5mm_monochrome_25_noSA.txt');
filepath7 = strcat(datpath,'energy_sim_amorphous1_80GeV_1.5mm_monochrome_5_noSA.txt');
filepath8 = strcat(datpath,'energy_sim_amorphous1_80GeV_1.5mm_monochrome_10_noSA.txt');
filepath9 = strcat(datpath,'energy_sim_amorphous1_80GeV_1.5mm_monochrome_40_noSA.txt');

nrg_80GeV_sim_5 = load(filepath2) * 6.2415091E9;
nrg_80GeV_sim_10 = load(filepath3) * 6.2415091E9;
nrg_80GeV_sim_25 = load(filepath4) * 6.2415091E9;
nrg_80GeV_sim_40 = load(filepath5) * 6.2415091E9;
nrg_80GeV_sim_25_noSA = load(filepath6) * 6.2415091E9;
nrg_80GeV_sim_5_noSA = load(filepath7) * 6.2415091E9;
nrg_80GeV_sim_10_noSA = load(filepath8) * 6.2415091E9;
nrg_80GeV_sim_40_noSA = load(filepath9) * 6.2415091E9;

[counts_sim_80GeV_mono_5, ~] = hist(nrg_80GeV_sim_5(nrg_80GeV_sim_5 < 80 & nrg_80GeV_sim_5 > 0), energy);
counts_sim_80GeV_mono_5_norm = energy .* counts_sim_80GeV_mono_5/3E+06;

[counts_sim_80GeV_mono_5_noSA, ~] = hist(nrg_80GeV_sim_5_noSA(nrg_80GeV_sim_5_noSA < 80 & nrg_80GeV_sim_5_noSA > 0), energy);
counts_sim_80GeV_mono_5_norm_noSA = energy .* counts_sim_80GeV_mono_5_noSA/3E+06;

[counts_sim_80GeV_mono_10, ~] = hist(nrg_80GeV_sim_10(nrg_80GeV_sim_10 < 80 & nrg_80GeV_sim_10 > 0), energy);
counts_sim_80GeV_mono_10_norm = energy .* counts_sim_80GeV_mono_10/5E+06;

[counts_sim_80GeV_mono_10_noSA, ~] = hist(nrg_80GeV_sim_10_noSA(nrg_80GeV_sim_10_noSA < 80 & nrg_80GeV_sim_10_noSA > 0), energy);
counts_sim_80GeV_mono_10_noSA_norm = energy .* counts_sim_80GeV_mono_10_noSA/1E+06;

[counts_sim_80GeV_mono_25, ~] = hist(nrg_80GeV_sim_25(nrg_80GeV_sim_25 < 80 & nrg_80GeV_sim_25 > 0), energy);
counts_sim_80GeV_mono_25_norm = energy .* counts_sim_80GeV_mono_25/5E+06;

[counts_sim_80GeV_mono_25_noSA, ~] = hist(nrg_80GeV_sim_25_noSA(nrg_80GeV_sim_25_noSA < 80 & nrg_80GeV_sim_25_noSA > 0), energy);
counts_sim_80GeV_mono_25_noSA_norm = energy .* counts_sim_80GeV_mono_25_noSA/1E+06;

[counts_sim_80GeV_mono_40, ~] = hist(nrg_80GeV_sim_40(nrg_80GeV_sim_40 < 80 & nrg_80GeV_sim_40 > 0), energy);
counts_sim_80GeV_mono_40_norm = energy .* counts_sim_80GeV_mono_40/10E+06;

[counts_sim_80GeV_mono_40_noSA, ~] = hist(nrg_80GeV_sim_40_noSA(nrg_80GeV_sim_40_noSA < 80 & nrg_80GeV_sim_40_noSA > 0), energy);
counts_sim_80GeV_mono_40_noSA_norm = energy .* counts_sim_80GeV_mono_40_noSA/10E+06;


f = figure;
hold on
box on
title('80GeV e- ; 1.5mm C','fontsize',22,'interpreter','latex')

plot(energy, counts_sim_80GeV_mono_5_norm/1.5, '-','linewidth',2.5)
plot(energy, counts_sim_80GeV_mono_10_norm/1.5, '-','linewidth',2.5)
plot(energy, counts_sim_80GeV_mono_25_norm/1.5, '-','linewidth',2.5)
plot(energy, counts_sim_80GeV_mono_40_norm/1.5, '-','linewidth',2.5)

% plot(energy, counts_sim_80GeV_mono_5_norm_noSA * max(counts_sim_80GeV_mono_5_norm/1.5)/max(counts_sim_80GeV_mono_5_norm_noSA),':','linewidth',2,'color',   [      0    0.4470    0.7410])
% plot(energy, counts_sim_80GeV_mono_10_noSA_norm * max(counts_sim_80GeV_mono_10_norm/1.5)/max(counts_sim_80GeV_mono_10_noSA_norm),':','linewidth',2,'color',[ 0.8500    0.3250    0.0980])
% plot(energy, counts_sim_80GeV_mono_25_noSA_norm * max(counts_sim_80GeV_mono_25_norm/1.5)/max(counts_sim_80GeV_mono_25_noSA_norm),':','linewidth',2,'color',[ 0.9290    0.6940    0.1250])
% plot(energy, counts_sim_80GeV_mono_40_noSA_norm * max(counts_sim_80GeV_mono_40_norm/1.5)/max(counts_sim_80GeV_mono_40_noSA_norm),':','linewidth',2,'color',[ 0.4940    0.1840    0.5560])

legend({'5 GeV','10 GeV','25 GeV','40 GeV'},'fontsize',22,'interpreter','latex');
set(gca, 'FontSize', 14)
xlabel('Energy [GeV]','fontsize',22,'interpreter','latex');ylabel('dP/dE [1/mm]','fontsize',22,'interpreter','latex')
box on
grid on

set(f, 'Units','centimeters','PaperUnits','centimeters', 'PaperSize',[18, 12],'PaperPosition',[0, 0, 18, 12],'Position',[0 0 18 12])
% print('../../figures/tracking_res.pdf', '-dpdf','-r600','-painters')
