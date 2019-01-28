clc; clear all; close all;
datpath = '/home/christian/Documents/cern2018/simdata/';
ulim = 80;
energy = linspace(0, ulim, 80);
%% load
filepath1 = strcat(datpath,'energy_sim_monochrome_80GeV_1.5mm_5GeV.txt');
filepath2 = strcat(datpath,'energy_sim_monochrome_80GeV_1.5mm_15GeV.txt');
filepath3 = strcat(datpath,'energy_sim_monochrome_80GeV_1.5mm_25GeV.txt');
filepath4 = strcat(datpath,'energy_sim_monochrome_80GeV_1.5mm_35GeV.txt');
filepath5 = strcat(datpath,'energy_sim_monochrome_80GeV_1.5mm_40GeV.txt');
filepath6 = strcat(datpath,'energy_sim_monochrome_80GeV_1.5mm_50GeV.txt');
filepath7 = strcat(datpath,'energy_sim_monochrome_80GeV_1.5mm_60GeV.txt');
filepath8 = strcat(datpath,'energy_sim_monochrome_80GeV_1.5mm_25GeV_inf.txt');
filepath9 = strcat(datpath,'energy_sim_monochrome_80GeV_1.5mm_35GeV_inf.txt');
filepath10 = strcat(datpath,'energy_sim_monochrome_80GeV_1.5mm_40GeV_inf.txt');
filepath11 = strcat(datpath,'energy_sim_monochrome_80GeV_1.5mm_5GeV_inf.txt');
filepath12 = strcat(datpath,'energy_sim_monochrome_80GeV_1.5mm_15GeV_inf.txt');

nrg_80GeV_sim_5 = load(filepath1) * 6.2415091E9;
nrg_80GeV_sim_15 = load(filepath2) * 6.2415091E9;
nrg_80GeV_sim_25 = load(filepath3) * 6.2415091E9;
nrg_80GeV_sim_35 = load(filepath4) * 6.2415091E9;
nrg_80GeV_sim_40 = load(filepath5) * 6.2415091E9;
nrg_80GeV_sim_50 = load(filepath6) * 6.2415091E9;
nrg_80GeV_sim_60 = load(filepath7) * 6.2415091E9;
nrg_80GeV_sim_5_inf = load(filepath11) * 6.2415091E9;
nrg_80GeV_sim_15_inf = load(filepath12) * 6.2415091E9;
nrg_80GeV_sim_25_inf = load(filepath8) * 6.2415091E9;
nrg_80GeV_sim_35_inf = load(filepath9) * 6.2415091E9;
nrg_80GeV_sim_40_inf = load(filepath10) * 6.2415091E9;

%% make hist
[counts_sim_80GeV_mono_5_inf, ~] = hist(nrg_80GeV_sim_5_inf(nrg_80GeV_sim_5_inf < ulim & nrg_80GeV_sim_5_inf > 0), energy);
counts_sim_80GeV_mono_5_norm_inf = energy .* counts_sim_80GeV_mono_5_inf/1E+05;

[counts_sim_80GeV_mono_5, ~] = hist(nrg_80GeV_sim_5(nrg_80GeV_sim_5 < ulim & nrg_80GeV_sim_5 > 0), energy);
counts_sim_80GeV_mono_5_norm = energy .* counts_sim_80GeV_mono_5/1E+05;

[counts_sim_80GeV_mono_15, ~] = hist(nrg_80GeV_sim_15(nrg_80GeV_sim_15 < ulim & nrg_80GeV_sim_15 > 0), energy);
counts_sim_80GeV_mono_15_norm = energy .* counts_sim_80GeV_mono_15/1E+05;

[counts_sim_80GeV_mono_15_inf, ~] = hist(nrg_80GeV_sim_15_inf(nrg_80GeV_sim_15_inf < ulim & nrg_80GeV_sim_15_inf > 0), energy);
counts_sim_80GeV_mono_15_norm_inf = energy .* counts_sim_80GeV_mono_15_inf/1E+05;

[counts_sim_80GeV_mono_25, ~] = hist(nrg_80GeV_sim_25(nrg_80GeV_sim_25 < ulim & nrg_80GeV_sim_25 > 0), energy);
counts_sim_80GeV_mono_25_norm = energy .* counts_sim_80GeV_mono_25/1E+06;

[counts_sim_80GeV_mono_25_inf, ~] = hist(nrg_80GeV_sim_25_inf(nrg_80GeV_sim_25_inf < ulim & nrg_80GeV_sim_25_inf > 0), energy);
counts_sim_80GeV_mono_25_norm_inf = energy .* counts_sim_80GeV_mono_25_inf/1E+06;

[counts_sim_80GeV_mono_35, ~] = hist(nrg_80GeV_sim_35(nrg_80GeV_sim_35 < ulim & nrg_80GeV_sim_35 > 0), energy);
counts_sim_80GeV_mono_35_norm = energy .* counts_sim_80GeV_mono_35/1E+06;

[counts_sim_80GeV_mono_35_inf, ~] = hist(nrg_80GeV_sim_35_inf(nrg_80GeV_sim_35_inf < ulim & nrg_80GeV_sim_35_inf > 0), energy);
counts_sim_80GeV_mono_35_norm_inf = energy .* counts_sim_80GeV_mono_35_inf/1E+06;

[counts_sim_80GeV_mono_40, ~] = hist(nrg_80GeV_sim_40(nrg_80GeV_sim_40 < ulim & nrg_80GeV_sim_40 > 0), energy);
counts_sim_80GeV_mono_40_norm = energy .* counts_sim_80GeV_mono_40/1E+06;

[counts_sim_80GeV_mono_40_inf, ~] = hist(nrg_80GeV_sim_40_inf(nrg_80GeV_sim_40_inf < ulim & nrg_80GeV_sim_40_inf > 0), energy);
counts_sim_80GeV_mono_40_norm_inf = energy .* counts_sim_80GeV_mono_40_inf/1E+06;

[counts_sim_80GeV_mono_50, ~] = hist(nrg_80GeV_sim_50(nrg_80GeV_sim_50 < ulim & nrg_80GeV_sim_50 > 0), energy);
counts_sim_80GeV_mono_50_norm = energy .* counts_sim_80GeV_mono_50/1E+06;

[counts_sim_80GeV_mono_60, ~] = hist(nrg_80GeV_sim_60(nrg_80GeV_sim_60 < ulim & nrg_80GeV_sim_60 > 0), energy);
counts_sim_80GeV_mono_60_norm = energy .* counts_sim_80GeV_mono_60/1E+06;

%% plot
colors =    [       0    0.4470    0.7410
                0.8500    0.3250    0.0980
                0.9290    0.6940    0.1250
                0.4940    0.1840    0.5560
                0.4660    0.6740    0.1880
                0.3010    0.7450    0.9330
                0.6350    0.0780    0.1840
            ];

% f = figure;
% hold on
% box on
% title('80GeV e- ; 1.5mm C | peaks unscaled','fontsize',22,'interpreter','latex')
% 
% plot(energy, counts_sim_80GeV_mono_5_norm, '-','linewidth',1.5)
% plot(energy, counts_sim_80GeV_mono_15_norm, '-','linewidth',1.5)
% plot(energy, counts_sim_80GeV_mono_25_norm, '-','linewidth',1.5)
% plot(energy, counts_sim_80GeV_mono_40_norm, '-','linewidth',1.5)
% plot(energy, counts_sim_80GeV_mono_5_norm_inf, '--','linewidth',1.5,'color',colors(1,:))
% plot(energy, counts_sim_80GeV_mono_15_norm_inf, '--','linewidth',1.5,'color',colors(2,:))
% plot(energy, counts_sim_80GeV_mono_25_norm_inf, '--','linewidth',1.5,'color',colors(3,:))
% plot(energy, counts_sim_80GeV_mono_40_norm_inf, '--','linewidth',1.5,'color',colors(4,:))
% 
% 
% legend({'5 GeV','15 GeV','25 GeV','40 GeV'});
% set(gca, 'FontSize', 14)
% xlabel('Energy [GeV]','fontsize',22,'interpreter','latex');ylabel('dP/dE [1/mm]','fontsize',22,'interpreter','latex')
% box on
% grid on
% 
% set(f, 'Units','centimeters','PaperUnits','centimeters', 'PaperSize',[18, 12],'PaperPosition',[0, 0, 18, 12],'Position',[0 0 18 12])
% print('/home/christian/Dropbox/Cern2018Experiment/figures/monochrome_80GeV_1.5mm.pdf', '-dpdf','-r600','-painters')

f = figure;
hold on
box on
plot(energy, counts_sim_80GeV_mono_5_norm_inf, '-','linewidth',1.5)
plot(energy, counts_sim_80GeV_mono_15_norm_inf, '-.','linewidth',1.5)
plot(energy, counts_sim_80GeV_mono_25_norm_inf, '--','linewidth',1.5)
plot(energy, counts_sim_80GeV_mono_40_norm_inf, ':','linewidth',1.5)
legend({'5 GeV','15 GeV','25 GeV','40 GeV'});
set(gca, 'FontSize', 14)
xlabel('Energy [GeV]','fontsize',22,'interpreter','latex');ylabel('dP/dE [1/mm]','fontsize',22,'interpreter','latex')
box on
grid on

set(f, 'Units','centimeters','PaperUnits','centimeters', 'PaperSize',[18, 12],'PaperPosition',[0, 0, 18, 12],'Position',[0 0 18 12])
print('/home/christian/Dropbox/Cern2018Experiment/figures/monochrome_80GeV_1.5mm.pdf', '-dpdf','-r600','-painters')
