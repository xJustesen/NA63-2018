clear all; close all;
global datpath;
datpath = '/home/christian/Documents/cern2018/simdata/';
altpath = '/home/christian/Dropbox/Cern2018Experiment/spectre/';
E20 = linspace(0,20,20);
E40 = linspace(0,40,30);
E80 = linspace(0,80,80);

%% 20GeV
% DATA 1mm
[counts_dat_20GeV_aligned_norm_tot,counts_dat_20GeV_aligned_norm_tot_err] = spectrum(E20, 'energy_','.txt',[103, 104, 112:114],876217+467348+286232+258194+74259);
[counts_dat_20GeV_amorph_norm_tot,counts_dat_20GeV_amorph_norm_tot_err] = spectrum(E20, 'energy_', '.txt', [109,115], 2764454 + 82089);
[counts_dat_20GeV_bg_norm_tot,counts_dat_20GeV_bg_norm_tot_err] = spectrum(E20, 'energy_', '.txt', [105:108, 110, 111], 324000 + 590172 + 624602 + 734446 + 1415716 + 1224254);

% DATA 1.5mm
[counts_dat_20GeV_amorph_1_5mm_norm_tot,counts_dat_20GeV_amorph_1_5mm_norm_tot_err] = spectrum(E20, 'energy_', '.txt', 66:69 ,1197107 + 432597 + 634614 + 386867);
[counts_dat_20GeV_bg_1_5mm_norm_tot,counts_dat_20GeV_bg_1_5mm_norm_tot_err] = spectrum(E20, 'energy_', '.txt', [60,65] ,812417 + 2530699);
[counts_dat_20GeV_aligned_norm_tot_1_5mm,counts_dat_20GeV_aligned_norm_tot_1_5mm_err] = spectrum(E20, 'energy_', '.txt',61:64,174257 + 474970 + 876508 + 792574);

% SIM 1mm
counts_sim_20GeV_amorph_bg_norm = spectrum(E20, 'energy_sim_amorphous', '_20GeV.txt', [], 500e7, altpath);
counts_sim_20GeV_amorph_norm = spectrum(E20, 'energy_sim_amorphous', '_20GeV_no_background.txt', 1:5,5e6);
counts_sim_20GeV_bg_norm = spectrum(E20, 'energy_sim_background', '_20GeV_16012019.txt', 1:5, 5e7);
counts_sim_20GeV_aligned_norm_woshot = spectrum(E20, 'energy_sim_aligned', '_20GeV_woshot.txt', 1:5, 5e7);
counts_sim_20GeV_aligned_norm_worr = spectrum(E20, 'energy_sim_aligned', '_20GeV_worr.txt', 1:5, 5e7);
counts_sim_20GeV_aligned_norm = spectrum(E20, 'energy_sim_aligned', '_20GeV.txt', 1:5, 5e7);

plot(E20,counts_sim_20GeV_amorph_bg_norm)

% SIM 1.5mm
counts_sim_20GeV_amorph_1_5mm_norm = spectrum(E20, 'energy_sim_amorphous', '_20GeV_1.5mm_no_background.txt', 1:5, 5e6);
counts_sim_20GeV_bg_norm_1_5mm = spectrum(E20, 'energy_sim_background', '_20GeV_1.5mm_16012019.txt', 1:5, 5e6);
counts_sim_20GeV_amorph_bg_norm_1_5mm = spectrum(E20, 'energy_sim_amorphous', '_20GeV_1.5mm_16012019.txt', 1:5, 5e6);
counts_sim_20GeV_aligned_1_5mm_norm_woshot = spectrum(E20, 'energy_sim_aligned', '_20GeV_woshot_1.5mm.txt', 1:5, 5e6);
counts_sim_20GeV_aligned_1_5mm_norm_worr = spectrum(E20, 'energy_sim_aligned', '_20GeV_worr_1.5mm.txt', 1:5, 5e6);
counts_sim_20GeV_aligned_1_5mm_norm = spectrum(E20, 'energy_sim_aligned', '_20GeV_1.5mm.txt', 1:5, 5e6);
%% 40GeV
% DATA 1mm
energy = linspace(0, 40, 40);
[counts_dat_40GeV_bg_norm,counts_dat_40GeV_bg_norm_err] = spectrum(E40, 'energy_', '.txt', 74, 3029506);
[counts_dat_40GeV_amorph_norm_tot,counts_dat_40GeV_amorph_norm_tot_err] = spectrum(E40, 'energy_', '.txt', [73 75:78], 1290988 + 1361162 + 1447462 + 715126 + 1456319);
[counts_dat_40GeV_aligned_norm_tot,counts_dat_40GeV_aligned_norm_tot_err] = spectrum(E40, 'energy_', '.txt', [71 72 79:81], 142959 + 473324 + 460625 + 1288624 + 1275493);

% DATA 1.5mm
[counts_dat_40GeV_amorph_norm_tot_1_5mm,counts_dat_40GeV_amorph_norm_tot_1_5mm_err] = spectrum(E40, 'energy_', '.txt', [32:34 39:41 43], 1449257 + 529097 + 724698 + 134167 + 692475 + 1694966 + 496471);
[counts_dat_40GeV_bg_norm_1_5mm,counts_dat_40GeV_bg_norm_1_5mm_err] = spectrum(E40, 'energy_', '.txt', 31, 2771767);
[counts_dat_40GeV_aligned_norm_tot_1_5mm,counts_dat_40GeV_aligned_norm_tot_1_5mm_err] = spectrum(E40, 'energy_', '.txt', [30 35:38], 172307 + 435890 + 538900 + 363630 + 209144);

% SIM 1mm
counts_sim_40GeV_amorph_bg_norm = spectrum(E40, 'energy_sim_amorphous', '_40GeV_16012019.txt', 1:5, 5*4e6);
counts_sim_40GeV_bg_norm = spectrum(E40, 'energy_sim_background', '_40GeV_16012019.txt', 1:5, 5*4e6);
counts_sim_40GeV_amorph_norm = spectrum(E40, 'energy_sim_amorphous', '_40GeV_no_background.txt', 1:5, 5e6);
counts_sim_40GeV_aligned_norm_woshot = spectrum(E40, 'energy_sim_aligned', '_40GeV_woshot.txt', 1:5, 5e7);
counts_sim_40GeV_aligned_norm_worr = spectrum(E40, 'energy_sim_aligned', '_40GeV_worr.txt', 1:5, 5e7);
counts_sim_40GeV_aligned_norm = spectrum(E40, 'energy_sim_aligned', '_40GeV.txt', 1:5, 5e7);

% SIM 1.5mm
counts_sim_40GeV_amorph_1_5mm_norm = spectrum(E40, 'energy_sim_amorphous', '_40GeV_1.5mm_no_background.txt', 1:5, 5e6);
counts_sim_40GeV_bg_norm_1_5mm = spectrum(E40, 'energy_sim_background', '_40GeV_1.5mm_16012019.txt', 1:5, 5e6);
counts_sim_40GeV_amorph_bg_norm_1_5mm = spectrum(E40, 'energy_sim_amorphous', '_40GeV_1.5mm_16012019.txt', 1:5, 5e6);
counts_sim_40GeV_aligned_1_5mm_norm_woshot = spectrum(E40, 'energy_sim_aligned', '_40GeV_woshot_1.5mm.txt', 1:5, 5e6);
counts_sim_40GeV_aligned_1_5mm_norm_worr = spectrum(E40, 'energy_sim_aligned', '_40GeV_worr_1.5mm.txt', 1:5, 5e6);
counts_sim_40GeV_aligned_1_5mm_norm = spectrum(E40, 'energy_sim_aligned', '_40GeV_1.5mm.txt', 1:5, 5e6);
%% 80GeV
% DAT 1mm
[counts_dat_80GeV_bg_norm,counts_dat_80GeV_bg_norm_err] = spectrum(E80, 'energy_', '.txt', [85 90 91], 1911215 + 1892132 + 1219189);
[counts_dat_80GeV_amorph_norm_tot,counts_dat_80GeV_amorph_norm_tot_err] = spectrum(E80, 'energy_', '.txt', 86:89, 1719461 + 1172521 + 1210722 + 538281);
[counts_dat_80GeV_aligned_norm_tot,counts_dat_80GeV_aligned_norm_tot_err] = spectrum(E80, 'energy_', '.txt', [84 92:95], 1134032 + 462880 + 943921 + 1833415 + 35424);

% SIM 1mm
counts_sim_80GeV_amorph_bg_norm = spectrum(E80, 'energy_sim_amorphous', '_80GeV_16012019.txt', 1:10, 15e7);
counts_sim_80GeV_amorph_norm = spectrum(E80, 'energy_sim_amorphous', '_80GeV_no_background_04012019.txt', 1:5, 5e6);
counts_sim_80GeV_bg_norm = spectrum(E80, 'energy_sim_background', '_80GeV_16012019.txt', 1:10, 15e7);
counts_sim_80GeV_aligned_norm_woshot = spectrum(E80, 'energy_sim_aligned', '_80GeV_woshot.txt', 1:5, 5e6);
counts_sim_80GeV_aligned_norm_worr = spectrum(E80, 'energy_sim_aligned', '_80GeV_worr.txt', 1:5, 5e6);
counts_sim_80GeV_aligned_norm = spectrum(E80, 'energy_sim_aligned', '_80GeV.txt', 1:5, 5e6);

% DATA 1.5mm
[counts_dat_80GeV_bg_norm_1_5mm_tot, counts_dat_80GeV_bg_norm_1_5mm_tot_err] = spectrum(E80, 'energy_', '.txt', [48 54 55], 2497897 + 312302 + 216860);
[counts_dat_80GeV_amorph_norm_1_5mm_tot,counts_dat_80GeV_amorph_norm_1_5mm_tot_err] = spectrum(E80, 'energy_', '.txt', [49:53 56 57], 873246 + 434680 + 847524 + 182889 + 18846 + 392613 + 495068);
[counts_dat_80GeV_aligned_norm_1_5mm_tot,counts_dat_80GeV_aligned_norm_1_5mm_tot_err] = spectrum(E80, 'energy_', '.txt', [46, 47], 431667 + 942149);

% SIM 1.5mm
counts_sim_80GeV_amorph_bg_1_5mm_norm = spectrum(energy, 'energy_sim_amorphous', '_80GeV_1.5mm_16012019.txt', 1:5, 5e6);
counts_sim_80GeV_amorph_1_5mm_norm = spectrum(energy, 'energy_sim_amorphous', '_80GeV_1.5mm_no_background.txt', 1:5, 5e6);
counts_sim_80GeV_bg_1_5mm_norm = spectrum(energy, 'energy_sim_background', '_80GeV_1.5mm_16012019.txt', 1:5, 5e6);
counts_sim_80GeV_aligned_1_5mm_norm_woshot = spectrum(energy, 'energy_sim_aligned', '_80GeV_woshot_1.5mm.txt', 1:5, 5e6);
counts_sim_80GeV_aligned_1_5mm_norm_worr = spectrum(energy, 'energy_sim_aligned', '_80GeV_worr_1.5mm.txt', 1:5, 5e6);
counts_sim_80GeV_aligned_1_5mm_norm = spectrum(energy, 'energy_sim_aligned', '_80GeV_1.5mm.txt', 1:5, 5e6);

%% Calibration factors
% 20 GeV
dat = counts_dat_20GeV_amorph_norm_tot - counts_dat_20GeV_bg_norm_tot;
err = counts_dat_20GeV_amorph_norm_tot_err+counts_dat_20GeV_bg_norm_tot_err;
sim = counts_sim_20GeV_amorph_bg_norm - counts_sim_20GeV_bg_norm;
eff_20 = callibrate(dat, err, sim)

% 40 GeV
dat = counts_dat_40GeV_amorph_norm_tot - counts_dat_40GeV_bg_norm;
err = counts_dat_40GeV_amorph_norm_tot_err+counts_dat_40GeV_bg_norm_err;
sim = counts_sim_40GeV_amorph_bg_norm - counts_sim_40GeV_bg_norm;
eff_40 = callibrate(dat,err,sim)

% 80 GeV
dat = counts_dat_80GeV_amorph_norm_tot - counts_dat_80GeV_bg_norm;
sim = counts_sim_80GeV_amorph_bg_norm - counts_sim_80GeV_bg_norm;
err = counts_dat_80GeV_amorph_norm_tot_err+counts_dat_80GeV_bg_norm_err;
eff_80 = callibrate(dat, err, sim)

%% plot
colors = [        0    0.4470    0.7410
             0.8500    0.3250    0.0980
             0.9290    0.6940    0.1250
             0.4940    0.1840    0.5560
             0.4660    0.6740    0.1880
             0.3010    0.7450    0.9330
             0.6350    0.0780    0.1840
         ];

f = figure;
[ha, ~] = tight_subplot(3,1,[.12 .04],[.05 .05],[.07 .07]);

energy = linspace(0, 20, 20);
axes(ha(1));
hold on
box on
title('20GeV e- ; 1.0mm C','fontsize',22,'interpreter','latex')
errorbar(E20, counts_dat_20GeV_amorph_norm_tot,counts_dat_20GeV_amorph_norm_tot_err,'^','MarkerFaceColor',colors(1,:))
errorbar(E20, counts_dat_20GeV_bg_norm_tot,counts_dat_20GeV_bg_norm_tot_err,'p','MarkerFaceColor',colors(2,:))
errorbar(E20, counts_dat_20GeV_amorph_norm_tot - counts_dat_20GeV_bg_norm_tot,counts_dat_20GeV_amorph_norm_tot_err+counts_dat_20GeV_bg_norm_tot_err,'s','MarkerFaceColor',colors(3,:))
plot(E20, eff_20*counts_sim_20GeV_amorph_bg_norm,'--','linewidth',2.5,'color',colors(1,:))
plot(E20, eff_20*counts_sim_20GeV_bg_norm,':','linewidth',2.5,'color',colors(2,:))
plot(E20, eff_20 * (counts_sim_20GeV_amorph_bg_norm - counts_sim_20GeV_bg_norm), '-','linewidth',2.5,'color',colors(3,:))
xlabel('Energy [GeV]','fontsize',22,'interpreter','latex');ylabel('dP/dE [1/mm]','fontsize',22,'interpreter','latex');
legend({'Amorph run','Background run','Amorph - background'},'location','northwest','fontsize',18,'interpreter','latex')
xticklabels('auto'); yticklabels('auto')
grid on
set(gca, 'FontSize', 18)
ax = gca;
ax.YAxis.Exponent = -3;


axes(ha(2));
hold on
box on
title('40GeV e- ; 1.0mm C','fontsize',22,'interpreter','latex')
errorbar(E40, counts_dat_40GeV_amorph_norm_tot, counts_dat_40GeV_amorph_norm_tot_err,'^','MarkerFaceColor',colors(1,:))
errorbar(E40, counts_dat_40GeV_bg_norm, counts_dat_40GeV_bg_norm_err,'p','MarkerFaceColor',colors(2,:))
errorbar(E40, counts_dat_40GeV_amorph_norm_tot - counts_dat_40GeV_bg_norm,counts_dat_40GeV_amorph_norm_tot_err+counts_dat_40GeV_bg_norm_err,'s','MarkerFaceColor',colors(3,:))
plot(E40, eff_40*counts_sim_40GeV_bg_norm,':','linewidth',2.5,'color',colors(2,:))
plot(E40, eff_40*counts_sim_40GeV_amorph_bg_norm,'--','linewidth',2.5,'color',colors(1,:))
plot(E40, eff_40*(counts_sim_40GeV_amorph_bg_norm - counts_sim_40GeV_bg_norm),'-','linewidth',2.5,'color',colors(3,:))
xlabel('Energy [GeV]','fontsize',22,'interpreter','latex');ylabel('dP/dE [1/mm]','fontsize',22,'interpreter','latex');
legend({'Amorph run','Background run','Amorph - background'},'fontsize',18,'interpreter','latex')
xticklabels('auto'); yticklabels('auto')
grid on
set(gca, 'FontSize', 18)
ax = gca;
ax.YAxis.Exponent = -3;


axes(ha(3));
hold on
box on
title('80GeV e- ; 1.0mm C','fontsize',22,'interpreter','latex')
errorbar(E80, counts_dat_80GeV_amorph_norm_tot, counts_dat_80GeV_amorph_norm_tot_err ,'^','MarkerFaceColor',colors(1,:))
errorbar(E80, counts_dat_80GeV_bg_norm, counts_dat_80GeV_bg_norm_err,'p','MarkerFaceColor',colors(2,:))
errorbar(E80, counts_dat_80GeV_amorph_norm_tot - counts_dat_80GeV_bg_norm, counts_dat_80GeV_amorph_norm_tot_err+counts_dat_80GeV_bg_norm_err,'s','MarkerFaceColor',colors(3,:))
plot(E80,eff_80*counts_sim_80GeV_bg_norm ,':','linewidth',2.5,'color',colors(2,:))
plot(E80, eff_80*counts_sim_80GeV_amorph_bg_norm, '--','linewidth',2.5,'color',colors(1,:))
plot(E80, eff_80* (counts_sim_80GeV_amorph_bg_norm - counts_sim_80GeV_bg_norm),'-','linewidth',2.5,'color',colors(3,:))
legend({'Amorph run','Background run','Amorph - background'},'fontsize',18,'interpreter','latex')
xlabel('Energy [GeV]','fontsize',22,'interpreter','latex');ylabel('dP/dE [1/mm]','fontsize',22,'interpreter','latex');
xticklabels('auto'); yticklabels('auto')
grid on
set(gca, 'FontSize', 18)
ax = gca;
ax.YAxis.Exponent = -3;
 

f = figure;
[ha, ~] = tight_subplot(1,3,[.12 .04],[.1 .15],[.1 .07]);

axes(ha(2));
hold on
box on
title('b)','Interpreter','latex')
errorbar(E40, counts_dat_40GeV_aligned_norm_tot - counts_dat_40GeV_bg_norm,counts_dat_40GeV_aligned_norm_tot_err+counts_dat_40GeV_bg_norm_err,'o','MarkerFaceColor',colors(1,:))
plot(E40, -eff_40 * counts_sim_40GeV_bg_norm + eff_40 * counts_sim_40GeV_aligned_norm,'-','linewidth',2.5)
plot(E40, -eff_40 * counts_sim_40GeV_bg_norm + eff_40 * counts_sim_40GeV_aligned_norm_woshot,'--','linewidth',2.5)
plot(E40, -eff_40 * counts_sim_40GeV_bg_norm + eff_40 * counts_sim_40GeV_aligned_norm_worr,':','linewidth',2.5)
xlabel('Energy [GeV]','fontsize',22,'interpreter','latex');
axpos = get(gca,'Position');
legend({'Data','BKC','BKCnoSchott','BKCnoRR'},'Interpreter','latex','location','northoutside','orientation','horizontal')
set(gca, 'Position', axpos)
set(gca, 'FontSize', 18)
ax = gca;
ax.FontSize = 18;
ax.YAxis.Exponent = -3;
xticklabels('auto'); yticklabels('auto')
grid on

axes(ha(3));
hold on
box on
title('c)','Interpreter','latex')
errorbar(E80, counts_dat_80GeV_aligned_norm_tot - counts_dat_80GeV_bg_norm,counts_dat_80GeV_aligned_norm_tot_err+counts_dat_80GeV_bg_norm_err,'o','MarkerFaceColor',colors(1,:))
plot(E80, -eff_80 * counts_sim_80GeV_bg_norm + eff_80 * counts_sim_80GeV_aligned_norm,'-','linewidth',2.5)
plot(E80, -eff_80 * counts_sim_80GeV_bg_norm + eff_80 * counts_sim_80GeV_aligned_norm_woshot,'--','linewidth',2.5)
plot(E80, -eff_80 * counts_sim_80GeV_bg_norm + eff_80 * counts_sim_80GeV_aligned_norm_worr,':','linewidth',2.5)
xlabel('Energy [GeV]','fontsize',22,'interpreter','latex');
xticklabels('auto'); yticklabels('auto')
xlim([0, 80])
set(gca, 'FontSize', 18)
ax = gca;
ax.YAxis.Exponent = -3;
grid on

axes(ha(1));
title('a)','Interpreter','latex')
hold on
box on
errorbar(E20, counts_dat_20GeV_aligned_norm_tot - counts_dat_20GeV_bg_norm_tot,counts_dat_20GeV_aligned_norm_tot_err+counts_dat_20GeV_bg_norm_tot_err,'o','MarkerFaceColor',colors(1,:))
plot(E20, -eff_20 * counts_sim_20GeV_bg_norm + eff_20 * counts_sim_20GeV_aligned_norm,'-','linewidth',2.5)
plot(E20, -eff_20 * counts_sim_20GeV_bg_norm + eff_20 * counts_sim_20GeV_aligned_norm_woshot,'--','linewidth',2.5)
plot(E20, -eff_20 * counts_sim_20GeV_bg_norm + eff_20 * counts_sim_20GeV_aligned_norm_worr,':','linewidth',2.5)
ylabel('dP/dE [1/mm]','fontsize',22,'interpreter','latex');xlabel('Energy [GeV]','fontsize',22,'interpreter','latex');
set(gca, 'FontSize', 18)
xticklabels('auto'); yticklabels('auto')
grid on
ax = gca;
ax.YAxis.Exponent = -3;

set(f, 'Units','centimeters','PaperUnits','centimeters', 'PaperSize',[36, 18],'PaperPosition',[0, 0, 36, 18],'Position',[0 0 36 18])

%% funktioner
function k = callibrate(dat, err, sim)
    dat = dat(2:end);
    err = err(2:end);
    sim = sim(2:end);
    chi2 = @(eff) sum(((dat - eff .* sim).^2)./(2 * err.^2)); % sum of square errors
    k = barebonesSA(0.8, 1e6, 1e4, 10, 1, 2, chi2); % find minimum with sim. annealeing
end

function [s, err] = spectrum(energy, filePrefix, fileSuffix, fileNum, M, altpath)
    global datpath;
    s = [];
    
    if exist('altpath','var')
        file = strcat(altpath, filePrefix,fileSuffix);
        s = load(file) * 6.2415091E9;
    else
        for i = fileNum
            file = strcat(datpath, filePrefix,num2str(i),fileSuffix);
            nrg = load(file) * 6.2415091E9;
            s = [s; nrg];
        end
    end
    
    s = hist(s(s < max(energy) & s > min(energy)), energy);
    err = energy.*sqrt(s)/M;
    s = energy.*s/M;
end