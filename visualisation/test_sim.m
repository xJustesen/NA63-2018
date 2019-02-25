clear all; close all

%% test energy
global datpath;
datpath = "~/Dropbox/Cern2018Experiment/spectre/";

E40 = linspace(0, 40, 40);
fprintf('Loading CERN data\n');
[counts_dat_40GeV_bg_norm,counts_dat_40GeV_bg_norm_err] = spectrum(E40, 'energy_', '.txt', 74, 3029506);
[counts_dat_40GeV_amorph_norm_tot,counts_dat_40GeV_amorph_norm_tot_err] = spectrum(E40, 'energy_', '.txt', [73 75:78], 1290988 + 1361162 + 1447462 + 715126 + 1456319);
[counts_dat_40GeV_aligned_norm_tot,counts_dat_40GeV_aligned_norm_tot_err] = spectrum(E40, 'energy_', '.txt', [71 72 79:81], 142959 + 473324 + 460625 + 1288624 + 1275493);

E80 = linspace(0, 80, 60);
fprintf('Loading SIM data\n');
counts_sim_80GeV_amorph_bg_norm = spectrum(E80, 'energy_sim_amorphous', '_80GeV.txt', [], 25e7, altpath);
counts_sim_80GeV_bg_norm = spectrum(E80, 'energy_sim_background', '_80GeV.txt', [], 25e7, altpath);
counts_sim_40GeV_aligned_norm_woshot = spectrum(E40, 'energy_sim_aligned', '_40GeV_woshot_test.txt', [], 1e6);
counts_sim_40GeV_aligned_norm_worr = spectrum(E40, 'energy_sim_aligned', '_40GeV_worr_test.txt', [], 1e6);
counts_sim_40GeV_aligned_norm_worr2 = spectrum(E40, 'energy_sim_aligned', '_40GeV_worr_test2.txt', [], 1e6);
counts_sim_40GeV_aligned_norm = spectrum(E40, 'energy_sim_aligned', '_40GeV_test.txt', [], 1e6);

dat = counts_dat_40GeV_amorph_norm_tot - counts_dat_40GeV_bg_norm;
err = counts_dat_40GeV_amorph_norm_tot_err+counts_dat_40GeV_bg_norm_err;
sim = counts_sim_40GeV_amorph_bg_norm - counts_sim_40GeV_bg_norm;
eff_40 = callibrate(dat,err,sim)

fprintf('Plotting data\n');
colors = [        0    0.4470    0.7410
             0.8500    0.3250    0.0980
             0.9290    0.6940    0.1250
             0.4940    0.1840    0.5560
             0.4660    0.6740    0.1880
             0.3010    0.7450    0.9330
             0.6350    0.0780    0.1840
         ];


figure
hold on
title('Amorph data + sim')
errorbar(E40, counts_dat_40GeV_amorph_norm_tot,counts_dat_40GeV_amorph_norm_tot_err,'s','MarkerFaceColor','auto')
errorbar(E40, counts_dat_40GeV_bg_norm,counts_dat_40GeV_bg_norm_err,'o','MarkerFaceColor','auto')
errorbar(E40, counts_dat_40GeV_amorph_norm_tot - counts_dat_40GeV_bg_norm,counts_dat_40GeV_amorph_norm_tot_err+counts_dat_40GeV_bg_norm_err,'^','MarkerFaceColor','auto')
plot(E40, counts_sim_40GeV_amorph_bg_norm,'-','linewidth',1.5,'color',colors(1,:))
plot(E40, counts_sim_40GeV_bg_norm,'-','linewidth',1.5,'color',colors(2,:))
plot(E40, (counts_sim_40GeV_amorph_bg_norm - counts_sim_40GeV_bg_norm),'-','linewidth',1.5,'color',colors(3,:))
legend('Amorph data','BG data','Amorph sim','BG sim')

figure
hold on
title('aligned data + sim')
errorbar(E40, counts_dat_40GeV_aligned_norm_tot - counts_dat_40GeV_bg_norm,counts_dat_40GeV_aligned_norm_tot_err + counts_dat_40GeV_bg_norm_err,'o','MarkerFaceColor','auto')
plot(E40, eff_40 * (counts_sim_40GeV_aligned_norm_woshot - counts_sim_40GeV_bg_norm),'-','linewidth',1.5,'color',colors(3,:))
plot(E40, eff_40 * (counts_sim_40GeV_aligned_norm_worr - counts_sim_40GeV_bg_norm),'-','linewidth',1.5,'color',colors(2,:))
plot(E40, eff_40 * (counts_sim_40GeV_aligned_norm - counts_sim_40GeV_bg_norm), '-','linewidth',1.5,'color',colors(1,:))
legend('Aligned data','BKCnoSchott','BKCnoRR','BKC')

figure
hold on
title('2x vs 1x')
plot(E40, eff_40 * (counts_sim_40GeV_aligned_norm_worr),'-','linewidth',1.5)
plot(E40, eff_40 * (counts_sim_40GeV_aligned_norm_worr2 )/4,'-','linewidth',1.5)


%% funktioner
function k = callibrate(dat, err, sim)
    dat = dat(2:end);
    err = err(2:end);
    sim = sim(2:end);
    chi2 = @(eff) sum(((dat - eff .* sim).^2)./(2 * err.^2)); % chi-squared
%     k = @() barebonesSA(0.8, 1e6, 1e4, 10, 1, 2, chi2); % slow, global
    k = fminsearch(chi2, 0.9); % fast, local
end

function [s, err] = spectrum(energy, filePrefix, fileSuffix, fileNum, N)
    global datpath;
    s = [];

    if isempty(fileNum)
        file = strcat(datpath, filePrefix,fileSuffix);
        s = load(file) * 6.2415091E9;
    else
        for i = fileNum
            file = strcat(datpath, filePrefix,num2str(i),fileSuffix);
            nrg = load(file) * 6.2415091E9;
            s = [s; nrg];
        end
    end
    
    s = hist(s(s < max(energy) & s > min(energy)), energy);
    err = energy .* sqrt(s) / N;
    s = energy .* s / N;
end