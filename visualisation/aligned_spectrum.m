clear all; close all
global datpath;
datpath = "~/Dropbox/Cern2018Experiment/spectre/";
% bremstrahlung = 1.00*1e-3*alpha*16/3*r0^2*(1-I_40_full(:,1)./(E)+(I_40_full(:,1)/(E)).^2)*(7*6*log(183*6^(-1/3)))*Na*di_density/di_m;

E = 40000;  
Na = 6.022*10^23;
alpha = 1/137;
r0 = 2.8*10^-15;
di_density = 3.520*10^6;
di_m = 12.0107;

%% LOAD DATA and MAKE HISTOGRAMS

% 40 GeV
% data
E40 = linspace(0, 40, 40);
counts_dat_40GeV_bg_norm = spectrum('energy_', '.txt', 74);
data_bg_40GeV_1mm = hist(counts_dat_40GeV_bg_norm(counts_dat_40GeV_bg_norm < max(E40) & counts_dat_40GeV_bg_norm > min(E40)), E40);
counts_dat_40GeV_amorph_norm_tot = spectrum('energy_', '.txt', [73 75:78]);
data_amorph_40GeV_1mm = hist(counts_dat_40GeV_amorph_norm_tot(counts_dat_40GeV_amorph_norm_tot < max(E40) & counts_dat_40GeV_amorph_norm_tot > min(E40)), E40);
data_amorph_40GeV_1mm_err = E40 .* sqrt(data_amorph_40GeV_1mm/(1290988 + 1361162 + 1447462 + 715126 + 1456319)^2 + data_bg_40GeV_1mm/(3029506)^2);
counts_dat_40GeV_aligned_norm_tot = spectrum('energy_', '.txt', [71 72 79:81]);%, 142959 + 473324 + 460625 + 1288624 + 1275493);
data_align_40GeV_1mm = hist(counts_dat_40GeV_aligned_norm_tot(counts_dat_40GeV_aligned_norm_tot < max(E40) & counts_dat_40GeV_aligned_norm_tot > min(E40)), E40);
data_align_40GeV_1mm_err = E40 .* sqrt(data_align_40GeV_1mm/(142959+473324+460625+1288624+1275493)^2 + data_bg_40GeV_1mm/(3029506)^2);
data_align_40GeV_1mm = E40 .* (data_align_40GeV_1mm / (142959 + 473324 + 460625 + 1288624 + 1275493) - data_bg_40GeV_1mm/3029506);

% sim
counts_sim_40GeV_amorph_bg_norm = spectrum('energy_sim_amorphous', '_40GeV.txt', []);
counts_sim_40GeV_bg_norm = spectrum('energy_sim_background', '_40GeV.txt', []);
counts_sim_40GeV_aligned_norm = spectrum('energy_sim_aligned', '_40GeV_test.txt', []);
sim_bg_40GeV_1mm = hist(counts_sim_40GeV_amorph_bg_norm(counts_sim_40GeV_amorph_bg_norm < max(E40) & counts_sim_40GeV_amorph_bg_norm > min(E40)), E40);
sim_align_40GeV_1mm = hist(counts_sim_40GeV_aligned_norm(counts_sim_40GeV_aligned_norm < max(E40) & counts_sim_40GeV_aligned_norm > min(E40)), E40);
sim_align_40GeV_1mm = E40 .* (sim_align_40GeV_1mm / 5e6 - sim_bg_40GeV_1mm / 1e8);

[eff_40, ~, ~, ~, data_amorph_40GeV_1mm, sim_amorph_40GeV_1mm] = calibrate(E40, counts_sim_40GeV_amorph_bg_norm,counts_sim_40GeV_bg_norm,...
                          counts_dat_40GeV_amorph_norm_tot,counts_dat_40GeV_bg_norm,...
                          1e8, 1290988 + 1361162 + 1447462 + 715126 + 1456319, 3029506);
fprintf('\t 40 GeV data+sim loaded; sorted into histograms\n')

% 20GeV
% data
E20 = linspace(0, 20, 20);
counts_dat_20GeV_aligned_norm_tot = spectrum('energy_','.txt',[103, 104, 112:114]);
counts_dat_20GeV_amorph_norm_tot = spectrum('energy_', '.txt', [109,115]);
counts_dat_20GeV_bg_norm_tot = spectrum('energy_', '.txt', [105:108, 110, 111]);
data_bg_20GeV_1mm = hist(counts_dat_20GeV_bg_norm_tot(counts_dat_20GeV_bg_norm_tot < max(E20) & counts_dat_20GeV_bg_norm_tot > min(E20)), E20);
data_amorph_20GeV_1mm = hist(counts_dat_20GeV_amorph_norm_tot(counts_dat_20GeV_amorph_norm_tot < max(E20) & counts_dat_20GeV_amorph_norm_tot > min(E20)), E20);
data_amorph_20GeV_1mm_err = E20 .* sqrt(data_amorph_20GeV_1mm/(2764454+82089)^2 + data_bg_20GeV_1mm/(324000+590172+624602+734446+1415716+1224254)^2);
data_align_20GeV_1mm = hist(counts_dat_20GeV_aligned_norm_tot(counts_dat_20GeV_aligned_norm_tot < max(E20) & counts_dat_20GeV_aligned_norm_tot > min(E20)), E20);
data_align_20GeV_1mm_err = E20 .* sqrt(data_align_20GeV_1mm/(876217+467348+286232+258194+74259)^2 + data_bg_20GeV_1mm/(324000+590172+624602+734446+1415716+1224254)^2);
data_align_20GeV_1mm = E20 .* data_align_20GeV_1mm/(876217+467348+286232+258194+74259);

% sim
counts_sim_20GeV_amorph_bg_norm = spectrum('energy_sim_amorphous', '_20GeV.txt', []);
counts_sim_20GeV_bg_norm = spectrum('energy_sim_background', '_20GeV.txt', []);
sim_bg_20GeV_1mm = hist(counts_sim_20GeV_bg_norm(counts_sim_20GeV_bg_norm < max(E20) & counts_sim_20GeV_bg_norm > min(E20)), E20);
counts_sim_20GeV_aligned_RR = spectrum('energy_sim_aligned', '_20GeV_test.txt',  []);
sim_align_20GeV_1mm_RR = hist(counts_sim_20GeV_aligned_RR(counts_sim_20GeV_aligned_RR < max(E20) & counts_sim_20GeV_aligned_RR > min(E20)), E20);
sim_align_20GeV_1mm_RR = E20 .* (sim_align_20GeV_1mm_RR / 5e6 - sim_bg_20GeV_1mm / 1e8);

[eff_20, ~, ~, ~, data_amorph_20GeV_1mm, sim_amorph_20GeV_1mm] = calibrate(E20, counts_sim_20GeV_amorph_bg_norm,...
                    counts_sim_20GeV_bg_norm, counts_dat_20GeV_amorph_norm_tot,...
                    counts_dat_20GeV_bg_norm_tot,...
                    1e8, 2764454 + 82089, 324000 + 590172 + 624602 + 734446 + 1415716 + 1224254);


fprintf('\t 20 GeV data+sim loaded; sorted into histogram\n')

% 80 GeV
% data
E80 = linspace(0, 80, 40);
counts_dat_80GeV_bg_norm = spectrum('energy_', '.txt', [85 90 91]);
data_bg_80GeV_1mm = hist(counts_dat_80GeV_bg_norm(counts_dat_80GeV_bg_norm < max(E80) & counts_dat_80GeV_bg_norm > min(E80)), E80);
counts_dat_80GeV_amorph_norm_tot = spectrum('energy_', '.txt', 86:89);
data_amorph_80GeV_1mm = hist(counts_dat_80GeV_amorph_norm_tot(counts_dat_80GeV_amorph_norm_tot < max(E80) & counts_dat_80GeV_amorph_norm_tot > min(E20)), E80);
data_amorph_80GeV_1mm_err = E80 .* sqrt(data_amorph_80GeV_1mm/(1719461+1172521+1210722+538281)^2 + data_bg_80GeV_1mm/(1911215+1892132+1219189)^2);
counts_dat_80GeV_aligned_norm_tot = spectrum('energy_', '.txt', [84 92:95]);
data_align_80GeV_1mm = hist(counts_dat_80GeV_aligned_norm_tot(counts_dat_80GeV_aligned_norm_tot < max(E80) & counts_dat_80GeV_aligned_norm_tot > min(E80)), E80);
data_align_80GeV_1mm_err = E80 .* sqrt(data_align_80GeV_1mm/(1134032+462880+943921+1833415+35424)^2 + data_bg_80GeV_1mm/(1911215+1892132+1219189)^2);
data_align_80GeV_1mm = E80 .* data_align_80GeV_1mm / (1134032 + 462880 + 943921 + 1833415 + 35424);

% sim
counts_sim_80GeV_amorph_bg_norm = spectrum('energy_sim_amorphous', '_80GeV.txt', []);
counts_sim_80GeV_bg_norm = spectrum('energy_sim_background', '_80GeV.txt', []);
counts_sim_80GeV_aligned_noSH = spectrum('energy_sim_aligned', '_80GeV_woshot_test.txt', []);
sim_bg_80GeV_1mm = hist(counts_sim_80GeV_amorph_bg_norm(counts_sim_80GeV_amorph_bg_norm < max(E80) & counts_sim_80GeV_amorph_bg_norm > min(E80)), E80);
sim_align_80GeV_1mm_noSH = hist(counts_sim_80GeV_aligned_noSH(counts_sim_80GeV_aligned_noSH < max(E80) & counts_sim_80GeV_aligned_noSH > min(E80)), E80);
sim_align_80GeV_1mm_noSH = E80 .* (sim_align_80GeV_1mm_noSH / 5e6 - sim_bg_80GeV_1mm / 1e8);
counts_sim_80GeV_aligned_noRR = spectrum('energy_sim_aligned', '_80GeV_worr_test.txt',  []);
sim_align_80GeV_1mm_noRR = hist(counts_sim_80GeV_aligned_noRR(counts_sim_80GeV_aligned_noRR < max(E80) & counts_sim_80GeV_aligned_noRR > min(E80)), E80);
sim_align_80GeV_1mm_noRR = E80 .* (sim_align_80GeV_1mm_noRR / 5e6 - sim_bg_80GeV_1mm / 1e8);
counts_sim_80GeV_aligned_RR = spectrum('energy_sim_aligned', '_80GeV_test.txt',  []);
sim_align_80GeV_1mm_RR = hist(counts_sim_80GeV_aligned_RR(counts_sim_80GeV_aligned_RR < max(E80) & counts_sim_80GeV_aligned_RR > min(E80)), E80);
sim_align_80GeV_1mm_RR = E80 .* (sim_align_80GeV_1mm_RR / 5e6 - sim_bg_80GeV_1mm / 1e8);

[eff_80, ~, ~, ~, data_amorph_80GeV_1mm, sim_amorph_80GeV_1mm] = calibrate(E80, counts_sim_80GeV_amorph_bg_norm,counts_sim_80GeV_bg_norm,...
                          counts_dat_80GeV_amorph_norm_tot,counts_dat_80GeV_bg_norm,...
                          1e8, 1719461 + 1172521 + 1210722 + 538281, 1911215 + 1892132 + 1219189);
fprintf('\t 80 GeV data+sim loaded; sorted into histogram\n')

%% CALCULATE ENHANCEMENT
% 40 GeV
I_40_RR  = load(strcat(datpath,'sum_angles40GeV1mm.txt'));
bremstrahlung_40GeV_1mm = 1.00*1e-3*alpha*16/3*r0^2*(1-I_40_RR(:,1)./(E)+(I_40_RR(:,1)/(E)).^2)*(7*6*log(183*6^(-1/3)))*Na*di_density/di_m;
enhance_dat_40GeV_1mm = (data_align_40GeV_1mm)./(data_amorph_40GeV_1mm);
enhance_err_40GeV_1mm = abs(enhance_dat_40GeV_1mm) .* sqrt ((data_align_40GeV_1mm_err ./ data_align_40GeV_1mm).^2 +  (data_amorph_40GeV_1mm_err ./ data_amorph_40GeV_1mm).^2);

% 20 GeV
E = 20000;
I_20_RR  = load(strcat(datpath,'sum_angles20GeV1mmRR.txt'));
bremstrahlung_20GeV_1mm = 1.00*1e-3*alpha*16/3*r0^2*(1-I_20_RR(:,1)./(E)+(I_20_RR(:,1)/(E)).^2)*(7*6*log(183*6^(-1/3)))*Na*di_density/di_m;

enhance_dat_20GeV_1mm = (data_align_20GeV_1mm)./(data_amorph_20GeV_1mm);
enhance_err_20GeV_1mm = abs(enhance_dat_20GeV_1mm) .* sqrt ((data_align_20GeV_1mm_err ./ data_align_20GeV_1mm).^2 +  (data_amorph_20GeV_1mm_err ./ data_amorph_20GeV_1mm).^2);

% 80 GeV
I_80_RR  = load(strcat(datpath,'sum_angles80GeV1mmRR.txt'));
I_80_noSH  = load(strcat(datpath,'sum_angles80GeV1mmnoSH.txt'));
I_80_noRR  = load(strcat(datpath,'sum_angles80GeV1mmnoRR.txt'));
E = 80000;
bremstrahlung_80GeV_1mm = 1.00*1e-3*alpha*16/3*r0^2*(1-I_80_RR(:,1)./(E)+(I_80_RR(:,1)/(E)).^2)*(7*6*log(183*6^(-1/3)))*Na*di_density/di_m;

enhance_dat_80GeV_1mm = (data_align_80GeV_1mm)./(data_amorph_80GeV_1mm);
enhance_err_80GeV_1mm = abs(enhance_dat_80GeV_1mm) .* sqrt ((data_align_80GeV_1mm_err ./ data_align_80GeV_1mm).^2 +  (data_amorph_80GeV_1mm_err ./ data_amorph_80GeV_1mm).^2);


%% PLOT
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
box on
title('40GeV e- ; 1.0mm amorphous C','fontsize',14,'interpreter','latex')
errorbar(E40, data_amorph_40GeV_1mm,data_amorph_40GeV_1mm_err,'s','MarkerFaceColor','auto')
plot(E40, eff_40*sim_amorph_40GeV_1mm,'-','linewidth',1.5,'color',colors(2,:))
legend({'Data','Sim', 'Sim nofac'},'fontsize',14,'interpreter','latex')
xlabel('Energy [GeV]','fontsize',14,'interpreter','latex');ylabel('dP/dE [1/mm]','fontsize',14,'interpreter','latex');
xticklabels('auto'); yticklabels('auto')
grid on
ylim([0,1e-3])
ax = gca;
ax.YAxis.Exponent = -3;

f = figure;
hold on
box on
grid on
set(gca, 'FontSize', 14)
plot(I_40_RR(:,1), I_40_RR(:,2),'-','linewidth',1.5,'color',colors(2,:))
xlabel('E [GeV]','fontsize',14,'interpreter','latex')
ylabel('$dP/d\hbar\omega$','fontsize',14,'interpreter','latex')
set(f, 'Units','centimeters','PaperUnits','centimeters', 'PaperSize',[18, 12],'PaperPosition',[0, 0, 18, 12],'Position',[0 0 18 12])


figure
hold on
errorbar(E40, enhance_dat_40GeV_1mm,enhance_err_40GeV_1mm,'o','MarkerFaceColor','auto')
plot(I_40_RR(:,1),  I_40_RR(:,2) ./ bremstrahlung_40GeV_1mm,'linewidth',1.5)
plot(E40, sim_align_40GeV_1mm./sim_amorph_40GeV_1mm,'linewidth',1.5)
legend('Data','Theory','Sim')
xticklabels('auto'); yticklabels('auto')
xlim([0, 40])
ylim([0, 100]);
grid on
box on
title('enhancement 40GeV')

figure
hold on
errorbar(E20, enhance_dat_20GeV_1mm,enhance_err_20GeV_1mm,'o','MarkerFaceColor','auto')
plot(I_20_RR(:,1),  I_20_RR(:,2) ./ bremstrahlung_20GeV_1mm,'linewidth',1.5)
plot(E20, sim_align_20GeV_1mm_RR ./ sim_amorph_20GeV_1mm,'linewidth',1.5)
legend('Data','Theory','Sim')
xticklabels('auto'); yticklabels('auto')
xlim([0, 20])
ylim([0, 100]);
grid on
box on
title('enhancement 20GeV')

figure
hold on
errorbar(E80, enhance_dat_80GeV_1mm,enhance_err_80GeV_1mm,'o','MarkerFaceColor','auto')
plot(I_80_RR(:,1), I_80_RR(:,2) ./ bremstrahlung_80GeV_1mm,'linewidth',1.5)
plot(I_80_noRR(:,1), I_80_noRR(:,2) ./ bremstrahlung_80GeV_1mm,'linewidth',1.5)
plot(I_80_noSH(:,1), I_80_noSH(:,2) ./ bremstrahlung_80GeV_1mm,'linewidth',1.5)
plot(E80, sim_align_80GeV_1mm_RR./sim_amorph_80GeV_1mm,'--','linewidth',1.5,'color',colors(2,:))
plot(E80, sim_align_80GeV_1mm_noRR./sim_amorph_80GeV_1mm,'--','linewidth',1.5,'color',colors(3,:))
plot(E80, sim_align_80GeV_1mm_noSH./sim_amorph_80GeV_1mm,'--','linewidth',1.5,'color',colors(4,:))
legend('Data','LL','noRR','noSH')
xticklabels('auto'); yticklabels('auto')
xlim([0, 80])
ylim([0, 100]);
grid on
box on
title('enhancement 80GeV')


%% funktioner
function [minimum,val, k, f, dat, sim, dat_err] = calibrate(E, amorphsim, bgsim, amorphdat, bgdat, Nsim, Ndat1, Ndat2)
    amorph = hist(amorphsim(amorphsim < max(E) & amorphsim > min(E)), E);
    bg = hist(bgsim(bgsim < max(E) & bgsim > min(E)), E);
    sim = E .* (amorph - bg)/Nsim;
    sim_err = E/Nsim .* sqrt(amorph + bg);

    amorph = hist(amorphdat(amorphdat < max(E) & amorphdat > min(E)), E);
    bg = hist(bgdat(bgdat < max(E) & bgdat > min(E)), E);
    dat = E .* (amorph/Ndat1 - bg/Ndat2);
    dat_err = E .* sqrt( amorph/Ndat1 + bg/Ndat2);

    chi2 = @(eff) sum(((dat(sim_err > 0) - eff .* sim(sim_err > 0)).^2)./(sim_err(sim_err > 0).^2 + dat_err(sim_err > 0).^2)); % chi-squared wo error
    
    f = zeros(100,1);
    k = linspace(0,2,100);
    for i = 1:100
        f(i) = chi2(k(i));
    end


    minimum = fminsearch(chi2, 0.5); % fast, local
    val = chi2(minimum);    
end

function s = spectrum(filePrefix, fileSuffix, fileNum)
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
end