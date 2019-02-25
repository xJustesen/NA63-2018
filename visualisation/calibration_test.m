%% DATA
%% preprocess
clc;
global datpath;
datpath = "~/Dropbox/Cern2018Experiment/spectre/";

counts_sim_80GeV_amorph_bg_norm = spectrum('energy_sim_amorphous', '_80GeV.txt', []);
counts_sim_80GeV_bg_norm = spectrum('energy_sim_background', '_80GeV.txt', []);
counts_dat_80GeV_bg_norm = spectrum('energy_', '.txt', [85 90 91]);
counts_dat_80GeV_amorph_norm_tot = spectrum('energy_', '.txt', 86:89);
E1 = linspace(0,80,40);
E2 = linspace(0,80,80);
E3 = linspace(0,80,160);
Nsim = 1e8;
Namorph = 1719461 + 1172521 + 1210722 + 538281;
Nbg = 1911215 + 1892132 + 1219189;
colors = [        0    0.4470    0.7410
             0.8500    0.3250    0.0980
             0.9290    0.6940    0.1250
             0.4940    0.1840    0.5560
             0.4660    0.6740    0.1880
             0.3010    0.7450    0.9330
             0.6350    0.0780    0.1840
         ];
bins = 2:160;
%% calibrate
min = zeros(length(bins),1);
for i = 1:length(bins)
    E = linspace(0,80,bins(i));
    min(i) = calibrate(E,counts_sim_80GeV_amorph_bg_norm,counts_sim_80GeV_bg_norm,counts_dat_80GeV_amorph_norm_tot,counts_dat_80GeV_bg_norm,Nsim,Namorph,Nbg);
    fprintf('iter %d of %d\n',i,length(bins))
end
[min1,val1,k1,f1,dat1,sim1] = calibrate(E1,counts_sim_80GeV_amorph_bg_norm,counts_sim_80GeV_bg_norm,counts_dat_80GeV_amorph_norm_tot,counts_dat_80GeV_bg_norm,Nsim,Namorph,Nbg);
[min2,val2,k2,f2,dat2,sim2] = calibrate(E2,counts_sim_80GeV_amorph_bg_norm,counts_sim_80GeV_bg_norm,counts_dat_80GeV_amorph_norm_tot,counts_dat_80GeV_bg_norm,Nsim,Namorph,Nbg);
[min3,val3,k3,f3,dat3,sim3] = calibrate(E3,counts_sim_80GeV_amorph_bg_norm,counts_sim_80GeV_bg_norm,counts_dat_80GeV_amorph_norm_tot,counts_dat_80GeV_bg_norm,Nsim,Namorph,Nbg);
min1
min2
min3


%% plot
figure
hold on
plot(bins, min,'linewidth',1.5)
xlabel('nbins');ylabel('k');
box on; grid on;

figure
hold on
plot(E1, dat1);
plot(E1, sim1);

figure
hold on
plot(E2, dat2);
plot(E2, sim2);

figure
hold on
plot(E3, dat3);
plot(E3, sim3);

figure
hold on
p(1) = plot(k1, f1,'linewidth',1.5);
p(2) = plot(k2, f2,'linewidth',1.5);
p(3) = plot(k3, f3,'linewidth',1.5);
plot(min1,val1,'o','color',colors(1,:),'markersize',5,'markerfacecolor',colors(1,:))
plot(min2,val2,'o','color',colors(2,:),'markersize',5,'markerfacecolor',colors(2,:))
plot(min3,val3,'o','color',colors(3,:),'markersize',5,'markerfacecolor',colors(3,:))
ax = gca;xlim = ax.XLim;ylim = ax.YLim;
annotation('textarrow',[0.3 min1/xlim(2) + 0.02],[0.2 val1/ylim(2) + 0.095],'String',num2str(min1))
annotation('textarrow',[0.6 min2/xlim(2) + 0.04],[0.5 val2/ylim(2) + 0.1],'String',num2str(min2))
annotation('textarrow',[0.3 min3/xlim(2) + 0.025],[0.5 val3/ylim(2) + 0.1],'String',num2str(min3))
legend(p,{'nbins = 40','nbins = 100','nbins = 160'})
xlabel('k');ylabel('chi2')
box on; grid on;

%% function defs
function [minimum,val, k, f, dat, sim] = calibrate(E, amorphsim, bgsim, amorphdat, bgdat, Nsim, Ndat1, Ndat2)
    amorph = hist(amorphsim(amorphsim < max(E) & amorphsim > min(E)), E);
    bg = hist(bgsim(bgsim < max(E) & bgsim > min(E)), E);
    sim = E .* (amorph - bg)/Nsim;
    sim_err = E/Nsim .* sqrt(amorph + bg);

    amorph = hist(amorphdat(amorphdat < max(E) & amorphdat > min(E)), E);
    bg = hist(bgdat(bgdat < max(E) & bgdat > min(E)), E);
    dat = E .* (amorph/Ndat1 - bg/Ndat2);
    dat_err = E .* sqrt( amorph/Ndat1 + bg/Ndat2);

    chi2 = @(eff) sum(((dat(sim_err > 0) - eff .* sim(sim_err > 0)).^2)./(sim_err(sim_err > 0).^2 + dat_err(sim_err > 0).^2)); % chi-squared wo error
%     chi2 = @(eff) norm(dat - eff .* sim);
%     chi2 = @(eff) trapz(E, dat - eff .* sim);

    
    
    f = zeros(100,1);
    k = linspace(0,2,100);
    for i = 1:100
        f(i) = chi2(k(i));
    end


%     minimum = fzero(chi2, 0.5); % fast, local
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