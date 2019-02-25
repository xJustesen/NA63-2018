clear all; close all; clc;
global datpath simpath;
datpath = '/home/christian/Dropbox/Cern2018Experiment/spectre/cuts/';
simpath = '/home/christian/Dropbox/Cern2018Experiment/grendel/spectre/cuts/';

E = 80000;  
Em = 0.511;
gamma = E/Em;
Na = 6.022*10^23;
alpha = 1/137;
r0 = 2.8*10^-15;
di_density = 3.520*10^6;
di_m = 12.0107;
di_Z = 6;
di_d = 3.57; % Å
e2 = 14.4; % e² i eVÅ

crit_angle_20 = sqrt(4 * di_Z * e2 / (20000 * di_d * 1e-6)) % se Allan's NIMB artikel
crit_angle_40 = crit_angle_20 * sqrt(1/2)
crit_angle_80 = crit_angle_20 * sqrt(1/4)

E20 = linspace(0,20,20);
E40 = linspace(0,40,40);
E80 = linspace(0,80,60);

types = {'aligned','amorph','bg'};

fprintf('Loading data\n')

nevents = load('/home/christian/Dropbox/speciale/code/src/events_run_cuts.txt');
dataruns = nevents(:,1);
events = nevents(:,2);
%% 20 GeV data+sim
% DATA 1mm
runs_20GeV_1mm.(types{1}) = [103, 104, 112:114];
runs_20GeV_1mm.(types{2}) = [109,115];
runs_20GeV_1mm.(types{3}) = [105:108, 110, 111];

events_20GeV_1mm = zeros(1,3); % [aligned, amorph, bg]

for i = 1:length(types)
   runs = runs_20GeV_1mm.(types{i});
   for j = 1:length(runs)
       events_20GeV_1mm(i) = events_20GeV_1mm(i) +  events(dataruns == runs(j));
   end
end

counts_dat_20GeV_aligned_norm_tot = spectrum('energy_','.txt',[103, 104, 112:114]);
counts_dat_20GeV_amorph_norm_tot = spectrum('energy_', '.txt', [109,115]);
counts_dat_20GeV_bg_norm_tot = spectrum('energy_', '.txt', [105:108, 110, 111]);

data_bg_20GeV_1mm = hist(counts_dat_20GeV_bg_norm_tot(counts_dat_20GeV_bg_norm_tot < max(E20) & counts_dat_20GeV_bg_norm_tot > min(E20)), E20);

data_align_20GeV_1mm = hist(counts_dat_20GeV_aligned_norm_tot(counts_dat_20GeV_aligned_norm_tot < max(E20) & counts_dat_20GeV_aligned_norm_tot > min(E20)), E20);
data_align_20GeV_1mm_err = E20 .* sqrt(data_align_20GeV_1mm/events_20GeV_1mm(1)^2 + data_bg_20GeV_1mm/events_20GeV_1mm(3)^2);
data_align_20GeV_1mm = E20 .* (data_align_20GeV_1mm/events_20GeV_1mm(1) - data_bg_20GeV_1mm/events_20GeV_1mm(3));

% DATA 1.5mm
runs_20GeV_1_5mm.(types{1}) = [61, 62, 63, 64];
runs_20GeV_1_5mm.(types{2}) = [66, 67, 68, 69];
runs_20GeV_1_5mm.(types{3}) = [60, 65];

events_20GeV_1_5mm = zeros(1,3); % [aligned, amorph, bg]

for i = 1:length(types)
   runs = runs_20GeV_1_5mm.(types{i});
   for j = 1:length(runs)
       events_20GeV_1_5mm(i) = events_20GeV_1_5mm(i) +  events(dataruns == runs(j));
   end
end
counts_dat_20GeV_amorph_1_5mm_norm_tot = spectrum('energy_', '.txt', 66:69 );
counts_dat_20GeV_bg_1_5mm_norm_tot = spectrum('energy_', '.txt', [60,65] );
counts_dat_20GeV_aligned_norm_tot_1_5mm = spectrum('energy_', '.txt',61:64);

data_bg_20GeV_1_5mm = hist(counts_dat_20GeV_bg_1_5mm_norm_tot(counts_dat_20GeV_bg_1_5mm_norm_tot < max(E20) & counts_dat_20GeV_bg_1_5mm_norm_tot > min(E20)), E20);

data_align_20GeV_1_5mm = hist(counts_dat_20GeV_aligned_norm_tot_1_5mm(counts_dat_20GeV_aligned_norm_tot_1_5mm < max(E20) & counts_dat_20GeV_aligned_norm_tot_1_5mm > min(E20)), E20);
data_align_20GeV_1_5mm_err = E20 .* sqrt(data_align_20GeV_1_5mm/(events_20GeV_1_5mm(1))^2 + data_bg_20GeV_1_5mm/(events_20GeV_1_5mm(3))^2);
data_align_20GeV_1_5mm = E20 .* data_align_20GeV_1_5mm/(events_20GeV_1_5mm(1));

fprintf('\t 20 GeV CERN data loaded\n')

% SIM 1mm
nevents_amorph_20GeV_1mm = load(strcat(simpath, 'events_run_cuts_sim_amorphous_20GeV.txt'));
nevents_bg_20GeV_1mm = load(strcat(simpath, 'events_run_cuts_sim_background_20GeV.txt'));

counts_sim_20GeV_amorph_bg_norm = spectrum('energy_sim_amorphous', '_20GeV.txt', []);
counts_sim_20GeV_bg_norm = spectrum('energy_sim_background', '_20GeV.txt', []);
sim_bg_20GeV_1mm = hist(counts_sim_20GeV_bg_norm(counts_sim_20GeV_bg_norm < max(E20) & counts_sim_20GeV_bg_norm > min(E20)), E20);


% SIM 1.5mm
nevents_amorph_20GeV_1_5mm = load(strcat(simpath, 'events_run_cuts_sim_amorphous_20GeV_1.5mm.txt'));
nevents_bg_20GeV_1_5mm = load(strcat(simpath, 'events_run_cuts_sim_background_20GeV_1.5mm.txt'));

counts_sim_20GeV_bg_norm_1_5mm = spectrum('energy_sim_background', '_20GeV_1.5mm.txt', []);
counts_sim_20GeV_amorph_bg_norm_1_5mm = spectrum('energy_sim_amorphous', '_20GeV_1.5mm.txt', []);
sim_bg_20GeV_1_5mm = hist(counts_sim_20GeV_bg_norm_1_5mm(counts_sim_20GeV_bg_norm_1_5mm < max(E20) & counts_sim_20GeV_bg_norm_1_5mm > min(E20)), E20);

fprintf('\t 20 GeV SIM data loaded\n')

%% 40 GeV data+sim
% DATA 1mm
runs_40GeV_1mm.(types{1}) = [71, 72, 79, 80, 81];
runs_40GeV_1mm.(types{2}) = [73, 75, 76, 77, 78];
runs_40GeV_1mm.(types{3}) = 74;

events_40GeV_1mm = zeros(1,3); % [aligned, amorph, bg]

for i = 1:length(types)
   runs = runs_40GeV_1mm.(types{i});
   for j = 1:length(runs)
       events_40GeV_1mm(i) = events_40GeV_1mm(i) +  events(dataruns == runs(j));
   end
end

counts_dat_40GeV_bg_norm = spectrum('energy_', '.txt', 74);
counts_dat_40GeV_amorph_norm_tot = spectrum('energy_', '.txt', [73 75:78]);
counts_dat_40GeV_aligned_norm_tot = spectrum('energy_', '.txt', [71 72 79:81]);%, 142959 + 473324 + 460625 + 1288624 + 1275493);

data_bg_40GeV_1mm = hist(counts_dat_40GeV_bg_norm(counts_dat_40GeV_bg_norm < max(E40) & counts_dat_40GeV_bg_norm > min(E40)), E40);

data_align_40GeV_1mm = hist(counts_dat_40GeV_aligned_norm_tot(counts_dat_40GeV_aligned_norm_tot < max(E40) & counts_dat_40GeV_aligned_norm_tot > min(E40)), E40);
data_align_40GeV_1mm_err = E40 .* sqrt(data_align_40GeV_1mm/events_40GeV_1mm(1)^2 + data_bg_40GeV_1mm/events_40GeV_1mm(3)^2);
data_align_40GeV_1mm = E40 .* (data_align_40GeV_1mm / events_40GeV_1mm(1) - data_bg_40GeV_1mm/events_40GeV_1mm(3));

% DATA 1.5mm
runs_40GeV_1_5mm.(types{1}) = [30, 35, 36, 37, 38];
runs_40GeV_1_5mm.(types{2}) = [32, 33, 34, 39, 40, 41, 43];
runs_40GeV_1_5mm.(types{3}) = 31;

events_40GeV_1_5mm = zeros(1,3); % [aligned, amorph, bg]

for i = 1:length(types)
   runs = runs_40GeV_1_5mm.(types{i});
   for j = 1:length(runs)
       events_40GeV_1_5mm(i) = events_40GeV_1_5mm(i) +  events(dataruns == runs(j));
   end
end
counts_dat_40GeV_amorph_norm_tot_1_5mm = spectrum('energy_', '.txt', [32:34 39:41 43]);
counts_dat_40GeV_bg_norm_1_5mm = spectrum('energy_', '.txt', 31);
counts_dat_40GeV_aligned_norm_tot_1_5mm = spectrum('energy_', '.txt', [30 35:38]);

data_bg_40GeV_1_5mm = hist(counts_dat_40GeV_bg_norm_1_5mm(counts_dat_40GeV_bg_norm_1_5mm < max(E40) & counts_dat_40GeV_bg_norm_1_5mm > min(E40)), E40);

data_align_40GeV_1_5mm = hist(counts_dat_40GeV_aligned_norm_tot_1_5mm(counts_dat_40GeV_aligned_norm_tot_1_5mm < max(E40) & counts_dat_40GeV_aligned_norm_tot_1_5mm > min(E40)), E40);
data_align_40GeV_1_5mm_err = E40 .* sqrt(data_align_40GeV_1_5mm/events_40GeV_1_5mm(1)^2 + data_bg_40GeV_1_5mm/events_40GeV_1_5mm(3)^2);
data_align_40GeV_1_5mm = E40 .* (data_align_40GeV_1_5mm/events_40GeV_1_5mm(1)  - data_bg_40GeV_1_5mm/events_40GeV_1_5mm(3));

fprintf('\t 40 GeV CERN data loaded\n')

% SIM 1mm
nevents_amorph_40GeV_1mm = load(strcat(simpath, 'events_run_cuts_sim_amorphous_40GeV.txt'));
nevents_bg_40GeV_1mm = load(strcat(simpath, 'events_run_cuts_sim_background_40GeV.txt'));

counts_sim_40GeV_amorph_bg_norm = spectrum('energy_sim_amorphous', '_40GeV.txt', []);
counts_sim_40GeV_bg_norm = spectrum('energy_sim_background', '_40GeV.txt', []);

sim_bg_40GeV_1mm = hist(counts_sim_40GeV_bg_norm(counts_sim_40GeV_bg_norm < max(E40) & counts_sim_40GeV_bg_norm > min(E40)), E40);


% SIM 1.5mm
nevents_amorph_40GeV_1_5mm = load(strcat(simpath, 'events_run_cuts_sim_amorphous_40GeV_1.5mm.txt'));
nevents_bg_40GeV_1_5mm = load(strcat(simpath, 'events_run_cuts_sim_background_40GeV_1.5mm.txt'));

counts_sim_40GeV_bg_norm_1_5mm = spectrum('energy_sim_background', '_40GeV_1.5mm.txt', []);
counts_sim_40GeV_amorph_bg_norm_1_5mm = spectrum('energy_sim_amorphous', '_40GeV_1.5mm.txt', []);

sim_bg_40GeV_1_5mm = hist(counts_sim_40GeV_bg_norm_1_5mm(counts_sim_40GeV_bg_norm_1_5mm < max(E40) & counts_sim_40GeV_bg_norm_1_5mm > min(E40)), E40);

fprintf('\t 40 GeV SIM data loaded\n')
 
%% 80 GeV data+sim
% DAT 1mm
runs_80GeV_1mm.(types{1}) = [84, 92, 93, 94, 95];
runs_80GeV_1mm.(types{2}) = [86, 87, 88, 89];
runs_80GeV_1mm.(types{3}) = [85, 90, 91];

events_80GeV_1mm = zeros(1,3); % [aligned, amorph, bg]

for i = 1:length(types)
   runs = runs_80GeV_1mm.(types{i});
   for j = 1:length(runs)
       events_80GeV_1mm(i) = events_80GeV_1mm(i) +  events(dataruns == runs(j));
   end
end

counts_dat_80GeV_bg_norm = spectrum('energy_', '.txt', [85 90 91]);
counts_dat_80GeV_amorph_norm_tot = spectrum('energy_', '.txt', 86:89);
counts_dat_80GeV_aligned_norm_tot = spectrum('energy_', '.txt', [84 92:95]);

data_bg_80GeV_1mm = hist(counts_dat_80GeV_bg_norm(counts_dat_80GeV_bg_norm < max(E80) & counts_dat_80GeV_bg_norm > min(E80)), E80);

data_align_80GeV_1mm = hist(counts_dat_80GeV_aligned_norm_tot(counts_dat_80GeV_aligned_norm_tot < max(E80) & counts_dat_80GeV_aligned_norm_tot > min(E80)), E80);
data_align_80GeV_1mm_err = E80 .* sqrt(data_align_80GeV_1mm/events_80GeV_1mm(1)^2 + data_bg_80GeV_1mm/events_80GeV_1mm(3)^2);
data_align_80GeV_1mm = E80 .* (data_align_80GeV_1mm / events_80GeV_1mm(1) - data_bg_80GeV_1mm/events_80GeV_1mm(3));

% DATA 1.5mm
runs_80GeV_1_5mm.(types{1}) = [46, 47];
runs_80GeV_1_5mm.(types{2}) = [49, 50, 51, 52, 53, 56, 57];
runs_80GeV_1_5mm.(types{3}) = [48, 54, 55];

events_80GeV_1_5mm = zeros(1,3); % [aligned, amorph, bg]

for i = 1:length(types)
   runs = runs_80GeV_1_5mm.(types{i});
   for j = 1:length(runs)
       events_80GeV_1_5mm(i) = events_80GeV_1_5mm(i) +  events(dataruns == runs(j));
   end
end

counts_dat_80GeV_bg_norm_1_5mm_tot = spectrum('energy_', '.txt', [48 54 55]);
counts_dat_80GeV_amorph_norm_1_5mm_tot = spectrum('energy_', '.txt', [49:53 56 57]);
counts_dat_80GeV_aligned_norm_1_5mm_tot = spectrum('energy_', '.txt', [46, 47]);

data_bg_80GeV_1_5mm = hist(counts_dat_80GeV_bg_norm_1_5mm_tot(counts_dat_80GeV_bg_norm_1_5mm_tot < max(E80) & counts_dat_80GeV_bg_norm_1_5mm_tot > min(E80)), E80);

data_align_80GeV_1_5mm = hist(counts_dat_80GeV_aligned_norm_1_5mm_tot(counts_dat_80GeV_aligned_norm_1_5mm_tot < max(E80) & counts_dat_80GeV_aligned_norm_1_5mm_tot > min(E80)), E80);
data_align_80GeV_1_5mm_err = E80 .* sqrt(data_align_80GeV_1_5mm/events_80GeV_1_5mm(1)^2 + data_bg_80GeV_1_5mm/events_80GeV_1_5mm(3)^2);
data_align_80GeV_1_5mm = E80 .* (data_align_80GeV_1_5mm / events_80GeV_1_5mm(1) - data_bg_80GeV_1_5mm/events_80GeV_1_5mm(3));
fprintf('\t 80 GeV CERN data loaded\n')

% SIM 1mm
nevents_amorph_80GeV_1mm = load(strcat(simpath, 'events_run_cuts_sim_amorphous_80GeV.txt'));
nevents_bg_80GeV_1mm = load(strcat(simpath, 'events_run_cuts_sim_background_80GeV.txt'));

counts_sim_80GeV_amorph_bg_norm = spectrum('energy_sim_amorphous', '_80GeV.txt', []);
counts_sim_80GeV_bg_norm = spectrum('energy_sim_background', '_80GeV.txt', []);

sim_bg_80GeV_1mm = hist(counts_sim_80GeV_bg_norm(counts_sim_80GeV_bg_norm < max(E80) & counts_sim_80GeV_bg_norm > min(E80)), E80);


% SIM 1.5mm
nevents_amorph_80GeV_1_5mm = load(strcat(simpath, 'events_run_cuts_sim_amorphous_80GeV_1.5mm.txt'));
nevents_bg_80GeV_1_5mm = load(strcat(simpath, 'events_run_cuts_sim_background_80GeV_1.5mm.txt'));

counts_sim_80GeV_amorph_bg_1_5mm_norm = spectrum('energy_sim_amorphous', '_80GeV_1.5mm.txt', []);
counts_sim_80GeV_bg_1_5mm_norm = spectrum('energy_sim_background', '_80GeV_1.5mm.txt', []);

sim_bg_80GeV_1_5mm = hist(counts_sim_80GeV_bg_1_5mm_norm(counts_sim_80GeV_bg_1_5mm_norm < max(E80) & counts_sim_80GeV_bg_1_5mm_norm > min(E80)), E80);

fprintf('\t 80 GeV SIM data loaded\n')

%% Kallibreringsfaktorer
fprintf('Calculating callibration factors\n')
% 20 GeV
%1mm
[eff_20, ~, ~, ~, data_amorph_20GeV_1mm, sim_amorph_20GeV_1mm, data_amorph_20GeV_1mm_err] = calibrate(E20,...
                    counts_sim_20GeV_amorph_bg_norm,counts_sim_20GeV_bg_norm,...
                    counts_dat_20GeV_amorph_norm_tot,counts_dat_20GeV_bg_norm_tot,...
                    nevents_amorph_20GeV_1mm, nevents_bg_20GeV_1mm,...
                    events_20GeV_1mm(2), events_20GeV_1mm(3));
                
%1.5mm
[eff_20_1_5mm, ~, ~, ~, data_amorph_20GeV_1_5mm, sim_amorph_20GeV_1_5mm, data_amorph_20GeV_1_5mm_err] = calibrate(E20,...
                    counts_sim_20GeV_amorph_bg_norm_1_5mm,counts_sim_20GeV_bg_norm_1_5mm,...
                    counts_dat_20GeV_amorph_1_5mm_norm_tot,counts_dat_20GeV_bg_1_5mm_norm_tot,...
                    nevents_amorph_20GeV_1_5mm, nevents_bg_20GeV_1_5mm,...
                    events_20GeV_1_5mm(2), events_20GeV_1_5mm(3));

% 40 GeV
%1mm
[eff_40, ~, ~, ~, data_amorph_40GeV_1mm, sim_amorph_40GeV_1mm, data_amorph_40GeV_1mm_err] = calibrate(E40, ...
                          counts_sim_40GeV_amorph_bg_norm,counts_sim_40GeV_bg_norm,...
                          counts_dat_40GeV_amorph_norm_tot,counts_dat_40GeV_bg_norm,...
                          nevents_amorph_40GeV_1mm,nevents_bg_40GeV_1mm,...
                          events_40GeV_1mm(2), events_40GeV_1mm(3));

%1.5mm
[eff_40_1_5mm, ~, ~, ~, data_amorph_40GeV_1_5mm, sim_amorph_40GeV_1_5mm, data_amorph_40GeV_1_5mm_err] = calibrate(E40, ...
                          counts_sim_40GeV_amorph_bg_norm_1_5mm,counts_sim_40GeV_bg_norm_1_5mm,...
                          counts_dat_40GeV_amorph_norm_tot_1_5mm,counts_dat_40GeV_bg_norm_1_5mm,...
                          nevents_amorph_40GeV_1_5mm, nevents_bg_40GeV_1_5mm,...
                          events_40GeV_1_5mm(2), events_40GeV_1_5mm(3));

% 80 GeV
%1mm
[eff_80, ~, ~, ~, data_amorph_80GeV_1mm, sim_amorph_80GeV_1mm, data_amorph_80GeV_1mm_err] = calibrate(E80,...
                          counts_sim_80GeV_amorph_bg_norm,counts_sim_80GeV_bg_norm,...
                          counts_dat_80GeV_amorph_norm_tot,counts_dat_80GeV_bg_norm,...
                          nevents_amorph_80GeV_1mm,nevents_bg_80GeV_1mm,...
                          events_80GeV_1mm(2), events_80GeV_1mm(3));
%1.5mm
[eff_80_1_5mm, ~, ~, ~, data_amorph_80GeV_1_5mm, sim_amorph_80GeV_1_5mm, data_amorph_80GeV_1_5mm_err] = calibrate(E80,...
                          counts_sim_80GeV_amorph_bg_1_5mm_norm,counts_sim_80GeV_bg_1_5mm_norm,...
                          counts_dat_80GeV_amorph_norm_1_5mm_tot,counts_dat_80GeV_bg_norm_1_5mm_tot,...
                          nevents_amorph_80GeV_1_5mm,nevents_bg_80GeV_1_5mm, ...
                          events_80GeV_1_5mm(2), events_80GeV_1_5mm(3));
                      
% fprintf('\n\t20 GeV\t40 GeV\t80 GeV\n%1.1fmm \t%1.3f \t%1.3f \t%1.3f\n%1.1fmm \t%1.3f \t%1.3f \t%1.3f\n',1, eff_20,eff_40,eff_80,1.5, eff_20_1_5mm,eff_40_1_5mm,eff_80_1_5mm);

%% CALCULATE ENHANCEMENT
datpath = '/home/christian/Dropbox/Cern2018Experiment/spectre/cuts/';

% 40 GeV
% 1mm
I_40_RR_1mm  = load(strcat(datpath,'sum_angles40GeV1mmRRcut.txt'));
bremstrahlung_40GeV_1mm_RR = 1.00*1e-3*alpha*16/3*r0^2*(1-I_40_RR_1mm(:,1)./(E)+(I_40_RR_1mm(:,1)/(E)).^2)*(7*6*log(183*6^(-1/3)))*Na*di_density/di_m;
enhance_dat_40GeV_1mm = (data_align_40GeV_1mm)./(data_amorph_40GeV_1mm);
enhance_err_40GeV_1mm = abs(enhance_dat_40GeV_1mm) .* sqrt ((data_align_40GeV_1mm_err ./ data_align_40GeV_1mm).^2 +  (data_amorph_40GeV_1mm_err ./ data_amorph_40GeV_1mm).^2);

% 1.5mm
I_40_RR_1_5mm  = load(strcat(datpath,'sum_angles40GeV1.5mmRRcut.txt'));
bremstrahlung_40GeV_1_5mm_RR = 1.50*1e-3*alpha*16/3*r0^2*(1-I_40_RR_1_5mm(:,1)./(E)+(I_40_RR_1_5mm(:,1)/(E)).^2)*(7*6*log(183*6^(-1/3)))*Na*di_density/di_m;
enhance_dat_40GeV_1_5mm = (data_align_40GeV_1_5mm)./(data_amorph_40GeV_1_5mm);
enhance_err_40GeV_1_5mm = abs(enhance_dat_40GeV_1_5mm) .* sqrt ((data_align_40GeV_1_5mm_err ./ data_align_40GeV_1_5mm).^2 +  (data_amorph_40GeV_1_5mm_err ./ data_amorph_40GeV_1_5mm).^2);


% 20 GeV
E = 20000;
% 1mm
I_20_RR_1mm  = load(strcat(datpath,'sum_angles20GeV1mmRRcut.txt'));
bremstrahlung_20GeV_1mm = 1.00*1e-3*alpha*16/3*r0^2*(1-I_20_RR_1mm(:,1)./(E)+(I_20_RR_1mm(:,1)/(E)).^2)*(7*6*log(183*6^(-1/3)))*Na*di_density/di_m;
enhance_dat_20GeV_1mm = (data_align_20GeV_1mm)./(data_amorph_20GeV_1mm);
enhance_err_20GeV_1mm = abs(enhance_dat_20GeV_1mm) .* sqrt ((data_align_20GeV_1mm_err ./ data_align_20GeV_1mm).^2 +  (data_amorph_20GeV_1mm_err ./ data_amorph_20GeV_1mm).^2);

%1.5mm
I_20_RR_1_5mm  = load(strcat(datpath,'sum_angles20GeV1.5mmRRcut.txt'));
bremstrahlung_20GeV_1_5mm_RR = 1.50*1e-3*alpha*16/3*r0^2*(1-I_20_RR_1_5mm(:,1)./(E)+(I_20_RR_1_5mm(:,1)/(E)).^2)*(7*6*log(183*6^(-1/3)))*Na*di_density/di_m;
enhance_dat_20GeV_1_5mm = (data_align_20GeV_1_5mm)./(data_amorph_20GeV_1_5mm);
enhance_err_20GeV_1_5mm = abs(enhance_dat_20GeV_1_5mm) .* sqrt ((data_align_20GeV_1_5mm_err ./ data_align_20GeV_1_5mm).^2 +  (data_amorph_20GeV_1_5mm_err ./ data_amorph_20GeV_1_5mm).^2);

% 80 GeV
E = 80000;

%1mm
I_80_RR_1mm  = load(strcat(datpath,'sum_angles80GeV1mmRRcut.txt'));
bremstrahlung_80GeV_1mm = 1.00*1e-3*alpha*16/3*r0^2*(1-I_80_RR_1mm(:,1)./(E)+(I_80_RR_1mm(:,1)/(E)).^2)*(7*6*log(183*6^(-1/3)))*Na*di_density/di_m;
enhance_dat_80GeV_1mm = (data_align_80GeV_1mm)./(data_amorph_80GeV_1mm);
enhance_err_80GeV_1mm = abs(enhance_dat_80GeV_1mm) .* sqrt ((data_align_80GeV_1mm_err ./ data_align_80GeV_1mm).^2 +  (data_amorph_80GeV_1mm_err ./ data_amorph_80GeV_1mm).^2);

%1.5mm
I_80_RR_1_5mm  = load(strcat(datpath,'sum_angles80GeV1.5mmRRcut.txt'));
bremstrahlung_80GeV_1_5mm_RR = 1.50*1e-3*alpha*16/3*r0^2*(1-I_80_RR_1_5mm(:,1)./(E)+(I_80_RR_1_5mm(:,1)/(E)).^2)*(7*6*log(183*6^(-1/3)))*Na*di_density/di_m;
enhance_dat_80GeV_1_5mm = (data_align_80GeV_1_5mm)./(data_amorph_80GeV_1_5mm);
enhance_err_80GeV_1_5mm = abs(enhance_dat_80GeV_1_5mm) .* sqrt ((data_align_80GeV_1_5mm_err ./ data_align_80GeV_1_5mm).^2 +  (data_amorph_80GeV_1_5mm_err ./ data_amorph_80GeV_1_5mm).^2);

fprintf('Finished calculating enhancement\n')

%% farver til plot
fprintf('Plotting data\n');
colors = [        0    0.4470    0.7410
             0.8500    0.3250    0.0980
             0.9290    0.6940    0.1250
             0.4940    0.1840    0.5560
             0.4660    0.6740    0.1880
             0.3010    0.7450    0.9330
             0.6350    0.0780    0.1840
         ];

%% amorph plot 1mm
f = figure;
[ha, ~] = tight_subplot(3,1,[.12 .04],[.05 .05],[.07 .07]);

axes(ha(1));
hold on
box on
title('20GeV e- ; 1.0mm amorphous C','fontsize',22,'interpreter','latex')
errorbar(E20, data_amorph_20GeV_1mm,data_amorph_20GeV_1mm_err,'s','MarkerFaceColor','auto')
plot(E20, eff_20 * sim_amorph_20GeV_1mm, '-','linewidth',2.5,'color',colors(2,:))
ylabel('dP/dE [1/mm]','fontsize',22,'interpreter','latex')
legend({'Data','sim'},'fontsize',18,'interpreter','latex')
xticklabels('auto'); yticklabels('auto')
grid on
set(gca, 'FontSize', 18)
ax = gca;
ax.YAxis.Exponent = -3;

axes(ha(2));
hold on
box on
title('40GeV e- ; 1.0mm amorphous C','fontsize',22,'interpreter','latex')
errorbar(E40, data_amorph_40GeV_1mm,data_amorph_40GeV_1mm_err,'s','MarkerFaceColor','auto')
plot(E40, eff_40*sim_amorph_40GeV_1mm,'-','linewidth',2.5,'color',colors(2,:))
xlabel('Energy [GeV]','fontsize',22,'interpreter','latex')
legend({'Data','Sim'},'fontsize',18,'interpreter','latex')
xticklabels('auto'); yticklabels('auto')
grid on
set(gca, 'FontSize', 18)
ax = gca;
ax.YAxis.Exponent = -3;

axes(ha(3));
hold on
box on
title('80GeV e- ; 1.0mm amorphous C','fontsize',22,'interpreter','latex')
errorbar(E80, data_amorph_80GeV_1mm,data_amorph_80GeV_1mm_err,'s','MarkerFaceColor','auto')
plot(E80, eff_80* sim_amorph_80GeV_1mm,'-','linewidth',2.5,'color',colors(2,:))
legend({'Data','Sim'},'fontsize',18,'interpreter','latex')
xlabel('Energy [GeV]','fontsize',22,'interpreter','latex')
xticklabels('auto'); yticklabels('auto')
grid on
set(gca, 'FontSize', 18)
ax = gca;
ax.YAxis.Exponent = -3;

set(f, 'Units','centimeters','PaperUnits','centimeters', 'PaperSize',[18, 36],'PaperPosition',[0, 0, 18, 36],'Position',[0 0 18, 36])

%% amorph plot 1.5mm
f = figure;
[ha, ~] = tight_subplot(3,1,[.12 .04],[.05 .05],[.07 .07]);

axes(ha(1));
hold on
box on
title('20GeV e- ; 1.5mm amorphous C','fontsize',22,'interpreter','latex')
errorbar(E20, data_amorph_20GeV_1_5mm,data_amorph_20GeV_1_5mm_err,'s','MarkerFaceColor','auto')
plot(E20, eff_20_1_5mm * sim_amorph_20GeV_1_5mm, '-','linewidth',2.5,'color',colors(2,:))
xlabel('Energy [GeV]','fontsize',22,'interpreter','latex');ylabel('dP/dE [1/mm]','fontsize',22,'interpreter','latex');
legend({'Data','Sim'},'fontsize',18,'interpreter','latex')
xticklabels('auto'); yticklabels('auto')
grid on
set(gca, 'FontSize', 18)
ax = gca;
ax.YAxis.Exponent = -3;

axes(ha(2));
hold on
box on
title('40GeV e- ; 1.5mm amorphous C','fontsize',22,'interpreter','latex')
errorbar(E40, data_amorph_40GeV_1_5mm,data_amorph_40GeV_1_5mm_err,'s','MarkerFaceColor','auto')
plot(E40, eff_40_1_5mm * sim_amorph_40GeV_1_5mm, '-','linewidth',2.5,'color',colors(2,:))
xlabel('Energy [GeV]','fontsize',22,'interpreter','latex');ylabel('dP/dE [1/mm]','fontsize',22,'interpreter','latex');
legend({'Data','Sim'},'fontsize',18,'interpreter','latex')
xticklabels('auto'); yticklabels('auto')
grid on
set(gca, 'FontSize', 18)
ax = gca;
ax.YAxis.Exponent = -3;

axes(ha(3));
hold on
box on
title('80GeV e- ; 1.5mm amorphous C','fontsize',22,'interpreter','latex')
errorbar(E80, data_amorph_80GeV_1_5mm,data_amorph_80GeV_1_5mm_err,'s','MarkerFaceColor','auto')
plot(E80, eff_80_1_5mm * sim_amorph_80GeV_1_5mm, '-','linewidth',2.5,'color',colors(2,:))
legend({'Data','Sim'},'fontsize',18,'interpreter','latex')
xlabel('Energy [GeV]','fontsize',22,'interpreter','latex');ylabel('dP/dE [1/mm]','fontsize',22,'interpreter','latex');
xticklabels('auto'); yticklabels('auto')
grid on
set(gca, 'FontSize', 18)
ax = gca;
ax.YAxis.Exponent = -3;

set(f, 'Units','centimeters','PaperUnits','centimeters', 'PaperSize',[18, 36],'PaperPosition',[0, 0, 18, 36],'Position',[0 0 18, 36])
 
%% aligned plot 1mm
% f = figure;
% [ha, ~] = tight_subplot(3,1,[.12 .04],[.05 .05],[.07 .07]);
% axes(ha(1));
% title('20GeV e- ; 1.0mm aligned C','fontsize',22,'interpreter','latex')
% hold on
% box on
% errorbar(E20, data_align_20GeV_1mm,data_align_20GeV_1mm_err,'o','MarkerFaceColor','auto')
% plot(E20, sim_align_20GeV_1mm_RR,'-','linewidth',1.5)
% ylabel('dP/dE [1/mm]','fontsize',22,'interpreter','latex');xlabel('Energy [GeV]','fontsize',22,'interpreter','latex');
% set(gca, 'FontSize', 18)
% xticklabels('auto'); yticklabels('auto')
% grid on
% ax = gca;
% ax.YAxis.Exponent = -3;
% legend('Data','Sim')

% axes(ha(2));
% hold on
% box on
% title('b)','Interpreter','latex')
% errorbar(E40, data_align_40GeV_1mm,data_align_40GeV_1mm_err,'o','MarkerFaceColor','auto')
% plot(E40, eff_40 * sim_align_40GeV_1mm_RR,'-','linewidth',2.5)
% plot(E40, eff_40 * sim_align_40GeV_1mm_noRR,'-','linewidth',2.5)
% plot(E40, eff_40 * sim_align_40GeV_1mm_noSH,'-','linewidth',2.5)
% xlabel('Energy [GeV]','fontsize',22,'interpreter','latex');
% axpos = get(gca,'Position');
% % legend({'Data','BKC','BKCnoSchott','BKCnoRR'},'Interpreter','latex','location','northoutside','orientation','horizontal')
% set(gca, 'Position', axpos)
% set(gca, 'FontSize', 18)
% ax = gca;
% ax.FontSize = 18;
% ax.YAxis.Exponent = -3;
% xticklabels('auto'); yticklabels('auto')
% grid on
% legend('Data','LL','noRR','noSH')

% axes(ha(3));
% hold on
% box on
% title('c)','Interpreter','latex')
% errorbar(E80, data_align_80GeV_1mm,data_align_80GeV_1mm_err,'o','MarkerFaceColor','auto')
% plot(E80, eff_80 * sim_align_80GeV_1mm_RR,'-','linewidth',2.5)
% plot(E80, eff_80 * sim_align_80GeV_1mm_noSH,'-','linewidth',2.5)
% plot(E80, eff_80 * sim_align_80GeV_1mm_noRR,'-','linewidth',2.5)
% xlabel('Energy [GeV]','fontsize',22,'interpreter','latex');
% xticklabels('auto'); yticklabels('auto')
% xlim([0, 80])
% set(gca, 'FontSize', 18)
% ax = gca;
% ax.YAxis.Exponent = -3;
% grid on
% legend('Data','LL','noRR','noSH')

% set(f, 'Units','centimeters','PaperUnits','centimeters', 'PaperSize',[18, 36],'PaperPosition',[0, 0, 18, 36],'Position',[0 0 18, 36])

%% aligned plot 1.5mm
% f = figure;
% [ha, ~] = tight_subplot(3,1,[.12 .04],[.05 .05],[.07 .07]);
% 
% axes(ha(1));
% title('a)','Interpreter','latex')
% hold on
% box on
% errorbar(E20, data_align_20GeV_1_5mm,data_align_20GeV_1_5mm_err,'^','MarkerFaceColor','auto')
% plot(E20, eff_20_1_5mm * sim_align_20GeV_1_5mm_RR,'-','linewidth',2.5)
% ylabel('dP/dE [1/mm]','fontsize',22,'interpreter','latex');xlabel('Energy [GeV]','fontsize',22,'interpreter','latex');
% set(gca, 'FontSize', 18)
% xticklabels('auto'); yticklabels('auto')
% grid on
% ax = gca;
% ax.YAxis.Exponent = -3;
% 
% axes(ha(2));
% hold on
% box on
% title('b)','Interpreter','latex')
% errorbar(E40, data_align_40GeV_1_5mm,data_align_40GeV_1_5mm_err,'^','MarkerFaceColor','auto')
% plot(E40, eff_40_1_5mm * sim_align_40GeV_1_5mm_RR,'-','linewidth',2.5)
% xlabel('Energy [GeV]','fontsize',22,'interpreter','latex');
% ax = gca;
% ax.FontSize = 18;
% ax.YAxis.Exponent = -3;
% xticklabels('auto'); yticklabels('auto')
% grid on
% 
% axes(ha(3));
% hold on
% box on
% title('c)','Interpreter','latex')
% errorbar(E80, data_align_80GeV_1_5mm,data_align_80GeV_1_5mm_err,'^','MarkerFaceColor','auto')
% plot(E80, eff_80_1_5mm * sim_align_80GeV_1_5mm_RR,'-','linewidth',2.5)
% xlabel('Energy [GeV]','fontsize',22,'interpreter','latex');
% xticklabels('auto'); yticklabels('auto')
% xlim([0, 80])
% set(gca, 'FontSize', 18)
% ax = gca;
% ax.YAxis.Exponent = -3;
% grid on
% 
% set(f, 'Units','centimeters','PaperUnits','centimeters', 'PaperSize',[18 36],'PaperPosition',[0, 0, 18 36],'Position',[0 0 18 36])

%% enhancement
f = figure;
[ha, ~] = tight_subplot(3,1,[.12 .04],[.05 .05],[.07 .07]);

axes(ha(1))
hold on
errorbar(E20, enhance_dat_20GeV_1mm,enhance_err_20GeV_1mm,'o','MarkerFaceColor','auto')
plot(I_20_RR_1mm(:,1),  I_20_RR_1mm(:,2) ./ bremstrahlung_20GeV_1mm,':','linewidth',1.5)
legend('Data','Theory', 'Sim')
xticklabels('auto'); yticklabels('auto')
ylim([-20, 150]);
grid on
box on
title('enhancement 20GeV 1mm')

axes(ha(2))
hold on
errorbar(E40, enhance_dat_40GeV_1mm,enhance_err_40GeV_1mm,'o','MarkerFaceColor','auto')
plot(I_40_RR_1mm(:,1),  I_40_RR_1mm(:,2) ./ bremstrahlung_40GeV_1mm_RR,':','linewidth',1.5)
legend('Data','LL','noRR','noSH')
xticklabels('auto'); yticklabels('auto')
xlim([0, 40])
ylim([0, 100]);
grid on
box on
title('enhancement 40GeV 1mm')

axes(ha(3))
hold on
errorbar(E80, enhance_dat_80GeV_1mm,enhance_err_80GeV_1mm,'o','MarkerFaceColor','auto')
plot(I_80_RR_1mm(:,1), I_80_RR_1mm(:,2) ./ bremstrahlung_80GeV_1mm,':','linewidth',1.5)
legend('Data','LL','noRR','noSH')
xticklabels('auto'); yticklabels('auto')
xlim([0, 80])
ylim([0, 100]);
grid on
box on
title('enhancement 80GeV 1mm')

set(f, 'Units','centimeters','PaperUnits','centimeters', 'PaperSize',[18, 36],'PaperPosition',[0, 0, 18, 36],'Position',[0 0 18, 36])

f = figure;
[ha, ~] = tight_subplot(3,1,[.12 .04],[.05 .05],[.07 .07]);
axes(ha(1))
hold on
errorbar(E20, enhance_dat_20GeV_1_5mm,enhance_err_20GeV_1_5mm,'o','MarkerFaceColor','auto')
plot(I_20_RR_1_5mm(:,1),  I_20_RR_1_5mm(:,2) ./ bremstrahlung_20GeV_1_5mm_RR,':','linewidth',1.5)
legend('Data','Theory')
xticklabels('auto'); yticklabels('auto')
xlim([0, 20])
ylim([0, 100]);
grid on
box on
title('enhancement 20GeV 1.5mm')
% 
axes(ha(2))
hold on
errorbar(E40, enhance_dat_40GeV_1_5mm,enhance_err_40GeV_1_5mm,'o','MarkerFaceColor','auto')
plot(I_40_RR_1_5mm(:,1),  I_40_RR_1_5mm(:,2) ./ bremstrahlung_40GeV_1_5mm_RR,':','linewidth',1.5)
legend('Data','LL')
xticklabels('auto'); yticklabels('auto')
xlim([0, 40])
ylim([0, 100]);
grid on
box on
title('enhancement 40GeV 1.5mm')

axes(ha(3))
hold on
errorbar(E80, enhance_dat_80GeV_1_5mm,enhance_err_80GeV_1_5mm,'o','MarkerFaceColor','auto')
plot(I_80_RR_1_5mm(:,1), I_80_RR_1_5mm(:,2) ./ bremstrahlung_80GeV_1_5mm_RR,':','linewidth',1.5)
legend('Data','LL')
xticklabels('auto'); yticklabels('auto')
xlim([0, 80])
ylim([0, 100]);
grid on
box on
title('enhancement 80GeV 1.5mm')
set(f, 'Units','centimeters','PaperUnits','centimeters', 'PaperSize',[18, 36],'PaperPosition',[0, 0, 18, 36],'Position',[0 0 18, 36])

%% FUNKTIONER
function [minimum,val, k, f, dat, sim, dat_err] = calibrate(E, amorphsim, bgsim, amorphdat, bgdat, NsimAmorph, NsimBg, Ndat1, Ndat2)
    amorph = hist(amorphsim(amorphsim < max(E) & amorphsim > min(E)), E);
    bg = hist(bgsim(bgsim < max(E) & bgsim > min(E)), E);
    sim = E .* (amorph/NsimAmorph - bg/NsimBg);

    amorph = hist(amorphdat(amorphdat < max(E) & amorphdat > min(E)), E);
    bg = hist(bgdat(bgdat < max(E) & bgdat > min(E)), E);
    dat = E .* (amorph/Ndat1 - bg/Ndat2);
    dat_err = E .* sqrt(amorph/Ndat1^2 + bg/Ndat2^2);

    chi2 = @(eff) norm(dat(dat > 0) - eff .* sim(dat > 0));

    f = zeros(100,1);
    k = linspace(0,2,100);
    for i = 1:100
        f(i) = chi2(k(i));
    end

    minimum = fminsearch(chi2, 0.5); % fast, local
    val = chi2(minimum);    
end

function s = spectrum(filePrefix, fileSuffix, fileNum)
    global datpath simpath;
    s = [];

    if isempty(fileNum)
        file = strcat(simpath, filePrefix,fileSuffix);
        s = load(file) * 6.2415091E9;
    else
        for i = fileNum
            file = strcat(datpath, filePrefix,num2str(i),fileSuffix);
            nrg = load(file) * 6.2415091E9;
            s = [s; nrg];
        end
    end
end
