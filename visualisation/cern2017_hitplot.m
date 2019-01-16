clear all; clc; close all
datpath = '/home/christian/Documents/cern2018/simdata/';

%% load data
input_pos   = load('/home/christian/Dropbox/Cern2018Experiment/kode/beam_position_4mm.txt');
input_ang   = load('/home/christian/Dropbox/Cern2018Experiment/kode/beam_direction_4mm.txt'); 
sampled     = load('/home/christian/Documents/cern2018/simdata/2017_hit_data.txt');

xpos = input_pos(:,3);
ypos = input_pos(:,4);
xang = input_ang(:,3);
yang = input_ang(:,4);

[xpc, xp] = histcounts(sampled(:,1), 2000);
[ypc, yp] = histcounts(sampled(:,2), 2000);
[xac, xa] = histcounts(sampled(:,3), 2000);
[yac, ya] = histcounts(sampled(:,4), 2000);

figure
hold on
plot(xp(1:end-1),xpc/max(xpc),'-')
plot(xpos, input_pos(:,1)/max(input_pos(:,1)))

figure
hold on
plot(yp(1:end-1),ypc/max(ypc),'-')
plot(ypos, input_pos(:,2)/max(input_pos(:,2)))

figure
hold on
plot(xa(1:end-1),xac/max(xac),'-')
plot(xang * 1e-6, input_ang(:,1)/max(input_ang(:,1)))

figure
hold on
plot(ya(1:end-1),yac/max(yac),'-')
plot(yang * 1e-6, input_ang(:,2)/max(input_ang(:,2)))

%%
energy = linspace(0,50,50);

filepath = strcat(datpath,'energy_sim_amorphous_50GeV_622mm.txt');
nrg_50GeV_sim_amorph = load(filepath) * 6.2415091E9;
[counts_sim_50GeV_amorph, ~] = hist(nrg_50GeV_sim_amorph(nrg_50GeV_sim_amorph < 50 & nrg_50GeV_sim_amorph > 0), energy);
counts_sim_50GeV_amorph_norm = energy.*counts_sim_50GeV_amorph/1e6;

filepath1 = strcat(datpath, 'energy_44.txt');
nrg_50GeV_dat1 = load(filepath1) * 6.2415091E9;
[counts_dat_50GeV_amorph_1, ~] = hist(nrg_50GeV_dat1(nrg_50GeV_dat1 < 50 & nrg_50GeV_dat1 > 0), energy);
counts_dat_50GeV_amorph_norm_1 = energy.*counts_dat_50GeV_amorph_1/201678;

filepath2 = strcat(datpath, 'energy_82_2017.txt');
nrg_50GeV_dat2 = load(filepath2) * 6.2415091E9;
[counts_dat_50GeV_bg_1, ~] = hist(nrg_50GeV_dat2(nrg_50GeV_dat2 < 50 & nrg_50GeV_dat2 > 0), energy);
counts_dat_50GeV_bg_norm_1 = energy.*counts_dat_50GeV_bg_1/274742;


delta_50 = @(k) norm(((counts_dat_50GeV_amorph_norm_1 - counts_dat_50GeV_bg_norm_1) - k .* counts_sim_50GeV_amorph_norm));
eff_50 = fminsearch(delta_50, 0.5)

figure
hold on
plot(energy, counts_dat_50GeV_amorph_norm_1 - counts_dat_50GeV_bg_norm_1,'o')
plot(energy, eff_50.*counts_sim_50GeV_amorph_norm)
