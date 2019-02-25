clear all; clc; close all
datpath = '/home/christian/Documents/cern2018/simdata/';

test = load('../../../Cern2018Experiment/spectre/sum_angles2mmfullsubstinobfac.txt');

trapz(test(2:end,1), test(2:end,2)./test(2:end,1))

%% load data
input_ang   = load('../beamParameters/2017beam_direction_50GeV_2mm.txt');
input_pos   = load('../beamParameters/2017beam_position_50GeV_2mm.txt'); 
sampled     = loaddat('hits_coord_data_50GeV_2mm_aligned_2017');

xpos = input_pos(:,3);
ypos = input_pos(:,4);
xang = input_ang(:,3);
yang = input_ang(:,4);

[xpc, xp] = histcounts(sampled.plane0(:,1), 1000);
[ypc, yp] = histcounts(sampled.plane0(:,2), 1000);

figure
hold on
plot(xp(1:end-1),xpc/max(xpc),'-')
plot(xpos, input_pos(:,1)/max(input_pos(:,1)))

figure
hold on
plot(yp(1:end-1),ypc/max(ypc),'-')
plot(ypos, input_pos(:,2)/max(input_pos(:,2)))

%%
energy = linspace(0,50,50);

% data
filepath1 = strcat('/home/christian/Dropbox/Cern2018Experiment/spectre/energy_2017_50GeV_2mm.txt');
nrg_50GeV_dat1 = load(filepath1) * 6.2415091E9;
[counts_dat_50GeV_amorph_1, ~] = hist(nrg_50GeV_dat1(nrg_50GeV_dat1 < 50 & nrg_50GeV_dat1 > 0), energy);
counts_dat_50GeV_aligned = energy.*counts_dat_50GeV_amorph_1/215529;

filepath1 = strcat('/home/christian/Dropbox/Cern2018Experiment/spectre/energy_2017_50GeV_2mm_bg.txt');
nrg_50GeV_dat1 = load(filepath1) * 6.2415091E9;
[counts_dat_50GeV_amorph_1, ~] = hist(nrg_50GeV_dat1(nrg_50GeV_dat1 < 50 & nrg_50GeV_dat1 > 0), energy);
counts_dat_50GeV_bg = energy.*counts_dat_50GeV_amorph_1/2263230;

filepath1 = strcat('/home/christian/Dropbox/Cern2018Experiment/spectre/energy_2017_50GeV_2mm_amorph.txt');
nrg_50GeV_dat1 = load(filepath1) * 6.2415091E9;
[counts_dat_50GeV_amorph_1, ~] = hist(nrg_50GeV_dat1(nrg_50GeV_dat1 < 50 & nrg_50GeV_dat1 > 0), energy);
counts_dat_50GeV_amorphous = energy.*counts_dat_50GeV_amorph_1/471216;

filepath1 = strcat('/home/christian/Dropbox/Cern2018Experiment/spectre/cut_amorf_2mm_momenta.txt');
nrg_50GeV_dat1 = load(filepath1);
[counts_dat_50GeV_amorph_1, ~] = hist(nrg_50GeV_dat1(nrg_50GeV_dat1 < 50 & nrg_50GeV_dat1 > 0), energy);
counts_dat_50GeV_amorphous_flohr = energy.*counts_dat_50GeV_amorph_1/1.8e6;

%sim 
filepath1 = strcat('/home/christian/Dropbox/Cern2018Experiment/spectre/energy_50GeV_2mm_aligned_2017_sim.txt');
nrg_50GeV_dat1 = load(filepath1) * 6.2415091E9;
[counts_dat_50GeV_amorph_1, ~] = hist(nrg_50GeV_dat1(nrg_50GeV_dat1 < 50 & nrg_50GeV_dat1 > 0), energy);
counts_dat_50GeV_aligned_sim = energy.*counts_dat_50GeV_amorph_1/1e6;

filepath1 = strcat('/home/christian/Dropbox/Cern2018Experiment/spectre/energy_50GeV_2mm_amorphous_2017_sim.txt');
nrg_50GeV_dat1 = load(filepath1) * 6.2415091E9;
[counts_dat_50GeV_amorph_1, ~] = hist(nrg_50GeV_dat1(nrg_50GeV_dat1 < 50 & nrg_50GeV_dat1 > 0), energy);
counts_dat_50GeV_amorphous_sim = energy.*counts_dat_50GeV_amorph_1/1e6;

filepath1 = strcat('/home/christian/Dropbox/Cern2018Experiment/spectre/energy_50GeV_2mm_background_2017_sim.txt');
nrg_50GeV_dat1 = load(filepath1) * 6.2415091E9;
[counts_dat_50GeV_amorph_1, ~] = hist(nrg_50GeV_dat1(nrg_50GeV_dat1 < 50 & nrg_50GeV_dat1 > 0), energy);
counts_dat_50GeV_background_sim = energy.*counts_dat_50GeV_amorph_1/1e6;

filepath1 = strcat('/home/christian/Dropbox/Cern2018Experiment/spectre/cutnobfac_substi_align_2mm_momenta_sim.txt');
nrg_50GeV_dat1 = load(filepath1);
[counts_dat_50GeV_amorph_1, ~] = hist(nrg_50GeV_dat1(nrg_50GeV_dat1 < 50 & nrg_50GeV_dat1 > 0), energy);
counts_dat_50GeV_aligned_sim_flohr_wbg = energy.*counts_dat_50GeV_amorph_1/1e6;

filepath1 = strcat('/home/christian/Dropbox/Cern2018Experiment/spectre/cutnobfac_amorf_2mm_momenta_sim.txt');
nrg_50GeV_dat1 = load(filepath1);
[counts_dat_50GeV_amorph_1, ~] = hist(nrg_50GeV_dat1(nrg_50GeV_dat1 < 50 & nrg_50GeV_dat1 > 0), energy);
counts_dat_50GeV_amorph_sim_flohr_wbg = energy.*counts_dat_50GeV_amorph_1/41e6;

figure
hold on
plot(energy, 0.8 * counts_dat_50GeV_background_sim)
plot(energy, 0.8 * counts_dat_50GeV_amorphous_sim);
plot(energy, 0.8 * counts_dat_50GeV_amorph_sim_flohr_wbg);
plot(energy, counts_dat_50GeV_amorphous,'x')
plot(energy, counts_dat_50GeV_amorphous_flohr,'x')
plot(energy, counts_dat_50GeV_bg,'o')


figure
hold on
plot(energy, counts_dat_50GeV_aligned ,'x')
plot(energy, counts_dat_50GeV_aligned - counts_dat_50GeV_bg,'o')
plot(energy, 0.8 * counts_dat_50GeV_aligned_sim,'-')
plot(energy, 0.8 * counts_dat_50GeV_aligned_sim_flohr_wbg,'--')

%% functions
function hitdat = loaddat(filename)
    datpath = '/home/christian/Documents/cern2018/simdata/';
    filepath = strcat(datpath,filename,'.txt');
    formatSpec = '%f %f';
    fileID = fopen(filepath,'r');
    datblocks = 6;
    for i = 1:datblocks
        field = strcat('plane',num2str(i-1));
        data = textscan(fileID,formatSpec,'HeaderLines',1);
        hitdat.(field) = [data{1},data{2}];
    end
end
