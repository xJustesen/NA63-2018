clc; clear all; close all;
datpath = '/home/christian/Documents/cern2018/simdata/';

%% DATA 20GeV 1mm e- amorphous
energy = linspace(0, 20, 20);

filepath1 = strcat(datpath, 'energy_109.txt');
filepath2 = strcat(datpath, 'energy_115.txt');

nrg_20GeV_dat1 = load(filepath1) * 6.2415091E9;
nrg_20GeV_dat2 = load(filepath2) * 6.2415091E9;
nrg_20GeV_dat_amorph_tot = [nrg_20GeV_dat1; nrg_20GeV_dat2];

NEvents_20GeV_amorph_tot = 2764454 + 82089;

[counts_dat_20GeV_amorph_1, ~] = hist(nrg_20GeV_dat1(nrg_20GeV_dat1 < 20 & nrg_20GeV_dat1 > 0), energy);
[counts_dat_20GeV_amorph_2, ~] = hist(nrg_20GeV_dat2(nrg_20GeV_dat2 < 20 & nrg_20GeV_dat2 > 0), energy);

counts_dat_20GeV_amorph_norm_1 = energy.*counts_dat_20GeV_amorph_1/2764454;
counts_dat_20GeV_amorph_norm_2 = energy.*counts_dat_20GeV_amorph_2/82089;

[counts_dat_20GeV_amorph_tot, ~] = hist(nrg_20GeV_dat_amorph_tot(nrg_20GeV_dat_amorph_tot < 20 & nrg_20GeV_dat_amorph_tot > 0), energy);
counts_dat_20GeV_amorph_norm_tot = energy.*counts_dat_20GeV_amorph_tot/NEvents_20GeV_amorph_tot;

%% DATA 20GeV 1mm e- background
filepath1 = strcat(datpath, 'energy_105.txt');
filepath2 = strcat(datpath, 'energy_106.txt');
filepath3 = strcat(datpath, 'energy_107.txt');
filepath4 = strcat(datpath, 'energy_108.txt');
filepath5 = strcat(datpath, 'energy_110.txt');
filepath6 = strcat(datpath, 'energy_111.txt');

nrg_20GeV_dat1 = load(filepath1) * 6.2415091E9;
nrg_20GeV_dat2 = load(filepath2) * 6.2415091E9;
nrg_20GeV_dat3 = load(filepath3) * 6.2415091E9;
nrg_20GeV_dat4 = load(filepath4) * 6.2415091E9;
nrg_20GeV_dat5 = load(filepath5) * 6.2415091E9;
nrg_20GeV_dat6 = load(filepath6) * 6.2415091E9;

nrg_20GeV_dat_bg_tot = [nrg_20GeV_dat1; nrg_20GeV_dat2; nrg_20GeV_dat3; nrg_20GeV_dat4; nrg_20GeV_dat5; nrg_20GeV_dat6];
NEvents_20GeV_dat_bg_tot = 324000 + 590172 + 624602 + 734446 + 1415716 + 1224254;

[counts_dat_20GeV_bg_tot, ~] = hist(nrg_20GeV_dat_bg_tot(nrg_20GeV_dat_bg_tot < 20 & nrg_20GeV_dat_bg_tot > 0), energy);
counts_dat_20GeV_bg_norm_tot = energy.*counts_dat_20GeV_bg_tot/NEvents_20GeV_dat_bg_tot;

[counts_dat_20GeV_bg_105, ~] = hist(nrg_20GeV_dat1(nrg_20GeV_dat1 < 20 & nrg_20GeV_dat1 > 0), energy);
counts_dat_20GeV_bg_norm_105 = energy.*counts_dat_20GeV_bg_105/324000;

%% DATA 20GeV 1mm e- aligned.
filepath1 = strcat(datpath,'energy_103.txt');
filepath2 = strcat(datpath,'energy_104.txt');
filepath3 = strcat(datpath,'energy_112.txt');
filepath4 = strcat(datpath,'energy_113.txt');
filepath5 = strcat(datpath,'energy_114.txt');

nrg_20GeV_dat1 = load(filepath1) * 6.2415091E9;
nrg_20GeV_dat2 = load(filepath2) * 6.2415091E9;
nrg_20GeV_dat3 = load(filepath3) * 6.2415091E9;
nrg_20GeV_dat4 = load(filepath4) * 6.2415091E9;
nrg_20GeV_dat5 = load(filepath5) * 6.2415091E9;

NEvents_20GeV_aligned_tot = 876217 + 467348 + 286232 + 258194 + 74259;
nrg_20GeV_dat_aligned_tot = [nrg_20GeV_dat1; nrg_20GeV_dat2; nrg_20GeV_dat3; nrg_20GeV_dat4; nrg_20GeV_dat5];
[counts_dat_20GeV_aligned_tot, ~] = hist(nrg_20GeV_dat_aligned_tot(nrg_20GeV_dat_aligned_tot < 20 & nrg_20GeV_dat_aligned_tot > 0), energy);

counts_dat_20GeV_aligned_norm_tot = energy.*counts_dat_20GeV_aligned_tot/NEvents_20GeV_aligned_tot;

%% DATA 20GeV 1.5mm e- Amorph
filepath1 = strcat(datpath, 'energy_66.txt');
filepath2 = strcat(datpath, 'energy_67.txt');
filepath3 = strcat(datpath, 'energy_68.txt');
filepath4 = strcat(datpath, 'energy_69.txt');

nrg_20GeV_dat1 = load(filepath1) * 6.2415091E9;
nrg_20GeV_dat2 = load(filepath2) * 6.2415091E9;
nrg_20GeV_dat3 = load(filepath3) * 6.2415091E9;
nrg_20GeV_dat4 = load(filepath4) * 6.2415091E9;

nrg_20GeV_dat_amorph_1_5mm_tot = [nrg_20GeV_dat1; nrg_20GeV_dat2; nrg_20GeV_dat3; nrg_20GeV_dat4];
NEvents_20GeV_dat_amorph_1_5mm_tot = 1197107 + 432597 + 634614 + 386867;

[counts_dat_20GeV_amorph_1_5mm_tot, ~] = hist(nrg_20GeV_dat_amorph_1_5mm_tot(nrg_20GeV_dat_amorph_1_5mm_tot < 20 & nrg_20GeV_dat_amorph_1_5mm_tot > 0), energy);
counts_dat_20GeV_amorph_1_5mm_norm_tot = energy.*counts_dat_20GeV_amorph_1_5mm_tot/NEvents_20GeV_dat_amorph_1_5mm_tot;

[counts_dat_20GeV_amorph_1_5mm_1, ~] = hist(nrg_20GeV_dat1(nrg_20GeV_dat1 < 20 & nrg_20GeV_dat1 > 0), energy);
[counts_dat_20GeV_amorph_1_5mm_2, ~] = hist(nrg_20GeV_dat2(nrg_20GeV_dat2 < 10000 & nrg_20GeV_dat2 > 0), energy);
[counts_dat_20GeV_amorph_1_5mm_3, ~] = hist(nrg_20GeV_dat3(nrg_20GeV_dat3 < 20 & nrg_20GeV_dat3 > 0), energy);
[counts_dat_20GeV_amorph_1_5mm_4, ~] = hist(nrg_20GeV_dat4(nrg_20GeV_dat4 < 20 & nrg_20GeV_dat4 > 0), energy);

counts_dat_20GeV_amorph_1_5mm_1_norm = energy.* counts_dat_20GeV_amorph_1_5mm_1/1197107;
counts_dat_20GeV_amorph_1_5mm_2_norm = energy.* counts_dat_20GeV_amorph_1_5mm_2/432597;
counts_dat_20GeV_amorph_1_5mm_3_norm = energy.* counts_dat_20GeV_amorph_1_5mm_3/634614;
counts_dat_20GeV_amorph_1_5mm_4_norm = energy.* counts_dat_20GeV_amorph_1_5mm_4/386867;

%% DATA 20GeV 1.5mm e- Background
filepath1 = strcat(datpath, 'energy_60.txt');
filepath2 = strcat(datpath, 'energy_65.txt');

nrg_20GeV_dat1 = load(filepath1) * 6.2415091E9;
nrg_20GeV_dat2 = load(filepath2) * 6.2415091E9;

nrg_20GeV_dat_bg_1_5mm_tot = [nrg_20GeV_dat1; nrg_20GeV_dat2];
NEvents_20GeV_dat_bg_1_5mm_tot = 812417 + 2530699;

[counts_dat_20GeV_bg_1_5mm_tot, ~] = hist(nrg_20GeV_dat_bg_1_5mm_tot(nrg_20GeV_dat_bg_1_5mm_tot < 20 & nrg_20GeV_dat_bg_1_5mm_tot > 0), energy);
counts_dat_20GeV_bg_1_5mm_norm_tot = energy .* counts_dat_20GeV_bg_1_5mm_tot/NEvents_20GeV_dat_bg_1_5mm_tot;

[counts_dat_20GeV_bg_1_5mm_1, ~] = hist(nrg_20GeV_dat1(nrg_20GeV_dat1 < 20 & nrg_20GeV_dat1 > 0), energy);
[counts_dat_20GeV_bg_1_5mm_2, ~] = hist(nrg_20GeV_dat2(nrg_20GeV_dat2 < 20 & nrg_20GeV_dat2 > 0), energy);

counts_dat_20GeV_bg_1_5mm_norm_1 = counts_dat_20GeV_bg_1_5mm_1/812417;
counts_dat_20GeV_bg_1_5mm_norm_2 = counts_dat_20GeV_bg_1_5mm_2/2530699;

%% DATA 20GeV 1.5mm e- aligned.
filepath1 = strcat(datpath,'energy_61.txt');
filepath2 = strcat(datpath,'energy_62.txt');
filepath3 = strcat(datpath,'energy_63.txt');
filepath4 = strcat(datpath,'energy_64.txt');

nrg_20GeV_dat1 = load(filepath1) * 6.2415091E9;
nrg_20GeV_dat2 = load(filepath2) * 6.2415091E9;
nrg_20GeV_dat3 = load(filepath3) * 6.2415091E9;
nrg_20GeV_dat4 = load(filepath4) * 6.2415091E9;

NEvents_20GeV_aligned_tot_1_5mm = 174257 + 474970 + 876508 + 792574;
nrg_20GeV_dat_aligned_tot_1_5mm = [nrg_20GeV_dat1; nrg_20GeV_dat2; nrg_20GeV_dat3; nrg_20GeV_dat4];
[counts_dat_20GeV_aligned_tot_1_5mm, ~] = hist(nrg_20GeV_dat_aligned_tot_1_5mm(nrg_20GeV_dat_aligned_tot_1_5mm < 20 & nrg_20GeV_dat_aligned_tot_1_5mm > 0), energy);

[counts_dat_20GeV_aligned_1_1_5mm, ~] = hist(nrg_20GeV_dat1(nrg_20GeV_dat1 < 20 & nrg_20GeV_dat1 > 0), energy);
[counts_dat_20GeV_aligned_2_1_5mm, ~] = hist(nrg_20GeV_dat2(nrg_20GeV_dat2 < 20 & nrg_20GeV_dat2 > 0), energy);
[counts_dat_20GeV_aligned_3_1_5mm, ~] = hist(nrg_20GeV_dat3(nrg_20GeV_dat3 < 20 & nrg_20GeV_dat3 > 0), energy);
[counts_dat_20GeV_aligned_4_1_5mm, ~] = hist(nrg_20GeV_dat4(nrg_20GeV_dat4 < 20 & nrg_20GeV_dat4 > 0), energy);

counts_dat_20GeV_aligned_norm_tot_1_5mm = energy.*counts_dat_20GeV_aligned_tot_1_5mm/NEvents_20GeV_aligned_tot_1_5mm;

%% SIM 20GeV 1mm e- amorph + backogrund
filepath1 = strcat(datpath,'energy_sim_amorphous1_20GeV_16012019.txt'); 
filepath2 = strcat(datpath,'energy_sim_amorphous2_20GeV_16012019.txt'); 
filepath3 = strcat(datpath,'energy_sim_amorphous3_20GeV_16012019.txt'); 
filepath4 =strcat(datpath,'energy_sim_amorphous4_20GeV_16012019.txt'); 
filepath5 =strcat(datpath,'energy_sim_amorphous5_20GeV_16012019.txt');

nrg_20GeV_sim1 = load(filepath1) * 6.2415091E9;
nrg_20GeV_sim2 = load(filepath2) * 6.2415091E9; 
nrg_20GeV_sim3 = load(filepath3) * 6.2415091E9; 
nrg_20GeV_sim4 = load(filepath4) * 6.2415091E9;
nrg_20GeV_sim5 = load(filepath5) * 6.2415091E9;

nrg_20GeV_sim_amorph_bg = [nrg_20GeV_sim1; nrg_20GeV_sim2; nrg_20GeV_sim3; nrg_20GeV_sim4; nrg_20GeV_sim5];
[counts_sim_20GeV_amorph_bg, ~] = hist(nrg_20GeV_sim_amorph_bg(nrg_20GeV_sim_amorph_bg < 20 & nrg_20GeV_sim_amorph_bg > 0), energy); 
counts_sim_20GeV_amorph_bg_norm = energy.*counts_sim_20GeV_amorph_bg/5e6;
 
%% SIM 20GeV 1mm e- amorph
filepath1 = strcat(datpath,'energy_sim_amorphous1_20GeV_no_background.txt');
filepath2 = strcat(datpath,'energy_sim_amorphous2_20GeV_no_background.txt');
filepath3 = strcat(datpath,'energy_sim_amorphous3_20GeV_no_background.txt');
filepath4 = strcat(datpath,'energy_sim_amorphous4_20GeV_no_background.txt');
filepath5 = strcat(datpath,'energy_sim_amorphous5_20GeV_no_background.txt');

nrg_20GeV_sim1 = load(filepath1) * 6.2415091E9;
nrg_20GeV_sim2 = load(filepath2) * 6.2415091E9;
nrg_20GeV_sim3 = load(filepath3) * 6.2415091E9;
nrg_20GeV_sim4 = load(filepath4) * 6.2415091E9;
nrg_20GeV_sim5 = load(filepath5) * 6.2415091E9;

nrg_20GeV_sim_amorph = [nrg_20GeV_sim1; nrg_20GeV_sim2; nrg_20GeV_sim3; nrg_20GeV_sim4; nrg_20GeV_sim5];
[counts_sim_20GeV_amorph, ~] = hist(nrg_20GeV_sim_amorph(nrg_20GeV_sim_amorph < 20 & nrg_20GeV_sim_amorph > 0), energy);
counts_sim_20GeV_amorph_norm = energy.*counts_sim_20GeV_amorph/5e6;

%% SIM 20GeV 1mm e- background
filepath1 = strcat(datpath,'energy_sim_background1_20GeV_16012019.txt');
filepath2 = strcat(datpath,'energy_sim_background2_20GeV_16012019.txt');
filepath3 = strcat(datpath,'energy_sim_background3_20GeV_16012019.txt');
filepath4 = strcat(datpath,'energy_sim_background4_20GeV_16012019.txt');
filepath5 = strcat(datpath,'energy_sim_background5_20GeV_16012019.txt');

nrg_20GeV_sim1 = load(filepath1) * 6.2415091E9;
nrg_20GeV_sim2 = load(filepath2) * 6.2415091E9;
nrg_20GeV_sim3 = load(filepath3) * 6.2415091E9;
nrg_20GeV_sim4 = load(filepath4) * 6.2415091E9;
nrg_20GeV_sim5 = load(filepath5) * 6.2415091E9;

nrg_20GeV_sim_bg = [nrg_20GeV_sim1; nrg_20GeV_sim2; nrg_20GeV_sim3; nrg_20GeV_sim4; nrg_20GeV_sim5];
[counts_sim_20GeV_bg, ~] = hist(nrg_20GeV_sim_bg(nrg_20GeV_sim_bg < 20 & nrg_20GeV_sim_bg > 0), energy);
counts_sim_20GeV_bg_norm = energy.*counts_sim_20GeV_bg/5e6;

%% SIM 20GeV 1mm e- aligned w/o schot
filepath1 = strcat(datpath,'energy_sim_aligned1_20GeV_woshot.txt');
filepath2 = strcat(datpath,'energy_sim_aligned2_20GeV_woshot.txt');
filepath3 = strcat(datpath,'energy_sim_aligned3_20GeV_woshot.txt');
filepath4 = strcat(datpath,'energy_sim_aligned4_20GeV_woshot.txt');
filepath5 = strcat(datpath,'energy_sim_aligned5_20GeV_woshot.txt');

nrg_20GeV_sim1 = load(filepath1) * 6.2415091E9;
nrg_20GeV_sim2 = load(filepath2) * 6.2415091E9;
nrg_20GeV_sim3 = load(filepath3) * 6.2415091E9;
nrg_20GeV_sim4 = load(filepath4) * 6.2415091E9;
nrg_20GeV_sim5 = load(filepath5) * 6.2415091E9;

nrg_20GeV_sim_aligned = [nrg_20GeV_sim1; nrg_20GeV_sim2; nrg_20GeV_sim3; nrg_20GeV_sim4; nrg_20GeV_sim5];
[counts_sim_20GeV_aligned_woshot, ~] = hist(nrg_20GeV_sim_aligned(nrg_20GeV_sim_aligned < 20 & nrg_20GeV_sim_aligned > 0), energy);
counts_sim_20GeV_aligned_norm_woshot = energy.*counts_sim_20GeV_aligned_woshot/5e6;

%% SIM 20GeV 1mm e- aligned w/o RR
filepath1 = strcat(datpath,'energy_sim_aligned1_20GeV_worr.txt');
filepath2 = strcat(datpath,'energy_sim_aligned2_20GeV_worr.txt');
filepath3 = strcat(datpath,'energy_sim_aligned3_20GeV_worr.txt');
filepath4 = strcat(datpath,'energy_sim_aligned4_20GeV_worr.txt');
filepath5 = strcat(datpath,'energy_sim_aligned5_20GeV_worr.txt');

nrg_20GeV_sim1 = load(filepath1) * 6.2415091E9;
nrg_20GeV_sim2 = load(filepath2) * 6.2415091E9;
nrg_20GeV_sim3 = load(filepath3) * 6.2415091E9;
nrg_20GeV_sim4 = load(filepath4) * 6.2415091E9;
nrg_20GeV_sim5 = load(filepath5) * 6.2415091E9;

nrg_20GeV_sim_aligned = [nrg_20GeV_sim1; nrg_20GeV_sim2; nrg_20GeV_sim3; nrg_20GeV_sim4; nrg_20GeV_sim5];
[counts_sim_20GeV_aligned_worr, ~] = hist(nrg_20GeV_sim_aligned(nrg_20GeV_sim_aligned < 20 & nrg_20GeV_sim_aligned > 0), energy);
counts_sim_20GeV_aligned_norm_worr = energy.*counts_sim_20GeV_aligned_worr/5e6;

%% SIM 20GeV 1mm e- aligned full LL
filepath1 = strcat(datpath,'energy_sim_aligned1_20GeV.txt');
filepath2 = strcat(datpath,'energy_sim_aligned2_20GeV.txt');
filepath3 = strcat(datpath,'energy_sim_aligned3_20GeV.txt');
filepath4 = strcat(datpath,'energy_sim_aligned4_20GeV.txt');
filepath5 = strcat(datpath,'energy_sim_aligned5_20GeV.txt');

nrg_20GeV_sim1 = load(filepath1) * 6.2415091E9;
nrg_20GeV_sim2 = load(filepath2) * 6.2415091E9;
nrg_20GeV_sim3 = load(filepath3) * 6.2415091E9;
nrg_20GeV_sim4 = load(filepath4) * 6.2415091E9;
nrg_20GeV_sim5 = load(filepath5) * 6.2415091E9;

nrg_20GeV_sim_aligned = [nrg_20GeV_sim1; nrg_20GeV_sim2; nrg_20GeV_sim3; nrg_20GeV_sim4; nrg_20GeV_sim5];
[counts_sim_20GeV_aligned, ~] = hist(nrg_20GeV_sim_aligned(nrg_20GeV_sim_aligned < 20 & nrg_20GeV_sim_aligned > 0), energy);
counts_sim_20GeV_aligned_norm = energy.*counts_sim_20GeV_aligned/5e6;

%% SIM 20GeV 1.5mm e- amorph
filepath1 = strcat(datpath,'energy_sim_amorphous1_20GeV_1.5mm_no_background.txt');
filepath2 = strcat(datpath,'energy_sim_amorphous2_20GeV_1.5mm_no_background.txt');
filepath3 = strcat(datpath,'energy_sim_amorphous3_20GeV_1.5mm_no_background.txt');
filepath4 = strcat(datpath,'energy_sim_amorphous4_20GeV_1.5mm_no_background.txt');
filepath5 = strcat(datpath,'energy_sim_amorphous5_20GeV_1.5mm_no_background.txt');

nrg_20GeV_sim1 = load(filepath1) * 6.2415091E9;
nrg_20GeV_sim2 = load(filepath2) * 6.2415091E9;
nrg_20GeV_sim3 = load(filepath3) * 6.2415091E9;
nrg_20GeV_sim4 = load(filepath4) * 6.2415091E9;
nrg_20GeV_sim5 = load(filepath5) * 6.2415091E9;

nrg_20GeV_sim_amorph_1_5mm = [nrg_20GeV_sim1; nrg_20GeV_sim2; nrg_20GeV_sim3; nrg_20GeV_sim4; nrg_20GeV_sim5];
[counts_sim_20GeV_amorph_1_5mm, ~] = hist(nrg_20GeV_sim_amorph_1_5mm(nrg_20GeV_sim_amorph_1_5mm < 20 & nrg_20GeV_sim_amorph_1_5mm > 0), energy);
counts_sim_20GeV_amorph_1_5mm_norm = energy.*counts_sim_20GeV_amorph_1_5mm/5e6;

%% SIM 20GeV 1.5mm e- background
filepath1 = strcat(datpath,'energy_sim_background1_20GeV_1.5mm_16012019.txt');
filepath2 = strcat(datpath,'energy_sim_background2_20GeV_1.5mm_16012019.txt');
filepath3 = strcat(datpath,'energy_sim_background3_20GeV_1.5mm_16012019.txt');
filepath4 = strcat(datpath,'energy_sim_background4_20GeV_1.5mm_16012019.txt');
filepath5 = strcat(datpath,'energy_sim_background5_20GeV_1.5mm_16012019.txt');

nrg_20GeV_sim1 = load(filepath1) * 6.2415091E9;
nrg_20GeV_sim2 = load(filepath2) * 6.2415091E9;
nrg_20GeV_sim3 = load(filepath3) * 6.2415091E9;
nrg_20GeV_sim4 = load(filepath4) * 6.2415091E9;
nrg_20GeV_sim5 = load(filepath5) * 6.2415091E9;

nrg_20GeV_sim_bg_1_5mm = [nrg_20GeV_sim1; nrg_20GeV_sim2; nrg_20GeV_sim3; nrg_20GeV_sim4; nrg_20GeV_sim5];
[counts_sim_20GeV_bg_1_5mm, ~] = hist(nrg_20GeV_sim_bg_1_5mm(nrg_20GeV_sim_bg_1_5mm < 20 & nrg_20GeV_sim_bg_1_5mm > 0), energy);
counts_sim_20GeV_bg_norm_1_5mm = energy.*counts_sim_20GeV_bg_1_5mm/5e6;

%% SIM 20GeV 1.5mm e- amorph + backogrund
filepath1 = strcat(datpath,'energy_sim_amorphous1_20GeV_1.5mm_16012019.txt');
filepath2 = strcat(datpath,'energy_sim_amorphous2_20GeV_1.5mm_16012019.txt');
filepath3 = strcat(datpath,'energy_sim_amorphous3_20GeV_1.5mm_16012019.txt');
filepath4 = strcat(datpath,'energy_sim_amorphous4_20GeV_1.5mm_16012019.txt');
filepath5 = strcat(datpath,'energy_sim_amorphous5_20GeV_1.5mm_16012019.txt');

nrg_20GeV_sim1 = load(filepath1) * 6.2415091E9;
nrg_20GeV_sim2 = load(filepath2) * 6.2415091E9;
nrg_20GeV_sim3 = load(filepath3) * 6.2415091E9;
nrg_20GeV_sim4 = load(filepath4) * 6.2415091E9;
nrg_20GeV_sim5 = load(filepath5) * 6.2415091E9;

nrg_20GeV_sim_amorph_bg_1_5mm = [nrg_20GeV_sim1; nrg_20GeV_sim2; nrg_20GeV_sim3; nrg_20GeV_sim4; nrg_20GeV_sim5];
[counts_sim_20GeV_amorph_bg_1_5mm, ~] = hist(nrg_20GeV_sim_amorph_bg_1_5mm(nrg_20GeV_sim_amorph_bg_1_5mm < 20 & nrg_20GeV_sim_amorph_bg_1_5mm > 0), energy);
counts_sim_20GeV_amorph_bg_norm_1_5mm = energy.*counts_sim_20GeV_amorph_bg_1_5mm/5e6;

%% SIM 20GeV 1.5mm e- aligned w/o schot
filepath1 = strcat(datpath,'energy_sim_aligned1_20GeV_woshot_1.5mm.txt');
filepath2 = strcat(datpath,'energy_sim_aligned2_20GeV_woshot_1.5mm.txt');
filepath3 = strcat(datpath,'energy_sim_aligned3_20GeV_woshot_1.5mm.txt');
filepath4 = strcat(datpath,'energy_sim_aligned4_20GeV_woshot_1.5mm.txt');
filepath5 = strcat(datpath,'energy_sim_aligned5_20GeV_woshot_1.5mm.txt');

nrg_20GeV_sim1 = load(filepath1) * 6.2415091E9;
nrg_20GeV_sim2 = load(filepath2) * 6.2415091E9;
nrg_20GeV_sim3 = load(filepath3) * 6.2415091E9;
nrg_20GeV_sim4 = load(filepath4) * 6.2415091E9;
nrg_20GeV_sim5 = load(filepath5) * 6.2415091E9;

nrg_20GeV_sim_aligned = [nrg_20GeV_sim1; nrg_20GeV_sim2; nrg_20GeV_sim3; nrg_20GeV_sim4; nrg_20GeV_sim5];
[counts_sim_20GeV_aligned_1_5mm_woshot, ~] = hist(nrg_20GeV_sim_aligned(nrg_20GeV_sim_aligned < 20 & nrg_20GeV_sim_aligned > 0), energy);
counts_sim_20GeV_aligned_1_5mm_norm_woshot = energy.*counts_sim_20GeV_aligned_1_5mm_woshot/5e6;

%% SIM 20GeV 1.5mm e- aligned w/o RR
filepath1 = strcat(datpath,'energy_sim_aligned1_20GeV_worr_1.5mm.txt');
filepath2 = strcat(datpath,'energy_sim_aligned2_20GeV_worr_1.5mm.txt');
filepath3 = strcat(datpath,'energy_sim_aligned3_20GeV_worr_1.5mm.txt');
filepath4 = strcat(datpath,'energy_sim_aligned4_20GeV_worr_1.5mm.txt');
filepath5 = strcat(datpath,'energy_sim_aligned5_20GeV_worr_1.5mm.txt');

nrg_20GeV_sim1 = load(filepath1) * 6.2415091E9;
nrg_20GeV_sim2 = load(filepath2) * 6.2415091E9;
nrg_20GeV_sim3 = load(filepath3) * 6.2415091E9;
nrg_20GeV_sim4 = load(filepath4) * 6.2415091E9;
nrg_20GeV_sim5 = load(filepath5) * 6.2415091E9;

nrg_20GeV_sim_aligned = [nrg_20GeV_sim1; nrg_20GeV_sim2; nrg_20GeV_sim3; nrg_20GeV_sim4; nrg_20GeV_sim5];
[counts_sim_20GeV_aligned_1_5mm_worr, ~] = hist(nrg_20GeV_sim_aligned(nrg_20GeV_sim_aligned < 20 & nrg_20GeV_sim_aligned > 0), energy);
counts_sim_20GeV_aligned_1_5mm_norm_worr = energy.*counts_sim_20GeV_aligned_1_5mm_worr/5e6;

%% SIM 20GeV 1.5mm e- aligned full LL
filepath1 = strcat(datpath,'energy_sim_aligned1_20GeV_1.5mm.txt');
filepath2 = strcat(datpath,'energy_sim_aligned2_20GeV_1.5mm.txt');
filepath3 = strcat(datpath,'energy_sim_aligned3_20GeV_1.5mm.txt');
filepath4 = strcat(datpath,'energy_sim_aligned4_20GeV_1.5mm.txt');
filepath5 = strcat(datpath,'energy_sim_aligned5_20GeV_1.5mm.txt');

nrg_20GeV_sim1 = load(filepath1) * 6.2415091E9;
nrg_20GeV_sim2 = load(filepath2) * 6.2415091E9;
nrg_20GeV_sim3 = load(filepath3) * 6.2415091E9;
nrg_20GeV_sim4 = load(filepath4) * 6.2415091E9;
nrg_20GeV_sim5 = load(filepath5) * 6.2415091E9;

nrg_20GeV_sim_aligned = [nrg_20GeV_sim1; nrg_20GeV_sim2; nrg_20GeV_sim3; nrg_20GeV_sim4; nrg_20GeV_sim5];
[counts_sim_20GeV_aligned_1_5mm, ~] = hist(nrg_20GeV_sim_aligned(nrg_20GeV_sim_aligned < 20 & nrg_20GeV_sim_aligned > 0), energy);
counts_sim_20GeV_aligned_1_5mm_norm = energy.*counts_sim_20GeV_aligned_1_5mm/5e6;

%% DATA 40GeV 1mm e- background
energy = linspace(0, 40, 40);

filepath1 = strcat(datpath,'energy_74.txt');

nrg_40GeV_dat_bg_1 = load(filepath1) * 6.2415091E9;

[counts_dat_40GeV_bg_1, ~] = hist(nrg_40GeV_dat_bg_1(nrg_40GeV_dat_bg_1 < 40 & nrg_40GeV_dat_bg_1 > 0), energy);
nwocuts = 3029506;
counts_dat_40GeV_bg_norm = energy .* counts_dat_40GeV_bg_1/nwocuts;

%% DATA 40GeV 1mm e- amoprh.
filepath1 = strcat(datpath,'energy_73.txt');
filepath2 = strcat(datpath,'energy_75.txt');
filepath3 = strcat(datpath,'energy_76.txt');
filepath4 = strcat(datpath,'energy_77.txt');
filepath5 = strcat(datpath,'energy_78.txt');

nrg_40GeV_dat1 = load(filepath1) * 6.2415091E9;
nrg_40GeV_dat2 = load(filepath2) * 6.2415091E9;
nrg_40GeV_dat3 = load(filepath3) * 6.2415091E9;
nrg_40GeV_dat4 = load(filepath4) * 6.2415091E9;
nrg_40GeV_dat5 = load(filepath5) * 6.2415091E9;

NEvents_40GeV_amorph_tot = 1290988 + 1361162 + 1447462 + 715126 + 1456319;
nrg_40GeV_dat_amorph_tot = [nrg_40GeV_dat1; nrg_40GeV_dat2; nrg_40GeV_dat3; nrg_40GeV_dat4; nrg_40GeV_dat5];
[counts_dat_40GeV_amorph_tot, ~] = hist(nrg_40GeV_dat_amorph_tot(nrg_40GeV_dat_amorph_tot < 40 & nrg_40GeV_dat_amorph_tot > 0), energy);

[counts_dat_40GeV_amorph_1, ~] = hist(nrg_40GeV_dat1(nrg_40GeV_dat1 < 40 & nrg_40GeV_dat1 > 0), energy);
[counts_dat_40GeV_amorph_2, ~] = hist(nrg_40GeV_dat2(nrg_40GeV_dat2 < 40 & nrg_40GeV_dat2 > 0), energy);
[counts_dat_40GeV_amorph_3, ~] = hist(nrg_40GeV_dat3(nrg_40GeV_dat3 < 40 & nrg_40GeV_dat3 > 0), energy);
[counts_dat_40GeV_amorph_4, ~] = hist(nrg_40GeV_dat4(nrg_40GeV_dat4 < 40 & nrg_40GeV_dat4 > 0), energy);
[counts_dat_40GeV_amorph_5, ~] = hist(nrg_40GeV_dat5(nrg_40GeV_dat5 < 40 & nrg_40GeV_dat5 > 0), energy);

counts_dat_40GeV_amorph_norm_tot = energy .* counts_dat_40GeV_amorph_tot/NEvents_40GeV_amorph_tot;

counts_dat_40GeV_amorph_norm_1 = energy .* counts_dat_40GeV_amorph_1/1412665;
counts_dat_40GeV_amorph_norm_2 = energy .* counts_dat_40GeV_amorph_2/1488348;
counts_dat_40GeV_amorph_norm_3 = energy .* counts_dat_40GeV_amorph_3/1581029;
counts_dat_40GeV_amorph_norm_4 = energy .* counts_dat_40GeV_amorph_4/715126;
counts_dat_40GeV_amorph_norm_5 = energy .* counts_dat_40GeV_amorph_5/1588513;

%% DATA 40GeV 1mm e- aligned
filepath1 = strcat(datpath,'energy_71.txt');
filepath2 = strcat(datpath,'energy_72.txt');
filepath3 = strcat(datpath,'energy_79.txt');
filepath4 = strcat(datpath,'energy_80.txt');
filepath5 = strcat(datpath,'energy_81.txt');

nrg_40GeV_dat1 = load(filepath1) * 6.2415091E9;
nrg_40GeV_dat2 = load(filepath2) * 6.2415091E9;
nrg_40GeV_dat3 = load(filepath3) * 6.2415091E9;
nrg_40GeV_dat4 = load(filepath4) * 6.2415091E9;
nrg_40GeV_dat5 = load(filepath5) * 6.2415091E9;

nrg_40GeV_dat_tot = [nrg_40GeV_dat1; nrg_40GeV_dat2; nrg_40GeV_dat3; nrg_40GeV_dat4; nrg_40GeV_dat5];
NEvents_40GeV_aligned_tot = 142959 + 473324 + 460625 + 1288624 + 1275493;
[counts_dat_40GeV_aligned_tot, ~] = hist(nrg_40GeV_dat_tot(nrg_40GeV_dat_tot < 40 & nrg_40GeV_dat_tot > 0), energy);

[counts_dat_40GeV_aligned_1, ~] = hist(nrg_40GeV_dat1(nrg_40GeV_dat1 < 40 & nrg_40GeV_dat1 > 0), energy);
[counts_dat_40GeV_aligned_2, ~] = hist(nrg_40GeV_dat2(nrg_40GeV_dat2 < 40 & nrg_40GeV_dat2 > 0), energy);
[counts_dat_40GeV_aligned_3, ~] = hist(nrg_40GeV_dat3(nrg_40GeV_dat3 < 40 & nrg_40GeV_dat3 > 0), energy);
[counts_dat_40GeV_aligned_4, ~] = hist(nrg_40GeV_dat4(nrg_40GeV_dat4 < 40 & nrg_40GeV_dat4 > 0), energy);
[counts_dat_40GeV_aligned_5, ~] = hist(nrg_40GeV_dat5(nrg_40GeV_dat5 < 40 & nrg_40GeV_dat5 > 0), energy);

counts_dat_40GeV_aligned_norm_1 = energy .* counts_dat_40GeV_aligned_1/92534;
counts_dat_40GeV_aligned_norm_2 = energy .* counts_dat_40GeV_aligned_2/520239;
counts_dat_40GeV_aligned_norm_3 = energy .* counts_dat_40GeV_aligned_3/504200;
counts_dat_40GeV_aligned_norm_4 = energy .* counts_dat_40GeV_aligned_4/1409253;
counts_dat_40GeV_aligned_norm_5 = energy .* counts_dat_40GeV_aligned_5/1395830;

counts_dat_40GeV_aligned_norm_tot = energy .* counts_dat_40GeV_aligned_tot/NEvents_40GeV_aligned_tot;

%% DATA 40GeV 1.5mm e- amoprh.
filepath1 = strcat(datpath,'energy_32.txt');
filepath2 = strcat(datpath,'energy_33.txt');
filepath3 = strcat(datpath,'energy_34.txt');
filepath4 = strcat(datpath,'energy_39.txt');
filepath5 = strcat(datpath,'energy_40.txt');
filepath6 = strcat(datpath,'energy_41.txt');
filepath7 = strcat(datpath,'energy_43.txt');

nrg_40GeV_dat1 = load(filepath1) * 6.2415091E9;
nrg_40GeV_dat2 = load(filepath2) * 6.2415091E9;
nrg_40GeV_dat3 = load(filepath3) * 6.2415091E9;
nrg_40GeV_dat4 = load(filepath4) * 6.2415091E9;
nrg_40GeV_dat5 = load(filepath5) * 6.2415091E9;
nrg_40GeV_dat6 = load(filepath6) * 6.2415091E9;
nrg_40GeV_dat7 = load(filepath7) * 6.2415091E9;

NEvents_40GeV_amorph_tot_1_5mm = 1449257 + 529097 + 724698 + 134167 + 692475 + 1694966 + 496471;
nrg_40GeV_dat_amorph_tot_1_5mm = [nrg_40GeV_dat1; nrg_40GeV_dat2; nrg_40GeV_dat3; nrg_40GeV_dat4; nrg_40GeV_dat5; nrg_40GeV_dat6; nrg_40GeV_dat7];
[counts_dat_40GeV_amorph_tot_1_5mm, ~] = hist(nrg_40GeV_dat_amorph_tot_1_5mm(nrg_40GeV_dat_amorph_tot_1_5mm < 40 & nrg_40GeV_dat_amorph_tot_1_5mm > 0), energy);

[counts_dat_40GeV_amorph_1_1_5mm, ~] = hist(nrg_40GeV_dat1(nrg_40GeV_dat1 < 40 & nrg_40GeV_dat1 > 0), energy);
[counts_dat_40GeV_amorph_2_1_5mm, ~] = hist(nrg_40GeV_dat2(nrg_40GeV_dat2 < 40 & nrg_40GeV_dat2 > 0), energy);
[counts_dat_40GeV_amorph_3_1_5mm, ~] = hist(nrg_40GeV_dat3(nrg_40GeV_dat3 < 40 & nrg_40GeV_dat3 > 0), energy);
[counts_dat_40GeV_amorph_4_1_5mm, ~] = hist(nrg_40GeV_dat4(nrg_40GeV_dat4 < 40 & nrg_40GeV_dat4 > 0), energy);
[counts_dat_40GeV_amorph_5_1_5mm, ~] = hist(nrg_40GeV_dat5(nrg_40GeV_dat5 < 40 & nrg_40GeV_dat5 > 0), energy);
[counts_dat_40GeV_amorph_6_1_5mm, ~] = hist(nrg_40GeV_dat6(nrg_40GeV_dat6 < 40 & nrg_40GeV_dat6 > 0), energy);
[counts_dat_40GeV_amorph_7_1_5mm, ~] = hist(nrg_40GeV_dat7(nrg_40GeV_dat7 < 40 & nrg_40GeV_dat7 > 0), energy);

counts_dat_40GeV_amorph_norm_tot_1_5mm = energy.*counts_dat_40GeV_amorph_tot_1_5mm/NEvents_40GeV_amorph_tot_1_5mm;

counts_dat_40GeV_amorph_1_1_5mm_norm = energy .* counts_dat_40GeV_amorph_1_1_5mm/1449257;
counts_dat_40GeV_amorph_2_1_5mm_norm = energy .* counts_dat_40GeV_amorph_2_1_5mm/529097;
counts_dat_40GeV_amorph_3_1_5mm_norm = energy .* counts_dat_40GeV_amorph_3_1_5mm/724698;
counts_dat_40GeV_amorph_4_1_5mm_norm = energy .* counts_dat_40GeV_amorph_4_1_5mm/134167;
counts_dat_40GeV_amorph_5_1_5mm_norm = energy .* counts_dat_40GeV_amorph_5_1_5mm/692475;
counts_dat_40GeV_amorph_6_1_5mm_norm = energy .* counts_dat_40GeV_amorph_6_1_5mm/1694966;
counts_dat_40GeV_amorph_7_1_5mm_norm = energy .* counts_dat_40GeV_amorph_7_1_5mm/496471;

%% DATA 40GeV 1.5mm e- background
energy = linspace(0, 40, 40);

filepath1 = strcat(datpath,'energy_31.txt');
nrg_40GeV_dat_bg_1_5mm = load(filepath1) * 6.2415091E9;
[counts_dat_40GeV_bg_1_5mm, ~] = hist(nrg_40GeV_dat_bg_1_5mm(nrg_40GeV_dat_bg_1_5mm < 40 & nrg_40GeV_dat_bg_1_5mm > 0), energy);
counts_dat_40GeV_bg_norm_1_5mm = energy.*counts_dat_40GeV_bg_1_5mm/2771767;

%% DATA 40GeV 1.5mm e- aligned.
filepath1 = strcat(datpath,'energy_30.txt');
filepath2 = strcat(datpath,'energy_35.txt');
filepath3 = strcat(datpath,'energy_36.txt');
filepath4 = strcat(datpath,'energy_37.txt');
filepath5 = strcat(datpath,'energy_38.txt');

nrg_40GeV_dat1 = load(filepath1) * 6.2415091E9;
nrg_40GeV_dat2 = load(filepath2) * 6.2415091E9;
nrg_40GeV_dat3 = load(filepath3) * 6.2415091E9;
nrg_40GeV_dat4 = load(filepath4) * 6.2415091E9;
nrg_40GeV_dat5 = load(filepath5) * 6.2415091E9;

NEvents_40GeV_aligned_tot_1_5mm = 172307 + 435890 + 538900 + 363630 + 209144;
nrg_40GeV_dat_aligned_tot_1_5mm = [nrg_40GeV_dat1; nrg_40GeV_dat2; nrg_40GeV_dat3; nrg_40GeV_dat4; nrg_40GeV_dat5];
[counts_dat_40GeV_aligned_tot_1_5mm, ~] = hist(nrg_40GeV_dat_aligned_tot_1_5mm(nrg_40GeV_dat_aligned_tot_1_5mm < 40 & nrg_40GeV_dat_aligned_tot_1_5mm > 0), energy);

[counts_dat_40GeV_aligned_1_1_5mm, ~] = hist(nrg_40GeV_dat1(nrg_40GeV_dat1 < 40 & nrg_40GeV_dat1 > 0), energy);
[counts_dat_40GeV_aligned_2_1_5mm, ~] = hist(nrg_40GeV_dat2(nrg_40GeV_dat2 < 40 & nrg_40GeV_dat2 > 0), energy);
[counts_dat_40GeV_aligned_3_1_5mm, ~] = hist(nrg_40GeV_dat3(nrg_40GeV_dat3 < 40 & nrg_40GeV_dat3 > 0), energy);
[counts_dat_40GeV_aligned_4_1_5mm, ~] = hist(nrg_40GeV_dat4(nrg_40GeV_dat4 < 40 & nrg_40GeV_dat4 > 0), energy);
[counts_dat_40GeV_aligned_5_1_5mm, ~] = hist(nrg_40GeV_dat5(nrg_40GeV_dat5 < 40 & nrg_40GeV_dat5 > 0), energy);

counts_dat_40GeV_aligned_norm_tot_1_5mm = energy.*counts_dat_40GeV_aligned_tot_1_5mm/NEvents_40GeV_aligned_tot_1_5mm;

%% SIM 40GeV 1mm e- amorph + backogrund
filepath1 = strcat(datpath,'energy_sim_amorphous1_40GeV_16012019.txt');
filepath2 = strcat(datpath,'energy_sim_amorphous2_40GeV_16012019.txt');
filepath3 = strcat(datpath,'energy_sim_amorphous3_40GeV_16012019.txt');
filepath4 = strcat(datpath,'energy_sim_amorphous4_40GeV_16012019.txt');
filepath5 = strcat(datpath,'energy_sim_amorphous5_40GeV_16012019.txt');

nrg_40GeV_sim1 = load(filepath1) * 6.2415091E9;
nrg_40GeV_sim2 = load(filepath2) * 6.2415091E9;
nrg_40GeV_sim3 = load(filepath3) * 6.2415091E9;
nrg_40GeV_sim4 = load(filepath4) * 6.2415091E9;
nrg_40GeV_sim5 = load(filepath5) * 6.2415091E9;

nrg_40GeV_sim_amorph_bg = [nrg_40GeV_sim1; nrg_40GeV_sim2; nrg_40GeV_sim3; nrg_40GeV_sim4; nrg_40GeV_sim5];
[counts_sim_40GeV_amorph_bg, ~] = hist(nrg_40GeV_sim_amorph_bg(nrg_40GeV_sim_amorph_bg < 40 & nrg_40GeV_sim_amorph_bg > 0), energy);
counts_sim_40GeV_amorph_bg_norm = energy .* counts_sim_40GeV_amorph_bg/5e6;

%% SIM 40GeV 1mm e- background
filepath1 = strcat(datpath,'energy_sim_background1_40GeV_16012019.txt');
filepath2 = strcat(datpath,'energy_sim_background2_40GeV_16012019.txt');
filepath3 = strcat(datpath,'energy_sim_background3_40GeV_16012019.txt');
filepath4 = strcat(datpath,'energy_sim_background4_40GeV_16012019.txt');
filepath5 = strcat(datpath,'energy_sim_background5_40GeV_16012019.txt');

nrg_40GeV_sim1 = load(filepath1) * 6.2415091E9;
nrg_40GeV_sim2 = load(filepath2) * 6.2415091E9;
nrg_40GeV_sim3 = load(filepath3) * 6.2415091E9;
nrg_40GeV_sim4 = load(filepath4) * 6.2415091E9;
nrg_40GeV_sim5 = load(filepath5) * 6.2415091E9;

nrg_40GeV_sim_bg = [nrg_40GeV_sim1; nrg_40GeV_sim2; nrg_40GeV_sim3; nrg_40GeV_sim4; nrg_40GeV_sim5];
[counts_sim_40GeV_bg, ~] = hist(nrg_40GeV_sim_bg(nrg_40GeV_sim_bg < 40 & nrg_40GeV_sim_bg > 0), energy);
counts_sim_40GeV_bg_norm = energy .* counts_sim_40GeV_bg/5e6;

%% SIM 40GeV 1mm e- amorph
filepath1 = strcat(datpath,'energy_sim_amorphous1_40GeV_no_background.txt');
filepath2 = strcat(datpath,'energy_sim_amorphous2_40GeV_no_background.txt');
filepath3 = strcat(datpath,'energy_sim_amorphous3_40GeV_no_background.txt');
filepath4 = strcat(datpath,'energy_sim_amorphous4_40GeV_no_background.txt');
filepath5 = strcat(datpath,'energy_sim_amorphous5_40GeV_no_background.txt');

nrg_40GeV_sim1 = load(filepath1) * 6.2415091E9;
nrg_40GeV_sim2 = load(filepath2) * 6.2415091E9;
nrg_40GeV_sim3 = load(filepath3) * 6.2415091E9;
nrg_40GeV_sim4 = load(filepath4) * 6.2415091E9;
nrg_40GeV_sim5 = load(filepath5) * 6.2415091E9;

nrg_40GeV_sim_amorphous = [nrg_40GeV_sim1; nrg_40GeV_sim2; nrg_40GeV_sim3; nrg_40GeV_sim4; nrg_40GeV_sim5];
[counts_sim_40GeV_amorph, ~] = hist(nrg_40GeV_sim_amorphous(nrg_40GeV_sim_amorphous < 40 & nrg_40GeV_sim_amorphous > 0), energy);
counts_sim_40GeV_amorph_norm = energy .* counts_sim_40GeV_amorph/5e6;

%% SIM 40GeV 1mm e- aligned w/o schot
filepath1 = strcat(datpath,'energy_sim_aligned1_40GeV_woshot.txt');
filepath2 = strcat(datpath,'energy_sim_aligned2_40GeV_woshot.txt');
filepath3 = strcat(datpath,'energy_sim_aligned3_40GeV_woshot.txt');
filepath4 = strcat(datpath,'energy_sim_aligned4_40GeV_woshot.txt');
filepath5 = strcat(datpath,'energy_sim_aligned5_40GeV_woshot.txt');

nrg_40GeV_sim1 = load(filepath1) * 6.2415091E9;
nrg_40GeV_sim2 = load(filepath2) * 6.2415091E9;
nrg_40GeV_sim3 = load(filepath3) * 6.2415091E9;
nrg_40GeV_sim4 = load(filepath4) * 6.2415091E9;
nrg_40GeV_sim5 = load(filepath5) * 6.2415091E9;

nrg_40GeV_sim_aligned = [nrg_40GeV_sim1; nrg_40GeV_sim2; nrg_40GeV_sim3; nrg_40GeV_sim4; nrg_40GeV_sim5];
[counts_sim_40GeV_aligned_woshot, ~] = hist(nrg_40GeV_sim_aligned(nrg_40GeV_sim_aligned < 40 & nrg_40GeV_sim_aligned > 0), energy);
counts_sim_40GeV_aligned_norm_woshot = energy.*counts_sim_40GeV_aligned_woshot/5e6;

%% SIM 40GeV 1mm e- aligned w/o RR
filepath1 = strcat(datpath,'energy_sim_aligned1_40GeV_worr.txt');
filepath2 = strcat(datpath,'energy_sim_aligned2_40GeV_worr.txt');
filepath3 = strcat(datpath,'energy_sim_aligned3_40GeV_worr.txt');
filepath4 = strcat(datpath,'energy_sim_aligned4_40GeV_worr.txt');
filepath5 = strcat(datpath,'energy_sim_aligned5_40GeV_worr.txt');

nrg_40GeV_sim1 = load(filepath1) * 6.2415091E9;
nrg_40GeV_sim2 = load(filepath2) * 6.2415091E9;
nrg_40GeV_sim3 = load(filepath3) * 6.2415091E9;
nrg_40GeV_sim4 = load(filepath4) * 6.2415091E9;
nrg_40GeV_sim5 = load(filepath5) * 6.2415091E9;

nrg_40GeV_sim_aligned = [nrg_40GeV_sim1; nrg_40GeV_sim2; nrg_40GeV_sim3; nrg_40GeV_sim4; nrg_40GeV_sim5];
[counts_sim_40GeV_aligned_worr, ~] = hist(nrg_40GeV_sim_aligned(nrg_40GeV_sim_aligned < 40 & nrg_40GeV_sim_aligned > 0), energy);
counts_sim_40GeV_aligned_norm_worr = energy.*counts_sim_40GeV_aligned_worr/5e6;

%% SIM 40GeV 1mm e- aligned full LL
filepath1 = strcat(datpath,'energy_sim_aligned1_40GeV.txt');
filepath2 = strcat(datpath,'energy_sim_aligned2_40GeV.txt');
filepath3 = strcat(datpath,'energy_sim_aligned3_40GeV.txt');
filepath4 = strcat(datpath,'energy_sim_aligned4_40GeV.txt');
filepath5 = strcat(datpath,'energy_sim_aligned5_40GeV.txt');

nrg_40GeV_sim1 = load(filepath1) * 6.2415091E9;
nrg_40GeV_sim2 = load(filepath2) * 6.2415091E9;
nrg_40GeV_sim3 = load(filepath3) * 6.2415091E9;
nrg_40GeV_sim4 = load(filepath4) * 6.2415091E9;
nrg_40GeV_sim5 = load(filepath5) * 6.2415091E9;
size(nrg_40GeV_sim1)
[counts_sim_40GeV_aligned_1, ~] = hist(nrg_40GeV_sim1(nrg_40GeV_sim1 < 40 & nrg_40GeV_sim1 > 0), energy);
counts_sim_40GeV_aligned_norm_1 = energy.*counts_sim_40GeV_aligned_1/1e6;


nrg_40GeV_sim_aligned = [nrg_40GeV_sim1; nrg_40GeV_sim2; nrg_40GeV_sim3; nrg_40GeV_sim4; nrg_40GeV_sim5];
[counts_sim_40GeV_aligned, ~] = hist(nrg_40GeV_sim_aligned(nrg_40GeV_sim_aligned < 40 & nrg_40GeV_sim_aligned > 0), energy);
counts_sim_40GeV_aligned_norm = energy.*counts_sim_40GeV_aligned/5e6;

%% SIM 40GeV 1.5mm e- amorph
filepath1 = strcat(datpath,'energy_sim_amorphous1_40GeV_1.5mm_no_background.txt');
filepath2 = strcat(datpath,'energy_sim_amorphous2_40GeV_1.5mm_no_background.txt');
filepath3 = strcat(datpath,'energy_sim_amorphous3_40GeV_1.5mm_no_background.txt');
filepath4 = strcat(datpath,'energy_sim_amorphous4_40GeV_1.5mm_no_background.txt');
filepath5 = strcat(datpath,'energy_sim_amorphous5_40GeV_1.5mm_no_background.txt');

nrg_40GeV_sim1 = load(filepath1) * 6.2415091E9;
nrg_40GeV_sim2 = load(filepath2) * 6.2415091E9;
nrg_40GeV_sim3 = load(filepath3) * 6.2415091E9;
nrg_40GeV_sim4 = load(filepath4) * 6.2415091E9;
nrg_40GeV_sim5 = load(filepath5) * 6.2415091E9;

nrg_40GeV_sim_amorph = [nrg_40GeV_sim1; nrg_40GeV_sim2; nrg_40GeV_sim3; nrg_40GeV_sim4; nrg_40GeV_sim5];
[counts_sim_40GeV_amorph_1_5mm, ~] = hist(nrg_40GeV_sim_amorph(nrg_40GeV_sim_amorph < 40 & nrg_40GeV_sim_amorph > 0), energy);
counts_sim_40GeV_amorph_1_5mm_norm = energy.*counts_sim_40GeV_amorph_1_5mm/5e6;

%% SIM 40GeV 1.5mm e- background
filepath1 = strcat(datpath,'energy_sim_background1_40GeV_1.5mm_16012019.txt');
filepath2 = strcat(datpath,'energy_sim_background2_40GeV_1.5mm_16012019.txt');
filepath3 = strcat(datpath,'energy_sim_background3_40GeV_1.5mm_16012019.txt');
filepath4 = strcat(datpath,'energy_sim_background4_40GeV_1.5mm_16012019.txt');
filepath5 = strcat(datpath,'energy_sim_background5_40GeV_1.5mm_16012019.txt');

nrg_40GeV_sim1 = load(filepath1) * 6.2415091E9;
nrg_40GeV_sim2 = load(filepath2) * 6.2415091E9;
nrg_40GeV_sim3 = load(filepath3) * 6.2415091E9;
nrg_40GeV_sim4 = load(filepath4) * 6.2415091E9;
nrg_40GeV_sim5 = load(filepath5) * 6.2415091E9;

nrg_40GeV_sim_bg_1_5mm = [nrg_40GeV_sim1; nrg_40GeV_sim2; nrg_40GeV_sim3; nrg_40GeV_sim4; nrg_40GeV_sim5];
[counts_sim_40GeV_bg_1_5mm, ~] = hist(nrg_40GeV_sim_bg_1_5mm(nrg_40GeV_sim_bg_1_5mm < 40 & nrg_40GeV_sim_bg_1_5mm > 0), energy);
counts_sim_40GeV_bg_norm_1_5mm = energy.*counts_sim_40GeV_bg_1_5mm/5e6;

%% SIM 40GeV 1.5mm e- amorph + backogrund
filepath1 = strcat(datpath,'energy_sim_amorphous1_40GeV_1.5mm_16012019.txt');
filepath2 = strcat(datpath,'energy_sim_amorphous2_40GeV_1.5mm_16012019.txt');
filepath3 = strcat(datpath,'energy_sim_amorphous3_40GeV_1.5mm_16012019.txt');
filepath4 = strcat(datpath,'energy_sim_amorphous4_40GeV_1.5mm_16012019.txt');
filepath5 = strcat(datpath,'energy_sim_amorphous5_40GeV_1.5mm_16012019.txt');

nrg_40GeV_sim1 = load(filepath1) * 6.2415091E9;
nrg_40GeV_sim2 = load(filepath2) * 6.2415091E9;
nrg_40GeV_sim3 = load(filepath3) * 6.2415091E9;
nrg_40GeV_sim4 = load(filepath4) * 6.2415091E9;
nrg_40GeV_sim5 = load(filepath5) * 6.2415091E9;

nrg_40GeV_sim_amorph_bg_1_5mm = [nrg_40GeV_sim1; nrg_40GeV_sim2; nrg_40GeV_sim3; nrg_40GeV_sim4; nrg_40GeV_sim5];
[counts_sim_40GeV_amorph_bg_1_5mm, ~] = hist(nrg_40GeV_sim_amorph_bg_1_5mm(nrg_40GeV_sim_amorph_bg_1_5mm < 40 & nrg_40GeV_sim_amorph_bg_1_5mm > 0), energy);
counts_sim_40GeV_amorph_bg_norm_1_5mm = energy.*counts_sim_40GeV_amorph_bg_1_5mm/5e6;

%% SIM 40GeV 1.5mm e- aligned w/o schot
filepath1 = strcat(datpath,'energy_sim_aligned1_40GeV_woshot_1.5mm.txt');
filepath2 = strcat(datpath,'energy_sim_aligned2_40GeV_woshot_1.5mm.txt');
filepath3 = strcat(datpath,'energy_sim_aligned3_40GeV_woshot_1.5mm.txt');
filepath4 = strcat(datpath,'energy_sim_aligned4_40GeV_woshot_1.5mm.txt');
filepath5 = strcat(datpath,'energy_sim_aligned5_40GeV_woshot_1.5mm.txt');

nrg_40GeV_sim1 = load(filepath1) * 6.2415091E9;
nrg_40GeV_sim2 = load(filepath2) * 6.2415091E9;
nrg_40GeV_sim3 = load(filepath3) * 6.2415091E9;
nrg_40GeV_sim4 = load(filepath4) * 6.2415091E9;
nrg_40GeV_sim5 = load(filepath5) * 6.2415091E9;

nrg_40GeV_sim_aligned = [nrg_40GeV_sim1; nrg_40GeV_sim2; nrg_40GeV_sim3; nrg_40GeV_sim4; nrg_40GeV_sim5];
[counts_sim_40GeV_aligned_1_5mm_woshot, ~] = hist(nrg_40GeV_sim_aligned(nrg_40GeV_sim_aligned < 40 & nrg_40GeV_sim_aligned > 0), energy);
counts_sim_40GeV_aligned_1_5mm_norm_woshot = energy.*counts_sim_40GeV_aligned_1_5mm_woshot/5e6;

%% SIM 40GeV 1.5mm e- aligned w/o RR
filepath1 = strcat(datpath,'energy_sim_aligned1_40GeV_worr_1.5mm.txt');
filepath2 = strcat(datpath,'energy_sim_aligned2_40GeV_worr_1.5mm.txt');
filepath3 = strcat(datpath,'energy_sim_aligned3_40GeV_worr_1.5mm.txt');
filepath4 = strcat(datpath,'energy_sim_aligned4_40GeV_worr_1.5mm.txt');
filepath5 = strcat(datpath,'energy_sim_aligned5_40GeV_worr_1.5mm.txt');

nrg_40GeV_sim1 = load(filepath1) * 6.2415091E9;
nrg_40GeV_sim2 = load(filepath2) * 6.2415091E9;
nrg_40GeV_sim3 = load(filepath3) * 6.2415091E9;
nrg_40GeV_sim4 = load(filepath4) * 6.2415091E9;
nrg_40GeV_sim5 = load(filepath5) * 6.2415091E9;

nrg_40GeV_sim_aligned = [nrg_40GeV_sim1; nrg_40GeV_sim2; nrg_40GeV_sim3; nrg_40GeV_sim4; nrg_40GeV_sim5];
[counts_sim_40GeV_aligned_1_5mm_worr, ~] = hist(nrg_40GeV_sim_aligned(nrg_40GeV_sim_aligned < 40 & nrg_40GeV_sim_aligned > 0), energy);
counts_sim_40GeV_aligned_1_5mm_norm_worr = energy.*counts_sim_40GeV_aligned_1_5mm_worr/5e6;

%% SIM 40GeV 1.5mm e- aligned full LL
filepath1 = strcat(datpath,'energy_sim_aligned1_40GeV_1.5mm.txt');
filepath2 = strcat(datpath,'energy_sim_aligned2_40GeV_1.5mm.txt');
filepath3 = strcat(datpath,'energy_sim_aligned3_40GeV_1.5mm.txt');
filepath4 = strcat(datpath,'energy_sim_aligned4_40GeV_1.5mm.txt');
filepath5 = strcat(datpath,'energy_sim_aligned5_40GeV_1.5mm.txt');

nrg_40GeV_sim1 = load(filepath1) * 6.2415091E9;
nrg_40GeV_sim2 = load(filepath2) * 6.2415091E9;
nrg_40GeV_sim3 = load(filepath3) * 6.2415091E9;
nrg_40GeV_sim4 = load(filepath4) * 6.2415091E9;
nrg_40GeV_sim5 = load(filepath5) * 6.2415091E9;

nrg_40GeV_sim_aligned = [nrg_40GeV_sim1; nrg_40GeV_sim2; nrg_40GeV_sim3; nrg_40GeV_sim4; nrg_40GeV_sim5];
[counts_sim_40GeV_aligned_1_5mm, ~] = hist(nrg_40GeV_sim_aligned(nrg_40GeV_sim_aligned < 40 & nrg_40GeV_sim_aligned > 0), energy);
counts_sim_40GeV_aligned_1_5mm_norm = energy.*counts_sim_40GeV_aligned_1_5mm/5e6;

%% DATA 80GeV 1mm e- background
energy = linspace(0, 80, 40);

filepath1 = strcat(datpath,'energy_85.txt');
filepath2 = strcat(datpath,'energy_90.txt');
filepath3 = strcat(datpath,'energy_91.txt');

NEvents_80GeV_dat_bg_tot = 1911215 + 1892132 + 1219189;
nrg_80GeV_dat_bg1 = load(filepath1) * 6.2415091E9;
nrg_80GeV_dat_bg2 = load(filepath2) * 6.2415091E9;
nrg_80GeV_dat_bg3 = load(filepath3) * 6.2415091E9;

nrg_80GeV_dat_bg_tot = [nrg_80GeV_dat_bg1; nrg_80GeV_dat_bg2; nrg_80GeV_dat_bg3];
[counts_dat_80GeV_bg_tot, ~] = hist(nrg_80GeV_dat_bg_tot(nrg_80GeV_dat_bg_tot < 80 & nrg_80GeV_dat_bg_tot > 0), energy);
counts_dat_80GeV_bg_norm = energy.*counts_dat_80GeV_bg_tot/NEvents_80GeV_dat_bg_tot;

%% DATA 80GeV 1mm e- amoprh.
filepath1 = strcat(datpath,'energy_86.txt');
filepath2 = strcat(datpath,'energy_87.txt');
filepath3 = strcat(datpath,'energy_88.txt');
filepath4 = strcat(datpath,'energy_89.txt');

nrg_80GeV_dat1 = load(filepath1) * 6.2415091E9;
nrg_80GeV_dat2 = load(filepath2) * 6.2415091E9;
nrg_80GeV_dat3 = load(filepath3) * 6.2415091E9;
nrg_80GeV_dat4 = load(filepath4) * 6.2415091E9;

NEvents_80GeV_amorph_tot = 1719461 + 1172521 + 1210722 + 538281;
nrg_80GeV_dat_amorph_tot = [nrg_80GeV_dat1; nrg_80GeV_dat2; nrg_80GeV_dat3; nrg_80GeV_dat4];
[counts_dat_80GeV_amorph_tot, ~] = hist(nrg_80GeV_dat_amorph_tot(nrg_80GeV_dat_amorph_tot < 80 & nrg_80GeV_dat_amorph_tot > 0), energy);

[counts_dat_80GeV_amorph_1, ~] = hist(nrg_80GeV_dat1(nrg_80GeV_dat1 < 80 & nrg_80GeV_dat1 > 0), energy);
[counts_dat_80GeV_amorph_2, ~] = hist(nrg_80GeV_dat2(nrg_80GeV_dat2 < 80 & nrg_80GeV_dat2 > 0), energy);
[counts_dat_80GeV_amorph_3, ~] = hist(nrg_80GeV_dat3(nrg_80GeV_dat3 < 80 & nrg_80GeV_dat3 > 0), energy);
[counts_dat_80GeV_amorph_4, ~] = hist(nrg_80GeV_dat4(nrg_80GeV_dat4 < 80 & nrg_80GeV_dat4 > 0), energy);

counts_dat_80GeV_amorph_norm_tot = energy.*counts_dat_80GeV_amorph_tot/NEvents_80GeV_amorph_tot;

counts_dat_80GeV_amorph_norm_1 = energy.*counts_dat_80GeV_amorph_1/1837156;
counts_dat_80GeV_amorph_norm_2 = energy.*counts_dat_80GeV_amorph_2/1254664;
counts_dat_80GeV_amorph_norm_3 = energy.*counts_dat_80GeV_amorph_3/1295201;
counts_dat_80GeV_amorph_norm_4 = energy.*counts_dat_80GeV_amorph_4/577897;

%% DATA 80GeV 1mm e- Aligned
filepath1 = strcat(datpath,'energy_84.txt');
filepath2 = strcat(datpath,'energy_92.txt');
filepath3 = strcat(datpath,'energy_93.txt');
filepath4 = strcat(datpath,'energy_94.txt');
filepath5 = strcat(datpath,'energy_95.txt');

nrg_80GeV_dat1 = load(filepath1) * 6.2415091E9;
nrg_80GeV_dat2 = load(filepath2) * 6.2415091E9;
nrg_80GeV_dat3 = load(filepath3) * 6.2415091E9;
nrg_80GeV_dat4 = load(filepath4) * 6.2415091E9;
nrg_80GeV_dat5 = load(filepath5) * 6.2415091E9;

NEvents_80GeV_aligned_tot = 1134032 + 462880 + 943921 + 1833415 + 35424;
nrg_80GeV_dat_aligned_tot = [nrg_80GeV_dat1; nrg_80GeV_dat2; nrg_80GeV_dat3; nrg_80GeV_dat4; nrg_80GeV_dat5];
[counts_dat_80GeV_aligned_tot, ~] = hist(nrg_80GeV_dat_aligned_tot(nrg_80GeV_dat_aligned_tot < 80 & nrg_80GeV_dat_aligned_tot > 0), energy);

[counts_dat_80GeV_aligned_1, ~] = hist(nrg_80GeV_dat1(nrg_80GeV_dat1 < 80 & nrg_80GeV_dat1 > 0), energy);
[counts_dat_80GeV_aligned_2, ~] = hist(nrg_80GeV_dat2(nrg_80GeV_dat2 < 80 & nrg_80GeV_dat2 > 0), energy);
[counts_dat_80GeV_aligned_3, ~] = hist(nrg_80GeV_dat3(nrg_80GeV_dat3 < 80 & nrg_80GeV_dat3 > 0), energy);
[counts_dat_80GeV_aligned_4, ~] = hist(nrg_80GeV_dat4(nrg_80GeV_dat4 < 80 & nrg_80GeV_dat4 > 0), energy);
[counts_dat_80GeV_aligned_5, ~] = hist(nrg_80GeV_dat5(nrg_80GeV_dat5 < 80 & nrg_80GeV_dat5 > 0), energy);

counts_dat_80GeV_aligned_norm_tot = energy.*counts_dat_80GeV_aligned_tot/NEvents_80GeV_aligned_tot;
% 

%% SIM 80GeV 1mm e- amorph + backogrund
filepath1 = strcat(datpath,'energy_sim_amorphous1_80GeV_16012019.txt');
filepath2 = strcat(datpath,'energy_sim_amorphous2_80GeV_16012019.txt');
filepath3 = strcat(datpath,'energy_sim_amorphous3_80GeV_16012019.txt');
filepath4 = strcat(datpath,'energy_sim_amorphous4_80GeV_16012019.txt');
filepath5 = strcat(datpath,'energy_sim_amorphous5_80GeV_16012019.txt');

nrg_80GeV_sim1 = load(filepath1) * 6.2415091E9;
nrg_80GeV_sim2 = load(filepath2) * 6.2415091E9;
nrg_80GeV_sim3 = load(filepath3) * 6.2415091E9;
nrg_80GeV_sim4 = load(filepath4) * 6.2415091E9;
nrg_80GeV_sim5 = load(filepath5) * 6.2415091E9;

nrg_80GeV_sim_amorph_bg = [nrg_80GeV_sim1; nrg_80GeV_sim2; nrg_80GeV_sim3; nrg_80GeV_sim4; nrg_80GeV_sim5];
[counts_sim_80GeV_amorph_bg, ~] = hist(nrg_80GeV_sim_amorph_bg(nrg_80GeV_sim_amorph_bg < 80 & nrg_80GeV_sim_amorph_bg > 0), energy);
counts_sim_80GeV_amorph_bg_norm = energy.*counts_sim_80GeV_amorph_bg/5e6;

%% SIM 80GeV 1mm e- amorph
filepath1 = strcat(datpath,'energy_sim_amorphous1_80GeV_no_background_04012019.txt');
filepath2 = strcat(datpath,'energy_sim_amorphous2_80GeV_no_background_04012019.txt');
filepath3 = strcat(datpath,'energy_sim_amorphous3_80GeV_no_background_04012019.txt');
filepath4 = strcat(datpath,'energy_sim_amorphous4_80GeV_no_background_04012019.txt');
filepath5 = strcat(datpath,'energy_sim_amorphous5_80GeV_no_background_04012019.txt');

nrg_80GeV_sim1 = load(filepath1) * 6.2415091E9;
nrg_80GeV_sim2 = load(filepath2) * 6.2415091E9;
nrg_80GeV_sim3 = load(filepath3) * 6.2415091E9;
nrg_80GeV_sim4 = load(filepath4) * 6.2415091E9;
nrg_80GeV_sim5 = load(filepath5) * 6.2415091E9;

nrg_80GeV_sim_amorph = [nrg_80GeV_sim1; nrg_80GeV_sim2; nrg_80GeV_sim3; nrg_80GeV_sim4; nrg_80GeV_sim5];
[counts_sim_80GeV_amorph, ~] = hist(nrg_80GeV_sim_amorph(nrg_80GeV_sim_amorph < 80 & nrg_80GeV_sim_amorph > 0), energy);
counts_sim_80GeV_amorph_norm = energy.*counts_sim_80GeV_amorph/5e6;

%% SIM 80GeV 1mm e- background
filepath1 = strcat(datpath,'energy_sim_background1_80GeV_16012019.txt');
filepath2 = strcat(datpath,'energy_sim_background2_80GeV_16012019.txt');
filepath3 = strcat(datpath,'energy_sim_background3_80GeV_16012019.txt');
filepath4 = strcat(datpath,'energy_sim_background4_80GeV_16012019.txt');
filepath5 = strcat(datpath,'energy_sim_background5_80GeV_16012019.txt');

nrg_80GeV_sim1 = load(filepath1) * 6.2415091E9;
nrg_80GeV_sim2 = load(filepath2) * 6.2415091E9;
nrg_80GeV_sim3 = load(filepath3) * 6.2415091E9;
nrg_80GeV_sim4 = load(filepath4) * 6.2415091E9;
nrg_80GeV_sim5 = load(filepath5) * 6.2415091E9;

nrg_80GeV_sim_bg = [nrg_80GeV_sim1; nrg_80GeV_sim2; nrg_80GeV_sim3; nrg_80GeV_sim4; nrg_80GeV_sim5];
[counts_sim_80GeV_bg, ~] = hist(nrg_80GeV_sim_bg(nrg_80GeV_sim_bg < 80 & nrg_80GeV_sim_bg > 0), energy);
counts_sim_80GeV_bg_norm = energy.*counts_sim_80GeV_bg/5e6;

%% SIM 80GeV 1mm e- aligned w/o schot
filepath1 = strcat(datpath,'energy_sim_aligned1_80GeV_woshot.txt');
filepath2 = strcat(datpath,'energy_sim_aligned2_80GeV_woshot.txt');
filepath3 = strcat(datpath,'energy_sim_aligned3_80GeV_woshot.txt');
filepath4 = strcat(datpath,'energy_sim_aligned4_80GeV_woshot.txt');
filepath5 = strcat(datpath,'energy_sim_aligned5_80GeV_woshot.txt');

nrg_80GeV_sim1 = load(filepath1) * 6.2415091E9;
nrg_80GeV_sim2 = load(filepath2) * 6.2415091E9;
nrg_80GeV_sim3 = load(filepath3) * 6.2415091E9;
nrg_80GeV_sim4 = load(filepath4) * 6.2415091E9;
nrg_80GeV_sim5 = load(filepath5) * 6.2415091E9;

nrg_80GeV_sim_aligned = [nrg_80GeV_sim1; nrg_80GeV_sim2; nrg_80GeV_sim3; nrg_80GeV_sim4; nrg_80GeV_sim5];
[counts_sim_80GeV_aligned_woshot, ~] = hist(nrg_80GeV_sim_aligned(nrg_80GeV_sim_aligned < 80 & nrg_80GeV_sim_aligned > 0), energy);
counts_sim_80GeV_aligned_norm_woshot = energy.*counts_sim_80GeV_aligned_woshot/5e6;

%% SIM 80GeV 1mm e- aligned w/o RR
filepath1 = strcat(datpath,'energy_sim_aligned1_80GeV_worr.txt');
filepath2 = strcat(datpath,'energy_sim_aligned2_80GeV_worr.txt');
filepath3 = strcat(datpath,'energy_sim_aligned3_80GeV_worr.txt');
filepath4 = strcat(datpath,'energy_sim_aligned4_80GeV_worr.txt');
filepath5 = strcat(datpath,'energy_sim_aligned5_80GeV_worr.txt');

nrg_80GeV_sim1 = load(filepath1) * 6.2415091E9;
nrg_80GeV_sim2 = load(filepath2) * 6.2415091E9;
nrg_80GeV_sim3 = load(filepath3) * 6.2415091E9;
nrg_80GeV_sim4 = load(filepath4) * 6.2415091E9;
nrg_80GeV_sim5 = load(filepath5) * 6.2415091E9;

nrg_80GeV_sim_aligned = [nrg_80GeV_sim1; nrg_80GeV_sim2; nrg_80GeV_sim3; nrg_80GeV_sim4; nrg_80GeV_sim5];
[counts_sim_80GeV_aligned_worr, ~] = hist(nrg_80GeV_sim_aligned(nrg_80GeV_sim_aligned < 80 & nrg_80GeV_sim_aligned > 0), energy);
counts_sim_80GeV_aligned_norm_worr = energy.*counts_sim_80GeV_aligned_worr/5e6;

%% SIM 80GeV 1mm e- aligned full LL
filepath1 = strcat(datpath,'energy_sim_aligned1_80GeV.txt');
filepath2 = strcat(datpath,'energy_sim_aligned2_80GeV.txt');
filepath3 = strcat(datpath,'energy_sim_aligned3_80GeV.txt');
filepath4 = strcat(datpath,'energy_sim_aligned4_80GeV.txt');
filepath5 = strcat(datpath,'energy_sim_aligned5_80GeV.txt');

nrg_80GeV_sim1 = load(filepath1) * 6.2415091E9;
nrg_80GeV_sim2 = load(filepath2) * 6.2415091E9;
nrg_80GeV_sim3 = load(filepath3) * 6.2415091E9;
nrg_80GeV_sim4 = load(filepath4) * 6.2415091E9;
nrg_80GeV_sim5 = load(filepath5) * 6.2415091E9;

nrg_80GeV_sim_aligned = [nrg_80GeV_sim1; nrg_80GeV_sim2; nrg_80GeV_sim3; nrg_80GeV_sim4; nrg_80GeV_sim5];
[counts_sim_80GeV_aligned, ~] = hist(nrg_80GeV_sim_aligned(nrg_80GeV_sim_aligned < 80 & nrg_80GeV_sim_aligned > 0), energy);
counts_sim_80GeV_aligned_norm = energy.*counts_sim_80GeV_aligned/5e6;

%% DATA 80GeV 1.5mm e- background
filepath1 = strcat(datpath,'energy_48.txt');
filepath2 = strcat(datpath,'energy_54.txt');
filepath3 = strcat(datpath,'energy_55.txt');

NEvents_80GeV_dat_bg_1_5mm_tot = 2497897 + 312302 + 216860;

nrg_80GeV_dat_bg1 = load(filepath1) * 6.2415091E9;
nrg_80GeV_dat_bg2 = load(filepath2) * 6.2415091E9;
nrg_80GeV_dat_bg3 = load(filepath3) * 6.2415091E9;

nrg_80GeV_dat_bg_1_5mm_tot = [nrg_80GeV_dat_bg1; nrg_80GeV_dat_bg2; nrg_80GeV_dat_bg3];
[counts_dat_80GeV_bg_1_5mm_tot, ~] = hist(nrg_80GeV_dat_bg_1_5mm_tot(nrg_80GeV_dat_bg_1_5mm_tot < 80 & nrg_80GeV_dat_bg_1_5mm_tot > 0), energy);
counts_dat_80GeV_bg_norm_1_5mm_tot = energy.*counts_dat_80GeV_bg_1_5mm_tot/NEvents_80GeV_dat_bg_1_5mm_tot;

%% DATA 80GeV 1.5mm e- amoprh.
filepath1 = strcat(datpath,'energy_49.txt');
filepath2 = strcat(datpath,'energy_50.txt');
filepath3 = strcat(datpath,'energy_51.txt');
filepath4 = strcat(datpath,'energy_52.txt');
filepath5 = strcat(datpath,'energy_53.txt');
filepath6 = strcat(datpath,'energy_56.txt');
filepath7 = strcat(datpath,'energy_57.txt');

nrg_80GeV_dat1 = load(filepath1) * 6.2415091E9;
nrg_80GeV_dat2 = load(filepath2) * 6.2415091E9;
nrg_80GeV_dat3 = load(filepath3) * 6.2415091E9;
nrg_80GeV_dat4 = load(filepath4) * 6.2415091E9;
nrg_80GeV_dat5 = load(filepath5) * 6.2415091E9;
nrg_80GeV_dat6 = load(filepath6) * 6.2415091E9;
nrg_80GeV_dat7 = load(filepath7) * 6.2415091E9;

NEvents_80GeV_amorph_1_5mm_tot = 873246 + 434680 + 847524 + 182889 + 18846 + 392613 + 495068;
nrg_80GeV_dat_amorph_1_5mm_tot = [nrg_80GeV_dat1; nrg_80GeV_dat2; nrg_80GeV_dat3; nrg_80GeV_dat4; nrg_80GeV_dat5; nrg_80GeV_dat6; nrg_80GeV_dat7];
[counts_dat_80GeV_amorph_1_5mm_tot, ~] = hist(nrg_80GeV_dat_amorph_1_5mm_tot(nrg_80GeV_dat_amorph_1_5mm_tot < 80 & nrg_80GeV_dat_amorph_1_5mm_tot > 0), energy);

[counts_dat_80GeV_amorph_1, ~] = hist(nrg_80GeV_dat1(nrg_80GeV_dat1 < 80 & nrg_80GeV_dat1 > 0), energy);
[counts_dat_80GeV_amorph_2, ~] = hist(nrg_80GeV_dat2(nrg_80GeV_dat2 < 80 & nrg_80GeV_dat2 > 0), energy);
[counts_dat_80GeV_amorph_3, ~] = hist(nrg_80GeV_dat3(nrg_80GeV_dat3 < 80 & nrg_80GeV_dat3 > 0), energy);
[counts_dat_80GeV_amorph_4, ~] = hist(nrg_80GeV_dat4(nrg_80GeV_dat4 < 80 & nrg_80GeV_dat4 > 0), energy);
[counts_dat_80GeV_amorph_5, ~] = hist(nrg_80GeV_dat5(nrg_80GeV_dat5 < 80 & nrg_80GeV_dat5 > 0), energy);
[counts_dat_80GeV_amorph_6, ~] = hist(nrg_80GeV_dat6(nrg_80GeV_dat6 < 80 & nrg_80GeV_dat6 > 0), energy);
[counts_dat_80GeV_amorph_7, ~] = hist(nrg_80GeV_dat7(nrg_80GeV_dat7 < 80 & nrg_80GeV_dat7 > 0), energy);

counts_dat_80GeV_amorph_norm_1_5mm_tot = energy.*counts_dat_80GeV_amorph_1_5mm_tot/NEvents_80GeV_amorph_1_5mm_tot;

counts_dat_80GeV_amorph_norm_1_5mm_1 = energy.*counts_dat_80GeV_amorph_1/928547;
counts_dat_80GeV_amorph_norm_1_5mm_2 = energy.*counts_dat_80GeV_amorph_2/460187;
counts_dat_80GeV_amorph_norm_1_5mm_3 = energy.*counts_dat_80GeV_amorph_3/894524;
counts_dat_80GeV_amorph_norm_1_5mm_4 = energy.*counts_dat_80GeV_amorph_4/192621;
counts_dat_80GeV_amorph_norm_1_5mm_5 = energy.*counts_dat_80GeV_amorph_5/19870;
counts_dat_80GeV_amorph_norm_1_5mm_6 = energy.*counts_dat_80GeV_amorph_6/418065;
counts_dat_80GeV_amorph_norm_1_5mm_7 = energy.*counts_dat_80GeV_amorph_7/526675;

%% DATA 80GeV 1.5mm e- aligned
filepath1 = strcat(datpath,'energy_46.txt');
filepath2 = strcat(datpath,'energy_47.txt');

nrg_80GeV_dat1 = load(filepath1) * 6.2415091E9;
nrg_80GeV_dat2 = load(filepath2) * 6.2415091E9;

NEvents_80GeV_aligned_1_5mm_tot = 431667 + 942149;
nrg_80GeV_dat_align_1_5mm_tot = [nrg_80GeV_dat1; nrg_80GeV_dat2];
[counts_dat_80GeV_aligned_1_5mm_tot, ~] = hist(nrg_80GeV_dat_align_1_5mm_tot(nrg_80GeV_dat_align_1_5mm_tot < 80 & nrg_80GeV_dat_align_1_5mm_tot > 0), energy);

counts_dat_80GeV_aligned_norm_1_5mm_tot = energy.*counts_dat_80GeV_aligned_1_5mm_tot/NEvents_80GeV_aligned_1_5mm_tot;

%% SIM 80GeV 1.5mm e- amorph + backogrund
filepath1 = strcat(datpath,'energy_sim_amorphous1_80GeV_1.5mm_16012019.txt');
filepath2 = strcat(datpath,'energy_sim_amorphous2_80GeV_1.5mm_16012019.txt');
filepath3 = strcat(datpath,'energy_sim_amorphous3_80GeV_1.5mm_16012019.txt');
filepath4 = strcat(datpath,'energy_sim_amorphous4_80GeV_1.5mm_16012019.txt');
filepath5 = strcat(datpath,'energy_sim_amorphous5_80GeV_1.5mm_16012019.txt');

nrg_80GeV_sim1 = load(filepath1) * 6.2415091E9;
nrg_80GeV_sim2 = load(filepath2) * 6.2415091E9;
nrg_80GeV_sim3 = load(filepath3) * 6.2415091E9;
nrg_80GeV_sim4 = load(filepath4) * 6.2415091E9;
nrg_80GeV_sim5 = load(filepath5) * 6.2415091E9;

nrg_80GeV_sim_amorph = [nrg_80GeV_sim1; nrg_80GeV_sim2; nrg_80GeV_sim3; nrg_80GeV_sim4; nrg_80GeV_sim5];
[counts_sim_80GeV_amorph_bg_1_5mm, ~] = hist(nrg_80GeV_sim_amorph(nrg_80GeV_sim_amorph < 80 & nrg_80GeV_sim_amorph > 0), energy);
counts_sim_80GeV_amorph_bg_1_5mm_norm = energy.*counts_sim_80GeV_amorph_bg_1_5mm/5e6;

%% SIM 80GeV 1.5mm e- amorph
filepath1 = strcat(datpath,'energy_sim_amorphous1_80GeV_1.5mm_no_background.txt');
filepath2 = strcat(datpath,'energy_sim_amorphous2_80GeV_1.5mm_no_background.txt');
filepath3 = strcat(datpath,'energy_sim_amorphous3_80GeV_1.5mm_no_background.txt');
filepath4 = strcat(datpath,'energy_sim_amorphous4_80GeV_1.5mm_no_background.txt');
filepath5 = strcat(datpath,'energy_sim_amorphous5_80GeV_1.5mm_no_background.txt');

nrg_80GeV_sim1 = load(filepath1) * 6.2415091E9;
nrg_80GeV_sim2 = load(filepath2) * 6.2415091E9;
nrg_80GeV_sim3 = load(filepath3) * 6.2415091E9;
nrg_80GeV_sim4 = load(filepath4) * 6.2415091E9;
nrg_80GeV_sim5 = load(filepath5) * 6.2415091E9;

nrg_80GeV_sim_amorph = [nrg_80GeV_sim1; nrg_80GeV_sim2; nrg_80GeV_sim3; nrg_80GeV_sim4; nrg_80GeV_sim5];
[counts_sim_80GeV_amorph_1_5mm, ~] = hist(nrg_80GeV_sim_amorph(nrg_80GeV_sim_amorph < 80 & nrg_80GeV_sim_amorph > 0), energy);
counts_sim_80GeV_amorph_1_5mm_norm = energy.*counts_sim_80GeV_amorph_1_5mm/5e6;

%% SIM 80GeV 1.5mm e- background
filepath1 = strcat(datpath,'energy_sim_background1_80GeV_1.5mm_16012019.txt');
filepath2 = strcat(datpath,'energy_sim_background2_80GeV_1.5mm_16012019.txt');
filepath3 = strcat(datpath,'energy_sim_background3_80GeV_1.5mm_16012019.txt');
filepath4 = strcat(datpath,'energy_sim_background4_80GeV_1.5mm_16012019.txt');
filepath5 = strcat(datpath,'energy_sim_background5_80GeV_1.5mm_16012019.txt');

nrg_80GeV_sim1 = load(filepath1) * 6.2415091E9;
nrg_80GeV_sim2 = load(filepath2) * 6.2415091E9;
nrg_80GeV_sim3 = load(filepath3) * 6.2415091E9;
nrg_80GeV_sim4 = load(filepath4) * 6.2415091E9;
nrg_80GeV_sim5 = load(filepath5) * 6.2415091E9;

nrg_80GeV_sim_bg = [nrg_80GeV_sim1; nrg_80GeV_sim2; nrg_80GeV_sim3; nrg_80GeV_sim4; nrg_80GeV_sim5];
[counts_sim_80GeV_bg_1_5mm, ~] = hist(nrg_80GeV_sim_bg(nrg_80GeV_sim_bg < 80 & nrg_80GeV_sim_bg > 0), energy);
counts_sim_80GeV_bg_1_5mm_norm = energy.*counts_sim_80GeV_bg_1_5mm/5e6;

%% SIM 80GeV 1.5mm e- aligned w/o schot
filepath1 = strcat(datpath,'energy_sim_aligned1_80GeV_woshot_1.5mm.txt');
filepath2 = strcat(datpath,'energy_sim_aligned2_80GeV_woshot_1.5mm.txt');
filepath3 = strcat(datpath,'energy_sim_aligned3_80GeV_woshot_1.5mm.txt');
filepath4 = strcat(datpath,'energy_sim_aligned4_80GeV_woshot_1.5mm.txt');
filepath5 = strcat(datpath,'energy_sim_aligned5_80GeV_woshot_1.5mm.txt');

nrg_80GeV_sim1 = load(filepath1) * 6.2415091E9;
nrg_80GeV_sim2 = load(filepath2) * 6.2415091E9;
nrg_80GeV_sim3 = load(filepath3) * 6.2415091E9;
nrg_80GeV_sim4 = load(filepath4) * 6.2415091E9;
nrg_80GeV_sim5 = load(filepath5) * 6.2415091E9;

nrg_80GeV_sim_aligned = [nrg_80GeV_sim1; nrg_80GeV_sim2; nrg_80GeV_sim3; nrg_80GeV_sim4; nrg_80GeV_sim5];
[counts_sim_80GeV_aligned_1_5mm_woshot, ~] = hist(nrg_80GeV_sim_aligned(nrg_80GeV_sim_aligned < 80 & nrg_80GeV_sim_aligned > 0), energy);
counts_sim_80GeV_aligned_1_5mm_norm_woshot = energy.*counts_sim_80GeV_aligned_1_5mm_woshot/5e6;

%% SIM 80GeV 1.5mm e- aligned w/o RR
filepath1 = strcat(datpath,'energy_sim_aligned1_80GeV_worr_1.5mm.txt');
filepath2 = strcat(datpath,'energy_sim_aligned2_80GeV_worr_1.5mm.txt');
filepath3 = strcat(datpath,'energy_sim_aligned3_80GeV_worr_1.5mm.txt');
filepath4 = strcat(datpath,'energy_sim_aligned4_80GeV_worr_1.5mm.txt');
filepath5 = strcat(datpath,'energy_sim_aligned5_80GeV_worr_1.5mm.txt');

nrg_80GeV_sim1 = load(filepath1) * 6.2415091E9;
nrg_80GeV_sim2 = load(filepath2) * 6.2415091E9;
nrg_80GeV_sim3 = load(filepath3) * 6.2415091E9;
nrg_80GeV_sim4 = load(filepath4) * 6.2415091E9;
nrg_80GeV_sim5 = load(filepath5) * 6.2415091E9;

nrg_80GeV_sim_aligned = [nrg_80GeV_sim1; nrg_80GeV_sim2; nrg_80GeV_sim3; nrg_80GeV_sim4; nrg_80GeV_sim5];
[counts_sim_80GeV_aligned_1_5mm_worr, ~] = hist(nrg_80GeV_sim_aligned(nrg_80GeV_sim_aligned < 80 & nrg_80GeV_sim_aligned > 0), energy);
counts_sim_80GeV_aligned_1_5mm_norm_worr = energy.*counts_sim_80GeV_aligned_1_5mm_worr/5e6;

%% SIM 80GeV 1.5mm e- aligned full LL
filepath1 = strcat(datpath,'energy_sim_aligned1_80GeV_1.5mm.txt');
filepath2 = strcat(datpath,'energy_sim_aligned2_80GeV_1.5mm.txt');
filepath3 = strcat(datpath,'energy_sim_aligned3_80GeV_1.5mm.txt');
filepath4 = strcat(datpath,'energy_sim_aligned4_80GeV_1.5mm.txt');
filepath5 = strcat(datpath,'energy_sim_aligned5_80GeV_1.5mm.txt');

nrg_80GeV_sim1 = load(filepath1) * 6.2415091E9;
nrg_80GeV_sim2 = load(filepath2) * 6.2415091E9;
nrg_80GeV_sim3 = load(filepath3) * 6.2415091E9;
nrg_80GeV_sim4 = load(filepath4) * 6.2415091E9;
nrg_80GeV_sim5 = load(filepath5) * 6.2415091E9;

nrg_80GeV_sim_aligned = [nrg_80GeV_sim1; nrg_80GeV_sim2; nrg_80GeV_sim3; nrg_80GeV_sim4; nrg_80GeV_sim5];
[counts_sim_80GeV_aligned_1_5mm, ~] = hist(nrg_80GeV_sim_aligned(nrg_80GeV_sim_aligned < 80 & nrg_80GeV_sim_aligned > 0), energy);
counts_sim_80GeV_aligned_1_5mm_norm = energy.*counts_sim_80GeV_aligned_1_5mm/5e6;
%  

%% Calibration factors
delta_40 = @(k) norm(counts_dat_40GeV_amorph_norm_tot - counts_dat_40GeV_bg_norm - k .* ( (counts_sim_40GeV_amorph_bg_norm - counts_sim_40GeV_bg_norm)));
delta_40_1_5mm = @(k) norm(counts_dat_40GeV_amorph_norm_tot_1_5mm - counts_dat_40GeV_bg_norm_1_5mm - k .* ( (counts_sim_40GeV_amorph_bg_norm_1_5mm - counts_sim_40GeV_bg_norm_1_5mm)));
eff_40 = fminsearch(delta_40, 0.9)
eff_40_1_5mm = fminsearch(delta_40_1_5mm, 0.9)

% eff_40 = 1;
% eff_40_1_5mm = 1;
% 
delta_20 = @(k) norm((counts_dat_20GeV_amorph_norm_tot - counts_dat_20GeV_bg_norm_tot - k .* (counts_sim_20GeV_amorph_bg_norm - counts_sim_20GeV_bg_norm)));
delta_20_1_5mm = @(k) norm(counts_dat_20GeV_amorph_1_5mm_norm_tot - counts_dat_20GeV_bg_1_5mm_norm_tot - k .* (counts_sim_20GeV_amorph_bg_norm_1_5mm - counts_sim_20GeV_bg_norm_1_5mm));
eff_20 = fminsearch(delta_20, 0.9)
eff_20_1_5mm = fminsearch(delta_20_1_5mm, 0.9)

% eff_20 = 1;
% eff_20_1_5mm = 1;

delta_80 = @(k) norm(counts_dat_80GeV_amorph_norm_tot - counts_dat_80GeV_bg_norm - k .* ( (counts_sim_80GeV_amorph_bg_norm - counts_sim_80GeV_bg_norm)));
delta_80_1_5mm = @(k) norm(counts_dat_80GeV_amorph_norm_1_5mm_tot - counts_dat_80GeV_bg_norm_1_5mm_tot - k .* ( (counts_sim_80GeV_amorph_bg_norm - counts_sim_80GeV_bg_norm)));
eff_80 = fminsearch(delta_80, 0.9)
eff_80_1_5mm = fminsearch(delta_80_1_5mm, 0.9)

% eff_80 = 1;
% eff_80_1_5mm = 1;


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

[ha, ~] = tight_subplot(3,2,[.12 .04],[.05 .05],[.07 .07]);

energy = linspace(0, 20, 20);

axes(ha(1));
hold on
box on
title('20GeV e- ; 1mm C','fontsize',22,'interpreter','latex')
errorbar(energy, counts_dat_20GeV_amorph_norm_tot - counts_dat_20GeV_bg_norm_tot,energy.*sqrt(counts_dat_20GeV_amorph_tot)/NEvents_20GeV_amorph_tot +  energy.*sqrt(counts_dat_20GeV_bg_tot)/NEvents_20GeV_dat_bg_tot,'s','markerfacecolor',colors(1,:))
plot(energy, eff_20*counts_sim_20GeV_amorph_norm, '-','linewidth',2.5)
xlabel('Energy [GeV]','fontsize',22,'interpreter','latex');ylabel('dP/dE [1/mm]','fontsize',22,'interpreter','latex');
legend({'Amorphous - background data','Amorphous only simulation'},'location','northwest','fontsize',16,'interpreter','latex')
xticklabels('auto'); yticklabels('auto')
set(gca, 'FontSize', 18)
ylim([-1.8e-4, 4e-4])
grid on

axes(ha(2));
hold on
box on
title('20GeV e- ; 1.5mm C','fontsize',22,'interpreter','latex')
errorbar(energy, counts_dat_20GeV_amorph_1_5mm_norm_tot/1.5 - counts_dat_20GeV_bg_1_5mm_norm_tot/1.5, energy.*sqrt(counts_dat_20GeV_amorph_1_5mm_tot)/NEvents_20GeV_dat_amorph_1_5mm_tot/1.5 +  energy.*sqrt(counts_dat_20GeV_bg_1_5mm_tot)/NEvents_20GeV_dat_bg_1_5mm_tot/1.5,'s','markerfacecolor',colors(1,:))
plot(energy, eff_20_1_5mm*counts_sim_20GeV_amorph_1_5mm_norm/1.5, '-','linewidth',2.5)
xlabel('Energy [GeV]','fontsize',22,'interpreter','latex');ylabel('dP/dE [1/mm]','fontsize',22,'interpreter','latex')
legend({'Amorphous - background data','Amorphous only simulation'},'location','northwest','fontsize',16,'interpreter','latex')
xticklabels('auto'); yticklabels('auto')
set(gca, 'FontSize', 18)
energy = linspace(0, 40, 40);
ylim([-1.8e-4, 4e-4])
grid on
ax = gca;
ax.YAxisLocation = 'right';

axes(ha(4));
hold on
box on
title('40GeV e- ; 1.5mm C','fontsize',22,'interpreter','latex')
errorbar(energy, counts_dat_40GeV_amorph_norm_tot_1_5mm/1.5 - counts_dat_40GeV_bg_norm_1_5mm/1.5,energy.*sqrt(counts_dat_40GeV_amorph_tot_1_5mm)/NEvents_40GeV_amorph_tot_1_5mm/1.5 +  energy.*sqrt(counts_dat_40GeV_bg_1_5mm)/2771767/1.5,'s','markerfacecolor',colors(1,:))
plot(energy, eff_40_1_5mm*counts_sim_40GeV_amorph_1_5mm_norm/1.5, '-','linewidth',2.5)
xlabel('Energy [GeV]','fontsize',22,'interpreter','latex');ylabel('dP/dE [1/mm]','fontsize',22,'interpreter','latex')
legend({'Amorphous - background data','Amorphous only simulation'},'location','northwest','fontsize',16,'interpreter','latex')
xticklabels('auto'); yticklabels('auto')
set(gca, 'FontSize', 18)
ylim([-1.8e-4, 6e-4])
grid on
ax = gca;
ax.YAxisLocation = 'right';

axes(ha(3));
hold on
box on
title('40GeV e- ; 1mm C','fontsize',22,'interpreter','latex')
errorbar(energy, counts_dat_40GeV_amorph_norm_tot - counts_dat_40GeV_bg_norm,energy.*sqrt(counts_dat_40GeV_amorph_tot)/NEvents_40GeV_amorph_tot +  energy.*sqrt(counts_dat_40GeV_bg_1)/3029506,'s','markerfacecolor',colors(1,:))
plot(energy, eff_40*counts_sim_40GeV_amorph_norm, '-','linewidth',2.5)
xlabel('Energy [GeV]','fontsize',22,'interpreter','latex');ylabel('dP/dE [1/mm]','fontsize',22,'interpreter','latex')
legend({'Amorphous - background data','Amorphous only simulation'},'location','northwest','fontsize',16,'interpreter','latex')
xticklabels('auto'); yticklabels('auto')
set(gca, 'FontSize', 18)
ylim([-1.8e-4, 6e-4])
grid on

energy = linspace(0, 80, 40);

axes(ha(5));

hold on
box on
title('80GeV e- ; 1mm C','fontsize',22,'interpreter','latex')
errorbar(energy, counts_dat_80GeV_amorph_norm_tot - counts_dat_80GeV_bg_norm,energy.*sqrt(counts_dat_80GeV_amorph_tot)/NEvents_80GeV_amorph_tot +  energy.*sqrt(counts_dat_80GeV_bg_tot)/NEvents_80GeV_dat_bg_tot,'s','markerfacecolor',colors(1,:))
plot(energy, eff_80*counts_sim_80GeV_amorph_norm, '-','linewidth',2.5)
legend({'Amorphous - background data','Amorphous only simulation'},'location','northwest','fontsize',16,'interpreter','latex')
xlabel('Energy [GeV]','fontsize',22,'interpreter','latex');ylabel('dP/dE [1/mm]','fontsize',22,'interpreter','latex')
xticklabels('auto'); yticklabels('auto')
set(gca, 'FontSize', 18)
ylim([-2e-4, 12e-4])
grid on

axes(ha(6));
hold on
box on
title('80GeV e- ; 1.5mm C','fontsize',22,'interpreter','latex')
errorbar(energy, counts_dat_80GeV_amorph_norm_1_5mm_tot/1.5 - counts_dat_80GeV_bg_norm_1_5mm_tot/1.5,energy.*sqrt(counts_dat_80GeV_amorph_1_5mm_tot)/NEvents_80GeV_amorph_1_5mm_tot/1.5 +  energy.*sqrt(counts_dat_80GeV_bg_1_5mm_tot)/NEvents_80GeV_dat_bg_1_5mm_tot/1.5,'s','markerfacecolor',colors(1,:))
plot(energy, eff_80*counts_sim_80GeV_amorph_1_5mm_norm/1.5, '-','linewidth',2.5)
xlabel('Energy [GeV]','fontsize',22,'interpreter','latex');ylabel('dP/dE [1/mm]','fontsize',22,'interpreter','latex')
legend({'Amorphous - background data','Amorphous only simulation'},'location','northwest','fontsize',16,'interpreter','latex')
xticklabels('auto'); yticklabels('auto')
set(gca, 'FontSize', 18)
ylim([-2e-4, 12e-4])
grid on
ax = gca;
ax.YAxisLocation = 'right';

% 
f = figure;
[ha, ~] = tight_subplot(3,2,[.12 .04],[.05 .05],[.07 .07]);

energy = linspace(0, 20, 20);

axes(ha(2));
hold on
box on
title('20GeV e- ; 1.5mm C','fontsize',22,'interpreter','latex')
errorbar(energy, counts_dat_20GeV_amorph_1_5mm_norm_tot/1.5, energy.*sqrt(counts_dat_20GeV_amorph_1_5mm_tot)/1.5/NEvents_20GeV_dat_amorph_1_5mm_tot,'^','MarkerFaceColor',colors(1,:))
errorbar(energy, counts_dat_20GeV_bg_1_5mm_norm_tot/1.5, energy.*sqrt(counts_dat_20GeV_bg_1_5mm_tot)/1.5/NEvents_20GeV_dat_bg_tot,'p','MarkerFaceColor',colors(2,:))
errorbar(energy, counts_dat_20GeV_amorph_1_5mm_norm_tot/1.5 - counts_dat_20GeV_bg_1_5mm_norm_tot/1.5, energy.*sqrt(counts_dat_20GeV_amorph_1_5mm_tot)/1.5/NEvents_20GeV_dat_amorph_1_5mm_tot + energy.*sqrt(counts_dat_20GeV_bg_1_5mm_tot)/1.5/NEvents_20GeV_dat_bg_tot,'s','MarkerFaceColor',colors(3,:))
plot(energy, eff_20_1_5mm * counts_sim_20GeV_bg_norm_1_5mm/1.5,':','linewidth',2.5,'color',colors(2,:))
plot(energy, eff_20_1_5mm * counts_sim_20GeV_amorph_bg_norm_1_5mm/1.5,'--','linewidth',2.5,'color',colors(1,:))
plot(energy, eff_20_1_5mm * counts_sim_20GeV_amorph_bg_norm_1_5mm/1.5 - eff_20_1_5mm * counts_sim_20GeV_bg_norm_1_5mm/1.5,'-','linewidth',2.5,'color',colors(3,:))
xlabel('Energy [GeV]','fontsize',22,'interpreter','latex');ylabel('dP/dE [1/mm]','fontsize',22,'interpreter','latex');
legend({'Amorph run','Background run','Amorph - background'},'location','northwest','fontsize',18,'interpreter','latex')
xticklabels('auto'); yticklabels('auto')
ylim([-2e-4, 12e-4])
grid on
set(gca, 'FontSize', 18)
ax = gca;
ax.YAxis.Exponent = -3;
ax.YAxisLocation = 'right';

axes(ha(1));
hold on
box on
title('20GeV e- ; 1.0mm C','fontsize',22,'interpreter','latex')
errorbar(energy, counts_dat_20GeV_amorph_norm_tot,energy.*sqrt(counts_dat_20GeV_amorph_tot)/NEvents_20GeV_amorph_tot,'^','MarkerFaceColor',colors(1,:))
errorbar(energy, counts_dat_20GeV_bg_norm_tot,energy.*sqrt(counts_dat_20GeV_bg_tot)/NEvents_20GeV_dat_bg_tot,'p','MarkerFaceColor',colors(2,:))
errorbar(energy, counts_dat_20GeV_amorph_norm_tot-counts_dat_20GeV_bg_norm_tot,energy.*sqrt(counts_dat_20GeV_amorph_tot)/NEvents_20GeV_amorph_tot+energy.*sqrt(counts_dat_20GeV_bg_tot)/NEvents_20GeV_dat_bg_tot,'s','MarkerFaceColor',colors(3,:))
plot(energy, eff_20*counts_sim_20GeV_amorph_bg_norm,'--','linewidth',2.5,'color',colors(1,:))
plot(energy, eff_20*counts_sim_20GeV_bg_norm,':','linewidth',2.5,'color',colors(2,:))
plot(energy, -eff_20*counts_sim_20GeV_bg_norm + eff_20*counts_sim_20GeV_amorph_bg_norm, '-','linewidth',2.5,'color',colors(3,:))
xlabel('Energy [GeV]','fontsize',22,'interpreter','latex');ylabel('dP/dE [1/mm]','fontsize',22,'interpreter','latex');
legend({'Amorph run','Background run','Amorph - background'},'location','northwest','fontsize',18,'interpreter','latex')
xticklabels('auto'); yticklabels('auto')
ylim([-2e-4, 15e-4])
grid on
set(gca, 'FontSize', 18)
ax = gca;
ax.YAxis.Exponent = -3;

energy = linspace(0, 40, 40);

axes(ha(4));
hold on
box on
title('40GeV e- ; 1.5mm C','fontsize',22,'interpreter','latex')
errorbar(energy, counts_dat_40GeV_amorph_norm_tot_1_5mm/1.5, energy.*sqrt(counts_dat_40GeV_amorph_tot_1_5mm)/NEvents_40GeV_amorph_tot_1_5mm/1.5 ,'^','MarkerFaceColor',colors(1,:))
errorbar(energy, counts_dat_40GeV_bg_norm_1_5mm/1.5, energy.*sqrt(counts_dat_40GeV_bg_1_5mm)/3023132/1.5,'p','MarkerFaceColor',colors(2,:))
errorbar(energy, counts_dat_40GeV_amorph_norm_tot_1_5mm/1.5 - counts_dat_40GeV_bg_norm_1_5mm/1.5,energy.*sqrt(counts_dat_40GeV_amorph_tot_1_5mm)/NEvents_40GeV_amorph_tot_1_5mm/1.5 +  energy.*sqrt(counts_dat_40GeV_bg_1_5mm)/3023132/1.5,'s','MarkerFaceColor',colors(3,:))
plot(energy, eff_40_1_5mm * counts_sim_40GeV_bg_norm_1_5mm/1.5,':','linewidth',2.5,'color',colors(2,:))
plot(energy, eff_40_1_5mm * counts_sim_40GeV_amorph_bg_norm_1_5mm/1.5,'--','linewidth',2.5,'color',colors(1,:))
plot(energy, -eff_40_1_5mm * counts_sim_40GeV_bg_norm_1_5mm/1.5 + eff_40_1_5mm * counts_sim_40GeV_amorph_bg_norm_1_5mm/1.5,'-','linewidth',2.5,'color',colors(3,:))
xlabel('Energy [GeV]','fontsize',22,'interpreter','latex');ylabel('dP/dE [1/mm]','fontsize',22,'interpreter','latex');
legend({'Amorph run','Background run','Amorph - background'},'fontsize',18,'interpreter','latex')
xticklabels('auto'); yticklabels('auto')
ylim([-0.5e-4, 1.5e-3])
grid on
set(gca, 'FontSize', 18)
ax = gca;
ax.YAxis.Exponent = -3;
ax.YAxisLocation = 'right';

axes(ha(3));
hold on
box on
title('40GeV e- ; 1.0mm C','fontsize',22,'interpreter','latex')
errorbar(energy, counts_dat_40GeV_amorph_norm_tot, energy.*sqrt(counts_dat_40GeV_amorph_tot)/NEvents_40GeV_amorph_tot ,'^','MarkerFaceColor',colors(1,:))
errorbar(energy, counts_dat_40GeV_bg_norm, energy.*sqrt(counts_dat_40GeV_bg_1)/3029506,'p','MarkerFaceColor',colors(2,:))
errorbar(energy, counts_dat_40GeV_amorph_norm_tot - counts_dat_40GeV_bg_norm,energy.*sqrt(counts_dat_40GeV_amorph_tot)/NEvents_40GeV_amorph_tot +  energy.*sqrt(counts_dat_40GeV_bg_1)/3029506,'s','MarkerFaceColor',colors(3,:))
plot(energy, eff_40*counts_sim_40GeV_bg_norm,':','linewidth',2.5,'color',colors(2,:))
plot(energy, eff_40*counts_sim_40GeV_amorph_bg_norm,'--','linewidth',2.5,'color',colors(1,:))
plot(energy, eff_40*counts_sim_40GeV_amorph_bg_norm - eff_40*counts_sim_40GeV_bg_norm,'-','linewidth',2.5,'color',colors(3,:))
xlabel('Energy [GeV]','fontsize',22,'interpreter','latex');ylabel('dP/dE [1/mm]','fontsize',22,'interpreter','latex');
legend({'Amorph run','Background run','Amorph - background'},'fontsize',18,'interpreter','latex')
xticklabels('auto'); yticklabels('auto')
ylim([-0.5e-4, 3e-3])
grid on
set(gca, 'FontSize', 18)
ax = gca;
ax.YAxis.Exponent = -3;

energy = linspace(0, 80, 40);

axes(ha(5));
hold on
box on
title('80GeV e- ; 1.0mm C','fontsize',22,'interpreter','latex')
errorbar(energy, counts_dat_80GeV_amorph_norm_tot, energy.*sqrt(counts_dat_80GeV_amorph_tot)/NEvents_80GeV_amorph_tot ,'^','MarkerFaceColor',colors(1,:))
errorbar(energy, counts_dat_80GeV_bg_norm, energy.*sqrt(counts_dat_80GeV_amorph_tot)/NEvents_80GeV_amorph_tot,'p','MarkerFaceColor',colors(2,:))
errorbar(energy, counts_dat_80GeV_amorph_norm_tot - counts_dat_80GeV_bg_norm, energy.*sqrt(counts_dat_80GeV_amorph_tot)/NEvents_80GeV_amorph_tot,'s','MarkerFaceColor',colors(3,:))
plot(energy,eff_80*counts_sim_80GeV_bg_norm ,':','linewidth',2.5,'color',colors(2,:))
plot(energy, eff_80*counts_sim_80GeV_amorph_bg_norm, '--','linewidth',2.5,'color',colors(1,:))
plot(energy, eff_80*counts_sim_80GeV_amorph_bg_norm - eff_80*counts_sim_80GeV_bg_norm,'-','linewidth',2.5,'color',colors(3,:))
legend({'Amorph run','Background run','Amorph - background'},'fontsize',18,'interpreter','latex')
xlabel('Energy [GeV]','fontsize',22,'interpreter','latex');ylabel('dP/dE [1/mm]','fontsize',22,'interpreter','latex');
xticklabels('auto'); yticklabels('auto')
ylim([-1e-4, 3e-3])
grid on
set(gca, 'FontSize', 18)
ax = gca;
ax.YAxis.Exponent = -3;

axes(ha(6));
hold on
box on
title('80GeV e- ; 1.5mm C','fontsize',22,'interpreter','latex')
errorbar(energy, counts_dat_80GeV_amorph_norm_1_5mm_tot/1.5, energy.*sqrt(counts_dat_80GeV_amorph_1_5mm_tot)/NEvents_80GeV_amorph_1_5mm_tot/1.5 ,'^','MarkerFaceColor',colors(1,:))
errorbar(energy, counts_dat_80GeV_bg_norm_1_5mm_tot/1.5, energy.*sqrt(counts_dat_80GeV_bg_1_5mm_tot)/NEvents_80GeV_dat_bg_1_5mm_tot/1.5,'p','MarkerFaceColor',colors(2,:))
errorbar(energy, counts_dat_80GeV_amorph_norm_1_5mm_tot/1.5 - counts_dat_80GeV_bg_norm_1_5mm_tot/1.5,energy.*sqrt(counts_dat_80GeV_amorph_1_5mm_tot)/NEvents_80GeV_amorph_1_5mm_tot/1.5 +  energy.*sqrt(counts_dat_80GeV_bg_1_5mm_tot)/NEvents_80GeV_dat_bg_1_5mm_tot/1.5,'s','MarkerFaceColor',colors(3,:))
plot(energy, eff_80_1_5mm * counts_sim_80GeV_bg_1_5mm_norm/1.5 ,':','linewidth',2.5,'color',colors(2,:))
plot(energy, eff_80_1_5mm * counts_sim_80GeV_amorph_bg_1_5mm_norm/1.5,'--','linewidth',2.5,'color',colors(1,:))
plot(energy, -eff_80_1_5mm * counts_sim_80GeV_bg_1_5mm_norm/1.5 + eff_80_1_5mm * counts_sim_80GeV_amorph_bg_1_5mm_norm/1.5,'-','linewidth',2.5,'color',colors(3,:))
xlabel('Energy [GeV]','fontsize',22,'interpreter','latex');ylabel('dP/dE [1/mm]','fontsize',22,'interpreter','latex');
legend({'Amorph run','Background run','Amorph - background'},'fontsize',18,'interpreter','latex')
xticklabels('auto'); yticklabels('auto')
ylim([-1e-4, 20e-4])
grid on
set(gca, 'FontSize', 18)
ax = gca;
ax.YAxisLocation = 'right';
ax.YAxis.Exponent = -3;
 

f = figure;
[ha, ~] = tight_subplot(3,2,[.12 .04],[.05 .05],[.07 .07]);

energy = linspace(0, 40, 40);

% figure
axes(ha(4));
hold on
box on
title('40GeV e- ; 1.5mm C','Interpreter','latex','fontsize',22,'interpreter','latex')
errorbar(energy, counts_dat_40GeV_aligned_norm_tot_1_5mm/1.5 - counts_dat_40GeV_bg_norm_1_5mm/1.5,energy.*sqrt(counts_dat_40GeV_aligned_tot_1_5mm)/NEvents_40GeV_aligned_tot_1_5mm/1.5 +   energy.*sqrt(counts_dat_40GeV_bg_1_5mm)/3023132/1.5,'o','MarkerFaceColor',colors(1,:))
plot(energy, -eff_40_1_5mm * counts_sim_40GeV_bg_norm_1_5mm/1.5 + eff_40_1_5mm * counts_sim_40GeV_aligned_1_5mm_norm/1.5,'-','linewidth',2.5)
plot(energy, -eff_40_1_5mm * counts_sim_40GeV_bg_norm_1_5mm/1.5 + eff_40_1_5mm * counts_sim_40GeV_aligned_1_5mm_norm_woshot/1.5,'--','linewidth',2.5)
plot(energy, -eff_40_1_5mm * counts_sim_40GeV_bg_norm_1_5mm/1.5 + eff_40_1_5mm * counts_sim_40GeV_aligned_1_5mm_norm_worr/1.5,':','linewidth',2.5)
xlabel('Energy [GeV]','fontsize',22,'interpreter','latex');ylabel('dP/dE [1/mm]','fontsize',22,'interpreter','latex')
legend({'Aligned data','Full','No Schott','No RR'},'Interpreter','latex')
set(gca, 'FontSize', 18)
xticklabels('auto'); yticklabels('auto')
ylim([-9E-5,5e-3]); 
grid on
ax = gca;
ax.YAxisLocation = 'right';
ax.YAxis.Exponent = -3;


% figure
axes(ha(3));
hold on
box on
title('40GeV e- ; 1mm C','Interpreter','latex')
errorbar(energy, counts_dat_40GeV_aligned_norm_tot - counts_dat_40GeV_bg_norm,energy.*sqrt(counts_dat_40GeV_aligned_tot)/NEvents_40GeV_aligned_tot +  energy.*sqrt(counts_dat_40GeV_bg_1)/3029506,'o','MarkerFaceColor',colors(1,:))
plot(energy, -eff_40 * counts_sim_40GeV_bg_norm + eff_40 * counts_sim_40GeV_aligned_norm,'-','linewidth',2.5)
plot(energy, -eff_40 * counts_sim_40GeV_bg_norm + eff_40 * counts_sim_40GeV_aligned_norm_woshot,'--','linewidth',2.5)
plot(energy, -eff_40 * counts_sim_40GeV_bg_norm + eff_40 * counts_sim_40GeV_aligned_norm_worr,':','linewidth',2.5)
xlabel('Energy [GeV]','fontsize',22,'interpreter','latex');ylabel('dP/dE [1/mm]','fontsize',22,'interpreter','latex')
legend({'Aligned data','Full','No Schott','No RR'},'Interpreter','latex')
set(gca, 'FontSize', 18)
ax = gca;
ax.FontSize = 18;
ax.YAxis.Exponent = -3;
xticklabels('auto'); yticklabels('auto')
ylim([-1.5E-4,10e-3]); 
grid on

energy = linspace(0, 80, 40);

% figure;
axes(ha(5));
hold on
box on
title('80GeV e- ; 1mm C','Interpreter','latex')
errorbar(energy, counts_dat_80GeV_aligned_norm_tot - counts_dat_80GeV_bg_norm,energy.*sqrt(counts_dat_80GeV_aligned_tot)/NEvents_80GeV_aligned_tot +   energy.*sqrt(counts_dat_80GeV_bg_tot)/NEvents_80GeV_dat_bg_tot,'o','MarkerFaceColor',colors(1,:))
plot(energy, -eff_80 * counts_sim_80GeV_bg_norm + eff_80 * counts_sim_80GeV_aligned_norm,'-','linewidth',2.5)
plot(energy, -eff_80 * counts_sim_80GeV_bg_norm + eff_80 * counts_sim_80GeV_aligned_norm_woshot,'--','linewidth',2.5)
plot(energy, -eff_80 * counts_sim_80GeV_bg_norm + eff_80 * counts_sim_80GeV_aligned_norm_worr,':','linewidth',2.5)
xlabel('Energy [GeV]','fontsize',22,'interpreter','latex');ylabel('dP/dE [1/mm]','fontsize',22,'interpreter','latex')
legend({'Aligned data','Full','No Schott','No RR'},'Interpreter','latex')
xticklabels('auto'); yticklabels('auto')
ylim([-6E-5,20e-3]); 
set(gca, 'FontSize', 18)
ax = gca;
ax.YAxis.Exponent = -3;
grid on

% figure
axes(ha(6));
hold on
box on
title('80GeV e- ; 1.5mm C','Interpreter','latex')
errorbar(energy, counts_dat_80GeV_aligned_norm_1_5mm_tot/1.5 - counts_dat_80GeV_bg_norm_1_5mm_tot/1.5,energy.*sqrt(counts_dat_80GeV_aligned_1_5mm_tot)/NEvents_80GeV_aligned_1_5mm_tot/1.5 +   energy.*sqrt(counts_dat_80GeV_bg_1_5mm_tot)/NEvents_80GeV_dat_bg_1_5mm_tot/1.5,'o','MarkerFaceColor',colors(1,:))
plot(energy, -eff_80_1_5mm * counts_sim_80GeV_bg_1_5mm_norm/1.5 + eff_80_1_5mm * counts_sim_80GeV_aligned_1_5mm_norm/1.5,'-','linewidth',2.5)
plot(energy, -eff_80_1_5mm * counts_sim_80GeV_bg_1_5mm_norm/1.5 + eff_80_1_5mm * counts_sim_80GeV_aligned_1_5mm_norm_woshot/1.5,'--','linewidth',2.5)
plot(energy, -eff_80_1_5mm * counts_sim_80GeV_bg_1_5mm_norm/1.5 + eff_80_1_5mm * counts_sim_80GeV_aligned_1_5mm_norm_worr/1.5,':','linewidth',2.5)
xlabel('Energy [GeV]','fontsize',22,'interpreter','latex');ylabel('dP/dE [1/mm]','fontsize',22,'interpreter','latex')
legend({'Aligned data','Full','No Schott','No RR'},'Interpreter','latex')
xticklabels('auto'); yticklabels('auto')
ylim([0,20e-3]); 
set(gca, 'FontSize', 18)
grid on
ax = gca;
ax.YAxisLocation = 'right';
ax.YAxis.Exponent = -3;

energy = linspace(0, 20, 20);

% figure
axes(ha(2));
hold on
box on
title('20GeV e- ; 1.5mm C','Interpreter','latex')
errorbar(energy, counts_dat_20GeV_aligned_norm_tot_1_5mm/1.5 - counts_dat_20GeV_bg_1_5mm_norm_tot/1.5,energy.*sqrt(counts_dat_20GeV_aligned_tot_1_5mm)/NEvents_20GeV_aligned_tot_1_5mm/1.5 +   energy.*sqrt(counts_dat_20GeV_bg_1_5mm_tot)/NEvents_20GeV_dat_bg_1_5mm_tot/1.5,'o','MarkerFaceColor',colors(1,:))
plot(energy, -eff_20_1_5mm * counts_sim_20GeV_bg_norm_1_5mm/1.5 + eff_20_1_5mm * counts_sim_20GeV_aligned_1_5mm_norm/1.5,'-','linewidth',2.5)
plot(energy, -eff_20_1_5mm * counts_sim_20GeV_bg_norm_1_5mm/1.5 + eff_20_1_5mm * counts_sim_20GeV_aligned_1_5mm_norm_woshot/1.5,'--','linewidth',2.5)
plot(energy, -eff_20_1_5mm * counts_sim_20GeV_bg_norm_1_5mm/1.5 + eff_20_1_5mm * counts_sim_20GeV_aligned_1_5mm_norm_worr/1.5,':','linewidth',2.5)
xlabel('Energy [GeV]','fontsize',22,'interpreter','latex');ylabel('dP/dE [1/mm]','fontsize',22,'interpreter','latex')
legend({'Aligned data','Full','No Schott','No RR'},'Interpreter','latex')
set(gca, 'FontSize', 18)
xticklabels('auto'); yticklabels('auto')
ylim([-5E-5,10e-4]); 
grid on
ax = gca;
ax.YAxisLocation = 'right';
ax.YAxis.Exponent = -3;

% figure
axes(ha(1));
hold on
box on
title('20GeV e- ; 1mm C','Interpreter','latex')
errorbar(energy, counts_dat_20GeV_aligned_norm_tot - counts_dat_20GeV_bg_norm_tot,energy.*sqrt(counts_dat_20GeV_aligned_tot_1_5mm)/NEvents_20GeV_aligned_tot_1_5mm/1.5 +   energy.*sqrt(counts_dat_20GeV_bg_1_5mm_tot)/NEvents_20GeV_dat_bg_1_5mm_tot/1.5,'o','MarkerFaceColor',colors(1,:))
plot(energy, -eff_20 * counts_sim_20GeV_bg_norm + eff_20 * counts_sim_20GeV_aligned_norm,'-','linewidth',2.5)
plot(energy, -eff_20 * counts_sim_20GeV_bg_norm + eff_20 * counts_sim_20GeV_aligned_norm_woshot,'--','linewidth',2.5)
plot(energy, -eff_20 * counts_sim_20GeV_bg_norm + eff_20 * counts_sim_20GeV_aligned_norm_worr,':','linewidth',2.5)
xlabel('Energy [GeV]','fontsize',22,'interpreter','latex');ylabel('dP/dE [1/mm]','fontsize',22,'interpreter','latex')
legend({'Aligned data','Full','No Schott','No RR'},'Interpreter','latex')
set(gca, 'FontSize', 18)
xticklabels('auto'); yticklabels('auto')
ylim([-9E-5,10e-4]); 
grid on
ax = gca;
ax.YAxis.Exponent = -3;


% figure;
% axes(ha(5));
% hold on
% box on
% title('Enhancement e- ; 1mm C','Interpreter','latex')
% energy = linspace(0, 20, 20);
% plot(energy, (counts_dat_20GeV_aligned_norm_tot - counts_dat_20GeV_bg_norm_tot)./(counts_dat_20GeV_amorph_norm_tot - counts_dat_20GeV_bg_norm_tot))
% figure;
% energy = linspace(0, 40, 40);
% plot(energy, (counts_dat_40GeV_aligned_norm_tot - counts_dat_40GeV_bg_norm)./(counts_dat_40GeV_amorph_norm_tot - counts_dat_40GeV_bg_norm))
% figure;
% energy = linspace(0, 80, 40);
% plot(energy, (counts_dat_80GeV_aligned_norm_tot - counts_dat_80GeV_bg_norm)./(counts_dat_80GeV_amorph_norm_tot - counts_dat_80GeV_bg_norm))
