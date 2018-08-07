close all; clc; clear all
datpath = '/home/christian/Dropbox/speciale/data/';
filepath = strcat(datpath,'zclosepos_73.txt');
zpos = load(filepath);

[counts, centers] = hist(zpos*1e-6,300);

figure
stairs(centers, counts);
xlabel('z [meter]')
ylabel('Counts')

