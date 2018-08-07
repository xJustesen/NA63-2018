clear all; clc; close all;
datpath = '/home/christian/Dropbox/speciale/data/';
filepath = strcat(datpath,'angles_73.txt');

hold on
ang_mat = load(filepath);
ang_in = ang_mat(:,1);
ang_out = ang_mat(:,2);
defl = ang_in - ang_out;
bins = linspace(-0.02, 0.02, 750);

[counts_o, centers_o] = hist(ang_out,bins);
[counts_i, centers_i] = hist(ang_in,bins);

figure(1)
hold on
plot(bins, counts_i);
plot(bins, counts_o);
xlabel('ang (rad)'); ylabel('counts');
legend('Angle in','Angle out');
