clear all; close all; clc;

data1 = load("/home/christian/Desktop/SimulationTestOutput/dat_amorphous_40GeV_1.5mm.txt") * 6.2415091E9;
data2 = load("/home/christian/Dropbox/Cern2018Experiment/spectre/energy_39.txt") * 6.2415091E9;

E40 = linspace(0, 40, 40);

[n, x] = hist(data1, E40);
[nn, xx] = hist(data2, E40); 

figure
hold on
plot(x, n,'-')
plot(xx, nn,'o')
