close all; clc; clear all
datpath = '/home/christian/Dropbox/speciale/data/';
% filepath1 = strcat(datpath,'zclosepos_sim_amorphous1_40GeV.txt');
filepath2 = strcat(datpath,'zclosepos_39.txt');

% zpos1 = load(filepath1);
zpos2 = load(filepath2);

% [counts1, centers1] = hist(zpos1*1e-6,300);
[counts2, centers2] = hist(zpos2*1e-6,300);

f = figure;
hold on
box on
errorbar(centers2, counts2, sqrt(counts2));
% plot(centers1, max(counts2)/max(counts1) * counts1,'linewidth',2);
legend('Data','Sim')
xlabel('z [meter]')
ylabel('counts')
set(f, 'Units','centimeters','PaperUnits','centimeters', 'PaperSize',[18, 12],'PaperPosition',[0, 0, 18, 12],'Position',[0 0 18 12])
% print(f,'../../figures/zpos_40GeV.pdf', '-dpdf','-r600','-painters')

