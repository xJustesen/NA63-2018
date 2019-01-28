clear all; close all
datpath = '/home/christian/Dropbox/Cern2018Experiment/spectre/';

%% AVG INTENSITY SPECTRUM
I_80_woschot = load(strcat(datpath,'sum_angles1.0mm80GeVbaiernoSHbfac.txt'));
I_80_worr    = load(strcat(datpath,'sum_angles1.0mm80GeVbaiernoRR.txt'));
I_80_full    = load(strcat(datpath,'sum_angles1.0mm80GeVbaierRRbfac.txt'));


I_40_woschot = load(strcat(datpath,'sum_angles1.0mm40GeVbaiernoSHbfac.txt'));
I_40_worr    = load(strcat(datpath,'sum_angles1.0mm40GeVbaiernoRR.txt'));
I_40_full    = load(strcat(datpath,'sum_angles1.0mm40GeVbaierRRbfac.txt'));

I_20_woschot = load(strcat(datpath,'sum_angles1.0mm20GeVbaiernoSHbfac.txt'));
I_20_worr    = load(strcat(datpath,'sum_angles1.0mm20GeVbaiernoRR.txt'));
I_20_full    = load(strcat(datpath,'sum_angles1.0mm20GeVbaierRRbfac.txt'));


colors =    [       0    0.4470    0.7410
                0.8500    0.3250    0.0980
                0.9290    0.6940    0.1250
                0.4940    0.1840    0.5560
                0.4660    0.6740    0.1880
                0.3010    0.7450    0.9330
                0.6350    0.0780    0.1840
            ];

f = figure;

hold on
box on
grid on
set(gca, 'FontSize', 14)

plot(I_20_full(:,1), I_20_full(:,2),'-','linewidth',1.5,'color',colors(1,:))
plot(I_20_woschot(:,1), I_20_woschot(:,2),'x','linewidth',1.5,'color',colors(1,:))
plot(I_20_worr(:,1), I_20_worr(:,2),'--','linewidth',1.5,'color',colors(1,:))

plot(I_40_full(:,1), I_40_full(:,2),'-','linewidth',1.5,'color',colors(2,:))
plot(I_40_woschot(:,1), I_40_woschot(:,2),'x','linewidth',1.5,'color',colors(2,:))
plot(I_40_worr(:,1), I_40_worr(:,2),'--','linewidth',1.5,'color',colors(2,:))

plot(I_80_full(:,1), I_80_full(:,2),'-','linewidth',1.5,'color',colors(3,:))
plot(I_80_woschot(:,1), I_80_woschot(:,2),'x','linewidth',1.5,'color',colors(3,:))
plot(I_80_worr(:,1), I_80_worr(:,2),'--','linewidth',1.5,'color',colors(3,:))

xlabel('E [GeV]','fontsize',22,'interpreter','latex')
ylabel('$dP/d\hbar\omega$','fontsize',22,'interpreter','latex')

% hack to manually change the color of legend entries
L(1) = plot(nan, nan, 'k-','linewidth',1.5);
L(2) = plot(nan, nan, 'kx','linewidth',1.5);
L(3) = plot(nan, nan, 'k--','linewidth',1.5);

legend(L,'BKC','BKCnoSchott','BKCnoRR','location','northoutside','orientation','horizontal')

axes('Position',[.55 .45 .3 .3])
box on
grid on
hold on

plot(I_20_full(:,1), I_20_full(:,2),'-','linewidth',1.5,'color',colors(1,:))
plot(I_20_woschot(:,1), I_20_woschot(:,2),'x','linewidth',1.5,'color',colors(1,:))
plot(I_20_worr(:,1), I_20_worr(:,2),'--','linewidth',1.5,'color',colors(1,:))

plot(I_40_full(:,1), I_40_full(:,2),'-','linewidth',1.5,'color',colors(2,:))
plot(I_40_woschot(:,1), I_40_woschot(:,2),'x','linewidth',1.5,'color',colors(2,:))
plot(I_40_worr(:,1), I_40_worr(:,2),'--','linewidth',1.5,'color',colors(2,:))

plot(I_80_full(:,1), I_80_full(:,2),'-','linewidth',1.5,'color',colors(3,:))
plot(I_80_woschot(:,1), I_80_woschot(:,2),'x','linewidth',1.5,'color',colors(3,:))
plot(I_80_worr(:,1), I_80_worr(:,2),'--','linewidth',1.5,'color',colors(3,:))

xlim([0, 20])
ylim([0.1, 0.6])

set(f, 'Units','centimeters','PaperUnits','centimeters', 'PaperSize',[18, 12],'PaperPosition',[0, 0, 18, 12],'Position',[0 0 18 12])
print('/home/christian/Dropbox/Cern2018Experiment/figures/theory_compare', '-dpdf','-r600','-painters')


%% ENERGY vs ANGLE SPECTRUM
dat = load('../crystalSimulations/sum_initials40GeV_full_peak.txt');
datsim = load('/media/christian/Elements/DataAnalysis/radiation/spectres/sum_initials80GeV_1.5mm_full_peak_alt_woshot.txt');
datsim = datsim';

angles = dat(4:5, :)';

dat(1:5, :) = [];
dat = dat';

% i = dat(:,3); a1 = angles(1:100,1); a2 = a1; isim = datsim(:,3);
% i = reshape(i, 100, 100); isim = reshape(isim, 100, 100);
% 
% datsim2 = load('/home/christian/Dropbox/speciale/data/photon_angles.txt');
% N = hist3(datsim2,'Nbins',[100,100]);
% 
% f = figure;
% imagesc(dat)
% 
% f = figure;
% spy(i-isim)
% 
% 
% [A1, A2] = meshgrid(a1, a2);
% 
% f = figure
% hold on
% contourf(A1, A2, i, 1000, 'edgecolor','none');
% colorbar
% xlabel('$\theta\quad [\mathrm{rad}]$','fontsize',15,'interpreter','latex')
% ylabel('$\phi\quad [\mathrm{rad}]$','fontsize',15,'interpreter','latex')
% title('DATA (specific energy)','fontsize',25,'interpreter','latex')
% axis equal
% set(f, 'Units','centimeters','PaperUnits','centimeters', 'PaperSize',[18, 12],'PaperPosition',[0, 0, 18, 12],'Position',[0 0 18 12])
% print('../../figures/angle_spec.pdf', '-dpdf','-r600','-painters')
% 
% f = figure
% hold on
% contourf(A1, A2, N, 1000, 'edgecolor','none');
% colorbar
% xlabel('$\theta\quad [\mathrm{rad}]$','fontsize',15,'interpreter','latex')
% ylabel('$\phi\quad [\mathrm{rad}]$','fontsize',15,'interpreter','latex')
% title('SIMULATION (across all energies)','fontsize',25,'interpreter','latex')
% axis equal
% set(f, 'Units','centimeters','PaperUnits','centimeters', 'PaperSize',[18, 12],'PaperPosition',[0, 0, 18, 12],'Position',[0 0 18 12])
% % print('../../figures/angle_spec.pdf', '-dpdf','-r600','-painters')

