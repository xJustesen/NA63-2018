close all; clear all

colors =    [       0    0.4470    0.7410
                0.8500    0.3250    0.0980
                0.9290    0.6940    0.1250
                0.4940    0.1840    0.5560
                0.4660    0.6740    0.1880
                0.3010    0.7450    0.9330
                0.6350    0.0780    0.1840
            ];

datpath = '/home/christian/Dropbox/speciale/data/';

dat.a = load(strcat(datpath, 'M1M2_dist_71.txt'));
dat.b = load(strcat(datpath, 'M2M3_dist_71.txt'));
dat.c = load(strcat(datpath, 'M5M6_dist_71.txt'));
dat.d = load(strcat(datpath, 'Match_d_foil_71.txt'));
dat.e = load(strcat(datpath, 'Match_d_MM_71.txt'));
dat.aa = load(strcat(datpath, 'M1M2_dist_sim_amorphous1_40GeV.txt'));
dat.bb = load(strcat(datpath, 'M2M3_dist_sim_amorphous1_40GeV.txt'));
dat.cc = load(strcat(datpath, 'M5M6_dist_sim_amorphous1_40GeV.txt'));
dat.dd = load(strcat(datpath, 'Match_d_foil_sim_amorphous1_40GeV.txt'));
dat.ee = load(strcat(datpath, 'Match_d_MM_sim_amorphous1_40GeV.txt'));

defl = load(strcat(datpath, "yz_defl_71.txt"));

angles_x = load(strcat(datpath, 'no_photons_x_alt46.txt'));
[~,idx] = sort(angles_x(:,1));
angles_x = angles_x(idx, :);
angx = angles_x(:,1); angxc = angles_x(:,2);

angles_y = load(strcat(datpath, 'no_photons_y_alt46.txt'));
angy = angles_y(:,1); angyc = angles_y(:,2);
[~,idx] = sort(angles_y(:,1));
angles_y = angles_y(idx, :);

stdev_x = 0.00013093;
stdev_y = 7.6515e-05;
xdir = 0.00016617;
ydir = -3.8038e-05;

xfilter = angx < 150e-6 & angx > 100e-6;
yfilter = angy < -50e-6 & angy > -100e-6;

fitx = fit(angx(xfilter), angxc(xfilter),'poly2');
fity = fit(angy(yfilter), angyc(yfilter),'poly2');

axisx = -fitx.p2/(2*fitx.p1)
axisy = -fity.p2/(2*fity.p1)
beamx = xdir
beamy = ydir
relative_x_dir = xdir - axisx
relative_y_dir = ydir - axisx

f = figure;
subplot(2,1,1)
grid on
title('a) x scan','fontsize',22,'interpreter','latex')
ax = gca;
% ax.YAxis.Exponent = -6;
ax.XAxis.Exponent = -6;
hold on
% ylim([0, 1e4]);
yl = ylim;
fill([xdir + stdev_x, xdir + stdev_x, xdir - stdev_x, xdir - stdev_x],[yl(1), yl(2), yl(2), yl(1)],'g','FaceAlpha',0.2,'edgecolor','none')
p(1) = plot(angles_x(:,1), angles_x(:,2),'-','linewidth',2.5);
p(2) = plot([xdir, xdir], yl,'--','linewidth',2.5,'color',colors(1,:));
plot(axisx, fitx(axisx),'o','MarkerFaceColor',colors(3,:),'MarkerSize',7)
legend(p,{'scan','$\mu_x$'},'fontsize',22,'interpreter','latex')
ylabel('Number of photons','fontsize',22,'interpreter','latex')
xlim([-200e-6, 600e-6]);
box on
set(gca, 'FontSize', 18)
subplot(2,1,2)
grid on
title('b) y scan','fontsize',22,'interpreter','latex')
ax = gca;
% ax.YAxis.Exponent = -6;
ax.XAxis.Exponent = -6;
hold on
% ylim([0, 1e4]);
yl = ylim;
fill([ydir + stdev_y, ydir + stdev_y, ydir - stdev_y, ydir - stdev_y],[yl(1), yl(2), yl(2), yl(1)],'g','FaceAlpha',0.2,'edgecolor','none')
p(1) = plot(angles_y(:,1), angles_y(:,2),'-','linewidth',2.5);
p(2) = plot([ydir, ydir], yl,'--','linewidth',2.5,'color',colors(1,:));
plot(axisy, fity(axisy),'o','MarkerFaceColor',colors(3,:),'MarkerSize',7)
legend(p, {'scan','$\mu_y$'},'fontsize',22,'interpreter','latex')
xlabel('$\langle d\theta \rangle$ [rad]','fontsize',22,'interpreter','latex');ylabel('Number of photons','fontsize',22,'interpreter','latex')
xlim([-250e-6, 250e-6]);
box on
set(gca, 'FontSize', 18)
set(f, 'Units','centimeters','PaperUnits','centimeters', 'PaperSize',[18, 12],'PaperPosition',[0, 0, 18, 12],'Position',[0 0 18 12])
% print('../../figures/80GeV_1.5mm_axis_v2.pdf', '-dpdf','-r600','-painters')


[defl_c, defl] = histcounts(defl, 'numbins', 1000);

f = figure;
subplot(1,2,1)
grid on
hold on
plot(defl(2:end), defl_c/max(defl_c),'-','linewidth',2.5)
xlabel('Angle [rad]','fontsize',22,'interpreter','latex');ylabel('Normalized counts','fontsize',22,'interpreter','latex')
xlim(1e-3 * [-3.5, 3.5])
box on
set(gca, 'FontSize', 18)
title('a) yz-deflection angle','fontsize',22,'interpreter','latex')

lims = [2000, 2000, 2000, 100, 100];
% lims = [inf, inf, inf, inf, inf];
fields = fieldnames(dat);

for i = 1:length(fieldnames(dat))/2
   
   field1 = fields{i};
   field2 = fields{i + length(fields)/2};
    
   [counts1, edges1] = histcounts(dat.(field1)(dat.(field1) < lims(i)),200);
   [counts2, edges2] = histcounts(dat.(field2)(dat.(field2) < lims(i)),200);
   
   counts.(field1) = counts1;
   counts.(field2) = counts2;
   edges.(field1) = edges1(2:end);
   edges.(field2) = edges2(2:end);
   
end

subplot(1,2,2)
grid on
hold on
box on
set(gca, 'FontSize', 18)
plot(edges.a, counts.a/max(counts.a),'-','linewidth',2.5)
plot(edges.b, counts.b/max(counts.b),'-','linewidth',2.5)
plot(edges.c, counts.c/max(counts.c),'-','linewidth',2.5)
% plot(edges.aa, counts.aa/max(counts.aa),'r-','MarkerFaceColor','r')
% plot(edges.bb, counts.bb/max(counts.bb),'g-','MarkerFaceColor','g')
% plot(edges.cc, counts.cc/max(counts.cc),'b-','MarkerFaceColor','b')

ylabel('Normalized counts','fontsize',22,'interpreter','latex');xlabel('Distance $\left[\mu \mathrm{m}\right]$','fontsize',22,'interpreter','latex')
legend({'M3 distance','M4 distance','MIMOSA Magnet distance'},'fontsize',14,'interpreter','latex')
title('b) matching distances','fontsize',22,'interpreter','latex')

set(f, 'Units','centimeters','PaperUnits','centimeters', 'PaperSize',[36, 12],'PaperPosition',[0, 0, 36, 12],'Position',[0 0 36 12])
% print(f, '../../figures/yz_MM_match_conditions.pdf', '-dpdf','-r600','-painters')

f = figure;
subplot(1,2,1)
hold on
plot(edges.d, counts.d/max(counts.d),'o-','LineWidth',2.5)
grid on
xlabel('Distance $\left[\mu \mathrm{m}\right]$','fontsize',22,'interpreter','latex');ylabel('Normalized counts','fontsize',22,'interpreter','latex')
box on
set(gca, 'FontSize', 18)
title('a) converter foil','fontsize',22,'interpreter','latex');
subplot(1,2,2)
grid on
hold on
plot(edges.e, counts.e/max(counts.e),'o-','LineWidth',2.5)
xlabel('Distance $\left[\mu \mathrm{m}\right]$','fontsize',22,'interpreter','latex');ylabel('Normalized counts','fontsize',22,'interpreter','latex')
box on
set(gca, 'FontSize', 18)
title('b) MIMOSA magnet','fontsize',22,'interpreter','latex');
set(f, 'Units','centimeters','PaperUnits','centimeters', 'PaperSize',[36, 12],'PaperPosition',[0, 0, 36, 12],'Position',[0 0 36 12])
% print(f,'../../figures/magnet_foil_match_dist.pdf', '-dpdf','-r600','-painters')
