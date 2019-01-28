close all; clear all

colors =    [       0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840
    ];

datpath = '/home/christian/Documents/cern2018/simdata/';



% dat = load(strcat(datpath, 'M5M6_dist_115.txt'));
% histogram(dat,200);
%

% dat.b = load(strcat(datpath, 'M2M3_dist_18.txt'));
% dat.c = load(strcat(datpath, 'M5M6_dist_18.txt'));
% dat.d = load(strcat(datpath, 'Match_d_foil_18.txt'));
% dat.e = load(strcat(datpath, 'Match_d_MM_18.txt'));
% dat.aa = load(strcat(datpath, 'M1M2_dist_sim_alignment.txt'));
% dat.bb = load(strcat(datpath, 'M2M3_dist_sim_alignment.txt'));
% dat.cc = load(strcat(datpath, 'M5M6_dist_sim_alignment.txt'));
% dat.dd = load(strcat(datpath, 'Match_d_foil_sim_aligned1_80GeV_1.5mm.txt'));
% dat.ee = load(strcat(datpath, 'Match_d_MM_sim_aligned1_80GeV_1.5mm.txt'));
%
% defl = load(strcat(datpath, "yz_defl_18.txt"));
% deflsim = load(strcat(datpath, "yz_defl_sim_alignment.txt"));


axisdir_1mm = [1.7300e-04	-6.1449e-05; -3.4474e-05	-1.1023e-04; 1.1261e-04	-1.5039e-04];
axisdir_1_5mm = [1.2060e-04	-3.7260e-05; -5.7066e-05	-3.7708e-05; 1.2077e-04	-8.0643e-05];
beamdir.A = [2.002e-06	-3.003e-05; 3.4034e-05	-2.6026e-05; 0.00015816 	-8.2082e-05];   % 1.0mm
beamdir.B = [1.8018e-05	-3.8038e-05; 6.006e-06	-1.001e-05; 0.00016617	-3.8038e-05];       % 1.5mm
stdev = [0.00030436	0.00013943	0.00028226	0.00013263; 0.00019894	0.00011562	0.00017684	0.00011222; 0.00012242	7.9916e-05	0.00013093	7.6515e-05];

relaxisdir_1mm = beamdir.A - axisdir_1mm
relaxisdir_1_5mm = beamdir.B - axisdir_1_5mm

run_no = [103, 61, 72, 35, 84, 46];

run_no = 35;


for run = run_no
    
    angles_x = load(strcat(datpath, 'angles_mean_x_',num2str(run),'.txt'));
    [~,idx] = sort(angles_x(:,1));
    angles_x = angles_x(idx, :);
    angx = angles_x(:,1); angxc = angles_x(:,2);
    
    angles_y = load(strcat(datpath, 'angles_mean_y_',num2str(run),'.txt'));
    angy = angles_y(:,1); angyc = angles_y(:,2);
    [~,idx] = sort(angles_y(:,1));
    angles_y = angles_y(idx, :);
    
    photonscountsx = load(strcat(datpath, 'no_photons_x_alt',num2str(run),'.txt'));
    [~,idx] = sort(photonscountsx(:,1));
    photonscountsx = photonscountsx(idx, :);
    angx_alt = photonscountsx(:,1); angxc_alt = photonscountsx(:,2);
    
    photonscountsy = load(strcat(datpath, 'no_photons_y_alt',num2str(run),'.txt'));
    [~,idx] = sort(photonscountsy(:,1));
    photonscountsy = photonscountsy(idx, :);
    angy_alt = photonscountsy(:,1); angyc_alt = photonscountsy(:,2);
    
    
    stdev_x = 0.00012242;
    stdev_y = 7.9916e-05;
    
    if run == 103
        xdir = beamdir.A(1,1);
        ydir = beamdir.A(1,2);
        stdev_x = stdev(1,1);
        stdev_y = stdev(1,2);
        xfilter = angx < 150e-6 & angx > 100e-6;
        yfilter = angy < -50e-6 & angy > -100e-6;
        E = 20;
        d = 1;
    elseif run == 61
        xdir = beamdir.B(1,1);
        ydir = beamdir.B(1,2);
        stdev_x = stdev(1,3);
        stdev_y = stdev(1,4);
        xfilter = angx < 150e-6 & angx > 100e-6;
        yfilter = angy < -50e-6 & angy > -100e-6;
        E = 20; d = 1.5;
    elseif run == 72
        xdir = beamdir.A(2,1);
        ydir = beamdir.A(2,2);
        stdev_x = stdev(2,1);
        stdev_y = stdev(2,2);
        xfilter = angx < 3e-05 & angx > -8.7e-5;
        yfilter = angy < -3.73e-05 & angy > -0.0001813;
        E = 40; d = 1;
    elseif run == 35
        xdir = beamdir.B(2,1);
        ydir = beamdir.B(2,2);
        stdev_x = stdev(2,3);
        stdev_y = stdev(2,4);
        xfilter = angx < -3e-06 & angx > -9.7e-5;
        yfilter = angy < -2.71e-6 & angy > -9.529e-05;
        E = 40; d = 1.5;
    elseif run == 84
        xdir = beamdir.A(3,1);
        ydir = beamdir.A(3,2);
        stdev_x = stdev(3,1);
        stdev_y = stdev(3,2);
        xfilter = angx < 150e-6 & angx > 100e-6;
        yfilter = angy < -9.92e-05 & angy > -0.0001952;
        E = 80; d = 1;
    elseif run == 46
        xdir = beamdir.B(3,1);
        ydir = beamdir.B(3,2);
        stdev_x = stdev(3,3);
        stdev_y = stdev(3,4);
        xfilter = angx < 150e-6 & angx > 100e-6;
        yfilter = angy < -50e-6 & angy > -100e-6;
        E = 80; d = 1.5;
    end
    
    
    fitx = fit(angx(xfilter), angxc(xfilter),'poly2');
    fity = fit(angy(yfilter), angyc(yfilter),'poly2');
    
    axisx = -fitx.p2/(2*fitx.p1);
    axisy = -fity.p2/(2*fity.p1);
    beamx = xdir;
    beamy = ydir;
    relative_x_dir = xdir - axisx;
    relative_y_dir = ydir - axisy;
    
    fitx_alt = fit(angx_alt(xfilter), angxc_alt(xfilter),'poly2');
    fity_alt = fit(angy_alt(yfilter), angyc_alt(yfilter),'poly2');
    
    axisx_alt = -fitx_alt.p2/(2*fitx_alt.p1);
    axisy_alt = -fity_alt.p2/(2*fity_alt.p1);
    relative_x_dir_alt = xdir - axisx_alt;
    relative_y_dir_alt = ydir - axisy_alt;
    
    formatSpec = ['\nrun %i\tEnergy %i GeV \td %f mm \n\tbeam x axis: %e\t crystal x axis (defl): %e\t crystal x axis (counts): %e\t rel. x dir (defl): %e\t rel. x. dir (counts): %e\n' ...
        '\tbeam y axis: %e\t crystal y axis (defl):%e\t crystal y axis (counts): %e\t rel. y dir (defl): %e\t rel. y dir (counts): %e\n'
        ];
    fprintf(formatSpec, run, E, d, xdir, axisx, axisx_alt, relative_x_dir, relative_x_dir_alt, ydir, axisy, axisy_alt, relative_y_dir, relative_y_dir_alt)
    
    
    f = figure;
    subplot(2,2,1)
    grid on
    title('a)','fontsize',22,'interpreter','latex')
    ax = gca;
    ax.YAxis.Exponent = -6;
    ax.XAxis.Exponent = -6;
    hold on
    ylim([20e-6, 120e-6]);
    yl = ylim;
    fill([xdir + stdev_x/2, xdir + stdev_x/2, xdir - stdev_x/2, xdir - stdev_x/2],[yl(1), yl(2), yl(2), yl(1)],'g','FaceAlpha',0.2,'edgecolor','none')
    p(1) = plot(angles_x(:,1), angles_x(:,2),'-','linewidth',2.5);
    p(2) = plot([xdir, xdir], yl,'--','linewidth',2.5,'color',colors(1,:));
    plot(axisx, fitx(axisx),'o','MarkerFaceColor',colors(3,:),'MarkerSize',7)
    legend(p,{'scan','$\mu_\mathrm{hor.}$'},'fontsize',22,'interpreter','latex')
    ylabel('$\langle d\phi_\mathrm{hor.} \rangle$ [rad]','fontsize',22,'interpreter','latex')
    box on
    set(gca, 'FontSize', 18)
    subplot(2,2,2)
    grid on
    title('b)','fontsize',22,'interpreter','latex')
    ax = gca;
    ax.XAxis.Exponent = -6;
    ax.YAxisLocation = 'right';
    hold on
    ylim([0, max(angxc_alt)]);
    yl = ylim;
    fill([xdir + stdev_x/2, xdir + stdev_x/2, xdir - stdev_x/2, xdir - stdev_x/2],[yl(1), yl(2), yl(2), yl(1)],'g','FaceAlpha',0.2,'edgecolor','none')
    p(1) = plot(angx_alt, angxc_alt,'-','linewidth',2.5);
    p(2) = plot([xdir, xdir], yl,'--','linewidth',2.5,'color',colors(1,:));
    plot(axisx_alt, fitx_alt(axisx_alt),'o','MarkerFaceColor',colors(3,:),'MarkerSize',7)
    legend(p,{'scan','$\mu_\mathrm{hor.}$'},'fontsize',22,'interpreter','latex')
    ylabel('$\left(\frac{N_\gamma}{N_E}\right)_\mathrm{cut}$','fontsize',22,'interpreter','latex')
    box on
    set(gca, 'FontSize', 18)
    subplot(2,2,3)
    grid on
    title('c)','fontsize',22,'interpreter','latex')
    ax = gca;
    ax.YAxis.Exponent = -6;
    ax.XAxis.Exponent = -6;
    hold on
    ylim([20e-6, 150e-6]);
    yl = ylim;
    fill([ydir + stdev_y/2, ydir + stdev_y/2, ydir - stdev_y/2, ydir - stdev_y/2],[yl(1), yl(2), yl(2), yl(1)],'g','FaceAlpha',0.2,'edgecolor','none')
    p(1) = plot(angles_y(:,1), angles_y(:,2),'-','linewidth',2.5);
    p(2) = plot([ydir, ydir], yl,'--','linewidth',2.5,'color',colors(1,:));
    plot(axisy, fity(axisy),'o','MarkerFaceColor',colors(3,:),'MarkerSize',7)
    legend(p, {'scan','$\mu_\mathrm{vert.}$'},'fontsize',22,'interpreter','latex')
    xlabel('$\langle d\theta \rangle$ [rad]','fontsize',22,'interpreter','latex');
    ylabel('$\langle d\phi_\mathrm{vert.} \rangle$ [rad]','fontsize',22,'interpreter','latex')
    box on
    set(gca, 'FontSize', 18)
    subplot(2,2,4)
    grid on
    title('d)','fontsize',22,'interpreter','latex')
    ax = gca;
    ax.XAxis.Exponent = -6;
    ax.YAxisLocation = 'right';
    hold on
    ylim([0, max(angyc_alt)]);
    yl = ylim;
    fill([ydir + stdev_y/2, ydir + stdev_y/2, ydir - stdev_y/2, ydir - stdev_y/2],[yl(1), yl(2), yl(2), yl(1)],'g','FaceAlpha',0.2,'edgecolor','none')
    p(1) = plot(angy_alt, angyc_alt,'-','linewidth',2.5);
    p(2) = plot([ydir, ydir], yl,'--','linewidth',2.5,'color',colors(1,:));
    plot(axisy_alt, fity_alt(axisy_alt),'o','MarkerFaceColor',colors(3,:),'MarkerSize',7)
    legend(p, {'scan','$\mu_\mathrm{vert.}$'},'fontsize',22,'interpreter','latex')
    xlabel('$\langle d\theta \rangle$ [rad]','fontsize',22,'interpreter','latex');
    ylabel('$\left(\frac{N_\gamma}{N_E}\right)_\mathrm{cut}$','fontsize',22,'interpreter','latex')
    box on
    set(gca, 'FontSize', 18)
    set(f, 'Units','centimeters','PaperUnits','centimeters', 'PaperSize',[36, 18],'PaperPosition',[0, 0, 36, 18],'Position',[0 0 36 18])
    print('/home/christian/Dropbox/Cern2018Experiment/figures/40GeV_1.5mm_axis_run35.pdf', '-dpdf','-r600','-painters')
    
end
%
% histogram(defl, 1000);
% [defl_c_sim, defl_sim] = histcounts(deflsim, 'numbins', 1000);
%
% f = figure;
% subplot(1,2,1)
% grid on
% hold on
% plot(defl(2:end), defl_c/max(defl_c),'o','linewidth',1.5)
% plot(defl_sim(2:end), defl_c_sim/max(defl_c_sim),'-','linewidth',2.5,'color',colors(1,:))
% xlabel('Angle [rad]','fontsize',22,'interpreter','latex');ylabel('Normalized counts','fontsize',22,'interpreter','latex')
% xlim(1e-3 * [-3.5, 3.5])
% box on
% set(gca, 'FontSize', 18)
% title('a) yz-deflection angle','fontsize',22,'interpreter','latex')
%
% lims = [2000, 100, 100, 100, 100];
% % lims = [inf, inf, inf, inf, inf];
% fields = fieldnames(dat);
%
% for i = 1:length(fieldnames(dat))/2
%
%    field1 = fields{i};
%    field2 = fields{i + length(fields)/2};
%
%    [counts1, edges1] = histcounts(dat.(field1)(dat.(field1) < lims(i)),200);
%    [counts2, edges2] = histcounts(dat.(field2)(dat.(field2) < lims(i)),200);
%
%    counts.(field1) = counts1;
%    counts.(field2) = counts2;
%    edges.(field1) = edges1(2:end);
%    edges.(field2) = edges2(2:end);
%
% end
%
% subplot(1,2,2)
% grid on
% hold on
% box on
% set(gca, 'FontSize', 18)
% plot(edges.a, counts.a/max(counts.a),'o','linewidth',1.5)
% plot(edges.b, counts.b/max(counts.b),'o','linewidth',1.5)
% plot(edges.c, counts.c/max(counts.c),'o','linewidth',1.5)
% plot(edges.aa, counts.aa/max(counts.aa),'-','linewidth',2.5,'color',colors(1,:))
% plot(edges.bb, counts.bb/max(counts.bb),'-','linewidth',2.5,'color',colors(2,:))
% plot(edges.cc, counts.cc/max(counts.cc),'-','linewidth',2.5,'color',colors(3,:))
%
% ylabel('Normalized counts','fontsize',22,'interpreter','latex');xlabel('Distance $\left[\mu \mathrm{m}\right]$','fontsize',22,'interpreter','latex')
% legend({'M3 distance','M4 distance','MIMOSA Magnet distance'},'fontsize',14,'interpreter','latex')
% title('b) matching distances','fontsize',22,'interpreter','latex')
%
% set(f, 'Units','centimeters','PaperUnits','centimeters', 'PaperSize',[36, 12],'PaperPosition',[0, 0, 36, 12],'Position',[0 0 36 12])
% print(f, '../../figures/yz_MM_match_conditions.pdf', '-dpdf','-r600','-painters')

% f = figure;
% subplot(1,2,1)
% hold on
% plot(edges.d, counts.d/max(counts.d),'o-','LineWidth',2.5)
% grid on
% xlabel('Distance $\left[\mu \mathrm{m}\right]$','fontsize',22,'interpreter','latex');ylabel('Normalized counts','fontsize',22,'interpreter','latex')
% box on
% set(gca, 'FontSize', 18)
% title('a) converter foil','fontsize',22,'interpreter','latex');
% subplot(1,2,2)
% grid on
% hold on
% plot(edges.e, counts.e/max(counts.e),'o-','LineWidth',2.5)
% xlabel('Distance $\left[\mu \mathrm{m}\right]$','fontsize',22,'interpreter','latex');ylabel('Normalized counts','fontsize',22,'interpreter','latex')
% box on
% set(gca, 'FontSize', 18)
% title('b) MIMOSA magnet','fontsize',22,'interpreter','latex');
% set(f, 'Units','centimeters','PaperUnits','centimeters', 'PaperSize',[36, 12],'PaperPosition',[0, 0, 36, 12],'Position',[0 0 36 12])
% print(f,'../../figures/magnet_foil_match_dist.pdf', '-dpdf','-r600','-painters')
