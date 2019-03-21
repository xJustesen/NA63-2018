clear all

colors =    [   0         0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840
    ];

datpath = '/home/christian/Documents/cern2018/simdata/';

axisdir_1mm = [1.7300e-04	-6.1449e-05; -3.4474e-05	-1.1023e-04; 1.1261e-04	-1.5039e-04];
axisdir_1_5mm = [1.2060e-04	-3.7260e-05; -5.7066e-05	-3.7708e-05; 1.2077e-04	-8.0643e-05];
beamdir.A = [2.002e-06	-3.003e-05; 3.4034e-05	-2.6026e-05; 0.00015816 	-8.2082e-05];   % 1.0mm
beamdir.B = [1.8018e-05	-3.8038e-05; 6.006e-06	-1.001e-05; 0.00016617	-3.8038e-05];       % 1.5mm
stdev = [0.00030436	0.00013943	0.00028226	0.00013263; 0.00019894	0.00011562	0.00017684	0.00011222; 0.00012242	7.9916e-05	0.00013093	7.6515e-05];

relaxisdir_1mm = beamdir.A - axisdir_1mm
relaxisdir_1_5mm = beamdir.B - axisdir_1_5mm

run_no = [103, 61];%, 72, 35, 84, 46];

% run_no = '111_2017';

run = 0;
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
    elseif run == 35 || run == 43 || run == 32
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
    
end
f = figure;

subplot(2,2,1)
grid on
title('a)','fontsize',22,'interpreter','latex')
ax = gca;
ax.YAxis.Exponent = -6;
ax.XAxis.Exponent = -6;
hold on
ylim([20e-6, 120e-6]);
p(1) = plot(angles_x(:,1), angles_x(:,2),'-','linewidth',2.5);
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
p(1) = plot(angx_alt, angxc_alt,'-','linewidth',2.5);
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
p(1) = plot(angles_y(:,1), angles_y(:,2),'-','linewidth',2.5);
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
p(1) = plot(angy_alt, angyc_alt,'-','linewidth',2.5);
legend(p, {'scan','$\mu_\mathrm{vert.}$'},'fontsize',22,'interpreter','latex')
xlabel('$\langle d\theta \rangle$ [rad]','fontsize',22,'interpreter','latex');
ylabel('$\left(\frac{N_\gamma}{N_E}\right)_\mathrm{cut}$','fontsize',22,'interpreter','latex')
box on
set(gca, 'FontSize', 18)

set(f, 'Units','centimeters','PaperUnits','centimeters', 'PaperSize',[36, 18],'PaperPosition',[0, 0, 36, 18],'Position',[0 0 36 18])
