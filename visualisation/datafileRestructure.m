clear all; clc; close all;
datapath = '/home/christian/Documents/cern2018/code/src/';
% load 2017 data
ang = load("/home/christian/Dropbox/Cern2018Experiment/spectre/beam_direction_4mm.txt");
pos = load("../beamParameters/2017beam_position_50GeV_2mm.txt");

% posx = load(strcat(datapath,'xpos.txt'));
% posy = load(strcat(datapath,'ypos.txt'));
% dx = load(strcat(datapath,'dx.txt'));
% dy = load(strcat(datapath,'dy.txt'));

% [n, edges] = histcounts(posx, 200);

% restrucutre to suit 2018 analysis program
% x = pos(:,3);
% y = pos(:,4);
% xc = pos(:,1);
% yc = pos(:,2);

ax = ang(:,3) * 1e-6;
ay = ang(:,4) * 1e-6;
axc = ang(:,1);
ayc = ang(:,2);

yfit_fakes = fit(ay(ayc < 1.5 * mean(sort(ayc))), ayc(ayc < 1.5 * mean(sort(ayc))),'poly2');
xfit_fakes = fit(ax(axc < 1.5 * mean(sort(axc))), axc(axc < 1.5 * mean(sort(axc))),'poly2');


figure
hold on
plot(ay(ayc < 1.5 * mean(sort(ayc))), ayc(ayc < 1.5 * mean(sort(ayc))))
plot(ax(axc < 1.5 * mean(sort(axc))), axc(axc < 1.5 * mean(sort(axc))))

xfit = fit(ax, axc - xfit_fakes(ax), 'gauss1');
yfit = fit(ay, ayc - yfit_fakes(ay), 'gauss1');

x = linspace(min(ax), max(ax), 1000);
counts_x_interp = interp1(ax, axc, x);
% Find the half max value.
halfMax = (min(xfit(x)) + max(xfit(x))) / 2;
% Find where the data first drops below half the max.
xl = x(find(xfit(x) >= halfMax, 1, 'first'));
xr = x(find(xfit(x) >= halfMax, 1, 'last'));
fwhmx = xr - xl;

y = linspace(min(ay), max(ay), 1000);
counts_y_interp = interp1(ay, ayc, y);
halfMax = (min(yfit(y)) + max(yfit(y))) / 2;
% Find where the data first drops below half the max.
yl = y(find(yfit(y) >= halfMax, 1, 'first'));
yr = y(find(yfit(y) >= halfMax, 1, 'last'));
fwhmy = yr - yl;

norm_cumsum_x = cumsum(counts_x_interp)/max(cumsum(counts_x_interp));
lb_x = max(x(norm_cumsum_x < 0.25));
ub_x = min(x(norm_cumsum_x > 0.75));

norm_cumsum_y = cumsum(counts_y_interp)/max(cumsum(counts_y_interp));
lb_y = max(y(norm_cumsum_y < 0.25));
ub_y = min(y(norm_cumsum_y > 0.75));

disp(['xmean = ',num2str(x(xfit(x) == max(xfit(x)))),' ; fwhmx = ', num2str(fwhmx), ' ; stddev = ', num2str(fwhmx/(2*sqrt(2*log(2))))]);
disp(['left limit = ',num2str(lb_x),' ; right limit = ',num2str(ub_x)]);
disp(['relative left limit = ', num2str(abs(lb_x - x(xfit(x) == max(xfit(x))))),' ; relative right limit = ', num2str(abs(ub_x - x(xfit(x) == max(xfit(x)))))])

disp(newline)
disp(['ymean = ',num2str(y(yfit(y) == max(yfit(y)))),' ; fwhmy = ', num2str(fwhmy), ' ; stddev = ', num2str(fwhmy/(2*sqrt(2*log(2))))]);
disp(['left limit = ',num2str(lb_y),' ; right limit = ',num2str(ub_y)]);
disp(['relative left limit = ', num2str(abs(lb_y - y(yfit(y) == max(yfit(y))))),' ; relative right limit = ', num2str(abs(ub_y - y(yfit(y) == max(yfit(y)))))])

% figure
% hold on
% plot(x, xc/max(xc));
% plot(edges(1:end-1), n/max(n));
% save '/home/christian/Documents/cern2018/code/beamParameters/2017_50GeV_2mm_beam.txt' 'data' -ascii

% figure
% hold on
% plot(x, xc)
% plot(y, yc)
% 
figure
hold on
plot(ax, xfit(ax))
plot(ay, yfit(ay))
plot(ax, axc)
plot(ay, ayc)

% save in 2018 compatible datafiles
% save '../beamParameters/angle_xweight_2017_50GeV_2mm.txt' 'axc' -ascii
% save '../beamParameters/angle_yweight_2017_50GeV_2mm.txt' 'ayc' -ascii
% save '../beamParameters/xweight_2017_50GeV_2mm.txt' 'xc' -ascii
% save '../beamParameters/yweight_2017_50GeV_2mm.txt' 'yc' -ascii
% save '../beamParameters/xdat_2017_50GeV_2mm.txt' 'x' -ascii
% save '../beamParameters/ydat_2017_50GeV_2mm.txt' 'y' -ascii
% save '../beamParameters/anglexdat_2017_50GeV_2mm.txt' 'ax' -ascii
% save '../beamParameters/angleydat_2017_50GeV_2mm.txt' 'ay' -ascii
