clear all; close all; clear all;
datpath = '/home/christian/Documents/cern2018/simdata/';

runs_dat = 44;
sim = 1:1;
angles = linspace(-2e-3, 2e-3, 200);
angx_tot = [];
angy_tot = [];
posx_tot = [];
posy_tot = [];

for i = runs_dat
    % Load data
    field = strcat('run_',num2str(i));
    filepath = strcat(datpath,'beam_divergence_',num2str(i),'.txt');
    dat = load(filepath);
    
    angx = dat(:,1);
    angy = dat(:,2);
    
    % Save angles and positions in matrix
    angx_tot = [angx_tot; angx];
    angy_tot = [angy_tot; angy];
end

% Make histogram for combined runs
[ang_counts_x_tot, ~] = hist(angx_tot(angx_tot > angles(1) & angx_tot < angles(end)), angles);
[ang_counts_y_tot, ~] = hist(angy_tot(angy_tot > angles(1) & angy_tot < angles(end)), angles);
ang_counts_x_tot = ang_counts_x_tot.';
ang_counts_y_tot = ang_counts_y_tot.';

angles = angles.';

yfit_fakes = fit(angles(ang_counts_y_tot < 1.5 * median(sort(ang_counts_y_tot))), ang_counts_y_tot(ang_counts_y_tot < 1.5 * median(sort(ang_counts_y_tot))),'poly2');
xfit_fakes = fit(angles(ang_counts_x_tot < 1.5 * median(sort(ang_counts_x_tot))), ang_counts_x_tot(ang_counts_x_tot < 1.5 * median(sort(ang_counts_x_tot))),'poly2');

input_ang   = load('/home/christian/Dropbox/Cern2018Experiment/kode/beam_direction_4mm.txt'); 

xang = input_ang(:,3);
yang = input_ang(:,4);

xw = (ang_counts_x_tot - xfit_fakes(angles))/max((ang_counts_x_tot - xfit_fakes(angles)));
yw = (ang_counts_y_tot - yfit_fakes(angles))/max((ang_counts_x_tot - xfit_fakes(angles)));

xw(xw < 1e-03) = 0;
yw(yw < 1e-03) = 0;


f = figure;
title('x')
hold on
plot(angles, xw,'-');
plot(xang * 1e-6, input_ang(:,1)/max(input_ang(:,1)))

f = figure
title('y')
hold on
plot(angles, yw,'-');
plot(yang * 1e-6, input_ang(:,2)/max(input_ang(:,2)))
