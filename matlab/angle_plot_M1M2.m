clc; close all;
datpath = '/home/christian/Dropbox/speciale/data/';

runs_dat = [49, 50, 51, 52, 53, 56, 57];
angles = linspace(-2e-3, 2e-3, 200);
angx_tot = [];
angy_tot = [];

for i = runs_dat
    % Load data
    field = strcat('run_',num2str(i));
    filepath = strcat(datpath,['angles_M1M2_',num2str(i),'.txt']);
    ang = load(filepath);
    angx = ang(:,1);
    angy = ang(:,2);
    
    % Save angles in matrix
    angx_tot((1 + numel(angx_tot)):numel(angx)+(numel(angx_tot))) = angx;
    angy_tot((1 + numel(angy_tot)):numel(angy)+(numel(angy_tot))) = angy;

    % Make histogram for each run
    [counts_x, ~] = hist(angx(angx >= min(angles) & angx <= max(angles)), angles);
    [counts_y, ~] = hist(angy(angy >= min(angles) & angy <= max(angles)), angles);
       
    % Save histogram for each run in struct
    Counts_x.(field) = counts_x;
    Counts_y.(field) = counts_y;
end

% Make histogram for combined runs
[counts_x_tot, ~] = hist(angx_tot(angx_tot > angles(1) & angx_tot < angles(end)), angles);
[counts_y_tot, ~] = hist(angy_tot(angy_tot > angles(1) & angy_tot < angles(end)), angles);

counts_x_tot = counts_x_tot.';
counts_y_tot = counts_y_tot.';
angles = angles.';
save('../cpp/angles.txt','angles','-ascii');
save('../cpp/angles_xweight_80GeV_beam_params_1.5mm.txt','counts_x_tot','-ascii');
save('../cpp/angles_yweight_80GeV_beam_params_1.5mm.txt','counts_y_tot','-ascii');


xfit = fit(angles, counts_x_tot, 'gauss1')
yfit = fit(angles, counts_y_tot, 'gauss1')
    
x = linspace(min(angles), max(angles), 1000);
% Find the half max value.
halfMax = (min(xfit(x)) + max(xfit(x))) / 2;
% Find where the data first drops below half the max.
xl = x(find(xfit(x) >= halfMax, 1, 'first'));
xr = x(find(xfit(x) >= halfMax, 1, 'last'));
fwhmx = xr - xl;
    
y = linspace(min(angles), max(angles), 1000);
halfMax = (min(yfit(y)) + max(yfit(y))) / 2;
% Find where the data first drops below half the max.
yl = y(find(yfit(y) >= halfMax, 1, 'first'));
yr = y(find(yfit(y) >= halfMax, 1, 'last'));
fwhmy = yr - yl;
disp('')
disp(['xmean = ',num2str(x(xfit(x) == max(xfit(x)))),' ; fwhmx = ', num2str(fwhmx), ' ; stddev = ', num2str(fwhmx/(2*sqrt(2*log(2))))]);
disp(['ymean = ',num2str(y(yfit(y) == max(yfit(y)))),' ; fwhmy = ', num2str(fwhmy), ' ; stddev = ', num2str(fwhmy/(2*sqrt(2*log(2))))]);
    

f = figure(1)
hold on
box on
grid on
stairs(angles, counts_x_tot);
stairs(angles, counts_y_tot);
xlabel('ang (rad)'); ylabel('counts');
legend('x','y')
set(f, 'Units','centimeters','PaperUnits','centimeters', 'PaperSize',[18, 12],'PaperPosition',[0, 0, 18, 12],'Position',[0 0 18 12])
% print('../../figures/20GeV_angular_distro', '-dpdf','-r600','-painters')


