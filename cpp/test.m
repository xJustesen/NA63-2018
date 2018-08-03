close all; clear all;

% intensity = load('sum_angles1mm40GeVelec.txt');
% energy = linspace(0, 40, length(intensity));
% intensity_interp = load('interp_intensity.txt');
% energy_interp = linspace(0, 40, length(intensity_interp));
% 
% figure
% hold on
% plot(energy_interp, intensity_interp, 'o')
% plot(energy, intensity, 'x-','linewidth',2)
% 
% filepath = 'distro_draws.txt';
% formatSpec = '%f %f';
% fileID = fopen(filepath,'r');
% datblocks = 212;
% 
% intensity = loaddat(fileID, formatSpec, datblocks);
% angles = linspace(0, 100,length(intensity.distro0));
% 
% dat = load("sum_initials1mm40GeVelec.txt");
% dat(1:5,:) = [];
% angles2 = linspace(0, 100, size(dat, 2));
% 
% figure
% hold on
% plot(angles, intensity.distro0,'-')
% plot(angles2, dat(1,:),'.')

dat = load('sum_initials1mm40GeVelec.txt');
theta = find_angles(dat);
energy = linspace(0, 40, 100);

intensity = meshgrid(dat(10, 1:5:end));
[tx, ty] = meshgrid(theta.x(1:5:end), theta.y(1:5:end));

surf(tx, ty, intensity)
xlabel('tx');ylabel('ty');zlabel(['Intensity for E = ', num2str(energy(5)), 'GeV'])

function theta = find_angles(dat)
    V = dat(1:3, :);
%     theta.x = acos(V(1, :));
theta.x = sqrt(2 * (1 - V(1,:)));
theta.y = sqrt(2 * (1 - V(2,:)));
    theta.y = acos(V(2, :));
end

function dat = loaddat(fileID, formatSpec, datblocks)
    for i = 1:datblocks
        field = strcat('distro',num2str(i-1));
        datastruct = textscan(fileID,formatSpec,'HeaderLines',1);
        dat.(field) = [datastruct{1}];
    end
end