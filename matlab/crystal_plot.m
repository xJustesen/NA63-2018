clc; clear all; close all
datpath = '/home/christian/Dropbox/speciale/data/';
crystal_coords1 = load(strcat(datpath,'crystal_image_71.txt'));
crystal_coords2 = load(strcat(datpath,'crystal_image_72.txt'));
crystal_coords3 = load(strcat(datpath,'crystal_image_79.txt'));
crystal_coords4 = load(strcat(datpath,'crystal_image_80.txt'));
crystal_coords5 = load(strcat(datpath,'crystal_image_81.txt'));

crystal_coords = [crystal_coords1; crystal_coords2; crystal_coords3; crystal_coords4; crystal_coords5];

crystal_coords = load(strcat(datpath,'crystal_image_80.txt'));

figure
hold on
plot(crystal_coords(:,1), crystal_coords(:,2),'b.','markersize',.1)
xlabel('xpos (\mum)'); ylabel('ypos (\mum)');
xlim([-10000, 10000])
ylim([-5000, 5000])
box on
axis equal
grid on
