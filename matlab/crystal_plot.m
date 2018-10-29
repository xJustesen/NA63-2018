clc; clear all; close all
datpath = '/home/christian/Dropbox/speciale/data/';

colors =    [       0    0.4470    0.7410
                0.8500    0.3250    0.0980
                0.9290    0.6940    0.1250
                0.4940    0.1840    0.5560
                0.4660    0.6740    0.1880
                0.3010    0.7450    0.9330
                0.6350    0.0780    0.1840
            ];

runs = [73, 75, 76, 77, 78];
crystal_coords = [];
sz = 0;
for i = runs
%     dat = load(strcat(datpath,'crystal_image_',num2str(i),'.txt'));
    dat = load(strcat(datpath,'crystal_image_30.txt'));
    crystal_coords = [crystal_coords; dat];
end

xvert = [-1.124e+04, -1.14e+04, 9611, 9709, -1.124e+04];
yvert = [-6178, 4206, 4672, -5606, -6178];

f = figure;
hold on
plot(xvert, yvert,'k')
plot(crystal_coords(:,1), crystal_coords(:,2),'o','markersize',3,'MarkerFaceColor',colors(1,:))
xlabel('xpos $\left[\mu\mathrm{m}\right]$','fontsize',22,'interpreter','latex'); ylabel('ypos $\left[\mu\mathrm{m}\right]$','fontsize',22,'interpreter','latex');
xlim([-12000, 11000])
ylim([-5300, 5300])
set(gca, 'FontSize', 18)

box on
axis equal
grid on
set(f, 'Units','centimeters','PaperUnits','centimeters', 'PaperSize',[18, 12],'PaperPosition',[0, 0, 18, 12],'Position',[0 0 18 12])
print('../../figures/crystal_plot.pdf', '-dpdf','-r600','-painters')

function hitdat = loaddat(fileID, formatSpec, datblocks)
    for i = 1:datblocks
        field = strcat('plane',num2str(i-1));
        data = textscan(fileID,formatSpec,'HeaderLines',1);
        hitdat.(field) = [data{1},data{2}];
    end
end