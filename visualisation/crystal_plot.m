clc; clear all;
datpath = '/home/christian/Documents/cern2018/simdata/';

colors =    [       0    0.4470    0.7410
                0.8500    0.3250    0.0980
                0.9290    0.6940    0.1250
                0.4940    0.1840    0.5560
                0.4660    0.6740    0.1880
                0.3010    0.7450    0.9330
                0.6350    0.0780    0.1840
            ];
xvert = [-1.124e+04, -1.14e+04, 9611, 9709, -1.124e+04];
yvert = [-6178, 4206, 4672, -5606, -6178];

runs = [46, 84];
crystal_coords = [];
sz = 0;
tstring = {'1.5mm crystal transformed', '1.0mm crystal transformed'};
for i = 1:length(runs)
    run = runs(i);
    dat = load(strcat(datpath,'crystal_image_',num2str(run),'.txt'));
    
    dat = transform(dat);
    
    
    f = figure;
    title(tstring{i});
    hold on
    plot(xvert, yvert,'k')
    plot(dat(:,1), dat(:,2),'o','markersize',3,'MarkerFaceColor',colors(1,:))
    xlabel('xpos $\left[\mu\mathrm{m}\right]$','fontsize',22,'interpreter','latex'); ylabel('ypos $\left[\mu\mathrm{m}\right]$','fontsize',22,'interpreter','latex');
    xlim([-12000, 11000])
    ylim([-5300, 5300])
    set(gca, 'FontSize', 18)
    box on
    axis equal
    grid on

end


% set(f, 'Units','centimeters','PaperUnits','centimeters', 'PaperSize',[18, 12],'PaperPosition',[0, 0, 18, 12],'Position',[0 0 18 12])
% print('../../presentation/figures/crystal_plot.svg', '-dsvg','-r600','-painters')

function coord = transform(coord)
    coord(:,1) = (2.060e6 - 1.832e6) * 1.6617e-5 + coord(:,1);
    coord(:,2) = (2.060e6 - 1.832e6) * -3.8038e-5 + coord(:,2);
end