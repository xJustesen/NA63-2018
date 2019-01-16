close all;

datpath = '/home/christian/Dropbox/speciale/data/';
colors =    [       0    0.4470    0.7410
                0.8500    0.3250    0.0980
                0.9290    0.6940    0.1250
                0.4940    0.1840    0.5560
                0.4660    0.6740    0.1880
                0.3010    0.7450    0.9330
                0.6350    0.0780    0.1840
            ];
for i = 2:5
    dat = load(strcat(datpath, "plane", num2str(i),"_105_dist.txt"));
    dat2 = load(strcat(datpath, "plane", num2str(i),"_103_dist.txt"));
    dat3 = load(strcat(datpath, "plane", num2str(i),"_4_dist.txt"));

    [counts, edges] = histcounts(dat(dat < 100),1000);
    [COUNTS, EDGES] = histcounts(dat2(dat2 < 100),1000);
    [C, E] = histcounts(dat2(dat3 < 100),1000);
    
    [counts2, edges2] = histcounts(dat,250);
    [COUNTS2, EDGES2] = histcounts(dat2,250);
    [C2, E2] = histcounts(dat3,250);


    f = figure(1);
    hold on
    plot(edges(2:end),counts/324000 - C/130042,'-','linewidth',2.5,'color',colors(i+1,:));
    plot(EDGES(2:end),COUNTS/876217 - C/130042,'x','linewidth',2.5,'color',colors(i+1,:));
    grid on
    box on
    
    f = figure(2);
    hold on
    plot(edges2(2:end),counts2/324000 - C2/130042,'-','linewidth',2.5,'color',colors(i+1,:));
    plot(EDGES2(2:end),COUNTS2/876217 - C2/130042,'x','linewidth',2.5,'color',colors(i+1,:));
    grid on
    box on

%     axes('position',[.55 .55 .3 .3])
%     edges = edges(2:end);
%     filter = edges < 100 & edges > 0;
%     
%     box on
%     grid on
%     hold on
%     plot(edges(filter), counts(filter),'-','linewidth',2.5,'color',colors(1,:))
%     axis tight
%     set(gca, 'FontSize', 18)


end

xlabel('d $[\mu\mathrm{m}]$','fontsize',22,'interpreter','latex');
ylabel('Counts','fontsize',22,'interpreter','latex');
set(gca, 'fontsize',18)