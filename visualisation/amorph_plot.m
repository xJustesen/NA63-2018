clear all; clc;
global datpath counts;
datpath = "~/Dropbox/Cern2018Experiment/spectre/";

%% 40 GeV
E40 = linspace(0, 40, 40);
N_40GeV_1mm = [1290988 1361162 1447462 715126 1456319];
runs_40GeV_1mm = [73 75:78];
makeHist(runs_40GeV_1mm, E40, N_40GeV_1mm)

N_40GeV_1mm_bg = 3029506;
bg_run_40GeV_1mm = 74;
counts.bg40 = spectrum('energy_', '.txt', bg_run_40GeV_1mm, E40) / N_40GeV_1mm_bg;
amorph_40GeV_tot = spectrum('energy_', '.txt', runs_40GeV_1mm, E40);
amorph_40GeV_tot = amorph_40GeV_tot / sum(N_40GeV_1mm);

amorph_40GeV_tot_alt = spectrum('energy_', '.txt', runs_40GeV_1mm(1:end-1), E40);
amorph_40GeV_tot_alt = amorph_40GeV_tot_alt / sum(N_40GeV_1mm(1:end-1));


%% 80 GeV
E80 = linspace(0, 80, 60);
N_80GeV_1mm = [1719461 1172521 1210722 538281];
runs_80GeV_1mm = 86:89;
makeHist(runs_80GeV_1mm, E80, N_80GeV_1mm)

N_80GeV_1mm_bg = [1911215 1892132 1219189];
bg_run_80GeV_1mm = [85 90 91];
counts.bg80 = spectrum('energy_', '.txt', bg_run_80GeV_1mm, E80) / sum(N_80GeV_1mm_bg);
amorph_80GeV_tot = spectrum('energy_', '.txt', runs_80GeV_1mm, E80);
amorph_80GeV_tot = amorph_80GeV_tot / sum(N_80GeV_1mm);

%% 20 GeV
E20 = linspace(0, 20, 20);
N_20GeV_1mm = [2764454 82089];
runs_20GeV_1mm = [109,115];
makeHist(runs_20GeV_1mm, E20, N_20GeV_1mm)

N_20GeV_1mm_bg = [324000 590172 624602 734446 1415716 1224254];
bg_run_20GeV_1mm = [105:108, 110, 111];
counts.bg20 = spectrum('energy_', '.txt', bg_run_20GeV_1mm, E20) / sum(N_20GeV_1mm_bg);
amorph_20GeV_tot = spectrum('energy_', '.txt', runs_20GeV_1mm, E20);
amorph_20GeV_tot = amorph_20GeV_tot / sum(N_20GeV_1mm);

%% plot

f = figure;
subplot(1,2,1)
hold on
plot(E40, counts.run73)
plot(E40, counts.run75)
plot(E40, counts.run76)
plot(E40, counts.run77)
plot(E40, counts.run78)
plot(E40, amorph_40GeV_tot,'k-','linewidth',1.5)
plot(E40, counts.bg40,'k-.','linewidth',1.5)
plot(E40, (amorph_40GeV_tot-counts.bg40),'k-.','linewidth',1.5)
legend('73','75','76','77','78','tot','background','tot - bg')
grid on; box on;

subplot(1,2,2)
hold on
plot(E40, E40 .* counts.run73)
plot(E40, E40 .* counts.run75)
plot(E40, E40 .* counts.run76)
plot(E40, E40 .* counts.run77)
plot(E40, E40 .* counts.run78)
plot(E40, E40 .* amorph_40GeV_tot,'k-','linewidth',1.5)
plot(E40, E40 .* counts.bg40,'k--','linewidth',1.5)
plot(E40, E40 .* (amorph_40GeV_tot-counts.bg40),'k-.','linewidth',1.5)
legend('73','75','76','77','78','tot','background','tot - bg')
grid on; box on;

set(f, 'Units','centimeters','PaperUnits','centimeters', 'PaperSize',[36, 18],'PaperPosition',[0, 0, 36, 18],'Position',[0 0 36, 18])

f = figure;
subplot(1,2,1)
hold on
plot(E80, counts.run86)
plot(E80, counts.run87)
plot(E80, counts.run88)
plot(E80, counts.run89)
plot(E80, amorph_80GeV_tot,'k-','linewidth',1.5)
plot(E80, counts.bg80,'k--','linewidth',1.5)
plot(E80, (amorph_80GeV_tot-counts.bg80),'k-.','linewidth',1.5)
legend('86','87','88','89','tot','background','tot - bg')
grid on; box on;

subplot(1,2,2)
hold on
plot(E80, E80.*counts.run86)
plot(E80, E80.*counts.run87)
plot(E80, E80.*counts.run88)
plot(E80, E80.*counts.run89)
plot(E80, E80.*amorph_80GeV_tot,'k-','linewidth',1.5)
plot(E80, E80.*counts.bg80,'k--','linewidth',1.5)
plot(E80, E80.*(amorph_80GeV_tot-counts.bg80),'k-.','linewidth',1.5)
legend('86','87','88','89','tot','background','tot - bg')
grid on; box on;

set(f, 'Units','centimeters','PaperUnits','centimeters', 'PaperSize',[36, 18],'PaperPosition',[0, 0, 36, 18],'Position',[0 0 36, 18])


f = figure;
subplot(1,2,1)
hold on
plot(E20, counts.run109)
plot(E20, counts.run115)
plot(E20, amorph_20GeV_tot,'k-','linewidth',1.5)
plot(E20, counts.bg20,'k--','linewidth',1.5)
plot(E20, (amorph_20GeV_tot-counts.bg20),'k-.','linewidth',1.5)
legend('109','115','tot','background','tot - bg')
grid on; box on;

subplot(1,2,2)
hold on
plot(E20, E20.*counts.run109)
plot(E20, E20.*counts.run115)
plot(E20, E20.*amorph_20GeV_tot,'k-','linewidth',1.5)
plot(E20, E20.*counts.bg20,'k--','linewidth',1.5)
plot(E20, E20.*(amorph_20GeV_tot-counts.bg20),'k-.','linewidth',1.5)
legend('109','115','tot','background','tot - bg')
set(f, 'Units','centimeters','PaperUnits','centimeters', 'PaperSize',[36, 18],'PaperPosition',[0, 0, 36, 18],'Position',[0 0 36, 18])
grid on; box on;

%% func.
function s = spectrum(filePrefix, fileSuffix, fileNum, E)
global datpath;
s = [];

if isempty(fileNum)
    file = strcat(datpath, filePrefix,fileSuffix);
    s = load(file) * 6.2415091E9;
else
    for i = fileNum
        file = strcat(datpath, filePrefix,num2str(i),fileSuffix);
        nrg = load(file) * 6.2415091E9;
        s = [s; nrg];
    end
end
    s = hist(s(s > min(E) & s < max(E)), E);
end

function makeHist(runs, E, N)
global counts;
for i = 1:length(runs)
   run = runs(i);
   name = strcat('run',num2str(run));
   counts.(name) = spectrum('energy_', '.txt', run, E) / N(i);
end
end