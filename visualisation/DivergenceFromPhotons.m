clear all; clc; close all
bins = linspace(-4e-3, 4e-3, 150);
%% 20 GeV 1mm
data = load('../beamParameters/crystal_image_amorph_dat_20GeV_1mm.txt');
[y, x] = hist(data(:,1), bins);
fit_object = fit(x.', y.', 'gauss1');
figure
hold on
bar(x, y)
plot(fit_object)
title('20GeV 1mm ; x')
sigx1mm20GeV = fit_object.c1/sqrt(2);
mux1mm20GeV = fit_object.b1;
[y, x] = hist(data(:,2), bins);
fit_object = fit(x.', y.', 'gauss1');
figure
hold on
bar(x, y)
plot(fit_object)
title('20GeV 1mm; y')
sigy1mm20GeV = fit_object.c1/sqrt(2);
muy1mm20GeV = fit_object.b1;

%% 20 GeV 1.5mm
data = load('../beamParameters/crystal_image_amorph_dat_20GeV_1.5mm.txt');
[y, x] = hist(data(:,1), bins);
fit_object = fit(x.', y.', 'gauss1');
figure
hold on
bar(x, y)
plot(fit_object)
title('20GeV 1.5mm ; x')
sigx1_5mm20GeV = fit_object.c1/sqrt(2);
mux1_5mm20GeV = fit_object.b1;
[y, x] = hist(data(:,2), bins);
fit_object = fit(x.', y.', 'gauss1');
figure
hold on
bar(x, y)
plot(fit_object)
title('20GeV 1.5mm; y')
sigy1_5mm20GeV = fit_object.c1/sqrt(2);
muy1_5mm20GeV = fit_object.b1;

%% 40 GeV 1mm
data = load('../beamParameters/crystal_image_amorph_data_40GeV_1mm.txt');
[y, x] = hist(data(:,1), bins);
fit_object = fit(x.', y.', 'gauss1');
figure
hold on
bar(x, y)
plot(fit_object)
title('40GeV 1mm ; x')
sigx1mm40GeV = fit_object.c1/sqrt(2);
mux1mm40GeV = fit_object.b1;
[y, x] = hist(data(:,2), bins);
fit_object = fit(x.', y.', 'gauss1');
figure
hold on
bar(x, y)
plot(fit_object)
title('40GeV 1mm ; y')
sigy1mm40GeV = fit_object.c1/sqrt(2);
muy1mm40GeV = fit_object.b1;

%% 40 GeV 1.5mm
data = load('../beamParameters/crystal_image_amorph_data_40GeV_1.5mm.txt');
[y, x] = hist(data(:,1), bins);
fit_object = fit(x.', y.', 'gauss1');
figure
hold on
bar(x, y)
plot(fit_object)
title('40GeV 1.5mm ; x')
sigx1_5mm40GeV = fit_object.c1/sqrt(2);
mux1_5mm40GeV = fit_object.b1;
[y, x] = hist(data(:,2), bins);
fit_object = fit(x.', y.', 'gauss1');
figure
hold on
bar(x, y)
plot(fit_object)
title('40GeV 1.5mm ; y')
sigy1_5mm40GeV = fit_object.c1/sqrt(2);
muy1_5mm40GeV = fit_object.b1;

%% 80 GeV 1mm
data = load('../beamParameters/crystal_image_amorph_data_80GeV_1mm.txt');
[y, x] = hist(data(:,1), bins);
fit_object = fit(x.', y.', 'gauss1','StartpoinT',[max(y), median(x), 200e-6]);
figure
hold on
bar(x, y)
plot(fit_object)
title('80GeV 1mm ; x')
sigx1mm80GeV = fit_object.c1/sqrt(2);
mux1mm80GeV = fit_object.b1;

[y, x] = hist(data(:,2), bins);
fit_object = fit(x.', y.', 'gauss1');
figure
hold on
bar(x, y)
plot(fit_object)
title('80GeV 1mm ; y')
sigy1mm80GeV = fit_object.c1/sqrt(2);
muy1mm80GeV = fit_object.b1;

%% 80 GeV 1.5mm
data = load('../beamParameters/crystal_image_amorph_data_80GeV_1.5mm.txt');
[y, x] = hist(data(:,1), bins);
fit_object = fit(x.', y.', 'gauss1','StartpoinT',[max(y), median(x), 200e-6]);
figure
hold on
bar(x, y)
plot(fit_object)
title('80GeV 1.5mm ; x')
sigx1_5mm80GeV = fit_object.c1/sqrt(2);
mux1_5mm80GeV = fit_object.b1;

[y, x] = hist(data(:,2), bins);
fit_object = fit(x.', y.', 'gauss1');
figure
hold on
bar(x, y)
plot(fit_object)
title('80GeV 1.5mm ; y')
sigy1_5mm80GeV = fit_object.c1/sqrt(2);
muy1_5mm80GeV = fit_object.b1;

%%
limit = @(sig, mu) mu + [erfinv(-0.95 * sqrt(2) * sig); erfinv(0.95 * sqrt(2) * sig)];

xlim20GeV1mm = limit(sigx1mm20GeV, mux1mm20GeV);
ylim20GeV1mm = limit(sigy1mm20GeV, muy1mm20GeV);
xlim40GeV1mm = limit(sigx1mm40GeV, mux1mm40GeV);
ylim40GeV1mm = limit(sigy1mm40GeV, muy1mm40GeV);
xlim80GeV1mm = limit(sigx1mm80GeV, mux1mm80GeV);
ylim80GeV1mm = limit(sigy1mm80GeV, muy1mm80GeV);
xlim20GeV1_5mm = limit(sigx1_5mm20GeV, mux1_5mm20GeV);
ylim20GeV1_5mm = limit(sigy1_5mm20GeV, muy1_5mm20GeV);
xlim40GeV1_5mm = limit(sigx1_5mm40GeV, mux1_5mm40GeV);
ylim40GeV1_5mm = limit(sigy1_5mm40GeV, muy1_5mm40GeV);
xlim80GeV1_5mm = limit(sigx1_5mm80GeV, mux1_5mm80GeV);
ylim80GeV1_5mm = limit(sigy1_5mm80GeV, muy1_5mm80GeV);


fprintf('\n\txmu\t\tymu\t\txmu\t\tymu\n');
fprintf('\t1.0mm\t\t1.0mm\t\t1.5mm\t\t1.5mm\n');
fprintf('20GeV\t%e\t%e\t%e\t%e\n',mux1mm20GeV,muy1mm20GeV,mux1_5mm20GeV,muy1_5mm20GeV);
fprintf('40GeV\t%e\t%e\t%e\t%e\n',mux1mm40GeV,muy1mm40GeV,mux1_5mm40GeV,muy1_5mm40GeV);
fprintf('80GeV\t%e\t%e\t%e\t%e\n',mux1mm80GeV,muy1mm80GeV,mux1_5mm80GeV,muy1_5mm80GeV);

fprintf('\n\txsig\t\tysig\t\txsig\t\tysig\n');
fprintf('\t1.0mm\t\t1.0mm\t\t1.5mm\t\t1.5mm\n');
fprintf('20GeV\t%e\t%e\t%e\t%e\n',sigx1mm20GeV,sigy1mm20GeV,sigx1_5mm20GeV,sigy1_5mm20GeV);
fprintf('40GeV\t%e\t%e\t%e\t%e\n',sigx1mm40GeV,sigy1mm40GeV,sigx1_5mm40GeV,sigy1_5mm40GeV);
fprintf('80GeV\t%e\t%e\t%e\t%e\n',sigx1mm80GeV,sigy1mm80GeV,sigx1_5mm80GeV,sigy1_5mm80GeV);

fprintf('\nlimits:\n\t20GeV\t\t\t\t40GeV\t\t\t\t80GeV\n');
fprintf('\t1.0mm\t\t1.5mm\t\t1.0mm\t\t1.5mm\t\t1.0mm\t\t1.5mm\n');
fprintf('\t%e\t%e\t%e\t%e\t%e\t%e\n',xlim20GeV1mm(1),xlim20GeV1_5mm(1),xlim40GeV1mm(1),xlim40GeV1_5mm(1),xlim80GeV1mm(1),xlim80GeV1_5mm(1));
fprintf('\t%e\t%e\t%e\t%e\t%e\t%e\n',xlim20GeV1mm(2),xlim20GeV1_5mm(2),xlim40GeV1mm(2),xlim40GeV1_5mm(2),xlim80GeV1mm(2),xlim80GeV1_5mm(2));
fprintf('\t%e\t%e\t%e\t%e\t%e\t%e\n',ylim20GeV1mm(1),ylim20GeV1_5mm(1),ylim40GeV1mm(1),ylim40GeV1_5mm(1),ylim80GeV1mm(1),ylim80GeV1_5mm(1));
fprintf('\t%e\t%e\t%e\t%e\t%e\t%e\n',ylim20GeV1mm(2),ylim20GeV1_5mm(2),ylim40GeV1mm(2),ylim40GeV1_5mm(2),ylim80GeV1mm(2),ylim80GeV1_5mm(2));




