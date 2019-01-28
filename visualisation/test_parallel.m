DATAPATH = "~/Dropbox/Cern2018Experiment/spectre/";
FILENAME = strcat(DATAPATH,"energy_parallel_test.txt");

DATA = load(FILENAME) * 6.2415091E9;;
ENERGY = linspace(0, 80, 80);

SPECTRUM = hist(DATA(DATA > ENERGY(1) & DATA < ENERGY(end)), ENERGY);

plot(ENERGY, SPECTRUM.*ENERGY)