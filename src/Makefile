CC	:= g++
CXXFLAGS = -I -std=c++14 -Ofast -fopenmp -Wall -larmadillo -I/home/christian/Root/include
LDLIBS = -lm -march=native -larmadillo -fopenmp $(shell /home/christian/Documents/root/bin/root-config --libs)
DATPATH := /media/christian/Elements/CERN2018NA63
LIBPATH := ./lib
CONFIGPATH := ./config

data_20GeV_1mm : main
	# 20 GeV ; Aligned ; 1mm C
	./main $(CONFIGPATH)/aligned_20GeV_1.0mm_dat_config.txt "/home/christian/Documents/cern2018/root data/run103_01.root"
	./main $(CONFIGPATH)/aligned_20GeV_1.0mm_dat_config.txt "/home/christian/Documents/cern2018/root data/run104_01.root"
	./main $(CONFIGPATH)/aligned_20GeV_1.0mm_dat_config.txt "/home/christian/Documents/cern2018/root data/run112_01.root"
	./main $(CONFIGPATH)/aligned_20GeV_1.0mm_dat_config.txt "/home/christian/Documents/cern2018/root data/run113_01.root"
	./main $(CONFIGPATH)/aligned_20GeV_1.0mm_dat_config.txt "/home/christian/Documents/cern2018/root data/run114_01.root"

data_20GeV_1.5mm : main
	# 20 GeV ; Aligned ; 1.5mm C
	./main $(CONFIGPATH)/aligned_20GeV_1.5mm_dat_config.txt "/home/christian/Documents/cern2018/root data/run61_01.root"
	./main $(CONFIGPATH)/aligned_20GeV_1.5mm_dat_config.txt "/home/christian/Documents/cern2018/root data/run62_01.root"
	./main $(CONFIGPATH)/aligned_20GeV_1.5mm_dat_config.txt "/home/christian/Documents/cern2018/root data/run63_01.root"
	./main $(CONFIGPATH)/aligned_20GeV_1.5mm_dat_config.txt "/home/christian/Documents/cern2018/root data/run64_01.root"

data_80GeV_1.5mm : main
	# 80 GeV ; Aligned ; 1.5mm C
	./main $(CONFIGPATH)/aligned_80GeV_1.5mm_dat_config.txt $(DATPATH)/run46_01.root
	./main $(CONFIGPATH)/aligned_80GeV_1.5mm_dat_config.txt $(DATPATH)/run62_01.root


tester : test
	./test $(CONFIGPATH)/amorph_40GeV_1.0mm_sim_config.txt

main : main.cc $(LIBPATH)/simulator.cc $(LIBPATH)/analyser.cc $(LIBPATH)/pre_processor.cc $(LIBPATH)/run_simulation.cc $(LIBPATH)/run_data_analysis.cc $(LIBPATH)/auxillary_functions.cc

test : test.cc $(LIBPATH)/simulator.cc $(LIBPATH)/analyser.cc $(LIBPATH)/pre_processor.cc $(LIBPATH)/run_simulation.cc $(LIBPATH)/run_data_analysis.cc $(LIBPATH)/auxillary_functions.cc

clean:
	rm -f main test
