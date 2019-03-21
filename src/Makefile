CC	:= g++
CXXFLAGS = -I /home/christian/Downloads/boost_1_67_0/ -std=c++14 -Ofast -fopenmp -Wall -larmadillo 
LDLIBS = -lm -march=native -larmadillo -fopenmp $(shell /home/christian/Documents/root/bin/root-config --libs)
DATPATH := /media/christian/Elements/CERN2018NA63
LIBPATH := ./lib
CONFIGPATH := ./config

data_20GeV_1mm : main
	# 20 GeV ; Aligned ; 1mm C
	./main $(CONFIGPATH)/data/aligned_20GeV_1.0mm_dat_config.txt "/home/christian/Documents/cern2018/root data/run103_01.root"
	./main $(CONFIGPATH)/data/aligned_20GeV_1.0mm_dat_config.txt "/home/christian/Documents/cern2018/root data/run104_01.root"
	./main $(CONFIGPATH)/data/aligned_20GeV_1.0mm_dat_config.txt "/home/christian/Documents/cern2018/root data/run112_01.root"
	./main $(CONFIGPATH)/data/aligned_20GeV_1.0mm_dat_config.txt "/home/christian/Documents/cern2018/root data/run113_01.root"
	./main $(CONFIGPATH)/data/aligned_20GeV_1.0mm_dat_config.txt "/home/christian/Documents/cern2018/root data/run114_01.root"
	# 20 GeV ; Amorphous ; 1mm C
	./main $(CONFIGPATH)/data/amorph_20GeV_1.0mm_dat_config.txt "/home/christian/Documents/cern2018/root data/run109_01.root"
	./main $(CONFIGPATH)/data/amorph_20GeV_1.0mm_dat_config.txt "/home/christian/Documents/cern2018/root data/run115_01.root"
	# # 20 GeV ; Background ; 1mm C
	./main $(CONFIGPATH)/data/background_20GeV_1.0mm_dat_config.txt "/home/christian/Documents/cern2018/root data/run105_01.root"
	./main $(CONFIGPATH)/data/background_20GeV_1.0mm_dat_config.txt "/home/christian/Documents/cern2018/root data/run106_01.root"
	./main $(CONFIGPATH)/data/background_20GeV_1.0mm_dat_config.txt "/home/christian/Documents/cern2018/root data/run107_01.root"
	./main $(CONFIGPATH)/data/background_20GeV_1.0mm_dat_config.txt "/home/christian/Documents/cern2018/root data/run108_01.root"
	./main $(CONFIGPATH)/data/background_20GeV_1.0mm_dat_config.txt "/home/christian/Documents/cern2018/root data/run110_01.root"
	./main $(CONFIGPATH)/data/background_20GeV_1.0mm_dat_config.txt "/home/christian/Documents/cern2018/root data/run111_01.root"

data_40GeV_1mm : main
	# 40 GeV ; Aligned ; 1.0mm C
	./main $(CONFIGPATH)/data/aligned_40GeV_1.0mm_dat_config.txt $(DATPATH)/run71_01.root
	./main $(CONFIGPATH)/data/aligned_40GeV_1.0mm_dat_config.txt $(DATPATH)/run72_01.root
	./main $(CONFIGPATH)/data/aligned_40GeV_1.0mm_dat_config.txt $(DATPATH)/run79_01.root
	./main $(CONFIGPATH)/data/aligned_40GeV_1.0mm_dat_config.txt $(DATPATH)/run80_01.root
	./main $(CONFIGPATH)/data/aligned_40GeV_1.0mm_dat_config.txt $(DATPATH)/run81_01.root
	# 40 GeV ; Amorphous ; 1.0mm C
	./main $(CONFIGPATH)/data/amorph_40GeV_1.0mm_dat_config.txt $(DATPATH)/run73_01.root
	./main $(CONFIGPATH)/data/amorph_40GeV_1.0mm_dat_config.txt $(DATPATH)/run75_01.root
	./main $(CONFIGPATH)/data/amorph_40GeV_1.0mm_dat_config.txt $(DATPATH)/run76_01.root
	./main $(CONFIGPATH)/data/amorph_40GeV_1.0mm_dat_config.txt $(DATPATH)/run77_01.root
	./main $(CONFIGPATH)/data/amorph_40GeV_1.0mm_dat_config.txt $(DATPATH)/run78_01.root
	# 20 GeV ; Background ; 1.0mm C
	./main $(CONFIGPATH)/data/background_40GeV_1.0mm_dat_config.txt $(DATPATH)/run74_01.root

data_80GeV_1mm : main
	# 80 GeV ; Aligned ; 1.0mm C
	./main $(CONFIGPATH)/data/aligned_80GeV_1.0mm_dat_config.txt $(DATPATH)/run84_01.root
	./main $(CONFIGPATH)/data/aligned_80GeV_1.0mm_dat_config.txt $(DATPATH)/run92_01.root
	./main $(CONFIGPATH)/data/aligned_80GeV_1.0mm_dat_config.txt $(DATPATH)/run93_01.root
	./main $(CONFIGPATH)/data/aligned_80GeV_1.0mm_dat_config.txt $(DATPATH)/run94_01.root
	./main $(CONFIGPATH)/data/aligned_80GeV_1.0mm_dat_config.txt $(DATPATH)/run95_01.root
	# 80 GeV ; Amorphous ; 1.0mm C
	./main $(CONFIGPATH)/data/amorph_80GeV_1.0mm_dat_config.txt $(DATPATH)/run86_01.root
	./main $(CONFIGPATH)/data/amorph_80GeV_1.0mm_dat_config.txt $(DATPATH)/run87_01.root
	./main $(CONFIGPATH)/data/amorph_80GeV_1.0mm_dat_config.txt $(DATPATH)/run88_01.root
	./main $(CONFIGPATH)/data/amorph_80GeV_1.0mm_dat_config.txt $(DATPATH)/run89_01.root
	# 80 GeV ; Background ; 1.0mm C
	./main $(CONFIGPATH)/data/background_80GeV_1.0mm_dat_config.txt $(DATPATH)/run85_01.root
	./main $(CONFIGPATH)/data/background_80GeV_1.0mm_dat_config.txt $(DATPATH)/run90_01.root
	./main $(CONFIGPATH)/data/background_80GeV_1.0mm_dat_config.txt $(DATPATH)/run91_01.root

data_20GeV_1.5mm : main
	# 20 GeV ; Aligned ; 1.5mm C
	./main $(CONFIGPATH)/data/aligned_20GeV_1.5mm_dat_config.txt "/home/christian/Documents/cern2018/root data/run61_01.root"
	./main $(CONFIGPATH)/data/aligned_20GeV_1.5mm_dat_config.txt "/home/christian/Documents/cern2018/root data/run62_01.root"
	./main $(CONFIGPATH)/data/aligned_20GeV_1.5mm_dat_config.txt "/home/christian/Documents/cern2018/root data/run63_01.root"
	./main $(CONFIGPATH)/data/aligned_20GeV_1.5mm_dat_config.txt "/home/christian/Documents/cern2018/root data/run64_01.root"
	# 20 GeV ; Amorphous ; 1.5mm C
	./main $(CONFIGPATH)/data/amorph_20GeV_1.5mm_dat_config.txt "/home/christian/Documents/cern2018/root data/run66_01.root"
	./main $(CONFIGPATH)/data/amorph_20GeV_1.5mm_dat_config.txt "/home/christian/Documents/cern2018/root data/run67_01.root"
	./main $(CONFIGPATH)/data/amorph_20GeV_1.5mm_dat_config.txt "/home/christian/Documents/cern2018/root data/run68_01.root"
	./main $(CONFIGPATH)/data/amorph_20GeV_1.5mm_dat_config.txt "/home/christian/Documents/cern2018/root data/run69_01.root"
	# 20 GeV ; Background ; 1.5mm C
	./main $(CONFIGPATH)/data/background_20GeV_1.5mm_dat_config.txt "/home/christian/Documents/cern2018/root data/run60_01.root"
	./main $(CONFIGPATH)/data/background_20GeV_1.5mm_dat_config.txt "/home/christian/Documents/cern2018/root data/run65_01.root"

data_40GeV_1.5mm : main
	# 40 GeV ; Aligned ; 1.5mm C
	./main $(CONFIGPATH)/data/aligned_40GeV_1.5mm_dat_config.txt $(DATPATH)/run30_01.root
	./main $(CONFIGPATH)/data/aligned_40GeV_1.5mm_dat_config.txt $(DATPATH)/run35_01.root
	./main $(CONFIGPATH)/data/aligned_40GeV_1.5mm_dat_config.txt $(DATPATH)/run36_01.root
	./main $(CONFIGPATH)/data/aligned_40GeV_1.5mm_dat_config.txt $(DATPATH)/run37_01.root
	./main $(CONFIGPATH)/data/aligned_40GeV_1.5mm_dat_config.txt $(DATPATH)/run38_01.root
	# 40 GeV ; Amorphous ; 1.5mm C
	./main $(CONFIGPATH)/data/amorph_40GeV_1.5mm_dat_config.txt $(DATPATH)/run86_01.root
	./main $(CONFIGPATH)/data/amorph_40GeV_1.5mm_dat_config.txt $(DATPATH)/run87_01.root
	./main $(CONFIGPATH)/data/amorph_40GeV_1.5mm_dat_config.txt $(DATPATH)/run88_01.root
	./main $(CONFIGPATH)/data/amorph_40GeV_1.5mm_dat_config.txt $(DATPATH)/run89_01.root
	# 40 GeV ; Background ; 1.5mm C
	./main $(CONFIGPATH)/data/background_40GeV_1.5mm_dat_config.txt $(DATPATH)/run31_01.root

data_80GeV_1.5mm : main
	# 80 GeV ; Aligned ; 1.5mm C
	./main $(CONFIGPATH)/data/aligned_80GeV_1.5mm_dat_config.txt $(DATPATH)/run46_01.root
	./main $(CONFIGPATH)/data/aligned_80GeV_1.5mm_dat_config.txt $(DATPATH)/run62_01.root
	# 80 GeV ; Amorph ; 1.5mm C
	./main $(CONFIGPATH)/data/amorph_80GeV_1.5mm_dat_config.txt $(DATPATH)/run49_01.root
	./main $(CONFIGPATH)/data/amorph_80GeV_1.5mm_dat_config.txt $(DATPATH)/run50_01.root
	./main $(CONFIGPATH)/data/amorph_80GeV_1.5mm_dat_config.txt $(DATPATH)/run51_01.root
	./main $(CONFIGPATH)/data/amorph_80GeV_1.5mm_dat_config.txt $(DATPATH)/run52_01.root
	./main $(CONFIGPATH)/data/amorph_80GeV_1.5mm_dat_config.txt $(DATPATH)/run53_01.root
	./main $(CONFIGPATH)/data/amorph_80GeV_1.5mm_dat_config.txt $(DATPATH)/run56_01.root
	./main $(CONFIGPATH)/data/amorph_80GeV_1.5mm_dat_config.txt $(DATPATH)/run57_01.root
	# 80 GeV ; Amorph ; 1.5mm C
	./main $(CONFIGPATH)/data/background_80GeV_1.5mm_dat_config.txt $(DATPATH)/run48_01.root
	./main $(CONFIGPATH)/data/background_80GeV_1.5mm_dat_config.txt $(DATPATH)/run54_01.root
	./main $(CONFIGPATH)/data/background_80GeV_1.5mm_dat_config.txt $(DATPATH)/run55_01.root

sim_aligned_20GeV_1mm : main
	./main $(CONFIGPATH)/sim/aligned_20GeV_1.0mm_sim_config_stochastic.txt

sim_amorph_40GeV_1.5mm : main
	./main $(CONFIGPATH)/sim/amorph_40GeV_1.5mm_sim_config.txt

sim_monochrome_80GeV_1.5mm : main
	./main $(CONFIGPATH)/sim/monochrome_80GeV_1.5mm_sim_config.txt

main : main.cc $(LIBPATH)/simulator.cc $(LIBPATH)/analyser.cc $(LIBPATH)/pre_processor.cc $(LIBPATH)/run_simulation.cc $(LIBPATH)/run_data_analysis.cc $(LIBPATH)/auxillary_functions.cc

test : test.cc

clean:
	rm -f main test
