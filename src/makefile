CC	:= g++
CXXFLAGS += -std=c++11 -Ofast -fopenmp -Wall -larmadillo -I/home/christian/Root/include
LDLIBS += -lm -march=native -larmadillo -fopenmp $(shell /home/christian/Documents/root/bin/root-config --libs)
DATPATH := /home/christian/Documents/CERN2018NA63
PARAMPATH := /home/christian/Dropbox/speciale/code/beamParameters

all : outputdata

outputdata : main
	# # 20 GeV ; Amorphous ; 1mm C
	# ./main $(DATPATH)/run109_01.root 109
	# ./main $(DATPATH)/run115_01.root 115
	# #
	# # # 20 GeV ; Background ; 1mm C
	# ./main $(DATPATH)/run105_01.root 105
	# ./main $(DATPATH)/run106_01.root 106
	# ./main $(DATPATH)/run107_01.root 107
	# ./main $(DATPATH)/run108_01.root 108
	# ./main $(DATPATH)/run110_01.root 110
	# ./main $(DATPATH)/run111_01.root 111
	#
	# # 20 GeV ; Amorphous ; 1.5mm C
	# ./main $(DATPATH)/run66_01.root 66
	# ./main $(DATPATH)/run67_01.root 67
	# ./main $(DATPATH)/run68_01.root 68
	# ./main $(DATPATH)/run69_01.root 69
	#
	# # 20 GeV ; Background ; 1.5mm C
	# ./main $(DATPATH)/run60_01.root 60
	# ./main $(DATPATH)/run65_01.root 65
	#
	# # 40 GeV ; Amorphous ; 1mm C
	# ./main $(DATPATH)/run73_01.root 73
	# ./main $(DATPATH)/run75_01.root 75
	# ./main $(DATPATH)/run76_01.root 76
	# ./main $(DATPATH)/run77_01.root 77
	# ./main $(DATPATH)/run78_01.root 78

	# 40 GeV ; Aligned ; 1mm C
	# ./main $(DATPATH)/run71_01.root 71
	# ./main $(DATPATH)/run72_01.root 72
	# ./main $(DATPATH)/run79_01.root 79
	# ./main $(DATPATH)/run80_01.root 80
	# ./main $(DATPATH)/run81_01.root 81

	# 40 GeV ; Background ; 1mm C
	# ./main $(DATPATH)/run74_01.root 74
	#
	# # 40 GeV ; Amorphous ; 1.5mm C
	# ./main $(DATPATH)/run32_01.root 32
	# ./main $(DATPATH)/run33_01.root 33
	# ./main $(DATPATH)/run34_01.root 34
	# ./main $(DATPATH)/run39_01.root 39
	# ./main $(DATPATH)/run40_01.root 40
	# ./main $(DATPATH)/run41_01.root 41
	# ./main $(DATPATH)/run43_01.root 43
	#
	# # # 40 GeV ; Background ; 1.5mm C
	# ./main $(DATPATH)/run31_01.root 31
	#
	# # 80 GeV ; Amorphous ; 1mm C
	# ./main $(DATPATH)/run86_01.root 86
	# ./main $(DATPATH)/run87_01.root 87
	# ./main $(DATPATH)/run88_01.root 88
	# ./main $(DATPATH)/run89_01.root 89
	#
	# # 80 GeV ; Background ; 1mm C
	# ./main $(DATPATH)/run85_01.root 85
	# ./main $(DATPATH)/run90_01.root 90
	# ./main $(DATPATH)/run91_01.root 91
	#
	# # 80 GeV ; Amorphous ; 1.5mm C
	# ./main $(DATPATH)/run49_01.root 49
	# ./main $(DATPATH)/run50_01.root 50
	# ./main $(DATPATH)/run51_01.root 51
	# ./main $(DATPATH)/run52_01.root 52
	./main $(DATPATH)/run53_01.root 53 $(PARAMPATH)/1mmC_cuts.txt
	# ./main $(DATPATH)/run56_01.root 56
	# ./main $(DATPATH)/run57_01.root 57
	#
	# # 80 GeV ; Background ; 1.5mm C
	# ./main $(DATPATH)/run48_01.root 48
	# ./main $(DATPATH)/run54_01.root 54
	# ./main $(DATPATH)/run55_01.root 55

main : main.cc analyser.cc

outputsim : test_sim

	# ./test_sim amorphous sim_amorphous1_80GeV_1.5mm_no_background 80GeV_beam_params_1.5mm.txt
	# ./test_sim amorphous sim_amorphous2_80GeV_1.5mm_no_background 80GeV_beam_params_1.5mm.txt
	# ./test_sim amorphous sim_amorphous3_80GeV_1.5mm_no_background 80GeV_beam_params_1.5mm.txt
	# ./test_sim amorphous sim_amorphous4_80GeV_1.5mm_no_background 80GeV_beam_params_1.5mm.txt
	# ./test_sim amorphous sim_amorphous5_80GeV_1.5mm_no_background 80GeV_beam_params_1.5mm.txt
	# ./test_sim amorphous sim_amorphous1_80GeV_1.5mm 80GeV_beam_params_1.5mm.txt
	# ./test_sim amorphous sim_amorphous2_80GeV_1.5mm 80GeV_beam_params_1.5mm.txt
	# ./test_sim amorphous sim_amorphous3_80GeV_1.5mm 80GeV_beam_params_1.5mm.txt
	# ./test_sim amorphous sim_amorphous4_80GeV_1.5mm 80GeV_beam_params_1.5mm.txt
	# ./test_sim amorphous sim_amorphous5_80GeV_1.5mm 80GeV_beam_params_1.5mm.txt
	# ./test_sim background sim_background1_80GeV_1.5mm 80GeV_beam_params_1.5mm.txt
	# ./test_sim background sim_background2_80GeV_1.5mm 80GeV_beam_params_1.5mm.txt
	# ./test_sim background sim_background3_80GeV_1.5mm 80GeV_beam_params_1.5mm.txt
	# ./test_sim background sim_background4_80GeV_1.5mm 80GeV_beam_params_1.5mm.txt
	# ./test_sim background sim_background5_80GeV_1.5mm 80GeV_beam_params_1.5mm.txt
	# # #
	# ./test_sim amorphous sim_amorphous1_80GeV_no_background 80GeV_beam_params.txt
	# ./test_sim amorphous sim_amorphous2_80GeV_no_background 80GeV_beam_params.txt
	# ./test_sim amorphous sim_amorphous3_80GeV_no_background 80GeV_beam_params.txt
	# ./test_sim amorphous sim_amorphous4_80GeV_no_background 80GeV_beam_params.txt
	# ./test_sim amorphous sim_amorphous5_80GeV_no_background 80GeV_beam_params.txt
	# ./test_sim amorphous sim_amorphous1_80GeV 80GeV_beam_params.txt
	# ./test_sim amorphous sim_amorphous2_80GeV 80GeV_beam_params.txt
	# ./test_sim amorphous sim_amorphous3_80GeV 80GeV_beam_params.txt
	# ./test_sim amorphous sim_amorphous4_80GeV 80GeV_beam_params.txt
	# ./test_sim amorphous sim_amorphous5_80GeV 80GeV_beam_params.txt
	# ./test_sim background sim_background1_80GeV 80GeV_beam_params.txt
	# ./test_sim background sim_background2_80GeV 80GeV_beam_params.txt
	# ./test_sim background sim_background3_80GeV 80GeV_beam_params.txt
	# ./test_sim background sim_background4_80GeV 80GeV_beam_params.txt
	# ./test_sim background sim_background5_80GeV 80GeV_beam_params.txt

	# ./test_sim amorphous sim_amorphous1_20GeV_no_background 20GeV_beam_params.txt
	# ./test_sim amorphous sim_amorphous2_20GeV_no_background 20GeV_beam_params.txt
	# ./test_sim amorphous sim_amorphous3_20GeV_no_background 20GeV_beam_params.txt
	# ./test_sim amorphous sim_amorphous4_20GeV_no_background 20GeV_beam_params.txt
	# ./test_sim amorphous sim_amorphous5_20GeV_no_background 20GeV_beam_params.txt
	# ./test_sim amorphous sim_amorphous1_20GeV 20GeV_beam_params.txt
	# ./test_sim amorphous sim_amorphous2_20GeV 20GeV_beam_params.txt
	# ./test_sim amorphous sim_amorphous3_20GeV 20GeV_beam_params.txt
	# ./test_sim amorphous sim_amorphous4_20GeV 20GeV_beam_params.txt
	# ./test_sim amorphous sim_amorphous5_20GeV 20GeV_beam_params.txt
	# ./test_sim background sim_background1_20GeV 20GeV_beam_params.txt
	# ./test_sim background sim_background2_20GeV 20GeV_beam_params.txt
	# ./test_sim background sim_background3_20GeV 20GeV_beam_params.txt
	# ./test_sim background sim_background4_20GeV 20GeV_beam_params.txt
	# ./test_sim background sim_background5_20GeV 20GeV_beam_params.txt
	#
	# ./test_sim amorphous sim_amorphous1_20GeV_1.5mm_no_background 20GeV_beam_params_1.5mm.txt
	# ./test_sim amorphous sim_amorphous2_20GeV_1.5mm_no_background 20GeV_beam_params_1.5mm.txt
	# ./test_sim amorphous sim_amorphous3_20GeV_1.5mm_no_background 20GeV_beam_params_1.5mm.txt
	# ./test_sim amorphous sim_amorphous4_20GeV_1.5mm_no_background 20GeV_beam_params_1.5mm.txt
	# ./test_sim amorphous sim_amorphous5_20GeV_1.5mm_no_background 20GeV_beam_params_1.5mm.txt
	# ./test_sim amorphous sim_amorphous1_20GeV_1.5mm 20GeV_beam_params_1.5mm.txt
	# ./test_sim amorphous sim_amorphous2_20GeV_1.5mm 20GeV_beam_params_1.5mm.txt
	# ./test_sim amorphous sim_amorphous3_20GeV_1.5mm 20GeV_beam_params_1.5mm.txt
	# ./test_sim amorphous sim_amorphous4_20GeV_1.5mm 20GeV_beam_params_1.5mm.txt
	# ./test_sim amorphous sim_amorphous5_20GeV_1.5mm 20GeV_beam_params_1.5mm.txt
	# ./test_sim background sim_background1_20GeV_1.5mm 20GeV_beam_params_1.5mm.txt
	# ./test_sim background sim_background2_20GeV_1.5mm 20GeV_beam_params_1.5mm.txt
	# ./test_sim background sim_background3_20GeV_1.5mm 20GeV_beam_params_1.5mm.txt
	# ./test_sim background sim_background4_20GeV_1.5mm 20GeV_beam_params_1.5mm.txt
	# ./test_sim background sim_background5_20GeV_1.5mm 20GeV_beam_params_1.5mm.txt

	./test_sim amorphous sim_amorphous1_40GeV_no_background 40GeV_beam_params.txt
	# ./test_sim amorphous sim_amorphous2_40GeV_no_background 40GeV_beam_params.txt
	# ./test_sim amorphous sim_amorphous3_40GeV_no_background 40GeV_beam_params.txt
	# ./test_sim amorphous sim_amorphous4_40GeV_no_background 40GeV_beam_params.txt
	# ./test_sim amorphous sim_amorphous5_40GeV_no_background 40GeV_beam_params.txt
	# ./test_sim amorphous sim_amorphous1_40GeV 40GeV_beam_params.txt
	# ./test_sim amorphous sim_amorphous2_40GeV 40GeV_beam_params.txt
	# ./test_sim amorphous sim_amorphous3_40GeV 40GeV_beam_params.txt
	# ./test_sim amorphous sim_amorphous4_40GeV 40GeV_beam_params.txt
	# ./test_sim amorphous sim_amorphous5_40GeV 40GeV_beam_params.txt
	# ./test_sim background sim_background1_40GeV 40GeV_beam_params.txt
	# ./test_sim background sim_background2_40GeV 40GeV_beam_params.txt
	# ./test_sim background sim_background3_40GeV 40GeV_beam_params.txt
	# ./test_sim background sim_background4_40GeV 40GeV_beam_params.txt
	# ./test_sim background sim_background5_40GeV 40GeV_beam_params.txt
	#
	# ./test_sim amorphous sim_amorphous1_40GeV_1.5mm_no_background 40GeV_beam_params_1.5mm.txt
	# ./test_sim amorphous sim_amorphous2_40GeV_1.5mm_no_background 40GeV_beam_params_1.5mm.txt
	# ./test_sim amorphous sim_amorphous3_40GeV_1.5mm_no_background 40GeV_beam_params_1.5mm.txt
	# ./test_sim amorphous sim_amorphous4_40GeV_1.5mm_no_background 40GeV_beam_params_1.5mm.txt
	# ./test_sim amorphous sim_amorphous5_40GeV_1.5mm_no_background 40GeV_beam_params_1.5mm.txt
	# ./test_sim amorphous sim_amorphous1_40GeV_1.5mm 40GeV_beam_params_1.5mm.txt
	# ./test_sim amorphous sim_amorphous2_40GeV_1.5mm 40GeV_beam_params_1.5mm.txt
	# ./test_sim amorphous sim_amorphous3_40GeV_1.5mm 40GeV_beam_params_1.5mm.txt
	# ./test_sim amorphous sim_amorphous4_40GeV_1.5mm 40GeV_beam_params_1.5mm.txt
	# ./test_sim amorphous sim_amorphous5_40GeV_1.5mm 40GeV_beam_params_1.5mm.txt
	# ./test_sim background sim_background1_40GeV_1.5mm 40GeV_beam_params_1.5mm.txt
	# ./test_sim background sim_background2_40GeV_1.5mm 40GeV_beam_params_1.5mm.txt
	# ./test_sim background sim_background3_40GeV_1.5mm 40GeV_beam_params_1.5mm.txt
	# ./test_sim background sim_background4_40GeV_1.5mm 40GeV_beam_params_1.5mm.txt
	# ./test_sim background sim_background5_40GeV_1.5mm 40GeV_beam_params_1.5mm.txt

test_sim : test_sim.cc simulator.cc analyser.cc

clean:
	rm -f main main.o test_sim test_sim.o
