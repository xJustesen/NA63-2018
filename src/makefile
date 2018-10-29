CC	:= g++
CXXFLAGS += -std=c++14 -Ofast -fopenmp -Wall -larmadillo -I/home/christian/Root/include
LDLIBS += -lm -march=native -larmadillo -fopenmp $(shell /home/christian/Documents/root/bin/root-config --libs)
DATPATH := /home/christian/Documents/CERN2018NA63
PARAMPATH := /home/christian/Dropbox/speciale/code/beamParameters
SPECPATHEXT := /media/christian/Elements/DataAnalysis/radiation/spectres
SPECPATHINT := /home/christian/Dropbox/speciale/code/crystalSimulations

all : outputdata

outputdata : main
	# # 20 GeV ; Amorphous ; 1mm C
	# ./main $(DATPATH)/run109_01.root 109
	# ./main $(DATPATH)/run115_01.root 115
	# #
	# 20 GeV ; Aligned ; 1mm C
	# ./main $(DATPATH)/run103_01.root 103 $(PARAMPATH)/20GeV_beam_params.txt
	# ./main $(DATPATH)/run104_01.root 104 $(PARAMPATH)/20GeV_beam_params.txt
	# ./main $(DATPATH)/run112_01.root 112 $(PARAMPATH)/20GeV_beam_params.txt
	# ./main $(DATPATH)/run113_01.root 113 $(PARAMPATH)/20GeV_beam_params.txt
	# ./main $(DATPATH)/run114_01.root 114 $(PARAMPATH)/20GeV_beam_params.txt

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

	# 20 GeV ; Aligned ; 1.5mm C
	# ./main $(DATPATH)/run61_01.root 61 $(PARAMPATH)/20GeV_beam_params_1.5mm.txt
	# ./main $(DATPATH)/run62_01.root 62 $(PARAMPATH)/20GeV_beam_params_1.5mm.txt
	# ./main $(DATPATH)/run63_01.root 63 $(PARAMPATH)/20GeV_beam_params_1.5mm.txt
	# ./main $(DATPATH)/run64_01.root 64 $(PARAMPATH)/20GeV_beam_params_1.5mm.txt


	# # 20 GeV ; Background ; 1.5mm C
	# ./main $(DATPATH)/run60_01.root 60
	# ./main $(DATPATH)/run65_01.root 65
	# #
	#
	# 40 GeV ; Aligned ; 1mm C
	# ./main $(DATPATH)/run71_01.root 71 $(PARAMPATH)/40GeV_beam_params.txt #$(PARAMPATH)/1mmC_cuts.txt
	# ./main $(DATPATH)/run72_01.root 72 $(PARAMPATH)/40GeV_beam_params.txt #$(PARAMPATH)/1mmC_cuts.txt
	# ./main $(DATPATH)/run79_01.root 79 $(PARAMPATH)/40GeV_beam_params.txt #$(PARAMPATH)/1mmC_cuts.txt
	# ./main $(DATPATH)/run80_01.root 80 $(PARAMPATH)/40GeV_beam_params.txt #$(PARAMPATH)/1mmC_cuts.txt
	# ./main $(DATPATH)/run81_01.root 81 $(PARAMPATH)/40GeV_beam_params.txt #$(PARAMPATH)/1mmC_cuts.txt
	#
	# # 40 GeV ; Amorphous ; 1mm C
	# ./main $(DATPATH)/run73_01.root 73 #$(PARAMPATH)/1mmC_cuts.txt
	# ./main $(DATPATH)/run75_01.root 75
	# ./main $(DATPATH)/run76_01.root 76
	# ./main $(DATPATH)/run77_01.root 77
	# ./main $(DATPATH)/run78_01.root 78
	#
	# # 40 GeV ; Background ; 1mm C
	# ./main $(DATPATH)/run74_01.root 74

	# 40 GeV ; Aligned ; 1.5mm C
	# ./main $(DATPATH)/run30_01.root 30  $(PARAMPATH)/40GeV_beam_params_1.5mm.txt
	# ./main $(DATPATH)/run35_01.root 35  $(PARAMPATH)/40GeV_beam_params_1.5mm.txt
	# ./main $(DATPATH)/run36_01.root 36  $(PARAMPATH)/40GeV_beam_params_1.5mm.txt
	# ./main $(DATPATH)/run37_01.root 37  $(PARAMPATH)/40GeV_beam_params_1.5mm.txt
	# ./main $(DATPATH)/run38_01.root 38  $(PARAMPATH)/40GeV_beam_params_1.5mm.txt
	#
	# #
	# # 40 GeV ; Amorphous ; 1.5mm C
	# ./main $(DATPATH)/run32_01.root 32
	# ./main $(DATPATH)/run33_01.root 33
	# ./main $(DATPATH)/run34_01.root 34
	# ./main $(DATPATH)/run39_01.root 39
	# ./main $(DATPATH)/run40_01.root 40
	# ./main $(DATPATH)/run41_01.root 41
	# ./main $(DATPATH)/run43_01.root 43
	# #
	# # # 40 GeV ; Background ; 1.5mm C
	# ./main $(DATPATH)/run31_01.root 31
	# #
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

	# # 80 GeV ; Aligned ; 1mm C
	# ./main $(DATPATH)/run84_01.root 84 $(PARAMPATH)/80GeV_beam_params.txt
	# ./main $(DATPATH)/run92_01.root 92 $(PARAMPATH)/80GeV_beam_params.txt
	# ./main $(DATPATH)/run93_01.root 93 $(PARAMPATH)/80GeV_beam_params.txt
	# ./main $(DATPATH)/run94_01.root 94 $(PARAMPATH)/80GeV_beam_params.txt
	# ./main $(DATPATH)/run95_01.root 95 $(PARAMPATH)/80GeV_beam_params.txt

	# 80 GeV ; Aligned ; 1.5mm C
	./main $(DATPATH)/run46_01.root 46  $(PARAMPATH)/80GeV_beam_params_1.5mm.txt
	# ./main $(DATPATH)/run47_01.root 47  $(PARAMPATH)/80GeV_beam_params_1.5mm.txt

	# # # 80 GeV ; Amorphous ; 1.5mm C
	# ./main $(DATPATH)/run49_01.root 49
	# ./main $(DATPATH)/run50_01.root 50
	# ./main $(DATPATH)/run51_01.root 51
	# ./main $(DATPATH)/run52_01.root 52
	# ./main $(DATPATH)/run53_01.root 53
	# ./main $(DATPATH)/run56_01.root 56
	# ./main $(DATPATH)/run57_01.root 57
	# #
	# # # 80 GeV ; Background ; 1.5mm C
	# ./main $(DATPATH)/run48_01.root 48
	# ./main $(DATPATH)/run54_01.root 54
	# ./main $(DATPATH)/run55_01.root 55
	#
	# # Alignment
	# ./main $(DATPATH)/run18_01.root 18

	# test run
	# ./main $(DATPATH)/run4_01.root 4

main : main.cc analyser.cc

outputsim : test_sim

	# ./test_sim amorphous sim_amorphous1_80GeV_1.5mm_no_background 80GeV_beam_params_1.5mm.txt
	# ./test_sim amorphous sim_amorphous2_80GeV_1.5mm_no_background 80GeV_beam_params_1.5mm.txt
	# ./test_sim amorphous sim_amorphous3_80GeV_1.5mm_no_background 80GeV_beam_params_1.5mm.txt
	# ./test_sim amorphous sim_amorphous4_80GeV_1.5mm_no_background 80GeV_beam_params_1.5mm.txt
	# ./test_sim amorphous sim_amorphous5_80GeV_1.5mm_no_background 80GeV_beam_params_1.5mm.txt

	./test_sim amorphous sim_amorphous1_80GeV_1.5mm 80GeV_beam_params_1.5mm.txt
	./test_sim amorphous sim_amorphous2_80GeV_1.5mm 80GeV_beam_params_1.5mm.txt
	./test_sim amorphous sim_amorphous3_80GeV_1.5mm 80GeV_beam_params_1.5mm.txt
	./test_sim amorphous sim_amorphous4_80GeV_1.5mm 80GeV_beam_params_1.5mm.txt
	./test_sim amorphous sim_amorphous5_80GeV_1.5mm 80GeV_beam_params_1.5mm.txt

	./test_sim background sim_background1_80GeV_1.5mm 80GeV_beam_params_1.5mm.txt
	./test_sim background sim_background2_80GeV_1.5mm 80GeV_beam_params_1.5mm.txt
	./test_sim background sim_background3_80GeV_1.5mm 80GeV_beam_params_1.5mm.txt
	./test_sim background sim_background4_80GeV_1.5mm 80GeV_beam_params_1.5mm.txt
	./test_sim background sim_background5_80GeV_1.5mm 80GeV_beam_params_1.5mm.txt

	# ./test_sim aligned sim_aligned1_80GeV_woshot_1.5mm 80GeV_beam_params_1.5mm.txt $(SPECPATHINT)/sum_initials80GeV_1.5mm_full_peak_alt_woshot.txt $(SPECPATHINT)/sum_angles80GeV_1.5mm_full_peak_alt_woshot.txt
	# ./test_sim aligned sim_aligned2_80GeV_woshot_1.5mm 80GeV_beam_params_1.5mm.txt $(SPECPATHINT)/sum_initials80GeV_1.5mm_full_peak_alt_woshot.txt $(SPECPATHINT)/sum_angles80GeV_1.5mm_full_peak_alt_woshot.txt
	# ./test_sim aligned sim_aligned3_80GeV_woshot_1.5mm 80GeV_beam_params_1.5mm.txt $(SPECPATHINT)/sum_initials80GeV_1.5mm_full_peak_alt_woshot.txt $(SPECPATHINT)/sum_angles80GeV_1.5mm_full_peak_alt_woshot.txt
	# ./test_sim aligned sim_aligned4_80GeV_woshot_1.5mm 80GeV_beam_params_1.5mm.txt $(SPECPATHINT)/sum_initials80GeV_1.5mm_full_peak_alt_woshot.txt $(SPECPATHINT)/sum_angles80GeV_1.5mm_full_peak_alt_woshot.txt
	# ./test_sim aligned sim_aligned5_80GeV_woshot_1.5mm 80GeV_beam_params_1.5mm.txt $(SPECPATHINT)/sum_initials80GeV_1.5mm_full_peak_alt_woshot.txt $(SPECPATHINT)/sum_angles80GeV_1.5mm_full_peak_alt_woshot.txt
	#
	# ./test_sim aligned sim_aligned1_80GeV_worr_1.5mm 80GeV_beam_params_1.5mm.txt $(SPECPATHINT)/sum_initials80GeV_1.5mm_full_peak_alt_worr.txt $(SPECPATHINT)/sum_angles80GeV_1.5mm_full_peak_alt_worr.txt
	# ./test_sim aligned sim_aligned2_80GeV_worr_1.5mm 80GeV_beam_params_1.5mm.txt $(SPECPATHINT)/sum_initials80GeV_1.5mm_full_peak_alt_worr.txt $(SPECPATHINT)/sum_angles80GeV_1.5mm_full_peak_alt_worr.txt
	# ./test_sim aligned sim_aligned3_80GeV_worr_1.5mm 80GeV_beam_params_1.5mm.txt $(SPECPATHINT)/sum_initials80GeV_1.5mm_full_peak_alt_worr.txt $(SPECPATHINT)/sum_angles80GeV_1.5mm_full_peak_alt_worr.txt
	# ./test_sim aligned sim_aligned4_80GeV_worr_1.5mm 80GeV_beam_params_1.5mm.txt $(SPECPATHINT)/sum_initials80GeV_1.5mm_full_peak_alt_worr.txt $(SPECPATHINT)/sum_angles80GeV_1.5mm_full_peak_alt_worr.txt
	# ./test_sim aligned sim_aligned5_80GeV_worr_1.5mm 80GeV_beam_params_1.5mm.txt $(SPECPATHINT)/sum_initials80GeV_1.5mm_full_peak_alt_worr.txt $(SPECPATHINT)/sum_angles80GeV_1.5mm_full_peak_alt_worr.txt
	#
	# ./test_sim aligned sim_aligned1_80GeV_1.5mm 80GeV_beam_params_1.5mm.txt $(SPECPATHINT)/sum_initials80GeV_1.5mm_full_peak_alt.txt $(SPECPATHINT)/sum_angles80GeV_1.5mm_full_peak_alt.txt
	# ./test_sim aligned sim_aligned2_80GeV_1.5mm 80GeV_beam_params_1.5mm.txt $(SPECPATHINT)/sum_initials80GeV_1.5mm_full_peak_alt.txt $(SPECPATHINT)/sum_angles80GeV_1.5mm_full_peak_alt.txt
	# ./test_sim aligned sim_aligned3_80GeV_1.5mm 80GeV_beam_params_1.5mm.txt $(SPECPATHINT)/sum_initials80GeV_1.5mm_full_peak_alt.txt $(SPECPATHINT)/sum_angles80GeV_1.5mm_full_peak_alt.txt
	# ./test_sim aligned sim_aligned4_80GeV_1.5mm 80GeV_beam_params_1.5mm.txt $(SPECPATHINT)/sum_initials80GeV_1.5mm_full_peak_alt.txt $(SPECPATHINT)/sum_angles80GeV_1.5mm_full_peak_alt.txt
	# ./test_sim aligned sim_aligned5_80GeV_1.5mm 80GeV_beam_params_1.5mm.txt $(SPECPATHINT)/sum_initials80GeV_1.5mm_full_peak_alt.txt $(SPECPATHINT)/sum_angles80GeV_1.5mm_full_peak_alt.txt

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
	#
	# ./test_sim background sim_background1_80GeV 80GeV_beam_params.txt
	# ./test_sim background sim_background2_80GeV 80GeV_beam_params.txt
	# ./test_sim background sim_background3_80GeV 80GeV_beam_params.txt
	# ./test_sim background sim_background4_80GeV 80GeV_beam_params.txt
	# ./test_sim background sim_background5_80GeV 80GeV_beam_params.txt

	# ./test_sim aligned sim_aligned1_80GeV_woshot 80GeV_beam_params.txt $(SPECPATHINT)/sum_initials80GeV_1mm_full_peak_alt_woshot.txt $(SPECPATHINT)/sum_angles80GeV_1mm_full_peak_alt_woshot.txt
	# ./test_sim aligned sim_aligned2_80GeV_woshot 80GeV_beam_params.txt $(SPECPATHINT)/sum_initials80GeV_1mm_full_peak_alt_woshot.txt $(SPECPATHINT)/sum_angles80GeV_1mm_full_peak_alt_woshot.txt
	# ./test_sim aligned sim_aligned3_80GeV_woshot 80GeV_beam_params.txt $(SPECPATHINT)/sum_initials80GeV_1mm_full_peak_alt_woshot.txt $(SPECPATHINT)/sum_angles80GeV_1mm_full_peak_alt_woshot.txt
	# ./test_sim aligned sim_aligned4_80GeV_woshot 80GeV_beam_params.txt $(SPECPATHINT)/sum_initials80GeV_1mm_full_peak_alt_woshot.txt $(SPECPATHINT)/sum_angles80GeV_1mm_full_peak_alt_woshot.txt
	# ./test_sim aligned sim_aligned5_80GeV_woshot 80GeV_beam_params.txt $(SPECPATHINT)/sum_initials80GeV_1mm_full_peak_alt_woshot.txt $(SPECPATHINT)/sum_angles80GeV_1mm_full_peak_alt_woshot.txt
	#
	# ./test_sim aligned sim_aligned1_80GeV_worr 80GeV_beam_params.txt $(SPECPATHINT)/sum_initials80GeV_1mm_full_peak_alt_worr.txt $(SPECPATHINT)/sum_angles80GeV_1mm_full_peak_alt_worr.txt
	# ./test_sim aligned sim_aligned2_80GeV_worr 80GeV_beam_params.txt $(SPECPATHINT)/sum_initials80GeV_1mm_full_peak_alt_worr.txt $(SPECPATHINT)/sum_angles80GeV_1mm_full_peak_alt_worr.txt
	# ./test_sim aligned sim_aligned3_80GeV_worr 80GeV_beam_params.txt $(SPECPATHINT)/sum_initials80GeV_1mm_full_peak_alt_worr.txt $(SPECPATHINT)/sum_angles80GeV_1mm_full_peak_alt_worr.txt
	# ./test_sim aligned sim_aligned4_80GeV_worr 80GeV_beam_params.txt $(SPECPATHINT)/sum_initials80GeV_1mm_full_peak_alt_worr.txt $(SPECPATHINT)/sum_angles80GeV_1mm_full_peak_alt_worr.txt
	# ./test_sim aligned sim_aligned5_80GeV_worr 80GeV_beam_params.txt $(SPECPATHINT)/sum_initials80GeV_1mm_full_peak_alt_worr.txt $(SPECPATHINT)/sum_angles80GeV_1mm_full_peak_alt_worr.txt
	#
	# ./test_sim aligned sim_aligned1_80GeV 80GeV_beam_params.txt $(SPECPATHINT)/sum_initials80GeV_1mm_full_peak_alt.txt $(SPECPATHINT)/sum_angles80GeV_1mm_full_peak_alt.txt
	# ./test_sim aligned sim_aligned2_80GeV 80GeV_beam_params.txt $(SPECPATHINT)/sum_initials80GeV_1mm_full_peak_alt.txt $(SPECPATHINT)/sum_angles80GeV_1mm_full_peak_alt.txt
	# ./test_sim aligned sim_aligned3_80GeV 80GeV_beam_params.txt $(SPECPATHINT)/sum_initials80GeV_1mm_full_peak_alt.txt $(SPECPATHINT)/sum_angles80GeV_1mm_full_peak_alt.txt
	# ./test_sim aligned sim_aligned4_80GeV 80GeV_beam_params.txt $(SPECPATHINT)/sum_initials80GeV_1mm_full_peak_alt.txt $(SPECPATHINT)/sum_angles80GeV_1mm_full_peak_alt.txt
	# ./test_sim aligned sim_aligned5_80GeV 80GeV_beam_params.txt $(SPECPATHINT)/sum_initials80GeV_1mm_full_peak_alt.txt $(SPECPATHINT)/sum_angles80GeV_1mm_full_peak_alt.txt

	# ./test_sim amorphous sim_amorphous1_20GeV_no_background 20GeV_beam_params.txt
	# ./test_sim amorphous sim_amorphous2_20GeV_no_background 20GeV_beam_params.txt
	# ./test_sim amorphous sim_amorphous3_20GeV_no_background 20GeV_beam_params.txt
	# ./test_sim amorphous sim_amorphous4_20GeV_no_background 20GeV_beam_params.txt
	# ./test_sim amorphous sim_amorphous5_20GeV_no_background 20GeV_beam_params.txt

	./test_sim amorphous sim_amorphous1_20GeV 20GeV_beam_params.txt
	./test_sim amorphous sim_amorphous2_20GeV 20GeV_beam_params.txt
	./test_sim amorphous sim_amorphous3_20GeV 20GeV_beam_params.txt
	./test_sim amorphous sim_amorphous4_20GeV 20GeV_beam_params.txt
	./test_sim amorphous sim_amorphous5_20GeV 20GeV_beam_params.txt

	./test_sim background sim_background1_20GeV 20GeV_beam_params.txt
	./test_sim background sim_background2_20GeV 20GeV_beam_params.txt
	./test_sim background sim_background3_20GeV 20GeV_beam_params.txt
	./test_sim background sim_background4_20GeV 20GeV_beam_params.txt
	./test_sim background sim_background5_20GeV 20GeV_beam_params.txt

	# ./test_sim aligned sim_aligned1_20GeV_woshot 20GeV_beam_params.txt $(SPECPATHINT)/sum_initials20GeV_1mm_full_peak_alt_woshot.txt $(SPECPATHINT)/sum_angles20GeV_1mm_full_peak_alt_woshot.txt
	# ./test_sim aligned sim_aligned2_20GeV_woshot 20GeV_beam_params.txt $(SPECPATHINT)/sum_initials20GeV_1mm_full_peak_alt_woshot.txt $(SPECPATHINT)/sum_angles20GeV_1mm_full_peak_alt_woshot.txt
	# ./test_sim aligned sim_aligned3_20GeV_woshot 20GeV_beam_params.txt $(SPECPATHINT)/sum_initials20GeV_1mm_full_peak_alt_woshot.txt $(SPECPATHINT)/sum_angles20GeV_1mm_full_peak_alt_woshot.txt
	# ./test_sim aligned sim_aligned4_20GeV_woshot 20GeV_beam_params.txt $(SPECPATHINT)/sum_initials20GeV_1mm_full_peak_alt_woshot.txt $(SPECPATHINT)/sum_angles20GeV_1mm_full_peak_alt_woshot.txt
	# ./test_sim aligned sim_aligned5_20GeV_woshot 20GeV_beam_params.txt $(SPECPATHINT)/sum_initials20GeV_1mm_full_peak_alt_woshot.txt $(SPECPATHINT)/sum_angles20GeV_1mm_full_peak_alt_woshot.txt
	#
	# ./test_sim aligned sim_aligned1_20GeV_worr 20GeV_beam_params.txt $(SPECPATHINT)/sum_initials20GeV_1mm_full_peak_alt_worr.txt $(SPECPATHINT)/sum_angles20GeV_1mm_full_peak_alt_worr.txt
	# ./test_sim aligned sim_aligned2_20GeV_worr 20GeV_beam_params.txt $(SPECPATHINT)/sum_initials20GeV_1mm_full_peak_alt_worr.txt $(SPECPATHINT)/sum_angles20GeV_1mm_full_peak_alt_worr.txt
	# ./test_sim aligned sim_aligned3_20GeV_worr 20GeV_beam_params.txt $(SPECPATHINT)/sum_initials20GeV_1mm_full_peak_alt_worr.txt $(SPECPATHINT)/sum_angles20GeV_1mm_full_peak_alt_worr.txt
	# ./test_sim aligned sim_aligned4_20GeV_worr 20GeV_beam_params.txt $(SPECPATHINT)/sum_initials20GeV_1mm_full_peak_alt_worr.txt $(SPECPATHINT)/sum_angles20GeV_1mm_full_peak_alt_worr.txt
	# ./test_sim aligned sim_aligned5_20GeV_worr 20GeV_beam_params.txt $(SPECPATHINT)/sum_initials20GeV_1mm_full_peak_alt_worr.txt $(SPECPATHINT)/sum_angles20GeV_1mm_full_peak_alt_worr.txt
	#
	# ./test_sim aligned sim_aligned1_20GeV 20GeV_beam_params.txt $(SPECPATHINT)/sum_initials20GeV_1mm_full_peak_alt.txt $(SPECPATHINT)/sum_angles20GeV_1mm_full_peak_alt.txt
	# ./test_sim aligned sim_aligned2_20GeV 20GeV_beam_params.txt $(SPECPATHINT)/sum_initials20GeV_1mm_full_peak_alt.txt $(SPECPATHINT)/sum_angles20GeV_1mm_full_peak_alt.txt
	# ./test_sim aligned sim_aligned3_20GeV 20GeV_beam_params.txt $(SPECPATHINT)/sum_initials20GeV_1mm_full_peak_alt.txt $(SPECPATHINT)/sum_angles20GeV_1mm_full_peak_alt.txt
	# ./test_sim aligned sim_aligned4_20GeV 20GeV_beam_params.txt $(SPECPATHINT)/sum_initials20GeV_1mm_full_peak_alt.txt $(SPECPATHINT)/sum_angles20GeV_1mm_full_peak_alt.txt
	# ./test_sim aligned sim_aligned5_20GeV 20GeV_beam_params.txt $(SPECPATHINT)/sum_initials20GeV_1mm_full_peak_alt.txt $(SPECPATHINT)/sum_angles20GeV_1mm_full_peak_alt.txt

	# ./test_sim amorphous sim_amorphous1_20GeV_1.5mm_no_background 20GeV_beam_params_1.5mm.txt
	# ./test_sim amorphous sim_amorphous2_20GeV_1.5mm_no_background 20GeV_beam_params_1.5mm.txt
	# ./test_sim amorphous sim_amorphous3_20GeV_1.5mm_no_background 20GeV_beam_params_1.5mm.txt
	# ./test_sim amorphous sim_amorphous4_20GeV_1.5mm_no_background 20GeV_beam_params_1.5mm.txt
	# ./test_sim amorphous sim_amorphous5_20GeV_1.5mm_no_background 20GeV_beam_params_1.5mm.txt

	./test_sim amorphous sim_amorphous1_20GeV_1.5mm 20GeV_beam_params_1.5mm.txt
	./test_sim amorphous sim_amorphous2_20GeV_1.5mm 20GeV_beam_params_1.5mm.txt
	./test_sim amorphous sim_amorphous3_20GeV_1.5mm 20GeV_beam_params_1.5mm.txt
	./test_sim amorphous sim_amorphous4_20GeV_1.5mm 20GeV_beam_params_1.5mm.txt
	./test_sim amorphous sim_amorphous5_20GeV_1.5mm 20GeV_beam_params_1.5mm.txt

	./test_sim background sim_background1_20GeV_1.5mm 20GeV_beam_params_1.5mm.txt
	./test_sim background sim_background2_20GeV_1.5mm 20GeV_beam_params_1.5mm.txt
	./test_sim background sim_background3_20GeV_1.5mm 20GeV_beam_params_1.5mm.txt
	./test_sim background sim_background4_20GeV_1.5mm 20GeV_beam_params_1.5mm.txt
	./test_sim background sim_background5_20GeV_1.5mm 20GeV_beam_params_1.5mm.txt
	#

		# ./test_sim aligned sim_aligned1_20GeV_woshot_1.5mm 20GeV_beam_params_1.5mm.txt $(SPECPATHINT)/sum_initials20GeV_1.5mm_full_peak_alt_woshot.txt $(SPECPATHINT)/sum_angles20GeV_1.5mm_full_peak_alt_woshot.txt
		# ./test_sim aligned sim_aligned2_20GeV_woshot_1.5mm 20GeV_beam_params_1.5mm.txt $(SPECPATHINT)/sum_initials20GeV_1.5mm_full_peak_alt_woshot.txt $(SPECPATHINT)/sum_angles20GeV_1.5mm_full_peak_alt_woshot.txt
		# ./test_sim aligned sim_aligned3_20GeV_woshot_1.5mm 20GeV_beam_params_1.5mm.txt $(SPECPATHINT)/sum_initials20GeV_1.5mm_full_peak_alt_woshot.txt $(SPECPATHINT)/sum_angles20GeV_1.5mm_full_peak_alt_woshot.txt
		# ./test_sim aligned sim_aligned4_20GeV_woshot_1.5mm 20GeV_beam_params_1.5mm.txt $(SPECPATHINT)/sum_initials20GeV_1.5mm_full_peak_alt_woshot.txt $(SPECPATHINT)/sum_angles20GeV_1.5mm_full_peak_alt_woshot.txt
		# ./test_sim aligned sim_aligned5_20GeV_woshot_1.5mm 20GeV_beam_params_1.5mm.txt $(SPECPATHINT)/sum_initials20GeV_1.5mm_full_peak_alt_woshot.txt $(SPECPATHINT)/sum_angles20GeV_1.5mm_full_peak_alt_woshot.txt
		#
		# ./test_sim aligned sim_aligned1_20GeV_worr_1.5mm 20GeV_beam_params_1.5mm.txt $(SPECPATHINT)/sum_initials20GeV_1.5mm_full_peak_alt_worr.txt $(SPECPATHINT)/sum_angles20GeV_1.5mm_full_peak_alt_worr.txt
		# ./test_sim aligned sim_aligned2_20GeV_worr_1.5mm 20GeV_beam_params_1.5mm.txt $(SPECPATHINT)/sum_initials20GeV_1.5mm_full_peak_alt_worr.txt $(SPECPATHINT)/sum_angles20GeV_1.5mm_full_peak_alt_worr.txt
		# ./test_sim aligned sim_aligned3_20GeV_worr_1.5mm 20GeV_beam_params_1.5mm.txt $(SPECPATHINT)/sum_initials20GeV_1.5mm_full_peak_alt_worr.txt $(SPECPATHINT)/sum_angles20GeV_1.5mm_full_peak_alt_worr.txt
		# ./test_sim aligned sim_aligned4_20GeV_worr_1.5mm 20GeV_beam_params_1.5mm.txt $(SPECPATHINT)/sum_initials20GeV_1.5mm_full_peak_alt_worr.txt $(SPECPATHINT)/sum_angles20GeV_1.5mm_full_peak_alt_worr.txt
		# ./test_sim aligned sim_aligned5_20GeV_worr_1.5mm 20GeV_beam_params_1.5mm.txt $(SPECPATHINT)/sum_initials20GeV_1.5mm_full_peak_alt_worr.txt $(SPECPATHINT)/sum_angles20GeV_1.5mm_full_peak_alt_worr.txt
		#
		# ./test_sim aligned sim_aligned1_20GeV_1.5mm 20GeV_beam_params_1.5mm.txt $(SPECPATHINT)/sum_initials20GeV_1.5mm_full_peak_alt.txt $(SPECPATHINT)/sum_angles20GeV_1.5mm_full_peak_alt.txt
		# ./test_sim aligned sim_aligned2_20GeV_1.5mm 20GeV_beam_params_1.5mm.txt $(SPECPATHINT)/sum_initials20GeV_1.5mm_full_peak_alt.txt $(SPECPATHINT)/sum_angles20GeV_1.5mm_full_peak_alt.txt
		# ./test_sim aligned sim_aligned3_20GeV_1.5mm 20GeV_beam_params_1.5mm.txt $(SPECPATHINT)/sum_initials20GeV_1.5mm_full_peak_alt.txt $(SPECPATHINT)/sum_angles20GeV_1.5mm_full_peak_alt.txt
		# ./test_sim aligned sim_aligned4_20GeV_1.5mm 20GeV_beam_params_1.5mm.txt $(SPECPATHINT)/sum_initials20GeV_1.5mm_full_peak_alt.txt $(SPECPATHINT)/sum_angles20GeV_1.5mm_full_peak_alt.txt
		# ./test_sim aligned sim_aligned5_20GeV_1.5mm 20GeV_beam_params_1.5mm.txt $(SPECPATHINT)/sum_initials20GeV_1.5mm_full_peak_alt.txt $(SPECPATHINT)/sum_angles20GeV_1.5mm_full_peak_alt.txt
	#
	# ./test_sim amorphous sim_amorphous1_40GeV_no_background 40GeV_beam_params.txt
	# ./test_sim amorphous sim_amorphous2_40GeV_no_background 40GeV_beam_params.txt
	# ./test_sim amorphous sim_amorphous3_40GeV_no_background 40GeV_beam_params.txt
	# ./test_sim amorphous sim_amorphous4_40GeV_no_background 40GeV_beam_params.txt
	# ./test_sim amorphous sim_amorphous5_40GeV_no_background 40GeV_beam_params.txt

	./test_sim amorphous sim_amorphous1_40GeV 40GeV_beam_params.txt
	./test_sim amorphous sim_amorphous2_40GeV 40GeV_beam_params.txt
	./test_sim amorphous sim_amorphous3_40GeV 40GeV_beam_params.txt
	./test_sim amorphous sim_amorphous4_40GeV 40GeV_beam_params.txt
	./test_sim amorphous sim_amorphous5_40GeV 40GeV_beam_params.txt

	./test_sim background sim_background1_40GeV 40GeV_beam_params.txt
	./test_sim background sim_background2_40GeV 40GeV_beam_params.txt
	./test_sim background sim_background3_40GeV 40GeV_beam_params.txt
	./test_sim background sim_background4_40GeV 40GeV_beam_params.txt
	./test_sim background sim_background5_40GeV 40GeV_beam_params.txt

	# ./test_sim aligned sim_aligned1_40GeV_woshot 40GeV_beam_params.txt $(SPECPATHINT)/sum_initials40GeV_1mm_full_peak_alt_woshot.txt $(SPECPATHINT)/sum_angles40GeV_1mm_full_peak_alt_woshot.txt
	# ./test_sim aligned sim_aligned2_40GeV_woshot 40GeV_beam_params.txt $(SPECPATHINT)/sum_initials40GeV_1mm_full_peak_alt_woshot.txt $(SPECPATHINT)/sum_angles40GeV_1mm_full_peak_alt_woshot.txt
	# ./test_sim aligned sim_aligned3_40GeV_woshot 40GeV_beam_params.txt $(SPECPATHINT)/sum_initials40GeV_1mm_full_peak_alt_woshot.txt $(SPECPATHINT)/sum_angles40GeV_1mm_full_peak_alt_woshot.txt
	# ./test_sim aligned sim_aligned4_40GeV_woshot 40GeV_beam_params.txt $(SPECPATHINT)/sum_initials40GeV_1mm_full_peak_alt_woshot.txt $(SPECPATHINT)/sum_angles40GeV_1mm_full_peak_alt_woshot.txt
	# ./test_sim aligned sim_aligned5_40GeV_woshot 40GeV_beam_params.txt $(SPECPATHINT)/sum_initials40GeV_1mm_full_peak_alt_woshot.txt $(SPECPATHINT)/sum_angles40GeV_1mm_full_peak_alt_woshot.txt
	#
	# ./test_sim aligned sim_aligned1_40GeV_worr 40GeV_beam_params.txt $(SPECPATHINT)/sum_initials40GeV_1mm_full_peak_alt_worr.txt $(SPECPATHINT)/sum_angles40GeV_1mm_full_peak_alt_worr.txt
	# ./test_sim aligned sim_aligned2_40GeV_worr 40GeV_beam_params.txt $(SPECPATHINT)/sum_initials40GeV_1mm_full_peak_alt_worr.txt $(SPECPATHINT)/sum_angles40GeV_1mm_full_peak_alt_worr.txt
	# ./test_sim aligned sim_aligned3_40GeV_worr 40GeV_beam_params.txt $(SPECPATHINT)/sum_initials40GeV_1mm_full_peak_alt_worr.txt $(SPECPATHINT)/sum_angles40GeV_1mm_full_peak_alt_worr.txt
	# ./test_sim aligned sim_aligned4_40GeV_worr 40GeV_beam_params.txt $(SPECPATHINT)/sum_initials40GeV_1mm_full_peak_alt_worr.txt $(SPECPATHINT)/sum_angles40GeV_1mm_full_peak_alt_worr.txt
	# ./test_sim aligned sim_aligned5_40GeV_worr 40GeV_beam_params.txt $(SPECPATHINT)/sum_initials40GeV_1mm_full_peak_alt_worr.txt $(SPECPATHINT)/sum_angles40GeV_1mm_full_peak_alt_worr.txt
	# #
	# ./test_sim aligned sim_aligned1_40GeV 40GeV_beam_params.txt $(SPECPATHINT)/sum_initials40GeV_1mm_full_peak_alt.txt $(SPECPATHINT)/sum_angles40GeV_1mm_full_peak_alt.txt
	# ./test_sim aligned sim_aligned2_40GeV 40GeV_beam_params.txt $(SPECPATHINT)/sum_initials40GeV_1mm_full_peak_alt.txt $(SPECPATHINT)/sum_angles40GeV_1mm_full_peak_alt.txt
	# ./test_sim aligned sim_aligned3_40GeV 40GeV_beam_params.txt $(SPECPATHINT)/sum_initials40GeV_1mm_full_peak_alt.txt $(SPECPATHINT)/sum_angles40GeV_1mm_full_peak_alt.txt
	# ./test_sim aligned sim_aligned4_40GeV 40GeV_beam_params.txt $(SPECPATHINT)/sum_initials40GeV_1mm_full_peak_alt.txt $(SPECPATHINT)/sum_angles40GeV_1mm_full_peak_alt.txt
	# ./test_sim aligned sim_aligned5_40GeV 40GeV_beam_params.txt $(SPECPATHINT)/sum_initials40GeV_1mm_full_peak_alt.txt $(SPECPATHINT)/sum_angles40GeV_1mm_full_peak_alt.txt
	#
	# ./test_sim amorphous sim_amorphous1_40GeV_1.5mm_no_background 40GeV_beam_params_1.5mm.txt
	# ./test_sim amorphous sim_amorphous2_40GeV_1.5mm_no_background 40GeV_beam_params_1.5mm.txt
	# ./test_sim amorphous sim_amorphous3_40GeV_1.5mm_no_background 40GeV_beam_params_1.5mm.txt
	# ./test_sim amorphous sim_amorphous4_40GeV_1.5mm_no_background 40GeV_beam_params_1.5mm.txt
	# ./test_sim amorphous sim_amorphous5_40GeV_1.5mm_no_background 40GeV_beam_params_1.5mm.txt

	./test_sim amorphous sim_amorphous1_40GeV_1.5mm 40GeV_beam_params_1.5mm.txt
	./test_sim amorphous sim_amorphous2_40GeV_1.5mm 40GeV_beam_params_1.5mm.txt
	./test_sim amorphous sim_amorphous3_40GeV_1.5mm 40GeV_beam_params_1.5mm.txt
	./test_sim amorphous sim_amorphous4_40GeV_1.5mm 40GeV_beam_params_1.5mm.txt
	./test_sim amorphous sim_amorphous5_40GeV_1.5mm 40GeV_beam_params_1.5mm.txt

	./test_sim background sim_background1_40GeV_1.5mm 40GeV_beam_params_1.5mm.txt
	./test_sim background sim_background2_40GeV_1.5mm 40GeV_beam_params_1.5mm.txt
	./test_sim background sim_background3_40GeV_1.5mm 40GeV_beam_params_1.5mm.txt
	./test_sim background sim_background4_40GeV_1.5mm 40GeV_beam_params_1.5mm.txt
	./test_sim background sim_background5_40GeV_1.5mm 40GeV_beam_params_1.5mm.txt

	# ./test_sim aligned sim_aligned1_40GeV_woshot_1.5mm 40GeV_beam_params_1.5mm.txt $(SPECPATHINT)/sum_initials40GeV_1.5mm_full_peak_alt_woshot.txt $(SPECPATHINT)/sum_angles40GeV_1.5mm_full_peak_alt_woshot.txt
	# ./test_sim aligned sim_aligned2_40GeV_woshot_1.5mm 40GeV_beam_params_1.5mm.txt $(SPECPATHINT)/sum_initials40GeV_1.5mm_full_peak_alt_woshot.txt $(SPECPATHINT)/sum_angles40GeV_1.5mm_full_peak_alt_woshot.txt
	# ./test_sim aligned sim_aligned3_40GeV_woshot_1.5mm 40GeV_beam_params_1.5mm.txt $(SPECPATHINT)/sum_initials40GeV_1.5mm_full_peak_alt_woshot.txt $(SPECPATHINT)/sum_angles40GeV_1.5mm_full_peak_alt_woshot.txt
	# ./test_sim aligned sim_aligned4_40GeV_woshot_1.5mm 40GeV_beam_params_1.5mm.txt $(SPECPATHINT)/sum_initials40GeV_1.5mm_full_peak_alt_woshot.txt $(SPECPATHINT)/sum_angles40GeV_1.5mm_full_peak_alt_woshot.txt
	# ./test_sim aligned sim_aligned5_40GeV_woshot_1.5mm 40GeV_beam_params_1.5mm.txt $(SPECPATHINT)/sum_initials40GeV_1.5mm_full_peak_alt_woshot.txt $(SPECPATHINT)/sum_angles40GeV_1.5mm_full_peak_alt_woshot.txt
	#
	# ./test_sim aligned sim_aligned1_40GeV_worr_1.5mm 40GeV_beam_params_1.5mm.txt $(SPECPATHINT)/sum_initials40GeV_1.5mm_full_peak_alt_worr.txt $(SPECPATHINT)/sum_angles40GeV_1.5mm_full_peak_alt_worr.txt
	# ./test_sim aligned sim_aligned2_40GeV_worr_1.5mm 40GeV_beam_params_1.5mm.txt $(SPECPATHINT)/sum_initials40GeV_1.5mm_full_peak_alt_worr.txt $(SPECPATHINT)/sum_angles40GeV_1.5mm_full_peak_alt_worr.txt
	# ./test_sim aligned sim_aligned3_40GeV_worr_1.5mm 40GeV_beam_params_1.5mm.txt $(SPECPATHINT)/sum_initials40GeV_1.5mm_full_peak_alt_worr.txt $(SPECPATHINT)/sum_angles40GeV_1.5mm_full_peak_alt_worr.txt
	# ./test_sim aligned sim_aligned4_40GeV_worr_1.5mm 40GeV_beam_params_1.5mm.txt $(SPECPATHINT)/sum_initials40GeV_1.5mm_full_peak_alt_worr.txt $(SPECPATHINT)/sum_angles40GeV_1.5mm_full_peak_alt_worr.txt
	# ./test_sim aligned sim_aligned5_40GeV_worr_1.5mm 40GeV_beam_params_1.5mm.txt $(SPECPATHINT)/sum_initials40GeV_1.5mm_full_peak_alt_worr.txt $(SPECPATHINT)/sum_angles40GeV_1.5mm_full_peak_alt_worr.txt
	# #
	# ./test_sim aligned sim_aligned1_40GeV_1.5mm 40GeV_beam_params_1.5mm.txt $(SPECPATHINT)/sum_initials40GeV_1.5mm_full_peak_alt.txt $(SPECPATHINT)/sum_angles40GeV_1.5mm_full_peak_alt.txt
	# ./test_sim aligned sim_aligned2_40GeV_1.5mm 40GeV_beam_params_1.5mm.txt $(SPECPATHINT)/sum_initials40GeV_1.5mm_full_peak_alt.txt $(SPECPATHINT)/sum_angles40GeV_1.5mm_full_peak_alt.txt
	# ./test_sim aligned sim_aligned3_40GeV_1.5mm 40GeV_beam_params_1.5mm.txt $(SPECPATHINT)/sum_initials40GeV_1.5mm_full_peak_alt.txt $(SPECPATHINT)/sum_angles40GeV_1.5mm_full_peak_alt.txt
	# ./test_sim aligned sim_aligned4_40GeV_1.5mm 40GeV_beam_params_1.5mm.txt $(SPECPATHINT)/sum_initials40GeV_1.5mm_full_peak_alt.txt $(SPECPATHINT)/sum_angles40GeV_1.5mm_full_peak_alt.txt
	# ./test_sim aligned sim_aligned5_40GeV_1.5mm 40GeV_beam_params_1.5mm.txt $(SPECPATHINT)/sum_initials40GeV_1.5mm_full_peak_alt.txt $(SPECPATHINT)/sum_angles40GeV_1.5mm_full_peak_alt.txt

	# ./test_sim alignmentrun sim_alignment alignment_beam_params.txt

test_sim : test_sim.cc simulator.cc analyser.cc

clean:
	rm -f main main.o test_sim test_sim.o
