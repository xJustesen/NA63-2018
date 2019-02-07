CC	:= g++
CXXFLAGS := -std=c++14 -Ofast -fopenmp -Wall -larmadillo -I/home/christian/Root/include
LDLIBS += -lm -march=native -larmadillo -fopenmp $(shell /home/christian/Documents/root/bin/root-config --libs)
DATPATH := /media/christian/Elements/CERN2018NA63
PARAMPATH := /home/christian/Dropbox/speciale/code/beamParameters
SPECPATHINT := /home/christian/Dropbox/Cern2018Experiment/spectre

data_all : data_20GeV_1mm data_20GeV_1.5mm data_40GeV_1mm data_40GeV_1.5mm data_80GeV_1mm data_80GeV_1.5mm

sim_all :  sim_40GeV_1mm sim_20GeV_1mm sim_20GeV_1.5mm#  sim_40GeV_1.5mm sim_80GeV_1mm sim_80GeV_1.5mm

data_20GeV_1mm : main_data
	# # 20 GeV ; Amorphous ; 1mm C
	./main_data $(DATPATH)/run109_01.root 109
	./main_data $(DATPATH)/run115_01.root 115

	# # 20 GeV ; Aligned ; 1mm C
	./main_data $(DATPATH)/run103_01.root 103
	./main_data $(DATPATH)/run104_01.root 104
	./main_data $(DATPATH)/run112_01.root 112
	./main_data $(DATPATH)/run113_01.root 113
	./main_data $(DATPATH)/run114_01.root 114
	#
	# # 20 GeV ; Background ; 1mm C
	./main_data $(DATPATH)/run105_01.root 105
	./main_data $(DATPATH)/run106_01.root 106
	./main_data $(DATPATH)/run107_01.root 107
	./main_data $(DATPATH)/run108_01.root 108
	./main_data $(DATPATH)/run110_01.root 110
	./main_data $(DATPATH)/run111_01.root 111

data_20GeV_1.5mm : main_data
	# 20 GeV ; Amorphous ; 1.5mm C
	./main_data $(DATPATH)/run66_01.root 66
	./main_data $(DATPATH)/run67_01.root 67
	./main_data $(DATPATH)/run68_01.root 68
	./main_data $(DATPATH)/run69_01.root 69

	# 20 GeV ; Aligned ; 1.5mm C
	./main_data $(DATPATH)/run61_01.root 61
	./main_data $(DATPATH)/run62_01.root 62
	./main_data $(DATPATH)/run63_01.root 63
	./main_data $(DATPATH)/run64_01.root 64

	# 20 GeV ; Background ; 1.5mm C
	./main_data $(DATPATH)/run60_01.root 60
	./main_data $(DATPATH)/run65_01.root 65

data_40GeV_1mm : main_data
	# 40 GeV ; Aligned ; 1mm C
	./main_data $(DATPATH)/run71_01.root 71
	./main_data $(DATPATH)/run72_01.root 72
	./main_data $(DATPATH)/run79_01.root 79
	./main_data $(DATPATH)/run80_01.root 80
	./main_data $(DATPATH)/run81_01.root 81

	# 40 GeV ; Amorphous ; 1mm C
	./main_data $(DATPATH)/run73_01.root 73
	./main_data $(DATPATH)/run75_01.root 75
	./main_data $(DATPATH)/run76_01.root 76
	./main_data $(DATPATH)/run77_01.root 77
	./main_data $(DATPATH)/run78_01.root 78

	# 40 GeV ; Background ; 1mm C
	./main_data $(DATPATH)/run74_01.root 74

data_40GeV_1.5mm : main_data
	# 40 GeV ; Aligned ; 1.5mm C
	./main_data $(DATPATH)/run30_01.root 30
	./main_data $(DATPATH)/run35_01.root 35
	./main_data $(DATPATH)/run36_01.root 36
	./main_data $(DATPATH)/run37_01.root 37
	./main_data $(DATPATH)/run38_01.root 38

	# 40 GeV ; Amorphous ; 1.5mm C
	./main_data $(DATPATH)/run32_01.root 32
	./main_data $(DATPATH)/run33_01.root 33
	./main_data $(DATPATH)/run34_01.root 34
	./main_data $(DATPATH)/run39_01.root 39
	./main_data $(DATPATH)/run40_01.root 40
	./main_data $(DATPATH)/run41_01.root 41
	./main_data $(DATPATH)/run43_01.root 43

	# 40 GeV ; Background ; 1.5mm C
	./main_data $(DATPATH)/run31_01.root 31

data_80GeV_1mm : main_data
	# 80 GeV ; Amorphous ; 1mm C
	./main_data $(DATPATH)/run86_01.root 86
	./main_data $(DATPATH)/run87_01.root 87
	./main_data $(DATPATH)/run88_01.root 88
	./main_data $(DATPATH)/run89_01.root 89

	# 80 GeV ; Background ; 1mm C
	./main_data $(DATPATH)/run85_01.root 85
	./main_data $(DATPATH)/run90_01.root 90
	./main_data $(DATPATH)/run91_01.root 91

	# 80 GeV ; Aligned ; 1mm C
	./main_data $(DATPATH)/run84_01.root 84
	./main_data $(DATPATH)/run92_01.root 92
	./main_data $(DATPATH)/run93_01.root 93
	./main_data $(DATPATH)/run94_01.root 94
	./main_data $(DATPATH)/run95_01.root 95

data_80GeV_1.5mm : main_data
	# 80 GeV ; Aligned ; 1.5mm C
	./main_data $(DATPATH)/run46_01.root 46
	./main_data $(DATPATH)/run47_01.root 47

	# 80 GeV ; Amorphous ; 1.5mm C
	./main_data $(DATPATH)/run49_01.root 49
	./main_data $(DATPATH)/run50_01.root 50
	./main_data $(DATPATH)/run51_01.root 51
	./main_data $(DATPATH)/run52_01.root 52
	./main_data $(DATPATH)/run53_01.root 53
	./main_data $(DATPATH)/run56_01.root 56
	./main_data $(DATPATH)/run57_01.root 57

	# 80 GeV ; Background ; 1.5mm C
	./main_data $(DATPATH)/run48_01.root 48
	./main_data $(DATPATH)/run54_01.root 54
	./main_data $(DATPATH)/run55_01.root 55

data_2017 : main_data
	# ./main_data /media/christian/Elements/dataCERN2017v2/tafconverted/run40/run44_01.root 2017_50GeV_2mm # aligned 2mm Si 50 GeV
	# ./main_data /media/christian/Elements/dataCERN2017v2/tafconverted/run32/run44_01.root 2017_50GeV_2mm_bg # bg 2mm Si 50 GeV
	./main_data /media/christian/Elements/dataCERN2017v2/tafconverted/run47/run44_01.root 2017_50GeV_2mm_amorph # amorph 2mm Si 50 GeV

data_alignmnet : main_data
	# # Alignment
	./main_data $(DATPATH)/run18_01.root 18

data_fakehits : main_data
	# Fake hits run
	./main_data $(DATPATH)/run4_01.root 4

main_data : main_data.cc analyser.cc

n = 4
aligned_40GeV : main_sim
	n=$(n); \
	while [ $${n} -gt 0 ] ; do \
		echo $$n ; \
		./main_sim aligned sim_aligned_40GeV_worr_test 40GeV_beam_params.txt $(SPECPATHINT)/sum_angles1.0mm40GeVbaiernoRR.txt ; \
		./main_sim aligned sim_aligned_40GeV_woshot_test 40GeV_beam_params.txt $(SPECPATHINT)/sum_angles1.0mm40GeVbaiernoSHbfac.txt ; \
		./main_sim aligned sim_aligned_40GeV_test 40GeV_beam_params.txt $(SPECPATHINT)/sum_angles1.0mm40GeVbaierRRbfac.txt; \
		n=`expr $$n - 1`; \
	done; \
	true

amorph_sim_1mm : main_sim
	n=$(n); \
	while [ $${n} -gt 0 ] ; do \
		echo $$n ; \
		./main_sim amorphous sim_amorphous_20GeV_test 20GeV_beam_params.txt ; \
		./main_sim background sim_background_20GeV_test 20GeV_beam_params.txt ; \
		./main_sim amorphous sim_amorphous_40GeV_test 40GeV_beam_params.txt ; \
		./main_sim background sim_background_40GeV_test 40GeV_beam_params.txt ; \
		./main_sim amorphous sim_amorphous_80GeV_test 80GeV_beam_params.txt ; \
		./main_sim background sim_background_80GeV_test 80GeV_beam_params.txt ; \
		n=`expr $$n - 1`; \
	done; \
	true

aligned_20GeV : main_sim
	n=$(n); \
	while [ $${n} -gt 0 ] ; do \
		echo $$n ; \
		./main_sim aligned sim_aligned_20GeV_worr_test 20GeV_beam_params.txt $(SPECPATHINT)/sum_angles1.0mm20GeVbaiernoRR.txt ; \
		./main_sim aligned sim_aligned_20GeV_woshot_test 20GeV_beam_params.txt $(SPECPATHINT)/sum_angles1.0mm20GeVbaiernoSHbfac.txt ; \
		./main_sim aligned sim_aligned_20GeV_test 20GeV_beam_params.txt $(SPECPATHINT)/sum_angles1.0mm20GeVbaierRRbfac.txt; \
		n=`expr $$n - 1`; \
	done; \
	true

aligned_80GeV : main_sim
	n=$(n); \
	while [ $${n} -gt 0 ] ; do \
		echo $$n ; \
		./main_sim aligned sim_aligned_80GeV_worr_test 80GeV_beam_params.txt $(SPECPATHINT)/sum_angles1.0mm80GeVbaiernoRR.txt ; \
		./main_sim aligned sim_aligned_80GeV_woshot_test 80GeV_beam_params.txt $(SPECPATHINT)/sum_angles1.0mm80GeVbaiernoSHbfac.txt ; \
		./main_sim aligned sim_aligned_80GeV_test 80GeV_beam_params.txt $(SPECPATHINT)/sum_angles1.0mm80GeVbaierRRbfac.txt; \
		n=`expr $$n - 1`; \
	done; \
	true

amorph_20GeV : main_sim
	n=$(n); \
	while [ $${n} -gt 0 ] ; do \
		echo $$n ; \
		./main_sim amorphous sim_amorphous_20GeV_test 20GeV_beam_params.txt ; \
		./main_sim background sim_background_20GeV_test 20GeV_beam_params.txt ; \
		n=`expr $$n - 1`; \
	done; \
	true

amorph_40GeV : main_sim
	n=$(n); \
	while [ $${n} -gt 0 ] ; do \
		echo $$n ; \
		./main_sim amorphous sim_amorphous_40GeV_test 40GeV_beam_params.txt ; \
		./main_sim background sim_background_40GeV_test 40GeV_beam_params.txt ; \
		n=`expr $$n - 1`; \
	done; \
	true

amorph_80GeV : main_sim
	n=$(n); \
	while [ $${n} -gt 0 ] ; do \
		echo $$n ; \
		./main_sim amorphous sim_amorphous_80GeV_test 80GeV_beam_params.txt ; \
		./main_sim background sim_background_80GeV_test 80GeV_beam_params.txt ; \
		n=`expr $$n - 1`; \
	done; \
	true

sim_alignment : main_sim
	./main_sim alignmentrun sim_alignment alignment_beam_params.txt

sim_monochrome : main_sim
	./main_sim amorphous sim_monochrome_80GeV_1.5mm_25GeV 80GeV_beam_params_1.5mm.txt

sim_2017 : main_sim
	./main_sim amorphous sim_amorphous_50GeV_622mm flohr_beam.txt

main_sim : main_sim.cc simulator.cc analyser.cc

clean:
	rm -f main_sim main_sim.o main_data main_data.o
