CC	:= g++
CXXFLAGS = -I -std=c++14 -Ofast -fopenmp -Wall -larmadillo -I/home/christian/Root/include
LDLIBS = -lm -march=native -larmadillo -fopenmp $(shell /home/christian/Documents/root/bin/root-config --libs)
DATPATH := /media/christian/Elements/CERN2018NA63
DATPATH17 := /media/christian/Elements/dataCERN2017v2/tafconverted
PARAMPATH := /home/christian/Dropbox/speciale/code/beamParameters
SPECPATHINT := /home/christian/Dropbox/Cern2018Experiment/spectre

data_all : data_80GeV_1mm data_80GeV_1.5mm #data_20GeV_1mm data_20GeV_1.5mm data_40GeV_1mm data_40GeV_1.5mm

sim_all :  sim_40GeV_1mm sim_20GeV_1mm sim_20GeV_1.5mm#  sim_40GeV_1.5mm sim_80GeV_1mm sim_80GeV_1.5mm

test : test.o
	$(CC) $(CXXFLAGS) -o test test.cc $(LDLIBS)

data_20GeV_1mm : main_data
	# 20 GeV ; Amorphous ; 1mm C
	# ./main_data $(DATPATH)/run109_01.root 109 $(PARAMPATH)/20GeV_1.0mmC_cuts.txt
	# ./main_data $(DATPATH)/run115_01.root 115 $(PARAMPATH)/20GeV_1.0mmC_cuts.txt

	# 20 GeV ; Aligned ; 1mm C
	./main_data $(DATPATH)/run103_01.root 103 $(PARAMPATH)/20GeV_1.0mmC_cuts.txt
	# ./main_data $(DATPATH)/run104_01.root 104 $(PARAMPATH)/20GeV_1.0mmC_cuts.txt
	# ./main_data $(DATPATH)/run112_01.root 112 $(PARAMPATH)/20GeV_1.0mmC_cuts.txt
	# ./main_data $(DATPATH)/run113_01.root 113 $(PARAMPATH)/20GeV_1.0mmC_cuts.txt
	# ./main_data $(DATPATH)/run114_01.root 114 $(PARAMPATH)/20GeV_1.0mmC_cuts.txt
	#
	# 20 GeV ; Background ; 1mm C
	# ./main_data $(DATPATH)/run105_01.root 105 $(PARAMPATH)/20GeV_1.0mmC_cuts.txt
	# ./main_data $(DATPATH)/run106_01.root 106 $(PARAMPATH)/20GeV_1.0mmC_cuts.txt
	# ./main_data $(DATPATH)/run107_01.root 107 $(PARAMPATH)/20GeV_1.0mmC_cuts.txt
	# ./main_data $(DATPATH)/run108_01.root 108 $(PARAMPATH)/20GeV_1.0mmC_cuts.txt
	# ./main_data $(DATPATH)/run110_01.root 110 $(PARAMPATH)/20GeV_1.0mmC_cuts.txt
	# ./main_data $(DATPATH)/run111_01.root 111 $(PARAMPATH)/20GeV_1.0mmC_cuts.txt

data_20GeV_1.5mm : main_data
	# 20 GeV ; Amorphous ; 1.5mm C
	# ./main_data $(DATPATH)/run66_01.root 66 $(PARAMPATH)/20GeV_1.5mmC_cuts.txt
	# ./main_data $(DATPATH)/run67_01.root 67 $(PARAMPATH)/20GeV_1.5mmC_cuts.txt
	# ./main_data $(DATPATH)/run68_01.root 68 $(PARAMPATH)/20GeV_1.5mmC_cuts.txt
	# ./main_data $(DATPATH)/run69_01.root 69 $(PARAMPATH)/20GeV_1.5mmC_cuts.txt

	# 20 GeV ; Aligned ; 1.5mm C
	# ./main_data $(DATPATH)/run61_01.root 61 $(PARAMPATH)/20GeV_1.5mmC_cuts.txt
	./main_data $(DATPATH)/run62_01.root 62 $(PARAMPATH)/20GeV_1.5mmC_cuts.txt
	./main_data $(DATPATH)/run63_01.root 63 $(PARAMPATH)/20GeV_1.5mmC_cuts.txt
	./main_data $(DATPATH)/run64_01.root 64 $(PARAMPATH)/20GeV_1.5mmC_cuts.txt

	# 20 GeV ; Background ; 1.5mm C
	./main_data $(DATPATH)/run60_01.root 60 $(PARAMPATH)/20GeV_1.5mmC_cuts.txt
	./main_data $(DATPATH)/run65_01.root 65 $(PARAMPATH)/20GeV_1.5mmC_cuts.txt

data_40GeV_1mm : main_data
	# 40 GeV ; Aligned ; 1mm C
	./main_data $(DATPATH)/run71_01.root 71 $(PARAMPATH)/40GeV_1.0mmC_cuts.txt
	./main_data $(DATPATH)/run72_01.root 72 $(PARAMPATH)/40GeV_1.0mmC_cuts.txt
	./main_data $(DATPATH)/run79_01.root 79 $(PARAMPATH)/40GeV_1.0mmC_cuts.txt
	./main_data $(DATPATH)/run80_01.root 80 $(PARAMPATH)/40GeV_1.0mmC_cuts.txt
	./main_data $(DATPATH)/run81_01.root 81 $(PARAMPATH)/40GeV_1.0mmC_cuts.txt

	# 40 GeV ; Amorphous ; 1mm C
	./main_data $(DATPATH)/run73_01.root 73 $(PARAMPATH)/40GeV_1.0mmC_cuts.txt
	./main_data $(DATPATH)/run75_01.root 75 $(PARAMPATH)/40GeV_1.0mmC_cuts.txt
	./main_data $(DATPATH)/run76_01.root 76 $(PARAMPATH)/40GeV_1.0mmC_cuts.txt
	./main_data $(DATPATH)/run77_01.root 77 $(PARAMPATH)/40GeV_1.0mmC_cuts.txt
	./main_data $(DATPATH)/run78_01.root 78 $(PARAMPATH)/40GeV_1.0mmC_cuts.txt

	# 40 GeV ; Background ; 1mm C
	./main_data $(DATPATH)/run74_01.root 74 $(PARAMPATH)/40GeV_1.0mmC_cuts.txt

data_40GeV_1.5mm : main_data
	# 40 GeV ; Aligned ; 1.5mm C
	./main_data $(DATPATH)/run30_01.root 30 $(PARAMPATH)/40GeV_1.5mmC_cuts.txt
	./main_data $(DATPATH)/run35_01.root 35 $(PARAMPATH)/40GeV_1.5mmC_cuts.txt
	./main_data $(DATPATH)/run36_01.root 36 $(PARAMPATH)/40GeV_1.5mmC_cuts.txt
	./main_data $(DATPATH)/run37_01.root 37 $(PARAMPATH)/40GeV_1.5mmC_cuts.txt
	./main_data $(DATPATH)/run38_01.root 38 $(PARAMPATH)/40GeV_1.5mmC_cuts.txt

	# 40 GeV ; Amorphous ; 1.5mm C
	./main_data $(DATPATH)/run32_01.root 32 $(PARAMPATH)/40GeV_1.5mmC_cuts.txt
	./main_data $(DATPATH)/run33_01.root 33 $(PARAMPATH)/40GeV_1.5mmC_cuts.txt
	./main_data $(DATPATH)/run34_01.root 34 $(PARAMPATH)/40GeV_1.5mmC_cuts.txt
	./main_data $(DATPATH)/run39_01.root 39 $(PARAMPATH)/40GeV_1.5mmC_cuts.txt
	./main_data $(DATPATH)/run40_01.root 40 $(PARAMPATH)/40GeV_1.5mmC_cuts.txt
	./main_data $(DATPATH)/run41_01.root 41 $(PARAMPATH)/40GeV_1.5mmC_cuts.txt
	./main_data $(DATPATH)/run43_01.root 43 $(PARAMPATH)/40GeV_1.5mmC_cuts.txt

	# 40 GeV ; Background ; 1.5mm C
	./main_data $(DATPATH)/run31_01.root 31 $(PARAMPATH)/40GeV_1.5mmC_cuts.txt

data_80GeV_1mm : main_data
	# 80 GeV ; Amorphous ; 1mm C
	# ./main_data $(DATPATH)/run86_01.root 86 $(PARAMPATH)/80GeV_1.0mmC_cuts.txt
	# ./main_data $(DATPATH)/run87_01.root 87 $(PARAMPATH)/80GeV_1.0mmC_cuts.txt
	# ./main_data $(DATPATH)/run88_01.root 88 $(PARAMPATH)/80GeV_1.0mmC_cuts.txt
	# ./main_data $(DATPATH)/run89_01.root 89 $(PARAMPATH)/80GeV_1.0mmC_cuts.txt

	# 80 GeV ; Background ; 1mm C
	./main_data $(DATPATH)/run85_01.root 85 $(PARAMPATH)/80GeV_1.0mmC_cuts.txt
	./main_data $(DATPATH)/run90_01.root 90 $(PARAMPATH)/80GeV_1.0mmC_cuts.txt
	./main_data $(DATPATH)/run91_01.root 91 $(PARAMPATH)/80GeV_1.0mmC_cuts.txt

	# 80 GeV ; Aligned ; 1mm C
	# ./main_data $(DATPATH)/run84_01.root 84 $(PARAMPATH)/80GeV_1.0mmC_cuts.txt
	# ./main_data $(DATPATH)/run92_01.root 92 $(PARAMPATH)/80GeV_1.0mmC_cuts.txt
	# ./main_data $(DATPATH)/run93_01.root 93 $(PARAMPATH)/80GeV_1.0mmC_cuts.txt
	# ./main_data $(DATPATH)/run94_01.root 94 $(PARAMPATH)/80GeV_1.0mmC_cuts.txt
	# ./main_data $(DATPATH)/run95_01.root 95 $(PARAMPATH)/80GeV_1.0mmC_cuts.txt

data_80GeV_1.5mm : main_data
	# 80 GeV ; Aligned ; 1.5mm C
	# ./main_data $(DATPATH)/run46_01.root 46 $(PARAMPATH)/80GeV_1.5mmC_cuts.txt
	# ./main_data $(DATPATH)/run47_01.root 47 $(PARAMPATH)/80GeV_1.5mmC_cuts.txt
	#
	# # 80 GeV ; Amorphous ; 1.5mm C
	# ./main_data $(DATPATH)/run49_01.root 49 $(PARAMPATH)/80GeV_1.5mmC_cuts.txt
	# ./main_data $(DATPATH)/run50_01.root 50 $(PARAMPATH)/80GeV_1.5mmC_cuts.txt
	# ./main_data $(DATPATH)/run51_01.root 51 $(PARAMPATH)/80GeV_1.5mmC_cuts.txt
	# ./main_data $(DATPATH)/run52_01.root 52 $(PARAMPATH)/80GeV_1.5mmC_cuts.txt
	# ./main_data $(DATPATH)/run53_01.root 53 $(PARAMPATH)/80GeV_1.5mmC_cuts.txt
	# ./main_data $(DATPATH)/run56_01.root 56 $(PARAMPATH)/80GeV_1.5mmC_cuts.txt
	# ./main_data $(DATPATH)/run57_01.root 57 $(PARAMPATH)/80GeV_1.5mmC_cuts.txt
	#
	# # 80 GeV ; Background ; 1.5mm C
	# ./main_data $(DATPATH)/run48_01.root 48 $(PARAMPATH)/80GeV_1.5mmC_cuts.txt
	# ./main_data $(DATPATH)/run54_01.root 54 $(PARAMPATH)/80GeV_1.5mmC_cuts.txt
	# ./main_data $(DATPATH)/run55_01.root 55 $(PARAMPATH)/80GeV_1.5mmC_cuts.txt

data_2017 : main_data
	# ./main_data /media/christian/Elements/dataCERN2017v2/tafconverted/run40/run44_01.root 2017_50GeV_2mm # aligned 2mm Si 50 GeV
	# ./main_data /media/christian/Elements/dataCERN2017v2/tafconverted/run32/run44_01.root 2017_50GeV_2mm_bg # bg 2mm Si 50 GeV
	# ./main_data /media/christian/Elements/dataCERN2017v2/tafconverted/run47/run44_01.root 2017_50GeV_2mm_amorph # amorph 2mm Si 50 GeV

data_alignmnet : main_data
	# # Alignment
	./main_data $(DATPATH)/run18_01.root 18

data_fakehits : main_data
	# Fake hits run
	./main_data $(DATPATH)/run4_01.root 4

data_findaxis : main_data
	# ./main_data $(DATPATH)/run35_01.root 35 $(PARAMPATH)/40GeV_beam_params_1.5mm.txt
	# ./main_data $(DATPATH)/run32_01.root 32 $(PARAMPATH)/40GeV_beam_params_1.5mm.txt
	./main_data $(DATPATH17)/run111/run44_01.root 111_2017 $(PARAMPATH)/2017_50GeV_4mm.txt # aligned 6mm Si 50 GeV

main_data : main_data.cc analyser.cc

n = 5
aligned_20GeV_1_5mm : main_sim
	n=$(n); \
	while [ $${n} -gt 0 ] ; do \
		echo $$n ; \
		./main_sim aligned sim_aligned_20GeV_1.5mm_test 20GeV_beam_params_1.5mm.txt $(SPECPATHINT)/sum_angles20GeV1.5mmRR.txt; \
		n=`expr $$n - 1`; \
	done; \
	true

aligned_40GeV_1_5mm : main_sim
	n=$(n); \
	while [ $${n} -gt 0 ] ; do \
		echo $$n ; \
		./main_sim aligned sim_aligned_40GeV_1.5mm_test 40GeV_beam_params_1.5mm.txt $(SPECPATHINT)/sum_angles40GeV1.5mmRR.txt; \
		n=`expr $$n - 1`; \
	done; \
	true

aligned_80GeV_1_5mm : main_sim
	n=$(n); \
	while [ $${n} -gt 0 ] ; do \
		echo $$n ; \
		./main_sim aligned sim_aligned_80GeV_1.5mm_test 80GeV_beam_params_1.5mm.txt $(SPECPATHINT)/sum_angles80GeV1.5mmRR.txt; \
		n=`expr $$n - 1`; \
	done; \
	true

aligned_20GeV : main_sim
	n=$(n); \
	while [ $${n} -gt 0 ] ; do \
		echo $$n ; \
		./main_sim aligned sim_aligned_20GeV 20GeV_beam_params.txt $(SPECPATHINT)/cuts/sum_angles20GeV1mmRRcut.txt $(PARAMPATH)/20GeV_1.0mmC_cuts.txt; \
		n=`expr $$n - 1`; \
	done; \
	true

aligned_40GeV : main_sim
	n=$(n); \
	while [ $${n} -gt 0 ] ; do \
		echo $$n ; \
		./main_sim aligned sim_aligned_40GeV_worr_test 40GeV_beam_params.txt $(SPECPATHINT)/sum_angles40GeV1mmnoRR.txt ; \
		./main_sim aligned sim_aligned_40GeV_woshot_test 40GeV_beam_params.txt $(SPECPATHINT)/sum_angles40GeV1mmnoSH.txt ; \
		./main_sim aligned sim_aligned_40GeV_test 40GeV_beam_params.txt $(SPECPATHINT)/sum_angles40GeV1mmRR.txt; \
		n=`expr $$n - 1`; \
	done; \
	true

aligned_80GeV : main_sim
	n=$(n); \
	while [ $${n} -gt 0 ] ; do \
		echo $$n ; \
		./main_sim aligned sim_aligned_80GeV_worr_test 80GeV_beam_params.txt $(SPECPATHINT)/sum_angles80GeV1mmnoRR.txt ; \
		./main_sim aligned sim_aligned_80GeV_woshot_test 80GeV_beam_params.txt $(SPECPATHINT)/sum_angles80GeV1mmnoSH.txt ; \
		./main_sim aligned sim_aligned_80GeV_test 80GeV_beam_params.txt $(SPECPATHINT)/sum_angles80GeV1mmRR.txt; \
		n=`expr $$n - 1`; \
	done; \
	true

amorph_20GeV : main_sim
	n=$(n); \
	while [ $${n} -gt 0 ] ; do \
		echo $$n ; \
		./main_sim amorphous sim_amorphous_20GeV 20GeV_beam_params.txt $(SPECPATHINT)/cuts/sum_angles20GeV1mmRRcut.txt $(PARAMPATH)/20GeV_1.0mmC_cuts.txt; \
		./main_sim background sim_background_20GeV 20GeV_beam_params.txt $(SPECPATHINT)/cuts/sum_angles20GeV1mmRRcut.txt $(PARAMPATH)/20GeV_1.0mmC_cuts.txt; \
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
	rm -f main_sim main_sim.o main_data main_data.o tester tester.o
