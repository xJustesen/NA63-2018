#ifndef PREPROCESSOR_H
#define PREPROCESSOR_H

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include "auxillary_functions.h"

using namespace std;

struct CONFIG
{
  double BeamEnergy, CrystalThickness, cut_lb_x, cut_ub_x, cut_lb_y, cut_ub_y;
  std::string runno, BeamProfile, CrystalType, DataPath, TheorySpectrum, OutputDataFilename, OutputEventFilename, InputDatafile;
  int IncludeBG, NEvents, Simulation;
};

class PreProcessor
{
public:
  vector<double> z_;
  CONFIG config_;

  PreProcessor(string config_file);

private:
  vector<string> legal_input_;
  void InitializeInputVariables(std::string file_name);
  void InitializeInputVariablesHelper(std::string, std::string);
  int SearchLegalInput(std::string key);
  void LoadConfigFile(const string &configfile);
};

#endif