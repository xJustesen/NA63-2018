#ifndef PREPROCESSOR_H
#define PREPROCESSOR_H

#include <string>
#include <vector>
#include <iostream>
#include <fstream>

using namespace std;

typedef struct
{
    double BeamEnergy, CrystalThickness, cut_lb_x, cut_ub_x, cut_lb_y, cut_ub_y;
    std::string runno, BeamProfile, CrystalType, DataPath, TheorySpectrum, OutputFilename, InputDatafile;
    int IncludeBG, NEvents, Simulation;
} CONFIG;

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
    int SearchList(std::vector<std::string> list, std::string key);
};

#endif