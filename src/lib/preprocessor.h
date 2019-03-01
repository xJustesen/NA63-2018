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

class preprocessor
{
  public:
    vector<double> z;
    CONFIG Config;
    preprocessor(string ConfigFile);

  private:
    vector<string> LegalInput;
    void InitializeInputVariables(std::string filename);
    void InitializeInputVariablesHelper(std::string, std::string);
    int SearchList(std::vector<std::string> List, std::string Key);
};

#endif