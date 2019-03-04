#ifndef RUNDATAANALYSIS_H
#define RUNDATAANALYSIS_H

#include "analyser.h"
#include "pre_processor.h"

using namespace std;

class RunDataAnalysis
{
  public:
    RunDataAnalysis(PreProcessor *initialize,  const char *input_data_file);
    vector<vector<double>> GetEnergies();

  private:
    vector<vector<double>> energies_;
};

#endif