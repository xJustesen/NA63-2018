#ifndef RUNDATAANALYSIS_H
#define RUNDATAANALYSIS_H

#include "analyser.h"
#include "preprocessor.h"

using namespace std;

class RunDataAnalysis
{
  public:
    RunDataAnalysis(preprocessor *initialize,  const char *InputDataFile);
    vector<vector<double>> get_energies();

  private:
    vector<vector<double>> energies;
};

#endif