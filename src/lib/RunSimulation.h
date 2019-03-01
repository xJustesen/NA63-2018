#ifndef RUNSIMULATION_H
#define RUNSIMULATION_H

#include "simulator.h"
#include "analyser.h"
#include "preprocessor.h"

using namespace std;

class RunSimulation
{
  public:
    RunSimulation(preprocessor *initialize);
    vector<vector<double>> get_energies();

  private:
    vector<vector<double>> energies;
};

#endif

