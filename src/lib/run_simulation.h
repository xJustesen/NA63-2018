#ifndef RUNSIMULATION_H
#define RUNSIMULATION_H

#include "simulator.h"
#include "analyser.h"
#include "pre_processor.h"

using namespace std;

class RunSimulation
{
  public:
    RunSimulation(PreProcessor *initialize);
    vector<vector<double>> GetEnergies();

  private:
    vector<vector<double>> energies_;
};

#endif

