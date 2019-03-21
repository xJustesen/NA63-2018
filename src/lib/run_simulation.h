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
  int GetEvents();

private:
  vector<vector<double>> energies_;
  int number_of_events_;
};

#endif
