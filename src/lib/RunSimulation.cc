#include "RunSimulation.h"

RunSimulation::RunSimulation(preprocessor *initialize)
{
    simulator *_simulate = new simulator(initialize->Config.NEvents,
                                         initialize->z,
                                         initialize->Config.CrystalType,
                                         initialize->Config.BeamEnergy,
                                         initialize->Config.CrystalThickness,
                                         initialize->Config.BeamProfile,
                                         initialize->Config.TheorySpectrum,
                                         initialize->Config.IncludeBG);
    _simulate->propagate_particles();
    cout << "\nData simulated, starting analysis\n";

    /* Analysis */
    const char *DummyROOTFilename = "";
    analyser *_analyse = new analyser(initialize->z,
                                      DummyROOTFilename,
                                      initialize->Config.runno,
                                      initialize->Config.BeamProfile);
    _analyse->hitcoords = _simulate->mimosas;
    delete _simulate;
    _analyse->Nevents = initialize->Config.NEvents;
    _analyse->construct_tracks(initialize->Config.cut_lb_x,
                               initialize->Config.cut_ub_x,
                               initialize->Config.cut_lb_y,
                               initialize->Config.cut_ub_y);
    _analyse->hitcoords.clear();         // clears elements from vector (size = 0, capacity unchanged)
    _analyse->hitcoords.shrink_to_fit(); // frees memory, capacity = 0
    cout << "Pairing tracks and calculating energies\n";
    _analyse->pair_tracks();
    energies = _analyse->get_energies();
}

vector<vector<double>> RunSimulation::get_energies()
{
    return energies;
}
