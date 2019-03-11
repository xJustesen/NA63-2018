#include "run_simulation.h"

RunSimulation::RunSimulation(PreProcessor *initialize)
{
    simulator *Simulate = new simulator(initialize->config_.NEvents,
                                        initialize->z_,
                                        initialize->config_.CrystalType,
                                        initialize->config_.BeamEnergy,
                                        initialize->config_.CrystalThickness,
                                        initialize->config_.BeamProfile,
                                        initialize->config_.TheorySpectrum,
                                        initialize->config_.IncludeBG);
    Simulate->PropagateParticles();
    cout << "\nData simulated, starting analysis\n";
    /* Analysis */
    const char *dummy_root_file_name = "";
    analyser *Analyse = new analyser(initialize->z_,
                                     dummy_root_file_name,
                                     initialize->config_.runno,
                                     initialize->config_.BeamProfile);
    Analyse->hits_container_ = Simulate->mimosas_;
    delete Simulate;
    Analyse->total_events_ = initialize->config_.NEvents;
    Analyse->ConstructTracks(initialize->config_.cut_lb_x,
                             initialize->config_.cut_ub_x,
                             initialize->config_.cut_lb_y,
                             initialize->config_.cut_ub_y);
    Analyse->hits_container_.clear();         // clears elements from vector (size = 0, capacity unchanged)
    Analyse->hits_container_.shrink_to_fit(); // frees memory, capacity = 0
    cout << "Pairing tracks and calculating energies\n";
    Analyse->PairTracks();
    energies_ = Analyse->GetEnergies();
}

vector<vector<double>> RunSimulation::GetEnergies()
{
    return energies_;
}
