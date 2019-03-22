#include "run_simulation.h"

// Run simulation.
// Input:
//      initialize  Instance of PreProcessor class
RunSimulation::RunSimulation(PreProcessor *initialize)
{
    // Initalize simulator class
    simulator *Simulate = new simulator(initialize->config_.NEvents,
                                        initialize->z_,
                                        initialize->config_.CrystalType,
                                        initialize->config_.BeamEnergy,
                                        initialize->config_.CrystalThickness,
                                        initialize->config_.BeamProfile,
                                        initialize->config_.TheorySpectrum,
                                        initialize->config_.IncludeBG);
    // Generate simulated data set
    Simulate->GenerateSyntheticData();
    cout << "\nData simulated, starting analysis\n";
    // Initialize analyser class
    const char *dummy_root_file_name = "";
    analyser *Analyse = new analyser(initialize->z_,
                                     dummy_root_file_name,
                                     initialize->config_.runno,
                                     initialize->config_.BeamProfile);
    Analyse->hits_container_ = Simulate->mimosas_;
    delete Simulate; // delete since no longer needed, everything lives in analyser class now
    Analyse->total_events_ = initialize->config_.NEvents;
    // Analyse data set
    Analyse->ConstructTracks(initialize->config_.cut_lb_x,
                             initialize->config_.cut_ub_x,
                             initialize->config_.cut_lb_y,
                             initialize->config_.cut_ub_y);
    Analyse->hits_container_.clear();         // clears elements from vector (size = 0, capacity unchanged)
    Analyse->hits_container_.shrink_to_fit(); // frees memory, capacity = 0
    cout << "Pairing tracks and calculating energies\n";
    Analyse->PairTracks();
    // Retrieve energies and events
    energies_ = Analyse->GetEnergies();
    number_of_events_ = Analyse->GetEventsWithinCut();
}

// Return energies
vector<vector<double>> RunSimulation::GetEnergies()
{
    return energies_;
}

// Returns events
int RunSimulation::GetEvents()
{
    return number_of_events_;
}
