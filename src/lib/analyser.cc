#include "analyser.h"

// Constructor for analyser-class. "analysis_configuration.json" is read and member parameters assigned corresponding
// values.
// Input:
//      z           : container with z-coordinates of mimosa detectors
//      name        : name of file which containts ROOT data.
//      run         : identifier for which type run we analyse (fx. amorphous_40GeV_1mm).
//      beamparams  : File containing beam parameters (direction and size).
analyser::analyser(vector<double> z, const char *name, string run, string beamparams)
{
    LoadConfigFile("/home/christian/Dropbox/speciale/code/src/config/analysis_configuration.json");
    M1_M2_zcoord_ = {detector_zcoord_[0], detector_zcoord_[1], detector_zcoord_[2]};                               // z-coords for M1-M3 PLANES
    M2_M3_zcoord_ = {detector_zcoord_[1], detector_zcoord_[2], detector_zcoord_[3]};                               // z-coords for M2-M4 planes
    M3_M4_zcoord_ = {detector_zcoord_[2], detector_zcoord_[3], (detector_zcoord_[4] + detector_zcoord_[3]) / 2.0}; // z-coords of M3-MM planes
    M6_M5_zcoord_ = {detector_zcoord_[5], detector_zcoord_[4], (detector_zcoord_[4] + detector_zcoord_[3]) / 2.0}; // z-coords for M5-MM planes
    run_number_ = run;
    file_name_ = name;
    beam_parameters_file_name_ = beamparams;
    hits_container_.resize(6);
    pixel_grids_.resize(6);
    hot_pixels_.resize(6);
    T.set_size(3, 3, 4);
    interplanar_distance_.resize(4);
    cout.setf(ios::fixed, ios::floatfield);
    cout.precision(3);
}
// Load configuration from configfile and assign class member data-structures the corresponding values
// Input:
//          kconfigfile : name of config file including path and .JSON extension
void analyser::LoadConfigFile(const string &kconfigfile)
{
    pt::ptree tree;                   // Create empty property tree object
    pt::read_json(kconfigfile, tree); // Parse the JSON into the property tree.
    // Assign parameters
    row_count_ = tree.get<int>("Detector.NumberOfRows");
    column_count_ = tree.get<int>("Detector.NumberOfColumns");
    detector_xcoord_min_ = tree.get<double>("Detector.XMin");
    detector_xcoord_max_ = tree.get<double>("Detector.XMax");
    detector_ycoord_min_ = tree.get<double>("Detector.YMin");
    detector_ycoord_max_ = tree.get<double>("Detector.YMax");
    M1_M2_proj_lim_ = tree.get<double>("ProjectionLimit.M1M2");
    M2_M3_proj_lim_ = tree.get<double>("ProjectionLimit.M2M3");
    M6_M5_proj_lim_ = tree.get<double>("ProjectionLimit.M6M5");
    foil_paired_tracks_lim_ = tree.get<double>("ProjectionLimit.Foil");
    T_matrix_convergence_tol_ = tree.get<double>("MatrixConvergenceTolerance");
    data_path_ = tree.get<string>("DataPath");
    for (auto i : as_vector<double>(tree, "AlignmentRadius"))
    {
        alignment_radius_lim_.push_back(i);
    }
    for (auto i : as_vector<double>(tree, "DetectorZCoordinates"))
    {
        detector_zcoord_.push_back(i);
    }
}

vector<vector<double>> analyser::GetEnergies()
{
    return energies;
}

int analyser::GetEventsWithinCut()
{
    return events_within_cut_;
}

void analyser::UpdatePixelgrids(int plane, vector<vector<double>> pixel_grid)
{
    pixel_grids_[plane] = pixel_grid;
}

void analyser::UpdateHitCoordinates(int plane, vector<vector<vector<double>>> hit_coordinates)
{
    hits_container_[plane] = hit_coordinates;
}

void analyser::UpdateHotPixels(int plane, vector<int> hot_pixel)
{
    hot_pixels_[plane] = hot_pixel;
}

// Append data to file name in "data_path_"
//  Input:
//          name : string with name of file (without extensions and path)
//          data : container with data to save, 1 column
void analyser::PrintVector(string name, vector<double> data)
{
    string file_name = data_path_ + "/" + name + ".txt";
    ofstream output(file_name, ios_base::app);
    for (size_t i = 0; i < data.size(); i++)
    {
        output << data[i] << "\n";
    }
}

// Append data to file name in "data_path_"
//  Input:
//          name : string with name of file (without extensions and path)
//          data : container with data to save, 2 column
void analyser::PrintVector(string name, vector<vector<double>> data)
{
    string file_name = data_path_ + "/" + name + ".txt";
    ofstream output(file_name, ios_base::app);
    for (size_t i = 0; i < data.size(); i++)
    {
        for (size_t j = 0; j < data[i].size(); j++)
        {
            output << data[i][j] << "\t";
        }
        output << "\n";
    }
}

// Append data to file name in "data_path_"
//  Input:
//          name : string with name of file (without extensions and path)
//          data : container with data to save, 3 column
void analyser::PrintVector(string name, vector<vector<vector<double>>> data)
{
    string file_name = data_path_ + "/" + name + ".txt";
    ofstream output(file_name);
    for (size_t i = 0; i < data.size(); i++)
    {
        for (size_t j = 0; j < data[i].size(); j++)
        {
            for (size_t k = 0; k < data[i][j].size(); k++)
            {
                output << data[i][j][k] << "\t";
            }
            output << "\n";
        }
    }
}

// Save pixelgrid data, ie. number of hits in pixels sum over all events, in "data_path_"
// Input:
//          name : string with name of file (without extensions and path)
void analyser::PrintPixels(string name)
{
    string file_name = data_path_ + "/pixeldata_run" + name + ".txt";
    ofstream output(file_name);
    for (size_t i = 0; i < pixel_grids_.size(); i++)
    {
        output << "## PLANE\t" << i << "\tPIXEL DATA\n";
        for (size_t j = 0; j < pixel_grids_[i].size(); j++)
        {
            output << pixel_grids_[i][j][1] << "\n";
        }
    }
}

// Save hotpixel data, ie. list of hotpixels in "data_path_"
// Input:
//          name : string with name of file (without extensions and path)
void analyser::PrintHotPixels(string name)
{
    for (size_t i = 0; i < hot_pixels_.size(); i++)
    {
        string file_name = data_path_ + "/hotpixels_run" + name + "_plane_" + to_string(i) + ".txt";
        ofstream output(file_name);
        for (size_t j = 0; j < hot_pixels_[i].size(); j++)
        {
            output << hot_pixels_[i][j] << "\n";
        }
    }
}

// Save disances between all pairs of points in the six different MIMOSA detectors in "data_path_"
// Input:
//          name : string with name of file (without extensions and path)
void analyser::PrintInterdistance(string name)
{
    string file_name = data_path_ + "/interdistance_data_" + name + ".txt";
    ofstream output(file_name);
    for (size_t i = 0; i < interplanar_distance_.size(); i++)
    {
        output << "## DATABLOCK\t" << i << "\tINTERDISTANCE DATA\n";
        for (size_t j = 0; j < interplanar_distance_[i].size(); j++)
        {
            output << interplanar_distance_[i][j] << "\n";
        }
        output << "\n\n";
    }
}

// Save x+y coordinates of hits in the six different MIMOSA detectors in "data_path_"
// Input:
//          name : string with name of file (without extensions and path)
void analyser::PrintHits(string name)
{
    string file_name = data_path_ + "/hits_coord_data_" + name + ".txt";
    ofstream output(file_name);
    for (size_t i = 0; i < hits_container_.size(); i++) // plane
    {
        output << "## PLANE\t" << i << "\tHIT DATA\n";
        for (size_t j = 0; j < hits_container_[i].size(); j++) // events
        {
            for (size_t k = 0; k < hits_container_[i][j].size(); k++) // hits
            {
                output << hits_container_[i][j][k][0] << ' ' << hits_container_[i][j][k][1] << '\n';
            }
        }
    }
}

// Save energy of emitted photons on disk
// Input:
//          file_name : string with name of file including file-extensions (.txt typically) and path
void analyser::PrintEnergy(string file_name)
{
    ofstream output;
    output.open(file_name, ios_base::app);
    for (size_t i = 0; i < energies.size(); i++)
    {
        for (size_t j = 0; j < energies[i].size(); j++)
        {
            output << energies[i][j] << "\n";
        }
    }
    output.close();
}

// Save incoming (column 1) and outgoing (column 2) angle of particle in MM  in "data_path_"
// Input:
//          name : string with name of file excluding file-extension and path
void analyser::PrintSlope(string name)
{
    string file_name = data_path_ + "/angles_" + name + ".txt";
    ofstream output(file_name);
    for (size_t i = 0; i < MM_slopes_.size(); i++)
    {
        output << MM_slopes_[i][0] << " " << MM_slopes_[i][1] << "\n";
    }
}

// Save angle of particle through M1+M2 in "data_path_"
// Input:
//          name : string with name of file excluding file-extension and path
void analyser::PrintM1M2Slope(string name)
{
    string file_name = data_path_ + "/angles_M1M2_" + name + ".txt";
    ofstream output(file_name);
    for (int i = 0; i < total_events_; i++)
    {
        for (size_t j = 0; j < M1_M2_slopes_[i].size(); j++)
        {
            output << M1_M2_slopes_[i][j][0] << "\t" << M1_M2_slopes_[i][j][1] << "\t" << M1_M2_slopes_[i][j][2] << "\t" << M1_M2_slopes_[i][j][3] << "\n";
        }
    }
}

// Calculate the slope of a particle path relative to detectors using small-angle approximation
// Input:
//          x : 2 element arma::vec with x coordinates of particle in plane upstream (0) and downstream (1)
//          y : 2 element arma::vec with y coordinates of particle in plane upstream (0) and downstream (1)
double analyser::CalculateSlope(vec x, vec y)
{
    return (y(1) - y(0)) / (x(1) - x(0));
}

// Calculate the slope of a particle path relative to detectors using small-angle approximation
// Input:
//          x0 : x coordinates of particle in plane upstream
//          x1 : x coordinates of particle in plane downstream
//          y0 : y coordinates of particle in plane upstream
//          y1 : y coordinates of particle in plane downstream
double analyser::CalculateSlope(double x0, double y0, double x1, double y1)
{
    return (y1 - y0) / (x1 - x0);
}

// Calculate direction through M1+M2 of particles that emit a photon.
// Calls member-function "BeamDivergencePhotons"
void analyser::PhotonTracksDivergence()
{
    string file_name = data_path_ + "/beam_direction_photons_" + run_number_ + ".txt";
    ofstream output;
    output.open(file_name, ios_base::app);
    for (size_t i = 0; i < paired_tracks_.size(); i++)
    {
        for (size_t j = 0; j < paired_tracks_[i].size(); j++)
        {
            double x_M1 = paired_tracks_[i][j][0][0][0],
                   y_M1 = paired_tracks_[i][j][0][0][1],
                   x_M2 = paired_tracks_[i][j][0][1][0],
                   y_M2 = paired_tracks_[i][j][0][1][1];
            vector<double> x_coordinates = {x_M1, x_M2};
            vector<double> y_coordinates = {y_M1, y_M2};
            vector<double> divergence = BeamDivergencePhotons(x_coordinates, y_coordinates);
            output << divergence[0] << "\t" << divergence[1] << "\n";
        }
    }
    output.close();
}

// Save M2 coordinates of particle that emit a photons in "data_path_"
// Input:
//          name : string with name of file excluding file-extension and path
void analyser::ImageCrystal()
{
    string file_name = data_path_ + "/crystal_image_" + run_number_ + ".txt";
    ofstream output;
    output.open(file_name);
    for (size_t i = 0; i < paired_tracks_.size(); i++)
    {
        for (size_t j = 0; j < paired_tracks_[i].size(); j++)
        {
            output << paired_tracks_[i][j][0][1][0] << "\t" << paired_tracks_[i][j][0][1][1] << "\n";
        }
    }
    output.close();
}

// Calculate distances between hits in a plane and save in "data_path_"
// Input:
//          plane : number of plane, from 0 to 5
void analyser::PrintInterplanarDistance(int plane)
{
    string plane_str = to_string(plane);
    string file_name = data_path_ + "/plane" + plane_str + "_" + run_number_ + "_dist.txt";
    ofstream output;
    output.open(file_name);
    double x1, x2, y1, y2;
    for (int i = 0; i < total_events_; i++)
    {
        for (size_t j = 0; j < hits_container_[plane][i].size(); j++)
        {
            for (size_t k = j + 1; k < hits_container_[plane][i].size(); k++)
            {
                x1 = hits_container_[plane][i][j][0];
                x2 = hits_container_[plane][i][k][0];
                y1 = hits_container_[plane][i][j][1];
                y2 = hits_container_[plane][i][k][1];
                output << sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2)) << "\n";
            }
        }
    }
    output.close();
}

// Calculate combined energy of matched pairs of particles, ie. photon energies
// Input:
//          paired_tracks : container with paired tracks
// Calls "CalculateSlope" member function.
double analyser::CalculatePairEnergy(vector<vector<vector<double>>> paired_tracks)
{
    double length = 0.15, charge = 1.6021766e-19, field_strength = 0.12, c = 299792458; // SI units
    double energy = 0;
    for (int i = 0; i < 2; i++)
    {
        vec M3M4_x = {paired_tracks[i][2][0], paired_tracks[i][3][0]};  // M3, M4 x coord
        vec M3_M4_zcoord_ = {detector_zcoord_[2], detector_zcoord_[3]}; // M3, M4 z coord
        vec M5M6_x = {paired_tracks[i][6][0], paired_tracks[i][7][0]};  // M5, M6 x coord
        vec M5M6_z = {detector_zcoord_[4], detector_zcoord_[5]};        // M5 M6 z coord
        double m = CalculateSlope(M3_M4_zcoord_, M3M4_x);
        double n = CalculateSlope(M5M6_z, M5M6_x);
        vector<double> slope = {m, n};
        MM_slopes_.push_back(slope);
        double ang = n - m;                                        // calculate deflection angle to find energy
        energy += charge * c * length * field_strength / abs(ang); // Joule
    }
    return energy;
}

// Method for pairing tracks in "tracks_" vector. Ensure that "ConstructTracks" method is used
// before running this method.
// Calls "CalculateEuclideanDistance" and "CalculatePairEnergy" member functions
void analyser::PairTracks(void)
{
    int pairs = 0;
    paired_tracks_.resize(tracks_.size());
    energies.resize(tracks_.size());
    vector<int> bad_events(tracks_.size());
    // Here we find pairs. We only want a track to pair exactly one other track.
    // To ensure this we iterate first from i = 0 -> number of tracks, then j = i + 1 -> number of tracks.
    // If multiple pairs are found we erase them from the tracks array.
    int photons_total_count = 0;
    for (size_t i = 0; i < tracks_.size(); i++)
    {
        int photons_event_count = 0;
        for (size_t j = 0; j < tracks_[i].size(); j++)
        {
            int matched_tracks = 0;
            for (size_t k = j + 1; k < tracks_[i].size(); k++)
            {
                double M3_distance = CalculateEuclideanDistance(tracks_[i][j][2][0], tracks_[i][j][2][1], tracks_[i][k][2][0], tracks_[i][k][2][1]);
                if (M3_distance < foil_paired_tracks_lim_)
                {
                    photons_event_count++;
                    matched_tracks++;
                    pairs++;
                    paired_tracks_[i].push_back({tracks_[i][j], tracks_[i][k]});
                    energies[i].push_back(CalculatePairEnergy({tracks_[i][j], tracks_[i][k]}));
                    // Remove additional matched track since we do not want it to potentially match any other track //
                    if (matched_tracks > 1)
                    {
                        tracks_[i].erase(tracks_[i].begin() + k);
                        k--;
                    }
                }
            } // tracks #2 done
        }     // tracks #1 done
        if (photons_event_count > 1)
        {
            bad_events[photons_total_count] = i;
            photons_total_count++;
        }
    }
    bad_events.resize(photons_total_count);
    // Discard photons from bad events //
    int erased_events = 0;
    int photons_removed = 0;
    for (size_t i = 0; i < bad_events.size(); i++)
    {
        photons_removed += energies[bad_events[i] - erased_events].size();
        energies.erase(energies.begin() + bad_events[i] - erased_events);
        paired_tracks_.erase(paired_tracks_.begin() + bad_events[i] - erased_events);
        erased_events++;
    }
    int N_photons_usable = 0;
    for (size_t i = 0; i < energies.size(); i++)
    {
        N_photons_usable += energies[i].size();
    }
    cout << "no. of events removed :\t\t\t" << erased_events << "\n";
    cout << "no. of photons before removal :\t\t" << pairs << "\nno. of photons after removal :\t\t" << N_photons_usable << "\n";
}

// Calculate M1-M2 direction of particles and save in M1_M2_slopes_
// Input:
//      event   : Which event is used
//      M1_hit  : Which hit in M1 is used
//      M2_hit  : Which hit in M2 is used
// Calls "CalculateSlope" member function
void analyser::BeamDivergence(int event, int M1_hit, int M2_hit)
{
    vec M1_M2_x = {hits_container_[0][event][M1_hit][0], hits_container_[1][event][M2_hit][0]};
    vec M1_M2_y = {hits_container_[0][event][M1_hit][1], hits_container_[1][event][M2_hit][1]};
    vec vec_M1_M2_z = {M1_M2_zcoord_[0], M1_M2_zcoord_[1]};
    vector<double> ang_pos;
    double ang_x = CalculateSlope(vec_M1_M2_z, M1_M2_x);
    double ang_y = CalculateSlope(vec_M1_M2_z, M1_M2_y);
    ang_pos = {ang_x, ang_y, hits_container_[0][event][M1_hit][0], hits_container_[0][event][M1_hit][1]};
    M1_M2_slopes_[event].push_back(ang_pos);
}

// Calculate M1-M2 direction of particles and returns angle in x and angle in y
// Input:
//      x  : Container with x-coordinate of upstream (0) and downstream (1) particle
//      y  : Container with y-coordinate of upstream (0) and downstream (1) particle
// Calls "CalculateSlope" member function
vector<double> analyser::BeamDivergencePhotons(vector<double> x, vector<double> y)
{
    vec vec_M1_M2_z = {M1_M2_zcoord_[0], M1_M2_zcoord_[1]};
    vector<double> ang;
    double ang_x = CalculateSlope(vec_M1_M2_z, x);
    double ang_y = CalculateSlope(vec_M1_M2_z, y);
    ang = {ang_x, ang_y};
    return ang;
}

// This method constructs particle tracks from M1 to M6. Tracks are saved in private memeber "tracks_".
// Input:
//      M1_M2_slope_lb_x:   The lower bound in x-entry angle into crystal. Defaults to -1E+17
//      M1_M2_slope_ub_x:   The upper bound in x-entry angle into crystal. Defaults to +1E+17
//      M1_M2_slope_lb_y:   The lower bound in y-entry angle into crystal. Defaults to -1E+17
//      M1_M2_slope_ub_y:   The upper bound in y-entry angle into crystal. Defaults to +1E+17
// Calls "RectilinearProjection", "BeamDivergence" and "CalculateEuclideanDistance" member functions.
void analyser::ConstructTracks(double M1_M2_slope_lb_x, double M1_M2_slope_ub_x, double M1_M2_slope_lb_y, double M1_M2_slope_ub_y)
{
    int M1_M4_tot_tracks = 0, M1_M6_tot_tracks = 0, M1_M3_tot_tracks = 0, M5_M6_tot_tracks = 0;
    events_within_cut_ = 0;
    tracks_.clear();
    tracks_.resize(total_events_);
    M1_M2_slopes_.resize(total_events_);
    int progress = 0;
    double start_time = omp_get_wtime();
    double T = 0;
#pragma omp parallel for
    for (int i = 0; i < total_events_; i++)
    {
        vector<vector<vector<double>>> M1_MM_tracks(hits_container_[0][i].size() * hits_container_[1][i].size() * hits_container_[2][i].size() * hits_container_[3][i].size());
        int M1_M4_tracks = 0;
        vector<vector<double>> M1_M6_track;
        bool included = false;
        // Here we construct tracks.Each loop iterates over hits in a mimosa detector(0->5).Conditions to find tracks checked periodically.All acceptable tracks are saved, which
        // means bad tracks are saved.These should later be discarded in "PairTracks" function
        for (size_t j = 0; j < hits_container_[0][i].size(); j++)
        {
            for (size_t k = 0; k < hits_container_[1][i].size(); k++)
            {
                vector<double> M1_M2_Proj = RectilinearProjection(hits_container_[0][i][j], hits_container_[1][i][k], M1_M2_zcoord_);
                BeamDivergence(i, j, k);
                // Only pick out certain angles of entry into crystal //
                if (M1_M2_slopes_[i].back()[0] > M1_M2_slope_lb_x and M1_M2_slopes_[i].back()[0] < M1_M2_slope_ub_x and M1_M2_slopes_[i].back()[1] > M1_M2_slope_lb_y and M1_M2_slopes_[i].back()[1] < M1_M2_slope_ub_y)
                {
                    if (!included)
                    {
#pragma omp atomic
                        events_within_cut_++;
                        included = true;
                    }
                    for (size_t l = 0; l < hits_container_[2][i].size(); l++)
                    {
                        double M1M2_d = CalculateEuclideanDistance(M1_M2_Proj[0], M1_M2_Proj[1], hits_container_[2][i][l][0], hits_container_[2][i][l][1]);
                        if (M1M2_d < M1_M2_proj_lim_)
                        {
                            vector<double> M2M3_Proj = RectilinearProjection(hits_container_[1][i][k], hits_container_[2][i][l], M2_M3_zcoord_);
#pragma omp atomic
                            M1_M3_tot_tracks++;
                            for (size_t m = 0; m < hits_container_[3][i].size(); m++)
                            {
                                double M2M3_d = CalculateEuclideanDistance(M2M3_Proj[0], M2M3_Proj[1], hits_container_[3][i][m][0], hits_container_[3][i][m][1]);
                                if (M2M3_d < M2_M3_proj_lim_)
                                {
                                    vector<double> M3M4_Proj = RectilinearProjection(hits_container_[2][i][l], hits_container_[3][i][m], M3_M4_zcoord_);
                                    M1_MM_tracks[M1_M4_tracks] = {hits_container_[0][i][j], hits_container_[1][i][k], hits_container_[2][i][l], hits_container_[3][i][m], M3M4_Proj};
                                    M1_M4_tracks++;
#pragma omp atomic
                                    M1_M4_tot_tracks++;
                                } // if statement for check of M2M3 limit done
                            }     // hits in M4 done
                        }         // if statement for check of M1M2 limit done
                    }             // hits in M3 done
                }                 // end angle cut
            }                     // hits in M2 done
        }                         // hits in M1 done
        M1_MM_tracks.resize(M1_M4_tracks);
        // Iterate over M6 -> M6 tracks //
        for (size_t n = 0; n < hits_container_[5][i].size(); n++)
        {
            vector<vector<double>> M1_M6_track;
            double shortest_distance = INFINITY;
            bool track_found = false;
            for (size_t o = 0; o < hits_container_[4][i].size(); o++)
            {
                vector<double> M6M5_Proj = RectilinearProjection(hits_container_[5][i][n], hits_container_[4][i][o], M6_M5_zcoord_);
                M5_M6_tot_tracks++;
                // For every M5-M6 combination, pick an M3M4 projection //
                for (size_t l = 0; l < M1_MM_tracks.size(); l++)
                {
                    vector<double> M3M4_Proj = M1_MM_tracks[l].back();
                    double M6M5_d = CalculateEuclideanDistance(M6M5_Proj[0], M6M5_Proj[1], M3M4_Proj[0], M3M4_Proj[1]);
                    if (M6M5_d < shortest_distance and M6M5_d < M6_M5_proj_lim_)
                    {
                        shortest_distance = M6M5_d;
                        vec vec_M5M6_z = {M6_M5_zcoord_[0], M6_M5_zcoord_[1]};
                        vec M5M6_y = {hits_container_[5][i][n][1], hits_container_[4][i][o][1]};
                        track_found = true;
                        M1_M6_track = {M1_MM_tracks[l][0], M1_MM_tracks[l][1], M1_MM_tracks[l][2], M1_MM_tracks[l][3], M3M4_Proj, M6M5_Proj, hits_container_[4][i][o], hits_container_[5][i][n]};
                    }
                }
            } // hits in M6 done

            if (track_found)
            {
#pragma omp critical
                {
                    tracks_[i].push_back(M1_M6_track);
                    M1_M6_tot_tracks++;
                }
            } // if statement for push_back of tracks
        }     // hits in M5 done
#pragma omp atomic
        progress++;

        // Report progress in terminal //
        if ((progress + 1) % (total_events_ / 10 + 1) == 0)
        {
            double dt = omp_get_wtime() - start_time;
            T += dt;
            cout << "Progress :\t" << floor(100 * double(progress) / (double)total_events_) << "%"
                 << "\ttime used :\t" << dt << "\ttotal time elapsed :\t" << T << "\ttime remaining :\t" << dt * (double)total_events_ / (total_events_ / 10 + 1) - T << "\n";
            start_time = omp_get_wtime();
        }
    } // events done
    cout << "\nTotal tracks from M1 -> M3:\t" << M1_M3_tot_tracks << "\n";
    cout << "Total tracks from M1 -> M4:\t" << M1_M4_tot_tracks << "\n";
    cout << "Total tracks from M6 -> M5:\t" << M5_M6_tot_tracks << "\n";
    cout << "Total tracks from M1 -> M6:\t" << M1_M6_tot_tracks << "\n";
    cout << "Number of events with tracks within cut:\t" << events_within_cut_ << "\n\n";
}
// This method determines the x+y direction of the axis using the normalized (to number of events) number
// of emitted photons as a function of the incoming angle of an electron. A 30murad window and
// 2urad increment is used in a the scan of entry angles. The number of photons is saved in "data_path_".
// Calls "RectilinearProjection", "CalculateSlope", "CalculateEuclideanDistance" and "SaveCutData" member functions.
void analyser::FindAxisCounts(void)
{
    vector<double> beam_parameters;
    double parameter;
    ifstream beam_parameters_file(beam_parameters_file_name_);
    if (beam_parameters_file.is_open())
    {
        while (true)
        {
            beam_parameters_file >> parameter;
            if (beam_parameters_file.eof())
            {
                break;
            }
            beam_parameters.push_back(parameter);
        }
        beam_parameters_file.close();
    }
    else
    {
        cout << "Unable to open file containing beam parameters";
    }
    vec vec_M1_M2_z = {M1_M2_zcoord_[0], M1_M2_zcoord_[1]};
    vec vec_M2_M3_z = {M2_M3_zcoord_[0], M2_M3_zcoord_[1]};
    double beam_x_direction = beam_parameters[2];
    double beam_y_direction = beam_parameters[3];
    double scan_x_range = 500E-6;
    double scan_y_range = 500E-6;
    double scan_window_size = 30E-6;
    double scan_increment = 2.0E-6;
    double ycut_lb = beam_y_direction - scan_y_range - scan_window_size / 2.0;
    double ycut_ub = beam_y_direction - scan_y_range + scan_window_size / 2.0;
    double xcut_lb = beam_x_direction - scan_x_range - scan_window_size / 2.0;
    double xcut_ub = beam_x_direction - scan_x_range + scan_window_size / 2.0;
    int no_cuts, no_cuts_x = 2 * scan_x_range / scan_increment, no_cuts_y = 2 * scan_y_range / scan_increment;
    vector<vector<double>> x_cuts(no_cuts_x);
    vector<vector<double>> y_cuts(no_cuts_y);
    for (size_t i = 0; i < x_cuts.size(); i++)
    {
        x_cuts[i].resize(2);
        x_cuts[i][0] = xcut_lb;
        x_cuts[i][1] = xcut_ub;
        xcut_lb += scan_increment;
        xcut_ub += scan_increment;
    }
    for (size_t i = 0; i < y_cuts.size(); i++)
    {
        y_cuts[i].resize(2);
        y_cuts[i][0] = ycut_lb;
        y_cuts[i][1] = ycut_ub;
        ycut_lb += scan_increment;
        ycut_ub += scan_increment;
    }
    string file_name;
    vector<vector<double>> photons_in_cut;
    // Iterate over cuts //
    for (int direction = 0; direction < 2; direction++)
    {
        vector<vector<double>> cuts;
        if (direction == 0)
        {
            file_name = data_path_ + "/no_photons_x_counts_" + (string)run_number_ + ".txt";
            no_cuts = no_cuts_x;
            cuts = x_cuts;
            cout << "\nMaking cuts in x-direction; total number of cuts:\t" << no_cuts << "\n\n";
        }
        else
        {
            file_name = data_path_ + "/no_photons_y_counts_" + (string)run_number_ + ".txt";
            no_cuts = no_cuts_y;
            cuts = y_cuts;
            cout << "\nMaking cuts in y-direction; total number of cuts:\t" << no_cuts << "\n\n";
        }
        photons_in_cut.resize(no_cuts);
        int progress = 0;
        double start_time = omp_get_wtime();
        double T = 0;
#pragma omp parallel for
        for (size_t cut = 0; cut < cuts.size(); cut++)
        {
            double M1_M2_slope_lb_x, M1_M2_slope_ub_x, M1_M2_slope_lb_y, M1_M2_slope_ub_y;
            double events_within_cut = 0;
            double photons_within_cut = 0;
            if (direction == 0)
            {
                M1_M2_slope_lb_x = cuts[cut][0];
                M1_M2_slope_ub_x = cuts[cut][1];
                M1_M2_slope_lb_y = (-1.0) * numeric_limits<double>::infinity();
                M1_M2_slope_ub_y = numeric_limits<double>::infinity();
            }
            else
            {
                M1_M2_slope_lb_x = (-1.0) * numeric_limits<double>::infinity();
                M1_M2_slope_ub_x = numeric_limits<double>::infinity();
                M1_M2_slope_lb_y = cuts[cut][0];
                M1_M2_slope_ub_y = cuts[cut][1];
            }
            // Determine the full tracks within cut //
            for (int i = 0; i < total_events_; i++)
            {
                vector<vector<vector<double>>> tracks;
                vector<vector<vector<double>>> M1_MM_tracks(hits_container_[0][i].size() * hits_container_[1][i].size() * hits_container_[2][i].size() * hits_container_[3][i].size());
                int M1_M4_tracks = 0;
                vector<vector<double>> M1_M6_track;
                bool included = false;
                for (size_t j = 0; j < hits_container_[0][i].size(); j++)
                {
                    for (size_t k = 0; k < hits_container_[1][i].size(); k++)
                    {
                        vector<double> M1_M2_Proj = RectilinearProjection(hits_container_[0][i][j], hits_container_[1][i][k], M1_M2_zcoord_);
                        vec M1_M2_x = {hits_container_[0][i][j][0], hits_container_[1][i][k][0]};
                        vec M1_M2_y = {hits_container_[0][i][j][1], hits_container_[1][i][k][1]};
                        vec vec_M1_M2_z = {M1_M2_zcoord_[0], M1_M2_zcoord_[1]};
                        double ang_x = CalculateSlope(vec_M1_M2_z, M1_M2_x);
                        double ang_y = CalculateSlope(vec_M1_M2_z, M1_M2_y);
                        // Only pick out certain angles of entry into crystal //
                        if (ang_x > M1_M2_slope_lb_x and ang_x < M1_M2_slope_ub_x and ang_y > M1_M2_slope_lb_y and ang_y < M1_M2_slope_ub_y)
                        {
                            if (!included)
                            {
#pragma omp critical
                                {
                                    events_within_cut++;
                                    included = true;
                                }
                            }
                            for (size_t l = 0; l < hits_container_[2][i].size(); l++)
                            {
                                double M1M2_d = CalculateEuclideanDistance(M1_M2_Proj[0], M1_M2_Proj[1], hits_container_[2][i][l][0], hits_container_[2][i][l][1]);
                                if (M1M2_d < M1_M2_proj_lim_)
                                {
                                    vector<double> M2M3_Proj = RectilinearProjection(hits_container_[1][i][k], hits_container_[2][i][l], M2_M3_zcoord_);
                                    for (size_t m = 0; m < hits_container_[3][i].size(); m++)
                                    {
                                        double M2M3_d = CalculateEuclideanDistance(M2M3_Proj[0], M2M3_Proj[1], hits_container_[3][i][m][0], hits_container_[3][i][m][1]);
                                        if (M2M3_d < M2_M3_proj_lim_)
                                        {
                                            vector<double> M3M4_Proj = RectilinearProjection(hits_container_[2][i][l], hits_container_[3][i][m], M3_M4_zcoord_);
                                            M1_MM_tracks[M1_M4_tracks] = {hits_container_[0][i][j], hits_container_[1][i][k], hits_container_[2][i][l], hits_container_[3][i][m], M3M4_Proj};
                                            M1_M4_tracks++;
                                        } // if statement for check of M2M3 limit done
                                    }     // hits in M4 done
                                }         // if statement for check of M1M2 limit done
                            }             // hits in M3 done
                        }                 // end angle cut
                    }                     // hits in M2 done
                }                         // hits in M1 done
                M1_MM_tracks.resize(M1_M4_tracks);
                // Iterate over M6 -> M6 tracks //
                for (size_t n = 0; n < hits_container_[5][i].size(); n++)
                {
                    vector<vector<double>> M1_M6_track;
                    double shortest_distance = INFINITY;
                    bool track_found = false;
                    for (size_t o = 0; o < hits_container_[4][i].size(); o++)
                    {
                        vector<double> M6M5_Proj = RectilinearProjection(hits_container_[5][i][n], hits_container_[4][i][o], M6_M5_zcoord_);
                        // For every M5-M6 combination, pick an M3M4 projection //
                        for (size_t l = 0; l < M1_MM_tracks.size(); l++)
                        {
                            vector<double> M3M4_Proj = M1_MM_tracks[l].back();
                            double M6M5_d = CalculateEuclideanDistance(M6M5_Proj[0], M6M5_Proj[1], M3M4_Proj[0], M3M4_Proj[1]);
                            if (M6M5_d < shortest_distance and M6M5_d < M6_M5_proj_lim_)
                            {
                                shortest_distance = M6M5_d;
                                track_found = true;
                                M1_M6_track = {M1_MM_tracks[l][0], M1_MM_tracks[l][1], M1_MM_tracks[l][2], M1_MM_tracks[l][3], M3M4_Proj, M6M5_Proj, hits_container_[4][i][o], hits_container_[5][i][n]};
                            }
                        }
                    } // hits in M6 done
                    if (track_found)
                    {
                        tracks.push_back(M1_M6_track);
                    }
                } // hits in M5 done
                // pair tracks
                for (size_t j = 0; j < tracks.size(); j++)
                {
                    int matched_tracks_counter = 0;
                    for (size_t k = j + 1; k < tracks.size(); k++)
                    {
                        double M3_distance = CalculateEuclideanDistance(tracks[j][2][0], tracks[j][2][1], tracks[k][2][0], tracks[k][2][1]);
                        if (M3_distance < foil_paired_tracks_lim_)
                        {
                            photons_within_cut++;
                            matched_tracks_counter++;
                            if (matched_tracks_counter > 1)
                            {
                                tracks.erase(tracks.begin() + k);
                                k--;
                            }
                        }
                    } // tracks #2 done
                }     // tracks #1 done
            }         // events done
            double photons_per_event;
            (events_within_cut != 0) ? photons_per_event = photons_within_cut / events_within_cut : photons_per_event = 0.0;
#pragma omp critical
            {
                photons_in_cut[cut] = {cuts[cut][0] + scan_window_size / 2.0, photons_per_event};
                progress++;
            }
            if ((progress + 1) % (no_cuts / 10 + 1) == 0)
            {
                double dt = omp_get_wtime() - start_time;
                T += dt;
                cout << "Progress :\t" << floor(100 * double(progress) / (double)no_cuts) << "%"
                     << "\ttime used :\t" << dt << "\ttotal time elapsed :\t" << T << "\ttime remaining :\t" << dt * (double)no_cuts / (no_cuts / 10 + 1) - T << "\n";
                start_time = omp_get_wtime();
            }
        } // cuts done
        SaveCutData(photons_in_cut, file_name);
    } // direction done
}

// This method determines the x+y direction of the axis using the deflection of an electron after it
// has traversed the crystal as a function of the incoming angle of an electron. A 30murad window and
// 2urad increment is used in a the scan of entry angles. The deflection is saved in "data_path_".
// Calls "RectilinearProjection", "CalculateSlope", "CalculateEuclideanDistance" and "SaveCutData" member functions.
void analyser::FindAxisDeflection(void)
{
    vector<double> beam_parameters;
    double parameter;
    ifstream beam_parameters_file(beam_parameters_file_name_);
    if (beam_parameters_file.is_open())
    {
        while (true)
        {
            beam_parameters_file >> parameter;
            if (beam_parameters_file.eof())
            {
                break;
            }
            beam_parameters.push_back(parameter);
        }
        beam_parameters_file.close();
    }
    else
    {
        cout << "Unable to open file containing beam parameters";
    }
    vec vec_M1_M2_z = {M1_M2_zcoord_[0], M1_M2_zcoord_[1]};
    vec vec_M2_M3_z = {M2_M3_zcoord_[0], M2_M3_zcoord_[1]};
    double beam_x_direction = beam_parameters[2];
    double beam_y_direction = beam_parameters[3];
    double scan_x_range = 500e-6;
    double scan_y_range = 500e-6;
    double scan_window_size = 30E-6;
    double scan_increment = 2.0E-6;
    double ycut_lb = beam_y_direction - scan_y_range - scan_window_size / 2.0;
    double ycut_ub = beam_y_direction - scan_y_range + scan_window_size / 2.0;
    double xcut_lb = beam_x_direction - scan_x_range - scan_window_size / 2.0;
    double xcut_ub = beam_x_direction - scan_x_range + scan_window_size / 2.0;
    int no_cuts, no_cuts_x = 2 * scan_x_range / scan_increment, no_cuts_y = 2 * scan_y_range / scan_increment;
    vector<vector<double>> x_cuts(no_cuts_x);
    vector<vector<double>> y_cuts(no_cuts_y);

    for (size_t i = 0; i < x_cuts.size(); i++)
    {
        x_cuts[i].resize(2);
        x_cuts[i][0] = xcut_lb;
        x_cuts[i][1] = xcut_ub;
        xcut_lb += scan_increment;
        xcut_ub += scan_increment;
    }
    for (size_t i = 0; i < y_cuts.size(); i++)
    {
        y_cuts[i].resize(2);
        y_cuts[i][0] = ycut_lb;
        y_cuts[i][1] = ycut_ub;
        ycut_lb += scan_increment;
        ycut_ub += scan_increment;
    }
    string file_name;
    vector<vector<double>> mean_deflection;
    for (int direction = 0; direction < 2; direction++)
    {
        vector<vector<double>> cuts;
        if (direction == 0)
        {
            file_name = data_path_ + "/no_photons_x_defl_" + (string)run_number_ + ".txt";
            no_cuts = no_cuts_x;
            cuts = x_cuts;
            cout << "Making cuts in x-direction; total number of cuts:\t" << no_cuts << "\n\n";
        }
        else
        {
            file_name = data_path_ + "/no_photons_y_defl_" + (string)run_number_ + ".txt";
            no_cuts = no_cuts_y;
            cuts = y_cuts;
            cout << "\nMaking cuts in y-direction; total number of cuts:\t" << no_cuts << "\n\n";
        }
        mean_deflection.resize(no_cuts);
        int progress = 0;
        double start_time = omp_get_wtime();
        double T = 0;
#pragma omp parallel for
        for (size_t cut = 0; cut < cuts.size(); cut++)
        {
            double M1_M2_slope_lb_x, M1_M2_slope_ub_x, M1_M2_slope_lb_y, M1_M2_slope_ub_y;
            double deflection_sum = 0.0;
            if (direction == 0)
            {
                M1_M2_slope_lb_x = cuts[cut][0];
                M1_M2_slope_ub_x = cuts[cut][1];
                M1_M2_slope_lb_y = (-1.0) * numeric_limits<double>::infinity();
                M1_M2_slope_ub_y = numeric_limits<double>::infinity();
            }
            else
            {
                M1_M2_slope_lb_x = (-1.0) * numeric_limits<double>::infinity();
                M1_M2_slope_ub_x = numeric_limits<double>::infinity();
                M1_M2_slope_lb_y = cuts[cut][0];
                M1_M2_slope_ub_y = cuts[cut][1];
            }
            int photons_in_cut = 0;
            for (int i = 0; i < total_events_; i++)
            {
                for (size_t j = 0; j < hits_container_[0][i].size(); j++)
                {
                    for (size_t k = 0; k < hits_container_[1][i].size(); k++)
                    {
                        vector<double> M1_M2_Proj = RectilinearProjection(hits_container_[0][i][j], hits_container_[1][i][k], M1_M2_zcoord_);
                        vec M1_M2_x = {hits_container_[0][i][j][0], hits_container_[1][i][k][0]};
                        vec M1_M2_y = {hits_container_[0][i][j][1], hits_container_[1][i][k][1]};
                        double M1_M2_ang_x = CalculateSlope(vec_M1_M2_z, M1_M2_x);
                        double M1_M2_ang_y = CalculateSlope(vec_M1_M2_z, M1_M2_y);
                        // Only pick out certain angles of entry into crystal //
                        if (M1_M2_ang_x > M1_M2_slope_lb_x and M1_M2_ang_x < M1_M2_slope_ub_x and M1_M2_ang_y > M1_M2_slope_lb_y and M1_M2_ang_y < M1_M2_slope_ub_y)
                        {
                            for (size_t l = 0; l < hits_container_[2][i].size(); l++)
                            {
                                double M1M2_d = CalculateEuclideanDistance(M1_M2_Proj[0], M1_M2_Proj[1], hits_container_[2][i][l][0], hits_container_[2][i][l][1]);
                                if (M1M2_d < M1_M2_proj_lim_)
                                {
                                    vec M2_M3_x = {hits_container_[1][i][k][0], hits_container_[2][i][l][0]};
                                    vec M2_M3_y = {hits_container_[1][i][k][1], hits_container_[2][i][l][1]};
                                    double M2_M3_ang_x = CalculateSlope(vec_M2_M3_z, M2_M3_x);
                                    double M2_M3_ang_y = CalculateSlope(vec_M2_M3_z, M2_M3_y);
                                    // Calculate sum of delta angles //
                                    (direction == 0) ? deflection_sum += abs(M2_M3_ang_x - M1_M2_ang_x) : deflection_sum += abs(M2_M3_ang_y - M1_M2_ang_y);
                                    photons_in_cut++;
                                }
                            } // hits in M3 done
                        }     // end cut
                    }         // hits in M2 done
                }             // hits in M1 done
            }                 // events done
            // Calculate mean of delta angles //
            deflection_sum /= photons_in_cut;
            mean_deflection[cut] = {cuts[cut][0] + scan_window_size / 2.0, deflection_sum};
#pragma omp atomic
            progress++;
            if ((progress + 1) % (no_cuts / 10 + 1) == 0)
            {
                double dt = omp_get_wtime() - start_time;
                T += dt;
                cout << "Progress :\t" << floor(100 * double(progress) / (double)no_cuts) << "%"
                     << "\ttime used :\t" << dt << "\ttotal time elapsed :\t" << T << "\ttime remaining :\t" << dt * (double)no_cuts / (no_cuts / 10 + 1) - T << "\n";
                start_time = omp_get_wtime();
            }
        } // cuts done
        SaveCutData(mean_deflection, file_name);
    }
}

// This method is called when running either "FindAxisCounts" or "FindAxisDeflection" and is used to save
// the data recorded when running these methods. It checks if a data-file exists, and if it does, it averages
// the previously recorded and newly recorded normalized counts / deflection for a given entry angle, and then
// overwrites the data on disk.
//
// Input:
//      cut_data:   container with normalized counts / deflection for different entry angles
//      file_name:  name of file to save on disk
void analyser::SaveCutData(vector<vector<double>> cut_data, string file_name)
{
    double dummy, number_of_photons;
    int cut_count = 0;
    vector<double> input_data(cut_data.size());
    ifstream input;
    input.open(file_name);
    if (input.is_open())
    {
        while (true)
        {
            input >> dummy >> number_of_photons;
            if (input.eof())
            {
                break;
            }
            input_data[cut_count] = number_of_photons;
            cut_count++;
        }
        input.close();
    }
    for (size_t i = 0; i < cut_data.size(); i++)
    {
        cut_data[i][1] += input_data[i];
        cut_data[i][1] /= 2;
    }
    ofstream output;
    output.open(file_name);
    for (size_t i = 0; i < cut_data.size(); i++)
    {
        output << cut_data[i][0] << "\t" << cut_data[i][1] << "\n";
    }
    output.close();
}

// Calculate distance between observed and projected hit in a detector across all events. Used for alignment.
//
// Input:
//      Hits0: container with hits in detector "0", detector furthest upstream
//      Hits1: container with hits in detector "1", middle detector
//      Hits2: container with hits in detector "2", detector furthest downstream (the one we project to)
// Calls "RectilinearProjection" and "CalculateEuclideanDistance" member functions
vector<double> analyser::CalculateInterdistance(vector<vector<vector<double>>> Hits0, vector<vector<vector<double>>> Hits1, vector<vector<vector<double>>> Hits2, vector<double> z)
{
    vector<double> interdistance;
    for (int i = 0; i < total_events_; i++)
    {
        vector<vector<double>> hits0 = Hits0[i];
        vector<vector<double>> hits1 = Hits1[i];
        vector<vector<double>> hits2 = Hits2[i];
        for (size_t j = 0; j < hits0.size(); j++)
        {
            for (size_t k = 0; k < hits1.size(); k++)
            {
                vector<double> hitp = RectilinearProjection(hits0[j], hits1[k], z);
                for (size_t l = 0; l < hits2.size(); l++)
                {
                    double distance = CalculateEuclideanDistance(hitp[0], hitp[1], hits2[l][0], hits2[l][1]);
                    interdistance.push_back(distance);
                }
            }
        }
    }
    return interdistance;
}

// Calclate euclidean distance between two points (A and B) in a plane.
// Input:
//          Ax: x-coordinate of point A
//          Ay: y-coordinate of point A
//          Bx: x-coordinate of point B
//          Bx: y-coordinate of point B
double analyser::CalculateEuclideanDistance(double Ax, double Ay, double Bx, double By)
{
    return sqrt((Bx - Ax) * (Bx - Ax) + (By - Ay) * (By - Ay));
}

// Calculate rectilinear projection (point C) using two points A and B
// Input:
//      point_A    : x+y coordinates of point A
//      point_B    : x+y  coorindates of point B
//      z          : z-coordinate of point A, B and C
vector<double> analyser::RectilinearProjection(vector<double> point_A, vector<double> point_B, vector<double> z)
{
    vector<double> projection(2);
    projection[0] = (z[2] - z[0]) * (point_A[0] - point_B[0]) / (z[0] - z[1]) + point_A[0];
    projection[1] = (z[2] - z[0]) * (point_A[1] - point_B[1]) / (z[0] - z[1]) + point_A[1];
    return projection;
}

// Align coordinates in a plane by applying an affine transformation.
// Input:
//          mat_Tot : Alignment matrix
//          Hits0   : Container with x+y coordinates of hits in upstream plane (already aligned)
//          Hits1   : Container with x+y coordinates of hits in middle plane (already aligned)
//          Hits2   : Container with x+y coordinates of hits in downsstream plane (not aligned)
//          z       : z-coordinate of planes
// Calls "ConstructTMatrix", "AdjustCoordinates" and "CheckMatrixConvergence" member functions
void analyser::AlignPlane(mat &mat_Tot, vector<vector<vector<double>>> Hits0, vector<vector<vector<double>>> Hits1, vector<vector<vector<double>>> &Hits2, vector<double> z)
{
    vector<double> alignment_radius_lim;
    if (z.back() == detector_zcoord_[2] or z.back() == detector_zcoord_[3])
    {
        alignment_radius_lim = {alignment_radius_lim_[0], alignment_radius_lim_[1], alignment_radius_lim_[2]};
    }
    else
    {
        alignment_radius_lim = alignment_radius_lim_;
    }
    for (size_t j = 0; j < alignment_radius_lim.size(); j++)
    {
        double dr_crit = alignment_radius_lim[j];
        bool convergence;
        mat mat_temp_T = zeros<mat>(3, 3);
        do
        {
            mat mat_T = ConstructTMatrix(Hits0, Hits1, Hits2, dr_crit, z);
            mat_Tot = mat_T * mat_Tot;
            AdjustCoordinates(Hits2, mat_T);
            convergence = CheckMatrixConvergence(mat_T, mat_temp_T);
            mat_temp_T = mat_T;
        } while (!convergence);
    }
}

// Calculate the affine transformation matrix
// Input:
//          Hits0       : Container with x+y coordinates of hits in upstream plane (already aligned)
//          Hits1       : Container with x+y coordinates of hits in middle plane (already aligned)
//          Hits2       : Container with x+y coordinates of hits in downsstream plane (not aligned)
//          dr_limit    : Convergence criteria, ie. maximum allowed distance between Hits2 and projections using Hits0 and Hits1
//          z           : z-coordinate of planes
// Calls "RectilinearProjection", "CalculateEuclideanDistance" and "SaveHit" member functions
mat analyser::ConstructTMatrix(vector<vector<vector<double>>> Hits0, vector<vector<vector<double>>> Hits1, vector<vector<vector<double>>> Hits2, double dr_limit, vector<double> z)
{
    int rowno = 0;
    mat mat_expected = zeros<mat>(Hits1.size() * 200, 3), mat_observed = zeros<mat>(Hits1.size() * 200, 3);
    for (size_t k = 0; k < Hits0.size(); k++)
    {
        vector<vector<double>> hits0 = Hits0[k];
        vector<vector<double>> hits1 = Hits1[k];
        vector<vector<double>> hits2 = Hits2[k];
        for (size_t l = 0; l < hits0.size(); l++)
        {
            for (size_t m = 0; m < hits1.size(); m++)
            {
                vector<double> hitp = RectilinearProjection(hits0[l], hits1[m], z);
                for (size_t n = 0; n < hits2.size(); n++)
                {
                    double dr = CalculateEuclideanDistance(hitp[0], hitp[1], hits2[n][0], hits2[n][1]);
                    if (dr < dr_limit)
                    {
                        SaveHit(hitp, hits2[n], mat_expected, mat_observed, rowno);
                        rowno++;
                    }
                }
            }
        }
    }
    mat_expected.resize(rowno, 3);
    mat_observed.resize(rowno, 3);
    return mat_expected.t() * mat_observed * (mat_observed.t() * mat_observed).i(); // formula derived in Tobias' thesis
}

// Apply affine transformation.
// Input:
//          hits:   Container with x+y coordinates
//          mat_T:  Affine transformation matrix
void analyser::AdjustCoordinates(vector<vector<vector<double>>> &hits, mat mat_T)
{
    for (size_t i = 0; i < hits.size(); i++)
    {
        for (size_t j = 0; j < hits[i].size(); j++)
        {
            vec hit = {{hits[i][j][0]}, {hits[i][j][1]}, {1}};
            vec temp = mat_T * hit;
            hits[i][j][0] = temp[0];
            hits[i][j][1] = temp[1];
        }
    }
}

// Save the projected and closest observed hit in a detector
// Input:
//          projected_hit   : x+y coordinates of projected hit
//          observed_hit    : x+y coordinates of observed hit
//          mat_expected    : matrix with expected hits in rows
//          mat_observed    : matrix with observed hits in rows
//          row             : row number
void analyser::SaveHit(vector<double> projected_hit, vector<double> observed_hit, mat &mat_expected, mat &mat_observed, int row)
{
    rowvec vec_expected = {projected_hit[0], projected_hit[1], 1};
    rowvec vec_observed = {observed_hit[0], observed_hit[1], 1};
    mat_expected.row(row) = vec_expected;
    mat_observed.row(row) = vec_observed;
}

// Determine if two matrices are essentially the same
// Input:
//          A   : Matrix A
//          B   : Matrix B
bool analyser::CheckMatrixConvergence(mat A, mat B)
{
    return sum(abs(vectorise(A) - vectorise(B))) < T_matrix_convergence_tol_;
}

// Fill container with interplanar distances
// Calls "CalculateInterdistance" member function
void analyser::FillInterplanarDistanceContainer(void)
{
    for (size_t i = 0; i < 4; i++)
    {
        vector<double> zproj = {detector_zcoord_[i], detector_zcoord_[i + 1], detector_zcoord_[i + 2]};
        interplanar_distance_[i] = CalculateInterdistance(hits_container_[i], hits_container_[i + 1], hits_container_[i + 2], zproj);
    }
}

// Align hits in detectors with transformation matrix already determined and fill interdistance container
void analyser::AlignWithTMatrix(void)
{
    for (size_t i = 0; i < 4; i++)
    {
        vector<double> zproj = {detector_zcoord_[i], detector_zcoord_[i + 1], detector_zcoord_[i + 2]};
        mat mat_T = T.slice(i);
        for (size_t j = 0; j < hits_container_[i + 2].size(); j++)
        {
            for (size_t k = 0; k < hits_container_[i + 2][j].size(); k++)
            {
                vec hit = {{hits_container_[i + 2][j][k][0]}, {hits_container_[i + 2][j][k][1]}, {1}};
                vec temp = mat_T * hit;
                hits_container_[i + 2][j][k][0] = temp[0];
                hits_container_[i + 2][j][k][1] = temp[1];
            }
        }
        cout << "Aligned plane " << i + 2 << "\n";
    }
}

// Aligns hits in detecors by calculating and applying transformation matrix
// Calls "AlignPlane" member function
void analyser::AlignWithoutTMatrix(void)
{
    for (int i = 0; i < 4; i++)
    {
        mat mat_T = eye<mat>(3, 3);
        vector<double> zproj = {detector_zcoord_[i], detector_zcoord_[i + 1], detector_zcoord_[i + 2]};
        AlignPlane(mat_T, hits_container_[i], hits_container_[i + 1], hits_container_[i + 2], zproj);
        cout << "Aligned plane " << i + 2 << "\n";
        T.slice(i) = mat_T;
    }
    T.save(data_path_ + "/Align/alignment_matrix.txt", arma_ascii);
}

// Make pixel grid, ie. number of hits in a pixel.
// Input:
//      pixel_grid : Contains 2 columns {pixel number, hits in pixel}
// This function does not fill the pixel grid!
void analyser::MakeGrid(vector<vector<double>> &pixel_grid)
{
    pixel_grid.resize(column_count_ * row_count_);
    for (size_t i = 0; i < pixel_grid.size(); i++)
    {
        pixel_grid[i] = {(double)i, (double)0};
    }
}

// Save coordinates of hits in a given detector and fill pixel grid
// Input:
//          hit_coordinates : Container where coordinates are stored
//          pixel_grid      : Contains 2 columns {pixel number, hits in pixel}
//          plane           : detector number
// Calls "Coord2Pixel" member function
void analyser::ExtractHitData(vector<vector<vector<double>>> &hit_coordinates, vector<vector<double>> &pixel_grid, int plane)
{
    hit_coordinates.resize(total_events_);
    for (int i = 0; i < total_events_; i++)
    {
        vector<vector<double>> event_data = events_container_[i];
        for (size_t j = 0; j < event_data.size(); j++)
        {
            vector<double> hit_data = event_data[j];
            if (hit_data[2] == plane)
            {
                int pixel = Coord2Pixel(hit_data[0], hit_data[1]);
                pixel_grid[pixel][1] += 1;
                hit_coordinates[i].push_back({hit_data[0], hit_data[1]});
            }
        }
    }
}

// Determine pixel number from coordinates
// Input:
//          hit_x_coordinate : x coordinate of hit
//          hit_y_coordinate : y coordinate of hit
int analyser::Coord2Pixel(double hit_x_coordinate, double hit_y_coordinate)
{
    double dx = (detector_xcoord_max_ - detector_xcoord_min_) / column_count_;
    double dy = (detector_ycoord_max_ - detector_ycoord_min_) / row_count_;
    int colno = (hit_x_coordinate + detector_xcoord_max_) / dx;
    int rowno = (hit_y_coordinate + detector_ycoord_max_) / dy;
    return rowno * column_count_ + colno;
}

// Count total number of hits across all events
// Input:
//          count           : Total number of hits
//          hit_coordinates : hits in events
void analyser::CountHits(int &count, vector<vector<vector<double>>> hit_coordinates)
{
    count = 0;
    for (size_t i = 0; i < hit_coordinates.size(); i++)
    {
        count += hit_coordinates[i].size();
    }
}

// Sort-function to sort vector by 2nd column, descending
// Input:
//          kiterator_begin : Iterator to beginning of vector
//          kiterator_end : Iterator to end of vector
bool analyser::SortVectorDescending(const vector<double> &kiterator_begin, const vector<double> &kiterator_end)
{
    return kiterator_begin[1] > kiterator_end[1];
}

// Locate hot pixels and save pixel number in "hot_pixels"
// Input:
//          pixelgrid   : contains pixels and total hits in a pixel
//          hot_pixels  : empty container, is filled with which pixels are hot pixels
// Uses "SortVectorDescending" member function.
void analyser::LocateHotPixels(vector<vector<double>> pixelgrid, vector<int> &hot_pixels)
{
    // We need to identify hot pixels.We do this by sorting the pixelgrid by number of pixels, descendind,
    // and then saving the pixelno.of the hotpixels.
    sort(pixelgrid.begin(), pixelgrid.end(), SortVectorDescending);
    double lim = 4E-04 * total_events_;
    for (size_t i = 0; i < pixelgrid.size(); i++)
    {
        if (pixelgrid[i][1] > lim)
        {
            hot_pixels.push_back(pixelgrid[i][0]);
        }
        else
        {
            break;
        }
    }
}

// Remove the number of counts in pixelgrid for hotpixels, and erase hitcoordinates from hotpixels.
// Input:
//          hit_coordinates : Container with Hits in detectorh
//          pixel_grid      : Pixel grid for detector
//          hot_pixels      : Container with hot pixels in detector
void analyser::RemoveHotPixels(vector<vector<vector<double>>> &hit_coordinates, vector<vector<double>> &pixel_grid, vector<int> hot_pixels)
{
    for (size_t i = 0; i < hot_pixels.size(); i++)
    {
        int hot_pixel = hot_pixels[i];
        pixel_grid[hot_pixel][1] = 0;
        for (size_t j = 0; j < hit_coordinates.size(); j++)
        {
            for (size_t k = 0; k < hit_coordinates[j].size(); k++)
            {
                double hit_x_coordinate = hit_coordinates[j][k][0];
                double hit_y_coordinate = hit_coordinates[j][k][1];
                int pixel = Coord2Pixel(hit_x_coordinate, hit_y_coordinate);
                if (pixel == hot_pixel)
                {
                    hit_coordinates[j].erase(hit_coordinates[j].begin() + k);
                    k--;
                }
            }
        }
    }
}

// Return vector of linearly spaced doubles
// Input:
//      min : minimum double in vector
//      max : maximum double in vector
//      N   : number of doubles in vector
vector<double> analyser::Linspace(double min, double max, int N)
{
    vector<double> range;
    range.resize(N);
    double delta = (max - min) / (N - 1); // step-size
    for (int i = 0; i < N - 1; i++)
    {
        range[i] = min + i * delta;
    }
    range.back() = max; // Set last entry to "max". This ensures that we get a range from min -> max
    return range;
}

// Extract raw data from root file, do pre-processing and save usable data in vectors
void analyser::ExtractRootData(void)
{
    // Open root file and obtain relevant data //
    TFile *f = TFile::Open(file_name_);
    TTree *T1 = (TTree *)f->Get("T");
    TLeaf *Hpk = (TLeaf *)T1->GetLeaf("fAHits.Hpk");
    TLeaf *Hu = (TLeaf *)T1->GetLeaf("fAHits.Hu");
    TLeaf *Hv = (TLeaf *)T1->GetLeaf("fAHits.Hv");
    TLeaf *fAHitsN = (TLeaf *)T1->GetLeaf("fAHitsN");
    TLeaf *ENumberOfTriggers = (TLeaf *)T1->GetLeaf("fHeader.ENumberOfTriggers");
    int total_events = T1->GetEntries();
    cout << "\nTotal number of events :\t" << total_events << "\n";
    events_container_.resize(total_events);
    total_events_ = 0;
    // Sort data into std::vector //
    for (int i = 0; i < total_events; i++)
    {
        T1->GetEntry(i);
        int hits_in_event = fAHitsN->GetValue();
        int triggers_in_event = ENumberOfTriggers->GetValue();
        if (triggers_in_event == 2 && hits_in_event > 0)
        {
            vector<vector<double>> event_data(hits_in_event);
            for (int j = 0; j < hits_in_event; j++)
            {
                vector<double> hit_data(3);
                hit_data[0] = Hu->GetValue(j);                           // xcoord
                hit_data[1] = Hv->GetValue(j);                           // ycoord
                hit_data[2] = static_cast<double>(Hpk->GetValue(j) - 1); // plane
                event_data[j] = hit_data;
            }
            events_container_[total_events_] = event_data;
            total_events_++;
        }
    }
    events_container_.resize(total_events_);
    cout << "Number of usable events :\t" << total_events_ << "\n\n";
}

// Check if point lies within boundaries of a polygon
// Input:
//          vertices            : number of vertices
//          vertex_x_coordinate : x coordinate of vertices
//          vertex_y_coordinate : y coordinate of vertices
//          x_coordinate        : x coordinate of point
//          y_coordinate        : y coordinate of point
int analyser::IsInsidePolygon(int vertices, vector<double> vertex_x_coordinate, vector<double> vertex_y_coordinate, double x_coordinate, double y_coordinate)
{
    int i, j, c = 0;
    for (i = 0, j = vertices - 1; i < vertices; j = i++)
    {
        if (((vertex_y_coordinate[i] > y_coordinate) != (vertex_y_coordinate[j] > y_coordinate)) && (x_coordinate < (vertex_x_coordinate[j] - vertex_x_coordinate[i]) * (y_coordinate - vertex_y_coordinate[i]) / (vertex_y_coordinate[j] - vertex_y_coordinate[i]) + vertex_x_coordinate[i]))
        {
            c = !c;
        }
    }
    return c;
}
