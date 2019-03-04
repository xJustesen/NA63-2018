#include "analyser.h"

analyser::analyser(vector<double> z, const char *name, string run, string beamparams)
{
    /* Paramters for MIMOSA detectors */
    row_count_ = 576;
    column_count_ = 1152;
    detector_xcoord_min_ = -11000, detector_xcoord_max_ = 11000, detector_ycoord_min_ = -5500, detector_ycoord_max_ = 5500;
    detector_zcoord_ = z;                                                                                          // z-coords of planes
    M1_M2_zcoord_ = {detector_zcoord_[0], detector_zcoord_[1], detector_zcoord_[2]};                               // z-coords for M1-M3 PLANES
    M2_M3_zcoord_ = {detector_zcoord_[1], detector_zcoord_[2], detector_zcoord_[3]};                               // z-coords for M2-M4 planes
    M3_M4_zcoord_ = {detector_zcoord_[2], detector_zcoord_[3], (detector_zcoord_[4] + detector_zcoord_[3]) / 2.0}; // z-coords of M3-MM planes
    M6_M5_zcoord_ = {detector_zcoord_[5], detector_zcoord_[4], (detector_zcoord_[4] + detector_zcoord_[3]) / 2.0}; // z-coords for M5-MM planes

    /* Conditions for track construction + pairing */
    M1_M2_proj_lim_ = 1400.0;       // acceptable distance between M1->M2 projected point and observed hit in M3 [micro-meter]
    M2_M3_proj_lim_ = 400.0;        // acceptable distance between M2->M3 projected point and observed hit in M4 [micro-meter]
    M6_M5_proj_lim_ = 200.0;        // acceptable distance in MM bewteen M1->MM and M6->MM arm of track [micro-meter]
    MM_paired_tracks_lim_ = 55.0e6; // matching ditance in MM between electron/positron paired track [micro-meter]
    foil_paired_tracks_lim_ = 55.0; // matching ditance in foil between electron/positron paired track [micro-meter]

    /* Alignment conditions */
    alignment_radius_lim_ = {5000, 1300, 600, 200, 60}; // radius for acceptable hits when aligning [micro-meter]
    T_matrix_convergence_tol_ = 1E-9;                   // convergence tolerance

    /* Info for data storage on disk */
    run_number_ = run;
    file_name_ = name;
    beam_parameters_file_name_ = beamparams;
    data_path_ = "/home/christian/Documents/cern2018/simdata";
    is_data_run_ = false;

    /* Vectors for data storage in RAM */
    hits_container_.resize(6);
    pixel_grids_.resize(6);
    hot_pixels_.resize(6);
    T.set_size(3, 3, 4);
    interplanar_distance_.resize(4);

    cout.setf(ios::fixed, ios::floatfield);
    cout.precision(3);
}

vector<vector<double>> analyser::GetEnergies()
{
    return energies;
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

void analyser::PrintVector(string name, vector<double> data)
{
    string file_name = data_path_ + "/" + name + ".txt";
    ofstream output(file_name);

    for (size_t i = 0; i < data.size(); i++)
    {
        output << data[i] << "\n";
    }
}

void analyser::PrintVector(string name, vector<vector<double>> data)
{
    string file_name = data_path_ + "/" + name + ".txt";
    ofstream output(file_name);

    for (size_t i = 0; i < data.size(); i++)
    {
        for (size_t j = 0; j < data[i].size(); j++)
        {
            output << data[i][j] << "\t";
        }

        output << "\n";
    }
}

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

void analyser::PrintHits(string name)
{
    string file_name = data_path_ + "/hits_coord_data_" + name + ".txt";
    ofstream output(file_name);

    for (size_t i = 0; i < hits_container_.size(); i++)
    {
        output << "## PLANE\t" << i << "\tHIT DATA\n";

        for (size_t j = 0; j < hits_container_[i].size(); j++)
        {
            for (size_t k = 0; k < hits_container_[i][j].size(); k++)
            {
                output << hits_container_[i][j][k][0] << ' ' << hits_container_[i][j][k][1] << '\n';
            }
        }
    }
}

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
}

void analyser::PrintSlope(string name)
{
    string file_name = data_path_ + "/angles_" + name + ".txt";
    ofstream output(file_name);

    for (size_t i = 0; i < MM_slopes_.size(); i++)
    {
        output << MM_slopes_[i][0] << " " << MM_slopes_[i][1] << "\n";
    }
}

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

double analyser::CalculateSlope(vec x, vec y)
{
    return (y(1) - y(0)) / (x(1) - x(0));
}

void analyser::ImageCrystal(string name)
{
    string file_name = data_path_ + "/crystal_image_" + name + ".txt";
    ofstream output(file_name);

    for (size_t i = 0; i < paired_tracks_.size(); i++)
    {

        for (size_t j = 0; j < paired_tracks_[i].size(); j++)
        {

            output << paired_tracks_[i][j][0][1][0] << "\t" << paired_tracks_[i][j][0][1][1] << "\n";
        }
    }
}

/* Calculate distances between hits in plane */
void analyser::PrintInterplanarDistance(int plane, string name)
{
    string plane_str = to_string(plane);
    string file_name = data_path_ + "/plane" + plane_str + "_" + name + "_dist.txt";
    ofstream output(file_name);

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
}

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

/* Pair tracks determined by "ConstructTracks" method */
void analyser::PairTracks(void)
{
    int pairs = 0;

    paired_tracks_.resize(tracks_.size());
    energies.resize(tracks_.size());
    vector<int> bad_events(tracks_.size());

    /*
           Here we find pairs. We only want a track to pair exactly one other track.
           To ensure this we iterate first from i = 0 -> number of tracks, then j = i + 1 -> number of tracks.
           If multiple pairs are found we erase them from the tracks array.
         */

    int photons_total_count = 0;

    for (int i = 0; i < total_events_; i++)
    {

        int photons_event_count = 0;

        for (size_t j = 0; j < tracks_[i].size(); j++)
        {

            int matched_tracks = 0;

            for (size_t k = j + 1; k < tracks_[i].size(); k++)
            {
                double distance = CalculateDistance(tracks_[i][j][4][0], tracks_[i][j][4][1], tracks_[i][k][4][0], tracks_[i][k][4][1]);

                if (distance < MM_paired_tracks_lim_)
                {
                    double M3_distance = CalculateDistance(tracks_[i][j][2][0], tracks_[i][j][2][1], tracks_[i][k][2][0], tracks_[i][k][2][1]);

                    if (M3_distance < foil_paired_tracks_lim_)
                    {
                        photons_event_count++;
                        matched_tracks++;
                        pairs++;

                        paired_tracks_[i].push_back({tracks_[i][j], tracks_[i][k]});
                        energies[i].push_back(CalculatePairEnergy({tracks_[i][j], tracks_[i][k]}));

                        /* Remove additional matched track since we do not want it to potentially match any other track */
                        if (matched_tracks > 1)
                        {
                            tracks_[i].erase(tracks_[i].begin() + k);
                            k--;
                        }
                    }
                }

            } // tracks #2 done

        } // tracks #1 done

        if (photons_event_count > 1)
        {
            bad_events[photons_total_count] = i;
            photons_total_count++;
        }
    }

    bad_events.resize(photons_total_count);

    /* Discard photons from bad events */
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

void analyser::BeamDivergence(int event_index, int hit_indx_0, int hit_indx_1)
{
    vec M1_M2_x = {hits_container_[0][event_index][hit_indx_0][0], hits_container_[1][event_index][hit_indx_1][0]};
    vec M1_M2_y = {hits_container_[0][event_index][hit_indx_0][1], hits_container_[1][event_index][hit_indx_1][1]};
    vec vec_M1_M2_z = {M1_M2_zcoord_[0], M1_M2_zcoord_[1]};

    vector<double> ang_pos;

    double ang_x = CalculateSlope(vec_M1_M2_z, M1_M2_x);
    double ang_y = CalculateSlope(vec_M1_M2_z, M1_M2_y);

    ang_pos = {ang_x, ang_y, hits_container_[0][event_index][hit_indx_0][0], hits_container_[0][event_index][hit_indx_0][1]};
    M1_M2_slopes_[event_index].push_back(ang_pos);
}

/* Construct particle tracks from M1 -> M6 */
void analyser::ConstructTracks(double M1_M2_slope_lb_x, double M1_M2_slope_ub_x, double M1_M2_slope_lb_y, double M1_M2_slope_ub_y)
{
    int M1_M4_tot_tracks = 0, M1_M6_tot_tracks = 0, M1_M3_tot_tracks = 0, M5_M6_tot_tracks = 0;
    int events_within_cut = 0;

    tracks_.clear();
    tracks_.resize(total_events_);
    M1_M2_slopes_.resize(total_events_);

    vector<double> M1_M2_dist(10 * total_events_);
    vector<double> M2_M3_dist(10 * total_events_);
    vector<double> M5_M6_dist(10 * total_events_);
    vector<double> yz_defl(10 * total_events_);

    int progress = 0;
    int M1_M2_count = 0;
    int M2_M3_count = 0;
    int M5_M6_count = 0;
    int yz_defl_count = 0;

    double start_time = omp_get_wtime();
    double T = 0;

#pragma omp parallel for
    for (int i = 0; i < total_events_; i++)
    {
        vector<vector<vector<double>>> M1_MM_tracks(hits_container_[0][i].size() * hits_container_[1][i].size() * hits_container_[2][i].size() * hits_container_[3][i].size());
        int M1_M4_tracks = 0;
        vector<vector<double>> M1_M6_track;
        bool included = false;

        /*
        Here we construct tracks. Each loop iterates over hits in a mimosa detector (0 -> 5).
        Conditions to find tracks checked periodically. All acceptable tracks are saved, which
        means bad tracks are saved. These should later be discarded in "PairTracks" function
        */
        for (size_t j = 0; j < hits_container_[0][i].size(); j++)
        {
            for (size_t k = 0; k < hits_container_[1][i].size(); k++)
            {
                vector<double> M1_M2_Proj = RectilinearProjection(hits_container_[0][i][j], hits_container_[1][i][k], M1_M2_zcoord_);
                BeamDivergence(i, j, k);
                /* Only pick out certain angles of entry into crystal */
                if (M1_M2_slopes_[i].back()[0] > M1_M2_slope_lb_x and M1_M2_slopes_[i].back()[0] < M1_M2_slope_ub_x and M1_M2_slopes_[i].back()[1] > M1_M2_slope_lb_y and M1_M2_slopes_[i].back()[1] < M1_M2_slope_ub_y)
                {
                    if (!included)
                    {
#pragma omp atomic
                        events_within_cut++;
                        included = true;
                    }
                    for (size_t l = 0; l < hits_container_[2][i].size(); l++)
                    {
                        double M1M2_d = CalculateDistance(M1_M2_Proj[0], M1_M2_Proj[1], hits_container_[2][i][l][0], hits_container_[2][i][l][1]);
                        M1_M2_dist[M1_M2_count] = M1M2_d;
                        M1_M2_count++;
                        if (M1M2_d < M1_M2_proj_lim_)
                        {
                            vector<double> M2M3_Proj = RectilinearProjection(hits_container_[1][i][k], hits_container_[2][i][l], M2_M3_zcoord_);
#pragma omp atomic
                            M1_M3_tot_tracks++;
                            for (size_t m = 0; m < hits_container_[3][i].size(); m++)
                            {
                                double M2M3_d = CalculateDistance(M2M3_Proj[0], M2M3_Proj[1], hits_container_[3][i][m][0], hits_container_[3][i][m][1]);
                                M2_M3_dist[M2_M3_count] = M2M3_d;
                                M2_M3_count++;
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
        /* Iterate over M6 -> M6 tracks */
        for (size_t n = 0; n < hits_container_[5][i].size(); n++)
        {
            vector<vector<double>> M1_M6_track;
            double shortest_distance = INFINITY;
            bool track_found = false;
            for (size_t o = 0; o < hits_container_[4][i].size(); o++)
            {
                vector<double> M6M5_Proj = RectilinearProjection(hits_container_[5][i][n], hits_container_[4][i][o], M6_M5_zcoord_);
                M5_M6_tot_tracks++;
                /* For every M5-M6 combination, pick an M3M4 projection */
                for (size_t l = 0; l < M1_MM_tracks.size(); l++)
                {
                    vector<double> M3M4_Proj = M1_MM_tracks[l].back();
                    double M6M5_d = CalculateDistance(M6M5_Proj[0], M6M5_Proj[1], M3M4_Proj[0], M3M4_Proj[1]);
                    M5_M6_dist[M5_M6_count] = M6M5_d;
                    M5_M6_count++;
                    vec vec_M3M4_z = {M3_M4_zcoord_[0], M3_M4_zcoord_[1]};
                    vec M3M4_y = {M1_MM_tracks[l][2][1], M1_MM_tracks[l][3][1]};
                    double M3M4_slope = CalculateSlope(vec_M3M4_z, M3M4_y);
                    if (M6M5_d < shortest_distance and M6M5_d < M6_M5_proj_lim_)
                    {
                        shortest_distance = M6M5_d;
                        vec vec_M5M6_z = {M6_M5_zcoord_[0], M6_M5_zcoord_[1]};
                        vec M5M6_y = {hits_container_[5][i][n][1], hits_container_[4][i][o][1]};
                        double M5M6_slope = CalculateSlope(vec_M5M6_z, M5M6_y);
#pragma omp critical
                        {
                            yz_defl[yz_defl_count] = M5M6_slope - M3M4_slope;
                            yz_defl_count++;
                        }
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

        /* Report progress in terminal */
        if ((progress + 1) % (total_events_ / 10 + 1) == 0)
        {
            double dt = omp_get_wtime() - start_time;
            T += dt;
            cout << "Progress :\t" << floor(100 * double(progress) / (double)total_events_) << "%"
                 << "\ttime used :\t" << dt << "\ttotal time elapsed :\t" << T << "\ttime remaining :\t" << dt * (double)total_events_ / (total_events_ / 10 + 1) - T << "\n";
            start_time = omp_get_wtime();
        }
    } // events done
    string eventfile_name;
    if (is_data_run_)
    {
        eventfile_name = "events_run_cuts.txt";
        fstream event_number_stream(eventfile_name, fstream::in | fstream::out | fstream::app);

        if (event_number_stream.is_open())
        {
            event_number_stream << run_number_ << "\t" << events_within_cut << "\n";
            event_number_stream.close();
        }
        else
            cerr << "Cannot open file:\t" << eventfile_name << "\n\n";
    }
    else
    {
        eventfile_name = "events_run_cuts_" + run_number_ + ".txt";
        ifstream input(eventfile_name);
        int events = 0;

        if (input.is_open())
        {
            while (true)
            {
                input >> events;
                if (input.eof())
                    break;
            }

            input.close();
        }
        else
            cerr << "Unable to open file containing no. of events\n";

        ofstream output(eventfile_name);
        output << events_within_cut + events;
        output.close();
    }

    string name0 = "M1M2_dist_" + (string)run_number_;
    string name1 = "M2M3_dist_" + (string)run_number_;
    string name2 = "M5M6_dist_" + (string)run_number_;
    string name3 = "yz_defl_" + (string)run_number_;
    string name4 = "beam_divergence_" + (string)run_number_;

    M1_M2_dist.resize(M1_M2_count);
    M2_M3_dist.resize(M2_M3_count);
    M5_M6_dist.resize(M5_M6_count);
    yz_defl.resize(yz_defl_count);

    cout << "\nTotal tracks from M1 -> M3:\t" << M1_M3_tot_tracks << "\n";
    cout << "Total tracks from M1 -> M4:\t" << M1_M4_tot_tracks << "\n";
    cout << "Total tracks from M6 -> M5:\t" << M5_M6_tot_tracks << "\n";
    cout << "Total tracks from M1 -> M6:\t" << M1_M6_tot_tracks << "\n";
    cout << "Number of events with tracks within cut:\t" << events_within_cut << "\n\n";
}

void analyser::FindAxisCounts(void)
{
    /* Load beam-parameters from file */
    vector<double> params;
    double param;
    string file = (string)beam_parameters_file_name_;
    ifstream paramters(file);

    if (paramters.is_open())
    {
        while (true)
        {
            paramters >> param;
            cout << param << "\n";
            if (paramters.eof())
                break;

            params.push_back(param);
        }

        paramters.close();
    }
    else
        cout << "Unable to open file containing beam parameters";

    vec vec_M1_M2_z = {M1_M2_zcoord_[0], M1_M2_zcoord_[1]};
    vec vec_M2_M3_z = {M2_M3_zcoord_[0], M2_M3_zcoord_[1]};

    double beam_x_direction = params[2];
    double beam_y_direction = params[3];
    double scan_x_range = 500e-6;
    double scan_y_range = 500e-6;
    double scan_window_size = 30E-6;
    double scan_increment = 2.0E-6;

    double ycut_lb = beam_y_direction - scan_y_range - scan_window_size;
    double ycut_ub = beam_y_direction + scan_y_range + scan_window_size;
    double xcut_lb = beam_x_direction - scan_x_range - scan_window_size;
    double xcut_ub = beam_x_direction + scan_x_range + scan_window_size;

    int no_cuts, no_cuts_x = 2 * scan_x_range / scan_increment, no_cuts_y = 2 * scan_y_range / scan_increment;
    string file_name;

    vector<vector<vector<vector<double>>>> tracks_alt;
    tracks_alt.resize(total_events_);

    double M1_M2_slope_lb_x, M1_M2_slope_ub_x, M1_M2_slope_lb_y, M1_M2_slope_ub_y;
    vector<vector<double>> photons_in_cut_x;
    photons_in_cut_x.resize(no_cuts_x);
    cout << "no_cuts_y:\t" << no_cuts_y << "\tno_cuts_x:\t" << no_cuts_x << "\n";
    vector<vector<double>> photons_in_cut_y;
    photons_in_cut_y.resize(no_cuts_y);

    for (int direction = 0; direction < 2; direction++)
    {
        if (direction == 0)
        {
            file_name = "no_photons_x_alt" + (string)run_number_;
            no_cuts = no_cuts_x;
        }
        else
        {
            file_name = "no_photons_y_alt" + (string)run_number_;
            no_cuts = no_cuts_y;
        }

        int cut_count = 0;

        // #pragma omp parallel for
        for (int cut = 0; cut < no_cuts; cut++)
        {
            int events_within_cut = 0;
            int photons_event_count = 0;

            if (direction == 0)
            {
                M1_M2_slope_lb_x = xcut_lb;
                M1_M2_slope_ub_x = xcut_ub;
                M1_M2_slope_lb_y = -1e+17;
                M1_M2_slope_ub_y = 1e+17;
            }
            else
            {
                M1_M2_slope_lb_x = -1e+17;
                M1_M2_slope_ub_x = 1e+17;
                M1_M2_slope_lb_y = ycut_lb;
                M1_M2_slope_ub_y = ycut_ub;
            }

            for (int i = 0; i < total_events_; i++)
            {
                int track_count = 0;

                tracks_alt[i].resize(10 * hits_container_[5][i].size());

                vector<vector<vector<double>>> M1_MM_tracks;
                vector<vector<double>> M1_M6_track;

                for (size_t j = 0; j < hits_container_[0][i].size(); j++)
                {
                    for (size_t k = 0; k < hits_container_[1][i].size(); k++)
                    {
                        vector<double> M1_M2_Proj = RectilinearProjection(hits_container_[0][i][j], hits_container_[1][i][k], M1_M2_zcoord_);

                        vec M1_M2_x = {hits_container_[0][i][j][0], hits_container_[1][i][k][0]};
                        vec M1_M2_y = {hits_container_[0][i][j][1], hits_container_[1][i][k][1]};
                        double M1_M2_ang_x = CalculateSlope(vec_M1_M2_z, M1_M2_x);
                        double M1_M2_ang_y = CalculateSlope(vec_M1_M2_z, M1_M2_y);

                        /* Only pick out certain angles of entry into crystal */
                        if (M1_M2_ang_x > M1_M2_slope_lb_x and M1_M2_ang_x < M1_M2_slope_ub_x and M1_M2_ang_y > M1_M2_slope_lb_y and M1_M2_ang_y < M1_M2_slope_ub_y)
                        {
                            events_within_cut++;

                            for (size_t l = 0; l < hits_container_[2][i].size(); l++)
                            {
                                double M1M2_d = CalculateDistance(M1_M2_Proj[0], M1_M2_Proj[1], hits_container_[2][i][l][0], hits_container_[2][i][l][1]);

                                if (M1M2_d < M1_M2_proj_lim_)
                                {
                                    vector<double> M2M3_Proj = RectilinearProjection(hits_container_[1][i][k], hits_container_[2][i][l], M2_M3_zcoord_);

                                    for (size_t m = 0; m < hits_container_[3][i].size(); m++)
                                    {
                                        double M2M3_d = CalculateDistance(M2M3_Proj[0], M2M3_Proj[1], hits_container_[3][i][m][0], hits_container_[3][i][m][1]);

                                        if (M2M3_d < M2_M3_proj_lim_)
                                        {
                                            vector<double> M3M4_Proj = RectilinearProjection(hits_container_[2][i][l], hits_container_[3][i][m], M3_M4_zcoord_);
                                            M1_MM_tracks.push_back({hits_container_[0][i][j], hits_container_[1][i][k], hits_container_[2][i][l], hits_container_[3][i][m], M3M4_Proj});
                                        }
                                    }
                                }

                            } // hits in M3 done

                        } // end cut

                    } // hits in M2 done

                } // hits in M1 done

                /* Iterate over M5 -> M6 tracks */
                for (size_t n = 0; n < hits_container_[5][i].size(); n++)
                {
                    vector<vector<double>> M1_M6_track;
                    double shortest_distance = INFINITY;
                    bool track_found = false;

                    for (size_t o = 0; o < hits_container_[4][i].size(); o++)
                    {
                        vector<double> M6M5_Proj = RectilinearProjection(hits_container_[5][i][n], hits_container_[4][i][o], M6_M5_zcoord_);

                        /* For every M5-M6 combination, pick an M3M4 projection */
                        for (size_t l = 0; l < M1_MM_tracks.size(); l++)
                        {
                            vector<double> M3M4_Proj = M1_MM_tracks[l].back();
                            double M6M5_d = CalculateDistance(M6M5_Proj[0], M6M5_Proj[1], M3M4_Proj[0], M3M4_Proj[1]);

                            vec vec_M3M4_z = {M3_M4_zcoord_[0], M3_M4_zcoord_[1]};
                            vec M3M4_y = {M1_MM_tracks[l][2][1], M1_MM_tracks[l][3][1]};

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
                        tracks_alt[i][track_count] = M1_M6_track;
                        track_count++;
                    }

                } // hits in M5 done

                tracks_alt[i].resize(track_count);

            } // events done

            /* Pair tracks */
            for (size_t i = 0; i < tracks_alt.size(); i++)
            {

                for (size_t j = 0; j < tracks_alt[i].size(); j++)
                {

                    int matchedtracks = 0;

                    for (size_t k = j + 1; k < tracks_alt[i].size(); k++)
                    {
                        double dist = CalculateDistance(tracks_alt[i][j][4][0], tracks_alt[i][j][4][1], tracks_alt[i][k][4][0], tracks_alt[i][k][4][1]);

                        if (dist < MM_paired_tracks_lim_)
                        {
                            double M3_dist = CalculateDistance(tracks_alt[i][j][2][0], tracks_alt[i][j][2][1], tracks_alt[i][k][2][0], tracks_alt[i][k][2][1]);

                            if (M3_dist < foil_paired_tracks_lim_)
                            {
                                photons_event_count += 1;

                                /* Remove additional matched track since we do not want it to potentially match any other track */
                                if (matchedtracks > 1)
                                {
                                    tracks_alt[i].erase(tracks_alt[i].begin() + k);
                                    k--;
                                }
                            }
                        }
                    }
                }
            }

            /* Calculate number of photons within cut */
            if (direction == 0)
            {
                photons_in_cut_x[cut] = {xcut_lb + scan_increment / 2.0, (double)photons_event_count / events_within_cut};
                xcut_lb += scan_increment;
                xcut_ub += scan_increment;
            }
            else
            {
                photons_in_cut_y[cut] = {ycut_lb + scan_window_size / 2.0, (double)photons_event_count / events_within_cut};
                ycut_lb += scan_increment;
                ycut_ub += scan_increment;
            }

            cut_count++;

            if ((cut_count + 1) % (no_cuts / 10 + 1) == 0)
                cout << "Progress: \t" << 100 * cut_count / no_cuts << "%\n";

        } // cuts done

        if (direction == 0)
            PrintVector(file_name, photons_in_cut_x);
        else
            PrintVector(file_name, photons_in_cut_y);
    }
}

void analyser::FindAxisDeflection(void)
{
    /* Load beam-parameters from file */
    vector<double> params;
    double param;
    string file = (string)beam_parameters_file_name_;
    ifstream paramters(file);

    if (paramters.is_open())
    {
        while (true)
        {
            paramters >> param;
            cout << param << "\n";
            if (paramters.eof())
                break;

            params.push_back(param);
        }

        paramters.close();
    }
    else
        cout << "Unable to open file containing beam parameters" << endl;

    vec vec_M1_M2_z = {M1_M2_zcoord_[0], M1_M2_zcoord_[1]};
    vec vec_M2_M3_z = {M2_M3_zcoord_[0], M2_M3_zcoord_[1]};

    double beam_x_direction = params[2];
    double beam_y_direction = params[3];
    double scan_x_range = 500e-6;
    double scan_y_range = 500e-6;
    double scan_window_size = 30E-6;
    double scan_increment = 2.0E-6;

    double ycut_lb = beam_y_direction - scan_y_range - scan_window_size;
    double ycut_ub = beam_y_direction + scan_y_range + scan_window_size;
    double xcut_lb = beam_x_direction - scan_x_range - scan_window_size;
    double xcut_ub = beam_x_direction + scan_x_range + scan_window_size;

    int no_cuts, no_cuts_x = 2 * scan_x_range / scan_increment, no_cuts_y = 2 * scan_y_range / scan_increment;
    string file_name;

    vector<vector<double>> mean_deflection_x;
    mean_deflection_x.resize(no_cuts_x);
    vector<vector<double>> mean_deflection_y;
    mean_deflection_y.resize(no_cuts_y);
    cout << "no_cuts_y:\t" << no_cuts_y << "\tno_cuts_x:\t" << no_cuts_x << "\n";

    double M1_M2_slope_lb_x, M1_M2_slope_ub_x, M1_M2_slope_lb_y, M1_M2_slope_ub_y, deflection_sum_x = 0, deflection_sum_y = 0;

    for (int direction = 0; direction < 2; direction++)
    {
        if (direction == 0)
        {
            file_name = "angles_mean_x_" + (string)run_number_;
            no_cuts = no_cuts_x;
        }
        else
        {
            file_name = "angles_mean_y_" + (string)run_number_;
            no_cuts = no_cuts_y;
        }

        int cut_count = 0;

        // #pragma omp parallel for
        for (int cut = 0; cut < no_cuts; cut++)
        {
            if (direction == 0)
            {
                M1_M2_slope_lb_x = xcut_lb;
                M1_M2_slope_ub_x = xcut_ub;
                M1_M2_slope_lb_y = -1e+17;
                M1_M2_slope_ub_y = 1e+17;
                deflection_sum_x = 0;
            }
            else
            {
                M1_M2_slope_lb_x = -1e+17;
                M1_M2_slope_ub_x = 1e+17;
                M1_M2_slope_lb_y = ycut_lb;
                M1_M2_slope_ub_y = ycut_ub;
                deflection_sum_y = 0;
            }

            int count = 0;

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

                        /* Only pick out certain angles of entry into crystal */
                        if (M1_M2_ang_x > M1_M2_slope_lb_x and M1_M2_ang_x < M1_M2_slope_ub_x and M1_M2_ang_y > M1_M2_slope_lb_y and M1_M2_ang_y < M1_M2_slope_ub_y)
                        {
                            for (size_t l = 0; l < hits_container_[2][i].size(); l++)
                            {
                                double M1M2_d = CalculateDistance(M1_M2_Proj[0], M1_M2_Proj[1], hits_container_[2][i][l][0], hits_container_[2][i][l][1]);

                                if (M1M2_d < M1_M2_proj_lim_)
                                {
                                    vec M2_M3_x = {hits_container_[1][i][k][0], hits_container_[2][i][l][0]};
                                    vec M2_M3_y = {hits_container_[1][i][k][1], hits_container_[2][i][l][1]};
                                    double M2_M3_ang_x = CalculateSlope(vec_M2_M3_z, M2_M3_x);
                                    double M2_M3_ang_y = CalculateSlope(vec_M2_M3_z, M2_M3_y);

                                    /* Calculate sum of delta angles */
                                    if (direction == 0)
                                        deflection_sum_x += abs(M2_M3_ang_x - M1_M2_ang_x);
                                    else
                                        deflection_sum_y += abs(M2_M3_ang_y - M1_M2_ang_y);

                                    count++;
                                }
                            } // hits in M3 done
                        }     // end cut
                    }         // hits in M2 done
                }             // hits in M1 done
            }                 // events done

            /* Calculate mean of delta angles */
            if (direction == 0)
            {
                deflection_sum_x /= count;
                mean_deflection_x[cut] = {xcut_lb + scan_increment / 2.0, deflection_sum_x};
                xcut_lb += scan_increment;
                xcut_ub += scan_increment;
            }
            else
            {
                deflection_sum_y /= count;
                mean_deflection_y[cut] = {ycut_lb + scan_window_size / 2.0, deflection_sum_y};
                ycut_lb += scan_increment;
                ycut_ub += scan_increment;
            }

            cut_count++;

            if ((cut_count + 1) % (no_cuts / 10 + 1) == 0)
                cout << "Progress: \t" << 100 * cut_count / no_cuts << "%\n";
        } // cuts done

        if (direction == 0)
            PrintVector(file_name, mean_deflection_x);
        else
            PrintVector(file_name, mean_deflection_y);
    }
}

/* Calculate distances between projected and observes hits */
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
                    double distance = CalculateDistance(hitp[0], hitp[1], hits2[l][0], hits2[l][1]);
                    interdistance.push_back(distance);
                }
            }
        }
    }
    return interdistance;
}

double analyser::CalculateDistance(double x0, double y0, double x1, double y1)
{
    return sqrt((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0)); // pythagoras
}

vector<double> analyser::RectilinearProjection(vector<double> hit0, vector<double> hit1, vector<double> z)
{
    vector<double> projection(2);

    projection[0] = (z[2] - z[0]) * (hit0[0] - hit1[0]) / (z[0] - z[1]) + hit0[0];
    projection[1] = (z[2] - z[0]) * (hit0[1] - hit1[1]) / (z[0] - z[1]) + hit0[1];

    return projection;
}

/* Align detector relative to M1-M2 */
void analyser::AlignPlane(mat &mat_Tot, vector<vector<vector<double>>> Hits0, vector<vector<vector<double>>> Hits1, vector<vector<vector<double>>> &Hits2, vector<double> z)
{
    vector<double> alignment_radius_lim;
    if (z.back() == detector_zcoord_[2] or z.back() == detector_zcoord_[3])
    {
        alignment_radius_lim = {alignment_radius_lim_[0], alignment_radius_lim_[1], alignment_radius_lim_[2]};
    }
    else
        alignment_radius_lim = alignment_radius_lim_;
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

/* Calculate matrix for detector alignment */
mat analyser::ConstructTMatrix(vector<vector<vector<double>>> Hits0, vector<vector<vector<double>>> Hits1, vector<vector<vector<double>>> Hits2, double dr_crit, vector<double> z)
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
                    double dr = CalculateDistance(hitp[0], hitp[1], hits2[n][0], hits2[n][1]);

                    if (dr < dr_crit)
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

/* Transform coordinates using transform matrix */
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

/* Allocate hitcoordinates as a row in a matrix */
void analyser::SaveHit(vector<double> projected_hit, vector<double> observed_hit, mat &mat_expected, mat &mat_observed, int indx)
{
    rowvec vec_expected = {projected_hit[0], projected_hit[1], 1};
    rowvec vec_observed = {observed_hit[0], observed_hit[1], 1};
    mat_expected.row(indx) = vec_expected;
    mat_observed.row(indx) = vec_observed;
}

/* Check similarity between two matrices */
bool analyser::CheckMatrixConvergence(mat A, mat B)
{
    return sum(abs(vectorise(A) - vectorise(B))) < T_matrix_convergence_tol_;
}

void analyser::FillInterplanarDistanceContainer(void)
{
    for (size_t i = 0; i < 4; i++)
    {
        vector<double> zproj = {detector_zcoord_[i], detector_zcoord_[i + 1], detector_zcoord_[i + 2]};
        interplanar_distance_[i] = CalculateInterdistance(hits_container_[i], hits_container_[i + 1], hits_container_[i + 2], zproj);
    }
}

/* Align planes with T matrix */
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

        interplanar_distance_[i] = CalculateInterdistance(hits_container_[i], hits_container_[i + 1], hits_container_[i + 2], zproj);
        cout << "Aligned plane " << i + 2 << "\n";
    }
}

/* Align planes without T matrix */
void analyser::AlignWithoutTMatrix(void)
{
    for (int i = 0; i < 4; i++)
    {
        mat mat_T = eye<mat>(3, 3);
        vector<double> zproj = {detector_zcoord_[i], detector_zcoord_[i + 1], detector_zcoord_[i + 2]};
        AlignPlane(mat_T, hits_container_[i], hits_container_[i + 1], hits_container_[i + 2], zproj);
        interplanar_distance_[i] = CalculateInterdistance(hits_container_[i], hits_container_[i + 1], hits_container_[i + 2], zproj);
        cout << "Aligned plane " << i + 2 << "\n";
        T.slice(i) = mat_T;
    }

    T.save(data_path_ + "/Align/alignment_matrix.txt", arma_ascii);
}

/* Construct pixel grid */
void analyser::MakeGrid(vector<vector<double>> &pixel_grid)
{
    pixel_grid.resize(column_count_ * row_count_);

    for (size_t i = 0; i < pixel_grid.size(); i++)
    {
        pixel_grid[i] = {(double)i, (double)0};
    }
}

/* Save coordinates of hits based on which detector recorded it */
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

/* Determine pixel no from coordinates */
int analyser::Coord2Pixel(double hit_x_coordinate, double hit_y_coordinate)
{
    double dx = (detector_xcoord_max_ - detector_xcoord_min_) / column_count_;
    double dy = (detector_ycoord_max_ - detector_ycoord_min_) / row_count_;
    int colno = (hit_x_coordinate + detector_xcoord_max_) / dx;
    int rowno = (hit_y_coordinate + detector_ycoord_max_) / dy;

    return rowno * column_count_ + colno;
}

/* Count total number of hits across all events */
void analyser::CountHits(int &count, vector<vector<vector<double>>> hit_coordinates)
{
    count = 0;
    for (size_t i = 0; i < hit_coordinates.size(); i++)
    {
        count += hit_coordinates[i].size();
    }
}

/*  Criteria to sort vector by 2nd column, descending */
bool analyser::SortVectorDescending_(const vector<double> &p1, const vector<double> &p2)
{
    return p1[1] > p2[1];
}

void analyser::LocateHotPixels(vector<vector<double>> pixelgrid, vector<int> &hot_pixels_, int i)
{
    /*
        We need to identify hot pixels. We do this by sorting the pixelgrid by number of pixels, descendind,
        and then saving the pixelno. of the hotpixels.
    */

    sort(pixelgrid.begin(), pixelgrid.end(), SortVectorDescending_);
    double lim = 4E-04 * total_events_;

    for (size_t i = 0; i < pixelgrid.size(); i++)
    {
        if (pixelgrid[i][1] > lim)
        {
            hot_pixels_.push_back(pixelgrid[i][0]);
        }
        else
            break;
    }
}

/* We remove both the number of counts in pixelgrid for hotpixels, and erase hitcoordinates from hotpixels. */
void analyser::RemoveHotPixels(vector<vector<vector<double>>> &hit_coordinates, vector<vector<double>> &pixel_grid, vector<int> hot_pixels_)
{
    for (size_t i = 0; i < hot_pixels_.size(); i++)
    {
        int hotpix = hot_pixels_[i];
        pixel_grid[hotpix][1] = 0;
        for (size_t j = 0; j < hit_coordinates.size(); j++)
        {
            for (size_t k = 0; k < hit_coordinates[j].size(); k++)
            {
                double hit_x_coordinate = hit_coordinates[j][k][0];
                double hit_y_coordinate = hit_coordinates[j][k][1];
                int pixel = Coord2Pixel(hit_x_coordinate, hit_y_coordinate);
                if (pixel == hotpix)
                {
                    hit_coordinates[j].erase(hit_coordinates[j].begin() + k);
                    k--;
                }
            }
        }
    }
}

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

/* Extract raw data from root file, do pre-processing and save usable data in std::vectors */
void analyser::ExtractRootData(void)
{
    is_data_run_ = true;
    /* Open root file and obtain relevant data */
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

    /* Sort data into std::vector */
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

    /* Save no. of events in run on disk */
    string eventfile_name = "events_run.txt";
    fstream event_number_stream(eventfile_name, fstream::in | fstream::out | fstream::app);

    if (event_number_stream.is_open())
    {
        event_number_stream << run_number_ << "\t" << total_events_ << "\n";
        event_number_stream.close();
    }
    else
        cerr << "Cannot open file:\t" << eventfile_name << "\n\n";
}

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
