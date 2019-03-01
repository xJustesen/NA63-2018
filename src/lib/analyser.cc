#include "analyser.h"

analyser::analyser(vector<double> z, const char *name, string run, string beamparams)
{
    /* Paramters for MIMOSA detectors */
    nrows = 576;
    ncols = 1152;
    xmin = -11000, xmax = 11000, ymin = -5500, ymax = 5500;
    zplanes = z;                                                        // z-coords of planes
    M1M2_z = {zplanes[0], zplanes[1], zplanes[2]};                      // z-coords for M1-M3 PLANES
    M2M3_z = {zplanes[1], zplanes[2], zplanes[3]};                      // z-coords for M2-M4 planes
    M3M4_z = {zplanes[2], zplanes[3], (zplanes[4] + zplanes[3]) / 2.0}; // z-coords of M3-MM planes
    M6M5_z = {zplanes[5], zplanes[4], (zplanes[4] + zplanes[3]) / 2.0}; // z-coords for M5-MM planes

    /* Conditions for track construction + pairing */
    M1M2_d_lim = 1400.0;         // acceptable distance between M1->M2 projected point and observed hit in M3 [micro-meter]
    M2M3_d_lim = 400.0;          // acceptable distance between M2->M3 projected point and observed hit in M4 [micro-meter]
    M6M5_d_lim = 200.0;          // acceptable distance in MM bewteen M1->MM and M6->MM arm of track [micro-meter]
    Match_d = 55.0e6;            // matching ditance in MM between electron/positron paired track [micro-meter]
    Match_d_foil = 55.0;         // matching ditance in foil between electron/positron paired track [micro-meter]
    yz_defl_lim = 3.50E-3 * 1e9; // angle of acceptable y-deflction [rad]

    /* Alignment conditions */
    dr_crit_list = {5000, 1300, 600, 200, 60}; // radius for acceptable hits when aligning [micro-meter]
    tol = 1E-9;                                // convergence tolerance

    /* Info for data storage on disk */
    runno = run;
    filename = name;
    paramsfile = beamparams;
    DATPATH = "/home/christian/Documents/cern2018/simdata";
    data = false;

    /* Vectors for data storage in RAM */
    hitcoords.resize(6);
    pixelgrids.resize(6);
    hotpixels.resize(6);
    T.set_size(3, 3, 4);
    distarray.resize(4);

    cout.setf(ios::fixed, ios::floatfield);
    cout.precision(3);
}

vector<vector<double>> analyser::get_energies()
{
    return energies;
}

void analyser::update_pixelgrids(int plane, vector<vector<double>> pixelgrid)
{
    pixelgrids[plane] = pixelgrid;
}

void analyser::update_hitcoords(int plane, vector<vector<vector<double>>> hitcoord)
{
    hitcoords[plane] = hitcoord;
}

void analyser::update_hotpixels(int plane, vector<int> hotpixel)
{
    hotpixels[plane] = hotpixel;
}

void analyser::save_vector(string name, vector<double> data)
{
    string filename = DATPATH + "/" + name + ".txt";
    ofstream output(filename);

    for (size_t i = 0; i < data.size(); i++)
    {
        output << data[i] << "\n";
    }
}

void analyser::save_vector(string name, vector<vector<double>> data)
{
    string filename = DATPATH + "/" + name + ".txt";
    ofstream output(filename);

    for (size_t i = 0; i < data.size(); i++)
    {
        for (size_t j = 0; j < data[i].size(); j++)
        {
            output << data[i][j] << "\t";
        }

        output << "\n";
    }
}

void analyser::save_vector(string name, vector<vector<vector<double>>> data)
{
    string filename = DATPATH + "/" + name + ".txt";
    ofstream output(filename);

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

void analyser::print_pixels(string name)
{
    string filename = DATPATH + "/pixeldata_run" + name + ".txt";
    ofstream output(filename);

    for (size_t i = 0; i < pixelgrids.size(); i++)
    {
        output << "## PLANE\t" << i << "\tPIXEL DATA\n";

        for (size_t j = 0; j < pixelgrids[i].size(); j++)
        {
            output << pixelgrids[i][j][1] << "\n";
        }
    }
}

void analyser::print_hotpixels(string name)
{
    for (size_t i = 0; i < hotpixels.size(); i++)
    {
        string filename = DATPATH + "/hotpixels_run" + name + "_plane_" + to_string(i) + ".txt";
        ofstream output(filename);

        for (size_t j = 0; j < hotpixels[i].size(); j++)
        {
            output << hotpixels[i][j] << "\n";
        }
    }
}

void analyser::print_interdistance(string name)
{
    string filename = DATPATH + "/interdistance_data_" + name + ".txt";
    ;
    ofstream output(filename);

    for (size_t i = 0; i < distarray.size(); i++)
    {
        output << "## DATABLOCK\t" << i << "\tINTERDISTANCE DATA\n";

        for (size_t j = 0; j < distarray[i].size(); j++)
        {
            output << distarray[i][j] << "\n";
        }

        output << "\n\n";
    }
}

void analyser::print_hits(string name)
{
    string filename = DATPATH + "/hits_coord_data_" + name + ".txt";
    ofstream output(filename);

    for (size_t i = 0; i < hitcoords.size(); i++)
    {
        output << "## PLANE\t" << i << "\tHIT DATA\n";

        for (size_t j = 0; j < hitcoords[i].size(); j++)
        {
            for (size_t k = 0; k < hitcoords[i][j].size(); k++)
            {
                output << hitcoords[i][j][k][0] << ' ' << hitcoords[i][j][k][1] << '\n';
            }
        }
    }
}

void analyser::print_energy(string filename)
{
    ofstream output;
    output.open(filename, ios_base::app);

    for (size_t i = 0; i < energies.size(); i++)
    {
        for (size_t j = 0; j < energies[i].size(); j++)
        {
            output << energies[i][j] << "\n";
        }
    }
}

void analyser::print_slope(string name)
{
    string filename = DATPATH + "/angles_" + name + ".txt";
    ofstream output(filename);

    for (size_t i = 0; i < slopes.size(); i++)
    {
        output << slopes[i][0] << " " << slopes[i][1] << "\n";
    }
}

void analyser::print_M1M2_slope(string name)
{
    string filename = DATPATH + "/angles_M1M2_" + name + ".txt";
    ofstream output(filename);

    for (int i = 0; i < Nevents; i++)
    {
        for (size_t j = 0; j < M1M2_slopes[i].size(); j++)
        {
            output << M1M2_slopes[i][j][0] << "\t" << M1M2_slopes[i][j][1] << "\t" << M1M2_slopes[i][j][2] << "\t" << M1M2_slopes[i][j][3] << "\n";
        }
    }
}

void analyser::print_zpos(string name)
{
    string filename = DATPATH + "/zclosepos_" + name + ".txt";
    ofstream output(filename);

    for (size_t i = 0; i < zclosepos.size(); i++)
    {
        for (size_t j = 0; j < zclosepos[i].size(); j++)
        {
            output << zclosepos[i][j] << "\n";
        }
    }
}

double analyser::calc_ang(double m, double n)
{
    return n - m;
}

double analyser::calc_slope(vec x, vec y)
{
    return (y(1) - y(0)) / (x(1) - x(0));
}

mat analyser::lines3d_nearestpoints(vec A, vec B, vec C, vec D)
{
    vec d1 = B - A;
    vec d2 = D - C;
    vec n = cross(d1, d2);
    vec n1 = cross(d1, n);
    vec n2 = cross(d2, n);
    vec c1 = A + dot(C - A, n2) / dot(d1, n2) * d1;
    vec c2 = C + dot(A - C, n1) / dot(d2, n1) * d2;

    mat closepos(3, 2);
    closepos.col(0) = c1;
    closepos.col(1) = c2;

    return closepos;
}

void analyser::image_crystal(string name)
{
    string filename = DATPATH + "/crystal_image_" + name + ".txt";
    ofstream output(filename);

    for (size_t i = 0; i < paired_tracks.size(); i++)
    { // events

        for (size_t j = 0; j < paired_tracks[i].size(); j++)
        { // tracks

            output << paired_tracks[i][j][0][1][0] << "\t" << paired_tracks[i][j][0][1][1] << "\n";
        }
    }
}

/* Calculate distances between hits in plane */
void analyser::print_planar_dist(int plane, string name)
{
    string plane_str = to_string(plane);
    string filename = DATPATH + "/plane" + plane_str + "_" + name + "_dist.txt";
    ofstream output(filename);

    double x1, x2, y1, y2;

    for (int i = 0; i < Nevents; i++)
    {
        for (size_t j = 0; j < hitcoords[plane][i].size(); j++)
        {
            for (size_t k = j + 1; k < hitcoords[plane][i].size(); k++)
            {
                x1 = hitcoords[plane][i][j][0];
                x2 = hitcoords[plane][i][k][0];
                y1 = hitcoords[plane][i][j][1];
                y2 = hitcoords[plane][i][k][1];

                output << sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2)) << "\n";
            }
        }
    }
}

double analyser::calc_pair_energy(vector<vector<vector<double>>> pairedtracks)
{
    double L = 0.15, q = 1.6021766e-19, B = 0.12, c = 299792458; // SI units
    double energy = 0;

    /*
           The energy of the photon corresponds to the total energy of each particle in the pair.
           The energy of each particle is calculated seperately in this loop
         */

    for (int i = 0; i < 2; i++)
    {
        vec M3M4_x = {pairedtracks[i][2][0], pairedtracks[i][3][0]}; // M3, M4 x coord
        vec M3M4_z = {zplanes[2], zplanes[3]};                       // M3, M4 z coord
        vec M5M6_x = {pairedtracks[i][6][0], pairedtracks[i][7][0]}; // M5, M6 x coord
        vec M5M6_z = {zplanes[4], zplanes[5]};                       // M5 M6 z coord

        double m = calc_slope(M3M4_z, M3M4_x);
        double n = calc_slope(M5M6_z, M5M6_x);
        vector<double> slope = {m, n};
        slopes.push_back(slope);

        double ang = calc_ang(m, n);        // calculate deflection angle to find energy
        energy += q * c * L * B / abs(ang); // Joule
    }

    return energy;
}

/* Pair tracks determined by "construct_tracks" method */
void analyser::pair_tracks(void)
{
    int pairs = 0;

    paired_tracks.resize(tracks.size());
    energies.resize(tracks.size());
    zclosepos.resize(tracks.size());
    vector<int> bad_events(tracks.size());

    /*
           Here we find pairs. We only want a track to pair exactly one other track.
           To ensure this we iterate first from i = 0 -> number of tracks, then j = i + 1 -> number of tracks.
           If multiple pairs are found we erase them from the tracks array.
         */

    int N_photons = 0;

    for (int i = 0; i < Nevents; i++)
    { // no. of events

        int no_photons = 0;

        for (size_t j = 0; j < tracks[i].size(); j++)
        { // tracks in event

            int matchedtracks = 0;

            for (size_t k = j + 1; k < tracks[i].size(); k++)
            {
                double dist = calc_dist(tracks[i][j][4][0], tracks[i][j][4][1], tracks[i][k][4][0], tracks[i][k][4][1]);

                if (dist < Match_d)
                {
                    double M3_dist = calc_dist(tracks[i][j][2][0], tracks[i][j][2][1], tracks[i][k][2][0], tracks[i][k][2][1]);

                    if (M3_dist < Match_d_foil)
                    {
                        vec P1 = {tracks[i][j][6][0], tracks[i][j][6][1], M6M5_z[0]}, Q1 = {tracks[i][j][7][0], tracks[i][j][7][1], M6M5_z[1]};
                        vec P2 = {tracks[i][k][2][0], tracks[i][k][2][1], M3M4_z[0]}, Q2 = {tracks[i][k][3][0], tracks[i][k][3][1], M6M5_z[1]};
                        mat closepos = lines3d_nearestpoints(P1, Q1, P2, Q2);
                        double zclosepospair = (closepos(2, 0) + closepos(2, 1)) / 2;
                        zclosepos[i].push_back(zclosepospair);

                        if (1 /*zclosepospair > 8e6 && zclosepospair < 1e7*/)
                        {
                            no_photons++;
                            matchedtracks++;
                            pairs++;

                            paired_tracks[i].push_back({tracks[i][j], tracks[i][k]});
                            energies[i].push_back(calc_pair_energy({tracks[i][j], tracks[i][k]}));

                            /* Remove additional matched track since we do not want it to potentially match any other track */
                            if (matchedtracks > 1)
                            {
                                tracks[i].erase(tracks[i].begin() + k);
                                k--;
                            }
                        }
                    }
                }

            } // tracks #2 done

        } // tracks #1 done

        if (no_photons > 1)
        {
            bad_events[N_photons] = i;
            N_photons++;
        }

    } // events done

    bad_events.resize(N_photons);

    /* Discard photons from bad events */
    int erased_events = 0;
    int photons_removed = 0;
    for (size_t i = 0; i < bad_events.size(); i++)
    {
        photons_removed += energies[bad_events[i] - erased_events].size();
        energies.erase(energies.begin() + bad_events[i] - erased_events);
        paired_tracks.erase(paired_tracks.begin() + bad_events[i] - erased_events);
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

void analyser::beam_divergence(int eventno, int hit_indx_0, int hit_indx_1, string name)
{
    string str1 = "full";
    string str2 = "partial";

    vec M1M2_x = {hitcoords[0][eventno][hit_indx_0][0], hitcoords[1][eventno][hit_indx_1][0]};
    vec M1M2_y = {hitcoords[0][eventno][hit_indx_0][1], hitcoords[1][eventno][hit_indx_1][1]};
    vec vec_M1M2_z = {M1M2_z[0], M1M2_z[1]};

    vector<double> ang_pos;

    double ang_x = calc_slope(vec_M1M2_z, M1M2_x);
    double ang_y = calc_slope(vec_M1M2_z, M1M2_y);

    if (name.compare(str1) == 0)
    {
        ang_pos = {ang_x, ang_y, hitcoords[0][eventno][hit_indx_0][0], hitcoords[0][eventno][hit_indx_0][1]};
        M1M2_slopes[eventno].push_back(ang_pos);
    }

    if (name.compare(str2) == 0)
    {
        ang_pos = {ang_x, ang_y, hitcoords[0][eventno][hit_indx_0][0], hitcoords[0][eventno][hit_indx_0][1]};
        divergence[eventno].push_back(ang_pos);
    }
}

/* Construct particle tracks from M1 -> M6 */
void analyser::construct_tracks(double M1M2_slope_lb_x, double M1M2_slope_ub_x, double M1M2_slope_lb_y, double M1M2_slope_ub_y)
{
    int M1_M4_tot_tracks = 0, M1_M6_tot_tracks = 0, M1_M3_tot_tracks = 0, M5_M6_tot_tracks = 0;
    int Nevents_within_cut = 0;

    divergence.resize(Nevents);

    string str1 = "full", str2 = "partial";

    tracks.clear();
    tracks.resize(Nevents);
    M1M2_slopes.resize(Nevents);

    vector<double> M1M2_dist(10 * Nevents);
    vector<double> M2M3_dist(10 * Nevents);
    vector<double> M5M6_dist(10 * Nevents);
    vector<double> yz_defl(10 * Nevents);

    int prog = 0;
    int M1M2_count = 0;
    int M2M3_count = 0;
    int M5M6_count = 0;
    int yz_defl_count = 0;

    double start_time = omp_get_wtime();
    double T = 0;

#pragma omp parallel for
    for (int i = 0; i < Nevents; i++)
    {
        vector<vector<vector<double>>> M1_MM_tracks(hitcoords[0][i].size() * hitcoords[1][i].size() * hitcoords[2][i].size() * hitcoords[3][i].size());
        int M1_M4_tracks = 0;
        vector<vector<double>> M1_M6_track;
        bool included = false;
        M1M2_slopes[i].resize(hitcoords[0][i].size() * hitcoords[1][i].size());

        /*
                   Here we construct tracks. Each loop iterates over hits in a mimosa detector (0 -> 5).
                   Conditions to find tracks checked periodically. All acceptable tracks are saved, which
                   means bad tracks are saved. These should later be discarded in "pair_tracks" function
                 */
        for (size_t j = 0; j < hitcoords[0][i].size(); j++)
        {
            for (size_t k = 0; k < hitcoords[1][i].size(); k++)
            {
                vector<double> M1M2_Proj = rect_project(hitcoords[0][i][j], hitcoords[1][i][k], M1M2_z);
                beam_divergence(i, j, k, str1);
                /* Only pick out certain angles of entry into crystal */
                if (M1M2_slopes[i].back()[0] > M1M2_slope_lb_x and M1M2_slopes[i].back()[0] < M1M2_slope_ub_x and M1M2_slopes[i].back()[1] > M1M2_slope_lb_y and M1M2_slopes[i].back()[1] < M1M2_slope_ub_y)
                {
                    if (!included)
                    {
#pragma omp atomic
                        Nevents_within_cut++;
                        included = true;
                    }
                    for (size_t l = 0; l < hitcoords[2][i].size(); l++)
                    {
                        double M1M2_d = calc_dist(M1M2_Proj[0], M1M2_Proj[1], hitcoords[2][i][l][0], hitcoords[2][i][l][1]);
                        M1M2_dist[M1M2_count] = M1M2_d;
                        M1M2_count++;
                        if (M1M2_d < M1M2_d_lim)
                        {
                            vector<double> M2M3_Proj = rect_project(hitcoords[1][i][k], hitcoords[2][i][l], M2M3_z);
#pragma omp atomic
                            M1_M3_tot_tracks++;
                            for (size_t m = 0; m < hitcoords[3][i].size(); m++)
                            {
                                double M2M3_d = calc_dist(M2M3_Proj[0], M2M3_Proj[1], hitcoords[3][i][m][0], hitcoords[3][i][m][1]);
                                M2M3_dist[M2M3_count] = M2M3_d;
                                M2M3_count++;
                                if (M2M3_d < M2M3_d_lim)
                                {
                                    vector<double> M3M4_Proj = rect_project(hitcoords[2][i][l], hitcoords[3][i][m], M3M4_z);
                                    M1_MM_tracks[M1_M4_tracks] = {hitcoords[0][i][j], hitcoords[1][i][k], hitcoords[2][i][l], hitcoords[3][i][m], M3M4_Proj};
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
        for (size_t n = 0; n < hitcoords[5][i].size(); n++)
        {
            vector<vector<double>> M1_M6_track;
            double shortdist = INFINITY;
            bool trackfound = false;
            for (size_t o = 0; o < hitcoords[4][i].size(); o++)
            {
                vector<double> M6M5_Proj = rect_project(hitcoords[5][i][n], hitcoords[4][i][o], M6M5_z);
                M5_M6_tot_tracks++;
                /* For every M5-M6 combination, pick an M3M4 projection */
                for (size_t l = 0; l < M1_MM_tracks.size(); l++)
                {
                    vector<double> M3M4_Proj = M1_MM_tracks[l].back();
                    double M6M5_d = calc_dist(M6M5_Proj[0], M6M5_Proj[1], M3M4_Proj[0], M3M4_Proj[1]);
                    M5M6_dist[M5M6_count] = M6M5_d;
                    M5M6_count++;
                    vec vec_M3M4_z = {M3M4_z[0], M3M4_z[1]};
                    vec M3M4_y = {M1_MM_tracks[l][2][1], M1_MM_tracks[l][3][1]};
                    double M3M4_slope = calc_slope(vec_M3M4_z, M3M4_y);
                    if (M6M5_d < shortdist and M6M5_d < M6M5_d_lim)
                    {
                        shortdist = M6M5_d;
                        vec vec_M5M6_z = {M6M5_z[0], M6M5_z[1]};
                        vec M5M6_y = {hitcoords[5][i][n][1], hitcoords[4][i][o][1]};
                        double M5M6_slope = calc_slope(vec_M5M6_z, M5M6_y);
#pragma omp critical
                        {
                            yz_defl[yz_defl_count] = M5M6_slope - M3M4_slope;
                            yz_defl_count++;
                        }
                        if (abs(M5M6_slope - M3M4_slope) < yz_defl_lim)
                        {
                            trackfound = true;
                            M1_M6_track = {M1_MM_tracks[l][0], M1_MM_tracks[l][1], M1_MM_tracks[l][2], M1_MM_tracks[l][3], M3M4_Proj, M6M5_Proj, hitcoords[4][i][o], hitcoords[5][i][n]};
                        }
                    }
                }
            } // hits in M6 done

            if (trackfound)
            {
#pragma omp critical
                {
                    tracks[i].push_back(M1_M6_track);
                    M1_M6_tot_tracks++;
                }
            } // if statement for push_back of tracks
        }     // hits in M5 done

#pragma omp atomic
        prog++;

        /* Report progress in terminal */
        if ((prog + 1) % (Nevents / 10 + 1) == 0)
        {
            double dt = omp_get_wtime() - start_time;
            T += dt;
            cout << "Progress :\t" << floor(100 * double(prog) / (double)Nevents) << "%"
                 << "\ttime used :\t" << dt << "\ttotal time elapsed :\t" << T << "\ttime remaining :\t" << dt * (double)Nevents / (Nevents / 10 + 1) - T << "\n";
            start_time = omp_get_wtime();
        }
    } // events done
    string eventfile_name;
    if (data)
    {
        eventfile_name = "events_run_cuts.txt";
        fstream event_number_stream(eventfile_name, fstream::in | fstream::out | fstream::app);

        if (event_number_stream.is_open())
        {
            event_number_stream << runno << "\t" << Nevents_within_cut << "\n";
            event_number_stream.close();
        }
        else
            cerr << "Cannot open file:\t" << eventfile_name << "\n\n";
    }
    else
    {
        eventfile_name = "events_run_cuts_" + runno + ".txt";
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
        output << Nevents_within_cut + events;
        output.close();
    }

    string name0 = "M1M2_dist_" + (string)runno;
    string name1 = "M2M3_dist_" + (string)runno;
    string name2 = "M5M6_dist_" + (string)runno;
    string name3 = "yz_defl_" + (string)runno;
    string name4 = "beam_divergence_" + (string)runno;

    M1M2_dist.resize(M1M2_count);
    M2M3_dist.resize(M2M3_count);
    M5M6_dist.resize(M5M6_count);
    yz_defl.resize(yz_defl_count);

    cout << "\nTotal tracks from M1 -> M3:\t" << M1_M3_tot_tracks << "\n";
    cout << "Total tracks from M1 -> M4:\t" << M1_M4_tot_tracks << "\n";
    cout << "Total tracks from M6 -> M5:\t" << M5_M6_tot_tracks << "\n";
    cout << "Total tracks from M1 -> M6:\t" << M1_M6_tot_tracks << "\n";
    cout << "Number of events with tracks within cut:\t" << Nevents_within_cut << "\n\n";
}

void analyser::find_axis_alt(void)
{
    /* Load beam-parameters from file */
    vector<double> params;
    double param;
    string file = (string)paramsfile;
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

    vec vec_M1M2_z = {M1M2_z[0], M1M2_z[1]};
    vec vec_M2M3_z = {M2M3_z[0], M2M3_z[1]};

    double xmean = params[2];
    double ymean = params[3];
    double deltax = 500e-6;
    double deltay = 500e-6;
    double theta = 30E-6;
    double dtheta = 2.0E-6;

    double ycut_lb = ymean - deltay - theta;
    double ycut_ub = ymean + deltay + theta;
    double xcut_lb = xmean - deltax - theta;
    double xcut_ub = xmean + deltax + theta;

    int no_cuts, no_cuts_x = 2 * deltax / dtheta, no_cuts_y = 2 * deltay / dtheta;
    string filename;

    vector<vector<vector<vector<double>>>> tracks_alt;
    tracks_alt.resize(Nevents);

    double M1M2_slope_lb_x, M1M2_slope_ub_x, M1M2_slope_lb_y, M1M2_slope_ub_y;
    vector<vector<double>> photons_in_cut_x;
    photons_in_cut_x.resize(no_cuts_x);
    cout << "no_cuts_y:\t" << no_cuts_y << "\tno_cuts_x:\t" << no_cuts_x << "\n";
    vector<vector<double>> photons_in_cut_y;
    photons_in_cut_y.resize(no_cuts_y);

    for (int dir = 0; dir < 2; dir++)
    {
        if (dir == 0)
        {
            filename = "no_photons_x_alt" + (string)runno;
            no_cuts = no_cuts_x;
        }
        else
        {
            filename = "no_photons_y_alt" + (string)runno;
            no_cuts = no_cuts_y;
        }

        int cutno = 0;

        // #pragma omp parallel for
        for (int cut = 0; cut < no_cuts; cut++)
        {
            int Nevents_within_cut = 0;
            int no_photons = 0;

            if (dir == 0)
            {
                M1M2_slope_lb_x = xcut_lb;
                M1M2_slope_ub_x = xcut_ub;
                M1M2_slope_lb_y = -1e+17;
                M1M2_slope_ub_y = 1e+17;
            }
            else
            {
                M1M2_slope_lb_x = -1e+17;
                M1M2_slope_ub_x = 1e+17;
                M1M2_slope_lb_y = ycut_lb;
                M1M2_slope_ub_y = ycut_ub;
            }

            for (int i = 0; i < Nevents; i++)
            {
                int track_count = 0;

                tracks_alt[i].resize(10 * hitcoords[5][i].size());

                vector<vector<vector<double>>> M1_MM_tracks;
                vector<vector<double>> M1_M6_track;

                for (size_t j = 0; j < hitcoords[0][i].size(); j++)
                {
                    for (size_t k = 0; k < hitcoords[1][i].size(); k++)
                    {
                        vector<double> M1M2_Proj = rect_project(hitcoords[0][i][j], hitcoords[1][i][k], M1M2_z);

                        vec M1M2_x = {hitcoords[0][i][j][0], hitcoords[1][i][k][0]};
                        vec M1M2_y = {hitcoords[0][i][j][1], hitcoords[1][i][k][1]};
                        double ang_x_M1M2 = calc_slope(vec_M1M2_z, M1M2_x);
                        double ang_y_M1M2 = calc_slope(vec_M1M2_z, M1M2_y);

                        /* Only pick out certain angles of entry into crystal */
                        if (ang_x_M1M2 > M1M2_slope_lb_x and ang_x_M1M2 < M1M2_slope_ub_x and ang_y_M1M2 > M1M2_slope_lb_y and ang_y_M1M2 < M1M2_slope_ub_y)
                        {
                            Nevents_within_cut++;

                            for (size_t l = 0; l < hitcoords[2][i].size(); l++)
                            {
                                double M1M2_d = calc_dist(M1M2_Proj[0], M1M2_Proj[1], hitcoords[2][i][l][0], hitcoords[2][i][l][1]);

                                if (M1M2_d < M1M2_d_lim)
                                {
                                    vector<double> M2M3_Proj = rect_project(hitcoords[1][i][k], hitcoords[2][i][l], M2M3_z);

                                    for (size_t m = 0; m < hitcoords[3][i].size(); m++)
                                    {
                                        double M2M3_d = calc_dist(M2M3_Proj[0], M2M3_Proj[1], hitcoords[3][i][m][0], hitcoords[3][i][m][1]);

                                        if (M2M3_d < M2M3_d_lim)
                                        {
                                            vector<double> M3M4_Proj = rect_project(hitcoords[2][i][l], hitcoords[3][i][m], M3M4_z);
                                            M1_MM_tracks.push_back({hitcoords[0][i][j], hitcoords[1][i][k], hitcoords[2][i][l], hitcoords[3][i][m], M3M4_Proj});
                                        }
                                    }
                                }

                            } // hits in M3 done

                        } // end cut

                    } // hits in M2 done

                } // hits in M1 done

                /* Iterate over M5 -> M6 tracks */
                for (size_t n = 0; n < hitcoords[5][i].size(); n++)
                {
                    vector<vector<double>> M1_M6_track;
                    double shortdist = INFINITY;
                    bool trackfound = false;

                    for (size_t o = 0; o < hitcoords[4][i].size(); o++)
                    {
                        vector<double> M6M5_Proj = rect_project(hitcoords[5][i][n], hitcoords[4][i][o], M6M5_z);

                        /* For every M5-M6 combination, pick an M3M4 projection */
                        for (size_t l = 0; l < M1_MM_tracks.size(); l++)
                        {
                            vector<double> M3M4_Proj = M1_MM_tracks[l].back();
                            double M6M5_d = calc_dist(M6M5_Proj[0], M6M5_Proj[1], M3M4_Proj[0], M3M4_Proj[1]);

                            vec vec_M3M4_z = {M3M4_z[0], M3M4_z[1]};
                            vec M3M4_y = {M1_MM_tracks[l][2][1], M1_MM_tracks[l][3][1]};
                            // double M3M4_slope = calc_slope(vec_M3M4_z, M3M4_y);

                            if (M6M5_d < shortdist and M6M5_d < M6M5_d_lim)
                            {
                                shortdist = M6M5_d;
                                vec vec_M5M6_z = {M6M5_z[0], M6M5_z[1]};
                                vec M5M6_y = {hitcoords[5][i][n][1], hitcoords[4][i][o][1]};
                                // double M5M6_slope = calc_slope(vec_M5M6_z, M5M6_y);

                                if (1 /*abs(M5M6_slope - M3M4_slope) < yz_defl_lim*/)
                                {
                                    trackfound = true;
                                    M1_M6_track = {M1_MM_tracks[l][0], M1_MM_tracks[l][1], M1_MM_tracks[l][2], M1_MM_tracks[l][3], M3M4_Proj, M6M5_Proj, hitcoords[4][i][o], hitcoords[5][i][n]};
                                }
                            }
                        }

                    } // hits in M6 done

                    if (trackfound)
                    {
                        tracks_alt[i][track_count] = M1_M6_track;
                        track_count++;
                    }

                } // hits in M5 done

                tracks_alt[i].resize(track_count);

            } // events done

            /* Pair tracks */
            for (size_t i = 0; i < tracks_alt.size(); i++)
            { // no. of events

                for (size_t j = 0; j < tracks_alt[i].size(); j++)
                { // tracks in event

                    int matchedtracks = 0;

                    for (size_t k = j + 1; k < tracks_alt[i].size(); k++)
                    {
                        double dist = calc_dist(tracks_alt[i][j][4][0], tracks_alt[i][j][4][1], tracks_alt[i][k][4][0], tracks_alt[i][k][4][1]);

                        if (dist < Match_d)
                        {
                            double M3_dist = calc_dist(tracks_alt[i][j][2][0], tracks_alt[i][j][2][1], tracks_alt[i][k][2][0], tracks_alt[i][k][2][1]);

                            if (M3_dist < Match_d_foil)
                            {
                                no_photons += 1;

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
            if (dir == 0)
            {
                photons_in_cut_x[cut] = {xcut_lb + dtheta / 2.0, (double)no_photons / Nevents_within_cut};
                xcut_lb += dtheta;
                xcut_ub += dtheta;
            }
            else
            {
                photons_in_cut_y[cut] = {ycut_lb + theta / 2.0, (double)no_photons / Nevents_within_cut};
                ycut_lb += dtheta;
                ycut_ub += dtheta;
            }

            cutno++;

            if ((cutno + 1) % (no_cuts / 10 + 1) == 0)
                cout << "Progress: \t" << 100 * cutno / no_cuts << "%\n";

        } // cuts done

        if (dir == 0)
            save_vector(filename, photons_in_cut_x);
        else
            save_vector(filename, photons_in_cut_y);
    }
}

void analyser::find_axis(void)
{
    /* Load beam-parameters from file */
    vector<double> params;
    double param;
    string file = (string)paramsfile;
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

    vec vec_M1M2_z = {M1M2_z[0], M1M2_z[1]};
    vec vec_M2M3_z = {M2M3_z[0], M2M3_z[1]};

    double xmean = params[2];
    double ymean = params[3];
    double deltax = 500e-6;
    double deltay = 500e-6;
    double theta = 30E-6;
    double dtheta = 2.0E-6;

    double ycut_lb = ymean - deltay - theta;
    double ycut_ub = ymean + deltay + theta;
    double xcut_lb = xmean - deltax - theta;
    double xcut_ub = xmean + deltax + theta;

    int no_cuts, no_cuts_x = 2 * deltax / dtheta, no_cuts_y = 2 * deltay / dtheta;
    string filename;

    vector<vector<double>> mean_anglesx;
    mean_anglesx.resize(no_cuts_x);
    vector<vector<double>> mean_anglesy;
    mean_anglesy.resize(no_cuts_y);
    cout << "no_cuts_y:\t" << no_cuts_y << "\tno_cuts_x:\t" << no_cuts_x << "\n";

    // int Nevents = hitcoords[0].size();

    double M1M2_slope_lb_x, M1M2_slope_ub_x, M1M2_slope_lb_y, M1M2_slope_ub_y, dv_x = 0, dv_y = 0;

    for (int dir = 0; dir < 2; dir++)
    {
        if (dir == 0)
        {
            filename = "angles_mean_x_" + (string)runno;
            no_cuts = no_cuts_x;
        }
        else
        {
            filename = "angles_mean_y_" + (string)runno;
            no_cuts = no_cuts_y;
        }

        int cutno = 0;

        // #pragma omp parallel for
        for (int cut = 0; cut < no_cuts; cut++)
        {
            if (dir == 0)
            {
                M1M2_slope_lb_x = xcut_lb;
                M1M2_slope_ub_x = xcut_ub;
                M1M2_slope_lb_y = -1e+17;
                M1M2_slope_ub_y = 1e+17;
                dv_x = 0;
            }
            else
            {
                M1M2_slope_lb_x = -1e+17;
                M1M2_slope_ub_x = 1e+17;
                M1M2_slope_lb_y = ycut_lb;
                M1M2_slope_ub_y = ycut_ub;
                dv_y = 0;
            }

            int count = 0;

            for (int i = 0; i < Nevents; i++)
            {
                for (size_t j = 0; j < hitcoords[0][i].size(); j++)
                {
                    for (size_t k = 0; k < hitcoords[1][i].size(); k++)
                    {
                        vector<double> M1M2_Proj = rect_project(hitcoords[0][i][j], hitcoords[1][i][k], M1M2_z);

                        vec M1M2_x = {hitcoords[0][i][j][0], hitcoords[1][i][k][0]};
                        vec M1M2_y = {hitcoords[0][i][j][1], hitcoords[1][i][k][1]};
                        double ang_x_M1M2 = calc_slope(vec_M1M2_z, M1M2_x);
                        double ang_y_M1M2 = calc_slope(vec_M1M2_z, M1M2_y);

                        /* Only pick out certain angles of entry into crystal */
                        if (ang_x_M1M2 > M1M2_slope_lb_x and ang_x_M1M2 < M1M2_slope_ub_x and ang_y_M1M2 > M1M2_slope_lb_y and ang_y_M1M2 < M1M2_slope_ub_y)
                        {
                            for (size_t l = 0; l < hitcoords[2][i].size(); l++)
                            {
                                double M1M2_d = calc_dist(M1M2_Proj[0], M1M2_Proj[1], hitcoords[2][i][l][0], hitcoords[2][i][l][1]);

                                if (M1M2_d < M1M2_d_lim)
                                {
                                    vec M2M3_x = {hitcoords[1][i][k][0], hitcoords[2][i][l][0]};
                                    vec M2M3_y = {hitcoords[1][i][k][1], hitcoords[2][i][l][1]};
                                    double ang_x_M2M3 = calc_slope(vec_M2M3_z, M2M3_x);
                                    double ang_y_M2M3 = calc_slope(vec_M2M3_z, M2M3_y);

                                    /* Calculate sum of delta angles */
                                    if (dir == 0)
                                        dv_x += abs(ang_x_M2M3 - ang_x_M1M2);
                                    else
                                        dv_y += abs(ang_y_M2M3 - ang_y_M1M2);

                                    count++;
                                }
                            } // hits in M3 done
                        }     // end cut
                    }         // hits in M2 done
                }             // hits in M1 done
            }                 // events done

            /* Calculate mean of delta angles */
            if (dir == 0)
            {
                dv_x /= count;
                mean_anglesx[cut] = {xcut_lb + dtheta / 2.0, dv_x};
                xcut_lb += dtheta;
                xcut_ub += dtheta;
            }
            else
            {
                dv_y /= count;
                mean_anglesy[cut] = {ycut_lb + theta / 2.0, dv_y};
                ycut_lb += dtheta;
                ycut_ub += dtheta;
            }

            cutno++;

            if ((cutno + 1) % (no_cuts / 10 + 1) == 0)
                cout << "Progress: \t" << 100 * cutno / no_cuts << "%\n";
        } // cuts done

        if (dir == 0)
            save_vector(filename, mean_anglesx);
        else
            save_vector(filename, mean_anglesy);
    }
}

/* Calculate distances between projected and observes hits */
vector<double> analyser::calc_interdistance(vector<vector<vector<double>>> Hits0, vector<vector<vector<double>>> Hits1, vector<vector<vector<double>>> Hits2, vector<double> z)
{
    vector<double> distvec;

    for (int i = 0; i < Nevents; i++)
    {
        vector<vector<double>> hits0 = Hits0[i];
        vector<vector<double>> hits1 = Hits1[i];
        vector<vector<double>> hits2 = Hits2[i];

        for (size_t j = 0; j < hits0.size(); j++)
        {
            for (size_t k = 0; k < hits1.size(); k++)
            {
                vector<double> hitp = rect_project(hits0[j], hits1[k], z);

                for (size_t l = 0; l < hits2.size(); l++)
                {
                    double newdist = calc_dist(hitp[0], hitp[1], hits2[l][0], hits2[l][1]);
                    distvec.push_back(newdist);
                }
            }
        }
    }
    return distvec;
}

double analyser::calc_dist(double x0, double y0, double x1, double y1)
{
    return sqrt((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0)); // pythagoras
}

vector<double> analyser::rect_project(vector<double> hit0, vector<double> hit1, vector<double> z)
{
    vector<double> proj(2);

    proj[0] = (z[2] - z[0]) * (hit0[0] - hit1[0]) / (z[0] - z[1]) + hit0[0];
    proj[1] = (z[2] - z[0]) * (hit0[1] - hit1[1]) / (z[0] - z[1]) + hit0[1];

    return proj;
}

/* Align detector relative to M1-M2 */
void analyser::align_plane(mat &mat_Tot, vector<vector<vector<double>>> Hits0, vector<vector<vector<double>>> Hits1, vector<vector<vector<double>>> &Hits2, vector<double> z)
{
    vector<double> dr;

    if (z.back() == zplanes[2] or z.back() == zplanes[3])
    {
        dr = {dr_crit_list[0], dr_crit_list[1], dr_crit_list[2]};
    }
    else
        dr = dr_crit_list;

    for (size_t j = 0; j < dr.size(); j++)
    {
        double dr_crit = dr[j];
        bool convergence;
        mat mat_temp_T = zeros<mat>(3, 3);

        do
        {
            mat mat_T = construct_T_mat(Hits0, Hits1, Hits2, dr_crit, z);
            mat_Tot = mat_T * mat_Tot;
            adjust_coordinates(Hits2, mat_T);
            convergence = check_convergence(mat_T, mat_temp_T);
            mat_temp_T = mat_T;

        } while (!convergence);
    }
}

/* Calculate matrix for detector alignment */
mat analyser::construct_T_mat(vector<vector<vector<double>>> Hits0, vector<vector<vector<double>>> Hits1, vector<vector<vector<double>>> Hits2, double dr_crit, vector<double> z)
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
                vector<double> hitp = rect_project(hits0[l], hits1[m], z);

                for (size_t n = 0; n < hits2.size(); n++)
                {
                    double dr = calc_dist(hitp[0], hitp[1], hits2[n][0], hits2[n][1]);

                    if (dr < dr_crit)
                    {
                        save_hit(hitp, hits2[n], mat_expected, mat_observed, rowno);
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
void analyser::adjust_coordinates(vector<vector<vector<double>>> &Hits, mat mat_T)
{
    for (size_t i = 0; i < Hits.size(); i++)
    {
        for (size_t j = 0; j < Hits[i].size(); j++)
        {
            vec hit = {{Hits[i][j][0]}, {Hits[i][j][1]}, {1}};
            vec temp = mat_T * hit;
            Hits[i][j][0] = temp[0];
            Hits[i][j][1] = temp[1];
        }
    }
}

/* Allocate hitcoordinates as a row in a matrix */
void analyser::save_hit(vector<double> hitp, vector<double> hit, mat &mat_expected, mat &mat_observed, int indx)
{
    rowvec vec_expected = {hitp[0], hitp[1], 1};
    rowvec vec_observed = {hit[0], hit[1], 1};
    mat_expected.row(indx) = vec_expected;
    mat_observed.row(indx) = vec_observed;
}

/* Check similarity between two matrices */
bool analyser::check_convergence(mat A, mat B)
{
    return sum(abs(vectorise(A) - vectorise(B))) < tol;
}

void analyser::construct_distarray(void)
{
    for (size_t i = 0; i < 4; i++)
    {
        vector<double> zproj = {zplanes[i], zplanes[i + 1], zplanes[i + 2]};
        distarray[i] = calc_interdistance(hitcoords[i], hitcoords[i + 1], hitcoords[i + 2], zproj);
    }
}

/* Align planes with T matrix */
void analyser::align_w_T(void)
{
    for (size_t i = 0; i < 4; i++)
    {
        vector<double> zproj = {zplanes[i], zplanes[i + 1], zplanes[i + 2]};
        mat mat_T = T.slice(i);

        for (size_t j = 0; j < hitcoords[i + 2].size(); j++)
        {
            for (size_t k = 0; k < hitcoords[i + 2][j].size(); k++)
            {
                vec hit = {{hitcoords[i + 2][j][k][0]}, {hitcoords[i + 2][j][k][1]}, {1}};
                vec temp = mat_T * hit;
                hitcoords[i + 2][j][k][0] = temp[0];
                hitcoords[i + 2][j][k][1] = temp[1];
            }
        }

        distarray[i] = calc_interdistance(hitcoords[i], hitcoords[i + 1], hitcoords[i + 2], zproj);
        cout << "Aligned plane " << i + 2 << "\n";
    }
}

/* Align planes without T matrix */
void analyser::align_wo_T(void)
{
    for (int i = 0; i < 4; i++)
    {
        mat mat_T = eye<mat>(3, 3);
        vector<double> zproj = {zplanes[i], zplanes[i + 1], zplanes[i + 2]};
        align_plane(mat_T, hitcoords[i], hitcoords[i + 1], hitcoords[i + 2], zproj);
        distarray[i] = calc_interdistance(hitcoords[i], hitcoords[i + 1], hitcoords[i + 2], zproj);
        cout << "Aligned plane " << i + 2 << "\n";
        T.slice(i) = mat_T;
    }

    T.save(DATPATH + "/Align/alignment_matrix.txt", arma_ascii);
}

/* Construct pixel grid */
void analyser::make_grid(vector<vector<double>> &pixelgrid, vector<double> &xgrid, vector<double> &ygrid)
{
    pixelgrid.resize(ncols * nrows);

    for (size_t i = 0; i < pixelgrid.size(); i++)
    {
        pixelgrid[i] = {(double)i, (double)0};
    }
}

/* Save coordinates of hits based on which detector recorded it */
void analyser::extract_hit_data(vector<vector<vector<double>>> &hitcoord, vector<vector<double>> &pixelgrid, int plane)
{
    hitcoord.resize(Nevents);

    for (int i = 0; i < Nevents; i++)
    {
        vector<vector<double>> eventdata = Events[i];

        for (size_t j = 0; j < eventdata.size(); j++)
        {
            vector<double> hitdata = eventdata[j];

            if (hitdata[2] == plane)
            {
                int pixel = coord2pixel(hitdata[0], hitdata[1]);
                pixelgrid[pixel][1] += 1;
                hitcoord[i].push_back({hitdata[0], hitdata[1]});
            }
        }
    }
}

/* Determine pixel no from coordinates */
int analyser::coord2pixel(double xhit, double yhit)
{
    double dx = (xmax - xmin) / ncols;
    double dy = (ymax - ymin) / nrows;
    int colno = (xhit + xmax) / dx;
    int rowno = (yhit + ymax) / dy;

    return rowno * ncols + colno;
}

/* Count total number of hits across all events */
void analyser::count_hits(int &count, vector<vector<vector<double>>> hitcoord)
{
    count = 0;

    for (size_t i = 0; i < hitcoord.size(); i++)
    {
        for (size_t j = 0; j < hitcoord[i].size(); j++)
        {
            count++;
        }
    }
}

/*  Criteria to sort vector by 2nd column, descending */
bool analyser::sortFunc(const vector<double> &p1, const vector<double> &p2)
{
    return p1[1] > p2[1];
}

void analyser::locate_hot_pixels(vector<vector<double>> pixelgrid, vector<int> &hotpixels, int i)
{
    /*
           We need to identify hot pixels. We do this by sorting the pixelgrid by number of pixels, descendind,
           and then saving the pixelno. of the hotpixels.
         */

    sort(pixelgrid.begin(), pixelgrid.end(), sortFunc);
    double lim = 4E-04 * Nevents;

    for (size_t i = 0; i < pixelgrid.size(); i++)
    {
        if (pixelgrid[i][1] > lim)
        {
            hotpixels.push_back(pixelgrid[i][0]);
        }
        else
            break;
    }
}

/* We remove both the number of counts in pixelgrid for hotpixels, and erase hitcoordinates from hotpixels. */
void analyser::remove_hot_pixels(vector<vector<vector<double>>> &hitcoord, vector<vector<double>> &pixelgrid, vector<int> hotpixels)
{
    for (size_t i = 0; i < hotpixels.size(); i++)
    {
        int hotpix = hotpixels[i];
        pixelgrid[hotpix][1] = 0;

        for (size_t j = 0; j < hitcoord.size(); j++)
        {
            for (size_t k = 0; k < hitcoord[j].size(); k++)
            {
                double xhit = hitcoord[j][k][0];
                double yhit = hitcoord[j][k][1];
                int pixel = coord2pixel(xhit, yhit);

                if (pixel == hotpix)
                {
                    hitcoord[j].erase(hitcoord[j].begin() + k);
                    k--;
                }
            }
        }
    }
}

vector<double> analyser::linspace(double min, double max, int N)
{
    vector<double> range;
    range.resize(N);

    /* Determine step-size */
    double delta = (max - min) / (N - 1);

    /* Iterate over range, excluding the last entry */
    for (int i = 0; i < N - 1; i++)
    {
        range[i] = min + i * delta;
    }

    /* Set last entry to "max". This ensures that we get a range from min -> max */
    range.back() = max;
    return range;
}

/* Extract raw data from root file, do pre-processing and save usable data in std::vectors */
void analyser::extract_root_data(void)
{
    data = true;
    /* Open root file and obtain relevant data */
    TFile *f = TFile::Open(filename);
    TTree *T1 = (TTree *)f->Get("T");
    TLeaf *Hpk = (TLeaf *)T1->GetLeaf("fAHits.Hpk");
    TLeaf *Hu = (TLeaf *)T1->GetLeaf("fAHits.Hu");
    TLeaf *Hv = (TLeaf *)T1->GetLeaf("fAHits.Hv");
    TLeaf *fAHitsN = (TLeaf *)T1->GetLeaf("fAHitsN");
    TLeaf *ENumberOfTriggers = (TLeaf *)T1->GetLeaf("fHeader.ENumberOfTriggers");

    int N_events = T1->GetEntries();
    cout << "\nTotal number of events :\t" << N_events << "\n";

    Events.resize(N_events);
    Nevents = 0;

    /* Sort data into std::vector */
    for (int i = 0; i < N_events; i++)
    {
        T1->GetEntry(i);
        int N_hits = fAHitsN->GetValue();
        int N_triggers = ENumberOfTriggers->GetValue();

        if (N_triggers == 2 && N_hits > 0)
        {
            vector<vector<double>> EventData(N_hits);

            for (int j = 0; j < N_hits; j++)
            {
                vector<double> HitData(3);
                HitData[0] = Hu->GetValue(j);                           // xcoord
                HitData[1] = Hv->GetValue(j);                           // ycoord
                HitData[2] = static_cast<double>(Hpk->GetValue(j) - 1); // plane
                EventData[j] = HitData;
            }

            Events[Nevents] = EventData;
            Nevents++;
        }
    }

    Events.resize(Nevents);
    cout << "Number of usable events :\t" << Nevents << "\n\n";

    /* Save no. of events in run on disk */
    string eventfile_name = "events_run.txt";
    fstream event_number_stream(eventfile_name, fstream::in | fstream::out | fstream::app);

    if (event_number_stream.is_open())
    {
        event_number_stream << runno << "\t" << Nevents << "\n";
        event_number_stream.close();
    }
    else
        cerr << "Cannot open file:\t" << eventfile_name << "\n\n";
}

int analyser::isInside(int nvert, vector<double> vertx, vector<double> verty, double testx, double testy)
{
    int i, j, c = 0;

    for (i = 0, j = nvert - 1; i < nvert; j = i++)
    {
        if (((verty[i] > testy) != (verty[j] > testy)) && (testx < (vertx[j] - vertx[i]) * (testy - verty[i]) / (verty[j] - verty[i]) + vertx[i]))
        {
            c = !c;
        }
    }

    return c;
}
