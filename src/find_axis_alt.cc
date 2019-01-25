void analyser::find_axis_alt (void) {

        /* Load beam-parameters from file */
        vector<double> params;
        double param;
        string file = (string)paramsfile;
        ifstream paramters (file);

        if (paramters.is_open()) {

                while (true) {

                        paramters >> param;
                        cerr << param << "\n";
                        if (paramters.eof()) break;

                        params.push_back(param);

                }

                paramters.close();

        } else cout << "Unable to open file containing beam parameters";


        vec vec_M1M2_z = {M1M2_z[0], M1M2_z[1]};
        vec vec_M2M3_z = {M2M3_z[0], M2M3_z[1]};

        double xmean = params[2];
        double ymean = params[3];
        double deltax = params[4]; // FWHM x
        double deltay = params[5]; // FWHM y
        double theta = 30E-6;
        double dtheta = 2.0E-6;

        double ycut_lb = ymean - deltay - theta;
        double ycut_ub = ymean - deltay + theta;
        double xcut_lb = xmean - deltax - theta;
        double xcut_ub = xmean - deltax + theta;

        int no_cuts = 0;
        string filename;

        vector<vector<vector<vector<double> > > > tracks_alt;
        tracks_alt.resize(Nevents);

        cerr << 738 << "\n";
        double M1M2_slope_lb_x, M1M2_slope_ub_x, M1M2_slope_lb_y, M1M2_slope_ub_y;
        vector<vector<double> > photons_in_cut_x; photons_in_cut_x.resize(2 * deltax / dtheta);
        cerr << "y:\t" << 2 * deltay / dtheta << "\t x:\t" << 2 * deltax / dtheta << "\n";
        vector<vector<double> > photons_in_cut_y; photons_in_cut_y.resize(2 * deltay / dtheta);
        cerr << 743 << "\n";
        for (int dir = 0; dir < 2; dir++) {

                if (dir == 0) {

                        filename = "no_photons_x_alt" + (string)runno;
                        no_cuts = 2 * deltax / dtheta;

                } else {

                        filename = "no_photons_y_alt" + (string)runno;
                        no_cuts = 2 * deltay / dtheta;

                }

                int cutno = 0;

                // #pragma omp parallel for
                for (int cut = 0; cut < no_cuts; cut++) {

                        if (dir == 0) {

                                M1M2_slope_lb_x = xcut_lb;
                                M1M2_slope_ub_x = xcut_ub;
                                M1M2_slope_lb_y = -1e+17;
                                M1M2_slope_ub_y = 1e+17;

                        } else {

                                M1M2_slope_lb_x = -1e+17;
                                M1M2_slope_ub_x = 1e+17;
                                M1M2_slope_lb_y = ycut_lb;
                                M1M2_slope_ub_y = ycut_ub;

                        }

                        for (int i = 0; i < Nevents; i++) {

                                int track_count = 0;

                                tracks_alt[i].resize(10 * hitcoords[5][i].size());

                                vector<vector<vector<double> > > M1_MM_tracks;
                                vector<vector<double> > M1_M6_track;

                                for (size_t j = 0; j < hitcoords[0][i].size(); j++) {

                                        for (size_t k = 0; k < hitcoords[1][i].size(); k++) {

                                                vector<double> M1M2_Proj = rect_project(hitcoords[0][i][j], hitcoords[1][i][k], M1M2_z);

                                                vec M1M2_x = {hitcoords[0][i][j][0], hitcoords[1][i][k][0]};
                                                vec M1M2_y = {hitcoords[0][i][j][1], hitcoords[1][i][k][1]};
                                                double ang_x_M1M2 = calc_slope(vec_M1M2_z, M1M2_x);
                                                double ang_y_M1M2 = calc_slope(vec_M1M2_z, M1M2_y);

                                                /* Only pick out certain angles of entry into crystal */
                                                if (ang_x_M1M2 > M1M2_slope_lb_x and ang_x_M1M2 < M1M2_slope_ub_x and ang_y_M1M2 > M1M2_slope_lb_y and ang_y_M1M2 < M1M2_slope_ub_y) {

                                                        for (size_t l = 0; l < hitcoords[2][i].size(); l++) {

                                                                double M1M2_d = calc_dist(M1M2_Proj[0], M1M2_Proj[1], hitcoords[2][i][l][0], hitcoords[2][i][l][1]);

                                                                if (M1M2_d < M1M2_d_lim) {

                                                                        vector<double> M2M3_Proj = rect_project(hitcoords[1][i][k], hitcoords[2][i][l], M2M3_z);

                                                                        for (size_t m = 0; m < hitcoords[3][i].size(); m++) {
                                                                                double M2M3_d = calc_dist(M2M3_Proj[0], M2M3_Proj[1], hitcoords[3][i][m][0], hitcoords[3][i][m][1]);

                                                                                if (M2M3_d < M2M3_d_lim) {

                                                                                        vector<double> M3M4_Proj = rect_project(hitcoords[2][i][l], hitcoords[3][i][m], M3M4_z);
                                                                                        M1_MM_tracks.push_back({hitcoords[0][i][j], hitcoords[1][i][k], hitcoords[2][i][l], hitcoords[3][i][m], M3M4_Proj});
                                                                                }

                                                                        }

                                                                }

                                                        } // hits in M3 done

                                                } // end cut

                                        } // hits in M2 done

                                } // hits in M1 done

                                /* Iterate over M6 -> M6 tracks */
                                for (size_t n = 0; n < hitcoords[5][i].size(); n++) {

                                        vector<vector<double> > M1_M6_track;
                                        double shortdist = INFINITY;
                                        bool trackfound = false;

                                        for (size_t o = 0; o < hitcoords[4][i].size(); o++) {

                                                vector<double> M6M5_Proj = rect_project(hitcoords[5][i][n], hitcoords[4][i][o], M6M5_z);

                                                /* For every M5-M6 combination, pick an M3M4 projection */
                                                for (size_t l = 0; l < M1_MM_tracks.size(); l++) {

                                                        vector<double> M3M4_Proj = M1_MM_tracks[l].back();
                                                        double M6M5_d = calc_dist(M6M5_Proj[0], M6M5_Proj[1], M3M4_Proj[0], M3M4_Proj[1]);

                                                        vec vec_M3M4_z = {M3M4_z[0], M3M4_z[1]};
                                                        vec M3M4_y = {M1_MM_tracks[l][2][1], M1_MM_tracks[l][3][1]};
                                                        double M3M4_slope = calc_slope(vec_M3M4_z, M3M4_y);

                                                        if (M6M5_d < shortdist and M6M5_d < M6M5_d_lim) {

                                                                shortdist = M6M5_d;
                                                                vec vec_M5M6_z = {M6M5_z[0], M6M5_z[1]};
                                                                vec M5M6_y = {hitcoords[5][i][n][1], hitcoords[4][i][o][1]};
                                                                double M5M6_slope = calc_slope(vec_M5M6_z, M5M6_y);

                                                                if (abs(M5M6_slope - M3M4_slope) < yz_defl_lim) {

                                                                        trackfound = true;
                                                                        M1_M6_track = {M1_MM_tracks[l][0], M1_MM_tracks[l][1], M1_MM_tracks[l][2], M1_MM_tracks[l][3], M3M4_Proj, M6M5_Proj, hitcoords[4][i][o], hitcoords[5][i][n]};

                                                                }

                                                        }

                                                }

                                        } // hits in M6 done

                                        if (trackfound) {

                                                tracks_alt[i][track_count] = M1_M6_track;
                                                track_count++;

                                        } // if statement for push_back of tracks

                                } // hits in M5 done

                                tracks_alt[i].resize(track_count);

                        } // events done

                        /* Pair tracks */
                        double no_photons = 0.0;
                        double energy_tot = 0;

                        for (size_t i = 0; i < tracks_alt.size(); i++) { // no. of events

                                for (size_t j = 0; j < tracks_alt[i].size(); j++) { // tracks in event

                                        int matchedtracks = 0;

                                        for (size_t k = j+1; k < tracks_alt[i].size(); k++) {

                                                double dist = calc_dist(tracks_alt[i][j][4][0], tracks_alt[i][j][4][1], tracks_alt[i][k][4][0], tracks_alt[i][k][4][1]);

                                                if (dist < Match_d) {

                                                        double M3_dist = calc_dist(tracks_alt[i][j][2][0], tracks_alt[i][j][2][1], tracks_alt[i][k][2][0], tracks_alt[i][k][2][1]);

                                                        if (M3_dist < Match_d_foil) {

                                                                double energy = calc_pair_energy({tracks_alt[i][j], tracks_alt[i][k]});
                                                                energy_tot += energy;
                                                                no_photons += 1;

                                                                /* Remove additional matched track since we do not want it to potentially match any other track */
                                                                if (matchedtracks > 1) {

                                                                        tracks_alt[i].erase(tracks_alt[i].begin() + k);
                                                                        k--;

                                                                }

                                                        }

                                                }

                                        }

                                }

                        }

                        /* Calculate number of photons within cut */
                        if (dir == 0) {

                                cerr << energy_tot*6.2415091E+9 << "\t" << "\t" << (double)no_photons << "\t" <<  6.2415091E+9*energy_tot/(double)no_photons << "\n";

                                photons_in_cut_x[cut] = {xcut_lb + dtheta/2.0, 6.2415091E+9 * energy_tot/(double)no_photons};
                                xcut_lb += dtheta;
                                xcut_ub += dtheta;

                        } else {

                                photons_in_cut_y[cut] = {ycut_lb + theta/2.0, 6.2415091E+9 *  energy_tot/(double)no_photons};
                                ycut_lb += dtheta;
                                ycut_ub += dtheta;

                        }

                        cutno++;

                        if ( (cutno + 1) % (no_cuts/10 + 1) == 0) cerr << "Progress: \t" << 100 * cutno/no_cuts << "%\n";

                } // cuts done

                if (dir == 0) save_vector(filename, photons_in_cut_x);
                else save_vector(filename, photons_in_cut_y);

        }

}
