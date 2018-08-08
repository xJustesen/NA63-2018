#include "analyser.h"

analyser::analyser(vector<double> z, char const *name, char const *run) {
  /* Paramters for MIMOSA detectors */
  nrows = 576;
  ncols = 1152;
  xmin = -11000, xmax = 11000, ymin = -5500, ymax = 5500;
  zplanes = z;  // z-coords of planes
  M1M2_z = {zplanes[0], zplanes[1], zplanes[2]};  // z-coords for M1-M3 PLANES
  M2M3_z = {zplanes[1], zplanes[2], zplanes[3]};  // z-coords for M2-M4 planes
  M3M4_z = {zplanes[2], zplanes[3], (zplanes[4]+zplanes[3])/2.0}; // z-coords of M3-MM planes
  M6M5_z = {zplanes[4], zplanes[5], (zplanes[4]+zplanes[3])/2.0}; // z-coords for M5-MM planes

  /* Conditions for track construction + pairing */
  M1M2_d_lim = 2500.0;  // acceptable distance between M1->M2 projected point and observed hit in M3 [micro-meter]
  M2M3_d_lim = 1000.0; // acceptable distance between M2->M3 projected point and observed hit in M4 [micro-meter]
  M6M5_d_lim = 500.0; // acceptable distance in MM bewteen M1->MM and M6->MM arm of track [micro-meter]
  Match_d = 250.0; // matching ditance in MM between electron/positron paired track [micro-meter]
  Match_d_foil = 250.0; // matching ditance in foil between electron/positron paired track [micro-meter]
  yz_defl_lim = 5.0e-3; // angle of acceptable y-deflction [rad]

  /* Alignment conditions */
  dr_crit_list = {3500, 2500, 1500, 500, 250, 150};  // radius for acceptable hits when aligning [micro-meter]
  tol = 1e-6; // convergence tolerance

  /* Info for data storage on disk */
  runno = run;
  filename = name;
  DATPATH = "/home/christian/Dropbox/speciale/data";

  /* Vectors for data storage in RAM */
  hitcoords.resize(6);
  pixelgrids.resize(6);
  hotpixels.resize(6);
  T.set_size(3, 3, 4);
  distarray.resize(4);
}

void analyser::update_pixelgrids(int plane, vector<vector<double>> pixelgrid) {
  pixelgrids[plane] = pixelgrid;
}

void analyser::update_hitcoords(int plane, vector<vector<vector<double>>> hitcoord) {
  hitcoords[plane] = hitcoord;
}

void analyser::update_hotpixels(int plane, vector<int> hotpixel) {
  hotpixels[plane] = hotpixel;
}

void analyser::print_pixels(void) {
  string filename = DATPATH + "/pixeldata_run" + runno + ".txt";
  ofstream output (filename);
  for (size_t i = 0; i < pixelgrids.size(); i++) {
    output << "## PLANE\t" << i << "\tPIXEL DATA\n";
    for (size_t j = 0; j < pixelgrids[i].size(); j++) {
      output << pixelgrids[i][j][1] << "\n";
    }
  }
}

void analyser::print_hotpixels(void) {
  for (size_t i = 0; i < hotpixels.size(); i++) {
    string filename = DATPATH + "/hotpixels_run" + runno +  "_plane_" + to_string(i) + ".txt";
    ofstream output (filename);
    for (size_t j = 0; j < hotpixels[i].size(); j++) {
      output << hotpixels[i][j] << "\n";
    }
  }
}

void analyser::print_interdistance(void){
  string filename = DATPATH + "/interdistance_data_" + runno + ".txt";;
  ofstream output (filename);
  for (size_t i = 0; i < distarray.size(); i++) {
    output << "## DATABLOCK\t" << i << "\tINTERDISTANCE DATA\n";
    for (size_t j = 0; j < distarray[i].size(); j++) {
      output << distarray[i][j] << "\n";
    }
    output << "\n\n";
  }
}

void analyser::print_hits(void){
  string filename = DATPATH + "/hits_coord_data_" + runno + ".txt";
  ofstream output (filename);
  for (size_t i = 0; i < hitcoords.size(); i++) {
    output << "## PLANE\t" << i << "\tHIT DATA\n";
    for (size_t j = 0; j < hitcoords[i].size(); j++) {
      for (size_t k = 0; k < hitcoords[i][j].size(); k++) {
      output << hitcoords[i][j][k][0] << ' ' << hitcoords[i][j][k][1] << '\n';
      }
    }
  }
}

void analyser::print_energy(void) {
  string filename = DATPATH + "/energy_" + runno + ".txt";
  ofstream output (filename);
  for (size_t i = 0; i < energies.size(); i++) {
    output << energies[i] << "\n";
  }
}

void analyser::print_slope(void) {
  string filename = DATPATH + "/angles_" + runno + ".txt";
  ofstream output (filename);
  for (size_t i = 0; i < slopes.size(); i++) {
    output << slopes[i][0] << " " << slopes[i][1] << "\n";
  }
}

void analyser::print_M1M2_slope(void) {
  string filename1 = DATPATH + "/angles_M1M2_" + runno + ".txt";
  string filename2 = DATPATH + "/beam_divergence_"  + runno + ".txt";
  ofstream output1 (filename1);
  ofstream output2 (filename2);
  for (int i = 0; i < Nevents; i++) {
    for (size_t j = 0; j < M1M2_slopes[i][0].size(); j++) {
      output1 << M1M2_slopes[i][0][j] << "\t" << M1M2_slopes[i][1][j] << "\n";
    }
    for (size_t j = 0; j < divergence[i][0].size(); j++) {
      output2 << divergence[i][0][j] << "\t" << divergence[i][1][j] << "\n";
    }
  }
}

void analyser::print_zpos(void) {
  string filename = DATPATH + "/zclosepos_" + runno + ".txt";
  ofstream output (filename);
  for (size_t i = 0; i < zclosepos.size(); i++) {
    output << zclosepos[i] << "\n";
  }
}

double analyser::calc_ang(double m, double n) {
  return n - m;
}

double analyser::calc_slope(vec x, vec y) {
    return (y(1) - y(0))/(x(1) - x(0));
}

mat analyser::lines3d_nearestpoints(vec A, vec B, vec C, vec D) {
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

void analyser::image_crystal(void) {
  string filename = DATPATH + "/crystal_image_" + runno + ".txt";
  ofstream output (filename);
  for (size_t i = 0; i < paired_tracks.size(); i++) {
    output << paired_tracks[i][0][1][0] << "\t" << paired_tracks[i][0][1][1] << "\n" << paired_tracks[i][1][1][0] << "\t" << paired_tracks[i][1][1][1] << "\n";
  }
}

double analyser::calc_pair_energy(vector<vector<vector<double>>> pairedtracks) {
  double L = 0.15, q = 1.6021766e-19, B = 0.12, c = 299792458; // SI units
  double energy = 0;
  /*
  The energy of the photon corresponds to the total energy of each particle in the pair.
  The energy of each particle is calculated seperately in this loop
  */
  for (int i = 0; i < 2; i++) {
    vec M4MM_x = {pairedtracks[i][3][0], pairedtracks[i][4][0]};
    vec M4MM_z = {zplanes[3], (zplanes[4]+zplanes[3])/2};
    vec M5M6_x = {pairedtracks[i][6][0], pairedtracks[i][7][0]};
    vec M5M6_z = {zplanes[4], zplanes[5]};

    double m = calc_slope(M4MM_z, M4MM_x);
    double n = calc_slope(M5M6_z, M5M6_x);
    vector<double> slope = {m, n};
    slopes.push_back(slope);

    double ang = calc_ang(m, n);  // calculate deflection angle to find energy
    energy += q*c*L*B/abs(ang); // Joule
  }

  return energy;
}

void analyser::pair_tracks(void) {
  /* Pair tracks determined by "construct_tracks" method */
  double dist;
  int pairs = 0, matchedtrack = 0, tot_matched_tracks = 0, erased_tracks = 0;

  paired_tracks.resize(tracks.size()*200);
  energies.resize(tracks.size()*200);
  zclosepos.resize(tracks.size()*200);

  /*
  Here we find pairs. We only want a track to pair exactly one other track.
  To ensure this we iterate first from i = 0 -> number of tracks, then j = i + 1 -> number of tracks.
  If multiple pairs are found we erase them from the tracks array.
  */
  for (size_t i = 0; i < tracks.size(); i++) {
    int matchedtracks = 0;
    for (size_t j = i+1; j < tracks.size(); j++) {
      if (tracks[i][8][0] == tracks[j][8][0]) {
        dist = calc_dist(tracks[i][4][0], tracks[i][4][1], tracks[j][4][0], tracks[j][4][1]);
        if (dist < Match_d) {
          double M3_dist = calc_dist(tracks[i][2][0], tracks[i][2][1], tracks[j][2][0], tracks[j][2][1]);
          if (M3_dist < Match_d_foil) {
            matchedtracks++;
            tot_matched_tracks++;
            if (matchedtracks > 1) {
              tracks.erase(tracks.begin() + j);
              erased_tracks++;
              j--;
            }
            else if (matchedtracks == 1) {
              matchedtrack = j;
            }
          }
        }
      }
    }
    /* Only save pair if there is a single match */
    if (matchedtracks == 1) {
      vec P1 = {tracks[i][6][0], tracks[i][6][1], M6M5_z[0]}, Q1 = {tracks[i][7][0], tracks[i][7][1], M6M5_z[1]};
      vec P2 = {tracks[matchedtrack][2][0], tracks[matchedtrack][2][1], M3M4_z[0]}, Q2 = {tracks[matchedtrack][3][0], tracks[matchedtrack][3][1], M6M5_z[1]};
      mat closepos = lines3d_nearestpoints(P1, Q1, P2, Q2);
      double zclosepospair = (closepos(2, 0) + closepos(2, 1))/2;
      if (zclosepospair > 8e6 && zclosepospair < 1e7){
        paired_tracks[pairs] = {tracks[i], tracks[matchedtrack]};
        energies[pairs] = calc_pair_energy(paired_tracks[pairs]);
        zclosepos[pairs] = zclosepospair;
        pairs++;
      }
    }
  }

  paired_tracks.resize(pairs);
  energies.resize(pairs);
  zclosepos.resize(pairs);

  cerr << "Number of paired tracks : " << pairs << "\n";
}

void analyser::beam_divergence(int eventno, int x_indx, int y_indx, string name) {
  string str1 = "partial";
  string str2 = "full";

  vec M1M2_x = {hitcoords[0][eventno][x_indx][0], hitcoords[1][eventno][y_indx][0]};
  vec M1M2_y = {hitcoords[0][eventno][x_indx][1], hitcoords[1][eventno][y_indx][1]};
  vec vec_M1M2_z = {M1M2_z[0], M1M2_z[1]};

  double ang_x = calc_slope(vec_M1M2_z, M1M2_x);
  double ang_y = calc_slope(vec_M1M2_z, M1M2_y);

  if (name.compare(str2) == 0) {
    divergence[eventno][0].push_back(ang_x);
    divergence[eventno][1].push_back(ang_y);
  }
  if (name.compare(str1) == 0) {
    M1M2_slopes[eventno][0].push_back(ang_x);
    M1M2_slopes[eventno][1].push_back(ang_y);
  }
}

void analyser::construct_tracks(void) {
  /* Construct particle tracks from M1 -> M6 */
  int M1_MM_tot_tracks = 0, M1_M6_tot_tracks = 0;

  double M1M2_slope_lb = -INFINITY, M1M2_slope_ub = INFINITY;

  M1M2_slopes.resize(Nevents);
  divergence.resize(Nevents);

  string str1 = "full", str2 = "partial";

  for (int i = 0; i < Nevents; i++) {
    vector<vector<double>> M1_M6_track;

    M1M2_slopes[i].resize(2);
    divergence[i].resize(2);

    /*
    Here we construct tracks. Each loop iterates over hits in a mimosa detector (0 -> 5).
    Conditions to find tracks checked periodically. All acceptable tracks are saved, which
    means bad tracks are saved. These should later be discarded in "pair_tracks" function
    */
    for (size_t j = 0; j < hitcoords[0][i].size(); j++) {
      for (size_t k = 0; k < hitcoords[1][i].size(); k++) {
        vector<double> M1M2_Proj = rect_project(hitcoords[0][i][j], hitcoords[1][i][k], M1M2_z);
        beam_divergence(i, j, k, str1);
        for (size_t l = 0; l < hitcoords[2][i].size(); l++) {
          double M1M2_d = calc_dist(M1M2_Proj[0], M1M2_Proj[1], hitcoords[2][i][l][0], hitcoords[2][i][l][1]);
          if (M1M2_d < M1M2_d_lim){
            vector<double> M2M3_Proj = rect_project(hitcoords[1][i][k], hitcoords[2][i][l], M2M3_z);
            for (size_t m = 0; m < hitcoords[3][i].size(); m++) {
              double M2M3_d = calc_dist(M2M3_Proj[0], M2M3_Proj[1], hitcoords[3][i][m][0], hitcoords[3][i][m][1]);
              vec vec_M3M4_z = {M3M4_z[0], M3M4_z[1]};
              vec M3M4_y = {hitcoords[2][i][l][1], hitcoords[3][i][m][1]};
              double M3M4_slope = calc_slope(vec_M3M4_z, M3M4_y);
              beam_divergence(i, j, k, str2);
              if (M1M2_slopes[i][0][j] > M1M2_slope_lb and M1M2_slopes[i][0][j] < M1M2_slope_ub) {
                if (M2M3_d < M2M3_d_lim) {
                  vector<double> M3M4_Proj = rect_project(hitcoords[2][i][l], hitcoords[3][i][m], M3M4_z);
                  M1_MM_tot_tracks++;
                  for (size_t n = 0; n < hitcoords[4][i].size(); n++) {
                    double shortdist = INFINITY;
                    bool trackfound = false;
                    for (size_t o = 0; o < hitcoords[5][i].size(); o++) {
                      vector<double> M6M5_Proj = rect_project(hitcoords[4][i][n], hitcoords[5][i][o], M6M5_z);
                      double M6M5_d = calc_dist(M6M5_Proj[0], M6M5_Proj[1], M3M4_Proj[0], M3M4_Proj[1]);
                      if (M6M5_d < shortdist and M6M5_d < M6M5_d_lim) {
                        vec vec_M5M6_z = {M6M5_z[0], M6M5_z[1]};
                        vec M5M6_y = {hitcoords[4][i][n][1], hitcoords[5][i][o][1]};
                        double M5M6_slope = calc_slope(vec_M5M6_z, M5M6_y);
                        if (abs(M5M6_slope - M3M4_slope) < yz_defl_lim) {
                          shortdist = M6M5_d;
                          trackfound = true;
                          M1_M6_track = {hitcoords[0][i][j], hitcoords[1][i][k], hitcoords[2][i][l], hitcoords[3][i][m], M3M4_Proj, M6M5_Proj, hitcoords[4][i][n], hitcoords[5][i][o], {(double)i}};
                        }
                      }
                    } // hits in M6 done
                    if (trackfound) {
                      tracks.push_back(M1_M6_track);
                      M1_M6_tot_tracks++;
                    }
                  } // hits in M5 done
                }
              }
            } // hits in M4 done
          }
        } // hits in M3 done
      } // hits in M2 done
    } // hits in M1 done
  } // events done

  cerr << "\nTotal tracks from M1 -> MM : " << M1_MM_tot_tracks << "\n";
  cerr << "Total tracks from M1 -> M6 : " << M1_M6_tot_tracks << "\n";
}

void analyser::calc_interdistance(vector<vector<vector<double>>> Hits0, vector<vector<vector<double>>> Hits1, vector<vector<vector<double>>> Hits2, vector<double> z, vector<double> &distvec) {
  /* Calculate distances between projected and observes hits */
  distvec.resize(500*Nevents);
  int indx = 0;
  for (int i = 0; i < Nevents; i++) {
    vector<vector<double>> hits0 = Hits0[i];
    vector<vector<double>> hits1 = Hits1[i];
    vector<vector<double>> hits2 = Hits2[i];
    for (size_t j = 0; j < hits0.size(); j++) {
      for (size_t k = 0; k < hits1.size(); k++) {
        vector<double> hitp = rect_project(hits0[j], hits1[k], z);
        for (size_t l = 0; l < hits2.size(); l++) {
          double newdist = calc_dist(hitp[0], hitp[1], hits2[l][0], hits2[l][1]);
          distvec[indx] = newdist;
          indx++;
        }
      }
    }
  }

  distvec.resize(indx);
}

double analyser::calc_dist(double x0, double y0, double x1,double y1) {
  return sqrt(pow(x1 - x0, 2) + pow(y1 - y0, 2)); // pythagoras
}

vector<double> analyser::rect_project(vector<double> hit0, vector<double> hit1, vector<double> z) {
  vector<double> proj(2);
  proj[0] = ((hit0[0]-hit1[0])/(z[0]-z[1]))*(z[2]-z[0]) + hit0[0];
  proj[1] = ((hit0[1]-hit1[1])/(z[0]-z[1]))*(z[2]-z[0]) + hit0[1];

  return proj;
}

void analyser::align_plane(mat &mat_Tot, vector<vector<vector<double>>>  Hits0, vector<vector<vector<double>>>  Hits1, vector<vector<vector<double>>>  &Hits2, vector<double> z) {
  /* Align detector relative to M1-M2 */
  for (size_t j = 0; j < dr_crit_list.size(); j++) {
    double dr_crit = dr_crit_list[j];
    mat mat_temp_T = zeros<mat>(3,3);

    do {
      mat mat_T = construct_T_mat(Hits0, Hits1, Hits2, dr_crit, z);
      mat_Tot = mat_T * mat_Tot;
      adjust_coordinates(Hits2, mat_T);
      bool convergence = check_convergence(mat_T, mat_temp_T);
      mat_temp_T = mat_T;
    } while (!convergence);
  }
}

mat analyser::construct_T_mat(vector<vector<vector<double>>> Hits0, vector<vector<vector<double>>> Hits1, vector<vector<vector<double>>> Hits2, double dr_crit,  vector<double> z) {
  /* Calculate matrix for detector alignment */
  int rowno = 0;
  mat mat_expected = zeros<mat>(Hits1.size()*200, 3), mat_observed = zeros<mat>(Hits1.size()*200,3);
  for (size_t k = 0; k < Hits0.size(); k++) {
    vector<vector<double>> hits0 = Hits0[k];
    vector<vector<double>> hits1 = Hits1[k];
    vector<vector<double>> hits2 = Hits2[k];
    for (size_t l = 0; l < hits0.size(); l++) {
      for (size_t m = 0; m < hits1.size(); m++) {
        vector<double> hitp = rect_project(hits0[l], hits1[m], z);
        for (size_t n = 0; n < hits2.size(); n++) {
          double dr = calc_dist(hitp[0], hitp[1], hits2[n][0], hits2[n][1]);
          if (dr < dr_crit){
            save_hit(hitp, hits2[n], mat_expected, mat_observed, rowno);
            rowno++;
          }
        }
      }
    }
  }

  mat_expected.resize(rowno, 3);
  mat_observed.resize(rowno, 3);

  return mat_expected.t() * mat_observed * (mat_observed.t() * mat_observed).i();  // formula derived in Tobias' thesis
}

void analyser::adjust_coordinates(vector<vector<vector<double>>> &Hits, mat mat_T) {
  /* Transform coordinates */
  for (size_t i = 0; i < Hits.size(); i++) {
    for (size_t j = 0; j < Hits[i].size(); j++) {
      vec hit = {{Hits[i][j][0]}, {Hits[i][j][1]}, {1}};
      vec temp = mat_T*hit;
      Hits[i][j][0] = temp[0];
      Hits[i][j][1] = temp[1];
    }
  }
}

void analyser::save_hit(vector<double> hitp, vector<double> hit, mat &mat_expected, mat &mat_observed, int indx) {
    /* Allocate hitoordinates as a row in a matrix */
    rowvec vec_expected = {hitp[0], hitp[1], 1};
    rowvec vec_observed = {hit[0], hit[1], 1};
    mat_expected.row(indx) = vec_expected;
    mat_observed.row(indx) = vec_observed;
}

bool analyser::check_convergence(mat A, mat B) {
  /* Check similarity between inputted matrices */
  return sum(abs(vectorise(A) - vectorise(B))) < tol;
}

void analyser::construct_distarray(void) {
  for (size_t i = 0; i < 4; i++) {
    vector<double> zproj = {zplanes[i], zplanes[i+1], zplanes[i+2]}, dist;
    dist.reserve(Nevents*200);

    calc_interdistance(hitcoords[i], hitcoords[i+1], hitcoords[i+2], zproj, dist);
    distarray[i] = dist;
  }
}

void analyser::align_w_T(void) {
  /* Align planes with T matrix */
  for (size_t i = 0; i < 4; i++) {
    vector<double> zproj = {zplanes[i], zplanes[i+1], zplanes[i+2]}, dist;
    mat mat_T = T.slice(i);
    dist.reserve(Nevents*200);
    for (size_t j = 0; j < hitcoords[i+2].size(); j++) {
      for (size_t k = 0; k < hitcoords[i+2][j].size(); k++) {
        vec hit = {{hitcoords[i+2][j][k][0]}, {hitcoords[i+2][j][k][1]}, {1}};
        vec temp = mat_T*hit;
        hitcoords[i+2][j][k][0] = temp[0];
        hitcoords[i+2][j][k][1] = temp[1];
      }
    }
    cerr << "Aligned plane " << i+2 << "\n";
    calc_interdistance(hitcoords[i], hitcoords[i+1], hitcoords[i+2], zproj, dist);
    distarray[i] = dist;
  }
}

void analyser::align_wo_T(void) {
  /* Align planes without T matrix */
  for (int i = 0; i < 4; i++) {
    mat mat_T = eye<mat>(3, 3);
    vector<double> zproj = {zplanes[i], zplanes[i+1], zplanes[i+2]}, dist;
    align_plane(mat_T, hitcoords[i], hitcoords[i+1], hitcoords[i+2], zproj);
    calc_interdistance(hitcoords[i], hitcoords[i+1], hitcoords[i+2], zproj, dist);
    distarray[i] = dist;
    cerr << "Aligned plane " << i+2 << "\n";
    T.slice(i) = mat_T;
  }

  T.save(DATPATH + "/Align/alignment_matrix.txt", arma_ascii);
}

void analyser::make_grid(vector<vector<double>> &pixelgrid, vector<double> &xgrid, vector<double> &ygrid) {
  /* Construct pixel grid */
  pixelgrid.resize(ncols*nrows);

  for (size_t i = 0; i < pixelgrid.size(); i++){
    pixelgrid[i] = {(double)i, (double)0};
  }

}

void analyser::extract_hit_data(vector<vector<vector<double>>> &hitcoord, vector<vector<double>> &pixelgrid, int plane) {
  /* Save coordinates of hits based on which detector recorded it */
  hitcoord.resize(Nevents);

  for (int i = 0; i < Nevents; i++) {
    vector<vector<double>> eventdata = Events[i];
    for (size_t j = 0; j < eventdata.size(); j++) {
      vector<double> hitdata = eventdata[j];
      if (hitdata[3] == plane) {
        int pixel = coord2pixel(hitdata[0], hitdata[1]);
        pixelgrid[pixel][1] += 1;
        hitcoord[i].push_back({hitdata[0], hitdata[1]});
      }
    }
  }
}

int analyser::coord2pixel(double xhit, double yhit) {
  /* Determine pixelno from coordinates */
  double dx = (xmax - xmin) / ncols;
  double dy = (ymax - ymin) / nrows;
  int colno = (xhit + xmax)/dx;
  int rowno = (yhit + ymax)/dy;
  return rowno * ncols + colno;
}

void analyser::count_hits(int &count, vector<vector<vector<double>>> hitcoord) {
  /* Count total number of hits across all events */
  count = 0;
  for (size_t i = 0; i < hitcoord.size(); i++) {
    for (size_t j = 0; j < hitcoord[i].size(); j++) {
      count++;
    }
  }
}

bool analyser::sortFunc(const vector<double> &p1, const vector<double> &p2) {
 /*  Criteria to sort vector by 2nd column, descending */
 return p1[1] > p2[1];
}

void analyser::locate_hot_pixels(vector<vector<double>> pixelgrid, vector<int> &hotpixels, int i) {

  /*
  We need to identify hot pixels. We do this by sorting the pixelgrid by number of pixels, descendind,
  and then saving the pixelno. of the hotpixels.
  */
  sort(pixelgrid.begin(), pixelgrid.end(), sortFunc);
  double lim = 4E-04 * Nevents;
  for (size_t i = 0; i < pixelgrid.size(); i++) {
    if (pixelgrid[i][1] > lim) {
      hotpixels.push_back(pixelgrid[i][0]);
    } else break;
  }
}

void analyser::remove_hot_pixels(vector<vector<vector<double>>> &hitcoord, vector<vector<double>> &pixelgrid, vector<int> hotpixels) {
  /* We remove both the number of counts in pixelgrid for hotpixels, and erase hitcoordinates from hotpixels. */
  for (size_t i = 0; i < hotpixels.size(); i++) {
    int hotpix = hotpixels[i];
    pixelgrid[hotpix][1] = 0;
    for (size_t j = 0; j < hitcoord.size(); j++) {
      for (size_t k = 0; k < hitcoord[j].size(); k++) {
        double xhit = hitcoord[j][k][0];
        double yhit = hitcoord[j][k][1];
        int pixel = coord2pixel(xhit, yhit);
        if (pixel == hotpix) {
          hitcoord[j].erase(hitcoord[j].begin() + k);
          k--;
        }
      }
    }
  }
}


void analyser::extract_root_data(void) {
  int N_hits_tot = 0;
  int N_hits_max = 0;
  int N_hits_min = 100;

  /* Open root file and obtain relevant data */
  TFile *f = TFile::Open(filename);
  TTree* T1 = (TTree *)f->Get("T");
  TLeaf *Hpk = (TLeaf *)T1->GetLeaf("fAHits.Hpk");
  TLeaf *Hu = (TLeaf *)T1->GetLeaf("fAHits.Hu");
  TLeaf *Hv = (TLeaf *)T1->GetLeaf("fAHits.Hv");
  TLeaf *fAHitsN = (TLeaf *)T1->GetLeaf("fAHitsN");
  TLeaf *ENumberOfTriggers = (TLeaf *)T1->GetLeaf("fHeader.ENumberOfTriggers");

  int N_events = T1->GetEntries();
  cerr << "\nTotal number of events :\t" << N_events << "\n";

  /* Extract data from each event */
  for (int i  = 0; i < N_events; i++) {
    T1->GetEntry(i);
    int N_hits = fAHitsN->GetValue();
    int N_triggers = ENumberOfTriggers->GetValue();
    if (N_triggers == 2){
      vector<vector<double>> EventData(N_hits);
      int N_hits_0 = 0;
      for (int j = 0; j < N_hits; j++) {
        vector<double> HitData(4);
        HitData[0] = Hu->GetValue(j);
        HitData[1] = Hv->GetValue(j);
        HitData[2] = 1.0;
        HitData[3] = static_cast<double>(Hpk->GetValue(j) -1);

        if (HitData[3] == 0) {
          N_hits_tot += 1;
          N_hits_0 += 1;
          if (N_hits_0 > N_hits_max) {
            N_hits_max = N_hits_0;
          }
          if (N_hits_0 < N_hits_min) {
            N_hits_min = N_hits_0;
          }
        }

        EventData[j] = HitData;
      }
      Events.push_back(EventData);
    }
  }

  Nevents = Events.size();
  cerr << "Number of usable events :\t" << Nevents << "\t average hits pr. event in plane0 :\t" << (double)N_hits_tot/(double)Nevents << "\t max hits :\t" << N_hits_max << "\t min hits :\t" << N_hits_min << "\n\n";
}
