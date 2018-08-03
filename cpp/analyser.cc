#include "analyser.h"

analyser::analyser(vector<double> z, char const *name, char const *run){
  /* Assign values to class attributes */
  nrows = 576;  // no. of rows in mimosa detector
  ncols = 1152; // no. of coloumns in mimosa detector
  xmin = -11000, xmax = 11000, ymin = -5500, ymax = 5500; // dimension of grid for Mimosa detectors [micro-meter]
  runno = run;
  M1M2_d_lim = 2500.0;  // acceptable distance between M1->M2 projected point and observed hit in M3 [micro-meter]
  M2M3_d_lim = 1000.0; // acceptable distance between M2->M3 projected point and observed hit in M4 [micro-meter]
  M6M5_d_lim = 500.0; // acceptable distance in MM bewteen M1->MM and M6->MM arm of track [micro-meter]
  Match_d = 250.0; // matching ditance in MM between electron/positron paired track [micro-meter]
  Match_d_foil = 250.0; // matching ditance in foil between electron/positron paired track [micro-meter]
  yz_defl_lim = 5.0e-3; // angle of acceptable y-deflction [rad]
  dr_crit_list = {3500, 2500, 1500, 500, 250, 150};  // radius for acceptable hits when aligning [micro-meter]
  tol = 1e-6; // convergence tolerance for (absolute) difference between alignment matrices for individual plane when aligning
  zplanes = z;  // z-coords of planes
  M1M2_z = {zplanes[0], zplanes[1], zplanes[2]};  // z-coords for M1-M3 PLANES
  M2M3_z = {zplanes[1], zplanes[2], zplanes[3]};  // z-coords for M2-M4 planes
  M3M4_z = {zplanes[2], zplanes[3], (zplanes[4]+zplanes[3])/2.0}; // z-coords of M3-MM planes
  M6M5_z = {zplanes[4], zplanes[5], (zplanes[4]+zplanes[3])/2.0}; // z-coords for M5-MM planes
  filename = name; // name of data file (root file)
  DATPATH = "/home/christian/Dropbox/speciale/data";  // path to save results

  /* Pre-allocate size of vectors */
  hitcoords.resize(6);
  pixelgrids.resize(6);
  hotpixels.resize(6);
  T.set_size(3, 3, 4);
  distarray.resize(4);
}

void analyser::update_pixelgrids(int plane, vector<vector<double>> pixelgrid){
  pixelgrids[plane] = pixelgrid;
}

void analyser::update_hitcoords(int plane, vector<vector<vector<double>>> hitcoord){
  hitcoords[plane] = hitcoord;
}

void analyser::update_hotpixels(int plane, vector<int> hotpixel){
  hotpixels[plane] = hotpixel;
}

void analyser::print_pixels(void){
  string filename = DATPATH + "/pixeldata_run" + runno + ".txt";
  ofstream output (filename);
  for (size_t i = 0; i < pixelgrids.size(); i++){
    output << "## PLANE\t" << i << "\tPIXEL DATA\n";  // header
    for (size_t j = 0; j < pixelgrids[i].size(); j++){
      output << pixelgrids[i][j][1] << "\n";
    }
  }
}

void analyser::print_hotpixels(void){
  string filename = DATPATH + "/hotpixels_run" + runno + ".txt";
  ofstream output (filename);
  for (size_t i = 0; i < hotpixels.size(); i++){
    output << "## PLANE\t" << i << "\tPIXEL DATA\n";  // header
    for (size_t j = 0; j < hotpixels[i].size(); j++){
      output << hotpixels[i][j] << "\n";
    }
  }
}

void analyser::print_interdistance(void){
  string filename = DATPATH + "/interdistance_data_" + runno + ".txt";;
  ofstream output (filename);
  for (size_t i = 0; i < distarray.size(); i++){
    output << "## DATABLOCK\t" << i << "\tINTERDISTANCE DATA\n";  // header
    for (size_t j = 0; j < distarray[i].size(); j++){
      output << distarray[i][j] << "\n";
    }
    output << "\n\n";
  }
}

void analyser::print_hits(void){
  string filename = DATPATH + "/hits_coord_data_" + runno + ".txt";
  ofstream output (filename);
  for (size_t i = 0; i < hitcoords.size(); i++){
    output << "## PLANE\t" << i << "\tHIT DATA\n";  // header
    for (size_t j = 0; j < hitcoords[i].size(); j++) {
      for (size_t k = 0; k < hitcoords[i][j].size(); k++){
      output << hitcoords[i][j][k][0] << ' ' << hitcoords[i][j][k][1] << '\n';
      }
    }
  }
}

void analyser::print_energy(void){
  string filename = DATPATH + "/energy_" + runno + ".txt";
  ofstream output (filename);
  for (size_t i = 0; i < energies.size(); i++){
    output << energies[i] << "\n";
  }
}

void analyser::print_slope(void){
  string filename = DATPATH + "/angles_" + runno + ".txt";
  ofstream output (filename);
  for (size_t i = 0; i < slopes.size(); i++){
    output << slopes[i][0] << " " << slopes[i][1] << "\n";
  }
}

void analyser::print_M1M2_slope(void){
  string filename1 = DATPATH + "/angles_M1M2_" + runno + ".txt";
  string filename2 = DATPATH + "/beam_divergence_"  + runno + ".txt";
  ofstream output1 (filename1);
  ofstream output2 (filename2);
  for (int i = 0; i < Nevents; i++){
    for (size_t j = 0; j < M1M2_slopes[i][0].size(); j++){
      output1 << M1M2_slopes[i][0][j] << "\t" << M1M2_slopes[i][1][j] << "\n";
    }
    for (size_t j = 0; j < divergence[i][0].size(); j++){
      output2 << divergence[i][0][j] << "\t" << divergence[i][1][j] << "\n";
    }
  }
}

void analyser::print_zpos(void){
  string filename = DATPATH + "/zclosepos_" + runno + ".txt";
  ofstream output (filename);
  for (size_t i = 0; i < zclosepos.size(); i++){
    output << zclosepos[i] << "\n";
  }
}

double analyser::calc_ang(double m, double n){
  return n - m;  // small angle approx. => angle is the difference between the slopes
}

double analyser::calc_slope(vec x, vec y){
    return (y(1) - y(0))/(x(1) - x(0));
}

mat analyser::lines3d_nearestpoints(vec A, vec B, vec C, vec D){
  /* Calculate the closest approach of two lines in 3D
        A,B,C,D : vector from origo to points A and B on line 1 and C and D on line 2
  */
  mat closepos(3, 2);
  vec d1 = B - A;  // directional vector for line 1
  vec d2 = D - C;  // directional vector for line 2
  vec n = cross(d1, d2);  // used to calculae n1, n2
  vec n1 = cross(d1, n);  // vector perpindicular to plane formed by translating line 1 along n, used to calculate points on line 1 nearest line 2
  vec n2 = cross(d2, n);  // vector perpindicular to plane formed by translating line 2 along n, used to calculatee points on line 2 nearest line 1
  vec c1 = A + dot(C - A, n2) / dot(d1, n2) * d1; // points on line 2 nearest line 1
  vec c2 = C + dot(A - C, n1) / dot(d2, n1) * d2; // points on line 1 nearest line 2
  closepos.col(0) = c1;
  closepos.col(1) = c2;
  return closepos;
}

void analyser::image_crystal(void){
  /* This function reconstructs an image of the crystal using the tracks obtained from "pair tracks" function */
  string filename = DATPATH + "/crystal_image_" + runno + ".txt";
  ofstream output (filename);
  for (size_t i = 0; i < paired_tracks.size(); i++){
    output << paired_tracks[i][0][1][0] << "\t" << paired_tracks[i][0][1][1] << "\n" << paired_tracks[i][1][1][0] << "\t" << paired_tracks[i][1][1][1] << "\n";
  }
}

void analyser::calc_pair_energy(double &energy,vector<vector<vector<double>>> pairedtracks){
  double L = 0.15, q = 1.6021766e-19, B = 0.12, c = 299792458; // SI units q = charge of positron, B = mag. field strength, c = speed of light in vac. L = length traveled through B field
  energy = 0;
  for (int i = 0; i < 2; i++){
    /* Define arma-vectors since calc_slope takes arma-vectors as input  */
    vec M4MM_x = {pairedtracks[i][3][0], pairedtracks[i][4][0]};
    vec M4MM_z = {zplanes[3], (zplanes[4]+zplanes[3])/2};
    vec M5M6_x = {pairedtracks[i][6][0], pairedtracks[i][7][0]};
    vec M5M6_z = {zplanes[4], zplanes[5]};
    double m = calc_slope(M4MM_z, M4MM_x);  // slope of linesegmet
    double n = calc_slope(M5M6_z, M5M6_x);
    vector<double> slope = {m, n};
    slopes.push_back(slope);  // push-back is slow, but not a huge issue
    double ang = calc_ang(m, n);  // calculate deflection angle to find energy
    energy += q*c*L*B/abs(ang); // calculate energy as sum of both particles energy in Joule
  }
}

void analyser::pair_tracks(void){
  double dist, nrg;
  int pairs = 0, matchedtrack = 0, tot_matched_tracks = 0, erased_tracks = 0;
  // vector<vector<vector<vector<double>>>> paired_tracks;
  paired_tracks.resize(tracks.size()*200);
  energies.resize(tracks.size()*200);
  zclosepos.resize(tracks.size()*200);

  for (size_t i = 0; i < tracks.size(); i++){
    int matchedtracks = 0;  // for indivdual track, zero matchedtracks counter since we only want 1 matched track for each track

    for (size_t j = i+1; j < tracks.size(); j++){  // iterate from i+1'th track. This ensures that we don't iterate over discarded tracks, without having to remove the i'th track from tracks.
      if (tracks[i][8][0] == tracks[j][8][0]){  // ensure we only look at tracks from the same event
        dist = calc_dist(tracks[i][4][0], tracks[i][4][1], tracks[j][4][0], tracks[j][4][1]);  // calculate distance between projected hit in mimosa-magnet for the two tracks
        if (dist < Match_d){
          double M3_dist = calc_dist(tracks[i][2][0], tracks[i][2][1], tracks[j][2][0], tracks[j][2][1]);
          if (M3_dist < Match_d_foil){  // don't proceed if tracks are far in M3 (same as foil) coordinates
            matchedtracks++;  // increase matchedtracks counter
            tot_matched_tracks++; // update total number of matchedtracks
            if (matchedtracks > 1){ // remove additional tracks if more than one is found
              tracks.erase(tracks.begin() + j); // erase additional track
              erased_tracks++;
              j--; // decrease iterator since tracks vector now has fewer elements than initially
            }
            else if (matchedtracks == 1){ // if only one is found then save the track number (j)
              matchedtrack = j;
            }
          }
        }
      }
    } // end for (tracks j)

    if (matchedtracks == 1){  // if only 1 track is found, then accept the track (j) as a pair (for i) and calculate energy
      vec P1 = {tracks[i][6][0], tracks[i][6][1], M6M5_z[0]}, Q1 = {tracks[i][7][0], tracks[i][7][1], M6M5_z[1]};
      vec P2 = {tracks[matchedtrack][2][0], tracks[matchedtrack][2][1], M3M4_z[0]}, Q2 = {tracks[matchedtrack][3][0], tracks[matchedtrack][3][1], M6M5_z[1]};
      mat closepos = lines3d_nearestpoints(P1, Q1, P2, Q2);
      double zclosepospair = (closepos(2, 0) + closepos(2, 1))/2;
      if (/*zclosepospair > 8e6 && zclosepospair < 1e7*/ true){ // throw away tracks if they don't meet near the middle of Mimosa magnet. This is a relatively slow calculation (requires calculating 4 cross-products), so perhaps there is a better way.... a brighter tomorrow.
        paired_tracks[pairs] = {tracks[i], tracks[matchedtrack]};
        calc_pair_energy(nrg, paired_tracks[pairs]);
        energies[pairs] = nrg;
        zclosepos[pairs] = zclosepospair; // this is saved so we can plot the closest position of our acceptated tracks
        pairs++;  // increase paired-tracks counter, used to later resize vectors as well as reported to terminal
      }
    }

  } // end for (tracks i)

  /* Re-size vectors (removes zeros) */
  paired_tracks.resize(pairs);
  energies.resize(pairs);
  zclosepos.resize(pairs);

  /* Report progress to terminal */
  cerr << "Number of paired tracks : " << pairs << "\n";
}

void analyser::beam_divergence(int eventno, int x_indx, int y_indx, string name){
  /* Calculate beam divergence in x- and y-direction */
  string str1 = "partial";
  string str2 = "full";
  vec M1M2_x = {hitcoords[0][eventno][x_indx][0], hitcoords[1][eventno][y_indx][0]};  // x-coords of M1-M2 track
  vec M1M2_y = {hitcoords[0][eventno][x_indx][1], hitcoords[1][eventno][y_indx][1]};  // y-coords of M1-M2 track
  vec vec_M1M2_z = {M1M2_z[0], M1M2_z[1]};  // z-coords of M1-M2 track
  double ang_x = calc_slope(vec_M1M2_z, M1M2_x); // we use this to enable analysis of tracks with specific entry angle
  double ang_y = calc_slope(vec_M1M2_z, M1M2_y);
  if (name.compare(str2) == 0){
    divergence[eventno][0].push_back(ang_x); // push_back is slow, but that is presently not a big issue
    divergence[eventno][1].push_back(ang_y);
  }
  if (name.compare(str1) == 0){
    M1M2_slopes[eventno][0].push_back(ang_x); // push_back is slow, but that is presently not a big issue
    M1M2_slopes[eventno][1].push_back(ang_y);
  }
}

void analyser::construct_tracks(void){
  /* Set counters to 0 and pre-allocate space for collecing slopes + beam divergence  */
  int M1_MM_tot_tracks = 0, M1_M6_tot_tracks = 0; // for counting number of tracks
  double M1M2_slope_lb = -INFINITY, M1M2_slope_ub = INFINITY; // upper and lower bounds for data cuts
  M1M2_slopes.resize(Nevents);
  divergence.resize(Nevents);
  string str1 = "full", str2 = "partial";

    /* Construct tracks for each event and save in "tracks" array */
    for (int i = 0; i < Nevents; i++){
      vector<vector<double>> M1_M6_track; // for saving final track
      M1M2_slopes[i].resize(2);
      divergence[i].resize(2);
      for (size_t j = 0; j < hitcoords[0][i].size(); j++){  // iterate over hits in M1 for i'th event to construct M1-M2 track
        for (size_t k = 0; k < hitcoords[1][i].size(); k++){  // iterate over hits in M2 for i'th Event to construct M1-M2 tracj
          vector<double> M1M2_Proj = rect_project(hitcoords[0][i][j], hitcoords[1][i][k], M1M2_z);  // project M1-M2 track into M3
          beam_divergence(i, j, k, str1); // calculate beam divergence so we can check the beam profile
          for (size_t l = 0; l < hitcoords[2][i].size(); l++){  // iterate over hits in M3 for i'th Event to construct M2-M3 track
            double M1M2_d = calc_dist(M1M2_Proj[0], M1M2_Proj[1], hitcoords[2][i][l][0], hitcoords[2][i][l][1]);
            if (M1M2_d < M1M2_d_lim){ // if the distance between projected and observed hit is small enough, begin constructing M1-M2 -> M3-M4 track
              vector<double> M2M3_Proj = rect_project(hitcoords[1][i][k], hitcoords[2][i][l], M2M3_z); // project into M4
              for (size_t m = 0; m < hitcoords[3][i].size(); m++){  // iterate over hits in M4 for i'th Event to construct M3-M4 track
                double M2M3_d = calc_dist(M2M3_Proj[0], M2M3_Proj[1], hitcoords[3][i][m][0], hitcoords[3][i][m][1]);
                vec vec_M3M4_z = {M3M4_z[0], M3M4_z[1]};  // z-coords of M3-M4 track
                vec M3M4_y = {hitcoords[2][i][l][1], hitcoords[3][i][m][1]};  // y-coords of M3-M4 track
                double M3M4_slope = calc_slope(vec_M3M4_z, M3M4_y); // slope in y-direction of M3-M4 track, calculated so we can later determine the total y-deflection of a track as it goes through the Mimosa magnet
                beam_divergence(i, j, k, str2); // calculate beam divergence so we can isolate tracks with given angle
                if (M1M2_slopes[i][0][j] > M1M2_slope_lb and M1M2_slopes[i][0][j] < M1M2_slope_ub){
                  if (M2M3_d < M2M3_d_lim){ // proceed from M1-M4 -> MM
                    vector<double> M3M4_Proj = rect_project(hitcoords[2][i][l], hitcoords[3][i][m], M3M4_z);  // project into MM
                    M1_MM_tot_tracks++;
                    for (size_t n = 0; n < hitcoords[4][i].size(); n++){  // iterate over hits in M5 so we can make M6-MM tracks
                      double shortdist = INFINITY;  // reset this for every hit (possible track), since we want the best track for every M6-M5 pair
                      bool trackfound = false;  // initially there is no track found, this is changed to "true" if a track is found
                      for (size_t o = 0; o < hitcoords[5][i].size(); o++){  // iterate over hits in M6
                        vector<double> M6M5_Proj = rect_project(hitcoords[4][i][n], hitcoords[5][i][o], M6M5_z);  // project into MM
                        double M6M5_d = calc_dist(M6M5_Proj[0], M6M5_Proj[1], M3M4_Proj[0], M3M4_Proj[1]);  // calculate distance in xy-plane between M1->MM and M6->MM projection
                        if (M6M5_d < shortdist and M6M5_d < M6M5_d_lim){  // if we have the best hit that also corresponds well with the projection, proceed
                          vec vec_M5M6_z = {M6M5_z[0], M6M5_z[1]};  // z-coords of M5-M6 track, used to calculate y-slope
                          vec M5M6_y = {hitcoords[4][i][n][1], hitcoords[5][i][o][1]};  // y-coords of M5-M6 track, used to calculate y-slope
                          double M5M6_slope = calc_slope(vec_M5M6_z, M5M6_y); // slope in y-direction of M5-M6 track, used to ensure minimal y-deflection
                          if (abs(M5M6_slope - M3M4_slope) < yz_defl_lim){ // if y-deflection is too large dont proceed, since manget only deflects in x-dir.
                            shortdist = M6M5_d; // update shortest distance so we only find the closest hit
                            trackfound = true;  // change boolean since track is found
                            M1_M6_track = {hitcoords[0][i][j], hitcoords[1][i][k], hitcoords[2][i][l], hitcoords[3][i][m], M3M4_Proj, M6M5_Proj, hitcoords[4][i][n], hitcoords[5][i][o], {(double)i}}; // save track: { M1(x,y), M2(x, y), M3(x, y), M4(x, y) MM_34(x, y) MM_56(x, y), M5(x, y) M6(x, y), eventno} where fx. M1(x,y) means (x, y) coords. in M1-plane
                          }
                        }
                      } // hits in M6 done
                      if (trackfound){
                        tracks.push_back(M1_M6_track); // update 'tracks' vector
                        M1_M6_tot_tracks++; // increase count of tracks
                      } // saving tracks done
                    } // hits in M5 done
                  } // M1-M4 -> MM done
                } // cuts done
              } // hits in M4 done
            } // if (distance good) -> construct M1M2->M3M4 done
          } // hits in M3 done
        } // hits in M2 done
      } // hits in M1 done
    } // events done

/* Report progress to terminal */
cerr << "\nTotal tracks from M1 -> MM : " << M1_MM_tot_tracks << "\n";
cerr << "Total tracks from M1 -> M6 : " << M1_M6_tot_tracks << "\n";
}

void analyser::calc_interdistance(vector<vector<vector<double>>> Hits0, vector<vector<vector<double>>> Hits1, vector<vector<vector<double>>> Hits2, vector<double> z, vector<double> &distvec){
  /* Calculate distances between projected and observes hits */
  distvec.resize(500*Nevents);
  int indx = 0; // this index counts the number of entries in "distvec"

  /* Project hits into neighbor plane */
  for (int i = 0; i < Nevents; i++){
    vector<vector<double>> hits0 = Hits0[i];
    vector<vector<double>> hits1 = Hits1[i];
    vector<vector<double>> hits2 = Hits2[i];
    for (size_t j = 0; j < hits0.size(); j++){

      /* Project hits */
      for (size_t k = 0; k < hits1.size(); k++){
        vector<double> hitp = rect_project(hits0[j], hits1[k], z);

        /* Calculate distances and save in "distvec" array */
        for (size_t l = 0; l < hits2.size(); l++){
          double newdist = calc_dist(hitp[0], hitp[1], hits2[l][0], hits2[l][1]); // calculate distance between projected and observed hit.
          distvec[indx] = newdist;
          indx++;
        }
      }
    }
  }

  /* Resize array to save space */
  distvec.resize(indx);
}

double analyser::calc_dist(double x0, double y0, double x1,double y1){
  return sqrt(pow(x1 - x0, 2) + pow(y1 - y0, 2));
}

vector<double> analyser::rect_project(vector<double> hit0, vector<double> hit1, vector<double> z){
  vector<double> proj(2);
  proj[0] = ((hit0[0]-hit1[0])/(z[0]-z[1]))*(z[2]-z[0]) + hit0[0];  // rectilinear projection, x coord
  proj[1] = ((hit0[1]-hit1[1])/(z[0]-z[1]))*(z[2]-z[0]) + hit0[1];  // rectilinear projection, z coord
  return proj;
}

void analyser::align_plane(mat &mat_Tot, vector<vector<vector<double>>>  Hits0, vector<vector<vector<double>>>  Hits1, vector<vector<vector<double>>>  &Hits2, vector<double> z){
  /* Iterate over allowed distance between projected and observed points */
  for (size_t j = 0; j < dr_crit_list.size(); j++) {
    double dr_crit = dr_crit_list[j];
    mat mat_temp_T = zeros<mat>(3,3);
    bool convergence = false;
    mat mat_T;

    /* While not convergent, construfct and update the tranformation matrix T */
    while (!convergence) {
      construct_T_mat(Hits0, Hits1, Hits2, mat_T, dr_crit, z);
      mat_Tot = mat_T * mat_Tot;  // the total transformation matrix, in the sense that this is the matrix we should multiply onto the observed hits in a plane
      adjust_coordinates(Hits2, mat_T); // adjust the observed hits using the transformation matrix we have just calculated. This must be done here, since it affects how mat_T is calculated, and therefore if convergence is reached.
      convergence = check_convergence(mat_T, mat_temp_T); // check convergence, ie. has the tranformation matrix changed significantly this iteration
      mat_temp_T = mat_T;
    }

  }
}

void analyser::construct_T_mat(vector<vector<vector<double>>> Hits0, vector<vector<vector<double>>> Hits1, vector<vector<vector<double>>> Hits2, mat &mat_T, double dr_crit,  vector<double> z){
  /* This function is almost identical to "calc_interdistance" */
  int rowno = 0; // number of rows in matrices mat_expected, mat_observed
  mat mat_expected = zeros<mat>(Hits1.size()*200, 3), mat_observed = zeros<mat>(Hits1.size()*200,3); // allocate arbitrarily large to avoid segmentation fault
  for (size_t k = 0; k < Hits0.size(); k++){
    vector<vector<double>> hits0 = Hits0[k];
    vector<vector<double>> hits1 = Hits1[k];
    vector<vector<double>> hits2 = Hits2[k];
    for (size_t l = 0; l < hits0.size(); l++){
      for (size_t m = 0; m < hits1.size(); m++){
        vector<double> hitp = rect_project(hits0[l], hits1[m], z); // project into final plane
        for (size_t n = 0; n < hits2.size(); n++){
          double dr = calc_dist(hitp[0], hitp[1], hits2[n][0], hits2[n][1]);
          if (dr < dr_crit){  // if the distance between the observed hit and the projected hit is small enough, we save both in matrices, and use the matrices to construct the tranformation matrix T
            save_hit(hitp, hits2[n], mat_expected, mat_observed, rowno);
            rowno++;
          }
        }
      }
    }
  }
  mat_expected.resize(rowno, 3);
  mat_observed.resize(rowno, 3);
  mat_T = mat_expected.t() * mat_observed * (mat_observed.t() * mat_observed).i();  // formula derived in Tobias' thesis
}

void analyser::adjust_coordinates(vector<vector<vector<double>>> &Hits, mat mat_T){
  for (size_t i = 0; i < Hits.size(); i++){
    for (size_t j = 0; j < Hits[i].size(); j++){
      vec hit = {{Hits[i][j][0]}, {Hits[i][j][1]}, {1}};  // use arma-vector to store observed hit, with z-coord = 1, so we can use arma-library for matrix-multiplication
      vec temp = mat_T*hit;
      Hits[i][j][0] = temp[0];  // update "Hits" with the adjusted hit's x-coord
      Hits[i][j][1] = temp[1];  // update "Hits" witht eh adjusted hit's y-coord
    }
  }
}

void analyser::save_hit(vector<double> hitp, vector<double> hit, mat &mat_expected, mat &mat_observed, int indx){
    /* This function simply allocates a hit's coordinates as a row in a matrix. This is done so we can later use these matrices to calculate the coordinate-tranformation matrix which ultimately aligns our planes */
    rowvec vec_expected = {hitp[0], hitp[1], 1};
    rowvec vec_observed = {hit[0], hit[1], 1};
    mat_expected.row(indx) = vec_expected;
    mat_observed.row(indx) = vec_observed;
}

bool analyser::check_convergence(mat mat_T, mat mat_temp_T){
  cerr << sum(abs(vectorise(mat_T) - vectorise(mat_temp_T))) << "\n";
  return sum(abs(vectorise(mat_T) - vectorise(mat_temp_T))) < tol; // this is used to determine whether we have converged for a given "dr" or not. Function is called in "align_plane" function
}

void analyser::construct_distarray(void){
  for (size_t i = 0; i < 4; i++){ // iterate over the planes we wish to align
    vector<double> zproj = {zplanes[i], zplanes[i+1], zplanes[i+2]}, dist;  // zproj are used to calculate distance between projected and observed points
    dist.reserve(Nevents*200); // arbitrarily large to avoid segmentation fault
    calc_interdistance(hitcoords[i], hitcoords[i+1], hitcoords[i+2], zproj, dist); // calculate distance between projected and observed hits in order to check if alignment is done correctly. This distance should be close to zero
    distarray[i] = dist;
  }
}

void analyser::align_w_T(void){
  /* This function aligns planes with the T matrix already determined. Use this for all runs other than the alignment run */
  for (size_t i = 0; i < 4; i++){ // iterate over the planes we wish to align
    vector<double> zproj = {zplanes[i], zplanes[i+1], zplanes[i+2]}, dist;  // zproj are used to calculate distance between projected and observed points
    mat mat_T = T.slice(i); // the alignment matrices are saved in an arma-cube. A "slice" of a cube is a matrix
    dist.reserve(Nevents*200); // arbitrarily large to avoid segmentation fault
    for (size_t j = 0; j < hitcoords[i+2].size(); j++){ // we only want to align planes M3,M4,M5,M6, ie. the 2,3,4,5 entry in hitcoords, therefore we iterate over i+2 with i running from 0 to 3
      for (size_t k = 0; k < hitcoords[i+2][j].size(); k++){  // align hit-by-hit
        vec hit = {{hitcoords[i+2][j][k][0]}, {hitcoords[i+2][j][k][1]}, {1}};
        vec temp = mat_T*hit; // align
        hitcoords[i+2][j][k][0] = temp[0];  // update hit's x-coord
        hitcoords[i+2][j][k][1] = temp[1];  // update hit's y-coord
      }
    }
    cerr << "Aligned plane " << i+2 << "\n";  // report progress to terminal
    calc_interdistance(hitcoords[i], hitcoords[i+1], hitcoords[i+2], zproj, dist); // calculate distance between projected and observed hits in order to check if alignment is done correctly. This distance should be close to zero
    distarray[i] = dist;
  }
}

void analyser::align_wo_T(void){
  /* This function aligns planes without the T matrix determined. Use this only for the alignement run */
  for (int i = 0; i < 4; i++){
    mat mat_T = eye<mat>(3, 3);
    vector<double> zproj = {zplanes[i], zplanes[i+1], zplanes[i+2]}, dist;  // z-coordinates used for projection and calculating distance between projected and observed points
    align_plane(mat_T, hitcoords[i], hitcoords[i+1], hitcoords[i+2], zproj);
    calc_interdistance(hitcoords[i], hitcoords[i+1], hitcoords[i+2], zproj, dist); // calculate distance between projected and observed hits in order to check if alignment is done correctly. This distance should be close to zero
    distarray[i] = dist;
    cerr << "Aligned plane " << i+2 << "\n";  // report progress to terminal
    T.slice(i) = mat_T; // save transformation matrix in arma-cube so we can align planes for data-runs without running alignment algorithm every time
  }
  T.save(DATPATH + "/Align/alignment_matrix.txt", arma_ascii);  // save arma-cube in txt-file
}

void analyser::make_grid(vector<vector<double>> &pixelgrid, vector<double> &xgrid, vector<double> &ygrid){
  pixelgrid.resize(ncols*nrows);  // contains data for each pixel in grid
  /* Construct grid */
  for (size_t i = 0; i < pixelgrid.size(); i++){
    pixelgrid[i] = {(double)i, (double)0};  // first entry is pixel number (starting top-left and counting rightward, like a matrix), second entry is number of hits in that pixel
  }
}

void analyser::extract_hit_data(vector<vector<vector<double>>> &hitcoord, vector<vector<double>> &pixelgrid, int plane){
  // EXTRACT COORDINATES FOR HITS
  hitcoord.resize(Nevents);
  for (int i = 0; i < Nevents; i++) {
    vector<vector<double>> eventdata = Events[i]; // we will use this to iterate over hits in event
    for (size_t j = 0; j < eventdata.size(); j++) {
      vector<double> hitdata = eventdata[j];  // contains x-, y-coord and plane for a given hit, we use this to fill pixelgrid and hitcoord vectors
      if (hitdata[3] == plane) {
        int pixel;
        coord2pixel(hitdata[0], hitdata[1], pixel); // finds the appropriate pixel given a set of x-y coordinates, used to fill pixelgrid vector
        pixelgrid[pixel][1] += 1; // we now have a hit in the pixel found using coord2pixel, so we increase the hitcount by 1. This pixelgrid is later saved in another vector "pixelgrids" to have indivdual pixelgrid for each plane
        hitcoord[i].push_back({hitdata[0], hitdata[1]});  // we push_back the x,y coordinate to the i'th entry. This keeps track of which plane has which hits
      } // end if-statement for plane
    } // end hits
  } // end events
}

void analyser::coord2pixel(double xhit, double yhit, int &pixel){
  double dx = (xmax - xmin) / ncols;  // the width of a column
  double dy = (ymax - ymin) / nrows;  // the width of a row
  int colno = (xhit + xmax)/dx;  // the column number for a given xhit. Casting as int rounds down, which we want since we count from 0 -> ncols.
  int rowno = (yhit + ymax)/dy; // the row number for a given xhit.
  pixel = rowno * ncols + colno;  // the pixelno given this row and column with the 0'th pixel as the top-right corner, and counting rightward (like a matrix)
}

void analyser::count_hits(int &count, vector<vector<vector<double>>> hitcoord){
  count = 0;  // we use this to output the total number of hits in an event to the terminal. Not used for data-analysis
  for (size_t i = 0; i < hitcoord.size(); i++){
    for (size_t j = 0; j < hitcoord[i].size(); j++){
      count++;
    }
  }
}

bool analyser::sortFunc(const vector<double> &p1, const vector<double> &p2) {
 return p1[1] > p2[1];  // sorts an std::vector by 2nd column, descending
}

void analyser::locate_hot_pixels(vector<vector<double>> pixelgrid, vector<int> &hotpixels, int i){
  sort(pixelgrid.begin(), pixelgrid.end(), sortFunc); // this sorts the pixelgrid by number of pixels, descending
  double lim = 4E-04 * Nevents;
  for (size_t i = 0; i < pixelgrid.size(); i++){
    if (pixelgrid[i][1] > lim){
      hotpixels.push_back(pixelgrid[i][0]);
    }
    else {
      break;
    }
  }
}

void analyser::remove_hot_pixels(vector<vector<vector<double>>> &hitcoord, vector<vector<double>> &pixelgrid, vector<int> hotpixels){
  for (size_t i = 0; i < hotpixels.size(); i++){  // iterate over number of hotpixels in order to remove them
    int hotpix = hotpixels[i];
    pixelgrid[hotpix][1] = 0;
    for (size_t j = 0; j < hitcoord.size(); j++){ // iterate over all hits in order to locate the hits in hotpixels
      for (size_t k = 0; k < hitcoord[j].size(); k++){  // iterate over hitdata for a given hit
        double xhit = hitcoord[j][k][0];
        double yhit = hitcoord[j][k][1];
        int pixel;
        coord2pixel(xhit, yhit, pixel); // here we figure out what pixel the hit is in
        if (pixel == hotpix){
          hitcoord[j].erase(hitcoord[j].begin() + k); // erase the hit if it is in a hotpixel
          k--;  // decrease k since we now have one entry less in hitcoord[j] vector
        }
      }
    }
  }
}


void analyser::extract_root_data(void){
  // EXTRACT DATA FROM ROOT FILE AND SAVE IN 'EVENTS'
  int N_hits_tot = 0;
  int N_hits_max = 0;
  int N_hits_min = 100;
  TFile *f = TFile::Open(filename);
  TTree* T1 = (TTree *)f->Get("T");
  TLeaf *Hpk = (TLeaf *)T1->GetLeaf("fAHits.Hpk");
  TLeaf *Hu = (TLeaf *)T1->GetLeaf("fAHits.Hu");
  TLeaf *Hv = (TLeaf *)T1->GetLeaf("fAHits.Hv");
  TLeaf *fAHitsN = (TLeaf *)T1->GetLeaf("fAHitsN");
  TLeaf *ENumberOfTriggers = (TLeaf *)T1->GetLeaf("fHeader.ENumberOfTriggers");
  int N_events = T1->GetEntries();  // no. of  scintilator detections, ie events. This is used extensively since we always iterate over indivdual events
  cerr << "\nTotal number of events :\t" << N_events << "\n";
  for (int i  = 0; i < N_events; i++) {
    T1->GetEntry(i);
    int N_hits = fAHitsN->GetValue();  // no. of mimosa detetctions, ie hits, in a given event
    int N_triggers = ENumberOfTriggers->GetValue(); // N_triggers = 2 -> 1, we extract this since we are only interested in events with 1 trigger
    if (N_triggers == 2){
      vector<vector<double>> EventData(N_hits); // vector 'EventData' with data for each hit in each event collected in a vector, later collected in Events vector which will contain data for all events
      int N_hits_0 = 0;
      for (int j = 0; j < N_hits; j++) {
        vector<double> HitData(4);
        HitData[0] = Hu->GetValue(j);  // x-coord.
        HitData[1] = Hv->GetValue(j);  // y-coord.
        HitData[2] = 1.0;  // z-coord.
        HitData[3] = static_cast<double>(Hpk->GetValue(j) -1);  // plane
        if (HitData[3] == 0) {
          N_hits_tot += 1;
          N_hits_0 += 1;
        }

        if (N_hits_0 > N_hits_max) {
          N_hits_max = N_hits_0;
        }

        if (N_hits_0 < N_hits_min) {
          N_hits_min = N_hits_0;
        }

        EventData[j] = HitData;
      }
      Events.push_back(EventData);  // push_back is slow, but it is not a huge issue here
    }
  }
  Nevents = Events.size();
  cerr << "Number of usable events :\t" << Nevents << "\t average hits pr. event in plane0 :\t" << (double)N_hits_tot/(double)Nevents << "\t max hits :\t" << N_hits_max << "\t min hits :\t" << N_hits_min << "\n\n"; // report progress to terminal
}
