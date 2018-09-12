#ifndef ANALYSER_H
#define ANALYSER_H

// EXT. LIBRARYS
#include <iostream>
#include <vector>
#include <algorithm>
#include <math.h>
#include <armadillo>
#include <fstream>
#include <omp.h>
#include <string>

// ROOT LIBRARYS
#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2.h"
#include "TH2F.h"
#include "TFile.h"
#include "TVersionCheck.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TObject.h"
#include "TStorage.h"
#include "TSpectrum.h"
#include "TVirtualFitter.h"
#include "TMath.h"
#include "TRandom.h"
#include "TF1.h"
#include "TGraph.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"

using namespace std;
using namespace arma;

class analyser {

  public:

    cube T; // armadillo cube for collecting alignment matrices (rank-3 tensor)
    int Nevents;  // total number of events
    vector<vector<vector<double>>> Events;  // containts root data
    vector<vector<vector<vector<double>>>> hitcoords;

    analyser(vector<double> z, char const *name, char const *runno);  // constructor
    void make_grid(vector<vector<double>> &pixelgrid, vector<double> &xgrid, vector<double> &ygrid);  // constructs a vector with pixeldata, except number of hits in pixel
    void extract_root_data(void); // extracts data from from root file and saves in a vector "Events"
    void extract_hit_data(vector<vector<vector<double>>> &hitcoord, vector<vector<double>> &pixelgrid, int plane);  // extracts and stores data for each hit in hitcoord, and fills pixelgrid with no. of hits in pixel
    void count_hits(int &count, vector<vector<vector<double>>> hitcoord);  // counts total number of hits in a plane
    void locate_hot_pixels(vector<vector<double>> pixelgrid, vector<int> &hotpixels, int i); // locates hot pixels
    void remove_hot_pixels(vector<vector<vector<double>>> &hitcoord, vector<vector<double>> &pixelgrid, vector<int> hotpixels);  // removes hot pixels
    void align_wo_T(void);  // determines alignment matrix by aligning planes
    void align_w_T(void);  // aligns planes using alignment matrix
    void construct_tracks(double M1M2_slope_lb_x, double M1M2_slope_ub_x, double M1M2_slope_lb_y, double M1M2_slope_ub_y);  // construct M1 -> M6 track
    void pair_tracks(void); // pair electron/positron tracks
    void update_pixelgrids(int plane, vector<vector<double>> pixelgrid);
    void update_hitcoords(int plane, vector<vector<vector<double>>> hitcoord);
    void update_hotpixels(int plane, vector<int> hotpixel);
    void print_pixels(string name); // saves pixeldata
    void print_hotpixels(string name); // prints hotpixels
    void print_hits(string name); // saves hitcoords
    void print_interdistance(string name); // saves distances
    void print_energy(string name);  // saves energy
    void print_slope(string name); // saves angles of a track's incoming + outgoing angle in the Mimosa magnet
    void print_M1M2_slope(string name); // saves M1-M2 angles
    void print_zpos(string name);  // saves z-position of closest approach
    void construct_distarray(void);
    void image_crystal(string name);
    void print_planar_dist(int plane, string name);

  private:

    string DATPATH; // directory to store data
    char const *runno;  // used to name output files
    int ncols, nrows; // pixel columns/rows of Mimosa-26
    double xmin, xmax, ymin, ymax; // maximum x,y coordinates in detector
    char const *filename; // name of alignment file
    double tol; // tolerance for convergence of T matrix
    double M1M2_d_lim;
    double M2M3_d_lim;
    double M6M5_d_lim;
    double Match_d;
    double Match_d_foil;
    double yz_defl_lim;
    vector<vector<vector<vector<vector<double>>>>> paired_tracks;
    vector<vector<vector<double>>> divergence;
    vector<vector<vector<double>>> M1M2_slopes;
    vector<vector<vector<double>>> pixelgrids;
    vector<vector<vector<vector<double>>>> tracks;
    vector<vector<double>> distarray;  // vector with distances between projected and observed hits
    vector<vector<double>> slopes;
    vector<vector<int>> hotpixels;
    vector<double> dr_crit_list; // std vector with descending values of dr for alignment
    vector<double> M1M2_z;
    vector<double> M2M3_z;
    vector<double> M3M4_z;
    vector<double> M6M5_z;
    vector<vector<double>> energies;
    vector<vector<double>> zclosepos;
    vector<double> zplanes;

    static bool sortFunc(const vector<double> &p1, const vector<double> &p2); // sort a vector<vector<>>, descending
    int coord2pixel(double xhit, double yhit);  // converts coordinates to pixelno
    void align_plane(mat &mat_Tot, vector<vector<vector<double>>>  Hits0, vector<vector<vector<double>>>  Hits1, vector<vector<vector<double>>> &Hits2, vector<double> z); // aligns plane '2' using '0' and '1'
    bool check_convergence(mat mat_T, mat mat_temp_T); // checks if T-matrix has changed
    void adjust_coordinates(vector<vector<vector<double>>> &Hits, mat mat_T); // adjusts coordinates for alignment
    void save_hit(vector<double> hitp, vector<double> hit, mat &mat_expected, mat &mat_observed, int indx);  // stores expected hit (projection) and observed hit
    mat construct_T_mat(vector<vector<vector<double>>> Hits0, vector<vector<vector<double>>> Hits1, vector<vector<vector<double>>> Hits2, double dr_crit, vector<double> z);  // uses all hits in 2 planes to project to third, and then stores acceptable observed hits in third planes with corresponding projectd hit
    vector<double> rect_project(vector<double> hit0, vector<double> hit1, vector<double> z);  // updates projected hit 'proj' by rectilinear projection using hit0, hit1 and coordinates z
    double calc_ang(double m, double n); // calculate angle between intersecting lines
    double calc_slope(vec x, vec y); // calculate slope of line
    double calc_dist(double x0, double y0, double x1,double y1); // calculates distance between two points in a plane
    double calc_pair_energy(vector<vector<vector<double>>> pairedtracks);
    vector<double> calc_interdistance(vector<vector<vector<double>>> Hits0, vector<vector<vector<double>>> Hits1, vector<vector<vector<double>>> Hits2, vector<double> z); // calculates the distance from every projected hit into a plane to every observed hit
    mat lines3d_nearestpoints(vec A, vec B, vec C, vec D); // determine the 2 closest points on 2 lines
    void beam_divergence(int, int, int, string);
    void save_vector(string name, vector<double> data);

};

#endif
