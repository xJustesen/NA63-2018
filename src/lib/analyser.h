#ifndef ANALYSER_H
#define ANALYSER_H

#include <math.h>
#include <omp.h>
#include <algorithm>
#include <armadillo>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <limits>
#include "TBranch.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2.h"
#include "TH2F.h"
#include "TLeaf.h"
#include "TMath.h"
#include "TObject.h"
#include "TROOT.h"
#include "TRandom.h"
#include "TSpectrum.h"
#include "TStorage.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"
#include "TVersionCheck.h"
#include "TVirtualFitter.h"
#include "auxillary_functions.h"

using namespace std;
using namespace arma;

class analyser
{
public:
  cube T;            // armadillo cube for collecting alignment matrices (rank-3 tensor)
  int total_events_; // total number of events
  int events_within_cut_;
  vector<vector<vector<double>>> events_container_; // containts root data
  vector<vector<vector<vector<double>>>> hits_container_;

  analyser(vector<double> z, const char *name, string runno, string beamparams);                                            // constructor
  void MakeGrid(vector<vector<double>> &pixelgrid);                                                                         // constructs a vector with pixeldata, except number of hits in pixel
  void ExtractRootData(void);                                                                                               // extracts data from from root file and saves in a vector "Events"
  void ExtractHitData(vector<vector<vector<double>>> &hitcoord, vector<vector<double>> &pixelgrid, int plane);              // extracts and stores data for each hit in hitcoord, and fills pixelgrid with no. of hits in pixel
  void CountHits(int &count, vector<vector<vector<double>>> hitcoord);                                                      // counts total number of hits in a plane
  void LocateHotPixels(vector<vector<double>> pixelgrid, vector<int> &hotpixels);                                           // locates hot pixels
  void RemoveHotPixels(vector<vector<vector<double>>> &hitcoord, vector<vector<double>> &pixelgrid, vector<int> hotpixels); // removes hot pixels
  void AlignWithoutTMatrix(void);                                                                                           // determines alignment matrix by aligning planes
  void AlignWithTMatrix(void);                                                                                              // aligns planes using alignment matrix
  void ConstructTracks(double M1M2_slope_lb_x, double M1M2_slope_ub_x, double M1M2_slope_lb_y, double M1M2_slope_ub_y);     // construct M1 -> M6 track
  void PairTracks(void);                                                                                                    // pair electron/positron tracks
  void UpdatePixelgrids(int plane, vector<vector<double>> pixelgrid);
  void UpdateHitCoordinates(int plane, vector<vector<vector<double>>> hitcoord);
  void UpdateHotPixels(int plane, vector<int> hotpixel);
  void PrintPixels(string name);        // saves pixeldata
  void PrintHotPixels(string name);     // prints hotpixels
  void PrintHits(string name);          // saves hitcoords
  void PrintInterdistance(string name); // saves distances
  void PrintEnergy(string name);        // saves energy
  void PrintSlope(string name);         // saves angles of a track's incoming + outgoing angle in the Mimosa magnet
  void PrintM1M2Slope(string name);     // saves M1-M2 angles
  void FillInterplanarDistanceContainer(void);
  void ImageCrystal(void);
  void PhotonTracksDivergence(void);
  void PrintInterplanarDistance(int plane);
  void FindAxisDeflection(void);
  void FindAxisCounts(void);
  vector<vector<double>> GetEnergies(void);
  int GetEventsWithinCut(void);
  vector<double> BeamDivergencePhotons(vector<double> x, vector<double> y);

private:
  string data_path_;                                                                             // directory to store data
  string run_number_;                                                                            // used to name output files
  int column_count_, row_count_;                                                                 // pixel columns/rows of Mimosa-26
  double detector_xcoord_min_, detector_xcoord_max_, detector_ycoord_min_, detector_ycoord_max_; // maximum x,y coordinates in detector
  const char *file_name_;                                                                        // name of ROOT file
  string beam_parameters_file_name_;
  double T_matrix_convergence_tol_; // tolerance for convergence of T matrix
  double M1_M2_proj_lim_;
  double M2_M3_proj_lim_;
  double M6_M5_proj_lim_;
  double MM_paired_tracks_lim_;
  double foil_paired_tracks_lim_;
  vector<vector<vector<vector<vector<double>>>>> paired_tracks_;
  vector<vector<vector<double>>> M1_M2_slopes_;
  vector<vector<vector<double>>> pixel_grids_;
  vector<vector<vector<vector<double>>>> tracks_;
  vector<vector<double>> interplanar_distance_; // vector with distances between projected and observed hits
  vector<vector<double>> MM_slopes_;
  vector<vector<int>> hot_pixels_;
  vector<double> alignment_radius_lim_; // std vector with descending values of dr for alignment
  vector<double> M1_M2_zcoord_;
  vector<double> M2_M3_zcoord_;
  vector<double> M3_M4_zcoord_;
  vector<double> M6_M5_zcoord_;
  vector<vector<double>> energies;
  vector<double> detector_zcoord_;

  static bool SortVectorDescending(const vector<double> &p1, const vector<double> &p2);                                                                                     // sort a vector<vector<double>>, descending
  int Coord2Pixel(double xhit, double yhit);                                                                                                                                // converts coordinates to pixelno
  void AlignPlane(mat &mat_Tot, vector<vector<vector<double>>> Hits0, vector<vector<vector<double>>> Hits1, vector<vector<vector<double>>> &Hits2, vector<double> z);       // aligns plane '2' using '0' and '1'
  bool CheckMatrixConvergence(mat mat_T, mat mat_temp_T);                                                                                                                   // checks if T-matrix has changed
  void AdjustCoordinates(vector<vector<vector<double>>> &Hits, mat mat_T);                                                                                                  // adjusts coordinates for alignment
  void SaveHit(vector<double> hitp, vector<double> hit, mat &mat_expected, mat &mat_observed, int indx);                                                                    // stores expected hit (projection) and observed hit
  mat ConstructTMatrix(vector<vector<vector<double>>> Hits0, vector<vector<vector<double>>> Hits1, vector<vector<vector<double>>> Hits2, double dr_crit, vector<double> z); // uses all hits in 2 planes to project to third, and then stores acceptable observed hits in third planes with corresponding projectd hit
  vector<double> RectilinearProjection(vector<double> hit0, vector<double> hit1, vector<double> z);                                                                         // updates projected hit 'proj' by rectilinear projection using hit0, hit1 and coordinates z
  double CalculateSlope(vec x, vec y);
  double CalculateSlope(double x0, double y0, double x1, double y1);             // calculate slope of line
  double CalculateEuclideanDistance(double x0, double y0, double x1, double y1); // calculates distance between two points in a plane
  double CalculatePairEnergy(vector<vector<vector<double>>> pairedtracks);
  vector<double> CalculateInterdistance(vector<vector<vector<double>>> Hits0, vector<vector<vector<double>>> Hits1, vector<vector<vector<double>>> Hits2, vector<double> z); // calculates the distance from every projected hit into a plane to every observed hit
  void BeamDivergence(int, int, int);
  void PrintVector(string name, vector<double> data);
  void PrintVector(string name, vector<vector<double>> data);
  void PrintVector(string name, vector<vector<vector<double>>> data);
  vector<double> Linspace(double min, double max, int N);
  void SaveCutData(vector<vector<double>> photons_in_cut, string file_name);
  int IsInsidePolygon(int nvert, vector<double> vertx, vector<double> verty, double testx, double testy);
  void LoadConfigFile(const string &configfile);
};

#endif
