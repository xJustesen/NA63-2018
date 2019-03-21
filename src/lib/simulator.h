#ifndef SIMULATOR_H
#define SIMULATOR_H

// EXT. LIBRARYS
#include <assert.h>
#include <math.h>
#include <omp.h>
#include <algorithm>
#include <armadillo>
#include <chrono>
#include <fstream>
#include <functional>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>

#include "auxillary_functions.h"

using namespace std;
using namespace arma;

struct DETECTOR
{
    vector<double> x;
    vector<double> y;
    double z;
    double resolution;
    double accuracy;
    double fake_hit_rate;
    int number;
};

class simulator
{
  public:
    double beam_energy_;                             // energy of particle
    double crystal_thicknes_;                        // size of crystal
    vector<vector<vector<vector<double>>>> mimosas_; // equivalent to the "hitcoords" vector in "analyser" class

    simulator(int N, vector<double> z, string run, double BeamEnergy, double CrystalThickness, string filename, string TheorySpec, int bg);
    void GenerateSyntheticData(void); // propagate particles through the experiment
    void PrintHits(void);
    void LinearInterpolation(vector<double> x_coordinate, vector<double> y_coordinate, vector<double> xi, vector<double> &yi);

  private:
    vector<vector<vector<double>>> photons_;   // stores photon data
    vector<vector<vector<double>>> particles_; // stores particles data
    vector<vector<int>> hot_pixels_;
    vector<double> detector_z_coordinates_; // z-coordinates of detectors
    vector<double> intensity_sum_;
    vector<double> intensity_sum_interp_;
    vector<double> emitted_energies_interp_; // interpolated values of "emitted_energies_"
    vector<double> emitted_energies_;        // summed distribution of simulated emitted energies for aligned crystal
    vector<double> x_angle_weight_;
    vector<double> y_angle_weight_;
    vector<double> x_coordinate_weight_;
    vector<double> y_coordinate_weight_;
    vector<double> x_coordinate_;
    vector<double> y_coordinate_;
    vector<double> angles_;
    vector<string> LegalInput; // list of accepted KEYs in config file
    cube alignment_matrix;
    double mean_entry_angle_x_; // mean entry angle of incoming particles
    double dev_entry_angle_x_;  // standard deviation of entry angles of incoming particles
    double mean_entry_angle_y_; // mean entry angle of incoming particles
    double dev_entry_angle_y_;  // standard deviation of entry angles of incoming particles
    double charge_;             // charge of positron in Coulombergy)
    double c_;                  // speed of light in vac. in m/s
    double electron_mass_;      // mass of electron/positron in kg
    double foil_thickness_;     // thickness of converter foil
    double X0_Si_amorph_;       // 9.370E+04,
    double X0_C_gem_;
    double X0_mimosa_;
    double X0_He_;
    double X0_air_;
    double X0_Ta_;
    double X0_tape_;
    double X0_mylar_;
    double intensity_integral_;
    double mimosa_resolution_;
    int include_background_radiation_;
    int photons_on_foil_;
    int total_photon_conversions_;
    int total_events_; // number of simulated events
    int photons_from_crystal_;
    string data_path_; // directory to store data
    string angle_spectrum_;
    string initial_spectrum_;
    string beam_spatial_distribution_;
    string crystal_type_; // name of crystal
    mt19937_64 global_generator_;
    normal_distribution<double> normal_distribution_;
    uniform_real_distribution<double> uniform_real_distribution_;
    DETECTOR M1;
    DETECTOR M2;
    DETECTOR M3;
    DETECTOR M4;
    DETECTOR M5;
    DETECTOR M6;

    void LoadBeamParameters(void);                                                      // generates beam profile
    void MultipleScattering(double X0, double z, mt19937_64, vector<vector<double>> &); // projection taking multiple scattering into account
    void MimosaMagnet(vector<vector<double>> &local_particles);                         // used to caluclate deflection in mimosa magnet
    void ConverterFoil(double X0, mt19937_64, int N_slices, vector<vector<double>> &photons, vector<vector<double>> &particles);
    void MimosaDetector(DETECTOR detector, int eventno, int &, mt19937_64, vector<vector<double>>); // adds a MIMOSA detector, ie simulates a detection
    void ElectronicEnergyDistribution(double &r);                                                   // analytical solution to the equation CDF(x) = r, where CDF is the cumulative distribution function for the fractional energy distribution of a produced electron/positron pair
    void MbplMagnet(vector<vector<double>> &);                                                      // clears particle-vector
    void AmorphMaterial(int &emitted, double X0, double z, double L, mt19937_64, vector<vector<double>> &, vector<vector<double>> &, int no_slices = 1);
    static double PhotonicEnergyDistribution(vec x, double randno, double E, double norm);
    void BorsellinoOpeningAngle(double E1, double E2, double E_phot, double &phi1, double &phi2); // the approximated Borsellino opening angle of e-/e+ pair
    void AlignedCrystal(int &emitted, int no_slices, mt19937_64, vector<vector<double>> &local_photons, vector<vector<double>> &local_particles);
    void AddPhotons(int &emitted, double X0, double d, mt19937_64, vector<vector<double>> &, vector<vector<double>>);
    int BinarySearch(vector<int> numbers, int low, int high, int val);
    int BinarySearch(vector<double> numbers, int low, int high, double value);
    int Coord2Pixel(double xhit, double yhit);
    void LoadHotpixels(void);
    void LoadDoubles(string filename, vector<double> &data);
    void LoadDoubles(string filename, vector<double> &data0, vector<double> &data1);
    void LoadDoubles(string filename, vector<double> &data0, vector<double> &data1, vector<double> &data2, vector<double> &data3);
    void LoadInt(string filename, vector<int> &data);
    void SimplexUpdate(vector<vec> simplex, vec fs, vec &centroid, int &ihigh, int &ilow);
    double SimplexSize(vector<vec> simplex);
    void SimplexReduce(vector<vec> &simplex, int ilow);
    vec SimplexExpand(vec highest, vec centroid);
    vec SimplexReflect(vec highest, vec centroid);
    vec SimplexContract(vec highest, vec centroid);
    int SelectMember(vector<double> weigths); // makes a random draw from a weighted list of numbers
    vector<double> Linspace(double min, double max, int N);
    int IsInsidePolygon(int nvert, vector<double> vertx, vector<double> verty, double testx, double testy); // check if point (testx, testy) is inside boundaries of polygon defined by verticecs (vertx, verty)
    double CalculateDistance(double x0, double y0, double x1, double y1);
    void ProjectPhotons(vector<vector<double>> &particletype, double zcoord); // rectilinear-projection hit into plane
    void SaveVector(string name, vector<double> data);
    void SaveVector(string name, vector<vector<double>> data);
    void SaveVector(string name, vector<vector<vector<double>>> data);
    double TrapezoidalIntegrator(vector<double> x_coordinate, vector<double> y_coordinate);
    void InitializeInputVariables(string filename);
    void InitializeInputVariablesHelper(string Key, string Value);
    int SearchList(vector<string> List, string Key);
    void LoadConfigFile(const string &configfile);
    void AssignDetectorParameters(string config, pt::ptree tree, DETECTOR &detector);

    template <typename lambda>
    void SimplexInitiate(vector<vec> simplex, lambda F, vec &fs)
    {
        for (size_t i = 0; i < simplex.size(); i++)
        {
            fs(i) = F(simplex[i]);
        }
    }

    /* Nelder-Mead simplex algorithm. Tested on Himmelblau's function and Rosenbrock function. In principle able to optimize n-dimensional problems. */
    template <typename lambda>
    vec SimplexNelderMead(lambda F, vector<vec> simplex, double simplex_size_goal)
    {
        vec fs;
        fs.resize(simplex.size());
        vec centroid = zeros<vec>(simplex.size() - 1);
        SimplexInitiate(simplex, F, fs);
        int ilow, ihigh;

        while (SimplexSize(simplex) > simplex_size_goal)
        {
            SimplexUpdate(simplex, fs, centroid, ihigh, ilow); // update simplex with new fs values, this updates centroid, ihigh, ilow.
            vec r = SimplexReflect(simplex[ihigh], centroid);  // reflection
            double fr = F(r);

            /* try expansion */
            if (fr < fs(ilow))
            {
                vec e = SimplexExpand(simplex[ihigh], centroid);
                double fe = F(e);

                if (fe < fr)
                { // accept expansion

                    simplex[ihigh] = e;
                    fs(ihigh) = fe;
                }

                else
                { // reject expansion and accept reflection

                    simplex[ihigh] = r;
                    fs(ihigh) = fr;
                }
            }

            else
            {
                /*if reflection is too poor, reflect in other direction */
                if (fr < fs(ihigh))
                { // accept reflection

                    simplex[ihigh] = r;
                    fs(ihigh) = fr;
                }

                else
                { // reject reflection, try contraction

                    vec c = SimplexContract(simplex[ihigh], centroid);
                    double fc = F(c);

                    if (fc < fs(ihigh))
                    { // accept contraction

                        simplex[ihigh] = c;
                        fs(ihigh) = fc;
                    }
                    else
                    { // reject contraction and reduce

                        SimplexReduce(simplex, ilow);
                        SimplexInitiate(simplex, F, fs);
                    }
                }
            }

        } // end while

        return simplex[ilow];
    }
};

#endif
