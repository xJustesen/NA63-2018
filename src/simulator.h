#ifndef SIMULATOR_H
#define SIMULATOR_H

// EXT. LIBRARYS
#include <iostream>
#include <vector>
#include <algorithm>
#include <math.h>
#include <armadillo>
#include <fstream>
#include <omp.h>
#include <random>
#include <functional>
#include <numeric>
#include <assert.h>

using namespace std;
using namespace arma;

class simulator {

  public:

    double E; // energy of particle
    double d_c; // size of crystal
    vector<vector<vector<vector<double>>>> mimosas; // equivalent to the "hitcoords" vector in "analyser" class

    /* Public methods */
    simulator(int N, vector<double> z, char const *run, vector<double> params, char const *filename, char const *s1, char const *s2);
    void propagate_particles(void);  // propaget particles through the experiment
    void print_hits(void);
    void print_energy(string name);

  private:

    vector<vector<vector<double>>> photons; // stores photon data
    vector<vector<vector<double>>> particles; // stores particles data
    vector<vector<double>> intensities;  // intensity at different points in space
    vector<vector<double>> angles;  // angles used to determine emission angle of photon
    vector<vector<double>> photon_angles;
    vector<vector<int>> hotpixels;
    vector<double> energies;  // photon energies used in simulation. Perfect detection should yield exactly these energies
    vector<double> zplanes; // z-coordinates of detectors
    vector<double> intensities_sum; // summed distribution of simulated energy intensities for aligned crystal
    vector<double> intensity_sum;
    vector<double> energies_interp; // interpolated values of "emitted_energies"
    vector<double> emitted_energies;  // summed distribution of simulated emitted energies for aligned crystal
    double mean_entry_angle_x; // mean entry angle of incoming particles
    double dev_entry_angle_x; // standard deviation of entry angles of incoming particles
    double mean_entry_angle_y; // mean entry angle of incoming particles
    double dev_entry_angle_y; // standard deviation of entry angles of incoming particles
    double q; // charge of positron in Coulombergy)
    double c; // speed of light in vac. in m/s
    double m; // mass of electron/positron in kg
    double d_f; // thickness of converter foil
    double X0_Si_amorph; // 9.370E+04,
    double X0_C_amorph;
    double X0_C_gem;
    double X0_Mimosa;
    double X0_He;
    double X0_air;
    double X0_Ta;
    double X0_tape;
    double X0_Mylar;
    double I_integral;
    double mimosa_res;
    int no_photons;
    int conversions;
    int Nevents;  // number of simulated events
    string DATPATH; // directory to store data
    string angle_spec;
    string initial_spec;
    string beam_spatial_distro;
    char const *name;

    /* Physics methods */
    void generate_beam(void);  // generates beam profile
    void SA_mult_scat(int i, double X0, double z); // projection taking multiple scattering into account
    void mimosa_magnet(int eventno);  // used to caluclate deflection in mimosa magnet
    void converter_foil(int eventno, double d_f, double X0, int N_slices = 1);
    void mimosa_detector(int planeno, int eventno, int&); // adds a MIMOSA detector, ie simulates a detection
    void pair_production(int eventno, double d_f, double X0, int N_slices = 1);  // adds a converter foil to calculate pair production
    bool pair_produced(double l, double X0_f); // determines whether a photon creates an electron/positron pair. R is a random double between 0 and 1
    void fractional_electron_energy(double &x); // caluclates the fractional (in terms of photon energy) energy of a pair-produced electron/positron
    void electronic_energy_distribution(double &r); // analytical solution to the equation CDF(x) = r, where CDF is the cumulative distribution function for the fractional energy distribution of a produced electron/positron pair
    void MBPL_magnet(int eventno); // clears particle-vector
    void amorph_material(int eventno, int &emitted, double X0, double z, double L, int no_slices = 1);
    void amorph_crystal(int eventno, int &emitted, double X0, double d_c, int no_slices = 1);  // simulate amorphous crystal
    bool photon_emitted_amorph(double l, double X0, double Epart);  // calculate wether or not photon is emitted (false: no emission)
    static double photonic_energy_distribution(vec x, double randno, double E, double norm);
    void Borsellino(double E1, double E2, double E_phot, double &phi1, double &phi2);  // the approximated Borsellino opening angle of e-/e+ pair
    void make_intensity_distro(vector<vector<double>> &intensity, vector<vector<double>> &angles);
    bool photon_emitted_aligned(double l);
    void aligned_crystal(int eventno, int &emitted, int no_slices);

    /* Numerical methods */
    int binarySearch(vector<int> numbers, int low, int high, int val);
    int coord2pixel(double xhit, double yhit);
    void load_hotpixels(void);
    void load_doubles(string filename, vector<double> &data);
    void load_int(string filename, vector<int> &data);
    void simplex_update(vector<vec> simplex, vec fs, vec &centroid, int &ihigh, int &ilow);
    double simplex_size(vector<vec> simplex);
    void simplex_reduce(vector<vec> &simplex, int ilow);
    vec simplex_expand(vec highest, vec centroid);
    vec simplex_reflect(vec highest, vec centroid);
    vec simplex_contract(vec highest, vec centroid);
    // void simplex_initiate(vector<vec> simplex, function<double(vec)> F, vec &fs);
    // vec simplex_NM(function<double(vec)> F, vector<vec> simplex, double simplex_size_goal = 1e-6);
    double VoseAliasMethod_draw(vector<int> Alias, vector<double> Prob);
    void VoseAliasMethod_table(vector<double> distro, vector<int> &Alias, vector<double> &Prob);
    int select_member(vector<double> numbers, vector<double> weigths); // makes a random draw from a weighted list of numbers
    vector<double> linspace(double min, double max, int N);
    void linterp(vector<double> x, vector<double> y, vector<double> xi, vector<double> &yi);
    int isInside(int nvert, vector<double> vertx, vector<double> verty, double testx, double testy); // check if point (testx, testy) is inside boundaries of polygon defined by verticecs (vertx, verty)
    double calc_dist(double x0, double y0, double x1,double y1);
    void project_particle(vector<vector<vector<double>>> &particletype, double zcoord, int eventno); // rectilinear-projection hit into plane
    void save_vector(string name, vector<vector<double>> data);

    template<typename lambda>
    void simplex_initiate(vector<vec> simplex, lambda F, vec &fs) {

      for (size_t i = 0; i < simplex.size(); i++) {

        fs(i) = F(simplex[i]);

      }

    }

    /* Nelder-Mead simplex algorithm. Tested on Himmelblau's function and Rosenbrock function. In principle able to optimize n-dimensional problems. */
    template<typename lambda>
    vec simplex_NM(lambda F, vector<vec> simplex, double simplex_size_goal) {

      vec fs; fs.resize(simplex.size());
      vec centroid = zeros<vec>(simplex.size() - 1);
      simplex_initiate(simplex, F, fs);
      int ilow, ihigh;

      while (simplex_size(simplex) > simplex_size_goal) {

        simplex_update(simplex, fs, centroid, ihigh, ilow); // update simplex with new fs values, this updates centroid, ihigh, ilow.
        vec r = simplex_reflect(simplex[ihigh], centroid);  // reflection
        double fr = F(r);

        /* try expansion */
        if (fr < fs(ilow)) {

          vec e = simplex_expand(simplex[ihigh], centroid);
          double fe = F(e);

          if (fe < fr ){ // accept expansion

            simplex[ihigh] = e;
            fs(ihigh) = fe;

          }

          else { // reject expansion and accept reflection

            simplex[ihigh] = r;
            fs(ihigh) = fr;

          }

        }

        else {

           /*if reflection is too poor, reflect in other direction */
          if (fr < fs(ihigh)) { // accept reflection

            simplex[ihigh] = r;
            fs(ihigh) = fr;

          }

          else { // reject reflection, try contraction

            vec c = simplex_contract(simplex[ihigh], centroid);
            double fc = F(c);

            if (fc < fs(ihigh)){ // accept contraction

              simplex[ihigh] = c;
              fs(ihigh) = fc;

            }
            else { // reject contraction and reduce

              simplex_reduce(simplex, ilow);
              simplex_initiate(simplex, F, fs);

            }

          }

        }

      } // end while

      return simplex[ilow];

    }

};

#endif
