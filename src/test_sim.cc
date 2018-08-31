#include "simulator.h"
#include "analyser.h"
#include "globalvars.h"

int main(int argc, char const *argv[]) {

  /* Define parameters for constructor */
  int N = 1E+06; // number of events to simulate
  vector<double> z = {0, 1832.3E+03, 8913E+03, 8989E+03, 9196.2E+03, 9273.7E+03}; // z-coordinates for planes (micro-meters)

  char const *name = "";
  char const *sim_type = argv[1];
  char const *runno = argv[2];
  char const *filename = argv[3];

  double cut;
  vector<double> cuts;

  if (argc == 5) {

    ifstream filename (argv[4]);

    if (filename.is_open()) {

      while (true) {

        filename >> cut;

        if (filename.eof()) break;

        cuts.push_back(cut);

      }

      filename.close();

    } else cout << "Unable to open file containing angles for data cuts";

  } else cuts = {-INFINITY, INFINITY};

  /* Load beam-parameters from file */
  vector<double> params;
  double param;
  string file = "../beamParameters/" + (string)filename;
  ifstream paramsfile (file);

  if (paramsfile.is_open()) {

    while (true){

      paramsfile >> param;

      if (paramsfile.eof()) break;

      params.push_back(param);

    }

    paramsfile.close();

  } else cout << "Unable to open file containing beam parameters";

  /* Initialise analyse and simulate classs */
  analyser analyse(z, name, runno);
  simulator simulate(N, z, sim_type, params, filename);


  /* Simulate experiment */
  simulate.propagate_particles();
  cerr << "\nData simulated, starting analysis\n";

  /* Analyse simulated data */
  analyse.hitcoords = simulate.mimosas;
  analyse.Nevents = N;

  for (size_t i = 0; i < cuts.size()/2; i++) {

    double cut_ub = cuts[2*i];
    double cut_lb = cuts[2*i + 1];

    string name;

    if (argc == 5) name = (string)argv[2] + "_cut_" + to_string(cut_ub);
    else name = argv[2];

    cerr << "\nEntry angles : {" << cut_ub << ", " << cut_lb << "} rad";
    analyse.construct_tracks(cut_ub, cut_lb);
    analyse.pair_tracks();
    /* Save data to txt-files */
    analyse.image_crystal(name);
    analyse.print_energy(name);
    analyse.print_hits(name);
    analyse.print_hotpixels(name);
    analyse.print_interdistance(name);
    analyse.print_slope(name);
    analyse.print_pixels(name);
    analyse.print_zpos(name);
    analyse.print_M1M2_slope(name);

    for (int i = 0; i < 6; i++) {

      analyse.print_planar_dist(i, "sim");

    }

    simulate.print_energy(name);

  }

  return 0;

}
