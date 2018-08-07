#include "simulator.h"
#include "analyser.h"
#include "globalvars.h"

int main(int argc, char const *argv[]){
  /* Define parameters for constructor */
  int N = 1E+06; // number of events to simulate
  vector<double> z = {0, 1832.3E+03, 8913E+03, 8989E+03, 9196.2E+03, 9273.7E+03}; // z-coordinates for planes (micro-meters)
  char const *name = "";
  char const *sim_type = argv[1];
  char const *runno = argv[2];
  char const *filename = argv[3];

  /* Load beam-parameters from file */
  // params = {Ebeam, d_crystal, mean_entry_angle_x, mean_entry_angle_y, dev_entry_angle_x, dev_entry_angle_y};
  vector<double> params;
  double param;
  string file = "../beamParameters/" + (string)filename;
  ifstream paramsfile (file);
  if (paramsfile.is_open()){
    while (true){
      paramsfile >> param;
      if (paramsfile.eof()) break;
      params.push_back(param);
    }
    paramsfile.close();
  }
  else cout << "Unable to open file containing beam parameters";


  /* Initialise analyse and simulate classs */
  analyser analyse(z, name, runno);
  simulator simulate(N, z, sim_type, params, filename);

  /* Simulate experiment */
  simulate.propagate_particles();
  cerr << "\nData simulated, starting analysis\n";

  /* Analyse simulated data */
  analyse.hitcoords = simulate.mimosas;
  analyse.Nevents = N;
  analyse.construct_tracks();
  analyse.pair_tracks();
  // analyse.construct_distarray();
  cerr << "\nFinished building + pairing tracks, saving energies\n";

  /* Save data */
  analyse.print_energy();
  analyse.print_hits();
  analyse.print_energy();
  analyse.print_zpos();
  analyse.print_slope();
  analyse.print_interdistance();
  analyse.print_M1M2_slope();
  analyse.image_crystal();
  simulate.print_energy();
  simulate.print_hits();
  return 0;
}
