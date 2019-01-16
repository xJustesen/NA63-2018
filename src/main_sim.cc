#include "simulator.h"
#include "analyser.h"
#include "globalvars.h"

int main(int argc, char const *argv[]) {

        /* Define parameters for constructor */
        int N = 1E+06; // number of events to simulate
        vector<double> z = {0, 1832.3E+03, 8913E+03, 8989E+03, 9196.2E+03, 9273.7E+03}; // z-coordinates for planes (micro-meters) 2018
        int bg = 1; // include background?

        char const *name = "";
        char const *sim_type = argv[1];
        char const *runno = argv[2];
        char const *filename = argv[3];
        char const *spectrum1;

        double cut;
        vector<double> cuts;

        cerr << "\nArguments : " << argc << "\n";
        for (int i = 0; i < argc; i++) cerr << argv[i] << "\n";
        cerr << "\n";

        if (argc < 4) {

                spectrum1 = "";

        } else {

                spectrum1 = argv[4];

        }

        if (argc > 6) {

                ifstream filename (argv[6]);

                if (filename.is_open()) {

                        while (true) {

                                filename >> cut;

                                if (filename.eof()) break;

                                cuts.push_back(cut);

                        }

                        filename.close();

                } else cout << "Unable to open file containing angles for data cuts";

        } else cuts = {-INFINITY, INFINITY, -INFINITY, INFINITY};

        /* Load beam-parameters from file */
        vector<double> params;
        double param;
        string file = "../beamParameters/" + (string)filename;
        ifstream paramsfile (file);

        if (paramsfile.is_open()) {

                while (true) {

                        paramsfile >> param;

                        if (paramsfile.eof()) break;

                        params.push_back(param);

                }

                paramsfile.close();

        } else cout << "Unable to open file containing beam parameters";

        /* Set limits for data cut */
        double cut_lb_x = cuts[0];
        double cut_ub_x = cuts[1];
        double cut_lb_y = cuts[2];
        double cut_ub_y = cuts[3];

        string Name;
        if (argc > 6) {

                Name = (string)argv[2] + "_cuts_halfSumPeak";

        } else Name = argv[2];

        cerr << "\nEntry angles x: {" << cut_lb_x << ", " << cut_ub_x << "} rad";
        cerr << "\nEntry angles y: {" << cut_lb_y << ", " << cut_ub_y << "} rad";


        /* Initialise simulate classs */
        simulator* simulate = new simulator(N, z, sim_type, params, filename, spectrum1, bg);

        /* Simulate experiment */
        simulate->propagate_particles();
        cerr << "\nData simulated, starting analysis\n";

        /* Analyse simulated data */
        // simulate->print_energy(Name);
        analyser* analyse = new analyser(z, name, runno, filename);
        analyse->hitcoords = simulate->mimosas;
        delete simulate; // delete simulate instace since no longer needed
        analyse->Nevents = N;
        analyse->construct_tracks(cut_lb_x, cut_ub_x, cut_lb_y, cut_ub_y);
        analyse->hitcoords.clear();   // clears elements from vector (size = 0, capacity unchanged)
        analyse->hitcoords.shrink_to_fit(); // frees memory, capacity = 0
        analyse->pair_tracks();
        // analyse->construct_distarray();
        /* Save data to txt-files */
        // for (int i = 0; i < 6; i++) {
        //
        //   analyse->print_planar_dist(i, Name);
        //
        // }

        // analyse->image_crystal(Name);
        // analyse->print_energy(Name);
        // analyse->print_hits(Name);
        // analyse->print_hotpixels(Name);
        // analyse->print_interdistance(Name);
        // analyse->print_slope(Name);
        // analyse->print_pixels(Name);
        // analyse->print_zpos(Name);
        // analyse->print_M1M2_slope(Name);

        return 0;

}
