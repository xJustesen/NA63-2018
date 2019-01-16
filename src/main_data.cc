#include "analyser.h"

void load_doubles(string filename, vector<double> &data) {

        ifstream datafile (filename);
        double val;

        if (datafile.is_open()) {

                while (true) {

                        datafile >> val;

                        if (datafile.eof()) break;

                        data.push_back(val);
                }

                datafile.close();

        } else cout << "Unable to open file: " << filename << "\n";

}


int main(int argc, char const *argv[]) {

        vector<double> z = {0, 1832.3e3, 8913e3, 8989e3, 9196.2e3, 9273.7e3}; // (micro-meters) 2018 experiment
        // vector<double> z = {0, 1757E+03, 9277E+03, 9355E+03, 9535E+03, 9613E+03}; // 2017 experiment

        vector<double> cuts;
        string name;

        cerr << "\nmain has " << argc << " arguments :\n";
        for (int i = 0; i < argc; i++) cerr << argv[i] << "\n";
        cerr << "\n";

        if (argc > 4) {

                string filename = argv[3];
                load_doubles(argv[3], cuts);
                name = (string)argv[2] + "_cuts_halfSumPeak";

        } else {

                cuts = {-1E+17, 1E+17, -1E+17, 1E+17};
                name = argv[2];

        }

        double cut_lb_x = cuts[0];
        double cut_ub_x = cuts[1];
        double cut_lb_y = cuts[2];
        double cut_ub_y = cuts[3];

        /* Initialise class "analyser" */
        analyser* Data = new analyser(z, argv[1], argv[2], argv[3]);

        /* Extract data from root file */
        cerr << "\n \t\t Inconsequential ROOT warnings: \n\n";
        Data->extract_root_data();
        cerr << "Root data extracted, remvoing hot pixels\n\n";

        /* Analyse data for each plane */
        #pragma omp parallel for
        for (size_t i = 0; i < 6; i++) {

                int nhits0, nhits1;
                vector<vector<vector<double> > > hitcoord;
                vector<vector<double> > pixelgrid;
                vector<int> hotpixel;
                vector<double> xgrid, ygrid;

                /* Make (x, y) coordinates for grid of pixels */
                Data->make_grid(pixelgrid, xgrid, ygrid);

                /* Construct and fill pixelgrids & extract coordinates of hit */
                Data->extract_hit_data(hitcoord, pixelgrid, i);
                Data->count_hits(nhits0, hitcoord);

                /* Identify and remove hot pixels */
                Data->locate_hot_pixels(pixelgrid, hotpixel, nhits0);
                Data->remove_hot_pixels(hitcoord, pixelgrid, hotpixel);
                Data->count_hits(nhits1, hitcoord);

                /* Update variables */
                Data->update_hitcoords(i, hitcoord);
                Data->update_pixelgrids(i, pixelgrid);
                Data->update_hotpixels(i, hotpixel);

                Data->print_planar_dist(i, name);

                cerr << "Plane " << i << ": Hits before removal: " << nhits0 <<  ";\thits after removal: " << nhits1 << ";\thits removed: " << nhits0 - nhits1 << ";\tPixels removed: " << hotpixel.size() << ";\tmax hits allowed per pixel: " << floor(4e-4 * Data->Nevents) <<"\n";

        }

        cerr << "\nHot pixels removed, beginning alignment\n\n";

        /* Align planes with T-matrices already determined */
        Data->T.load("/home/christian/Documents/cern2018/alignment_matrix.txt", arma_ascii);
        // Data->T.load("/home/christian/Dropbox/Cern2018Experiment/alignment_matrices_1.txt",arma_ascii);
        Data->align_w_T();
        //
        /* Determine T-matrices & align planes */
        // Data->align_wo_T(); // only use this for alignment run

        cerr << "\nAlignment finished, building tracks\n\n";

        /* Construct & pair particles tracks & determine energy */
        cerr << "Entry angles x: {" << cut_lb_x << ", " << cut_ub_x << "} rad\n";
        cerr << "Entry angles y: {" << cut_lb_y << ", " << cut_ub_y << "} rad\n\n";

        // Data->find_axis_alt();

        Data->construct_tracks(cut_lb_x, cut_ub_x, cut_lb_y, cut_ub_y);
        Data->pair_tracks();
        // Data->construct_distarray();
        // /* Save data to txt-files */
        // Data->image_crystal(name);
        Data->print_energy(name);
        // Data->print_hits(name);
        // Data->print_hotpixels(name);
        // Data->print_interdistance(name);
        // Data->print_slope(name);
        // Data->print_pixels(name);
        // Data->print_zpos(name);
        // Data->print_M1M2_slope(name);
        //
        // for (int i = 0; i < 6; i++) {
        //
        //   Data->print_planar_dist(i, "data");
        //
        // }
        //
        return 0;

}
