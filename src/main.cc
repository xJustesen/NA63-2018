#include "analyser.h"

void load_doubles(string filename, vector<double> &data) {
  ifstream datafile (filename);
  double val;
  if (datafile.is_open()){
    while (true){
      datafile >> val;
      if (datafile.eof()) break;
      data.push_back(val);
    }
    datafile.close();
  }
  else cout << "Unable to open file: " << filename << "\n";
}


int main(int argc, char const *argv[]){
  vector<double> z = {0, 1832.3e3, 8913e3, 8989e3, 9196.2e3, 9273.7e3}; // (micro-meters)
  vector<double> cuts;
  cerr << argc << "\n";
  if (argc == 4) {
    string filename = argv[3];
    load_doubles(argv[3], cuts);
  } else {
    cuts = {-INFINITY, INFINITY};
  }
  cerr << cuts[0] << "\t" << cuts[1] << "\n";
  /* Initialise class "analyser" */
  analyser Data(z, argv[1], argv[2]);

  /* Extract data from root file */
  cerr << "\n \t\t Inconsequential ROOT warnings: \n\n";
  Data.extract_root_data();
  cerr << "Root data extracted, remvoing hot pixels\n\n";
  /* Analyse data for each plane */
  #pragma omp parallel for
  for (size_t i = 0; i < 6; i++){
    int nhits0, nhits1;
    vector<vector<vector<double>>> hitcoord;
    vector<vector<double>> pixelgrid;
    vector<int> hotpixel;
    vector<double> xgrid, ygrid;

    /* Make (x, y) coordinates for grid of pixels */
    Data.make_grid(pixelgrid, xgrid, ygrid);

    /* Construct and fill pixelgrids & extract coordinates of hit */
    Data.extract_hit_data(hitcoord, pixelgrid, i);
    Data.count_hits(nhits0, hitcoord);

    /* Identify and remove hot pixels */
    Data.locate_hot_pixels(pixelgrid, hotpixel, nhits0);
    Data.remove_hot_pixels(hitcoord, pixelgrid, hotpixel);
    Data.count_hits(nhits1, hitcoord);

    /* Update variables */
    Data.update_hitcoords(i, hitcoord);
    Data.update_pixelgrids(i, pixelgrid);
    Data.update_hotpixels(i, hotpixel);

    cerr << "Plane " << i << ": Hits before removal: " << nhits0 <<  ";\thits after removal: " << nhits1 << ";\thits removed: " << nhits0 - nhits1 << ";\tPixels removed: " << hotpixel.size() << ";\tmax hits allowed per pixel: " << floor(4e-4 * Data.Nevents) <<"\n";
  }

  cerr << "\nHot pixels removed, beginning alignment\n\n";

  /* Align planes with T-matrices already determined */
  Data.T.load("/home/christian/Dropbox/speciale/data/Align/alignment_matrix.txt", arma_ascii);
  Data.align_w_T();

  /* Determine T-matrices & align planes */
  // Data.align_wo_T(); // only use this for alignment run

  cerr << "\nAlignment finished, building tracks\n";

  /* Construct & pair particles tracks & determine energy */
  for (size_t i = 0; i < cuts.size()/2; i++) {
    string name;
    if (argc == 4) {
      name = (string)argv[2] + "_cut_" + to_string(i);
    } else {
      name = argv[2];
    }
    double cut_ub = cuts[2*i];
    double cut_lb = cuts[2*i + 1];
    Data.construct_tracks(cut_ub, cut_lb);
    Data.pair_tracks();
    /* Save data to txt-files */
    Data.image_crystal(name);
    Data.print_energy(name);
    Data.print_hits(name);
    Data.print_hotpixels(name);
    Data.print_interdistance(name);
    Data.print_slope(name);
    Data.print_pixels(name);
    Data.print_zpos(name);
    Data.print_M1M2_slope(name);
  }


  /* Terminate main */
  return 0;
}
