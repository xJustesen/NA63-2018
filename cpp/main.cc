#include "analyser.h"

int main(int argc, char const *argv[]){
  /* List of z-coordinates of planes, used in constructor of analyser class */
  vector<double> z = {0, 1832.3e3, 8913e3, 8989e3, 9196.2e3, 9273.7e3}; // (micro-meters)

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

    /* Report status to terminal */
    cerr << "Plane " << i << ": Hits before removal: " << nhits0 <<  ";\thits after removal: " << nhits1 << ";\thits removed: " << nhits0 - nhits1 << ";\tPixels removed: " << hotpixel.size() << ";\tmax hits allowed per pixel: " << floor(4e-4 * Data.Nevents) <<"\n";
  }

  /* Report progress */
  cerr << "\nHot pixels removed, beginning alignment\n\n";

  /* Align planes with T-matrices already determined */
  Data.T.load("/home/christian/Dropbox/speciale/data/Align/alignment_matrix.txt", arma_ascii);
  Data.align_w_T();

  /* Determine T-matrices & align planes */
  // Data.align_wo_T(); // only use this for alignment run

  /* Report progress */
  cerr << "\nAlignment finished, building tracks\n";

  /* Construct & pair particles tracks & determine energy */
  Data.construct_tracks();
  Data.pair_tracks();
  /* Save data to txt-files */
  Data.image_crystal();
  Data.print_energy();
  Data.print_hits();
  Data.print_hotpixels();
  // Data.construct_distarray();
  Data.print_interdistance();
  Data.print_slope();
  Data.print_pixels();
  Data.print_zpos();
  Data.print_M1M2_slope();

  /* Terminate main */
  return 0;
}
