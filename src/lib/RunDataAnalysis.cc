#include "RunDataAnalysis.h"

RunDataAnalysis::RunDataAnalysis(preprocessor *initialize, const char *InputDataFile)
{
    analyser *_analyse = new analyser(initialize->z,
                                      InputDataFile,
                                      initialize->Config.runno,
                                      initialize->Config.BeamProfile);
    /* Extract data from root file */
    cout << "\n \t\t Inconsequential ROOT warnings: \n\n";
    _analyse->extract_root_data();
    cout << "Root data extracted, remvoing hot pixels\n\n";

    /* Analyse data for each plane */
#pragma omp parallel for
    for (size_t i = 0; i < 6; i++)
    {
        int nhits0, nhits1;
        vector<vector<vector<double>>> hitcoord;
        vector<vector<double>> pixelgrid;
        vector<int> hotpixel;
        vector<double> xgrid, ygrid;

        /* Make (x, y) coordinates for grid of pixels */
        _analyse->make_grid(pixelgrid, xgrid, ygrid);

        /* Construct and fill pixelgrids & extract coordinates of hit */
        _analyse->extract_hit_data(hitcoord, pixelgrid, i);
        _analyse->count_hits(nhits0, hitcoord);

        /* Identify and remove hot pixels */
        _analyse->locate_hot_pixels(pixelgrid, hotpixel, nhits0);
        _analyse->remove_hot_pixels(hitcoord, pixelgrid, hotpixel);
        _analyse->count_hits(nhits1, hitcoord);

        /* Update variables */
        _analyse->update_hitcoords(i, hitcoord);
        _analyse->update_pixelgrids(i, pixelgrid);
        _analyse->update_hotpixels(i, hotpixel);

        cout << "Plane " << i << ": Hits before removal: " << nhits0 << ";\thits after removal: " << nhits1 << ";\thits removed: " << nhits0 - nhits1 << ";\tPixels removed: " << hotpixel.size() << ";\tmax hits allowed per pixel: " << floor(4e-4 * _analyse->Nevents) << "\n";
    }
    cout << "\nHot pixels removed, beginning alignment\n\n";

    /* Align planes with T-matrices already determined */
    _analyse->T.load("/home/christian/Documents/cern2018/alignment_matrix.txt", arma_ascii); // 2018
    _analyse->align_w_T();

    /* Build tracks and determine energy */
    cout << "\nAlignment finished, building tracks\n\n";
    _analyse->construct_tracks(initialize->Config.cut_lb_x,
                               initialize->Config.cut_ub_x,
                               initialize->Config.cut_lb_y,
                               initialize->Config.cut_ub_y);
    _analyse->hitcoords.clear();         // clears elements from vector (size = 0, capacity unchanged)
    _analyse->hitcoords.shrink_to_fit(); // frees memory, capacity = 0
    cout << "Pairing tracks and calculating energies\n";
    _analyse->pair_tracks();
    energies = _analyse->get_energies();
}

vector<vector<double>> RunDataAnalysis::get_energies()
{
    return energies;
}