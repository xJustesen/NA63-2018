#include "run_data_analysis.h"

RunDataAnalysis::RunDataAnalysis(PreProcessor *initialize, const char *input_data_file)
{
    analyser *Analyse = new analyser(initialize->z_,
                                     input_data_file,
                                     initialize->config_.runno,
                                     initialize->config_.BeamProfile);
    cout << "\n \t\t Inconsequential ROOT warnings: \n\n";
    Analyse->ExtractRootData();
    cout << "Root data extracted, remvoing hot pixels\n\n";
    /* Analyse data for each plane */
#pragma omp parallel for
    for (size_t i = 0; i < 6; i++)
    {
        int hits_count_w_hot_pixels, hits_count_wo_hot_pixels;
        vector<vector<vector<double>>> hit_coordinates;
        vector<vector<double>> pixel_grid;
        vector<int> hot_pixel;
        Analyse->MakeGrid(pixel_grid);
        Analyse->ExtractHitData(hit_coordinates, pixel_grid, i);
        Analyse->CountHits(hits_count_w_hot_pixels, hit_coordinates);
        Analyse->LocateHotPixels(pixel_grid, hot_pixel, hits_count_w_hot_pixels);
        Analyse->RemoveHotPixels(hit_coordinates, pixel_grid, hot_pixel);
        Analyse->CountHits(hits_count_wo_hot_pixels, hit_coordinates);
        Analyse->UpdateHitCoordinates(i, hit_coordinates);
        Analyse->UpdatePixelgrids(i, pixel_grid);
        Analyse->UpdateHotPixels(i, hot_pixel);
        cout << "Plane " << i << ": Hits before removal: " << hits_count_w_hot_pixels << ";\thits after removal: " << hits_count_wo_hot_pixels << ";\thits removed: " << hits_count_w_hot_pixels - hits_count_wo_hot_pixels << ";\tPixels removed: " << hot_pixel.size() << ";\tmax hits allowed per pixel: " << floor(4e-4 * Analyse->total_events_) << "\n";
    }
    cout << "\nHot pixels removed, beginning alignment\n\n";
    /* Align planes with T-matrices already determined */
    Analyse->T.load("/home/christian/Documents/cern2018/alignment_matrix.txt", arma_ascii); // 2018
    Analyse->AlignWithTMatrix();
    // /* Build tracks and determine energy */
    // cout << "\nAlignment finished, building tracks\n\n";
    // Analyse->ConstructTracks(initialize->config_.cut_lb_x,
    //                            initialize->config_.cut_ub_x,
    //                            initialize->config_.cut_lb_y,
    //                            initialize->config_.cut_ub_y);
    // Analyse->hits_container_.clear();         // clears elements from vector (size = 0, capacity unchanged)
    // Analyse->hits_container_.shrink_to_fit(); // frees memory, capacity = 0
    // cout << "Pairing tracks and calculating energies\n";
    // Analyse->PairTracks();
    // energies_ = Analyse->GetEnergies();
    cout << "\nDetermining axis location using number of photons\n";
    Analyse->FindAxisCounts();
    // cout << "\nDetermining axis location using mean deflection of electrons\n\n";
    // Analyse->FindAxisDeflection();
}

vector<vector<double>> RunDataAnalysis::GetEnergies()
{
    return energies_;
}