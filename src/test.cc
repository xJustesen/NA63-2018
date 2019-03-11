#include "./lib/analyser.h"
#include "./lib/simulator.h"
#include "./lib/auxillary_functions.h"
#include "./lib/pre_processor.h"
#include "./lib/run_simulation.h"
#include "./lib/run_data_analysis.h"

int main(int argc, char const *argv[])
{
    /* Check user input */
    if (argc == 1)
    {
        cerr << "No input given main, remember to input config file\n";
        return 0;
    }
    else if (argc > 3)
    {
        cerr << "Illegal input given main:\n";
        for (int i = 1; i < argc; i++)
        {
            cerr << argv[i] << "\n";
        }
        cerr << "Only legal input for main is a config.txt file and run no.\n";
        return 0;
    }
    string config_file = argv[1];
    /* Pre-processing */
    PreProcessor *Initializer = new PreProcessor(config_file);
    simulator *Simulate = new simulator(Initializer->config_.NEvents,
                                         Initializer->z_,
                                         Initializer->config_.CrystalType,
                                         Initializer->config_.BeamEnergy,
                                         Initializer->config_.CrystalThickness,
                                         Initializer->config_.BeamProfile,
                                         Initializer->config_.TheorySpectrum,
                                         Initializer->config_.IncludeBG);
    vector<double> x, y, xi, yi;
    x = y = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    xi = {1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0};
    Simulate->LinearInterpolation(x, y, xi, yi);
    return 0;
}
