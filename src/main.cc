#include "./lib/analyser.h"
#include "./lib/simulator.h"
#include "./lib/auxfunctions.h"
#include "./lib/preprocessor.h"
#include "./lib/RunSimulation.h"
#include "./lib/RunDataAnalysis.h"

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
    string ConfigFile = argv[1];

    /* Pre-processing */
    preprocessor *initialize = new preprocessor(ConfigFile);
    vector<vector<double>> energies;

    /* Processing */
    if (initialize->Config.Simulation)
    {
        RunSimulation *simulation = new RunSimulation(initialize);
        energies = simulation->get_energies();
    }
    else if (!initialize->Config.Simulation)
    {
        if (argc < 3)
        {
            cerr << "\nMissing main input argument for data analysis.\nExpected inputs are:\n\tconfig.txt\n\tdata.root\n";
            cerr << "Input given:\n";
            for (int i = 1; i < argc; i++)
            {
                cerr << "\t" << argv[i] << "\n";
            }
            return 0;
        }
        cout << "InputDataFile: " << argv[2] << "\n";
        const char *InputDataFile = argv[2];
        RunDataAnalysis *analyser = new RunDataAnalysis(initialize, InputDataFile);
        energies = analyser->get_energies();
    }

    /* Post-processing */
    cout << "\nAnalysis finished, saving to file\n";
    print_vector(initialize->Config.OutputFilename, energies);
    cout << "\nCompleted\n";
    return 0;
}
