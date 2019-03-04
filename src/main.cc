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
    vector<vector<double>> energies;

    /* Processing */
    if (Initializer->config_.Simulation)
    {
        RunSimulation *Simulator = new RunSimulation(Initializer);
        energies = Simulator->GetEnergies();
    }
    else if (!Initializer->config_.Simulation)
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
        RunDataAnalysis *Analyser = new RunDataAnalysis(Initializer, InputDataFile);
        energies = Analyser->GetEnergies();
    }

    /* Post-processing */
    cout << "\nAnalysis finished, saving to file\n";
    PrintVector(Initializer->config_.OutputFilename, energies);
    cout << "\nCompleted\n";
    return 0;
}
