#include "main.h"

int main(int argc, char const *argv[])
{
    /* Check if user has provided the proper input */
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
        cerr << "Only legal input for main is a config.txt file and run no.\nTerminating program\n";
        return 0;
    }
    /* Load configuration file and assign values to config struct */
    string config_file = argv[1];
    PreProcessor *Initializer = new PreProcessor(config_file);
    const char *InputDataFile = "";
    /* Check that user has provided ROOT-file if Simulation flag is 0 */
    if (!Initializer->config_.Simulation)
    {
        if (argc < 3)
        {
            cerr << "\nERROR: Simulation flag is set to 0 in configuration file but no ROOT-datafile is given as input\n";
            cerr << "Input given:\n";
            for (int i = 1; i < argc; i++)
            {
                cerr << "\t" << argv[i] << "\n";
            }
            cerr << "Terminating program\n";
            return 0;
        }
        cout << "InputDataFile: " << argv[2] << "\n";
        InputDataFile = argv[2];
    }
    /* Processing */
    vector<vector<double>> energies;
    int number_of_events;
    if (Initializer->config_.Simulation)
    {
        RunSimulation *Simulator = new RunSimulation(Initializer);
        energies = Simulator->GetEnergies();
        number_of_events = Simulator->GetEvents();
    }
    else if (!Initializer->config_.Simulation)
    {
        RunDataAnalysis *Analyser = new RunDataAnalysis(Initializer, InputDataFile);
        energies = Analyser->GetEnergies();
        number_of_events = Analyser->GetEvents();
    }
    /* Post-processing */
    cout << "\nAnalysis finished, saving to file\n";
    PrintVector(Initializer->config_.OutputDataFilename, energies);
    PrintIntSum(Initializer->config_.OutputEventFilename, number_of_events);
    cout << "\nCompleted\n";
    return 0;
}
