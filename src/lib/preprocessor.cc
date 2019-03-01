#include "preprocessor.h"

preprocessor::preprocessor(string ConfigFile)
{
    /* Define parameters for constructor */
    LegalInput = {
        "BeamEnergy",
        "CrystalThickness",
        "BeamProfile",
        "cut_lb_x",
        "cut_ub_x",
        "cut_lb_y",
        "cut_ub_y",
        "OutputFilename",
        "CrystalType",
        "runno",
        "DataPath",
        "TheorySepctrum",
        "IncludeBG",
        "NEvents",
        "Simulation"};

    /* Assign default values to angular cuts */
    Config.cut_lb_x = Config.cut_lb_y = (-1.0) * 1E+17;
    Config.cut_ub_x = Config.cut_ub_y = 1E+17;

    /* Initialize simulation input parameters */
    cout << "\nConfig-struct contains:\n";
    try
    {
        InitializeInputVariables(ConfigFile);
    }
    catch (const char *msg)
    {
        cerr << msg << endl;
    };
    z = {0, 1832.3E+03, 8913E+03, 8989E+03, 9196.2E+03, 9273.7E+03};
    cout << "Entry angles x : {" << Config.cut_lb_x << ", " << Config.cut_ub_x << "} rad\n";
    cout << "Entry angles y : {" << Config.cut_lb_y << ", " << Config.cut_ub_y << "} rad\n";
    cout << "z : {" << z[0] << ", " << z[1] << ", " << z[2] << ", " << z[3] << ", " << z[4] << ", " << z[5] << "}\n";
}

void preprocessor::InitializeInputVariables(std::string filename)
{
    std::ifstream ConfigFile(filename);
    std::string Key;
    std::string Value;

    bool EmptyFile = true;
    if (ConfigFile.is_open())
    {
        while (true)
        {
            ConfigFile >> Key >> Value;
            if (ConfigFile.eof())
            {
                if (EmptyFile)
                {
                    throw "Config file appears to be empty!\n";
                }
                break;
            }
            EmptyFile = false;
            InitializeInputVariablesHelper(Key, Value);
        }
    }
    else
    {
        throw "Unable to open config file!\n";
    }
}

void preprocessor::InitializeInputVariablesHelper(std::string Key, std::string Value)
{
    int Index = SearchList(LegalInput, Key);
    switch (Index)
    {
    case 0:
        Config.BeamEnergy = std::stod(Value);
        std::cout << "BeamEnergy : " << Config.BeamEnergy << "\n";
        break;
    case 1:
        Config.CrystalThickness = std::stod(Value);
        std::cout << "CrystalThickness : " << Config.CrystalThickness << "\n";
        break;
    case 2:
        Config.BeamProfile = Value;
        std::cout << "BeamProfile : " << Config.BeamProfile << "\n";
        break;
    case 3:
        Config.cut_lb_x = std::stod(Value);
        std::cout << "cut_lb_x : " << Config.cut_lb_x << "\n";
        break;
    case 4:
        Config.cut_ub_x = std::stod(Value);
        std::cout << "cut_ub_x : " << Config.cut_ub_x << "\n";
        break;
    case 5:
        Config.cut_lb_y = std::stod(Value);
        std::cout << "cut_lb_y : " << Config.cut_lb_y << "\n";
        break;
    case 6:
        Config.cut_ub_y = std::stod(Value);
        std::cout << "cut_ub_y : " << Config.cut_ub_y << "\n";
        break;
    case 7:
        Config.OutputFilename = Value;
        std::cout << "OutputFilename : " << Config.OutputFilename << "\n";
        break;
    case 8:
        Config.CrystalType = Value;
        std::cout << "CrystalType : " << Config.CrystalType << "\n";
        break;
    case 9:
        Config.runno = Value;
        std::cout << "runno : " << Config.runno << "\n";
        break;
    case 10:
        Config.DataPath = Value;
        std::cout << "DataPath : " << Config.DataPath << "\n";
        break;
    case 11:
        Config.TheorySpectrum = Value;
        std::cout << "TheorySpectrum : " << Config.TheorySpectrum << "\n";
        break;
    case 12:
        Config.IncludeBG = stoi(Value);
        std::cout << "IncludeBG : " << Config.IncludeBG << "\n";
        break;
    case 13:
        Config.NEvents = stoi(Value);
        std::cout << "NEvents : " << Config.NEvents << "\n";
        break;
    case 14:
        Config.Simulation = stoi(Value);
        std::cout << "Simulation : " << Config.Simulation << "\n";
        break;
    case -1:
        std::cerr << "\nIllegal variable name:\t'" << Key << "'\tin config file\n";
        std::cerr << "Legal variable names are:\n";
        for (std::string Key : LegalInput)
        {
            std::cerr << Key << "\n";
        }
        std::cerr << std::endl;
        break;
    }
}

int preprocessor::SearchList(std::vector<std::string> List, std::string Key)
{
    for (size_t i = 0; i < List.size(); i++)
    {
        if (Key == List[i])
        {
            return i;
        }
    }
    return -1;
}