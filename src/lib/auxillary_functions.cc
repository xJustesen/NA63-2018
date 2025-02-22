#include "auxillary_functions.h"

// Load list of doubles
// Input:
//          file_name   Name (and path + extension) of file from which to load data
//          data        vector where data is saved
void LoadDoubles(std::string file_name, std::vector<double> &data)
{
    std::ifstream data_file(file_name);
    double val;
    if (data_file.is_open())
    {
        while (true)
        {
            data_file >> val;

            if (data_file.eof())
                break;

            data.push_back(val);
        }
        data_file.close();
    }
    else
    {
        std::cout << "Unable to open file: " << file_name << "\n";
    }
}

// Print list of doubles
// Input:
//          file_name   Name (and path + extension) of file from which to load data
//          input       vector with data
void PrintVector(std::string file_name, std::vector<std::vector<double>> input)
{
    std::ofstream output;
    output.open(file_name, std::ios_base::app);
    for (size_t i = 0; i < input.size(); i++)
    {
        for (size_t j = 0; j < input[i].size(); j++)
        {
            output << input[i][j] << "\n";
        }
    }
    output.close();
}

// Print list of ints, if file exists then sum all entries in file, add new value and overwrite file
// Input:
//          file_name       Name (and path + extension) of file from which to load data
//          input_number    Inputted number
void PrintIntSum(std::string file_name, int input_number)
{
    std::ifstream input(file_name);
    int previous_number = 0;
    if (input.is_open())
    {
        while (true)
        {
            input >> previous_number;
            if (input.eof())
            {
                break;
            }
        }
        input.close();
    }
    else
    {
        std::cerr << "Unable to open file containing no. of events\n";
    }
    std::ofstream output(file_name);
    output << input_number + previous_number;
    output.close();
}