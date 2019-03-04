#include "auxillary_functions.h"

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
        std::cout << "Unable to open file: " << file_name << "\n";
}

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
}
