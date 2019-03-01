#include "auxfunctions.h"

void load_doubles(std::string filename, std::vector<double> &data) {
    std::ifstream datafile(filename);
    double val;

    if (datafile.is_open()) {
        while (true) {
            datafile >> val;

            if (datafile.eof()) break;

            data.push_back(val);
        }

        datafile.close();

    } else
        std::cout << "Unable to open file: " << filename << "\n";
}

void print_vector(std::string filename, std::vector<std::vector<double>> input)
{
    std::ofstream output;
    output.open(filename, std::ios_base::app);

    for (size_t i = 0; i < input.size(); i++)
    {
        for (size_t j = 0; j < input[i].size(); j++)
        {
            output << input[i][j] << "\n";
        }
    }
}
