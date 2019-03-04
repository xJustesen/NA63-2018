#ifndef AUXFUNCTIONS_H
#define AUXFUNCTIONS_H

#include <string>
#include <vector>
#include <fstream>
#include <iostream>

void LoadDoubles(std::string file_name, std::vector<double> &data);
void PrintVector(std::string file_name, std::vector<std::vector<double>> input);

#endif