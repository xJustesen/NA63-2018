#ifndef AUXFUNCTIONS_H
#define AUXFUNCTIONS_H

#include <string>
#include <vector>
#include <fstream>
#include <iostream>

void load_doubles(std::string filename, std::vector<double> &data);
void print_vector(std::string filename, std::vector<std::vector<double>> input);

#endif