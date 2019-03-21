#ifndef AUXFUNCTIONS_H
#define AUXFUNCTIONS_H

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/foreach.hpp>

namespace pt = boost::property_tree;

void LoadDoubles(std::string file_name, std::vector<double> &data);
void PrintVector(std::string file_name, std::vector<std::vector<double>> input);
void PrintIntSum(std::string file_name, int input_number);

template <typename T>
std::vector<T> as_vector(pt::ptree const &tree, pt::ptree::key_type const &key)
{
    std::vector<T> r;
    for (auto &item : tree.get_child(key))
        r.push_back(item.second.get_value<T>());
    return r;
}

#endif