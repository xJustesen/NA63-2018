#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/foreach.hpp>
#include <string>
#include <set>
#include <exception>
#include <iostream>
#include <vector>

using namespace std;
namespace pt = boost::property_tree;

template <typename T>
std::vector<T> as_vector(pt::ptree const &tree, pt::ptree::key_type const &key)
{
    std::vector<T> r;
    for (auto &item : tree.get_child(key))
        r.push_back(item.second.get_value<T>());
    return r;
}

struct DETECTOR
{
    vector<double> x;
    vector<double> y;
    double z;
    double resoultion;
    double accuracy;
    double fake_hit_rate;
};

struct SCATTERINGLENGTH
{
    double Si;
    double C;
    double He;
    double Air;
    double Ta;
    double Tape;
    double Mylar;
};

struct CONFIG
{
    DETECTOR M1;
    DETECTOR M2;
    DETECTOR M3;
    DETECTOR M4;
    DETECTOR M5;
    DETECTOR M6;
    double c;
    double q;
    double e_m;
    SCATTERINGLENGTH X0;

    void load(const string &configfile);
    void AssignDetectorParameters(string config, pt::ptree tree, DETECTOR &detector);
};

void CONFIG::AssignDetectorParameters(string config, pt::ptree tree, DETECTOR &detector)
{
    for (auto i : as_vector<double>(tree, config + ".x"))
    {
        detector.x.push_back(i);
    }
    for (auto i : as_vector<double>(tree, config + ".y"))
    {
        detector.y.push_back(i);
    }
    detector.resoultion = tree.get<double>(config + ".Resolution");
    detector.accuracy = tree.get<double>(config + ".Accuracy");
    detector.z = tree.get<double>(config + ".z");
    detector.fake_hit_rate = tree.get<double>(config + ".FakeHitrate");
}

void CONFIG::load(const string &configfile)
{
    // Create empty property tree object
    pt::ptree tree;

    // Parse the JSON into the property tree.
    pt::read_json(configfile, tree);

    // Find physical constants
    c = tree.get<double>("PhysicalConstants.SpeedOfLight");
    q = tree.get<double>("PhysicalConstants.ElementaryCharge");
    e_m = tree.get<double>("PhysicalConstants.ElectronMass");
    X0.He = tree.get<double>("PhysicalConstants.ScatteringLength.He");
    X0.Si = tree.get<double>("PhysicalConstants.ScatteringLength.Si");
    X0.C = tree.get<double>("PhysicalConstants.ScatteringLength.C");
    X0.Air = tree.get<double>("PhysicalConstants.ScatteringLength.Air");
    X0.Ta = tree.get<double>("PhysicalConstants.ScatteringLength.Ta");
    X0.Tape = tree.get<double>("PhysicalConstants.ScatteringLength.Tape");
    X0.Mylar = tree.get<double>("PhysicalConstants.ScatteringLength.Mylar");

    // Find the detector parameters
    AssignDetectorParameters("Detectors.M1", tree, M1);
    AssignDetectorParameters("Detectors.M2", tree, M2);
    AssignDetectorParameters("Detectors.M3", tree, M3);
    AssignDetectorParameters("Detectors.M4", tree, M4);
    AssignDetectorParameters("Detectors.M5", tree, M5);
    AssignDetectorParameters("Detectors.M6", tree, M6);
}

int main(int argc, char const *argv[])
{
    try
    {
        CONFIG config;
        config.load("./config/experimental_configuration.json");
        cout << config.M1.x[0] << "\n";
    }
    catch (const std::exception &e)
    {
        std::cerr << "ERROR: " << e.what() << '\n';
    }

    return 0;
}
