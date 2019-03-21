#include "simulator.h"

simulator::simulator(int N, vector<double> z, string run, double BeamEnergy, double CrystalThickness, string filename, string TheorySpectrum, int bg)
{
    // Simulation paramters
    total_events_ = N;                             // number of simulated events
    detector_z_coordinates_ = z;                   // z-coordinates of mimosa detectors
    crystal_type_ = run;                           // suffix for produced data-files
    beam_spatial_distribution_ = (string)filename; // suffix of file containing beam's spaital distro.
    assert("amorphous" == crystal_type_ or "aligned" == crystal_type_ or "background" == crystal_type_ or "alignmentrun" == crystal_type_);
    beam_energy_ = BeamEnergy;            // GeV
    crystal_thicknes_ = CrystalThickness; // micrometer
    total_photon_conversions_ = 0;
    photons_on_foil_ = 0;
    photons_from_crystal_ = 0;
    include_background_radiation_ = bg;
    data_path_ = "/home/christian/Desktop/SimulationTestOutput"; // path to save results
    // Load experimental paramters
    string config_file_path = "./config/experimental_configuration.json";
    LoadConfigFile(config_file_path);
    // Prepare vectors for photon, particle and mimosa data
    particles_.resize(total_events_);
    photons_.resize(total_events_);
    mimosas_.resize(6);
    for (int i = 0; i < 6; i++)
    {
        mimosas_[i].resize(5 * total_events_);
    }
    // Load aligned spectrum
    if (crystal_type_ == "aligned")
    {
        angle_spectrum_ = TheorySpectrum;
        string aligned_crystal_sim = angle_spectrum_;
        LoadDoubles(aligned_crystal_sim, emitted_energies_, intensity_sum_);
        // Do linear interpolation to obtain better energy resoultion
        emitted_energies_interp_ = Linspace(100.0 * 0.5109989461E-03, emitted_energies_.back(), 100 * intensity_sum_.size());
        LinearInterpolation(emitted_energies_, intensity_sum_, emitted_energies_interp_, intensity_sum_interp_);
        for (size_t i = 1; i < intensity_sum_interp_.size(); i++)
        {
            intensity_sum_interp_[i] /= emitted_energies_interp_[i];
        }
        // Calculate integral of energy vs. intensity spectrum. Used to determine photon emission probability
        intensity_integral_ = TrapezoidalIntegrator(emitted_energies_interp_, intensity_sum_interp_);
    }
    else
    {
        initial_spectrum_ = " ";
        angle_spectrum_ = " ";
    }
    // Load alignment of detectors and hot pixels
    alignment_matrix.load("/home/christian/Documents/cern2018/alignment_matrix.txt", arma_ascii);
    hot_pixels_.resize(6);
    LoadHotpixels();
    // Make probability distributions and seed rng for Monte-Carlo methods
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    uniform_real_distribution_.param(uniform_real_distribution<double>(0.0, 1.0).param());
    normal_distribution_.param(normal_distribution<double>(0.0, 1.0).param());
    global_generator_.seed(seed);
    srand(time(NULL));
    // Report progress and loaded variables to terminal
    cout.setf(ios::fixed, ios::floatfield);
    cout.precision(3);
    cout << "\nFinished initializing simulator class.";
    cout << "\nNumber of events:\t" << N;
    cout << "\nCrystal thickness: \t" << crystal_thicknes_ << " micron ";
    cout << "\nBeam energy: \t\t" << beam_energy_ << " GeV";
    if (crystal_type_ == "aligned")
    {
        cout << "\nEmission % / slice: \t" << intensity_integral_ / 300.0;
    }
    cout << "\nConversion rate: \t" << 200.0 * (7.0) / (9.0 * X0_Ta_) * 100 << "%\n";
}

void simulator::AssignDetectorParameters(string config, pt::ptree tree, DETECTOR &detector)
{
    for (auto i : as_vector<double>(tree, config + ".x"))
    {
        detector.x.push_back(i);
    }
    for (auto i : as_vector<double>(tree, config + ".y"))
    {
        detector.y.push_back(i);
    }
    detector.resolution = tree.get<double>(config + ".Resolution");
    detector.accuracy = tree.get<double>(config + ".Accuracy");
    detector.z = tree.get<double>(config + ".z");
    detector.fake_hit_rate = tree.get<double>(config + ".FakeHitrate");
    detector.number = tree.get<int>(config + ".Number");
}

void simulator::LoadConfigFile(const string &configfile)
{
    // Create empty property tree object
    pt::ptree tree;

    // Parse the JSON into the property tree.
    pt::read_json(configfile, tree);

    // Find physical constants
    c_ = tree.get<double>("PhysicalConstants.SpeedOfLight");
    charge_ = tree.get<double>("PhysicalConstants.ElementaryCharge");
    electron_mass_ = tree.get<double>("PhysicalConstants.ElectronMass");
    X0_He_ = tree.get<double>("PhysicalConstants.ScatteringLength.He");
    X0_Si_amorph_ = tree.get<double>("PhysicalConstants.ScatteringLength.Si");
    X0_C_gem_ = tree.get<double>("PhysicalConstants.ScatteringLength.C");
    X0_air_ = tree.get<double>("PhysicalConstants.ScatteringLength.Air");
    X0_Ta_ = tree.get<double>("PhysicalConstants.ScatteringLength.Ta");
    X0_tape_ = tree.get<double>("PhysicalConstants.ScatteringLength.Tape");
    X0_mylar_ = tree.get<double>("PhysicalConstants.ScatteringLength.Mylar");
    foil_thickness_ = tree.get<double>("FoilThickness");

    // Find the detector parameters
    AssignDetectorParameters("Detectors.M1", tree, M1);
    AssignDetectorParameters("Detectors.M2", tree, M2);
    AssignDetectorParameters("Detectors.M3", tree, M3);
    AssignDetectorParameters("Detectors.M4", tree, M4);
    AssignDetectorParameters("Detectors.M5", tree, M5);
    AssignDetectorParameters("Detectors.M6", tree, M6);
}

double simulator::TrapezoidalIntegrator(vector<double> x_coordinate, vector<double> y_coordinate)
{
    double integral = 0.0;
    for (size_t i = 1; i < y_coordinate.size(); i++)
    {
        double dx = x_coordinate[i] - x_coordinate[i - 1];
        integral += 0.5 * dx * (y_coordinate[i - 1] + y_coordinate[i]);
    }
    return integral;
}

int simulator::Coord2Pixel(double x_coordinate, double y_coordinate)
{
    int nrows = 576;
    int ncols = 1152;
    double xmin = -11000, xmax = 11000, ymin = -5500, ymax = 5500;
    // Determine pixelno from coordinates
    double dx = (xmax - xmin) / ncols;
    double dy = (ymax - ymin) / nrows;
    int colno = (x_coordinate + xmax) / dx;
    int rowno = (y_coordinate + ymax) / dy;
    return rowno * ncols + colno;
}

void simulator::LoadHotpixels(void)
{
    for (int i = 0; i < 6; i++)
    {
        string file_name = "/home/christian/Documents/cern2018/simdata/hotpixels_run53_plane_" + to_string(i) + ".txt";
        LoadInt(file_name, hot_pixels_[i]);
        sort(hot_pixels_[i].begin(), hot_pixels_[i].end());
    }
}

void simulator::LoadDoubles(string file_name, vector<double> &data)
{
    ifstream data_file(file_name);
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
        cerr << "Unable to open file: " << file_name << "\n";
    }
}

void simulator::LoadDoubles(string file_name, vector<double> &data0, vector<double> &data1)
{
    ifstream data_file(file_name);
    double val0;
    double val1;
    if (data_file.is_open())
    {
        while (true)
        {
            data_file >> val0 >> val1;

            if (data_file.eof())
                break;

            data0.push_back(val0);
            data1.push_back(val1);
        }

        data_file.close();
    }
    else
    {
        cerr << "Unable to open file: " << file_name << "\n";
    }
}

void simulator::LoadDoubles(string file_name, vector<double> &data0, vector<double> &data1, vector<double> &data2, vector<double> &data3)
{
    ifstream data_file(file_name);
    double val0, val1, val2, val3;
    if (data_file.is_open())
    {
        while (true)
        {
            data_file >> val0 >> val1 >> val2 >> val3;

            if (data_file.eof())
                break;

            data0.push_back(val0);
            data1.push_back(val1);
            data2.push_back(val2);
            data3.push_back(val3);
        }

        data_file.close();
    }
    else
    {
        cerr << "Unable to open file: " << file_name << "\n";
    }
}

void simulator::LoadInt(string file_name, vector<int> &data)
{
    ifstream data_file(file_name);
    int val;
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
        cerr << "Unable to open file: " << file_name << "\n";
    }
}

void simulator::PrintHits(void)
{
    string file_name = data_path_ + "/simulated_hits_coord_data" + crystal_type_ + ".txt";
    ofstream output(file_name);
    for (size_t i = 0; i < mimosas_.size(); i++)
    {
        output << "## PLANE\t" << i << "\tHIT DATA\n"; // block header
        for (size_t j = 0; j < mimosas_[i].size(); j++)
        {
            for (size_t k = 0; k < mimosas_[i][j].size(); k++)
            {
                output << mimosas_[i][j][k][0] << ' ' << mimosas_[i][j][k][1] << '\n';
            }
        }
    }
}

// Generate a beam-profile using measured data and store hits in Events
void simulator::LoadBeamParameters(void)
{
    // Open the text file for reading
    string file_name1 = "../beamParameters/angle_xweight_" + beam_spatial_distribution_;
    string file_name2 = "../beamParameters/angle_yweight_" + beam_spatial_distribution_;
    string file_name3 = "../beamParameters/xweight_" + beam_spatial_distribution_;
    string file_name4 = "../beamParameters/yweight_" + beam_spatial_distribution_;
    string file_name5 = "../beamParameters/angles.txt";
    string file_name6 = "../beamParameters/xdat_" + beam_spatial_distribution_;
    string file_name7 = "../beamParameters/ydat_" + beam_spatial_distribution_;
    LoadDoubles(file_name1, x_angle_weight_);
    LoadDoubles(file_name2, y_angle_weight_);
    LoadDoubles(file_name3, x_coordinate_weight_);
    LoadDoubles(file_name4, y_coordinate_weight_);
    LoadDoubles(file_name5, angles_);
    LoadDoubles(file_name6, x_coordinate_);
    LoadDoubles(file_name7, y_coordinate_);
}

// Propagates particles through the experiment. All radiation lengths are taken from PDG
void simulator::GenerateSyntheticData(void)
{
    double start_time = omp_get_wtime();
    cout << "\nLoading beam parameters";
    LoadBeamParameters();
    cout << "\nBeam parameters succesfully loaded\n";
    int total_emitted = 0, total_detected = 0;
    double total_time = 0;
    int progress = 0;
    cout << "\nSimulating events\n";
    mt19937_64 private_generator = global_generator_;
    arma_rng::set_seed_random();
#pragma omp parallel for
    for (int i = 0; i < total_events_; i++)
    {
        vector<vector<double>> local_photons;
        int x_slope_index = SelectMember(x_angle_weight_);
        int y_slope_index = SelectMember(y_angle_weight_);
        int x_position_index = SelectMember(x_coordinate_weight_);
        int y_position_index = SelectMember(y_coordinate_weight_);
        double x_slope = angles_[x_slope_index];
        double y_slope = angles_[y_slope_index];
        double x_coordinate = x_coordinate_[x_position_index];
        double y_coordinate = y_coordinate_[y_position_index];
        vector<vector<double>> local_particles(1);
        local_particles[0] = {x_coordinate, y_coordinate, detector_z_coordinates_[0], 0, charge_, beam_energy_, x_slope, y_slope};
        int emitted_event = 0;
        if (include_background_radiation_)
        {
            AddPhotons(emitted_event, X0_Si_amorph_, 3.0E+03, private_generator, local_photons, local_particles); // add photons to match measured background spectrum
        }
        // MIMOSA 1 detector
        if (include_background_radiation_)
        {
            AmorphMaterial(emitted_event, X0_tape_, 50.0, 50.0, private_generator, local_photons, local_particles);
        }
        else
        {
            MultipleScattering(X0_tape_, 50.0, private_generator, local_particles);
        }
        if (include_background_radiation_)
        {
            AmorphMaterial(emitted_event, X0_Si_amorph_, 100.0, 50.0, private_generator, local_photons, local_particles);
        }
        else
        {
            MultipleScattering(X0_Si_amorph_, 100.0, private_generator, local_particles);
        }
        MimosaDetector(M1, i, total_detected, private_generator, local_particles);
        // Helium between M1 and M2
        if (include_background_radiation_)
        {
            AmorphMaterial(emitted_event, X0_He_, detector_z_coordinates_[1] - 50.0, detector_z_coordinates_[1] - 150.0, private_generator, local_photons, local_particles);
        }
        else
        {
            MultipleScattering(X0_He_, detector_z_coordinates_[1] - 50.0, private_generator, local_particles);
        }
        // MIMOSA 2 detector
        if (include_background_radiation_)
        {
            AmorphMaterial(emitted_event, X0_tape_, detector_z_coordinates_[1], 50.0, private_generator, local_photons, local_particles);
        }
        else
        {
            MultipleScattering(X0_tape_, detector_z_coordinates_[1], private_generator, local_particles);
        }
        if (include_background_radiation_)
        {
            AmorphMaterial(emitted_event, X0_Si_amorph_, detector_z_coordinates_[1] + 50.0, 50.0, private_generator, local_photons, local_particles);
        }
        else
        {
            MultipleScattering(X0_Si_amorph_, detector_z_coordinates_[1] + 50.0, private_generator, local_particles);
        }
        MimosaDetector(M2, i, total_detected, private_generator, local_particles);
        // Last mylar window He encasing
        if (include_background_radiation_)
        {
            AmorphMaterial(emitted_event, X0_mylar_, detector_z_coordinates_[1] + 100.0, 50.0, private_generator, local_photons, local_particles);
        }
        else
        {
            MultipleScattering(X0_mylar_, detector_z_coordinates_[1] + 100.0, private_generator, local_particles);
        }
        // Air
        if (include_background_radiation_)
        {
            AmorphMaterial(emitted_event, X0_air_, 2060E+03, 2060E+03 - (detector_z_coordinates_[1] + 100.0), private_generator, local_photons, local_particles);
        }
        else
        {
            MultipleScattering(X0_air_, 2060E+03, private_generator, local_particles);
        }
        // Traverse crystal
        if ("amorphous" == crystal_type_)
        { // if amorphous
            AmorphMaterial(emitted_event, X0_C_gem_, 2060E+03 + crystal_thicknes_, crystal_thicknes_, private_generator, local_photons, local_particles, 300);
        }
        else if ("aligned" == crystal_type_)
        { // if aligned

            AlignedCrystal(emitted_event, 300, private_generator, local_photons, local_particles);
        }
        // Air
        if (include_background_radiation_)
        {
            AmorphMaterial(emitted_event, X0_air_, 2310E+03, 2310E+03 - (2060E+03 + crystal_thicknes_), private_generator, local_photons, local_particles);
        }
        else
        {
            MultipleScattering(X0_air_, 2310.0E+03, private_generator, local_particles);
        }
        // Traverse Scintilators
        if (include_background_radiation_)
        {
            AmorphMaterial(emitted_event, X0_Si_amorph_, 2310E+03 + 1.0E+03, 1.0E+03, private_generator, local_photons, local_particles, 300);
        }
        else
        {
            MultipleScattering(X0_Si_amorph_, 2310E+03 + 1.0E+03, private_generator, local_particles);
        }
        // Air between S2 and vacuum tube
        if (include_background_radiation_)
        {
            AmorphMaterial(emitted_event, X0_air_, 2310E+03 + 1.0E+03 + 0.5E+06, 0.5E+06, private_generator, local_photons, local_particles);
        }
        else
        {
            MultipleScattering(X0_air_, 2310E+03 + 1.0E+03 + 0.5E+06, private_generator, local_particles);
        }
        // 1st vacuum tube window
        if (include_background_radiation_)
        {
            AmorphMaterial(emitted_event, X0_mylar_, 2310E+03 + 1.0E+03 + 0.5E+06 + 120.0, 120.0, private_generator, local_photons, local_particles);
        }
        else
        {
            MultipleScattering(X0_mylar_, 2310E+03 + 1.0E+03 + 0.5E+06 + 120.0, private_generator, local_particles);
        }
        if (include_background_radiation_)
        {
            AmorphMaterial(emitted_event, X0_tape_, 2310E+03 + 1.0E+03 + 0.5E+06 + 120.0 + 100.0, 100.0, private_generator, local_photons, local_particles);
        }
        else
        {
            MultipleScattering(X0_tape_, 2310E+03 + 1.0E+03 + 0.5E+06 + 120.0 + 100.0, private_generator, local_particles);
        }
        // If not alignment run
        if ("alignmentrun" != crystal_type_)
        {
            MbplMagnet(local_particles);
        }
        else
        {
            // Air intead of vacuum?
            if (include_background_radiation_)
            {
                AmorphMaterial(emitted_event, X0_air_ * 1e14, 2.8112e+006 + 5.14e+006, 5.14e+006, private_generator, local_photons, local_particles);
            }
            else
            {
                MultipleScattering(X0_air_ * 1e14, 2.8112e+006 + 5.14e+006, private_generator, local_particles);
            }
            // 2nd vacuum tube window
            if (include_background_radiation_)
            {
                AmorphMaterial(emitted_event, X0_mylar_, 2.8112e+006 + 5.14e+006 + 120.0, 120.0, private_generator, local_photons, local_particles);
            }
            else
            {
                MultipleScattering(X0_mylar_, 2.8112e+006 + 5.14e+006 + 120.0, private_generator, local_particles);
            }
            if (include_background_radiation_)
            {
                AmorphMaterial(emitted_event, X0_tape_, 2.8112e+006 + 5.14e+006 + 120.0 + 100.0, 100.0, private_generator, local_photons, local_particles);
            }
            else
            {
                MultipleScattering(X0_tape_, 2.8112e+006 + 5.14e+006 + 120.0 + 100.0, private_generator, local_particles);
            }
            // Air between vacuum tube and M3
            if (include_background_radiation_)
            {
                AmorphMaterial(emitted_event, X0_air_, detector_z_coordinates_[2] - foil_thickness_, detector_z_coordinates_[2] - foil_thickness_ - (2.8112e+006 + 5.14e+006 + 120.0 + 100.0), private_generator, local_photons, local_particles);
            }
            else
            {
                MultipleScattering(X0_air_, detector_z_coordinates_[2] - foil_thickness_, private_generator, local_particles);
            }
        }
        // Traverse converter foil
        ProjectPhotons(local_photons, detector_z_coordinates_[2] - foil_thickness_);
        ConverterFoil(X0_Ta_, private_generator, 300, local_photons, local_particles);
        MultipleScattering(X0_Ta_, detector_z_coordinates_[2], private_generator, local_particles);
        // MIMOSA 3 detector
        if (include_background_radiation_)
        {
            AmorphMaterial(emitted_event, X0_tape_, detector_z_coordinates_[2] + 50.0, 50.0, private_generator, local_photons, local_particles);
        }
        else
        {
            MultipleScattering(X0_tape_, detector_z_coordinates_[2] + 50.0, private_generator, local_particles);
        }
        if (include_background_radiation_)
        {
            AmorphMaterial(emitted_event, X0_Si_amorph_, detector_z_coordinates_[2] + 100.0, 50.0, private_generator, local_photons, local_particles);
        }
        else
        {
            MultipleScattering(X0_Si_amorph_, detector_z_coordinates_[2] + 100.0, private_generator, local_particles);
        }
        MimosaDetector(M3, i, total_detected, private_generator, local_particles);
        // Air between M3 and M4
        if (include_background_radiation_)
        {
            AmorphMaterial(emitted_event, X0_air_, detector_z_coordinates_[3], detector_z_coordinates_[3] - detector_z_coordinates_[2] - 100.0, private_generator, local_photons, local_particles);
        }
        else
        {
            MultipleScattering(X0_air_, detector_z_coordinates_[3], private_generator, local_particles);
        }
        // MIMOSA 4 detector
        if (include_background_radiation_)
        {
            AmorphMaterial(emitted_event, X0_tape_, detector_z_coordinates_[3] + 50.0, 50.0, private_generator, local_photons, local_particles);
        }
        else
        {
            MultipleScattering(X0_tape_, detector_z_coordinates_[3] + 50.0, private_generator, local_particles);
        }
        if (include_background_radiation_)
        {
            AmorphMaterial(emitted_event, X0_Si_amorph_, detector_z_coordinates_[3] + 100.0, 50.0, private_generator, local_photons, local_particles);
        }
        else
        {
            MultipleScattering(X0_Si_amorph_, detector_z_coordinates_[3] + 100.0, private_generator, local_particles);
        }
        MimosaDetector(M4, i, total_detected, private_generator, local_particles);
        // If not alignment run
        if ("alignmentrun" != beam_spatial_distribution_)
        {
            // Air between M4 and middle of MIMOSA magnet
            if (include_background_radiation_)
            {
                AmorphMaterial(emitted_event, X0_air_, (detector_z_coordinates_[4] + detector_z_coordinates_[3]) / 2.0, (detector_z_coordinates_[4] - detector_z_coordinates_[3] - 100.0) / 2.0, private_generator, local_photons, local_particles);
            }
            else
            {
                MultipleScattering(X0_air_, (detector_z_coordinates_[4] + detector_z_coordinates_[3]) / 2.0, private_generator, local_particles);
            }
            // MIMOSA magnet
            MimosaMagnet(local_particles);
            // Air between middle of MIMOSA magnet and M5
            if (include_background_radiation_)
            {
                AmorphMaterial(emitted_event, X0_air_, detector_z_coordinates_[4], (detector_z_coordinates_[4] - detector_z_coordinates_[3]) / 2.0, private_generator, local_photons, local_particles);
            }
            else
            {
                MultipleScattering(X0_air_, detector_z_coordinates_[4], private_generator, local_particles);
            }
        }
        else
        {
            if (include_background_radiation_)
            {
                AmorphMaterial(emitted_event, X0_air_, detector_z_coordinates_[4], detector_z_coordinates_[4] - detector_z_coordinates_[3] - 100.0, private_generator, local_photons, local_particles);
            }
            else
            {
                MultipleScattering(X0_air_, detector_z_coordinates_[4], private_generator, local_particles);
            }
        }
        // MIMOSA 5 detector
        if (include_background_radiation_)
        {
            AmorphMaterial(emitted_event, X0_tape_, detector_z_coordinates_[4] + 50.0, 50.0, private_generator, local_photons, local_particles);
        }
        else
        {
            MultipleScattering(X0_tape_, detector_z_coordinates_[4] + 50.0, private_generator, local_particles);
        }
        if (include_background_radiation_)
        {
            AmorphMaterial(emitted_event, X0_Si_amorph_, detector_z_coordinates_[4] + 100.0, 50.0, private_generator, local_photons, local_particles);
        }
        else
        {
            MultipleScattering(X0_Si_amorph_, detector_z_coordinates_[4] + 100.0, private_generator, local_particles);
        }
        MimosaDetector(M5, i, total_detected, private_generator, local_particles);
        // Air between MIMOSA 5 and MIMOSA 6 detectors
        if (include_background_radiation_)
        {
            AmorphMaterial(emitted_event, X0_air_, detector_z_coordinates_[5], detector_z_coordinates_[5] - detector_z_coordinates_[4] - 100.0, private_generator, local_photons, local_particles);
        }
        else
        {
            MultipleScattering(X0_air_, detector_z_coordinates_[5], private_generator, local_particles);
        }
        // MIMOSA 6 detector
        if (include_background_radiation_)
        {
            AmorphMaterial(emitted_event, X0_tape_, detector_z_coordinates_[5] + 50.0, 50.0, private_generator, local_photons, local_particles);
        }
        else
        {
            MultipleScattering(X0_tape_, detector_z_coordinates_[5] + 50.0, private_generator, local_particles);
        }
        if (include_background_radiation_)
        {
            AmorphMaterial(emitted_event, X0_Si_amorph_, detector_z_coordinates_[5] + 100.0, 50.0, private_generator, local_photons, local_particles);
        }
        else
        {
            MultipleScattering(X0_Si_amorph_, detector_z_coordinates_[5] + 100.0, private_generator, local_particles);
        }
        MimosaDetector(M6, i, total_detected, private_generator, local_particles);
#pragma omp critical
        {
            total_emitted += emitted_event;
            progress++;
            photons_[progress] = local_photons;
            particles_[progress] = local_particles;
        }
        // Report progress in terminal
        if ((progress + 1) % (total_events_ / 10 + 1) == 0)
        {
            double dt = omp_get_wtime() - start_time;
            total_time += dt;
            cout << "Progress :\t" << floor(100 * double(progress) / (double)total_events_) << "%"
                 << "\ttime used :\t" << dt << "\ttotal time elapsed :\t" << total_time << "\ttime remaining :\t" << dt * (double)total_events_ / (total_events_ / 10 + 1) - total_time << "\n";
            start_time = omp_get_wtime();
        }
    }
    // Report results to terminal
    cout << "\n Finished constructing tracks\n";
    cout << "\nTotal photons emitted: " << total_emitted << "\n";
    cout << "Total photons emitted by crystal: " << photons_from_crystal_ << "\n";
    cout << "Signal-to-Noise Ratio (SNR): " << double(photons_from_crystal_) / double(total_emitted - photons_from_crystal_) << "\n";
    cout << "Total photons incident on foil: " << photons_on_foil_ << "\n";
    cout << "Total conversions: " << total_photon_conversions_ << "\n";
}

// Simulate a particle travelling through an amorphous crystal.
void simulator::AddPhotons(int &emitted_event, double X0, double d, mt19937_64 generator, vector<vector<double>> &local_photons, vector<vector<double>> local_particles)
{
    double minimum_energy = 2 * 0.5109989461E-03; // 2 * electron mass in GeV
    for (size_t i = 0; i < local_particles.size(); i++)
    {
        double particle_energy = local_particles[i][5];
        if (particle_energy < minimum_energy)
        {
            break;
        }
        // Determine if photon is emitted
        double random_number = (double)rand() / RAND_MAX;
        double emission_probability = d / X0 * ((4.0 / 3.0) * log(beam_energy_ / minimum_energy) - (4.0 / 3.0) * (beam_energy_ - minimum_energy) / beam_energy_ + (3.0 / 4.0) * (1.0 - minimum_energy * minimum_energy / (2.0 * beam_energy_ * beam_energy_))); // BH-cross section * crystal length
        if (random_number < emission_probability)
        {
            // Determine photon energy
            double random_number_2 = (double)rand() / RAND_MAX;
            particle_energy = local_particles[i][5];
            double norm = 4.0 / 3.0 * log(particle_energy / minimum_energy) - 4.0 / (3.0 * particle_energy) * (particle_energy - minimum_energy) + 1.0 / (2.0 * particle_energy * particle_energy) * (particle_energy * particle_energy - minimum_energy * minimum_energy);
            function<double(vec)> energy = [random_number_2, particle_energy, norm](vec x) {
                return PhotonicEnergyDistribution(x, random_number_2, particle_energy, norm);
            }; // make lambda-function in order to use same random_number during iteration
            vec simplex_corner_1 = {0.1};
            vec simplex_corner_2 = {20.0};
            vector<vec> initial_simplex = {simplex_corner_1, simplex_corner_2};      // initial simplex for Nelder-Mead. The initial guess is hugely important for convergence
            vec photon_energy = SimplexNelderMead(energy, initial_simplex, 1.0E-08); // solve for energy using Nelder-Mead simplex.
            // Determine direction of photon
            double gamma = local_particles[i][5] / (5.109989461E-4);
            uniform_real_distribution<double> deflection_angle(-1.0 / gamma, 1.0 / gamma);
            // Update particles and photons vectors
            vector<double> photon = local_particles[i];
            photon[4] = 0.0; // charge
            photon[5] = photon_energy(0);
            double dx = deflection_angle(generator);
            double dy = sqrt(1.0 / (gamma * gamma) - dx * dx);
            photon[6] += dx;
            photon[7] += dy;
            local_photons.push_back(photon);
            emitted_event++;
        }
    }
}

void simulator::AmorphMaterial(int &emitted_event, double X0, double z, double length, mt19937_64 generator, vector<vector<double>> &local_photons, vector<vector<double>> &local_particles, int slices)
{
    double minimum_energy = 2.0 * 0.5109989461E-03;
    double slice_length = length / (double)slices;
    for (size_t i = 0; i < local_particles.size(); i++)
    {
        for (int j = 0; j < slices; j++)
        {
            double particle_energy = local_particles[i][5];
            // Multiple scattering in C crystal slice. Calculate x/y displacement and deflection independently
            for (int k = 0; k < 2; k++)
            {
                double z1 = randn();
                double z2 = randn();
                double theta0 = 0.0136 / particle_energy * sqrt(slice_length / X0) * (1.0 + 0.038 * log(slice_length / X0));
                double dy = slice_length * theta0 * (z1 / sqrt(12) + z2 / 2.0);
                double dtheta0 = z2 * theta0;
                local_particles[i][k] += dy + slice_length * local_particles[i][6 + k]; // total (x/y) displacement
                local_particles[i][6 + k] += dtheta0;                                   // x_slope/y_slope
            }
            local_particles[i][2] += slice_length;
            // Proceed only if photon is emitted
            double random_number = (double)rand() / RAND_MAX;                                                                                                                                                                                                                    // random double between 0, 1
            bool emission = random_number < slice_length / X0 * ((4.0 / 3.0) * log(beam_energy_ / minimum_energy) - (4.0 / 3.0) * (beam_energy_ - minimum_energy) / beam_energy_ + (3.0 / 4.0) * (1.0 - minimum_energy * minimum_energy / (2.0 * beam_energy_ * beam_energy_))); // BH-cross section * crystal length
            if (emission and particle_energy > minimum_energy)
            {
                if (X0 == X0_C_gem_)
                {
#pragma omp atomic
                    photons_from_crystal_++;
                }
                // Determine photon energy
                double random_number_2 = (double)rand() / RAND_MAX;
                double norm = 4.0 / 3.0 * log(particle_energy / minimum_energy) - 4.0 / (3.0 * particle_energy) * (particle_energy - minimum_energy) + 1.0 / (2.0 * particle_energy * particle_energy) * (particle_energy * particle_energy - minimum_energy * minimum_energy);
                vec simplex_corner_1 = {0.001}, simplex_corner_2 = {40.0};
                vector<vec> initial_simplex = {simplex_corner_1, simplex_corner_2};
                vec photon_energy = SimplexNelderMead([random_number_2, particle_energy, norm, minimum_energy](vec x) {
                    return x(0) > 0 ? abs((4.0 / 3.0 * log(x(0) / minimum_energy) - 4.0 / 3.0 * (x(0) - minimum_energy) / particle_energy + 1.0 / (2.0 * particle_energy * particle_energy) * (x(0) * x(0) - minimum_energy * minimum_energy)) - random_number_2 * norm) : 1E+17;
                },
                                                      initial_simplex, 1.0E-08);
                // Determine direction of photon
                photon_energy(0) = 20.0;
                double gamma = local_particles[i][5] / (5.109989461E-4);
                uniform_real_distribution<double> deflection_angle(-1.0 / gamma, 1.0 / gamma);
                double x_deflection_angle = deflection_angle(generator);
                double y_deflection_angle = sqrt(1.0 / (gamma * gamma) - x_deflection_angle * x_deflection_angle);
                // Update particles and photons vectors
                vector<double> photon = local_particles[i];
                photon[4] = 0.0; // charge
                photon[5] = photon_energy(0);
                photon[6] += x_deflection_angle;
                photon[7] += y_deflection_angle;
                local_particles[i][5] -= photon_energy(0);
                local_photons.push_back(photon);
#pragma omp atomic
                emitted_event++;
            }
        } // slices end
    }     // events end
}

void simulator::ConverterFoil(double X0, mt19937_64 generator, int slices, vector<vector<double>> &local_photons, vector<vector<double>> &local_particles)
{
    double slice_thickness = foil_thickness_ / (double)slices;
    uniform_real_distribution<double> rotation(0.0, 2.0 * M_PI);
    double conversion_probability = slice_thickness * (7.0) / (9.0 * X0);
    for (size_t i = 0; i < local_photons.size(); i++)
    {
#pragma omp atomic
        photons_on_foil_ += 1;
        for (int j = 0; j < slices; j++)
        {
            // Proceed only if a conversion happens
            double random_number = (double)rand() / RAND_MAX;
            if (random_number < conversion_probability)
            {
                // Calculate energy gained by e+/e- pair
                double photon_energy = local_photons[i][5];
                double electron_energy = (double)rand() / RAND_MAX; // random number for the inverse transform sampling in "ElectronicEnergyDistribution"
                ElectronicEnergyDistribution(electron_energy);      // the fractional electron energy, ie E_e-/ E_photon
                electron_energy *= photon_energy;
                double positron_energy = photon_energy - electron_energy; // energy conservation
                // Calculate deflection of e+/e- pair
                double electron_defl;
                double positron_defl;
                BorsellinoOpeningAngle(electron_energy, positron_energy, photon_energy, electron_defl, positron_defl); // deflection angle based on approximated Borsellino distribution
                double emission_angle = rotation(generator);
                // Add e+/e- pair to "particles" array
                vector<double> electron = local_photons[i];
                electron[2] += slice_thickness * (j + 1);
                electron[4] = (-1.0) * charge_;
                electron[5] = electron_energy;
                electron[6] += cos(emission_angle) * electron_defl;
                electron[7] += sin(emission_angle) * electron_defl;
                vector<double> positron = local_photons[i];
                positron[2] += slice_thickness * (j + 1);
                positron[4] = charge_;
                positron[5] = positron_energy;
                positron[6] += cos(emission_angle) * positron_defl;
                positron[7] += sin(emission_angle) * positron_defl;
                local_particles.push_back(electron);
                local_particles.push_back(positron);
#pragma omp atomic
                total_photon_conversions_++;
                break; // break since photon no longer exists
            }
        }
    }
}

void simulator::MbplMagnet(vector<vector<double>> &local_particles)
{
    local_particles.clear(); // the MPBL magnet removes all particles
}

// Simulate a particle travelling through an aligned crystal. Multiple scattering is taken into account in the "propagte particles" function.
void simulator::AlignedCrystal(int &emitted_event, int slices, mt19937_64 generator, vector<vector<double>> &local_photons, vector<vector<double>> &local_particles)
{
    double slice_thickness = crystal_thicknes_ / (double)slices;
    for (size_t i = 0; i < local_particles.size(); i++)
    {
        for (int j = 0; j < slices; j++)
        {
            double particle_energy = local_particles[i][5];
            // Multiple scattering in C crystal slice. Calculate x/y displacement and deflection independently
            for (int k = 0; k < 2; k++)
            {
                double z1 = randn();
                double z2 = randn();
                double theta0 = 0.0136 / particle_energy * sqrt(slice_thickness / X0_C_gem_) * (1.0 + 0.038 * log(slice_thickness / X0_C_gem_));
                double dy = slice_thickness * theta0 * (z1 / sqrt(12) + z2 / 2.0);
                double dtheta0 = z2 * theta0;
                local_particles[i][k] += dy + slice_thickness * local_particles[i][6 + k]; // total (x/y) displacement
                local_particles[i][6 + k] += dtheta0;                                      // x_slope/y_slope
            }
            local_particles[i][2] += slice_thickness;
            // Proceed only if photon is emitted
            double random_number = (double)rand() / RAND_MAX;
            bool emission = random_number < (intensity_integral_ / (double)slices);
            if (emission)
            {
#pragma omp atomic
                photons_from_crystal_++;
                // Determine photon energy
                int E_indx = SelectMember(intensity_sum_interp_);
                // Determine direction of emission of photon
                double gamma = local_particles[i][5] / (5.109989461E-4);
                uniform_real_distribution<double> deflection_angle(-1.0 / gamma, 1.0 / gamma);
                double dx = deflection_angle(generator);
                double dy = sqrt(1.0 / (gamma * gamma) - dx * dx);
                // Add photon to "photons" vector
                vector<double> photon = local_particles[i];
                photon[4] = 0.0; // charge
                photon[5] = emitted_energies_interp_[E_indx];
                photon[6] += dx;
                photon[7] += dy;
                local_photons.push_back(photon);
#pragma omp atomic
                emitted_event++;
            }
        } // end slices
    }     // end events
}

void simulator::MimosaDetector(DETECTOR detector, int event, int &detections, mt19937_64 generator, vector<vector<double>> local_particles)
{
    // Take into account detector resoultion by adding a random (Dx, Dy) displacement drawn from gaussian distro with stdev = 3.5mu, mean = 0mu
    normal_distribution<double> distribution(0.0, detector.accuracy);
    // Simulate real hit detection
    for (size_t i = 0; i < local_particles.size(); i++)
    {
        // Check if particle is within physical boundaries of detector.
        bool inside_bounds = IsInsidePolygon(4, detector.x, detector.y, local_particles[i][0], local_particles[i][1]);
        // inside_bounds = true;
        if (inside_bounds)
        {
            double Dx = distribution(generator);
            double Dy = distribution(generator);
            // Combine neighbor-hits into single hit if they are closer than detectors resoultion
            if (mimosas_[detector.number][event].size() > 0)
            {
                double dist = CalculateDistance(local_particles[i][0] + Dx, local_particles[i][1] + Dy, mimosas_[detector.number][event][0][0], mimosas_[detector.number][event][0][1]);
                // Save hit in "mimosas_" array. This array has the same structure and contains equivalent data as the "hitcoords" array in the "analyser" class
                if (dist > detector.resolution)
                {
                    mimosas_[detector.number][event].push_back({local_particles[i][0] + Dx, local_particles[i][1] + Dy});
                }
                else
                {
                    // The distance is shorter than detector resolution, so we combine two hits. This is equivalent to moving the previously recorded hit by half the seperation
                    double dx = mimosas_[detector.number][event][0][0] - local_particles[i][0] - Dx;
                    double dy = mimosas_[detector.number][event][0][1] - local_particles[i][1] - Dy;
                    mimosas_[detector.number][event][0][0] += dx / 2.0;
                    mimosas_[detector.number][event][0][1] += dy / 2.0;
                }
            }
            else
            {
                // If there are no previous hits in detector simply add the hit to "mimoas" array
                mimosas_[detector.number][event].push_back({local_particles[i][0] + Dx, local_particles[i][1] + Dy});
            }
        }
    }
    // Simulate fake hit detection
    uniform_real_distribution<double> fake_hit_x(detector.x[0], detector.x[4]);
    uniform_real_distribution<double> fake_hit_y(detector.y[0], detector.y[4]);
    if ((double)rand() / RAND_MAX < detector.fake_hit_rate)
    {
        double x_coordinate = fake_hit_x(generator);
        double y_coordinate = fake_hit_y(generator);
        vector<double> hit = {x_coordinate, y_coordinate};
        mimosas_[detector.number][event].push_back(hit);
    }
}

double simulator::CalculateDistance(double x0, double y0, double x1, double y1)
{
    return sqrt(pow(x1 - x0, 2) + pow(y1 - y0, 2));
}

void simulator::ProjectPhotons(vector<vector<double>> &local_photons, double z_coordinate)
{
    for (size_t i = 0; i < local_photons.size(); i++)
    {
        double dx = (z_coordinate - local_photons[i][2]) * local_photons[i][6];
        double dy = (z_coordinate - local_photons[i][2]) * local_photons[i][7];
        local_photons[i][0] += dx;
        local_photons[i][1] += dy;
        local_photons[i][2] = z_coordinate;
    }
}

void simulator::MultipleScattering(double X0, double z_coordinate, mt19937_64 generator, vector<vector<double>> &local_particles)
{
    // Iterate over paticles traversing medium
    for (size_t i = 0; i < local_particles.size(); i++)
    {
        double length = z_coordinate - local_particles[i][2];
        // Calculate x/y displacement and deflection independently
        for (int j = 0; j < 2; j++)
        {
            double random_number_1 = randn();
            double random_number_2 = randn();
            double theta0 = 0.0136 / local_particles[i][5] * sqrt(length / X0) * (1.0 + 0.038 * log(length / X0));
            double displacement = length * theta0 * (random_number_1 / sqrt(12) + random_number_2 / 2.0);
            double direction_change = random_number_2 * theta0;
            local_particles[i][j] += displacement + length * local_particles[i][6 + j]; // total (x/y) displacement
            local_particles[i][6 + j] += direction_change;                              // x_slope/y_slope
        }
        local_particles[i][2] = z_coordinate;
    }
}

// Calculate the deflection of a charged particle traversing the mimosa magnet. Magnet only deflects in x-direction.
void simulator::MimosaMagnet(vector<vector<double>> &local_particles)
{
    double length = 0.15;         // length traveled through Mimosa magnet in m
    double field_strength = 0.12; // strength of magnetic field in T
    for (size_t hitno = 0; hitno < local_particles.size(); hitno++)
    {
        local_particles[hitno][6] += (local_particles[hitno][4] * length * field_strength * c_) / (local_particles[hitno][5] * 1.6021766E-10); // update direction of particle
    }
}

void simulator::BorsellinoOpeningAngle(double energy_particle_1, double energy_particle_2, double energy_photon, double &deflection_particle_1, double &deflection_particle_2)
{
    double u1 = (double)rand() / RAND_MAX, u2 = (double)rand() / RAND_MAX;
    double nu0 = energy_photon * electron_mass_ * c_ / (energy_particle_1 * energy_particle_2) * 1.6021766E+10;
    double nu = nu0 * (sqrt(1.0 / (1.0 - u1) - 1.0) + 0.70 * u2);
    deflection_particle_1 = nu / (energy_particle_1 / energy_particle_2 + (1.0 - nu * nu / 2.0)); // deflection of particle 1
    deflection_particle_2 = nu - deflection_particle_1;                                           // deflection of particle 2
}

// Analytical solution to CDF(x) = random_number (see Flohr's thesis) integrated and solved 3rd order poly. equation (2.58)
void simulator::ElectronicEnergyDistribution(double &random_number)
{
    random_number = 0.25 * (1.5874 * pow(sqrt(196.0 * random_number * random_number - 196.0 * random_number + 81.0) + 14 * random_number - 7.0, 1.0 / 3.0) - 5.03968 / pow(sqrt(196.0 * random_number * random_number - 196.0 * random_number + 81.0) + 14.0 * random_number - 7.0, 1.0 / 3.0) + 2.0);
}

double simulator::PhotonicEnergyDistribution(vec energy, double random_number, double beam_energy, double norm)
{
    double minimum_energy = 2.0 * 0.5109989461E-3; // 2 * electron mass GeV
    return energy(0) > 0 ? abs((4.0 / 3.0 * log(energy(0) / minimum_energy) - 4.0 / 3.0 * (energy(0) - minimum_energy) / beam_energy + 1.0 / (2.0 * beam_energy * beam_energy) * (energy(0) * energy(0) - minimum_energy * minimum_energy)) - random_number * norm) : 1E+17;
}

int simulator::BinarySearch(vector<int> numbers, int low, int high, int value)
{
    if (high >= low)
    {
        int mid = (high + low) / 2;
        if (numbers[mid] == value)
        {
            return mid;
        }
        if (numbers[mid] > value)
        {
            return BinarySearch(numbers, low, mid - 1, value);
        }
        return BinarySearch(numbers, mid + 1, high, value);
    }
    return -1; // returned only if "value" is not in "numbers"
}

int simulator::BinarySearch(vector<double> numbers, int low, int high, double value)
{
    if (high >= low)
    {
        int mid = (high + low) / 2;
        if (abs(numbers[mid] - value) < 1e-05)
        {
            return mid;
        }
        if (numbers[mid] > value)
        {
            return BinarySearch(numbers, low, mid - 1, value);
        }
        return BinarySearch(numbers, mid + 1, high, value);
    }
    return -1; // returned only if "value" is not in "numbers"
}

// This function takes a test point (x_coordinate, y_coordinate) and checks if it is within an vertices-sided polygon defined by the vertices vertex_x_coordinate and vertex_y_coordinate. Function is a slightly modified version of the one posted on https://wrf.ecse.rpi.edu//Research/Short_Notes/pnpoly.html
int simulator::IsInsidePolygon(int vertices, vector<double> vertex_x_coordinate, vector<double> vertex_y_coordinate, double x_coordinate, double y_coordinate)
{
    int i, j, c = 0;
    for (i = 0, j = vertices - 1; i < vertices; j = i++)
    {
        if (((vertex_y_coordinate[i] > y_coordinate) != (vertex_y_coordinate[j] > y_coordinate)) && (x_coordinate < (vertex_x_coordinate[j] - vertex_x_coordinate[i]) * (y_coordinate - vertex_y_coordinate[i]) / (vertex_y_coordinate[j] - vertex_y_coordinate[i]) + vertex_x_coordinate[i]))
        {
            c = !c;
        }
    }
    return c;
}

// This function selects a member from a weighted ("weights") list and returns the indx
int simulator::SelectMember(vector<double> weights)
{
    // Calculate cumulative sum of weights to represent the probability that a number will be picked
    vector<double> weights_cumsum;
    weights_cumsum.resize(weights.size());
    partial_sum(weights.begin(), weights.end(), weights_cumsum.begin());
    uniform_real_distribution<double> distribution(0.0, weights_cumsum.back());
    // Generate random number, used to determine index
    double r = distribution(global_generator_);
    // Binary search to find which number to pick based on weight
    int low = 0, high = weights.size() - 1;
    while (high >= low)
    {
        int guess = (low + high) / 2;
        if (weights_cumsum[guess] < r)
        {
            low = guess + 1;
        }
        else if (weights_cumsum[guess] - weights[guess] > r)
        {
            high = guess - 1;
        }
        else
        {
            return guess; // return the index for selected member in "numbers" array
        }
    }
    return -1; // if this is returned, somehting has gone wrong
}

vector<double> simulator::Linspace(double min, double max, int N)
{
    vector<double> range;
    range.resize(N);
    double delta = (max - min) / (N - 1);
    for (int i = 0; i < N - 1; i++)
    {
        range[i] = min + i * delta;
    }
    range.back() = max; // Set last entry to "max". This ensures that we get a range from min -> max
    return range;
}

void simulator::LinearInterpolation(vector<double> x_coordinates, vector<double> y_coordinates, vector<double> x_interpolated_coordinates, vector<double> &y_interpolated_coordinates)
{
    y_interpolated_coordinates.resize(x_interpolated_coordinates.size());
    // Iterate over interpolated values
    for (size_t j = 0; j < x_interpolated_coordinates.size(); j++)
    {
        // Find the interval over which to interpolate with binary search
        int low = 0;
        int high = x_coordinates.size() - 1;
        while (high - low > 1)
        {
            int guess = (high + low) / 2;
            if (x_interpolated_coordinates[j] > x_coordinates[guess])
                low = guess;
            else
                high = guess;
        }
        // Calculate slope and do interpolation
        double dydx = (y_coordinates[low + 1] - y_coordinates[low]) / (x_coordinates[low + 1] - x_coordinates[low]);
        y_interpolated_coordinates[j] = y_coordinates[low] + dydx * (x_interpolated_coordinates[j] - x_coordinates[low]);
    }
}

// Reflect highest point against centroid
vec simulator::SimplexReflect(vec highest, vec centroid)
{
    return 2.0 * centroid - highest;
}

// Reflect highest point against centroid, then double distance to centroid
vec simulator::SimplexExpand(vec highest, vec centroid)
{
    return 3.0 * centroid - 2.0 * highest;
}

// Highest point halves its distance from centroid
vec simulator::SimplexContract(vec highest, vec centroid)
{
    return 0.5 * (centroid + highest);
}

// All points, except lowest, halves their distance to lowest point
void simulator::SimplexReduce(vector<vec> &simplex, int low_index)
{
    for (size_t i = 0; i < simplex.size(); i++)
    {
        if ((int)i != low_index)
        {
            simplex[i] = 0.5 * (simplex[i] + simplex[low_index]);
        }
    }
}

// Calculate size of the simplex-polygon
double simulator::SimplexSize(vector<vec> simplex)
{
    double current_size = norm(simplex[0] - simplex.back(), 2);
    for (size_t i = 1; i < simplex.size(); i++)
    {
        double new_size = norm(simplex[i] - simplex[i - 1], 2);
        if (current_size < new_size)
        {
            current_size = new_size;
        }
    }
    return current_size;
}

// Update high_index, low_index and centroid.
void simulator::SimplexUpdate(vector<vec> simplex, vec function, vec &centroid, int &high_index, int &low_index)
{
    high_index = 0;
    double function_high_value = function(0);
    low_index = 0;
    double function_low_value = function(0);
    for (size_t i = 1; i < simplex.size(); i++)
    {
        if (function(i) > function_high_value)
        {
            function_high_value = function(i);
            high_index = i;
        }
        if (function(i) < function_low_value)
        {
            function_low_value = function(i);
            low_index = i;
        }
    }
    vec s;
    s.resize(simplex[0].size());
    for (int i = 0; i < (int)simplex.size(); i++)
    {
        if (i != high_index)
        {
            s += simplex[i];
        }
    }
    for (int i = 0; i < (int)s.size(); i++)
    {
        centroid(i) = s(i) / (simplex.size() - 1);
    }
}

void simulator::SaveVector(string name, vector<double> data)
{
    string file_name = data_path_ + "/" + name + ".txt";
    ofstream output(file_name);
    for (size_t i = 0; i < data.size(); i++)
    {
        output << data[i] << "\n";
    }
}

void simulator::SaveVector(string name, vector<vector<double>> data)
{
    string file_name = data_path_ + "/" + name + ".txt";
    ofstream output(file_name);
    for (size_t i = 0; i < data.size(); i++)
    {
        for (size_t j = 0; j < data[i].size(); j++)
        {
            output << data[i][j] << "\t";
        }
        output << "\n";
    }
}

void simulator::SaveVector(string name, vector<vector<vector<double>>> data)
{
    string file_name = data_path_ + "/" + name + ".txt";
    ofstream output(file_name);
    for (size_t i = 0; i < data.size(); i++)
    {
        for (size_t j = 0; j < data[i].size(); j++)
        {
            for (size_t k = 0; k < data[i][j].size(); k++)
            {
                output << data[i][j][k] << "\t";
            }
            output << "\n";
        }
    }
}
