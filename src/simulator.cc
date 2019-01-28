#include "simulator.h"
#include <fstream>
/* Look up global variables (globalvars.h) */
extern default_random_engine generator;
extern normal_distribution<double> R_normal;
extern uniform_real_distribution<double> R;

struct point
{
        double angx;
        double angy;
        double posx;
        double posy;
};

simulator::simulator(int N, vector<double> z, char const *run, vector<double> params, char const *filename, char const* s1, int bg) {

        /* Simulation paramters */
        Nevents = N; // number of simulated events
        zplanes = z; // z-coordinates of mimosa detectors
        name = run; // suffix for produced data-files
        beam_spatial_distro = (string) filename; // suffix of file containing beam's spaital distro.
        assert (!strcmp("amorphous", name) or !strcmp("aligned", name) or !strcmp("background", name) or !strcmp("alignmentrun", name));
        E = params[0]; // energy of beam [GeV]
        d_c = params[1]; // thickness of crystal target [micrometer]
        conversions = 0;
        no_photons = 0;
        emitted_crystal = 0;
        bg_include = bg;
        DATPATH = "/home/christian/Documents/cern2018/simdata"; // path to save results

        /* Physics constants */
        q = -1.6021766208E-19; // charge of electron in Coulomb
        c = 299792458; // speed of light in vac. in m/s
        m = 9.10938356E-31; // electron/positron mass in kg

        /* Radiation length of different materials. Used to calculate SAMS of particle. Units micro-meter (1e4 micrometer = 1 cm)*/
        X0_Si_amorph = 9.370E+04;
        X0_C_amorph = 21.35E+04;
        X0_C_gem = 12.13E+04;
        X0_Mimosa = 9.370E+04;
        X0_He = 5.671E+09;
        X0_air = 3.039E+08;
        X0_Ta = 0.4094E+04;
        X0_tape = 19.63E+04;
        X0_Mylar =  28.54E+04;

        /* Prepare vectors for photon, particle and mimosa data */
        particles.resize(Nevents);
        photons.resize(Nevents);
        mimosas.resize(6);

        for (int i = 0; i < 6; i++) {

                mimosas[i].resize(5*Nevents);

        }

        if (!strcmp("aligned", name)) {

                angle_spec = s1;

                string aligned_crystal_sim = angle_spec; //sum_angles40GeV_full_peak.txt";
                load_doubles(aligned_crystal_sim, emitted_energies, intensity_sum);

                /* Calculate integral of energy vs. intensity spectrum. Used to determine photon emission probability */
                I_integral = 0;

                for (size_t i = 1; i < intensity_sum.size(); i++) {
                        double dE = (emitted_energies[i] - emitted_energies[i-1]);
                        double emitted_energy = 0.5 * (emitted_energies[i] + emitted_energies[i-1]);
                        I_integral += (0.5 * dE * (intensity_sum[i-1] + intensity_sum[i])) / emitted_energy;
                }

                cerr << "\nIntegral = " << I_integral << "\n";

                /* Do linear interpolation to obtain better energy resoultion */
                energies_interp = linspace( 2*0.5109989461E-03, emitted_energies.back(), 4*intensity_sum.size());

                linterp(emitted_energies, intensity_sum, energies_interp, intensity_sum_interp);

                for (size_t i = 1; i < intensity_sum_interp.size(); i++) {

                        intensity_sum_interp[i] /= energies_interp[i];

                }

        } else {

                initial_spec = " ";
                angle_spec = " ";

        }

        hotpixels.resize(6);
        load_hotpixels();


        /* Make probability distributions and seed rng */
        R.param(uniform_real_distribution<double>(0.0, 1.0).param());
        R_normal.param(normal_distribution<double>(0.0, 1.0).param());
        generator.seed(random_device{} ());

        /* Report progress to terminal */
        cerr << "\nFinished initializing simulator class.\n";
        cerr << "\nNumber of events:\t" << N;
        cerr << "\nCrystal thickness: \t" << d_c << " micron ";
        cerr << "\nBeam energy: \t\t" << E << " GeV\n\n";

}

int simulator::coord2pixel(double xhit, double yhit) {

        int nrows = 576;
        int ncols = 1152;
        double xmin = -11000, xmax = 11000, ymin = -5500, ymax = 5500;

        /* Determine pixelno from coordinates */
        double dx = (xmax - xmin) / ncols;
        double dy = (ymax - ymin) / nrows;
        int colno = (xhit + xmax)/dx;
        int rowno = (yhit + ymax)/dy;

        return rowno * ncols + colno;

}

void simulator::load_hotpixels(void) {

        for (int i = 0; i < 6; i++) {

                string filename = "/home/christian/Documents/cern2018/simdata/hotpixels_run53_plane_" + to_string(i) + ".txt";
                load_int(filename, hotpixels[i]);
                sort(hotpixels[i].begin(), hotpixels[i].end());

        }

}

void simulator::load_doubles(string filename, vector<double> &data) {

        ifstream datafile (filename);
        double val;

        if (datafile.is_open()) {

                while (true) {

                        datafile >> val;

                        if (datafile.eof()) break;

                        data.push_back(val);

                }

                datafile.close();

        } else cout << "Unable to open file: " << filename << "\n";

}

void simulator::load_doubles(string filename, vector<double> &data0, vector<double> &data1) {

        ifstream datafile (filename);
        double val0;
        double val1;

        if (datafile.is_open()) {

                while (true) {

                        datafile >> val0 >> val1;

                        if (datafile.eof()) break;

                        data0.push_back(val0);
                        data1.push_back(val1);

                }

                datafile.close();

        } else cout << "Unable to open file: " << filename << "\n";

}

void simulator::load_doubles(string filename, vector<double> &data0, vector<double> &data1, vector<double> &data2, vector<double> &data3) {

        ifstream datafile (filename);
        double val0, val1, val2, val3;

        if (datafile.is_open()) {

                while (true) {

                        datafile >> val0 >> val1 >> val2 >> val3;

                        if (datafile.eof()) break;

                        data0.push_back(val0);
                        data1.push_back(val1);
                        data2.push_back(val2);
                        data3.push_back(val3);


                }

                datafile.close();

        } else cout << "Unable to open file: " << filename << "\n";

}


void simulator::load_int(string filename, vector<int> &data) {

        ifstream datafile (filename);
        int val;

        if (datafile.is_open()) {

                while (true) {

                        datafile >> val;

                        if (datafile.eof()) break;

                        data.push_back(val);

                }

                datafile.close();

        } else cout << "Unable to open file: " << filename << "\n";

}

void simulator::print_hits(void) {

        string filename = DATPATH + "/simulated_hits_coord_data" + name + ".txt";
        ofstream output (filename);

        for (size_t i = 0; i < mimosas.size(); i++) { // planes

                output << "## PLANE\t" << i << "\tHIT DATA\n"; // block header

                for (size_t j = 0; j < mimosas[i].size(); j++) { // events

                        for (size_t k = 0; k < mimosas[i][j].size(); k++) { // hits

                                output << mimosas[i][j][k][0] << ' ' << mimosas[i][j][k][1] << '\n';

                        } // end hits

                } // end events

        } // end planes

}

void simulator::print_energy(string name) {

        string filename = DATPATH + "/photon_energy_sim" + name + ".txt";
        ofstream output (filename);

        for (size_t i = 0; i < energies.size(); i++) {

                output << energies[i] << "\n";

        }

}

/* Generate a beam-profile using measured data and store hits in Events */
void simulator::generate_beam(void) {

        /* A collection of points structures + decleration of variables for temporarily storing poins */
        vector<double> xaw;
        vector<double> yaw;
        vector<double> xpw;
        vector<double> ypw;
        vector<double> ax;
        vector<double> ay;
        vector<double> x;
        vector<double> y;
        vector<double> a;

        /* Open the text file for reading */
        string filename1 = "../beamParameters/angle_xweight_" + beam_spatial_distro;
        string filename2 = "../beamParameters/angle_yweight_" + beam_spatial_distro;
        string filename3 = "../beamParameters/xweight_" + beam_spatial_distro;
        string filename4 = "../beamParameters/yweight_" + beam_spatial_distro;
        string filename5 = "../beamParameters/angles.txt";
        string filename6 = "../beamParameters/xdat_" + beam_spatial_distro;
        string filename7 = "../beamParameters/ydat_" + beam_spatial_distro;

        load_doubles(filename1, xaw);
        load_doubles(filename2, yaw);
        load_doubles(filename3, xpw);
        load_doubles(filename4, ypw);
        load_doubles(filename5, a);
        load_doubles(filename6, x);
        load_doubles(filename7, y);

        int nohits = 1;

        vector<vector<double> > hits(Nevents);


        for (int i = 0; i < Nevents; i++) {

                particles[i].resize(nohits);

                for (int j = 0; j < nohits; j++) {

                        int xslope_indx = select_member(xaw);
                        int yslope_indx = select_member(yaw);
                        int xpos_indx = select_member(xpw);
                        int ypos_indx = select_member(ypw);
                        double xslope = a[xslope_indx];
                        double yslope = a[yslope_indx];
                        double xcoord = x[xpos_indx];
                        double ycoord = y[ypos_indx];

                        vector<double> hitdata = {xcoord, ycoord, zplanes[0], 0, q, E, xslope, yslope}; // x, y, z, plane, charge, energy of particle, xslope, yslope

                        particles[i][j] = hitdata;

                        hits[i] = {xcoord, ycoord, xslope, yslope};

                }

        }

        // save_vector("2017_hit_data",hits);
}

/* Propagates particles through the experiment. All radiation lengths are taken from PDG  */
void simulator::propagate_particles(void) {

        generate_beam();
        int total_emitted = 0, total_detected = 0;

        double d_f = 200.0; // thickness of Tantalum foil (micro-meter)

        double start_time = omp_get_wtime();
        double T = 0;
        int prog = 0;

        #pragma omp parallel for
        for (int i = 0; i < Nevents; i++) {

                default_random_engine private_generator;
                private_generator.seed(random_device{} ());


                int emitted = 0;

                if (bg_include) add_photons(i, emitted, X0_Si_amorph, 3.0E+03,1,private_generator); // add photons to match measured background spectrum

                /* MIMOSA 1 detector */
                if (bg_include) amorph_material(i, emitted, X0_tape, 50.0, 50.0,1,private_generator);
                SA_mult_scat(i, X0_tape, 50.0,private_generator);
                if (bg_include) amorph_material(i, emitted, X0_Mimosa, 100.0, 50.0,1,private_generator);
                SA_mult_scat(i, X0_Mimosa, 100.0,private_generator);
                mimosa_detector(0, i, total_detected,private_generator);

                /* Helium between M1 and M2 */
                if (bg_include) amorph_material(i, emitted, X0_He, zplanes[1] - 50.0, zplanes[1] - 150.0,1,private_generator);
                SA_mult_scat(i,  X0_He, zplanes[1] - 50.0,private_generator);

                /* MIMOSA 2 detector */
                if (bg_include) amorph_material(i, emitted, X0_tape, zplanes[1], 50.0,1,private_generator);
                SA_mult_scat(i, X0_tape, zplanes[1],private_generator);
                if (bg_include) amorph_material(i, emitted, X0_Mimosa, zplanes[1] + 50.0, 50.0,1,private_generator);
                SA_mult_scat(i, X0_Mimosa, zplanes[1] + 50.0,private_generator);
                mimosa_detector(1, i, total_detected,private_generator);

                /* Last mylar window He encasing */
                if (bg_include) amorph_material(i, emitted, X0_Mylar, zplanes[1] + 100.0, 50.0,1,private_generator);
                SA_mult_scat(i, X0_Mylar, zplanes[1] + 100.0,private_generator);

                /* Air */
                if (bg_include) amorph_material(i, emitted, X0_air, 2060E+03, 2060E+03 - (zplanes[1] + 100.0),1,private_generator);
                SA_mult_scat(i, X0_air, 2060E+03,private_generator);

                /* Traverse crystal */
                if (!strcmp("amorphous", name)) { // if amorphous

                        amorph_material(i, emitted, X0_C_amorph, 2060E+03 + d_c, d_c, 300,private_generator);
                        SA_mult_scat(i,  X0_C_amorph, 2060E+03 + d_c,private_generator);

                } else if (!strcmp("aligned", name)) { // if aligned

                        aligned_crystal(i, emitted, d_c/1000.0 * 300,private_generator);
                        // SA_mult_scat(i, X0_C_gem, 2060E+03);

                }

                /* Air */
                if (bg_include) amorph_material(i, emitted, X0_air, 2310E+03, 2310E+03 - (2060E+03 + d_c),1,private_generator);
                SA_mult_scat(i, X0_air, 2310.0E+03,private_generator);

                /* Traverse Scintilators */
                if (bg_include) amorph_material(i, emitted, X0_Si_amorph, 2310E+03 + 1.0E+03, 1.0E+03, 300,private_generator);
                SA_mult_scat(i, X0_Si_amorph, 2310E+03 + 1.0E+03,private_generator);

                /* Air between S2 and vacuum tube */
                if (bg_include) amorph_material(i, emitted, X0_air, 2310E+03 + 1.0E+03 + 0.5E+06, 0.5E+06,1,private_generator);
                SA_mult_scat(i, X0_air,  2310E+03 + 1.0E+03 + 0.5E+06,private_generator);

                /* 1st vacuum tube window */
                if (bg_include) amorph_material(i, emitted, X0_Mylar, 2310E+03 + 1.0E+03  + 0.5E+06 + 120.0, 120.0,1,private_generator);
                SA_mult_scat(i, X0_Mylar, 2310E+03 + 1.0E+03  + 0.5E+06 + 120.0,private_generator);
                if (bg_include) amorph_material(i, emitted, X0_tape, 2310E+03 + 1.0E+03  + 0.5E+06 + 120.0 + 100.0, 100.0,1,private_generator);
                SA_mult_scat(i, X0_tape,  2310E+03 + 1.0E+03  + 0.5E+06 + 120.0 + 100.0,private_generator);

                /* If not alignment run */
                if (strcmp("alignmentrun", name) != 0) {

                        MBPL_magnet(i);

                } else {

                        /* Air intead of vacuum? */
                        if (bg_include) amorph_material(i, emitted, X0_air*1e14, 2.8112e+006 + 5.14e+006, 5.14e+006,1,private_generator);
                        SA_mult_scat(i,  X0_air*1e14, 2.8112e+006 + 5.14e+006,private_generator);

                        /* 2nd vacuum tube window */
                        if (bg_include) amorph_material(i, emitted, X0_Mylar, 2.8112e+006 + 5.14e+006 + 120.0, 120.0,1,private_generator);
                        SA_mult_scat(i, X0_Mylar, 2.8112e+006 + 5.14e+006 + 120.0,private_generator);
                        if (bg_include) amorph_material(i, emitted, X0_tape, 2.8112e+006 + 5.14e+006 + 120.0 + 100.0, 100.0,1,private_generator);
                        SA_mult_scat(i, X0_tape, 2.8112e+006 + 5.14e+006 + 120.0 + 100.0,private_generator);

                        /* Air between vacuum tube and M3 */
                        if (bg_include) amorph_material(i, emitted, X0_air, zplanes[2] - d_f, zplanes[2] - d_f - ( 2.8112e+006 + 5.14e+006 + 120.0 + 100.0),1,private_generator);
                        SA_mult_scat(i, X0_air, zplanes[2] - d_f, private_generator);

                }

                /* Traverse converter foil */
                project_particle(photons, zplanes[2] - d_f, i);
                converter_foil(i, d_f, X0_Ta, 300,private_generator);
                SA_mult_scat(i, X0_Ta, zplanes[2],private_generator);

                /* MIMOSA 3 detector */
                if (bg_include) amorph_material(i, emitted, X0_tape, zplanes[2] + 50.0, 50.0,1,private_generator);
                SA_mult_scat(i, X0_tape, zplanes[2] + 50.0,private_generator);
                if (bg_include) amorph_material(i, emitted, X0_Mimosa, zplanes[2] + 100.0, 50.0,1,private_generator);
                SA_mult_scat(i, X0_Mimosa, zplanes[2] + 100.0,private_generator);
                mimosa_detector(2, i, total_detected,private_generator);

                /* Air between M3 and M4 */
                if (bg_include) amorph_material(i, emitted, X0_air, zplanes[3], zplanes[3] - zplanes[2] - 100.0,1,private_generator);
                SA_mult_scat(i, X0_air, zplanes[3],private_generator);

                /* MIMOSA 4 detector */
                if (bg_include) amorph_material(i, emitted, X0_tape, zplanes[3] + 50.0, 50.0,1,private_generator);
                SA_mult_scat(i, X0_tape, zplanes[3] + 50.0,private_generator);
                if (bg_include) amorph_material(i, emitted, X0_Mimosa, zplanes[3] + 100.0, 50.0,1,private_generator);
                SA_mult_scat(i, X0_Mimosa, zplanes[3] + 100.0,private_generator);
                mimosa_detector(3, i, total_detected,private_generator);

                /* If not alignment run */
                if (strcmp("alignmentrun", name) != 0) {

                        /* Air between M4 and middle of MIMOSA magnet */
                        if (bg_include) amorph_material(i, emitted, X0_air, (zplanes[4] + zplanes[3])/2.0, (zplanes[4] - zplanes[3] - 100.0)/2.0,1,private_generator);
                        SA_mult_scat(i, X0_air, (zplanes[4] + zplanes[3])/2.0,private_generator);

                        /* MIMOSA magnet */
                        mimosa_magnet(i);

                        /* Air between middle of MIMOSA magnet and M5 */
                        if (bg_include) amorph_material(i, emitted, X0_air, zplanes[4], (zplanes[4] - zplanes[3])/2.0,1,private_generator);
                        SA_mult_scat(i, X0_air, zplanes[4],private_generator);

                } else {

                        if (bg_include) amorph_material(i, emitted, X0_air, zplanes[4], zplanes[4] - zplanes[3] - 100.0,1,private_generator);
                        SA_mult_scat(i, X0_air, zplanes[4],private_generator);

                }

                /* MIMOSA 5 detector */
                if (bg_include) amorph_material(i, emitted, X0_tape, zplanes[4] + 50.0, 50.0,1,private_generator);
                SA_mult_scat(i, X0_tape, zplanes[4] + 50.0,private_generator);
                if (bg_include) amorph_material(i, emitted, X0_Mimosa, zplanes[4] + 100.0, 50.0,1,private_generator);
                SA_mult_scat(i, X0_Mimosa, zplanes[4] + 100.0,private_generator);
                mimosa_detector(4, i, total_detected,private_generator);

                /* Air between MIMOSA 5 and MIMOSA 6 detectors */
                if (bg_include) amorph_material(i, emitted, X0_air, zplanes[5], zplanes[5] - zplanes[4] - 100.0,1,private_generator);
                SA_mult_scat(i, X0_air, zplanes[5], private_generator);

                /* MIMOSA 6 detector */
                if (bg_include) amorph_material(i, emitted, X0_tape, zplanes[5] + 50.0, 50.0,1,private_generator);
                SA_mult_scat(i, X0_tape, zplanes[5] + 50.0,private_generator);
                if (bg_include) amorph_material(i, emitted, X0_Mimosa, zplanes[5] + 100.0, 50.0,1,private_generator);
                SA_mult_scat(i, X0_Mimosa, zplanes[5] + 100.0,private_generator);
                mimosa_detector(5, i, total_detected,private_generator);

                #pragma omp critical
                {
                        total_emitted += emitted;
                        prog++;
                }
                /* Report progress in terminal */
                if ( (prog+1) % (Nevents/10 + 1) == 0) {

                        double dt = omp_get_wtime() - start_time;
                        T += dt;
                        cerr << "Progress :\t" << floor(100 * double(prog)/(double)Nevents) << "%" << "\ttime used :\t" << dt << "\ttotal time elapsed :\t" << T << "\ttime remaining :\t" << dt * (double)Nevents/(Nevents/10 + 1) - T << "\n";
                        start_time = omp_get_wtime();

                }

        }

        /* Report results to terminal */
        cerr << "\nTotal photons emitted : " << total_emitted << "\n";
        cerr << "Total photons emitted by crystal : " << emitted_crystal << "\n";
        cerr << "Signal-to-Noise Ratio (SNR) : " << double(emitted_crystal) / double(total_emitted - emitted_crystal) << "\n";
        cerr << "Total photons incident on foil : " << no_photons << "\n";
        cerr << "Total conversions : " << conversions << "\n";

}

/* Simulate a particle travelling through an amorphous crystal. */
void simulator::add_photons(int eventno, int &emitted, double X0, double d, int no_slices, default_random_engine generator) {

        normal_distribution<double> distribution(0.0, 1.0);
        double dl = d/(double)no_slices;
        double Emin = 2*0.5109989461E-03; // 2 * electron mass in GeV

        for (size_t hitno = 0; hitno < particles[eventno].size(); hitno++) {

                for (int i = 0; i < no_slices; i++) {

                        double Epart = particles[eventno][hitno][5];

                        if (Epart < Emin) break;

                        /* Determine if photon is emitted */
                        bool emission = photon_emitted_amorph(dl, X0, Epart,generator);

                        /* If photon is emitted update "photons" array and particles array */
                        if (emission) {
                                emitted++;

                                /* Determine photon energy */
                                double randno = R(generator);
                                Epart = particles[eventno][hitno][5];
                                double norm = 4.0/3.0 * log(Epart/Emin) - 4.0/(3.0*Epart) * (Epart - Emin) + 1.0/(2.0*Epart*Epart) * (Epart*Epart - Emin*Emin);
                                function<double(vec)> energy =  [randno, Epart, norm] (vec x) {
                                                                       return photonic_energy_distribution(x, randno, Epart, norm);
                                                               };                                                                      // make lambda-function in order to use same randno during iteration
                                vec sc1 = {0.1}; vec sc2 = {20.0}; vector<vec> initial_simplex = {sc1, sc2}; // initial simplex for Nelder-Mead. The initial guess is hugely important for convergence
                                vec photon_energy = simplex_NM(energy, initial_simplex, 1.0E-08); // solve for energy using Nelder-Mead simplex.

                                /* Determine direction of photon */
                                double gamma = particles[eventno][hitno][5]/(5.109989461E-4);
                                uniform_real_distribution<double> defl_angle(-1.0/gamma, 1.0/gamma);

                                /* Update particles and photons vectors */
                                vector<double> photon = particles[eventno][hitno];
                                photon[4] = 0.0; // charge
                                photon[5] = photon_energy(0);
                                double dx = defl_angle(generator);
                                double dy = sqrt(1.0/(gamma*gamma) - dx*dx);
                                photon[6] += dx;
                                photon[7] += dy;
                                photons[eventno].push_back(photon);

                        }

                }

        }

}

void simulator::amorph_material(int eventno, int &emitted, double X0, double z, double L, int N_slices, default_random_engine generator) {
        double Emin = 2.0*0.5109989461E-03;

        double dl = L/(double)N_slices;
        vector<double> dy_part = {0.0, 0.0};
        vector<double> dtheta_part = {0.0, 0.0};
        // bool inside;

        for (size_t hitno = 0; hitno < particles[eventno].size(); hitno++) {

                double L_frac = (z - particles[eventno][hitno][2])/L;
                int no_slices = L_frac * N_slices;

                for (int i = 0; i < no_slices; i++) {

                        double Epart = particles[eventno][hitno][5];

                        /* Proceed only if photon is emitted */
                        bool emission = photon_emitted_amorph(dl, X0, Epart,generator);

                        if (emission and Epart > Emin) {

                                emitted++;
                                if (X0 == X0_C_amorph or X0 == X0_C_gem) emitted_crystal++;

                                /* Determine photon energy */
                                double randno = R(generator);
                                double norm = 4.0/3.0 * log(Epart/Emin) - 4.0/(3.0*Epart) * (Epart - Emin) + 1.0/(2.0*Epart*Epart) * (Epart*Epart - Emin*Emin);

                                vec sc1 = {0.001}, sc2 = {40.0};
                                vector<vec> initial_simplex = {sc1, sc2};

                                vec photon_energy = simplex_NM([randno, Epart, norm, Emin] (vec x) {
                                        return x(0) > 0 ? abs((4.0/3.0 * log(x(0)/Emin) - 4.0/3.0 * (x(0) - Emin)/Epart + 1.0/(2.0*Epart*Epart) * (x(0)*x(0) - Emin*Emin)) - randno*norm) : 1E+17;
                                }, initial_simplex, 1.0E-08);

                                /* Determine direction of photon */
                                double gamma = particles[eventno][hitno][5]/(5.109989461E-4);
                                uniform_real_distribution<double> defl_angle(-1.0/gamma, 1.0/gamma);

                                /* Update particles and photons vectors */
                                vector<double> photon = particles[eventno][hitno];
                                photon[0] += dy_part[0];
                                photon[1] += dy_part[1];
                                photon[2] += (1 + i) * dl;           // z-val
                                photon[4] = 0.0;           // charge
                                photon[5] = photon_energy(0);
                                double dx = defl_angle(generator);
                                double dy = sqrt(1.0/(gamma*gamma) - dx*dx);
                                photon[6] += dx;
                                photon[7] += dy;
                                particles[eventno][hitno][5] -= photon_energy(0);
                                photons[eventno].push_back(photon);

                        }

                }           // slices end

        } // events end

}

void simulator::converter_foil(int eventno, double d_f, double X0, int no_slices, default_random_engine generator) {

        double Emin = 2*0.5109989461E-03; // 2 * electron mass in GeV
        double dl = d_f/(double)no_slices;
        uniform_real_distribution<double> rotation(0.0, 2.0*M_PI);

        for (size_t i = 0; i < photons[eventno].size(); i++) {
                no_photons += 1;

                for (int j = 0; j < no_slices; j++) {

                        /* Proceed only if a conversion happens */
                        double randno = R(generator);
                        bool multiple_conversion = false;

                        if (randno < d_f/(double)no_slices * (7.0)/(9.0*X0) and photons[eventno][i][5] > Emin) {

                                if (!multiple_conversion) {

                                        conversions++;
                                        multiple_conversion = true;

                                }

                                /* Calculate energy gained by e+/e- pair */
                                double photon_energy = photons[eventno][i][5];
                                double electron_energy = R(generator); // random number for the inverse transform sampling in "electronic_energy_distribution"
                                electronic_energy_distribution(electron_energy); // the fractional electron energy, ie E_e-/ E_photon
                                electron_energy *= photon_energy;
                                double positron_energy = photon_energy - electron_energy; // energy conservation

                                /* Calculate deflection of e+/e- pair */
                                double electron_defl;
                                double positron_defl;
                                Borsellino(electron_energy, positron_energy, photon_energy, electron_defl, positron_defl,generator); // deflection angle based on approximated Borsellino distribution

                                double r1 = rotation(generator);

                                /* Add e+/e- pair to "particles" array */
                                vector<double> electron = photons[eventno][i];
                                electron[2] += dl * (j + 1);
                                electron[4] = (-1.0)*q;
                                electron[5] = electron_energy;
                                electron[6] += cos(r1)*electron_defl;
                                electron[7] += sin(r1)*electron_defl;

                                vector<double> positron = photons[eventno][i];
                                positron[2] += dl * (j + 1);
                                positron[4] = q;
                                positron[5] = positron_energy;
                                positron[6] += cos(r1)*positron_defl;
                                positron[7] += sin(r1)*positron_defl;

                                particles[eventno].push_back(positron);
                                particles[eventno].push_back(electron);

                                photons[eventno].erase(photons[eventno].begin() + i); // remove i'th photon since it converted

                                i -= 1; // subtract 1 from iterator since we removed element
                                break;

                        }

                }

        }

}

void simulator::MBPL_magnet(int eventno) {

        particles[eventno].clear(); // the MPBL magnet removes all particles

}

void simulator::make_intensity_distro(vector<vector<double> > &intensities, vector<vector<double> > &angles) {

        /* Load data of spatial distribution of energy into matrix */
        string file = initial_spec;
        mat intensity_angles; intensity_angles.load(file);

        /* Obtain emission angles */
        angles.resize(intensity_angles.n_cols);

        for (size_t i = 0; i < intensity_angles.n_cols; i++) {

                angles[i].resize(2);
                angles[i][0] = intensity_angles(3, i);
                angles[i][1] = intensity_angles(4, i);

        }

        /* Define vector with intensities */
        intensity_angles.shed_rows(0, 4);
        intensities.resize(energies_interp.size()); // contains angle/intensity spectrum for every energy
        for (size_t i = 0; i < intensities.size(); i++) intensities[i].resize(intensity_angles.n_cols);

        vector<double> intensities_i(intensity_angles.n_rows);

        /* Calculate spatial intensity distro. for every energy */
        for (size_t j = 0; j < energies_interp.size(); j++) {

                for (size_t i = 0; i < intensity_angles.n_cols; i++) {

                        /* Extract column i from intensity_spatial matrix */
                        for (size_t k = 0; k < intensity_angles.n_rows; k++) {

                                intensities_i[k] = intensity_angles(k, i); // angle/intensity spectrum for direction k

                        }

                        /* Do linear interpolation to increase energy/intensity resoultion */
                        linterp(emitted_energies, intensities_i, energies_interp, intensities_i);

                        /* Fill "intensities" vector. The j'th entry corresponds to given energy, the i'th entry is the intensity at the angle corresponding to the i'th column in "sum_initials1mm40GeVelec" vector */
                        intensities[j][i] = intensities_i[j]; // j : energy, i : angle/intensity spectrum for given direction

                }

        }

        // save_vector("intensity_distro", intensities);

}

bool simulator::photon_emitted_aligned(double no_slices, default_random_engine generator) {

        double randno = R(generator);

        return randno < I_integral/(double)no_slices;

}

/* Simulate a particle travelling through an aligned crystal. Multiple scattering is taken into account in the "propagte particles" function. */
void simulator::aligned_crystal(int eventno, int &emitted, int no_slices, default_random_engine generator) {

        double Emin = 2*0.5109989461E-03; // 2 * electron mass in GeV

        double dl = d_c/(double)no_slices;

        for (size_t hitno = 0; hitno < particles[eventno].size(); hitno++) {

                for (int i = 0; i < no_slices; i++) {

                        double Epart = particles[eventno][hitno][5];

                        /* Multiple scattering in C crystal slice. Calculate x/y displacement and deflection independently */
                        for (int k = 0; k < 2; k++) {

                                double z1 = R_normal(generator);
                                double z2 = R_normal(generator);
                                double theta0 = 0.0136/Epart * sqrt(dl/X0_C_gem) * (1.0 + 0.038 * log(dl/X0_C_gem));
                                double dy = dl * theta0 * (z1/sqrt(12) +  z2/2.0);
                                double dtheta0 = z2 * theta0;
                                particles[eventno][hitno][k] += dy + dl*particles[eventno][hitno][6+k];         // total (x/y) displacement
                                particles[eventno][hitno][6+k] += dtheta0;         // xslope/yslope

                        }

                        particles[eventno][hitno][2] += dl;

                        /* Proceed only if photon is emitted */
                        bool emission = photon_emitted_aligned(no_slices,generator);

                        if (emission and Epart > Emin) {

                                emitted++;
                                emitted_crystal++;

                                /* Determine which energy photon is emitted with */
                                int E_indx = select_member(intensity_sum_interp);

                                /* Determine where in space photon is emitted */
                                double gamma = particles[eventno][hitno][5]/(5.109989461E-4);
                                uniform_real_distribution<double> defl_angle(-1.0/gamma, 1.0/gamma);

                                double dx = defl_angle(generator);
                                double dy = sqrt(1.0/(gamma*gamma) - dx*dx);

                                /* Add photon to "photons" vector */
                                vector<double> photon = particles[eventno][hitno];
                                photon[4] = 0.0;         // charge
                                photon[5] = energies_interp[E_indx];
                                photon[6] += dx;
                                photon[7] += dy;
                                particles[eventno][hitno][5] -= energies_interp[E_indx];         // energy lost by particle as it emits photon
                                photons[eventno].push_back(photon);

                        }

                }         // end slices

        } // end events

}

void simulator::mimosa_detector(int planeno, int eventno, int &detections, default_random_engine generator) {

        /* Define mimosa parameters */
        double x1 = 0.0, x2 = 0.0, x3 = 0.0, x4 = 0.0, y1 = 0.0, y2 = 0.0, y3 = 0.0, y4 = 0.0, fakehitprob = 0;

        switch (planeno) {

        case 0:
                x1 = -1.048E+04;
                x2 = -1.052E+04;
                x3 = 1.053E+04;
                x4 = 1.059E+04;

                y1 = -0.5097E+04;
                y2 = 0.5290E+04;
                y3 = 0.5272E+04;
                y4 = -0.5235E+04;

                fakehitprob = 1;

                break;

        case 1:
                x1 = -1.048E+04;
                x2 = -1.052E+04;
                x3 = 1.053E+04;
                x4 = 1.059E+04;

                y1 = -0.5097E+04;
                y2 = 0.5290E+04;
                y3 = 0.5272E+04;
                y4 = -0.5235E+04;

                fakehitprob = 0.95;

                break;

        case 2:
                x1 = -1.154E+04;
                x2 =  -1.162E+04;
                x3 = 0.9258E+04;
                x4 = 0.9356E+04;

                y1 = -0.6043E+04;
                y2 = 0.4317E+04;
                y3 = 0.4742E+04;
                y4 = -0.5635E+04;

                fakehitprob = 0.43;

                break;

        case 3:
                x1 = -1.129E+04;
                x2 = -1.139E+04;
                x3 = 0.9605E+04;
                x4 = 0.9711E+04;

                y1 = -0.6198E+04;
                y2 = 0.4216E+04;
                y3 = 0.4659E+04;
                y4 = -0.5742E+04;

                fakehitprob = 0.43;

                break;

        case 4:
                x1 = -1.177E+04;
                x2 = -1.187E+04;
                x3 = 0.9459E+04;
                x4 = 0.9579E+04;

                y1 = -0.5652E+04;
                y2 = 0.4868E+04;
                y3 = 0.5380E+04;
                y4 = -0.5121E+04;

                fakehitprob = 0.80;

                break;

        case 5:
                x1 = -1.288E+04;
                x2 = -1.29E+04;
                x3 = 0.8563E+04;
                x4 = 0.8601E+04;

                y1 = -0.5410E+04;
                y2 = 0.5137E+04;
                y3 = 0.5446E+04;
                y4 = -0.5064E+04;

                fakehitprob =0.52;

                break;
        }

        vector<double> xbounds = {x1, x2, x3, x4}; // xcoords of corners
        vector<double> ybounds = {y1, y2, y3, y4}; // ycoords of corners

        mimosa_res = 35.0; // spatial resoultion of mimosa detector (micro-m) (smaller = better res.)

        /* Take into account detector resoultion by adding a random (Dx, Dy) displacement between -4 -> 4 micro-meter */
        normal_distribution<double> distribution(0.0, 3.5);

        /* Simulate real hit detection */
        for (size_t hitno = 0; hitno < particles[eventno].size(); hitno++) {

                double Dx = distribution(generator);
                double Dy = distribution(generator);

                /* Check if particle is within physical boundaries of detector.*/
                bool inside_bounds = isInside(4, xbounds, ybounds, particles[eventno][hitno][0], particles[eventno][hitno][1]);

                if (inside_bounds) {

                        /* Proceed only if particle does not hit hot pixel */
                        int pixel = coord2pixel(particles[eventno][hitno][0], particles[eventno][hitno][1]);
                        int high = hotpixels[planeno].size() - 1;
                        int low = 0;
                        int hotpix = binarySearch(hotpixels[planeno], low, high, pixel);

                        if (hotpix == -1) {

                                /* Combine neighbor-hits into single hit if they are closer than detectors resoultion */
                                if (mimosas[planeno][eventno].size() > 0) {

                                        double dist = calc_dist(particles[eventno][hitno][0] + Dx, particles[eventno][hitno][1] + Dy, mimosas[planeno][eventno][0][0], mimosas[planeno][eventno][0][1]);

                                        /* Save hit in "mimosas" array. This array has the same structure and contains equivalent data as the "hitcoords" array in the "analyser" class */
                                        if (dist > mimosa_res) {

                                                mimosas[planeno][eventno].push_back({particles[eventno][hitno][0] + 0*Dx, particles[eventno][hitno][1] + 0*Dy});

                                        } else {
                                                /* The distance is shorter than detector resolution, so we combine two hits. This is equivalent to moving the previously recorded hit by half the seperation */
                                                double dx = mimosas[planeno][eventno][0][0] - particles[eventno][hitno][0] - Dx;
                                                double dy = mimosas[planeno][eventno][0][1] - particles[eventno][hitno][1] - Dy;
                                                mimosas[planeno][eventno][0][0] += dx/2.0;
                                                mimosas[planeno][eventno][0][1] += dy/2.0;

                                        }

                                } else {

                                        /* If there are no previous hits in detector simply add the hit to "mimoas" array */
                                        mimosas[planeno][eventno].push_back({particles[eventno][hitno][0] + Dx, particles[eventno][hitno][1] + Dy});

                                }


                        }

                }

        }

        /* Simulate fake hit detection */
        uniform_real_distribution<double> fake_hit_x(x1, x4);
        uniform_real_distribution<double> fake_hit_y(y1, y2);

        if (R(generator) < fakehitprob) {

                double xcoord = fake_hit_x(generator);
                double ycoord = fake_hit_y(generator);

                vector<double> hit = {xcoord, ycoord};

                mimosas[planeno][eventno].push_back(hit);

        }

}

double simulator::calc_dist(double x0, double y0, double x1,double y1) {

        return sqrt(pow(x1 - x0, 2) + pow(y1 - y0, 2));

}

void simulator::project_particle(vector<vector<vector<double> > > &particletype, double zcoord, int eventno) {

        for (size_t hitno = 0; hitno < particletype[eventno].size(); hitno++) {

                double dx = (zcoord - particletype[eventno][hitno][2]) *  particletype[eventno][hitno][6];
                double dy = (zcoord - particletype[eventno][hitno][2]) *  particletype[eventno][hitno][7];
                particletype[eventno][hitno][0] += dx;
                particletype[eventno][hitno][1] += dy;
                particletype[eventno][hitno][2] = zcoord;

        }

}

/* (i,j) : (eventno, hitno). The following formulae taken from PDG ch. 32, 2014 (particles through matter) */
void simulator::SA_mult_scat(int i, double X0, double zcoord, default_random_engine generator) {

        /* Iterate over paticles traversing medium */
        for (size_t j = 0; j < particles[i].size(); j++) {

                double l = zcoord - particles[i][j][2];

                /* Calculate x/y displacement and deflection independently */
                for (int k = 0; k < 2; k++) {

                        double z1 = R_normal(generator);
                        double z2 = R_normal(generator);
                        double theta0 = 0.0136/particles[i][j][5] * sqrt(l/X0) * (1.0 + 0.038 * log(l/X0));
                        double dy = l * theta0 * (z1/sqrt(12) +  z2/2.0);
                        double dtheta0 = z2 * theta0;
                        particles[i][j][k] += dy + l*particles[i][j][6+k]; // total (x/y) displacement
                        particles[i][j][6+k] += dtheta0; // xslope/yslope

                }

                particles[i][j][2] = zcoord;

        }

}

/* Calculate the deflection of a charged particle traversing the mimosa magnet. Magnet only deflects in x-direction. */
void simulator::mimosa_magnet(int eventno) {

        double L = 0.15; // length traveled through Mimosa magnet in m
        double B = 0.12; // strength of magnetic field in T

        for (size_t hitno = 0; hitno < particles[eventno].size(); hitno++) {

                particles[eventno][hitno][6] += (particles[eventno][hitno][4]*L*B*c)/(particles[eventno][hitno][5] * 1.6021766E-10); // update direction of particle

        }

}

/* Approximated (SAA) Borsellino opening angle of produced electron (1) positron (2) pair */
void simulator::Borsellino(double E1, double E2, double E_phot, double &phi1, double &phi2, default_random_engine generator) {

        double u1 = R(generator), u2 = R(generator);
        double nu0 = E_phot*m*c/(E1*E2) * 1.6021766E+10;
        double nu = nu0*(sqrt(1.0/(1.0 - u1) - 1.0) + 0.70*u2);
        phi1 = nu/(E1/E2 + (1.0 - nu*nu/2.0)); // deflection of particle 1
        phi2 = nu - phi1; // deflection of particle 2

}

void simulator::pair_production(int eventno, double d_f, double X0, int N_slices, default_random_engine generator) {

        double dl = d_f/(double)N_slices;

        for (size_t i = 0; i < photons[eventno].size(); i++) {

                for (int j = 0; j < N_slices; j++) {

                        /* Proceed only if a conversion happens */
                        double randno = R(generator);

                        if (randno < d_f/(double)N_slices * (7.0)/(9.0*X0) and photons[eventno][i][5] < 2*0.5109989461E-03) {

                                /* Calculate energy gained by e+/e- pair */
                                double photon_energy = photons[eventno][i][5];
                                // energies.push_back(photon_energy);
                                double electron_energy = R(generator); // random number for the inverse transform sampling in "electronic_energy_distribution"
                                electronic_energy_distribution(electron_energy); // the fractional electron energy, ie E_e-/ E_photon
                                electron_energy *= photon_energy;
                                double positron_energy = photon_energy - electron_energy; // energy conservation

                                /* Calculate deflection of e+/e- pair */
                                double electron_defl;
                                double positron_defl;
                                Borsellino(electron_energy, positron_energy, photon_energy, electron_defl, positron_defl,generator); // deflection angle based on approximated Borsellino distribution

                                /* Add e+/e- pair to "particles" array */
                                vector<double> electron = photons[eventno][i];
                                electron[2] += dl * (j + 1);
                                electron[4] = (-1.0)*q;
                                electron[5] = electron_energy;
                                electron[6] += electron_defl;

                                vector<double> positron = photons[eventno][i];
                                positron[2] += dl * (j + 1);
                                positron[4] = q;
                                positron[5] = positron_energy;
                                positron[6] += positron_defl;

                                particles[eventno].push_back(positron);
                                particles[eventno].push_back(electron);

                                /* Remove photon since it has converted */
                                photons[eventno].erase(photons[eventno].begin() + i);

                                i -= 1;

                                break;

                        }

                }

        }

}

/* See http://pdg.lbl.gov/2014/reviews/rpp2014-rev-passage-particles-matter.pdf for formulae and constants */
bool simulator::pair_produced(double l, double X0_f, default_random_engine generator) {

        double randno = R(generator); // random double between 0, 1 used to determine if a pair is produced

        return randno < (7.0*l)/(9.0*X0_f);

}

/* Calcualte wheter or not a photon is emitted when traveling through an amorphous crystal. */
bool simulator::photon_emitted_amorph(double l, double X0, double E, default_random_engine generator) {

        double randno = R(generator); // random double between 0, 1 used to determine if a pair is produced
        double Emin = 2*0.5109989461E-03; // 2 * Electron mass in GeV must be the minimum energy of photon, since we cannot observe lower energy photons anyway
        double P = l/X0 * ( (4.0/3.0) * log(E/Emin) - (4.0/3.0) * (E - Emin)/E  + (3.0/4.0) * (1.0 - Emin*Emin/(2.0*E*E)));// BH-cross section * crystal length

        return randno < P;

}

void simulator::electronic_energy_distribution(double &r) {

        /* Analytical solution to CDF(x) = r (see Flohr's thesis) integrated and solved 3rd order poly. equation (2.58)
            r : random number */
        r = 0.25*(1.5874*pow(sqrt(196.0*r*r - 196.0*r + 81.0) + 14*r - 7.0, 1.0/3.0) - 5.03968/pow(sqrt(196.0*r*r - 196.0*r + 81.0) + 14.0*r - 7.0, 1.0/3.0) + 2.0);

}

double simulator::photonic_energy_distribution(vec x, double r, double Eb, double norm) {

        /* This function is minimized to determine the photon energy.
           The Bethe-Heitler cross section as presented in Flohr's thesis, should be solved using inverse transform sampling with a Nelder-Mead simplex solver. The equation is therefore normalized to ensure the probability distro. is betwen (0, 1)
            x   : solution, ie E_e- / E_photon
            r   : random number betwenn (0, 1)
            Eb  : beam energy
         */

        double Emin = 2.0 * 0.5109989461E-3; // 2 * electron mass GeV
        return x(0) > 0 ? abs((4.0/3.0 * log(x(0)/Emin) - 4.0/3.0 * (x(0) - Emin)/Eb + 1.0/(2.0*Eb*Eb) * (x(0)*x(0) - Emin*Emin)) - r*norm) : 1E+17;


}

int simulator::binarySearch(vector<int> numbers, int low, int high, int val) {

        if (high >= low) {

                int mid = (high + low)/2;
                if (numbers[mid] == val) return mid;
                if (numbers[mid] > val) return binarySearch(numbers, low, mid-1, val);
                return binarySearch(numbers, mid + 1, high, val);

        }

        return -1; // returned only if "val" is not in "numbers"

}



/* This function takes a test point (testx, testy) and checks if it is within an nvert-sided polygon defined by the vertices vertx and verty. Function is a slightly modified version of the one posted on https://wrf.ecse.rpi.edu//Research/Short_Notes/pnpoly.html */
int simulator::isInside(int nvert, vector<double> vertx, vector<double> verty, double testx, double testy) {

        int i, j, c = 0;

        for (i = 0, j = nvert-1; i < nvert; j = i++) {

                if ( ((verty[i]>testy) != (verty[j]>testy)) && (testx < (vertx[j]-vertx[i]) * (testy-verty[i]) / (verty[j]-verty[i]) + vertx[i]) ) {

                        c = !c;

                }

        }

        return c;

}

/* This function selects a member from a weighted ("weights") list and returns the indx */
int simulator::select_member(vector<double> weights) {

        /* Calculate cumulative sum of weights to represent the probability that a number will be picked */
        vector<double> weights_cumsum; weights_cumsum.resize(weights.size());
        partial_sum(weights.begin(), weights.end(), weights_cumsum.begin());
        uniform_real_distribution<double> distribution(0.0, weights_cumsum.back());

        /* Generate random number, used to determine index */
        double r = distribution(generator);

        /* Binary search to find which number to pick based on weight */
        int low = 0, high = weights.size() - 1;

        while (high >= low) {

                int guess = (low + high)/2;

                if (weights_cumsum[guess] < r) {

                        low = guess + 1;

                }

                else if (weights_cumsum[guess] - weights[guess] > r) {

                        high = guess - 1;

                }

                else {

                        return guess; // return the index for selected member in "numbers" array

                }

        }

        return -1; // if this is returned, somehting has gone wrong

}

vector<double> simulator::linspace(double min, double max, int N) {

        vector<double> range; range.resize(N);

        /* Determine step-size */
        double delta = (max-min)/(N-1);

        /* Iterate over range, excluding the last entry */
        for (int i = 0; i < N-1; i++) {

                range[i] = min + i * delta;

        }

        /* Set last entry to "max". This ensures that we get a range from min -> max */
        range.back() = max;
        return range;

}

void simulator::linterp(vector<double> x, vector<double> y, vector<double> xi, vector<double> &yi) {

        yi.resize(xi.size());

        //double dx = x[1] - x[0]; // x is monotonically increasing, so we can define dx outside loop. Otherwise calculate as dx = x[low+1] - x[low]

        for (size_t j = 0; j < xi.size(); j++) {

                /* Find the interval over which to interpolate with binary search */
                int low = 0;
                int high = x.size() - 1;

                while (high - low > 1) {

                        int guess = (high + low)/2;

                        if (xi[j] > x[guess]) {
                                low = guess;
                        } else high = guess;

                }

                /* Calculate slope and do interpolation */
                double dydx = (y[low + 1] - y[low])/(x[low + 1] - x[low]);
                yi[j] = y[low] + dydx * (xi[j] - x[low]);

        }

}

/* Draw from Alias table */
double simulator::VoseAliasMethod_draw(vector<int> Alias, vector<double> Prob) {

        /* Make fair roll of n-sided die */
        int n = Prob.size(); uniform_int_distribution<int> distribution(0, n);
        int i = distribution(generator);

        /* Draw double from uniform dist. on [0, 1] */
        double r = R(generator);

        /* Determine output based comparison between biasesd and fair coin-toss */
        return r < Prob[i] ? i : Alias[i];

}

/* Construct Alias table */
void simulator::VoseAliasMethod_table(vector<double> distro_in, vector<int> &Alias, vector<double> &Prob) {

        /* Allocate space for Alias and Prob arrays */
        vector<double> distro = distro_in;
        int n = distro.size();
        Alias.resize(n); Prob.resize(n);

        /* average probability, used for comparison with distro. value */
        double average = 1.0 / n;

        /* Fill "Small" and "Large" array using distro. */
        vector<double> Small, Large; // used to fill Alias and Prob

        for (int i = 0; i < n; i++) {

                if (distro[i] < average) {

                        Small.push_back(i);

                }
                else {

                        Large.push_back(i);

                }

        }

        /* Large may be emptied first due to innacuracy of doubles so use both criteria here. Mathematically Small should always be emptied first. */
        while (!Small.empty() and !Large.empty()) {

                /* Obtain index of small (l) and large (g) prob, used to populate Alias and Prob arrays */
                double l = Small.back(); Small.pop_back();
                double g = Large.back(); Large.pop_back();

                /* Populate the Prob and Alias arrays */
                Prob[l] = distro[l] * n;
                Alias[l] = g;
                distro[g] += distro[l] - average;

                /* Populate Large and Small array */
                if (distro[g] < average) {

                        Small.push_back(g);

                } else Large.push_back(g);

        } // end while

        /* Now empty both lists setting probabilities to 1.0. Ensures numerical stability */
        while (!Large.empty()) {

                double g = Large.back(); Large.pop_back();
                Prob[g] = 1.0;

        }
        while (!Small.empty()) {

                double l = Small.back(); Small.pop_back();
                Prob[l] = 1.0;

        }

}

/* Reflect highest point against centroid */
vec simulator::simplex_reflect(vec highest, vec centroid) {

        return 2.0 * centroid - highest;

}

/* Reflect highest point against centroid, then double distance to centroid */
vec simulator::simplex_expand(vec highest, vec centroid) {

        return 3.0 * centroid - 2.0 * highest;

}

/* Highest point halves its distance from centroid */
vec simulator::simplex_contract(vec highest, vec centroid) {

        return 0.5 * (centroid + highest);

}

/* All points, except lowest, halves their distance to lowest point */
void simulator::simplex_reduce(vector<vec> &simplex, int ilow) {

        for (size_t i = 0; i < simplex.size(); i++) {

                if ((int)i != ilow) {

                        simplex[i] = 0.5 * (simplex[i] + simplex[ilow]);

                }

        }

}

/* Calculate size of the simplex-polygon */
double simulator::simplex_size(vector<vec> simplex) {

        double s = norm(simplex[0] - simplex.back(), 2);

        for (size_t i = 1; i < simplex.size(); i++) {

                double d = norm(simplex[i] - simplex[i - 1], 2);

                if (s < d) s = d;

        }

        return s;

}

/* Update ihigh, ilow and centroid.  */
void simulator::simplex_update(vector<vec> simplex, vec fs, vec &centroid, int &ihigh, int &ilow) {

        ihigh = 0; double fhigh = fs(0);
        ilow = 0; double flow = fs(0);

        for (size_t i = 1; i < simplex.size(); i++) {

                if (fs(i) > fhigh) {

                        fhigh = fs(i); ihigh = i;

                }
                if (fs(i) < flow) {

                        flow = fs(i); ilow = i;

                }

        }

        vec s; s.resize(simplex[0].size());

        for (int i = 0; i < (int)simplex.size(); i++) {

                if (i != ihigh) {

                        s += simplex[i];

                }

        }

        for (int i = 0; i < (int)s.size(); i++) {

                centroid(i) = s(i) / (simplex.size() - 1);

        }

}

void simulator::save_vector(string name, vector<vector<double> > data) {

        string filename = DATPATH + "/" + name + ".txt";
        ofstream output (filename);

        for (size_t i = 0; i < data.size(); i++) {

                for (size_t j = 0; j < data[i].size(); j++) {

                        output << data[i][j] << "\t";

                }

                output << "\n";

        }

}
