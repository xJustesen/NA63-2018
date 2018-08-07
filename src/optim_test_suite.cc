void simualtor::amorph_material(int eventno, int &emitted, int no_slices, double X0, double d){
  vector<vector<double>> particledata = particles[eventno];
  for (size_t i = 0; i < particledata.size(); i++){
    amorph_material_helper(eventno, i, emitted, no_slices, X0, d);
  }
}

void simulator::amorph_material_helper(int eventno, int hitno, int &emitted, int no_slices, double X0, double d){
  for (int i = 0; i < no_slices; i++) {
    double Epart = particledata[hitno][5];
    if (Epart < 2.0*0.5109989461E-03){
      break;
    }

    /* Determine if photon is emitted */
    bool emission = photon_emitted_amorph(dl, X0, Epart);

    /* If photon is emitted update "photons" array and particles array */
    if (emission) {
      emitted++;

      /* Determine photon energy */
      double randno = R(generator);
      Epart = particledata[hitno][5];
      double norm = 4.0/3.0 * log(Epart/Emin) - 4.0/(3.0*Epart) * (Epart - Emin) + 1.0/(2.0*Epart*Epart) * (Epart*Epart - Emin*Emin);
      function<double(vec)> energy =  [randno, Epart, norm] (vec x) {return photonic_energy_distribution(x, randno, Epart, norm);};  // make lambda-function in order to use same randno during iteration
      vec sc1 = {0.1}; vec sc2 = {20.0}; vector<vec> initial_simplex = {sc1, sc2};  // initial simplex for Nelder-Mead. The initial guess is hugely important for convergence
      vec photon_energy = simplex_NM(energy, initial_simplex, 1.0E-09);  // solve for energy using Nelder-Mead simplex.
      energies.push_back(photon_energy(0));

      /* Determine direction of photon */
      double gamma = particledata[hitno][5]/(m*c*c);
      uniform_real_distribution<double> defl_angle(-1.0/gamma, 1.0/gamma);

      /* Update particles and photons vectors */
      vector<double> photon = particledata[hitno];  // LOCAL PHOTON VECTOR
      photon[4] = 0.0;  // charge
      photon[5] = photon_energy(0);
      double dx = defl_angle(generator);
      double dy = sqrt(1.0/(gamma*gamma) - dx*dx);
      photon[6] += dx;
      photon[7] += dy;
      particledata[hitno][5] -= photon_energy(0); // GLOBAL PARTILCES VECTOR
      photons[eventno].push_back(photon);

      /* Multiple scatter particle through slice */
      for (int j = 0; j < 2; j++){
        double z1 = distribution(generator);
        double z2 = distribution(generator);
        double z = particles[eventno][hitno][4]/q;
        double energy = particles[eventno][hitno][5] * 1.6021766E-10;  // energy in Joule
        double theta0 = 2.179E-12/(energy) * z * sqrt(l/X0) * (1 + 0.038*log(l/X0));
        double yplane = (z1*l*theta0)/sqrt(12) +  (z2*l*theta0)/2;
        double thetaplane = z2 * theta0;
        particles[eventno][hitno][j] += yplane + l * particles[eventno][hitno][6+j]; // GLOBAL particles vector
        particles[eventno][hitno][j+6] += thetaplane; // GLOBAL particles vector
      }

      /* Determine if a conversion happens */
      int indx = 0;
      while (indx < no_slices - i){
        bool conversion = pair_produced(d_f/(double)no_slices, X0);
        if (photons[eventno][i][5] < 2*0.5109989461E-03) conversion = false;
        if (conversion) {
          conversions++;
          /* Calculate energy gained by e+/e- pair */
          double electron_energy = R(generator);  // random number for the inverse transform sampling in "electronic_energy_distribution"
          electronic_energy_distribution(electron_energy); // the fractional electron energy, ie E_e-/ E_photon
          electron_energy *= photon_energy(0);
          double positron_energy = photon_energy(0) - electron_energy; // energy conservation
          photons[eventno][hitno][5] -= electron_energy + positron_energy;

          /* Calculate deflection of e+/e- pair */
          double electron_defl = Borsollini(electron_energy, positron_energy, photon_energy); // deflection angle based on approximated Borsollini distribution
          double positron_defl = electron_energy * electron_defl / positron_energy; // conservation of momentum

          /* Add e+/e- pair to "particles" array */
          vector<double> electron = photons[eventno][hitno];
          electron[4] = (-1.0)*q;
          electron[5] = electron_energy;
          electron[6] += electron_defl;

          vector<double> positron = photons[eventno][hitno];
          positron[4] = q;
          positron[5] = positron_energy;
          positron[6] += positron_defl;

          particles[eventno].push_back(positron); // GLOBAL particles vector
          particles[eventno].push_back(electron); // GLOBAL particles vector
          photons[eventno].erase(photons[eventno].begin() + hitno);

          /* Multiple scattering of electron + positron  through slice */
          double l = d/no_slices;
          int sz = particledata.size();
          normal_distribution<double> distribution(0.0, 1.0);

          for (int j = 1; j < 3; j++){
            for (int k = 0; k < 2; k++){
              double z1 = distribution(generator);
              double z2 = distribution(generator);
              double z = particles[eventno][sz - 1][4]/q;
              double energy = particles[eventno][sz - 1][5] * 1.6021766E-10;  // energy in Joule
              double theta0 = 2.179E-12/(energy) * z * sqrt(l/X0) * (1 + 0.038*log(l/X0));
              double yplane = (z1*l*theta0)/sqrt(12) +  (z2*l*theta0)/2;
              double thetaplane = z2 * theta0;
              particles[eventno][sz - j][k] += yplane + l * particles[eventno][sz - 1][6+k]; // GLOBAL particles vector
              particles[eventno][sz - j][k+6] += thetaplane; // GLOBAL particles vector
            }
          }
          amorph_material(eventno, sz - 1, emitted, no_slices - i, X0, ((double)no_slices - (double)i)/((double)no_slices) * d);
          amorph_material(eventno, sz - 2, emitted, no_slices - i, X0, ((double)no_slices - (double)i)/((double)no_slices) * d);
        }
      }
    }
    /* Multiple scatter particle through slice */
    for (int j = 0; j < 2; j++){
      double z1 = distribution(generator);
      double z2 = distribution(generator);
      double z = particles[eventno][hitno][4]/q;
      double energy = particles[eventno][hitno][5] * 1.6021766E-10;  // energy in Joule
      double theta0 = 2.179E-12/(energy) * z * sqrt(l/X0) * (1 + 0.038*log(l/X0));
      double yplane = (z1*l*theta0)/sqrt(12) +  (z2*l*theta0)/2;
      double thetaplane = z2 * theta0;
      particles[eventno][hitno][j] += yplane + l * particles[eventno][hitno][6+j]; // GLOBAL particles vector
      particles[eventno][hitno][j+6] += thetaplane; // GLOBAL particles vector
    }
  }
}
