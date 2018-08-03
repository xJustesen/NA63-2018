void simulator::make_intensity_distro(void){
  /* Load data of spatial distribution of energy into matrix, define vector for saving intensities */
  mat intensity_spatial; intensity_spatial.load("sum_initials1mm40GeVelec.txt");
  vector<vector<double>> intensities(energies_interp.size());
  for (size_t i = 0; i < intensities.size(); i++) intensities[i].resize(intensity_spatial.n_cols);
  // vector<double> intensities_i(intensity_spatial.n_rows);
  vec vec_emitted_energies = linspace<vec>(0.0, 40.0, intensity_spatial.n_rows);
  vec vec_energies_interp(energies_interp.size());
  for (size_t k = 0; k < energies_interp.size(); k++){
    vec_energies_interp(k) = energies_interp[k];
  }
  /* Calculate spatial intensity distro. for every energy */
  for (size_t j = 0; j < energies_interp.size(); j++){
    cerr << "Outer loop : Iteration no. " << j+1 << " of " << energies_interp.size() << "\n";
    for (size_t i = 0; i < intensity_spatial.n_cols; i++){
      /* Extract column i from intensity_spatial matrix */
      vec intensities_i = intensity_spatial.col(i);

      vec interp_intensities_i;
      interp1(vec_emitted_energies, intensities_i, vec_energies_interp, interp_intensities_i,"*linear");

      /* Fill "intensities" vector. The j'th entry corresponds to given energy, the i'th entry is the intensity at the angle corresponding to the i'th column in "sum_initials1mm40GeVelec" vector */
      intensities[j][i] = interp_intensities_i(j);
    }
  }
}
