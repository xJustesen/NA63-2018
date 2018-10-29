#ifndef SIMPLEX_H
#define SIMPLEX_H

#include <vector>
#include <armadillo>

using namespace std;
using namespace arma;

void simplex_update(vector<vec> simplex, vec fs, vec &centroid, int &ihigh, int &ilow);
double simplex_size(vector<vec> simplex);
void simplex_reduce(vector<vec> &simplex, int ilow);
vec simplex_expand(vec highest, vec centroid);
vec simplex_reflect(vec highest, vec centroid);
vec simplex_contract(vec highest, vec centroid);
void simplex_initiate(vector<vec> simplex, function<double(vec)> F, vec &fs);

template<typename lambda>
void simulator::simplex_initiate(vector<vec> simplex, lambda&& F, vec &fs) {

  for (size_t i = 0; i < simplex.size(); i++) {

    fs(i) = F(simplex[i]);

  }

}

/* Nelder-Mead simplex algorithm. Tested on Himmelblau's function and Rosenbrock function. In principle able to optimize n-dimensional problems. */
template<typename lambda>
vec simulator::simplex_NM(lambda&& F, vector<vec> simplex, double simplex_size_goal) {

  vec fs; fs.resize(simplex.size());
  vec centroid = zeros<vec>(simplex.size() - 1);
  simplex_initiate(simplex, F, fs);
  int ilow, ihigh;

  while (simplex_size(simplex) > simplex_size_goal) {

    simplex_update(simplex, fs, centroid, ihigh, ilow); // update simplex with new fs values, this updates centroid, ihigh, ilow.
    vec r = simplex_reflect(simplex[ihigh], centroid);  // reflection
    double fr = F(r);

    /* try expansion */
    if (fr < fs(ilow)) {

      vec e = simplex_expand(simplex[ihigh], centroid);
      double fe = F(e);

      if (fe < fr ){ // accept expansion

        simplex[ihigh] = e;
        fs(ihigh) = fe;

      }

      else { // reject expansion and accept reflection

        simplex[ihigh] = r;
        fs(ihigh) = fr;

      }

    }

    else {

       /*if reflection is too poor, reflect in other direction */
      if (fr < fs(ihigh)) { // accept reflection

        simplex[ihigh] = r;
        fs(ihigh) = fr;

      }

      else { // reject reflection, try contraction

        vec c = simplex_contract(simplex[ihigh], centroid);
        double fc = F(c);

        if (fc < fs(ihigh)){ // accept contraction

          simplex[ihigh] = c;
          fs(ihigh) = fc;

        }
        else { // reject contraction and reduce

          simplex_reduce(simplex, ilow);
          simplex_initiate(simplex, F, fs);

        }

      }

    }

  } // end while

  return simplex[ilow];

}

#endif
