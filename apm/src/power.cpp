#include "power.h"

#include <cmath>
// #include <algorithm>

namespace apm {

powerResults_t powerMethod(matrixSparse pi,
                           double s,
                           int t,
                           vInt k,
                           double alpha,
                           bool learningRate,
                           double r0) {
  powerResults_t results;   // structure for all the results
  const int n = pi.size();  // get the number of nodes

  vDouble zeta(t);  // eigenvalue
  vDouble psi(t);   // SCGF: log(zeta)

  vDouble r(n, r0);            // right eigenvector
  vDouble rTemp(n);            // temporary vector for the matrix-vector product
  matrixSparse pi_tilde = pi;  // driven process transition matrix

  double a = 1;  // learning rate

  if (!learningRate) {
    a = alpha;  // use a constant value
  }

  // construct the tilted matrix
  for (int i = 0; i < n; ++i) {
    pi_tilde.at(i).multiply(std::exp(s * k.at(i)));
  }

  // power method iterations
  for (int l = 0; l < t; ++l) {
    // matrix-vector product
    for (int i = 0; i < n; ++i) {
      rTemp.at(i) = pi_tilde.at(i).dot(r);
    }
    // find the largest value an estimate that as the eigen value
    zeta.at(l) = *std::max_element(rTemp.begin(), rTemp.end());
    psi.at(l) = std::log(zeta.at(l));

    // copy the temporary vector and divide by the eigenvalue
    // (maximum norm to 1)
    for (int i = 0; i < n; ++i) {
      r.at(i) = (1 - a) * r.at(i) + a * rTemp.at(i) / zeta.at(l);
    }

    // update the learning rate
    if (learningRate) {
      a = std::pow(l + 1, -alpha);
    }
  }

  // move quantities to return struct
  results.s = s;
  results.zeta = std::move(zeta);
  results.psi = std::move(psi);
  results.r = std::move(r);

  return results;
}

}  // namespace apm
