#include "rng.h"

#include <chrono>
#include <iostream>

namespace apm {
// class Rng

Rng::Rng(bool logging) {
  // initialize the RNG
  auto time = std::chrono::system_clock::now();
  unsigned int seed = time.time_since_epoch().count();
  this->seed = seed;
  if (logging) {
    std::cout << "[Info]\t[rng]\tseed: " << seed << '\n';
  }
  srand48(seed);
}

double Rng::real() {
  // int modulus = 1000000;
  // int r = rand();
  // r %= modulus;
  // double d = static_cast<double>(r) / modulus;
  // return d;
  return drand48();
}

int Rng::integer(int n) {
  return lrand48() % n;
}

int Rng::choice(vDouble &probabilities) {
  const double tol = 1e-6;
  double x = Rng::real();
  size_t n = probabilities.size();
  vDouble pSum(n);
  bool found = false;
  int index = 0;

  // create the cum sum
  pSum.at(0) = probabilities.at(0);
  for (size_t i = 1; i < n; ++i) {
    pSum.at(i) = pSum.at(i - 1) + probabilities.at(i);
  }
  // check the sum of probabilities
  if (std::abs(1 - pSum.at(n - 1)) > tol) {
    std::cout << "[Error]\t[rng]\tProbabilities do not sum to one but: "
              << pSum.at(n - 1)
              << ". Returning index -1.\n";
    return -1;
  }

  // see where x is located.
  while ((found == false) && (index < static_cast<int>(n))) {
    // continue if we need to go larger
    if (x > pSum.at(index)) {
      ++index;
    } else {
      found = true;
    }
  }

  return index;
}

}  // namespace apm
