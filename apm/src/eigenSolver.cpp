#include "eigenSolver.h"

#include <cmath>
#include <iostream>

namespace apm {

void printMatrix(double *a, int n) {
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      std::cout << a[fcm(i, j, n)] << ", ";
    }
    std::cout << "\n";
  }
}

void printMatrix(matrixDense a) {
  for (size_t i = 0; i < a.size(); i++) {
    for (size_t j = 0; j < a.at(i).size(); j++) {
      std::cout << a.at(i).at(j) << ", ";
    }
    std::cout << "\n";
  }
}

int dominatEigVec(double *a, int n, double *la, double *r, double *l) {
  double *wr = new double[n];
  double *wi = new double[n];
  double *vl = new double[n * n];
  double *vr = new double[n * n];

  // char Nchar = 'N';
  char Vchar = 'V';
  // int one = 1;
  int lwork = 6 * n;
  double *work = new double[lwork];
  int info;

  int i0 = 0;
  double lMax = 0;
  double lTemp = 0;

  dgeev_(&Vchar,
         &Vchar,
         &n,
         a,
         &n,
         wr,
         wi,
         vl,
         &n,
         vr,
         &n,
         work,
         &lwork,
         &info);

  // std::cout << "dgeev_ info: " << info << std::endl;

  // find the largest eigen value
  for (int j = 0; j < n; j++) {
    lTemp = wr[j] * wr[j] + wi[j] * wi[j];

    if (lTemp > lMax) {
      i0 = j;
      lMax = lTemp;
    }

    // std::cout << wr[j] << " + " << wi[j] << " i, (";
    // for (int i = 0; i < n; i++) {
    //   std::cout << vr[fcm(i, j, n)] << ", ";
    // }
    // std::cout << ")\n";
  }

  *la = wr[i0];
  if (l != nullptr) {
    for (int i = 0; i < n; i++) {
      r[i] = vr[fcm(i, i0, n)];
      l[i] = vl[fcm(i, i0, n)];
    }
  } else {
    for (int i = 0; i < n; i++) {
      r[i] = vr[fcm(i, i0, n)];
    }
  }

  delete[] work;
  delete[] wr;
  delete[] wi;
  delete[] vr;
  delete[] vl;

  return info;
}

int dominatEigVal(double *a, int n, double *la) {
  double *wr = new double[n];
  double *wi = new double[n];
  double *vl = nullptr;
  double *vr = nullptr;

  char Nchar = 'N';
  // char Vchar = 'V';
  int one = 1;
  int lwork = 6 * n;
  double *work = new double[lwork];
  int info;

  int i0 = 0;
  double lMax = 0;
  double lTemp = 0;

  dgeev_(&Nchar,
         &Nchar,
         &n,
         a,
         &n,
         wr,
         wi,
         vl,
         &one,
         vr,
         &n,
         work,
         &lwork,
         &info);

  // std::cout << "dgeev_ info: " << info << std::endl;

  // find the largest eigen value
  for (int j = 0; j < n; j++) {
    // lTemp = wr[j] * wr[j] + wi[j] * wi[j];
    lTemp = wr[j];

    if (lTemp > lMax) {
      i0 = j;
      lMax = lTemp;
    }

    // std::cout << wr[j] << " + " << wi[j] << " i, (";
    // for (int i = 0; i < n; i++) {
    //   std::cout << vr[fcm(i, j, n)] << ", ";
    // }
    // std::cout << ")\n";
  }

  if (wr[i0] <= 0) {
    for (int i = 0; i < n; i++) {
      std::cout << i << ": " << wr[i] << ", " << wi[0] << "\n";
    }
  }

  *la = wr[i0];

  delete[] work;
  delete[] wr;
  delete[] wi;

  return info;
}

}  // namespace apm