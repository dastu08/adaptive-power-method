#ifndef _EIGEN_SOLVER_H_
#define _EIGEN_SOLVER_H_

#include "general.h"

namespace apm {

// dgeev_ is a symbol in the LAPACK library files
extern "C" {
// https://netlib.org/lapack/explore-html/d9/d8e/group__double_g_eeigen_ga66e19253344358f5dee1e60502b9e96f.html
extern int dgeev_(char *JOBVL,
                  char *JOBVR,
                  int *N,
                  double *A,
                  int *LDA,
                  double *WR,
                  double *WI,
                  double *VL,
                  int *LDVL,
                  double *VR,
                  int *LDVR,
                  double *WORK,
                  int *LWORK,
                  int *INFO);
}

// Fortran 2D array indexing (column major)
inline int fcm(int i, int j, int n) {
  return n * j + i;
}

// print the n-by-n matrix a to stdout
void printMatrix(double *a, int n);
void printMatrix(matrixDense a);

// solve for the eigenvalues of a and
// return the dominant eigenvalue la and vector r
int dominatEigVec(double *a, int n, double *la, double *r, double *l = nullptr);

// compute the dominant eigenvalue la of a
int dominatEigVal(double *a, int n, double *la);

}  // namespace apm

#endif  // _EIGEN_SOLVER_H_
