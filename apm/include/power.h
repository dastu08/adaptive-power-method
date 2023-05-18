#ifndef _POWER_H_
#define _POWER_H_

#include "general.h"

namespace apm {

// time series of the power method values
typedef struct
{
  vDouble zeta;
  vDouble psi;
  vDouble r;
  double s;
} powerResults_t;

/* Power method to find the largest eigenvalue and vector

**Parameter**
  - pi: transition matrix of the unbiased random walk
  - s: s parameter of the tilted matrix
  - k: degree sequence of the graph
  - alpha: learning rate parameter
  - learningrate: flag to switch from a sequence a = l^(-alpha) if true to a
    constant value a = alpha if false
  - r0: initial condition for the eigenvector components

**Return**
  Struct containing the power method iteration values.

**Description**
  Initialize the right eigenvector with the initial condition.
  Iterate the matrix-vector products. Choose the vector normalization to the
  maximums norm to 1. Estimate the eigenvalue by the largest component of the
  vector after the matrix-vector product (before normalization). Finish after
  t iterations.
*/
powerResults_t powerMethod(matrixSparse pi,
                           double s,
                           int t,
                           vInt k,
                           double alpha,
                           bool learningRate,
                           double r0 = 1.);

}  // namespace apm

#endif  // _POWER_H_