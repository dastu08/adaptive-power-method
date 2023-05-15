#ifndef _APM_H_
#define _APM_H_

#include <functional>
#include <set>
#include <vector>

#include "general.h"
#include "rng.h"

namespace apm {

// results of the apm run
typedef struct {
  vInt x;
  vDouble zeta;
  vDouble zeta2;
  vDouble psi;
  vDouble psi2;
  vDouble r;
  vDouble cns;
  vDouble cns2;
  vDouble kns;
  vDouble kns2;
  vDouble psiEst;
  vDouble psiEst2;
  vDouble s;
  vDouble rUpdate;
} apmResults_t;

// results of the apm rate function run
typedef struct {
  vDouble s;
  vDouble psi;
  vDouble psi2;
  vDouble cns;
  vDouble cns2;
  vDouble kns;
  vDouble kns2;
  vDouble psiEst;
  vDouble psiEst2;
} apmEstimators_t;

/* Adaptive Power Method (APM)

**Parameter**
  - rng: random number generator object
  - pi: transition matrix of the unbiased random walk
  - s: s parameter of the driven process
  - t: number of time steps to run
  - k: degree sequence of the graph
  - alpha: learning rate parameter
  - learningrate: flag to switch from a sequence a = l^(-alpha) if true to a
    constant value a = alpha if false
  - r0: initial condition for the eigenvector components

**Return**
  Struct containing the time series of the APM quanties.

**Desciption**
  Initialize the right eigenvector with the initial condition `r0`. Choose a
  random starting node. In a loop generate a new node. Update a component of
  the eigenvector, update the eigenvalue, update the transition matrix.
  After t iterations stop.
*/
apmResults_t adaptivePowerMethod(Rng &rng,
                                 matrixSparse pi,
                                 double s,
                                 int t,
                                 vInt k,
                                 double alpha,
                                 bool learningRate = true,
                                 double r0 = 1.);

/* Adaptive Power Method (APM) with transfer learning

**Parameter**
  - rng: random number generator object
  - pi: transition matrix of the unbiased random walk
  - s: final s parameter of the driven process
  - epochs: number of epochs where the current value of s gets changed
  - t: number of time steps per epoch
  - k: degree sequence of the graph
  - alpha: learning rate parameter
  - learningrate: flag to switch from a sequence a = l^(-alpha) if true to a
    constant value a = alpha if false
  - r0: initial condition for the eigenvector components

**Return**
  Struct containing the time series of the APM quanties.

**Desciption**
  Initialize the right eigenvector with the initial condition `r0`. Choose a
  random starting node. In each epoch: in a loop generate a new node.
  Update a component of the eigenvector, update the eigenvalue, update the
  transition matrix. After t iterations change the value of s by s/(epochs-1).
*/
apmResults_t apmTransferLearning(Rng &rng,
                                 matrixSparse pi,
                                 double s,
                                 int epochs,
                                 int t,
                                 vInt k,
                                 double alpha,
                                 bool learningRate = true,
                                 double r0 = 1.);

/* Compute the rate function form the APM with transfer learning

**Parameter**
  - rng: random number generator object
  - pi: transition matrix of the unbiased random walk
  - s: final s parameter of the driven process
  - epochs: number of epochs where the current value of s gets changed
  - t: number of time steps per epoch
  - k: degree sequence of the graph
  - alpha: learning rate parameter
  - learningrate: flag to switch from a sequence a = l^(-alpha) if true to a
    constant value a = alpha if false
  - r0: initial condition for the eigenvector components

**Return**
  Struct containing the values of the APM quanties at the end of each epoch.

**Desciption**
  Run the APM with transfer learning, see `apmTransferLearning()` but save only
  the values at the end of each epoch.
*/
apmEstimators_t apmRateFunction(Rng &rng,
                                matrixSparse pi,
                                double s,
                                int epochs,
                                int t,
                                vInt k,
                                double alpha,
                                bool learningRate = true,
                                double r0 = 1.);

/* Learn the s parameter though the APM with transfer learning

**Parameter**
  - rng: random number generator object
  - pi: transition matrix of the unbiased random walk
  - c: wanted value of the observable of the driven process
  - epochs: number of epochs where the value of s gets changed
  - t: number of time steps per epoch
  - k: degree sequence of the graph
  - alpha: learning rate parameter for the APM steps
  - beta: learning rate parameter for the learning s in between epochs
  - b0: scaling factor in front of the learning rate with beta
  - learningrate: flag to switch from a sequence a = l^(-alpha) if true to a
    constant value a = alpha if false
  - r0: initial condition for the eigenvector components

**Return**
  Struct containing the time series of the APM quanties.

**Description**
  Initialize the right eigenvector with the initial condition `r0`. Choose a
  random starting node. In each epoch: in a loop generate a new node.
  Update a component of the eigenvector, update the eigenvalue, update the
  transition matrix. After t iterations change the value of s by
  `- b * (cCurrent - c)` where `b = b0 * (e+1)^(-beta)`.
*/
apmResults_t apmLearnS(Rng &rng,
                       matrixSparse pi,
                       double c,
                       int epochs,
                       int t,
                       vInt k,
                       double alpha,
                       double beta,
                       double b0,
                       bool learningRate = true,
                       double r0 = 1.);

// Add point wise the values in the structs.
void apmResultsAdd(apmResults_t &left, apmResults_t &right);
// Divide point wise the values in res by repeats.
void apmResultsDivide(apmResults_t &res, int repeats);

// Add point wise the values in the structs.
void apmRateFunctionAdd(apmEstimators_t &left, apmEstimators_t &right);
// Divide point wise the value in res by repeats.
void apmRateFunctionDivide(apmEstimators_t &res, int repeats);

}  // namespace apm

#endif  // _APM_H_