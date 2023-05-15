#include "apm.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>

#include "graph.h"
#include "rng.h"

namespace apm {

apmResults_t adaptivePowerMethod(Rng &rng,
                                 matrixSparse pi,
                                 double s,
                                 int t,
                                 vInt k,
                                 double alpha,
                                 bool learningRate,
                                 double r0) {
  apmResults_t results;     // structure for all the results
  const int n = pi.size();  // get the number of nodes

  vInt x(t);            // trajectory
  vDouble cns(t);       // observable estimator
  vDouble cns2(t);      // square of cns
  vDouble kns(t);       // rate function estimator
  vDouble kns2(t);      // square of kns
  vDouble zeta(t);      // eigenvalue
  vDouble zeta2(t);     // square of the eigenvalue
  vDouble psi(t);       // SCGF: log(zeta)
  vDouble psi2(t);      // square of psi
  vDouble psiEst(t);    // estimator of psi
  vDouble psiEst2(t);   // square of psiEst
  vDouble sList(t, s);  // s parameter

  vDouble r(n, r0);        // right eigenvector
  vDouble gamma(n, r0);    // normalization factors
  matrixSparse pi_s = pi;  // driven process transition matrix

  double a = 1;  // learning rate

  int i, j, idx;
  double rTemp;
  vDouble rUpdate(t);  // list of stochastica updates made to r

  if (!learningRate) {
    a = alpha;  // use a constant value
  }

  // initial values
  i = rng.integer(n);
  x.at(0) = i;
  cns.at(0) = k.at(i);
  kns.at(0) = 0;
  zeta.at(0) = 1;

#ifdef LOG_DEBUG
  std::cout << "[Debug]\t[apm]\t"
            << t << " iterations with starting node: "
            << i << '\n';
#endif

  for (int l = 1; l < t; ++l) {
    // choose the index of the neighbors of i
    idx = rng.choice(pi_s.at(i).value);
    // get the node number of the new choice
    j = pi_s.at(i).getIndex(idx);
    x.at(l) = j;

    // update the estimators
    cns.at(l) = cns.at(l - 1) + k.at(j);  // f(i) = k.at(j)
    kns.at(l) = kns.at(l - 1) - std::log(pi.at(i).get(j) / pi_s.at(i).get(j));

    // update the value of r at the last index i
    rUpdate.at(l) = std::exp(s * k.at(i)) * gamma.at(i) / zeta.at(l - 1) - r.at(i);
    rTemp = r.at(i) + a * rUpdate.at(l);
    r.at(i) = rTemp;

    // new eigenvalue estimate
    zeta.at(l) = *std::max_element(r.begin(), r.end());

    // compute the new normalization
    gamma.at(i) = pi.at(i).dot(r);

    // update the i-th row of pi_s
    for (size_t j = 0; j < pi_s.at(i).length(); ++j) {
      pi_s.at(i).value.at(j) = r.at(pi_s.at(i).getIndex(j)) / gamma.at(i) / k.at(i);
    }

    // update the old index
    i = j;
    // update the learning rate
    if (learningRate) {
      a = std::pow(l + 1, -alpha);
    }
  }

  for (size_t i = 0; i < zeta.size(); i++) {
    // compute the squares for the second moment average of multiple runs
    zeta2.at(i) = zeta.at(i) * zeta.at(i);
    psi.at(i) = std::log(zeta.at(i));
    psi2.at(i) = psi.at(i) * psi.at(i);

    // normalize the estimators
    cns.at(i) /= (i + 1);
    if (i > 0) {
      kns.at(i) /= i;
    }

    cns2.at(i) = cns.at(i) * cns.at(i);
    kns2.at(i) = kns.at(i) * kns.at(i);
    // perform the legendre transform to estimate psi
    psiEst.at(i) = cns.at(i) * s - kns.at(i);
    psiEst2.at(i) = psiEst.at(i) * psiEst.at(i);
  }

  // move quantities to return struct
  results.x = std::move(x);
  results.zeta = std::move(zeta);
  results.zeta2 = std::move(zeta2);
  results.psi = std::move(psi);
  results.psi2 = std::move(psi2);
  results.r = std::move(r);
  results.s = std::move(sList);
  results.cns = std::move(cns);
  results.cns2 = std::move(cns2);
  results.kns = std::move(kns);
  results.kns2 = std::move(kns2);
  results.psiEst = std::move(psiEst);
  results.psiEst2 = std::move(psiEst2);
  results.rUpdate = std::move(rUpdate);

  return results;
}

apmResults_t apmTransferLearning(Rng &rng,
                                 matrixSparse pi,
                                 double s,
                                 int epochs,
                                 int t,
                                 vInt k,
                                 double alpha,
                                 bool learningRate,
                                 double r0) {
  apmResults_t results;     // structure for all the results
  const int n = pi.size();  // get the number of nodes

  vInt x(epochs * t);           // trajectory
  vDouble cns(epochs * t);      // observable estimator
  vDouble cns2(epochs * t);     // square of cns
  vDouble kns(epochs * t);      // rate function estimator
  vDouble kns2(epochs * t);     // square of kns
  vDouble zeta(epochs * t);     // eigenvalue
  vDouble zeta2(epochs * t);    // square of the eigenvalue
  vDouble psi(epochs * t);      // SCGF: log(zeta)
  vDouble psi2(epochs * t);     // square of psi
  vDouble psiEst(epochs * t);   // estimator of psi
  vDouble psiEst2(epochs * t);  // square of psiEst
  vDouble sList(epochs * t);    // s parameter

  vDouble r(n, r0);        // right eigenvector
  vDouble gamma(n, r0);    // normalization factors
  matrixSparse pi_s = pi;  // driven process transition matrix

  // vDouble sList(epochs);            // list of s values per epoch
  double sStep = s / (epochs - 1);  // change in s between epochs

  double a = 1;  // learning rate

  int i, j, idx, lGlobal, lmin;
  double rTemp;
  vDouble rUpdate(epochs * t); 

  if (!learningRate) {
    a = alpha;  // use a constant value
  }

#ifdef LOG_DEBUG
  std::cout << "[Debug]\t[apm]\t"
            << epochs << " epochs with each "
            << t << " iterations to s = "
            << s << '\n';
#endif

  // initial values
  zeta.at(0) = 1;
  i = rng.integer(n);
  x.at(0) = i;
  sList.at(0) = 0;

  for (int e = 0; e < epochs; e++) {
    s = e * sStep;
    lGlobal = e * t;
#ifdef LOG_DEBUG
    std::cout << "[Debug]\t[apm]\tEpoch " << e
              << ", s = " << s
              << ", starting node " << i
              << '\n';
#endif

    if (e > 0) {
      lmin = 0;
      // cns.at(lGlobal) = 0;
      // kns.at(lGlobal) = 0;
      // psiEst.at(lGlobal) = cns.at(lGlobal) * s - kns.at(lGlobal);
      // keep the last zeta value from the epoch before
      // zeta.at(e * t) = zeta.at(e * t - 1);
    } else {
      lmin = 1;
    }
    cns.at(lGlobal) = k.at(i);
    kns.at(lGlobal) = 0;
    psiEst.at(lGlobal) = cns.at(lGlobal) * s - kns.at(lGlobal);

    // reset the learning rate
    if (!learningRate) {
      a = alpha;  // use a constant value
    } else {
      a = 1;
    }

    for (int l = lmin; l < t; ++l) {
      lGlobal = e * t + l;
      sList.at(lGlobal) = s;
      // choose the index of the neighbors of i
      idx = rng.choice(pi_s.at(i).value);
      // get the node number of the new choice
      j = pi_s.at(i).getIndex(idx);
      x.at(lGlobal) = j;

      if (l > 1) {
        cns.at(lGlobal) = cns.at(lGlobal - 1) + k.at(j);  // f(i) = k.at(j)
        kns.at(lGlobal) = kns.at(lGlobal - 1) - std::log(pi.at(i).get(j) / pi_s.at(i).get(j));
      } else {
        // reset the estimators in each epoch
        cns.at(lGlobal) += k.at(j);
        kns.at(lGlobal) -= std::log(pi.at(i).get(j) / pi_s.at(i).get(j));
      }
      psiEst.at(lGlobal) = cns.at(lGlobal) * s - kns.at(lGlobal);

      // update the value of r at the last index i
      rUpdate.at(lGlobal) = std::exp(s * k.at(i)) * gamma.at(i) / zeta.at(lGlobal - 1) - r.at(i);
      rTemp = r.at(i) + a * rUpdate.at(lGlobal);
      r.at(i) = rTemp;

      // new eigenvalue estimate
      zeta.at(lGlobal) = *std::max_element(r.begin(), r.end());

      // compute the new normalization
      gamma.at(i) = pi.at(i).dot(r);

      // update the i-th row of pi_s
      for (size_t j = 0; j < pi_s.at(i).length(); ++j) {
        pi_s.at(i).value.at(j) = r.at(pi_s.at(i).getIndex(j)) / gamma.at(i) / k.at(i);
      }

      // update the old index
      i = j;
      // update the learning rate
      if (learningRate) {
        a = std::pow(l + 1, -alpha);
      }
    }  // for during one epoch

    for (int l = 0; l < t; l++) {
      // normalize the estimators
      cns.at(e * t + l) /= (l + 1);
      kns.at(e * t + l) /= l;
      psiEst.at(e * t + l) /= (l + 1);
    }

  }  // for over epochs

  for (size_t i = 0; i < zeta.size(); i++) {
    // compute the squares for the second moment average of multiple runs
    zeta2.at(i) = zeta.at(i) * zeta.at(i);
    psi.at(i) = std::log(zeta.at(i));
    psi2.at(i) = psi.at(i) * psi.at(i);
    cns2.at(i) = cns.at(i) * cns.at(i);
    kns2.at(i) = kns.at(i) * kns.at(i);
    psiEst2.at(i) = psiEst.at(i) * psiEst.at(i);
  }

  results.x = std::move(x);
  results.zeta = std::move(zeta);
  results.zeta2 = std::move(zeta2);
  results.psi = std::move(psi);
  results.psi2 = std::move(psi2);
  results.r = std::move(r);
  results.s = std::move(sList);
  results.cns = std::move(cns);
  results.cns2 = std::move(cns2);
  results.kns = std::move(kns);
  results.kns2 = std::move(kns2);
  results.psiEst = std::move(psiEst);
  results.psiEst2 = std::move(psiEst2);
  results.rUpdate = std::move(rUpdate);

  return results;
}

apmEstimators_t apmRateFunction(Rng &rng,
                                matrixSparse pi,
                                double s,
                                int epochs,
                                int t,
                                vInt k,
                                double alpha,
                                bool learningRate,
                                double r0) {
  apmEstimators_t results;  // structure for all the results
  const int n = pi.size();  // get the number of nodes

  vDouble zeta(epochs * t);  // eigenvalue

  // points per epoch
  vDouble sList(epochs);
  vDouble psi(epochs);
  vDouble psi2(epochs);
  vDouble cns(epochs);
  vDouble cns2(epochs);
  vDouble kns(epochs);
  vDouble kns2(epochs);
  vDouble psiEst(epochs);
  vDouble psiEst2(epochs);

  vDouble r(n, r0);        // right eigenvector
  vDouble gamma(n, r0);    // normalization factors
  matrixSparse pi_s = pi;  // driven process transition matrix

  double a = 1;  // learning rate

  double sStep = s / (epochs - 1);  // change in s between epochs
  int i, j, idx, lmin, lGlobal;
  double rTemp;

  if (!learningRate) {
    a = alpha;  // use a constant value
  }

  // initial values
  i = rng.integer(n);
  cns.at(0) += k.at(i);
  zeta.at(0) = 1;

#ifdef LOG_DEBUG
  std::cout << "[Debug]\t[apm]\trate function: "
            << epochs << " epochs with each "
            << t << " iterations to s = "
            << s << " with starting node: "
            << i << '\n';
#endif

  for (int e = 0; e < epochs; e++) {
    s = e * sStep;

    if (e > 0) {
      lmin = 0;
      cns.at(e) = 0;
      kns.at(e) = 0;
    } else {
      lmin = 1;
      cns.at(0) = k.at(i);
      kns.at(0) = 0;
    }

    // reset the learning rate
    if (!learningRate) {
      a = alpha;  // use a constant value
    } else {
      a = 1;
    }

    for (int l = lmin; l < t; ++l) {
      lGlobal = e * t + l;
      // choose the index of the neighbors of i
      idx = rng.choice(pi_s.at(i).value);
      // get the node number of the new choice
      j = pi_s.at(i).getIndex(idx);

      cns.at(e) += k.at(j);  // f(i) = k.at(j)
      kns.at(e) -= std::log(pi.at(i).get(j) / pi_s.at(i).get(j));

      // update the value of r at the last index i
      rTemp = (1 - a) * r.at(i) + a * std::exp(s * k.at(i)) * gamma.at(i) / zeta.at(lGlobal - 1);
      r.at(i) = rTemp;

      // new eigenvalue estimate
      zeta.at(lGlobal) = *std::max_element(r.begin(), r.end());

      // compute the new normalization
      gamma.at(i) = pi.at(i).dot(r);

      // update the i-th row of pi_s
      for (size_t j = 0; j < pi_s.at(i).length(); ++j) {
        pi_s.at(i).value.at(j) = r.at(pi_s.at(i).getIndex(j)) / gamma.at(i) / k.at(i);
      }

      // update the old index
      i = j;
      // update the learning rate
      if (learningRate) {
        a = std::pow(l + 1, -alpha);
      }
    }  // for during one epoch

    psi.at(e) = std::log(zeta.at((e + 1) * t - 1));
    sList.at(e) = s;
  }  // for over epochs

  for (size_t e = 0; e < cns.size(); e++) {
    // normalize the estimators
    cns.at(e) /= t;
    kns.at(e) /= (t - 1);
    psiEst.at(e) = cns.at(e) * sList.at(e) - kns.at(e);

    // compute the second moments
    psi2.at(e) = psi.at(e) * psi.at(e);
    cns2.at(e) = cns.at(e) * cns.at(e);
    kns2.at(e) = kns.at(e) * kns.at(e);
    psiEst2.at(e) = psiEst.at(e) * psiEst.at(e);
  }

  results.s = std::move(sList);
  results.psi = std::move(psi);
  results.psi2 = std::move(psi2);
  results.cns = std::move(cns);
  results.cns2 = std::move(cns2);
  results.kns = std::move(kns);
  results.kns2 = std::move(kns2);
  results.psiEst = std::move(psiEst);
  results.psiEst2 = std::move(psiEst2);

  return results;
}

void apmResultsAdd(apmResults_t &left, apmResults_t &right) {
  size_t lenLeft = left.zeta.size();
  size_t lenRight = right.zeta.size();

  if (lenLeft == 0) {
    // copy if left side is empty
    left.zeta = right.zeta;
    left.zeta2 = right.zeta2;
    left.psi = right.psi;
    left.psi2 = right.psi2;
    left.r = right.r;
    left.s = right.s;
    left.cns = right.cns;
    left.cns2 = right.cns2;
    left.kns = right.kns;
    left.kns2 = right.kns2;
    left.psiEst = right.psiEst;
    left.psiEst2 = right.psiEst2;
    left.rUpdate = right.rUpdate;
  } else if (lenLeft == lenRight) {
    // add right to left
    for (size_t i = 0; i < lenLeft; i++) {
      left.s.at(i) = right.s.at(i);
      left.cns.at(i) += right.cns.at(i);
      left.cns2.at(i) += right.cns2.at(i);
      left.kns.at(i) += right.kns.at(i);
      left.kns2.at(i) += right.kns2.at(i);
      left.zeta.at(i) += right.zeta.at(i);
      left.zeta2.at(i) += right.zeta2.at(i);
      left.psi.at(i) += right.psi.at(i);
      left.psi2.at(i) += right.psi2.at(i);
      left.psiEst.at(i) += right.psiEst.at(i);
      left.psiEst2.at(i) += right.psiEst2.at(i);
      left.rUpdate.at(i) += right.rUpdate.at(i);
    }
    for (size_t i = 0; i < left.r.size(); i++) {
      left.r.at(i) += right.r.at(i);
    }

  } else {
    std::cout << "[Error]\t[apm]\tCannot add results ad the sizes disagree. Doing nothing!\n";
  }
}

void apmResultsDivide(apmResults_t &res, int repeats) {
  // skip if there is only 1 repeat
  if (repeats > 1) {
    for (size_t i = 0; i < res.zeta.size(); i++) {
      res.s.at(i) /= repeats;
      res.cns.at(i) /= repeats;
      res.cns2.at(i) /= repeats;
      res.kns.at(i) /= repeats;
      res.kns2.at(i) /= repeats;
      res.zeta.at(i) /= repeats;
      res.zeta2.at(i) /= repeats;
      res.psi.at(i) /= repeats;
      res.psi2.at(i) /= repeats;
      res.psiEst.at(i) /= repeats;
      res.psiEst2.at(i) /= repeats;
      res.rUpdate.at(i) /= repeats;
    }
    for (size_t i = 0; i < res.r.size(); i++) {
      res.r.at(i) /= repeats;
    }
  }
}

void apmRateFunctionAdd(apmEstimators_t &left, apmEstimators_t &right) {
  size_t lenLeft = left.psi.size();
  size_t lenRight = right.psi.size();

  if (lenLeft == 0) {
    // copy if left side is empty
    left.s = right.s;
    left.psi = right.psi;
    left.psi2 = right.psi2;
    left.cns = right.cns;
    left.cns2 = right.cns2;
    left.kns = right.kns;
    left.kns2 = right.kns2;
    left.psiEst = right.psiEst;
    left.psiEst2 = right.psiEst2;
  } else if (lenLeft == lenRight) {
    // add right to left
    // left.s = right.s;
    for (size_t i = 0; i < lenLeft; i++) {
      left.psi.at(i) += right.psi.at(i);
      left.psi2.at(i) += right.psi2.at(i);
      left.cns.at(i) += right.cns.at(i);
      left.cns2.at(i) += right.cns2.at(i);
      left.kns.at(i) += right.kns.at(i);
      left.kns2.at(i) += right.kns2.at(i);
      left.psiEst.at(i) += right.psiEst.at(i);
      left.psiEst2.at(i) += right.psiEst2.at(i);
    }

  } else {
    std::cout << "[Error]\t[apm]\tCannot add results ad the sizes disagree. Doing nothing!\n";
  }
}

void apmRateFunctionDivide(apmEstimators_t &res, int repeats) {
  // skip if there is only 1 repeat
  if (repeats > 1) {
    for (size_t i = 0; i < res.psi.size(); i++) {
      res.psi.at(i) /= repeats;
      res.psi2.at(i) /= repeats;
      res.cns.at(i) /= repeats;
      res.cns2.at(i) /= repeats;
      res.kns.at(i) /= repeats;
      res.kns2.at(i) /= repeats;
      res.psiEst.at(i) /= repeats;
      res.psiEst2.at(i) /= repeats;
    }
  }
}

apmResults_t apmLearnS(Rng &rng,
                       matrixSparse pi,
                       double c,
                       int epochs,
                       int t,
                       vInt k,
                       double alpha,
                       double beta,
                       double b0,
                       bool learningRate,
                       double r0) {
  apmResults_t results;     // structure for all the results
  const int n = pi.size();  // get the number of nodes

  vInt x(epochs * t);           // trajectory
  vDouble cns(epochs * t);      // observable estimator
  vDouble cns2(epochs * t);     // square of cns
  vDouble kns(epochs * t);      // rate function estimator
  vDouble kns2(epochs * t);     // square of kns
  vDouble zeta(epochs * t);     // eigenvalue
  vDouble zeta2(epochs * t);    // square of the eigenvalue
  vDouble psi(epochs * t);      // SCGF: log(zeta)
  vDouble psi2(epochs * t);     // square of psi
  vDouble psiEst(epochs * t);   // estimator of psi
  vDouble psiEst2(epochs * t);  // square of psiEst
  vDouble sList(epochs * t);    // s parameter

  vDouble r(n, r0);        // right eigenvector
  vDouble gamma(n, r0);    // normalization factors
  matrixSparse pi_s = pi;  // driven process transition matrix

  double s = 0;
  double cCurrent = 0;

  double a = 1;  // learning rate
  double b = 1;

  int i, j, idx, lGlobal, lmin;
  double rTemp;

  if (!learningRate) {
    a = alpha;  // use a constant value
  }

#ifdef LOG_DEBUG
  std::cout << "[Debug]\t[apm]\t"
            << epochs << " epochs with each "
            << t << " iterations to c = "
            << c << '\n';
#endif

  // initial values
  zeta.at(0) = 1;
  i = rng.integer(n);
  x.at(0) = i;
  sList.at(0) = 0;

  for (int e = 0; e < epochs; e++) {
    lGlobal = e * t;
    b = b0 * std::pow(e + 1, -beta);
    // b = beta / (e + 1);
#ifdef LOG_DEBUG
    std::cout << "[Debug]\t[apm]\tEpoch " << e
              << ", s = " << s
              << ", starting node " << i
              << ", c = " << cCurrent
              << ", b = " << b
              << '\n';
#endif

    if (e > 0) {
      lmin = 0;
    } else {
      lmin = 1;
    }

    cns.at(lGlobal) = k.at(i);
    kns.at(lGlobal) = 0;
    psiEst.at(lGlobal) = cns.at(lGlobal) * s - kns.at(lGlobal);

    // reset the learning rate
    if (!learningRate) {
      a = alpha;  // use a constant value
    } else {
      a = 1;
    }

    for (int l = lmin; l < t; ++l) {
      lGlobal = e * t + l;
      sList.at(lGlobal) = s;
      // choose the index of the neighbors of i
      idx = rng.choice(pi_s.at(i).value);
      // get the node number of the new choice
      j = pi_s.at(i).getIndex(idx);
      x.at(lGlobal) = j;

      if (l > 1) {
        cns.at(lGlobal) = cns.at(lGlobal - 1) + k.at(j);  // f(i) = k.at(j)
        kns.at(lGlobal) = kns.at(lGlobal - 1) - std::log(pi.at(i).get(j) / pi_s.at(i).get(j));
      } else {
        // reset the estimators in each epoch
        cns.at(lGlobal) += k.at(j);
        kns.at(lGlobal) -= std::log(pi.at(i).get(j) / pi_s.at(i).get(j));
      }
      psiEst.at(lGlobal) = cns.at(lGlobal) * s - kns.at(lGlobal);

      // update the value of r at the last index i
      rTemp = (1 - a) * r.at(i) + a * std::exp(s * k.at(i)) * gamma.at(i) / zeta.at(lGlobal - 1);
      r.at(i) = rTemp;

      // new eigenvalue estimate
      zeta.at(lGlobal) = *std::max_element(r.begin(), r.end());

      // compute the new normalization
      gamma.at(i) = pi.at(i).dot(r);

      // update the i-th row of pi_s
      for (size_t j = 0; j < pi_s.at(i).length(); ++j) {
        pi_s.at(i).value.at(j) = r.at(pi_s.at(i).getIndex(j)) / gamma.at(i) / k.at(i);
      }

      // update the old index
      i = j;
      // update the learning rate
      if (learningRate) {
        a = std::pow(l + 1, -alpha);
      }
    }  // for during one epoch

    for (int l = 0; l < t; l++) {
      // normalize the estimators
      cns.at(e * t + l) /= (l + 1);
      kns.at(e * t + l) /= l;
      psiEst.at(e * t + l) /= (l + 1);
    }

    cCurrent = cns.at(lGlobal);
    s -= b * (cCurrent - c);

  }  // for over epochs

  for (size_t i = 0; i < zeta.size(); i++) {
    // compute the squares for the second moment average of multiple runs
    zeta2.at(i) = zeta.at(i) * zeta.at(i);
    psi.at(i) = std::log(zeta.at(i));
    psi2.at(i) = psi.at(i) * psi.at(i);
    cns2.at(i) = cns.at(i) * cns.at(i);
    kns2.at(i) = kns.at(i) * kns.at(i);
    psiEst2.at(i) = psiEst.at(i) * psiEst.at(i);
  }

  results.x = std::move(x);
  results.zeta = std::move(zeta);
  results.zeta2 = std::move(zeta2);
  results.psi = std::move(psi);
  results.psi2 = std::move(psi2);
  results.r = std::move(r);
  results.s = std::move(sList);
  results.cns = std::move(cns);
  results.cns2 = std::move(cns2);
  results.kns = std::move(kns);
  results.kns2 = std::move(kns2);
  results.psiEst = std::move(psiEst);
  results.psiEst2 = std::move(psiEst2);

  return results;
}

}  // namespace apm