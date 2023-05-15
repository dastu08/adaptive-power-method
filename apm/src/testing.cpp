#include "testing.h"

#include <iostream>
#include <vector>

#include "apm.h"
#include "data2hdf5.h"
#include "general.h"
#include "graph.h"
#include "rng.h"

namespace apm {

void testVectorSparse() {
  std::cout << "*** Starting testVectorSparse\n";
  std::vector<int> idx{1, 5, 8, 9, 23};
  std::vector<double> val{1, 6, 4.5, 23, -3.4};
  VectorSparse row(idx, val);
  row.print();
  std::cout << "*** Finished testVectorSparse\n";
}

void testRng() {
  std::cout << "*** Starting testRng\n";
  Rng rng;
  double x = 0;
  int t = 10000;
  int larger = 0;
  int smaller = 0;
  for (int k = 0; k < 20; k++) {
    x = 0;
    for (int i = 0; i < t; i++) {
      x += rng.real();
    }
    x /= t;
    std::cout << x << '\n';
    if (x > 0.5) {
      larger++;
    } else {
      smaller++;
    }
  }
  std::cout << "smaller: " << smaller
            << ", larger: " << larger << '\n';
  std::cout << "*** Finished testRng\n";
}

void testGraph() {
  Rng rng;
  std::cout << "*** Starting testGraph\n";
  Graph graph(10, 2, rng);
  graph.generate()
      .findGiantComponent()
      .print()
      .printGiantComponentSet();

  matrixSparse adj = graph.adjacencyMatrix();
  std::cout << "adjacency matrix:\n";
  for (auto row : adj) {
    row.print();
  }
  std::cout << std::endl;

  graph.computeDegreeSequence()
      .printDegreeSeqence();

  matrixSparse pi = graph.transitionMatrix();
  std::cout << "transition matrix:\n";
  for (auto row : pi) {
    row.print();
  }
  std::cout << "*** Finished testGrap\n";
}

void testRngChoice() {
  std::cout << "*** Starting testRngChoice\n";
  Rng rng;
  vDouble p{0.05, 0.45, 0.2, 0.3};
  vDouble hist(p.size());
  int num;
  size_t t = 1000000;

  for (size_t i = 0; i < t; i++) {
    num = rng.choice(p);
    if (num < 0) {
      break;
    }

    // add the number to the histogram
    int j = 0;
    while (j != num) {
      j++;
    }
    hist.at(j) += 1;
  }

  for (size_t j = 0; j < hist.size(); ++j) {
    hist.at(j) /= t;
    std::cout << j << ": " << hist.at(j) << '\n';
  }

  std::cout << "*** Finished testRngChoice\n";
}

void testRngInteger() {
  std::cout << "*** Starting testRngInteger\n";
  Rng rng;
  int n = 10;
  // vDouble p{0.05, 0.45, 0.2, 0.3};
  vDouble hist(n);
  int num;
  size_t t = 1000000;

  for (size_t i = 0; i < t; i++) {
    num = rng.integer(n);
    if (num < 0) {
      break;
    }

    // add the number to the histogram
    int j = 0;
    while (j != num) {
      j++;
    }
    hist.at(j) += 1;
  }

  for (size_t j = 0; j < hist.size(); ++j) {
    hist.at(j) /= t;
    std::cout << j << ": " << hist.at(j) << '\n';
  }
  std::cout << "*** Finished testRngInteger\n";
}

void testApmSingle() {
  std::cout << "*** Starting testApmSingle\n";
  size_t t = 10000;
  double alpha = 0.1;
  double s = 1.;
  apmResults_t results;
  Rng rng;

  Graph graph(100, 3, rng);
  graph.generateGiantComponent();
  // .print()
  // .printGiantComponentSet()
  // .printDegreeSeqence();

  vvInt edgeList = graph.getEdgeList();
  matrixSparse pi = graph.transitionMatrix();
  vInt k = graph.getDegreeSequence();

  results = adaptivePowerMethod(rng, pi, s, t, k, alpha);
  std::cout << "Last eigenvalue estimate: " << results.zeta.at(t - 1)
            // << " (log: " << )"
            << ", Cns = " << results.cns.at(t - 1)
            << ", Kns = " << results.kns.at(t - 1)
            << ", s = " << results.s.at(t - 1)
            << ", Cns * s - Kns = " << results.cns.at(t - 1) * results.s.at(t - 1) - results.kns.at(t - 1)
            << std::endl;

  saveTrajectorySingle("data/test_apm_single.h5",
                       results, k, edgeList);

  std::cout << "*** Finished testApmSingle\n";
}

void testApmAveraged() {
  std::cout << "*** Starting testApmAveraged\n";
  size_t t = 10000;
  double alpha = 0.1;
  double s = 1.;
  double repeats = 10;
  apmResults_t results, resultsAveraged;
  Rng rng;

  Graph graph(100, 3, rng);
  graph.generateGiantComponent();

  vvInt edgeList = graph.getEdgeList();
  matrixSparse pi = graph.transitionMatrix();
  vInt k = graph.getDegreeSequence();

  for (int i = 0; i < repeats; i++) {
    std::cout << "Repetition " << i << '\n';
    results = adaptivePowerMethod(rng, pi, s, t, k, alpha);
    apmResultsAdd(resultsAveraged, results);
  }
  apmResultsDivide(resultsAveraged, repeats);

  saveTrajectoryAveraged("data/test_apm_averaged.h5",
                         resultsAveraged,
                         k,
                         edgeList);

  std::cout << "*** Finished testApmAveraged\n";
}

void testApmRateFunction() {
  std::cout << "*** Starting testApmRateFunction\n";
  size_t t = 10000;
  int epochs = 200;
  double alpha = 0.1;
  apmEstimators_t resultsLeft, resultsRight;
  Rng rng;

  Graph graph(100, 3, rng);
  graph.generateGiantComponent();

  vvInt edgeList = graph.getEdgeList();
  matrixSparse pi = graph.transitionMatrix();
  vInt k = graph.getDegreeSequence();

  resultsRight = apm::apmRateFunction(rng, pi, 1., epochs, t, k, alpha);
  resultsLeft = apm::apmRateFunction(rng, pi, -1., epochs, t, k, alpha);

  saveRateFunction("data/test_apm_ratefunc.h5",
                   resultsRight,
                   resultsLeft,
                   k,
                   edgeList);

  std::cout << "*** Finished testApmRateFunction\n";
}

}  // namespace apm
