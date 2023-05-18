#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "apm.h"
#include "data2hdf5.h"
#include "general.h"
#include "graph.h"
#include "parser.h"
#include "power.h"
#include "rng.h"
#include "testing.h"

using namespace apm;

// ./apm test
// ./apm single -g 100 3 -s 1. -t 10000 -f "data"
// ./apm averaging -g 100 3 -s 1. -t 10000 -r 10 -f "data"
// ./apm ratefunc -g 100 3 -s 1. -t 1000 -e 10 -f "data"
// ./apm learningrate -g 100 3 -s 1. -t 10000 -f "data" -r 10
// ./apm transfer -g 100 3 -s 1. -t 10000 -f "data" -e 4
// ./apm transavg -g 100 3 -s 1. -t 10000 -f "data" -e 4 -r 10
// ./apm power -g 100 3 -s 1. -t 10000 -f "data" -a 0.3
// ./apm learns -g 100 3 -c 3 -t 2000 -f "data" -a 0.1 -e 10
int main(int argc, char **argv) {
  herr_t status = 0;  // hdf5 status

  Rng rng;  // random number generator

  apmMode_t mode = NONE;

  vInt times;         // list of t values
  vDouble sVals;      // list of s values
  vDouble alphaVals;  // list of alpha values

  apmParam_t params;          // simulation parameters
  vInt danglingChains;        // list of dangling chain endpoints
  bool modeLooping = false;   // flag for looping everyting again
  size_t chainThreshold = 1;  // threshold for regenerating the graph

  size_t tEnd = 0;
  std::string learningrate = "";

  vString arg_mode_list;

  // check if we got a valid mode
  if (argc <= 1) {
    std::cout << "[Error]\t[main]\t no arguments where given.\n";
    return 1;
  }

#ifdef LOG_DEBUG
  std::cout << "[Debug]\t[main]\tParsing mode argument\n";
#endif
  arg_mode_list = split(argv[1], ',');

  // parse and check the command line arguments
  if (parseArgs(params, times, sVals, alphaVals, learningrate, argc, argv)) {
    return 1;
  }

  // setting up the graph
  Graph graph(params.n, params.kmean, rng);

  // generate the graph as long as the the number of dangling chains
  // is below the threshold
  do {
    graph.generateGiantComponent();
    // graph.printDegreeSeqence();
    danglingChains = graph.findDanglingChains();
  } while (danglingChains.size() < chainThreshold);

  vvInt edgeList = graph.getEdgeList();
  matrixSparse pi = graph.transitionMatrix();
  vInt k = graph.getDegreeSequence();

  // loop over different modes given
  for (std::string arg_mode : arg_mode_list) {
    if (arg_mode == "test") {
      mode = TEST;
    } else if (arg_mode == "single") {
      mode = SINGLE;
    } else if (arg_mode == "ratefunc") {
      mode = RATE_FUNC;
    } else if (arg_mode == "power") {
      mode = POWER;
    } else if (arg_mode == "transfer") {
      mode = TRANSFER;
    } else {
      mode = NONE;
    }

    if (mode == NONE) {
      std::cout << "[Error]\t[main] Could not determine a mode.\n";
      return 1;
    }
    params.arg_mode = arg_mode;

    // loop again if modeLoop is set
    do {
      // adapt the arg_mode if the learning rate needs to be in it
      if (learningrate == "both") {
        if (params.lr) {
          params.arg_mode = arg_mode + "-" + "lr";
        } else {
          params.arg_mode = arg_mode + "-" + "nolr";
        }
      }

      // loop over the list of parameter values
      for (size_t t : times) {
        for (double s : sVals) {
          for (double alpha : alphaVals) {
            // result variables
            apmResults_t results;
            apmResults_t resultsAveraged;
            powerResults_t resultsPower;
            apmEstimators_t resLeft;
            apmEstimators_t resRight;
            apmEstimators_t resLeftAvg;
            apmEstimators_t resRightAvg;

            // modify initial condition
            if (params.initialCondition) {
              params.r0 = 1. / pi.size();
              std::cout << "[Info]\t[main]\t Initial condition changed to 1/n.\n";
            }

            // save loop variables into struct
            params.t = t;
            params.s = s;
            params.alpha = alpha;
            constructOutName(params);

            // performing the desired action
            switch (mode) {
              case TEST:

                testVectorSparse();
                testRng();
                testGraph();
                testRngChoice();
                testRngInteger();
                testApmSingle();
                testApmAveraged();
                testApmRateFunction();
                status = 0;

                break;

              case SINGLE:
                // Run the APM bare
                std::cout << "[Info]\t[main]\tSINGLE"
                          << ", n=" << params.n
                          << ", kmean=" << params.kmean
                          << ", s=" << params.s
                          << ", t=" << params.t
                          << ", repeats=" << params.repeats
                          << ", alpha=" << params.alpha
                          << ", lr=" << params.lr
                          << ", r0=" << params.r0
                          << ", data_folder=" << params.data_folder
                          << '\n';

                for (int i = 0; i < params.repeats; i++) {
#ifdef LOG_VERBOSE
                  std::cout << "[Verbose]\t[main]\tRepetition " << i << "\t";
#endif
                  results = adaptivePowerMethod(rng,
                                                pi,
                                                params.s,
                                                params.t,
                                                k,
                                                params.alpha,
                                                params.lr,
                                                params.r0);
                  apmResultsAdd(resultsAveraged, results);
                }
                apmResultsDivide(resultsAveraged, params.repeats);

                status = saveTrajectoryAveraged(params.file_path,
                                                resultsAveraged,
                                                k,
                                                edgeList);

                tEnd = params.t - 1;
                std::cout << "[Info]\t[main]\tFinal values"
                          << ", zeta=" << resultsAveraged.zeta.at(tEnd)
                          << ", psi=" << resultsAveraged.psi.at(tEnd)
                          << ", cns=" << resultsAveraged.cns.at(tEnd)
                          << ", kns=" << resultsAveraged.kns.at(tEnd)
                          << ", s=" << resultsAveraged.s.at(tEnd)
                          << ", cns*s-kns=" << resultsAveraged.psiEst.at(tEnd)
                          << std::endl;

                break;

              case TRANSFER:
                // Run the APM with transfer learning
                std::cout << "[Info]\t[main]\tTRANSFER"
                          << ", n=" << params.n
                          << ", kmean=" << params.kmean
                          << ", s=" << params.s
                          << ", t=" << params.t
                          << ", epochs=" << params.epochs
                          << ", repeats=" << params.repeats
                          << ", alpha=" << params.alpha
                          << ", lr=" << params.lr
                          << ", r0=" << params.r0
                          << ", data_folder=" << params.data_folder
                          << '\n';

                for (int i = 0; i < params.repeats; i++) {
#ifdef LOG_VERBOSE
                  std::cout << "[Verbose]\t[main]\tRepetition " << i << "\t";
#endif
                  results = apmTransferLearning(rng,
                                                pi,
                                                params.s,
                                                params.epochs,
                                                params.t,
                                                k,
                                                params.alpha,
                                                params.lr,
                                                params.r0);
                  apmResultsAdd(resultsAveraged, results);
                }
                apmResultsDivide(resultsAveraged, params.repeats);

                status = saveTrajectoryAveraged(params.file_path,
                                                resultsAveraged,
                                                k,
                                                edgeList);

                tEnd = params.epochs * params.t - 1;
                std::cout << "[Info]\t[main]\tFinal values"
                          << ", zeta=" << resultsAveraged.zeta.at(tEnd)
                          << ", psi=" << resultsAveraged.psi.at(tEnd)
                          << ", cns=" << resultsAveraged.cns.at(tEnd)
                          << ", kns=" << resultsAveraged.kns.at(tEnd)
                          << ", s=" << resultsAveraged.s.at(tEnd)
                          << ", cns*s-kns=" << resultsAveraged.psiEst.at(tEnd)
                          << std::endl;
                break;

              case RATE_FUNC:
                // Run the APM with transfer learning and extract the rate function
                std::cout << "[Info]\t[main]\tRATE_FUNC"
                          << ", n=" << params.n
                          << ", kmean=" << params.kmean
                          << ", s=" << params.s
                          << ", t=" << params.t
                          << ", epochs=" << params.epochs
                          << ", repeats=" << params.repeats
                          << ", alpha=" << params.alpha
                          << ", lr=" << params.lr
                          << ", data_folder=" << params.data_folder
                          << '\n';

                // flip the sign of s to make it positive
                if (params.s < 0) {
                  params.s = -params.s;
                }

                for (int i = 0; i < params.repeats; i++) {
#ifdef LOG_VERBOSE
                  std::cout << "[Verbose]\t[main]\tRepetition " << i
                            << std::endl;
#endif
                  resLeft = apmRateFunction(rng,
                                            pi,
                                            -params.s,
                                            params.epochs,
                                            params.t,
                                            k,
                                            params.alpha,
                                            params.lr,
                                            params.r0);
                  resRight = apmRateFunction(rng,
                                             pi,
                                             params.s,
                                             params.epochs,
                                             params.t,
                                             k,
                                             params.alpha,
                                             params.lr,
                                             params.r0);

                  apmRateFunctionAdd(resLeftAvg, resLeft);
                  apmRateFunctionAdd(resRightAvg, resRight);
                }
                apmRateFunctionDivide(resLeftAvg, params.repeats);
                apmRateFunctionDivide(resRightAvg, params.repeats);

                status = saveRateFunction(params.file_path,
                                          resLeftAvg,
                                          resRightAvg,
                                          k,
                                          edgeList);
                break;

              case POWER:
                // Run the Power Method
                std::cout << "[Info]\t[main]\tPOWER"
                          << ", n=" << params.n
                          << ", kmean=" << params.kmean
                          << ", s=" << params.s
                          << ", t=" << params.t
                          << ", alpha=" << params.alpha
                          << ", lr=" << params.lr
                          << ", data_folder=" << params.data_folder
                          << '\n';

                resultsPower = powerMethod(pi,
                                           params.s,
                                           params.t,
                                           k,
                                           params.alpha,
                                           params.lr,
                                           params.r0);

                status = savePowerMethod(params.file_path,
                                         resultsPower,
                                         k,
                                         edgeList);

                tEnd = params.t - 1;
                std::cout << "[Info]\t[main]\tFinal values"
                          << ", zeta=" << resultsPower.zeta.at(tEnd)
                          << ", psi=" << resultsPower.psi.at(tEnd)
                          << ", s=" << resultsPower.s
                          << std::endl;

                break;

              default:
                std::cout << "[Error] Unkown mode. Doing nothing!\n";
                status = 0;
                break;
            }  // end switch mode
          }    // end for over different alpha values
        }      // end for over different s values
      }        // end for over different t values

      // loop again if we want both learning rate options and
      // we are just finished with the first one
      if (learningrate == "both") {
        if (params.lr) {
          params.lr = false;
          modeLooping = true;
        } else {
          modeLooping = false;
        }
      }
    } while (modeLooping);  // end while different modes for learning rate

  }  // end for over arg_mode_list

  return status;
}
