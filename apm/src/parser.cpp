#include "parser.h"

#include <iomanip>
#include <iostream>
#include <sstream>

namespace apm {

bool parseArgs(apmParam_t &params,
               vInt &times,
               vDouble &sVals,
               vDouble &alphaVals,
               std::string &learningrate,
               int argc,
               char **argv) {
  bool invalid_args = false;  // flag if the arguments are not valid
  int counter = 2;            // counter of the arguments
  std::string buff;

  // basic length check before going into the parsing
  if ((argc % 2 == 1)) {
    // parse the arguments from argv
    while (counter < argc) {
      buff = argv[counter];

      if (buff == "-g") {
        params.n = std::stoi(argv[++counter]);
        params.kmean = std::stod(argv[++counter]);
      } else if (buff == "-s") {
        sVals.push_back(std::stod(argv[++counter]));
      } else if (buff == "-t") {
        times.push_back(std::stoi(argv[++counter]));
      } else if (buff == "-f") {
        params.data_folder = std::string(argv[++counter]);
      } else if (buff == "-r") {
        params.repeats = std::stoi(argv[++counter]);
      } else if (buff == "-e") {
        params.epochs = std::stoi(argv[++counter]);
      } else if (buff == "-p") {
        params.prefix = std::string(argv[++counter]);
      } else if (buff == "-a") {
        alphaVals.push_back(std::stod(argv[++counter]));
        // } else if (buff == "-b") {
        //   params.beta = std::stod(argv[++counter]);
        // } else if (buff == "-c") {
        //   params.c = std::stod(argv[++counter]);
        //   sVals.push_back(params.c);
      } else if (buff == "-lr") {
        std::string boolString = argv[++counter];
        if (boolString == "false") {
          params.lr = false;
          learningrate = "false";
        } else if (boolString == "both") {
          params.lr = true;
          learningrate = "both";
        } else {
          params.lr = true;
          learningrate = "true";
        }
      } else if (buff == "-ic") {
        std::string boolString = argv[++counter];
        if (boolString == "true") {
          params.initialCondition = true;
          std::cout << "[Info]\t[main]\tUsing norm1 initial condition.\n";
        }
      }
      ++counter;
    }  // end while parsing commands from argv

    // check if the parameter values are okay

    // check graph size
    if (params.n < 1) {
      invalid_args = true;
      std::cout << "[Error]\t[main]\tn=" << params.n << " < 1.\n";
    }

    // check mean degree
    if (params.kmean < 0) {
      invalid_args = true;
      std::cout << "[Error]\t[main]\tkmean=" << params.kmean << " < 0\n";
    }

    // check the s values
    if (sVals.size() == 0) {
      invalid_args = true;
      std::cout << "[Error]\t[main]\tno value for s specified.\n";
    }

    // check number of iterations
    if (times.size() > 0) {
      for (int t : times) {
        if (t < 1) {
          invalid_args = true;
          std::cout << "[Error]\t[main]\tt=" << t << " < 1\n";
        }
      }
    } else {
      invalid_args = true;
      std::cout << "[Error]\t[main]\tno value of t specified.\n";
    }

    // check number of repeats
    if (params.repeats < 1) {
      invalid_args = true;
      std::cout << "[Error]\t[main]\trepeats=" << params.repeats << " < 1\n";
    }

    // check number of epochs
    if (params.epochs < 1) {
      invalid_args = true;
      std::cout << "[Error]\t[main]\tepochs=" << params.epochs << " < 1\n";
    }

    // check data folder
    if (params.data_folder == "") {
      invalid_args = true;
      std::cout << "[Error]\t[main]\tno data folder was specified.\n";
    }

    // check learning rate
    if (alphaVals.size() > 0) {
      for (double alpha : alphaVals) {
        if (alpha < 0) {
          invalid_args = true;
          std::cout << "[Error]\t[main]alpha=" << alpha << " < 0\n";
        }
      }
    } else {
      invalid_args = true;
      std::cout << "[Error]\t[main]\talpha not specified.\n";
    }

    // if (params.beta <= 0) {
    //   invalid_args = true;
    //   std::cout << "[Error]\t[main]\tbeta value is less than 0.\n";
    // }

  } else {
    invalid_args = true;
    std::cout << "[Error]\t[main]\tWrong number of command line arguments. Doing nothing!\n";
  }  // end if checking that the argument length is even

  return invalid_args;
}

void constructOutName(apmParam_t &params) {
  std::ostringstream oss;
  oss << params.data_folder << "/";
  if (params.prefix != "") {
    oss << params.prefix << "_";
  }
  oss << params.arg_mode
      << "_" << std::setw(3) << std::setfill('0') << params.n
      << "_" << params.kmean
      << "_s" << std::setw(3) << std::setfill('0') << std::internal
      << static_cast<int>(params.s * 100)
      << "_t" << std::setw(4) << std::setfill('0') << params.t;
  if (params.arg_mode != "power") {
    oss << "_r" << std::setw(3) << std::setfill('0') << params.repeats;
  }
  oss << "_a" << std::setw(3) << std::setfill('0')
      << static_cast<int>(params.alpha * 100);
  oss << ".h5";
  params.file_path = oss.str();
}

}  // namespace apm