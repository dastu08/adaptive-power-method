#ifndef _PARSER_H_
#define _PARSER_H_

#include <string>

#include "general.h"

namespace apm {

// Mode switch for performing different simulaitons
typedef enum {
  NONE,
  TEST,
  SINGLE,
  TRANSFER,
  RATE_FUNC,
  POWER
} apmMode_t;

// parameters for a simulation
typedef struct {
  int n = 1;
  int kmean = 0;
  int epochs = 2;
  int repeats = 1;
  size_t t = 0;
  size_t tEnd = 0;
  double s = 0;
  // double c = 0;
  double alpha = 0.1;
  // double beta = 1;
  bool lr = true;
  double r0 = 1.;
  bool initialCondition = false;

  std::string data_folder = "";
  std::string file_path = "";
  std::string arg_mode = "";
  std::string prefix = "";

} apmParam_t;

/* Parse the command line arguments and save them in a struct

**Parameter**
  - params: structs to hold the parameters in
  - times: vector of values of t
  - sVals: vector of values of s
  - alphaVals: vector of values of alpha
  - learningrate: string of the `lr` option
  - argc: argc from main
  - argv: argv from main

**Return**
  Value of `invalid_args`.

**Description**
  Go throught the list of the command line arguments given in `argv`. Start
  from the index 2 that is the first argument after the mode.
*/
bool parseArgs(apmParam_t &params,
               vInt &times,
               vDouble &sVals,
               vDouble &alphaVals,
               std::string &learningrate,
               int argc,
               char **argv);

/* Generate the file name from the parameters

**Parameter**
  - params: struct of the simulation parameters

**Description**
  Build the string of the output file name and save it into the struct.
  Example: "single_050_3_s100_t10000_r010_a010.h5"
*/
void constructOutName(apmParam_t &params);

}  // namespace apm

#endif  // _PARSER_H_