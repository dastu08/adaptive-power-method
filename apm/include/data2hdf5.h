#ifndef _DATA2HDF5_H_
#define _DATA2HDF5_H_

#include <string>

#include "apm.h"
#include "general.h"
#include "hdf5.h"
#include "power.h"

namespace apm {

/* Write a vector of double into a HDF5 file.

**Parameters**
  - vec: vector to save
  - name: name of the vector in the file, should start with a '/'
  - file_id: id of a writeable HDF5 file

**Return**
  HDF5 status type.

**Description**
  Create a one dimensional data set and write the vector as the data set.
*/
herr_t writeH5VectorDouble(vDouble &vec,
                           std::string name,
                           hid_t &file_id);

/* Write a vector of int into a HDF5 file.

**Parameters**
  - vec: vector to save
  - name: name of the vector in the file, should start with a '/'
  - file_id: id of a writeable HDF5 file

**Return**
  HDF5 status type.

**Description**
  Create a one dimensional data set and write the vector as the data set.
*/
herr_t writeH5VectorInt(vInt &vec,
                        std::string name,
                        hid_t &file_id);

/* Write an edge list into a HDF5 file.

**Parameters**
  - edgeList: edge list to save
  - name: prefix of the edge list names in the file, should start with a '/'
  - file_id: id of a writeable HDF5 file

**Return**
  HDF5 status type.

**Description**
  Create two one dimensional data sets.  Write the edge list as `<name>_left`
  and `<name>_right` into the file with each index corresponding to one edge.
*/
herr_t writeH5EdgeList(vvInt &edgeList,
                       std::string name,
                       hid_t &file_id);

/* Write a matreix of double into a HDF5 file.

**Parameters**
  - mat: matrix to save
  - name: name of the matrix in the file, should start with a '/'
  - file_id: id of a writeable HDF5 file

**Return**
  HDF5 status type.

**Description**
  Create a two dimensional data set and write the matrix as the data set.
*/
herr_t writeH5MatrixDouble(vvDouble &mat,
                           std::string name,
                           hid_t &file_id);

/* Write a matreix of integers into a HDF5 file.

**Parameters**
  - mat: matrix to save
  - name: name of the matrix in the file, should start with a '/'
  - file_id: id of a writeable HDF5 file

**Return**
  HDF5 status type.

**Description**
  Create a two dimensional data set and write the matrix as the data set.
*/
herr_t writeH5MatrixInt(vvInt &mat,
                        std::string name,
                        hid_t &file_id);

/* Save an APM trajectory into a HDF5 file.

**Parameters**
  - fileName: name of the file to save, can contain directories
  - results: struct of the APM results that are saved
  - k: degree sequence of the graph
  - edgeList: edge list of the graph

**Return**
  HDF5 status type.

**Description**
  Create the file specified by `fileName`. Save the time series of `s`, `x`,
  `zeta`, `psi`, `cns`, `kns`, `psiEst`. Save the final right eigenvector `r`.
  Save the graph information `k` and `edgeList`. Close the file.
*/
herr_t saveTrajectorySingle(std::string fileName,
                            apmResults_t &results,
                            vInt &k,
                            vvInt &edgeList);

/* Save an APM averaged trajectory into a HDF5 file.

**Parameters**
  - fileName: name of the file to save, can contain directories
  - results: struct of the APM results that are saved
  - k: degree sequence of the graph
  - edgeList: edge list of the graph

**Return**
  HDF5 status type.

**Description**
  Like `saveTrajectorySingle` but also saving the 2nd moments
  `psi2`, `zeta2`, `cns2`, `kns2` and `psiEst2`.
*/
herr_t saveTrajectoryAveraged(std::string fileName,
                              apmResults_t &results,
                              vInt &k,
                              vvInt &edgeList);

/* Save an APM two averaged trajectories into a HDF5 file.

**Parameters**
  - fileName: name of the file to save, can contain directories
  - results: struct of the APM results that are saved
  - resultsNoLr: struct of the APM results that are saved, post-fix `nolr`
  - k: degree sequence of the graph
  - edgeList: edge list of the graph

**Return**
  HDF5 status type.

**Description**
  Saving `psi`, `psi2`, `zeta` and `zeta2` both for the regular case and the
  nolr case. Post-fix is `nolr`. Save the graph information `k` and `edgeList`.
*/
herr_t saveTrajectoryAveragedLr(std::string fileName,
                                apmResults_t &results,
                                apmResults_t &resultsNoLr,
                                vInt &k,
                                vvInt &edgeList);

/* Save the data points of the rate funciton into a HDF5 file.

**Parameters**
  - fileName: name of the file to save, can contain directories
  - resultsLeft: struct of the APM rate function results for s < 0
  - resultsRight: struct of the APM rate function results for s > 0
  - k: degree sequence of the graph
  - edgeList: edge list of the graph

**Return**
  HDF5 status type.

**Description**
  Saving the data points of `s` `psi`, `psi2`, `cns`, `cns2`, `kns`, `kns2`,
  `psiEst` and `psiEst2` for the left and right side of the rate function.
  Use post-fixes `right` and `left`. Save the graph information `k` and
  `edgeList`.
*/
herr_t saveRateFunction(std::string fileName,
                        apmEstimators_t &resultsLeft,
                        apmEstimators_t &resultsRight,
                        vInt &k,
                        vvInt &edgeList);

/* Save the iteration results fo the power method into a HDF5 file.

**Parameter**
  - fileName: name of the file to save, can contain directories
  - results: struct of the power method results
  - k: degree sequence of the graph
  - edgeList: edge list of the graph

**Return**
  HDF5 status type.

**Description**
  Save the series of `zeta`, `psi`. Save the final right eigenvector `r`.
  Save the graph information `k` and `edgeList`.
*/
herr_t savePowerMethod(std::string fileName,
                       powerResults_t &results,
                       vInt &k,
                       vvInt &edgeList);

}  // namespace apm

#endif  // _DATA2HDF5_H_