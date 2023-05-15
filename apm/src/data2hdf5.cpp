#include "data2hdf5.h"

#include <iostream>

namespace apm {

herr_t writeH5VectorDouble(std::vector<double> &vec, std::string name, hid_t &file_id) {
  hsize_t dims[1];
  dims[0] = vec.size();
  hid_t dataspace_id, dataset_id;
  herr_t status;

  dataspace_id = H5Screate_simple(1, dims, NULL);

  dataset_id = H5Dcreate(file_id,
                         name.c_str(),
                         H5T_IEEE_F64BE,
                         dataspace_id,
                         H5P_DEFAULT,
                         H5P_DEFAULT,
                         H5P_DEFAULT);

  status = H5Dwrite(dataset_id,
                    H5T_NATIVE_DOUBLE,
                    H5S_ALL,
                    H5S_ALL,
                    H5P_DEFAULT,
                    vec.data());

  status = H5Dclose(dataset_id);
  status = H5Sclose(dataspace_id);

#ifdef LOG_DEBUG
  std::cout << "[Debug]\t[h5]\tWrote " << name << " to hdf5 file.\n";
#endif

  return status;
}

herr_t writeH5VectorInt(std::vector<int> &vec, std::string name, hid_t &file_id) {
  hsize_t dims[1];
  dims[0] = vec.size();
  hid_t dataspace_id, dataset_id;
  herr_t status;

  dataspace_id = H5Screate_simple(1, dims, NULL);

  dataset_id = H5Dcreate(file_id,
                         name.c_str(),
                         H5T_STD_I32BE,
                         dataspace_id,
                         H5P_DEFAULT,
                         H5P_DEFAULT,
                         H5P_DEFAULT);

  status = H5Dwrite(dataset_id,
                    H5T_NATIVE_INT,
                    H5S_ALL,
                    H5S_ALL,
                    H5P_DEFAULT,
                    vec.data());

  status = H5Dclose(dataset_id);
  status = H5Sclose(dataspace_id);

#ifdef LOG_DEBUG
  std::cout << "[Debug]\t[h5]\tWrote " << name << " to hdf5 file.\n";
#endif

  return status;
}

herr_t writeH5EdgeList(vvInt &edgeList, std::string name, hid_t &file_id) {
  hsize_t dims[2];
  dims[0] = edgeList.size();
  dims[1] = 1;
  hid_t dataspace_id, datasetLeft_id, datasetRight_id;
  herr_t status;
  std::string nameLeft = name + "_left";
  std::string nameRight = name + "_right";

  vInt left(edgeList.size());
  vInt right(edgeList.size());
  for (size_t i = 0; i < edgeList.size(); i++) {
    left.at(i) = edgeList.at(i).at(0);
    right.at(i) = edgeList.at(i).at(1);
  }

  dataspace_id = H5Screate_simple(1, dims, NULL);

  datasetLeft_id = H5Dcreate(file_id,
                             nameLeft.c_str(),
                             H5T_STD_I32BE,
                             dataspace_id,
                             H5P_DEFAULT,
                             H5P_DEFAULT,
                             H5P_DEFAULT);

  status = H5Dwrite(datasetLeft_id,
                    H5T_NATIVE_INT,
                    H5S_ALL,
                    H5S_ALL,
                    H5P_DEFAULT,
                    left.data());

  status = H5Dclose(datasetLeft_id);

  datasetRight_id = H5Dcreate(file_id,
                              nameRight.c_str(),
                              H5T_STD_I32BE,
                              dataspace_id,
                              H5P_DEFAULT,
                              H5P_DEFAULT,
                              H5P_DEFAULT);

  status = H5Dwrite(datasetRight_id,
                    H5T_NATIVE_INT,
                    H5S_ALL,
                    H5S_ALL,
                    H5P_DEFAULT,
                    right.data());

  status = H5Dclose(datasetRight_id);

  status = H5Sclose(dataspace_id);

#ifdef LOG_DEBUG
  std::cout << "[Debug]\t[h5]\tWrote " << nameLeft
            << " and " << nameRight
            << " to hdf5 file.\n";
#endif

  return status;
}

herr_t writeH5MatrixDouble(vvDouble &mat,
                           std::string name,
                           hid_t &file_id) {
  hsize_t dims[2];
  dims[0] = mat.size();
  dims[1] = mat.at(0).size();
  hid_t dataspace_id, dataset_id, datatype_id;
  herr_t status;
  vDouble buff(dims[0] * dims[1]);

  for (size_t i = 0; i < dims[0]; ++i) {
    for (size_t j = 0; j < dims[1]; ++j) {
      buff.at(i * dims[1] + j) = mat.at(i).at(j);
    }
  }

  datatype_id = H5Tcopy(H5T_IEEE_F64BE);

  dataspace_id = H5Screate_simple(2, dims, NULL);

  dataset_id = H5Dcreate(file_id,
                         name.c_str(),
                         datatype_id,
                         dataspace_id,
                         H5P_DEFAULT,
                         H5P_DEFAULT,
                         H5P_DEFAULT);

  status = H5Dwrite(dataset_id,
                    H5T_NATIVE_DOUBLE,
                    H5S_ALL,
                    H5S_ALL,
                    H5P_DEFAULT,
                    buff.data());

  status = H5Dclose(dataset_id);
  status = H5Sclose(dataspace_id);

#ifdef LOG_DEBUG
  std::cout << "[Debug]\t[h5]\tWrote " << name << " to hdf5 file.\n";
#endif

  return status;
}

herr_t writeH5MatrixInt(vvInt &mat,
                        std::string name,
                        hid_t &file_id) {
  hsize_t dims[2];
  dims[0] = mat.size();
  dims[1] = mat.at(0).size();
  hid_t dataspace_id, dataset_id, datatype_id;
  herr_t status;
  vInt buff(dims[0] * dims[1]);

  for (size_t i = 0; i < dims[0]; ++i) {
    for (size_t j = 0; j < dims[1]; ++j) {
      buff.at(i * dims[1] + j) = mat.at(i).at(j);
    }
  }

  datatype_id = H5Tcopy(H5T_STD_I32BE);

  dataspace_id = H5Screate_simple(2, dims, NULL);

  dataset_id = H5Dcreate(file_id,
                         name.c_str(),
                         datatype_id,
                         dataspace_id,
                         H5P_DEFAULT,
                         H5P_DEFAULT,
                         H5P_DEFAULT);

  status = H5Dwrite(dataset_id,
                    H5T_NATIVE_INT,
                    H5S_ALL,
                    H5S_ALL,
                    H5P_DEFAULT,
                    buff.data());

  status = H5Dclose(dataset_id);
  status = H5Sclose(dataspace_id);

#ifdef LOG_DEBUG
  std::cout << "[Debug]\t[h5]\tWrote " << name << " to hdf5 file.\n";
#endif

  return status;
}

herr_t saveTrajectorySingle(std::string fileName,
                            apmResults_t &results,
                            vInt &k,
                            vvInt &edgeList) {
  hid_t file_id;
  herr_t status = 1;

  file_id = H5Fcreate(fileName.c_str(),
                      H5F_ACC_TRUNC,
                      H5P_DEFAULT,
                      H5P_DEFAULT);

  if (file_id > 0) {
    status = writeH5VectorDouble(results.s, "/s", file_id);
    status = writeH5VectorInt(results.x, "/x", file_id);
    status = writeH5VectorDouble(results.zeta, "/zeta", file_id);
    status = writeH5VectorDouble(results.psi, "/psi", file_id);
    status = writeH5VectorDouble(results.cns, "/cns", file_id);
    status = writeH5VectorDouble(results.kns, "/kns", file_id);
    status = writeH5VectorDouble(results.psiEst, "/psiEst", file_id);
    status = writeH5VectorDouble(results.r, "/r", file_id);
    status = writeH5VectorInt(k, "/k", file_id);
    status = writeH5EdgeList(edgeList, "/edgeList", file_id);

    status = H5Fclose(file_id);

    std::cout << "[Info]\t[h5]\tSaved data into file: " << fileName << std::endl;
  } else {
    std::cout << "[Error]\t[h5]\tcould not open file: " << fileName << std::endl;
  }
  return status;
}

herr_t saveTrajectoryAveraged(std::string fileName,
                              apmResults_t &results,
                              vInt &k,
                              vvInt &edgeList) {
  hid_t file_id;
  herr_t status = 1;

  file_id = H5Fcreate(fileName.c_str(),
                      H5F_ACC_TRUNC,
                      H5P_DEFAULT,
                      H5P_DEFAULT);

  if (file_id > 0) {
    status = writeH5VectorDouble(results.s, "/s", file_id);
    status = writeH5VectorDouble(results.psi, "/psi", file_id);
    status = writeH5VectorDouble(results.psi2, "/psi2", file_id);
    status = writeH5VectorDouble(results.zeta, "/zeta", file_id);
    status = writeH5VectorDouble(results.zeta2, "/zeta2", file_id);
    status = writeH5VectorDouble(results.cns, "/cns", file_id);
    status = writeH5VectorDouble(results.cns2, "/cns2", file_id);
    status = writeH5VectorDouble(results.kns, "/kns", file_id);
    status = writeH5VectorDouble(results.kns2, "/kns2", file_id);
    status = writeH5VectorDouble(results.psiEst, "/psiEst", file_id);
    status = writeH5VectorDouble(results.psiEst2, "/psiEst2", file_id);
    status = writeH5VectorDouble(results.r, "/r", file_id);
    status = writeH5VectorInt(k, "/k", file_id);
    status = writeH5EdgeList(edgeList, "/edgeList", file_id);
    status = writeH5VectorDouble(results.rUpdate, "/rUpdate", file_id);

    status = H5Fclose(file_id);

    std::cout << "[Info]\t[h5]\tSaved data into file: " << fileName << std::endl;
  } else {
    std::cout << "[Error]\t[h5]\tcould not open file: " << fileName << std::endl;
  }
  return status;
}

herr_t saveTrajectoryAveragedLr(std::string fileName,
                                apmResults_t &results,
                                apmResults_t &resultsNoLr,
                                vInt &k,
                                vvInt &edgeList) {
  hid_t file_id;
  herr_t status = 1;

  file_id = H5Fcreate(fileName.c_str(),
                      H5F_ACC_TRUNC,
                      H5P_DEFAULT,
                      H5P_DEFAULT);

  if (file_id > 0) {
    status = writeH5VectorDouble(results.psi, "/psi", file_id);
    status = writeH5VectorDouble(results.psi2, "/psi2", file_id);
    status = writeH5VectorDouble(results.zeta, "/zeta", file_id);
    status = writeH5VectorDouble(results.zeta2, "/zeta2", file_id);

    status = writeH5VectorDouble(resultsNoLr.psi, "/psi_nolr", file_id);
    status = writeH5VectorDouble(resultsNoLr.psi2, "/psi2_nolr", file_id);
    status = writeH5VectorDouble(resultsNoLr.zeta, "/zeta_nolr", file_id);
    status = writeH5VectorDouble(resultsNoLr.zeta2, "/zeta2_nolr", file_id);

    status = writeH5VectorDouble(results.r, "/r", file_id);
    status = writeH5VectorInt(k, "/k", file_id);
    status = writeH5EdgeList(edgeList, "/edgeList", file_id);

    status = H5Fclose(file_id);

    std::cout << "[Info]\t[h5]\tSaved data into file: " << fileName << std::endl;
  } else {
    std::cout << "[Error]\t[h5]\tcould not open file: " << fileName << std::endl;
  }
  return status;
}

herr_t saveRateFunction(std::string fileName,
                        apmEstimators_t &resultsLeft,
                        apmEstimators_t &resultsRight,
                        vInt &k,
                        vvInt &edgeList) {
  hid_t file_id;
  herr_t status = 1;

  file_id = H5Fcreate(fileName.c_str(),
                      H5F_ACC_TRUNC,
                      H5P_DEFAULT,
                      H5P_DEFAULT);

  if (file_id > 0) {
    status = writeH5VectorDouble(resultsRight.s, "/s_right", file_id);
    status = writeH5VectorDouble(resultsRight.psi, "/psi_right", file_id);
    status = writeH5VectorDouble(resultsRight.psi2, "/psi2_right", file_id);
    status = writeH5VectorDouble(resultsRight.cns, "/cns_right", file_id);
    status = writeH5VectorDouble(resultsRight.cns2, "/cns2_right", file_id);
    status = writeH5VectorDouble(resultsRight.kns, "/kns_right", file_id);
    status = writeH5VectorDouble(resultsRight.kns2, "/kns2_right", file_id);
    status = writeH5VectorDouble(resultsRight.psiEst, "/psiEst_right", file_id);
    status = writeH5VectorDouble(resultsRight.psiEst2, "/psiEst2_right", file_id);

    status = writeH5VectorDouble(resultsLeft.s, "/s_left", file_id);
    status = writeH5VectorDouble(resultsLeft.psi, "/psi_left", file_id);
    status = writeH5VectorDouble(resultsLeft.psi2, "/psi2_left", file_id);
    status = writeH5VectorDouble(resultsLeft.cns, "/cns_left", file_id);
    status = writeH5VectorDouble(resultsLeft.cns2, "/cns2_left", file_id);
    status = writeH5VectorDouble(resultsLeft.kns, "/kns_left", file_id);
    status = writeH5VectorDouble(resultsLeft.kns2, "/kns2_left", file_id);
    status = writeH5VectorDouble(resultsLeft.psiEst, "/psiEst_left", file_id);
    status = writeH5VectorDouble(resultsLeft.psiEst2, "/psiEst2_left", file_id);

    status = writeH5VectorInt(k, "/k", file_id);
    status = writeH5EdgeList(edgeList, "/edgeList", file_id);

    status = H5Fclose(file_id);

    std::cout << "[Info]\t[h5]\tSaved data into file: " << fileName << std::endl;
  } else {
    std::cout << "[Error]\t[h5]\tcould not open file: " << fileName << std::endl;
  }

  return status;
}

herr_t savePowerMethod(std::string fileName,
                       powerResults_t &results,
                       vInt &k,
                       vvInt &edgeList) {
  hid_t file_id;
  herr_t status = 1;

  file_id = H5Fcreate(fileName.c_str(),
                      H5F_ACC_TRUNC,
                      H5P_DEFAULT,
                      H5P_DEFAULT);

  if (file_id > 0) {
    status = writeH5VectorDouble(results.zeta, "/zeta", file_id);
    status = writeH5VectorDouble(results.psi, "/psi", file_id);
    status = writeH5VectorDouble(results.r, "/r", file_id);
    status = writeH5VectorInt(k, "/k", file_id);
    status = writeH5EdgeList(edgeList, "/edgeList", file_id);

    status = H5Fclose(file_id);

    std::cout << "[Info]\t[h5]\tSaved data into file: " << fileName << std::endl;
  } else {
    std::cout << "[Error]\t[h5]\tcould not open file: " << fileName << std::endl;
  }
  return status;
}

}  // namespace apm