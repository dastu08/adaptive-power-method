set(HDF5_DIR "${PROJECT_SOURCE_DIR}/HDF5-1.14.0-Linux/HDF_Group/HDF5/1.14.0")
set(APM_DIR "${PROJECT_SOURCE_DIR}/apm")

add_library(apm STATIC
  apm.cpp 
  data2hdf5.cpp 
  general.cpp 
  graph.cpp 
  rng.cpp 
  testing.cpp 
  power.cpp 
  parser.cpp
  eigenSolver.cpp)

# local includes
target_include_directories(apm PUBLIC "${APM_DIR}/include")

# hdf5 library
target_include_directories(apm PUBLIC "${HDF5_DIR}/include")
target_link_directories(apm PUBLIC "${HDF5_DIR}/lib")
target_link_libraries(apm PUBLIC hdf5)
