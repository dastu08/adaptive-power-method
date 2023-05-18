cmake_minimum_required(VERSION 3.0.0)
project(adaptive-power-method VERSION 1.0.1)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -Wall")

# apm libary
set(APM_DIR "./apm")
add_subdirectory("${APM_DIR}/src" "${APM_DIR}/build")
link_libraries(apm)

# exectuables
add_executable(simapm main.cpp)
