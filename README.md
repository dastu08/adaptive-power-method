# Adaptive Power Method (APM)

This project is an implementation of the Adaptive Power Method which I use in my
master thesis [[1]()]. It simulates a random walk on an Erdős-Rény random graph
to compute the scaled culumant generating function (SCGF) and rate function for
the time-additive observable of the mean degree.

**Table of Contents**
- [Adaptive Power Method (APM)](#adaptive-power-method-apm)
    - [Project Structure](#project-structure)
    - [Prerequities](#prerequities)
    - [Compiling](#compiling)
    - [Usage](#usage)
    - [Data](#data)
    - [License](#license)
    - [Software versions](#software-versions)

---

## Project Structure
This project consists of a library `apm` in the `apm/` subfolder an the main
file `main.cpp` for the executable `simapm` defined in `CMakeLists.txt`. 

To just use the program check the [prerequisites](#prerequities) and then go to
the [compiling](#compiling) section. How to use the program is described in
[usage](#usage).

If you want to use the `apm` library in another project just add the following
lines to your `CMakeLists.txt` file before you define an executable with
`add_executable()`,
```cmake
set(APM_DIR "<path-to-apm-sub-folder>")
add_subdirectory("${APM_DIR}/src" "${APM_DIR}/build")
link_libraries(apm)
```

Have a look at `main.cpp` to see how you can use the `apm` function in your
project. 

## Prerequities
The data produced by the APM is saved in the HDF5 data format. To compile the
library `apm` you need the binary and header file of the HDF5 library. The
version used at the time of writing is `1.14.0`.

I recommend to build the HDF5 library with [CMake](www.cmake.org). Download
`CMake-hdf5-1.14.0.tar.gz` from [HDF5
Group](https://www.hdfgroup.org/downloads/hdf5/source-code/) and check the
`sha256` checksum.

In a desired folder run the following to build the HDF5 library.
```bash
tar xvzf CMake-hdf5-1.14.0.tar.gz 
cd CMake-hdf5-1.14.0/
sh build-unix.sh
```
The result is a file `HDF5-1.14.0-Linux.tar.gz`. Unpack this into the `apm`
folder.
```bash
cd apm/
tar xvzf HDF5-1.14.0-Linux.tar.gz
```
Check that the path set in `apm/src/CMakeLists.txt` under `APM_DIR` matches the
path of the HDF5 files.


## Compiling

Make sure you have compiled the HDF5 library as described under
[prerequisites](#prerequities). Also you need the programs `CMake`, `Ninja` and
`g++` (from GCC). The version that I used when starting this project can be
found under [software versions](#software-versions).

When you are ready, open a terminal inside the source directory where the top most
`CMakeLists.txt` file is located. The following commands create a directory for
the build files and build the project in `Release` mode. The binaries are
produced in the build folder.

```bash
mkdir build/
cd build/
cmake -G Ninja -DCMAKE_BUILD_TYPE:STRING=Release ..
cmake --build .
```

## Usage

The different modes of the APM are invoked by using the correct command line
parameters. The first argument is a mode that is a string of the list `single`,
`transfer`, `ratefunc` or `power`. After that follows a list of `key`-`value`
pairs. The order of the key-value pair does not matter. Important is that the
key has to start with a hyphen `-`. The most important key is `-g` which has 2
values `-g 100 3` the first being the number of nodes for the graph and the
second the average node degree. All other pairs have only one value.

Below you find a list of the possible keys and a description of what they do.

| Key         | Description                                                                                                             |
| ----------- | ----------------------------------------------------------------------------------------------------------------------- |
| `-g 100 3`  | Graph size is 100 nodes with mean degree 3                                                                              |
| `-t 10000`  | run 10000 time steps                                                                                                    |
| `-s -0.5`   | value for the control parameter s                                                                                       |
| `-f ./data` | Directory where the output file will be written to. WATCH OUT that the the directory MUST EXIST before you run the APM. |
| `-a 0.1`    | Value of the exponent in the learning rate                                                                              |
| `-r 100`    | Repeat the APM 100 times and average the obtained quantities                                                            |
| `-e 10`     | Run the APM for 10 epochs if the mode used transfer learning (`transfer` or `ratefunc`)                                 |


The following example show the most common use cases and which arguments they
need. It is assumed that you have compiled the executable `simapm` as described
in the section [compiling](#compiling).

> You need to create the data folder `data` before running the commands.
```bash
# run the APM once
./build/simapm single -g 50 3 -s 1 -t 10000 -a 0.1 -f "./data" -p "demo" 
# run the APM 100 and average
./build/simapm single -g 50 3 -s 1 -t 10000 -a 0.1 -f "./data" -p "demo"  -r 100
# run the APM with transfer learning 
./build/simapm transfer -g 50 3 -s 1 -t 500 -a 0.1 -f "./data" -p "demo"  -r 100 -e 10
# run the APM to compute the rate function
./build/simapm ratefunc -g 50 3 -s 1 -t 1000 -a 0.1 -f "./data" -p "demo" -r 100 -e 50
# run the APM and the power method for the same graph
./build/simapm "power,single" -g 50 3 -s 1 -t 10000 -r 100 -a 0.1 -f "./data" -p "demo"
```

## Data
The data produced by `simapm` is saved in HDF5 files ending in `.h5`. You can
inspect these files with the HDF5 tool `h5ls <file>`. A typical output of the
file created in the example above is:

```bash
h5ls ./data/demo_single_050_3_s-100_t10000_r001_a010.h5
cns                      Dataset {10000}
cns2                     Dataset {10000}
edgeList_left            Dataset {142}
edgeList_right           Dataset {142}
k                        Dataset {50}
kns                      Dataset {10000}
kns2                     Dataset {10000}
psi                      Dataset {10000}
psi2                     Dataset {10000}
psiEst                   Dataset {10000}
psiEst2                  Dataset {10000}
r                        Dataset {50}
rUpdate                  Dataset {10000}
s                        Dataset {10000}
zeta                     Dataset {10000}
zeta2                    Dataset {10000}
```
All entries in the file are arrays whose length is given in the brackets `{}`.
The ones with a length of `10000` are computed along the random which we
specified with `-t 10000`. Those include `zeta`, `psi`, `cns`, `kns`, `psiEst`,
`s` and `rUpdate`. The keys ending with a `2` are just the second moments of the
quantities when you use for example `-r 10` where the number of repeats for the
averaging is larger than 1.

| Quantity                         | Description                                         |
| -------------------------------- | --------------------------------------------------- |
| `s`                              | control parameter of the APM                        |
| `zeta`                           | dominant eigenvalue                                 |
| `psi`                            | SCGF (log of `zeta`)                                |
| `cns`                            | time-additive obserable of the mean degree          |
| `kns`                            | time-additive obserable of the rate function        |
| `psiEst`                         | SCGF computed point wise from `cns * s - kns`       |
| `rUpdate`                        | time series of the changes of the components of `r` |
| `k`                              | node degrees                                        |
| `r`                              | final right eigenvector                             |
| `edgeList_left`/`edgeList_right` | pairs describing the edges of the graph             |

> You can easily load the data into python with the `h5py` package. 

## License
This projected is licensed under the GPL-3.0 license found in the file
[LICENSE](./LICENSE) by me David Stuhrmann.


## Software versions
| Software | Version |
| -------- | ------- |
| CMake    | 3.22.1  |
| GCC      | 11.3.0  |
| HDF5     | 1.14.0  |
| Ninja    | 1.10.1  |