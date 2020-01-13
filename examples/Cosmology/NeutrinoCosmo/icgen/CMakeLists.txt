cmake_minimum_required (VERSION 2.6)
project (icgen)
add_library(create_grf create_grf.cpp)
add_library(read_transfer read_transfer.cpp)
add_library(sampler sampler.cpp)
add_library(infini_lpt infini_lpt.cpp)
add_executable(icgen icgen.cpp)

# find_package (Eigen3 3.3 REQUIRED NO_MODULE)

set(CMAKE_LEGACY_CYGWIN_WIN32, 0)
set(CMAKE_CXX_FLAGS "-I/usr/local/include -L/usr/local/lib -I/usr/include/eigen3 -I/usr/include/hdf5/serial -L/usr/lib/x86_64-linux-gnu/hdf5/serial -std=c++11")
target_link_libraries(icgen fftw3 fftw3_omp fftw3f fftw3f_omp m hdf5_cpp hdf5)
target_link_libraries(icgen create_grf read_transfer sampler infini_lpt)
