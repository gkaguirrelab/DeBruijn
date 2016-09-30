cmake_minimum_required (VERSION 2.6)

project (DeBruijn CXX)

find_package (PkgConfig REQUIRED)
pkg_check_modules (FFTW REQUIRED fftw3)

add_executable (debruijn debruijn.cpp)
include_directories (${FFTW_INCLUDE_DIRS})
target_link_libraries (debruijn ${FFTW_LIBRARIES})

install (TARGETS debruijn RUNTIME DESTINATION bin)
