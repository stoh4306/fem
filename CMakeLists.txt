cmake_minimum_required(VERSION 3.7)
project (fem LANGUAGES Fortran)
add_executable(fem fem.f90 fem-matrix.f90)
