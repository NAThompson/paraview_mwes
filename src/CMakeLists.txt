##============================================================================
##  Copyright (c) Kitware, Inc.
##  All rights reserved.
##  See LICENSE.txt for details.
##
##  This software is distributed WITHOUT ANY WARRANTY; without even
##  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
##  PURPOSE.  See the above copyright notice for more information.
##============================================================================
cmake_minimum_required(VERSION 3.8...3.15 FATAL_ERROR)
project(ParaviewTutorial CXX)

find_package(VTKm REQUIRED QUIET)

include_directories(/usr/local/include)

add_executable(WaveletTransform WaveletTransform.cxx)
target_link_libraries(WaveletTransform PRIVATE vtkm_filter)
set_property(TARGET WaveletTransform PROPERTY CXX_STANDARD 17)

add_executable(KdV KdV.cxx)
target_link_libraries(KdV PRIVATE vtkm_filter)
set_property(TARGET KdV PROPERTY CXX_STANDARD 17)

add_executable(JacobiTheta JacobiTheta.cxx)
target_link_libraries(JacobiTheta PRIVATE vtkm_filter)
set_property(TARGET JacobiTheta PROPERTY CXX_STANDARD 17)

add_executable(Fredholm Fredholm.cxx)
target_link_libraries(Fredholm PRIVATE vtkm_filter)
set_property(TARGET Fredholm PROPERTY CXX_STANDARD 17)

add_executable(EulerSpiral EulerSpiral.cxx)
target_link_libraries(EulerSpiral PRIVATE vtkm_filter)
set_property(TARGET EulerSpiral PROPERTY CXX_STANDARD 17)

add_executable(GrayScott GrayScott.cxx)
target_link_libraries(GrayScott PRIVATE vtkm_filter)
set_property(TARGET GrayScott PROPERTY CXX_STANDARD 17)

add_executable(ScatteredData ScatteredData.cxx)
target_link_libraries(ScatteredData PRIVATE vtkm_filter)
set_property(TARGET ScatteredData PROPERTY CXX_STANDARD 17)

add_executable(Hydrogen Hydrogen.cxx)
target_link_libraries(Hydrogen PRIVATE vtkm_filter)
set_property(TARGET Hydrogen PROPERTY CXX_STANDARD 17)
