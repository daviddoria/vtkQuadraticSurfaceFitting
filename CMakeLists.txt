cmake_minimum_required(VERSION 2.6)

PROJECT(vtkQuadricSurfaceFittingExample)

FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

include_directories(${CMAKE_SOURCE_DIR}/../vtkBestFitPlane)

ADD_EXECUTABLE(vtkQuadricSurfaceFittingExample vtkQuadricSurfaceFittingExample.cxx
vtkQuadricSurfaceFitting.cxx
${CMAKE_SOURCE_DIR}/../vtkBestFitPlane/vtkBestFitPlane.cxx
${CMAKE_SOURCE_DIR}/../vtkBestFitPlane/vtkCenterOfMass.cxx)
TARGET_LINK_LIBRARIES(vtkQuadricSurfaceFittingExample vtkHybrid )
