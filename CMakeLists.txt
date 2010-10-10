cmake_minimum_required(VERSION 2.6)

PROJECT(vtkQuadricSurfaceFittingExample)

FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

ADD_EXECUTABLE(vtkQuadricSurfaceFittingExample vtkQuadricSurfaceFittingExample.cxx
vtkQuadricSurfaceFitting.cxx)
TARGET_LINK_LIBRARIES(vtkQuadricSurfaceFittingExample vtkHybrid )
