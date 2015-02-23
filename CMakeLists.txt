PROJECT(ViennaPermon C CXX)

include(ExternalProject)
include(FindPackageHandleStandardArgs)

# CMake 2.8.2 is required as a minimum due to the
# ExternalProject_Add mechanism supporting GIT repositories
cmake_minimum_required(VERSION 2.8.2)

OPTION(BUILD_TESTS   "Build Tests"  ON)

#list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

#set default build type
IF (NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE "Release")
ENDIF()

#set compiler flags
SET(CMAKE_CXX_FLAGS "-Wall -pedantic")
SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g ")
SET(CMAKE_CXX_FLAGS_RELEASE "-O3 ")
SET(CMAKE_CXX_FLAGS_DEBUG  "-O0 -g")


INCLUDE_DIRECTORIES(".")
INCLUDE_DIRECTORIES("include/")

INCLUDE(cmake/GetViennaGrid.cmake)
MESSAGE(STATUS "ViennaGrid Include Path: " ${VIENNAGRID_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${VIENNAGRID_INCLUDE_DIRS})

IF(BUILD_TESTS)
  #INCLUDE(CTest)
  enable_testing()
  ADD_SUBDIRECTORY(tests)
ENDIF(BUILD_TESTS)

##Compatibility with Qt-Creator
file( GLOB_RECURSE QtCreatorCompatibility_SRC
    viennamini/*.hpp viennamini/*.h)
add_custom_target( QtCreatorCompatibility SOURCES ${QtCreatorCompatibility_SRC} )
