cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
include ($ENV{VIENNAUTILSPATH}/cmake/vienna_build_common.cmake NO_POLICY_SCOPE)

vienna_build_setup_project(ViennaPermon viennautils viennagrid viennamesh)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_BINARY_DIR}/install/lib/cmake" "${PROJECT_SOURCE_DIR}/cmake")

# CMake 2.8.2 is required as a minimum due to the
# ExternalProject_Add mechanism supporting GIT repositories
cmake_minimum_required(VERSION 2.8.2)

OPTION(BUILD_TESTS   "Build Tests"  ON)
OPTION(BUILD_SHARED_LIBS  "Build all libraries as shared (dynamic) libraries" ON )

#list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

#set default build type
IF (NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE "Release")
ENDIF()

#set compiler flags
SET(CMAKE_CXX_FLAGS "-Wall -Wno-long-long -std=c++11")
SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g ")
SET(CMAKE_CXX_FLAGS_RELEASE "-O3 ")
SET(CMAKE_CXX_FLAGS_DEBUG  "-O0 -g")


FIND_PACKAGE(MPI REQUIRED)
IF(MPI_FOUND)
  SET(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_CXX_LINK_FLAGS})
  INCLUDE_DIRECTORIES(${MPI_CXX_INCLUDE_PATH})
  SET(LIBRARIES ${LIBRARIES} ${MPI_CXX_LIBRARIES})
ENDIF(MPI_FOUND)


FIND_PACKAGE(LIBMESH REQUIRED)
INCLUDE_DIRECTORIES(${LIBMESH_INCLUDE_DIR})
SET(LIBRARIES ${LIBRARIES} ${LIBMESH_LIBRARIES})


vienna_build_add_subprojects()

INCLUDE_DIRECTORIES(".")
INCLUDE_DIRECTORIES("include/")

# ViennaGrid
SET(LIBRARIES ${LIBRARIES} viennagridpp)

# ViennaMesh
SET(LIBRARIES ${LIBRARIES} viennameshpp)

IF(BUILD_TESTS)
  enable_testing()
  ADD_SUBDIRECTORY(tests)
ENDIF(BUILD_TESTS)

#build the application
ADD_EXECUTABLE(viennapermon_sim     src/viennapermon_sim.cpp)
TARGET_LINK_LIBRARIES(viennapermon_sim ${LIBRARIES})

#Compatibility with Qt-Creator
file( GLOB_RECURSE QtCreatorCompatibility_SRC
    viennamini/*.hpp viennamini/*.h)
add_custom_target( QtCreatorCompatibility SOURCES ${QtCreatorCompatibility_SRC} )
