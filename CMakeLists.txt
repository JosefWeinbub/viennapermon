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

#FIND_PACKAGE(VTK REQUIRED)
#INCLUDE_DIRECTORIES(${VTK_INCLUDE_DIRS})
#SET(LIBRARIES ${LIBRARIES} ${VTK_LIBRARIES})

#IF(VTK_FOUND)
#  message(STATUS "Found VTK version " ${VTK_MAJOR_VERSION}.${VTK_MINOR_VERSION})
#  message(STATUS "VTK include path: ${VTK_INCLUDE_DIRS}")
##  message(STATUS "VTK library path: ${VTK_LIBRARY_DIRS}")
##  message(STATUS "VTK libraries: ${VTK_LIBRARIES}")
##  message(STATUS "VTK available kits: ${VTK_KITS}")
#ENDIF(VTK_FOUND)

foreach(DEF ${LIBMESH_DEFINITIONS})
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${DEF}")
endforeach(DEF)


# Configure subproject: ViennaMesh
option(ENABLE_PLUGIN_TRIANGLE "Enable netgen plugin" ON)
option(ENABLE_PLUGIN_IO "Enable io plugin" ON)
option(ENABLE_PLUGIN_TETGEN "Enable tetgen plugin" ON)
option(ENABLE_PLUGIN_NETGEN "Enable netgen plugin" ON)
option(ENABLE_PLUGIN_CGAL "Enable cgal plugin" OFF)
option(ENABLE_PLUGIN_POISSON "Enable poisson plugin" OFF)
option(ENABLE_PLUGIN_METIS "Enable metis plugin" ON)
option(ENABLE_PLUGIN_TDR "Enable TDR plugin" OFF)
option(ENABLE_PLUGIN_SYMMETRY "Enable symmetry plugin" OFF)
option(ENABLE_PLUGIN_STATISTICS "Enable symmetry plugin" ON)
option(ENABLE_PLUGIN_MESH_HEALING "Enable mesh healing plugin" ON)
option(BUILD_EXAMPLES "Build examples" OFF)

vienna_build_add_subprojects()

INCLUDE_DIRECTORIES(".")
INCLUDE_DIRECTORIES("include/")

# ViennaGrid
SET(LIBRARIES ${LIBRARIES} viennagridpp)

# ViennaMesh
SET(LIBRARIES ${LIBRARIES} viennameshpp)


IF(BUILD_TESTS)
  enable_testing()
#  ADD_SUBDIRECTORY(tests)
ENDIF(BUILD_TESTS)

#build the application
ADD_EXECUTABLE(viennapermon_sim     src/viennapermon_sim.cpp)
TARGET_LINK_LIBRARIES(viennapermon_sim ${LIBRARIES})

ADD_EXECUTABLE(devel_libmetis_parmesh     src/devel_libmetis_parmesh.cpp)
TARGET_LINK_LIBRARIES(devel_libmetis_parmesh ${LIBRARIES})



#Compatibility with Qt-Creator
file( GLOB_RECURSE QtCreatorCompatibility_SRC
    viennamini/*.hpp viennamini/*.h)
add_custom_target( QtCreatorCompatibility SOURCES ${QtCreatorCompatibility_SRC} )
