PROJECT(ViennaPermon C CXX)

include(ExternalProject)
include(FindPackageHandleStandardArgs)

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
SET(CMAKE_CXX_FLAGS "-Wall -pedantic -Wno-long-long")
SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g ")
SET(CMAKE_CXX_FLAGS_RELEASE "-O3 ")
SET(CMAKE_CXX_FLAGS_DEBUG  "-O0 -g")


FIND_PACKAGE(MPI REQUIRED)
IF(MPI_FOUND)
  SET(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_CXX_LINK_FLAGS})
  INCLUDE_DIRECTORIES(${MPI_CXX_INCLUDE_PATH})
  SET(LIBRARIES ${LIBRARIES} ${MPI_CXX_LIBRARIES})
ENDIF(MPI_FOUND)


INCLUDE_DIRECTORIES(".")
INCLUDE_DIRECTORIES("include/")

INCLUDE(cmake/GetViennaGrid.cmake)
MESSAGE(STATUS "ViennaGrid Include Path: " ${VIENNAGRID_INCLUDE_DIRS})
MESSAGE(STATUS "ViennaGrid Libraries:    " ${VIENNAGRID_LIBRARIES})
INCLUDE_DIRECTORIES(${VIENNAGRID_INCLUDE_DIRS})
SET(LIBRARIES ${LIBRARIES} ${VIENNAGRID_LIBRARIES})


INCLUDE(cmake/GetViennaMesh.cmake)
MESSAGE(STATUS "ViennaMesh Include Path: " ${VIENNAMESH_INCLUDE_DIRS})
MESSAGE(STATUS "ViennaMesh Libraries:    " ${VIENNAMESH_LIBRARIES})

set(VIENNAMESH_PLUGIN_DIRECTORIES "${CMAKE_CURRENT_BINARY_DIR}/viennamesh/plugins")
MESSAGE(STATUS "ViennaGrid Plugins Directory:    ${VIENNAMESH_PLUGIN_DIRECTORIES}" )
add_definitions( "-DVIENNAMESH_DEFAULT_PLUGIN_DIRECTORY=\"${VIENNAMESH_PLUGIN_DIRECTORIES}\"" )

INCLUDE_DIRECTORIES(${VIENNAMESH_INCLUDE_DIRS})
SET(LIBRARIES ${LIBRARIES} ${VIENNAMESH_LIBRARIES})





IF(BUILD_TESTS)
  enable_testing()
  ADD_SUBDIRECTORY(tests)
ENDIF(BUILD_TESTS)

#build the application
ADD_EXECUTABLE(viennapermon_sim     src/viennapermon_sim.cpp)
TARGET_LINK_LIBRARIES(viennapermon_sim ${LIBRARIES})
ADD_DEPENDENCIES(viennapermon_sim viennagrid)
ADD_DEPENDENCIES(viennapermon_sim viennamesh)

#Compatibility with Qt-Creator
file( GLOB_RECURSE QtCreatorCompatibility_SRC
    viennamini/*.hpp viennamini/*.h)
add_custom_target( QtCreatorCompatibility SOURCES ${QtCreatorCompatibility_SRC} )
