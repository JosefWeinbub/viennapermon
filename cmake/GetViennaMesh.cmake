
IF(DOWNLOAD)
  ExternalProject_Add(viennamesh
    PREFIX viennamesh
    GIT_REPOSITORY https://github.com/viennamesh/viennamesh-dev.git
    GIT_TAG next
    BINARY_DIR "${CMAKE_CURRENT_BINARY_DIR}/viennamesh"
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
  )
  ExternalProject_Get_Property(viennamesh SOURCE_DIR)
ELSE(DOWNLOAD)
  ExternalProject_Add(viennamesh
    PREFIX viennamesh
    SOURCE_DIR $ENV{VIENNAMESHPATH}
    BINARY_DIR "${CMAKE_CURRENT_BINARY_DIR}/viennamesh"
    CMAKE_ARGS -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
      -DBUILD_TESTING=OFF -DBUILD_EXAMPLES=OFF
    INSTALL_COMMAND ""
  )
  ExternalProject_Get_Property(viennamesh SOURCE_DIR)
  ExternalProject_Get_Property(viennamesh BINARY_DIR)
ENDIF(DOWNLOAD)


SET(VIENNAMESH_INCLUDE_DIRS ${SOURCE_DIR}/include ${SOURCE_DIR}/external/pugixml-1.5/src)

IF(BUILD_SHARED_LIBS)
  set(LIBSUFFIX ".so")
ELSE(BUILD_SHARED_LIBS)
  set(LIBSUFFIX ".a")
ENDIF(BUILD_SHARED_LIBS)

SET(VIENNAMESH_LIBRARIES ${BINARY_DIR}/src/${CMAKE_FIND_LIBRARY_PREFIXES}viennamesh${LIBSUFFIX})
