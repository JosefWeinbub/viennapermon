
IF(DOWNLOAD)
  ExternalProject_Add(viennagrid
    PREFIX viennagrid
    GIT_REPOSITORY https://github.com/viennagrid/viennagrid-dev.git
    GIT_TAG next
    BINARY_DIR "${CMAKE_CURRENT_BINARY_DIR}/viennagrid"
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
  )
  ExternalProject_Get_Property(viennagrid SOURCE_DIR)
ELSE(DOWNLOAD)
  ExternalProject_Add(viennagrid
    PREFIX viennagrid
    SOURCE_DIR $ENV{VIENNAGRIDPATH}
    BINARY_DIR "${CMAKE_CURRENT_BINARY_DIR}/viennagrid"
    CMAKE_ARGS -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
      -DBUILD_TESTING=OFF -DBUILD_EXAMPLES=OFF
    INSTALL_COMMAND ""
  )
  ExternalProject_Get_Property(viennagrid SOURCE_DIR)
  ExternalProject_Get_Property(viennagrid BINARY_DIR)
ENDIF(DOWNLOAD)


SET(VIENNAGRID_INCLUDE_DIRS ${SOURCE_DIR}/include)

IF(BUILD_SHARED_LIBS)
  set(LIBSUFFIX ".so")
ELSE(BUILD_SHARED_LIBS)
  set(LIBSUFFIX ".a")
ENDIF(BUILD_SHARED_LIBS)

SET(VIENNAGRID_LIBRARIES ${BINARY_DIR}/src/${CMAKE_FIND_LIBRARY_PREFIXES}viennagrid${LIBSUFFIX})
