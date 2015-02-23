

ADD_EXECUTABLE(viennagrid_mpi-test     viennagrid_mpi.cpp)
TARGET_LINK_LIBRARIES(viennagrid_mpi-test ${LIBRARIES})
ADD_DEPENDENCIES(viennagrid_mpi-test viennagrid)
ADD_TEST(viennagrid_mpi ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 viennagrid_mpi-test)
