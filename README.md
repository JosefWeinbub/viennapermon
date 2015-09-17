# viennapermon

This is the development repository for the joint effort of interfacing the Vienna* libraries with Permon* libraries.
The goal is to build a large-scale distributed FETI application by making use of synergy effects provided by
different free open source scientific software infrastructures.

Dependencies:
- CMake
- Metis
- MPI
- LibMesh
- ViennaUtils: https://github.com/viennautils/viennautils-dev  branch: master
- ViennaGrid:  https://github.com/viennagrid/viennagrid-dev    branch: FlorianRudolf/feature-viennagrid-3
- ViennaMesh:  https://github.com/viennamesh/viennamesh-dev    branch: FlorianRudolf/dev-viennamesh-2


Environment variables:
export VIENNAUTILSPATH=/path/to/your/viennautils-dev
export VIENNAGRIDPATH=/path/to/your/viennagrid-dev
export VIENNAMESHPATH=/path/to/your/viennamesh-dev


Configure:
cmake -DLIBMESH_DIR=/path/to/libmesh/instal/directory/ ..





libMesh Build Instructions: (tested with 0.9.4)
  ./configure --prefix=/install/path/ --disable-mpi --disable-fortran --disable-tetgen


