

#include <iostream>
#include <mpi.h>

#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/metis_partitioner.h"


#include "libmesh/exodusII_io.h"
//#include "libmesh/vtk_io.h"

int main(int argc, char* argv[])
{
  libMesh::LibMeshInit init (argc, argv);
  {
  std::cout << "Testing build_cube -> parallel mesh mechanism .." << std::endl;
  libMesh::Mesh libmesh(init.comm());
  int ps = 15;
  int dim = 3;
  libMesh::Real halfwidth = dim > 1 ? 1. : 0.;
  libMesh::Real halfheight = dim > 2 ? 1. : 0.;
  libMesh::MeshTools::Generation::build_cube (libmesh,
                                     ps,
                                     (dim>1) ? ps : 0,
                                     (dim>2) ? ps : 0,
                                     -1., 1.,
                                     -halfwidth, halfwidth,
                                     -halfheight, halfheight,
                                     (dim==1)    ? libMesh::EDGE2 :
                                     ((dim == 2) ? libMesh::QUAD4 : libMesh::HEX8));

  libmesh.print_info();
  }
  {
  std::cout << std::endl;
  std::cout << "Testing mesh reader -> parallel mesh mechanism .." << std::endl;
  libMesh::Mesh mesh(init.comm());
  mesh.read("../data/miscellaneous_ex9.exo");
  mesh.print_info();
  }

  return 0;
}
