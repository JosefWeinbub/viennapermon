
// System includes
#include <iostream>

// ViennaGrid includes
#include "viennagridpp/core.hpp"
#include "viennagridpp/io/netgen_reader.hpp"

// Local includes
#include "viennagrid2libmesh.hpp"

int main(int argc, char* argv[])
{
  //  ----
  // Setup ViennaGrid mesh; read a netgen mesh 
  //  ----
  viennagrid::mesh   vgrid_mesh;
  viennagrid::io::netgen_reader reader;
  reader(vgrid_mesh, "../../tests/data/cube48.mesh");

  //  ----
  // Instantiate libMesh data structure and transfer ViennGrid mesh to it
  //  ----
  libMesh::SerialMesh libmesh;
  viennagrid2libmesh(vgrid_mesh, libmesh);

  //  ----
  // Write the libMesh data structure to a gmsh mesh output for verification
  //  ----
  libMesh::GmshIO io(libmesh);
  io.write("libmesh_output.msh");
}
