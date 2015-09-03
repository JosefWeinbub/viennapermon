
// System includes
#include <iostream>
#include <mpi.h>

// ViennaGrid includes
#include "viennagridpp/core.hpp"
#include "viennagridpp/io/netgen_reader.hpp"

// libMesh includes
#include "libmesh/serial_mesh.h"
#include "libmesh/cell_tet4.h"
#include "libmesh/gmsh_io.h"

int main(int argc, char* argv[])
{
  //  ----
  // Setup test ViennaGrid mesh
  //  ----
  typedef viennagrid::mesh                                          MeshType;
  typedef viennagrid::result_of::vertex_range<MeshType>::type       VertexRange;
  typedef viennagrid::result_of::iterator<VertexRange>::type        VertexIterator;
  typedef viennagrid::result_of::cell_range<MeshType>::type       CellRange;
  typedef viennagrid::result_of::iterator<CellRange>::type          CellIterator;
  typedef viennagrid::result_of::element<MeshType>::type            CellType;
  typedef viennagrid::result_of::vertex_range<CellType>::type       VertexOnCellRange;
  typedef viennagrid::result_of::iterator<VertexOnCellRange>::type  VertexOnCellIterator;

  MeshType   vgrid_mesh;
  viennagrid::io::netgen_reader reader;
  reader(vgrid_mesh, "../../tests/data/cube48.mesh");

  //  ----
  // Transfer to libMesh mesh data structure
  //  ----
  libMesh::LibMeshInit init (argc, argv);
  libMesh::SerialMesh libmesh(init.comm());
  libmesh.set_mesh_dimension(3);

  viennagrid::mesh_point_accessor mesh_point_accessor(vgrid_mesh);

  unsigned int i = 0;
  std::map<unsigned int, unsigned int> nodetrans;

  VertexRange vgrid_vertices(vgrid_mesh);
  for (VertexIterator vit = vgrid_vertices.begin();  //STL-like iteration
                      vit != vgrid_vertices.end();
                    ++vit)
  {
    libmesh.add_point(libMesh::Point(mesh_point_accessor.get( *vit )[0], mesh_point_accessor.get( *vit )[1], mesh_point_accessor.get( *vit )[2]), i);
    nodetrans[(*vit).id()] = i++;
  }

  CellRange cells_mesh(vgrid_mesh);
  for (CellIterator cit = cells_mesh.begin();
                    cit != cells_mesh.end();
                  ++cit)
  {
    libMesh::Elem* elem = new libMesh::Tet4();

    VertexOnCellRange vertex_on_cells(*cit);
    int j = 0;
    for(VertexOnCellIterator vocit  = vertex_on_cells.begin();
                             vocit != vertex_on_cells.end();
                           ++vocit)
    { 
      elem->set_node(j) = libmesh.node_ptr(libMesh::cast_int<libMesh::dof_id_type>( nodetrans[ (*vocit).id() ] ));

      j++;
    }
    elem->set_id((*cit).id());
    libmesh.add_elem(elem);
 }

  libMesh::GmshIO io(libmesh);
  io.write("libmesh_output.msh");



}
