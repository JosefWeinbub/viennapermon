#ifndef VIENNAGRID_2_LIBMESH_HPP
#define VIENNAGRID_2_LIBMESH_HPP

// libMesh includes
#include "libmesh/serial_mesh.h"
#include "libmesh/cell_tet4.h"
#include "libmesh/gmsh_io.h"

void viennagrid2libmesh(viennagrid::mesh const& vgrid_mesh, libMesh::SerialMesh& libmesh)
{
  typedef viennagrid::mesh                                          MeshType;
  typedef viennagrid::result_of::vertex_range<MeshType>::type       VertexRange;
  typedef viennagrid::result_of::iterator<VertexRange>::type        VertexIterator;
  typedef viennagrid::result_of::cell_range<MeshType>::type         CellRange;
  typedef viennagrid::result_of::iterator<CellRange>::type          CellIterator;
  typedef viennagrid::result_of::element<MeshType>::type            CellType;
  typedef viennagrid::result_of::vertex_range<CellType>::type       VertexOnCellRange;
  typedef viennagrid::result_of::iterator<VertexOnCellRange>::type  VertexOnCellIterator;

  libmesh.clear();
  libmesh.set_mesh_dimension(3);

  // transfer geometry
  //
  unsigned int i = 0;
  viennagrid::mesh_point_accessor mesh_point_accessor(vgrid_mesh);
  std::map<viennagrid::element_id, unsigned int> nodetrans;

  VertexRange vgrid_vertices(vgrid_mesh);
  for (VertexIterator vit = vgrid_vertices.begin();  //STL-like iteration
                      vit != vgrid_vertices.end();
                    ++vit)
  {
    libmesh.add_point(libMesh::Point(mesh_point_accessor.get( *vit )[0], mesh_point_accessor.get( *vit )[1], mesh_point_accessor.get( *vit )[2]), i);
    nodetrans[(*vit).id()] = i++;
  }

  // transfer topology
  //
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
    elem->set_id((*cit).id().index());
    libmesh.add_elem(elem);
 }
}



#endif
