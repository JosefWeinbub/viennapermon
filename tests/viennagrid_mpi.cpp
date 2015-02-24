
// System includes
#include <iostream>
#include <mpi.h>

// ViennaGrid includes
#include "viennagrid/core.hpp"
#include "viennagrid/io/vtk_reader.hpp"
#include "viennagrid/io/vtk_writer.hpp"

// Local includes
#include "viennagrid_mpi.hpp"

int main(int argc, char* argv[])
{
  int mpi_size, mpi_rank;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

  if(mpi_size != 2)
  {
    std::cout << "This application must be executed with 2 MPI processes - aborting!" << std::endl;
    MPI_Finalize();
    return 1;
  }

  typedef viennagrid::mesh_hierarchy_t MeshHierarchyType;
  typedef viennagrid::result_of::mesh<MeshHierarchyType>::type MeshType;

  viennagrid_serialized_mesh_hierarchy serialized_mesh;
  viennagrid_serialized_mesh_hierarchy_make(&serialized_mesh);

  if(mpi_rank == 0)
  {
    viennagrid::io::vtk_reader<MeshType> reader;
    MeshHierarchyType mesh_hierarchy;
    MeshType mesh = mesh_hierarchy.root();
    reader(mesh, "../../tests/data/tets_with_data_main.pvd");
    mesh_hierarchy.serialize(serialized_mesh, true);

    viennagrid::mpi::serializedmesh_operator sermeshop(serialized_mesh);

    MPI_Send(&sermeshop.get_blocklen_size(), 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
    MPI_Send(sermeshop.get_blocklen(), sermeshop.get_blocklen_size(), MPI_INT, 1, 0, MPI_COMM_WORLD);
    MPI_Send(&sermeshop.get_mesh(), 1, sermeshop.get_mpi_type(), 1, 0, MPI_COMM_WORLD);
  }
  else
  {
    int blocklen_size = 0;
    MPI_Recv(&blocklen_size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    int* blocklen = new int[blocklen_size];
    MPI_Recv(blocklen, blocklen_size, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    viennagrid::mpi::serializedmesh_operator sermeshop(serialized_mesh, blocklen, blocklen_size);
//    MPI_Recv(&sermeshop.get_mesh(), 1, sermeshop.get_mpi_type(), 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

//    MeshHierarchyType mesh_hierarchy;
//    mesh_hierarchy.deserialize(sermeshop.get_mesh());
//    MeshType mesh = mesh_hierarchy.root();
//    viennagrid::io::vtk_writer<MeshType> writer;
//    writer(mesh, "deserialized");
  }

  viennagrid_serialized_mesh_hierarchy_delete(serialized_mesh);




//  viennagrid_serialized_mesh_hierarchy serialized_mesh;
//  viennagrid_serialized_mesh_hierarchy_make(&serialized_mesh);

//--> mpi send the serialized_mesh

//backend/api.h:218
//size(points)=vertex_count*geometric_dimension
//size(hole_points_offsets)=hole_point_element_count+1
//size(hole_points)=hole_points_offsets[hole_point_element_count]
//size(cell_element_tags)=cell_count
//size(cell_vertex_offsets)=cell_count+1
//size(cell_vertices)=cell_vertex_offsets[cell_count]
//size(cell_parents)=cell_count
//size(cell_region_offsets)=cell_count+1
//size(cell_regions)=cell_region_offsets[cell_count]
//size(mesh_parents)=mesh_count
//size(mesh_cell_count)=mesh_count
//size(mesh_cells)=mesh_count
//size(mesh_cells[i])=mesh_cell_count[i]
//size(region_ids)=region_count
//size(region_names)=region_count
//size(region_names[i])=null_terminated (aka. strlen(region_names[i])+1)

  MPI_Finalize();
  return 0;
}
