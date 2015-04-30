
// System includes
#include <iostream>
#include <mpi.h>

// ViennaGrid includes
#include "viennagridpp/core.hpp"
#include "viennagridpp/io/vtk_reader.hpp"
#include "viennagridpp/io/vtk_writer.hpp"

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

  if(mpi_rank == 0)
  {
    // Test 1
    // Transmit a mesh hierarchy
    {
      MeshHierarchyType mesh_hierarchy;
      MeshType mesh = mesh_hierarchy.root();
      viennagrid::io::vtk_reader<MeshType> reader;
      reader(mesh, "../../data/tets_with_data_main.pvd");
      viennagrid::mpi::send(mesh_hierarchy, 1, MPI_COMM_WORLD);
    }
    // Test 2
    // Transmit a serialized mesh
    {
      MeshHierarchyType mesh_hierarchy;
      MeshType mesh = mesh_hierarchy.root();
      viennagrid::io::vtk_reader<MeshType> reader;
      reader(mesh, "../../data/tets_with_data_main.pvd");
      viennagrid_serialized_mesh_hierarchy serialized_mesh;
      viennagrid_serialized_mesh_hierarchy_make(&serialized_mesh);
      mesh_hierarchy.serialize(serialized_mesh, true);
      viennagrid::mpi::send(serialized_mesh, 1, MPI_COMM_WORLD);
      viennagrid_serialized_mesh_hierarchy_delete(serialized_mesh);
    }
  }
  else
  {
    // Test 1
    // Receive a mesh hierarchy
    {
      MeshHierarchyType mesh_hierarchy;
      viennagrid::mpi::recv(mesh_hierarchy, 0, MPI_COMM_WORLD);
      viennagrid::io::vtk_writer<MeshType> writer;
      MeshType mesh = mesh_hierarchy.root();
      writer(mesh, "received_hierarchy");
    }
    // Test 2
    // Receive a serialized mesh
    {
      MeshHierarchyType mesh_hierarchy;
      viennagrid_serialized_mesh_hierarchy serialized_mesh;
      viennagrid_serialized_mesh_hierarchy_make(&serialized_mesh);
      viennagrid::mpi::recv(serialized_mesh, 0, MPI_COMM_WORLD);
      mesh_hierarchy.deserialize(serialized_mesh);
      viennagrid::io::vtk_writer<MeshType> writer;
      MeshType mesh = mesh_hierarchy.root();
      writer(mesh, "received_serialized_mesh");
      viennagrid_serialized_mesh_hierarchy_delete(serialized_mesh);
    }
  }


  MPI_Finalize();
  return 0;
}
