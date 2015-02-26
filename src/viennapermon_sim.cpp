
// System includes
#include <iostream>
#include <mpi.h>

// ViennaGrid includes
#include "viennagrid/core.hpp"
#include "viennagrid/io/vtk_reader.hpp"
#include "viennagrid/io/vtk_writer.hpp"

// Local includes
#include "viennagrid_mpi.hpp"


void process_local(viennagrid::mesh_hierarchy_t& mesh_hierarchy_local)
{
  // TODO Perform the local FEM-based simulation
  // Every MPI process executes this bit for its own mesh
  //
}

int main(int argc, char* argv[])
{
  int mpi_size, mpi_rank;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

  typedef viennagrid::mesh_hierarchy_t MeshHierarchyType;
  typedef viennagrid::result_of::mesh<MeshHierarchyType>::type MeshType;

  if(argc != 2)
  {
    if(mpi_rank == 0)
      std::cout << "Missing parameter - Usage: " << argv[0] << " input_meshfile.pvd" << std::endl;
    MPI_Finalize();
    return 1;
  }

  if(mpi_size < 2)
  {
    std::cout << "This application must be executed with at least 2 MPI processes - aborting!" << std::endl;
    MPI_Finalize();
    return 1;
  }

  // Master process
  if(mpi_rank == 0)
  {
    // Load the input mesh
    //
    MeshHierarchyType mesh_hierarchy;
    MeshType mesh = mesh_hierarchy.root();
    viennagrid::io::vtk_reader<MeshType> reader;
    reader(mesh, argv[1]);

    // Partition the input mesh - one partition for each MPI process
    //
    std::vector<MeshHierarchyType*> mesh_partitions(mpi_size);
    // ... TODO populate 'mesh_partitions' via a METIS partitioning and ViennaMesh postprocessing
    // for testing, we use/distribute the input mesh to all processes
    for(std::size_t i = 0; i < mesh_partitions.size(); i++)
      mesh_partitions[i] = &mesh_hierarchy;

    // Distribute the partitions
    // partition-0 is processed by the master process (this process, mpi_rank=0)
    // partitions-1..n are processed by the worker processes (mpi_rank=1..mpi_size-1)
    //
    for(int target = 1; target < mpi_size; target++)
    {
      // transfer this partition
      viennagrid::mpi::send(*mesh_partitions[target], target, MPI_COMM_WORLD);
    }

    // Master simulates the first partition
    //
    process_local(*mesh_partitions.front());
  }
  // Worker processes
  else
  {
    // Receive the local mesh which is to be processed/simulated
    //
    MeshHierarchyType mesh_hierarchy;
    viennagrid::mpi::recv(mesh_hierarchy, 0, MPI_COMM_WORLD);

    // Process/simulate the local mesh
    //
    process_local(mesh_hierarchy);
  }

  MPI_Finalize();
  return 0;
}
