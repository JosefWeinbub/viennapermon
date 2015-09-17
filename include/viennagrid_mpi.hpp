#ifndef VIENNAGRID_MPI_HPP
#define VIENNAGRID_MPI_HPP

// System includes
#include <mpi.h>



namespace viennagrid {
namespace mpi {

#define MAX_REGION_NAME_LENGTH 100

enum tags{
  SERIALIZED_MESH_SIZE,
  SERIALIZED_MESH
};


/**
 * Sends a ViennaGrid serialized mesh over the MPI in a blocking manner
 * @param serialized_mesh A ViennaGrid serialized mesh (which is a pointer) to be transmitted via MPI
 * @param destination The destination MPI process ID
 * @param comm The MPI communicator
 * @return void
 *
 */
void send(void * serialized_mesh, viennagrid_int serialized_mesh_size, int destination, MPI_Comm comm)
{
  MPI_Send(&serialized_mesh_size, 1, MPI_INT, destination, SERIALIZED_MESH_SIZE, MPI_COMM_WORLD);
  MPI_Send(serialized_mesh, serialized_mesh_size, MPI_BYTE, destination, SERIALIZED_MESH, MPI_COMM_WORLD);
}

/**
 * Sends a ViennaGrid mesh hierarchy over the MPI in a blocking manner
 * @param mesh_hierarchy A ViennaGrid mesh hierarchy to be transmitted via MPI
 * @param destination The destination MPI process ID
 * @param comm The MPI communicator
 * @return void
 *
 */
void send(viennagrid::mesh & mesh, int destination, MPI_Comm comm)
{
  void * serialized_mesh;
  viennagrid_int serialized_mesh_size;

  mesh.serialize(&serialized_mesh, &serialized_mesh_size);

  send(serialized_mesh, serialized_mesh_size, destination, comm);

  viennagrid_delete(&serialized_mesh);
}

/**
 * Receives a ViennaGrid serialized mesh over the MPI in a blocking manner
 * @param serialized_mesh A ViennaGrid serialized mesh (which is a pointer) to be populated via MPI
 * @param source The source MPI process ID
 * @param comm The MPI communicator
 * @return void
 *
 */
void recv(void ** serialized_mesh, viennagrid_int * serialized_mesh_size, int source, MPI_Comm comm)
{
  MPI_Recv(serialized_mesh_size, 1, MPI_INT, source, SERIALIZED_MESH_SIZE, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  viennagrid_new(*serialized_mesh_size, serialized_mesh);
  MPI_Recv(*serialized_mesh, *serialized_mesh_size, MPI_BYTE, source, SERIALIZED_MESH, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

/**
 * Receives a ViennaGrid mesh hierarchy over the MPI in a blocking manner
 * @param serialized_mesh A ViennaGrid mesh hierarchy to be populated via MPI
 * @param source The source MPI process ID
 * @param comm The MPI communicator
 * @return void
 *
 */
void recv(viennagrid::mesh & mesh, int source, MPI_Comm comm)
{
  viennagrid_int serialized_mesh_size;
  void * serialized_mesh;

  recv(&serialized_mesh, &serialized_mesh_size, source, comm);

  mesh.deserialize(serialized_mesh, serialized_mesh_size);

  viennagrid_delete(&serialized_mesh);
}


} // mpi
} // viennagrid

#endif
