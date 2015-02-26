#ifndef VIENNAGRID_MPI_HPP
#define VIENNAGRID_MPI_HPP

// System includes
#include <mpi.h>



namespace viennagrid {
namespace mpi {

#define MAX_REGION_NAME_LENGTH 100

enum tags{
  SCALAR_VALUES,
  HOLE_POINTS_OFFSETS,
  CELL_VERTEX_OFFSETS,
  CELL_VERTICES,
  CELL_PARENTS,
  CELL_REGION_OFFSETS,
  CELL_REGIONS,
  MESH_PARENTS,
  MESH_CELL_COUNT,
  MESH_CELLS,
  REGION_IDS,
  REGION_NAMES,
  POINTS,
  HOLE_POINTS,
  CELL_ELEMENT_TAGS
};

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

/**
 * Sends a ViennaGrid serialized mesh over the MPI in a blocking manner
 * @param serialized_mesh A ViennaGrid serialized mesh (which is a pointer) to be transmitted via MPI
 * @param destination The destination MPI process ID
 * @param comm The MPI communicator
 * @return void
 *
 */
void send(viennagrid_serialized_mesh_hierarchy serialized_mesh, int destination, MPI_Comm comm)
{
  std::vector<int> scalar_values(8);
  scalar_values[0] = serialized_mesh->geometric_dimension;
  scalar_values[1] = serialized_mesh->vertex_count;
  scalar_values[2] = serialized_mesh->hole_point_element_count;
  scalar_values[3] = serialized_mesh->cell_count;
  scalar_values[4] = serialized_mesh->cell_dimension;
  scalar_values[5] = serialized_mesh->mesh_count;
  scalar_values[6] = serialized_mesh->region_count;
  scalar_values[7] = serialized_mesh->data_owned;
  MPI_Send(&scalar_values.front(), scalar_values.size(), MPI_INT, destination, SCALAR_VALUES, MPI_COMM_WORLD);

  MPI_Send(serialized_mesh->cell_vertex_offsets, serialized_mesh->cell_count+1, MPI_INT, destination, CELL_VERTEX_OFFSETS, MPI_COMM_WORLD);
  MPI_Send(serialized_mesh->cell_vertices,       serialized_mesh->cell_vertex_offsets[serialized_mesh->cell_count], MPI_INT, destination, CELL_VERTICES, MPI_COMM_WORLD);
  MPI_Send(serialized_mesh->cell_parents,        serialized_mesh->cell_count, MPI_INT, destination, CELL_PARENTS, MPI_COMM_WORLD);
  MPI_Send(serialized_mesh->cell_region_offsets, serialized_mesh->cell_count+1, MPI_INT, destination, CELL_REGION_OFFSETS, MPI_COMM_WORLD);
  MPI_Send(serialized_mesh->cell_regions,        serialized_mesh->cell_region_offsets[serialized_mesh->cell_count], MPI_INT, destination, CELL_REGIONS, MPI_COMM_WORLD);
  MPI_Send(serialized_mesh->mesh_parents,        serialized_mesh->mesh_count, MPI_INT, destination, MESH_PARENTS, MPI_COMM_WORLD);
  MPI_Send(serialized_mesh->mesh_cell_count,     serialized_mesh->mesh_count, MPI_INT, destination, MESH_CELL_COUNT, MPI_COMM_WORLD);
  for(int i = 0; i < serialized_mesh->mesh_count; i++)
    MPI_Send(serialized_mesh->mesh_cells[i],     serialized_mesh->mesh_cell_count[i], MPI_INT, destination, MESH_CELLS, MPI_COMM_WORLD);
  MPI_Send(serialized_mesh->region_ids,          serialized_mesh->region_count, MPI_INT, destination, REGION_IDS, MPI_COMM_WORLD);
  for(int i = 0; i < serialized_mesh->region_count; i++)
  {
    std::string temp_name(serialized_mesh->region_names[i]);
    MPI_Send(&temp_name[0], temp_name.size()+1, MPI_CHAR, destination, REGION_NAMES, MPI_COMM_WORLD);
  }

  MPI_Send(serialized_mesh->points,      serialized_mesh->vertex_count*serialized_mesh->geometric_dimension, MPI_DOUBLE, destination, POINTS, MPI_COMM_WORLD);

  if(serialized_mesh->hole_point_element_count > 0)
  {
    MPI_Send(serialized_mesh->hole_points_offsets, serialized_mesh->hole_point_element_count+1, MPI_INT, destination, HOLE_POINTS_OFFSETS, MPI_COMM_WORLD);
    MPI_Send(serialized_mesh->hole_points, serialized_mesh->hole_points_offsets[serialized_mesh->hole_point_element_count], MPI_DOUBLE, destination, HOLE_POINTS, MPI_COMM_WORLD);
  }

  MPI_Send(serialized_mesh->cell_element_tags, serialized_mesh->cell_count, MPI_UNSIGNED_CHAR, destination, CELL_ELEMENT_TAGS, MPI_COMM_WORLD);
}

/**
 * Sends a ViennaGrid mesh hierarchy over the MPI in a blocking manner
 * @param mesh_hierarchy A ViennaGrid mesh hierarchy to be transmitted via MPI
 * @param destination The destination MPI process ID
 * @param comm The MPI communicator
 * @return void
 *
 */
void send(viennagrid::mesh_hierarchy_t& mesh_hierarchy, int destination, MPI_Comm comm)
{
  viennagrid_serialized_mesh_hierarchy serialized_mesh;
  viennagrid_serialized_mesh_hierarchy_make(&serialized_mesh);
  mesh_hierarchy.serialize(serialized_mesh, true);

  send(serialized_mesh, destination, comm);

  viennagrid_serialized_mesh_hierarchy_delete(serialized_mesh);
}

/**
 * Receives a ViennaGrid serialized mesh over the MPI in a blocking manner
 * @param serialized_mesh A ViennaGrid serialized mesh (which is a pointer) to be populated via MPI
 * @param source The source MPI process ID
 * @param comm The MPI communicator
 * @return void
 *
 */
void recv(viennagrid_serialized_mesh_hierarchy serialized_mesh, int source, MPI_Comm comm)
{
  std::vector<int> scalar_values(8);
  MPI_Recv(&scalar_values.front(), scalar_values.size(), MPI_INT, source, SCALAR_VALUES, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  serialized_mesh->geometric_dimension      = scalar_values[0];
  serialized_mesh->vertex_count             = scalar_values[1];
  serialized_mesh->hole_point_element_count = scalar_values[2];
  serialized_mesh->cell_count               = scalar_values[3];
  serialized_mesh->cell_dimension           = scalar_values[4];
  serialized_mesh->mesh_count               = scalar_values[5];
  serialized_mesh->region_count             = scalar_values[6];
  serialized_mesh->data_owned               = scalar_values[7];

  viennagrid_new((serialized_mesh->cell_count+1)*sizeof(int), (void**)&(serialized_mesh->cell_vertex_offsets));
  MPI_Recv(serialized_mesh->cell_vertex_offsets, serialized_mesh->cell_count+1, MPI_INT, source, CELL_VERTEX_OFFSETS, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  viennagrid_new((serialized_mesh->cell_vertex_offsets[serialized_mesh->cell_count])*sizeof(int), (void**)&(serialized_mesh->cell_vertices));
  MPI_Recv(serialized_mesh->cell_vertices, serialized_mesh->cell_vertex_offsets[serialized_mesh->cell_count], MPI_INT, source, CELL_VERTICES, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  viennagrid_new((serialized_mesh->cell_count)*sizeof(int), (void**)&(serialized_mesh->cell_parents));
  MPI_Recv(serialized_mesh->cell_parents, serialized_mesh->cell_count, MPI_INT, source, CELL_PARENTS, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  viennagrid_new((serialized_mesh->cell_count+1)*sizeof(int), (void**)&(serialized_mesh->cell_region_offsets));
  MPI_Recv(serialized_mesh->cell_region_offsets, serialized_mesh->cell_count+1, MPI_INT, source, CELL_REGION_OFFSETS, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  viennagrid_new((serialized_mesh->cell_region_offsets[serialized_mesh->cell_count])*sizeof(int), (void**)&(serialized_mesh->cell_regions));
  MPI_Recv(serialized_mesh->cell_regions, serialized_mesh->cell_region_offsets[serialized_mesh->cell_count], MPI_INT, source, CELL_REGIONS, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  viennagrid_new((serialized_mesh->mesh_count)*sizeof(int), (void**)&(serialized_mesh->mesh_parents));
  MPI_Recv(serialized_mesh->mesh_parents, serialized_mesh->mesh_count, MPI_INT, source, MESH_PARENTS, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  viennagrid_new((serialized_mesh->mesh_count)*sizeof(int), (void**)&(serialized_mesh->mesh_cell_count));
  MPI_Recv(serialized_mesh->mesh_cell_count, serialized_mesh->mesh_count, MPI_INT, source, MESH_CELL_COUNT, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  viennagrid_new((serialized_mesh->mesh_count)*sizeof(int*), (void**)&(serialized_mesh->mesh_cells));
  for(int i = 0; i < serialized_mesh->mesh_count; i++)
  {
    viennagrid_new((serialized_mesh->mesh_cell_count[i])*sizeof(int), (void**)&(serialized_mesh->mesh_cells[i]));
    MPI_Recv(serialized_mesh->mesh_cells[i], serialized_mesh->mesh_cell_count[i], MPI_INT, source, MESH_CELLS, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

  viennagrid_new((serialized_mesh->region_count)*sizeof(int), (void**)&(serialized_mesh->region_ids));
  MPI_Recv(serialized_mesh->region_ids, serialized_mesh->region_count, MPI_INT, source, REGION_IDS, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  viennagrid_new((serialized_mesh->region_count)*sizeof(const char *), (void**)&(serialized_mesh->region_names));
  for(int i = 0; i < serialized_mesh->region_count; i++)
  {
    char* temp_name;
    viennagrid_new((MAX_REGION_NAME_LENGTH)*sizeof(char), (void**)&temp_name);
    MPI_Status status;
    MPI_Recv(temp_name, MAX_REGION_NAME_LENGTH, MPI_CHAR, source, REGION_NAMES, MPI_COMM_WORLD, &status);
    int name_length = 0;
    MPI_Get_count(&status, MPI_CHAR, &name_length);
    char* temp_name_cut;
    viennagrid_new((name_length)*sizeof(char), (void**)&temp_name_cut);
    std::copy(temp_name, temp_name+name_length, temp_name_cut);
    serialized_mesh->region_names[i] = (const char *)temp_name_cut;
    viennagrid_delete((void**)&temp_name);
  }

  viennagrid_new((serialized_mesh->vertex_count*serialized_mesh->geometric_dimension)*sizeof(double), (void**)&(serialized_mesh->points));
  MPI_Recv(serialized_mesh->points, serialized_mesh->vertex_count*serialized_mesh->geometric_dimension, MPI_DOUBLE, source, POINTS, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  if(serialized_mesh->hole_point_element_count > 0)
  {
    viennagrid_new((serialized_mesh->hole_point_element_count+1)*sizeof(int), (void**)&(serialized_mesh->hole_points_offsets));
    MPI_Recv(serialized_mesh->hole_points_offsets, serialized_mesh->hole_point_element_count+1, MPI_INT, source, HOLE_POINTS_OFFSETS, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    viennagrid_new((serialized_mesh->hole_points_offsets[serialized_mesh->hole_point_element_count])*sizeof(double), (void**)&(serialized_mesh->hole_points));
    MPI_Recv(serialized_mesh->hole_points, serialized_mesh->hole_points_offsets[serialized_mesh->hole_point_element_count], MPI_DOUBLE, source, HOLE_POINTS, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

  viennagrid_new((serialized_mesh->cell_count)*sizeof(unsigned char), (void**)&(serialized_mesh->cell_element_tags));
  MPI_Recv(serialized_mesh->cell_element_tags, serialized_mesh->cell_count, MPI_UNSIGNED_CHAR, source, CELL_ELEMENT_TAGS, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

/**
 * Receives a ViennaGrid mesh hierarchy over the MPI in a blocking manner
 * @param serialized_mesh A ViennaGrid mesh hierarchy to be populated via MPI
 * @param source The source MPI process ID
 * @param comm The MPI communicator
 * @return void
 *
 */
void recv(viennagrid::mesh_hierarchy_t& mesh_hierarchy, int source, MPI_Comm comm)
{
  viennagrid_serialized_mesh_hierarchy serialized_mesh;
  viennagrid_serialized_mesh_hierarchy_make(&serialized_mesh);

  recv(serialized_mesh, source, comm);

  mesh_hierarchy.deserialize(serialized_mesh);
  viennagrid_serialized_mesh_hierarchy_delete(serialized_mesh);
}


} // mpi
} // viennagrid

#endif
