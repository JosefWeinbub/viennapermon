#ifndef VIENNAGRID_MPI_HPP
#define VIENNAGRID_MPI_HPP

// System includes
#include <mpi.h>



namespace viennagrid {
namespace mpi {

struct serializedmesh_operator
{
public:

  serializedmesh_operator(viennagrid_serialized_mesh_hierarchy serialized_mesh, MPI_Datatype* mpi_datatype, int* blocklen, int blocklen_size) :
    serialized_mesh_(serialized_mesh), commited_(false), delete_blocklen_(false)
  {
    typesize_ = blocklen_size;
    register_datatype(mpi_datatype, blocklen, blocklen_size);
  }

  serializedmesh_operator(viennagrid_serialized_mesh_hierarchy serialized_mesh, MPI_Datatype* mpi_datatype) :
    serialized_mesh_(serialized_mesh), commited_(false), delete_blocklen_(false)
  {
//    if(serialized_mesh->cell_count == 0) return;

    int mesh_cells_size = 0;
    for(int i = 0; i < serialized_mesh->mesh_count; i++)
      mesh_cells_size += serialized_mesh->mesh_cell_count[i];

    int int_len =
      1+ // geometric_dimension
      1+ // vertex_count
      1+ // hole_point_element_count
      serialized_mesh->hole_point_element_count+1+ // size(hole_points_offsets)
      1+ // cell_count
      1+ // cell_dimension
      serialized_mesh->cell_count+1+ // size(cell_vertex_offsets)
      serialized_mesh->cell_vertex_offsets[serialized_mesh->cell_count]+ // size(cell_vertices)
      serialized_mesh->cell_count+ // size(cell_parents)
      serialized_mesh->cell_count+1+ // size(cell_region_offsets)
      serialized_mesh->cell_region_offsets[serialized_mesh->cell_count]+ // size(cell_regions)
      1+ // mesh_count
      serialized_mesh->mesh_count+ // size(mesh_parents)
      serialized_mesh->mesh_count+ // size(mesh_cell_count)
      mesh_cells_size+ // size(mesh_cells*)
      1+ // region_count
      serialized_mesh->region_count+// size(region_ids)
      1; // data_owned

    int double_len = serialized_mesh->vertex_count*serialized_mesh->geometric_dimension; // size(points)
    if(serialized_mesh->hole_point_element_count > 0)
      double_len += serialized_mesh->hole_points_offsets[serialized_mesh->hole_point_element_count]; // size(hole_points)

    int unsigned_char_len = serialized_mesh->cell_count; // size(cell_element_tags)

    int char_len = 0;
    for(int i = 0; i < serialized_mesh->region_count; i++)
      char_len += strlen(serialized_mesh->region_names[i])+1;

//    typesize_  = 4;
//    blocklen_ = new int[typesize_];
//    delete_blocklen_ = true;
//    blocklen_[0] = int_len;
//    blocklen_[1] = double_len;
//    blocklen_[2] = unsigned_char_len;
//    blocklen_[3] = char_len;

    typesize_  = 1;
    blocklen_ = new int[typesize_];
    delete_blocklen_ = true;
    blocklen_[0] = double_len;


    register_datatype(mpi_datatype, blocklen_, typesize_);
  }

  ~serializedmesh_operator()
  {
//    if(commited_)
//    {
//      MPI_Type_free(&get_mpi_type());
//      delete disp_;
//      if(delete_blocklen_) delete blocklen_;
//    }
  }

private:
  void register_datatype(MPI_Datatype* mpi_datatype, int* blocklen, int blocklen_size)
  {
    std::cout << blocklen_size << " blocklen size " << std::endl;
//    std::cout << blocklen[0] << " ints" << std::endl;
//    std::cout << blocklen[1] << " doubles" << std::endl;
//    std::cout << blocklen[2] << " unsigned char" << std::endl;
//    std::cout << blocklen[3] << " char" << std::endl;

    std::cout << blocklen[0] << " doubles" << std::endl;

//    MPI_Datatype  types[]    = {MPI_INT, MPI_DOUBLE, MPI_UNSIGNED_CHAR, MPI_CHAR};
//    disp_ = new MPI_Aint[blocklen_size];
//    MPI_Address(&(serialized_mesh_->geometric_dimension),   &disp_[0]);
//    MPI_Address(&(serialized_mesh_->points),                &disp_[1]);
//    MPI_Address(&(serialized_mesh_->cell_element_tags),     &disp_[2]);
//    MPI_Address(&(serialized_mesh_->region_names),          &disp_[3]);

    MPI_Datatype  types[]    = {MPI_DOUBLE};
    disp_ = new MPI_Aint[blocklen_size];
    MPI_Address(&(serialized_mesh_->points),                &disp_[0]);


    for(int i = blocklen_size-1; i >= 0; i--)
      disp_[i] -= disp_[0];

    MPI_Type_struct(blocklen_size, blocklen, disp_, types, mpi_datatype);
    MPI_Type_commit(mpi_datatype);
    commited_ = true;
  }


public:
//  MPI_Datatype& get_mpi_type()
//  {
//    return MPI_VIENNAGRID_SERIALIZEDMESH_;
//  }

  viennagrid_serialized_mesh_hierarchy& get_mesh()
  {
    return serialized_mesh_;
  }

  int* get_blocklen()
  {
    return blocklen_;
  }

  int& get_blocklen_size()
  {
    return typesize_;
  }


  viennagrid_serialized_mesh_hierarchy& serialized_mesh_;
//  MPI_Datatype  MPI_VIENNAGRID_SERIALIZEDMESH_;
  MPI_Aint*     disp_;
  bool          commited_;
  int           typesize_;
  int*          blocklen_;
  bool          delete_blocklen_;
};


} // mpi
} // viennagrid

#endif
