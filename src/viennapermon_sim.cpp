
// System includes
#include <iostream>
#include <mpi.h>

// ViennaGrid includes
#include "viennagrid/viennagrid.hpp"
#include "viennagrid/io/vtk_reader.hpp"
#include "viennagrid/io/vtk_writer.hpp"

// ViennaMesh includes
#include "viennameshpp/algorithm_pipeline.hpp"

// Local includes
#include "viennagrid_mpi.hpp"
#include "viennagrid2libmesh.hpp"



void process_local(viennamesh::context_handle & context,
                   int mpi_rank,
                   viennagrid::mesh mesh_local)
{
  viennamesh::algorithm_handle mark_hull_regions = context.make_algorithm("mark_hull_regions");
  mark_hull_regions.set_input( "mesh", mesh_local.internal() );
  {
    viennamesh::LoggingStack s("mark_hull_regions");
    mark_hull_regions.run();
  }

  // generate local volume mesh
  //
  viennamesh::algorithm_handle mesher = context.make_algorithm("netgen_make_mesh");
  mesher.set_default_source(mark_hull_regions);
  {
    viennamesh::LoggingStack s("netgen_make_mesh");
    mesher.run();
  }

  viennamesh::algorithm_handle mesh_writer = context.make_algorithm("mesh_writer");
  mesh_writer.set_default_source(mesher);
  mesh_writer.set_input( "filename", "mpi_mesh_" + boost::lexical_cast<std::string>(mpi_rank) + "_mesh.vtu" );
  {
    viennamesh::LoggingStack s("mesh_writer");
    mesh_writer.run();
  }

//  libMesh::SerialMesh libmesh;
//  viennagrid2libmesh(mesher.get_output<viennagrid_mesh>("mesh")(), libmesh);
//  libmesh.print_info();



}

int main(int argc, char* argv[])
{
//  libMesh::LibMeshInit init (argc, argv);

  int mpi_size, mpi_rank;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

  typedef viennagrid::mesh MeshType;

  if(argc != 2)
  {
    if(mpi_rank == 0)
      std::cout << "Missing parameter - Usage: " << argv[0] << " input_meshfile.poly" << std::endl;
    MPI_Finalize();
    return 1;
  }



  if(mpi_size < 2)
  {
    std::cout << "This application must be executed with at least 2 MPI processes - aborting!" << std::endl;
    MPI_Finalize();
    return 1;
  }

  viennamesh::context_handle context;

  // Master process
  if(mpi_rank == 0)
  {
    std::string filename = argv[1];

    viennamesh::algorithm_handle mesh_reader = context.make_algorithm("plc_reader");
    mesh_reader.set_input( "filename", filename );
    {
      viennamesh::LoggingStack s("plc_reader");
      mesh_reader.run();
    }

    viennamesh::algorithm_handle mesher = context.make_algorithm("tetgen_make_mesh");
    mesher.set_default_source(mesh_reader);
    mesher.set_input("cell_size", 10.0);
    {
      viennamesh::LoggingStack s("tetgen_make_mesh");
      mesher.run();
    }

    viennamesh::algorithm_handle metis_partitioning = context.make_algorithm("metis_mesh_partitioning");
    metis_partitioning.set_default_source(mesher);
    metis_partitioning.set_input( "region_count", mpi_size );
    {
      viennamesh::LoggingStack s("metis_mesh_partitioning");
      metis_partitioning.run();
    }

    viennamesh::algorithm_handle extract_boundary = context.make_algorithm("extract_boundary");
    extract_boundary.set_default_source(metis_partitioning);
    extract_boundary.run();

    viennamesh::algorithm_handle extract_plc_geometry = context.make_algorithm("extract_plc_geometry");
    extract_plc_geometry.set_default_source(extract_boundary);
    extract_plc_geometry.set_input("coplanar_tolerance", 1e-2);
    extract_plc_geometry.set_input("colinear_tolerance", 1e-2);
    {
      viennamesh::LoggingStack s("extract_plc_geometry");
      extract_plc_geometry.run();
    }

    viennamesh::algorithm_handle triangle_make_hull = context.make_algorithm("triangle_make_hull");
    triangle_make_hull.set_default_source(extract_plc_geometry);
    triangle_make_hull.set_input("cell_size", 1.0);
    {
      viennamesh::LoggingStack s("triangle_make_hull");
      triangle_make_hull.run();
    }

    viennamesh::algorithm_handle mark_hull_regions = context.make_algorithm("mark_hull_regions");
    mark_hull_regions.set_default_source(triangle_make_hull);
    {
      viennamesh::LoggingStack s("mark_hull_regions");
      mark_hull_regions.run();
    }

    viennamesh::algorithm_handle split_mesh = context.make_algorithm("split_mesh");
    split_mesh.set_default_source(mark_hull_regions);
    {
      viennamesh::LoggingStack s("split_mesh");
      split_mesh.run();
    }

    // Distribute the partitions
    // partition-0 is processed by the master process (this process, mpi_rank=0)
    // partitions-1..n are processed by the worker processes (mpi_rank=1..mpi_size-1)
    //

//     viennamesh::data_handle<viennagrid_mesh> meshes = split_mesh.get_output<viennagrid_mesh>( "mesh" );
    for(int target = 1; target < mpi_size; target++)
    {
      viennagrid::mesh mesh = split_mesh.get_output<viennagrid_mesh>("mesh[" + boost::lexical_cast<std::string>(target) + "]")();
      viennagrid::mpi::send(mesh, target, MPI_COMM_WORLD);
    }

    // Master simulates the first partition
    //

//     viennagrid::mesh mesh = meshes(0);
    viennagrid::mesh mesh = split_mesh.get_output<viennagrid_mesh>("mesh[0]")();
    process_local(context, mpi_rank, mesh);
  }
  // Worker processes
  else
  {
    // Receive the local mesh which is to be processed/simulated
    //
    MeshType mesh;
    viennagrid::mpi::recv(mesh, 0, MPI_COMM_WORLD);

    // Process/simulate the local mesh
    //
    process_local(context, mpi_rank, mesh);
  }

  MPI_Finalize();
  return 0;
}
