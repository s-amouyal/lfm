#ifndef MESH_PARTITIONING_H
#define MESH_PARTITIONING_H

#include <vector>
#include <metis.h>
#include "elements.h"
#include "mpi_env.h"
#include "gmsh_reader.h"

//***************************************************************************************************
// Separates the boundary cells from original mesh
// 		- returns a new t_mesh containing only the interior cells
//***************************************************************************************************
gmsh_mesh set_interior_mesh( gmsh_mesh        *gmsh_local,
							 MPI_env           mpi_env,
							 std::vector<int> &old2new,
							 std::vector<int> &old2new_bnd,
							 std::vector<int> &int2old );

//***************************************************************************************************
// Partition mesh based on number of parts
// 		- only library with direct contact with libmetis
//***************************************************************************************************
void metis_partition_mesh( gmsh_mesh *gmsh,						// Original mesh
						   int nParts,							// Number of metis partitions
						   std::vector<int> &partitions);		// (Output) Maps cells to partitions

void metis_partition_mesh_new( gmsh_mesh *gmsh,						// Original mesh
						   int nParts,							// Number of metis partitions
						   std::vector<int> &partitions);		// (Output) Maps cells to partitions


//***************************************************************************************************
// Initialize the parameters, indices and connectivity of every mesh partition
// 		- Returns a vector of t_mesh of size (nParts+1): nParts partitions + boundary hpath
//***************************************************************************************************
std::vector<gmsh_mesh> init_submesh_partitions( gmsh_mesh       *gmsh,
												gmsh_mesh       *gmsh_int,
												int              nParts,
												std::vector<int> partitions,
												std::vector<std::vector<int>> &int2old );

//***************************************************************************************************
// Partition mesh based on number of parts
// 		- only library with direct contact with libmetis
//***************************************************************************************************
gmsh_mesh partition_mesh_metis_mpi( gmsh_mesh        *gmsh,					// Original mesh
									MPI_env           mpi_env,				// MPI environment
									std::vector<int> &mpi2global );			// Index mapping

//***************************************************************************************************
// Partition mesh based on number of parts
// 		- only library with direct contact with libmetis
//***************************************************************************************************
std::vector<gmsh_mesh> partition_local_mesh( gmsh_mesh                     *gmsh_local,
													gmsh_mesh                     *gmsh_int,
													int                            nSubparts,
													std::vector<int>               old2parts,
													MPI_env		                   mpi_env,
													std::vector<int>               tmp2old,
													std::vector<std::vector<int>> &int2old );

//***************************************************************************************************
// Initialize the parameters, indices and connectivity of every mesh partition
// 		- Returns a vector of t_mesh of size (nParts+1): nParts partitions + boundary hpath
//***************************************************************************************************
gmsh_mesh init_submesh_partitions_mpi( gmsh_mesh        *gmsh,
									   MPI_env           mpi_env,
									   std::vector<int>  partitions,
									   std::vector<int> &old2parts );

//***************************************************************************************************
//***************************************************************************************************
void update_mesh_mpi_connectivity( gmsh_mesh *mesh_local,
								   gmsh_mesh *mesh_bnd );

//***************************************************************************************************
// Update periodic's neighbors IDs after boundary Hpath
// 	- only serial for now
//***************************************************************************************************
void update_periodic_connectivity( gmsh_mesh *mesh_local,
								   gmsh_mesh *mesh_bnd );

//***************************************************************************************************
//***************************************************************************************************
gmsh_mesh partition_mesh_mpi_efficient( MPI_env mpi_env,
										gmsh_mesh *gmsh_domain );

//***************************************************************************************************
//***************************************************************************************************
gmsh_mesh node_local_mesh_partitioning( MPI_env          mpi_env,
										gmsh_mesh       *gmsh_domain,
										std::vector<int> partitions_mesh );

//***************************************************************************************************
//***************************************************************************************************
struct t_mesh_size get_mesh_size( gmsh_mesh *mesh );
struct t_mesh_size get_mesh_size( gmsh_mesh *mesh, int part_num, std::vector<int> partitions );

//***************************************************************************************************
//***************************************************************************************************
gmsh_mesh get_mesh_partition( MPI_env mpi_env,
							  gmsh_mesh *mesh_domain,
							  int part_num,
							  std::vector<int> partitions,
							  std::vector<int> &old2new_faces );

//***************************************************************************************************
//***************************************************************************************************
gmsh_mesh allocate_mesh( struct t_mesh_size mesh_size );

//***************************************************************************************************
//***************************************************************************************************
void update_mpi_connectivity( MPI_env    mpi_env,
							  gmsh_mesh *gmsh_domain,
							  gmsh_mesh *gmsh_local );

void update_mpi_connectivity_faces( MPI_env    mpi_env,
									gmsh_mesh *gmsh_domain,
									gmsh_mesh *gmsh_local );
#endif






























