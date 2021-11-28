#ifndef GMSH_READER_H
#define GMSH_READER_H

#include <vector>
#include "elements.h"
#include "mpi_env.h"

using namespace std;

//***************************************************************************************************
// Functions
//***************************************************************************************************

void read_boundaries( ifstream &fin, unsigned int &num_boundaries  );
void read_entities  ( ifstream &fin, unsigned int num_boundaries, vector<vector<vector<int>>> &entity2boundary );
void skip_lines     ( ifstream *fin, int num_skip_lines );
void get_periodic_boundary_conditions( MPI_env mpi_env,
									   ifstream &fin,
									   vector<vector<vector<int>>> &entity2boundary,
									   std::vector<vector<int>>    &periodic_nodes);

// New gmsh
void read_nodes( ifstream &fin,
				 vector<vector<vector<int>>> &entity2boundary,
				 std::vector<std::vector<int>> &nodes_bound,
				 gmsh_mesh *mesh );

void read_cells( ifstream  &fin,
				 int       &n_neigh,
				 gmsh_mesh *mesh,
				 std::vector<std::vector<gmsh_neighbor>> &gmsh_neighs );

void init_faces( int n_neigh,
				 std::vector<std::vector<int>> nodes_bound,
				 std::vector<vector<int>>     &periodic_nodes,
				 std::vector<int>             &periodic_faces,
				 gmsh_mesh                    *mesh );

void reorder_cell_local_nodes( int        n_neigh,
							   gmsh_mesh *mesh );

void set_cells_faces( int        n_neigh,
					  gmsh_mesh *mesh );

void reorder_cell_local_faces( int        n_neigh,
							   gmsh_mesh *mesh );

void set_boundary_conditions( int n_neigh,
							  gmsh_mesh *mesh,
							  std::vector<std::vector<int>> nodes_bound,
							  std::vector<int> &periodic_faces,
							  std::vector<std::vector<gmsh_neighbor>> &gmsh_neighs );

void read_gmsh_file( MPI_env mpi_env, std::string filename, gmsh_mesh *mesh );

void set_cell_neighbors( gmsh_mesh *mesh, std::vector<std::vector<gmsh_neighbor>> gmsh_neighs );
void compute_mesh_parameters( MPI_env mpi_env, gmsh_mesh *mesh );

void init_mesh_size_mpi_type( MPI_env &mpi_env );

MPI_Datatype get_mesh_mpi_type( MPI_env &mpi_env, gmsh_mesh *mesh );

MPI_Datatype get_mpitype_meshnode ( gmsh_mesh *mesh );
MPI_Datatype get_mpitype_meshface ( gmsh_mesh *mesh );
MPI_Datatype get_mpitype_meshcell ( gmsh_mesh *mesh );
MPI_Datatype get_mpitype_meshneigh( gmsh_mesh *mesh );

template < typename T >
T* from_2Dvec_to_1Darray( int nrows, int ncols, std::vector<std::vector<T>> vecT );

template < typename T >
std::vector<std::vector<T>> from_1Darray_to_2Dvec( int nrows, int ncols, T *arrayT );

#endif








































