#ifndef MESH_PREPROCESSOR_H
#define MESH_PREPROCESSOR_H

#include <vector>
#include "elements.h"
#include "mpi_env.h"
#include "gmsh_reader.h"

using namespace std;


void preprocess( std::string filename, std::string file_base, int nSubparts );
void gather_boundary_hpaths( gmsh_mesh gmsh_local, MPI_env mpi_env, int nSubparts, std::vector<int> mpi2global, std::vector<std::vector<std::vector<int>>> all_hpaths );
void gather_interior_hpaths( std::vector<gmsh_mesh> gmsh_parts, MPI_env mpi_env, std::vector<std::vector<int>> int2old, std::vector<int> mpi2global, std::vector<std::vector<std::vector<int>>> all_hpaths );
void update_mesh_cells( gmsh_mesh *mesh, std::vector<std::vector<int>> local2hpath );
void output_mesh_files( std::string file_base, std::vector<gmsh_mesh> gmsh_parts );
void update_node_ordering( MPI_env mpi_env, int ipart, gmsh_mesh *mesh );
std::vector<std::vector<int>> get_local2hpath_mapping( std::vector<gmsh_mesh> gmsh_parts );
void print_hpaths( std::string filebase, std::vector<gmsh_mesh> gmsh_parts );
void print_mesh_cells ( std::string file_base, gmsh_mesh *gmsh );
void tecplot_cell_ids( char* filename, gmsh_mesh gmsh );

#endif
