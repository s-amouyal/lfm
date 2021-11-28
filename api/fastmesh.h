/**
 * Copyright (C) Huawei Technologies Co., Ltd. 2020. ALL RIGHTS RESERVED.
 *
 * See file LICENSE for terms.
 *
 *
 * (lib)FastMesh
 * =============
 *
 * This is a prototype for a new memory layout for representing a mesh, intended
 * to facilitate faster iterative solvers. This is achieved by enhanced memory
 * proximity and a layout which is friendly to both SIMD instructions and
 * network offloading APIs (i.e. RDMA).
 */

#ifndef NEIGHBOR_H_
#define NEIGHBOR_H_

#include <math.h>
#include <vector>
#include <memory>
#include <ostream>
#include <iostream>
#include <functional>

#include "gmsh_reader.h"
#include "elements.h"
#include "mpi_env.h"
#include "cfdv0_solver.h"

namespace fastmesh {

//***************************************************************************************************
// List of solvers
//***************************************************************************************************
typedef enum fastmesh_solvers {
    SOLVER_CAAFOAM = 0,
} fastmesh_solver_t;

//***************************************************************************************************
// The Mesh class focuses on data exchange between Submesh instances (while the Submesh focuses on
// solvers). This means the mesh is responsible for moving the contents of boundary cells between
// Submeshes - after every iteration. This exchange can be either within the same host or over the
// network.
//***************************************************************************************************
template <typename PRECISION, unsigned DIM_CNT>
class Mesh {

public:
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    /**< Create an empty mesh */
    Mesh() {};

    /**< Clone an existing mesh */
    Mesh( Mesh<PRECISION, DIM_CNT> *original) {};

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Geometrical and Solver objects
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    std::vector<gmsh_mesh> gmsh_parts;
    std::vector<var_faced_solver> submesh_cfd;

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // MPI environment
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MPI_env mpi_env;
    std::vector<int> mpi_neighbors;

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Functions
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    void initialize( std::string filebase );

	void initialize_solver( fastmesh_solver_t solver );
	void output_solution( PRECISION time, PRECISION PRINT_FILE_FREQUENCY, std::string file_base );

	void solve();

	// Wrap-up simulation, including MPI environmnet
	void finalize(){
		mpi_env.finalize();
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// RK variables (temporary)
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	int num_RK_steps = 1;

	void read_mesh_files( MPI_env mpi_env, std::string file_base );
	void read_mesh_faces( ifstream &fin, int n_faces, gmsh_mesh *mesh );
	void read_mesh_cells( ifstream &fin, int n_cells, gmsh_mesh *mesh );
	void read_mesh_nodes( ifstream &fin, int n_nodes, gmsh_mesh *mesh );
	void set_params( int ipart, gmsh_mesh *mesh, std::vector<gmsh_mesh> mesh_parts );
	void calculate_v_neigh( int ipart, gmsh_mesh *mesh, std::vector<gmsh_mesh> mesh_parts, PRECISION domain[DIM_CNT] );
};

} // namespace fastmesh

#endif // NEIGHBOR_H_








































