#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <initializer_list>
#include <array>

#include "math.h"
#include "api/cfdv0_solver.h"
#include "api/mesh_writer.h"

using namespace std;

//***************************************************************************************************
//											Templates
//***************************************************************************************************

// Restart simulation from t_start
template void CFDv0_solver<float , 2, 3> :: restart( float  t_start, float  PRINT_FILE_FREQUENCY, int submesh );
template void CFDv0_solver<float , 2, 4> :: restart( float  t_start, float  PRINT_FILE_FREQUENCY, int submesh );
template void CFDv0_solver<double, 2, 3> :: restart( double t_start, double PRINT_FILE_FREQUENCY ,int submesh );
template void CFDv0_solver<double, 2, 4> :: restart( double t_start, double PRINT_FILE_FREQUENCY, int submesh );

// Prepare time-step
template void CFDv0_solver<float , 2, 3> :: prepare_for_timestep( int rk_step );
template void CFDv0_solver<float , 2, 4> :: prepare_for_timestep( int rk_step );
template void CFDv0_solver<double, 2, 3> :: prepare_for_timestep( int rk_step );
template void CFDv0_solver<double, 2, 4> :: prepare_for_timestep( int rk_step );

// Calculate Gradients
template void CFDv0_solver<float , 2, 3> :: calc_gradients( MPI_env &mpi_env );
template void CFDv0_solver<float , 2, 4> :: calc_gradients( MPI_env &mpi_env );
template void CFDv0_solver<double, 2, 3> :: calc_gradients( MPI_env &mpi_env );
template void CFDv0_solver<double, 2, 4> :: calc_gradients( MPI_env &mpi_env );

// AD
template void CFDv0_solver<float , 2, 3> :: calc_AD();
template void CFDv0_solver<float , 2, 4> :: calc_AD();
template void CFDv0_solver<double, 2, 3> :: calc_AD();
template void CFDv0_solver<double, 2, 4> :: calc_AD();

// Viscosity
template void CFDv0_solver<float , 2, 3> :: calc_VIS( MPI_env &mpi_env );
template void CFDv0_solver<float , 2, 4> :: calc_VIS( MPI_env &mpi_env );
template void CFDv0_solver<double, 2, 3> :: calc_VIS( MPI_env &mpi_env );
template void CFDv0_solver<double, 2, 4> :: calc_VIS( MPI_env &mpi_env );

// One RK step
template void CFDv0_solver<float , 2, 3> :: one_rk_step_v1( int rk_step, float  dt, MPI_env &mpi_env );
template void CFDv0_solver<float , 2, 4> :: one_rk_step_v1( int rk_step, float  dt, MPI_env &mpi_env );
template void CFDv0_solver<double, 2, 3> :: one_rk_step_v1( int rk_step, double dt, MPI_env &mpi_env );
template void CFDv0_solver<double, 2, 4> :: one_rk_step_v1( int rk_step, double dt, MPI_env &mpi_env );

// Set boundary conditions (BC)
template void CFDv0_solver<float , 2, 3> :: set_boundary_conditions();
template void CFDv0_solver<float , 2, 4> :: set_boundary_conditions();
template void CFDv0_solver<double, 2, 3> :: set_boundary_conditions();
template void CFDv0_solver<double, 2, 4> :: set_boundary_conditions();

// Exchange ghost cells
template void CFDv0_solver<float , 2, 3> :: exchange_ghost_cells( MPI_env &mpi_env );
template void CFDv0_solver<float , 2, 4> :: exchange_ghost_cells( MPI_env &mpi_env );
template void CFDv0_solver<double, 2, 3> :: exchange_ghost_cells( MPI_env &mpi_env );
template void CFDv0_solver<double, 2, 4> :: exchange_ghost_cells( MPI_env &mpi_env );

// Print submesh info (rho, vel, p)
template void CFDv0_solver<float , 2, 3> :: gather_info( MPI_env &mpi_env, struct t_domain_info<float ,2> &domain_info );
template void CFDv0_solver<float , 2, 4> :: gather_info( MPI_env &mpi_env, struct t_domain_info<float ,2> &domain_info );
template void CFDv0_solver<double, 2, 3> :: gather_info( MPI_env &mpi_env, struct t_domain_info<double,2> &domain_info );
template void CFDv0_solver<double, 2, 4> :: gather_info( MPI_env &mpi_env, struct t_domain_info<double,2> &domain_info );

//***************************************************************************************************
// gmsh
//***************************************************************************************************
template void CFDv0_solver<float , 2, 3> :: allocate( gmsh_mesh *mesh, std::vector<t_solution_vars<float ,2>*> &cells_ptr );
template void CFDv0_solver<float , 2, 4> :: allocate( gmsh_mesh *mesh, std::vector<t_solution_vars<float ,2>*> &cells_ptr );
template void CFDv0_solver<double, 2, 3> :: allocate( gmsh_mesh *mesh, std::vector<t_solution_vars<double,2>*> &cells_ptr );
template void CFDv0_solver<double, 2, 4> :: allocate( gmsh_mesh *mesh, std::vector<t_solution_vars<double,2>*> &cells_ptr );

template void CFDv0_solver<float , 2, 3> :: allocate_ghost_cells( MPI_env &mpi_env, gmsh_mesh *mesh, std::vector<std::vector<t_solution_vars<float ,2>*>> &ghost_ptr, std::vector<std::vector<t_solution_vars<float ,2>*>> &ghost_bnd );
template void CFDv0_solver<float , 2, 4> :: allocate_ghost_cells( MPI_env &mpi_env, gmsh_mesh *mesh, std::vector<std::vector<t_solution_vars<float ,2>*>> &ghost_ptr, std::vector<std::vector<t_solution_vars<float ,2>*>> &ghost_bnd );
template void CFDv0_solver<double, 2, 3> :: allocate_ghost_cells( MPI_env &mpi_env, gmsh_mesh *mesh, std::vector<std::vector<t_solution_vars<double,2>*>> &ghost_ptr, std::vector<std::vector<t_solution_vars<double,2>*>> &ghost_bnd );
template void CFDv0_solver<double, 2, 4> :: allocate_ghost_cells( MPI_env &mpi_env, gmsh_mesh *mesh, std::vector<std::vector<t_solution_vars<double,2>*>> &ghost_ptr, std::vector<std::vector<t_solution_vars<double,2>*>> &ghost_bnd );

template void CFDv0_solver<float , 2, 3> :: assign_pointers( gmsh_mesh *mesh, std::vector<std::vector<t_solution_vars<float ,2>*>> cells_addr, std::vector<std::vector<t_solution_vars<float ,2>*>> ghost_addr, std::vector<std::vector<t_solution_vars<float ,2>*>> ghost_bnd_addr );
template void CFDv0_solver<float , 2, 4> :: assign_pointers( gmsh_mesh *mesh, std::vector<std::vector<t_solution_vars<float ,2>*>> cells_addr, std::vector<std::vector<t_solution_vars<float ,2>*>> ghost_addr, std::vector<std::vector<t_solution_vars<float ,2>*>> ghost_bnd_addr );
template void CFDv0_solver<double, 2, 3> :: assign_pointers( gmsh_mesh *mesh, std::vector<std::vector<t_solution_vars<double,2>*>> cells_addr, std::vector<std::vector<t_solution_vars<double,2>*>> ghost_addr, std::vector<std::vector<t_solution_vars<double,2>*>> ghost_bnd_addr );
template void CFDv0_solver<double, 2, 4> :: assign_pointers( gmsh_mesh *mesh, std::vector<std::vector<t_solution_vars<double,2>*>> cells_addr, std::vector<std::vector<t_solution_vars<double,2>*>> ghost_addr, std::vector<std::vector<t_solution_vars<double,2>*>> ghost_bnd_addr );

template void CFDv0_solver<float , 2, 3> :: initialize( MPI_env &mpi_env, gmsh_mesh *mesh );
template void CFDv0_solver<float , 2, 4> :: initialize( MPI_env &mpi_env, gmsh_mesh *mesh );
template void CFDv0_solver<double, 2, 3> :: initialize( MPI_env &mpi_env, gmsh_mesh *mesh );
template void CFDv0_solver<double, 2, 4> :: initialize( MPI_env &mpi_env, gmsh_mesh *mesh );

template void CFDv0_solver<float , 2, 3> :: output_cfd( float  time, std::string filename, gmsh_mesh *mesh );
template void CFDv0_solver<float , 2, 4> :: output_cfd( float  time, std::string filename, gmsh_mesh *mesh );
template void CFDv0_solver<double, 2, 3> :: output_cfd( double time, std::string filename, gmsh_mesh *mesh );
template void CFDv0_solver<double, 2, 4> :: output_cfd( double time, std::string filename, gmsh_mesh *mesh );

//template float  CFDv0_solver<float , 2, 3> :: compute_cfl( float  dt );
//template float  CFDv0_solver<float , 2, 4> :: compute_cfl( float  dt );
//template double CFDv0_solver<double, 2, 3> :: compute_cfl( double dt );
//template double CFDv0_solver<double, 2, 4> :: compute_cfl( double dt );

template float  CFDv0_solver<float , 2, 3> :: compute_dt( float  cflMax );
template float  CFDv0_solver<float , 2, 4> :: compute_dt( float  cflMax );
template double CFDv0_solver<double, 2, 3> :: compute_dt( double cflMax );
template double CFDv0_solver<double, 2, 4> :: compute_dt( double cflMax );

template void CFDv0_solver<float , 2, 3> :: halo_comm_wait( MPI_env &mpi_env );
template void CFDv0_solver<float , 2, 4> :: halo_comm_wait( MPI_env &mpi_env );
template void CFDv0_solver<double, 2, 3> :: halo_comm_wait( MPI_env &mpi_env );
template void CFDv0_solver<double, 2, 4> :: halo_comm_wait( MPI_env &mpi_env );

//***************************************************************************************************
//***************************************************************************************************
//									Public functions
//***************************************************************************************************
//***************************************************************************************************

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
void CFDv0_solver<PRECISION, DIM_CNT, FACE_CNT> :: allocate( gmsh_mesh *mesh,
															 std::vector<t_solution_vars<PRECISION,DIM_CNT>*> &cells_addr ){

	int myrank;
	MPI_Comm_rank( MPI_COMM_WORLD, &myrank );

	// Allocate
	cells_cfd .resize( mesh->cells.size() );
	cells_addr.resize( mesh->cells.size() );
	dummy_cell.resize( 1 );

	// Initialize Hpath cells
	for( size_t icell=0; icell < cells_cfd.size(); icell++ ){
		cells_cfd[icell].vars.id    = icell;
		cells_cfd[icell].vars.sm_id = this->sm_id;
		cells_cfd[icell].vars.is_ghost = false;

		for( unsigned idim=0; idim < DIM_CNT+2; idim++ ){
			cells_cfd[icell].vars.q_old[idim] = 0.;
			cells_cfd[icell].vars.q_new[idim] = 0.;
			cells_cfd[icell].vars.delta_q[idim] = 0.;
		}
		cells_cfd[icell].vars.q_old[0] = 100000 + icell;
		cells_cfd[icell].vars.q_old[1] = 200000 + icell;
		cells_cfd[icell].vars.q_old[2] = 300000 + icell;
		cells_cfd[icell].vars.q_old[3] = 400000 + icell;
	}

	// Dummy cell, used for NEIGHBOR_NONE:
	// 		- instead of having a special case for out-of-bounds neighbors, have a "dummy cell" with
	// 			all variables initialized to zero
	// 		- saves a condition statement or extra loop depending on solution chosen
	dummy_var.q_old[0] = 123456.;

	// Record local cells' memory addresses
	for( size_t icell=0; icell < cells_cfd.size(); icell++ ){
		cells_addr[icell] = &( cells_cfd[icell].vars );
	}

	// Print info

	if( myrank == 0 ){
		cout << "    - num_cells        =" << setw(10) << sm_id
									   	   << setw(10) << cells_cfd.size()
									   	   << endl;
		cout << "    - size of CFD cell =" << setw(10) << sizeof( cells_cfd[0] ) << endl;
		cout << "    - size of CFD vars =" << setw(10) << sizeof( cells_cfd[0].vars ) << endl;
	}
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Initialize the CFD solver: parameters, mpi type, initial/boundary conditions, etc...
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
void CFDv0_solver<PRECISION, DIM_CNT, FACE_CNT> :: initialize( MPI_env &mpi_env, gmsh_mesh *mesh ){

	// Initialize all parameters, mostly geometrical variables for the CFD solver
	init_params( mesh );

	// Initilize MPI environment: only for boundary submesh
	if( mesh->id == INDEX_BND_SUBMESH ){
		initialize_mpi_env( mpi_env, mesh );
	}

	// Setup initial conditions
	initial_condition init_case = CYLINDER_TUNNEL;
	set_init_conditions( init_case );

	if( mesh->id == INDEX_BND_SUBMESH ){
		// Setup boundary conditions, not always necessary
		set_boundary_conditions();
	}
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
void CFDv0_solver<PRECISION, DIM_CNT, FACE_CNT> :: initialize_mpi_env( MPI_env &mpi_env, gmsh_mesh *mesh ){

	// Get MPI neighbors
	std::vector<int > mpi_neighbors;
	std::vector<bool> registered_rank( mpi_env.size(), false );

	for( size_t icell=0; icell < cells_cfd.size(); icell++ ){
		for( unsigned iface=0; iface < FACE_CNT; iface++ ){

			// Isolate MPI/periodic neighbors
			auto *neigh = &( mesh->cell2neigh[icell][iface] );
			switch( neigh->type ){
				case CELL_REGULAR:
				case CELL_SUBMESH:
				case CELL_NONE:
					continue;
				case CELL_MPI:
				case CELL_PERIODIC:
				case CELL_PERIODIC_MPI:
					break;
			}

			// Filter out MPI neighbors we already registered
			if( registered_rank[neigh->sm] ){
				continue;
			}

			// Register neighbor
			mpi_neighbors.push_back( neigh->sm );
			registered_rank[neigh->sm] = true;
		}
	}

	mpi_env.set_mpi_neighbors( mpi_neighbors );

	// Should only happen for a serial case without periodic boundary conditions
	if( mpi_neighbors.size() == 0 )
		return;

	init_mpi_types( mpi_env );
	mpi_env.set_halo_comm_reqs( 2*mpi_env.mpi_neighbors.size() );
	mpi_env.set_halo_comm_type( MPI_TWOSIDED_NONB, cells_cfd, ghost_mpi );
}

//***************************************************************************************************
template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
void CFDv0_solver<PRECISION, DIM_CNT, FACE_CNT> :: halo_comm_wait( MPI_env &mpi_env ){

	// Wait for ghost cells
	mpi_env.halo_comm_wait();

	// Fix ghost cells
	for( size_t ineigh=0; ineigh < ghost_mpi.size(); ineigh++ ){
		for( size_t icell=0; icell < ghost_mpi[ineigh].size(); icell++ ){
			ghost_mpi[ineigh][icell].vars.sm_id    = 10000;
		}
	}
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// New method
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
void CFDv0_solver<PRECISION, DIM_CNT, FACE_CNT> ::
				allocate_ghost_cells( MPI_env &mpi_env,
									  gmsh_mesh *mesh,
									  std::vector<std::vector<t_solution_vars<PRECISION,DIM_CNT>*>> &ghost_ptr,
									  std::vector<std::vector<t_solution_vars<PRECISION,DIM_CNT>*>> &ghost_ptr_bnd ){

	if( mesh->id != INDEX_BND_SUBMESH )
		return;

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Initialize boundary ghost cells
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	boundaries.resize( 4 );
	for( size_t icell=0; icell < mesh->cells.size(); icell++ ){
		for( unsigned iface=0; iface < FACE_CNT; iface++ ){

			t_gmsh_neighbor *neigh = &( mesh->cell2neigh[icell][iface] );

			if( neigh->type != CELL_NONE )
				continue;

			if( neigh->bound > 0 && neigh->bound < 5 ){
				int ibound = neigh->bound - 1;

				boundaries[ibound].push_back( { (int) icell, (int) iface } );
				neigh->id = boundaries[ibound].size() -1 ;
			}
		}
	}
	ghost_bnd    .resize( 4 );
	ghost_ptr_bnd.resize( 4 );

	for( size_t ibound=0; ibound < ghost_ptr_bnd.size(); ibound++ ){
		ghost_bnd    [ibound].resize( boundaries[ibound].size() );
		ghost_ptr_bnd[ibound].resize( boundaries[ibound].size() );
		for( size_t icell=0; icell < ghost_bnd[ibound].size(); icell++ ){
			for( unsigned idim=0; idim < DIM_CNT+2; idim++ ){
				ghost_bnd[ibound][icell].vars.q_old[idim] = 0.;
				ghost_bnd[ibound][icell].vars.q_new[idim] = 0.;
			}
		}
	}

	for( size_t ibound=0; ibound < ghost_bnd.size(); ibound++ ){
		for( size_t icell=0; icell < ghost_bnd[ibound].size(); icell++ ){
			ghost_ptr_bnd[ibound][icell] = &( ghost_bnd[ibound][icell].vars );

			ghost_bnd[ibound][icell].vars.sm_id = 10000;
			if (ibound == 0)
				ghost_bnd[ibound][icell].vars.id = mesh->cells.size() + icell;
			else
				ghost_bnd[ibound][icell].vars.id = mesh->cells.size() + ibound*ghost_bnd[ibound-1].size() + icell;
		}
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Initialize MPI ghost cells
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	std::vector<int > mpi_neighbors;
	std::vector<bool> registered_rank( mpi_env.size(), false );

	for( size_t icell=0; icell < cells_cfd.size(); icell++ ){
		for( unsigned iface=0; iface < FACE_CNT; iface++ ){

			// Isolate MPI/periodic neighbors
			auto *neigh = &( mesh->cell2neigh[icell][iface] );
			switch( neigh->type ){
				case CELL_REGULAR:
				case CELL_SUBMESH:
				case CELL_NONE:
					continue;
				case CELL_MPI:
				case CELL_PERIODIC:
				case CELL_PERIODIC_MPI:
					break;
			}

			// Filter out MPI neighbors we already registered
			if( registered_rank[neigh->sm] ){
				continue;
			}

			// Register neighbor
			mpi_neighbors.push_back( neigh->sm );
			registered_rank[neigh->sm] = true;
		}
	}

	// Should only happen for a serial case without periodic boundary conditions
	// 		==> true for now but will never happen once we implement support for ALL BCs (wall, etc...)
	if( mpi_neighbors.size() == 0 )
		return;

	// Initialize mapping from rank ID to local MPI neighbor
	rank2local.resize( mpi_env.size(), -1 );
	for( size_t ineigh=0; ineigh < mpi_neighbors.size(); ineigh++ ){
		rank2local[mpi_neighbors[ineigh]] = ineigh;
	}

	// Count number of ghost cells from each neighbor
	int million = 1000;
	std::vector<MPI_Request> comm_reqs( 2*mpi_neighbors.size() );
	std::vector<int> ghost_size( mpi_neighbors.size(), -1 );
	for( size_t ineigh=0; ineigh < mpi_neighbors.size(); ineigh++ ){

		int neigh_rank = mpi_neighbors[ineigh];

		// Hpath cellls
		int send_tag = million*mpi_env.rank() + mpi_neighbors[ineigh];
		int recv_tag = mpi_env.rank() + million*mpi_neighbors[ineigh];

		int local_size = (int) cells_cfd.size();

		MPI_Isend( &local_size, 1, MPI_INT, neigh_rank, send_tag, MPI_COMM_WORLD, &comm_reqs[ineigh] );
		MPI_Irecv( &ghost_size[ineigh], 1, MPI_INT, neigh_rank, recv_tag, MPI_COMM_WORLD, &comm_reqs[ineigh+mpi_neighbors.size()] );
	}
	MPI_Waitall( comm_reqs.size(), &comm_reqs[0], MPI_STATUS_IGNORE );

	// Allocate ghost cells
	ghost_mpi.resize( mpi_neighbors.size() );
	ghost_ptr  .resize( mpi_neighbors.size() );
	for( size_t ineigh=0; ineigh < ghost_mpi.size(); ineigh++ ){
			ghost_mpi[ineigh].resize( ghost_size[ineigh] );
			ghost_ptr  [ineigh].resize( ghost_size[ineigh] );
	}

	// Store ghost cells' memory addresses
	for( size_t ineigh=0; ineigh < ghost_mpi.size(); ineigh++ ){
		for( size_t icell=0; icell < ghost_mpi[ineigh].size(); icell++ ){
			ghost_ptr[ineigh][icell] = &( ghost_mpi[ineigh][icell].vars );

			ghost_mpi[ineigh][icell].vars.sm_id = mpi_env.rank();
			ghost_mpi[ineigh][icell].vars.id    = icell;

			ghost_mpi[ineigh][icell].vars.q_old[0] = 0.0;
			ghost_mpi[ineigh][icell].vars.q_old[1] = 0.0;
			ghost_mpi[ineigh][icell].vars.q_old[2] = 0.0;
			ghost_mpi[ineigh][icell].vars.q_old[3] = 0.0;
		}
	}
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// New method
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
void CFDv0_solver<PRECISION, DIM_CNT, FACE_CNT> :: assign_pointers(
				gmsh_mesh                                                      *mesh,
				std::vector<std::vector<t_solution_vars<PRECISION, DIM_CNT>*>>  cells_ptr,
				std::vector<std::vector<t_solution_vars<PRECISION, DIM_CNT>*>>  ghost_ptr,
				std::vector<std::vector<t_solution_vars<PRECISION, DIM_CNT>*>>  ghost_ptr_bnd ){


	// Assign pointers
	for( size_t icell=0; icell < cells_cfd.size(); icell++ ){
		for( unsigned iface=0; iface < FACE_CNT; iface++ ){
			t_gmsh_neighbor *neigh = &( mesh->cell2neigh[icell][iface] );

			switch( neigh->type ){
				// Same submesh neighbor
				case CELL_PERIODIC:
				case CELL_REGULAR:
					cells_cfd[icell].neighs[iface] = cells_ptr[this->sm_id][neigh->id];
					break;
				// Different submesh, same MPI rank
				case CELL_SUBMESH:
					cells_cfd[icell].neighs[iface] = cells_ptr[neigh->sm][neigh->id];
					break;
				// All of the following neighbor types require ghost cells
				case CELL_NONE:
					cells_cfd[icell].neighs[iface] = ghost_ptr_bnd[neigh->bound-1][neigh->id];
					break;
				case CELL_MPI:
				case CELL_PERIODIC_MPI:
					cells_cfd[icell].neighs[iface] = ghost_ptr[rank2local[neigh->sm]][neigh->id];
					break;
				default:
					break;
			}
		}
	}
}

//***************************************************************************************************
//***************************************************************************************************
//									PRIVATE FUNCTIONS
//***************************************************************************************************
//***************************************************************************************************

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Initialize the CFD solver
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
void CFDv0_solver<PRECISION, DIM_CNT, FACE_CNT> :: init_params( gmsh_mesh *mesh ){

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Hpath cells
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	cells_cfd.set_start( 0 );
	cells_cfd.set_end  ( cells_cfd.size() );

	// Cell's center
	xc.resize( cells_cfd.size(), std::vector<PRECISION>( DIM_CNT ));
	for( size_t icell=0; icell < cells_cfd.size(); icell++ ){
		for( unsigned idim = 0; idim < DIM_CNT; idim++ ){
			xc[icell][idim] = mesh->cells[icell].xc[idim] ;
		}
	}

	// Compute face surface vector
	for( size_t icell=0; icell < cells_cfd.size(); icell++ ){
		for( unsigned iface=0; iface < FACE_CNT; iface++ ){
			int gface = mesh->cell2face[icell][iface];
			for( unsigned idim = 0; idim < DIM_CNT; idim++ ){
				cells_cfd[icell].S[iface][idim] = mesh->v_sign[icell][iface] * mesh->v_norm[gface][idim] ;
			}
		}
	}

	// Vector from face center to cell's center
	std::vector<std::vector<std::vector<PRECISION>>> dP( cells_cfd.size(), std::vector<std::vector<PRECISION>> ( FACE_CNT ) );
	for( size_t icell=0; icell < cells_cfd.size(); icell++ ){
		for( unsigned iface=0; iface < FACE_CNT; iface++ ){
			dP[icell][iface].resize( DIM_CNT, 0.0 );
		}
	}
	for( size_t icell=0; icell < cells_cfd.size(); icell++ ){
		for( unsigned iface=0; iface < FACE_CNT; iface++ ){
			int gface = mesh->cell2face[icell][iface];
			for( unsigned idim = 0; idim < DIM_CNT; idim++ ){
				dP[icell][iface][idim] = mesh->cells[icell].xc[idim] - mesh->faces[gface].xf[idim] ;
			}
		}
	}

	// Vector from cell center to neighbor's center
	for( size_t icell=0; icell < cells_cfd.size(); icell++ ){
		for( unsigned iface=0; iface < FACE_CNT; iface++ ){
			int gface = mesh->cell2face[icell][iface];
			auto *neigh = &( mesh->cell2neigh[icell][iface] );
			if( neigh->type == CELL_NONE ){
				PRECISION S_mag_sqrt = 0.0;
				PRECISION dP_dot_S = 0.0;
				for( unsigned idim = 0; idim < DIM_CNT; idim++ ){
					S_mag_sqrt += cells_cfd[icell].S[iface][idim] * cells_cfd[icell].S[iface][idim];
				}
				for( unsigned idim = 0; idim < DIM_CNT; idim++ ){
					dP_dot_S += dP[icell][iface][idim] * cells_cfd[icell].S[iface][idim];
				}
				for( unsigned idim = 0; idim < DIM_CNT; idim++ ){
					cells_cfd[icell].d[iface][idim] = 2.0 * abs( dP_dot_S ) * cells_cfd[icell].S[iface][idim] / S_mag_sqrt;
				}
			}else{
				for( unsigned idim = 0; idim < DIM_CNT; idim++ ){
					cells_cfd[icell].d[iface][idim] = mesh->v_sign[icell][iface] * mesh->v_neigh[gface][idim] ;
				}
			}
		}
	}

	// Vector from face center to neighbor's center
	std::vector<std::vector<std::vector<PRECISION>>> dN( cells_cfd.size(), std::vector<std::vector<PRECISION>> ( FACE_CNT ) );
	for( size_t icell=0; icell < cells_cfd.size(); icell++ ){
		for( unsigned iface=0; iface < FACE_CNT; iface++ ){
			dN[icell][iface].resize( DIM_CNT, 0.0 );
		}
	}
	for( size_t icell=0; icell < cells_cfd.size(); icell++ ){
		for( unsigned iface=0; iface < FACE_CNT; iface++ ){
			for( unsigned idim = 0; idim < DIM_CNT; idim++ ){
				dN[icell][iface][idim] = dP[icell][iface][idim] + cells_cfd[icell].d[iface][idim] ;
			}
		}
	}

	// beta_pos and beta_neg calculations
	for( size_t icell=0; icell < cells_cfd.size(); icell++ ){
		for( unsigned iface=0; iface < FACE_CNT; iface++ ){
			cells_cfd[icell].weight_linear[iface] = 0.0;
			PRECISION divider = 0.0;
			for( unsigned idim = 0; idim < DIM_CNT; idim++ ){
				cells_cfd[icell].weight_linear[iface] += cells_cfd[icell].S[iface][idim] * dN[icell][iface][idim] ;
				divider += ( - cells_cfd[icell].S[iface][idim] * dP[icell][iface][idim] + cells_cfd[icell].S[iface][idim] * dN[icell][iface][idim] );
			}
			cells_cfd[icell].weight_linear[iface] /= divider ;
		}
	}

	// Cell's volume
	for( size_t icell=0; icell < cells_cfd.size(); icell++ ){
		cells_cfd[icell].vars.vol_inv = ONE / mesh->cells_vol[icell];
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Constants
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	M  = 1.0;
	AoA = 0.0;
	Cp = 1004.7;
	Cv = Cp - Rgas;
	gamma = Cp / Cv;
	gamma_m_one = gamma - ONE;
	gamma_m_one_inv = ONE / gamma_m_one;
	Pr = 0.71;
	Pr_inv = ONE / Pr ;
	mu0 = 1.8e-5;
	T0 = 273.15;
	S = 110.4;
	L = 1.0;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Initial conditions: currently 2D isentropic vortex
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
void CFDv0_solver<PRECISION, DIM_CNT, FACE_CNT> :: set_init_conditions(initial_condition init_case){

	for( size_t icell=0; icell < cells_cfd.size(); icell++ ){
		for( unsigned idim = 0; idim < DIM_CNT+2; idim++ ){
			cells_cfd[icell].vars.q_old[idim] = ZERO;
			cells_cfd[icell].vars.q_new[idim] = ZERO;
			cells_cfd[icell].vars.delta_q[idim] = ZERO;
		}
	}

	vector<PRECISION> rho, u[DIM_CNT], E;

	rho.resize( cells_cfd.size() );
	for( unsigned idim = 0; idim < DIM_CNT; idim++ ){
		u[idim].resize( cells_cfd.size() );
	}
	E.resize( cells_cfd.size() );

	switch (init_case) {
		case ISENTROPIC_VORTEX: {
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Setup case/initial conditions: 2D ISENTROPIC VORTEX
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			PRECISION alpha = 0.0;
			PRECISION beta  = 5.0;
			PRECISION mesh_xc[2] = {0.0, 0.0};

			for( size_t icell=0; icell < cells_cfd.size(); icell++ ){
				PRECISION r = sqrt( pow( xc[icell][0] - mesh_xc[0], 2)
								  + pow( xc[icell][1] - mesh_xc[1], 2) );

				rho[icell] = pow(ONE - beta * beta * (gamma - ONE) * exp(ONE - r * r) / (8.0 * gamma * M_PI * M_PI), ONE / (gamma - ONE));

				u[0][icell] = M * cos(alpha * M_PI / 180.0) - beta * ( xc[icell][1] - mesh_xc[1]) * exp(0.5 * (ONE - r * r)) / (2.0 * M_PI);
				u[1][icell] = M * sin(alpha * M_PI / 180.0) + beta * ( xc[icell][0] - mesh_xc[0]) * exp(0.5 * (ONE - r * r)) / (2.0 * M_PI);

				PRECISION u_mag = 0.0;
				for( unsigned idim = 0; idim < DIM_CNT; idim++ ){
					u_mag += u[idim][icell] * u[idim][icell];
				}
				u_mag = sqrt(u_mag);

				PRECISION p = pow(rho[icell], gamma);

				E[icell] = p / (rho[icell] * (gamma - ONE)) + 0.5 * u_mag * u_mag;
			}
			break;
		}
		case DOUBLE_SHEAR_LAYER_1: {
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Setup case/initial conditions: 2D DOUBLE SHEAR LAYER
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			for( size_t icell=0; icell < cells_cfd.size(); icell++ ){
				if (xc[icell][1] >= 0.25 || xc[icell][1] < -0.25) {
					rho[icell]  =  1.0;
					u[0][icell] =  1.0;
				} else {
					rho[icell]  =  2.0;
					u[0][icell] = -1.0;
				}
				u[1][icell] =  0.01 * sin(2.0 * M_PI * 2.0 * xc[icell][0]);

				PRECISION u_mag = 0.0;
				for( unsigned idim = 0; idim < DIM_CNT; idim++ ){
					u_mag += u[idim][icell] * u[idim][icell];
				}
				u_mag = sqrt(u_mag);

				PRECISION p = 2.5;

				E[icell] = p / (rho[icell] * (gamma - ONE)) + 0.5 * u_mag * u_mag;
			}
			break;
		}
		case DOUBLE_SHEAR_LAYER_2: {
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Setup case/initial conditions: 2D DOUBLE SHEAR LAYER (Another version)
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			for( size_t icell=0; icell < cells_cfd.size(); icell++ ){
				rho[icell] = 1.0;

				if (xc[icell][1] <= 0.5){
					u[0][icell] = tanh(80.0 * (xc[icell][1] - 0.25));
				} else {
					u[0][icell] = tanh(80.0 * (0.75 - xc[icell][1]));
				}
				u[1][icell] = 0.05 * sin(2.0 * M_PI * (xc[icell][0] + 0.25));

				PRECISION u_mag = 0.0;
				for( unsigned idim = 0; idim < DIM_CNT; idim++ ){
					u_mag += u[idim][icell] * u[idim][icell];
				}
				u_mag = sqrt(u_mag);

				PRECISION p = pow(rho[icell], 1.4);

				E[icell] = p / (rho[icell] * (gamma - ONE)) + 0.5 * u_mag * u_mag;
			}
			break;
		}
		case TAYLOR_GREEN_VORTEX: {
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Setup case/initial conditions: TAYLOR GREEN VORTEX
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			PRECISION rho0 = 1.0;
			PRECISION u0 = 1.0;
			PRECISION p0 = 1.0;

			for( size_t icell=0; icell < cells_cfd.size(); icell++ ){
				rho[icell] = rho0;

				u[0][icell] =   u0 * cos(xc[icell][0] / L) * sin(xc[icell][1] / L);
				u[1][icell] = - u0 * sin(xc[icell][0] / L) * cos(xc[icell][1] / L);

				PRECISION u_mag = 0.0;
				for( unsigned idim = 0; idim < DIM_CNT; idim++ ){
					u_mag += u[idim][icell] * u[idim][icell];
				}
				u_mag = sqrt(u_mag);

				PRECISION p = p0 - rho0 / 4.0 * (cos(2.0 * xc[icell][0] / L) + cos(2.0 * xc[icell][1] / L));

				E[icell] = p / (rho[icell] * (gamma - ONE)) + 0.5 * u_mag * u_mag;
			}
			break;
		}
		case CYLINDER_TUNNEL: {
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Setup case/initial conditions: Cylinder in Tunnel
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			PRECISION rho0 = 1.0;
			PRECISION p0 = 1.0;

			for( size_t icell=0; icell < cells_cfd.size(); icell++ ){
				rho[icell] = rho0;

				u[0][icell] = 6.0 * (0.41 / 2.0 - xc[icell][1]) * (0.41 / 2.0 + xc[icell][1] ) / ( 0.41 * 0.41 );
				u[1][icell] = 0.0;

				PRECISION u_mag = 0.0;
				for( unsigned idim = 0; idim < DIM_CNT; idim++ ){
					u_mag += u[idim][icell] * u[idim][icell];
				}
				u_mag = sqrt(u_mag);

				PRECISION p = p0;

				E[icell] = p / (rho[icell] * (gamma - ONE)) + 0.5 * u_mag * u_mag;
			}
			break;
		}
	}

	for( size_t icell=0; icell < cells_cfd.size(); icell++ ){
		cells_cfd[icell].vars.q_new[0] = rho[icell];
		for( unsigned idim = 0; idim < DIM_CNT; idim++ ){
			cells_cfd[icell].vars.q_new[idim+1] = rho[icell] * u[idim][icell];
		}
		cells_cfd[icell].vars.q_new[DIM_CNT+1] = rho[icell] * E[icell];
	}
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Set boundary conditions
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
void CFDv0_solver<PRECISION, DIM_CNT, FACE_CNT> :: set_boundary_conditions(){

	// WALL
	for ( size_t icell = 0; icell < ghost_bnd[0].size(); icell++ ) {
		ghost_bnd[0][icell].vars.q_new[0] = cells_cfd[boundaries[0][icell][0]].vars.q_new[0];
		for ( size_t idim = 0; idim < DIM_CNT; idim++ ) {
			ghost_bnd[0][icell].vars.q_new[idim+1] = - cells_cfd[boundaries[0][icell][0]].vars.q_new[idim+1];
		}
		ghost_bnd[0][icell].vars.q_new[DIM_CNT+1] = cells_cfd[boundaries[0][icell][0]].vars.q_new[DIM_CNT+1];

		for ( size_t idim = 0; idim < DIM_CNT; idim++ ) {
			ghost_bnd[0][icell].vars.rho_grad[idim] = cells_cfd[boundaries[0][icell][0]].vars.rho_grad[idim];
			for ( size_t jdim = 0; jdim < DIM_CNT; jdim++ ) {
				ghost_bnd[0][icell].vars.rhoU_grad[idim][jdim] = - cells_cfd[boundaries[0][icell][0]].vars.rhoU_grad[idim][jdim];
				ghost_bnd[0][icell].vars.dudx[idim][jdim]      = - cells_cfd[boundaries[0][icell][0]].vars.dudx[idim][jdim];
				ghost_bnd[0][icell].vars.tauMC[idim][jdim]     = - cells_cfd[boundaries[0][icell][0]].vars.tauMC[idim][jdim];
			}
			ghost_bnd[0][icell].vars.rhoE_grad[idim] = cells_cfd[boundaries[0][icell][0]].vars.rhoE_grad[idim];
			ghost_bnd[0][icell].vars.Rpsi_grad[idim] = cells_cfd[boundaries[0][icell][0]].vars.Rpsi_grad[idim];
			ghost_bnd[0][icell].vars.c_grad[idim]    = cells_cfd[boundaries[0][icell][0]].vars.c_grad[idim];

			ghost_bnd[0][icell].vars.dTdx[idim]      = cells_cfd[boundaries[0][icell][0]].vars.dTdx[idim];
		}
	}

	// INLET
	for ( size_t icell = 0; icell < ghost_bnd[1].size(); icell++ ) {

		int gcell = boundaries[1][icell][0];

		PRECISION rho = 1.0;
		PRECISION Uvec[DIM_CNT];
		Uvec[0] = 6.0 * (0.41 / 2.0 - xc[gcell][1]) * (0.41 / 2.0 + xc[gcell][1] ) / ( 0.41 * 0.41 );
		Uvec[1] = 0.0;
		PRECISION Umag_sqrt = 0.0;
		for ( size_t idim = 0; idim < DIM_CNT; idim++ ) {
			Umag_sqrt += Uvec[idim] * Uvec[idim];
		}
		PRECISION p = 1.0;
		PRECISION E = p / (rho * (gamma - ONE)) + 0.5 * Umag_sqrt;

		ghost_bnd[1][icell].vars.q_new[0] = rho;
		for ( size_t idim = 0; idim < DIM_CNT; idim++ ) {
			ghost_bnd[1][icell].vars.q_new[idim+1] = rho * Uvec[idim];
		}
		ghost_bnd[1][icell].vars.q_new[DIM_CNT+1] = rho * E;

		for ( size_t idim = 0; idim < DIM_CNT; idim++ ) {
			ghost_bnd[1][icell].vars.rho_grad[idim] = 0.0;
			for ( size_t jdim = 0; jdim < DIM_CNT; jdim++ ) {
				ghost_bnd[1][icell].vars.rhoU_grad[idim][jdim] = 0.0;
				ghost_bnd[1][icell].vars.dudx[idim][jdim]      = 0.0;
				ghost_bnd[1][icell].vars.tauMC[idim][jdim]     = 0.0;
			}
			ghost_bnd[1][icell].vars.rhoE_grad[idim] = 0.0;
			ghost_bnd[1][icell].vars.Rpsi_grad[idim] = 0.0;
			ghost_bnd[1][icell].vars.c_grad[idim]    = 0.0;

			ghost_bnd[1][icell].vars.dTdx[idim]      = 0.0;
		}
	}

	// OUTLET
	for ( size_t icell = 0; icell < ghost_bnd[2].size(); icell++ ) {
		for ( size_t idim = 0; idim < DIM_CNT+2; idim++ ) {
			ghost_bnd[2][icell].vars.q_new[idim] = cells_cfd[boundaries[2][icell][0]].vars.q_new[idim];
		}

		for ( size_t idim = 0; idim < DIM_CNT; idim++ ) {
			ghost_bnd[2][icell].vars.rho_grad[idim] = 0.0;
			for ( size_t jdim = 0; jdim < DIM_CNT; jdim++ ) {
				ghost_bnd[2][icell].vars.rhoU_grad[idim][jdim] = 0.0;
				ghost_bnd[2][icell].vars.dudx[idim][jdim]      = 0.0;
				ghost_bnd[2][icell].vars.tauMC[idim][jdim]     = 0.0;
			}
			ghost_bnd[2][icell].vars.rhoE_grad[idim] = 0.0;
			ghost_bnd[2][icell].vars.Rpsi_grad[idim] = 0.0;
			ghost_bnd[2][icell].vars.c_grad[idim]    = 0.0;

			ghost_bnd[2][icell].vars.dTdx[idim]      = 0.0;
		}
	}

	// FAR-FIELD
	for ( size_t icell = 0; icell < ghost_bnd[3].size(); icell++ ) {

		PRECISION Udirection[DIM_CNT] = { (PRECISION) cos(AoA * M_PI / 180.0), (PRECISION) sin(AoA * M_PI / 180.0) };
		PRECISION T = 0.0;
		for ( size_t idim = 0; idim < DIM_CNT; idim++ ) {
			T += Udirection[idim] * cells_cfd[boundaries[3][icell][0]].S[idim][boundaries[3][icell][1]];
		}

		if ( T >= 0.0 ) {		// INLET

			PRECISION rho = 1.0;
			PRECISION Uvec[DIM_CNT] = {1.0, 0.0};
			PRECISION Umag_sqrt = 0.0;
			for ( size_t idim = 0; idim < DIM_CNT; idim++ ) {
				Umag_sqrt += Uvec[idim] * Uvec[idim];
			}
			PRECISION p = 1.0;
			PRECISION E = p / (rho * (gamma - ONE)) + 0.5 * Umag_sqrt;

			ghost_bnd[3][icell].vars.q_new[0] = rho;
			for ( size_t idim = 0; idim < DIM_CNT; idim++ ) {
				ghost_bnd[3][icell].vars.q_new[idim+1] = rho * Uvec[idim];
			}
			ghost_bnd[3][icell].vars.q_new[DIM_CNT+1] = rho * E;

		} else {				 // OUTLET

			for ( size_t idim = 0; idim < DIM_CNT+2; idim++ ) {
				ghost_bnd[3][icell].vars.q_new[idim] = cells_cfd[boundaries[3][icell][0]].vars.q_new[idim];
			}

		}

		for ( size_t idim = 0; idim < DIM_CNT; idim++ ) {
			ghost_bnd[3][icell].vars.rho_grad[idim] = 0.0;
			for ( size_t jdim = 0; jdim < DIM_CNT; jdim++ ) {
				ghost_bnd[3][icell].vars.rhoU_grad[idim][jdim] = 0.0;
				ghost_bnd[3][icell].vars.dudx[idim][jdim]      = 0.0;
				ghost_bnd[3][icell].vars.tauMC[idim][jdim]     = 0.0;
			}
			ghost_bnd[3][icell].vars.rhoE_grad[idim] = 0.0;
			ghost_bnd[3][icell].vars.Rpsi_grad[idim] = 0.0;
			ghost_bnd[3][icell].vars.c_grad[idim]    = 0.0;

			ghost_bnd[3][icell].vars.dTdx[idim]      = 0.0;
		}
	}
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Initialize MPI types
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
void CFDv0_solver<PRECISION, DIM_CNT, FACE_CNT> :: init_mpi_types( MPI_env &mpi_env ){

	if( mpi_env.size() == 1 ){
		return;
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// MPI Datatype for t_solution_vars
	// 		- 15 fields total
	// 		- without sm_id, id, and D2U
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	const int solvars_nfields = 15;
	MPI_Aint  solvars_disps[solvars_nfields];
	int       solvars_blocklens[] = { DIM_CNT+2, DIM_CNT+2, DIM_CNT+2, DIM_CNT, DIM_CNT*DIM_CNT,
									  DIM_CNT, DIM_CNT, DIM_CNT, DIM_CNT*DIM_CNT, DIM_CNT,
									  DIM_CNT*DIM_CNT, 1, 1, 1, 1 };

	MPI_Datatype mpi_precision;
	if( std::is_same< PRECISION, float >::value ){
		mpi_precision = MPI_FLOAT;
	}else if( std::is_same< PRECISION, double >::value ){
		mpi_precision = MPI_DOUBLE;
	}else{
		MPI_Abort( MPI_COMM_WORLD, 254 );
	}

	MPI_Datatype solvars_types[] = { mpi_precision, mpi_precision, mpi_precision, mpi_precision,
									 mpi_precision, mpi_precision, mpi_precision, mpi_precision,
									 mpi_precision, mpi_precision, mpi_precision,
									 mpi_precision, MPI_UNSIGNED, MPI_UNSIGNED_SHORT, MPI_LOGICAL };

	t_solution_vars<PRECISION, DIM_CNT> tmp_solvars;

	solvars_disps[ 0] = (char*) &tmp_solvars.q_old     - (char*) &tmp_solvars;
	solvars_disps[ 1] = (char*) &tmp_solvars.q_new     - (char*) &tmp_solvars;
	solvars_disps[ 2] = (char*) &tmp_solvars.delta_q   - (char*) &tmp_solvars;
	solvars_disps[ 3] = (char*) &tmp_solvars.rho_grad  - (char*) &tmp_solvars;
	solvars_disps[ 4] = (char*) &tmp_solvars.rhoU_grad - (char*) &tmp_solvars;
	solvars_disps[ 5] = (char*) &tmp_solvars.rhoE_grad - (char*) &tmp_solvars;
	solvars_disps[ 6] = (char*) &tmp_solvars.Rpsi_grad - (char*) &tmp_solvars;
	solvars_disps[ 7] = (char*) &tmp_solvars.c_grad    - (char*) &tmp_solvars;
	solvars_disps[ 8] = (char*) &tmp_solvars.dudx      - (char*) &tmp_solvars;
	solvars_disps[ 9] = (char*) &tmp_solvars.dTdx      - (char*) &tmp_solvars;
	solvars_disps[10] = (char*) &tmp_solvars.tauMC     - (char*) &tmp_solvars;
	solvars_disps[11] = (char*) &tmp_solvars.vol_inv   - (char*) &tmp_solvars;
	solvars_disps[12] = (char*) &tmp_solvars.id        - (char*) &tmp_solvars;
	solvars_disps[13] = (char*) &tmp_solvars.sm_id     - (char*) &tmp_solvars;
	solvars_disps[14] = (char*) &tmp_solvars.is_ghost  - (char*) &tmp_solvars;

	MPI_Datatype type_solvars;
	MPI_Type_create_struct ( solvars_nfields, solvars_blocklens, solvars_disps, solvars_types, &type_solvars );
	MPI_Type_commit( &type_solvars );

	// MPI Datatype for the CFDv0 cell
	const int cfdcell_nfields = 5;
	MPI_Aint  cfdcell_disps[cfdcell_nfields];
	int       cfdcell_blocklens[] = { 1,				// solution vars struct
									  FACE_CNT,			// array of pointers to neighbors' sol vars
									  FACE_CNT*DIM_CNT,	// S
									  FACE_CNT*DIM_CNT,	// d
									  FACE_CNT,			// weight_linear
									 };

	MPI_Datatype cfdcell_types[] = { type_solvars,			// t_solution_vars
									 mpi_precision,			// neighbors
									 mpi_precision,			// S
									 mpi_precision,			// d
									 mpi_precision,			// weight_linear
									};

	// Compute displacements in the CFDv0 cell
	CFDv0_cell<PRECISION, DIM_CNT, FACE_CNT> tmp_cfdcell;

	cfdcell_disps[ 0] = (char*) &tmp_cfdcell.vars          - (char*) &tmp_cfdcell;
	cfdcell_disps[ 1] = (char*) &tmp_cfdcell.neighs        - (char*) &tmp_cfdcell;
	cfdcell_disps[ 2] = (char*) &tmp_cfdcell.S             - (char*) &tmp_cfdcell;
	cfdcell_disps[ 3] = (char*) &tmp_cfdcell.d             - (char*) &tmp_cfdcell;
	cfdcell_disps[ 4] = (char*) &tmp_cfdcell.weight_linear - (char*) &tmp_cfdcell;

	MPI_Aint total_memory_disp = (char*) &cells_cfd[1] - (char*) &cells_cfd[0];

	MPI_Datatype tmp_type, type_cfdcell;
	MPI_Type_create_struct ( cfdcell_nfields, cfdcell_blocklens, cfdcell_disps, cfdcell_types, &tmp_type );
	MPI_Type_create_resized( tmp_type, 0, total_memory_disp, &type_cfdcell );
	MPI_Type_commit( &type_cfdcell );

	// Commit type to MPI env
	mpi_env.set_type_cfdcell( type_cfdcell );
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Restart simulation from t_start if t_start > 0.0
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
void CFDv0_solver<PRECISION, DIM_CNT, FACE_CNT> :: restart( PRECISION t_start, PRECISION PRINT_FILE_FREQUENCY, int submesh ){

	int myrank;
	MPI_Comm_rank( MPI_COMM_WORLD, &myrank );

	stringstream filename;
	if( myrank < 10 ){
		filename << "output/cfdv0_p00" << myrank;
	}else if( myrank < 100 ){
		filename << "output/cfdv0_p0"  << myrank;
	}else{
		filename << "output/cfdv0_p"   << myrank;
	}

	if( submesh < 10 ){
		filename << "_s0" << submesh << "_t" << t_start / PRINT_FILE_FREQUENCY << ".plt";
	}else if( submesh < 100 ){
		filename << "_s"  << submesh << "_t" << t_start / PRINT_FILE_FREQUENCY << ".plt";
	}

	ifstream fin;
	fin.open( filename.str() );

	for (int i = 0; i < 2; i++){
		fin.ignore( 1000, '\n' );
	}

	string dummy;
	int numNodes, numElements;
	fin >> dummy >> dummy >> dummy >> numNodes >> dummy >> dummy >> dummy >> numElements;
	
	for (int i = 0; i < 3; i++){
		fin.ignore(numeric_limits<streamsize>::max(), '\n');
	}

	vector<PRECISION> rho;
	vector<PRECISION> u;
	vector<PRECISION> v;
	vector<PRECISION> E;
	PRECISION var;

	for (int i = 0; i < numElements; i++){
		fin >> var;
		rho.push_back(var);
	}

	for (int i = 0; i < numElements; i++){
		fin >> var;
		u.push_back(var);
	}

	for (int i = 0; i < numElements; i++){
		fin >> var;
		v.push_back(var);
	}

	fin.ignore(numeric_limits<streamsize>::max(), '\n');
	fin.ignore(numeric_limits<streamsize>::max(), '\n');

	for (int i = 0; i < numElements; i++){
		fin >> var;
		E.push_back(var);
	}

	for( size_t icell=0; icell < cells_cfd.size(); icell++ ){
		cells_cfd[icell].vars.q_new[0] = rho[icell];
		cells_cfd[icell].vars.q_new[1] = rho[icell] * u[icell];
		cells_cfd[icell].vars.q_new[2] = rho[icell] * v[icell];
		cells_cfd[icell].vars.q_new[3] = rho[icell] * E[icell];
	}
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Make preparations for the time-step integration (q_old->q_new, etc...)
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
void CFDv0_solver<PRECISION, DIM_CNT, FACE_CNT> :: prepare_for_timestep( int rk_step ){

	// Copy new ----> old
	for( size_t icell=0; icell < cells_cfd.size(); icell++ ){
		for( unsigned idim=0; idim < DIM_CNT+2; idim++ ){
			cells_cfd[icell].vars.q_old[idim]   = cells_cfd[icell].vars.q_new  [idim];
			cells_cfd[icell].vars.delta_q[idim] *= Ak[rk_step];
		}

		// Initialize vars to zero
		init_array( cells_cfd[icell].vars.rho_grad  );
		init_array( cells_cfd[icell].vars.rhoU_grad );
		init_array( cells_cfd[icell].vars.rhoE_grad );
		init_array( cells_cfd[icell].vars.Rpsi_grad );
		init_array( cells_cfd[icell].vars.c_grad    );
		init_array( cells_cfd[icell].vars.dudx      );
		init_array( cells_cfd[icell].vars.dTdx      );
		init_array( cells_cfd[icell].vars.tauMC     );
	}

	if( sm_id != INDEX_BND_SUBMESH )
		return;

	// MPI ghost cells
	for( size_t ineigh=0; ineigh < ghost_mpi.size(); ineigh++ ){
		for( size_t icell=0; icell < ghost_mpi[ineigh].size(); icell++ ){
			//ghost_mpi[ineigh][icell].vars.sm_id    = 10000;
			//ghost_mpi[ineigh][icell].vars.is_ghost = true;
			for( unsigned idim=0; idim < DIM_CNT+2; idim++ ){
				ghost_mpi[ineigh][icell].vars.q_old[idim] = ghost_mpi[ineigh][icell].vars.q_new[idim];
			}
		}
	}

	// Boundary ghost cells
	for( size_t ibound=0; ibound < ghost_bnd.size(); ibound++ ){
		for( size_t icell=0; icell < ghost_bnd[ibound].size(); icell++ ){
			//ghost_bnd[ibound][icell].vars.sm_id    = 10000;
			//ghost_bnd[ineigh][icell].vars.is_ghost = true;
			for( unsigned idim=0; idim < DIM_CNT+2; idim++ ){
				ghost_bnd[ibound][icell].vars.q_old[idim] = ghost_bnd[ibound][icell].vars.q_new[idim];
			}
		}
	}
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Calculate needed gradients for cell to face interpolation schemes
// 	- only linear interpolation
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
void CFDv0_solver<PRECISION, DIM_CNT, FACE_CNT> :: calc_gradients( MPI_env &mpi_env ) {

	for( size_t icell=cells_cfd.start(); icell < cells_cfd.end(); icell++ ){
		auto *cell = &cells_cfd[icell];
		auto *vars = &cells_cfd[icell].vars;

		// Cell primitive vars
		PRECISION cell_Rpsi = compute_Rpsi( vars->q_old );
		PRECISION cell_c    = sqrt( gamma * cell_Rpsi );

		for( unsigned iface=0; iface < FACE_CNT; iface++ ){
			auto *neigh = cells_cfd[icell].neighs[iface];

			if( neigh->sm_id < sm_id )
				continue;
			if( neigh->sm_id == sm_id && neigh->id < icell )
				continue;

			// Neighbor primitive vars
			PRECISION adjc_Rpsi = compute_Rpsi( neigh->q_old );
			PRECISION adjc_c    = sqrt( gamma * adjc_Rpsi );

			// Fluxes
			PRECISION rho = interp_linear( cell->weight_linear[iface],	vars ->q_old[0],
																		neigh->q_old[0] );
			PRECISION rhoU[DIM_CNT];
			for( unsigned idim = 0; idim < DIM_CNT; idim++ ){
				rhoU[idim] = interp_linear( cell->weight_linear[iface], vars ->q_old[idim+1],
																		neigh->q_old[idim+1] );
			}
			PRECISION rhoE = interp_linear( cell->weight_linear[iface], vars ->q_old[DIM_CNT+1],
																		neigh->q_old[DIM_CNT+1] );

			PRECISION Rpsi = interp_linear( cell->weight_linear[iface], cell_Rpsi, adjc_Rpsi );
			PRECISION c    = interp_linear( cell->weight_linear[iface], cell_c,    adjc_c    );

			// Cell and neighbor (surface / vol );
			PRECISION cell_surf_over_vol[DIM_CNT];
			PRECISION adjc_surf_over_vol[DIM_CNT];
			for( unsigned idim = 0; idim < DIM_CNT; idim++ ){
				cell_surf_over_vol[idim] =   cell->S[iface][idim] * vars ->vol_inv;
				adjc_surf_over_vol[idim] = - cell->S[iface][idim] * neigh->vol_inv;
			}

			// Cell and neighbors' gradients
			for( unsigned idim = 0; idim < DIM_CNT; idim++ ){
				vars->rho_grad [idim] += rho  * cell_surf_over_vol[idim];
				for( unsigned jdim = 0; jdim < DIM_CNT; jdim++ ){
					vars->rhoU_grad[jdim][idim] += rhoU[jdim] * cell_surf_over_vol[idim];
				}
				vars->rhoE_grad[idim] += rhoE * cell_surf_over_vol[idim];
				vars->Rpsi_grad[idim] += Rpsi * cell_surf_over_vol[idim];
				vars->c_grad   [idim] += c    * cell_surf_over_vol[idim];

				neigh->rho_grad [idim] += rho  * adjc_surf_over_vol[idim];
				for( unsigned jdim = 0; jdim < DIM_CNT; jdim++ ){
					neigh->rhoU_grad[jdim][idim] += rhoU[jdim] * adjc_surf_over_vol[idim];
				}
				neigh->rhoE_grad[idim] += rhoE * adjc_surf_over_vol[idim];
				neigh->Rpsi_grad[idim] += Rpsi * adjc_surf_over_vol[idim];
				neigh->c_grad   [idim] += c    * adjc_surf_over_vol[idim];
			}
		}
	}
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Calculate D2U for artificial diffusion
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
void CFDv0_solver<PRECISION, DIM_CNT, FACE_CNT> :: calc_AD(){

	////
	//// Artificial Diffusion
	////
	//for( size_t icell=cells_cfd.start(); icell < cells_cfd.end(); icell++ ){
	//	auto *cell = &cells_cfd[icell];

	//	for( unsigned idim=0; idim < DIM_CNT+2; idim++ ){
	//		cell->vars.D2U[idim] = 0.0;

	//		for( unsigned iface=0; iface < FACE_CNT; iface++ ){
	//			auto *neigh = cells_cfd[icell].neighs[iface];

	//			cell->vars.D2U[idim] += neigh->q_old[idim];
	//		}
	//		cell->vars.D2U[idim] -= FACE_CNT * cell->vars.q_old[idim];
	//	}
	//}

}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Calculate Derivatives for Viscosity
// 		- interpolation: linear or midPoint
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
void CFDv0_solver<PRECISION, DIM_CNT, FACE_CNT> :: calc_VIS( MPI_env &mpi_env ){

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Viscosity
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	for( size_t icell=cells_cfd.start(); icell < cells_cfd.end(); icell++ ){
		auto *cell = &cells_cfd[icell];
		auto *vars = &cells_cfd[icell].vars;

		// Cell primitive vars
		PRECISION cell_Rpsi = compute_Rpsi( vars->q_old );

		for( unsigned iface=0; iface < FACE_CNT; iface++ ){
			auto *neigh = cells_cfd[icell].neighs[iface];

			if( neigh->sm_id < sm_id )
				continue;
			if( neigh->sm_id == sm_id && neigh->id < icell )
				continue;

			// Neighbor primitive vars
			PRECISION adjc_Rpsi = compute_Rpsi( neigh->q_old );

			// Fluxes
			PRECISION face_rho_inv = ONE / interp_linear( cell->weight_linear[iface], vars ->q_old[0],
																					  neigh->q_old[0] );
			PRECISION face_rhoU[DIM_CNT];
			for( unsigned idim = 0; idim < DIM_CNT; idim++ ){
				face_rhoU[idim] = interp_linear( cell->weight_linear[iface], vars ->q_old[idim+1],
																			 neigh->q_old[idim+1] );
			}
			PRECISION face_Rpsi = interp_linear( cell->weight_linear[iface], cell_Rpsi, adjc_Rpsi );
			PRECISION face_T    = face_Rpsi * Rgas_inv;

			PRECISION face_U[DIM_CNT];
			for( unsigned idim = 0; idim < DIM_CNT; idim++ ){
				face_U[idim] = face_rhoU[idim] * face_rho_inv;
			}

			// Cell and neighbor (surface vector / vol )
			PRECISION cell_surf_over_vol[DIM_CNT];
			PRECISION adjc_surf_over_vol[DIM_CNT];
			for( unsigned idim = 0; idim < DIM_CNT; idim++ ){
				cell_surf_over_vol[idim] =  cell->S[iface][idim] * vars ->vol_inv ;
				adjc_surf_over_vol[idim] = -cell->S[iface][idim] * neigh->vol_inv ;
			}

			for( unsigned idim = 0; idim < DIM_CNT; idim++ ){
				for( unsigned jdim = 0; jdim < DIM_CNT; jdim++ ){
					vars ->dudx[idim][jdim] += face_U[idim] * cell_surf_over_vol[jdim];
					neigh->dudx[idim][jdim] += face_U[idim] * adjc_surf_over_vol[jdim];
				}
			}

			for( unsigned idim=0; idim < DIM_CNT; idim++ ){
				vars ->dTdx[idim] += face_T * cell_surf_over_vol[idim] ;
				neigh->dTdx[idim] += face_T * adjc_surf_over_vol[idim] ;
			}

			for( unsigned idim = 0; idim < DIM_CNT; idim++ ){
				for( unsigned jdim = 0; jdim < DIM_CNT; jdim++ ){
					vars ->tauMC[idim][jdim] += mu0 * face_U[jdim] * cell_surf_over_vol[idim];
					neigh->tauMC[idim][jdim] += mu0 * face_U[jdim] * adjc_surf_over_vol[idim];
					if ( idim == jdim ) {
						vars ->tauMC[idim][jdim] -= mu0 * 2.0 / 3.0 * dot_product( face_U, cell_surf_over_vol );
						neigh->tauMC[idim][jdim] -= mu0 * 2.0 / 3.0 * dot_product( face_U, adjc_surf_over_vol );
					}
				}
			}
		}
	}
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Perform numerical scheme for the Hpath and regular cells
// 	- interpolation: linear (convection) or midPoint (viscous)
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
void CFDv0_solver<PRECISION, DIM_CNT, FACE_CNT> :: one_rk_step_v1( int rk_step, PRECISION dt, MPI_env &mpi_env ){

	for( size_t icell=cells_cfd.start(); icell < cells_cfd.end(); icell++ ){

		auto *cell  = &cells_cfd[icell].vars;
		auto *param = &cells_cfd[icell];

		// Cell primitive variables
		PRECISION cell_rho_inv = ONE / cell->q_old[0];
		PRECISION cell_Uvec[DIM_CNT];
		for( unsigned idim = 0; idim < DIM_CNT; idim++ ){
			cell_Uvec[idim] = cell->q_old[idim+1] * cell_rho_inv ;
		}
		PRECISION cell_Rpsi = compute_Rpsi( cell->q_old );
		PRECISION cell_T = cell_Rpsi * Rgas_inv;

		// Neighbors contribution
		for( unsigned iface=0; iface < FACE_CNT; iface++ ){
			// Fetch neighbor solution vars
			auto *neigh = param->neighs[iface];

			if( neigh->sm_id < sm_id )
				continue;
			if( neigh->sm_id == sm_id && neigh->id < icell )
				continue;

			PRECISION rhs[DIM_CNT+2];

			// Neighbor primitive variables
			PRECISION adjc_rho_inv = ONE / neigh->q_old[0];
			PRECISION adjc_Uvec[DIM_CNT];
			for( unsigned idim = 0; idim < DIM_CNT; idim++ ){
				adjc_Uvec[idim] = neigh->q_old[idim+1] * adjc_rho_inv ;
			}
			PRECISION adjc_Rpsi = compute_Rpsi( neigh->q_old );
			PRECISION adjc_T = adjc_Rpsi * Rgas_inv;

			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			// Convection
			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

			// Fluxes
			PRECISION rhoPos = interp_minmod( cell->q_old[0], neigh->q_old[0], cell ->rho_grad, param->d[iface], param->weight_linear[iface],   ONE );
			PRECISION rhoNeg = interp_minmod( cell->q_old[0], neigh->q_old[0], neigh->rho_grad, param->d[iface], param->weight_linear[iface], M_ONE );

			PRECISION rhoPos_inv = ONE / rhoPos;
			PRECISION rhoNeg_inv = ONE / rhoNeg;

			PRECISION rhoUPos[DIM_CNT];
			PRECISION rhoUNeg[DIM_CNT];
			for( unsigned idim = 0; idim < DIM_CNT; idim++ ){
				rhoUPos[idim] = interp_minmod( cell->q_old[idim+1], neigh->q_old[idim+1], cell ->rhoU_grad[idim], param->d[iface], param->weight_linear[iface],   ONE );
				rhoUNeg[idim] = interp_minmod( cell->q_old[idim+1], neigh->q_old[idim+1], neigh->rhoU_grad[idim], param->d[iface], param->weight_linear[iface], M_ONE );
			}
			PRECISION rhoEPos = interp_minmod( cell->q_old[DIM_CNT+1], neigh->q_old[DIM_CNT+1], cell ->rhoE_grad, param->d[iface], param->weight_linear[iface],   ONE );
			PRECISION rhoENeg = interp_minmod( cell->q_old[DIM_CNT+1], neigh->q_old[DIM_CNT+1], neigh->rhoE_grad, param->d[iface], param->weight_linear[iface], M_ONE );

			PRECISION RpsiPos = interp_minmod( cell_Rpsi, adjc_Rpsi, cell ->Rpsi_grad, param->d[iface], param->weight_linear[iface],   ONE );
			PRECISION RpsiNeg = interp_minmod( cell_Rpsi, adjc_Rpsi, neigh->Rpsi_grad, param->d[iface], param->weight_linear[iface], M_ONE );

			PRECISION pPos = rhoPos * RpsiPos;
			PRECISION pNeg = rhoNeg * RpsiNeg;

			PRECISION cP = sqrt(gamma * cell_Rpsi );
			PRECISION cN = sqrt(gamma * adjc_Rpsi );

			PRECISION cPos = interp_minmod( cP, cN, cell ->c_grad, param->d[iface], param->weight_linear[iface],   ONE );
			PRECISION cNeg = interp_minmod( cP, cN, neigh->c_grad, param->d[iface], param->weight_linear[iface], M_ONE );

			PRECISION phiPos = 0.0;
			PRECISION phiNeg = 0.0;
			PRECISION uPos[DIM_CNT] = {0.0};
			PRECISION uNeg[DIM_CNT] = {0.0};

			for( unsigned idim = 0; idim < DIM_CNT; idim++ ){
				uPos[idim] = rhoUPos[idim] * rhoPos_inv ;
				uNeg[idim] = rhoUNeg[idim] * rhoNeg_inv ;

				phiPos +=   uPos[idim] * param->S[iface][idim];
				phiNeg += - uNeg[idim] * param->S[iface][idim];
			}

			PRECISION S_mag = 0.;
			for( unsigned idim=0; idim < DIM_CNT; idim++ ){
				S_mag += pow( cells_cfd[icell].S[iface][idim], 2 );
			}
			S_mag = sqrt( S_mag );

			PRECISION psiPos = max( { phiPos + cPos * S_mag, - phiNeg + cNeg * S_mag, ZERO } );
			PRECISION psiNeg = min( { phiPos - cPos * S_mag, - phiNeg - cNeg * S_mag, ZERO } );

			// RHS contribution
			PRECISION a0 = ONE / (psiPos - psiNeg );
			PRECISION a1 = psiPos * psiNeg ;
			PRECISION aPos = psiPos * phiPos ;
			PRECISION aNeg = psiNeg * phiNeg ;

			rhs[0] = - (aPos*rhoPos + aNeg*rhoNeg + (rhoNeg - rhoPos) * a1 ) * a0;
			for( unsigned idim=0; idim < DIM_CNT; idim++ ){
				rhs[idim+1] = - (aPos*rhoUPos[idim] + aNeg*rhoUNeg[idim] + (rhoUNeg[idim] - rhoUPos[idim]) * a1 ) * a0 - (pPos * psiPos - pNeg * psiNeg) * a0 * param->S[iface][idim];
			}
			rhs[DIM_CNT+1] = - ( aPos*rhoEPos + aNeg*rhoENeg + (rhoENeg - rhoEPos) * a1 + (aPos*pPos + aNeg*pNeg) ) * a0;

			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			// Viscosity
			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

			// Derivatives
			PRECISION mu = mu0;

			PRECISION dudx[DIM_CNT][DIM_CNT];
			PRECISION dTdx[DIM_CNT];

			PRECISION d_mag = vector_mag( param->d[iface] );
			PRECISION dmag_inv = ONE / d_mag ;
			PRECISION d_norm[DIM_CNT];
			d_norm[0] = param->d[iface][0] * dmag_inv ;
			d_norm[1] = param->d[iface][1] * dmag_inv ;

			PRECISION dudx_avg[DIM_CNT][DIM_CNT];
			for( unsigned idim=0; idim < DIM_CNT; idim++ ){
			for( unsigned jdim=0; jdim < DIM_CNT; jdim++ ){
				dudx_avg[idim][jdim] = interp_linear( param->weight_linear[iface],
												cell->dudx[idim][jdim], neigh->dudx[idim][jdim] );
			}
			}

			for( unsigned idim=0; idim < DIM_CNT; idim++ ){
				PRECISION tmp = dot_product( dudx_avg[idim], d_norm ) - (adjc_Uvec[idim] - cell_Uvec[idim]) * dmag_inv;
				for( unsigned jdim=0; jdim < DIM_CNT; jdim++ ){
					dudx[idim][jdim] = dudx_avg[idim][jdim] - tmp * d_norm[jdim] ;
				}
			}

			PRECISION dTdx_avg[DIM_CNT];
			for( unsigned idim=0; idim < DIM_CNT; idim++ ){
				dTdx_avg[idim] = interp_linear( param->weight_linear[iface], cell->dTdx[idim], neigh->dTdx[idim] );
			}
			PRECISION tmp = dot_product(dTdx_avg, d_norm) - (adjc_T - cell_T) * dmag_inv;
			for( unsigned idim=0; idim < DIM_CNT; idim++ ){
				dTdx[idim] = dTdx_avg[idim] - tmp * d_norm[idim] ;
			}

			// Stress tensor
			PRECISION tau[DIM_CNT][DIM_CNT];
			tau[0][0] = 2.0 / 3.0 * mu * (2.0 * dudx[0][0] - dudx[1][1]);
			tau[0][1] = mu * (dudx[0][1] + dudx[1][0]);
			tau[1][0] = tau[0][1];
			tau[1][1] = 2.0 / 3.0 * mu * (2.0 * dudx[1][1] - dudx[0][0]);

			// Head conductivity
			PRECISION Cp = gamma * Rgas * gamma_m_one_inv;
			PRECISION k  = Cp * mu * Pr_inv;

			PRECISION q_x[DIM_CNT];
			for( unsigned idim=0; idim < DIM_CNT; idim++ ){
				q_x[idim] = - k * dTdx[idim];
			}

			// div(tau)
			PRECISION divTauMC[DIM_CNT];
			PRECISION tauMC[DIM_CNT][DIM_CNT];

			for( unsigned idim=0; idim < DIM_CNT; idim++ ){
			for( unsigned jdim=0; jdim < DIM_CNT; jdim++ ){
				tauMC[idim][jdim] = interp_linear( param->weight_linear[iface], cell ->tauMC[idim][jdim],
																				neigh->tauMC[idim][jdim] );
			}
			}
			for( unsigned idim=0; idim < DIM_CNT; idim++ ){
				divTauMC[idim] = dot_product( tauMC[idim], param->S[iface] );
			}

			// Laplacians
			PRECISION laplacianU[DIM_CNT];
			for( unsigned idim=0; idim < DIM_CNT; idim++ ){
				laplacianU[idim] = mu * dot_product( dudx[idim], param->S[iface] );
			}
			PRECISION laplacianT = - dot_product( q_x, param->S[iface] );

			// div(sigmaU)
			PRECISION sigmaU[DIM_CNT];
			for( unsigned idim=0; idim < DIM_CNT; idim++ ){
				sigmaU[idim] = dot_product( cell_Uvec, tau[idim] );
			}
			PRECISION divSigmaU = dot_product( sigmaU, param->S[iface] );

			// RHS contribution
			for( unsigned idim=0; idim < DIM_CNT; idim++ ){
				rhs[idim+1] += divTauMC[idim] + laplacianU[idim];
			}
			rhs[DIM_CNT+1] += divSigmaU + laplacianT;

			for( unsigned idim = 0; idim < DIM_CNT+2; idim++ ){
				cell ->delta_q[idim] += rhs[idim];
				neigh->delta_q[idim] -= rhs[idim];
			}
		}

		// 6- Compute next RK step
		for( unsigned idim = 0; idim < DIM_CNT+2; idim++ ){
			cell->q_new[idim] = cell->q_old[idim] + Bk[rk_step] * dt * cell->delta_q[idim] * cell->vol_inv;
		}
	}
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Exchange ghosts cells (Solal)
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
void CFDv0_solver<PRECISION, DIM_CNT, FACE_CNT> :: exchange_ghost_cells( MPI_env &mpi_env ){

	// Serial case
	if( mpi_env.is_serial() ){
		// TODO: Add flag periodic
		// 	- for now, we assume ONLY periodic BCs
		// 	- something like ---> if( mpi_env.periodic ){ do....}
		if ( !ghost_mpi.empty() && !ghost_mpi[0].empty() ) {
			for( size_t icell=0; icell < cells_cfd.size(); icell++ ){
				ghost_mpi[0][icell].vars = cells_cfd[icell].vars;
			}
		}
		return;
	}

	// Depending on framework used...
	switch( mpi_env.get_halo_comm_type() ){
		case MPI_TWOSIDED_NONB:
			mpi_env.isendrecv( cells_cfd, ghost_mpi );
			break;
		case MPI_TWOSIDED_PERS:
		case MPI_ONESIDED_BLCK:
		case MPI_ONESIDED_NONB:
		case MPI_ONESIDED_PERS:
			break;
		case MPI_NEIGHCOLL_BLCK:
			mpi_env.neighbor_alltoallw( &cells_cfd[0], &ghost_mpi[0][0] );
			break;
		case MPI_NEIGHCOLL_NONB:
			mpi_env.ineighbor_alltoallw( &cells_cfd[0], &ghost_mpi[0][0] );
			break;
		case MPI_NEIGHCOLL_PERS:
			mpi_env.neighbor_alltoallw_start();
			break;
		case MPI_NONE:
			cout << "MPI halo communication type has not been set. Terminate simulation." << endl;
			MPI_Abort( MPI_COMM_WORLD, 910 );
		default:
			break;
	}
}

//***************************************************************************************************
/*
template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
PRECISION CFDv0_solver<PRECISION, DIM_CNT, FACE_CNT> :: compute_cfl( PRECISION dt ){

	PRECISION cfl, cfl_max = 0.;
	for( size_t icell=0; icell < cells_cfd.size(); icell++ ){

		cfl = 0;
		PRECISION rho_inv = ONE / cells_cfd[icell].vars.q_old[0];

		for( unsigned iface=0; iface < FACE_CNT; iface++ ){
			PRECISION flux = 0.;
			for( unsigned idim=0; idim < DIM_CNT; idim++ ){
				flux += cells_cfd[icell].vars.q_old[idim+1] * cells_cfd[icell].S[iface][idim];
			}
			cfl += abs(flux);
		}
		cfl *= rho_inv * cells_cfd[icell].vars.vol_inv * dt;

		if( cfl_max < cfl )
			cfl_max = cfl;
	}

	return cfl_max;
}*/

//***************************************************************************************************
template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
PRECISION CFDv0_solver<PRECISION, DIM_CNT, FACE_CNT> :: compute_dt( PRECISION cflMax ){

	PRECISION sumFlux, dt, dt_min=100.0;
	for( size_t icell=0; icell < cells_cfd.size(); icell++ ){

		sumFlux = 0;
		PRECISION rho = cells_cfd[icell].vars.q_old[0];

		for( unsigned iface=0; iface < FACE_CNT; iface++ ){
			PRECISION sMag = 0.;
			for( unsigned idim=0; idim < DIM_CNT; idim++ ){
				sMag += cells_cfd[icell].S[iface][idim] * cells_cfd[icell].S[iface][idim];
			}
			sMag = sqrt( sMag );
			
			PRECISION flux = 0.;
			for( unsigned idim=0; idim < DIM_CNT; idim++ ){
				flux += cells_cfd[icell].vars.q_old[idim+1] * cells_cfd[icell].S[iface][idim] / sMag;
			}
			sumFlux += abs(flux);
		}
		PRECISION cells_cfd_vol = 1.0 / cells_cfd[icell].vars.vol_inv;
		dt = cflMax * rho * cells_cfd_vol / sumFlux;

		if( dt_min > dt )
			dt_min = dt;
	}

	return dt_min;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
void CFDv0_solver<PRECISION, DIM_CNT, FACE_CNT> :: gather_info(
									MPI_env									 &mpi_env,
									struct t_domain_info<PRECISION, DIM_CNT> &domain_info ){

	domain_info.rho_min = cells_cfd[0].vars.q_new[0];
	domain_info.rho_max = cells_cfd[0].vars.q_new[0];

	for( size_t icell=0; icell < cells_cfd.size(); icell++ ){

		if( domain_info.rho_min > cells_cfd[icell].vars.q_new[0] )
			domain_info.rho_min = cells_cfd[icell].vars.q_new[0];

		if( domain_info.rho_max < cells_cfd[icell].vars.q_new[0] )
			domain_info.rho_max = cells_cfd[icell].vars.q_new[0];
	}
}

//***************************************************************************************************
// Save the results
//***************************************************************************************************
template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
void CFDv0_solver<PRECISION, DIM_CNT, FACE_CNT> :: finalize() {

	ofstream output;
	output.open("Results.txt");
	for( size_t icell=0; icell < cells_cfd.size(); icell++ ) {

	  for( size_t idim=0; idim < DIM_CNT; idim++ ) {
	    output << cells_cfd[icell].x[idim] << "\t";
	  }

	  PRECISION rho = cells_cfd[icell].vars.q_new[0];
	  PRECISION u[DIM_CNT], u_mag = 0.0;
	  for( unsigned idim = 0; idim < DIM_CNT; idim++ ){
	    u[idim] = cells_cfd[icell].vars.q_new[idim+1] / cells_cfd[icell].vars.q_new[0];
	    u_mag += u[idim] * u[idim];
	  }
	  u_mag = sqrt(u_mag);
	  PRECISION E = cells_cfd[icell].vars.q_new[DIM_CNT+1] / cells_cfd[icell].vars.q_new[0];
	  PRECISION p = rho * (gamma - ONE) * (E - 0.5 * u_mag * u_mag);

	  output << rho << "\t";

	  for( size_t idim=0; idim < DIM_CNT; idim++ ) {
	    output << u[idim] << "\t";
	  }

	  output << u_mag << "\t" << E << "\t" << p << "\n";
	}

	for( size_t icell=0; icell < cells_cfd.size(); icell++ ) {

	  for( size_t idim=0; idim < DIM_CNT; idim++ ) {
	    output << cells_cfd[icell].x[idim] << "\t";
	  }

	  PRECISION rho = cells_cfd[icell].vars.q_new[0];
	  PRECISION u[DIM_CNT], u_mag = 0.0;
	  for( unsigned idim = 0; idim < DIM_CNT; idim++ ){
	    u[idim] = cells_cfd[icell].vars.q_new[idim+1] / cells_cfd[icell].vars.q_new[0];
	    u_mag += u[idim] * u[idim];
	  }
	  u_mag = sqrt(u_mag);
	  PRECISION E = cells_cfd[icell].vars.q_new[DIM_CNT+1] / cells_cfd[icell].vars.q_new[0];
	  PRECISION p = rho * (gamma - ONE) * (E - 0.5 * u_mag * u_mag);

	  output << rho << "\t";

	  for( size_t idim=0; idim < DIM_CNT; idim++ ) {
	    output << u[idim] << "\t";
	  }

	  output << u_mag << "\t" << E << "\t" << p << "\n";
	}
	output.close();
}

//***************************************************************************************************
template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
void CFDv0_solver<PRECISION, DIM_CNT, FACE_CNT> :: output_cfd( PRECISION time, std::string filename, gmsh_mesh *mesh ){

	MeshWriter mesh_writer;

	mesh_writer.output_cfd( time, filename, mesh, this->cells_cfd );
}































