#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <initializer_list>

#include "math.h"
#include "api/cfdv1_templates.h"
#include "api/cfdv1_solver.h"
#include "api/mesh_writer.h"

using namespace std;

//***************************************************************************************************
// Allocate cells and record addresses for neighbors' connectivity
//***************************************************************************************************
template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
void CFDv1_solver<PRECISION, DIM_CNT, FACE_CNT> :: allocate( MPI_env mpi_env,
															 gmsh_mesh *mesh,
															 std::vector<CFDv1_soln_vars<PRECISION,DIM_CNT>*> &addr_soln,
															 std::vector<CFDv1_conv_vars<PRECISION,DIM_CNT>*> &addr_conv,
															 std::vector<CFDv1_visc_vars<PRECISION,DIM_CNT>*> &addr_visc ){

	// Allocate
	cells_soln  .resize( mesh->cells.size() );
	cells_conv  .resize( mesh->cells.size() );
	cells_visc  .resize( mesh->cells.size() );
	cells_params.resize( mesh->cells.size() );

	addr_soln.resize( mesh->cells.size() );
	addr_visc.resize( mesh->cells.size() );
	addr_conv.resize( mesh->cells.size() );

	// Initialize Hpath cells
	for( size_t icell=0; icell < mesh->cells.size(); icell++ ){
		for( unsigned idim=0; idim < DIM_CNT+2; idim++ ){
			cells_soln[icell].q_old  [idim] = 0.;
			cells_soln[icell].q_new  [idim] = 0.;
			cells_soln[icell].delta_q[idim] = 0.;
		}
	}

	// Record local cells' memory addresses
	for( size_t icell=0; icell < mesh->cells.size(); icell++ ){
		addr_soln[icell] = &( cells_soln[icell] );
		addr_conv[icell] = &( cells_conv[icell] );
		addr_visc[icell] = &( cells_visc[icell] );
	}

	// Print info
	if( mpi_env.is_master() ){
		cout << "local cell number =" << setw(10) << sm_id
									  << setw(10) << mesh->cells.size()
									  << endl;
		cout << "    - size of CFD solution   cell =" << setw(10) << sizeof( cells_soln[0] ) << endl;
		cout << "    - size of CFD convection cell =" << setw(10) << sizeof( cells_conv[0] ) << endl;
		cout << "    - size of CFD viscous    cell =" << setw(10) << sizeof( cells_visc[0] ) << endl;
	}
}

//***************************************************************************************************
//***************************************************************************************************
template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
void CFDv1_solver<PRECISION, DIM_CNT, FACE_CNT> ::
				allocate_ghost_cells( MPI_env &mpi_env,
									  gmsh_mesh *mesh,
									  std::vector<std::vector<CFDv1_soln_vars<PRECISION,DIM_CNT>*>> &ptr_mpi_soln,
									  std::vector<std::vector<CFDv1_conv_vars<PRECISION,DIM_CNT>*>> &ptr_mpi_conv,
									  std::vector<std::vector<CFDv1_visc_vars<PRECISION,DIM_CNT>*>> &ptr_mpi_visc,
									  std::vector<std::vector<CFDv1_soln_vars<PRECISION,DIM_CNT>*>> &ptr_bnd_soln,
									  std::vector<std::vector<CFDv1_conv_vars<PRECISION,DIM_CNT>*>> &ptr_bnd_conv,
									  std::vector<std::vector<CFDv1_visc_vars<PRECISION,DIM_CNT>*>> &ptr_bnd_visc ){

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
	ghost_bnd_soln.resize( 4 );
	ghost_bnd_conv.resize( 4 );
	ghost_bnd_visc.resize( 4 );

	ptr_bnd_soln.resize( 4 );
	ptr_bnd_conv.resize( 4 );
	ptr_bnd_visc.resize( 4 );

	for( size_t ibound=0; ibound < ptr_bnd_soln.size(); ibound++ ){
		ghost_bnd_soln[ibound].resize( boundaries[ibound].size() );
		ghost_bnd_conv[ibound].resize( boundaries[ibound].size() );
		ghost_bnd_visc[ibound].resize( boundaries[ibound].size() );

		ptr_bnd_soln[ibound].resize( boundaries[ibound].size() );
		ptr_bnd_conv[ibound].resize( boundaries[ibound].size() );
		ptr_bnd_visc[ibound].resize( boundaries[ibound].size() );

		for( size_t icell=0; icell < ghost_bnd_soln[ibound].size(); icell++ ){
			for( unsigned idim=0; idim < DIM_CNT+2; idim++ ){
				ghost_bnd_soln[ibound][icell].q_old[idim] = 0.;
				ghost_bnd_soln[ibound][icell].q_new[idim] = 0.;
			}
		}
	}

	for( size_t ibound=0; ibound < ghost_bnd_soln.size(); ibound++ ){
		for( size_t icell=0; icell < ghost_bnd_soln[ibound].size(); icell++ ){
			ptr_bnd_soln[ibound][icell] = &( ghost_bnd_soln[ibound][icell] );
			ptr_bnd_conv[ibound][icell] = &( ghost_bnd_conv[ibound][icell] );
			ptr_bnd_visc[ibound][icell] = &( ghost_bnd_visc[ibound][icell] );

			ghost_bnd_conv[ibound][icell].sm_id = 10000;
			ghost_bnd_visc[ibound][icell].sm_id = 10000;

			if (ibound == 0){
				ghost_bnd_conv[ibound][icell].id = mesh->cells.size() + icell;
				ghost_bnd_visc[ibound][icell].id = mesh->cells.size() + icell;
			}else{
				ghost_bnd_conv[ibound][icell].id = mesh->cells.size() + ibound*ghost_bnd_conv[ibound-1].size() + icell;
				ghost_bnd_visc[ibound][icell].id = mesh->cells.size() + ibound*ghost_bnd_visc[ibound-1].size() + icell;
			}
		}
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Initialize MPI ghost cells
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	std::vector<int > mpi_neighbors;
	std::vector<bool> registered_rank( mpi_env.size(), false );

	for( size_t icell=0; icell < mesh->cells.size(); icell++ ){
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

		int local_size = (int) mesh->cells.size();

		MPI_Isend( &local_size, 1, MPI_INT, neigh_rank, send_tag, MPI_COMM_WORLD, &comm_reqs[ineigh] );
		MPI_Irecv( &ghost_size[ineigh], 1, MPI_INT, neigh_rank, recv_tag, MPI_COMM_WORLD, &comm_reqs[ineigh+mpi_neighbors.size()] );
	}
	MPI_Waitall( comm_reqs.size(), &comm_reqs[0], MPI_STATUS_IGNORE );

	// Allocate ghost cells
	ghost_mpi_soln.resize( mpi_neighbors.size() );
	ghost_mpi_conv.resize( mpi_neighbors.size() );
	ghost_mpi_visc.resize( mpi_neighbors.size() );

	ptr_mpi_soln.resize( mpi_neighbors.size() );
	ptr_mpi_conv.resize( mpi_neighbors.size() );
	ptr_mpi_visc.resize( mpi_neighbors.size() );

	for( size_t ineigh=0; ineigh < ghost_mpi_soln.size(); ineigh++ ){
			ghost_mpi_soln[ineigh].resize( ghost_size[ineigh] );
			ghost_mpi_conv[ineigh].resize( ghost_size[ineigh] );
			ghost_mpi_visc[ineigh].resize( ghost_size[ineigh] );

			ptr_mpi_soln[ineigh].resize( ghost_size[ineigh] );
			ptr_mpi_conv[ineigh].resize( ghost_size[ineigh] );
			ptr_mpi_visc[ineigh].resize( ghost_size[ineigh] );
	}

	// Store ghost cells' memory addresses
	for( size_t ineigh=0; ineigh < ghost_mpi_soln.size(); ineigh++ ){
		for( size_t icell=0; icell < ghost_mpi_soln[ineigh].size(); icell++ ){
			ptr_mpi_soln[ineigh][icell] = &( ghost_mpi_soln[ineigh][icell] );
			ptr_mpi_conv[ineigh][icell] = &( ghost_mpi_conv[ineigh][icell] );
			ptr_mpi_visc[ineigh][icell] = &( ghost_mpi_visc[ineigh][icell] );

			ghost_mpi_conv[ineigh][icell].sm_id = 10000; // TODO: change to MAX_UNSIGNED_SHORT
			ghost_mpi_visc[ineigh][icell].sm_id = 10000; // TODO: change to MAX_UNSIGNED_SHORT

			for( unsigned idim=0; idim < DIM_CNT; idim++ ){
				ghost_mpi_soln[ineigh][icell].q_old  [idim] = 0.0;
				ghost_mpi_soln[ineigh][icell].q_new  [idim] = 0.0;
				ghost_mpi_soln[ineigh][icell].delta_q[idim] = 0.0;
			}
		}
	}
}


//***************************************************************************************************
//***************************************************************************************************
template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
void CFDv1_solver<PRECISION, DIM_CNT, FACE_CNT> :: assign_pointers(
				gmsh_mesh *mesh,
				std::vector<std::vector<CFDv1_soln_vars<PRECISION, DIM_CNT>*>> ptr_cell_soln,
				std::vector<std::vector<CFDv1_conv_vars<PRECISION, DIM_CNT>*>> ptr_cell_conv,
				std::vector<std::vector<CFDv1_visc_vars<PRECISION, DIM_CNT>*>> ptr_cell_visc,
				std::vector<std::vector<CFDv1_soln_vars<PRECISION, DIM_CNT>*>> ptr_mpi_soln,
				std::vector<std::vector<CFDv1_conv_vars<PRECISION, DIM_CNT>*>> ptr_mpi_conv,
				std::vector<std::vector<CFDv1_visc_vars<PRECISION, DIM_CNT>*>> ptr_mpi_visc,
				std::vector<std::vector<CFDv1_soln_vars<PRECISION, DIM_CNT>*>> ptr_bnd_soln,
				std::vector<std::vector<CFDv1_conv_vars<PRECISION, DIM_CNT>*>> ptr_bnd_conv,
				std::vector<std::vector<CFDv1_visc_vars<PRECISION, DIM_CNT>*>> ptr_bnd_visc ){


	// Assign pointers
	for( size_t icell=0; icell < mesh->cells.size(); icell++ ){
		for( unsigned iface=0; iface < FACE_CNT; iface++ ){
			t_gmsh_neighbor *neigh = &( mesh->cell2neigh[icell][iface] );

			switch( neigh->type ){
				// Same submesh neighbor
				case CELL_PERIODIC:
				case CELL_REGULAR:
					neigh_soln[icell].neighs[iface] = ptr_cell_soln[this->sm_id][neigh->id];
					neigh_conv[icell].neighs[iface] = ptr_cell_conv[this->sm_id][neigh->id];
					neigh_visc[icell].neighs[iface] = ptr_cell_visc[this->sm_id][neigh->id];
					break;
				// Different submesh, same MPI rank
				case CELL_SUBMESH:
					neigh_soln[icell].neighs[iface] = ptr_cell_soln[neigh->sm][neigh->id];
					neigh_conv[icell].neighs[iface] = ptr_cell_conv[neigh->sm][neigh->id];
					neigh_visc[icell].neighs[iface] = ptr_cell_visc[neigh->sm][neigh->id];
					break;
				// All of the following neighbor types require ghost cells
				case CELL_NONE:
					neigh_soln[icell].neighs[iface] = ptr_bnd_soln[neigh->bound-1][neigh->id];
					neigh_conv[icell].neighs[iface] = ptr_bnd_conv[neigh->bound-1][neigh->id];
					neigh_visc[icell].neighs[iface] = ptr_bnd_visc[neigh->bound-1][neigh->id];
					break;
				case CELL_MPI:
				case CELL_PERIODIC_MPI:
					neigh_soln[icell].neighs[iface] = ptr_mpi_soln[rank2local[neigh->sm]][neigh->id];
					neigh_conv[icell].neighs[iface] = ptr_mpi_conv[rank2local[neigh->sm]][neigh->id];
					neigh_visc[icell].neighs[iface] = ptr_mpi_visc[rank2local[neigh->sm]][neigh->id];
					break;
				default:
					break;
			}
		}
	}
}

//***************************************************************************************************
//***************************************************************************************************
template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
void CFDv1_solver<PRECISION, DIM_CNT, FACE_CNT> :: initialize( MPI_env &mpi_env, gmsh_mesh *mesh ){

	// Initialize all parameters, mostly geometrical variables for the CFD solver
	init_params( mesh );

	// Initilize MPI environment: only for boundary submesh
	if( mesh->id == INDEX_BND_SUBMESH ){
		initialize_mpi_env( mpi_env, mesh );
	}

	//// Setup initial conditions
	//initial_condition init_case = CYLINDER_TUNNEL;
	//set_init_conditions( init_case );

	//if( mesh->id == INDEX_BND_SUBMESH ){
	//	// Setup boundary conditions, not always necessary
	//	set_boundary_conditions();
	//}
}

//***************************************************************************************************
//***************************************************************************************************
//									PRIVATE FUNCTIONS
//***************************************************************************************************
//***************************************************************************************************

//***************************************************************************************************
// Initialize the CFD solver
//***************************************************************************************************
template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
void CFDv1_solver<PRECISION, DIM_CNT, FACE_CNT> :: init_params( gmsh_mesh *mesh ){

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Hpath cells
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	// Cell's center
	xc.resize( mesh->cells.size(), std::vector<PRECISION>( DIM_CNT ));
	for( size_t icell=0; icell < mesh->cells.size(); icell++ ){
		for( unsigned idim = 0; idim < DIM_CNT; idim++ ){
			xc[icell][idim] = mesh->cells[icell].xc[idim] ;
		}
	}

	// Compute face surface vector
	for( size_t icell=0; icell < mesh->cells.size(); icell++ ){
		for( unsigned iface=0; iface < FACE_CNT; iface++ ){
			int gface = mesh->cell2face[icell][iface];
			for( unsigned idim = 0; idim < DIM_CNT; idim++ ){
				cells_params[icell].v_surf[iface][idim] = mesh->v_sign[icell][iface] * mesh->v_norm[gface][idim] ;
			}
		}
	}

	// Vector from face center to cell's center
	std::vector<std::vector<std::vector<PRECISION>>> dP( mesh->cells.size(), std::vector<std::vector<PRECISION>> ( FACE_CNT ) );
	for( size_t icell=0; icell < mesh->cells.size(); icell++ ){
		for( unsigned iface=0; iface < FACE_CNT; iface++ ){
			dP[icell][iface].resize( DIM_CNT, 0.0 );
		}
	}
	for( size_t icell=0; icell < mesh->cells.size(); icell++ ){
		for( unsigned iface=0; iface < FACE_CNT; iface++ ){
			int gface = mesh->cell2face[icell][iface];
			for( unsigned idim = 0; idim < DIM_CNT; idim++ ){
				dP[icell][iface][idim] = mesh->cells[icell].xc[idim] - mesh->faces[gface].xf[idim] ;
			}
		}
	}

	// Vector from cell center to neighbor's center
	for( size_t icell=0; icell < mesh->cells.size(); icell++ ){
		for( unsigned iface=0; iface < FACE_CNT; iface++ ){
			int gface = mesh->cell2face[icell][iface];
			auto *neigh = &( mesh->cell2neigh[icell][iface] );
			if( neigh->type == CELL_NONE ){
				PRECISION S_mag_sqrt = 0.0;
				PRECISION dP_dot_S = 0.0;
				for( unsigned idim = 0; idim < DIM_CNT; idim++ ){
					S_mag_sqrt += cells_params[icell].v_surf[iface][idim] * cells_params[icell].v_surf[iface][idim];
				}
				for( unsigned idim = 0; idim < DIM_CNT; idim++ ){
					dP_dot_S += dP[icell][iface][idim] * cells_params[icell].v_surf[iface][idim];
				}
				for( unsigned idim = 0; idim < DIM_CNT; idim++ ){
					cells_params[icell].v_neigh[iface][idim] = 2.0 * abs( dP_dot_S ) * cells_params[icell].v_surf[iface][idim] / S_mag_sqrt;
				}
			}else{
				for( unsigned idim = 0; idim < DIM_CNT; idim++ ){
					cells_params[icell].v_neigh[iface][idim] = mesh->v_sign[icell][iface] * mesh->v_neigh[gface][idim] ;
				}
			}
		}
	}

	// Vector from face center to neighbor's center
	std::vector<std::vector<std::vector<PRECISION>>> dN( mesh->cells.size(), std::vector<std::vector<PRECISION>> ( FACE_CNT ) );
	for( size_t icell=0; icell < mesh->cells.size(); icell++ ){
		for( unsigned iface=0; iface < FACE_CNT; iface++ ){
			dN[icell][iface].resize( DIM_CNT, 0.0 );
		}
	}
	for( size_t icell=0; icell < mesh->cells.size(); icell++ ){
		for( unsigned iface=0; iface < FACE_CNT; iface++ ){
			for( unsigned idim = 0; idim < DIM_CNT; idim++ ){
				dN[icell][iface][idim] = dP[icell][iface][idim] + cells_params[icell].v_neigh[iface][idim] ;
			}
		}
	}

	// beta_pos and beta_neg calculations
	for( size_t icell=0; icell < mesh->cells.size(); icell++ ){
		for( unsigned iface=0; iface < FACE_CNT; iface++ ){
			cells_params[icell].weight[iface] = 0.0;
			PRECISION divider = 0.0;
			for( unsigned idim = 0; idim < DIM_CNT; idim++ ){
				cells_params[icell].weight[iface] += cells_params[icell].v_surf[iface][idim] * dN[icell][iface][idim] ;
				divider += ( - cells_params[icell].v_surf[iface][idim] * dP[icell][iface][idim] + cells_params[icell].v_surf[iface][idim] * dN[icell][iface][idim] );
			}
			cells_params[icell].weight[iface] /= divider ;
		}
	}

	// Cell's volume
	for( size_t icell=0; icell < mesh->cells.size(); icell++ ){
		cells_params[icell].vol_inv = ONE / mesh->cells_vol[icell];
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

//***************************************************************************************************
//***************************************************************************************************
template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
void CFDv1_solver<PRECISION, DIM_CNT, FACE_CNT> :: initialize_mpi_env( MPI_env &mpi_env, gmsh_mesh *mesh ){

	// Get MPI neighbors
	std::vector<int > mpi_neighbors;
	std::vector<bool> registered_rank( mpi_env.size(), false );

	for( size_t icell=0; icell < mesh->cells.size(); icell++ ){
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

	//init_mpi_types( mpi_env );
	//mpi_env.set_halo_comm_reqs( 2*mpi_env.mpi_neighbors.size() );
	//mpi_env.set_halo_comm_type( MPI_TWOSIDED_NONB, cells_cfd, ghost_mpi );
}

////***************************************************************************************************
//// Initialize MPI types
////***************************************************************************************************
//template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
//void CFDv1_solver<PRECISION, DIM_CNT, FACE_CNT> :: init_mpi_types( MPI_env &mpi_env ){
//
//	if( mpi_env.size() == 1 ){
//		return;
//	}
//
//	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//	// Solution vars
//	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//	const int solvars_nfields = 3;
//	MPI_Aint  solvars_disps[solvars_nfields];
//	int       solvars_blocklens[] = { DIM_CNT+2, DIM_CNT+2, DIM_CNT+2 };
//
//	MPI_Datatype mpi_precision;
//	if( std::is_same< PRECISION, float >::value ){
//		mpi_precision = MPI_FLOAT;
//	}else if( std::is_same< PRECISION, double >::value ){
//		mpi_precision = MPI_DOUBLE;
//	}else{
//		MPI_Abort( MPI_COMM_WORLD, 254 );
//	}
//
//	MPI_Datatype solvars_types[] = { mpi_precision };
//
//	CFDv1_soln_vars<PRECISION, DIM_CNT> tmp_soln;
//
//	solvars_disps[0] = (char*) &tmp_soln.q_old   - (char*) &tmp_soln;
//	solvars_disps[1] = (char*) &tmp_soln.q_new   - (char*) &tmp_soln;
//	solvars_disps[2] = (char*) &tmp_soln.delta_q - (char*) &tmp_soln;
//
//	MPI_Datatype type_solvars;
//	MPI_Type_create_struct( solvars_nfields, solvars_blocklens, solvars_disps, solvars_types, &type_solvars );
//	MPI_Type_commit( &type_solvars );
//
//	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//	// Convection variables
//	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//	const int solvars_nfields = 7;
//	MPI_Aint  solvars_disps[solvars_nfields];
//	int       solvars_blocklens[] = { 1, 1, DIM_CNT, DIM_CNT*DIM_CNT, DIM_CNT, DIM_CNT, DIM_CNT };
//
//	MPI_Datatype mpi_precision;
//	if( std::is_same< PRECISION, float >::value ){
//		mpi_precision = MPI_FLOAT;
//	}else if( std::is_same< PRECISION, double >::value ){
//		mpi_precision = MPI_DOUBLE;
//	}else{
//		MPI_Abort( MPI_COMM_WORLD, 254 );
//	}
//
//	MPI_Datatype solvars_types[] = { mpi_precision };
//
//	CFDv1_conv_vars<PRECISION, DIM_CNT> tmp_conv;
//
//	solvars_disps[0] = (char*) &tmp_conv.sm_id     - (char*) &tmp_conv;
//	solvars_disps[1] = (char*) &tmp_conv.id        - (char*) &tmp_conv;
//	solvars_disps[2] = (char*) &tmp_conv.rho_grad  - (char*) &tmp_conv;
//	solvars_disps[3] = (char*) &tmp_conv.rhoU_grad - (char*) &tmp_conv;
//	solvars_disps[4] = (char*) &tmp_conv.rhoE_grad - (char*) &tmp_conv;
//	solvars_disps[5] = (char*) &tmp_conv.Rpsi_grad - (char*) &tmp_conv;
//	solvars_disps[6] = (char*) &tmp_conv.c_grad    - (char*) &tmp_conv;
//
//	MPI_Datatype type_conv;
//	MPI_Type_create_struct( solvars_nfields, solvars_blocklens, solvars_disps, solvars_types, &type_conv );
//	MPI_Type_commit( &type_conv );
//
//	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//	// Viscous variables
//	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//	const int solvars_nfields = 5;
//	MPI_Aint  solvars_disps[solvars_nfields];
//	int       solvars_blocklens[] = { 1, 1, DIM_CNT*DIM_CNT, DIM_CNT, DIM_CNT*DIM_CNT };
//
//	MPI_Datatype mpi_precision;
//	if( std::is_same< PRECISION, float >::value ){
//		mpi_precision = MPI_FLOAT;
//	}else if( std::is_same< PRECISION, double >::value ){
//		mpi_precision = MPI_DOUBLE;
//	}else{
//		MPI_Abort( MPI_COMM_WORLD, 254 );
//	}
//
//	MPI_Datatype solvars_types[] = { mpi_precision };
//
//	CFDv1_visc_vars<PRECISION, DIM_CNT> tmp_visc;
//
//	solvars_disps[0] = (char*) &tmp_visc.sm_id - (char*) &tmp_visc;
//	solvars_disps[1] = (char*) &tmp_visc.id    - (char*) &tmp_visc;
//	solvars_disps[2] = (char*) &tmp_visc.dudx  - (char*) &tmp_visc;
//	solvars_disps[3] = (char*) &tmp_visc.dTdx  - (char*) &tmp_visc;
//	solvars_disps[4] = (char*) &tmp_visc.tauMC - (char*) &tmp_visc;
//
//	MPI_Datatype type_visc;
//	MPI_Type_create_struct( solvars_nfields, solvars_blocklens, solvars_disps, solvars_types, &type_visc );
//	MPI_Type_commit( &type_visc );
//
//	// Save types
//	mpi_env.set_type_soln( type_soln );
//	mpi_env.set_type_conv( type_conv );
//	mpi_env.set_type_visc( type_visc );
//}
//
////***************************************************************************************************
////***************************************************************************************************
//template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
//void CFDv1_solver<PRECISION, DIM_CNT, FACE_CNT> :: set_init_conditions(initial_condition init_case){
//
//	for( size_t icell=0; icell < cells_soln.size(); icell++ ){
//		for( unsigned idim = 0; idim < DIM_CNT+2; idim++ ){
//			cells_soln[icell].q_old  [idim] = ZERO;
//			cells_soln[icell].q_new  [idim] = ZERO;
//			cells_soln[icell].delta_q[idim] = ZERO;
//		}
//	}
//
//	std::vector<PRECISION> rho, u[DIM_CNT], E;
//
//	rho.resize( cells_soln.size() );
//	for( unsigned idim = 0; idim < DIM_CNT; idim++ ){
//		u[idim].resize( cells_soln.size() );
//	}
//	E.resize( cells_soln.size() );
//
//	switch (init_case) {
//		case ISENTROPIC_VORTEX: {
//		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//		// Setup case/initial conditions: 2D ISENTROPIC VORTEX
//		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//			PRECISION alpha = 0.0;
//			PRECISION beta  = 5.0;
//			PRECISION mesh_xc[2] = {0.0, 0.0};
//
//			for( size_t icell=0; icell < cells_soln.size(); icell++ ){
//				PRECISION r = sqrt( pow( xc[icell][0] - mesh_xc[0], 2)
//								  + pow( xc[icell][1] - mesh_xc[1], 2) );
//
//				rho[icell] = pow(ONE - beta * beta * (gamma - ONE) * exp(ONE - r * r) / (8.0 * gamma * M_PI * M_PI), ONE / (gamma - ONE));
//
//				u[0][icell] = M * cos(alpha * M_PI / 180.0) - beta * ( xc[icell][1] - mesh_xc[1]) * exp(0.5 * (ONE - r * r)) / (2.0 * M_PI);
//				u[1][icell] = M * sin(alpha * M_PI / 180.0) + beta * ( xc[icell][0] - mesh_xc[0]) * exp(0.5 * (ONE - r * r)) / (2.0 * M_PI);
//
//				PRECISION u_mag = 0.0;
//				for( unsigned idim = 0; idim < DIM_CNT; idim++ ){
//					u_mag += u[idim][icell] * u[idim][icell];
//				}
//				u_mag = sqrt(u_mag);
//
//				PRECISION p = pow(rho[icell], gamma);
//
//				E[icell] = p / (rho[icell] * (gamma - ONE)) + 0.5 * u_mag * u_mag;
//			}
//			break;
//		}
//		case DOUBLE_SHEAR_LAYER_1: {
//		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//		// Setup case/initial conditions: 2D DOUBLE SHEAR LAYER
//		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//			for( size_t icell=0; icell < cells_soln.size(); icell++ ){
//				if (xc[icell][1] >= 0.25 || xc[icell][1] < -0.25) {
//					rho[icell]  =  1.0;
//					u[0][icell] =  1.0;
//				} else {
//					rho[icell]  =  2.0;
//					u[0][icell] = -1.0;
//				}
//				u[1][icell] =  0.01 * sin(2.0 * M_PI * 2.0 * xc[icell][0]);
//
//				PRECISION u_mag = 0.0;
//				for( unsigned idim = 0; idim < DIM_CNT; idim++ ){
//					u_mag += u[idim][icell] * u[idim][icell];
//				}
//				u_mag = sqrt(u_mag);
//
//				PRECISION p = 2.5;
//
//				E[icell] = p / (rho[icell] * (gamma - ONE)) + 0.5 * u_mag * u_mag;
//			}
//			break;
//		}
//		case DOUBLE_SHEAR_LAYER_2: {
//		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//		// Setup case/initial conditions: 2D DOUBLE SHEAR LAYER (Another version)
//		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//			for( size_t icell=0; icell < cells_soln.size(); icell++ ){
//				rho[icell] = 1.0;
//
//				if (xc[icell][1] <= 0.5){
//					u[0][icell] = tanh(80.0 * (xc[icell][1] - 0.25));
//				} else {
//					u[0][icell] = tanh(80.0 * (0.75 - xc[icell][1]));
//				}
//				u[1][icell] = 0.05 * sin(2.0 * M_PI * (xc[icell][0] + 0.25));
//
//				PRECISION u_mag = 0.0;
//				for( unsigned idim = 0; idim < DIM_CNT; idim++ ){
//					u_mag += u[idim][icell] * u[idim][icell];
//				}
//				u_mag = sqrt(u_mag);
//
//				PRECISION p = pow(rho[icell], 1.4);
//
//				E[icell] = p / (rho[icell] * (gamma - ONE)) + 0.5 * u_mag * u_mag;
//			}
//			break;
//		}
//		case TAYLOR_GREEN_VORTEX: {
//		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//		// Setup case/initial conditions: TAYLOR GREEN VORTEX
//		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//			PRECISION rho0 = 1.0;
//			PRECISION u0 = 1.0;
//			PRECISION p0 = 1.0;
//
//			for( size_t icell=0; icell < cells_soln.size(); icell++ ){
//				rho[icell] = rho0;
//
//				u[0][icell] =   u0 * cos(xc[icell][0] / L) * sin(xc[icell][1] / L);
//				u[1][icell] = - u0 * sin(xc[icell][0] / L) * cos(xc[icell][1] / L);
//
//				PRECISION u_mag = 0.0;
//				for( unsigned idim = 0; idim < DIM_CNT; idim++ ){
//					u_mag += u[idim][icell] * u[idim][icell];
//				}
//				u_mag = sqrt(u_mag);
//
//				PRECISION p = p0 - rho0 / 4.0 * (cos(2.0 * xc[icell][0] / L) + cos(2.0 * xc[icell][1] / L));
//
//				E[icell] = p / (rho[icell] * (gamma - ONE)) + 0.5 * u_mag * u_mag;
//			}
//			break;
//		}
//		case CYLINDER_TUNNEL: {
//		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//		// Setup case/initial conditions: Cylinder in Tunnel
//		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//			PRECISION rho0 = 1.0;
//			PRECISION p0 = 1.0;
//
//			for( size_t icell=0; icell < cells_soln.size(); icell++ ){
//				rho[icell] = rho0;
//
//				u[0][icell] = 6.0 * (0.41 / 2.0 - xc[icell][1]) * (0.41 / 2.0 + xc[icell][1] ) / ( 0.41 * 0.41 );
//				u[1][icell] = 0.0;
//
//				PRECISION u_mag = 0.0;
//				for( unsigned idim = 0; idim < DIM_CNT; idim++ ){
//					u_mag += u[idim][icell] * u[idim][icell];
//				}
//				u_mag = sqrt(u_mag);
//
//				PRECISION p = p0;
//
//				E[icell] = p / (rho[icell] * (gamma - ONE)) + 0.5 * u_mag * u_mag;
//			}
//			break;
//		}
//	}
//
//	for( size_t icell=0; icell < cells_soln.size(); icell++ ){
//		cells_soln[icell].vars.q_new[0] = rho[icell];
//		for( unsigned idim = 0; idim < DIM_CNT; idim++ ){
//			cells_soln[icell].vars.q_new[idim+1] = rho[icell] * u[idim][icell];
//		}
//		cells_soln[icell].vars.q_new[DIM_CNT+1] = rho[icell] * E[icell];
//	}
//}
//
////***************************************************************************************************
//// Set boundary conditions
////***************************************************************************************************
//template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
//void CFDv1_solver<PRECISION, DIM_CNT, FACE_CNT> :: set_boundary_conditions(){
//
//	// WALL
//	for ( size_t icell = 0; icell < ghost_bnd_soln[0].size(); icell++ ) {
//		ghost_bnd_soln[0][icell].q_new[0] = cells_soln[boundaries[0][icell][0]].q_new[0];
//		for ( size_t idim = 0; idim < DIM_CNT; idim++ ) {
//			ghost_bnd_soln[0][icell].q_new[idim+1] = - cells_soln[boundaries[0][icell][0]].q_new[idim+1];
//		}
//		ghost_bnd_soln[0][icell].q_new[DIM_CNT+1] = cells_soln[boundaries[0][icell][0]].q_new[DIM_CNT+1];
//
//		for ( size_t idim = 0; idim < DIM_CNT; idim++ ) {
//			ghost_bnd_conv[0][icell].rho_grad[idim] = cells_soln[boundaries[0][icell][0]].rho_grad[idim];
//			for ( size_t jdim = 0; jdim < DIM_CNT; jdim++ ) {
//				ghost_bnd_conv[0][icell].rhoU_grad[idim][jdim] = - cells_soln[boundaries[0][icell][0]].rhoU_grad[idim][jdim];
//				ghost_bnd_visc[0][icell].dudx[idim][jdim]      = - cells_soln[boundaries[0][icell][0]].dudx[idim][jdim];
//				ghost_bnd_visc[0][icell].tauMC[idim][jdim]     = - cells_soln[boundaries[0][icell][0]].tauMC[idim][jdim];
//			}
//			ghost_bnd_conv[0][icell].rhoE_grad[idim] = cells_soln[boundaries[0][icell][0]].rhoE_grad[idim];
//			ghost_bnd_conv[0][icell].Rpsi_grad[idim] = cells_soln[boundaries[0][icell][0]].Rpsi_grad[idim];
//			ghost_bnd_conv[0][icell].c_grad[idim]    = cells_soln[boundaries[0][icell][0]].c_grad[idim];
//
//			ghost_bnd_visc[0][icell].dTdx[idim]      = cells_soln[boundaries[0][icell][0]].dTdx[idim];
//		}
//	}
//
//	// INLET
//	for ( size_t icell = 0; icell < ghost_bnd_soln[1].size(); icell++ ) {
//
//		int gcell = boundaries[1][icell][0];
//
//		PRECISION rho = 1.0;
//		PRECISION Uvec[DIM_CNT];
//		Uvec[0] = 6.0 * (0.41 / 2.0 - xc[gcell][1]) * (0.41 / 2.0 + xc[gcell][1] ) / ( 0.41 * 0.41 );
//		Uvec[1] = 0.0;
//		PRECISION Umag_sqrt = 0.0;
//		for ( size_t idim = 0; idim < DIM_CNT; idim++ ) {
//			Umag_sqrt += Uvec[idim] * Uvec[idim];
//		}
//		PRECISION p = 1.0;
//		PRECISION E = p / (rho * (gamma - ONE)) + 0.5 * Umag_sqrt;
//
//		ghost_bnd_soln[1][icell].q_new[0] = rho;
//		for ( size_t idim = 0; idim < DIM_CNT; idim++ ) {
//			ghost_bnd_soln[1][icell].q_new[idim+1] = rho * Uvec[idim];
//		}
//		ghost_bnd_soln[1][icell].q_new[DIM_CNT+1] = rho * E;
//
//		for ( size_t idim = 0; idim < DIM_CNT; idim++ ) {
//			ghost_bnd_conv[1][icell].rho_grad[idim] = 0.0;
//			for ( size_t jdim = 0; jdim < DIM_CNT; jdim++ ) {
//				ghost_bnd_conv[1][icell].rhoU_grad[idim][jdim] = 0.0;
//				ghost_bnd_visc[1][icell].dudx[idim][jdim]      = 0.0;
//				ghost_bnd_visc[1][icell].tauMC[idim][jdim]     = 0.0;
//			}
//			ghost_bnd_conv[1][icell].rhoE_grad[idim] = 0.0;
//			ghost_bnd_conv[1][icell].Rpsi_grad[idim] = 0.0;
//			ghost_bnd_conv[1][icell].c_grad[idim]    = 0.0;
//
//			ghost_bnd_visc[1][icell].dTdx[idim]      = 0.0;
//		}
//	}
//
//	// OUTLET
//	for ( size_t icell = 0; icell < ghost_bnd_soln[2].size(); icell++ ) {
//		for ( size_t idim = 0; idim < DIM_CNT+2; idim++ ) {
//			ghost_bnd_soln[2][icell].q_new[idim] = cells_soln[boundaries[2][icell][0]].q_new[idim];
//		}
//
//		for ( size_t idim = 0; idim < DIM_CNT; idim++ ) {
//			ghost_bnd_conv[2][icell].rho_grad[idim] = 0.0;
//			for ( size_t jdim = 0; jdim < DIM_CNT; jdim++ ) {
//				ghost_bnd_conv[2][icell].rhoU_grad[idim][jdim] = 0.0;
//				ghost_bnd_visc[2][icell].dudx[idim][jdim]      = 0.0;
//				ghost_bnd_visc[2][icell].tauMC[idim][jdim]     = 0.0;
//			}
//			ghost_bnd_conv[2][icell].rhoE_grad[idim] = 0.0;
//			ghost_bnd_conv[2][icell].Rpsi_grad[idim] = 0.0;
//			ghost_bnd_conv[2][icell].c_grad[idim]    = 0.0;
//
//			ghost_bnd_visc[2][icell].dTdx[idim]      = 0.0;
//		}
//	}
//
//	// FAR-FIELD
//	for ( size_t icell = 0; icell < ghost_bnd_soln[3].size(); icell++ ) {
//
//		PRECISION Udirection[DIM_CNT] = { (PRECISION) cos(AoA * M_PI / 180.0), (PRECISION) sin(AoA * M_PI / 180.0) };
//		PRECISION T = 0.0;
//		for ( size_t idim = 0; idim < DIM_CNT; idim++ ) {
//			T += Udirection[idim] * cells_soln[boundaries[3][icell][0]].S[idim][boundaries[3][icell][1]];
//		}
//
//		if ( T >= 0.0 ) {		// INLET
//
//			PRECISION rho = 1.0;
//			PRECISION Uvec[DIM_CNT] = {1.0, 0.0};
//			PRECISION Umag_sqrt = 0.0;
//			for ( size_t idim = 0; idim < DIM_CNT; idim++ ) {
//				Umag_sqrt += Uvec[idim] * Uvec[idim];
//			}
//			PRECISION p = 1.0;
//			PRECISION E = p / (rho * (gamma - ONE)) + 0.5 * Umag_sqrt;
//
//			ghost_bnd_soln[3][icell].q_new[0] = rho;
//			for ( size_t idim = 0; idim < DIM_CNT; idim++ ) {
//				ghost_bnd_soln[3][icell].q_new[idim+1] = rho * Uvec[idim];
//			}
//			ghost_bnd_soln[3][icell].q_new[DIM_CNT+1] = rho * E;
//
//		} else {				 // OUTLET
//
//			for ( size_t idim = 0; idim < DIM_CNT+2; idim++ ) {
//				ghost_bnd_soln[3][icell].q_new[idim] = cells_soln[boundaries[3][icell][0]].q_new[idim];
//			}
//
//		}
//
//		for ( size_t idim = 0; idim < DIM_CNT; idim++ ) {
//			ghost_bnd_conv[3][icell].rho_grad[idim] = 0.0;
//			for ( size_t jdim = 0; jdim < DIM_CNT; jdim++ ) {
//				ghost_bnd_conv[3][icell].rhoU_grad[idim][jdim] = 0.0;
//				ghost_bnd_visc[3][icell].dudx[idim][jdim]      = 0.0;
//				ghost_bnd_visc[3][icell].tauMC[idim][jdim]     = 0.0;
//			}
//			ghost_bnd_conv[3][icell].rhoE_grad[idim] = 0.0;
//			ghost_bnd_conv[3][icell].Rpsi_grad[idim] = 0.0;
//			ghost_bnd_conv[3][icell].c_grad[idim]    = 0.0;
//
//			ghost_bnd_visc[3][icell].dTdx[idim]      = 0.0;
//		}
//	}
//}






























