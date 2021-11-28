#include <iomanip>

#include "api/mpi_env.h"
#include "mpi-ext.h"
#include <cstddef>


#include "api/templates_mpienv.h"

//***************************************************************************************************
void MPI_env::initialize( bool is_solver_sim ){

	MPI_Init( NULL, NULL );
	MPI_Comm_size( MPI_COMM_WORLD, &nProcs );
	MPI_Comm_rank( MPI_COMM_WORLD, &myrank );

	if( myrank == 0 ){
		master = true;
	}else{
		master = false;
	}

	if( is_solver_sim ){
		return;
	}

	// Get node master
	char local_proc_name[MPI_MAX_PROCESSOR_NAME];
	int   name_length;
	MPI_Get_processor_name( local_proc_name, &name_length );

	std::string name( local_proc_name );

	char *global_proc_name = new char[ this->size() * MPI_MAX_PROCESSOR_NAME ];
	MPI_Allgather( &local_proc_name, name_length, MPI_CHAR, &global_proc_name[0], name_length, MPI_CHAR, MPI_COMM_WORLD );

	std::vector<std::string> all_proc_names( this->size(), "" );

	for( int irank=0; irank < this->size(); irank++ ){
		for( int ichar=0; ichar < name_length; ichar++ ){
			all_proc_names[irank] += global_proc_name[ichar + name_length*irank];
		}
	}
	int node_master_rank = this->size();
	std::vector<bool> is_node_local( this->size(), false );

	for( int irank=0; irank < this->size(); irank++ ){
		if( all_proc_names[this->rank()].compare( all_proc_names[irank] ) == 0 ){

			// Add node-local MPI ranks
			node_local_ranks.push_back( irank );
			is_node_local[irank] = true;

			// Identify node master
			if( node_master_rank > irank ){
				node_master_rank = irank;
			}
			node_master_rank = node_master_rank > irank ? irank : node_master_rank;
		}
	}

	if( myrank == node_master_rank ){
		node_master = true;
	}else{
		node_master = false;
	}

	// Prepare node-local communicator
	MPI_Comm node_comm;
	MPI_Comm_split( MPI_COMM_WORLD, node_master_rank, this->rank(), &node_comm );

	int node_size, node_rank;
	MPI_Comm_size( node_comm, &node_size );
	MPI_Comm_rank( node_comm, &node_rank );

	this->node_comm.set_comm( node_comm );
	this->node_comm.set_size( node_size );
	this->node_comm.set_rank( node_rank );
}

//***************************************************************************************************
void MPI_env::finalize(){
	MPI_Finalize();
}

//***************************************************************************************************
template <typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT>
void MPI_env :: set_halo_comm_type( enum t_mpi_halo_comm type_mpi_comm,
									fm_vector<CFDv0_cell<PRECISION, DIM_CNT, FACE_CNT>> cells_cfd,
									fm_vector<fm_vector<CFDv0_cell<PRECISION, DIM_CNT, FACE_CNT>>> ghost_cells ){

	halo_comm_type = type_mpi_comm;

	// Neighborhood collectives
	switch( halo_comm_type ){
		case MPI_NEIGHCOLL_BLCK:
		case MPI_NEIGHCOLL_NONB:

			// Init neighborhood graph communication
			init_graph( mpi_neighbors );

			// Hpath cells
			sendcounts = new int          [ mpi_neighbors.size() ];
			senddisps  = new MPI_Aint     [ mpi_neighbors.size() ];
			sendtypes  = new MPI_Datatype [ mpi_neighbors.size() ];

			recvcounts = new int          [ mpi_neighbors.size() ];
			recvdisps  = new MPI_Aint     [ mpi_neighbors.size() ];
			recvtypes  = new MPI_Datatype [ mpi_neighbors.size() ];

			for( size_t ineigh=0; ineigh < mpi_neighbors.size(); ineigh++ ){
				sendcounts[ineigh] = (int) cells_cfd.size();
				senddisps [ineigh] = 0;
				sendtypes [ineigh] = type_cfdcell;

				recvcounts[ineigh] = ghost_cells[ineigh].size();
				recvdisps [ineigh] = (char*) &(ghost_cells[ineigh][0]) - (char*) &(ghost_cells[0][0]);
				recvtypes [ineigh] = type_cfdcell;
			}

			break;

		case MPI_NEIGHCOLL_PERS:

			// Init neighborhood graph communication
			init_graph( mpi_neighbors );

			// Hpath cells
			sendcounts = new int          [ mpi_neighbors.size() ];
			senddisps  = new MPI_Aint     [ mpi_neighbors.size() ];
			sendtypes  = new MPI_Datatype [ mpi_neighbors.size() ];

			recvcounts = new int          [ mpi_neighbors.size() ];
			recvdisps  = new MPI_Aint     [ mpi_neighbors.size() ];
			recvtypes  = new MPI_Datatype [ mpi_neighbors.size() ];

			for( size_t ineigh=0; ineigh < mpi_neighbors.size(); ineigh++ ){
				sendcounts[ineigh] = (int) cells_cfd.size();
				senddisps [ineigh] = 0;
				sendtypes [ineigh] = type_cfdcell;

				recvcounts[ineigh] = ghost_cells[ineigh].size();
				recvdisps [ineigh] = (char*) &(ghost_cells[ineigh][0]) - (char*) &(ghost_cells[0][0]);
				recvtypes [ineigh] = type_cfdcell;
			}

			MPIX_Neighbor_alltoallw_init( &cells_cfd  [0],    sendcounts, senddisps, sendtypes,
										  &ghost_cells[0][0], recvcounts, recvdisps, recvtypes,
										  graph_comm(), MPI_INFO_NULL, &neighcoll_req );
			MPI_Start( &neighcoll_req   );
		default:
			break;
	}

}

//***************************************************************************************************
// MPI_Wait/Waitall
// 		- can be used with any non-blocking and persistent calls
// 		- when blocking calls are used, simply returns without doing anything
//***************************************************************************************************
void MPI_env:: halo_comm_wait(){

	// Serial simulations
	if( this->is_serial() ) return;

	// Depending on communication type...
	int num_mpi_neighbors = (int) comm_reqs.size() ;

	switch( halo_comm_type ){

		case MPI_TWOSIDED_BLCK:
		case MPI_ONESIDED_BLCK:
		case MPI_NEIGHCOLL_BLCK:
			return;

		case MPI_TWOSIDED_NONB:
		case MPI_TWOSIDED_PERS:
		case MPI_ONESIDED_NONB:
		case MPI_ONESIDED_PERS:
			MPI_Waitall( num_mpi_neighbors, &comm_reqs[0], MPI_STATUS_IGNORE );
			break;

		case MPI_NEIGHCOLL_NONB:
			MPI_Waitall( 1, &comm_reqs[0], MPI_STATUS_IGNORE );
			break;
		case MPI_NEIGHCOLL_PERS:
			MPI_Waitall( 1, &neighcoll_req, MPI_STATUS_IGNORE );
			break;
		case MPI_NONE:
			cout << "Error with halo communication type" << endl;
			MPI_Abort( MPI_COMM_WORLD, 928 );
	}
}

//***************************************************************************************************
// MPI_Wait/Waitall
// 		- can be used with any non-blocking and persistent calls
// 		- when blocking calls are used, simply returns without doing anything
//***************************************************************************************************
void MPI_env:: halo_comm_test( int &flag ){

	// Serial simulations
	if( this->is_serial() ) return;

	// Depending on communication type...
	switch( halo_comm_type ){

		case MPI_TWOSIDED_BLCK:
		case MPI_ONESIDED_BLCK:
		case MPI_NEIGHCOLL_BLCK:
			return;

		case MPI_TWOSIDED_NONB:
		case MPI_TWOSIDED_PERS:
		case MPI_ONESIDED_NONB:
		case MPI_ONESIDED_PERS:

			MPI_Testall( comm_reqs.size(), &comm_reqs[0], &flag, MPI_STATUS_IGNORE );
			break;

		case MPI_NEIGHCOLL_NONB:
			MPI_Waitall( 1, &comm_reqs[0], MPI_STATUS_IGNORE );
			break;
		case MPI_NEIGHCOLL_PERS:
			MPI_Waitall( 1, &neighcoll_req, MPI_STATUS_IGNORE );
			break;
		case MPI_NONE:
			cout << "Error with halo communication type" << endl;
			MPI_Abort( MPI_COMM_WORLD, 928 );
	}
}

//***************************************************************************************************
void MPI_env:: init_graph( const std::vector<int> mpi_neighbors ){

	// Initialize graph topology
	int  src_size  = (int) mpi_neighbors.size();			// Number of neighbors
	int* src_ids   = new int[ mpi_neighbors.size() ];
	int* src_wghts = new int[ mpi_neighbors.size() ];

	int  des_size  = (int) mpi_neighbors.size();			// Number of neighbors
	int* des_ids   = new int[ mpi_neighbors.size() ];
	int* des_wghts = new int[ mpi_neighbors.size() ];

	bool allow_reorder = true;

	for( size_t ineigh=0; ineigh < mpi_neighbors.size(); ineigh++ ){
		src_ids[ineigh] = mpi_neighbors[ineigh];
		des_ids[ineigh] = mpi_neighbors[ineigh];

		src_wghts[ineigh] = 1;
		des_wghts[ineigh] = 1;
	}

	MPI_Dist_graph_create_adjacent( MPI_COMM_WORLD,
									src_size, src_ids, src_wghts,
									des_size, des_ids, des_wghts,
									MPI_INFO_NULL, allow_reorder, &comm_graph );

	int* test_src_ids   = new int[ mpi_neighbors.size() ];
	int* test_src_wghts = new int[ mpi_neighbors.size() ];
	int* test_des_ids   = new int[ mpi_neighbors.size() ];
	int* test_des_wghts = new int[ mpi_neighbors.size() ];

	MPI_Dist_graph_neighbors( comm_graph, src_size, test_src_ids, test_src_wghts,
										  des_size, test_des_ids, test_des_wghts );

	delete[] src_ids;
	delete[] src_wghts;
	delete[] des_ids;
	delete[] des_wghts;

	delete[] test_src_ids;
	delete[] test_src_wghts;
	delete[] test_des_ids;
	delete[] test_des_wghts;
}

//***************************************************************************************************
//***************************************************************************************************
// 									CFD cells
//***************************************************************************************************
//***************************************************************************************************

// Isendrecv
template <typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT>
void MPI_env:: isendrecv(
						  fm_vector<CFDv0_cell<PRECISION, DIM_CNT, FACE_CNT>>  &cells_cfd,
				fm_vector<fm_vector<CFDv0_cell<PRECISION, DIM_CNT, FACE_CNT>>> &ghost_cells ){

	const int million = 1000;

	for( size_t ineigh=0; ineigh < mpi_neighbors.size(); ineigh++ ){

		int src  = mpi_neighbors[ineigh];
		int dest = mpi_neighbors[ineigh];

		int send_tag = million*rank() + mpi_neighbors[ineigh];
		int recv_tag = rank() + million*mpi_neighbors[ineigh];

		int local_size = (int) cells_cfd.size();
		int neigh_size = (int) ghost_cells[ineigh].size();

		MPI_Isend( &cells_cfd[0]          , local_size, type_cfdcell, dest, send_tag, MPI_COMM_WORLD, &comm_reqs[ineigh] );
		MPI_Irecv( &ghost_cells[ineigh][0], neigh_size, type_cfdcell, src , recv_tag, MPI_COMM_WORLD, &comm_reqs[ineigh+mpi_neighbors.size()] );
	}
}

// Neighbor_alltoallw
template <typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT>
void MPI_env:: neighbor_alltoallw( CFDv0_cell<PRECISION, DIM_CNT, FACE_CNT> *cells_cfd,
								   CFDv0_cell<PRECISION, DIM_CNT, FACE_CNT> *ghost_cells ){

	MPI_Neighbor_alltoallw( cells_cfd,   sendcounts, senddisps, sendtypes,
							ghost_cells, recvcounts, recvdisps, recvtypes, graph_comm() );
}

// Ineighbor_alltoallw
template <typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT>
void MPI_env:: ineighbor_alltoallw( CFDv0_cell<PRECISION, DIM_CNT, FACE_CNT> *cells_cfd,
								    CFDv0_cell<PRECISION, DIM_CNT, FACE_CNT> *ghost_cells ){

	MPI_Ineighbor_alltoallw( cells_cfd,   sendcounts, senddisps, sendtypes,
							 ghost_cells, recvcounts, recvdisps, recvtypes,
							 graph_comm(), &comm_reqs[0] );
}









































