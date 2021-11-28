#ifndef MPI_ENV_H
#define MPI_ENV_H

#include "mpi.h"
#include "elements.h"
#include "mpi_env.h"
#include "cfdv0_elements.h"

struct t_mpi_neighbors{
	int id_local;
	int id_global;
};

enum t_mpi_halo_comm {
	MPI_TWOSIDED_BLCK  = 0, // Plain cell
	MPI_TWOSIDED_NONB  = 1, // Neighbor is in a neighboring submesh
	MPI_TWOSIDED_PERS  = 2, // Neighbor is in a neighboring submesh
	MPI_ONESIDED_BLCK  = 3, // Plain cell
	MPI_ONESIDED_NONB  = 4, // Neighbor is in a neighboring submesh
	MPI_ONESIDED_PERS  = 5, // Neighbor is in a neighboring submesh
	MPI_NEIGHCOLL_BLCK = 6,
	MPI_NEIGHCOLL_NONB = 7,
	MPI_NEIGHCOLL_PERS = 8,
	MPI_NONE           = 9, // Error
};

class t_mpi_comm{

	public:
		// Constructors
		t_mpi_comm(){}
		t_mpi_comm( MPI_Comm tmp_comm, int tmp_size, int tmp_rank ){
			m_comm = tmp_comm;
			m_size = tmp_size;
			m_rank = tmp_rank;
		}

		// Communicator
		void set_comm( MPI_Comm tmp_comm ){
			m_comm = tmp_comm;
		}
		MPI_Comm self(){
			return m_comm;
		}

		// Size
		void set_size( int tmp_size ){
			m_size = tmp_size;
		}
		int size(){
			return m_size;
		}

		// Rank
		void set_rank( int tmp_rank ){
			m_rank = tmp_rank;
		}
		int rank(){
			return m_rank;
		}

		// Master
		int master(){
			return 0;
		}

	private:
		MPI_Comm m_comm;
		int m_size;
		int m_rank;
};

//***************************************************************************************************
// MPI environment class
// 		- handles all calls to MPI functions/routines
// 		- setup and wrapping of MPI environment
// 		- MPI datatypes, buffer communication, collectives...
//***************************************************************************************************
class MPI_env {

	public:

		// Initialize MPI environment
		void initialize( bool is_solver_sim );
		void finalize();

		//MPI_Datatype get_regular_mesh_MPI_datatype( t_mesh *mesh );
		//void test_mpi_exchange( t_mesh *mesh, MPI_Datatype type_node );

		int size(){
			return nProcs;
		}

		int rank(){
			return myrank;
		}

		int get_master(){
			return 0;
		}

		bool is_master(){
			return master;
		}

		int get_node_master(){
			return node_master_rank;
		}

		bool is_node_master(){
			return node_master;
		}

		bool is_serial(){
			return (nProcs == 1);
		}

		void set_cell_type( MPI_Datatype tmp_cell_type ){ type_cell = tmp_cell_type; };
		MPI_Datatype cell_type()      { return type_cell;       }
		MPI_Comm     graph_comm()     { return comm_graph;      }

		void barrier(void)          { MPI_Barrier( MPI_COMM_WORLD ); }
		void abort( int error_code ){ MPI_Abort  ( MPI_COMM_WORLD, error_code ); }
		void broadcast_partitions( std::vector<int> partitions );

		// Initialize MPI graph topology
		void init_graph( std::vector<int> mpi_neighbors );

		template <typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT>
		void set_halo_comm_type( enum t_mpi_halo_comm type_mpi_comm,
									fm_vector<CFDv0_cell<PRECISION, DIM_CNT, FACE_CNT>> cells_cfd,
									fm_vector<fm_vector<CFDv0_cell<PRECISION, DIM_CNT, FACE_CNT>>> ghost_cells );

		enum t_mpi_halo_comm get_halo_comm_type(){ return halo_comm_type; }

		void set_mpi_neighbors( std::vector<int> m_mpi_neighbors ){
			mpi_neighbors = m_mpi_neighbors;
		}

		void print_from_master( std::string print_statement ){
			barrier();
			if( is_master() ) cout << print_statement << endl;
		}

		void set_mesh_size_type( MPI_Datatype type_mesh_size ){ type_meshsize = type_mesh_size; }
		MPI_Datatype get_type_meshsize(){ return type_meshsize; }

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Wait functions
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		void halo_comm_wait();
		void halo_comm_test( int &flag );

		// Neighborhood collective
		void save_neighborhood_request( MPI_Request request_hpath, MPI_Request request_deadend ){
			neighcoll_req_hpath   = request_hpath   ;
			neighcoll_req_deadend = request_deadend ;
		}

		void neighbor_alltoallw_start  (){
			MPI_Start( &neighcoll_req );
		}

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// CFD Solver
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		// Two-sided
		template <typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT>
		void isendrecv(          fm_vector<CFDv0_cell<PRECISION, DIM_CNT, FACE_CNT>>  &cells_hpath,
					   fm_vector<fm_vector<CFDv0_cell<PRECISION, DIM_CNT, FACE_CNT>>> &ghost_hpath );

		// Neighborhood collective
		template <typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT>
		void  neighbor_alltoallw( CFDv0_cell<PRECISION, DIM_CNT, FACE_CNT> *cells_cfd,
								  CFDv0_cell<PRECISION, DIM_CNT, FACE_CNT> *ghost_cells );

		template <typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT>
		void ineighbor_alltoallw( CFDv0_cell<PRECISION, DIM_CNT, FACE_CNT> *cells_cfd,
								  CFDv0_cell<PRECISION, DIM_CNT, FACE_CNT> *ghost_cells );

		// New functions
		inline void set_type_cfdcell( MPI_Datatype t_cfdcell ){ type_cfdcell = t_cfdcell; }
		inline void set_halo_comm_reqs( int n ){ comm_reqs.resize( n ); };

		std::vector<int> mpi_neighbors;
		std::vector<MPI_Request> comm_reqs;
		std::vector<int> node_local_ranks;

		// Node communicator
		class t_mpi_comm node_comm;
	private:
		int  nProcs;
		int  myrank;
		bool master;
		int  node_master_rank;
		bool node_master;

		MPI_Datatype type_cell;
		MPI_Comm comm_graph;
		MPI_Request neighcoll_req_hpath;
		MPI_Request neighcoll_req_deadend;

		MPI_Datatype type_cfdcell;

		std::vector<MPI_Request> twosided_send_req_hpath, twosided_send_req_reg;
		std::vector<MPI_Request> twosided_recv_req_hpath, twosided_recv_req_reg;

		int          *sendcounts, *recvcounts;
		MPI_Aint     *senddisps , *recvdisps ;
		MPI_Datatype *sendtypes , *recvtypes ;

		int          *sendcounts_reg, *recvcounts_reg;
		MPI_Aint     *senddisps_reg , *recvdisps_reg ;
		MPI_Datatype *sendtypes_reg , *recvtypes_reg ;

		MPI_Request neighcoll_req;

		enum t_mpi_halo_comm halo_comm_type = MPI_NONE;

		MPI_Datatype type_meshsize;
};

#endif








































