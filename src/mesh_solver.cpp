#include <iostream>
#include <iomanip>
#include <fstream>
#include <ctime>
#include <sstream>
#include <array>
#include "api/timing.h"

#include "api/fastmesh.h"
#include "api/mesh_writer.h"
#include "api/mpi_env.h"

using namespace std;
using namespace fastmesh;

template void Mesh<double, 2> :: initialize_solver( fastmesh_solvers solver );
template void Mesh<double, 2> :: output_solution( double time, double PRINT_FILE_FREQUENCY, std::string file_base );
template void Mesh<double, 2> :: solve();

//***************************************************************************************************
// Initialize sovler
// 	- for now, only jacobi solver
// 	- next: extend for multiple solvers w/ function pointers
//***************************************************************************************************
template <typename PRECISION, unsigned DIM_CNT >
void Mesh<PRECISION, DIM_CNT> :: initialize_solver( fastmesh_solvers solver ){

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Initialize MPI environment
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	int m_count_3 = 0;
	int m_count_4 = 0;

	std::vector<std::vector<t_solution_vars<PRECISION,DIM_CNT>*>> cells_ptr( gmsh_parts.size() );
	std::vector<std::vector<t_solution_vars<PRECISION,DIM_CNT>*>> ghost_ptr;
	std::vector<std::vector<t_solution_vars<PRECISION,DIM_CNT>*>> ghost_bnd_ptr;

	std::vector<CFDv0_solver<PRECISION, DIM_CNT, 3>> tmp_solver_3( gmsh_parts.size() );
	std::vector<CFDv0_solver<PRECISION, DIM_CNT, 4>> tmp_solver_4( gmsh_parts.size() );

	for( size_t ipart=0; ipart < gmsh_parts.size(); ipart++ ){
		switch( gmsh_parts[ipart].faces_per_cell ){
			case 3:
				tmp_solver_3[m_count_3].sm_id = m_count_3;
				tmp_solver_3[m_count_3].allocate( &gmsh_parts[ipart], cells_ptr[ipart] );
				tmp_solver_3[m_count_3].allocate_ghost_cells( mpi_env, &gmsh_parts[ipart], ghost_ptr, ghost_bnd_ptr );

				m_count_3++;
				break;
			case 4:
				tmp_solver_4[m_count_4].sm_id = m_count_4;
				tmp_solver_4[m_count_4].allocate( &gmsh_parts[ipart], cells_ptr[ipart] );
				tmp_solver_4[m_count_4].allocate_ghost_cells( mpi_env, &gmsh_parts[ipart], ghost_ptr, ghost_bnd_ptr );

				m_count_4++;
				break;
		}
	}

	m_count_3 = 0;
	m_count_4 = 0;
	for( size_t ipart=0; ipart < gmsh_parts.size(); ipart++ ){
		switch( gmsh_parts[ipart].faces_per_cell ){
			case 3:
				tmp_solver_3[m_count_3].assign_pointers( &gmsh_parts[ipart], cells_ptr, ghost_ptr, ghost_bnd_ptr );
				tmp_solver_3[m_count_3].initialize( mpi_env, &gmsh_parts[ipart] );
				submesh_cfd.push_back( tmp_solver_3[m_count_3] );

				m_count_3++;
				break;
			case 4:
				tmp_solver_4[m_count_4].assign_pointers( &gmsh_parts[ipart], cells_ptr, ghost_ptr, ghost_bnd_ptr );
				tmp_solver_4[m_count_4].initialize( mpi_env, &gmsh_parts[ipart] );
				submesh_cfd.push_back( tmp_solver_4[m_count_4] );

				m_count_4++;
				break;
		}
	}

	std::string file_base = "output/";
	output_solution( 0.0, 1.0, file_base );
}

//***************************************************************************************************
//***************************************************************************************************
template < typename PRECISION, unsigned DIM_CNT >
void Mesh<PRECISION, DIM_CNT> :: solve(){

	if( mpi_env.is_master() ){
		cout << endl;
		cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
		cout << "CFD solver: caafoam" << endl;
		cout << "    - Flow type    : Euler equation" << endl;
		cout << "    - MPI framework: ";
		switch( mpi_env.get_halo_comm_type() ){
			case MPI_TWOSIDED_BLCK:
				cout << "two-sided, blocking" << endl;
				break;
			case MPI_TWOSIDED_NONB:
				cout << "two-sided, non-blocking" << endl;
				break;
			case MPI_TWOSIDED_PERS:
				cout << "two-sided, persistent" << endl;
				break;
			case MPI_NEIGHCOLL_BLCK:
				cout << "neighborhood collective, blocking" << endl;
				break;
			case MPI_NEIGHCOLL_NONB:
				cout << "neighborhood collective, non-blocking" << endl;
				break;
			case MPI_NEIGHCOLL_PERS:
				cout << "neighborhood collective, persistent" << endl;
				break;
			default:
				cout << "None: serial simulation" << endl;
				break;
		}
		cout << "    - MPI cores    :" << setw(10) << mpi_env.size() << endl;
		cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Initialization phase
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	// CPU timing
	uint64_t t0, t1, t2, t3, t5, t6, t10, t11;
	double *t_submesh = new double [submesh_cfd.size()];
	double t_init    = 0.;
	double t_solver  = 0.;
	double t_output  = 0.;
	double t_visc    = 0.;
	double t_ghost   = 0.;
	double t_wait    = 0.;
	double t_mpi     = 0.;
	double t_grads   = 0.;
	double t_wait0   = 0.;
	double t_int0    = 0.;
	double t_int1    = 0.;
	double t_bnd     = 0.;
	double t_int     = 0.;
	double t_barrier = 0.;

	cpu_timing timing;
	timing.calibrate_rdtsc();

	for( size_t ipart=0; ipart < submesh_cfd.size(); ipart++ ){
		t_submesh[ipart] = 0.;
	}

	bool MINMOD_INTERPOLATION_EXISTS = true;
	bool TIME_STEP_CONST = false;

	// Local variables
	int time_step     = 0;

	PRECISION dt      = 0.0001;
	PRECISION t_start = 0.;
	PRECISION t_end   = 60.;
	PRECISION time    = t_start;
	PRECISION cflMax  = 1.0;

	const int PRINT_INFO_FREQUENCY	= 1;
	PRECISION PRINT_FILE_FREQUENCY	= 0.025;
	PRECISION SAVE_FILE_TIME		= t_start + PRINT_FILE_FREQUENCY;	// If time step not constant
	const int SAVE_FILE_STEP		= PRINT_FILE_FREQUENCY / dt;		// If time step constant

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// CFD solver algorithm. For each time-step:
	// 		- prepapre solver for time-step (e.g. q_old = q_new, etc...)
	// 		- solve for boundary submesh
	// 		- start non-blocking MPI exchange of ghost cells
	// 		- solve for interior submeshes
	// 		- wait for MPI comm to be done
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	if ( t_start > 0. ) {
		int submesh = 0;
		for( auto all_sm = submesh_cfd.begin(); all_sm < submesh_cfd.end(); all_sm++ ){
			all_sm->restart( t_start, PRINT_FILE_FREQUENCY, submesh );
			submesh++;
		}
		if ( mpi_env.is_master() ){
			cout << endl;
			cout << "Simulation restarted successfully!" << endl;
			cout << endl;
		}
	}

	auto bnd_sm = submesh_cfd.begin() + INDEX_BND_SUBMESH;
	MPI_Barrier( MPI_COMM_WORLD );
	bnd_sm->exchange_ghost_cells( mpi_env );
	bnd_sm->set_boundary_conditions();
	mpi_env.halo_comm_wait();
	MPI_Barrier( MPI_COMM_WORLD );

	t0 = clock();

	while( time < t_end ){

		// For each RK-step...
		t5 = clock();
		for( int rk_step=0; rk_step < num_RK_steps; rk_step++ ){

			//cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
			//cout << "rk step =" << setw(10) << rk_step << endl;
			//cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;

			// Prepare for next iteration
			t2 = clock();
			for( auto all_sm = submesh_cfd.begin(); all_sm < submesh_cfd.end(); all_sm++ ){
				all_sm->prepare_for_timestep( rk_step );
			}
			t3 = clock(); t_init += (double) (t3-t2) / CLOCKS_PER_SEC;

			//t10 = clock();
			//mpi_env.barrier();
			//t11 = clock(); t_barrier += (double) (t11-t10) / CLOCKS_PER_SEC;

			// Wait for MPI comm of last interiation to finish
			t10 = clock();
			bnd_sm->halo_comm_wait( mpi_env );
			t11 = clock(); t_wait  += (double) (t11-t10) / CLOCKS_PER_SEC;

			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			// Boundary submesh: gradients & viscous first step
			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			if( MINMOD_INTERPOLATION_EXISTS ){
				t10 = clock();
				bnd_sm->calc_gradients( mpi_env );
				t11 = clock(); t_grads += (double) (t11-t10) / CLOCKS_PER_SEC;
				t_bnd += (double) (t11-t10) / CLOCKS_PER_SEC;
			}
			t10 = clock();
			bnd_sm->calc_VIS( mpi_env );
			t11 = clock(); t_visc += (double) (t11-t10) / CLOCKS_PER_SEC;
			t_bnd += (double) (t11-t10) / CLOCKS_PER_SEC;

			t10 = clock();
			bnd_sm->set_boundary_conditions();
			bnd_sm->exchange_ghost_cells( mpi_env );
			t11 = clock(); t_ghost += (double) (t11-t10) / CLOCKS_PER_SEC;

			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			// Interior submesh: gradients & viscous first step
			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			if( MINMOD_INTERPOLATION_EXISTS ){
				t10 = clock();
				for( auto int_sm = submesh_cfd.begin()+1; int_sm < submesh_cfd.end(); int_sm++ ){
					int_sm->calc_gradients( mpi_env );
				}
				t11 = clock(); t_grads += (double) (t11-t10) / CLOCKS_PER_SEC;
				t_int0 += (double) (t11-t10) / CLOCKS_PER_SEC;
				t_int  += (double) (t11-t10) / CLOCKS_PER_SEC;
			}
			t10 = clock();
			for( auto int_sm = submesh_cfd.begin()+1; int_sm < submesh_cfd.end(); int_sm++ ){
				int_sm->calc_VIS( mpi_env );
			}
			t11 = clock(); t_visc += (double) (t11-t10) / CLOCKS_PER_SEC;
			t_int0 += (double) (t11-t10) / CLOCKS_PER_SEC;
			t_int  += (double) (t11-t10) / CLOCKS_PER_SEC;

			//t10 = clock();
			//mpi_env.barrier();
			//t11 = clock(); t_barrier += (double) (t11-t10) / CLOCKS_PER_SEC;

			t10 = clock();
			bnd_sm->halo_comm_wait( mpi_env );
			t11 = clock(); t_wait0 += (double) (t11-t10) / CLOCKS_PER_SEC;

			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			// Boundary submesh: advance time-step (convection + viscous, step 2)
			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			t10 = clock();
			bnd_sm->one_rk_step_v1( rk_step, dt, mpi_env );
			t11 = clock(); t_solver += (double) (t11-t10) / CLOCKS_PER_SEC;

			t10 = clock();
			bnd_sm->exchange_ghost_cells( mpi_env );
			t11 = clock(); t_ghost += (double) (t11-t10) / CLOCKS_PER_SEC;

			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			// Interior submesh: advance time-step (convection + viscous, step 2)
			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			t10 = clock();
			for( auto int_sm = submesh_cfd.begin()+1; int_sm < submesh_cfd.end(); int_sm++ ){
				int_sm->one_rk_step_v1( rk_step, dt, mpi_env );
			}
			t11 = clock(); t_solver += (double) (t11-t10) / CLOCKS_PER_SEC;
			t_int1 += (double) (t11-t10) / CLOCKS_PER_SEC;
			t_int  += (double) (t11-t10) / CLOCKS_PER_SEC;

			// Wait for MPI comm of last interiation to finish
			t10 = clock();
			bnd_sm->halo_comm_wait( mpi_env );
			t11 = clock(); t_wait  += (double) (t11-t10) / CLOCKS_PER_SEC;
		}
		t6 = clock();

		// CFL
		if ( !TIME_STEP_CONST ) {
			PRECISION dt_local = 100.0;
			for( auto it_sm = submesh_cfd.begin(); it_sm < submesh_cfd.end(); it_sm++ ){

				PRECISION tmp_dt = it_sm->compute_dt( cflMax );
				if( dt_local > tmp_dt )
					dt_local = tmp_dt;
			}
			MPI_Datatype mpi_precision;
			if( sizeof(PRECISION) == sizeof(double) )
				mpi_precision = MPI_DOUBLE;
			else
				mpi_precision = MPI_FLOAT;

			MPI_Allreduce( &dt_local, &dt, 1, mpi_precision, MPI_MAX, MPI_COMM_WORLD );

			if ( dt == 100.0 )
				break;

			if ( time+dt > SAVE_FILE_TIME)
				dt = SAVE_FILE_TIME - time;
		}

		time += dt;
		time_step++;

		// Sim info
		if( time_step % PRINT_INFO_FREQUENCY == 0 ){
			std::vector<struct t_domain_info<PRECISION, DIM_CNT>> domain_info( submesh_cfd.size() );

			if( mpi_env.is_master() )
				cout << "timestep = " << setw(10) << "Time step: " << time_step
									<< setw(15) << "Time: " << time
									<< setw(15) << "dt: " << dt
									<< setw(15) << (double) (t6-t5) /  CLOCKS_PER_SEC
									<< endl;
			//int ind_sm = 0;
			//for( auto it_sm = submesh_cfd.begin(); it_sm < submesh_cfd.end(); it_sm++ ){
			//	it_sm->gather_info( mpi_env, domain_info[ind_sm++] );
			//}
		}

		// Output solution
		if ( TIME_STEP_CONST ) {
			if ( time_step % SAVE_FILE_STEP == 0 ) {
				t10 = clock();
				std::string file_base = "output/";

				output_solution( time, PRINT_FILE_FREQUENCY, file_base );
				t11 = clock(); t_output += (double) (t11-t10) / CLOCKS_PER_SEC;
			}
		} else {
			if ( time == SAVE_FILE_TIME ) {
				t10 = clock();
				std::string file_base = "output/";

				output_solution( time, PRINT_FILE_FREQUENCY, file_base );
				t11 = clock(); t_output += (double) (t11-t10) / CLOCKS_PER_SEC;

				SAVE_FILE_TIME += PRINT_FILE_FREQUENCY;
			}
		}
	}
	t1 = clock();

	double t_total = (double) (t1 - t0) / CLOCKS_PER_SEC;
	double tmp_compute;
	double *t_compute = new double[mpi_env.size()];
	double global_cpu = -100;
	double global_ghost, global_wait, global_output, global_visc, global_init;
	double global_grads, global_solver, global_compute;
	double global_wait0, global_int0, global_int1, global_int, global_bnd;
	for( size_t ism=0; ism < submesh_cfd.size(); ism++ ){
		t_solver += t_submesh[ism];
	}

	MPI_Allreduce( &t_total  , &global_cpu    , 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
	MPI_Allreduce( &t_ghost  , &global_ghost  , 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
	MPI_Allreduce( &t_wait   , &global_wait   , 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
	MPI_Allreduce( &t_output , &global_output , 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
	MPI_Allreduce( &t_visc   , &global_visc   , 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
	MPI_Allreduce( &t_init   , &global_init   , 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
	MPI_Allreduce( &t_grads  , &global_grads  , 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
	MPI_Allreduce( &t_solver , &global_solver , 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

	MPI_Allreduce( &t_wait0  , &global_wait0  , 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
	MPI_Allreduce( &t_int0   , &global_int0   , 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
	MPI_Allreduce( &t_int1   , &global_int1   , 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
	MPI_Allreduce( &t_int    , &global_int    , 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
	MPI_Allreduce( &t_bnd    , &global_bnd    , 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

	double min_barrier, max_barrier, avg_barrier;
	MPI_Allreduce( &t_barrier, &min_barrier, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD );
	MPI_Allreduce( &t_barrier, &max_barrier, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
	MPI_Allreduce( &t_barrier, &avg_barrier, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

	tmp_compute = t_init + t_grads + t_visc + t_grads;
	global_compute = global_solver + global_init + global_visc + global_grads;
	t_mpi = global_ghost + global_wait + global_wait0;

	MPI_Gather( &tmp_compute, 1, MPI_DOUBLE, &t_compute[0], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );

	if( mpi_env.is_master() ){
		if ( dt < 100.0 ){
			//double global_ideal0 = global_wait0 < global_int0 ? 0.0 : global_wait0 - global_int0;
			//double global_ideal1 = global_wait  < global_int1 ? 0.0 : global_wait  - global_int1;
			cout << endl;
			cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
			cout << "Compute =" << setw(15) << std::setprecision(5) << global_compute    / mpi_env.size() << endl;
			cout << "    - solver:  " << setw(15) << std::setprecision(5) << global_solver / mpi_env.size() << endl;
			cout << "    - init:    " << setw(15) << std::setprecision(5) << global_init   / mpi_env.size() << endl;
			cout << "    - visc:    " << setw(15) << std::setprecision(5) << global_visc   / mpi_env.size() << endl;
			cout << "    - grad:    " << setw(15) << std::setprecision(5) << global_grads  / mpi_env.size() << endl;
			cout << "    - interior:" << setw(15) << std::setprecision(5) << global_int    / mpi_env.size() << endl;
			cout << "    - boundary:" << setw(15) << std::setprecision(5) << global_bnd    / mpi_env.size() << endl;
			cout << "MPI comm =" << setw(15) << std::setprecision(5) << t_mpi / mpi_env.size() << endl;
			cout << "    - ghost   :" << setw(15) << std::setprecision(5) << global_ghost  / mpi_env.size() << endl;
			//cout << "    - wait0   :" << setw(15) << std::setprecision(5) << global_wait0  / mpi_env.size() << endl;
			//cout << "    - overlap0:" << setw(15) << std::setprecision(5) << global_int0   / mpi_env.size() << endl;
			//cout << "    - ideal0  :" << setw(15) << std::setprecision(5) << global_ideal0 / mpi_env.size() << endl;
			//cout << "    - wait1   :" << setw(15) << std::setprecision(5) << global_wait   / mpi_env.size() << endl;
			//cout << "    - overlap1:" << setw(15) << std::setprecision(5) << global_int1   / mpi_env.size() << endl;
			//cout << "    - ideal1  :" << setw(15) << std::setprecision(5) << global_ideal1 / mpi_env.size() << endl;
			cout << "    - avg barrier:" << setw(15) << std::setprecision(5) << avg_barrier / mpi_env.size() << endl;
			cout << "    - min barrier:" << setw(15) << std::setprecision(5) << min_barrier << endl;
			cout << "    - max barrier:" << setw(15) << std::setprecision(5) << max_barrier << endl;
			cout << "File output =" << setw(15) << std::setprecision(5) << t_output << endl;
			cout << endl;
			cout << "Total CPU-time   =" << setw(15) << std::setprecision(5) << global_cpu / mpi_env.size() << endl;
			//for( size_t ipart=0; ipart < submesh_cfd.size(); ipart++ ){
			//	cout << "    - submesh =" << setw(5) << std::setprecision(5) << t_submesh[ipart] << endl;
			//}
			//for( int irank=0; irank < mpi_env.size(); irank++ ){
			//	cout << "    - MPI     =" << setw(10) << irank << setw(15) << std::setprecision(5) << t_compute[irank] << endl;
			//}
			cout << endl;
			cout << "Simulation finished successfully." << endl;
			cout << endl;
		} else {
			cout << endl;
			cout << "Simulation aborted. Please check your setting and rerun the case." << endl;
			cout << endl;
		}
	}

	mpi_env.barrier();
}

//***************************************************************************************************
template < typename PRECISION, unsigned DIM_CNT>
void Mesh<PRECISION, DIM_CNT> :: output_solution( PRECISION time, PRECISION PRINT_FILE_FREQUENCY, std::string file_base ){

	int myrank;
	MPI_Comm_rank( MPI_COMM_WORLD, &myrank );

	std::stringstream tmp_filename;

	if( myrank == 0 ) cout << "output solution =" << setw(15) << time << endl;

	if( myrank < 10 ){
		tmp_filename << file_base << "cfdv0_p00" << myrank;
	}else if( myrank < 100 ){
		tmp_filename << file_base << "cfdv0_p0"  << myrank;
	}else{
		tmp_filename << file_base << "cfdv0_p"   << myrank;
	}

	std::string filename;

	int ipart=0;
	for( auto it_sm = submesh_cfd.begin(); it_sm < submesh_cfd.end(); it_sm++ ){
		std::stringstream tmp;
		if( ipart < 10 ){
			tmp << tmp_filename.str() << "_s0" << ipart << "_t" << time / PRINT_FILE_FREQUENCY << ".plt";
		}else if( ipart < 100 ){
			tmp << tmp_filename.str() << "_s"  << ipart << "_t" << time / PRINT_FILE_FREQUENCY << ".plt";
		}
		filename = tmp.str();
		it_sm->output_cfd( time, filename, &gmsh_parts[ipart] );
		ipart++;
	}
}













































