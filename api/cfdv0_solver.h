#ifndef CFDV0_H
#define CFDV0_H

#include <memory>
#include <functional>
#include <array>
#include <sstream>
#include "elements.h"
#include "cfdv0_elements.h"
#include "mpi_env.h"
#include "gmsh_reader.h"

//***************************************************************************************************
//***************************************************************************************************
//
//***************************************************************************************************
//***************************************************************************************************
template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
class CFDv0_solver{

	public:

		//
		int sm_id;

		// Runge-Kutta variables
		PRECISION ti, tf, dt;

		PRECISION Ak[1] = { 0.0 };
		PRECISION Bk[1] = { 1.0 };
		//PRECISION Ak[5] = { 0.0            , -0.4178904745   , -1.192151694643 , -1.697784692471 , -1.514183444257  };
		//PRECISION Bk[5] = { 0.1496590219993,  0.3792103129999,  0.8229550293869,  0.6994504559488,  0.1530572479681 };

		// Cells
		fm_vector<CFDv0_cell<PRECISION, DIM_CNT, FACE_CNT>> cells_cfd;

		// Dummy cell: used for NEIGHBOR_NONE access
		fm_vector<CFDv0_cell<PRECISION, DIM_CNT, FACE_CNT>> dummy_cell;
		t_solution_vars<PRECISION, DIM_CNT> dummy_var;

		std::vector<std::vector<PRECISION>> xc;

		// Ghost cells
		fm_vector<fm_vector<CFDv0_cell<PRECISION, DIM_CNT, FACE_CNT>>> ghost_mpi;
		fm_vector<fm_vector<CFDv0_cell<PRECISION, DIM_CNT, FACE_CNT>>> ghost_bnd;

		// Ghost boundary cells
		std::vector<std::vector<array<int, 2>>> boundaries;

		// Mapping from rank ID to local MPI neighbor ID
		std::vector<int> rank2local;

		// Constants
		PRECISION ZERO  =  0.0;
		PRECISION HALF  =  0.5;
		PRECISION ONE   =  1.0;
		PRECISION M_ONE = -1.0;

		PRECISION Pr, Re, gamma, L, M, gamma_m_one, Pr_inv, gamma_m_one_inv, AoA;
		PRECISION mu0, T0, S;
		PRECISION k, Cp, Cv;

		PRECISION rho_i, u_i[DIM_CNT], E_i, u_mag_i, p_i, T_i;		// Initial  conditions variables.
		PRECISION rho_b, u_b[DIM_CNT], E_b, u_mag_b, p_b, T_b;		// Boundary conditions variables.
		const PRECISION Runiversal = 8.31447 ;						// Universal gas constant. TODO: one header for all constants
		const PRECISION molWeight = 28.9647 ;						// Molar Weight of the air in g/mol.
		const PRECISION Rgas = Runiversal / molWeight * 1000.0 ;	// Specific gas constant.
		const PRECISION Rgas_inv = ONE / Rgas;

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Functions
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		void restart( PRECISION t_start, PRECISION PRINT_FILE_FREQUENCY, int submesh );
		void prepare_for_timestep( int rk_step );
		void calc_gradients      ( MPI_env &mpi_env );
		void calc_AD             ();
		void calc_VIS            ( MPI_env &mpi_env );
		void one_rk_step_v1      ( int rk_step, PRECISION dt, MPI_env &mpi_env );
		void exchange_ghost_cells( MPI_env &mpi_env );
		void gather_info         ( MPI_env &mpi_env, struct t_domain_info<PRECISION, DIM_CNT> &domain_info );

		void finalize();
		void set_boundary_conditions();

		// gmsh
		void allocate       ( gmsh_mesh *mesh, std::vector<            t_solution_vars<PRECISION, DIM_CNT>*> &cells_ptr );
		void allocate_ghost_cells( MPI_env &mpi_env,
								   gmsh_mesh *mesh,
								   std::vector<std::vector<t_solution_vars<PRECISION,DIM_CNT>*>> &ghost_ptr,
								   std::vector<std::vector<t_solution_vars<PRECISION,DIM_CNT>*>> &ghost_ptr_bnd );
		void assign_pointers( gmsh_mesh *mesh,
							  std::vector<std::vector<t_solution_vars<PRECISION, DIM_CNT>*>> cells_ptr,
							  std::vector<std::vector<t_solution_vars<PRECISION, DIM_CNT>*>> ghost_ptr,
							  std::vector<std::vector<t_solution_vars<PRECISION, DIM_CNT>*>> ghost_bnd );
		void init_params    ( gmsh_mesh *mesh );
		void initialize     ( MPI_env &mpi_env, gmsh_mesh *mesh );
		void output_cfd     ( PRECISION time, std::string filename, gmsh_mesh *mesh );
		void initialize_mpi_env( MPI_env &mpi_env, gmsh_mesh *mesh );

		//PRECISION compute_cfl( PRECISION dt );
		PRECISION compute_dt( PRECISION cflMax );

		void halo_comm_wait( MPI_env &mpi_env );
	private:

		MPI_Datatype type_cfdcell;

		void init_mpi_types( MPI_env &mpi_env );
		void set_init_conditions( initial_condition init_case );

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Fetch cell neighbor
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		inline const CFDv0_cell<PRECISION, DIM_CNT, FACE_CNT>* getCellNeighbor(
					        CFDv0_cell<PRECISION, DIM_CNT, FACE_CNT> *current,
							unsigned                                  neighbor_index ){
			if (neighbor_index < 2) {
				return current + ((2 * (int)neighbor_index) - 1);
			}
			return &cells_cfd[current->neighbors[neighbor_index-2]] ;
		}

		inline const t_solution_vars<PRECISION, DIM_CNT>* getCellVars(
					        CFDv0_cell<PRECISION, DIM_CNT, FACE_CNT> *current,
							unsigned                                  neighbor_index ){

			return current->neighs[neighbor_index];
		}

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Compute primitives variables from conservative vector
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		inline void compute_primitives( PRECISION q_sol[DIM_CNT+2], PRECISION &rho_inv,
										PRECISION Uvec[DIM_CNT], PRECISION &Umag_sqr,
										PRECISION &energy, PRECISION &p, PRECISION &Rpsi ){

			rho_inv = ONE / q_sol[0];

			Uvec[0]  = q_sol[1] * rho_inv;
			Umag_sqr = Uvec[0] * Uvec[0];
			for( unsigned idim = 1; idim < DIM_CNT; idim++ ){
				Uvec[idim] = q_sol[idim+1] * rho_inv ;
				Umag_sqr += Uvec[idim] * Uvec[idim];
			}
			energy = q_sol[DIM_CNT+1] * rho_inv ;
			p      = q_sol[0] * (gamma-1.) * (energy - 0.5 * Umag_sqr);
			Rpsi   = p * rho_inv ;
		}

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		inline PRECISION compute_Rpsi( PRECISION q_sol[DIM_CNT+2] ){

			PRECISION rho_inv  = ONE / q_sol[0];
			PRECISION rhoU_sqr = q_sol[1] * q_sol[1];
			for( unsigned idim = 1; idim < DIM_CNT; idim++ ){
				rhoU_sqr += q_sol[idim+1] * q_sol[idim+1];
			}

			PRECISION Rpsi = gamma_m_one * (q_sol[DIM_CNT+1] - 0.5 * rhoU_sqr * rho_inv) * rho_inv;
			return Rpsi;
		}

		inline PRECISION dot_product( PRECISION lvec[DIM_CNT], PRECISION rvec[DIM_CNT] ){

			PRECISION dot_prod = lvec[0] * rvec[0];
			for( unsigned idim=1; idim < DIM_CNT; idim++ ){
				dot_prod += lvec[idim] * rvec[idim];
			}
			return dot_prod;
		}

		inline PRECISION vector_mag( PRECISION vec[DIM_CNT] ){

			PRECISION vec_mag = vec[0] * vec[0];
			for( unsigned idim=1; idim < DIM_CNT; idim++ ){
				vec_mag += vec[idim] * vec[idim];
			}
			return sqrt( vec_mag );
		}

		inline void init_array( PRECISION array[DIM_CNT] ){
			for( unsigned idim=0; idim < DIM_CNT; idim++ ){
				array[idim] = 0.;
			}
		}
		inline void init_array( PRECISION array[DIM_CNT][DIM_CNT] ){
			for( unsigned idim=0; idim < DIM_CNT; idim++ ){
				for( unsigned jdim=0; jdim < DIM_CNT; jdim++ ){
					array[idim][jdim] = 0.;
				}
			}
		}

		inline void compute_tauMC( PRECISION dudx[DIM_CNT][DIM_CNT], PRECISION divU, PRECISION tauMC[DIM_CNT][DIM_CNT] ){

			PRECISION dev2[DIM_CNT][DIM_CNT];
			if( DIM_CNT == 2 ){
				dev2[0][0] = dudx[0][0] - 2.0 / 3.0 * divU;
				dev2[0][1] = dudx[1][0];
				dev2[1][0] = dudx[0][1];
				dev2[1][1] = dudx[1][1] - 2.0 / 3.0 * divU;
			}

			// tauMC
			for( unsigned idim = 0; idim < DIM_CNT; idim++ ){
			for( unsigned jdim = 0; jdim < DIM_CNT; jdim++ ){
				tauMC[idim][jdim] += mu0 * dev2[idim][jdim];
			}
			}
		}

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Linear interpolation scheme
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		inline PRECISION interp_linear( PRECISION weight, PRECISION cell_phi, PRECISION adjc_phi ){
			return weight*cell_phi + (ONE-weight)*adjc_phi;
		}

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Minmod interpolation scheme
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		// Get sign of variable a
		inline PRECISION sign(PRECISION a) {
			return a < 0.0 ? -1.0 : 1.0 ;
		}

		// Compute "r"
		inline PRECISION calc_r( PRECISION phiP, PRECISION phiN, PRECISION phiGrad[DIM_CNT], PRECISION d[DIM_CNT] ){

			PRECISION gradf  = phiN - phiP + 1.0e-30;
			PRECISION gradcf = ZERO;

			for (unsigned idim=0; idim<DIM_CNT; idim++) {
				gradcf += d[idim] * phiGrad[idim];
			}

			if (abs(gradcf) >= 1000.0 * abs(gradf)) {
				return 2.0 * 1000.0 * sign(gradcf) * sign(gradf) - ONE;
			} else {
				return 2.0 * (gradcf / gradf) - ONE;
			}
		}

		// Minmod interpolation
		inline PRECISION interp_minmod( PRECISION cell_phi, PRECISION adjc_phi,
										PRECISION grad_phi[DIM_CNT], PRECISION d[DIM_CNT],
										PRECISION weight_linear, PRECISION flux ){

			PRECISION r       = calc_r(cell_phi, adjc_phi, grad_phi, d);
			PRECISION limiter = max(min(r, ONE), ZERO);
			PRECISION weight  = limiter * weight_linear + ( ONE - limiter ) * (ONE + flux) * HALF;
			//PRECISION weight = limiter * weight_linear + ( ONE - limiter );
			//weight = HALF * (ONE - flux) + flux * weight;

			return weight * cell_phi + (ONE - weight) * adjc_phi;
		}
};

//***************************************************************************************************
//***************************************************************************************************
// var_faced_solver
// 		- C++ trick to store solver objects of different template parameters in same vector
// 		- Do we really need this? If we define our own vector, we could by pass all that.
//***************************************************************************************************
//***************************************************************************************************
struct var_faced_solver {
	template <class T>
	var_faced_solver(T&& t) {
		auto ptr = std::make_shared<std::remove_reference_t<T>>(std::forward<T>(t));

		halo_comm_wait = [ptr]( MPI_env &mpi_env ){
			return ptr->halo_comm_wait( mpi_env );
		};

		set_boundary_conditions = [ptr](){
			return ptr->set_boundary_conditions();
		};

		/*compute_cfl = [ptr]( double dt ){
			return ptr->compute_cfl( dt );
		};*/
		compute_dt = [ptr]( double cflMax ){
			return ptr->compute_dt( cflMax );
		};

		output_cfd = [ptr]( double time, std::string filename, gmsh_mesh *mesh ){
			return ptr->output_cfd( time, filename, mesh );
		};

		restart = [ptr]( double t_start, double PRINT_FILE_FREQUENCY, int submesh ){
			return ptr->restart( t_start, PRINT_FILE_FREQUENCY, submesh );
		};

		prepare_for_timestep = [ptr]( int rk_step ){
			return ptr->prepare_for_timestep( rk_step );
		};

		// Gradients
		calc_gradients = [ptr]( MPI_env &mpi_env ){
			return ptr->calc_gradients( mpi_env );
		};

		// Artificial Diffusion
		calc_AD = [ptr](){
			return ptr->calc_AD();
		};

		// Viscosity
		calc_VIS = [ptr]( MPI_env &mpi_env ){
			return ptr->calc_VIS( mpi_env );
		};

		one_rk_step_v1 = [ptr]( int rk_step, double dt, MPI_env &mpi_env ){
			return ptr->one_rk_step_v1( rk_step, dt, mpi_env );
		};

		exchange_ghost_cells = [ptr]( MPI_env &mpi_env ){
			return ptr->exchange_ghost_cells( mpi_env );
		};

		gather_info = [ptr]( MPI_env &mpi_env, struct t_domain_info<double, 2> &domain_info ){
			return ptr->gather_info( mpi_env, domain_info );
		};

		submesh_obj = &((typeof(T))ptr);
	}

	void* submesh_obj;

	std::function<void(MPI_env&)> halo_comm_wait;
	std::function<void()> set_boundary_conditions;
	std::function<void(double, std::string, gmsh_mesh*)> output_cfd;
	//std::function<double(double)> compute_cfl;
	std::function<double(double)> compute_dt;

	std::function<void(double,double,int)> restart;
	std::function<void(int)> prepare_for_timestep;
	std::function<void(MPI_env&)> calc_gradients;
	std::function<void()> calc_AD;
	std::function<void(MPI_env&)> calc_VIS;
	std::function<void(int,double,MPI_env&)> one_rk_step_v1;
	std::function<void(MPI_env&)> exchange_ghost_cells;
	std::function<void(MPI_env&, struct t_domain_info<double,2>& )> gather_info;
};

#endif






























