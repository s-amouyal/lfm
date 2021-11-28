#ifndef CFDV1_H
#define CFDV1_H

#include <vector>
#include <array>
#include "api/constants.h"
#include "api/cfdv1_elements.h"
#include "api/mpi_env.h"

//***************************************************************************************************
//***************************************************************************************************
template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
class CFDv1_solver{

	public:

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Cells and neighbors
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		// Cells
		std::vector<CFDv1_soln_vars<PRECISION,DIM_CNT>> cells_soln;
		std::vector<CFDv1_conv_vars<PRECISION,DIM_CNT>> cells_conv;
		std::vector<CFDv1_visc_vars<PRECISION,DIM_CNT>> cells_visc;
		std::vector<CFDv1_params   <PRECISION,DIM_CNT,FACE_CNT>> cells_params;

		// MPI ghost cells
		std::vector<std::vector<CFDv1_soln_vars<PRECISION,DIM_CNT>>> ghost_mpi_soln;
		std::vector<std::vector<CFDv1_conv_vars<PRECISION,DIM_CNT>>> ghost_mpi_conv;
		std::vector<std::vector<CFDv1_visc_vars<PRECISION,DIM_CNT>>> ghost_mpi_visc;

		// Boundary ghost cells
		std::vector<std::vector<array<int, 2>>> boundaries;
		std::vector<std::vector<CFDv1_soln_vars<PRECISION,DIM_CNT>>> ghost_bnd_soln;
		std::vector<std::vector<CFDv1_conv_vars<PRECISION,DIM_CNT>>> ghost_bnd_conv;
		std::vector<std::vector<CFDv1_visc_vars<PRECISION,DIM_CNT>>> ghost_bnd_visc;

		// Neighbors Pointers
		CFDv1_soln_vars<PRECISION,DIM_CNT>* neigh_soln[FACE_CNT];
		CFDv1_conv_vars<PRECISION,DIM_CNT>* neigh_conv[FACE_CNT];
		CFDv1_visc_vars<PRECISION,DIM_CNT>* neigh_visc[FACE_CNT];

		// Cell centers: used for initial conditions
		std::vector<std::vector<PRECISION>> xc;

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Solver variables
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		// Submesh ID
		int sm_id;

		// RK45 coefficients
		PRECISION Ak[1] = { 0.0 };
		PRECISION Bk[1] = { 1.0 };

		// Time
		PRECISION ti, tf, dt;

		// Fluid domain-constants
		PRECISION Pr, Re, gamma, L, M, gamma_m_one, Pr_inv, gamma_m_one_inv, AoA;
		PRECISION mu0, T0, S;
		PRECISION k, Cp, Cv;

		// Initial/Boundary conditions
		PRECISION rho_i, u_i[DIM_CNT], E_i, u_mag_i, p_i, T_i;
		PRECISION rho_b, u_b[DIM_CNT], E_b, u_mag_b, p_b, T_b;

		// Constants
		const PRECISION Runiversal = 8.31447 ;						// Universal gas constant.
		const PRECISION molWeight = 28.9647 ;						// Molar Weight of the air in g/mol.
		const PRECISION Rgas = Runiversal / molWeight * 1000.0 ;	// Specific gas constant.
		const PRECISION Rgas_inv = ONE / Rgas;						// Inverse specific gas constant

		// MPI variables
		std::vector<int> rank2local;

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Solver functions
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

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Solver preparation functions
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		void allocate( MPI_env mpi_env,
					   gmsh_mesh *mesh,
					   std::vector<CFDv1_soln_vars<PRECISION,DIM_CNT>*> &addr_soln,
					   std::vector<CFDv1_conv_vars<PRECISION,DIM_CNT>*> &addr_conv,
					   std::vector<CFDv1_visc_vars<PRECISION,DIM_CNT>*> &addr_visc );

		void allocate_ghost_cells( MPI_env &mpi_env,
								   gmsh_mesh *mesh,
								   std::vector<std::vector<CFDv1_soln_vars<PRECISION,DIM_CNT>*>> &ptr_mpi_soln,
								   std::vector<std::vector<CFDv1_conv_vars<PRECISION,DIM_CNT>*>> &ptr_mpi_conv,
								   std::vector<std::vector<CFDv1_visc_vars<PRECISION,DIM_CNT>*>> &ptr_mpi_visc,
								   std::vector<std::vector<CFDv1_soln_vars<PRECISION,DIM_CNT>*>> &ptr_bnd_soln,
								   std::vector<std::vector<CFDv1_conv_vars<PRECISION,DIM_CNT>*>> &ptr_bnd_conv,
								   std::vector<std::vector<CFDv1_visc_vars<PRECISION,DIM_CNT>*>> &ptr_bnd_visc );
		void assign_pointers( gmsh_mesh *mesh,
							  std::vector<std::vector<CFDv1_soln_vars<PRECISION, DIM_CNT>*>>  ptr_cell_soln,
							  std::vector<std::vector<CFDv1_conv_vars<PRECISION, DIM_CNT>*>>  ptr_cell_conv,
							  std::vector<std::vector<CFDv1_visc_vars<PRECISION, DIM_CNT>*>>  ptr_cell_visc,
							  std::vector<std::vector<CFDv1_soln_vars<PRECISION, DIM_CNT>*>>  ptr_mpi_soln,
							  std::vector<std::vector<CFDv1_conv_vars<PRECISION, DIM_CNT>*>>  ptr_mpi_conv,
							  std::vector<std::vector<CFDv1_visc_vars<PRECISION, DIM_CNT>*>>  ptr_mpi_visc,
							  std::vector<std::vector<CFDv1_soln_vars<PRECISION, DIM_CNT>*>>  ptr_bnd_soln,
							  std::vector<std::vector<CFDv1_conv_vars<PRECISION, DIM_CNT>*>>  ptr_bnd_conv,
							  std::vector<std::vector<CFDv1_visc_vars<PRECISION, DIM_CNT>*>>  ptr_bnd_visc );

		void init_params    ( gmsh_mesh *mesh );
		void initialize     ( MPI_env &mpi_env, gmsh_mesh *mesh );
		void output_cfd     ( PRECISION time, std::string filename, gmsh_mesh *mesh );
		void initialize_mpi_env( MPI_env &mpi_env, gmsh_mesh *mesh );

	private:

		MPI_Datatype type_soln, type_conv, type_visc;

		void init_mpi_types( MPI_env &mpi_env );
		void set_init_conditions( initial_condition init_case );

		//*******************************************************************************************
		// Inline functions
		//*******************************************************************************************

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Compute primitives variables from conservative vector
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		inline void compute_primitives( PRECISION q_sol[DIM_CNT+2], PRECISION &rho_inv,
										PRECISION Uvec[DIM_CNT], PRECISION &Umag_sqr, PRECISION &energy,
										PRECISION &p, PRECISION &Rpsi ){

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
};

////***************************************************************************************************
//// Compute Rpsi
////***************************************************************************************************
//inline PRECISION compute_Rpsi( PRECISION q_sol[DIM_CNT+2] ){
//
//	PRECISION rho_inv  = ONE / q_sol[0];
//	PRECISION rhoU_sqr = q_sol[1] * q_sol[1];
//	for( unsigned idim = 1; idim < DIM_CNT; idim++ ){
//		rhoU_sqr += q_sol[idim+1] * q_sol[idim+1];
//	}
//
//	PRECISION Rpsi = gamma_m_one * (q_sol[DIM_CNT+1] - 0.5 * rhoU_sqr * rho_inv) * rho_inv;
//	return Rpsi;
//}
//
////***************************************************************************************************
//// Dot product
////***************************************************************************************************
//inline PRECISION dot_product( PRECISION lvec[DIM_CNT], PRECISION rvec[DIM_CNT] ){
//
//	PRECISION dot_prod = lvec[0] * rvec[0];
//	for( unsigned idim=1; idim < DIM_CNT; idim++ ){
//		dot_prod += lvec[idim] * rvec[idim];
//	}
//	return dot_prod;
//}
//
////***************************************************************************************************
//// Vector magnitude
////***************************************************************************************************
//inline PRECISION vector_mag( PRECISION vec[DIM_CNT] ){
//
//	PRECISION vec_mag = vec[0] * vec[0];
//	for( unsigned idim=1; idim < DIM_CNT; idim++ ){
//		vec_mag += vec[idim] * vec[idim];
//	}
//	return sqrt( vec_mag );
//}
//
////***************************************************************************************************
//// Initialization functions
////***************************************************************************************************
//
//// Vector
//inline void init_array( PRECISION array[DIM_CNT] ){
//	for( unsigned idim=0; idim < DIM_CNT; idim++ ){
//		array[idim] = 0.;
//	}
//}
//
//// Matrix
//inline void init_array( PRECISION array[DIM_CNT][DIM_CNT] ){
//	for( unsigned idim=0; idim < DIM_CNT; idim++ ){
//		for( unsigned jdim=0; jdim < DIM_CNT; jdim++ ){
//			array[idim][jdim] = 0.;
//		}
//	}
//}
//
////***************************************************************************************************
//// Compute viscous stress tensor from velocity derivative tensor, divergence and tauMC
////***************************************************************************************************
//inline void compute_tauMC( PRECISION dudx[DIM_CNT][DIM_CNT], PRECISION divU, PRECISION tauMC[DIM_CNT][DIM_CNT] ){
//
//	PRECISION dev2[DIM_CNT][DIM_CNT];
//	if( DIM_CNT == 2 ){
//		dev2[0][0] = dudx[0][0] - 2.0 / 3.0 * divU;
//		dev2[0][1] = dudx[1][0];
//		dev2[1][0] = dudx[0][1];
//		dev2[1][1] = dudx[1][1] - 2.0 / 3.0 * divU;
//	}
//
//	// tauMC
//	for( unsigned idim = 0; idim < DIM_CNT; idim++ ){
//	for( unsigned jdim = 0; jdim < DIM_CNT; jdim++ ){
//		tauMC[idim][jdim] += mu0 * dev2[idim][jdim];
//	}
//	}
//}
//
////***************************************************************************************************
//// Linear interpolation scheme
////***************************************************************************************************
//inline PRECISION interp_linear( PRECISION weight, PRECISION cell_phi, PRECISION adjc_phi ){
//	return weight*cell_phi + (ONE-weight)*adjc_phi;
//}
//
////***************************************************************************************************
//// Minmod interpolation scheme
////***************************************************************************************************
//
//// Get sign of variable a
//inline PRECISION sign(PRECISION a) {
//	return a < 0.0 ? -1.0 : 1.0 ;
//}
//
//// Compute "r"
//inline PRECISION calc_r( PRECISION phiP, PRECISION phiN, PRECISION phiGrad[DIM_CNT], PRECISION d[DIM_CNT] ){
//
//	PRECISION gradf  = phiN - phiP + 1.0e-30;
//	PRECISION gradcf = ZERO;
//
//	for (unsigned idim=0; idim<DIM_CNT; idim++) {
//		gradcf += d[idim] * phiGrad[idim];
//	}
//
//	if (abs(gradcf) >= 1000.0 * abs(gradf)) {
//		return 2.0 * 1000.0 * sign(gradcf) * sign(gradf) - ONE;
//	} else {
//		return 2.0 * (gradcf / gradf) - ONE;
//	}
//}
//
//// Minmod interpolation
//inline PRECISION interp_minmod( PRECISION cell_phi, PRECISION adjc_phi, PRECISION grad_phi[DIM_CNT],
//								PRECISION d[DIM_CNT], PRECISION weight_linear, PRECISION flux ){
//
//	PRECISION r       = calc_r(cell_phi, adjc_phi, grad_phi, d);
//	PRECISION limiter = max(min(r, ONE), ZERO);
//	PRECISION weight  = limiter * weight_linear + ( ONE - limiter ) * (ONE + flux) * HALF;
//
//	return weight * cell_phi + (ONE - weight) * adjc_phi;
//}

#endif






























