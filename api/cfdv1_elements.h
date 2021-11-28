#ifndef CFDV1_ELEMENTS_H
#define CFDV1_ELEMENTS_H


//***************************************************************************************************
// Only contains solutions variables: q_old, q_new and delta_q
// 		- 2D, DP: 96 bytes
//***************************************************************************************************
template < typename PRECISION, unsigned DIM_CNT >
struct CFDv1_soln_vars{

	PRECISION q_old  [DIM_CNT+2];
	PRECISION q_new  [DIM_CNT+2];
	PRECISION delta_q[DIM_CNT+2];
};

//***************************************************************************************************
// Convective variables required by the minmod interpolation.
// 		- 2D, DP: 96 bytes
//***************************************************************************************************
template < typename PRECISION, unsigned DIM_CNT >
struct CFDv1_conv_vars{

	unsigned short int sm_id;
	unsigned       int id;

	PRECISION rho_grad [DIM_CNT];
	PRECISION rhoU_grad[DIM_CNT][DIM_CNT];
	PRECISION rhoE_grad[DIM_CNT];
	PRECISION Rpsi_grad[DIM_CNT];
	PRECISION c_grad   [DIM_CNT];
};

//***************************************************************************************************
// Viscous variables: used for all interpolation schemes
// 		- 2D, DP: 80 bytes
//***************************************************************************************************
template < typename PRECISION, unsigned DIM_CNT >
struct CFDv1_visc_vars{

	unsigned short int sm_id;
	unsigned       int id;

	PRECISION dudx[DIM_CNT][DIM_CNT];
	PRECISION dTdx[DIM_CNT];
	PRECISION tauMC[DIM_CNT][DIM_CNT];
};

//***************************************************************************************************
// Geometric parameters
//		- tri:  128 bytes
//		- quad: 168 bytes
//***************************************************************************************************
template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
struct CFDv1_params{

	PRECISION v_surf [FACE_CNT][DIM_CNT];
	PRECISION v_neigh[FACE_CNT][DIM_CNT];
	PRECISION weights[FACE_CNT];

	PRECISION vol_inv;
};


#endif






























