//***************************************************************************************************
// Function template declaration
//***************************************************************************************************

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// CFD cells
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//
template void MPI_env::set_halo_comm_type<float , 2, 3>(
										  enum t_mpi_halo_comm halo_comm_type,
										  fm_vector<CFDv0_cell<float ,2,3>>  cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<float ,2,3>>> ghost_hpath );
template void MPI_env::set_halo_comm_type<float , 2, 4>(
										  enum t_mpi_halo_comm halo_comm_type,
										  fm_vector<CFDv0_cell<float ,2,4>>  cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<float ,2,4>>> ghost_hpath );
template void MPI_env::set_halo_comm_type<double, 2, 3>(
										  enum t_mpi_halo_comm halo_comm_type,
										  fm_vector<CFDv0_cell<double,2,3>>  cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<double,2,3>>> ghost_hpath );
template void MPI_env::set_halo_comm_type<double, 2, 4>(
										  enum t_mpi_halo_comm halo_comm_type,
										  fm_vector<CFDv0_cell<double,2,4>>  cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<double,2,4>>> ghost_hpath );

// Send/recv, non-blocking
template void MPI_env::isendrecv<float , 2, 3>(
										  fm_vector<CFDv0_cell<float ,2,3>>  &cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<float ,2,3>>> &ghost_hpath );
template void MPI_env::isendrecv<float , 2, 4>(
										  fm_vector<CFDv0_cell<float ,2,4>>  &cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<float ,2,4>>> &ghost_hpath );
template void MPI_env::isendrecv<double, 2, 3>(
										  fm_vector<CFDv0_cell<double,2,3>>  &cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<double,2,3>>> &ghost_hpath );
template void MPI_env::isendrecv<double, 2, 4>(
										  fm_vector<CFDv0_cell<double,2,4>>  &cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<double,2,4>>> &ghost_hpath );

// Blocking
template void MPI_env::  neighbor_alltoallw<float , 2, 3>( CFDv0_cell<float ,2,3> *cells_hpath,
														   CFDv0_cell<float ,2,3> *ghost_hpath );
template void MPI_env::  neighbor_alltoallw<float , 2, 4>( CFDv0_cell<float ,2,4> *cells_hpath,
														   CFDv0_cell<float ,2,4> *ghost_hpath );
template void MPI_env::  neighbor_alltoallw<double, 2, 3>( CFDv0_cell<double,2,3> *cells_hpath,
														   CFDv0_cell<double,2,3> *ghost_hpath );
template void MPI_env::  neighbor_alltoallw<double, 2, 4>( CFDv0_cell<double,2,4> *cells_hpath,
														   CFDv0_cell<double,2,4> *ghost_hpath );

template void MPI_env:: ineighbor_alltoallw<float , 2, 3>( CFDv0_cell<float ,2,3> *cells_hpath,
														   CFDv0_cell<float ,2,3> *ghost_hpath );
template void MPI_env:: ineighbor_alltoallw<float , 2, 4>( CFDv0_cell<float ,2,4> *cells_hpath,
														   CFDv0_cell<float ,2,4> *ghost_hpath );
template void MPI_env:: ineighbor_alltoallw<double, 2, 3>( CFDv0_cell<double,2,3> *cells_hpath,
														   CFDv0_cell<double,2,3> *ghost_hpath );
template void MPI_env:: ineighbor_alltoallw<double, 2, 4>( CFDv0_cell<double,2,4> *cells_hpath,
														   CFDv0_cell<double,2,4> *ghost_hpath );

