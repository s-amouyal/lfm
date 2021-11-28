#ifndef CFDV1_SOLVER_H
#define CFDV1_SOLVER_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <initializer_list>
#include <vector>

#include "math.h"
#include "api/cfdv1_solver.h"
#include "api/mesh_writer.h"

using namespace std;

// Allocate
template void CFDv1_solver<float , 2, 3> :: allocate( MPI_env mpi_env, gmsh_mesh *mesh, std::vector<CFDv1_soln_vars<float ,2>*> &addr_soln, std::vector<CFDv1_conv_vars<float ,2>*> &addr_conv, std::vector<CFDv1_visc_vars<float ,2>*> &addr_visc );
template void CFDv1_solver<float , 2, 4> :: allocate( MPI_env mpi_env, gmsh_mesh *mesh, std::vector<CFDv1_soln_vars<float ,2>*> &addr_soln, std::vector<CFDv1_conv_vars<float ,2>*> &addr_conv, std::vector<CFDv1_visc_vars<float ,2>*> &addr_visc );
template void CFDv1_solver<double, 2, 3> :: allocate( MPI_env mpi_env, gmsh_mesh *mesh, std::vector<CFDv1_soln_vars<double,2>*> &addr_soln, std::vector<CFDv1_conv_vars<double,2>*> &addr_conv, std::vector<CFDv1_visc_vars<double,2>*> &addr_visc );
template void CFDv1_solver<double, 2, 4> :: allocate( MPI_env mpi_env, gmsh_mesh *mesh, std::vector<CFDv1_soln_vars<double,2>*> &addr_soln, std::vector<CFDv1_conv_vars<double,2>*> &addr_conv, std::vector<CFDv1_visc_vars<double,2>*> &addr_visc );

template void CFDv1_solver<float , 2, 3> :: allocate_ghost_cells( MPI_env &mpi_env, gmsh_mesh *mesh, std::vector<std::vector<CFDv1_soln_vars<float ,2>*>> &ptr_mpi_soln, std::vector<std::vector<CFDv1_conv_vars<float ,2>*>> &ptr_mpi_conv, std::vector<std::vector<CFDv1_visc_vars<float ,2>*>> &ptr_mpi_visc, std::vector<std::vector<CFDv1_soln_vars<float ,2>*>> &ptr_bnd_soln, std::vector<std::vector<CFDv1_conv_vars<float ,2>*>> &ptr_bnd_conv, std::vector<std::vector<CFDv1_visc_vars<float ,2>*>> &ptr_bnd_visc );
template void CFDv1_solver<float , 2, 4> :: allocate_ghost_cells( MPI_env &mpi_env, gmsh_mesh *mesh, std::vector<std::vector<CFDv1_soln_vars<float ,2>*>> &ptr_mpi_soln, std::vector<std::vector<CFDv1_conv_vars<float ,2>*>> &ptr_mpi_conv, std::vector<std::vector<CFDv1_visc_vars<float ,2>*>> &ptr_mpi_visc, std::vector<std::vector<CFDv1_soln_vars<float ,2>*>> &ptr_bnd_soln, std::vector<std::vector<CFDv1_conv_vars<float ,2>*>> &ptr_bnd_conv, std::vector<std::vector<CFDv1_visc_vars<float ,2>*>> &ptr_bnd_visc );
template void CFDv1_solver<double, 2, 3> :: allocate_ghost_cells( MPI_env &mpi_env, gmsh_mesh *mesh, std::vector<std::vector<CFDv1_soln_vars<double,2>*>> &ptr_mpi_soln, std::vector<std::vector<CFDv1_conv_vars<double,2>*>> &ptr_mpi_conv, std::vector<std::vector<CFDv1_visc_vars<double,2>*>> &ptr_mpi_visc, std::vector<std::vector<CFDv1_soln_vars<double,2>*>> &ptr_bnd_soln, std::vector<std::vector<CFDv1_conv_vars<double,2>*>> &ptr_bnd_conv, std::vector<std::vector<CFDv1_visc_vars<double,2>*>> &ptr_bnd_visc );
template void CFDv1_solver<double, 2, 4> :: allocate_ghost_cells( MPI_env &mpi_env, gmsh_mesh *mesh, std::vector<std::vector<CFDv1_soln_vars<double,2>*>> &ptr_mpi_soln, std::vector<std::vector<CFDv1_conv_vars<double,2>*>> &ptr_mpi_conv, std::vector<std::vector<CFDv1_visc_vars<double,2>*>> &ptr_mpi_visc, std::vector<std::vector<CFDv1_soln_vars<double,2>*>> &ptr_bnd_soln, std::vector<std::vector<CFDv1_conv_vars<double,2>*>> &ptr_bnd_conv, std::vector<std::vector<CFDv1_visc_vars<double,2>*>> &ptr_bnd_visc );
//
//template void CFDv1_solver<float , 2, 3> :: assign_pointers( gmsh_mesh *mesh, std::vector<std::vector<CFDv1_soln_vars<float , 2>*>>  ptr_cell_soln, std::vector<std::vector<CFDv1_conv_vars<float , 2>*>>  ptr_cell_conv, std::vector<std::vector<CFDv1_visc_vars<float , 2>*>>  ptr_cell_visc, std::vector<std::vector<CFDv1_soln_vars<float , 2>*>>  ptr_mpi_soln, std::vector<std::vector<CFDv1_conv_vars<float , 2>*>>  ptr_mpi_conv, std::vector<std::vector<CFDv1_visc_vars<float , 2>*>>  ptr_mpi_visc, std::vector<std::vector<CFDv1_soln_vars<float , 2>*>>  ptr_bnd_soln, std::vector<std::vector<CFDv1_conv_vars<float , 2>*>>  ptr_bnd_conv, std::vector<std::vector<CFDv1_visc_vars<float , 2>*>>  ptr_bnd_visc );
//template void CFDv1_solver<float , 2, 4> :: assign_pointers( gmsh_mesh *mesh, std::vector<std::vector<CFDv1_soln_vars<float , 2>*>>  ptr_cell_soln, std::vector<std::vector<CFDv1_conv_vars<float , 2>*>>  ptr_cell_conv, std::vector<std::vector<CFDv1_visc_vars<float , 2>*>>  ptr_cell_visc, std::vector<std::vector<CFDv1_soln_vars<float , 2>*>>  ptr_mpi_soln, std::vector<std::vector<CFDv1_conv_vars<float , 2>*>>  ptr_mpi_conv, std::vector<std::vector<CFDv1_visc_vars<float , 2>*>>  ptr_mpi_visc, std::vector<std::vector<CFDv1_soln_vars<float , 2>*>>  ptr_bnd_soln, std::vector<std::vector<CFDv1_conv_vars<float , 2>*>>  ptr_bnd_conv, std::vector<std::vector<CFDv1_visc_vars<float , 2>*>>  ptr_bnd_visc );
//template void CFDv1_solver<double, 2, 3> :: assign_pointers( gmsh_mesh *mesh, std::vector<std::vector<CFDv1_soln_vars<double, 2>*>>  ptr_cell_soln, std::vector<std::vector<CFDv1_conv_vars<double, 2>*>>  ptr_cell_conv, std::vector<std::vector<CFDv1_visc_vars<double, 2>*>>  ptr_cell_visc, std::vector<std::vector<CFDv1_soln_vars<double, 2>*>>  ptr_mpi_soln, std::vector<std::vector<CFDv1_conv_vars<double, 2>*>>  ptr_mpi_conv, std::vector<std::vector<CFDv1_visc_vars<double, 2>*>>  ptr_mpi_visc, std::vector<std::vector<CFDv1_soln_vars<double, 2>*>>  ptr_bnd_soln, std::vector<std::vector<CFDv1_conv_vars<double, 2>*>>  ptr_bnd_conv, std::vector<std::vector<CFDv1_visc_vars<double, 2>*>>  ptr_bnd_visc );
//template void CFDv1_solver<double, 2, 4> :: assign_pointers( gmsh_mesh *mesh, std::vector<std::vector<CFDv1_soln_vars<double, 2>*>>  ptr_cell_soln, std::vector<std::vector<CFDv1_conv_vars<double, 2>*>>  ptr_cell_conv, std::vector<std::vector<CFDv1_visc_vars<double, 2>*>>  ptr_cell_visc, std::vector<std::vector<CFDv1_soln_vars<double, 2>*>>  ptr_mpi_soln, std::vector<std::vector<CFDv1_conv_vars<double, 2>*>>  ptr_mpi_conv, std::vector<std::vector<CFDv1_visc_vars<double, 2>*>>  ptr_mpi_visc, std::vector<std::vector<CFDv1_soln_vars<double, 2>*>>  ptr_bnd_soln, std::vector<std::vector<CFDv1_conv_vars<double, 2>*>>  ptr_bnd_conv, std::vector<std::vector<CFDv1_visc_vars<double, 2>*>>  ptr_bnd_visc );
//
//template void CFDv1_solver<float , 2, 3> :: initialize( MPI_env &mpi_env, gmsh_mesh *mesh );
//template void CFDv1_solver<float , 2, 4> :: initialize( MPI_env &mpi_env, gmsh_mesh *mesh );
//template void CFDv1_solver<double, 2, 3> :: initialize( MPI_env &mpi_env, gmsh_mesh *mesh );
//template void CFDv1_solver<double, 2, 4> :: initialize( MPI_env &mpi_env, gmsh_mesh *mesh );

//// Restart simulation from t_start
//template void CFDv0_solver<float , 2, 3> :: restart( float  t_start, float  PRINT_FILE_FREQUENCY, int submesh );
//template void CFDv0_solver<float , 2, 4> :: restart( float  t_start, float  PRINT_FILE_FREQUENCY, int submesh );
//template void CFDv0_solver<double, 2, 3> :: restart( double t_start, double PRINT_FILE_FREQUENCY ,int submesh );
//template void CFDv0_solver<double, 2, 4> :: restart( double t_start, double PRINT_FILE_FREQUENCY, int submesh );
//
//// Prepare time-step
//template void CFDv0_solver<float , 2, 3> :: prepare_for_timestep( int rk_step );
//template void CFDv0_solver<float , 2, 4> :: prepare_for_timestep( int rk_step );
//template void CFDv0_solver<double, 2, 3> :: prepare_for_timestep( int rk_step );
//template void CFDv0_solver<double, 2, 4> :: prepare_for_timestep( int rk_step );
//
//// Calculate Gradients
//template void CFDv0_solver<float , 2, 3> :: calc_gradients( MPI_env &mpi_env );
//template void CFDv0_solver<float , 2, 4> :: calc_gradients( MPI_env &mpi_env );
//template void CFDv0_solver<double, 2, 3> :: calc_gradients( MPI_env &mpi_env );
//template void CFDv0_solver<double, 2, 4> :: calc_gradients( MPI_env &mpi_env );
//
//// AD
//template void CFDv0_solver<float , 2, 3> :: calc_AD();
//template void CFDv0_solver<float , 2, 4> :: calc_AD();
//template void CFDv0_solver<double, 2, 3> :: calc_AD();
//template void CFDv0_solver<double, 2, 4> :: calc_AD();
//
//// Viscosity
//template void CFDv0_solver<float , 2, 3> :: calc_VIS( MPI_env &mpi_env );
//template void CFDv0_solver<float , 2, 4> :: calc_VIS( MPI_env &mpi_env );
//template void CFDv0_solver<double, 2, 3> :: calc_VIS( MPI_env &mpi_env );
//template void CFDv0_solver<double, 2, 4> :: calc_VIS( MPI_env &mpi_env );
//
//// One RK step
//template void CFDv0_solver<float , 2, 3> :: one_rk_step_v1( int rk_step, float  dt, MPI_env &mpi_env );
//template void CFDv0_solver<float , 2, 4> :: one_rk_step_v1( int rk_step, float  dt, MPI_env &mpi_env );
//template void CFDv0_solver<double, 2, 3> :: one_rk_step_v1( int rk_step, double dt, MPI_env &mpi_env );
//template void CFDv0_solver<double, 2, 4> :: one_rk_step_v1( int rk_step, double dt, MPI_env &mpi_env );
//
//// Set boundary conditions (BC)
//template void CFDv0_solver<float , 2, 3> :: set_boundary_conditions();
//template void CFDv0_solver<float , 2, 4> :: set_boundary_conditions();
//template void CFDv0_solver<double, 2, 3> :: set_boundary_conditions();
//template void CFDv0_solver<double, 2, 4> :: set_boundary_conditions();
//
//// Exchange ghost cells
//template void CFDv0_solver<float , 2, 3> :: exchange_ghost_cells( MPI_env &mpi_env );
//template void CFDv0_solver<float , 2, 4> :: exchange_ghost_cells( MPI_env &mpi_env );
//template void CFDv0_solver<double, 2, 3> :: exchange_ghost_cells( MPI_env &mpi_env );
//template void CFDv0_solver<double, 2, 4> :: exchange_ghost_cells( MPI_env &mpi_env );
//
//// Print submesh info (rho, vel, p)
//template void CFDv0_solver<float , 2, 3> :: gather_info( MPI_env &mpi_env, struct t_domain_info<float ,2> &domain_info );
//template void CFDv0_solver<float , 2, 4> :: gather_info( MPI_env &mpi_env, struct t_domain_info<float ,2> &domain_info );
//template void CFDv0_solver<double, 2, 3> :: gather_info( MPI_env &mpi_env, struct t_domain_info<double,2> &domain_info );
//template void CFDv0_solver<double, 2, 4> :: gather_info( MPI_env &mpi_env, struct t_domain_info<double,2> &domain_info );
//
//template void CFDv0_solver<float , 2, 3> :: allocate_ghost_cells( MPI_env &mpi_env, gmsh_mesh *mesh, std::vector<std::vector<t_solution_vars<float ,2>*>> &ghost_ptr, std::vector<std::vector<t_solution_vars<float ,2>*>> &ghost_bnd );
//template void CFDv0_solver<float , 2, 4> :: allocate_ghost_cells( MPI_env &mpi_env, gmsh_mesh *mesh, std::vector<std::vector<t_solution_vars<float ,2>*>> &ghost_ptr, std::vector<std::vector<t_solution_vars<float ,2>*>> &ghost_bnd );
//template void CFDv0_solver<double, 2, 3> :: allocate_ghost_cells( MPI_env &mpi_env, gmsh_mesh *mesh, std::vector<std::vector<t_solution_vars<double,2>*>> &ghost_ptr, std::vector<std::vector<t_solution_vars<double,2>*>> &ghost_bnd );
//template void CFDv0_solver<double, 2, 4> :: allocate_ghost_cells( MPI_env &mpi_env, gmsh_mesh *mesh, std::vector<std::vector<t_solution_vars<double,2>*>> &ghost_ptr, std::vector<std::vector<t_solution_vars<double,2>*>> &ghost_bnd );
//
//template void CFDv0_solver<float , 2, 3> :: assign_pointers( gmsh_mesh *mesh, std::vector<std::vector<t_solution_vars<float ,2>*>> cells_addr, std::vector<std::vector<t_solution_vars<float ,2>*>> ghost_addr, std::vector<std::vector<t_solution_vars<float ,2>*>> ghost_bnd_addr );
//template void CFDv0_solver<float , 2, 4> :: assign_pointers( gmsh_mesh *mesh, std::vector<std::vector<t_solution_vars<float ,2>*>> cells_addr, std::vector<std::vector<t_solution_vars<float ,2>*>> ghost_addr, std::vector<std::vector<t_solution_vars<float ,2>*>> ghost_bnd_addr );
//template void CFDv0_solver<double, 2, 3> :: assign_pointers( gmsh_mesh *mesh, std::vector<std::vector<t_solution_vars<double,2>*>> cells_addr, std::vector<std::vector<t_solution_vars<double,2>*>> ghost_addr, std::vector<std::vector<t_solution_vars<double,2>*>> ghost_bnd_addr );
//template void CFDv0_solver<double, 2, 4> :: assign_pointers( gmsh_mesh *mesh, std::vector<std::vector<t_solution_vars<double,2>*>> cells_addr, std::vector<std::vector<t_solution_vars<double,2>*>> ghost_addr, std::vector<std::vector<t_solution_vars<double,2>*>> ghost_bnd_addr );
//
//
//template void CFDv0_solver<float , 2, 3> :: output_cfd( float  time, std::string filename, gmsh_mesh *mesh );
//template void CFDv0_solver<float , 2, 4> :: output_cfd( float  time, std::string filename, gmsh_mesh *mesh );
//template void CFDv0_solver<double, 2, 3> :: output_cfd( double time, std::string filename, gmsh_mesh *mesh );
//template void CFDv0_solver<double, 2, 4> :: output_cfd( double time, std::string filename, gmsh_mesh *mesh );
//
//template float  CFDv0_solver<float , 2, 3> :: compute_dt( float  cflMax );
//template float  CFDv0_solver<float , 2, 4> :: compute_dt( float  cflMax );
//template double CFDv0_solver<double, 2, 3> :: compute_dt( double cflMax );
//template double CFDv0_solver<double, 2, 4> :: compute_dt( double cflMax );
//
//template void CFDv0_solver<float , 2, 3> :: halo_comm_wait( MPI_env &mpi_env );
//template void CFDv0_solver<float , 2, 4> :: halo_comm_wait( MPI_env &mpi_env );
//template void CFDv0_solver<double, 2, 3> :: halo_comm_wait( MPI_env &mpi_env );
//template void CFDv0_solver<double, 2, 4> :: halo_comm_wait( MPI_env &mpi_env );

#endif
