#include <iostream>
#include <iomanip>
#include <math.h>
#include "mpi.h"

#include "api/hpath.h"
#include "api/elements.h"
#include "api/mesh_writer.h"

using namespace std;

const int PI = 4*atan(1.);

//***************************************************************************************************
// Find boundary Hpath
//***************************************************************************************************
void Hpath::get_hpath_boundary( gmsh_mesh* mesh ){

	int myrank;
	MPI_Comm_rank( MPI_COMM_WORLD, &myrank );

	std::vector<int>  tmp_hpath ( mesh->cells.size(), -1 );
	std::vector<bool> is_deadend( mesh->cells.size(), false );

	//std::vector<bool> bnd_neigh ( mesh->cells.size(), false );
	std::vector<int> bnd_neigh ( mesh->cells.size(), -1 );

	// Identify boundary cells and their neighbors
	for( size_t icell=0; icell < mesh->cells.size(); icell++ ){

		if( mesh->cells[icell].is_bcell == true ){

			for( int iface = 0; iface < mesh->faces_per_cell; iface++ ){
				if( mesh->cell2neigh[icell][iface].type != CELL_REGULAR ){
					continue;
				}
				int adjc_id = (int) mesh->cell2neigh[icell][iface].id ;

				//bnd_neigh[adjc_id] = true;
				bnd_neigh[adjc_id] = icell;
			}
		}
	}

	// Find starting cell
	const int NUM_STARTING_CELLS = 10;
	int num_starting_cells = 0;
	std::vector<int> starting_cells( NUM_STARTING_CELLS );

	for( size_t icell=0; icell < mesh->cells.size(); icell++ ){
		int found = false;

		for( int iface=0; iface < mesh->faces_per_cell; iface++ ){
			t_gmsh_neighbor *neigh = &( mesh->cell2neigh[icell][iface] );

			if( neigh->type != CELL_REGULAR ){
				starting_cells[num_starting_cells++] = icell;
				found = true;
				break;
			}

		}
		if( found ){
			break;
		}
	}

	// Special case for quads
	if( mesh->faces_per_cell == 4 ){
		for( size_t icell=0; icell < mesh->cells.size(); icell++ ){
			int found = false;

			for( int iface=0; iface < mesh->faces_per_cell; iface++ ){
				t_gmsh_neighbor *neigh = &( mesh->cell2neigh[icell][iface] );

				if( neigh->type != CELL_REGULAR ){

					int num_neighs_regular = 0;
					for( int jface=0; jface < mesh->faces_per_cell; jface++ ){
						t_gmsh_neighbor *neigh2 = &( mesh->cell2neigh[icell][jface] );

						if( neigh2->type == CELL_REGULAR ){
							if( mesh->cells[neigh2->id].is_bcell ){
								num_neighs_regular++;
							}
						}
					}

					if( num_neighs_regular > 1 ){
						starting_cells[num_starting_cells++] = icell;
						found = true;
						break;
					}
				}

			}
			if( found && num_starting_cells >= NUM_STARTING_CELLS ){
				break;
			}else{
				found = false;
			}
		}
	}

	if( starting_cells.size() == 0 ){
		cout << "Could not find any good starting cells for the boundary Hpath on rank " << myrank << ".";
		cout << "Terminating simulation prematurely." << endl;
		MPI_Abort( MPI_COMM_WORLD, 923 );
	}


	// Try different starting cells and find boundary Hpath
	MPI_Barrier( MPI_COMM_WORLD );

	bool success;
	for( size_t icell=0; icell < starting_cells.size(); icell++ ){
		success = find_boundary_hpath( mesh, starting_cells[icell], bnd_neigh );

		if( success ) break;
	}

	if( !success ){
		cout << "Failed boundary Hpath on rank " << myrank << ". Terminating simulation prematurely." << endl;
		MPI_Abort( MPI_COMM_WORLD, 924 );
	}
}

//***************************************************************************************************
bool Hpath::find_boundary_hpath( gmsh_mesh *mesh, int start_cell, std::vector<int> bnd_neigh ){

	int myrank;
	MPI_Comm_rank( MPI_COMM_WORLD, &myrank );

	std::vector<int>  tmp_hpath ( mesh->cells.size(), -1 );
	std::vector<bool> is_deadend( mesh->cells.size(), false );

	// Keep track of previously walked cells
	std::vector<bool> walked_cells( mesh->cells.size(), false );
	walked_cells[start_cell] = true;

	// Store first cell in hpath
	int num_bnd_cells = 0;
	tmp_hpath[num_bnd_cells++] = start_cell;

	int next_cell = get_next_bnd_hpath_cell( mesh, -1, start_cell, -1, walked_cells, bnd_neigh );
	int last_cell = -1;

	for( int iface=0; iface < mesh->faces_per_cell; iface++ ){
		int adjc_id = mesh->cell2neigh[start_cell][iface].id ;

		if( mesh->cell2neigh[start_cell][iface].type != CELL_REGULAR ){
			continue;
		}

		if( mesh->faces_per_cell == 3 ){
			if( adjc_id != next_cell ){
				last_cell = adjc_id;
			}
		}else if( mesh->faces_per_cell == 4 ){
			if( mesh->cells[adjc_id].is_bcell ){
				last_cell = adjc_id;
			}
		}
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Walk Hamiltonian path
	int this_cell = next_cell;
	int deadend_count  = 0;
	int nonhpath_cells = 0;
	std::vector<int>  tmp_deadend; tmp_deadend.reserve( mesh->cells.size() );

	do{
		// Store next cell
		tmp_hpath[num_bnd_cells++] = next_cell;

		if( next_cell < 0 || next_cell > (int) walked_cells.size() ){
			//for( int icell=0; icell < num_bnd_cells; icell++ ){
			//	cout << "hpath =" << setw(10) << icell
			//					  << setw(10) << tmp_hpath[icell]
			//					  << endl;
			//}
			return false;
			MPI_Abort( MPI_COMM_WORLD, 9034 );
		}
		walked_cells[next_cell] = true;

		// Get next cell
		this_cell = next_cell;

		int deadend_neighbor = get_deadend_neighbor_bnd( mesh, this_cell, last_cell, walked_cells );

		if( deadend_neighbor > -1 ){
			nonhpath_cells++;
			tmp_deadend.push_back( deadend_neighbor );
			deadend_count++;
			is_deadend[deadend_neighbor] = true;
		}

		next_cell = get_next_bnd_hpath_cell( mesh, start_cell, this_cell, deadend_neighbor, walked_cells, bnd_neigh);

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// If path fails prematurely, most likely reason is long strip of deadend cells
		if( next_cell == -1 && num_bnd_cells + nonhpath_cells < (int) mesh->cells.size() ){
			// Walk back Hpath out of deadend and change direction
			//cout << "-------------> failed hpath, rank =" << setw(10) << myrank << endl;

			bool found_exit = false;
			int walkedback_cells = 1;
			std::vector<bool> walked_cells_restart = walked_cells;
			int restart_next_cell = -1;

			// 1- Starting from last walked cell, walk back one cell
			// 2- set last walked cell as a deadend neighbor
			// 3- look for the next cell considering deadend neighbor
			// 4- If another path is found, take it
			//    else, go back to 1
			do{
				if( num_bnd_cells - walkedback_cells < 0 )
					break;
				int deadend_neighbor  = tmp_hpath[num_bnd_cells-walkedback_cells  ];
				int restart_this_cell = tmp_hpath[num_bnd_cells-walkedback_cells-1];

				if( restart_this_cell == -1 ) break;

				restart_next_cell = get_next_hpath_cell( mesh, restart_this_cell, deadend_neighbor, walked_cells );

				if( restart_next_cell != -1 ){
					found_exit = true;
				}else{
					walkedback_cells++;
				}
				if( walkedback_cells == 10 ){
					break;
				}
			}while( !found_exit );

			if( restart_next_cell == -1 ){
				break;
			}

			for( int icell=0; icell < walkedback_cells; icell++ ){
				int new_deadend_cell = tmp_hpath[num_bnd_cells-icell-1];

				walked_cells[new_deadend_cell] = false;
				is_deadend  [new_deadend_cell] = true;
				tmp_deadend.push_back(new_deadend_cell);
				tmp_hpath[num_bnd_cells-icell-1] = -1;

				nonhpath_cells++;
			}
			next_cell = restart_next_cell;
			num_bnd_cells -= walkedback_cells;
		}
	} while( next_cell != last_cell );

	if( next_cell != -1 ) walked_cells[next_cell] = true;

	// Check if boundary cells are missing and if so, add to deadends
	for( size_t icell=0; icell < mesh->cells.size(); icell++ ){
		if( !mesh->cells[icell].is_bcell )
			continue;

		if( walked_cells[icell] || is_deadend[icell] )
			continue;

		// If we reach this point, it means one of the boundary cells was not "walked"
		nonhpath_cells++;
		tmp_deadend.push_back( icell );
		deadend_count++;
	}

	tmp_hpath[num_bnd_cells++] = next_cell;
	tmp_hpath.resize( num_bnd_cells );

	mesh->hpath = tmp_hpath;
	mesh->deadend_cells = tmp_deadend;

	return true;
}

//***************************************************************************************************
// Find interior Hpath without splitting cells
//***************************************************************************************************
void Hpath::get_hpath_interior( gmsh_mesh *mesh ){

	int myrank, nProcs;
	MPI_Comm_size( MPI_COMM_WORLD, &nProcs );
	MPI_Comm_rank( MPI_COMM_WORLD, &myrank );

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// No cells have been walked yes
	std::vector<bool> walked_cells( mesh->cells.size(), false );

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Get a list of starting cells
	std::vector<int> starting_cells = get_starting_cells( mesh );

	std::vector<int> deadend_cells;
	deadend_cells.reserve( mesh->cells.size() );

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Walk Hamiltonian path
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	std::vector<int> hpath( mesh->cells.size() );
	int ind_trial = 0;
	int max_trials = ( (int) starting_cells.size() / 1 );

	// Keep track of failed Hpaths
	int total_hpath_cells;
	int best_hpath_cells = -1;
	int best_hpath_trial = -1;
	//int best_hpath_algo  = -1;
	double best_hpath_prct  = -1;

	// Try simple algorithm first...
	do{
		find_hpath_nodeadend_simple ( mesh, starting_cells[ind_trial], total_hpath_cells, hpath, deadend_cells );

		if( (int) mesh->cells.size() != (int) (hpath.size() + deadend_cells.size()) ){
			cout << "Error: some deadend cells are missing" << endl;
			//MeshWriter mesh_writer;
			//mesh_writer.tecplot_hpath_nosplit( "output/failed_hpath.plt", mesh, hpath );
			MPI_Abort( MPI_COMM_WORLD, 922 );
		}
		if( best_hpath_cells < total_hpath_cells ){
			best_hpath_cells = total_hpath_cells;
			best_hpath_trial = ind_trial;
			//best_hpath_algo  = 0;
			best_hpath_prct  = (total_hpath_cells) / (double) mesh->cells.size() * 100;
		}

		if( best_hpath_prct > 99 )
			break;
		if( ind_trial > 20 && best_hpath_prct > 96 )
			break;
		if( ind_trial > 40 )
			break;
	}while( ++ind_trial < max_trials );

	//// If unsuccessful, try complex algo
	//if( !successful_hpath ){
	//	do{ find_hpath_nodeadend_complex ( mesh, starting_cells[ind_trial], hpath, deadend_cells);
	//	}while( !successful_hpath && ++ind_trial < max_trials );
	//}

	// Find final Hpath
	find_hpath_nodeadend_simple( mesh, starting_cells[best_hpath_trial], total_hpath_cells, hpath, deadend_cells );

	mesh->hpath         = hpath;
	mesh->hpath.resize( mesh->cells.size() - deadend_cells.size() );
	mesh->deadend_cells = deadend_cells;

	// Print Hpath statistics
	MPI_Barrier( MPI_COMM_WORLD );

	struct value_loc{
		double val;
		int   rank;
	};
	double percent_hpath = ((double) best_hpath_cells) / (double) mesh->cells.size() * 100;

	struct value_loc hpath_in, hpath_min, hpath_max;

	hpath_in.val  = percent_hpath;
	hpath_in.rank = myrank;

	MPI_Reduce( &hpath_in, &hpath_min, 1, MPI_DOUBLE_INT, MPI_MINLOC, 0, MPI_COMM_WORLD );
	MPI_Reduce( &hpath_in, &hpath_max, 1, MPI_DOUBLE_INT, MPI_MAXLOC, 0, MPI_COMM_WORLD );

	double hpath_avg;
	MPI_Reduce( &hpath_in.val, &hpath_avg, 1, MPI_DOUBLE, MPI_SUM   , 0, MPI_COMM_WORLD );

	double percent_hpath_min;
	double percent_hpath_max;
	MPI_Reduce( &percent_hpath, &percent_hpath_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD );
	MPI_Reduce( &percent_hpath, &percent_hpath_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD );

	if( myrank == 0 ){
		cout << "Best  Hpath =" << setw(10) << hpath_max.val << ", rank =" << setw(10) << hpath_max.rank << endl;
		cout << "Worst Hpath =" << setw(10) << hpath_min.val << ", rank =" << setw(10) << hpath_min.rank << endl;
		cout << "Avg   Hpath =" << setw(10) << hpath_avg / (double) nProcs << endl;
		cout << endl;
	}

}

//***************************************************************************************************
void Hpath:: find_hpath_nodeadend_simple( gmsh_mesh			*mesh,
										  int				starting_cell,
										  int				&total_hpath_cells,
										  std::vector<int>	&hpath,
										  std::vector<int> &deadend_cells ){

	int myrank;
	MPI_Comm_rank( MPI_COMM_WORLD, &myrank );

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Find all boundary cells and their regular neighbors
	std::vector<bool> cell_is_boundary   ( mesh->cells.size(), false );
	std::vector<bool> cell_is_boundneigh ( mesh->cells.size(), false );

	for( size_t icell=0; icell < mesh->cells.size(); icell++ ){
		for( int iface=0; iface < mesh->faces_per_cell; iface++ ){
			if( mesh->cell2neigh[icell][iface].type != CELL_REGULAR ){
				cell_is_boundary[icell] = true;
				break;
			}
		}
	}

	t_gmsh_neighbor *neigh;
	for( size_t icell=0; icell < mesh->cells.size(); icell++ ){
		for( int iface=0; iface < mesh->faces_per_cell; iface++ ){
			neigh = &( mesh->cell2neigh[icell][iface] );

			if( neigh->type != CELL_REGULAR ){
				continue;
			}

			if( cell_is_boundary[neigh->id] ){
				cell_is_boundneigh[icell] = true;
				break;
			}
		}
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Walk Hamiltonian path
	int m_count = 0;
	int this_cell;
	int next_cell = starting_cell;

	bool update_hpath_algorithm = false;

	std::vector<int>  tmp_hpath   ( mesh->cells.size(), 0 );
	std::vector<bool> walked_cells( mesh->cells.size(), false );
	std::vector<int>  tmp_deadend; tmp_deadend.reserve( mesh->cells.size() );

	//int  index_update_boundary  = -1;
	bool first_time = true;

	int nonhpath_cells = 0;

	// Spiral algorithm
	std::vector<int>  hpath_neigh        ( mesh->cells.size(), -1 );
	std::vector<int>  hpath_neighofneigh ( mesh->cells.size(), -1 );

	do{

		this_cell = next_cell;

		// Mark cell as "walked"
		walked_cells[this_cell] = true;

		// Update hpath
		tmp_hpath[m_count++] = this_cell;

		if( first_time && tmp_hpath[0] != starting_cell ){
			cout << "------------------>Weird thing happened" << endl;
			first_time = false;
			//exit(0);
		}


		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// If one of the neighboring cell is a deadend, split this_cell and deadend_cell in two
		// Only store information for now, the mesh and hpath will be updated at the end
		int deadend_neighbor = get_deadend_neighbor( mesh, this_cell, walked_cells );

		if( deadend_neighbor > -1 ){
			nonhpath_cells++;

			tmp_deadend.push_back( deadend_neighbor );
		}

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Get the best next cell
		if( update_hpath_algorithm ){
			next_cell = get_next_hpath_cell( mesh, this_cell, deadend_neighbor, walked_cells );

			//next_cell = get_next_hpath_cell_spiral( mesh, hpath_neigh, hpath_neighofneigh,
			//										this_cell, deadend_neighbor, walked_cells,
			//										update_hpath_algorithm );
		}else{
			next_cell = get_next_hpath_cell_boundary( mesh,
													cell_is_boundary,
													cell_is_boundneigh,
													this_cell,
													deadend_neighbor,
													walked_cells,
													update_hpath_algorithm );
			if( update_hpath_algorithm ){
				//index_update_boundary = m_count;
			}
		}

		//next_cell = get_next_hpath_cell_spiral( mesh, hpath_neigh, hpath_neighofneigh,
		//										this_cell, deadend_neighbor, walked_cells,
		//										update_hpath_algorithm );

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// If path fails prematurely, most likely reason is long strip of deadend cells
		if( next_cell == -1 && m_count + nonhpath_cells < (int) mesh->cells.size() ){
			// Walk back Hpath out of deadend and change direction
			//cout << "-------------> failed hpath" << endl;

			bool found_exit = false;
			int walkedback_cells = 1;
			std::vector<bool> walked_cells_restart = walked_cells;
			int restart_next_cell = -1;

			// 1- Starting from last walked cell, walk back one cell
			// 2- set last walked cell as a deadend neighbor
			// 3- look for the next cell considering deadend neighbor
			// 4- If another path is found, take it
			//    else, go back to 1
			do{
				if( m_count - walkedback_cells < 0 )
					break;
				int deadend_neighbor  = tmp_hpath[m_count-walkedback_cells  ];
				int restart_this_cell = tmp_hpath[m_count-walkedback_cells-1];

				if( restart_this_cell == -1 ) break;

				restart_next_cell = get_next_hpath_cell( mesh, restart_this_cell, deadend_neighbor, walked_cells );

				if( restart_next_cell != -1 ){
					found_exit = true;
				}else{
					walkedback_cells++;
				}

				if( walkedback_cells == 10 ){
					break;
				}
			}while( !found_exit );

			if( restart_next_cell == -1 ){
				break;
			}

			for( int icell=1; icell < walkedback_cells-1; icell++ ){
				int new_deadend_cell = tmp_hpath[m_count-icell];
				walked_cells[new_deadend_cell] = false;
				tmp_deadend.push_back(new_deadend_cell);
				tmp_hpath[m_count-icell] = -1;

				nonhpath_cells++;
			}
			next_cell = restart_next_cell;
			m_count -= walkedback_cells;
		}

		if( next_cell != -1 ){
			walked_cells[next_cell] = true;
		}
	}while( next_cell != -1 );

	tmp_hpath.resize( m_count );
	hpath = tmp_hpath;

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Set deadend cells:
	// 		- most of deadends were (hopefully) recorder along with the Hpath algo
	// 		- it's possible that large chunks of the mesh were left out of the Hpath ==> record
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	std::vector<bool> is_deadend( mesh->cells.size(), true );
	deadend_cells.clear();
	deadend_cells.reserve( mesh->cells.size() );

	// Define the deadend cells as cells not walked by Hpath
	for( size_t icell=0; icell < hpath.size(); icell++ ){
		if( hpath[icell] == -1 ){
			cout << "Error with Hpath algorithm." << endl;
			MPI_Abort( MPI_COMM_WORLD, 921 );
		}
		is_deadend[hpath[icell]] = false;
	}

	// Walk the Hpath, and record deadends along the way
	int m_deadends = 0;
	std::vector<bool> recorded_deadend( mesh->cells.size(), false );

	for( size_t icell=0; icell < hpath.size(); icell++ ){
		int gcell = hpath[icell];
		for( int iface=0; iface < mesh->faces_per_cell; iface++ ){
			t_gmsh_neighbor *neigh = &( mesh->cell2neigh[gcell][iface] );

			if( neigh->type != CELL_REGULAR )
				continue;

			if( is_deadend[neigh->id] && !recorded_deadend[neigh->id] ){
				deadend_cells.push_back( neigh->id );
				recorded_deadend[neigh->id] = true;
				m_deadends++;
			}
		}
	}

	// Check if there are deadend cells not yet recorded (i.e. with no Hpath cells neighbor )
	for( size_t icell=0; icell < mesh->cells.size(); icell++ ){
		if( is_deadend[icell] && !recorded_deadend[icell] ){
			deadend_cells.push_back( icell );
			recorded_deadend[icell] = true;
			m_deadends++;
		}
	}
	deadend_cells.resize( m_deadends );

	total_hpath_cells = m_count;

	if( mesh->cells.size() == hpath.size() + deadend_cells.size() ) {
		return;
	}
	MPI_Abort( MPI_COMM_WORLD, 934 );
}

////***************************************************************************************************
//void Hpath:: find_hpath_complex( gmsh_mesh		*mesh,
//								 int				 starting_cell,
//						 		 gmsh_mesh		    *mesh_deadend,
//						 		 std::vector<int>	&hpath,
//						 		 int				&deadend_count,
//						 		 bool				&successful_hpath,
//						 		 int                &total_hpath_cells ){
//
//	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//	// Find all boundary cells and their regular neighbors
//	int num_priority_cells = 0;
//	std::vector<bool> cell_is_boundary   ( mesh->cells.size(), false );
//	std::vector<bool> cell_is_boundneigh ( mesh->cells.size(), false );
//
//	for( size_t icell=0; icell < mesh->cells.size(); icell++ ){
//		for( int iface=0; iface < mesh->faces_per_cell; iface++ ){
//			if( mesh->cells[icell].neighbors[iface].type != CELL_REGULAR ){
//				cell_is_boundary[icell] = true;
//				num_priority_cells++;
//				break;
//			}
//		}
//	}
//
//	t_neighbor *neigh;
//	for( size_t icell=0; icell < mesh->cells.size(); icell++ ){
//		for( int iface=0; iface < mesh->faces_per_cell; iface++ ){
//			neigh = &( mesh->cells[icell].neighbors[iface] );
//
//			if( neigh->type != CELL_REGULAR ){
//				continue;
//			}
//
//			if( cell_is_boundary[neigh->id] ){
//				num_priority_cells++;
//				cell_is_boundneigh[icell] = true;
//				break;
//			}
//		}
//	}
//
//	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//	// Walk Hamiltonian path
//	int m_count = 0;
//	int count_detours = 0;
//	int prev_cell;
//	int this_cell;
//	int next_cell = starting_cell;
//
//	int  index_update_boundary = -1;
//	bool update_hpath_algorithm = false;
//
//	std::vector<int>  tmp_hpath   ( mesh->cells.size(), 0 );
//	std::vector<bool> walked_cells( mesh->cells.size(), false );
//
//	std::vector<int>  hpath_neigh        ( mesh->cells.size(), -1 );
//	std::vector<int>  hpath_neighofneigh ( mesh->cells.size(), -1 );
//
//	hpath_algorithm hpath_algorithm_m = ALGO_BOUNDARY;
//
//	int num_boundhpath_cells = 0;
//	int index_hpath_bound  = -1;
//	int index_hpath_spiral = -1;
//
//	std::vector<int> record_bound_hpath_neigh       ( mesh->cells.size(), -1 );
//	std::vector<int> record_bound_hpath_neighofneigh( mesh->cells.size(), -1 );
//
//	std::vector<int> record_spiral_hpath_neigh       ( mesh->cells.size(), -1 );
//	std::vector<int> record_spiral_hpath_neighofneigh( mesh->cells.size(), -1 );
//
//	std::vector<int> record_center_hpath_neigh       ( mesh->cells.size(), -1 );
//	std::vector<int> record_center_hpath_neighofneigh( mesh->cells.size(), -1 );
//
//	do{
//
//		prev_cell = this_cell;
//		this_cell = next_cell;
//
//		// Mark cell as "walked"
//		walked_cells[this_cell] = true;
//
//		// Update hpath
//		tmp_hpath[m_count++] = this_cell;
//
//		// If one of the neighboring cell is a deadend, split this_cell and deadend_cell in two
//		// Only store information for now, the mesh and hpath will be updated at the end
//		int deadend_neighbor = get_deadend_neighbor( mesh, this_cell, walked_cells );
//
//		if( deadend_neighbor > -1 ){
//			deadend_count++;
//
//			// Split deadend cell in two
//			split_deadend_neighbor( mesh, prev_cell, this_cell, deadend_count, deadend_neighbor,
//																mesh_deadend );
//		}
//
//		switch( hpath_algorithm_m ){
//			case ALGO_BOUNDARY:
//				next_cell = get_next_hpath_cell_boundary( mesh,
//													  	  cell_is_boundary,
//													  	  cell_is_boundneigh,
//													  	  this_cell,
//													  	  deadend_neighbor,
//													  	  walked_cells,
//													  	  update_hpath_algorithm );
//				if( update_hpath_algorithm ){
//					hpath_algorithm_m    = ALGO_TRUE_SPIRAL;
//
//					index_hpath_bound = m_count;
//					record_bound_hpath_neigh        = hpath_neigh;
//					record_bound_hpath_neighofneigh = hpath_neighofneigh;
//				}
//
//				num_boundhpath_cells++;
//
//				break;
//			case ALGO_TRUE_SPIRAL:
//				next_cell = get_next_hpath_cell_spiral( mesh, hpath_neigh, hpath_neighofneigh,
//														this_cell, deadend_neighbor, walked_cells,
//														update_hpath_algorithm );
//
//				if( m_count > num_boundhpath_cells*5 ){
//					hpath_algorithm_m = ALGO_MESH_CENTER;
//
//					index_hpath_spiral = m_count;
//					record_spiral_hpath_neigh        = hpath_neigh;
//					record_spiral_hpath_neighofneigh = hpath_neighofneigh;
//				}
//				break;
//			case ALGO_MESH_CENTER:
//				next_cell = get_next_hpath_cell( mesh, this_cell, deadend_neighbor, walked_cells );
//				break;
//		}
//
//		record_hpath_neighbors( mesh, m_count, next_cell, deadend_neighbor, walked_cells,
//								hpath_neigh, hpath_neighofneigh );
//
//		// If path fails prematurely, most likely reason is long strip of deadend cells
//
//		if( next_cell == -1 && m_count+deadend_count < (int) mesh->cells.size() && count_detours < 1 ){
//			count_detours++;
//			//make_detour_new( mesh, mesh_deadend, m_count, deadend_count, deadend_neighbor,
//			//					 index_update_boundary, update_hpath_algorithm, walked_cells, tmp_hpath,
//			//					 next_cell );
//
//			make_detour_spiral( mesh, mesh_deadend, m_count, deadend_count, deadend_neighbor,
//								 index_update_boundary, update_hpath_algorithm, walked_cells, tmp_hpath,
//								 next_cell,
//								 hpath_algorithm_m, hpath_neigh, hpath_neighofneigh,
//								 index_hpath_bound,  record_bound_hpath_neigh, record_bound_hpath_neighofneigh,
//								 index_hpath_spiral, record_spiral_hpath_neigh, record_spiral_hpath_neighofneigh);
//
//			// Update this_cell
//			this_cell = tmp_hpath[m_count-1];
//		}
//
//		if( next_cell != -1 ){
//			walked_cells[next_cell] = true;
//		}
//	}while( next_cell != -1 );
//
//	hpath = tmp_hpath;
//	successful_hpath = ( m_count + deadend_count == (int) mesh->cells.size() ) ? true : false;
//	total_hpath_cells = m_count + deadend_count;
//}

//***************************************************************************************************
// Record the cells neighboring the hpath (and their neighbors) for the true-spiral algorithm
//***************************************************************************************************
void Hpath:: record_hpath_neighbors( gmsh_mesh *mesh,
									 int     m_count,
									 int     next_cell,
									 int     deadend_neighbor,
									 std::vector<bool>  walked_cells,
									 std::vector<int>  &hpath_neigh,
									 std::vector<int>  &hpath_neighofneigh ){


	if( next_cell == -1 ){
		return;
	}

	for( int iface=0; iface < mesh->faces_per_cell; iface++ ){
		t_gmsh_neighbor *neigh = &( mesh->cell2neigh[next_cell][iface] );

		// Skip walked cells, out-of-bounds cells and deadends
		if( walked_cells[neigh->id] ){
			continue;
		}
		if( neigh->type != CELL_REGULAR ){
			continue;
		}
		if( neigh->id == deadend_neighbor ){
			continue;
		}
		if( hpath_neigh[neigh->id] > -1 ){
			continue;
		}

		hpath_neigh[neigh->id] = m_count;

		for( int jface=0; jface < mesh->faces_per_cell; jface++ ){
			t_gmsh_neighbor *neigh_of_neigh = &( mesh->cell2neigh[neigh->id][jface] );

			// Skip walked cells, out-of-bounds cells and deadends
			if( walked_cells[neigh_of_neigh->id] ){
				continue;
			}
			if( neigh_of_neigh->type != CELL_REGULAR ){
				continue;
			}
			if( neigh_of_neigh->id == deadend_neighbor ){
				continue;
			}
			if( neigh_of_neigh->id == next_cell ){
				continue;
			}
			if( hpath_neighofneigh[neigh_of_neigh->id] > -1 ){
				continue;
			}
			hpath_neighofneigh[neigh_of_neigh->id] = m_count;
		}
	}
}

//***************************************************************************************************
// First part of hpath algorithm, go through all the boundary cells
//***************************************************************************************************
int Hpath::get_next_hpath_cell_boundary( gmsh_mesh* mesh,
									     std::vector<bool>  cell_is_boundary,
									     std::vector<bool>  cell_is_boundneigh,
									     int current_cell,
									     int deadend_neighbor,
									     std::vector<bool> walked_cells,
									     bool &update_hpath_algorithm ){

	t_gmsh_neighbor *neigh;

	// Priority 1: boundary cells
	for( int iface=0; iface < mesh->faces_per_cell; iface++ ){
		neigh = &( mesh->cell2neigh[current_cell][iface] );

		// Skip out-of-bounds cells
		if( neigh->type != CELL_REGULAR ){
			continue;
		}

		// Skip already walked cells
		if( walked_cells[neigh->id] ){
			continue;
		}

		// Skip deadend
		if( neigh->id == deadend_neighbor ){
			continue;
		}

		if( cell_is_boundary[neigh->id] && !walked_cells[neigh->id] ){
			return neigh->id;
		}
	}

	// Priority 2: boundary's neighbors
	for( int iface=0; iface < mesh->faces_per_cell; iface++ ){
		neigh = &( mesh->cell2neigh[current_cell][iface] );

		// Skip out-of-bounds cells
		if( neigh->type != CELL_REGULAR ){
			continue;
		}

		// Skip already walked cells
		if( walked_cells[neigh->id] ){
			continue;
		}

		// Skip deadend
		if( neigh->id == deadend_neighbor ){
			continue;
		}

		if( cell_is_boundneigh[neigh->id] && !walked_cells[neigh->id] ){
			return neigh->id;
		}
	}

	// Couldn't find a boundary or its neighbor, happens sometimes on concave meshes
	// 		-----> look for the neighbor's neighbors
	for( int iface=0; iface < mesh->faces_per_cell; iface++ ){
		neigh = &( mesh->cell2neigh[current_cell][iface] );

		// Skip already walked cells
		if( walked_cells[neigh->id] ){
			continue;
		}

		// Skip non-regular cells
		if( neigh->type != CELL_REGULAR ){
			continue;
		}

		// Skip deadend
		if( neigh->id == deadend_neighbor ){
			continue;
		}

		for( int jface=0; jface < mesh->faces_per_cell; jface++ ){
			t_gmsh_neighbor *neigh_of_neigh = &( mesh->cell2neigh[neigh->id][jface] );


			// Skip already walked cells
			if( walked_cells[neigh_of_neigh->id] ){
				continue;
			}

			// Skip non-regular cells
			if( neigh_of_neigh->type != CELL_REGULAR ){
				continue;
			}

			// Skip deadend
			if( neigh_of_neigh->id == deadend_neighbor ){
				continue;
			}

			if( cell_is_boundneigh[neigh_of_neigh->id] ){
				return neigh->id;
			}
		}

	}

	// If we've reached this point, the boundary cells and their neighbors are probably all walked
	for( int iface=0; iface < mesh->faces_per_cell; iface++ ){

		neigh = &( mesh->cell2neigh[current_cell][iface] );

		// Skip already walked cells
		if( walked_cells[neigh->id] ){
			continue;
		}

		// Skip non-regular cells
		if( neigh->type != CELL_REGULAR ){
			continue;
		}

		// Skip deadend
		if( neigh->id == deadend_neighbor ){
			continue;
		}

		update_hpath_algorithm = true;
		return neigh->id;
	}

	//int myrank;
	//MPI_Comm_rank( MPI_COMM_WORLD, &myrank );
	//cout << "Error while finding boundary hpath for partition:" << setw(10) << myrank << setw(10) << mesh->part_id << endl;
	//exit(0);
	return -1;
}

//***************************************************************************************************
// Second part of hpath algorithm: start with a "true" spiral.
// 		- regular spiral-like algorithm works fine except for concave/convex shapes
// 		- this algorithm does not use the mesh center
// 		- gives priority to cells with the oldest hpath cell as a neighbor
//***************************************************************************************************
int Hpath::get_next_hpath_cell_spiral( gmsh_mesh* mesh,
									   std::vector<int> hpath_neigh,
									   std::vector<int> hpath_neighofneigh,
									   int current_cell,
									   int deadend_neighbor,
									   std::vector<bool> walked_cells,
									   bool &update_hpath_algorithm ){

	int hpath_index = -2;
	int cell_id;

	// First check for cells neighboring an Hpath
	for( int iface=0; iface < mesh->faces_per_cell; iface++ ){
		t_gmsh_neighbor *neigh = &( mesh->cell2neigh[current_cell][iface] );

		// Skip out-of-bounds cells
		if( neigh->type != CELL_REGULAR ){
			continue;
		}

		// Skip already walked cells
		if( walked_cells[neigh->id] ){
			continue;
		}

		// Skip deadend
		if( neigh->id == deadend_neighbor ){
			continue;
		}

		if( hpath_index == -2 ){
			hpath_index = hpath_neigh[neigh->id];
			cell_id = neigh->id;
		}

		if( hpath_index < hpath_neigh[neigh->id] ){
			hpath_index = hpath_neigh[neigh->id];
			cell_id = neigh->id;
		}

		if( hpath_index < hpath_neighofneigh[neigh->id] ){
			hpath_index = hpath_neighofneigh[neigh->id];
			cell_id = neigh->id;
		}

		for( int jface=0; jface < mesh->faces_per_cell; jface++ ){
			t_gmsh_neighbor *neigh_of_neigh = &( mesh->cell2neigh[neigh->id][jface] );

			// Skip out-of-bounds cells
			if( neigh_of_neigh->type != CELL_REGULAR ){
				continue;
			}

			// Skip already walked cells
			if( walked_cells[neigh_of_neigh->id] ){
				continue;
			}

			// Skip deadend
			if( neigh_of_neigh->id == deadend_neighbor ){
				continue;
			}

			if( hpath_neigh[neigh_of_neigh->id] == -1 ){
				continue;
			}

			if( hpath_index < hpath_neigh[neigh_of_neigh->id] ){
				hpath_index = hpath_neigh[neigh_of_neigh->id];
				cell_id = neigh->id;
			}

			if( hpath_index < hpath_neighofneigh[neigh_of_neigh->id] ){
				hpath_index = hpath_neighofneigh[neigh_of_neigh->id];
				cell_id = neigh->id;
			}

		}
	}

	if( hpath_index != -2 ){
		return cell_id;
	}else{
		return -1;
	}
}

//***************************************************************************************************
int Hpath::get_next_hpath_cell( gmsh_mesh* mesh, int current_cell, int deadend_neighbor,
								std::vector<bool> walked_cells){

	int    next_cell = -1;
	double distance, max_distance = 0;

	// Compute the number of available cells
	std::vector<int> free_cells;
	free_cells.reserve( mesh->faces_per_cell );

	for( int iface=0; iface < mesh->faces_per_cell; iface++ ){

		t_gmsh_neighbor *neigh = &( mesh->cell2neigh[current_cell][iface] );

		// Skip cells not on this submesh
		if( neigh->type != CELL_REGULAR ){
			continue;
		}
		// Skip cells already walked
		if( walked_cells[neigh->id] ){
			continue;
		}
		free_cells.push_back( neigh->id );
	}

	// If only 1 cell is free, no further work required
	if( free_cells.size() == 1 ){
		return free_cells[0];
	}

	// If no cells are available, we've reached the end of the Hpath
	if( free_cells.size()== 0 ){
		return -1;
	}

	// If multiple free cells...
	for( size_t icell=0; icell < free_cells.size(); icell++ ){

		gmsh_cell *this_cell = &( mesh->cells[free_cells[icell]] );

		// Skip deadend neighbor: will be added after Hpath is done
		if( this_cell->id == deadend_neighbor ){
			continue;
		}

		distance = sqrt( pow( this_cell->xc[0]-mesh->x[0],2) + pow(this_cell->xc[1]-mesh->x[1],2));

		if( max_distance < distance ){
			max_distance = distance;
			next_cell    = this_cell->id;
		}
	}
	return next_cell;
}

//***************************************************************************************************
// Requirements for starting cell:
//		- cell's neighbors cannot be a irregular cell's neighbor (irregular cell ==> cell with 2
//					non-regular neighbors )
//		- cell's two neighbors must have at least two regular neighbors
//		- cell's two neighbors must have at least one non-boundary cell neighbor
//***************************************************************************************************
std::vector<int> Hpath:: get_starting_cells( gmsh_mesh *mesh ){

	// Find boundary cells
	int num_boundary_cells = 0;
	std::vector<int>  bound_cells( mesh->cells.size(), -1 );
	std::vector<bool> is_boundary( mesh->cells.size(), false );

	for( size_t icell=0; icell < mesh->cells.size(); icell++ ){
		for( int iface=0; iface < mesh->faces_per_cell; iface++ ){
			if( mesh->cell2neigh[icell][iface].type != CELL_REGULAR ){
				bound_cells[num_boundary_cells++ ] = mesh->cells[icell].id;
				is_boundary[mesh->cells[icell].id] = true;
				break;
			}
		}
	}

	// Find irregular cells (cells with only one regular neighbor)
	std::vector<int> irregular_neighs;
	std::vector<bool> is_irregular_cell_neigh( mesh->cells.size(), false );
	int tmp_irregular_neigh;
	for( int bcell=0; bcell < num_boundary_cells; bcell++ ){

		int num_regular_neighbors = 0;
		for( int iface=0; iface < mesh->faces_per_cell; iface++ ){
			t_gmsh_neighbor *neigh = &( mesh->cell2neigh[bound_cells[bcell]][iface] );

			if( neigh->type != CELL_REGULAR ){
				continue;
			}
			num_regular_neighbors++;

			tmp_irregular_neigh = neigh->id;
		}
		if( num_regular_neighbors < 2 ){
			irregular_neighs.push_back( tmp_irregular_neigh );
			is_irregular_cell_neigh[tmp_irregular_neigh] = true;
		}
	}

	// Find cell next to and not part of boundary hpath
	int num_starting_cells = 0;
	std::vector<int> tmp_starting_cells( mesh->cells.size(), -1 );

	for( int bcell=0; bcell < num_boundary_cells; bcell++ ){

		gmsh_cell *this_cell = &( mesh->cells[bound_cells[bcell]] );

		// Condition 0: no neighbor can be an irregular cell's neighbor
		bool valid_starting_cell = true;
		for( int iface=0; iface < mesh->faces_per_cell; iface++ ){
			t_gmsh_neighbor *neigh = &( mesh->cell2neigh[bound_cells[bcell]][iface] );

			if( neigh->type != CELL_REGULAR ){
				continue;
			}

			if( is_irregular_cell_neigh[neigh->id] ){
				valid_starting_cell = false;
			}
		}

		if( !valid_starting_cell ){
			continue;
		}

		// Condition 1: the two neighbors must have at least two regular cells other than the current
		// 				cell
		int num_valid_neighbors = 0;
		for( int iface=0; iface < mesh->faces_per_cell; iface++ ){
			t_gmsh_neighbor *neigh = &( mesh->cell2neigh[bound_cells[bcell]][iface] );

			if( neigh->type != CELL_REGULAR ){
				continue;
			}

			int num_regular_cells = 0;
			for( int jface=0; jface < mesh->faces_per_cell; jface++ ){
				t_gmsh_neighbor *neigh_of_neigh = &( mesh->cell2neigh[neigh->id][jface] );

				// Skip non-regular cells
				if( neigh_of_neigh->type != CELL_REGULAR ){
					continue;
				}

				// Skip current cell
				if( neigh_of_neigh->id == this_cell->id ){
					continue;
				}

				num_regular_cells++;
			}
			if( num_regular_cells > 1 ){
				num_valid_neighbors++;
			}
		}

		if( num_valid_neighbors < 2 ){
			continue;
		}

		// Condition 2: the two neighbors must have at least one interior cell
		num_valid_neighbors = 0;
		for( int iface=0; iface < mesh->faces_per_cell; iface++ ){
			t_gmsh_neighbor *neigh = &( mesh->cell2neigh[bound_cells[bcell]][iface] );

			if( neigh->type != CELL_REGULAR ){
				continue;
			}

			bool interior_cell = false;
			for( int jface=0; jface < mesh->faces_per_cell; jface++ ){
				t_gmsh_neighbor *neigh_of_neigh = &( mesh->cell2neigh[neigh->id][jface] );

				if( !is_boundary[neigh_of_neigh->id] ){
					interior_cell = true;
					break;
				}
			}

			if( interior_cell ){
				num_valid_neighbors++;
			}
		}

		if( num_valid_neighbors >= 2 ){
			tmp_starting_cells[num_starting_cells++] = this_cell->id;
		}
	}

	// Return final vector
	std::vector<int> starting_cells( num_starting_cells, -1 );

	for( size_t icell=0; icell < starting_cells.size(); icell++ ){
		starting_cells[icell] = tmp_starting_cells[icell];
	}

	if( starting_cells.size() == 0 ){
		cout << "Could not find any good starting cells. Terminating simulation prematurely." << endl;
		MPI_Abort( MPI_COMM_WORLD, 925 );
	}

	return starting_cells;

}

//***************************************************************************************************
int Hpath::get_deadend_neighbor( gmsh_mesh *mesh,
								 int current_cell,
								 std::vector<bool> walked_cells ){

	int neigh_cell;
	bool neigh_deadend = false;

	// If only one path is available, no choice anyway
	int cells_left = mesh->faces_per_cell;
	for( int iface=0; iface < mesh->faces_per_cell; iface++ ){

		if( mesh->cell2neigh[current_cell][iface].type != CELL_REGULAR ){
			cells_left--;
			continue;
		}

		neigh_cell = mesh->cell2neigh[current_cell][iface].id;

		if( walked_cells[neigh_cell] ){
			cells_left--;
		}
	}
	if( cells_left == 1 ){
		return -1;
	}

	// If more than one path is available, check if one of them is a deadend
	for( int iface=0; iface < mesh->faces_per_cell; iface++ ){

		// Skip out of bounds cells
		if( mesh->cell2neigh[current_cell][iface].type != CELL_REGULAR ){
			continue;
		}

		neigh_cell = mesh->cell2neigh[current_cell][iface].id;

		// Skip walked cells
		if( walked_cells[neigh_cell] ){
			continue;
		}

		// Check if cell is deadend
		neigh_deadend = is_neighbor_deadend( mesh, current_cell, neigh_cell, walked_cells );

		if( neigh_deadend ){
			break;
		}
	}

	return neigh_deadend ? neigh_cell : -1 ;
}

//***************************************************************************************************
int Hpath::get_deadend_neighbor_bnd( gmsh_mesh *mesh,
									 int current_cell,
									 int last_cell,
									 std::vector<bool> walked_cells ){

	int myrank;
	MPI_Comm_rank( MPI_COMM_WORLD, &myrank );

	int neigh_cell;
	bool neigh_deadend = false;

	// If only one path is available, no choice anyway
	int cells_left = mesh->faces_per_cell;
	for( int iface=0; iface < mesh->faces_per_cell; iface++ ){

		if( mesh->cell2neigh[current_cell][iface].type != CELL_REGULAR ){
			cells_left--;
			continue;
		}

		neigh_cell = mesh->cell2neigh[current_cell][iface].id;

		if( walked_cells[neigh_cell] ){
			cells_left--;
		}
	}
	if( cells_left == 1 ){
		return -1;
	}

	// If more than one path is available, check if one of them is a deadend
	for( int iface=0; iface < mesh->faces_per_cell; iface++ ){

		// Skip out of bounds cells
		if( mesh->cell2neigh[current_cell][iface].type != CELL_REGULAR ){
			continue;
		}

		neigh_cell = mesh->cell2neigh[current_cell][iface].id;

		// Skip walked cells
		if( walked_cells[neigh_cell] ){
			continue;
		}

		// Skip the last cell: not a deadend
		if( neigh_cell == last_cell ){
			continue;
		}

		// Check if cell is deadend
		neigh_deadend = is_neighbor_deadend( mesh, current_cell, neigh_cell, walked_cells );

		if( neigh_deadend ){
			break;
		}
	}

	return neigh_deadend ? neigh_cell : -1 ;
}

//***************************************************************************************************
bool Hpath::is_neighbor_deadend( gmsh_mesh *mesh, int current_cell, int adjc_cell,
								 std::vector<bool> walked_cells ){

	// Check if cell is a deadend
	int gcell;
	bool deadend = true;

	for( int jface=0; jface < mesh->faces_per_cell; jface++ ){
		gcell = mesh->cell2neigh[adjc_cell][jface].id;

		// Skip current cell
		if( gcell == current_cell ){
			continue;
		}
		// Skip out-of-bounds cell
		if( mesh->cell2neigh[adjc_cell][jface].type != CELL_REGULAR ){
			continue;
		}

		if( !walked_cells[gcell] ){
			deadend = false;
			break;
		}
	}

	return deadend;
}

//***************************************************************************************************
int Hpath::get_next_bnd_hpath_cell( gmsh_mesh        *mesh,
									int               starting_cell,
									int               current_cell,
									int               deadend_neighbor,
									std::vector<bool> walked_cells,
									std::vector<int>  bnd_neigh ){

	// Init vars
	int adjc_id;
	double distance, max_distance = 0;
	int next_cell = -1;
	int myrank;
	MPI_Comm_rank( MPI_COMM_WORLD, &myrank );

	std::vector<int> free_cells;

	// Count the number of available cells
	for( int iface=0; iface < mesh->faces_per_cell; iface++ ){

		if( mesh->cell2neigh[current_cell][iface].type != CELL_REGULAR ){
			continue;
		}

		adjc_id = mesh->cell2neigh[current_cell][iface].id ;

		if( walked_cells[adjc_id] )
			continue;
		if( adjc_id == deadend_neighbor )
			continue;

		free_cells.push_back( adjc_id );
	}

	// If only one cell is left, no choice anyway
	if( (int) free_cells.size() == 1 ){
		return free_cells[0];
	}

	// This point shouldn't be reached
	if( (int) free_cells.size() == 0 ){
		return -1;
		cout << "Error in generating boundary hpath. rank =" << setw(10) << myrank << endl;
		MPI_Abort( MPI_COMM_WORLD, 15 );
		//exit(0);
	}

	// If multiple cells are left, look for boundary cells
	for( int iface=0; iface < mesh->faces_per_cell; iface++ ){
		if( mesh->cell2neigh[current_cell][iface].type != CELL_REGULAR ){
			continue;
		}
		adjc_id = mesh->cell2neigh[current_cell][iface].id;

		if( walked_cells[adjc_id] )
			continue;
		if( adjc_id == deadend_neighbor )
			continue;
		if( mesh->cells[adjc_id].is_bcell )
			return adjc_id;
	}

	// If multiple cells are left, and no boundary cells, look for boundary cells' neighbors
	for( int iface=0; iface < mesh->faces_per_cell; iface++ ){
		if( mesh->cell2neigh[current_cell][iface].type != CELL_REGULAR ){
			continue;
		}
		adjc_id = mesh->cell2neigh[current_cell][iface].id;

		if( walked_cells[adjc_id] ){
			continue;
		}
		if( adjc_id == deadend_neighbor ){
			continue;
		}

		if( bnd_neigh[adjc_id] != -1 ){
			bool found_neigh = false;
			for( int jface=0; jface < mesh->faces_per_cell; jface++ ){
				int neigh_of_neigh = mesh->cell2neigh[adjc_id][jface].id;

				if( mesh->cells[neigh_of_neigh].is_bcell && !walked_cells[neigh_of_neigh] ){
					found_neigh = true;
					break;
				}
			}
			if( found_neigh ){
				return adjc_id;
			}
		}
	}

	// Check cells' neighbors' neighbors until a boundary cell is found
	for( int iface=0; iface < mesh->faces_per_cell; iface++ ){

		// Neighbor
		t_gmsh_neighbor *neigh = &( mesh->cell2neigh[current_cell][iface] );

		if( neigh->type != CELL_REGULAR )
			continue;
		if( walked_cells[neigh->id] )
			continue;
		if( neigh->id == deadend_neighbor )
			continue;

		for( int jface=0; jface < mesh->faces_per_cell; jface++ ){

			// Neighbor of neighbor
			t_gmsh_neighbor *neigh_of_neigh = &( mesh->cell2neigh[neigh->id][jface] );

			if( neigh_of_neigh->type != CELL_REGULAR )
				continue;
			if( neigh_of_neigh->id == current_cell )
				continue;
			if( walked_cells[neigh_of_neigh->id] && neigh_of_neigh->id != starting_cell )
				continue;

			if( mesh->cells[neigh_of_neigh->id].is_bcell ){
				return neigh->id;
			}
		}
	}

	// Check cells' neighbors' neighbors until a boundary cell is found
	for( int iface=0; iface < mesh->faces_per_cell; iface++ ){

		// Neighbor
		t_gmsh_neighbor *neigh = &( mesh->cell2neigh[current_cell][iface] );

		if( neigh->type != CELL_REGULAR )
			continue;
		if( walked_cells[neigh->id] )
			continue;
		if( neigh->id == deadend_neighbor )
			continue;

		for( int jface=0; jface < mesh->faces_per_cell; jface++ ){

			// Neighbor of neighbor
			t_gmsh_neighbor *neigh_of_neigh = &( mesh->cell2neigh[neigh->id][jface] );

			if( neigh_of_neigh->type != CELL_REGULAR )
				continue;
			if( neigh_of_neigh->id == current_cell )
				continue;
			if( walked_cells[neigh_of_neigh->id] && neigh_of_neigh->id != starting_cell )
				continue;

			for( int kface=0; kface < mesh->faces_per_cell; kface++ ){

				// Neighbor of neighbor of neighbor
				t_gmsh_neighbor *neigh_of_neigh_of_neigh = &( mesh->cell2neigh[neigh_of_neigh->id][kface] );

				if( neigh_of_neigh_of_neigh->type != CELL_REGULAR )
					continue;
				if( neigh_of_neigh_of_neigh->id == current_cell )
					continue;
				if( walked_cells[neigh_of_neigh_of_neigh->id] && neigh_of_neigh_of_neigh->id != starting_cell )
					continue;

				if( mesh->cells[neigh_of_neigh_of_neigh->id].is_bcell )
					return neigh->id;
			}
		}
	}

	return -1;

	// If no cells was found, look for bnd cell's neighbors' neighbors
	for( int iface=0; iface < mesh->faces_per_cell; iface++ ){

		// Neighbor
		t_gmsh_neighbor *neigh = &( mesh->cell2neigh[current_cell][iface] ) ;
		if( neigh->type != CELL_REGULAR ){
			continue;
		}

		if( walked_cells[neigh->id] )
			continue;

		if( neigh->id == deadend_neighbor ){
			continue;
		}

		for( int jface=0; jface < mesh->faces_per_cell; jface ++ ){

			// Neighbor of neighbor
			t_gmsh_neighbor *neigh_of_neigh = &( mesh->cell2neigh[neigh->id][jface] );

			if( neigh_of_neigh->type != CELL_REGULAR )
				continue;

			if( walked_cells[neigh_of_neigh->id] )
				continue;

			if( bnd_neigh[neigh_of_neigh->id] != -1 )
				return neigh->id;
		}
	}

	// Loop of cell faces
	for( int iface=0; iface < mesh->faces_per_cell; iface++ ){

		// Neighbor ID adjacent to iface
		adjc_id = mesh->cell2neigh[current_cell][iface].id;

		// Skip boundary cell
		if( adjc_id == deadend_neighbor ){
			continue;
		}

		// Skip cell if already part of Hamiltonian path
		if( walked_cells[adjc_id] ){
			continue;
		}

		// Distance between neighbor center and grid center
		distance = sqrt( pow( mesh->cells[adjc_id].xc[0] - mesh->x[0], 2)
					   + pow( mesh->cells[adjc_id].xc[1] - mesh->x[1], 2) );

		// Get cell furthest from grid center
		if( max_distance < distance ){
			max_distance = distance;
			next_cell = adjc_id;
		}
	}

	return next_cell;
}




















































