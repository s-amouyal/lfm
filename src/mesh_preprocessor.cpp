#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <sstream>

#include <sys/stat.h>

#include "api/mesh_preprocessor.h"
#include "api/gmsh_reader.h"
#include "api/hpath.h"
#include "api/mesh_partitioning.h"
#include "api/mesh_writer.h"


using namespace std;
std::vector<int> bnd_hpath2original;

//***************************************************************************************************
//***************************************************************************************************
void preprocess( std::string filename, std::string file_base, int nSubparts ){

	MPI_env mpi_env;
	MeshWriter mesh_writer;

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Initialize MPI environment
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	mpi_env.initialize( false );

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Read and initialize global mesh
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	gmsh_mesh mesh_domain;
	read_gmsh_file( mpi_env, filename, &mesh_domain );
	mpi_env.barrier();

	if( mpi_env.is_master() ){

		struct stat info;

		if( stat( "./output/", &info ) != 0 )
			mkdir("output", 0777);

		if( stat( "./mesh_files/", &info ) != 0 )
			mkdir("mesh_files", 0777);

		tecplot_cell_ids( (char*) "output/mesh_original_new.plt", mesh_domain );
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Partition mesh into Np MPI ranks:
	// 		- gmsh_local: mesh local to the current MPI rank. Used to build boundary Hpath and
	// 						separate interior cells
	// 		- gmsh_int  : interior mesh only, further decomposed into submesh partitions
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	std::vector<int> mpi2global;
	std::vector<std::vector<std::vector<int>>> all_hpaths;

	// New approach
	gmsh_mesh mesh_local = partition_mesh_mpi_efficient( mpi_env, &mesh_domain );
	print_mesh_cells( file_base, &mesh_local );

	// Finc boundary Hpath
	Hpath hpath;
	hpath.get_hpath_boundary( &mesh_local );

	// Separate interior from boundary cells
	std::vector<int> old2new_int( mesh_local.cells.size(), -1 );
	std::vector<int> old2new_bnd( mesh_local.cells.size(), -1 );
	std::vector<int> tmp2old;
	std::vector<std::vector<int>> int2old;

	gmsh_mesh mesh_int = set_interior_mesh( &mesh_local, mpi_env, old2new_int, old2new_bnd, tmp2old );

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Apply second level of partitioning on interior cells
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	std::vector<gmsh_mesh> mesh_parts;
	mesh_parts = partition_local_mesh( &mesh_local, &mesh_int, nSubparts, old2new_int,
											  mpi_env, tmp2old, int2old );

	for( size_t ipart=0; ipart < mesh_parts.size(); ipart++ ){

		gmsh_mesh *mesh = &mesh_parts[ipart];

		// Direction of the face normals
		for( size_t icell=0; icell < mesh->cells.size(); icell++ ){
			for( int iface=0; iface < mesh->faces_per_cell; iface++ ){
				int gface = mesh->cell2face[ icell ][ iface ];

				// face nodes
				int fn0 = mesh->face2node[ gface ][ 0 ];
				int fn1 = mesh->face2node[ gface ][ 1 ];

				// cell nodes
				int cn0 = mesh->cell2node[ icell ][ iface ];
				int cn1 = mesh->cell2node[ icell ][ (iface+1) % mesh->faces_per_cell ];

				// Check wether the surface normal vector on each face is pointing outwards
				//      ---> true if the node ordering from the cell's and face's perspective is the same
				bool connectivity_right = (fn0 == cn0 && fn1 == cn1) || (fn0 == cn1 && fn1 == cn0) ;

				if( ! connectivity_right ){
					cout << "Issue with cells' faces ordering in preprocessor 1." << endl;
					cout << "icell =" << setw(10) << icell
								  	  << setw(10) << fn0
								  	  << setw(10) << fn1
								  	  << setw(10) << cn0
								  	  << setw(10) << cn1
								  	  << endl;
					MPI_Abort( MPI_COMM_WORLD, 986 );
				}
			}
		}
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Find Hamiltonian paths
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	for( int ipart=0; ipart < nSubparts+1; ipart++ ){

		// Skip boundary submesh
		if( ipart == INDEX_BND_SUBMESH ) {
			for( size_t icell=0; icell < mesh_parts[ipart].hpath.size(); icell++ ){
				mesh_parts[ipart].hpath[icell] = icell;
			}
			for( size_t icell=0; icell < mesh_parts[ipart].deadend_cells.size(); icell++ )
				mesh_parts[ipart].deadend_cells[icell] = icell + mesh_parts[ipart].hpath.size();

			continue;
		}

		// New mesh
		mesh_parts[ipart].id = ipart;
		hpath.get_hpath_interior( &mesh_parts[ipart] );
	}

	// Gather_interior_hpaths( gmsh_parts, mpi_env, int2old, mpi2global, all_hpaths );
	print_hpaths( file_base, mesh_parts );

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Update cells/faces/nodes ordering based on Hpath
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	// Index mapping: from local part ID to local Hpath
	std::vector<std::vector<int>> local2hpath;
	local2hpath = get_local2hpath_mapping( mesh_parts );

	// Nodes/faces
	for( size_t ipart=0; ipart < mesh_parts.size(); ipart++ ){
		update_node_ordering( mpi_env, ipart, &mesh_parts[ipart] );
	}

	// Cells
	for( size_t ipart=0; ipart < mesh_parts.size(); ipart++ ){
		update_mesh_cells( &mesh_parts[ipart], local2hpath );
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Print mesh files
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	output_mesh_files( file_base, mesh_parts );
}

//***************************************************************************************************
// Set submesh geometric parameters
//***************************************************************************************************
void update_node_ordering( MPI_env mpi_env, int ipart, gmsh_mesh *mesh ){

	int myrank;
	MPI_Comm_rank( MPI_COMM_WORLD, &myrank );

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// 		Reorder nodes:
	// 	Not all computations traverse the mesh along the cells. Solution variables are also computed
	// 	on all nodes. Therefore, for better performance, it is advantageous to reorder the nodes
	// 	according to their place in the boundary and interior paths, in that order.
	// 	Note that no Hamiltonian paths can be simultaneously obtained for the nodes and the cells.
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	std::vector<bool  > not_counted( mesh->nodes.size(), true );
	std::vector<size_t> hpath_nodes( mesh->nodes.size(), 0 );
	std::vector<size_t> nodes2hpath( mesh->nodes.size(), 0 );
	int gnode, gcell;
	int m_count = 0;

	// Interior cells
	for( size_t icell=0; icell < mesh->cells.size(); icell++ ){
		if( (int) icell < (int) mesh->hpath.size() ){
			if( ipart == INDEX_BND_SUBMESH ){
				gcell = icell;
			}else{
				gcell = mesh->hpath[icell];
			}
		}else{
			gcell = mesh->deadend_cells[icell-mesh->hpath.size()];
		}
		if( gcell == -1 ){
			cout << "Error in Hpath, check algorithm." << endl;
			MPI_Abort( MPI_COMM_WORLD, 912 );
		}

		for( int inode=0; inode < mesh->nodes_per_cell; inode++ ){

			gnode = mesh->cell2node[gcell][inode];

			if( not_counted[gnode] ){
				not_counted[gnode] = false;
				hpath_nodes[m_count++] = gnode;
				nodes2hpath[gnode] = m_count-1;
			}
		}
	}

	// Reorder nodes
	std::vector<gmsh_node> tmp_nodes( mesh->nodes.size() );

	for( size_t inode=0; inode < mesh->nodes.size(); inode++ ){
		tmp_nodes[inode].id_global = mesh->nodes[inode].id_global;

		tmp_nodes[inode].xn[0] = mesh->nodes[inode].xn[0];
		tmp_nodes[inode].xn[1] = mesh->nodes[inode].xn[1];
		tmp_nodes[inode].xn[2] = mesh->nodes[inode].xn[2];
	}


	for( size_t inode=0; inode < mesh->nodes.size(); inode++ ){
		// TODO: tmp fix
		gnode = hpath_nodes[inode];
		//gnode = inode;

		mesh->nodes[inode].id = gnode;
		mesh->nodes[inode].xn[0]  = tmp_nodes[gnode].xn[0];
		mesh->nodes[inode].xn[1]  = tmp_nodes[gnode].xn[1];
		mesh->nodes[inode].xn[2]  = tmp_nodes[gnode].xn[2];

		mesh->nodes[inode].id_global = tmp_nodes[gnode].id_global;
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Update connectivity
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	// cell2node
	std::vector<std::vector<int>> tmp_cell2node( mesh->cells.size(), std::vector<int> (mesh->nodes_per_cell, -1 ) );

	for( size_t icell=0; icell < mesh->cells.size(); icell++ ){
		for( int inode=0; inode < (int) mesh->nodes_per_cell; inode++ ){
			tmp_cell2node[icell][inode] = mesh->cell2node[icell][inode];
		}
	}

	for( size_t icell=0; icell < mesh->cells.size(); icell++ ){
		for( int inode=0; inode < (int) mesh->nodes_per_cell; inode++ ){
			gnode = tmp_cell2node[icell][inode];
			//gnode = inode;

			mesh->cell2node[icell][inode] = nodes2hpath[gnode];
		}
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Reorder faces
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	int gface;
	m_count = 0;

	not_counted.clear();
	not_counted.resize( mesh->faces.size(), true );

	std::vector<int> faces2hpath( mesh->faces.size(), 0 );
	std::vector<int> hpath_faces( mesh->faces.size(), 0 );

	// Traverse interior hpath
	for( size_t icell=0; icell < mesh->cells.size(); icell++ ){

		// The boundary hpath has already been set
		if( icell < mesh->hpath.size() ){
			if( ipart == INDEX_BND_SUBMESH ){
				gcell = icell;
			}else{
				gcell = mesh->hpath[icell];
			}
		}else{
			gcell = mesh->deadend_cells[icell-mesh->hpath.size()];
		}

		for( int iface=0; iface < (int) mesh->nodes_per_cell; iface++ ){
			gface = mesh->cell2face[gcell][iface];

			if( not_counted[gface] ){
				not_counted[gface] = false;
				faces2hpath[gface] = m_count;

				hpath_faces[m_count++] = gface;
			}
		}
	}

	// Reorder faces according to the hpaths
	for( size_t iface=0; iface < mesh->faces.size(); iface++ ){
		gface = hpath_faces[iface];
		//gface = iface;

		mesh->faces[iface].id = gface;
	}

	// face2node
	std::vector<std::vector<int>> tmp_face2node( mesh->faces.size(), std::vector<int> (mesh->nodes_per_face, -1) );

	for( size_t iface=0; iface < mesh->faces.size(); iface++ ){
		for( int inode=0; inode < mesh->nodes_per_face; inode++ ){
			tmp_face2node[iface][inode] = mesh->face2node[iface][inode];
		}
	}

	for( size_t iface=0; iface < mesh->faces.size(); iface++ ){
		gface = hpath_faces[iface];
		//gface = iface;

		for( int inode=0; inode < mesh->nodes_per_face; inode++ ){
			gnode = tmp_face2node[gface][inode];

			mesh->face2node[iface][inode] = nodes2hpath[gnode];
			//mesh->face2node[iface][inode] = gnode;
		}
	}

	// cell2face
	std::vector< std::vector<int> > tmp_cell2face;
	tmp_cell2face.resize( mesh->cells.size(), std::vector<int>( mesh->faces_per_cell, 0 ) );

	for( size_t icell=0; icell < mesh->cells.size(); icell++ ){
		for( int iface=0; iface < mesh->faces_per_cell; iface++ ){
			tmp_cell2face[icell][iface] = mesh->cell2face[icell][iface];
		}
	}

	for( size_t icell=0; icell < mesh->cells.size(); icell++ ){
		for( int iface=0; iface < mesh->faces_per_cell; iface++ ){
			gface = tmp_cell2face[icell][iface];

			mesh->cell2face[icell][iface] = faces2hpath[gface];
			//mesh->cell2face[icell][iface] = gface;
		}
	}

	// face2cell
	std::vector< std::vector<int> > tmp_face2cell;
	tmp_face2cell.resize( mesh->faces.size(), std::vector<int>( mesh->cells_per_face, -1 ) );

	for( size_t iface=0; iface < mesh->faces.size(); iface++ ){
		for( int icell=0; icell < mesh->cells_per_face; icell++ ){
			mesh->face2cell[iface][icell] = -1;
		}
	}

	for( size_t icell=0; icell < mesh->cells.size(); icell++ ){
		for( int iface=0; iface < mesh->faces_per_cell; iface++ ){
			gface = mesh->cell2face[icell][iface];

			if( mesh->face2cell[gface][0] == -1 ){
				mesh->face2cell[gface][0] = icell;
			}else if( mesh->face2cell[gface][1] == -1 ){
				mesh->face2cell[gface][1] = icell;
			}
		}
	}

	// Direction of the face normals
	for( size_t icell=0; icell < mesh->cells.size(); icell++ ){
		for( int iface=0; iface < mesh->faces_per_cell; iface++ ){
			int gface = mesh->cell2face[ icell ][ iface ];

			// face nodes
			int fn0 = mesh->face2node[ gface ][ 0 ];
			int fn1 = mesh->face2node[ gface ][ 1 ];

			// cell nodes
			int cn0 = mesh->cell2node[ icell ][ iface ];
			int cn1 = mesh->cell2node[ icell ][ (iface+1) % mesh->faces_per_cell ];

			// Check wether the surface normal vector on each face is pointing outwards
			//      ---> true if the node ordering from the cell's and face's perspective is the same
			bool connectivity_right = (fn0 == cn0 && fn1 == cn1) || (fn0 == cn1 && fn1 == cn0) ;

			if( ! connectivity_right ){
				cout << "Issue with cells' faces ordering in preprocessor " << myrank << " ." << endl;
				cout << "icell =" << setw(10) << icell
								  << setw(10) << fn0
								  << setw(10) << fn1
								  << setw(10) << cn0
								  << setw(10) << cn1
								  << endl;
				MPI_Abort( MPI_COMM_WORLD, 906 );
			}
		}
	}
}

//***************************************************************************************************
void update_mesh_cells( gmsh_mesh *mesh, std::vector<std::vector<int>> local2hpath ){

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Initialization
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	int myrank;
	MPI_Comm_rank( MPI_COMM_WORLD, &myrank );

	// Hpath cells
	std::vector<gmsh_cell> updated_cells( mesh->cells.size() );
	std::vector<std::vector<t_gmsh_neighbor>> updated_neighs( mesh->cells.size() );

	for( size_t icell=0; icell < mesh->cells.size(); icell++ ){
		if( icell < mesh->hpath.size() ){
			updated_cells[icell].id = mesh->hpath[icell];
		}else{
			updated_cells[icell].id = mesh->deadend_cells[icell-mesh->hpath.size()];
		}
		updated_neighs[icell].resize( mesh->faces_per_cell );
	}

	// old2new: mapping of old to new cell IDs
	std::vector<int> old2new( mesh->cells.size(), -1 );
	for( size_t icell=0; icell < mesh->cells.size(); icell++ ){
		if( icell < mesh->hpath.size() ){
			old2new[mesh->hpath[icell]] = icell;
		}else{
			old2new[mesh->deadend_cells[icell-mesh->hpath.size()]] = icell;
		}
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Update cells
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	// Default values
	for( size_t icell=0; icell < updated_cells.size(); icell++ ){
		for( int iface=0; iface < mesh->faces_per_cell; iface++ ){
			updated_neighs[icell][iface].type = CELL_NONE;
			updated_neighs[icell][iface].id    = 10000;
			updated_neighs[icell][iface].sm    = 10000;
			updated_neighs[icell][iface].proc  = 10000;
			updated_neighs[icell][iface].bound = 0;
		}
	}

	// Traverse Hpath and setup neighbor connectivity
	for( size_t icell=0; icell < mesh->cells.size(); icell++ ){

		int cell_id = updated_cells[icell].id;

		for( int iface=0; iface < mesh->faces_per_cell; iface++ ){

			t_gmsh_neighbor* neigh_old = &( mesh->cell2neigh[cell_id][iface] );

			// Neighbor is outside the boundary
			switch( neigh_old->type ){
				case CELL_REGULAR:
					updated_neighs[icell][iface].id   = old2new[neigh_old->id];
					updated_neighs[icell][iface].type = CELL_REGULAR;
					updated_neighs[icell][iface].sm   = mesh->id;
					break;
				case CELL_SUBMESH:
					updated_neighs[icell][iface].id   = local2hpath[neigh_old->sm][neigh_old->id];
					updated_neighs[icell][iface].type = CELL_SUBMESH;
					updated_neighs[icell][iface].sm   = neigh_old->sm;
					break;
				case CELL_MPI:
					updated_neighs[icell][iface].id   = neigh_old->id;
					updated_neighs[icell][iface].type = CELL_MPI;
					updated_neighs[icell][iface].sm   = neigh_old->proc;
					break;
				case CELL_PERIODIC:
					updated_neighs[icell][iface].id   = neigh_old->id;
					updated_neighs[icell][iface].type = CELL_PERIODIC;
					updated_neighs[icell][iface].sm   = neigh_old->proc;
					break;
				case CELL_PERIODIC_MPI:
					updated_neighs[icell][iface].id   = neigh_old->id;
					updated_neighs[icell][iface].type = CELL_PERIODIC_MPI;
					updated_neighs[icell][iface].sm   = neigh_old->proc;
					break;
				case CELL_NONE:
					updated_neighs[icell][iface].id    = 10000000;
					updated_neighs[icell][iface].type  = CELL_NONE;
					updated_neighs[icell][iface].sm    = 10000000;
					updated_neighs[icell][iface].bound = neigh_old->bound;
					break;
			}
		}
	}
	mesh->cells = updated_cells;
	for( size_t icell=0; icell < mesh->cells.size(); icell++ ){
		for( int iface=0; iface < mesh->faces_per_cell; iface++ ){
			mesh->cell2neigh[icell][iface].type  = updated_neighs[icell][iface].type  ;
			mesh->cell2neigh[icell][iface].id    = updated_neighs[icell][iface].id    ;
			mesh->cell2neigh[icell][iface].sm    = updated_neighs[icell][iface].sm    ;
			mesh->cell2neigh[icell][iface].proc  = updated_neighs[icell][iface].proc  ;
			mesh->cell2neigh[icell][iface].bound = updated_neighs[icell][iface].bound ;
		}
	}

	if( updated_cells.size() != mesh->cells.size() ){
		cout << "Not all boundary cells part of boundary Hpath =" << setw(10) << myrank
																  << setw(10) << mesh->id
																  << setw(10) << updated_cells.size()
																  << setw(10) << mesh->cells.size()
																  << endl;
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Connectivity
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	std::vector<std::vector<int>> tmp_cell2node;
	std::vector<std::vector<int>> tmp_cell2face;
	std::vector<std::vector<int>> tmp_face2cell;
	std::vector<std::vector<int>> tmp_face2node;

	tmp_cell2node.resize( mesh->cells.size(), std::vector<int>( mesh->nodes_per_cell, -1 ) );
	tmp_cell2face.resize( mesh->cells.size(), std::vector<int>( mesh->faces_per_cell, -1 ) );
	tmp_face2cell.resize( mesh->faces.size(), std::vector<int>( mesh->cells_per_face, -1 ) );
	tmp_face2node.resize( mesh->faces.size(), std::vector<int>( mesh->nodes_per_face, -1 ) );

	// Hpath cells
	for( size_t icell=0; icell < mesh->cells.size(); icell++ ){
		for( int inode=0; inode < mesh->nodes_per_cell; inode++ ){
			tmp_cell2node[icell][inode] = mesh->cell2node[mesh->cells[icell].id][inode];
		}
	}
	for( size_t icell = 0; icell < mesh->cells.size(); icell++ ){
		for( int iface = 0; iface < mesh->faces_per_cell; iface++ ){
			tmp_cell2face[icell][iface] = mesh->cell2face[mesh->cells[icell].id][iface];
		}
	}

	for( size_t iface = 0; iface < mesh->faces.size(); iface++ ){
		for( int inode = 0; inode < mesh->nodes_per_face; inode++ ){
			tmp_face2node[iface][inode] = mesh->face2node[mesh->faces[iface].id][inode];
		}
	}

	// Update all
	mesh->cell2node = tmp_cell2node;
	mesh->cell2face = tmp_cell2face;
	mesh->face2cell = tmp_face2cell;
	//mesh->face2node = tmp_face2node;

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Periodic faces
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	// Initialize all faces
	if( mesh->id != INDEX_BND_SUBMESH ){
		for( size_t iface=0; iface < mesh->faces.size(); iface++ ){
			mesh->faces[iface].twin_face.type  = PERIODIC_NEIGH_NONE;
			mesh->faces[iface].twin_face.index = -1;
			mesh->faces[iface].twin_face.rank  = -1;
		}
		return;
	}

	for( size_t icell=0; icell < mesh->cells.size(); icell++ ){
		for( int iface=0; iface < mesh->faces_per_cell; iface++ ){

			int gface = mesh->cell2face[icell][iface];

			if( mesh->cell2neigh[icell][iface].type != CELL_PERIODIC &&
				mesh->cell2neigh[icell][iface].type != CELL_PERIODIC_MPI ){

				mesh->faces[gface].twin_face.type  = PERIODIC_NEIGH_NONE;
				mesh->faces[gface].twin_face.index = -1;
				mesh->faces[gface].twin_face.rank  = -1;

				continue;
			}

			if( mesh->cell2neigh[icell][iface].type == CELL_PERIODIC ){
				mesh->faces[gface].twin_face.type = PERIODIC_NEIGH_REGULAR;
			}else{
				mesh->faces[gface].twin_face.type = PERIODIC_NEIGH_MPI;
			}
		}
	}
}

//***************************************************************************************************
// Index mapping from local part to local Hpath cell id
//***************************************************************************************************
std::vector<std::vector<int>> get_local2hpath_mapping( std::vector<gmsh_mesh> gmsh_parts ){

	int myrank;
	MPI_Comm_rank( MPI_COMM_WORLD, &myrank );

	// Initialize vector
	std::vector<std::vector<int>> local2hpath( gmsh_parts.size() );

	for( size_t ipart=0; ipart < gmsh_parts.size(); ipart++ ){
		local2hpath[ipart].resize( gmsh_parts[ipart].cells.size(), -1 );
	}

	// Map cells
	for( size_t ipart=0; ipart < gmsh_parts.size(); ipart++ ){
		// Hpath
		for( size_t icell=0; icell < gmsh_parts[ipart].hpath.size(); icell++ ){
			if( ipart == 0 ){
				local2hpath[ipart][icell] = icell;
			}else{
				local2hpath[ipart][gmsh_parts[ipart].hpath[icell]] = icell;
			}
		}
		// Deadends
		for( size_t icell=0; icell < gmsh_parts[ipart].deadend_cells.size(); icell++ ){
			if( ipart == 0 ){
				local2hpath[ipart][gmsh_parts[ipart].deadend_cells[icell]] = gmsh_parts[ipart].hpath.size() + icell;
			}else{
				local2hpath[ipart][gmsh_parts[ipart].deadend_cells[icell]] = gmsh_parts[ipart].hpath.size() + icell;
			}
		}
	}

	return local2hpath;
}

//***************************************************************************************************
void gather_boundary_hpaths( gmsh_mesh gmsh_local, MPI_env mpi_env, int nSubparts,
							 std::vector<int> mpi2global,
							 std::vector<std::vector<std::vector<int>>> all_hpaths ){

	// Local boundary hpaths
	std::vector<int> local_hpaths( gmsh_local.hpath.size() + gmsh_local.deadend_cells.size(), -1 );
	bnd_hpath2original.resize    ( gmsh_local.hpath.size() + gmsh_local.deadend_cells.size(), -1 );

	for( size_t icell=0; icell < local_hpaths.size(); icell++ ){
		if( icell < gmsh_local.hpath.size() )
			local_hpaths[icell] = mpi2global[gmsh_local.hpath[icell]];
		else
			local_hpaths[icell] = mpi2global[gmsh_local.deadend_cells[icell-gmsh_local.hpath.size()]];

		bnd_hpath2original[icell] = local_hpaths[icell];
	}

	// Exchange size of boundary Hpath
	int local_hpath_size = (int) local_hpaths.size();
	std::vector<int> all_bnd_sizes( mpi_env.size(), -1 );

	MPI_Allgather( &local_hpath_size, 1, MPI_INT,
				   &all_bnd_sizes[0], 1, MPI_INT, MPI_COMM_WORLD );

	// Allocate all_hpaths vector
	if( mpi_env.is_master() ){
		all_hpaths.resize( mpi_env.size() );
		for( int irank=0; irank < mpi_env.size(); irank++ ){
			all_hpaths[irank].resize( nSubparts+1 );
			all_hpaths[irank][0].resize( all_bnd_sizes[irank], -1 );
		}
	}

	int total_bnd_sizes = 0;
	for( int irank=0; irank < mpi_env.size(); irank++ ){
		total_bnd_sizes += all_bnd_sizes[irank];
	}

	std::vector<int> tmp_hpaths( total_bnd_sizes, -1 );

	std::vector<int> disps( mpi_env.size(), -1 );
	for( int irank=0; irank < mpi_env.size(); irank++ ){
		if( irank == 0 )
			disps[irank] = 0;
		else
			disps[irank] = all_bnd_sizes[irank-1] + disps[irank-1];
	}

	MPI_Gatherv( &local_hpaths[0], local_hpaths.size(), MPI_INT,
				 &tmp_hpaths[0], &all_bnd_sizes[0], &disps[0], MPI_INT, 0, MPI_COMM_WORLD );

	if( mpi_env.is_master() ){
		int local_count = 0;
		int m_rank = 0;

		for( int icell=0; icell < total_bnd_sizes; icell++ ){
			if( local_count == all_bnd_sizes[m_rank] ){
				m_rank++;
				local_count = 0;
			}

			all_hpaths[m_rank][0][local_count] = tmp_hpaths[icell];

			local_count++;
		}
	}
}

//***************************************************************************************************
void gather_interior_hpaths( std::vector<gmsh_mesh> gmsh_parts, MPI_env mpi_env,
							 std::vector<std::vector<int>> int2old, std::vector<int> mpi2global,
							 std::vector<std::vector<std::vector<int>>> all_hpaths ){

	// Initialize single vector for all parts
	int total_local_cells = 0;
	for( size_t ipart=0; ipart < gmsh_parts.size(); ipart++ ){
		total_local_cells += gmsh_parts[ipart].cells.size();
	}

	std::vector<int> local_hpaths( total_local_cells, -1 );

	//
	int m_count = 0;
	for( size_t ipart=0; ipart < gmsh_parts.size(); ipart++ ){

		// Skip boundary submesh
		if( ipart == INDEX_BND_SUBMESH ){
			for( size_t icell=0; icell < gmsh_parts[ipart].cells.size(); icell++ ){
				local_hpaths[m_count++] = bnd_hpath2original[icell];
			}
			continue;
		}

		// Record Hpath/deadend cells
		for( size_t icell=0; icell < gmsh_parts[ipart].hpath.size(); icell++ ){
			local_hpaths[m_count++] = mpi2global[int2old[ipart][gmsh_parts[ipart].hpath[icell]]];
		}
		for( size_t icell=0; icell < gmsh_parts[ipart].deadend_cells.size(); icell++ ){
			local_hpaths[m_count++] = mpi2global[int2old[ipart][gmsh_parts[ipart].deadend_cells[icell]]];
		}
	}

	// Exchange size of Hpaths
	int local_hpaths_size = local_hpaths.size();
	std::vector<int> all_hpaths_sizes( mpi_env.size() , -1 );

	MPI_Gather( &local_hpaths_size, 1, MPI_INT,
				&all_hpaths_sizes[0], 1, MPI_INT, 0, MPI_COMM_WORLD );

	// Exchange Hpaths
	std::vector<int> recvcounts( mpi_env.size(), -1 );
	std::vector<int> disps     ( mpi_env.size(), -1 );

	if( mpi_env.is_master() ){
		for( int irank=0; irank < mpi_env.size(); irank++ ){
			recvcounts[irank] = all_hpaths_sizes[irank];

			if( irank == 0 )
				disps[irank] = 0;
			else
				disps[irank] = disps[irank-1] + all_hpaths_sizes[irank-1];
		}
	}

	int total_global_cells = -1;
	MPI_Allreduce( &total_local_cells, &total_global_cells, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
	std::vector<int> tmp_hpaths( total_global_cells, -1 );


	MPI_Gatherv( &local_hpaths[0], local_hpaths.size(), MPI_INT,
				 &tmp_hpaths[0], &recvcounts[0], &disps[0], MPI_INT, 0, MPI_COMM_WORLD );


	std::vector<int> cell2rank;

	if( mpi_env.is_master() ){

		cell2rank.resize( tmp_hpaths.size(), -1 );

		int m_rank = 0;
		int threshold = all_hpaths_sizes[0];
		for( size_t icell=0; icell < tmp_hpaths.size(); icell++ ){
			if( (int) icell >= threshold ){
				threshold += all_hpaths_sizes[++m_rank];
			}

			cell2rank[icell] = m_rank;
		}
	}

	if( mpi_env.is_master() ){
		ofstream fout;
		char* filename = (char*) "LidDrivenCavity_Hpath_0050K_mpi01_sm01_lines.txt";
		fout.open( filename );

		for( size_t icell=0; icell < tmp_hpaths.size(); icell++ ){
			fout << setw(10) << cell2rank[icell]
				 << setw(10) << icell
				 << setw(10) << tmp_hpaths[icell]
				 << endl;
		}
	}
}


//***************************************************************************************************
// Output mesh files for each processor
//***************************************************************************************************

//***************************************************************************************************
void output_mesh_files( std::string file_base, std::vector<gmsh_mesh> gmsh_parts ){

	int myrank;
	MPI_Comm_rank( MPI_COMM_WORLD, &myrank );

	std::stringstream tmp_filename;

	if( myrank < 10 ){
		tmp_filename << "mesh_files/" << file_base << "_p00" << myrank << ".lfm" ;
	}else if( myrank < 100 ){
		tmp_filename << "mesh_files/" << file_base << "_p0"  << myrank << ".lfm" ;
	}else{
		tmp_filename << "mesh_files/" << file_base << "_p"   << myrank << ".lfm" ;
	}

	std::string filename = tmp_filename.str();

	MeshWriter mesh_writer;

	mesh_writer.output_mesh_files( filename, gmsh_parts );
}

//***************************************************************************************************
void print_mesh_cells( std::string file_base, gmsh_mesh *gmsh ){

	int myrank;
	MPI_Comm_rank( MPI_COMM_WORLD, &myrank );

	std::string file_string;

	if( myrank < 10 ){
		file_string = "output/" + file_base + "_original_p00" + std::to_string(myrank) + ".plt";
	}else if( myrank < 100 ){
		file_string = "output/" + file_base + "_original_p0"  + std::to_string(myrank) + ".plt";
	}else if( myrank < 1000){
		file_string = "output/" + file_base + "_original_p"   + std::to_string(myrank) + ".plt";
	}
	char *filename = (char*) file_string.c_str();

	MeshWriter mesh_writer;
	mesh_writer.tecplot_cell_ids( filename, gmsh );
}

//***************************************************************************************************
void tecplot_cell_ids( char* filename, gmsh_mesh gmsh ){

	// Open file
	ofstream fout;
	fout.open( filename );

	// Tecplot header
	fout << "TITLE = \"Example: 2D Finite Element Data\"" << endl;
	fout << "VARIABLES = \"X\", \"Y\", \"cell_id\"" << endl;

	// Assume triangular elements for now
	if( gmsh.nodes_per_cell == 3 )
	fout << "ZONE NODES = " << setw(5) << gmsh.nodes.size()
		 << ", ELEMENTS = " << setw(5) << gmsh.cells.size()
		 << ", DATAPACKING=BLOCK"
		 << ", ZONETYPE=FETRIANGLE"
		 << ", VarLocation=([3]=CellCentered)"
		 << endl;
	if( gmsh.nodes_per_cell == 4 )
	fout << "ZONE NODES = " << setw(5) << gmsh.nodes.size()
		 << ", ELEMENTS = " << setw(5) << gmsh.cells.size()
		 << ", DATAPACKING=BLOCK"
		 << ", ZONETYPE=FEQUADRILATERAL"
		 << ", VarLocation=([3]=CellCentered)"
		 << endl;

	// Node coordinates
	for( auto inode=gmsh.nodes.begin(); inode < gmsh.nodes.end(); inode++ ){
		fout << setw(15) << inode->xn[0];
	}
	fout << endl;

	for( auto inode=gmsh.nodes.begin(); inode < gmsh.nodes.end(); inode++ ){
		fout << setw(15) << inode->xn[1];
	}
	fout << endl;

	// Cell IDs
	for( size_t icell = 0; icell < gmsh.cells.size(); icell++ ){
		fout << setw(15) << gmsh.cells[icell].id;
	}
	fout << endl;

	// Mesh connectivity
	int gcell;
	for( size_t icell = 0; icell < gmsh.cells.size(); icell++ ){
		gcell = gmsh.cells[icell].id;

		for( int inode = 0; inode < gmsh.nodes_per_cell; inode++ ){
			fout << setw(15) << gmsh.cell2node[gcell][inode] + 1;
		}
		fout << endl;
	}
	fout.close();
}

//***************************************************************************************************
void print_hpaths( std::string file_base, std::vector<gmsh_mesh> gmsh_parts ){

	std::string tmp_filename;

	int myrank;
	MPI_Comm_rank( MPI_COMM_WORLD, &myrank );

	if( myrank < 10 ){
		tmp_filename = "output/" + file_base + "_hpath_p00" + std::to_string(myrank);
	}else if( myrank < 100 ){
		tmp_filename = "output/" + file_base + "_hpath_p0"  + std::to_string(myrank);
	}else{
		tmp_filename = "output/" + file_base + "_hpath_p"   + std::to_string(myrank);
	}

	MeshWriter mesh_writer;
	for( size_t ipart=0; ipart < gmsh_parts.size(); ipart++ ){

		std::string filename;

		if( ipart < 10 ){
			filename = tmp_filename + "_00" + std::to_string(ipart) + ".plt";
		}else if( ipart < 100 ){
			filename = tmp_filename + "_0"  + std::to_string(ipart) + ".plt";
		}else{
			filename = tmp_filename + "_"   + std::to_string(ipart) + ".plt";
		}

		mesh_writer.tecplot_hpath_nosplit( filename, &gmsh_parts[ipart], gmsh_parts[ipart].hpath );
	}
}




















































