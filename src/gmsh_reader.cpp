#include <iostream>
#include <iomanip>
#include <fstream>
#include "math.h"

#include "api/gmsh_reader.h"

using namespace std;

const double PI = 4*atan(1.);

//***************************************************************************************************
void read_gmsh_file( MPI_env mpi_env, std::string filename, gmsh_mesh *mesh ){

	// To minimize memory consumption, have only 1 core per node read the mesh
	if( !mpi_env.is_node_master() )
		return;

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Open file
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	ifstream fin;

	fin.open( filename );

	if( fin.fail() ){
		cout << "Failed to open file = " << filename << endl;
		exit(0);
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Physical Names
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	unsigned int num_boundaries;
	read_boundaries( fin, num_boundaries );

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Entities
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	vector<vector<vector<int>>> entity2boundary;
	entity2boundary.resize(4);
	read_entities( fin, num_boundaries, entity2boundary );

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Nodes
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	std::vector<std::vector<int>> nodes_bound;
	read_nodes( fin, entity2boundary, nodes_bound, mesh );

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Cells
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	int n_neigh = -1;
	std::vector<std::vector<gmsh_neighbor>> gmsh_neighs;
	read_cells( fin, n_neigh, mesh, gmsh_neighs );

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Get periodic boundary conditions (if exist)
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	vector<vector<int>> periodic_nodes( mesh->nodes.size(), vector<int>() );
	get_periodic_boundary_conditions( mpi_env, fin, entity2boundary, periodic_nodes );

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Get all faces
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	std::vector<int> periodic_faces;

	init_faces              ( n_neigh, nodes_bound, periodic_nodes, periodic_faces, mesh );
	reorder_cell_local_nodes( n_neigh, mesh );
	set_cells_faces         ( n_neigh, mesh );
	reorder_cell_local_faces( n_neigh, mesh );

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Set boundary conditions
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	set_boundary_conditions( n_neigh, mesh, nodes_bound, periodic_faces, gmsh_neighs );

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Element type: for now, only uniform-element mesh is supported
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// TODO: wrong statement, use enums
	mesh->elmnt_type = n_neigh - 2;

	if( mesh->elmnt_type == 1 ){
		mesh->nodes_per_cell = 3;
		mesh->faces_per_cell = 3;
		mesh->cells_per_face = 2;
		mesh->nodes_per_face = 2;
	}else if( mesh->elmnt_type == 2 ){
		mesh->nodes_per_cell = 4;
		mesh->faces_per_cell = 4;
		mesh->cells_per_face = 2;
		mesh->nodes_per_face = 2;
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Compute cell neighbors and mesh paramters
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	set_cell_neighbors( mesh, gmsh_neighs );
	compute_mesh_parameters( mpi_env, mesh );

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Print mesh info
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	if( mpi_env.is_master() ){
		cout << endl;
		cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
		cout << "n_nodes = " << mesh->nodes.size() << endl;
		cout << "n_faces = " << mesh->faces.size() << endl;
		cout << "n_cells = " << mesh->cells.size() << endl << endl;

		if( mesh->elmnt_type == 1 ){
			cout << "elements = 2D triangles" << endl;
		}else if( mesh->elmnt_type == 2 ){
			cout << "elements = 2D quadrilateral" << endl;
		}
		cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl << endl;
		//cout << "xmin/xmax   = " << setw(15) << grid_xmin[0] << setw(15) << grid_xmax[0] << endl;
		//cout << "ymin/ymax   = " << setw(15) << grid_xmin[1] << setw(15) << grid_xmax[1] << endl << endl;

		cout << "grid center = " << setw(15) << mesh->x[0] << setw(15) << mesh->x[1] << endl;
		cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl << endl;
	}
}

//***************************************************************************************************
void skip_lines( ifstream *fin, int num_skip_lines ){

	for (int ii = 0; ii < num_skip_lines; ii++){
		fin->ignore( 1000, '\n' );
	}
}

//***************************************************************************************************
void read_boundaries( ifstream &fin, unsigned int &num_boundaries ){

	std::string tp;

	while(tp.compare("$PhysicalNames") != 0) {
		fin >> tp;
		skip_lines( &fin, 1 );
	}

	int n_boundaries;
	fin >> n_boundaries;
	num_boundaries = n_boundaries;

	skip_lines( &fin, n_boundaries+1 );
}

//***************************************************************************************************
void read_entities( ifstream &fin, unsigned int num_boundaries, vector<vector<vector<int>>> &entity2boundary ){

	std::string tp;
	double dummy_double;
	int dummy_int;

	while(tp.compare("$Entities") != 0) {
		fin >> tp;
		skip_lines( &fin, 1 );
	}

	int n_nodes, n_faces, n_cells, n_volumes;
	fin >> n_nodes >> n_faces >> n_cells >> n_volumes;

	entity2boundary[0].resize(n_nodes);
	entity2boundary[1].resize(n_faces);
	entity2boundary[2].resize(n_cells);
	entity2boundary[3].resize(n_volumes);

	for (int n=0; n<4; n++) {
		for (size_t i=0; i<entity2boundary[n].size(); i++) {
			for (int j=0; j<(3 + 3 * (n ? 1 : 0) + 1); j++) {
				fin >> dummy_double;
			}

			int num_entities, entity;
			fin >> num_entities;
			if (num_entities>0) {
				for (int j=0; j<num_entities; j++) {
					fin >> entity;
					if( entity == (int) num_boundaries )
						entity2boundary[n][i].push_back(0);
					else
						entity2boundary[n][i].push_back(entity);
				}

				if ( n > 0 ) {
					int j_final;
					fin >> j_final;
					for (int j=0; j<j_final; j++) {
						fin >> dummy_int;
						for (size_t j=0; j<entity2boundary[n][i].size(); j++) {
							if (entity2boundary[n][i][j]>0)
								entity2boundary[n-1][abs(dummy_int)-1].push_back(entity2boundary[n][i][j]);
						}
					}
				}
			} else if ( n > 0 ){
				entity2boundary[n][i].push_back(0);
				int j_final;
				fin >> j_final;
				for (int j=0; j<j_final; j++) {
					fin >> dummy_int;
				}
			}
		}
	}
}

//***************************************************************************************************
void read_nodes( ifstream &fin, vector<vector<vector<int>>> &entity2boundary,
				 std::vector<std::vector<int>> &nodes_bound, gmsh_mesh *mesh ){

	std::string tp;
	int dummy_int;
	int n_blocks;

	while(tp.compare("$Nodes") != 0) {
		fin >> tp;
		skip_lines( &fin, 1 );
	}

	int n_nodes, min_node, max_node;

	fin >> n_blocks >> n_nodes >> min_node >> max_node ;

	// Allocate nodes objects
	mesh->nodes.resize( n_nodes );
	nodes_bound.resize( n_nodes );

	int entity_dim, entity, local_block_nodes;
	int m_count = 0;

	// Loop over all the node blocks
	for( int iblocks=0; iblocks < n_blocks; iblocks++ ){
		fin >> entity_dim >> entity >> dummy_int >> local_block_nodes;

		// Get the global id of nodes of the current block
		for( int ilocal_node = 0; ilocal_node < local_block_nodes; ilocal_node++ ){
			fin >> mesh->nodes[m_count].id;

			for (size_t ibound=0; ibound<entity2boundary[entity_dim][entity-1].size(); ibound++) {
				nodes_bound[m_count].push_back( entity2boundary[entity_dim][entity-1][ibound] );
			}

			m_count ++;
		}
		m_count -= local_block_nodes;
		// Get x,y,z coordinates of nodes of the current block
		for( int ilocal_node = 0; ilocal_node < local_block_nodes; ilocal_node++ ){
			fin >> mesh->nodes[m_count].xn[0]
				>> mesh->nodes[m_count].xn[1]
				>> mesh->nodes[m_count].xn[2] ;

			m_count ++;
		}
	}
	skip_lines( &fin, 1 );

	// Rescale node IDs so they start at 0
	for( size_t inode=0; inode < mesh->nodes.size(); inode++ ){
		mesh->nodes[inode].id = mesh->nodes[inode].id - 1;
	}
}

//***************************************************************************************************
void read_cells( ifstream &fin, int &n_neigh, gmsh_mesh *mesh,
				 std::vector<std::vector<gmsh_neighbor>> &gmsh_neighs ){

	skip_lines( &fin, 2 );

	int n_cells, n_blocks, dummy_int;

	fin >> n_blocks >> n_cells >> dummy_int >> dummy_int;

	int block_type, local_block_cells;

	// Skip elements which are not cells
	int iblock;
	for( iblock = 0; iblock < n_blocks; iblock++ ){

		fin >> block_type >> dummy_int >> n_neigh >> local_block_cells;

		if( block_type == 0 ){
			fin >> dummy_int >> dummy_int;
		} else if( block_type == 1 ){
			for( int iline=0; iline < local_block_cells; iline++ ){
				fin >> dummy_int >> dummy_int >> dummy_int;
			}
		} else if( block_type > 1 ){
			// 2D Element, exit loop
			break;
		}

		n_cells -= local_block_cells;
	}

	// Read in cells
	n_neigh  += 1;
	n_blocks -= iblock;

	// Allocate
	mesh->cells    .resize( n_cells );
	gmsh_neighs    .resize( n_cells, std::vector<gmsh_neighbor>( n_neigh ) );
	mesh->cell2node.resize( n_cells, std::vector<int>( n_neigh, -1 ) );

	int cell_num = 0;

	while ( n_blocks > 0 ) {
		n_blocks -= 1;

		for( int icell=0; icell < local_block_cells; icell++ ){
			fin >> mesh->cells[cell_num].id;
			for( int inode=0; inode < n_neigh; inode++ ){
				fin >> mesh->cell2node[cell_num][inode];
			}

			cell_num += 1;
		}
		if (n_blocks > 0)
			fin >> block_type >> dummy_int >> dummy_int >> local_block_cells;
	}

	// Reset cells/nodes IDs so they start from 0
	int zeroth_id = mesh->cells[0].id;
	for( size_t icell=0; icell < mesh->cells.size(); icell++ ){
		mesh->cells[icell].id -= zeroth_id;
		for( int inode=0; inode < n_neigh; inode++ ){
			mesh->cell2node[icell][inode] -= 1;
		}
	}
}

//***************************************************************************************************
void get_periodic_boundary_conditions( MPI_env mpi_env, ifstream &fin,
									   vector<vector<vector<int>>>  &entity2boundary,
									   std::vector<std::vector<int>> &periodic_nodes){

	std::string tp;
	int dummy_int;

	skip_lines( &fin, 2 );

	fin >> tp;

	if (tp.compare("$Periodic") == 0) {
		if( mpi_env.is_master() ){
			cout << "Periodic boundary conditions discovered..." << endl;
		}

		int n_periodics;

		fin >> n_periodics;

		int entity_dim, entity1, entity2, local_block_cells, periodic_type;
		int periodic_blocks = 0;

		// Skip elements which are not lines
		for( int iblock = 0; iblock < n_periodics; iblock++ ){

			fin >> entity_dim >> entity1 >> entity2;

			skip_lines( &fin, 2 );

			if( entity_dim == 0 ){
				fin >> local_block_cells;
				for( int iline=0; iline < local_block_cells; iline++ ){
					fin >> dummy_int >> dummy_int;
				}
			} else if( entity_dim == 1 ){
				// Line, exit loop
				periodic_blocks = n_periodics - iblock;
				break;
			}
		}

		int temp1_int, temp2_int;

		for( int iblock = 0; iblock < periodic_blocks; iblock++ ){
			fin >> local_block_cells;

			periodic_type = entity2boundary[entity_dim][entity1-1][0];

			for( int iline=0; iline < local_block_cells; iline++ ){
				fin >> temp1_int >> temp2_int;

				periodic_nodes[temp1_int - 1].push_back(periodic_type);
				periodic_nodes[temp2_int - 1].push_back(periodic_type);

				periodic_nodes[temp1_int - 1].push_back(temp2_int - 1);
				periodic_nodes[temp2_int - 1].push_back(temp1_int - 1);
			}

			fin >> entity_dim >> entity1 >> entity2;

			skip_lines( &fin, 2 );
		}
	}

	fin.close();
}

//***************************************************************************************************
void init_faces( int n_neigh,
				 std::vector<std::vector<int>> nodes_bound,
				 std::vector<vector<int>>     &periodic_nodes,
				 std::vector<int>             &periodic_faces,
				 gmsh_mesh                    *mesh ){

	mesh->faces.resize( n_neigh*mesh->nodes.size() );

	// Connectivity vectors
	mesh->cell2face.resize( mesh->cells.size(), std::vector<int>( n_neigh, -1 ) );
	mesh->face2node.resize( mesh->faces.size(), std::vector<int>(       2, -1 ) );
	mesh->face2cell.resize( mesh->faces.size(), std::vector<int>(       2, -1 ) );

	//
	vector<vector<int>> nodes_neighbors( mesh->nodes.size(), vector<int>( 12, -1 ));
	vector<vector<int>> face_ids       ( mesh->nodes.size(), vector<int>( 12, -1 ));
	vector<int> periodic_faces_temp[6];

	int face_count = 0;
	int gnode0, gnode1, n0, n1;
	int gface;

	bool face_exists;

	for( size_t icell=0; icell < mesh->cells.size(); icell++ ){
		for( int inode=0; inode < n_neigh; inode++ ){
			n0 = inode;
			n1 = (inode+1) % n_neigh;

			gnode0 = mesh->cell2node[icell][n0];
			gnode1 = mesh->cell2node[icell][n1];

			if( mesh->cell2node[icell][n0] < mesh->cell2node[icell][n1] ){
				gnode0 = mesh->cell2node[icell][n0];
				gnode1 = mesh->cell2node[icell][n1];
			}else{
				gnode0 = mesh->cell2node[icell][n1];
				gnode1 = mesh->cell2node[icell][n0];
			}

			// Check if face already exists
			int local_node = -1;
			face_exists = false;
			for( size_t jnode=0; jnode < 12; jnode++ ){
				if( nodes_neighbors[gnode0][jnode] == -1 ){
					local_node = jnode;
					break;
				}
				if( nodes_neighbors[gnode0][jnode] == gnode1 ){
					face_exists = true;
					local_node  = jnode;
					break;
				}
			}

			if( !face_exists ){
				mesh->faces[face_count].id  = face_count;
				face_ids[gnode0][local_node] = face_count;
				nodes_neighbors[gnode0][local_node] = gnode1;

				// Save face
				if( mesh->cell2node[icell][n0] < mesh->cell2node[icell][n1] ){
					mesh->face2node[face_count][0] = mesh->cell2node[icell][n0];
					mesh->face2node[face_count][1] = mesh->cell2node[icell][n1];
				}else{
					mesh->face2node[face_count][0] = mesh->cell2node[icell][n1];
					mesh->face2node[face_count][1] = mesh->cell2node[icell][n0];
				}


				// Boundary
				for( size_t ii=0; ii < nodes_bound[gnode0].size(); ii++ ){
					for( size_t jj=0; jj < nodes_bound[gnode1].size(); jj++ ){
						if (nodes_bound[gnode0][ii] == nodes_bound[gnode1][jj] && nodes_bound[gnode0][ii]!=0) {
							mesh->faces[face_count].bound = nodes_bound[gnode0][ii];
							ii = nodes_bound[gnode0].size();
							jj = nodes_bound[gnode1].size();
						}
					}
				}

				// Periodic Boundary Conditions (if exist)
				int fnode0 = mesh->face2node[face_count][0];
				int fnode1 = mesh->face2node[face_count][1];
				if( !periodic_nodes[fnode0].empty() && !periodic_nodes[fnode1].empty()) {
					for( size_t ii=0; ii < periodic_nodes[fnode0].size()/2; ii++) {
						for( size_t jj=0; jj < periodic_nodes[fnode1].size()/2; jj++) {
							if( periodic_nodes[fnode0][ii*2] != periodic_nodes[fnode1][jj*2])
								continue;

							periodic_faces_temp[0].push_back(face_count);
							periodic_faces_temp[1].push_back(fnode0);
							periodic_faces_temp[2].push_back(fnode1);
							periodic_faces_temp[3].push_back(-1);
							periodic_faces_temp[4].push_back(periodic_nodes[fnode0][ii*2+1]);
							periodic_faces_temp[5].push_back(periodic_nodes[fnode1][jj*2+1]);

						}
					}
				}
			}

			if( !face_exists ){
				gface = face_count;
			}else{
				gface = face_ids[gnode0][local_node];
			}

			if( mesh->face2cell[gface][0] == -1 ){
				mesh->face2cell[gface][0] = mesh->cells[icell].id;
			}else{
				mesh->face2cell[gface][1] = mesh->cells[icell].id;
			}

			if( !face_exists ){
				face_count++;
			}
		}
	}
	mesh->faces.resize( face_count );

	for( int iface=0; iface < int(periodic_faces_temp[0].size()-1); iface++) {
		for( int jface=iface+1; jface<int(periodic_faces_temp[0].size()); jface++){

			if( periodic_faces_temp[3][jface] == -1) {
				int gface  = mesh->faces[periodic_faces_temp[0][jface]].id;
				int fnode0 = mesh->face2node[gface][0];
				int fnode1 = mesh->face2node[gface][1];

				if( (fnode0 == periodic_faces_temp[4][iface] && fnode1 == periodic_faces_temp[5][iface]) ||
					(fnode1 == periodic_faces_temp[4][iface] && fnode0 == periodic_faces_temp[5][iface])) {
						periodic_faces_temp[3][iface] = periodic_faces_temp[0][jface];
						periodic_faces_temp[3][jface] = periodic_faces_temp[0][iface];
				}
			}
		}
	}

	periodic_faces.resize( mesh->faces.size() );
	for ( size_t iface=0; iface < mesh->faces.size(); iface++ ) {
		periodic_faces[iface] = -1;
	}

	for (size_t iface=0; iface<periodic_faces_temp[0].size(); iface++) {
		periodic_faces[periodic_faces_temp[0][iface]] = periodic_faces_temp[3][iface];
		periodic_faces[periodic_faces_temp[3][iface]] = periodic_faces_temp[0][iface];
	}

	for( size_t iface=0; iface < mesh->faces.size(); iface++ ){
		mesh->faces[iface].twin_face.type    = PERIODIC_NEIGH_NONE;
		mesh->faces[iface].twin_face.rank    = -1;
		mesh->faces[iface].twin_face.index   = -1;
		mesh->faces[iface].twin_face.orig_id = -1;
	}

	int ind_face = 0;
	for( size_t iface=0; iface < periodic_faces.size(); iface++ ){
		if( periodic_faces[iface] == -1 )
			continue;

		mesh->faces[iface].twin_face.index   = periodic_faces[iface];
		mesh->faces[iface].twin_face.orig_id = periodic_faces[iface];

		ind_face++;
	}

	// Check periodic faces
	for( size_t iface=0; iface < mesh->faces.size(); iface++ ){
		if( mesh->faces[iface].twin_face.index == -1 )
			continue;

		int twin_face = mesh->faces[iface].twin_face.index;

		if( mesh->faces[twin_face].twin_face.index != (int) iface ){
			cout << "Issue with periodic faces, terminating simulation prematurely." << endl;
			MPI_Abort( MPI_COMM_WORLD, 943 );
		}
	}
}

//***************************************************************************************************
void set_cells_faces( int n_neigh, gmsh_mesh *mesh ){

	// Store faces ID in cells objects
	int jface;

	for( size_t iface=0; iface < mesh->faces.size(); iface++ ){

		// cell0 id
		int fcell0 = mesh->face2cell[iface][0];

		// skip if boundary cell
		if( fcell0 > -1 ) {
			jface = 0;
			while( jface < n_neigh && mesh->cell2face[fcell0][jface] != -1 ){
				jface++;
			}

			if( jface < n_neigh ){
				// store face id
				mesh->cell2face[fcell0][jface] = mesh->faces[iface].id;
			}
		}

		// cell1 id
		int fcell1 = mesh->face2cell[iface][1];

		// skip if boundary cell
		if( fcell1 > -1 ){
			jface = 0;
			while( jface < n_neigh && mesh->cell2face[fcell1][jface] != -1 ){
				jface++;
			}

			if( jface < n_neigh ){
				// store face id
				mesh->cell2face[fcell1][jface] = mesh->faces[iface].id;
			}
		}
	}
}

//***************************************************************************************************
// For each cell, reorder the local nodes numbering counter-clockwise
//***************************************************************************************************
void reorder_cell_local_nodes( int n_neigh, gmsh_mesh *mesh ){

	const int DIM_CNT = 2;

	std::vector<std::vector<float>> angles( mesh->cells.size(), std::vector<float>( n_neigh, 0.0 ) );

	for( size_t icell=0; icell < mesh->cells.size(); icell++ ){

		// Cell centroid
		double xc[DIM_CNT] = { 0.0 };
		for( int inode=0; inode < n_neigh; inode++ ){
			int gnode = mesh->cell2node[icell][inode];
			xc[0] += mesh->nodes[gnode].xn[0] ;
			xc[1] += mesh->nodes[gnode].xn[1] ;
		}
		xc[0] /= (double) n_neigh;
		xc[1] /= (double) n_neigh;

		// Angles for each node
		for( int inode=0; inode < n_neigh; inode++ ){
			int gnode = mesh->cell2node[icell][inode];
			angles[icell][inode] = atan2( ( mesh->nodes[gnode].xn[1] - xc[1] ) , ( mesh->nodes[gnode].xn[0] - xc[0] ));
		}

		// Convert to degrees
		for( int inode=0; inode < n_neigh; inode++ ){
			angles[icell][inode] *= 360.0 / (2.0*PI) ;
		}

		// Keep angles between 0 and 360 degrees
		for( int inode=0; inode < n_neigh; inode++ ){
			angles[icell][inode] = angles[icell][inode] < 0 ? angles[icell][inode]+360.0
															: angles[icell][inode] ;
		}

		// Reorder nodes counter-clockwise based on angles
		for( int inode=0; inode < n_neigh-1; inode++ ){
			for( int jnode=inode+1; jnode < n_neigh; jnode++ ){
				if( angles[icell][inode] > angles[icell][jnode] ){
					double tmp_angle = angles[icell][inode];

					angles[icell][inode] = angles[icell][jnode];
					angles[icell][jnode] = tmp_angle;

					int tmp_node = mesh->cell2node[icell][inode];

					mesh->cell2node[icell][inode] = mesh->cell2node[icell][jnode];
					mesh->cell2node[icell][jnode] = tmp_node;
				}
			}
		}
	}
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// For each cell, reorder the local nodes numbering counter-clockwise
void reorder_cell_local_faces( int n_neigh, gmsh_mesh *mesh ){

	int fnode0, fnode1;
	int cnode0, cnode1;

	int inode, gface;

	for( size_t icell=0; icell < mesh->cells.size(); icell++ ){
		for( int iface=0; iface < n_neigh; iface++ ){
			inode = iface;

			// Face nodes obtained from cell
			cnode0 = mesh->cell2node[icell][inode];
			cnode1 = mesh->cell2node[icell][(inode+1) % n_neigh];

			for( int jface=0; jface < n_neigh; jface++ ){

				// Skip same faces
				if( iface == jface ) continue;

				// Face nodes obtained from face
				gface = mesh->cell2face[icell][jface];
				if( gface == -1 ){
					cout << "wrong gface =" << setw(10) << mesh->cells[icell].id;
					for( int idim=0; idim < n_neigh; idim++ ){
						cout << setw(10) << mesh->cell2node[icell][idim];
					}
					cout << endl;
					MPI_Abort( MPI_COMM_WORLD, 48984 );
				}
				fnode0 = mesh->face2node[gface][0];
				fnode1 = mesh->face2node[gface][1];

				// Switch faces
				if( (cnode0==fnode0 && cnode1==fnode1) || (cnode0==fnode1 && cnode1==fnode0) ){

					int tmp_face = mesh->cell2face[icell][iface];

					mesh->cell2face[icell][iface] = mesh->cell2face[icell][jface];
					mesh->cell2face[icell][jface] = tmp_face;
				}
			}
		}
	}
}

//***************************************************************************************************
// Set boundary conditions types to all faces, including periodic faces (if exist)
//***************************************************************************************************
void set_boundary_conditions( int n_neigh,
							  gmsh_mesh *mesh,
							  std::vector<std::vector<int>> nodes_bound,
							  std::vector<int> &periodic_faces,
							  std::vector<std::vector<gmsh_neighbor>> &gmsh_neighs ){

	int face_id;

	for( size_t icell=0; icell < mesh->cells.size(); icell++ ){
		for( int iface=0; iface < n_neigh; iface++ ){

			face_id = mesh->cell2face[icell][iface];

			if( periodic_faces[face_id] != -1 ){
				gmsh_neighs[icell][iface].id = mesh->face2cell[periodic_faces[face_id]][0];
			}

			if( mesh->face2cell[face_id][1] != -1 ){

				gmsh_neighs[icell][iface].type = CELL_NONE;

				if( mesh->cells[icell].id != mesh->face2cell[face_id][0] ){
					gmsh_neighs[icell][iface].id = mesh->face2cell[face_id][0];
				}else{
					gmsh_neighs[icell][iface].id = mesh->face2cell[face_id][1];
				}
			}else{
				int fnode0 = mesh->face2node[face_id][0];
				int fnode1 = mesh->face2node[face_id][1];

				for( size_t ib=0; ib < nodes_bound[fnode0].size(); ib++ ) {
					for( size_t jb=0; jb < nodes_bound[fnode1].size(); jb++ ) {

						if( nodes_bound[fnode0][ib] == nodes_bound[fnode1][jb] ){
							gmsh_neighs[icell][iface].type = nodes_bound[fnode1][jb];
							ib = nodes_bound[fnode0].size();
							jb = nodes_bound[fnode1].size();
						}
					}
				}
			}
		}
	}
}

//***************************************************************************************************
//
//***************************************************************************************************
void set_cell_neighbors( gmsh_mesh *mesh, std::vector<std::vector<gmsh_neighbor>> gmsh_neighs ){

	// Allocate connectivity arrays
	std::vector<std::vector<int>> cell2neighTYPE( mesh->cells.size(), vector<int> ( mesh->faces_per_cell, -1 ) );
	std::vector<std::vector<int>> cell2neighID  ( mesh->cells.size(), vector<int> ( mesh->faces_per_cell, -1 ) );

	// Cells
	for( size_t icell=0; icell < mesh->cells.size(); icell++ ){
		for ( int iface=0; iface < mesh->faces_per_cell; iface++ ) {
			cell2neighTYPE[icell][iface] = gmsh_neighs[icell][iface].type;
			cell2neighID  [icell][iface] = gmsh_neighs[icell][iface].id;
		}
	}

	// Allocate neighbors
	mesh->cell2neigh.resize( mesh->cells.size() );
	for( size_t icell=0; icell < mesh->cells.size(); icell++ ){
		mesh->cell2neigh[icell].resize( mesh->faces_per_cell );
	}

	// Set neighbors
	for( size_t icell=0; icell < mesh->cells.size(); icell++ ){
		for( int iface=0; iface < mesh->faces_per_cell; iface++ ){

			// Get ID of neighbors
			int face_id = mesh->cell2face[icell][iface];
			int cell_id = mesh->cells[icell].id;

			int adjc_id = ( mesh->face2cell[face_id][0] != cell_id ? mesh->face2cell[face_id][0]
																   : mesh->face2cell[face_id][1] );
			int jface = 0;

			// Interior neighbor
			if( adjc_id > -1 ){
				mesh->cell2neigh[icell][iface].type = CELL_REGULAR;
				mesh->cell2neigh[icell][iface].id   = adjc_id;
				mesh->cell2neigh[icell][iface].sm   = 1000000;

			// Boundary neighbor
			}else{

				mesh->faces[face_id].is_bface = true;

				// Periodic BC
				if( cell2neighTYPE[icell][iface] != 0 && cell2neighID[icell][iface] != -1 ){
					adjc_id = cell2neighID[icell][iface];

					for( jface=0; jface < mesh->faces_per_cell; jface++ ){
						if( cell_id == cell2neighID[adjc_id][jface] )
							break;
					}

					mesh->cell2neigh[icell][iface].type = CELL_PERIODIC;
					mesh->cell2neigh[icell][iface].id   = adjc_id;
					mesh->cell2neigh[icell][iface].sm   = 100000;

				// Other BCs
				}else{
					mesh->cell2neigh[icell][iface].type = CELL_NONE;
					mesh->cell2neigh[icell][iface].id   = 100000;
					mesh->cell2neigh[icell][iface].sm   = 100000;
				}
				mesh->cell2neigh[icell][iface].bound = mesh->faces[face_id].bound;
			}
		}
	}
}

//***************************************************************************************************
//***************************************************************************************************
void compute_mesh_parameters( MPI_env mpi_env, gmsh_mesh *mesh ){

	// Mesh center
	double grid_xmin[2];
	double grid_xmax[2];

	grid_xmin[0] = mesh->nodes[0].xn[0]; grid_xmax[0] = mesh->nodes[0].xn[0];
	grid_xmin[1] = mesh->nodes[0].xn[1]; grid_xmax[1] = mesh->nodes[0].xn[1];

	for( size_t inode=0; inode < mesh->nodes.size(); inode++ ){

		if( grid_xmin[0] > mesh->nodes[inode].xn[0] )
			grid_xmin[0] = mesh->nodes[inode].xn[0];

		if( grid_xmin[1] > mesh->nodes[inode].xn[1] )
			grid_xmin[1] = mesh->nodes[inode].xn[1];

		if( grid_xmax[0] < mesh->nodes[inode].xn[0] )
			grid_xmax[0] = mesh->nodes[inode].xn[0];

		if( grid_xmax[1] < mesh->nodes[inode].xn[1] )
			grid_xmax[1] = mesh->nodes[inode].xn[1];
	}
	mesh->x[0] = (grid_xmax[0] + grid_xmin[0]) / 2;
	mesh->x[1] = (grid_xmax[1] + grid_xmin[1]) / 2;

	// Cell centers
	for( size_t icell=0; icell < mesh->cells.size(); icell++ ){

		mesh->cells[icell].xc[0] = 0;
		mesh->cells[icell].xc[1] = 0;

		for( int inode=0; inode < mesh->nodes_per_cell; inode++ ){
			int node_id = mesh->cell2node[icell][inode];

			mesh->cells[icell].xc[0] += mesh->nodes[node_id].xn[0];
			mesh->cells[icell].xc[1] += mesh->nodes[node_id].xn[1];
		}
		mesh->cells[icell].xc[0] /= mesh->nodes_per_cell;
		mesh->cells[icell].xc[1] /= mesh->nodes_per_cell;
	}

	// Mesh center
	mesh->x[0] = 0;
	mesh->x[1] = 0;
	for( size_t icell=0; icell < mesh->cells.size(); icell++ ){
		mesh->x[0] += mesh->cells[icell].xc[0];
		mesh->x[1] += mesh->cells[icell].xc[1];
	}
	mesh->x[0] = mesh->x[0] / (double) mesh->cells.size();
	mesh->x[1] = mesh->x[1] / (double) mesh->cells.size();
}

//***************************************************************************************************
void init_mesh_size_mpi_type( MPI_env &mpi_env ){

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	const int mesh_size_nfields = 7;
	MPI_Aint  mesh_size_disps[mesh_size_nfields];
	int       mesh_size_blocklens[] = { 1, 1, 1, 1, 1, 1, 1 };

	MPI_Datatype mesh_size_types[] = { MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT };

	struct t_mesh_size tmp_mesh_size;

	mesh_size_disps[0] = (char*) &tmp_mesh_size.n_nodes        - (char*) &tmp_mesh_size;
	mesh_size_disps[1] = (char*) &tmp_mesh_size.n_faces        - (char*) &tmp_mesh_size;
	mesh_size_disps[2] = (char*) &tmp_mesh_size.n_cells        - (char*) &tmp_mesh_size;
	mesh_size_disps[3] = (char*) &tmp_mesh_size.nodes_per_face - (char*) &tmp_mesh_size;
	mesh_size_disps[4] = (char*) &tmp_mesh_size.cells_per_face - (char*) &tmp_mesh_size;
	mesh_size_disps[5] = (char*) &tmp_mesh_size.nodes_per_cell - (char*) &tmp_mesh_size;
	mesh_size_disps[6] = (char*) &tmp_mesh_size.faces_per_cell - (char*) &tmp_mesh_size;

	MPI_Datatype type_mesh_size;
	MPI_Type_create_struct( mesh_size_nfields, mesh_size_blocklens, mesh_size_disps, mesh_size_types, &type_mesh_size );
	MPI_Type_commit( &type_mesh_size );

	mpi_env.set_mesh_size_type( type_mesh_size );
}

//***************************************************************************************************
MPI_Datatype get_mpitype_meshnode( gmsh_mesh *mesh ){

	//
	const int node_nfields = 3;
	MPI_Aint  node_disps[node_nfields];
	int       node_blocklens[] = { 1, 1, 3 };

	MPI_Datatype node_types[] = { MPI_INT, MPI_INT, MPI_DOUBLE };

	//
	gmsh_node tmp_node;

	node_disps[0] = (char*) &tmp_node.id        - (char*) &tmp_node;
	node_disps[1] = (char*) &tmp_node.id_global - (char*) &tmp_node;
	node_disps[2] = (char*) &tmp_node.xn        - (char*) &tmp_node;

	MPI_Aint total_memory_disp = (char*) &mesh->nodes[1] - (char*) &mesh->nodes[0];

	//
	MPI_Datatype tmp_type, type_node;
	MPI_Type_create_struct ( node_nfields, node_blocklens, node_disps, node_types, &tmp_type );
	MPI_Type_create_resized( tmp_type, 0, total_memory_disp, &type_node );
	MPI_Type_commit( &type_node );

	return type_node;
}

//***************************************************************************************************
MPI_Datatype get_mpitype_meshface( gmsh_mesh *mesh ){

	// fastmesh_periodic_face
	const int periodic_face_nfields = 5;
	MPI_Aint  periodic_face_disps[periodic_face_nfields];
	int       periodic_face_blocklens[] = { 1, 1, 1, 1, 1 };

	MPI_Datatype periodic_face_types[] = { MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT };

	struct fastmesh_periodic_face tmp_periodic_face;

	periodic_face_disps[0] = (char*) &tmp_periodic_face.type    - (char*) &tmp_periodic_face;
	periodic_face_disps[1] = (char*) &tmp_periodic_face.orig_id - (char*) &tmp_periodic_face;
	periodic_face_disps[2] = (char*) &tmp_periodic_face.index   - (char*) &tmp_periodic_face;
	periodic_face_disps[3] = (char*) &tmp_periodic_face.rank    - (char*) &tmp_periodic_face;
	periodic_face_disps[4] = (char*) &tmp_periodic_face.cell_id - (char*) &tmp_periodic_face;

	MPI_Datatype type_periodic_face;
	MPI_Type_create_struct ( periodic_face_nfields, periodic_face_blocklens, periodic_face_disps,
							 periodic_face_types, &type_periodic_face );
	MPI_Type_commit( &type_periodic_face );

	// MeshFace
	const int face_nfields = 7;
	MPI_Aint  face_disps[face_nfields];
	int       face_blocklens[] = { 1, 1, 3, 1, 1, 1, 1 };

	MPI_Datatype face_types[] = { MPI_INT, MPI_INT, MPI_DOUBLE, MPI_INT, MPI_LOGICAL, MPI_LOGICAL,
								  type_periodic_face  };

	gmsh_face tmp_face;

	face_disps[0] = (char*) &tmp_face.id          - (char*) &tmp_face;
	face_disps[1] = (char*) &tmp_face.id_global   - (char*) &tmp_face;
	face_disps[2] = (char*) &tmp_face.xf          - (char*) &tmp_face;
	face_disps[3] = (char*) &tmp_face.bound       - (char*) &tmp_face;
	face_disps[4] = (char*) &tmp_face.is_bface    - (char*) &tmp_face;
	face_disps[5] = (char*) &tmp_face.is_periodic - (char*) &tmp_face;
	face_disps[6] = (char*) &tmp_face.twin_face   - (char*) &tmp_face;

	MPI_Aint total_memory_disp = (char*) &mesh->faces[1] - (char*) &mesh->faces[0];

	MPI_Datatype tmp_type, type_face;
	MPI_Type_create_struct ( face_nfields, face_blocklens, face_disps, face_types, &tmp_type );
	MPI_Type_create_resized( tmp_type, 0, total_memory_disp, &type_face );
	MPI_Type_commit( &type_face );

	return type_face;
}

//***************************************************************************************************
MPI_Datatype get_mpitype_meshcell( gmsh_mesh *mesh ){

	// mesh cell
	const int cell_nfields = 4;
	MPI_Aint  cell_disps[cell_nfields];
	int       cell_blocklens[] = { 1, 1, 3, 1 };

	MPI_Datatype cell_types[] = { MPI_INT, MPI_INT, MPI_DOUBLE, MPI_LOGICAL };

	gmsh_cell tmp_cell;
	cell_disps[0] = (char*) &tmp_cell.id        - (char*) &tmp_cell;
	cell_disps[1] = (char*) &tmp_cell.id_global - (char*) &tmp_cell;
	cell_disps[2] = (char*) &tmp_cell.xc        - (char*) &tmp_cell;
	cell_disps[3] = (char*) &tmp_cell.is_bcell  - (char*) &tmp_cell;

	MPI_Aint total_memory_disp = (char*) &mesh->cells[1] - (char*) &mesh->cells[0];

	MPI_Datatype tmp_type, type_cell;
	MPI_Type_create_struct ( cell_nfields, cell_blocklens, cell_disps, cell_types, &tmp_type );
	MPI_Type_create_resized( tmp_type, 0, total_memory_disp, &type_cell );
	MPI_Type_commit( &type_cell );

	return type_cell;
}

//***************************************************************************************************
MPI_Datatype get_mpitype_meshneigh( gmsh_mesh *mesh ){

	// mesh neigh
	const int neigh_nfields = 5;
	MPI_Aint  neigh_disps[neigh_nfields];
	int       neigh_blocklens[] = { 1, 1, 1, 1, 1 };

	MPI_Datatype neigh_types[] = { MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT };

	t_gmsh_neighbor tmp_neigh;
	neigh_disps[0] = (char*) &tmp_neigh.type      - (char*) &tmp_neigh;
	neigh_disps[1] = (char*) &tmp_neigh.id        - (char*) &tmp_neigh;
	neigh_disps[2] = (char*) &tmp_neigh.sm        - (char*) &tmp_neigh;
	neigh_disps[3] = (char*) &tmp_neigh.proc      - (char*) &tmp_neigh;
	neigh_disps[4] = (char*) &tmp_neigh.bound     - (char*) &tmp_neigh;

	MPI_Datatype type_neigh;
	MPI_Type_create_struct ( neigh_nfields, neigh_blocklens, neigh_disps, neigh_types, &type_neigh );
	MPI_Type_commit( &type_neigh );

	return type_neigh;
}















































