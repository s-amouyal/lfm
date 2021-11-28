#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <map>
#include <math.h>

#include "api/mesh_partitioning.h"
#include "api/mesh_writer.h"
#include "api/hpath.h"

using namespace std;

std::vector<std::vector<int>> new2old;
std::vector<std::vector<int>> old2new_nodes;
std::vector<std::vector<int>> old2new_faces;

//***************************************************************************************************
// Initialize metis' mesh patitions
// 		- only function to interact directly with libmetis
//***************************************************************************************************
gmsh_mesh partition_mesh_metis_mpi( gmsh_mesh *gmsh_domain, MPI_env mpi_env, std::vector<int> &mpi2global  ){

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Perform first level of mesh partitioning: one partition per MPI rank
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	std::vector<int> partitions_mesh;
	if( mpi_env.is_master() ){
		metis_partition_mesh( gmsh_domain, mpi_env.size(), partitions_mesh );
	}else{
		partitions_mesh.resize( gmsh_domain->cells.size() );
	}
	MPI_Bcast( &partitions_mesh[0], partitions_mesh.size(), MPI_INT, mpi_env.get_master(), MPI_COMM_WORLD );
	mpi_env.barrier();

	// Build new local t_mesh
	std::vector<int> old2parts;
	gmsh_mesh gmsh_local;

	gmsh_local = init_submesh_partitions_mpi( gmsh_domain, mpi_env, partitions_mesh, old2parts );

	mpi2global.resize( gmsh_local.cells.size(), -1 );
	for( size_t icell=0; icell < old2parts.size(); icell++ ){
		if( old2parts[icell] == -1 )
			continue;

		mpi2global[old2parts[icell]] = icell;
	}

	return gmsh_local;
}

//***************************************************************************************************
// Initialize metis' mesh patitions
// 		- only function to interact directly with libmetis
// 		- used both for the MPI and submesh partitioning
//***************************************************************************************************
void metis_partition_mesh( gmsh_mesh *gmsh, int nParts, std::vector<int> &partitions ){

	// Initialize partitions vector and handle 1 partition case
	partitions.resize( gmsh->cells.size(), -1 );

	if( nParts == 1 ){
		for( size_t icell=0; icell < gmsh->cells.size(); icell++ ){
			partitions[icell] = 0;
		}
		return;
	}

	// Count the number of metis "edges": the total number of neighbors for all cells
	int total_neighbors = 0;

	for( size_t icell=0; icell < gmsh->cells.size(); icell++ ){
		for( int iface=0; iface < gmsh->faces_per_cell; iface++ ){
			if( gmsh->cell2neigh[icell][iface].type == CELL_REGULAR ){
				total_neighbors++;
			}
		}
	}

	// Setup metis variables
	idx_t nVertices = (int) gmsh->cells.size();
	idx_t nEdges    = total_neighbors;
	idx_t nWeights  = 1;
	idx_t m_nParts  = (idx_t) nParts;

	idx_t objval;
	idx_t *part = new idx_t [nVertices];
	idx_t options[METIS_NOPTIONS];

	// xadj
	idx_t *xadj = new idx_t [nVertices+1];
	for( int icell=0; icell < nVertices+1; icell++ ){
		xadj[icell] = -1;
	}

	xadj[0] = 0;
	idx_t m_count = 0;

	for( size_t icell=0; icell < gmsh->cells.size(); icell++ ){
		for( int iface=0; iface < gmsh->faces_per_cell; iface++ ){
			if( gmsh->cell2neigh[icell][iface].type == CELL_REGULAR ){
				m_count++;
			}

		}
		xadj[icell+1] = m_count;
	}

	// adjncy
	idx_t *adjncy = new idx_t [ 2*nEdges ];
	for( int icell=0; icell < 2*nEdges; icell++ ){
		adjncy[icell] = -1;
	}

	idx_t m_ind = 0;

	int myrank;
	MPI_Comm_rank( MPI_COMM_WORLD, &myrank );

	for( size_t icell=0; icell < gmsh->cells.size(); icell++ ){
		for( int iface=0; iface < gmsh->faces_per_cell; iface++ ){
			if( gmsh->cell2neigh[icell][iface].type == CELL_REGULAR ){
				adjncy[m_ind++] = gmsh->cell2neigh[icell][iface].id ;
			}
		}
	}
	for( int icell=0; icell < nEdges; icell++ ){
		adjncy[icell+nEdges] = adjncy[icell];
	}

	// Options
	METIS_SetDefaultOptions( options );
	options[METIS_OPTION_CONTIG] = 1;

	// Partition mesh
	METIS_PartGraphKway( &nVertices, &nWeights, xadj, adjncy, NULL, NULL, NULL, &m_nParts, NULL,
						 NULL, options, &objval, part );

	// Fill in output vector
	for( size_t icell=0; icell < (size_t) nVertices; icell++ ){
		partitions[icell] = (int) part[icell];
	}

	delete[] part;
	delete[] adjncy;
	delete[] xadj;
}

//***************************************************************************************************
// Initialize metis' mesh patitions
// 		- only function to interact directly with libmetis
// 		- used both for the MPI and submesh partitioning
//***************************************************************************************************
void metis_partition_mesh_new( gmsh_mesh *gmsh, int nParts, std::vector<int> &partitions ){

	// Initialize partitions vector and handle 1 partition case
	partitions.resize( gmsh->cells.size(), -1 );

	if( nParts == 1 ){
		for( size_t icell=0; icell < gmsh->cells.size(); icell++ ){
			partitions[icell] = 0;
		}
		return;
	}

	// Count the number of metis "edges": the total number of neighbors for all cells
	int total_neighbors = 0;

	for( size_t icell=0; icell < gmsh->cells.size(); icell++ ){
		for( int iface=0; iface < gmsh->faces_per_cell; iface++ ){
			if( gmsh->cell2neigh[icell][iface].type == CELL_REGULAR ){
				total_neighbors++;
			}
		}
	}

	// Setup metis variables
	idx_t nVertices = (int) gmsh->cells.size();
	idx_t nEdges    = total_neighbors;
	idx_t nWeights  = 1;
	idx_t m_nParts  = (idx_t) nParts;

	idx_t objval;
	idx_t *part = new idx_t [nVertices];
	idx_t options[METIS_NOPTIONS];

	// xadj
	idx_t *xadj = new idx_t [nVertices+1];
	for( int icell=0; icell < nVertices+1; icell++ ){
		xadj[icell] = -1;
	}

	xadj[0] = 0;
	idx_t m_count = 0;

	for( size_t icell=0; icell < gmsh->cells.size(); icell++ ){
		for( int iface=0; iface < gmsh->faces_per_cell; iface++ ){
			if( gmsh->cell2neigh[icell][iface].type == CELL_REGULAR ){
				m_count++;
			}

		}
		xadj[icell+1] = m_count;
	}

	// adjncy
	idx_t *adjncy = new idx_t [ 2*nEdges ];
	for( int icell=0; icell < 2*nEdges; icell++ ){
		adjncy[icell] = -1;
	}

	idx_t m_ind = 0;

	int myrank;
	MPI_Comm_rank( MPI_COMM_WORLD, &myrank );

	for( size_t icell=0; icell < gmsh->cells.size(); icell++ ){
		for( int iface=0; iface < gmsh->faces_per_cell; iface++ ){
			if( gmsh->cell2neigh[icell][iface].type == CELL_REGULAR ){
				adjncy[m_ind++] = gmsh->cell2neigh[icell][iface].id ;
			}
		}
	}
	for( int icell=0; icell < nEdges; icell++ ){
		adjncy[icell+nEdges] = adjncy[icell];
	}

	// Options
	METIS_SetDefaultOptions( options );
	options[METIS_OPTION_CONTIG] = 1;

	// Partition mesh
	METIS_PartGraphKway( &nVertices, &nWeights, xadj, adjncy, NULL, NULL, NULL, &m_nParts, NULL,
						 NULL, options, &objval, part );

	// Fill in output vector
	for( size_t icell=0; icell < (size_t) nVertices; icell++ ){
		partitions[icell] = (int) part[icell];
	}

	delete[] part;
	delete[] adjncy;
	delete[] xadj;
}

//***************************************************************************************************
// Initialize submesh partitions
//***************************************************************************************************
gmsh_mesh init_submesh_partitions_mpi( gmsh_mesh *gmsh, MPI_env mpi_env, std::vector<int> partitions,
									   std::vector<int> &old2parts ){

	// Compute the number of cells for each mesh
	int local_ncells = 0;
	for( size_t icell=0; icell < gmsh->cells.size(); icell++ ){
		if( partitions[icell] == mpi_env.rank() ){
			local_ncells++;
		}
	}
	std::vector<int> new2old( local_ncells );

	// Create one new mesh for each interior partition
	gmsh_mesh gmsh_local;

	old2parts.resize( gmsh->cells.size(), -1 );

	gmsh_local.cells.resize( local_ncells );

	gmsh_local.cell2neigh.resize( local_ncells );
	for( int icell=0; icell < local_ncells; icell++ ){
		gmsh_local.cell2neigh[icell].resize( gmsh->faces_per_cell );
	}

	gmsh_local.id             = mpi_env.rank();
	gmsh_local.faces_per_cell = gmsh->faces_per_cell;
	gmsh_local.nodes_per_cell = gmsh->nodes_per_cell;
	gmsh_local.cells_per_face = gmsh->cells_per_face;
	gmsh_local.nodes_per_face = gmsh->nodes_per_face;
	gmsh_local.elmnt_type     = gmsh->elmnt_type;

	gmsh_local.cell2node.resize( local_ncells, std::vector<int> (gmsh_local.nodes_per_cell, -1));
	gmsh_local.cell2face.resize( local_ncells, std::vector<int> (gmsh_local.faces_per_cell, -1));

	// Initialize id and new-old indices connectivity
	int local_ind = 0;
	for( size_t icell=0; icell < gmsh->cells.size(); icell++ ){
		if( partitions[icell] == mpi_env.rank() ){
			gmsh_local.cells[local_ind].id = local_ind;
			new2old[local_ind] = gmsh->cells[icell].id ;

			old2parts[icell] = local_ind++;
		}
	}
	std::vector<int> global_old2parts( gmsh->cells.size(), -1 );
	{
		std::vector<int> local_ind( mpi_env.size(), 0 );
		for( size_t icell=0; icell < gmsh->cells.size(); icell++ ){
			global_old2parts[icell] = local_ind[partitions[icell]]++;
		}
	}

	// Initialize cells
	int gcell;
	for( size_t icell=0; icell < gmsh_local.cells.size(); icell++ ){

		gcell = new2old[icell];

		gmsh_local.cells[icell].xc[0] = gmsh->cells[gcell].xc[0]  ;
		gmsh_local.cells[icell].xc[1] = gmsh->cells[gcell].xc[1]  ;

		for( int iface=0; iface < gmsh_local.faces_per_cell; iface++ ){

			t_gmsh_neighbor *neigh = &(gmsh->cell2neigh[gcell][iface]) ;

			// Cell is out of domain boundary
			if( neigh->type == CELL_NONE ){
				gmsh_local.cell2neigh[icell][iface].type = CELL_NONE;
				gmsh_local.cell2neigh[icell][iface].id   = -1;
				gmsh_local.cell2neigh[icell][iface].sm   = -1;
				gmsh_local.cell2neigh[icell][iface].proc = -1;
			// Cell belongs to the current MPI mesh
			}else if( partitions[neigh->id] == mpi_env.rank() && neigh->type != CELL_PERIODIC ){
				gmsh_local.cell2neigh[icell][iface].type = CELL_REGULAR;
				gmsh_local.cell2neigh[icell][iface].id   = old2parts[neigh->id];
				gmsh_local.cell2neigh[icell][iface].sm   = -1;
				gmsh_local.cell2neigh[icell][iface].proc = mpi_env.rank();
			// Cell is on a neighboring processor
			}else if( neigh->type != CELL_PERIODIC ){
				gmsh_local.cell2neigh[icell][iface].type = CELL_MPI;
				gmsh_local.cell2neigh[icell][iface].id   = global_old2parts[neigh->id];
				gmsh_local.cell2neigh[icell][iface].sm   = INDEX_BND_SUBMESH;
				gmsh_local.cell2neigh[icell][iface].proc = partitions[neigh->id];
			}

			if( neigh->type == CELL_PERIODIC ){
				if( partitions[neigh->id] == mpi_env.rank() ){
					gmsh_local.cell2neigh[icell][iface].type = CELL_PERIODIC;
					gmsh_local.cell2neigh[icell][iface].id   = old2parts[neigh->id];
					gmsh_local.cell2neigh[icell][iface].sm   = INDEX_BND_SUBMESH;
					gmsh_local.cell2neigh[icell][iface].proc = mpi_env.rank();
				}else{
					gmsh_local.cell2neigh[icell][iface].type = CELL_PERIODIC_MPI;
					gmsh_local.cell2neigh[icell][iface].id   = global_old2parts[neigh->id];
					gmsh_local.cell2neigh[icell][iface].sm   = INDEX_BND_SUBMESH;
					gmsh_local.cell2neigh[icell][iface].proc = partitions[neigh->id];
				}
			}
		}
	}

	// Initialize nodes
	int gnode;
	std::vector<bool> unregistered_nodes( gmsh->nodes.size(), true );

	int n_nodes = 0;

	for( size_t icell=0; icell < gmsh_local.cells.size(); icell++ ){
		for( int inode=0; inode < gmsh->nodes_per_cell; inode++ ){
			gcell = new2old[icell];
			gnode = gmsh->cell2node[gcell][inode];

			if( unregistered_nodes[gnode] ){
				unregistered_nodes[gnode] = false;
				n_nodes++;
			}
		}
	}

	gmsh_local.nodes.resize( n_nodes );

	for( size_t inode=0; inode < gmsh->nodes.size(); inode++ ){
		unregistered_nodes[inode] = true;
	}

	int node_ind = 0;
	std::vector<int> old2new_nodes( gmsh->nodes.size(), -1 );

	for( size_t icell=0; icell < gmsh_local.cells.size(); icell++ ){

		for( int inode=0; inode < gmsh->nodes_per_cell; inode++ ){
			gcell = new2old[icell];
			gnode = gmsh->cell2node[gcell][inode];

			if( unregistered_nodes[gnode] ){
				unregistered_nodes[gnode] = false;

				gmsh_local.nodes[node_ind].id = gmsh->nodes[gnode].id ;
				gmsh_local.nodes[node_ind].xn[0] = gmsh->nodes[gnode].xn[0]  ;
				gmsh_local.nodes[node_ind].xn[1] = gmsh->nodes[gnode].xn[1]  ;
				gmsh_local.nodes[node_ind].xn[2] = gmsh->nodes[gnode].xn[2]  ;

				gmsh_local.nodes[node_ind].id_global = gmsh->nodes[gnode].id;

				old2new_nodes[gnode] = node_ind;

				node_ind++;
			}
		}
	}

	for( size_t icell=0; icell < gmsh_local.cells.size(); icell++ ){
		for( int inode=0; inode < gmsh->nodes_per_cell; inode++ ){
			gcell = new2old[icell];
			gnode = gmsh->cell2node[gcell][inode];

			gmsh_local.cell2node[icell][inode] = old2new_nodes[gmsh->cell2node[gcell][inode]];
		}
	}

	// Initialize faces
	int gface;
	std::vector<bool> unregistered_faces( gmsh->faces.size(), true );

	int n_faces = 0;

	for( size_t icell=0; icell < gmsh_local.cells.size(); icell++ ){
		for( int iface=0; iface < gmsh->faces_per_cell; iface++ ){
			gcell = new2old[icell];
			gface = gmsh->cell2face[gcell][iface];

			if( unregistered_faces[gface] ){
				unregistered_faces[gface] = false;
				n_faces++;
			}
		}
	}

	gmsh_local.faces.resize( n_faces );

	for( size_t iface=0; iface < gmsh->faces.size(); iface++ ){
		unregistered_faces[iface] = true;
	}

	int face_ind = 0;
	std::vector<int> old2new_faces( gmsh->faces.size(), -1 );
	std::vector<int> new2old_faces( gmsh_local.faces.size(), -1 );

	for( size_t icell=0; icell < gmsh_local.cells.size(); icell++ ){

		for( int iface=0; iface < gmsh->faces_per_cell; iface++ ){
			gcell = new2old[icell];
			gface = gmsh->cell2face[gcell][iface];

			if( unregistered_faces[gface] ){
				unregistered_faces[gface] = false;

				gmsh_local.faces[face_ind].id = face_ind;
				old2new_faces[gface] = face_ind;
				new2old_faces[face_ind] = gface;

				face_ind++;
			}
		}
	}

	// cell2face
	for( size_t icell=0; icell < gmsh_local.cells.size(); icell++ ){
		for( int iface=0; iface < gmsh->faces_per_cell; iface++ ){
			gcell = new2old[icell];

			gmsh_local.cell2face[icell][iface] = old2new_faces[gmsh->cell2face[gcell][iface]];
		}
	}

	// face2cell
	gmsh_local.face2cell.resize( gmsh_local.faces.size(),
								std::vector<int> ( gmsh_local.cells_per_face, -1 ) );

	for( size_t icell=0; icell < gmsh_local.cells.size(); icell++ ){
		for( int iface=0; iface < gmsh_local.faces_per_cell; iface++ ){
			gface = gmsh_local.cell2face[icell][iface];

			// If both cells are registered, job is already done
			if( gmsh_local.face2cell[gface][0] != -1 &&
				gmsh_local.face2cell[gface][1] != -1 ){
				continue;
			}

			//
			if( gmsh_local.face2cell[gface][0] == -1 ){
				gmsh_local.face2cell[gface][0] = icell;
			}else{
				gmsh_local.face2cell[gface][1] = icell;
			}
		}
	}

	// face2node
	int n0, n1, gface_new;

	gmsh_local.face2node.resize( gmsh_local.faces.size(),
								 std::vector<int> ( gmsh_local.nodes_per_face, -1 ) );

	for( size_t icell=0; icell < gmsh_local.cells.size(); icell++ ){
		for( int iface=0; iface < gmsh->faces_per_cell; iface++ ){
			gcell = new2old[icell];

			n0 = gmsh->cell2node[gcell][ iface ];
			n1 = gmsh->cell2node[gcell][ (iface+1) % gmsh->faces_per_cell ];

			gface_new = gmsh_local.cell2face[icell][iface];

			gmsh_local.face2node[gface_new][0] = old2new_nodes[n0];
			gmsh_local.face2node[gface_new][1] = old2new_nodes[n1];
		}
	}

	// Compute mesh center
	gmsh_local.x[0] = 0.;
	gmsh_local.x[1] = 0.;

	for( size_t icell=0; icell < gmsh_local.cells.size(); icell++ ){
		gmsh_local.x[0] += gmsh_local.cells[icell].xc[0] ;
		gmsh_local.x[1] += gmsh_local.cells[icell].xc[1] ;
	}

	gmsh_local.x[0] = gmsh_local.x[0] / (double) gmsh_local.cells.size() ;
	gmsh_local.x[1] = gmsh_local.x[1] / (double) gmsh_local.cells.size() ;

	// Submesh connectivity
	if( mpi_env.rank() == 1 )
		cout << "Test 100, bcell" << endl;

	for( size_t icell=0; icell < gmsh_local.cells.size(); icell++ ){
		for( int iface=0; iface < gmsh_local.faces_per_cell; iface++ ){
			if( gmsh_local.cell2neigh[icell][iface].type != CELL_REGULAR ){
				if( mpi_env.rank() == 1 ){
					cout << "bcell =" << setw(10) << icell << endl;
				}
				gmsh_local.cells[icell].is_bcell = true;
				break;
			}
		}
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// MPI connectivity
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	std::vector<bool> unregistered_mpi_neigh( mpi_env.size(), true );
	for( size_t icell=0; icell < gmsh_local.cells.size(); icell++ ){
		for( int iface=0; iface < gmsh_local.faces_per_cell; iface++ ){
			t_gmsh_neighbor *neigh_cell = &( gmsh_local.cell2neigh[icell][iface] );
			if( neigh_cell->type == CELL_MPI  || neigh_cell->type == CELL_PERIODIC ||
				neigh_cell->type == CELL_NONE || neigh_cell->type == CELL_PERIODIC_MPI ){
				if( unregistered_mpi_neigh[neigh_cell->proc] ){
					unregistered_mpi_neigh[neigh_cell->proc] = false;

					gmsh_local.mpi_neighbors.push_back( neigh_cell->proc );
				}
			}
		}
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Periodic faces
	//		- need to gather new faces ID local to each MPI rank
	//		- Need to update twin periodic faces IDs when they are on different MPI ranks
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	std::vector<int> faces_per_rank( mpi_env.size(), -1 );
	int num_faces  = (int) gmsh_local.faces.size();
	MPI_Allgather( &num_faces, 1, MPI_INT, &faces_per_rank[0], 1, MPI_INT, MPI_COMM_WORLD );

	// Need temporary 1D array since MPI does not deal well with jagged arrays
	std::vector<int> tmp_faces;
	int total_faces = 0;
	for( int irank=0; irank < mpi_env.size(); irank++ ){
		total_faces += faces_per_rank[irank];
	}
	tmp_faces.resize( gmsh->faces.size() * mpi_env.size() );

	std::vector<int> disp( mpi_env.size(), 0 );
	for( int irank=1; irank < mpi_env.size(); irank++ ){
		disp[irank] = faces_per_rank[irank-1];
	}

	// Gather the old2new mapping of each MPI rank
	MPI_Barrier( MPI_COMM_WORLD );
	MPI_Allgather( &old2new_faces[0], gmsh->faces.size(), MPI_INT,
				   &tmp_faces[0],     gmsh->faces.size(), MPI_INT, MPI_COMM_WORLD );

	std::vector<std::vector<int>> old2new_faces_global( mpi_env.size(), std::vector<int>( gmsh->faces.size(), -1 ) );
	int ind = 0;
	for( int irank=0; irank < mpi_env.size(); irank++ ){
		for( size_t iface=0; iface < old2new_faces.size(); iface++ ){
			old2new_faces_global[irank][iface] = tmp_faces[ind++];
		}
	}

	for( size_t iface=0; iface < gmsh_local.faces.size(); iface++ ){
		gmsh_local.faces[iface].twin_face.type    = PERIODIC_NEIGH_NONE;
		gmsh_local.faces[iface].twin_face.rank    = -1;
		gmsh_local.faces[iface].twin_face.index   = -1;
		gmsh_local.faces[iface].twin_face.orig_id = -1;
	}

	for( size_t icell=0; icell < gmsh_local.cells.size(); icell++ ){
		for( int iface=0; iface < gmsh_local.faces_per_cell; iface++ ){

			t_gmsh_neighbor *neigh = &( gmsh_local.cell2neigh[icell][iface] ) ;

			if( neigh->type != CELL_PERIODIC && neigh->type != CELL_PERIODIC_MPI ){
				continue;
			}
			int neigh_proc = -1;
			int gface      = gmsh_local.cell2face[icell][iface];
			int old_face   = new2old_faces[gface];
			int twin_face  = gmsh->faces[old_face].twin_face.index;

			if( neigh->type == CELL_PERIODIC ){
				neigh_proc = mpi_env.rank();
			}
			if( neigh->type == CELL_PERIODIC_MPI ){
				neigh_proc = neigh->proc;
			}
			gmsh_local.faces[gface].twin_face.index = old2new_faces_global[neigh_proc][twin_face];
			gmsh_local.faces[gface].twin_face.rank  = neigh_proc;

			if( neigh->type == CELL_PERIODIC )
				gmsh_local.faces[gface].twin_face.type = PERIODIC_NEIGH_REGULAR;

			if( neigh->type == CELL_PERIODIC_MPI )
				gmsh_local.faces[gface].twin_face.type = PERIODIC_NEIGH_MPI;
		}
	}

	return gmsh_local;
}

//***************************************************************************************************
// Create new mesh only for all the Mesh's interior cells
// 		- only fill the new t_mesh with the information relevant for metis
// 		- n_cells, faces_per_cell, nodes_per_cell, cells.neigh_id
//***************************************************************************************************
gmsh_mesh set_interior_mesh( gmsh_mesh *gmsh_local, MPI_env mpi_env, std::vector<int > &old2new,
							 std::vector<int> &old2new_bnd, std::vector<int> &int2old ){

	int myrank;
	MPI_Comm_rank( MPI_COMM_WORLD, &myrank );

	// Create new interior mesh
	gmsh_mesh gmsh_int;

	int n_cells = gmsh_local->cells.size() - ( gmsh_local->hpath.size() + gmsh_local->deadend_cells.size() );
	gmsh_int.faces_per_cell = gmsh_local->faces_per_cell;
	gmsh_int.nodes_per_cell = gmsh_local->nodes_per_cell;

	gmsh_int.cells     .resize( n_cells );
	gmsh_int.cell2neigh.resize( n_cells );
	for( size_t icell=0; icell < gmsh_int.cells.size(); icell++ ){
		gmsh_int.cell2neigh[icell].resize( gmsh_int.faces_per_cell );
	}

	// Mark the boundary cells
	cout << "num cells =" << setw(10) << mpi_env.rank() << setw(10) << gmsh_local->cells.size() << endl;
	std::vector<bool> is_boundary( gmsh_local->cells.size(), false );
	for( size_t icell=0; icell < gmsh_local->hpath.size(); icell++ ){
		if( gmsh_local->hpath[icell] < 0 || gmsh_local->hpath[icell] >= (int) is_boundary.size() ){
			cout << "error with boundary Hpath =" << setw(10) << myrank
												  << setw(10) << gmsh_local->hpath[icell]
												  << setw(10) << is_boundary.size()
												  << endl;
			for( size_t jcell=0; jcell < gmsh_local->hpath.size(); jcell++ ){
				cout << "hpath =" << setw(10) << jcell
								  << setw(10) << gmsh_local->hpath[jcell]
								  << endl;
			}
			MPI_Abort( MPI_COMM_WORLD, 2873 );
		}
		is_boundary[gmsh_local->hpath[icell]] = true;
	}
	for( size_t icell=0; icell < gmsh_local->deadend_cells.size(); icell++ ){
		is_boundary[gmsh_local->deadend_cells[icell]] = true;
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Index mapping: from original and to new cell ordering (local mesh minus boundary path )
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	// Boundary cells: already have Hpath
	for( size_t icell=0; icell < gmsh_local->hpath.size(); icell++ ){
		old2new    [gmsh_local->hpath[icell]] = icell;
		old2new_bnd[gmsh_local->hpath[icell]] = icell;
	}

	// Interior cels
	int cell_int = 0;
	int2old.resize( gmsh_int.cells.size() );
	for( size_t icell=0; icell < gmsh_local->cells.size(); icell++ ){

		// Skip boundary cells
		if( is_boundary[icell] ){
			continue;
		}

		// Cell ID
		if( cell_int < 0 || cell_int > (int) gmsh_int.cells.size() ){
			cout << "Error in setting interior mesh, #102, rank =" << setw(10) << myrank
																   << setw(10) << cell_int
																   << setw(10) << gmsh_int.cells.size()
																   << setw(10) << gmsh_local->cells.size()
																   << endl;
			MPI_Abort( MPI_COMM_WORLD, 102 );
		}
		gmsh_int.cells[cell_int].id = cell_int;

		// Mesh mapping
		int2old[cell_int] = icell;
		old2new[icell] = cell_int++;
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Setup the cell neighbors connectivity for local boundary/interior meshes
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	// Boundary cells only
	for( size_t icell=0; icell < gmsh_local->cells.size(); icell++ ){

		// Skip interior cells
		if( !is_boundary[icell] ){
			continue;
		}

		for( int iface=0; iface < gmsh_local->faces_per_cell; iface++ ){
			t_gmsh_neighbor *neigh = &( gmsh_local->cell2neigh[icell][iface] );

			// Skip MPI and out-of-boundary neighbors
			if( neigh->type == CELL_MPI || neigh->type == CELL_NONE ||
				neigh->type == CELL_PERIODIC || neigh->type == CELL_PERIODIC_MPI ){
				continue;
			}

			if( is_boundary[neigh->id] ){
				neigh->type = CELL_REGULAR;
				neigh->id   = neigh->id;
				neigh->sm   = 0;
			}else{
				neigh->type = CELL_SUBMESH;
				neigh->id   = neigh->id;
				neigh->sm   = 1;
			}
		}
	}

	// Interior cells
	cell_int = 0;
	for( size_t icell=0; icell < gmsh_local->cells.size(); icell++ ){

		// Skip boundary cells
		if( is_boundary[icell] ){
			continue;
		}

		for( int iface=0; iface < gmsh_local->faces_per_cell; iface++ ){

			t_gmsh_neighbor *old_neigh = &( gmsh_local->cell2neigh[icell][iface] );

			if( old_neigh->id < 0 || old_neigh->id > (int) is_boundary.size() ){
				cout << "Error setting interior mesh, rank =" << setw(10) << myrank << endl;
				MPI_Abort( MPI_COMM_WORLD, 935 );
			}

			if( is_boundary[old_neigh->id] ){
				gmsh_int.cell2neigh[cell_int][iface].type = CELL_SUBMESH;
				gmsh_int.cell2neigh[cell_int][iface].id   = old_neigh->id;
				gmsh_int.cell2neigh[cell_int][iface].sm   = INDEX_BND_SUBMESH;
			}else{
				gmsh_int.cell2neigh[cell_int][iface].type = CELL_REGULAR;
				gmsh_int.cell2neigh[cell_int][iface].id   = old2new[old_neigh->id];
				gmsh_int.cell2neigh[cell_int][iface].sm   = 0;
			}
		}
		cell_int++;
	}

	return gmsh_int;
}

//***************************************************************************************************
// Initialize metis' mesh patitions
// 		- only function to interact directly with libmetis
//***************************************************************************************************
std::vector<gmsh_mesh> partition_local_mesh( gmsh_mesh        *gmsh_local,
													gmsh_mesh        *gmsh_int,
													int               nSubparts,
													std::vector<int>  old2new,
													MPI_env		      mpi_env,
													std::vector<int>  tmp2old,
													std::vector<std::vector<int>> &int2old ){

	int myrank;
	MPI_Comm_rank( MPI_COMM_WORLD, &myrank );

	// Partition interior mesh into Ns submeshes
	std::vector<int>  partitions_tmp;
	metis_partition_mesh( gmsh_int, nSubparts, partitions_tmp );

	// Identify partitions
	std::vector<int>  partitions_core;
	std::vector<bool> is_boundary( gmsh_local->cells.size(), false );

	for( size_t icell=0; icell < gmsh_local->hpath.size(); icell++ ){
		is_boundary[gmsh_local->hpath[icell]] = true;
	}
	for( size_t icell=0; icell < gmsh_local->deadend_cells.size(); icell++ ){
		is_boundary[gmsh_local->deadend_cells[icell]] = true;
	}

	partitions_core.resize( gmsh_local->cells.size(), -1 );
	for( size_t icell=0; icell < partitions_core.size(); icell++ ){
		if( is_boundary[icell] ){
			partitions_core[icell] = INDEX_BND_SUBMESH;
		}else{
			partitions_core[icell] = partitions_tmp[old2new[icell]] + 1;
		}
	}

	// Store submeshes into vector part
	std::vector<gmsh_mesh> gmsh_parts;

	gmsh_parts = init_submesh_partitions( gmsh_local, gmsh_int, nSubparts, partitions_core, int2old );

	// Fill boundary hpath (already found)
	gmsh_parts[INDEX_BND_SUBMESH].hpath         = gmsh_local->hpath ;
	gmsh_parts[INDEX_BND_SUBMESH].deadend_cells = gmsh_local->deadend_cells ;

	// Update MPI connectivity
	update_mesh_mpi_connectivity( gmsh_local, &gmsh_parts[0] );

	// Periodic boundary conditions: only serial for now
	update_periodic_connectivity( gmsh_local, &gmsh_parts[0] );

	return gmsh_parts;
}

//***************************************************************************************************
//***************************************************************************************************
void update_mesh_mpi_connectivity( gmsh_mesh* mesh_local, gmsh_mesh* mesh_bnd ){

	const int million = 10000;
	int nProcs;
	int myrank;
	MPI_Comm_size( MPI_COMM_WORLD, &nProcs );
	MPI_Comm_rank( MPI_COMM_WORLD, &myrank );

	// Index mapping: from original to Hpath cell ordering
	std::vector<int> old2hpath( mesh_local->cells.size(), -1 );
	for( size_t icell=0; icell < mesh_bnd->cells.size(); icell++ ){
		if( icell < mesh_bnd->hpath.size() )
			old2hpath[mesh_bnd->hpath[icell]] = icell;
		else
			old2hpath[mesh_bnd->deadend_cells[icell-mesh_bnd->hpath.size()]] = icell;
	}

	// Index mapping: from MPI rank to local mpi neighbors
	std::vector<int> proc2local( nProcs, -1 );
	for( size_t ineigh=0; ineigh < mesh_bnd->mpi_neighbors.size(); ineigh++ ){
		proc2local[mesh_bnd->mpi_neighbors[ineigh]] = ineigh;
	}

	// Share old2hpath vector with each of the MPI neighbors
	// 		==> let the update of mpi connectivity happen locally
	std::vector<std::vector<int>> old2hpath_mpi( mesh_bnd->mpi_neighbors.size() );

	std::vector<MPI_Request> send_req( mesh_bnd->mpi_neighbors.size() );
	std::vector<MPI_Request> recv_req( mesh_bnd->mpi_neighbors.size() );

	std::vector<int> neigh_size( mesh_bnd->mpi_neighbors );

	MPI_Status *mpi_status = new MPI_Status[ mesh_bnd->mpi_neighbors.size() ];
	memset( mpi_status, 0, sizeof(MPI_Status) * mesh_bnd->mpi_neighbors.size() );

	for( size_t ineigh=0; ineigh < mesh_bnd->mpi_neighbors.size(); ineigh++ ){

		// Prepare communication
		int send_tag = million*myrank + mesh_bnd->mpi_neighbors[ineigh];
		int recv_tag = myrank + million*mesh_bnd->mpi_neighbors[ineigh];

		int src  = mesh_bnd->mpi_neighbors[ineigh];
		int dest = mesh_bnd->mpi_neighbors[ineigh];

		int local_size = mesh_local->cells.size();

		MPI_Isend( &local_size        , 1, MPI_INT, dest, send_tag, MPI_COMM_WORLD, &send_req[ineigh] );
		MPI_Irecv( &neigh_size[ineigh], 1, MPI_INT, src , recv_tag, MPI_COMM_WORLD, &recv_req[ineigh] );
	}

	MPI_Waitall( send_req.size(), &send_req[0], MPI_STATUS_IGNORE );
	MPI_Waitall( recv_req.size(), &recv_req[0], MPI_STATUS_IGNORE );

	for( size_t ineigh=0; ineigh < mesh_bnd->mpi_neighbors.size(); ineigh++ ){

		// Prepare communication
		int send_tag = million*myrank + mesh_bnd->mpi_neighbors[ineigh];
		int recv_tag = myrank + million*mesh_bnd->mpi_neighbors[ineigh];

		int src  = mesh_bnd->mpi_neighbors[ineigh];
		int dest = mesh_bnd->mpi_neighbors[ineigh];

		int local_size = mesh_local->cells.size();


		// Allocate vector with neighbor's size
		old2hpath_mpi[ineigh].resize( neigh_size[ineigh] );

		// Send hpath array to neighbor
		MPI_Isend( &old2hpath            [0], local_size,         MPI_INT, dest, send_tag, MPI_COMM_WORLD, &send_req[ineigh] );
		MPI_Irecv( &old2hpath_mpi[ineigh][0], neigh_size[ineigh], MPI_INT, src , recv_tag, MPI_COMM_WORLD, &recv_req[ineigh] );
	}

	MPI_Waitall( send_req.size(), &send_req[0], mpi_status );
	MPI_Waitall( recv_req.size(), &recv_req[0], mpi_status );

	// Update local connectivity with mpi neighbors
	for( size_t icell=0; icell < mesh_bnd->cells.size(); icell++ ){
		for( int iface=0; iface < mesh_bnd->faces_per_cell; iface++ ){
			t_gmsh_neighbor* neigh = &( mesh_bnd->cell2neigh[icell][iface] );
			if( neigh->type == CELL_MPI || neigh->type == CELL_PERIODIC_MPI ){
				neigh->id = old2hpath_mpi[proc2local[neigh->proc]][neigh->id];
			}
		}
	}
	//MPI_Barrier( MPI_COMM_WORLD );

	delete[] mpi_status;
}

//***************************************************************************************************
// Update periodic's neighbors IDs after boundary Hpath
// 	- only serial for now
//***************************************************************************************************
void update_periodic_connectivity( gmsh_mesh *mesh_local, gmsh_mesh *mesh_bnd ){

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// CELL_PERIODIC: when periodic neighbor sits on same processor
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	// Index mapping: from original to Hpath cell ordering
	std::vector<int> old2hpath( mesh_local->cells.size(), -1 );
	for( size_t icell=0; icell < mesh_bnd->cells.size(); icell++ ){
		if( icell < mesh_bnd->hpath.size() )
			old2hpath[mesh_bnd->hpath[icell]] = icell;
		else
			old2hpath[mesh_bnd->deadend_cells[icell-mesh_bnd->hpath.size()]] = icell;
	}

	// Index mapping: from original to Hpath cell ordering
	int myrank;
	MPI_Comm_rank( MPI_COMM_WORLD, &myrank );
	for( size_t icell=0; icell < mesh_bnd->cells.size(); icell++ ){
		for( int iface=0; iface < mesh_bnd->faces_per_cell; iface++ ){
			t_gmsh_neighbor *neigh = &( mesh_bnd->cell2neigh[icell][iface] );

			if( neigh->type == CELL_PERIODIC ){
				neigh->id = old2hpath[neigh->id];
			}
		}
	}
}

//***************************************************************************************************
// Initialize submesh partitions
//***************************************************************************************************
std::vector<gmsh_mesh> init_submesh_partitions( gmsh_mesh *gmsh, gmsh_mesh *gmsh_int, int nParts,
												std::vector<int> partitions,
												std::vector<std::vector<int>> &int2old ){

	int myrank;
	MPI_Comm_rank( MPI_COMM_WORLD, &myrank );

	// Compute the number of cells for each mesh
	int max_cells = 0;
	std::vector<int> part_ncells( nParts+1, 0 );

	for( size_t icell=0; icell < gmsh->cells.size(); icell++ ){
		part_ncells[partitions[icell]]++;
	}

	for( int ipart=0; ipart < nParts+1; ipart++ ){
		if( max_cells < part_ncells[ipart] ){
			max_cells = part_ncells[ipart];
		}
	}
	new2old.resize( nParts+1, std::vector<int> ( max_cells, -1 ) );

	// Create one new mesh for each interior partition
	std::vector<gmsh_mesh> gmsh_parts( nParts+1 );
	std::vector<int>    local_ind ( nParts+1, 0 );
	std::vector<int> old2parts( gmsh->cells.size(), -1 );

	int part_ind;

	for( int ipart=0; ipart < nParts+1; ipart++ ){
		gmsh_parts[ipart].cells     .resize( part_ncells[ipart] );
		gmsh_parts[ipart].cell2neigh.resize( part_ncells[ipart] );
		for( int icell=0; icell < part_ncells[ipart]; icell++ ){
			gmsh_parts[ipart].cell2neigh[icell].resize( gmsh->faces_per_cell );
		}

		gmsh_parts[ipart].id             = ipart;
		gmsh_parts[ipart].faces_per_cell = gmsh->faces_per_cell;
		gmsh_parts[ipart].nodes_per_cell = gmsh->nodes_per_cell;
		gmsh_parts[ipart].cells_per_face = gmsh->cells_per_face;
		gmsh_parts[ipart].nodes_per_face = gmsh->nodes_per_face;
		gmsh_parts[ipart].elmnt_type     = gmsh->elmnt_type;
		gmsh_parts[ipart].mpi_neighbors  = gmsh->mpi_neighbors;

		gmsh_parts[ipart].cell2node.resize( gmsh_parts[ipart].cells.size(),
											std::vector<int> (gmsh_parts[ipart].nodes_per_cell, -1));
		gmsh_parts[ipart].cell2face.resize( gmsh_parts[ipart].cells.size(),
											std::vector<int> (gmsh_parts[ipart].faces_per_cell, -1));
	}

	// Initialize boundary cells separately since Hpath is already setup
	for( size_t icell=0; icell < gmsh->hpath.size(); icell++ ){
		gmsh_parts[INDEX_BND_SUBMESH].cells[icell].id = icell;

		new2old[INDEX_BND_SUBMESH][icell] = gmsh->hpath[icell];
		old2parts[gmsh->hpath[icell]] = icell;
	}
	for( size_t icell=0; icell < gmsh->deadend_cells.size(); icell++ ){
		int cell_id = icell + gmsh->hpath.size();

		gmsh_parts[INDEX_BND_SUBMESH].cells[cell_id].id = cell_id;

		new2old[INDEX_BND_SUBMESH][cell_id] = gmsh->deadend_cells[icell];
		old2parts[gmsh->deadend_cells[icell]] = cell_id;
	}

	std::vector<bool> is_boundary( gmsh->cells.size(), false );
	for( size_t icell=0; icell < gmsh->hpath.size(); icell++ ){
		is_boundary[gmsh->hpath[icell]] = true;
	}
	for( size_t icell=0; icell < gmsh->deadend_cells.size(); icell++ ){
		is_boundary[gmsh->deadend_cells[icell]] = true;
	}

	// Initialize id and new-old indices connectivity
	for( size_t icell=0; icell < gmsh->cells.size(); icell++ ){

		part_ind = partitions[icell];

		if( part_ind == INDEX_BND_SUBMESH ){
			continue;
		}

		gmsh_parts[part_ind].cells[local_ind[part_ind]].id = local_ind[part_ind];

		new2old[part_ind][local_ind[part_ind]] = gmsh->cells[icell].id ;
		old2parts [icell] = local_ind[part_ind];

		local_ind[part_ind]++;
	}
	int2old = new2old;

	// Update neighbor connectivity for interior mesh
	//		- TODO: should do this before
	for( size_t icell=0; icell < gmsh->cells.size(); icell++ ){

		// Skip boundary
		if( is_boundary[icell] ){
			continue;
		}

		for( int iface=0; iface < gmsh->faces_per_cell; iface++ ){
			t_gmsh_neighbor* neigh = &( gmsh->cell2neigh[icell][iface] );
			if( is_boundary[neigh->id] ){
				neigh->type = CELL_SUBMESH;
				neigh->sm   = INDEX_BND_SUBMESH;
			}
		}
	}

	// Initialize cells
	int gcell;
	for( int ipart=0; ipart < nParts+1; ipart++ ){
		for( size_t icell=0; icell < gmsh_parts[ipart].cells.size(); icell++ ){

			gcell = new2old[ipart][icell];

			gmsh_parts[ipart].cells[icell].xc[0] = gmsh->cells[gcell].xc[0]  ;
			gmsh_parts[ipart].cells[icell].xc[1] = gmsh->cells[gcell].xc[1]  ;

			for( int iface=0; iface < gmsh_parts[ipart].faces_per_cell; iface++ ){

				t_gmsh_neighbor *neigh_old = &( gmsh->cell2neigh[gcell][iface] ) ;
				t_gmsh_neighbor *neigh_new = &( gmsh_parts[ipart].cell2neigh[icell][iface] );

				switch( neigh_old->type ){
					case CELL_REGULAR:
					case CELL_SUBMESH:
						neigh_new->type = neigh_old->type;
						neigh_new->id   = old2parts[neigh_old->id] ;
						neigh_new->sm   = partitions[neigh_old->id];
						neigh_new->proc = myrank;
						break;
					case CELL_MPI:
					case CELL_PERIODIC:
						neigh_new->type = neigh_old->type;
						neigh_new->id   = neigh_old->id  ;
						neigh_new->sm   = INDEX_BND_SUBMESH;
						neigh_new->proc = neigh_old->proc;
						break;
					case CELL_PERIODIC_MPI:
						neigh_new->type = neigh_old->type;
						neigh_new->id   = neigh_old->id  ;
						neigh_new->sm   = INDEX_BND_SUBMESH;
						neigh_new->proc = neigh_old->proc;
						break;
					case CELL_NONE:
						neigh_new->type = neigh_old->type;
						neigh_new->id   = -1;
						neigh_new->sm   = -1;
						neigh_new->proc = -1;
						break;
					default:
						break;
				}
				neigh_new->bound = neigh_old->bound;

				if( neigh_old->type == CELL_MPI || neigh_old->type == CELL_PERIODIC_MPI || neigh_old->type == CELL_NONE ){
					continue;
				}

				if( partitions[neigh_old->id] != ipart ){
					neigh_new->type = CELL_SUBMESH;
				}
			}
		}
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Initialize nodes
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	int gnode;
	std::vector<std::vector<bool>> unregistered_nodes( nParts+1, std::vector<bool>(gmsh->nodes.size(),true));

	std::vector<int> n_nodes( nParts+1, 0 );
	for( int ipart=0; ipart < nParts+1; ipart++ ){

		for( size_t icell=0; icell < gmsh_parts[ipart].cells.size(); icell++ ){
			for( int inode=0; inode < gmsh->nodes_per_cell; inode++ ){
				gcell = new2old[ipart][icell];
				gnode = gmsh->cell2node[gcell][inode];

				if( unregistered_nodes[ipart][gnode] ){
					unregistered_nodes[ipart][gnode] = false;
					n_nodes[ipart]++;
				}
			}
		}
	}

	for( int ipart=0; ipart < nParts+1; ipart++ ){
		gmsh_parts[ipart].nodes.resize( n_nodes[ipart] );

		for( size_t inode=0; inode < gmsh->nodes.size(); inode++ ){
			unregistered_nodes[ipart][inode] = true;
		}
	}

	int node_ind;
	old2new_nodes.resize( nParts+1, std::vector<int>( gmsh->nodes.size(), -1 ) );
	for( int ipart=0; ipart < nParts+1; ipart++ ){

		node_ind = 0;
		for( size_t icell=0; icell < gmsh_parts[ipart].cells.size(); icell++ ){

			for( int inode=0; inode < gmsh->nodes_per_cell; inode++ ){
				gcell = new2old[ipart][icell];
				gnode = gmsh->cell2node[gcell][inode];

				if( unregistered_nodes[ipart][gnode] ){
					unregistered_nodes[ipart][gnode] = false;

					gmsh_parts[ipart].nodes[node_ind].id    = gmsh->nodes[gnode].id ;
					gmsh_parts[ipart].nodes[node_ind].xn[0] = gmsh->nodes[gnode].xn[0] ;
					gmsh_parts[ipart].nodes[node_ind].xn[1] = gmsh->nodes[gnode].xn[1] ;
					gmsh_parts[ipart].nodes[node_ind].xn[2] = gmsh->nodes[gnode].xn[2] ;

					old2new_nodes[ipart][gnode] = node_ind;

					node_ind++;
				}
			}
		}
	}
	for( int ipart=0; ipart < nParts+1; ipart++ ){
		for( size_t icell=0; icell < gmsh_parts[ipart].cells.size(); icell++ ){
			for( int inode=0; inode < gmsh->nodes_per_cell; inode++ ){
				gcell = new2old[ipart][icell];
				gnode = gmsh->cell2node[gcell][inode];

				gmsh_parts[ipart].cell2node[icell][inode] = old2new_nodes[ipart][gmsh->cell2node[gcell][inode]];
			}
		}
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Initialize faces
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	int gface;
	std::vector<std::vector<bool>> unregistered_faces( nParts+1, std::vector<bool>(gmsh->faces.size(),true));

	std::vector<int> n_faces( nParts+1, 0 );

	for( int ipart=0; ipart < nParts+1; ipart++ ){

		for( size_t icell=0; icell < gmsh_parts[ipart].cells.size(); icell++ ){
			for( int iface=0; iface < gmsh->faces_per_cell; iface++ ){
				gcell = new2old[ipart][icell];
				gface = gmsh->cell2face[gcell][iface];

				if( gface > (int) unregistered_faces[ipart].size() || gface < 0 ){
					cout << "Error with cell2face on rank " << myrank << ".";
					cout << "Terminating simulation prematurely." << endl;
					MPI_Abort( MPI_COMM_WORLD, 931 );
				}
				if( unregistered_faces[ipart][gface] ){
					unregistered_faces[ipart][gface] = false;
					n_faces[ipart]++;
				}
			}
		}
	}

	for( int ipart=0; ipart < nParts+1; ipart++ ){
		gmsh_parts[ipart].faces.resize( n_faces[ipart] );

		for( size_t iface=0; iface < gmsh->faces.size(); iface++ ){
			unregistered_faces[ipart][iface] = true;
		}
	}

	int face_ind;
	old2new_faces.resize( nParts+1, std::vector<int>( gmsh->faces.size(), -1 ) );

	MPI_Barrier( MPI_COMM_WORLD );

	for( int ipart=0; ipart < nParts+1; ipart++ ){

		face_ind = 0;
		for( size_t icell=0; icell < gmsh_parts[ipart].cells.size(); icell++ ){

			for( int iface=0; iface < gmsh->faces_per_cell; iface++ ){
				gcell = new2old[ipart][icell];
				gface = gmsh->cell2face[gcell][iface];

				if( unregistered_faces[ipart][gface] ){
					unregistered_faces[ipart][gface] = false;

					gmsh_parts[ipart].faces[face_ind].id = gmsh->faces[gface].id ;
					gmsh_parts[ipart].faces[face_ind].twin_face.type  = gmsh->faces[gface].twin_face.type ;
					gmsh_parts[ipart].faces[face_ind].twin_face.index = gmsh->faces[gface].twin_face.index;
					gmsh_parts[ipart].faces[face_ind].twin_face.rank  = gmsh->faces[gface].twin_face.rank ;

					old2new_faces[ipart][gface] = face_ind;

					face_ind++;
				}
			}
		}
	}

	for( int ipart=0; ipart < nParts+1; ipart++ ){
		for( size_t icell=0; icell < gmsh_parts[ipart].cells.size(); icell++ ){
			for( int iface=0; iface < gmsh->faces_per_cell; iface++ ){
				gcell = new2old[ipart][icell];
				gface = gmsh->cell2face[gcell][iface];

				gmsh_parts[ipart].cell2face[icell][iface] = old2new_faces[ipart][gmsh->cell2face[gcell][iface]];
			}
		}
	}

	// face2cell
	for( int ipart=0; ipart < nParts+1; ipart++ ){
		gmsh_parts[ipart].face2cell.resize( gmsh_parts[ipart].faces.size(),
									std::vector<int> ( gmsh_parts[ipart].cells_per_face, -1 ) );

		for( size_t icell=0; icell < gmsh_parts[ipart].cells.size(); icell++ ){
			for( int iface=0; iface < gmsh_parts[ipart].faces_per_cell; iface++ ){
				gface = gmsh_parts[ipart].cell2face[icell][iface];

				// If both cells are registered, job is already done
				if( gmsh_parts[ipart].face2cell[gface][0] != -1 &&
					gmsh_parts[ipart].face2cell[gface][1] != -1 ){
					continue;
				}

				//
				if( gmsh_parts[ipart].face2cell[gface][0] == -1 ){
					gmsh_parts[ipart].face2cell[gface][0] = icell;
				}else{
					gmsh_parts[ipart].face2cell[gface][1] = icell;
				}
			}
		}
	}

	// face2node
	int n0, n1, gface_new;

	for( int ipart=0; ipart < nParts+1; ipart++ ){
		gmsh_parts[ipart].face2node.resize( gmsh_parts[ipart].faces.size(),
									std::vector<int> ( gmsh_parts[ipart].nodes_per_face, -1 ) );

		for( size_t icell=0; icell < gmsh_parts[ipart].cells.size(); icell++ ){
			for( int iface=0; iface < gmsh->faces_per_cell; iface++ ){
				gcell = new2old[ipart][icell];

				n0 = gmsh->cell2node[gcell][ iface ];
				n1 = gmsh->cell2node[gcell][ (iface+1) % gmsh->faces_per_cell ];

				gface_new = gmsh_parts[ipart].cell2face[icell][iface];

				gmsh_parts[ipart].face2node[gface_new][0] = old2new_nodes[ipart][n0];
				gmsh_parts[ipart].face2node[gface_new][1] = old2new_nodes[ipart][n1];
			}
		}
	}

	// Compute mesh center
	for( int ipart=0; ipart < nParts+1; ipart++ ){
		gmsh_parts[ipart].x[0] = 0.;
		gmsh_parts[ipart].x[1] = 0.;

		for( size_t icell=0; icell < gmsh_parts[ipart].cells.size(); icell++ ){
			gmsh_parts[ipart].x[0] += gmsh_parts[ipart].cells[icell].xc[0] ;
			gmsh_parts[ipart].x[1] += gmsh_parts[ipart].cells[icell].xc[1] ;
		}

		gmsh_parts[ipart].x[0] = gmsh_parts[ipart].x[0] / (double) gmsh_parts[ipart].cells.size() ;
		gmsh_parts[ipart].x[1] = gmsh_parts[ipart].x[1] / (double) gmsh_parts[ipart].cells.size() ;
	}

	// Mark boundary cells
	for( int ipart=0; ipart < nParts+1; ipart++ ){
		for( size_t icell=0; icell < gmsh_parts[ipart].cells.size(); icell++ ){
			for( int iface=0; iface < gmsh_parts[ipart].faces_per_cell; iface++ ){
				if( gmsh_parts[ipart].cell2neigh[icell][iface].type == CELL_SUBMESH ){
					gmsh_parts[ipart].cells[icell].is_bcell = true;
				}
			}
		}
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Periodic boundary conditions: need to reorder MPI twin faces' connectivity
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	int nProcs;
	MPI_Comm_size( MPI_COMM_WORLD, &nProcs );

	gmsh_mesh *mesh_bnd = &( gmsh_parts[INDEX_BND_SUBMESH] );

	// First gather the number of faces of each MPI rank
	std::vector<int> faces_per_rank( nProcs, -1 );
	int num_faces  = (int) gmsh->faces.size();
	MPI_Allgather( &num_faces, 1, MPI_INT, &faces_per_rank[0], 1, MPI_INT, MPI_COMM_WORLD );

	// Need temporary 1D array since MPI does not deal well with jagged arrays
	std::vector<int> tmp_faces;
	int total_faces = 0;
	for( int irank=0; irank < nProcs; irank++ ){
		total_faces += faces_per_rank[irank];
	}
	tmp_faces.resize( total_faces, -1 );

	std::vector<int> disp( nProcs, 0 );
	for( int irank=1; irank < nProcs; irank++ ){
		disp[irank] = faces_per_rank[irank-1] + disp[irank-1];
	}

	// Gather the old2new mapping of each MPI rank
	MPI_Barrier( MPI_COMM_WORLD );
	std::vector<int> faces2hpath( old2new_faces[INDEX_BND_SUBMESH].size(), -1 );
	for( size_t iface=0; iface < old2new_faces[INDEX_BND_SUBMESH].size(); iface++ ){
		faces2hpath[iface] = old2new_faces[INDEX_BND_SUBMESH][iface];
	}

	MPI_Allgatherv( &faces2hpath[0], faces2hpath.size(), MPI_INT,
				    &tmp_faces[0],   &faces_per_rank[0], &disp[0], MPI_INT, MPI_COMM_WORLD );

	std::vector<std::vector<int>> global_faces2hpath( nProcs );
	for( int irank=0; irank < nProcs; irank++ ){
		global_faces2hpath[irank].resize( faces_per_rank[irank], -1 );
	}

	int ind = 0;
	for( int irank=0; irank < nProcs; irank++ ){
		for( size_t iface=0; iface < global_faces2hpath[irank].size(); iface++ ){
			global_faces2hpath[irank][iface] = tmp_faces[ind++];
		}
	}

	// Correct periodic twin faces accross MPI rank
	for( size_t icell=0; icell < mesh_bnd->cells.size(); icell++ ){
		for( int iface=0; iface < mesh_bnd->faces_per_cell; iface++ ){

			if( mesh_bnd->cell2neigh[icell][iface].type != CELL_PERIODIC &&
				mesh_bnd->cell2neigh[icell][iface].type != CELL_PERIODIC_MPI ){
				continue;
			}

			int gface       = mesh_bnd->cell2face[icell][iface];
			int neigh_rank  = mesh_bnd->cell2neigh[icell][iface].proc;
			int neigh_index = mesh_bnd->faces[gface].twin_face.index ;

			mesh_bnd->faces[gface].twin_face.index = global_faces2hpath[neigh_rank][neigh_index];
		}
	}

	return gmsh_parts;
}

//***************************************************************************************************
// Initialize metis' mesh patitions
// 		- only function to interact directly with libmetis
//***************************************************************************************************
gmsh_mesh partition_mesh_mpi_efficient( MPI_env mpi_env, gmsh_mesh *gmsh_domain ){

	gmsh_mesh gmsh_local;

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Perform first level of mesh partitioning: one partition per MPI rank
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	std::vector<int>       partitions_mesh;
	std::vector<gmsh_mesh> gmsh_node;

	init_mesh_size_mpi_type( mpi_env );
	if( mpi_env.is_node_master() )
		metis_partition_mesh_new( gmsh_domain, mpi_env.size(), partitions_mesh );

	gmsh_local = node_local_mesh_partitioning( mpi_env, gmsh_domain, partitions_mesh );

	// Temporary fix
	for( size_t icell=0; icell < gmsh_local.cells.size(); icell++ ){
		for( int iface=0; iface < gmsh_local.faces_per_cell; iface++ ){
			if( gmsh_local.cell2neigh[icell][iface].type != CELL_REGULAR ){
				gmsh_local.cells[icell].is_bcell = true;
			}
		}
	}

	// Check for bugs
	for( size_t icell=0; icell < gmsh_local.cells.size(); icell++ ){
		for( int iface=0; iface < gmsh_local.faces_per_cell; iface++ ){
			int gface = gmsh_local.cell2face[ icell ][ iface ];

			// face nodes
			int fn0 = gmsh_local.face2node[ gface ][ 0 ];
			int fn1 = gmsh_local.face2node[ gface ][ 1 ];

			// cell nodes
			int cn0 = gmsh_local.cell2node[ icell ][ iface ];
			int cn1 = gmsh_local.cell2node[ icell ][ (iface+1) % gmsh_local.faces_per_cell ];

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
				MPI_Abort( MPI_COMM_WORLD, 1906 );
			}
		}
	}

	return gmsh_local;
}

//***************************************************************************************************
//***************************************************************************************************
gmsh_mesh node_local_mesh_partitioning( MPI_env          mpi_env,
										gmsh_mesh       *gmsh_domain,
										std::vector<int> partitions ){

	//if( mpi_env.size() == 1 ){
	//	gmsh_mesh gmsh_local = *gmsh_domain;
	//	return gmsh_local;
	//}

	const int thousand = 1000;

	// Create a mesh for each node-local CPU
	struct t_mesh_size mesh_size;
	std::vector<int> old2new_faces;

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Share faces normal
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	int total_faces;
	if( mpi_env.is_node_master() ){
		total_faces = gmsh_domain->faces.size();
	}
	MPI_Bcast( &total_faces, 1, MPI_INT, mpi_env.node_comm.master(), mpi_env.node_comm.self() );

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Exchange mesh size information
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	if( mpi_env.is_node_master() ){
		for( int irank=0; irank < mpi_env.node_comm.size(); irank++ ){

			// Exchange mesh size
			int send_rank_world = mpi_env.node_local_ranks[irank];
			struct t_mesh_size tmp_mesh_size = get_mesh_size( gmsh_domain, send_rank_world, partitions );

			if( mpi_env.rank() == send_rank_world ){
				mesh_size = tmp_mesh_size;
				continue;
			}

			int send_rank_node = (mpi_env.node_local_ranks[irank] - mpi_env.rank());
			int send_tag       = thousand*mpi_env.node_comm.master() + send_rank_node;	// Works for up to 999 cores per rank

			MPI_Send( &tmp_mesh_size, 1, mpi_env.get_type_meshsize(), send_rank_node, send_tag, mpi_env.node_comm.self() );
		}
	}else{
		// Receive mesh size
		int recv_tag = mpi_env.node_comm.rank() + thousand*mpi_env.node_comm.master() ;
		MPI_Recv( &mesh_size, 1, mpi_env.get_type_meshsize(), mpi_env.node_comm.master(), recv_tag,
			  	  mpi_env.node_comm.self(), MPI_STATUS_IGNORE );
	}

	// Allocate mesh
	gmsh_mesh gmsh_local;
	gmsh_local = allocate_mesh( mesh_size );

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Send elements
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	MPI_Datatype type_node  = get_mpitype_meshnode ( gmsh_domain );
	MPI_Datatype type_face  = get_mpitype_meshface ( gmsh_domain );
	MPI_Datatype type_cell  = get_mpitype_meshcell ( gmsh_domain );
	MPI_Datatype type_neigh = get_mpitype_meshneigh( gmsh_domain );

	// Send nodes
	if( mpi_env.is_node_master() ){
		for( int irank=0; irank < mpi_env.node_comm.size(); irank++ ){

			// Exchange mesh size
			int send_rank_world = mpi_env.node_local_ranks[irank];
			int send_rank_node  = (mpi_env.node_local_ranks[irank] - mpi_env.rank());	// TODO: not 100% it works for all cases
			int send_tag        = thousand*mpi_env.node_comm.master() + send_rank_node;	// Works for up to 999 cores per rank
			int send_size;

			std::vector<int> tmp_old2new_faces;
			gmsh_mesh mesh_rank = get_mesh_partition( mpi_env, gmsh_domain, send_rank_world, partitions, tmp_old2new_faces );

			struct t_mesh_size mesh_size = get_mesh_size( &mesh_rank );

			if( mpi_env.rank() == send_rank_world ){
				gmsh_local = mesh_rank;
				old2new_faces = tmp_old2new_faces;

				if( mpi_env.is_serial() ){
					return gmsh_local;
				}

				continue;
			}

			// Send nodes
			send_tag = thousand*mpi_env.node_comm.master() + send_rank_node;
			MPI_Send( &mesh_rank.nodes[0], mesh_size.n_nodes, type_node, send_rank_node, send_tag,
					  mpi_env.node_comm.self() );

			// Send faces
			send_tag = thousand*mpi_env.node_comm.master() + send_rank_node + 1;
			MPI_Send( &mesh_rank.faces[0], mesh_size.n_faces, type_face, send_rank_node, send_tag,
					  mpi_env.node_comm.self() );

			// Send cells
			send_tag = thousand*mpi_env.node_comm.master() + send_rank_node + 2;
			MPI_Send( &mesh_rank.cells[0], mesh_size.n_cells, type_cell, send_rank_node, send_tag,
					  mpi_env.node_comm.self() );

			// Neighbors
			t_gmsh_neighbor *tmp_neighs = from_2Dvec_to_1Darray( mesh_size.n_cells,
																 mesh_size.faces_per_cell,
																 mesh_rank.cell2neigh );
			send_tag  = thousand*mpi_env.node_comm.master() + send_rank_node + 3;
			send_size = mesh_size.n_cells * mesh_size.faces_per_cell;
			MPI_Send( &tmp_neighs[0], send_size, type_neigh, send_rank_node, send_tag,
					  mpi_env.node_comm.self() );
			delete[] tmp_neighs;

			// cell2node
			send_size = mesh_size.n_cells * mesh_size.nodes_per_cell ;
			int *tmp_cell2node = from_2Dvec_to_1Darray( mesh_size.n_cells,
														mesh_size.nodes_per_cell,
														mesh_rank.cell2node );
			send_tag  = thousand*mpi_env.node_comm.master() + send_rank_node + 4;
			MPI_Send( &tmp_cell2node[0], send_size, MPI_INT, send_rank_node, send_tag,
					  mpi_env.node_comm.self() );
			delete[] tmp_cell2node;

			// cell2face
			send_size = mesh_size.n_cells * mesh_size.faces_per_cell ;
			int *tmp_cell2face = from_2Dvec_to_1Darray( mesh_size.n_cells,
														mesh_size.faces_per_cell,
														mesh_rank.cell2face );
			send_tag  = thousand*mpi_env.node_comm.master() + send_rank_node + 5;
			MPI_Send( &tmp_cell2face[0], send_size, MPI_INT, send_rank_node, send_tag,
					  mpi_env.node_comm.self() );
			delete[] tmp_cell2face;

			// face2node
			send_size = mesh_size.n_faces * mesh_size.nodes_per_face ;
			int *tmp_face2node = from_2Dvec_to_1Darray( mesh_size.n_faces,
														mesh_size.nodes_per_face,
														mesh_rank.face2node );
			send_tag  = thousand*mpi_env.node_comm.master() + send_rank_node + 6;
			MPI_Send( &tmp_face2node[0], send_size, MPI_INT, send_rank_node, send_tag,
					  mpi_env.node_comm.self() );
			delete[] tmp_face2node;

			// face2cell
			send_size = mesh_size.n_faces * mesh_size.cells_per_face ;
			int *tmp_face2cell = from_2Dvec_to_1Darray( mesh_size.n_faces,
														mesh_size.cells_per_face,
														mesh_rank.face2cell );
			send_tag  = thousand*mpi_env.node_comm.master() + send_rank_node + 7;
			MPI_Send( &tmp_face2cell[0], send_size, MPI_INT, send_rank_node, send_tag,
					  mpi_env.node_comm.self() );
			delete[] tmp_face2cell;

			// Old2new faces
			send_size = (int) tmp_old2new_faces.size();
			send_tag  = thousand*mpi_env.node_comm.master() + send_rank_node + 8;
			MPI_Send( &send_size, 1, MPI_INT, send_rank_node, send_tag, mpi_env.node_comm.self() );

			send_tag  = thousand*mpi_env.node_comm.master() + send_rank_node + 9;
			MPI_Send( &tmp_old2new_faces[0], send_size, MPI_INT, send_rank_node, send_tag,
					  mpi_env.node_comm.self() );

		}
	}else{
		int recv_tag, recv_size;

		// Receive nodes
		recv_tag = mpi_env.node_comm.rank() + thousand*mpi_env.node_comm.master();
		MPI_Recv( &gmsh_local.nodes[0], mesh_size.n_nodes, type_node, mpi_env.node_comm.master(),
			  	  recv_tag, mpi_env.node_comm.self(), MPI_STATUS_IGNORE );

		// Receive faces
		recv_tag = mpi_env.node_comm.rank() + thousand*mpi_env.node_comm.master() + 1;
		MPI_Recv( &gmsh_local.faces[0], mesh_size.n_faces, type_face, mpi_env.node_comm.master(),
			  	  recv_tag, mpi_env.node_comm.self(), MPI_STATUS_IGNORE );

		// Receive cells
		recv_tag = mpi_env.node_comm.rank() + thousand*mpi_env.node_comm.master() + 2;
		MPI_Recv( &gmsh_local.cells[0], mesh_size.n_cells, type_cell, mpi_env.node_comm.master(),
			  	  recv_tag, mpi_env.node_comm.self(), MPI_STATUS_IGNORE );

		// Connectivity vectors
		t_gmsh_neighbor *tmp_neighs = new t_gmsh_neighbor[mesh_size.n_cells*mesh_size.faces_per_cell];
		recv_tag  = mpi_env.node_comm.rank() + thousand*mpi_env.node_comm.master() + 3;
		recv_size = mesh_size.n_cells * mesh_size.faces_per_cell;
		MPI_Recv( &tmp_neighs[0], recv_size, type_neigh, mpi_env.node_comm.master(),
			  	  recv_tag, mpi_env.node_comm.self(), MPI_STATUS_IGNORE );

		gmsh_local.cell2neigh = from_1Darray_to_2Dvec( mesh_size.n_cells,
													   mesh_size.faces_per_cell,
													   tmp_neighs );
		delete[] tmp_neighs;

		// cell2node
		recv_size = mesh_size.n_cells * mesh_size.nodes_per_cell;
		int *tmp_cell2node = new int [ recv_size ];
		recv_tag = mpi_env.node_comm.rank() + thousand*mpi_env.node_comm.master() + 4;
		MPI_Recv( &tmp_cell2node[0], recv_size, MPI_INT, mpi_env.node_comm.master(),
			  	  recv_tag, mpi_env.node_comm.self(), MPI_STATUS_IGNORE );

		gmsh_local.cell2node = from_1Darray_to_2Dvec( mesh_size.n_cells,
													  mesh_size.nodes_per_cell,
													  tmp_cell2node );
		delete[] tmp_cell2node;

		// cell2face
		recv_size = mesh_size.n_cells * mesh_size.faces_per_cell;
		int *tmp_cell2face = new int [ recv_size ];
		recv_tag = mpi_env.node_comm.rank() + thousand*mpi_env.node_comm.master() + 5;
		MPI_Recv( &tmp_cell2face[0], recv_size, MPI_INT, mpi_env.node_comm.master(),
			  	  recv_tag, mpi_env.node_comm.self(), MPI_STATUS_IGNORE );

		gmsh_local.cell2face = from_1Darray_to_2Dvec( mesh_size.n_cells,
													  mesh_size.faces_per_cell,
													  tmp_cell2face );
		delete[] tmp_cell2face;

		// face2node
		recv_size = mesh_size.n_faces * mesh_size.nodes_per_face;
		int *tmp_face2node = new int [ recv_size ];
		recv_tag = mpi_env.node_comm.rank() + thousand*mpi_env.node_comm.master() + 6;
		MPI_Recv( &tmp_face2node[0], recv_size, MPI_INT, mpi_env.node_comm.master(),
			  	  recv_tag, mpi_env.node_comm.self(), MPI_STATUS_IGNORE );

		gmsh_local.face2node = from_1Darray_to_2Dvec( mesh_size.n_faces,
													  mesh_size.nodes_per_face,
													  tmp_face2node );
		delete[] tmp_face2node;

		// face2cell
		recv_size = mesh_size.n_faces * mesh_size.faces_per_cell;
		int *tmp_face2cell = new int [ recv_size ];
		recv_tag = mpi_env.node_comm.rank() + thousand*mpi_env.node_comm.master() + 7;
		MPI_Recv( &tmp_face2cell[0], recv_size, MPI_INT, mpi_env.node_comm.master(),
			  	  recv_tag, mpi_env.node_comm.self(), MPI_STATUS_IGNORE );

		gmsh_local.face2cell = from_1Darray_to_2Dvec( mesh_size.n_faces,
													  mesh_size.cells_per_face,
													  tmp_face2cell );
		delete[] tmp_face2cell;

		// Old2new faces
		int total_faces = -1;
		recv_tag = mpi_env.node_comm.rank() + thousand*mpi_env.node_comm.master() + 8;
		MPI_Recv( &total_faces, 1, MPI_INT, mpi_env.node_comm.master(), recv_tag,
				  mpi_env.node_comm.self(), MPI_STATUS_IGNORE );

		old2new_faces.resize( total_faces, -1 );
		recv_size = total_faces;
		recv_tag = mpi_env.node_comm.rank() + thousand*mpi_env.node_comm.master() + 9;

		MPI_Recv( &old2new_faces[0], recv_size, MPI_INT, mpi_env.node_comm.master(),
			  	  recv_tag, mpi_env.node_comm.self(), MPI_STATUS_IGNORE );
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// MPI connectivity
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	// Neighbors
	update_mpi_connectivity( mpi_env, gmsh_domain, &gmsh_local );

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Periodic faces
	//     - for each MPI neighbor, find the shared faces and update the numbering
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	update_mpi_connectivity_faces( mpi_env, gmsh_domain, &gmsh_local );

	//std::vector<gmsh_mesh> mesh_mpi( 1 );
	return gmsh_local;
}

//***************************************************************************************************
gmsh_mesh get_mesh_partition( MPI_env mpi_env, gmsh_mesh *mesh_domain, int part_num,
							  std::vector<int> partitions, std::vector<int> &old2new_faces ){

	// Compute the number of cells for each mesh
	int local_ncells = 0;
	for( size_t icell=0; icell < mesh_domain->cells.size(); icell++ ){
		if( partitions[icell] == part_num ){
			local_ncells++;
		}
	}
	std::vector<int> part2old( local_ncells );

	// Create one new mesh for each interior partition
	gmsh_mesh mesh_local;

	std::vector<int> old2part( mesh_domain->cells.size(), -1 );

	mesh_local.cells.resize( local_ncells );

	mesh_local.cell2neigh.resize( local_ncells );
	for( int icell=0; icell < local_ncells; icell++ ){
		mesh_local.cell2neigh[icell].resize( mesh_domain->faces_per_cell );
	}
	mesh_local.id             = part_num;
	mesh_local.faces_per_cell = mesh_domain->faces_per_cell;
	mesh_local.nodes_per_cell = mesh_domain->nodes_per_cell;
	mesh_local.cells_per_face = mesh_domain->cells_per_face;
	mesh_local.nodes_per_face = mesh_domain->nodes_per_face;
	mesh_local.elmnt_type     = mesh_domain->elmnt_type;

	mesh_local.cell2node.resize( local_ncells, std::vector<int> (mesh_local.nodes_per_cell, -1));
	mesh_local.cell2face.resize( local_ncells, std::vector<int> (mesh_local.faces_per_cell, -1));

	// Initialize id and new-old indices connectivity
	int local_ind = 0;
	for( size_t icell=0; icell < mesh_domain->cells.size(); icell++ ){
		if( partitions[icell] == part_num ){

			mesh_local.cells[local_ind].id        = local_ind;
			mesh_local.cells[local_ind].id_global = icell;

			part2old[local_ind] = mesh_domain->cells[icell].id ;

			old2part[icell] = local_ind++;
		}
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Initialize cells
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	// Allocate neighbors
	for( size_t icell=0; icell < mesh_local.cells.size(); icell++ ){
		mesh_local.cell2neigh[icell].resize( mesh_local.faces_per_cell );
		for( int iface=0; iface < mesh_local.faces_per_cell; iface++ ){
			mesh_local.cell2neigh[icell][iface].type  = CELL_NONE;
			mesh_local.cell2neigh[icell][iface].id    = -1;
			mesh_local.cell2neigh[icell][iface].sm    = -1;
			mesh_local.cell2neigh[icell][iface].proc  = -1;
			mesh_local.cell2neigh[icell][iface].bound =  0;
		}
	}

	// Initialize neighbors: only local for now
	for( size_t icell=0; icell < mesh_local.cells.size(); icell++ ){

		int gcell = part2old[icell];

		mesh_local.cells[icell].xc[0] = mesh_domain->cells[gcell].xc[0] ;
		mesh_local.cells[icell].xc[1] = mesh_domain->cells[gcell].xc[1] ;

		for( int iface=0; iface < mesh_local.faces_per_cell; iface++ ){

			t_gmsh_neighbor *neigh = &(mesh_domain->cell2neigh[gcell][iface]) ;

			// Cell is out of domain boundary
			if( neigh->type == CELL_NONE ){
				mesh_local.cell2neigh[icell][iface].type  = CELL_NONE;
				mesh_local.cell2neigh[icell][iface].id    = -1;
				mesh_local.cell2neigh[icell][iface].sm    = -1;
				mesh_local.cell2neigh[icell][iface].proc  = -1;
				mesh_local.cell2neigh[icell][iface].bound = neigh->bound;
			// Cell belongs to the current MPI mesh
			}else if( partitions[neigh->id] == part_num && neigh->type != CELL_PERIODIC ){
				mesh_local.cell2neigh[icell][iface].type = CELL_REGULAR;
				mesh_local.cell2neigh[icell][iface].id   = old2part[neigh->id];
				mesh_local.cell2neigh[icell][iface].sm   = -1;
				mesh_local.cell2neigh[icell][iface].proc = part_num;
			// Cell is on a neighboring processor
			}else{
				mesh_local.cell2neigh[icell][iface].type = CELL_MPI;
				mesh_local.cell2neigh[icell][iface].id   = neigh->id;
				mesh_local.cell2neigh[icell][iface].sm   = INDEX_BND_SUBMESH;
				mesh_local.cell2neigh[icell][iface].proc = partitions[neigh->id];
			}

			// Special case: periodic neighbors
			if( neigh->type == CELL_PERIODIC ){
				if( partitions[neigh->id] == part_num ){
					mesh_local.cell2neigh[icell][iface].type = CELL_PERIODIC;
					mesh_local.cell2neigh[icell][iface].id   = old2part[neigh->id];
					mesh_local.cell2neigh[icell][iface].sm   = INDEX_BND_SUBMESH;
					mesh_local.cell2neigh[icell][iface].proc = part_num;
				}else{
					mesh_local.cell2neigh[icell][iface].type = CELL_PERIODIC_MPI;
					mesh_local.cell2neigh[icell][iface].id   = neigh->id;
					mesh_local.cell2neigh[icell][iface].sm   = INDEX_BND_SUBMESH;
					mesh_local.cell2neigh[icell][iface].proc = partitions[neigh->id];
				}
			}
		}
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Initialize nodes
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	std::vector<bool> unregistered_nodes( mesh_domain->nodes.size(), true );

	int n_nodes = 0;

	for( size_t icell=0; icell < mesh_local.cells.size(); icell++ ){
		for( int inode=0; inode < mesh_domain->nodes_per_cell; inode++ ){
			int gcell = part2old[icell];
			int gnode = mesh_domain->cell2node[gcell][inode];

			if( unregistered_nodes[gnode] ){
				unregistered_nodes[gnode] = false;
				n_nodes++;
			}
		}
	}

	mesh_local.nodes.resize( n_nodes );

	for( size_t inode=0; inode < mesh_domain->nodes.size(); inode++ ){
		unregistered_nodes[inode] = true;
	}

	int node_ind = 0;
	std::vector<int> old2new_nodes( mesh_domain->nodes.size(), -1 );

	for( size_t icell=0; icell < mesh_local.cells.size(); icell++ ){

		for( int inode=0; inode < mesh_domain->nodes_per_cell; inode++ ){
			int gcell = part2old[icell];
			int gnode = mesh_domain->cell2node[gcell][inode];

			if( unregistered_nodes[gnode] ){
				unregistered_nodes[gnode] = false;

				mesh_local.nodes[node_ind].id = mesh_domain->nodes[gnode].id ;
				mesh_local.nodes[node_ind].xn[0] = mesh_domain->nodes[gnode].xn[0]  ;
				mesh_local.nodes[node_ind].xn[1] = mesh_domain->nodes[gnode].xn[1]  ;
				mesh_local.nodes[node_ind].xn[2] = mesh_domain->nodes[gnode].xn[2]  ;

				mesh_local.nodes[node_ind].id_global = mesh_domain->nodes[gnode].id;

				old2new_nodes[gnode] = node_ind;

				node_ind++;
			}
		}
	}

	for( size_t icell=0; icell < mesh_local.cells.size(); icell++ ){
		for( int inode=0; inode < mesh_domain->nodes_per_cell; inode++ ){
			int gcell = part2old[icell];

			mesh_local.cell2node[icell][inode] = old2new_nodes[mesh_domain->cell2node[gcell][inode]];
		}
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Initialize faces
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	std::vector<bool> unregistered_faces( mesh_domain->faces.size(), true );

	int n_faces = 0;

	for( size_t icell=0; icell < mesh_local.cells.size(); icell++ ){
		for( int iface=0; iface < mesh_domain->faces_per_cell; iface++ ){
			int gcell = part2old[icell];
			int gface = mesh_domain->cell2face[gcell][iface];

			if( unregistered_faces[gface] ){
				unregistered_faces[gface] = false;
				n_faces++;
			}
		}
	}

	mesh_local.faces.resize( n_faces );

	for( size_t iface=0; iface < mesh_domain->faces.size(); iface++ ){
		unregistered_faces[iface] = true;
	}

	int face_ind = 0;
	old2new_faces.resize( mesh_domain->faces.size(), -1 );
	std::vector<int> part2old_faces( mesh_local.faces.size(), -1 );

	for( size_t icell=0; icell < mesh_local.cells.size(); icell++ ){

		for( int iface=0; iface < mesh_domain->faces_per_cell; iface++ ){
			int gcell = part2old[icell];
			int gface = mesh_domain->cell2face[gcell][iface];

			if( unregistered_faces[gface] ){
				unregistered_faces[gface] = false;

				mesh_local.faces[face_ind].id        = face_ind;
				mesh_local.faces[face_ind].id_global = gface;
				old2new_faces[gface] = face_ind;
				part2old_faces[face_ind] = gface;

				face_ind++;
			}
		}
	}

	// cell2face
	for( size_t icell=0; icell < mesh_local.cells.size(); icell++ ){
		for( int iface=0; iface < mesh_domain->faces_per_cell; iface++ ){
			int gcell = part2old[icell];

			mesh_local.cell2face[icell][iface] = old2new_faces[mesh_domain->cell2face[gcell][iface]];
		}
	}

	// face2cell
	mesh_local.face2cell.resize( mesh_local.faces.size(),
								std::vector<int> ( mesh_local.cells_per_face, -1 ) );

	for( size_t icell=0; icell < mesh_local.cells.size(); icell++ ){
		for( int iface=0; iface < mesh_local.faces_per_cell; iface++ ){
			int gface = mesh_local.cell2face[icell][iface];

			// If both cells are registered, job is already done
			if( mesh_local.face2cell[gface][0] != -1 &&
				mesh_local.face2cell[gface][1] != -1 ){
				continue;
			}

			//
			if( mesh_local.face2cell[gface][0] == -1 ){
				mesh_local.face2cell[gface][0] = icell;
			}else{
				mesh_local.face2cell[gface][1] = icell;
			}
		}
	}

	// face2node
	mesh_local.face2node.resize( mesh_local.faces.size(), std::vector<int> ( mesh_local.nodes_per_face, -1 ) );

	for( size_t icell=0; icell < mesh_local.cells.size(); icell++ ){
		for( int iface=0; iface < mesh_domain->faces_per_cell; iface++ ){
			int gcell = part2old[icell];

			int n0 = mesh_domain->cell2node[gcell][ iface ];
			int n1 = mesh_domain->cell2node[gcell][ (iface+1) % mesh_domain->faces_per_cell ];

			int gface_new = mesh_local.cell2face[icell][iface];

			mesh_local.face2node[gface_new][0] = old2new_nodes[n0];
			mesh_local.face2node[gface_new][1] = old2new_nodes[n1];
		}
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Miscellaneous
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	// Compute mesh center
	mesh_local.x[0] = 0.;
	mesh_local.x[1] = 0.;

	for( size_t icell=0; icell < mesh_local.cells.size(); icell++ ){
		mesh_local.x[0] += mesh_local.cells[icell].xc[0] ;
		mesh_local.x[1] += mesh_local.cells[icell].xc[1] ;
	}

	mesh_local.x[0] = mesh_local.x[0] / (double) mesh_local.cells.size() ;
	mesh_local.x[1] = mesh_local.x[1] / (double) mesh_local.cells.size() ;

	// Submesh connectivity
	for( size_t icell=0; icell < mesh_local.cells.size(); icell++ ){
		for( int iface=0; iface < mesh_local.faces_per_cell; iface++ ){
			if( mesh_local.cell2neigh[icell][iface].type != CELL_REGULAR ){
				if( mpi_env.rank() == 1 )
					cout << "New bcell =" << setw(10) << icell << endl;
				mesh_local.cells[icell].is_bcell = true;
				break;
			}
		}
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// twin_face: only works for serial case
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	if( mpi_env.size() == 1 ){
		for( size_t iface=0; iface < mesh_local.faces.size(); iface++ ){
			mesh_local.faces[iface].twin_face.type    = PERIODIC_NEIGH_NONE;
			mesh_local.faces[iface].twin_face.rank    = -1;
			mesh_local.faces[iface].twin_face.index   = -1;
			mesh_local.faces[iface].twin_face.orig_id = -1;
		}

		for( size_t iface=0; iface < mesh_local.faces.size(); iface++ ){
			int gface = part2old_faces[iface];

			if( mesh_domain->faces[gface].twin_face.index != -1 ){
				mesh_local.faces[iface].twin_face.index = old2new_faces[mesh_domain->faces[gface].twin_face.index];
			}
			mesh_local.faces[iface].twin_face.rank  = 0;
		}
		for( size_t icell=0; icell < mesh_local.cells.size(); icell++ ){
			for( int iface=0; iface < mesh_domain->faces_per_cell; iface++ ){
				int gface = mesh_local.cell2face[icell][iface];

				if( mesh_local.cell2neigh[icell][iface].type == CELL_PERIODIC )
					mesh_local.faces[gface].twin_face.type = PERIODIC_NEIGH_REGULAR;
				if( mesh_local.cell2neigh[icell][iface].type == CELL_PERIODIC_MPI )
					mesh_local.faces[gface].twin_face.type = PERIODIC_NEIGH_MPI;
			}
		}
	}else{
		for( size_t iface=0; iface < mesh_local.faces.size(); iface++ ){
			mesh_local.faces[iface].twin_face.type    = PERIODIC_NEIGH_NONE;
			mesh_local.faces[iface].twin_face.rank    = -1;
			mesh_local.faces[iface].twin_face.index   = -1;
			mesh_local.faces[iface].twin_face.orig_id = -1;
		}

		for( size_t iface=0; iface < mesh_local.faces.size(); iface++ ){
			int gface = part2old_faces[iface];
			mesh_local.faces[iface].twin_face.index = mesh_domain->faces[gface].twin_face.index;
		}
		for( size_t icell=0; icell < mesh_local.cells.size(); icell++ ){
			for( int iface=0; iface < mesh_domain->faces_per_cell; iface++ ){
				int gface = mesh_local.cell2face[icell][iface];

				if( mesh_local.cell2neigh[icell][iface].type == CELL_PERIODIC )
					mesh_local.faces[gface].twin_face.type = PERIODIC_NEIGH_REGULAR;
				if( mesh_local.cell2neigh[icell][iface].type == CELL_PERIODIC_MPI )
					mesh_local.faces[gface].twin_face.type = PERIODIC_NEIGH_MPI;
			}
		}
	}

	return mesh_local;
}

//***************************************************************************************************
// Based on a given mesh object, return the mesh sizes (used to send mesh objects with MPI)
//***************************************************************************************************
struct t_mesh_size get_mesh_size( gmsh_mesh *mesh ){

	t_mesh_size mesh_size;

	mesh_size.n_nodes = (int) mesh->nodes.size();
	mesh_size.n_faces = (int) mesh->faces.size();
	mesh_size.n_cells = (int) mesh->cells.size();

	mesh_size.nodes_per_face = mesh->nodes_per_face;
	mesh_size.cells_per_face = mesh->cells_per_face;
	mesh_size.nodes_per_cell = mesh->nodes_per_cell;
	mesh_size.faces_per_cell = mesh->faces_per_cell;

	return mesh_size;
}

//***************************************************************************************************
// Based on the domain mesh, return the partition's mesh_size
//***************************************************************************************************
struct t_mesh_size get_mesh_size( gmsh_mesh *mesh_domain, int part_num, std::vector<int> partitions ){

	// Count cells
	int n_cells = 0;
	for( size_t icell=0; icell < mesh_domain->cells.size(); icell++ ){
		if( partitions[icell] == part_num ){
			n_cells++;
		}
	}

	// Partition to domain mesh index mapping
	std::vector<int> part2old( n_cells, -1 );

	int local_ind = 0;
	for( size_t icell=0; icell < mesh_domain->cells.size(); icell++ ){
		if( partitions[icell] == part_num ){
			part2old[local_ind++] = mesh_domain->cells[icell].id ;
		}
	}

	// Count nodes
	std::vector<bool> unregistered_nodes( mesh_domain->nodes.size(), true );

	int n_nodes = 0;

	for( int icell=0; icell < n_cells; icell++ ){
		for( int inode=0; inode < mesh_domain->nodes_per_cell; inode++ ){
			int gcell = part2old[icell];
			int gnode = mesh_domain->cell2node[gcell][inode];

			if( unregistered_nodes[gnode] ){
				unregistered_nodes[gnode] = false;
				n_nodes++;
			}
		}
	}

	// Count faces
	std::vector<bool> unregistered_faces( mesh_domain->faces.size(), true );

	int n_faces = 0;

	for( int icell=0; icell < n_cells; icell++ ){
		for( int iface=0; iface < mesh_domain->faces_per_cell; iface++ ){
			int gcell = part2old[icell];
			int gface = mesh_domain->cell2face[gcell][iface];

			if( unregistered_faces[gface] ){
				unregistered_faces[gface] = false;
				n_faces++;
			}
		}
	}

	// Initialize mesh_size
	t_mesh_size mesh_size;

	mesh_size.n_nodes = n_nodes;
	mesh_size.n_faces = n_faces;
	mesh_size.n_cells = n_cells;

	mesh_size.nodes_per_face = mesh_domain->nodes_per_face;
	mesh_size.cells_per_face = mesh_domain->cells_per_face;
	mesh_size.nodes_per_cell = mesh_domain->nodes_per_cell;
	mesh_size.faces_per_cell = mesh_domain->faces_per_cell;

	return mesh_size;
}

//***************************************************************************************************
gmsh_mesh allocate_mesh( struct t_mesh_size mesh_size ){

	gmsh_mesh mesh;

	mesh.nodes.resize( mesh_size.n_nodes );
	mesh.faces.resize( mesh_size.n_faces );
	mesh.cells.resize( mesh_size.n_cells );

	mesh.nodes_per_face = mesh_size.nodes_per_face;
	mesh.cells_per_face = mesh_size.cells_per_face;
	mesh.nodes_per_cell = mesh_size.nodes_per_cell;
	mesh.faces_per_cell = mesh_size.faces_per_cell;

	mesh.face2node.resize( mesh.faces.size(), std::vector<int>( mesh.nodes_per_face, -1 ));
	mesh.face2node.resize( mesh.faces.size(), std::vector<int>( mesh.cells_per_face, -1 ));
	mesh.cell2node.resize( mesh.cells.size(), std::vector<int>( mesh.nodes_per_cell, -1 ));
	mesh.cell2face.resize( mesh.cells.size(), std::vector<int>( mesh.faces_per_cell, -1 ));


	mesh.cell2neigh.resize( mesh.cells.size(), std::vector<t_gmsh_neighbor>( mesh.faces_per_cell ) );
	//mesh.cell2neigh.allocate[ mesh.cells.size()][ mesh.faces_per_cell ];

	return mesh;
}

//***************************************************************************************************
//***************************************************************************************************
void update_mpi_connectivity( MPI_env    mpi_env,
							  gmsh_mesh *gmsh_domain,
							  gmsh_mesh *gmsh_local ){

	const int thousand = 1000;

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Initialization
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	// Identify MPI neighbors
	std::vector<bool> unregistered_mpi_neigh( mpi_env.size(), true );
	for( size_t icell=0; icell < gmsh_local->cells.size(); icell++ ){
		for( int iface=0; iface < gmsh_local->faces_per_cell; iface++ ){
			t_gmsh_neighbor *neigh_cell = &( gmsh_local->cell2neigh[icell][iface] );
			if( neigh_cell->type == CELL_MPI  || neigh_cell->type == CELL_PERIODIC ||
				neigh_cell->type == CELL_PERIODIC_MPI ){

				if( neigh_cell->proc < 0 )
					continue;

				if( unregistered_mpi_neigh[neigh_cell->proc] ){
					unregistered_mpi_neigh[neigh_cell->proc] = false;

					gmsh_local->mpi_neighbors.push_back( neigh_cell->proc );
				}
			}
		}
	}

	// MPI rank: global 2 local mapping
	std::map<int,int> rank_global2local;
	for( size_t ineigh=0; ineigh < gmsh_local->mpi_neighbors.size(); ineigh++ ){
		rank_global2local.insert( pair<int,int>( gmsh_local->mpi_neighbors[ineigh], (int) ineigh) );
	}

	// Old2new cell index mapping
	std::map<int,int> global2local_cells;
	for( size_t icell=0; icell < gmsh_local->cells.size(); icell++ ){
		global2local_cells.insert( pair<int,int>( gmsh_local->cells[icell].id_global, (int) icell ) );
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Exchange number of cells to be updated between each set of neighbors
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	std::vector<int>  this_num_mpi_cells( gmsh_local->mpi_neighbors.size(),  0 );
	std::vector<int> neigh_num_mpi_cells( gmsh_local->mpi_neighbors.size(), -1 );
	for( size_t icell=0; icell < gmsh_local->cells.size(); icell++ ){
		for( int iface=0; iface < gmsh_local->faces_per_cell; iface++ ){
			t_gmsh_neighbor *neigh = &( gmsh_local->cell2neigh[icell][iface] );

			if( neigh->type != CELL_MPI && neigh->type != CELL_PERIODIC_MPI )
				continue;

			this_num_mpi_cells[rank_global2local[neigh->proc]]++;
		}
	}

	for( size_t ineigh=0; ineigh < gmsh_local->mpi_neighbors.size(); ineigh++ ){
		int send_rank = gmsh_local->mpi_neighbors[ineigh];
		int send_tag  = thousand*(mpi_env.rank()+1) + (send_rank+1);
		int send_size = 1;
		MPI_Send( &this_num_mpi_cells[ineigh], send_size, MPI_INT, send_rank, send_tag, MPI_COMM_WORLD );
	}

	for( size_t ineigh=0; ineigh < gmsh_local->mpi_neighbors.size(); ineigh++ ){
		int recv_rank = gmsh_local->mpi_neighbors[ineigh];
		int recv_tag  = thousand*(recv_rank+1) + (mpi_env.rank()+1);
		int recv_size = 1;
		MPI_Recv( &neigh_num_mpi_cells[ineigh], recv_size, MPI_INT, recv_rank, recv_tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Exchange MPI cells' IDs that need to be updated
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	// Send
	for( size_t ineigh=0; ineigh < gmsh_local->mpi_neighbors.size(); ineigh++ ){

		// Skip MPI neighbors for which we have no cells to share (happens for periodic conditiosn)
		if( this_num_mpi_cells[ineigh] == 0 )
			continue;

		int neigh_rank = gmsh_local->mpi_neighbors[ineigh];

		int ind = 0;
		std::vector<int> mpi_cell_neighs( gmsh_local->cells.size(), -1 );

		for( size_t icell=0; icell < gmsh_local->cells.size(); icell++ ){
			for( int iface=0; iface < gmsh_local->faces_per_cell; iface++ ){
				t_gmsh_neighbor *neigh = &( gmsh_local->cell2neigh[icell][iface] );

				if( neigh->type != CELL_MPI && neigh->type != CELL_PERIODIC_MPI )
					continue;

				if( neigh->proc != neigh_rank )
					continue;

				mpi_cell_neighs[ind++] = neigh->id;
			}
		}

		int send_rank = neigh_rank;
		int send_tag  = thousand*(mpi_env.rank()+1) + (send_rank+1);
		int send_size = this_num_mpi_cells[ineigh];

		MPI_Send( &mpi_cell_neighs[0], send_size, MPI_INT, send_rank, send_tag, MPI_COMM_WORLD );
	}

	// Receive
	std::vector<std::vector<int>> neigh_mpi_cells;
	neigh_mpi_cells.resize( gmsh_local->mpi_neighbors.size() );
	for( size_t ineigh=0; ineigh < gmsh_local->mpi_neighbors.size(); ineigh++ ){

		int recv_rank = gmsh_local->mpi_neighbors[ineigh];
		int recv_tag  = thousand*(recv_rank+1) + (mpi_env.rank()+1);
		int recv_size = neigh_num_mpi_cells[ineigh];

		if( recv_size == 0 )
			continue;

		neigh_mpi_cells[ineigh].resize( recv_size, -1 );
		MPI_Recv( &neigh_mpi_cells[ineigh][0], recv_size, MPI_INT, recv_rank, recv_tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Perform mapping and update mpi neighbors' connectivity
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	// Map and send
	for( size_t ineigh=0; ineigh < gmsh_local->mpi_neighbors.size(); ineigh++ ){
		int neigh_rank = gmsh_local->mpi_neighbors[ineigh];

		if( neigh_mpi_cells[ineigh].size() == 0 )
			continue;

		std::vector<int> send_newcells_id( neigh_mpi_cells[ineigh].size(), -1 );

		for( size_t icell=0; icell < neigh_mpi_cells[ineigh].size(); icell++ ){
			send_newcells_id[icell] = global2local_cells[neigh_mpi_cells[ineigh][icell]];
		}

		int send_rank = neigh_rank;
		int send_tag  = thousand*(mpi_env.rank()+1) + (send_rank+1);
		int send_size = send_newcells_id.size();

		MPI_Send( &send_newcells_id[0], send_size, MPI_INT, send_rank, send_tag, MPI_COMM_WORLD );
	}

	// Receive and update
	for( size_t ineigh=0; ineigh < gmsh_local->mpi_neighbors.size(); ineigh++ ){

		// Receive updated indices
		int neigh_rank = gmsh_local->mpi_neighbors[ineigh];

		if( neigh_mpi_cells[ineigh].size() == 0 )
			continue;

		int recv_rank = neigh_rank;
		int recv_tag  = thousand*(recv_rank+1) + (mpi_env.rank()+1);
		int recv_size = neigh_num_mpi_cells[ineigh];

		std::vector<int> recv_newcells_id( recv_size, -1 );

		MPI_Recv( &recv_newcells_id[0], recv_size, MPI_INT, recv_rank, recv_tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE );

		// Update MPI connectivity
		int ind=0;
		for( size_t icell=0; icell < gmsh_local->cells.size(); icell++ ){
			for( int iface=0; iface < gmsh_local->faces_per_cell; iface++ ){
				t_gmsh_neighbor *neigh = &( gmsh_local->cell2neigh[icell][iface] );

				if( neigh->type != CELL_MPI && neigh->type != CELL_PERIODIC_MPI )
					continue;

				if( neigh->proc != neigh_rank )
					continue;

				neigh->id = recv_newcells_id[ind++];
			}
		}
	}
}

//***************************************************************************************************
//***************************************************************************************************
void update_mpi_connectivity_faces( MPI_env    mpi_env,
									gmsh_mesh *gmsh_domain,
									gmsh_mesh *gmsh_local ){

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Initialization
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	const int thousand = 1000;

	// MPI rank: global 2 local mapping
	std::map<int,int> rank_global2local;
	for( size_t ineigh=0; ineigh < gmsh_local->mpi_neighbors.size(); ineigh++ ){
		rank_global2local.insert( pair<int,int>( gmsh_local->mpi_neighbors[ineigh], (int) ineigh) );
	}

	// Old2new cell index mapping
	std::map<int,int> global2local_faces;
	for( size_t iface=0; iface < gmsh_local->faces.size(); iface++ ){
		global2local_faces.insert( pair<int,int>( gmsh_local->faces[iface].id_global,
												  gmsh_local->faces[iface].id ));
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Exchange number of cells to be updated between each set of neighbors
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	std::vector<int>  this_num_mpi_cells( gmsh_local->mpi_neighbors.size(),  0 );
	std::vector<int> neigh_num_mpi_cells( gmsh_local->mpi_neighbors.size(), -1 );
	for( size_t icell=0; icell < gmsh_local->cells.size(); icell++ ){
		for( int iface=0; iface < gmsh_local->faces_per_cell; iface++ ){
			t_gmsh_neighbor *neigh = &( gmsh_local->cell2neigh[icell][iface] );

			if( neigh->type != CELL_PERIODIC_MPI && neigh->type != CELL_PERIODIC )
				continue;

			this_num_mpi_cells[rank_global2local[neigh->proc]]++;
		}
	}

	for( size_t ineigh=0; ineigh < gmsh_local->mpi_neighbors.size(); ineigh++ ){
		int send_rank = gmsh_local->mpi_neighbors[ineigh];
		int send_tag  = thousand*(mpi_env.rank()+1) + (send_rank+1);
		int send_size = 1;
		MPI_Send( &this_num_mpi_cells[ineigh], send_size, MPI_INT, send_rank, send_tag, MPI_COMM_WORLD );
	}

	for( size_t ineigh=0; ineigh < gmsh_local->mpi_neighbors.size(); ineigh++ ){
		int recv_rank = gmsh_local->mpi_neighbors[ineigh];
		int recv_tag  = thousand*(recv_rank+1) + (mpi_env.rank()+1);
		int recv_size = 1;
		MPI_Recv( &neigh_num_mpi_cells[ineigh], recv_size, MPI_INT, recv_rank, recv_tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Exchange MPI cells' IDs that need to be updated
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	// Send
	for( size_t ineigh=0; ineigh < gmsh_local->mpi_neighbors.size(); ineigh++ ){

		// Skip MPI neighbors for which we have no cells to share (happens for periodic conditiosn)
		if( this_num_mpi_cells[ineigh] == 0 )
			continue;

		int neigh_rank = gmsh_local->mpi_neighbors[ineigh];

		int ind = 0;
		std::vector<int> mpi_faces( gmsh_local->cells.size(), -1 );

		for( size_t icell=0; icell < gmsh_local->cells.size(); icell++ ){
			for( int iface=0; iface < gmsh_local->faces_per_cell; iface++ ){
				t_gmsh_neighbor *neigh = &( gmsh_local->cell2neigh[icell][iface] );

				if( neigh->type != CELL_PERIODIC_MPI && neigh->type != CELL_PERIODIC )
					continue;

				if( neigh->proc != neigh_rank )
					continue;

				int gface = gmsh_local->cell2face[icell][iface];
				mpi_faces[ind++] = gmsh_local->faces[gface].twin_face.index;
			}
		}

		int send_rank = neigh_rank;
		int send_tag  = thousand*(mpi_env.rank()+1) + (send_rank+1);
		int send_size = this_num_mpi_cells[ineigh];

		MPI_Send( &mpi_faces[0], send_size, MPI_INT, send_rank, send_tag, MPI_COMM_WORLD );
	}

	// Receive
	std::vector<std::vector<int>> neigh_mpi_faces;
	neigh_mpi_faces.resize( gmsh_local->mpi_neighbors.size() );
	for( size_t ineigh=0; ineigh < gmsh_local->mpi_neighbors.size(); ineigh++ ){

		int recv_rank = gmsh_local->mpi_neighbors[ineigh];
		int recv_tag  = thousand*(recv_rank+1) + (mpi_env.rank()+1);
		int recv_size = neigh_num_mpi_cells[ineigh];

		if( recv_size == 0 )
			continue;

		neigh_mpi_faces[ineigh].resize( recv_size, -1 );
		MPI_Recv( &neigh_mpi_faces[ineigh][0], recv_size, MPI_INT, recv_rank, recv_tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Perform mapping and update mpi neighbors' connectivity
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	// Map and send
	for( size_t ineigh=0; ineigh < gmsh_local->mpi_neighbors.size(); ineigh++ ){
		int neigh_rank = gmsh_local->mpi_neighbors[ineigh];

		if( neigh_mpi_faces[ineigh].size() == 0 )
			continue;

		std::vector<int> send_newfaces_id( neigh_mpi_faces[ineigh].size(), -1 );

		for( size_t icell=0; icell < neigh_mpi_faces[ineigh].size(); icell++ ){
			send_newfaces_id[icell] = global2local_faces[neigh_mpi_faces[ineigh][icell]];
		}

		int send_rank = neigh_rank;
		int send_tag  = thousand*(mpi_env.rank()+1) + (send_rank+1);
		int send_size = send_newfaces_id.size();

		MPI_Send( &send_newfaces_id[0], send_size, MPI_INT, send_rank, send_tag, MPI_COMM_WORLD );
	}

	// Receive and update
	for( size_t ineigh=0; ineigh < gmsh_local->mpi_neighbors.size(); ineigh++ ){

		// Receive updated indices
		int neigh_rank = gmsh_local->mpi_neighbors[ineigh];

		if( neigh_mpi_faces[ineigh].size() == 0 )
			continue;

		int recv_rank = neigh_rank;
		int recv_tag  = thousand*(recv_rank+1) + (mpi_env.rank()+1);
		int recv_size = neigh_num_mpi_cells[ineigh];

		std::vector<int> recv_newfaces_id( recv_size, -1 );

		MPI_Recv( &recv_newfaces_id[0], recv_size, MPI_INT, recv_rank, recv_tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE );

		// Update MPI connectivity
		int ind=0;
		for( size_t icell=0; icell < gmsh_local->cells.size(); icell++ ){
			for( int iface=0; iface < gmsh_local->faces_per_cell; iface++ ){
				t_gmsh_neighbor *neigh = &( gmsh_local->cell2neigh[icell][iface] );

				if( neigh->type != CELL_PERIODIC_MPI && neigh->type != CELL_PERIODIC )
					continue;

				if( neigh->proc != neigh_rank )
					continue;


				int gface = gmsh_local->cell2face[icell][iface];

				gmsh_local->faces[gface].twin_face.index   = recv_newfaces_id[ind++];
				gmsh_local->faces[gface].twin_face.rank    = neigh_rank;
				gmsh_local->faces[gface].twin_face.cell_id = neigh->id;
			}
		}
	}
}

//***************************************************************************************************
template < typename T >
T* from_2Dvec_to_1Darray( int nrows, int ncols, std::vector<std::vector<T>> vecT ){

	T* tmp_T = new T [ nrows * ncols ];

	for( int ii=0; ii < nrows; ii++ ){
		for( int jj=0; jj < ncols; jj++ ){
			tmp_T[ii*ncols+jj] = vecT[ii][jj];
		}
	}
	return tmp_T;
}

//***************************************************************************************************
template < typename T >
std::vector<std::vector<T>> from_1Darray_to_2Dvec( int nrows, int ncols, T *arrayT ){

	std::vector<std::vector<T>> vecT( nrows, std::vector<T>( ncols ) );

	for( int ii=0; ii < nrows; ii++ ){
		for( int jj=0; jj < ncols; jj++ ){
			vecT[ii][jj] = arrayT[ii*ncols+jj];
		}
	}
	return vecT;
}





















































