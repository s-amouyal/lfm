#include <iostream>
#include <iomanip>
#include <fstream>
#include <map>
#include <sstream>
#include "api/timing.h"

#include "api/fastmesh.h"
#include "api/mesh_writer.h"
#include "api/mpi_env.h"

using namespace std;
using namespace fastmesh;

template void Mesh<double, 2> :: initialize( std::string file_base );

//***************************************************************************************************
//									PUBLIC FUNCTION
//***************************************************************************************************

//***************************************************************************************************
// Initialize simulation:
//		- MPI framework
//		- read LFM mesh files
//		- setup mesh parameters
//***************************************************************************************************
template <typename PRECISION, unsigned DIM_CNT >
void Mesh<PRECISION, DIM_CNT> :: initialize( std::string file_base ){

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Initialize MPI environment
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	mpi_env.initialize( true );

	// Read mesh files (one per MPI rank)
	if( mpi_env.is_master() ) cout << "read mesh files..." << endl;

	read_mesh_files( mpi_env, file_base );

	for( size_t ipart=0; ipart < gmsh_parts.size(); ipart++ ){
		set_params( ipart, &gmsh_parts[ipart], gmsh_parts );
	}

	// TODO Solal: understand if all this below is necessary
	PRECISION domain[DIM_CNT] = {0.0};
	PRECISION x_max[DIM_CNT];
	PRECISION x_min[DIM_CNT];

	// Calcualte v_neigh vectors
	for( size_t ipart=0; ipart < gmsh_parts.size(); ipart++ ){

		// Skip non-boundary submeshes
		if( ipart == INDEX_BND_SUBMESH ){

			for( unsigned idim=0; idim < DIM_CNT; idim++ ){
				x_max[idim] = gmsh_parts[ipart].faces[0].xf[idim];
				x_min[idim] = gmsh_parts[ipart].faces[0].xf[idim];
			}

			for( size_t icell=0; icell < gmsh_parts[ipart].cells.size(); icell++ ){
				for( int iface=0; iface < gmsh_parts[ipart].faces_per_cell; iface++ ){

					int gface = gmsh_parts[ipart].cell2face[icell][iface];

					for( unsigned idim=0; idim < DIM_CNT; idim++ ){
						if (gmsh_parts[ipart].faces[gface].xf[idim] > x_max[idim])
							x_max[idim] = gmsh_parts[ipart].faces[gface].xf[idim];

						if (gmsh_parts[ipart].faces[gface].xf[idim] < x_min[idim])
							x_min[idim] = gmsh_parts[ipart].faces[gface].xf[idim];
					}
				}
			}

			for( unsigned idim=0; idim < DIM_CNT; idim++ ){
				domain[idim] = x_max[idim] - x_min[idim];
			}
		}
		calculate_v_neigh( ipart, &gmsh_parts[ipart], gmsh_parts, domain );
	}
}

//***************************************************************************************************
//***************************************************************************************************
//									PRIVATE FUNCTIONS
//***************************************************************************************************
//***************************************************************************************************

//***************************************************************************************************
template <typename PRECISION, unsigned DIM_CNT >
void Mesh<PRECISION, DIM_CNT> :: read_mesh_files( MPI_env mpi_env, std::string file_base ){

	// open file for each mpi rank
	ifstream fin;

	std::stringstream local_filename;
	if( mpi_env.rank() < 10 ){
		local_filename << file_base << "_p00" << mpi_env.rank() << ".lfm" ;
	}else if( mpi_env.rank() < 100 ){
		local_filename << file_base << "_p0"  << mpi_env.rank() << ".lfm" ;
	}else{
		local_filename << file_base << "_p"   << mpi_env.rank() << ".lfm" ;
	}
	fin.open( local_filename.str() );

	if( fin.fail() ){
		cout << "could not read lfm file properly." << endl << endl;
		MPI_Abort( MPI_COMM_WORLD, 924 );
	}

	// read header
	const int num_lfm_header_files = 18;
	for (int ii = 0; ii < num_lfm_header_files; ii++){
		fin.ignore( 1000, '\n' );
	}

	int n_submesh, total_nodes, total_faces, total_cells;

	fin >> n_submesh >> total_nodes >> total_faces >> total_cells;

	gmsh_parts.resize( n_submesh );

	if( mpi_env.is_master() ){
		cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
		cout << "Reading LFM file..." << endl;
		cout << "    - num submeshes =" << setw(10) << n_submesh   << endl;
		cout << "    - num nodes     =" << setw(10) << total_nodes << endl;
		cout << "    - num faces     =" << setw(10) << total_faces << endl;
		cout << "    - num cells     =" << setw(10) << total_cells << endl;
		cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
	}

	mpi_env.barrier();

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Read submesh
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	for( int ipart=0; ipart < n_submesh; ipart++ ){

		int sm_id, elmnt_type, n_nodes, n_faces, n_cells;

		// TODO: remove duplicate variables
		gmsh_parts[ipart].id = ipart;

		fin >> sm_id >> elmnt_type >> n_nodes >> n_faces >> n_cells;

		switch( elmnt_type ){
			case 3:
				gmsh_parts[ipart].nodes_per_cell = 3;
				gmsh_parts[ipart].faces_per_cell = 3;
				break;
			case 4:
				gmsh_parts[ipart].nodes_per_cell = 4;
				gmsh_parts[ipart].faces_per_cell = 4;
				break;
		}

		read_mesh_nodes( fin, n_nodes, &gmsh_parts[ipart] );
		read_mesh_faces( fin, n_faces, &gmsh_parts[ipart] );
		read_mesh_cells( fin, n_cells, &gmsh_parts[ipart] );
	}

	// Wrap up
	fin.close();
}

//***************************************************************************************************
//***************************************************************************************************
template <typename PRECISION, unsigned DIM_CNT >
void Mesh<PRECISION, DIM_CNT> :: read_mesh_nodes( ifstream &fin, int n_nodes, gmsh_mesh *mesh ){

	mesh->nodes.resize( n_nodes );

	for( size_t inode=0; inode < mesh->nodes.size(); inode++ ){
		fin >> mesh->nodes[inode].id
			>> mesh->nodes[inode].xn[0]
			>> mesh->nodes[inode].xn[1]
			>> mesh->nodes[inode].xn[2];
	}
}

//***************************************************************************************************
//***************************************************************************************************
template <typename PRECISION, unsigned DIM_CNT >
void Mesh<PRECISION, DIM_CNT> :: read_mesh_faces( ifstream &fin, int n_faces, gmsh_mesh *mesh ){

	int dummy;

	mesh->faces.resize( n_faces );
	mesh->face2node.resize( n_faces, std::vector<int> (2, -1) );
	mesh->face2cell.resize( n_faces, std::vector<int> (2, -1) );

	if( mesh->id == INDEX_BND_SUBMESH  ){
		for( size_t iface=0; iface < mesh->faces.size(); iface++ ){
			fin >> dummy
				>> mesh->face2node[iface][0]
				>> mesh->face2node[iface][1]
				>> mesh->face2cell[iface][0]
				>> mesh->face2cell[iface][1] ;

			mesh->faces[iface].id = dummy;
		}
	}else{
		for( size_t iface=0; iface < mesh->faces.size(); iface++ ){
			fin >> dummy
				>> mesh->face2node[iface][0]
				>> mesh->face2node[iface][1]
				>> mesh->face2cell[iface][0]
				>> mesh->face2cell[iface][1] ;
			mesh->faces[iface].id = dummy;
		}
	}
}

//***************************************************************************************************
//***************************************************************************************************
template <typename PRECISION, unsigned DIM_CNT >
void Mesh<PRECISION, DIM_CNT> :: read_mesh_cells( ifstream &fin, int n_cells, gmsh_mesh *mesh ){

	// Allocate vectors
	mesh->cells     .resize( n_cells );
	mesh->cell2node .resize( n_cells, std::vector<int>            ( mesh->nodes_per_cell, -1 ) );
	mesh->cell2face .resize( n_cells, std::vector<int>            ( mesh->faces_per_cell, -1 ) );
	mesh->cell2neigh.resize( n_cells, std::vector<t_gmsh_neighbor>( mesh->faces_per_cell) );

	// Read cells
	for( size_t icell=0; icell < mesh->cells.size(); icell++ ){

		int dummy;
		fin >> dummy;

		mesh->cells[icell].id = dummy;
		for( int inode=0; inode < mesh->nodes_per_cell; inode++ ){
			fin >> mesh->cell2node[icell][inode];
		}
		for( int iface=0; iface < mesh->faces_per_cell; iface++ ){
			fin >> mesh->cell2face[icell][iface];
		}

		for( int iface=0; iface < mesh->faces_per_cell; iface++ ){

			int neigh_type, neigh_sm, neigh_id, neigh_bound;
			fin >> neigh_type >> neigh_sm >> neigh_id >> neigh_bound;

			switch( neigh_type ){
				case CELL_REGULAR:
					mesh->cell2neigh[icell][iface].type = CELL_REGULAR;
					break;
				case CELL_SUBMESH:
					mesh->cell2neigh[icell][iface].type  = CELL_SUBMESH;
					break;
				case CELL_MPI:
					mesh->cell2neigh[icell][iface].type  = CELL_MPI;
					break;
				case CELL_PERIODIC:
					mesh->cell2neigh[icell][iface].type  = CELL_PERIODIC;
					break;
				case CELL_PERIODIC_MPI:
					mesh->cell2neigh[icell][iface].type  = CELL_PERIODIC_MPI;
					break;
				case CELL_NONE:
					mesh->cell2neigh[icell][iface].type  = CELL_NONE;
					break;

			}

			mesh->cell2neigh[icell][iface].sm    = neigh_sm;
			mesh->cell2neigh[icell][iface].id    = neigh_id;
			mesh->cell2neigh[icell][iface].bound = neigh_bound;
		}
	}
}

//***************************************************************************************************
//***************************************************************************************************
template <typename PRECISION, unsigned DIM_CNT >
void Mesh<PRECISION, DIM_CNT> :: set_params( int ipart, gmsh_mesh *mesh, std::vector<gmsh_mesh> mesh_parts ){

	mesh->id = ipart;

	int myrank;
	MPI_Comm_rank( MPI_COMM_WORLD, &myrank );

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Hpath cells
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	// Cell geometric parameters
	mesh->cells_vol.resize( mesh->cells.size(), 0. );
	for( size_t icell=0; icell < mesh->cells.size(); icell++ ){

		mesh->cells[icell].id = icell;
		mesh->cells[icell].xc[0] = 0.;
		mesh->cells[icell].xc[1] = 0.;
		mesh->cells[icell].xc[2] = 0.;
	}

	// Face center & Face area of Front and Back plains
	vector<vector<PRECISION>> fCtrs;
	vector<PRECISION> fArea;
	fCtrs.resize(mesh->cells.size());
	fArea.resize(mesh->cells.size());
	for( size_t i=0; i < fCtrs.size(); ++i) {
    	fCtrs[i].resize(3);
		for (int j=0; j<3; ++j) {
			fCtrs[i][j] = 0.0;
		}
	}

	for( size_t icell=0; icell < mesh->cells.size(); icell++ ){
		for( int inode=0; inode < mesh->faces_per_cell; inode++ ){
			int node_id = mesh->cell2node[icell][inode];

			fCtrs[icell][0] += mesh->nodes[node_id].xn[0];
			fCtrs[icell][1] += mesh->nodes[node_id].xn[1];
			fCtrs[icell][2] += mesh->nodes[node_id].xn[2];
		}
		fCtrs[icell][0] /= (PRECISION) mesh->faces_per_cell;
		fCtrs[icell][1] /= (PRECISION) mesh->faces_per_cell;
		fCtrs[icell][2] /= (PRECISION) mesh->faces_per_cell;

		if( mesh->faces_per_cell == 3 ){
			fArea[icell] = 0.0;
			for( int iface=0; iface < mesh->faces_per_cell; iface++ ){
				int fn0 = mesh->cell2node[icell][ iface ];
				int fn1 = mesh->cell2node[icell][ (iface+1) % mesh->faces_per_cell ];
				fArea[icell] += ( mesh->nodes[fn0].xn[0] * mesh->nodes[fn1].xn[1]
								- mesh->nodes[fn0].xn[1] * mesh->nodes[fn1].xn[0]);
			}
			fArea[icell] = 0.5 * abs( fArea[icell] );
		}else{
			PRECISION sumN[3]  = { 0.0 };
            PRECISION sumA     =   0.0;
            PRECISION sumAc[3] = { 0.0 };
			for (int inode = 0; inode < mesh->faces_per_cell; inode++) {
                int nextP = (inode == mesh->faces_per_cell-1 ? 0 : inode+1);
                vector<PRECISION> nextPoint;
                vector<PRECISION> thisPoint;
				nextPoint.resize(3);
				thisPoint.resize(3);

				int node_id = mesh->cell2node[icell][nextP];
				nextPoint[0] = mesh->nodes[node_id].xn[0];
				nextPoint[1] = mesh->nodes[node_id].xn[1];
				nextPoint[2] = mesh->nodes[node_id].xn[2];

				node_id = mesh->cell2node[icell][inode];
				thisPoint[0] = mesh->nodes[node_id].xn[0];
				thisPoint[1] = mesh->nodes[node_id].xn[1];
				thisPoint[2] = mesh->nodes[node_id].xn[2];

                PRECISION c[3] = {thisPoint[0] + nextPoint[0] + fCtrs[icell][0],
								  thisPoint[1] + nextPoint[1] + fCtrs[icell][1],
								  thisPoint[2] + nextPoint[2] + fCtrs[icell][2]};
                PRECISION n[3] = {(nextPoint[1] - thisPoint[1]) * (fCtrs[icell][2] - thisPoint[2]) - (fCtrs[icell][1] - thisPoint[1]) * (nextPoint[2] - thisPoint[2]),
								  (fCtrs[icell][0] - thisPoint[0]) * (nextPoint[2] - thisPoint[2]) - (nextPoint[0] - thisPoint[0]) * (fCtrs[icell][2] - thisPoint[2]),
								  (nextPoint[0] - thisPoint[0]) * (fCtrs[icell][1] - thisPoint[1]) - (fCtrs[icell][0] - thisPoint[0]) * (nextPoint[1] - thisPoint[1])};
                PRECISION a = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);

				sumN[0] += n[0];
				sumN[1] += n[1];
				sumN[2] += n[2];

                sumA += a;

                sumAc[0] += a*c[0];
				sumAc[1] += a*c[1];
				sumAc[2] += a*c[2];
            }

            fCtrs[icell][0] = (1.0/3.0)*sumAc[0] / sumA;
			fCtrs[icell][1] = (1.0/3.0)*sumAc[1] / sumA;
			fCtrs[icell][2] = (1.0/3.0)*sumAc[2] / sumA;
            fArea[icell] = 0.5*sqrt(sumN[0]*sumN[0] + sumN[1]*sumN[1] + sumN[2]*sumN[2]);
		}
	}

	// Face attributes
	mesh->faces_area.resize( mesh->faces.size() );
	for( size_t iface=0; iface < mesh->faces.size(); iface++ ){

		// Get global IDs of nodes
		int n0 = mesh->face2node[ iface ][ 0 ];
		int n1 = mesh->face2node[ iface ][ 1 ];

		// Face centroid
		mesh->faces[iface].xf[0] = 0.5 * ( mesh->nodes[n0].xn[0] + mesh->nodes[n1].xn[0] );
		mesh->faces[iface].xf[1] = 0.5 * ( mesh->nodes[n0].xn[1] + mesh->nodes[n1].xn[1] );

		// Face area
		mesh->faces_area[iface] = sqrt( pow(mesh->nodes[n1].xn[0]-mesh->nodes[n0].xn[0], 2)
									  + pow(mesh->nodes[n1].xn[1]-mesh->nodes[n0].xn[1], 2) );
	}

	// Compute tangeantial and normal face vectors
	mesh->v_tang.resize( mesh->faces.size(), std::vector<double>( DIM_CNT) );
	mesh->v_norm.resize( mesh->faces.size(), std::vector<double>( DIM_CNT) );
	for( size_t icell=0; icell < mesh->cells.size(); icell++ ){

		for( int inode=0; inode < mesh->faces_per_cell; inode++ ){

			int n0 = mesh->cell2node[icell][ inode ];
			int n1 = mesh->cell2node[icell][ (inode+1) % mesh->faces_per_cell ];

			int gface = mesh->cell2face[icell][ inode ];

			// Face tangeantial vector
			mesh->v_tang[gface][0] =  mesh->nodes[n1].xn[0] - mesh->nodes[n0].xn[0] ;
			mesh->v_tang[gface][1] =  mesh->nodes[n1].xn[1] - mesh->nodes[n0].xn[1] ;

			// Face normal vector
			mesh->v_norm[gface][0] =  mesh->v_tang[gface][1] ;
			mesh->v_norm[gface][1] = -mesh->v_tang[gface][0] ;
		}
	}

	mesh->v_sign.resize( mesh->cells.size(), std::vector<double>( mesh->faces_per_cell, 1.0 ) );
	for( size_t icell=0; icell < mesh->cells.size(); icell++ ){

		for( int inode=0; inode < mesh->faces_per_cell; inode++ ){

			int n0 = mesh->cell2node[icell][ inode ];
			int n1 = mesh->cell2node[icell][ (inode+1) % mesh->faces_per_cell ];

			int gface = mesh->cell2face[icell][ inode ];

			// Face tangeantial vector
			double v_tang[2];
			v_tang[0] =  mesh->nodes[n1].xn[0] - mesh->nodes[n0].xn[0] ;
			v_tang[1] =  mesh->nodes[n1].xn[1] - mesh->nodes[n0].xn[1] ;

			if( v_tang[0] != mesh->v_tang[gface][0] || v_tang[1] != mesh->v_tang[gface][1] ){
				mesh->v_sign[icell][inode] = -1.0;
			}
		}
	}

	// Cells centers & Cells volumes
	vector<vector<PRECISION>> cEst;
	cEst.resize(mesh->cells.size());
	for( size_t i=0; i < cEst.size(); ++i){
    	cEst[i].resize(DIM_CNT);
	}

	for( size_t icell=0; icell < mesh->cells.size(); icell++ ){
		cEst[icell][0] = 2.0 * fCtrs[icell][0];
		cEst[icell][1] = 2.0 * fCtrs[icell][1];

		for( int iface=0; iface < mesh->faces_per_cell; iface++ ){
			int gface = mesh->cell2face[icell][iface];

			cEst[icell][0] += mesh->faces[gface].xf[0];
			cEst[icell][1] += mesh->faces[gface].xf[1];
		}
		cEst[icell][0] /= (PRECISION) mesh->faces_per_cell + 2.0;
		cEst[icell][1] /= (PRECISION) mesh->faces_per_cell + 2.0;
	}

	mesh->cells_vol.resize( mesh->cells.size(), 0.0 );
	for( size_t icell=0; icell < mesh->cells.size(); icell++ ){
		PRECISION pyr3Vol = 0.0;
		PRECISION pc[2] = { 0.0 };

		for( int iface=0; iface < mesh->faces_per_cell; iface++ ){

			int gface = mesh->cell2face[icell][iface];

			pyr3Vol = mesh->v_sign[icell][iface] * mesh->v_norm[gface][0] * (mesh->faces[gface].xf[0] - cEst[icell][0]) +
					  mesh->v_sign[icell][iface] * mesh->v_norm[gface][1] * (mesh->faces[gface].xf[1] - cEst[icell][1]);

			pc[0] = (3.0/4.0)*mesh->faces[gface].xf[0] + (1.0/4.0)*cEst[icell][0];
			pc[1] = (3.0/4.0)*mesh->faces[gface].xf[1] + (1.0/4.0)*cEst[icell][1];

			mesh->cells[icell].xc[0] += pyr3Vol * pc[0];
			mesh->cells[icell].xc[1] += pyr3Vol * pc[1];

			mesh->cells_vol[icell] += pyr3Vol;
		}
		pyr3Vol = 0.5 * fArea[icell];										// 0.5 for extrude of 1.0 for a 2D mesh
		pc[0] = ((3.0/4.0)*fCtrs[icell][0] + (1.0/4.0)*cEst[icell][0]);
		pc[1] = ((3.0/4.0)*fCtrs[icell][1] + (1.0/4.0)*cEst[icell][1]);

		mesh->cells[icell].xc[0] += 2.0 * pyr3Vol * pc[0];
		mesh->cells[icell].xc[1] += 2.0 * pyr3Vol * pc[1];

		mesh->cells_vol[icell] += 2.0 * pyr3Vol;
	}

	for( size_t icell=0; icell < mesh->cells.size(); icell++ ){
		if (mesh->cells_vol[icell] > 1e-03)
        {
            mesh->cells[icell].xc[0] /= mesh->cells_vol[icell];
			mesh->cells[icell].xc[1] /= mesh->cells_vol[icell];
        }
        else
        {
            mesh->cells[icell].xc[0] = cEst[icell][0];
			mesh->cells[icell].xc[1] = cEst[icell][1];
        }
		mesh->cells_vol[icell] /= 3.0;
	}

	// MPI connectivity
	std::vector<bool> unregistered_mpi_neigh( mpi_env.size(), true );
	for( size_t icell=0; icell < mesh->cells.size(); icell++ ){
		for( int iface=0; iface < mesh->faces_per_cell; iface++ ){
			t_gmsh_neighbor *neigh_cell = &( mesh->cell2neigh[icell][iface] );
			if( neigh_cell->type == CELL_MPI  || neigh_cell->type == CELL_PERIODIC ||
				neigh_cell->type == CELL_PERIODIC_MPI ){

				if( neigh_cell->sm < 0 )
					continue;

				if( unregistered_mpi_neigh[neigh_cell->sm] ){
					unregistered_mpi_neigh[neigh_cell->sm] = false;

					mesh->mpi_neighbors.push_back( neigh_cell->sm );
				}
			}
		}
	}
}

//***************************************************************************************************
// Compute vector to neighbor's cell center
// 		- requires special work for MPI neighbors
//***************************************************************************************************
template <typename PRECISION, unsigned DIM_CNT >
void Mesh<PRECISION, DIM_CNT> :: calculate_v_neigh( int                    ipart,
													gmsh_mesh             *mesh,
													std::vector<gmsh_mesh> mesh_parts,
													PRECISION              domain[DIM_CNT] ){

	mesh->v_neigh.resize( mesh->faces.size(), std::vector<double>( DIM_CNT ));

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Local neighbors
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	PRECISION neigh_xc[DIM_CNT];
	for( size_t icell=0; icell < mesh->cells.size(); icell++ ){

		for( int iface=0; iface < mesh->faces_per_cell; iface++ ){

			int gface = mesh_parts[ipart].cell2face[icell][iface];
			t_gmsh_neighbor *neigh = &( mesh->cell2neigh[icell][iface] );

			switch( neigh->type ){
				case CELL_REGULAR:
					for( unsigned idim=0; idim < DIM_CNT; idim++ ){
						neigh_xc[idim] = mesh->cells[neigh->id].xc[idim];
					}
					break;
				case CELL_SUBMESH:
					for( unsigned idim=0; idim < DIM_CNT; idim++ ){
						neigh_xc[idim] = mesh_parts[neigh->sm].cells[neigh->id].xc[idim];
					}
					break;
				case CELL_NONE:
					for( unsigned idim=0; idim < DIM_CNT; idim++ ){
						neigh_xc[idim] = mesh->faces[gface].xf[idim] + (mesh->faces[gface].xf[idim] - mesh->cells[icell].xc[idim]);
					}
					break;
				case CELL_MPI:
				case CELL_PERIODIC:
				case CELL_PERIODIC_MPI:
					continue;
				default:
					break;
			}

			// Special case: periodic neighbors
			for( unsigned idim=0; idim < DIM_CNT; idim++ ){
				mesh->v_neigh[gface][idim] = neigh_xc[idim] - mesh->cells[icell].xc[idim];
			}
		}
	}

	if( mesh->id != INDEX_BND_SUBMESH )
		return;

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// MPI neighbors
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	// t_cell_center
	struct t_cell_center{
		int id;
		double x[3];
	};

	// MPI rank: global 2 local mapping
	std::map<int,int> rank_global2local;
	for( size_t ineigh=0; ineigh < mesh->mpi_neighbors.size(); ineigh++ ){
		rank_global2local.insert( pair<int,int>( mesh->mpi_neighbors[ineigh], (int) ineigh) );
	}

	// Count number of MPI/periodic neighbors
	std::vector<int> num_mpi_neighs( mesh->mpi_neighbors.size(), 0 );

	for( size_t icell=0; icell < mesh->cells.size(); icell++ ){
		for( int iface=0; iface < mesh->faces_per_cell; iface++ ){

			t_gmsh_neighbor *neigh = &( mesh->cell2neigh[icell][iface] );

			if( neigh->type != CELL_MPI && neigh->type != CELL_PERIODIC && neigh->type != CELL_PERIODIC_MPI )
				continue;

			num_mpi_neighs[rank_global2local[neigh->sm]]++;
		}
	}

	// Store cell centers
	std::vector<std::vector<struct t_cell_center>> send_cells_coord( mesh->mpi_neighbors.size() );
	for( size_t ineigh=0; ineigh < mesh->mpi_neighbors.size(); ineigh++ ){
		send_cells_coord[ineigh].resize( num_mpi_neighs[ineigh] );
	}

	const int center_nfields = 2;
	MPI_Aint  center_disps[center_nfields];
	int       center_blocklens[] = { 1, 3 };

	MPI_Datatype center_types[] = { MPI_INT, MPI_DOUBLE };

	struct t_cell_center tmp_center;
	center_disps[0] = (char*) &tmp_center.id - (char*) &tmp_center;
	center_disps[1] = (char*)  tmp_center.x  - (char*) &tmp_center;

	MPI_Datatype type_center;
	MPI_Type_create_struct ( center_nfields, center_blocklens, center_disps, center_types, &type_center );
	MPI_Type_commit( &type_center );

	//
	std::vector<int> ind( mesh->mpi_neighbors.size(), 0 );

	for( size_t icell=0; icell < mesh->cells.size(); icell++ ){
		for( int iface=0; iface < mesh->faces_per_cell; iface++ ){

			t_gmsh_neighbor *neigh = &( mesh->cell2neigh[icell][iface] );

			if( neigh->type != CELL_MPI && neigh->type != CELL_PERIODIC && neigh->type != CELL_PERIODIC_MPI )
				continue;

			int local_mpi_ind = rank_global2local[neigh->sm];

			send_cells_coord[local_mpi_ind][ind[local_mpi_ind]].id = icell;
			if( neigh->type == CELL_MPI ){
				for( unsigned idim=0; idim < DIM_CNT; idim++ ){
					send_cells_coord[local_mpi_ind][ind[local_mpi_ind]].x[idim] = mesh->cells[icell].xc[idim];
				}
			}else{
				// For periodic cells, send the distance between the cell and face center
				int gface = mesh->cell2face[icell][iface];
				for( unsigned idim=0; idim < DIM_CNT; idim++ ){
					send_cells_coord[local_mpi_ind][ind[local_mpi_ind]].x[idim] =
								mesh->cells[icell].xc[idim] - mesh->faces[gface].xf[idim];
				}
			}
			ind[local_mpi_ind]++;
		}
	}

	// Send cell centers to MPI neighbors
	int thousand = 1000;
	std::vector<MPI_Request> send_req( mesh->mpi_neighbors.size() );
	std::vector<MPI_Request> recv_req( mesh->mpi_neighbors.size() );

	for( size_t ineigh=0; ineigh < mesh->mpi_neighbors.size(); ineigh++ ){
		int send_rank = mesh->mpi_neighbors[ineigh];
		int send_tag  = thousand*(mpi_env.rank()+1) + (send_rank+1);
		int send_size = send_cells_coord[ineigh].size();
		MPI_Send( &send_cells_coord[ineigh][0], send_size, type_center, send_rank, send_tag, MPI_COMM_WORLD );
	}

	std::vector<std::vector<struct t_cell_center>> recv_cells_coord( mesh->mpi_neighbors.size() );
	for( size_t ineigh=0; ineigh < mesh->mpi_neighbors.size(); ineigh++ ){
		recv_cells_coord[ineigh].resize( num_mpi_neighs[ineigh] );
	}
	for( size_t ineigh=0; ineigh < mesh->mpi_neighbors.size(); ineigh++ ){
		int recv_rank = mesh->mpi_neighbors[ineigh];
		int recv_tag  = mpi_env.rank()+1 + thousand*(recv_rank+1);
		int recv_size = recv_cells_coord[ineigh].size();
		MPI_Recv( &recv_cells_coord[ineigh][0], recv_size, type_center, recv_rank, recv_tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
	}

	// Get neigh
	std::vector<std::map<int,int>> cell_global2local( mesh->mpi_neighbors.size() );
	for( size_t ineigh=0; ineigh < mesh->mpi_neighbors.size(); ineigh++ ){
		for( size_t icell=0; icell < recv_cells_coord[ineigh].size(); icell++ ){
			cell_global2local[ineigh].insert( pair<int,int>( recv_cells_coord[ineigh][icell].id, icell ));
		}
	}

	for( size_t icell=0; icell < mesh->cells.size(); icell++ ){
		for( int iface=0; iface < mesh->faces_per_cell; iface++ ){
			t_gmsh_neighbor *neigh = &( mesh->cell2neigh[icell][iface] );

			if( neigh->type != CELL_MPI && neigh->type != CELL_PERIODIC && neigh->type != CELL_PERIODIC_MPI )
				continue;

			int local_neigh = rank_global2local[neigh->sm];
			int local_cell  = cell_global2local[local_neigh][neigh->id];

			int gface = mesh->cell2face[icell][iface];

			if( neigh->type == CELL_MPI ){
				double neigh_xc[DIM_CNT];
				for( unsigned idim=0; idim < DIM_CNT; idim++ ){
					neigh_xc[idim] = recv_cells_coord[local_neigh][local_cell].x[idim];

					mesh->v_neigh[gface][idim] = neigh_xc[idim] - mesh->cells[icell].xc[idim];
				}
			}else{
				// Special treatment for periodic neighbors
				double neigh_dist_cell2face[DIM_CNT];
				for( unsigned idim=0; idim < DIM_CNT; idim++ ){
					neigh_dist_cell2face[idim] = recv_cells_coord[local_neigh][local_cell].x[idim];

					mesh->v_neigh[gface][idim] = (mesh->faces[gface].xf[idim] - mesh->cells[icell].xc[idim])
											   +  neigh_dist_cell2face[idim];
				}
			}
		}
	}
}

////***************************************************************************************************
//// Write Hpath for each submesh
////***************************************************************************************************
//template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
//void Submesh<PRECISION, DIM_CNT, FACE_CNT>::write_tecplot( char* filename, int myrank, unsigned ipart ){
//
//	MeshWriter mesh_writer;
//
//	std::stringstream tmp_filename;
//
//	if( myrank < 10 ){
//		tmp_filename << filename << "_p00" << myrank;
//	}else if( ipart < 100 ){
//		tmp_filename << filename << "_p0"  << myrank;
//	}else{
//		tmp_filename << filename << "_p"   << myrank;
//	}
//
//	if( ipart < 10 ){
//		tmp_filename << "_00" << ipart << ".plt";
//	}else if( ipart < 100 ){
//		tmp_filename << "_0"  << ipart << ".plt";
//	}else{
//		tmp_filename << "_"   << ipart << ".plt";
//	}
//
//	std::string filename_str = tmp_filename.str();
//	const char* m_filename = filename_str.c_str();
//	mesh_writer.tecplot_hpath<PRECISION, DIM_CNT, FACE_CNT>( (char*) m_filename, this );
//}










































