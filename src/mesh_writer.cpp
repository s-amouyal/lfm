#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <cmath>
#include "api/mesh_writer.h"

using namespace std;

template void MeshWriter::output_cfd<float , 2, 3>( float  time, std::string filename,
													gmsh_mesh *mesh,
													fm_vector<CFDv0_cell<float ,2,3>> cells_cfd );
template void MeshWriter::output_cfd<float , 2, 4>( float  time, std::string filename,
													gmsh_mesh *mesh,
													fm_vector<CFDv0_cell<float ,2,4>> cells_cfd );
template void MeshWriter::output_cfd<double, 2, 3>( double time, std::string filename,
													gmsh_mesh *mesh,
													fm_vector<CFDv0_cell<double,2,3>> cells_cfd );
template void MeshWriter::output_cfd<double, 2, 4>( double time, std::string filename,
													gmsh_mesh *mesh,
													fm_vector<CFDv0_cell<double,2,4>> cells_cfd );

//***************************************************************************************************
// Cell variables
//***************************************************************************************************
void MeshWriter::tecplot_cell_ids( char* filename, gmsh_mesh *gmsh ){

	// Open file
	ofstream fout;
	fout.open( filename );

	// Tecplot header
	fout << "TITLE = \"Example: 2D Finite Element Data\"" << endl;
	fout << "VARIABLES = \"X\", \"Y\", \"cell_id\"" << endl;

	// Assume triangular elements for now
	if( gmsh->nodes_per_cell == 3 )
	fout << "ZONE NODES = " << setw(5) << gmsh->nodes.size()
		 << ", ELEMENTS = " << setw(5) << gmsh->cells.size()
		 << ", DATAPACKING=BLOCK"
		 << ", ZONETYPE=FETRIANGLE"
		 << ", VarLocation=([3]=CellCentered)"
		 << endl;
	if( gmsh->nodes_per_cell == 4 )
	fout << "ZONE NODES = " << setw(5) << gmsh->nodes.size()
		 << ", ELEMENTS = " << setw(5) << gmsh->cells.size()
		 << ", DATAPACKING=BLOCK"
		 << ", ZONETYPE=FEQUADRILATERAL"
		 << ", VarLocation=([3]=CellCentered)"
		 << endl;

	// Node coordinates
	for( auto inode=gmsh->nodes.begin(); inode < gmsh->nodes.end(); inode++ ){
		fout << setw(15) << inode->xn[0];
	}
	fout << endl;

	for( auto inode=gmsh->nodes.begin(); inode < gmsh->nodes.end(); inode++ ){
		fout << setw(15) << inode->xn[1];
	}
	fout << endl;

	// Cell IDs
	for( size_t icell = 0; icell < gmsh->cells.size(); icell++ ){
		fout << setw(15) << icell;
	}
	fout << endl;

	// Mesh connectivity
	for( size_t icell = 0; icell < gmsh->cells.size(); icell++ ){
		//int gcell = gmsh->cells[icell].id;

		for( int inode = 0; inode < gmsh->nodes_per_cell; inode++ ){
			fout << setw(15) << gmsh->cell2node[icell][inode] + 1;
		}
		fout << endl;
	}
	fout.close();
}

//***************************************************************************************************
// Cell variables
//***************************************************************************************************
void MeshWriter::tecplot_node_ids( char* filename, gmsh_mesh *gmsh ){


	// Open file
	ofstream fout;
	fout.open( filename );

	// Tecplot header
	fout << "TITLE = \"Example: 2D Finite Element Data\"" << endl;
	fout << "VARIABLES = \"X\", \"Y\", \"index\"" << endl;

	// Assume triangular elements for now
	fout << "ZONE NODES = " << setw(5) << gmsh->nodes.size()
		 << ", ELEMENTS = " << setw(5) << gmsh->cells.size()
		 << ", DATAPACKING=BLOCK"
		 << ", ZONETYPE=FETRIANGLE"
		 << endl;

	// Node coordinates
	for( auto inode=gmsh->nodes.begin(); inode < gmsh->nodes.end(); inode++ ){
		fout << setw(15) << inode->xn[0];
	}
	fout << endl;

	for( auto inode=gmsh->nodes.begin(); inode < gmsh->nodes.end(); inode++ ){
		fout << setw(15) << inode->xn[1];
	}
	fout << endl;

	// Node IDs
	for( size_t inode=0; inode < gmsh->nodes.size(); inode++ ){
		fout << setw(15) << inode;
	}

	// Mesh connectivity
	for( size_t icell = 0; icell < gmsh->cells.size(); icell++ ){
		for( int inode = 0; inode < gmsh->nodes_per_cell; inode++ ){
			fout << setw(15) << gmsh->cell2node[icell][inode] + 1;
		}
		fout << endl;
	}
	fout.close();
}

//***************************************************************************************************
// Cell variables
//***************************************************************************************************
void MeshWriter::tecplot_partitions( char* filename, gmsh_mesh *gmsh, std::vector<int> partitions ){

	// Open file
	ofstream fout;
	fout.open( filename );

	// Tecplot header
	fout << "TITLE = \"Example: 2D Finite Element Data\"" << endl;
	fout << "VARIABLES = \"X\", \"Y\", \"cell_id\"" << endl;

	// Assume triangular elements for now
	fout << "ZONE NODES = " << setw(5) << gmsh->nodes.size()
		 << ", ELEMENTS = " << setw(5) << gmsh->cells.size()
		 << ", DATAPACKING=BLOCK"
		 << ", ZONETYPE=FETRIANGLE"
		 << ", VarLocation=([3]=CellCentered)"
		 << endl;

	// Node coordinates
	for( auto inode=gmsh->nodes.begin(); inode < gmsh->nodes.end(); inode++ ){
		fout << setw(15) << inode->xn[0];
	}
	fout << endl;

	for( auto inode=gmsh->nodes.begin(); inode < gmsh->nodes.end(); inode++ ){
		fout << setw(15) << inode->xn[1];
	}
	fout << endl;

	// Cell IDs
	for( size_t icell = 0; icell < gmsh->cells.size(); icell++ ){
		fout << setw(15) << partitions[icell];
	}
	fout << endl;

	// Mesh connectivity
	int gcell;
	for( size_t icell = 0; icell < gmsh->cells.size(); icell++ ){
		gcell = gmsh->cells[icell].id;

		for( int inode = 0; inode < 3; inode++ ){
			fout << setw(15) << gmsh->cell2node[gcell][inode] + 1;
		}
		fout << endl;
	}
	fout.close();
}

//***************************************************************************************************
// Cell variables
//***************************************************************************************************
void MeshWriter::tecplot_bnd_hpath( char* filename, gmsh_mesh *gmsh, std::vector<int> hpath_bnd ){

	// Open file
	ofstream fout;
	fout.open( filename );

	// Tecplot header
	fout << "TITLE = \"Example: 2D Finite Element Data\"" << endl;
	fout << "VARIABLES = \"X\", \"Y\", \"cell_id\"" << endl;

	// Assume triangular elements for now
	if( gmsh->nodes_per_cell == 3 )
	fout << "ZONE NODES = " << setw(5) << gmsh->nodes.size()
		 << ", ELEMENTS = " << setw(5) << gmsh->cells.size()
		 << ", DATAPACKING=BLOCK"
		 << ", ZONETYPE=FETRIANGLE"
		 << ", VarLocation=([3]=CellCentered)"
		 << endl;
	if( gmsh->nodes_per_cell == 4 )
	fout << "ZONE NODES = " << setw(5) << gmsh->nodes.size()
		 << ", ELEMENTS = " << setw(5) << gmsh->cells.size()
		 << ", DATAPACKING=BLOCK"
		 << ", ZONETYPE=FEQUADRILATERAL"
		 << ", VarLocation=([3]=CellCentered)"
		 << endl;

	// Node coordinates
	for( auto inode=gmsh->nodes.begin(); inode < gmsh->nodes.end(); inode++ ){
		fout << setw(15) << inode->xn[0];
	}
	fout << endl;

	for( auto inode=gmsh->nodes.begin(); inode < gmsh->nodes.end(); inode++ ){
		fout << setw(15) << inode->xn[1];
	}
	fout << endl;

	// Cell IDs
	std::vector<int> old2hpath( gmsh->cells.size(), -1 );
	for( size_t icell=0; icell < hpath_bnd.size(); icell++ ){
		old2hpath[hpath_bnd[icell]] = icell;
	}
	for( size_t icell = 0; icell < gmsh->cells.size(); icell++ ){
		fout << setw(15) << old2hpath[icell];
	}
	fout << endl;

	// Mesh connectivity
	int gcell;
	for( size_t icell = 0; icell < gmsh->cells.size(); icell++ ){
		gcell = gmsh->cells[icell].id;

		for( int inode = 0; inode < gmsh->nodes_per_cell; inode++ ){
			fout << setw(15) << gmsh->cell2node[gcell][inode] + 1;
		}
		fout << endl;
	}
	fout.close();
}

////***************************************************************************************************
//// Cell variables
////***************************************************************************************************
//template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT>
//void MeshWriter::tecplot_cell_vars( char* filename,
//									Submesh<PRECISION, DIM_CNT, FACE_CNT> *submesh,
//									std::vector<PRECISION> phi_bnd,
//									std::vector<PRECISION> phi_int ){
//
//	//// Open file
//	//ofstream fout;
//	//fout.open( filename );
//
//	//// Tecplot header
//	//fout << "TITLE = \"Example: 2D Finite Element Data\"" << endl;
//	//fout << "VARIABLES = \"X\", \"Y\", \"phi\"" << endl;
//
//	//// Assume triangular elements for now
//	//fout << "ZONE NODES = " << setw(5) << submesh->nodes.size()
//	//	 << ", ELEMENTS = " << setw(5) << submesh->cells.size() + submesh->boundary.size()
//	//	 << ", DATAPACKING=BLOCK"
//	//	 << ", ZONETYPE=FETRIANGLE"
//	//	 << ", VarLocation=([3]=CellCentered)"
//	//	 << endl;
//
//	//// Node coordinates
//	//for( auto inode=submesh->nodes.begin(); inode < submesh->nodes.end(); inode++ ){
//	//	fout << setw(15) << inode->x;
//	//}
//	//fout << endl;
//
//	//for( auto inode=submesh->nodes.begin(); inode < submesh->nodes.end(); inode++ ){
//	//	fout << setw(15) << inode->y;
//	//}
//	//fout << endl;
//
//	////// Cell IDs
//	////for( size_t bcell = 0; bcell < submesh->boundary.size(); bcell++ ){
//	////	fout << setw(15) << submesh->boundary[bcell].id;
//	////}
//	////for( size_t icell = 0; icell < submesh->cells.size(); icell++ ){
//	////	fout << setw(15) << submesh->cells[icell].id;
//	////}
//	////fout << endl;
//
//	//// Solution var
//	//for( size_t bcell = 0; bcell < submesh->boundary.size(); bcell++ ){
//	//	fout << setw(15) << phi_bnd[bcell];
//	//}
//	//for( size_t icell = 0; icell < submesh->cells.size(); icell++ ){
//	//	fout << setw(15) << phi_int[icell];
//	//}
//	//fout << endl;
//
//	//// Mesh connectivity
//	//int gcell;
//	//for( size_t bcell = 0; bcell < submesh->boundary.size(); bcell++ ){
//	//	gcell = submesh->boundary[bcell].id;
//
//	//	for( int inode = 0; inode < 3; inode++ ){
//	//		fout << setw(15) << submesh->link_cell2node[gcell][inode] + 1;
//	//	}
//	//	fout << endl;
//	//}
//	//for( size_t icell = 0; icell < submesh->cells.size(); icell++ ){
//	//	gcell = submesh->cells[icell].id;
//
//	//	for( int inode = 0; inode < 3; inode++ ){
//	//		fout << setw(15) << submesh->link_cell2node[gcell][inode] + 1;
//	//	}
//	//	fout << endl;
//	//}
//	//fout.close();
//}
//
////***************************************************************************************************
//// Print node IDs
////***************************************************************************************************
//template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT>
//void MeshWriter::tecplot_node_ids( char* filename, Submesh <PRECISION, DIM_CNT, FACE_CNT> *submesh ){
//
//	// Open file
//	ofstream fout;
//	fout.open( filename );
//
//	// Tecplot header
//	fout << "TITLE = \"Example: 2D Finite Element Data\"" << endl;
//	fout << "VARIABLES = \"X\", \"Y\", \"index\"" << endl;
//
//	// Assume triangular elements for now
//	fout << "ZONE NODES = " << setw(5) << submesh->nodes.size()
//		 << ", ELEMENTS = " << setw(5) << submesh->cells.size()
//		 << ", DATAPACKING=BLOCK"
//		 << ", ZONETYPE=FETRIANGLE"
//		 << endl;
//
//	// Node coordinates
//	for( auto inode=submesh->nodes.begin(); inode < submesh->nodes.end(); inode++ ){
//		fout << setw(15) << inode->x;
//	}
//	fout << endl;
//
//	for( auto inode=submesh->nodes.begin(); inode < submesh->nodes.end(); inode++ ){
//		fout << setw(15) << inode->y;
//	}
//	fout << endl;
//
//	// Node IDs
//	for( size_t inode=0; inode < submesh->nodes.size(); inode++ ){
//		fout << setw(15) << inode;
//		//fout << setw(15) << submesh->nodes[inode].id_global;
//	}
//
//	// Mesh connectivity
//	for( size_t icell = 0; icell < submesh->cells.size(); icell++ ){
//		for( unsigned inode = 0; inode < FACE_CNT; inode++ ){
//			fout << setw(15) << submesh->link_cell2node[icell][inode] + 1;
//		}
//		fout << endl;
//	}
//	fout.close();
//}

//***************************************************************************************************
void MeshWriter::tecplot_hpath( char* fileName, gmsh_mesh *gmsh, std::vector<int> hpath ){

	int myrank;
	MPI_Comm_rank( MPI_COMM_WORLD, &myrank );

	// Open file
	ofstream fout;
	fout.open( fileName );

	// Tecplot header
	fout << "TITLE = \"Example: 2D Finite Element Data\"" << endl;
	fout << "VARIABLES = \"X\", \"Y\", \"hpath\"" << endl;

	fout << "ZONE NODES = " << setw(5) << gmsh->nodes.size()
		 << ", ELEMENTS = " << setw(5) << gmsh->cells.size()
		 << ", DATAPACKING=BLOCK"
		 << ", ZONETYPE=FETRIANGLE"
		 << ", VarLocation=([3]=CellCentered)" << endl;

	// Grid and solution variables
	for( size_t inode=0; inode < gmsh->nodes.size(); inode++ ){
		fout << setw(15) << gmsh->nodes[inode].xn[0];
	}
	fout << endl;

	for( size_t inode=0; inode < gmsh->nodes.size(); inode++ ){
		fout << setw(15) << gmsh->nodes[inode].xn[1];
	}
	fout << endl;

	// Hpath
	std::vector<int> old2hpath( gmsh->cells.size(), -1 );
	for( size_t icell=0; icell < hpath.size(); icell++ ){
		if( hpath[icell] == -1 )
			break;
		old2hpath[hpath[icell]] = icell;
	}
	for( size_t icell = 0; icell < gmsh->cells.size(); icell++ ){
		fout << setw(15) << old2hpath[icell];
	}
	fout << endl;

	// Connectivity
	for( size_t icell=0; icell < gmsh->cells.size(); icell++ ){
		for( int inode=0; inode < gmsh->nodes_per_cell; inode++ ){
			fout << setw(15) << gmsh->cell2node[icell][inode] + 1;
		}
		fout << endl;
	}
	fout.close();
}

//***************************************************************************************************
void MeshWriter::tecplot_hpath_nosplit( std::string filename, gmsh_mesh *gmsh, std::vector<int> hpath ){

	int myrank;
	MPI_Comm_rank( MPI_COMM_WORLD, &myrank );

	// Initialize hpath + deadend cells
	std::vector<int> old2hpath( gmsh->cells.size(), -1 );
	for( size_t icell=0; icell < hpath.size(); icell++ ){
		if( hpath[icell] == -1 )
			break;
		old2hpath[hpath[icell]] = icell;
	}
	for( size_t icell=0; icell < gmsh->deadend_cells.size(); icell++ ){
		old2hpath[gmsh->deadend_cells[icell]] = icell + hpath.size();
		old2hpath[gmsh->deadend_cells[icell]] *= -1;
	}

	// Open file
	ofstream fout;
	fout.open( filename );

	// Tecplot header
	fout << "TITLE = \"Example: 2D Finite Element Data\"" << endl;
	fout << "VARIABLES = \"X\", \"Y\", \"hpath\"" << endl;

	// Assume triangular elements for now
	if( gmsh->nodes_per_cell == 3 )
	fout << "ZONE NODES = " << setw(5) << gmsh->nodes.size()
		 << ", ELEMENTS = " << setw(5) << gmsh->cells.size()
		 << ", DATAPACKING=BLOCK"
		 << ", ZONETYPE=FETRIANGLE"
		 << ", VarLocation=([3]=CellCentered)"
		 << endl;
	if( gmsh->nodes_per_cell == 4 )
	fout << "ZONE NODES = " << setw(5) << gmsh->nodes.size()
		 << ", ELEMENTS = " << setw(5) << gmsh->cells.size()
		 << ", DATAPACKING=BLOCK"
		 << ", ZONETYPE=FEQUADRILATERAL"
		 << ", VarLocation=([3]=CellCentered)"
		 << endl;


	// Grid and solution variables
	for( size_t inode=0; inode < gmsh->nodes.size(); inode++ ){
		fout << setw(15) << gmsh->nodes[inode].xn[0];
	}
	fout << endl;

	for( size_t inode=0; inode < gmsh->nodes.size(); inode++ ){
		fout << setw(15) << gmsh->nodes[inode].xn[1];
	}
	fout << endl;

	// Hpath
	for( size_t icell = 0; icell < gmsh->cells.size(); icell++ ){
		fout << setw(15) << old2hpath[icell];
	}
	fout << endl;

	// Connectivity
	for( size_t icell=0; icell < gmsh->cells.size(); icell++ ){
		for( int inode=0; inode < gmsh->nodes_per_cell; inode++ ){
			fout << setw(15) << gmsh->cell2node[icell][inode] + 1;
		}
		fout << endl;
	}
	fout.close();
}

//***************************************************************************************************
void MeshWriter::output_mesh_files( std::string filename, std::vector<gmsh_mesh> gmsh_parts ){


	// Open file
	ofstream fout;
	fout.open( filename );

	// Header
	fout << "#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
	fout << "# Post-processed gmsh file." << endl;
	fout << "# libfastmesh pre-processed mesh file." << endl;
	fout << "# Contains all local MPI info for elements, Hpath/deadend ordering, submeshes and connectivity" << endl;
	fout << "#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl << endl;
	fout << "# Format:" << endl;
	fout << "#" << endl;
	fout << "# n_submeshes, total_nodes, total_faces, total_cells" << endl;
	fout << "# For all submeshes{" << endl;
	fout << "# 	- submesh_id, elmnt_type, n_nodes, n_faces, n_hpath_cells, n_deadend_cells" << endl;
	fout << "#	- nodes: node_id, node_x, node_y, node_z" << endl;
	fout << "#	- faces: face_id, node1_id, node2_id, cell1_id, cell2_id" << endl;
	fout << "#	- cells: cell_id, node1_id, ..., nodeN_id," << endl;
	fout << "#	                  face1_id, ..., faceN_id," << endl;
	fout << "#			          neigh1_type, neigh1_sm, neigh1_id, ...," << endl;
	fout << "#				      neighN_type, neighN_sm, neighN_id" << endl;
	fout << "#" << endl;

	int myrank;
	MPI_Comm_rank( MPI_COMM_WORLD, &myrank );

	// Local mesh information
	int total_nodes = 0;
	int total_faces = 0;
	int total_cells = 0;

	for( size_t ipart=0; ipart < gmsh_parts.size(); ipart++ ){
		total_nodes += (int) gmsh_parts[ipart].nodes     .size() ;
		total_faces += (int) gmsh_parts[ipart].faces    .size() ;
		total_cells += (int) gmsh_parts[ipart].cells.size() ;
	}
	fout << setw(4) << gmsh_parts.size() << setw(10) << total_nodes
										 << setw(10) << total_faces
										 << setw(10) << total_cells
										 << endl;

	// Submeshes
	for( size_t ipart=0; ipart < gmsh_parts.size(); ipart++ ){
		fout << setw( 4) << ipart
			 << setw( 4) << gmsh_parts[ipart].nodes_per_cell
			 << setw(10) << gmsh_parts[ipart].nodes    .size()
			 << setw(10) << gmsh_parts[ipart].faces    .size()
			 << setw(10) << gmsh_parts[ipart].cells.size()
			 << endl;

		// Nodes
		for( size_t inode=0; inode < gmsh_parts[ipart].nodes.size(); inode++ ){
			fout << setw(10) << inode << scientific
				 << setw(18) << gmsh_parts[ipart].nodes[inode].xn[0]
				 << setw(18) << gmsh_parts[ipart].nodes[inode].xn[1]
				 << setw(18) << 0.0
				 << endl;
		}

		// Faces
		if( ipart == INDEX_BND_SUBMESH ){
			for( size_t iface=0; iface < gmsh_parts[ipart].faces.size(); iface++ ){
				fout << setw(10) << iface
				 	 << setw(10) << gmsh_parts[ipart].face2node[iface][0]
				 	 << setw(10) << gmsh_parts[ipart].face2node[iface][1]
				 	 << setw(10) << gmsh_parts[ipart].face2cell[iface][0]
				 	 << setw(10) << gmsh_parts[ipart].face2cell[iface][1]
				 	 << endl;
			}
		}else{
			for( size_t iface=0; iface < gmsh_parts[ipart].faces.size(); iface++ ){
				fout << setw(10) << iface
				 	 << setw(10) << gmsh_parts[ipart].face2node[iface][0]
				 	 << setw(10) << gmsh_parts[ipart].face2node[iface][1]
				 	 << setw(10) << gmsh_parts[ipart].face2cell[iface][0]
				 	 << setw(10) << gmsh_parts[ipart].face2cell[iface][1]
				 	 << endl;
			}
		}

		// Hpath cells
		for( size_t icell=0; icell < gmsh_parts[ipart].cells.size(); icell++ ){
			fout << setw(10) << icell;

			for( int inode=0; inode < gmsh_parts[ipart].nodes_per_cell; inode++ ){
				fout << setw(10) << gmsh_parts[ipart].cell2node[icell][inode];
			}
			for( int iface=0; iface < gmsh_parts[ipart].faces_per_cell; iface++ ){
				fout << setw(10) << gmsh_parts[ipart].cell2face[icell][iface];
			}
			for( int iface=0; iface < gmsh_parts[ipart].faces_per_cell; iface++ ){
				fout << setw( 4) << gmsh_parts[ipart].cell2neigh[icell][iface].type
					 << setw(10) << gmsh_parts[ipart].cell2neigh[icell][iface].sm
					 << setw(10) << gmsh_parts[ipart].cell2neigh[icell][iface].id
					 << setw(10) << gmsh_parts[ipart].cell2neigh[icell][iface].bound;;
			}
			fout << endl;
		}
	}

	fout.close();
}

//***************************************************************************************************
template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT>
void MeshWriter::output_cfd( PRECISION time, std::string  filename,
							 gmsh_mesh   *mesh,
							 fm_vector<CFDv0_cell<PRECISION, DIM_CNT, FACE_CNT>> cells_cfd ){

	// Open tecplot file
	ofstream fout;
	fout.open( filename );

	// Tecplot header
	fout << "TITLE = \"Example: 2D Finite Element Data\"" << endl;
	fout << "VARIABLES = \"X\", \"Y\", \"rho\", \"u\", \"v\", \"U_mag\", \"E\", \"p\"" << endl;

	// Assume triangular elements for now
	if( FACE_CNT == 3 ){
		fout << "ZONE NODES = " << setw(5) << mesh->nodes.size()
		 	 << ", ELEMENTS = " << setw(5) << cells_cfd.size()
		 	 << ", DATAPACKING=BLOCK"
			 << ", SOLUTIONTIME = " << setw(5) << time
		 	 << ", ZONETYPE=FETRIANGLE"
			 << ", VarLocation=([3,4,5,6,7,8]=CellCentered)"
		 	 << endl;
	}else{
		fout << "ZONE NODES = " << setw(5) << mesh->nodes.size()
		 	 << ", ELEMENTS = " << setw(5) << cells_cfd.size()
		 	 << ", DATAPACKING=BLOCK"
			 << ", SOLUTIONTIME = " << setw(5) << time
		 	 << ", ZONETYPE=FEQUADRILATERAL"
			 << ", VarLocation=([3,4,5,6,7,8]=CellCentered)"
		 	 << endl;
	}

	std::vector<PRECISION> r( cells_cfd.size(), 0.0 );
	std::vector<PRECISION> u( cells_cfd.size(), 0.0 );
	std::vector<PRECISION> v( cells_cfd.size(), 0.0 );
	std::vector<PRECISION> U( cells_cfd.size(), 0.0 );
	std::vector<PRECISION> E( cells_cfd.size(), 0.0 );
	std::vector<PRECISION> p( cells_cfd.size(), 0.0 );

	for( size_t icell=0; icell < cells_cfd.size(); icell++ ){
		r[icell] = cells_cfd[icell].vars.q_new[0];
		u[icell] = cells_cfd[icell].vars.q_new[1] / r[icell];
		v[icell] = cells_cfd[icell].vars.q_new[2] / r[icell];
		U[icell] = sqrt(u[icell] * u[icell] + v[icell] * v[icell]);
		E[icell] = cells_cfd[icell].vars.q_new[3] / r[icell];
		p[icell] = r[icell] * (1.4 - 1.0) * (E[icell] - 0.5 * U[icell] * U[icell]);
	}

	// Node coordinates
	for( auto inode=mesh->nodes.begin(); inode < mesh->nodes.end(); inode++ ){
		fout << setw(15) << inode->xn[0];
	}
	fout << endl;

	for( auto inode=mesh->nodes.begin(); inode < mesh->nodes.end(); inode++ ){
		fout << setw(15) << inode->xn[1];
	}
	fout << endl;

	// rho
	for( size_t icell=0; icell < cells_cfd.size(); icell++ ){
		fout << setw(15) << std::fixed << std::setprecision(10) << r[icell];
	}
	fout << endl;

	// u
	for( size_t icell=0; icell < cells_cfd.size(); icell++ ){
		fout << setw(15) << std::fixed << std::setprecision(10) << u[icell];
	}
	fout << endl;

	// v
	for( size_t icell=0; icell < cells_cfd.size(); icell++ ){
		fout << setw(15) << std::fixed << std::setprecision(10) << v[icell];
	}
	fout << endl;

	// U_mag
	for( size_t icell=0; icell < cells_cfd.size(); icell++ ){
		fout << setw(15) << std::fixed << std::setprecision(10) << U[icell];
	}
	fout << endl;

	// E
	for( size_t icell=0; icell < cells_cfd.size(); icell++ ){
		fout << setw(15) << std::fixed << std::setprecision(10) << E[icell];
	}
	fout << endl;

	// p
	for( size_t icell=0; icell < cells_cfd.size(); icell++ ){
		fout << setw(15) << std::fixed << std::setprecision(10) << p[icell];
	}
	fout << endl;

	// Mesh connectivity
	for( size_t icell = 0; icell < cells_cfd.size(); icell++ ){
		for( unsigned inode = 0; inode < FACE_CNT; inode++ ){
			fout << setw(15) << mesh->cell2node[icell][inode] + 1;
		}
		fout << endl;
	}
	fout.close();
}





















