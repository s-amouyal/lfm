#ifndef MESH_WRITER_H
#define MESH_WRITER_H

#include "elements.h"
#include "cfdv0_elements.h"
#include "gmsh_reader.h"

class MeshWriter{

	public:
		MeshWriter(){
		}

		// Mesh
		void tecplot_cell_ids  ( char* filename, gmsh_mesh *gmsh );
		void tecplot_node_ids  ( char* filename, gmsh_mesh *gmsh );
		void tecplot_partitions( char* filename, gmsh_mesh *gmsh, std::vector<int> partitions );
		void tecplot_hpath     ( char* fileName, gmsh_mesh *gmsh, std::vector<int> hpath );
		void tecplot_bnd_hpath ( char* filename, gmsh_mesh *gmsh, std::vector<int> hpath_bnd );
		void tecplot_hpath_nosplit( std::string filename, gmsh_mesh *gmsh, std::vector<int> hpath );

		// CFD
		template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT>
		void output_cfd( PRECISION time, std::string filename,
						 gmsh_mesh  *mesh,
						 fm_vector<CFDv0_cell<PRECISION, DIM_CNT, FACE_CNT>> cells_cfd );
		void output_mesh_files( std::string filename, std::vector<gmsh_mesh> gmsh_parts );
};

#endif






























