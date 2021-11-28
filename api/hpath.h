#ifndef HPATH_H
#define HPATH_H

#include "elements.h"
#include "gmsh_reader.h"

typedef enum hpath_algorithm_t {
	ALGO_BOUNDARY    = 0,
	ALGO_TRUE_SPIRAL = 1,
	ALGO_MESH_CENTER = 2,
} hpath_algorithm;

class Hpath{

	public:
		void get_hpath_boundary( gmsh_mesh *mesh );
		void get_hpath_interior( gmsh_mesh *mesh );
	private:

		bool find_boundary_hpath( gmsh_mesh *mesh, int starting_cell, std::vector<int> bnd_neigh );

		void find_hpath_complex( gmsh_mesh *mesh,
								 int starting_cell,
						 		 gmsh_mesh *mesh_deadend,
						 		 std::vector<int> &hpath,
						 		 int			  &deadend_count,
						 		 bool			  &successful_hpath,
						 		 int              &total_hpath_cells );

		int get_next_hpath_cell_spiral( gmsh_mesh* mesh,
										std::vector<int> hpath_neigh,
										std::vector<int> hpath_neighofneigh,
										int current_cell,
										int deadend_neighbor,
										std::vector<bool> walked_cells,
										bool &update_boundary );


		void record_hpath_neighbors( gmsh_mesh *mesh,
								 	 int     m_count,
								 	 int     next_cell,
								 	 int     deadend_neighbor,
								 	 std::vector<bool>  walked_cells,
								 	 std::vector<int>  &hpath_neigh,
								 	 std::vector<int>  &hpath_neighofneigh );

		std::vector<int> get_starting_cells( gmsh_mesh *mesh );

		void find_hpath_nodeadend_simple( gmsh_mesh *mesh,
										  int starting_cell,
										  int &total_hpath_cells,
										  std::vector<int> &hpath,
										  std::vector<int> &deadend_cells );
		bool is_neighbor_deadend ( gmsh_mesh *mesh, int current_cell, int adjc_cell, std::vector<bool> walked_cells );
		int  get_deadend_neighbor( gmsh_mesh *mesh, int current_cell, std::vector<bool> walked_cells );
		int  get_next_hpath_cell ( gmsh_mesh *mesh, int current_cell, int deadend_neighbor, std::vector<bool> walked_cells);
		int  get_next_hpath_cell_boundary( gmsh_mesh* mesh,
										   std::vector<bool>  cell_is_boundary,
										   std::vector<bool>  cell_is_boundneigh,
										   int current_cell,
										   int deadend_neighbor,
										   std::vector<bool> walked_cells,
										   bool &update_boundary );
		int get_next_bnd_hpath_cell( gmsh_mesh *mesh, int starting_cell, int current_cell,
									 int deadend_neighbor, std::vector<bool> walked_cells, std::vector<int> bnd_neigh );
		int get_deadend_neighbor_bnd( gmsh_mesh *mesh, int current_cell, int last_cell, std::vector<bool> walked_cells );

};
#endif


















































