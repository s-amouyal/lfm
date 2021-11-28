#include <cstddef>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <unistd.h>
#include <signal.h>
#include <getopt.h>

#include "api/fastmesh.h"

static void show_usage( char *argv[] );
int read_runtime_args( int argc, char *argv[], std::string &file_base, std::string &output_dir );

//***************************************************************************************************
int main( int argc, char *argv[] ){

	// Read run-time arguments
	std::string file_base, output_dir;
	if( read_runtime_args( argc, argv, file_base, output_dir ) != 0 ){
		return 0;
	}

	// Mesh object:
	//		- contains all geometrical information of the mesh and submeshes
	//		- includes hpath cells ordering, connectivity, nodes, etc...
	//std::cout << std::setprecision(20);
	fastmesh::Mesh<double, 2> mesh;

	mesh.initialize( file_base );
	mesh.initialize_solver( fastmesh::SOLVER_CAAFOAM );
	mesh.solve();

	//if( mesh.mpi_env.is_master() ) cout << "Simulation finished." << endl << endl;
	mesh.finalize();

	return 0;
}

//***************************************************************************************************
// Run-time arguments interface
//***************************************************************************************************
static void show_usage( char *argv[] ){
    std::cerr << "Usage: " << argv[0] << " <option(s)> ARGUMENT\n"
              << "Options:\n"
              << "\t-h,--help          \tShow this help message\n"
              << "\t-i,--input_mesh    \tinput LFM file (preprocessed)  \tSpecify input mesh file path \n"
              << "\t-o,--output_dir    \toutput LFM file path           \tSpecify output file base path\n"
              << std::endl;
}

//***************************************************************************************************
// Run-time arguments
//***************************************************************************************************
int read_runtime_args( int argc, char *argv[], std::string &file_base, std::string &output_dir ){

	// Check if number of arguments is correct
	if( argc < 2 ){
		show_usage( argv );

		return 1;
	}

	// Go over options
	while (1) {

		// Define required arguments
		int option_index = 0;
		static struct option long_options[] = {
			{"input_mesh", required_argument, 0, 0 },
			{"output_dir", required_argument, 0, 0 },
		};

		// Read argument
		int c = getopt_long(argc, argv, "i:o:", long_options, &option_index);

		if (c == -1) break;

		switch (c) {
			// Long option
			case 0:
            	if( option_index == 0 ){
            		file_base   = optarg;
            	}
            	if( option_index == 1 ){
            		output_dir = optarg;
            	}
            	break;

			// Short options (if needed)
			case 'i':
				file_base = optarg;
            	break;
			case 'o':
				output_dir = optarg;
            	break;
        }
	}

	// Not sure what this does :)
	if (optind < argc) {
		printf("non-option ARGV-elements: ");
   		while (optind < argc)
			printf("%s ", argv[optind++]);
		printf("\n");
	}

	return 0;
}







































































