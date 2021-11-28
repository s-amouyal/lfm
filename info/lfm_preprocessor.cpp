#include <cstddef>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <unistd.h>
#include <signal.h>
#include <getopt.h>
#include <cstring>

#include "api/mesh_preprocessor.h"

static void show_usage( char *argv[] ){

    std::cerr << "Usage: " << argv[0] << " <option(s)> ARGUMENT\n"
              << "Options:\n"
              << "\t-h,--help             Show this help message\n"
              << "\t-i,--input_raw_mesh   input gmsh file            \tSpecify input mesh file path \n"
              << "\t-o,--output_file_base output LFM file path       \tSpecify output file base path\n"
              << "\t-s,--num_submesh      number of interior submesh \tSpecify number of interior submeshes\n"
              << std::endl;
}

int read_runtime_args( int argc, char *argv[],
					   std::string &filename, std::string &mesh_files, int &nSubparts );

//***************************************************************************************************
int main( int argc, char *argv[] ){

	int nSubparts = 1;
	std::string filename;
	std::string mesh_files;

	// Read run-time arguments
	if( read_runtime_args( argc, argv, filename, mesh_files, nSubparts ) != 0 ){
		return 0;
	}

	// Mesh object:
	//		- contains all geometrical information of the mesh and submeshes
	//		- includes hpath cells ordering, connectivity, nodes, etc...
	//fastmesh::Mesh<double, 2> mesh;
	//mesh.preprocess( filename, mesh_files, nSubparts );
	preprocess( filename, mesh_files, nSubparts );

	MPI_Finalize();

	//if( mesh.mpi_env.is_master() ) cout << "Pre-processing successful" << endl;
	//mesh.finalize();

	return 0;
}


//***************************************************************************************************
// Run-time arguments
//***************************************************************************************************
int read_runtime_args( int argc, char *argv[],
					   std::string &filename, std::string &mesh_files, int &nSubparts ){

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
			{"input_raw_mesh",   required_argument, 0, 0 },
			{"output_file_base", required_argument, 0, 0 },
			{"num_submesh",      required_argument, 0, 0 },
		};

		// Read argument
		int c = getopt_long(argc, argv, "i:o:s:", long_options, &option_index);

		if (c == -1) break;

		switch (c) {
			// Long option
			case 0:
            	if( option_index == 0 ){
            		filename = optarg;
            	}
            	if( option_index == 1 ){
            		mesh_files = optarg;
            	}
            	if( option_index == 2 ){
            		const char* tmp = (const char*) optarg;
            		nSubparts = atoi( tmp );
            	}
            	break;

			// Short options (if needed)
			case 'i':
				filename = optarg;
            	break;
			case 'o':
				mesh_files = optarg;
            	break;
			case 's':
				nSubparts = atoi(optarg);
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





































































