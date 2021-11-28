# libFastMesh - Finite Volume CFD Solver

## Installation
Solver's requirements:
   - Metis 5.1.0 or higher

The path to Metis directory needs to be specified in **Makefile.am** in line 7:
```[bash]
DIR_METIS = $(HOME)/usr/metis-5.1.0
```

After that, the 'fastmesh' folder needs to be placed in directory where 'bin', 'include', 'lib' etc. will be created.

When everything ready, to compile the solver next commands need to be executed:
```[bash]
autoreconf -i
./configure --prefix=${INSTALL_DIR} CC=$(which mpicc) CXX=$(which mpic++)
make
make install
```

${INSTALL_DIR} is the location where the user wants to compile the solver.

For more detailed installation process, please read **INSTALL** document.
After installation, the ${INSTLL_DIR}/bin directory should be added to PATH.
## Execution
In the 'fastmesh/testCases/Cylinder_in_Tunnel' directory, there is an example test case that can be executed. In order to execute it, the user need to be in the directory where the input ***msh*** file located and run next commands:
First:
```[bash]
mpirun ... lfm_preprocessor -i INPUT_FILE -o OUTPUT_FILE -s NUM_SUB_MESHES
```

This is the preprocess step, in which the script prepare the mesh for the solver itself. The partitioning and the renumbering of the cells happens here.
If the user want to execute given example case on 4 MPI ranks with 1 submesh on each rank, the next command needs to be executed:
```[bash]
mpirun -np 4 lfm_preprocessor -i input/Cylinder_in_Tunnel_0066K.msh -o Cylinder_in_Tunnel_0066K -s 1
```

Second:
After the preprocess, the simulation ready for the run. The command to begin the simulation is:
```[bash]
mpirun -np 4 lfm_solve -i mesh_files/<OUTPUT_FILE> -o output/
```

It is important to have the same number of processors in both steps, in preprocess step and in solver step.
For example, if the user wants to execute previous preprocessed case,  next command needs to be executed:
```[bash]
mpirun -np 4 lfm_solve -i mesh_files/Cylinder_in_Tunnel_0066K -o output/
```

The results will be saved in output folder in TecPlot format (.plt) and can be presented in VisIt or paraview.

The initial time, final time, time step, CFLmax, etc., can be determined in "mesh_solver.cpp", lines 157-166.
The properties of the fluid can be setup in "cfd_v0.cpp", lines 559-571.
After each modification, the script must be compiled once more. Inside fastmesh directory execute the command:
```[bash]
make -j install
```

In future we plan to have an input file that will contain all options that the user will need to setup before the case in order to make the use of the solver more user-friendly.

### Note
Right now, the solver can handle only five types of boundary conditions: Wall, Inlet, Outlet, Far-Field and Periodic boundary condition. The order to create correct mesh file the order of these boundary conditions must remain as shown, i.e. in ***geo*** file the user must specify the entities to be like this:
```[bash]
Physical Curve("Wall", 1) = ...
Physical Curve("Inlet", 2) = ...
Physical Curve("Outlet", 3) = ...
Physical Curve("Far-field", 4) = ...
Physical Curve("Periodic1", 5) = ...
Physical Curve("Periodic2", 6) = ...
Physical Curve("Periodic3", 7) = ...
...
Physical Curve("PeriodicN", N+4) = ...
Physical Surface("Fluid", N+5) = ...
```
Not all boundary conditions must be implemented. If user has only Wall and Far-Field boundary conditions, the Inlet and Outlet lines can be deleted, but the number of the entety must remains the same. Number of periodic boundary conditions can be as much as user needs. The "Fluid" field must be setup in the end.

If everything setup correctly, the ***msh*** file in the end will contain next lines:
```[bash]
$PhysicalNames
6
1 1 "Cylinder"
1 2 "Inlet"
1 3 "Outlet"
1 4 "Far-Field"
1 5 "Periodic1"
1 6 "Periodic2"
1 7 "Periodic3"
...
1 N+4 "PeriodicN"
2 N+5 "Fluid"
$EndPhysicalNames
```

More general approach to handle boundary conditions will be implemented in next version of libFastMesh.

## Post-process
To make the post-process more user-friendly, the **Combine_Results.py** python script was created. The script combines the results for each saved time step from all ranks and all submeshes. To execute the script, the user needs to setup the initial number of the result's file, the final number and the time step in the script, lines 6-8.
