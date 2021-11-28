#
# Copyright (C) Huawei Technologies Co., Ltd. 2020. ALL RIGHTS RESERVED.
#
# See file LICENSE for terms.
#

DIR_METIS = $(HOME)/usr/metis-5.1.0
LIB_METIS = -L$(DIR_METIS)/lib -lmetis
INC_METIS = -I$(DIR_METIS)/include

LDFLAGS                      = -lm $(LIB_METIS)
CXXFLAGS                     = -O3 -g -Wall $(INC_METIS)

lib_LTLIBRARIES              = libfastmesh.la
libfastmesh_la_LDFLAGS       = -ldl -version-info $(SOVERSION)
libfastmesh_la_SOURCES       =  src/mesh_reader.cpp src/mesh_writer.cpp src/mesh_solver.cpp			\
								src/mpi_env.cpp src/mesh_partitioning.cpp src/hpath.cpp				\
								src/cfd_v0.cpp src/gmsh_reader.cpp src/mesh_preprocessor.cpp		\
								src/cfdv1_solver.cpp
libfastmesh_ladir            = $(includedir)
dist_libfastmesh_la_HEADERS  = api/constants.h api/fastmesh.h api/mesh_writer.h api/mpi_env.h		\
							   api/mesh_partitioning.h api/hpath.h api/cfdv0_elements.h				\
							   api/cfdv0_solver.h api/gmsh_reader.h api/mesh_preprocessor.h			\
							   api/cfdv1_elements.h api/cfdv1_solver.h api/cfdv1_templates.h

bin_PROGRAMS                 = lfm_preprocessor lfm_solve

# Pre-processor binary:
lfm_preprocessor_LDADD = libfastmesh.la
lfm_preprocessor_SOURCES = info/lfm_preprocessor.cpp

# Solver binary:
lfm_solve_LDADD = libfastmesh.la
lfm_solve_SOURCES = info/lfm_solve.cpp
