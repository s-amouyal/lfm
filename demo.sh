#!/bin/sh
set -e
autoreconf --install
./configure --prefix=`pwd`/build
make
make install
make dist
./build/bin/fastmesh_info
