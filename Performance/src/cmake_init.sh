#!/usr/bin/env bash

# source $HOME/usr/source_env.sh first before running this 

rm -rf build
mkdir -p build
cd build && \
    cmake -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpic++ -DCMAKE_CXX_FLAGS=-pg -DCMAKE_EXE_LINKER_FLAGS=-pg -DCMAKE_SHARED_LINKER_FLAGS=-pg -DCMAKE_BUILD_TYPE=RelWithDebInfo ..
