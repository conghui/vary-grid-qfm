#!/bin/bash

jobdir=$PWD
rootdir=../../
cd $rootdir; scons -j4 && \
  cd $jobdir; rm -rf ./VTId-3d-GPU.rsf; scons && scons VTId-3d-GPU.rsf test

