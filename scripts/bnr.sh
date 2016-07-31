#!/bin/bash

jobdir=$PWD
rootdir=../../
cd $rootdir; scons -j4 && \
  cd $jobdir; rm -rf dd-3d.rsf; scons

