#!/bin/bash

cd ..; scons -j4 && cd job; rm -rf dd-3d.rsf; scons
