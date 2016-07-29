#!/bin/bash

(cd ..; scons -j4) && rm -rf dd-3d.rsf; scons
