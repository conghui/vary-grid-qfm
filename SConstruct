# vim: fdm=marker fdl=0

import os
from nvcc import *

# compiler options
additional_includes = [os.environ['RSFROOT'] + '/include']
additional_libpath  = [os.environ['RSFROOT'] + '/lib']
additional_libs     = ['rsf', 'su']

c_compiler    = ['/usr/bin/gcc',  '-std=c99']
cxx_compiler  = ['g++', '-fopenmp']
cur_cflags    = ['-Wall', '-Wextra', '-fopenmp', '-g', '-O2', '-DNO_BLAS',]
cuda_cc       = ['nvcc', '--compiler-bindir', '/usr/bin/gcc', '-arch=sm_35', '-O2']

# set the sub directories (key, value), where value is the name of directory#{{{
# please make sure the source code are in src subdirectory
dirlist = [
   ('lib', 'lib'),
   ('bin', 'bin'),
   ('src', 'src'),
]
dirs = dict(dirlist)
#}}}
# set includes and libs#{{{
inc_path          = []
libpath           = ['#lib']
libs = []
for inc in additional_includes:
  inc_path += ['-isystem', inc]
for l in additional_libpath:
  libpath += [l]
for lp in additional_libs:
  libs += [lp]
#}}}
# setup environment#{{{
cudaenv     = Environment(
                  CC      = c_compiler,
                  CCFLAGS = cur_cflags + inc_path,
                  LIBPATH = libpath,
                  LIBS    = libs,
                  LINK    = [cuda_cc] +  '-Xcompiler -fopenmp'.split())
cudaenv.Tool('nvcc', toolpath='./nvcc.py')
#}}}

# compile#{{{
for d in dirlist:
  if not 'src' in d[1]:
    continue

  SConscript(d[1] + '/SConscript',
             variant_dir = d[1].replace('src', 'build'),
             duplicate = 0,
             exports = 'cudaenv dirs')
#}}}
