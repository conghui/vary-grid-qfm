# vim: fdm=marker fdl=0

import os
from nvcc import *

# compiler options
additional_includes = [os.environ['RSFROOT'] + '/include']
additional_libpath  = [os.environ['RSFROOT'] + '/lib']
additional_libs     = ['rsf', 'su']

c_compiler   = ['/usr/bin/gcc', '-std=c99', '-march=native']
cxx_compiler = ['/usr/bin/g++', '-fopenmp', '-march=native']
linker       = ['/usr/bin/g++', '-fopenmp', '-march=native']
cur_cflags   = ['-Wall', '-Wextra', '-fopenmp', '-O3', '-DNO_BLAS',]
cuda_cc      = ['nvcc', '--compiler-bindir', '/usr/bin/g++', '-arch=sm_35']
cuda_flags   = ['-O2']

# set the sub directories (key, value), where value is the name of directory#{{{
# please make sure the source code are in src subdirectory
dirlist = [
   ('lib', 'lib'),
   ('bin', 'bin'),
   ('q-elasgpu', 'src/q-elasgpu'),
   ('elasgpu', 'src/elasgpu'),
   ('elascpu', 'src/elascpu'),
   ('acouscpu', 'src/acoustic-3d-c'),
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
cudaenv.Replace(NVCC = cuda_cc)
cudaenv.Prepend(NVCCFLAGS = cuda_flags)
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
