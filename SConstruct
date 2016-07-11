# vim: fdm=marker fdl=0

# usage:
# scons [-j n]: compile in release mode
# scons debug=1 [-j n]: compile in debug mode
#
# author: heconghui@gmail.com

import os

# compiler options
compiler_set        = 'gnu' # intel, gnu, sw, swintel
debug_mode          = 1
additional_includes = []
additional_libpath  = []
additional_libs     = []

if compiler_set == 'gnu':#{{{
  c_compiler      = ['gcc',  '-fopenmp']
  cxx_compiler    = ['g++', '-fopenmp']
  linker          = cxx_compiler
  warn_flags      = ['-Wall', '-Wextra', '-Wno-write-strings']
  optimize_flags  = ['-O2']
  debug_flags     = ['-O0', '-g']
  other_flags     = ['-DNO_BLAS', '-DMPICH_IGNORE_CXX_SEEK']
  link_flags      = ['-O1', ]
#}}}
# set the sub directories (key, value), where value is the name of directory#{{{
# please make sure the source code are in src subdirectory
dirlist = [
   ('lib', 'lib'),
   ('bin', 'bin'),
   ('tool', 'src/tool'),
   ('rsf', 'src/rsf'),
]
dirs = dict(dirlist)
#}}}
# choose debug/release flags#{{{
is_debug_mode = ARGUMENTS.get('debug', debug_mode)
if int(is_debug_mode):
  print 'Debug mode'
  cur_cflags = debug_flags + warn_flags + other_flags
else:
  print 'Release mode'
  cur_cflags = optimize_flags + warn_flags + other_flags
#}}}
# set includes and libs#{{{
inc_path          = []
libpath           = ['#lib']
libs              = []
#}}}
# setup environment#{{{
env = Environment(CC      = c_compiler,
                  CXX     = cxx_compiler,
                  LINK    = linker + link_flags,
                  CCFLAGS = cur_cflags + inc_path,
                  LIBPATH = libpath,
                  LIBS    = libs,
                  ENV     = os.environ)
#}}}
# compile#{{{
for d in dirlist:
  if not 'src/' in d[1]:
    continue

  SConscript(d[1] + '/SConscript',
             variant_dir = d[1].replace('src', 'build'),
             duplicate = 0,
             exports = 'env dirs compiler_set c_compiler cxx_compiler')
#}}}
