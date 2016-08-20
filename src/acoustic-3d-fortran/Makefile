MAKE_DEPEND=yes
CC=gcc
FC=gfortran
include ${SEPINC}/SEP.top
RESDIR=Fig
#UF90FLAGS= -fast -W0 -openmp
#UF90FLAGS= -fast -W0 -openmp -C
#UF90FLAGS= -W0 -g -O0
UF90LIBS=-openmp  -lseptravel -lsepaux
RESULTSER=qnoq spectra timeQ qtot follow timeTot sampleCompare single
RESULTSCR=full

default: install

main.x:   sinc_bob.o

install: make_time.x
	cp make_time.x ../../bin/maketime3d

source.H: source.P
	Wavelet par=source.P >$@

times.H: vel.3d.H make_time.x
	Window3d < vel.3d.H >vsmall.H j1=8 j2=8 j3=8
	./make_time.x < vsmall.H >$@ j2=4 j3=4

vel.3d.H:
	#Vel n1=600 n2=1600 n3=200 d1=8 d2=8 d3=8  rc=2000 > $@
	#echo o3=1800 >>$@
	Vel n1=400 n2=400 n3=400 d1=32 d2=32 d3=32  rc=2000 > $@

oneshot.H: ../../data/oneshot.H
	ln -s ../../data/oneshot.H .

single.H: main.x vel.3d.H oneshot.H source.H times.H
	./main.x source=source.H data=oneshot.H vel=vel.3d.H img=$@  src_x=5750 src_z=0 src_y=10 verbose=1 >junk.H qfact=200 times=times.H writeSection=1 timeBlocks=20
	#./main.x source=source.H data=oneshot.H vel=0.H vel2=1.H img=$@  src_x=5750 src_z=0 src_y=10 verbose=1 >junk.H qfact=200 times=times.H writeSection=1 timeBlocks=20
	rm -f single.H


clean: jclean

dummy:
#single.H: main.x oneshot.H source.H vel.2d.H times.H
	#./main.x source=source.H data=oneshot.H vel=vel.2d.H img=$@  src_x=5750 src_z=0 verbose=1 >junk.H qfact=200 times=times.H

#single.H: main.x oneshot.H source.H vel.2d.H times.H
	#./main.x source=source.H data=oneshot.H vel=vel.2d.H img=$@  src_x=5750 src_z=0 verbose=1 >junk.H qfact=200 times=times.H
include ${SEPINC}/SEP.bottom
