MAKE_DEPEND=yes
include ${SEPINC}/SEP.top
RESDIR=Fig
#UF90FLAGS= -fast -W0
UF90FLAGS= -W0 -g -O0
UF90LIBS=-openmp  -lseptravel -lsepaux
RESULTSER=qnoq spectra timeQ qtot follow timeTot sampleCompare single
RESULTSCR=full

default: main.x

main.x:   sinc_bob.o

elf.shots.H: /data/elf_L7d/elf.shots.HH
	Cp  /data/elf_L7d/elf.shots.HH $@

source.H: source.P
	Wavelet par=source.P >$@

times.H: vel.3d.H make_time.x
	Window3d < vel.3d.H >vsmall.H j1=8 j2=8
	./make_time.x < vsmall.H >$@ j2=4

vel.3d.H:
	Vel n1=100 n2=100 n3=100 d1=0.0125 d2=0.0125 d3=0.0123 rc=2000 > $@

oneshot.H: elf.shots.H
	Window3d < elf.shots.H f3=230 n3=1  >x.H
	<x.H Interp o1out=0. d1out=.004 n1out=1500 >$@

test: main.x vel.3d.H oneshot.H source.H times.H
	./main.x source=source.H data=oneshot.H vel=vel.3d.H img=$@  src_x=5750 src_z=0 verbose=1 >junk.H qfact=200 times=times.H


clean: jclean

dummy:
#single.H: main.x oneshot.H source.H vel.2d.H times.H
	#./main.x source=source.H data=oneshot.H vel=vel.2d.H img=$@  src_x=5750 src_z=0 verbose=1 >junk.H qfact=200 times=times.H

#single.H: main.x oneshot.H source.H vel.2d.H times.H
	#./main.x source=source.H data=oneshot.H vel=vel.2d.H img=$@  src_x=5750 src_z=0 verbose=1 >junk.H qfact=200 times=times.H
include ${SEPINC}/SEP.bottom
