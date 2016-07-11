MAKE_DEPEND=yes
include ${SEPINC}/SEP.top
RESDIR=Fig
UF90FLAGS= -fast
UF90LIBS=-openmp  -lseptravel -lsepaux
RESULTSER=qnoq spectra timeQ qtot follow timeTot sampleCompare single
RESULTSCR=full


main.x:   sinc_bob.o

elf.shots.H: /data/elf_L7d/elf.shots.HH
	Cp  /data/elf_L7d/elf.shots.HH $@



oneshot.H: elf.shots.H
	Window3d < elf.shots.H f3=230 n3=1  >x.H
	<x.H Interp o1out=0. d1out=.004 n1out=1500 >$@
twoshot.H: elf.shots.H
	Window3d < elf.shots.H f3=100 j3=200 n3=2  >x.H
	<x.H Interp o1out=0. d1out=.004 n1out=1500 >$@

allshots.H: elf.shots.H
	<elf.shots.H Interp o1out=0. d1out=.004 n1out=1500 >$@


times.H: vel.2d.H make_time.x
	Window3d < vel.2d.H >vsmall.H j1=8 j2=8
	./make_time.x < vsmall.H >$@ j2=4

q.no.%.H: speedup.x
	./speedup.x qvalue=$* >$@ follow=n

q.yes.%.H: speedup.x
	./speedup.x qvalue=$* >$@ follow=y


${RESDIR}/timeQ.v:  q.no.150.H q.no.200.H q.no.250.H q.no.300.H q.no.350.H
	Cat3d  q.no.150.H q.no.200.H q.no.250.H q.no.300.H q.no.350.H |Window3d f2=1 >x.H
	Graph <x.H legend=y curvelabel="150:200:250:300:350" legendtitle="Q" label2="Speedup factor" out=$@ >a.V legendloc=tl max2=200 title=" " dash=0,1,2,3,4,5

${RESDIR}/timeTot.v: q.yes.150.H q.yes.200.H q.yes.250.H q.yes.300.H q.yes.350.H
	Cat3d q.yes.150.H q.yes.200.H q.yes.250.H q.yes.300.H q.yes.350.H |Window3d f2=1 >x.H
	Graph <x.H legend=y curvelabel="150:200:250:300:350" legendtitle="Q" label2="Speedup factor"  >a.V legendloc=tl max2=200 title=" " dash=0,1,2,3,4,5
	Cat3d  q.yes.150.H q.yes.200.H q.yes.250.H q.yes.300.H q.yes.350.H |Window3d f2=1>x.H
	Math file1=x.H exp="1./file1" >y.H
	Transp < y.H |Stack >z.H normalize=y
	Math file1=z.H exp="1./file1" >c.H
	echo o1=150 d1=50 >>c.H
	Graph <c.H label1="Q value" label2="Total speedup" title=" " >b.V  min2=0
	Vppen gridnum=2,1 a.V b.V >c.V  vpystyle=n out=$@


${RESDIR}/follow.v: speedup.x
	./speedup.x qvalue=2000 >a.H follow=y
	./speedup.x qvalue=2000 >b.H follow=y aqx=8000 aqy=4000
	Cat a.H b.H |Window3d f2=1 >x.H
	Graph <x.H legend=y curvelabel="Modeling:Migration" legendtitle="Type" label2="Speedup factor" out=$@ >a.V legendloc=tr title=" " dash=0,1,2,3,4,5 max2=200


${RESDIR}/qtot.v: q.no.150.H q.no.200.H q.no.250.H q.no.300.H q.no.350.H
	Cat3d q.no.150.H q.no.200.H q.no.250.H q.no.300.H q.no.350.H |Window3d f2=1 >x.H
	Math file1=x.H exp="1./file1" >y.H
	Transp < y.H |Stack >z.H normalize=y
	Math file1=z.H exp="1./file1" >c.H
	echo o1=150 d1=50 >>c.H
	Graph <c.H label1="Q value" label2="Total speedup" title=" " >a.V out=$@  min2=1


source.H: source.P
	Wavelet par=source.P >$@

elf.vel.H: /data/elf_L7d/elf.vel.HH
	Cp /data/elf_L7d/elf.vel.HH elf.vel.H

vel.2d.H: elf.vel.H
	Interp type=2 < elf.vel.H  o1out=0 o2out=0 d1out=.004 d2out=.004 n1out=1250 n2out=4300  >a.H maxsize=2000
	Scale dscale=1000. < a.H >$@
	echo o1=0 d1=4 d2=4 o2=0 >>$@

full.H: main.x allshots.H source.H vel.2d.H times.H
	./main.x source=source.H data=allshots.H vel=vel.2d.H img=$@  src_x=5750 src_z=0 verbose=1 >junk.H qfact=50 times=times.H

full_slow.H: main.x allshots.H source.H vel.2d.H times.H
	./main.x source=source.H data=allshots.H vel=vel.2d.H img=$@  src_x=5750 src_z=0 verbose=1 >junk.H qfact=50 times=times.H slow=1

single.H: main.x oneshot.H source.H vel.2d.H times.H
	./main.x source=source.H data=oneshot.H vel=vel.2d.H img=$@  src_x=5750 src_z=0 verbose=1 >junk.H qfact=200 times=times.H

single_slow.H: main.x oneshot.H source.H vel.2d.H times.H
	main.x source=source.H data=oneshot.H vel=vel.2d.H img=$@  src_x=5750 src_z=0 verbose=1 >junk.H qfact=2000 times=times.H slow=1

band.%.H: main.x oneshot.H source.H vel.2d.H times.H
	./main.x source=source.H data=oneshot.H vel=vel.2d.H img=a.H  src_x=5750 src_z=0 verbose=1 >junk.H qfact=200 times=times.H \
        slow=1 qtest=y writeSection=$* wfield=$@


spec.%.H: band.%.H
	Window3d < band.$*.H j1=4 j2=4 j3=8 >a.H
	 < a.H Transp maxsize=2000 plane=13 >b.H
	 < b.H Spectra j2=9999 j3=9999 >$@

wfield.q.full.H: oneshot.H source.H vel.2d.H main.x times.H
	./main.x source=source.H data=oneshot.H vel=vel.2d.H img=a.H  src_x=5750 src_z=0 verbose=1 >junk.H qfact=200 times=times.H \
        slow=1 qtest=y writeCheckpoint=7 wfield=$@

wfield.q.fast.H: oneshot.H source.H vel.2d.H main.x times.H
	./main.x source=source.H data=oneshot.H vel=vel.2d.H img=a.H  src_x=5750 src_z=0 verbose=1 >junk.H qfact=200 times=times.H \
        writeCheckpoint=7 wfield=$@


wfield.full.H: oneshot.H source.H vel.2d.H main.x times.H
	./main.x source=source.H data=oneshot.H vel=vel.2d.H img=a.H  src_x=5750 src_z=0 verbose=1 >junk.H qfact=30 times=times.H \
        slow=1  writeCheckpoint=7 wfield=$@



${RESDIR}/qnoq.v: wfield.q.full.H wfield.full.H
	Grey < wfield.full.H >a.V title=" " label1="Depth(m)" label2="X Position(m)"
	Grey < wfield.q.full.H >b.V title=" " label1="Depth(m)" label2="X Position(m)"
	Vppen vpstyle=n < a.V b.V >c.V out=$@ gridnum=2,1

${RESDIR}/sampleCompare.v: wfield.q.full.H wfield.q.fast.H
	Grey < wfield.q.full.H >a.V title=" " label1="Depth(m)" label2="X Position(m)"
	Grey < wfield.q.fast.H >b.V title=" " label1="Depth(m)" label2="X Position(m)"
	Vppen vpstyle=n < a.V b.V  >d.V out=$@ gridnum=2,1

${RESDIR}/spectra.v:  spec.2.H spec.12.H spec.22.H spec.32.H
	Cat spec.2.H spec.12.H spec.22.H spec.32.H axis=2 |Graph legend=1 curvelabel=".3:1.3:2.3:3.3" >a.V out=$@ legendloc=tr max1=120 title=" " dash=1,2,3,4,5,6 title=" "


${RESDIR}/spectra-pres.v:  spec.2.H spec.12.H spec.22.H spec.32.H
	Cat spec.2.H spec.12.H spec.22.H spec.32.H axis=2 |\
	 Graph legend=1 curvelabel=".3:1.3:2.3:3.3" >a.V out=$@ legendloc=tr max1=120 \
	  title=" " dash=1,2,3,4,5,6 title=" " titlefat=4 labelfat=4 axisfat=4 \
	   titlesz=8 labelsz=8 plotfat=5 legendtitle="Time" crowd=.6 axiscol=0 \
	  legendfat=4 legendsz=8


${RESDIR}/e-pres.v: wfield.q.full.H wfield.full.H
	Grey < wfield.full.H min2=3800 max2=8200 max1=3600 title=" " label1="Depth(m)" label2="X Position(m)" titlesz=8 axissz=8 axiscol=0 labelsz=8 labelfat=4 legendsz=8 axisfat=4 plotfat=5 legendfat=4 crowd=.6 axiscol=0 out=$@ >/dev/null

${RESDIR}/f-pres.v: wfield.q.full.H wfield.q.fast.H
	Grey < wfield.q.full.H  min2=3800 max2=8200 max1=3600 title=" " label1="Depth(m)" label2="X Position(m)" titlesz=8 axissz=8 axiscol=0 labelsz=8 labelfat=4 legendsz=8 axisfat=4 plotfat=5 legendfat=4 crowd=.6 axiscol=0 out=$@ >/dev/null

${RESDIR}/g-pres.v: wfield.q.full.H wfield.q.fast.H
	Grey < wfield.q.fast.H min2=3800 max2=8200 max1=3600 title=" " label1="Depth(m)" label2="X Position(m)"	titlesz=8 axissz=8 axiscol=0 labelsz=8 labelfat=4 legendsz=8 axisfat=4 plotfat=5 legendfat=4 crowd=.6 axiscol=0 out=$@ >/dev/null

${RESDIR}/a-pres.v:  q.no.150.H q.no.200.H q.no.250.H q.no.300.H q.no.350.H
	Cat3d  q.no.150.H q.no.200.H q.no.250.H q.no.300.H q.no.350.H |Window3d f2=1 >x.H
	Graph <x.H legend=y curvelabel="150:200:250:300:350" legendtitle="Q" \
	 titlesz=8 axissz=8 axiscol=0 labelsz=8 labelfat=4 legendsz=8 axisfat=4 plotfat=5 legendfat=4 crowd=.6 axiscol=0\
	 label2="Speedup factor" out=$@ >a.V legendloc=tl max2=200 title=" " dash=0,1,2,3,4,5

${RESDIR}/b-pres.v: q.yes.150.H q.yes.200.H q.yes.250.H q.yes.300.H q.yes.350.H
	Cat3d q.yes.150.H q.yes.200.H q.yes.250.H q.yes.300.H q.yes.350.H |Window3d f2=1 >x.H
	Graph <x.H legend=y curvelabel="150:200:250:300:350" legendtitle="Q" label2="Speedup factor"   \
		axissz=8 titlesz=8 axiscol=0 labelsz=8 labelfat=4 legendsz=8 axisfat=4 plotfat=5 legendfat=4 crowd=.6 axiscol=0\
	legendloc=tl max2=200 title=" " dash=0,1,2,3,4,5 out=$@


${RESDIR}/d-pres.v: q.yes.150.H q.yes.200.H q.yes.250.H q.yes.300.H q.yes.350.H
	Cat3d q.yes.150.H q.yes.200.H q.yes.250.H q.yes.300.H q.yes.350.H |Window3d f2=1 >x.H
	Graph <x.H legend=y curvelabel="150:200:250:300:350" \
	legendtitle="Q" label2="Speedup factor"  >a.V legendloc=tl max2=200 title=" " dash=0,1,2,3,4,5
	Cat3d  q.yes.150.H q.yes.200.H q.yes.250.H q.yes.300.H q.yes.350.H |Window3d f2=1>x.H
	Math file1=x.H exp="1./file1" >y.H
	Transp < y.H |Stack >z.H normalize=y
	Math file1=z.H exp="1./file1" >c.H
	echo o1=150 d1=50 >>c.H
	Graph <c.H label1="Q value" label2="Total speedup" title=" "  min2=0  \
		 axissz=8 titlesz=8 labelfat=4 axiscol=0  labelsz=8 legendsz=8 axisfat=4 plotfat=5 legendfat=4 crowd=.6 axiscol=0\
		 out=$@


${RESDIR}/c-pres.v: speedup.x
	./speedup.x qvalue=2000 >a.H follow=y
	./speedup.x qvalue=2000 >b.H follow=y aqx=8000 aqy=4000
	Cat a.H b.H |Window3d f2=1 >x.H
	Graph <x.H legend=y curvelabel="Modeling:Migration" \
		 axissz=8 titlesz=8 labelsz=8 axiscol=0 legendsz=8 axisfat=4 plotfat=5 legendfat=4 crowd=.6 axiscol\
		 legendtitle="Type" labelfat=4 label2="Speedup factor" out=$@ >a.V legendloc=tr title=" " dash=0,1,2,3,4,5 max2=200


${RESDIR}/%-pres.ps: ${RESDIR}/%-pres.v
	/opt/SEP/SEP8.0/bin/pstexpen Fig/$*-pres.v Fig/$*-pres.ps color=y fat=1 fatmult=1.5 invras=n force=

${RESDIR}/single.v:  single.H single_slow.H
	Grey < single_slow.H >a.V title=" " label1="Depth(m)" label2="X Position(m)"
	Grey < single.H >b.V title=" " label1="Depth(m)" label2="X Position(m)"
	Vppen vpstyle=n < b.V a.V  >d.V out=$@ gridnum=2,1

${RESDIR}/m-pres.v: times.H
	Window3d < times.H f3=80 n3=1 |Grey   title=" " label1="Depth(m)" label2="X Position(m)"  titlesz=8 axissz=8 axiscol=0 labelsz=8 labelfat=4 legendsz=8 axisfat=4 plotfat=5 legendfat=4 crowd=.6 axiscol=0 out=$@ \
	  color=j wantscalebar=y barlabel="Time(s)" >/dev/null allpos=y


${RESDIR}/h-pres.v: single.H
	Grey < single.H  title=" " label1="Depth(m)" label2="X Position(m)" max2=9000 min2=2300 titlesz=8 axissz=8 axiscol=0 labelsz=8 labelfat=4 legendsz=8 axisfat=4 plotfat=5 legendfat=4 crowd=.6 axiscol=0 out=$@ clip=7 >/dev/null

${RESDIR}/i-pres.v: single_slow.H
	Grey < single_slow.H  title=" " label1="Depth(m)" label2="X Position(m)" max2=9000 min2=2300 titlesz=8 axissz=8 axiscol=0 labelsz=8 labelfat=4 legendsz=8 axisfat=4 plotfat=5 legendfat=4 crowd=.6 axiscol=0 clip=120 out=$@ >/dev/null


${RESDIR}/full.v: diff.x full.H
	diff.x < full.H |Grey title=" " label1="Depth(m)" label2="X Position(m)" out=$@ >a.V max2=14000 min2=0

${RESDIR}/j-pres.v: diff.x full.H
	diff.x < full.H |Grey title=" " min2=0 max2=14000 titlesz=8 axissz=8 axiscol=0 labelsz=8 labelfat=4 legendsz=8 axisfat=4 plotfat=5 legendfat=4 crowd=.6 axiscol=0 label1="Depth(m)" label2="X Position(m)" out=$@ >a.V

${RESDIR}/k-pres.v: diff.x full.H
	diff.x < full.H |Window3d min1=600 max1=1200 min2=2000 max2=4000 | Grey titlesz=8 axissz=8 axiscol=0 labelsz=8 labelfat=4 legendsz=8 axisfat=4 plotfat=5 legendfat=4 crowd=.6 axiscol=0 title=" " label1="Depth(m)" label2="X Position(m)" out=$@ >a.V




clean: jclean

include ${SEPINC}/SEP.bottom
