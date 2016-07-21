

/tmp/Patch_join.730.clapp:   clapp@cees-node-064   Sun Apr  9 23:15:18 2006
		junkME=0
	n2=216  o2=1600.000000  d2=20.000000   label2="hx"   unit2="Undefined"
	n3=72  o3=0.000000  d3=25.000000   label3="hy"   unit3="Undefined"
	n4=1  o4=-5.000000  d4=40.000000   label4="freq"   unit4="Undefined"
	n5=360  o5=0.000000  d5=0.166667   label5="Frequency"   unit5="Undefined"
		gff="-1"
		hff="-1"
hff=-1 gff=-1  n6=1  n7=1  n8=1  n9=1 
	n1=672  o1=0.000000  d1=20.000000   label1="cmpy"   unit1="Undefined"
		hff="-1"
		gff="-1"
		esize=8
		sets next: in="/data/temp2/bob/ELF_elf.inv_freq.H@"
		data_format="xdr_float"
o4=0


/tmp/Patch_join.18640.clapp:   clapp@cees-node-006   Mon Apr 10 14:14:10 2006
		junkME=0
	n2=216  o2=1600.000000  d2=20.000000   label2="hx"   unit2="Undefined"
	n3=72  o3=0.000000  d3=25.000000   label3="hy"   unit3="Undefined"
	n4=1  o4=-5.000000  d4=40.000000   label4="freq"   unit4="Undefined"
	n5=360  o5=0.000000  d5=0.166667   label5="Frequency"   unit5="Undefined"
		gff="-1"
		hff="-1"
hff=-1 gff=-1  n6=1  n7=1  n8=1  n9=1 
	n1=672  o1=0.000000  d1=20.000000   label1="cmpy"   unit1="Undefined"
		hff="-1"
		gff="-1"
		esize=8
		sets next: in="/data/temp2/bob/ELF_elf.inv_freq.H@"
		data_format="xdr_float"
o4=0
o4=0


Window3d:   clapp@cees-opteron.stanford.edu   Fri Oct  9 14:41:01 2009
		sets next: in="/data/temp2/bob/ELF_test.H@"
#Window #f1=0  #j1=1 #n1=672
#Window #f2=120  #j2=1 #n2=1
#Window #f3=0  #j3=1 #n3=72
#Window #f4=0  #j4=1 #n4=1
#Window #f5=0  #j5=1 #n5=360
		junkME=0
	n2=1  o2=4000.000000  d2=20.000000   label2="hx"   unit2="Undefined"
	n3=72  o3=0.000000  d3=25.000000   label3="hy"   unit3="Undefined"
	n4=1  o4=0.000000  d4=40.000000   label4="freq"   unit4="Undefined"
	n5=360  o5=0.000000  d5=0.166667   label5="Frequency"   unit5="Undefined"
		gff="-1"
		hff="-1"
hff=-1 gff=-1  n6=1  n7=1  n8=1  n9=1 
	n1=672  o1=0.000000  d1=20.000000   label1="cmpy"   unit1="Undefined"
		hff="-1"
		gff="-1"
		esize=8
		data_format="xdr_float"



Transf:   clapp@cees-opteron.stanford.edu   Fri Oct  9 14:42:28 2009
		sets next: in="/data/temp2/bob/ELF_time.H@"
		From par: maxsize=50
		From par: is2d=no
		From par: shotmig=no
		From par: shift2=0
		From par: verb=no
		From par: kill_nan=no
		From par: max_amp=1e+07
		From par: wei=yes
		esize=4
	 	n1=1024 o1=0 d1=0.00585936 label1="Frequency" 

	 	n2=672 o2=0 d2=20 label2="cmpy" 

	 	n3=1 o3=4000 d3=20 label3="hx" 

	 	n4=72 o4=0 d4=25 label4="hy" 

	 	n5=1 o5=0 d5=40 label5="freq" 

		data_format="xdr_float"



Cp3d:   clapp@cees-opteron.stanford.edu   Fri Oct  9 14:44:15 2009
		sepstr_ready_out() : sets next in="stdin"
		copy_intag="time.H"
		data_format="xdr_float"



Cp:   bob@pompei   Sat Oct 10 19:49:19 2009
		sets next: in="/net/pompei/scr1/bob/TEST_dat.H@"
		copy_intag="claudio.HH"
		data_format="xdr_float"


Create3d:   bob@pompei   Sat Oct 10 19:49:46 2009
		sets next: in="/net/pompei/scr1/bob/TEST_dat.3d.H@"
		hff="dat.3d.H@@"
	 	n2=48384 o2=1 d2=1 label2="trace_in_bin" 

		sets next: in="/net/pompei/scr1/bob/TEST_dat.H@"
		data_format="xdr_float"



Headermath:   bob@pompei   Sat Oct 10 19:54:44 2009
		sets next: in="/net/pompei/scr1/bob/TEST_math.3d.H@"
		eqn1="cmp_x-h_x"
		eqn2="h_x*2"
		hff="math.3d.H@@"
		sets next: in="/net/pompei/scr1/bob/TEST_dat.H@"
		data_format="xdr_float"



Sort3d:   bob@pompei   Sat Oct 10 20:02:31 2009
		sets next: in="/net/pompei/scr1/bob/TEST_elf.sort.H@"
		Read  from param: #synch=0
		Read  from param: #compress_tr=1
		From par: max_size max_memory=200
		Read  from param: #nkeys=2
		Read  from param: #verb=1
		Read  from param: #debug=1
		Read  from param: #just_grid=0
		Read  from param: #preserve_nd=0
		Read  from param: #plane="none"
		hff="elf.sort.H@@"
		gff="elf.sort.H@@@@"
		sets next: in="/net/pompei/scr1/bob/TEST_dat.H@"
		same_record_number=0
	n3=1 n4=1 n5=1 n6=1 n7=1 n8=1
		data_format="xdr_float"



Stack3d:   bob@pompei   Sat Oct 10 20:02:56 2009
		sets next: in="/net/pompei/scr1/bob/TEST_elf.2d.H@"
		From par: maxsize=10
		From par: normalize=no
		From par: verb=no
		From par: axes=2
		From par: keep_headers=no
		junkME=0
		gff="-1"
		hff="-1"
hff=-1 gff=-1  n4=1  n5=1  n6=1  n7=1  n8=1  n9=1 
	n1=1024  o1=0.000000  d1=0.005859   label1="Frequency"   unit1="Undefined"
		hff="-1"
		gff="-1"
		esize=4
		data_format="xdr_float"



Cp:   bob@pompei   Sat Oct 10 21:10:28 2009
		sepstr_ready_out() : sets next in="stdin"
		copy_intag="elf.2d.H"
		data_format="xdr_float"



Cp:   bob@pompei.stanford.edu   Tue Feb  2 21:25:43 2010
		sets next: in="/net/server/bob/data/elf.shot.H@"
		copy_intag="elf.shots.H"
		data_format="xdr_float"


Cp3d:   bob@gertie   Fri May  1 10:42:51 2015
		sepstr_ready_out() : sets next in="stdin"
		copy_intag="/homes/sep/bob/Src/rtm/2d_v2/elf.shots.H"
		data_format="xdr_float"



Cp:   rice@storm.stanford.edu   Thu Jul 21 13:09:23 2016
		sets next: in="/scr1/rice/3D_elf.shots.H@"
		copy_intag="/data/elf_L7d/elf.shots.HH"
		data_format="xdr_float"
	n2=70  o2=100.000000  d2=50.000000   label2="offset"   unit2="unit3"
        n3=1 o3=0 d3=25
	n4=538  o4=0.000000  d4=25.000000   label3="s_x"   unit3="unit4"
       n5=1 o5=1900 d5=25
