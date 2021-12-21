#!/bin/bash

set -x -e

export EOS='/Users/shirose/optab/eos/FastChem/chem_output_table.h5'
export OPTAB='/Users/shirose/optab'
export DATABASE='/Users/shirose/database'
export MPIBIN='/opt/local/bin'
  

##### DIRECTORIES
export WORKDIR=`basename "$0" | sed -e 's/.sh//'`
\rm -rf $WORKDIR
mkdir -p $WORKDIR
mkdir -p $WORKDIR/input
mkdir -p $WORKDIR/output

##### CD TO WORKING DIRECTORY
cd $WORKDIR

##### SELECT A SINGLE LINE-SOURCE FOR EACH MOLECULAR ISOTOPOLOGUE
cat <<EOF > input/species_id.dat
	Species	Isotopologue		HITRAN		HITEMP			Exomol		
1	H2O	1H2-16O			0	HITRAN	1	HITEMP2010	0	POKAZATEL	0	BT2
1	H2O	1H2-18O			0	HITRAN	1	HITEMP2010	0	HotWat78		
1	H2O	1H2-17O			0	HITRAN	1	HITEMP2010	0	HotWat78		
1	H2O	1H-2H-16O		0	HITRAN	1	HITEMP2010	0	VTT		
1	H2O	1H-2H-18O		0 	HITRAN	1	HITEMP2010				
1	H2O	1H-2H-17O		0	HITRAN	1  	HITEMP2010				
1	H2O	2H2-16O			0	HITRAN	1	HITEMP2010				
2	CO2	12C-16O2		0	HITRAN	1	HITEMP2010	0	UCL-4000				
2	CO2	13C-16O2		0 	HITRAN	1	HITEMP2010				
2	CO2	16O-12C-18O		0	HITRAN	1	HITEMP2010				
2	CO2	16O-12C-17O		0	HITRAN	1	HITEMP2010				
2	CO2	16O-13C-18O		0	HITRAN	1	HITEMP2010				
2	CO2	16O-13C-17O		0	HITRAN	1	HITEMP2010				
2	CO2	12C-18O2		0	HITRAN	1	HITEMP2010				
2	CO2	17O-12C-18O		0	HITRAN	1	HITEMP2010				
2	CO2	12C-17O2		0	HITRAN	1	HITEMP2010				
2	CO2	13C-18O2		0	HITRAN	1	HITEMP2010				
2	CO2	18O-13C-17O		0	HITRAN  1	HITEMP2010				
2	CO2	13C-17O2		0	HITRAN	1 	HITEMP2010				
3	O3	16O3			1	HITRAN						
3	O3	16O-16O-18O		1	HITRAN						
3	O3	16O-18O-16O		1	HITRAN						
3	O3	16O-16O-17O		1	HITRAN						
3	O3	16O-17O-16O		1	HITRAN						
4	N2O	14N2-16O		0	HITRAN	1	HITEMP2019				
4	N2O	14N-15N-16O		0	HITRAN	1	HITEMP2019				
4	N2O	15N-14N-16O		0	HITRAN	1	HITEMP2019				
4	N2O	14N2-18O		0	HITRAN	1	HITEMP2019				
4	N2O	14N2-17O		0	HITRAN	1	HITEMP2019				
5	CO	12C-16O			0	HITRAN	1	HITEMP2019	0	Li2015		
5	CO	13C-16O			0 	HITRAN	1	HITEMP2019	0	Li2015		
5	CO	12C-18O			0	HITRAN	1	HITEMP2019	0	Li2015		
5	CO	12C-17O			0	HITRAN	1	HITEMP2019	0	Li2015		
5	CO	13C-18O			0	HITRAN	1	HITEMP2019	0	Li2015		
5	CO	13C-17O			0	HITRAN	1	HITEMP2019	0	Li2015		
6	CH4	12C-1H4			0	HITRAN	1	HITEMP2020	0	YT34to10	0	YT10to10
6	CH4	13C-1H4			0	HITRAN	1	HITEMP2020				
6	CH4	12C-1H3-2H		0	HITRAN	1	HITEMP2020				
6	CH4	13C-1H3-2H		0	HITRAN	1	HITEMP2020				
7	O2	16O2			1	HITRAN						
7	O2	16O-18O			1	HITRAN						
7	O2	16O-17O			1	HITRAN						
8	NO	14N-16O			0	HITRAN	1 	HITEMP2019	0	NOname		
8	NO	15N-16O			0	HITRAN	1	HITEMP2019	0	NOname		
8	NO	14N-18O			0	HITRAN	1	HITEMP2019	0	NOname		
8	NO	15N-18O		 						0	NOname		
8	NO	14N-17O								0	NOname		
8	NO	14N-18O								0	NOname		
9	SO2	32S-16O2		1	HITRAN				0	ExoAmes		
9	SO2	34S-16O2		1	HITRAN						
10	NO2	14N-16O2		0  	HITRAN	1 	HITEMP2019				
10	NO2	15N-16O2		  	HITRAN(no q)
11	NH3	14N-1H3			1	HITRAN				0	CoYuTe		0	BYTe
11	NH3	15N-1H3			1	HITRAN				0	BYTe-15		
12	HNO3	1H-14N-16O3		1	HITRAN				0	AIJS		
12	HNO3	1H-15N-16O3		1	HITRAN						
13	OH	16O-1H			0	HITRAN	1	HITEMP2010	0	MoLLIST		
13	OH	18O-1H			0	HITRAN	1	HITEMP2010				
13	OH	16O-2H			0	HITRAN	1	HITEMP2010				
14	HF	1H-19F		 	1	HITRAN				0	Coxon-Hajig		
14	HF	2H-19F			1	HITRAN						
15	HCl	1H-35Cl			1	HITRAN						
15	HCl	1H-37Cl			1	HITRAN						
15	HCl	2H-35Cl			1	HITRAN						
15	HCl	2H-37Cl			1	HITRAN						
16	HBr	1H-79Br			1	HITRAN						
16	HBr	1H-81Br			1	HITRAN						
16	HBr	2H-79Br			1	HITRAN						
16	HBr	2H-81Br			1 	HITRAN						
17	HI	1H-127I			1	HITRAN						
17	HI	2H-127I			1	HITRAN						
18	ClO	35Cl-16O		1 	HITRAN						
18	ClO	37Cl-16O		1	HITRAN						
19	OCS	16O-12C-32S		1	HITRAN						
19	OCS	16O-12C-34S		1	HITRAN						
19	OCS	16O-13C-32S		1	HITRAN						
19	OCS	16O-12C-33S		1	HITRAN						
19	OCS	18O-12C-32S		1	HITRAN						
20	H2CO	1H2-12C-16O		1	HITRAN				0	AYTY		
20	H2CO	1H2-13C-16O		1	HITRAN						
20	H2CO	1H2-12C-18O		1	HITRAN						
21	HOCl	1H-16O-35Cl		1	HITRAN						
21	HOCl	1H-16O-37Cl		1	HITRAN						
22	N2	14N2			1	HITRAN				0	WCCRMT		
22	N2	14N-15N			1	HITRAN						
23	HCN	1H-12C-14N		1	HITRAN				0	Harris		
23	HCN	1H-13C-14N		1	HITRAN				0	Larner		
23	HCN	1H-12C-15N		1	HITRAN						
24	CH3Cl	12C-1H3-35Cl		1	HITRAN				0	OYT		
24	CH3Cl	12C-1H3-37Cl		1	HITRAN				0	OYT		
25	H2O2	1H2-16O2		1	HITRAN				0	APTY		
26	C2H2	12C2-1H2		1	HITRAN				0	aCeTY		
26	C2H2	1H-12C-13C-1H		1	HITRAN						
26	C2H2	1H-12C-12C-2H		1	HITRAN						
27	C2H6	12C2-1H6		1 	HITRAN						
27	C2H6	12C-1H3-13C-1H3		1	HITRAN						
28	PH3	31P-1H3			1	HITRAN				0	SAlTY		
29	COF2	12C-16O-19F2		1	HITRAN						
29	COF2	13C-16O-19F2		1	HITRAN						
30	SF6	32S-19F6		1	HITRAN						
31	H2S	1H2-32S			1	HITRAN				0	AYT2		
31	H2S	1H2-34S			1	HITRAN						
31	H2S	1H2-33S			1	HITRAN						
32	HCOOH	1H-12C-16O-16O-1H	1	HITRAN						
33	HO2	1H-16O2			1	HITRAN						
34	OH+	16O-1H_p		 					0	MoLLIST		
35	ClONO2	35Cl-16O-14N-16O2	1	HITRAN						
35	ClONO2	37Cl-16O-14N-16O2	1	HITRAN						
36	NO+	14N-16O_p		1	HITRAN						
37	HOBr	1H-16O-79Br		1	HITRAN						
37	HOBr	1H-16O-81Br		1	HITRAN						
38	C2H4	12C2-1H4		1	HITRAN				0	MaYTY		
38	C2H4	12C-1H2-13C-1H2		1	HITRAN						
39	CH3OH	12C-1H3-16O-1H		1	HITRAN						
40	CH3Br	12C-1H3-79Br		1	HITRAN						
40	CH3Br	12C-1H3-81Br		1	HITRAN						
41	CH3CN	12C-1H3-12C-14N		1	HITRAN						
42	CF4	12C-19F4		1	HITRAN						
43	C4H2	12C4-1H2		1	HITRAN						
44	HC3N	1H-12C3-14N		1	HITRAN						
45	H2	1H2			1	HITRAN				0	RACPPK		
45	H2	1H-2H			1	HITRAN				0	ADJSAAM		
46	CS	12C-32S			1	HITRAN				0	JnK		
46	CS	12C-34S			1	HITRAN				0	JnK		
46	CS	13C-32S			1	HITRAN				0	JnK		
46	CS	12C-33S			1 	HITRAN				0	JnK		
46	CS	12C-36S								0	JnK		
46	CS	13C-33S								0	JnK		
46	CS	13C-34S								0	JnK		
46	CS	13C-36S								0	JnK		
47	SO3	32S-16O3		1	HITRAN				0	UYT2		
48	C2N2	12C2-14N2		1	HITRAN						
49	COCl2	12C-16O-35Cl2		1	HITRAN						
49	COCl2	12C-16O-35Cl-37Cl	1 	HITRAN						
50	VO	51V-16O			 					0	VOMYT		
51	AlO	27Al-16O							0	ATP		
51	AlO	26Al-16O							0	ATP		
51	AlO	27Al-17O							0	ATP		
51	AlO	27Al-18O							0	ATP		
52	MgO	24Mg-16O							0	LiTY		
52	MgO	25Mg-16O							0	LiTY		
52	MgO	26Mg-16O							0	LiTY		
52	MgO	24Mg-17O							0	LiTY		
52	MgO	24Mg-18O							0	LiTY		
53	TiO	49Ti-16O							0	Toto		0	Kurucz
53	TiO	48Ti-16O							0	Toto		0	Kurucz
53	TiO	47Ti-16O							0	Toto		0	Kurucz
53	TiO	50Ti-16O							0	Toto		0	Kurucz
53	TiO	46Ti-16O							0	Toto		0	Kurucz
54	SiO	28Si-16O							0	EBJT		
54	SiO	28Si-17O							0	EBJT		
54	SiO	28Si-18O							0	EBJT		
54	SiO	29Si-16O							0	EBJT		
54	SiO	30Si-16O							0	EBJT		
55	CaO	40Ca-16O							0	VBATHY		
56	YO	89Y-16O								0	SSYT		
57	NH	14N-1H								0	MoLLIST		
58	CH	12C-1H								0	MoLLIST		
58	CH	13C-1H								0	MoLLIST		
59	SiH	28Si-1H								0	SiGHTLY		
59	SiH	29Si-1H								0	SiGHTLY		
59	SiH	30Si-1H								0	SiGHTLY		
59	SiH	28Si-2H								0	SiGHTLY		
60	SH	32S-1H								0	GYT		0	SNaSH
60	SH	33S-1H								0	GYT		0	SNaSH
60	SH	34S-1H								0	GYT		0	SNaSH
60	SH	36S-1H								0	GYT		0	SNaSH
60	SH	32S-2H								0	GYT		0	SNaSH
61	PH	31P-1H								0	LaTY		
62	SiH2	28Si-1H2							0	CATS		
63	SiO2	28Si-16O2							00	OYT3		
64	MgH	24Mg-1H								0	Yadin		
64	MgH	25Mg-1H								0	Yadin		
64	MgH	26Mg-1H								0	Yadin		
65	NaH	23Na-1H								0	Rivlin		
65	NaH	23Na-2H								0	Rivlin		
66	AlH	27Al-1H								0	AlHambra		
66	AlH	27Al-2H								0	AlHambra		
66	AlH	26Al-1H								0	AlHambra		
67	CrH	52Cr-1H								0	MoLLIST		
68	CaH	40Ca-1H								0	Yadin		
69	BeH	9Be-1H								0	Darby-Lewis		
69	BeH	9Be-2H								0	Darby-Lewis		
69	BeH	9Be-3H								0	Darby-Lewis		
70	TiH	48Ti-1H								0	MoLLIST		
71	FeH	56Fe-1H								0	MoLLIST		
72	LiH	7Li-1H								0	CLT-LiH
73	ScH	45Sc-1H								0	LYT		
74	SiH4	28Si-1H4							00	OY2T		
75	CH3F	12C-1H3-19F							0	OYKYT		
76	AsH3	75As-1H3							0	CYT18		
77	P2H2	31P2-1H2							00	Trans		
77	P2H2	31P2-1H-2H							00	Trans		
78	PF3	31P-19F3							00	MCYTY		
79	CH3	12C-1H3								0	AYYJ		
80	PO	31P-16O								0	POPS		
81	PN	31P-14N								0	YYLT		
81	PN	31P-15N								0	YYLT		
82	KCl	39K-35Cl							0	Barton		
82	KCl	39K-37Cl							0	Barton		
82	KCl	41K-35Cl							0	Barton		
82	KCl	41K-37Cl							0	Barton		
83	NaCl	23Na-35Cl							0 	Barton		
83	NaCl	23Na-37Cl							0	Barton		
84	LiCl	7Li-35Cl							0	MoLLIST		
84	LiCl	7Li-37Cl							0	MoLLIST		
84	LiCl	6Li-35Cl							0	MoLLIST		
84	LiCl	6Li-37Cl							0	MoLLIST		
85	CN	12C-14N								0 	MoLLIST		
86	C2	12C2								0	8states		
86	C2	12C-13C								0	8states		
86	C2	13C2								0	8states		
87	CP	12C-31P								0	MoLLIST		
88	PS	31P-32S								0	POPS		
89	NS	14N-32S								0	SNaSH		
89	NS	15N-32S								0	SNaSH		
89	NS	14N-33S								0	SNaSH		
89	NS	14N-34S								0	SNaSH		
89	NS	14N-36S								0	SNaSH		
90	SiS	28Si-32S							0	UCTY		
90	SiS	28Si-33S							0	UCTY		
90	SiS	28Si-34S							0	UCTY		
90	SiS	28Si-36S							0	UCTY		
90	SiS	29Si-32S							0	UCTY		
90	SiS	29Si-34S							0	UCTY		
90	SiS	29Si-36S							0	UCTY		
90	SiS	30Si-32S							0	UCTY		
90	SiS	30Si-33S							0	UCTY		
90	SiS	30Si-34S							0	UCTY		
90	SiS	30Si-36S							0	UCTY		
90	SiS	29Si-33S							0	UCTY		
91	NaF	23Na-19F							0	MoLLIST		
92	AlCl	27Al-35Cl							0	MoLLIST		
92	AlCl	27Al-37Cl							0	MoLLIST		
93	AlF	27Al-19F							0	MoLLIST		
94	KF	39K-19F								0	MoLLIST		
94	KF	41K-19F								0	MoLLIST		
95	LiF	7Li-19F								0	MoLLIST		
95	LiF	6Li-19F								0	MoLLIST		
96	CaF	40Ca-19F							0	MoLLIST		
96	CaF	42Ca-19F							0	MoLLIST		
96	CaF	43Ca-19F							0	MoLLIST		
96	CaF	44Ca-19F							0	MoLLIST		
96	CaF	46Ca-19F							0 	MoLLIST		
96	CaF	48Ca-19F							0 	MoLLIST		
97	MgF	24Mg-19F							0	MoLLIST		
97	MgF	25Mg-19F							0	MoLLIST		
97	MgF	26Mg-19F							0	MoLLIST		
98	LiH+	7Li-1H_p							00	CLT(no pf)
99	H2+	1H-2H_p								0	ADJSAAM		0	LT
100	HeH+	4He-1H_p							0	ADJSAAM		0	Engel
100	HeH+	3He-1H_p							0	ADJSAAM		0	Engel
100	HeH+	4He-2H_p							0	ADJSAAM		0	Engel
100	HeH+	3He-2H_p							0 	ADJSAAM		0	Engel
101	H3+	1H3_p								0	MiZATeP		
101	H3+	1H2-2H_p							0	ST
999	dummy	dummy								0	dummy
EOF

##### SELECT OPACITY SOURCES TO BE CONSIDERED (1: SELECTED, 0: NOT SELECTED)
cat <<EOF > input/fort.5
&switches ! selection of opacity sources
line_molecules = 1           ! molecular lines
line_kurucz_gfpred = 1       ! Kurucz gfpred lines
line_kurucz_gfall = 1        ! Kurucz gfall lines
rayleigh_scattering_h2 = 1   ! Rayleigh scattering by H2
rayleigh_scattering_he = 1   ! Rayleigh scattering by He
rayleigh_scattering_h = 1    ! Rayleigh scattering by H
electron_scattering = 1      ! electron scattering
cia = 0                      ! Collision-induced absorption (EXPERIMENTAL)
photoion_h2 = 1              ! Photoionization by H2
photoion_topbase = 0         ! TOPbase photoionization
photoion_mathisen = 1        ! Mathisen photoionization
photoion_verner = 1          ! Verner photoionization
photoion_h_minus = 1         ! Photoionization by H-
brems_h_minus = 1            ! Bremsstrahlung by H-
brems_h2_minus = 1           ! Bremsstrahlung by H2-
brems_atomicions = 1         ! Bremsstrahlung by atomic ions
/
&cutoffs ! for line evaluation
cutoff0_voigt = 1d2 ! cutoff for Voigt profile [in wavelenth(cm^-1)]
cutoff0_gauss = 3d0 ! cutoff for Gaussian profile [in gaussian width]
delta_crit = 1d-4   ! criteria for discarding weak lines
delta_voigt = 1d0   ! criteria for adopting Voigt profile
/
&radtemp ! radiation temperature for 2-temp Planck-mean opacity
temp2 = 6000d0
/
&grid_log_const ! wavenumber grid
k_total = 100000 ! total number of grid points
grd_min = 1d0    ! min value of wavenumber grid
grd_max = 7d0    ! max value of wavenumber grid
/
&mpi_decomp ! total number of MPI processes = kprc * jprc * mprc * jprc
kprc = 1  ! number of processes in wavenumber grid (EXPERIMENTAL)
lprc = 1  ! number of processes in line loop (EXPERIMENTAL)
mprc = 1  ! number of processes in reading line-block loop (EXPERIMENTAL)
jprc = 8  ! number of processes in layer loop
/
&block_cyclic
iblock = 1
/
EOF

##### TOTAL MPI PROCS REQUIRED
kprc=`cat input/fort.5 | grep "kprc = " | awk '{print $3}'`
lprc=`cat input/fort.5 | grep "lprc = " | awk '{print $3}'`
mprc=`cat input/fort.5 | grep "mprc = " | awk '{print $3}'`
jprc=`cat input/fort.5 | grep "jprc = " | awk '{print $3}'`
nprc=$(($kprc * $lprc * $mprc * $jprc))

echo $kprc $lprc $mprc $jprc $nprc

##### USER-DEFINED EOS FILE
cp $EOS input/eos.h5

##### OPACITY DATA FILES
ln -sf $DATABASE/* input/

##### EXEC FILE
export EXEC=a.out
cp $OPTAB/src/$EXEC .
rsync -avu $OPTAB/src/*.F90 $OPTAB/src/Makefile src/

##### EXECUTION
time $MPIBIN/mpirun -np $nprc ./$EXEC
#afplay /System/Library/Sounds/Blow.aiff
#say $WORKDIR ended
