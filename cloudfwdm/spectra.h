c... spectral header file 
	integer no_line		! max no. of lines for all species
	integer no_mol		! max no. of molecules
	parameter (no_line= 1000, no_mol = 20)
	real*8 v0(no_line)
	real*4 gse(no_line)
	real*4 ist(no_line)
	real*4 wth(no_line)
	real*4 nth(no_line)
	real*4 qlg(3,no_line)
	real*4 delta(no_line)
	real*4 n1(no_line)
	real*4 gamma(no_line)
	real*4 n2(no_line)
	integer mol		! molecule mass number

	integer nmol		! no. of molecules
	integer ncnt(no_mol)	! no. of lines used for each molecule

