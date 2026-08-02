#####################################################
#
# Copy sacio.a from SACHOME/lib/ to current direction
#
# Here is a copy of sacio.a, but it could be not available on your system.
# Please copy a new one from SACHOME/lib/, and overwrite the one saved in this direction.
#
#####################################################

install:
	gfortran -o hilbertsac hilbertsac.f90 fftsub.f90 sacio.a
	gfortran -o evtlocate evtlocate.f90 timesub.f90 libsun.f libtau.f sacio.a -fopenmp
	chmod +x timcpot.tcl
	cp timcpot.tcl hilbertsac evtlocate ../bin

