## source this to set up for FMKMF

export FMKMF_F90="epcf90 -g "
export FMKMF_SPATH=../../lib:.

echo Rebuilding Dummy SDP toolkit with ${FMKMF_F90}
${FMKMF_F90} -c ../../lib/SDPToolkitSubstitute.f90
ar -rcs libSDPToolkitSubstitute.a SDPToolkitSubstitute.o

export FMKMF_LINKOPTS="-L. -lSDPToolkitSubstitute"



