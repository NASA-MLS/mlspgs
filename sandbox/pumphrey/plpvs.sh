export FMKMF_F90="f90 -g -B108 -YEXT_NAMES=LCS"

echo Rebuilding Dummy SDP toolkit with ${FMKMF_F90}
${FMKMF_F90} -c ../../lib/SDPToolkitSubstitute.f90
ar -rcs libSDPToolkitSubstitute.a SDPToolkitSubstitute.o

export FMKMF_F90="f90 -B108 -YEXT_NAMES=LCS"
export FMKMF_SPATH=../../lib:.:./plp

export FMKMF_LINKOPTS="-L. -lplplotdX -lf2c -lSDPToolkitSubstitute "



