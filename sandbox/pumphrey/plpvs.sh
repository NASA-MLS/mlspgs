## source this to set up for FMKMF with plplot.
export FMKMF_F90="f90 -B108 -YEXT_NAMES=LCS"
export FMKMF_SPATH=.:../../lib:./plp


echo Rebuilding Dummy SDP toolkit with ${FMKMF_F90}
${FMKMF_F90} -c nomod_SDPToolkit.f90
ar -rcs libSDPToolkit.a nomod_SDPToolkit.o

export FMKMF_LINKOPTS="-L. -lplplotdX -lf2c -lSDPToolkit " 



