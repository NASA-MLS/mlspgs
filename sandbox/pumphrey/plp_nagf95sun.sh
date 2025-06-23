export FMKMF_F90="nagf95 -I${TKROOT}/include  -I../../srclib"
TKROOT=/eosmls/ret4/toolkits/NAG/TOOLKIT5.2.7.3_NAG

echo Rebuilding Dummy SDP toolkit with ${FMKMF_F90}
${FMKMF_F90} -c ../../lib/SDPToolkitSubstitute.f90
ar -rcs libSDPToolkitSubstitute.a SDPToolkitSubstitute.o

#export FMKMF_F90=""
export FMKMF_SPATH=../../lib:./plp:../../lib/machines/NAG.Sun:.

export FMKMF_LINKOPTS="-L. -L/eosmls/ret4/toolkits/PLPLOT/plplot_lib/lib -lplplotdX -lSDPToolkitSubstitute "



