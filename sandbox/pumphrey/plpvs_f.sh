## source this to set up for FMKMF with plplot.
## source this to set up for FMKMF with plplot using F.
export FMKMF_F90="F -Xf90"
export FMKMF_SPATH=.:../../lib:./plp

export FMKMF_LINKOPTS="-lplplotdX -lf2c "
