F95=f95

.SUFFIXES: .o .f90
.f90.o:
	${F95} ${F95FLAGS} ${MODULEPATHPREFIX}${MODULEPATH} -c $<

MODULEPATHPREFIX=-I
MODULEPATH=${HOME}/mlspgs/lib/NAG.Sun
MLSPGS_LIB_PATH=${HOME}/mlspgs/lib/NAG.Sun
TOOLKIT_LIB=${PGSLIB}/libPGSTK.a $
HDF_LIB_PATH=/software/SRC/TOOLKIT-5.2.6v1.00.NAG/hdf/sun5/HDF4.1r3/lib
HDF_LIB=$(HDF_LIB_PATH)/libmfhdf.a $(HDF_LIB_PATH)/libdf.a \
	$(HDF_LIB_PATH)/libz.a $(HDF_LIB_PATH)/libjpeg.a /usr/lib/libnsl.a
MLSLIB= ${MLSPGS_LIB_PATH}/libmls.a
# F95FLAGS=-g -w=all -C=all
F95FLAGS=-O -w=all
OBJSL= MLSL1Utils.o MLSL1Common.o L0_sci_tbls.o EngTbls.o Interpolation.o \
 Calibration.o MLSPCF1.o InitPCFs.o Init_tables_module.o tree_checker.o \
 MLSL1Config.o MLSL1Rad.o OutputL1B.o Orbit.o OpenInit.o L0Utils.o SciUtils.o \
 EngUtils.o SortQualify.o Radiances.o Scan.o TkL1B.o L1BOutUtils.o \
 WriteMetaL1.o Close_Files.o

LDFLAGS= -unsharedf95 -Bstatic

all: MLSL1

OutputL1B.o : OutputL1B.f90
	${F95} ${F95FLAGS} ${MODULEPATHPREFIX}${MODULEPATH} -dusty -c $<

tree_checker.o: ${HOME}/mlspgs/srclib/tree_checker.f90
	${F95} ${F95FLAGS} ${MODULEPATHPREFIX}${MODULEPATH} -c \
	${HOME}/mlspgs/srclib/tree_checker.f90

MLSL1: $(OBJSL) MLSL1.o
	${F95} ${LDFLAGS} -o $@ ${OBJSL} MLSL1.o ${MLSLIB} ${TOOLKIT_LIB} ${HDF_LIB}

MLSL1.o: $(OBJSL)

clean:
	-rm -f *.o *.a *.M *.mod
