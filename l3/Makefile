objs   = MLSL3.o MLSMessageModule.o MLSStrings.o L3DMData.o dates_module.o MLSCF.o L3CF.o OpenInit.o PCFModule.o OutputClose.o
hdflib = /software/SRC/TOOLKIT-5.2.6v1.00.NAG/hdf/sun5/HDF4.1r3/lib
hdfeos = /software/SRC/TOOLKIT-5.2.6v1.00.NAG/hdfeos/lib/sun5

MLSL3 :	${objs}	\
		/duke/nakamura/L1BOA/EPH/PGS_TD_ADEOSIItoTAI.o  \
		${PGSLIB}/libPGSTK.a	\
		${hdfeos}/libhdfeos.a	\
                ${hdfeos}/libGctp.a	\
		${hdflib}/libmfhdf.a	\
		${hdflib}/libdf.a	\
		${hdflib}/libz.a	\
		${hdflib}/libjpeg.a
	f95 -o  MLSL3	\
		${objs}		\
		/duke/nakamura/L1BOA/EPH/PGS_TD_ADEOSIItoTAI.o	\
		${PGSLIB}/libPGSTK.a	\
		${hdfeos}/libhdfeos.a	\
                ${hdfeos}/libGctp.a	\
		${hdflib}/libmfhdf.a	\
		${hdflib}/libdf.a	\
		${hdflib}/libz.a	\
		${hdflib}/libjpeg.a	\
	-L/usr/lib -lnsl

MLSL3.o : MLSL3.f90
	f95 -c MLSL3.f90

MLSL3.o : MLSMessageModule.o MLSStrings.o L3DMData.o dates_module.o MLSCF.o L3CF.o OpenInit.o PCFModule.o OutputClose.o

SDPToolkit.o : SDPToolkit.f90
	f95 -c SDPToolkit.f90

MLSPCF.o : MLSPCF.f90
	f95 -c MLSPCF.f90

MLSMessageModule.o : SDPToolkit.o MLSMessageModule.f90
	f95 -c MLSMessageModule.f90

MLSCommon.o : MLSCommon.f90
	f95 -c MLSCommon.f90

MLSStrings.o : MLSStrings.f90
	f95 -c MLSStrings.f90

L2GPData.o : L2GPData.f90
	f95 -c L2GPData.f90

Hdf.o : Hdf.f90
	f95 -c Hdf.f90

OutputL2GP.o : OutputL2GP.f90
	f95 -c -mismatch_all OutputL2GP.f90

L3DMData.o : MLSCommon.o L2GPData.o Hdf.o OutputL2GP.o MLSPCF.o L3DMData.f90
	f95 -c -mismatch_all L3DMData.f90

dates_module.o : dates_module.f90
	f95 -c dates_module.f90

MLSCF.o : MLSCF.f90
	f95 -c MLSCF.f90

L3CF.o : L3CF.f90
	f95 -c L3CF.f90

OpenInit.o : SDPToolkit.o MLSPCF.o OpenInit.f90
	f95 -c OpenInit.f90

PCFModule.o : PCFModule.f90
	f95 -c PCFModule.f90

OutputClose.o : OutputClose.f90
	f95 -c OutputClose.f90
