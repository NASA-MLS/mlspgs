hdflib = /software/SRC/TOOLKIT-5.2.6v1.00.NAG/hdf/sun5/HDF4.1r3/lib
hdfeos = /software/SRC/TOOLKIT-5.2.6v1.00.NAG/hdfeos/lib/sun5
mlslib = /user2/nakamura/mlspgs/lib/NAG.Sun
l1lib = /user2/nakamura/mlspgs/l1/NAG.Sun

objlib = ${mlslib}/Hdf.o ${mlslib}/MLSCommon.o ${mlslib}/MLSMessageModule.o \
         ${mlslib}/SDPToolkit.o

objl1 = ${l1lib}/Orbit.o ${l1lib}/OutputL1B.o ${l1lib}/Scan.o ${l1lib}/TkL1B.o

objs   = sids_l1boa.o L1boa.o Read.o Sd.o Time.o

sids_l1boa :	${objs} ${objlib} ${objl1}	\
		/user2/nakamura/L1BOA/EPH/PGS_TD_ADEOSIItoTAI.o  \
		${PGSLIB}/libPGSTK.a	\
		${hdflib}/libmfhdf.a	\
		${hdflib}/libdf.a	\
		${hdflib}/libz.a	\
		${hdflib}/libjpeg.a	\
		${mlslib}/libmls.a
	f95 -o  sids_l1boa	\
		${objs}	${objlib} ${objl1}	\
		/user2/nakamura/L1BOA/EPH/PGS_TD_ADEOSIItoTAI.o	\
		${PGSLIB}/libPGSTK.a	\
		${hdflib}/libmfhdf.a	\
		${hdflib}/libdf.a	\
		${hdflib}/libz.a	\
		${hdflib}/libjpeg.a	\
		${mlslib}/libmls.a	\
	-L/usr/lib -lnsl

sids_l1boa.o : sids_l1boa.f90
	f95 -c -I${mlslib} -I${l1lib} sids_l1boa.f90

sids_l1boa.o : ${mlslib}/Hdf.o ${mlslib}/MLSCommon.o \
               ${mlslib}/MLSMessageModule.o ${mlslib}/SDPToolkit.o \
               ${l1lib}/Orbit.o ${l1lib}/OutputL1B.o L1boa.o Read.o Sd.o \
               Time.o

Read.o : ${mlslib}/MLSCommon.o ${mlslib}/SDPToolkit.o Read.f90
	f95 -c -I${mlslib} Read.f90

Sd.o : ${mlslib}/Hdf.o ${mlslib}/MLSMessageModule.o ${l1lib}/OutputL1B.o Sd.f90
	f95 -c -I${mlslib} -I${l1lib} -dusty Sd.f90

Time.o : ${mlslib}/MLSCommon.o ${mlslib}/MLSMessageModule.o ${l1lib}/Scan.o \
         ${l1lib}/TkL1B.o Time.f90
	f95 -c -I${mlslib} -I${l1lib} Time.f90

L1boa.o : ${mlslib}/MLSCommon.o ${mlslib}/MLSMessageModule.o ${l1lib}/Orbit.o \
          ${l1lib}/OutputL1B.o ${l1lib}/Scan.o ${l1lib}/TkL1B.o Sd.o L1boa.f90
	f95 -c -I${mlslib} -I${l1lib} L1boa.f90
