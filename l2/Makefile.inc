# Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
# U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

# $Id$

# Variables set before this file is included:
# E = Beginning of path to HDFEOS library
# F90 = Name of Fortran compiler
# FOPTS = Options for Fortran compiler, other than -c
# I = option to tell compiler path to includes and .mod files (usually -I)
# M = path to MLS library's .mod files and libmls.a
# MLS = path to MLS library's source files, to get tree_checker.f90
# PLATFORM = The "platform" part of the path to HDF and HDFEOS libraries,
#            e.g. linux, sgi or sun5.
# S = path to source files
# TOOLKIT = Beginning of path to HDFEOS's GCT library, and to HDF's libraries,
#           e.g. /software/SRC/TOOLKIT-5.2.6v1.00.NAG

# Any of these variables can be set by environment variables, or by
# command-line arguments of the "make" command, overriding the values
# set in the Makefile that includes this file.

L2OBJS = Construct.o ConstructQuantityTemplates.o \
       ConstructVectorTemplates.o dumper.o Fill.o get_3d_view.o \
       global_settings.o HGrid.o init_tables_module.o Join.o L2AUXData.o \
       MLSL2.o obtainclimatology.o obtaindao.o obtain_mlscf.o \
       obtainncep.o open_init.o OutputAndClose.o QuantityTemplates.o \
       ScanDivide.o tree_checker.o tree_walker.o VGrid.o VectorsModule.o

all: mlsl2

Construct.o: $S/Construct.f90 ConstructQuantityTemplates.o \
          ConstructVectorTemplates.o dumper.o HGrid.o init_tables_module.o \
          $M/MLSMessageModule.o $M/MLSSignalNomenclature.o QuantityTemplates.o \
          $M/toggles_core.o $M/trace_m.o $M/tree.o $M/tree_types.o \
          VectorsModule.o VGrid.o
	$(F90) -c $(FOPTS) $I$M $I$S $I$M $S/Construct.f90

ConstructQuantityTemplates.o: $S/ConstructQuantityTemplates.f90 \
          $M/declaration_table.o HGrid.o init_tables_module.o $M/lexer_core.o \
          $M/MLSMessageModule.o $M/MLSSignalNomenclature.o $M/output_m.o \
          QuantityTemplates.o $M/string_table.o $M/toggles_core.o \
          $M/trace_m.o $M/tree.o VGrid.o
	$(F90) -c $(FOPTS) $I$M $I$S $I$M $S/ConstructQuantityTemplates.f90

ConstructVectorTemplates.o: $S/ConstructVectorTemplates.f90 \
          $M/Allocate_Deallocate.o init_tables_module.o $M/MLSMessageModule.o \
          QuantityTemplates.o $M/toggles_core.o $M/trace_m.o $M/tree.o \
          VectorsModule.o
	$(F90) -c $(FOPTS) $I$M $I$S $I$M $S/ConstructVectorTemplates.f90

dumper.o: $S/dumper.f90 $M/dump_0.o HGrid.o init_tables_module.o \
          $M/output_m.o QuantityTemplates.o $M/string_table.o VGrid.o
	$(F90) -c $(FOPTS) $I$M $I$S $I$M $S/dumper.f90

Fill.o: $S/Fill.f90 dumper.o init_tables_module.o $M/lexer_core.o \
          $M/output_m.o $M/toggles_core.o $M/trace_m.o $M/tree.o \
          $M/tree_types.o VectorsModule.o
	$(F90) -c $(FOPTS) $I$M $I$S $I$M $S/Fill.f90

get_3d_view.o: $S/get_3d_view.f90
	$(F90) -c $(FOPTS) $I$M $S/get_3d_view.f90

global_settings.o: $S/global_settings.f90 init_tables_module.o \
          $M/toggles_core.o $M/trace_m.o $M/tree.o
	$(F90) -c $(FOPTS) $I$M $I$S $I$M $S/global_settings.f90

HGrid.o: $S/HGrid.f90 $S/addItemToDatabase.f9h $M/expr_m.o \
          init_tables_module.o $M/lexer_core.o $M/output_m.o \
          $M/string_table.o $M/toggles_core.o $M/trace_m.o $M/tree.o
	$(F90) -c $(FOPTS) $I$M $I$S $I$M $S/HGrid.f90

init_tables_module.o: $S/init_tables_module.f90 $M/intrinsic.o \
          $M/symbol_table.o $M/symbol_types.o $M/tree.o $M/tree_types.o
	$(F90) -c $(FOPTS) $I$M $I$S $I$M $S/init_tables_module.f90

Join.o: $S/Join.f90 init_tables_module.o L2AUXData.o $M/L2GPData.o \
          $M/lexer_core.o $M/MLSMessageModule.o $M/output_m.o \
          QuantityTemplates.o $M/tree.o $M/tree_types.o VectorsModule.o
	$(F90) -c $(FOPTS) $I$M $I$S $I$M $S/Join.f90

L2AUXData.o: $S/L2AUXData.f90 $S/addItemToDatabase.f9h \
          $M/Allocate_Deallocate.o $M/MLSMessageModule.o
	$(F90) -c $(FOPTS) $I$S $I$M $S/L2AUXData.f90

MatrixModule.o: $S/MatrixModule.f90 $S/addItemToDatabase.f9h \
          $M/Allocate_Deallocate.o dump_0.o $M/MLSMessageModule.o \
          $M/output_m.o $M/string_table.o
	$(F90) -c $(FOPTS) $I$M $I$S $I$M $S/MatrixModule.f90

MLSL2.o: $S/MLSL2.f90 $M/declaration_table.o init_tables_module.o \
          $M/lexer_core.o $M/machine.o open_init.o $M/parser.o \
          $M/string_table.o $M/toggles_core.o $M/tree.o $M/ tree_checker.o \
          tree_walker.o $M/units.o
	$(F90) -c $(FOPTS) $I$M $I$S $I$M $S/MLSL2.f90

obtainclimatology.o: $M/Allocate_Deallocate.o $S/obtainclimatology.f90 \
          $M/MLSMessageModule.o
	$(F90) -c $(FOPTS) $I$M $I$S $I$M $S/obtainclimatology.f90

obtaindao.o: $S/obtaindao.f90 $M/MLSMessageModule.o
	$(F90) -c $(FOPTS) $I$M $I$S $I$M $S/obtaindao.f90

obtain_mlscf.o: $S/obtain_mlscf.f90 $M/MLSMessageModule.o
	$(F90) -c $(FOPTS) $I$M $I$S $I$M $S/obtain_mlscf.f90

obtainncep.o: $S/obtainncep.f90 $M/MLSMessageModule.o
	$(F90) -c $(FOPTS) $I$M $I$S $I$M $S/obtainncep.f90

open_init.o: $S/open_init.f90 $M/MLSMessageModule.o $M/MLSSignalNomenclature.o \
          $M/string_table.o
	$(F90) -c $(FOPTS) $I$M $I$S $I$M $S/open_init.f90

OutputAndClose.o: $S/OutputAndClose.f90 init_tables_module.o \
          L2AUXData.o $M/L2GPData.o $M/lexer_core.o $M/MLSMessageModule.o \
          $M/output_m.o $M/string_table.o $M/tree.o \
          $M/tree_types.o
	$(F90) -c $(FOPTS) $I$M $I$S $I$M $S/OutputAndClose.f90

QuantityTemplates.o: $S/QuantityTemplates.f90 $S/addItemToDatabase.f9h \
          $M/Allocate_Deallocate.o $M/MLSMessageModule.o \
          $M/MLSSignalNomenclature.o
	$(F90) -c $(FOPTS) $I$M $I$S $I$M $S/QuantityTemplates.f90

ScanDivide.o: $S/ScanDivide.f90 dumper.o $M/expr_m.o init_tables_module.o \
          $M/MLSMessageModule.o $M/string_table.o $M/tree.o
	$(F90) -c $(FOPTS) $I$M $I$S $I$M $S/ScanDivide.f90

tree_checker.o: $(MLS)/tree_checker.f90 $M/declaration_table.o \
          init_tables_module.o $M/lexer_core.o $M/output_m.o \
          $M/string_table.o $M/toggles_core.o $M/trace_m.o $M/tree.o \
          $M/tree_types.o
	$(F90) -c $(FOPTS) $I$M $I$S $I$M $(MLS)/tree_checker.f90

tree_walker.o: $S/tree_walker.f90 Construct.o $M/declaration_table.o dumper.o \
          global_settings.o Fill.o init_tables_module.o Join.o L2AUXData.o \
          $M/L2GPData.o $M/lexer_core.o obtainclimatology.o \
          obtaindao.o obtainncep.o \
          open_init.o $M/output_m.o OutputAndClose.o QuantityTemplates.o \
          ScanDivide.o \
          $M/string_table.o $M/toggles_core.o $M/trace_m.o $M/tree.o \
          $M/tree_types.o VectorsModule.o
	$(F90) -c $(FOPTS) $I$M $I$S $I$M $S/tree_walker.f90

VectorsModule.o: $S/VectorsModule.f90 $S/addItemToDatabase.f9h \
          $M/Allocate_Deallocate.o dumper.o get_3d_view.o \
          $M/MLSMessageModule.o $M/MLSSignalNomenclature.o $M/output_m.o \
          QuantityTemplates.o
	$(F90) -c $(FOPTS) $I$M $I$S $I$M $S/VectorsModule.f90

VGrid.o: $S/VGrid.f90 $S/addItemToDatabase.f9h $M/expr_m.o \
          init_tables_module.o $M/lexer_core.o $M/output_m.o $M/string_table.o \
          $M/toggles_core.o $M/trace_m.o $M/tree.o
	$(F90) -c $(FOPTS) $I$M $I$S $I$M $S/VGrid.f90

# ========================================================================
# Temporary until the toolkit gets installed:
pgs_smf_generatestatusreport.o: $S/pgs_smf_generatestatusreport.f
	$(F90) -c $(FOPTS) $I$M $I$S $I$M $S/pgs_smf_generatestatusreport.f
PG =
#PG = pgs_smf_generatestatusreport.o
#PG = pgs_smf_generatestatusreport.o \
#     $(HOME)/mlspgs/SDPToolkit/src/IO/GEN/PGS_IO_Gen_OpenF.o \
#     $(HOME)/mlspgs/SDPToolkit/src/IO/GEN/PGS_IO_Gen_CloseF.o

GDapi.o: $(TOOLKIT)/hdfeos/src/GDapi.c
	cc -c -DNAGf90Fortran -I$(HOME)/mlspgs/SDPToolkit/include/ \
           -I$(TOOLKIT)/hdf/linux/HDF4.1r3/include/ \
           -I$(TOOLKIT)/hdf/linux/HDF4.1r3/include \
           -I$(TOOLKIT)/hdfeos/include/ $(TOOLKIT)/hdfeos/src/GDapi.c
#GD=GDapi.o
GD=
# ========================================================================

mlsl2: $(L2OBJS) $(GD) $(PG) $M/SDPToolkit$(SDP).o $M/libmls.a
	$(F90) -o mlsl2 $(L2OBJS) \
              $M/SDPToolkit$(SDP).o $(PG) $(GD) \
              -L$M -lmls \
              -L$(PGSTK) -lPGSTK \
              -L$(HDFEOS)/hdfeos/lib/$(PLATFORM) -lhdfeos \
              -L$(TOOLKIT)/hdfeos/lib/$(PLATFORM) -lGctp \
              -L$(TOOLKIT)/hdf/$(PLATFORM)/HDF4.1r3/lib \
                -lmfhdf -ldf -ljpeg -lz

#mlsl2: $(L2OBJS) $(GD) $(PG) $M/SDPToolkit$(SDP).o $M/libmls.a
#	$(F90) -o mlsl2 $(L2OBJS) \
#             $M/SDPToolkit$(SDP).o $(PG) $(GD) \
#             -L$M -lmls \
#             -L$(TOOLKIT)/lib/$(PLATFORM) -lPGSTK \
#             -L$(TOOLKIT)/hdfeos/lib/$(PLATFORM) -lhdfeos
#             -L$(TOOLKIT)/hdfeos/lib/$(PLATFORM) -lGctp \
#             -L$(TOOLKIT)/hdf/$(PLATFORM)/HDF4.1r3/lib \
#               -lmfhdf -ldf -ljpeg -lz

# ========================================================================

clean:
	rm -f *.o *.mod mlsl2

tar:
	rm -f parser.tgz
	tar hcvzf parser.tgz * --exclude \*.o --exclude \*.mod --exclude test

# ========================================================================

# $Log$
# Revision 2.6  2000/10/12 00:38:52  vsnyder
# Account for lib/cf_parser moving to lib
#
# Revision 2.5  2000/10/10 20:33:13  vsnyder
# Repairs to get rid of WVS pathnames
#
# Revision 2.4  2000/10/10 00:24:43  vsnyder
# Account for moving L2GPData.f90 and machine.f90 to lib
#
# Revision 2.3  2000/10/05 01:38:31  vsnyder
# Added CVS comments
#
