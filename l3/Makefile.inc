# Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
# U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

# $Id$

# Variables set before this file is included:
# A = Path to directory for application's source, to get init_tables_module.f90
# F90 = Name of Fortran compiler
# FOPTS = Options for Fortran compiler, other than -c
# I = option to tell compiler path to includes and .mod files (usually -I)
# M = path to MLS library's .mod files and libmls.a
# MLS = path to MLS library's source files
# PLATFORM = The "platform" part of the path to HDF and HDFEOS libraries,
#            e.g. linux, sgi or sun5.
# S = path to source files
# SDP = "Substitute" or empty -- used as: SDPToolkit$(SDP).o
# TOOLKIT = Beginning of path to HDFEOS's GCT library, and to HDF's libraries,
#           e.g. /software/SRC/TOOLKIT-5.2.6v1.00.NAG

# Object files needed from lib:

OBJLIB = $M/declaration_table.o $M/error_handler.o     $M/expr_m.o \
         $M/hash.o              $M/intrinsic.o         $M/io_stuff.o \
         $M/lexer_core.o        $M/lexer_m.o           $M/machine.o \
         $M/MLSCF.o             $M/MLSMessageModule.o  $M/output_m.o \
         $M/parser.o            $M/string_table.o      $M/symbol_table.o \
         $M/symbol_types.o      $M/table_dumper.o      $M/table_generator.o \
         $M/toggles_core.o      $M/trace_m.o           $M/tree.o \
         $M/tree_types.o        $M/units.o             $M/MLSStrings.o \
         $M/MLSPCF.o            $M/dates_module.o

# These are used to make tar, which is from the parent of $(MLS):
SRCLIB = lib/declaration_table.f90 lib/error_handler.f90     lib/expr_m.f90 \
         lib/hash.f90              lib/intrinsic.f90         lib/io_stuff.f90 \
         lib/lexer_core.f90        lib/lexer_m.f90           lib/machine.f90 \
         lib/MLSCF.f90             lib/MLSMessageModule.f90  lib/output_m.f90 \
         lib/parser.f90            lib/string_table.f90      lib/symbol_table.f90 \
         lib/symbol_types.f90      lib/table_dumper.f90      lib/table_generator.f90 \
         lib/toggles_core.f90      lib/trace_m.f90           lib/tree.f90 \
         lib/tree_types.f90        lib/units.f90             lib/MLSStrings.f90 \
         lib/MLSPCF.f90            lib/dates_module.f90

# Files that must be compiled here from sources in the application:
 
OBJ_L3 = init_tables_module.o L3CF.o L3DMData.o OpenInit.o OutputClose.o \
         PCFModule.o 

# These are used to make tar, which is from the parent of $(MLS):
SRC_L3 = l3/init_tables_module.f90 l3/L3CF.f90 l3/L3DMData.f90 \
         l3/OpenInit.f90 l3/OutputClose.f90 l3/PCFModule.f90


# Files that must be compiled here from sources in lib, because
# they depend on init_tables_module:

OBJ_SRC_LIB = getCF_m.o tree_checker.o

# These are used to make tar, which is from the parent of $(MLS):
SRC_LIB = lib/getCF_m.f90 lib/tree_checker.f90

all: MLSL3

init_tables_module.o: $A/init_tables_module.f90 $M/intrinsic.o \
                      $M/symbol_table.o $M/symbol_types.o $M/tree.o \
                      $M/tree_types.o
	$(F90) -c $(FOPTS) $I$M $A/init_tables_module.f90

tree_checker.o: $(MLS)/tree_checker.f90 $M/declaration_table.o \
                init_tables_module.o $M/lexer_core.o $M/output_m.o \
                $M/string_table.o $M/toggles_core.o $M/trace_m.o $M/tree.o \
                $M/tree_types.o
	$(F90) -c $(FOPTS) $I$M $(MLS)/tree_checker.f90

getCF_m.o: $(MLS)/getCF_m.f90 $M/declaration_table.o init_tables_module.o \
           $M/lexer_core.o $M/MLSCF.o $M/output_m.o $M/parser.o \
           $M/string_table.o $M/table_dumper.o $M/table_generator.o \
           $M/toggles_core.o tree_checker.o $M/tree.o $M/units.o
	$(F90) -c $(FOPTS) $I$M $(MLS)/getCF_m.f90

L3DMData.o : $M/MLSMessageModule.o $M/MLSCommon.o $M/L2GPData.o $M/Hdf.o \
             $M/MLSPCF.o
	$(F90) -c $(FOPTS) $(DUSTY) $I$M $A/L3DMData.f90

PCFModule.o : $M/MLSCommon.o $M/MLSMessageModule.o $M/MLSPCF.o
	$(F90) -c $(FOPTS) $I$M $A/PCFModule.f90

L3CF.o : L3DMData.o $M/MLSCommon.o $M/MLSCF.o $M/MLSStrings.o \
         $M/MLSMessageModule.o PCFModule.o $M/MLSPCF.o
	$(F90) -c $(FOPTS) $I$M $A/L3CF.f90

OpenInit.o : $A/OpenInit.f90 $M/MLSPCF.o $M/MLSMessageModule.o L3DMData.o \
             $M/dates_module.o $M/MLSCF.o L3CF.o getCF_m.o
	$(F90) -c $(FOPTS) $I$M $I$A $A/OpenInit.f90

OutputClose.o : $M/MLSCommon.o $M/MLSMessageModule.o L3DMData.o L3CF.o \
                $M/MLSStrings.o $M/MLSPCF.o OpenInit.o $M/MLSCF.o PCFModule.o
	$(F90) -c $(FOPTS) $I$M $I$A $A/OutputClose.f90

$M/libmls.a:
	cd $M; make libmls.a

# ========================================================================

OBJS = MLSL3.o $(OBJLIB) $(OBJ_L3) $(OBJ_SRC_LIB)

MLSL3.o: $S/MLSL3.f90 OpenInit.o $M/MLSCF.o L3CF.o $M/MLSMessageModule.o \
        L3DMData.o OutputClose.o
	$(F90) -c $(FOPTS) $I$S $I$M $I$A $S/MLSL3.f90

MLSL3: $(OBJS) $M/libmls.a
	$(F90) $(LDOPTS) -o MLSL3 $(OBJS) $M/SDPToolkit$(SDP).o \
              -L$(PGSTK) -lPGSTK \
              -L$(HDFEOS)/hdfeos/lib/$(PLATFORM) -lhdfeos \
              -L$(TOOLKIT)/hdfeos/lib/$(PLATFORM) -lGctp \
              -L$(TOOLKIT)/hdf/$(PLATFORM)/HDF4.1r3/lib \
                -lmfhdf -ldf -ljpeg -lz -lnsl

# =====     Additional Targets     =======================================

clean:
	rm -f *.o *.mod MLSL3 Log[RSU]* *.tgz

tar:
	rm -f $S/parser.tgz
	cd $(MLS)/.. ; \
	tar hcvzf lib/test_cf_parser/parser.tgz $(SRCLIB) $(SRC_L3) \
          $(SRC_LIB) lib/test_cf_parser/*

# $Log$
# Revision 1.1  2000/10/11 20:54:39  vsnyder
# Initial addition
#
