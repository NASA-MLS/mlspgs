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
         $M/tree_types.o        $M/units.o

# These are used to make tar, which is from the parent of $(MLS):
SRCLIB = lib/declaration_table.f90 lib/error_handler.f90     lib/expr_m.f90 \
         lib/hash.f90              lib/intrinsic.f90         lib/io_stuff.f90 \
         lib/lexer_core.f90        lib/lexer_m.f90           lib/machine.f90 \
         lib/MLSCF.f90             lib/MLSMessageModule.f90  lib/output_m.f90 \
         lib/parser.f90            lib/string_table.f90      lib/symbol_table.f90 \
         lib/symbol_types.f90      lib/table_dumper.f90      lib/table_generator.f90 \
         lib/toggles_core.f90      lib/trace_m.f90           lib/tree.f90 \
         lib/tree_types.f90        lib/units.f90

# Files that must be compiled here from sources in the application:
 
OBJ_L2 = init_tables_module.o

# These are used to make tar, which is from the parent of $(MLS):
SRC_L2 = l2/init_tables_module.f90

# Files that must be compiled here from sources in lib, because
# they depend on init_tables_module:

OBJ_SRC_LIB = getCF_m.o tree_checker.o

# These are used to make tar, which is from the parent of $(MLS):
SRC_LIB = lib/getCF_m.f90 lib/tree_checker.f90

all: test

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

$M/libmls.a:
	cd $M; make libmls.a

# =====     For test harness     =========================================

OBJS = test.o $(OBJLIB) $(OBJ_L2) $(OBJ_SRC_LIB)

test.o: $S/test.f90 getCF_m.o $M/machine.o $M/MLSCF.o $M/output_m.o \
        $M/toggles_core.o
	$(F90) -c $(FOPTS) $I$S $I$M $S/test.f90

test: $(OBJS) $M/libmls.a
	$(F90) $(LDOPTS) -o test $(OBJS) $M/SDPToolkit$(SDP).o \
              -L$(PGSTK) -lPGSTK \

out.test: test $S/in.test
	./test -a $S/in.test out.test

# =====     Additional Targets     =======================================

clean:
	rm -f *.o *.mod test out.test Log[RSU]* *.tgz

tar:
	rm -f $S/parser.tgz
	cd $(MLS)/.. ; \
	tar hcvzf lib/test_cf_parser/parser.tgz $(SRCLIB) $(SRC_L2) \
          $(SRC_LIB) lib/test_cf_parser/*

# $Log$
