# $Id$

# Notice the conspicuous absence of init_tables_module.o in the following
# list.  It is used by the test harness and in other applications of the
# library, but is not part of the library.

OBJLIB = declaration_table.o  error_handler.o expr_m.o        hash.o \
         io_stuff.o           lexer_core.o    lexer_m.o       machine.o \
         output_m.o           parser.o        string_table.o  symbol_table.o \
         symbol_types.o       toggles_core.o  trace_m.o       tree.o \
         tree_checker.o       tree_types.o    units.o

all: test

declaration_table.o: $S/declaration_table.f90 init_tables_module.o \
                     machine.o output_m.o string_table.o toggles_core.o
	$(FC) -c $(FOPTS) $I$S $S/declaration_table.f90

error_handler.o: $S/error_handler.f90 machine.o output_m.o
	$(FC) -c $(FOPTS) $I$S $S/error_handler.f90

expr_m.o: $S/expr_m.f90 declaration_table.o init_tables_module.o \
          string_table.o toggles_core.o trace_m.o tree.o tree_types.o
	$(FC) -c $(FOPTS) $I$S $S/expr_m.f90

hash.o: $S/hash.f90
	$(FC) -c $(FOPTS) $I$S $S/hash.f90

init_tables_module.o: $S/init_tables_module.f90 \
                      symbol_table.o symbol_types.o tree.o tree_types.o
	$(FC) -c $(FOPTS) $I$S $S/init_tables_module.f90

io_stuff.o: $S/io_stuff.f90
	$(FC) -c $(FOPTS) $I$S $S/io_stuff.f90

lexer_core.o: $S/lexer_core.f90 string_table.o symbol_table.o toggles_core.o
	$(FC) -c $(FOPTS) $I$S $S/lexer_core.f90

lexer_m.o: $S/lexer_m.f90 error_handler.o lexer_core.o output_m.o \
           string_table.o symbol_table.o symbol_types.o toggles_core.o
	$(FC) -c $(FOPTS) $I$S $S/lexer_m.f90

MLSCF.o: $S/MLSCF.f90
	$(FC) -c $(FOPTS) $I$S $S/MLSCF.f90

output_m.o: $S/output_m.f90
	$(FC) -c $(FOPTS) $I$S $S/output_m.f90

parser.o: $S/parser.f90 lexer_core.o lexer_m.o symbol_types.o tree.o \
          tree_types.o
	$(FC) -c $(FOPTS) $I$S $S/parser.f90

string_table.o: $S/string_table.f90 hash.o io_stuff.o machine.o output_m.o
	$(FC) -c $(FOPTS) $I$S $S/string_table.f90

symbol_table.o: $S/symbol_table.f90 machine.o output_m.o string_table.o \
                symbol_types.o toggles_core.o
	$(FC) -c $(FOPTS) $I$S $S/symbol_table.f90

symbol_types.o: $S/symbol_types.f90 string_table.o
	$(FC) -c $(FOPTS) $I$S $S/symbol_types.f90

trace_m.o: $S/trace_m.f90 lexer_core.o output_m.o tree.o
	$(FC) -c $(FOPTS) $I$S $S/trace_m.f90

tree_checker.o: $S/tree_checker.f90 declaration_table.o init_tables_module.o \
                lexer_core.o output_m.o string_table.o toggles_core.o \
                trace_m.o tree.o tree_types.o
	$(FC) -c $(FOPTS) $I$S $S/tree_checker.f90

toggles_core.o: $S/toggles_core.f90
	$(FC) -c $(FOPTS) $I$S $S/toggles_core.f90

tree.o: $S/tree.f90 error_handler.o machine.o output_m.o string_table.o \
        symbol_table.o symbol_types.o toggles_core.o tree_types.o
	$(FC) -c $(FOPTS) $I$S $S/tree.f90

tree_types.o: $S/tree_types.f90 string_table.o symbol_types.o
	$(FC) -c $(FOPTS) $I$S $S/tree_types.f90

units.o: $S/units.f90 declaration_table.o init_tables_module.o tree.o
	$(FC) -c $(FOPTS) $I$S $S/units.f90

# =====     The Real Targets     =========================================

libCF.a: $(OBJLIB)
	ar ru libCF.a $(OBJLIB)

clean:
	rm -f *.o *.a *.mod *.g90 $S/*.g90 test out.test

tar:
	rm -f parser.tgz
	tar hcvzf parser.tgz * --exclude \*.o --exclude \*.mod \
                               --exclude \*.g90 --exclude CVS \
                               --exclude \*.a --exclude test --exclude out.test

# =====     For test harness     =========================================

OBJS = init_tables_module.o MLSCF.o table_dumper.o \
       table_generator.o test.o

table_dumper.o: $S/table_dumper.f90 declaration_table.o init_tables_module.o \
                MLSCF.o output_m.o string_table.o
	$(FC) -c $(FOPTS) $I$S $S/table_dumper.f90

table_generator.o: $S/table_generator.f90 declaration_table.o expr_m.o \
                   init_tables_module.o lexer_core.o MLSCF.o output_m.o \
                   string_table.o toggles_core.o trace_m.o tree.o tree_types.o
	$(FC) -c $(FOPTS) $I$S $S/table_generator.f90

test.o: $S/test.f90 declaration_table.o init_tables_module.o lexer_core.o \
        machine.o MLSCF.o parser.o string_table.o tree_checker.o \
        table_generator.o toggles_core.o tree.o units.o \
        table_dumper.o
	$(FC) -c $(FOPTS) $I$S $S/test.f90

test: $(OBJS) libCF.a
	$(FC) $(LDOPTS) -o test $(OBJS) -L. -lCF

out.test: test $S/in.test
	./test -a <$S/in.test >out.test

# $Log$
# Revision 1.4  2000/07/28 18:30:18  vsnyder
# Added a comment to make it clear that init_tables_module is not part
# of the library.
#
# Revision 1.3  2000/07/12 17:51:27  vsnyder
# Added some more stuff to 'clean:'
#
# Revision 1.2  2000/07/06 18:03:17  vsnyder
# Clean up -- remove unused options and library searches
#
# Revision 1.1  2000/07/06 01:43:12  vsnyder
# Initial check-in
#
