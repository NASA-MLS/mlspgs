#	Makefile.h
# An include file for all mainline Makefiles
# We assemble here all the dependencies, including hierarchical ones,
# to avoid having to repeat them when, e.g., building in l2 what
# we already had to say about when and how to build ib lib
#
# Copyright 2012, by the California Institute of Technology. ALL
# RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
# commercial use must be negotiated with the Office of Technology Transfer
# at the California Institute of Technology.

# This software may be subject to U.S. export control laws. By accepting this
# software, the user agrees to comply with all applicable U.S. export laws and
# regulations. User has the responsibility to obtain export licenses, or other
# export authority as may be required before exporting such information to
# foreign countries or providing access to foreign persons.

# "$Id$"

# What is your preferred language?
# Choose one of Eng, Es, Fr, De, It
PREFERREDLANG=Eng

# Oh, man do we really need this? Only if lit_names.txt has changed
EVERY_MLSCONFG = $(shell ${REECHO} -d -dirn $(CONFDIR)/lib/machines -excl CVS -excl NAG.nogc)

# Two useful utilities
$(INSTALLDIR)/init_gen: $(UTILDIR)/init_gen.f90
	$(UTILDIR)/build_f90_in_misc.sh -p init_gen \
	-d $(INSTALLDIR) -t $(TESTSDIR) -M $(MAKE) \
	-C $(MLSCFILE) \
   $(UTILDIR)/init_gen.f90
$(INSTALLDIR)/lr: $(UTILDIR)/lr/*.f90
	$(UTILDIR)/build_f90_in_misc.sh -d $(INSTALLDIR) -t $(TESTSDIR) \
	-c $(MLSCONFG) -p lr -M $(MAKE) -O short_name=lr_custom -m lib \
	-C $(MLSCFILE) \
        $(UTILDIR)/lr/*.[Ff]9[0h] $(UTILDIR)/lr/*.grm \
        $(CONFDIR)/lib/symbol_table.f90 $(CONFDIR)/lib/tree.f90

# Used only by wvs-149.tex
$(INSTALLDIR)/wvs-149: $(DOCDIR)/*.f90 $(DOCDIR)/Math77/*.f
	$(UTILDIR)/build_f90_in_misc.sh -d $(INSTALLDIR) -t $(TESTSDIR) \
	-c $(MLSCONFG) -p wvs-149 -M $(MAKE)  \
	-C $(MLSCFILE) \
        $(DOCDIR)/wvs-149.f90 $(DOCDIR)/Math77/*.f

ifneq ($(short_name),doc)
ifndef CASCADE
#	l1

CalibWeightsFlags.o: calibweightsflags.mod
	$(UTILDIR)/mark_as_uptodate.sh -M $(MAKE) -t \
          -T CalibWeightsFlags.o calibweightsflags.mod
calibweightsflags.mod:
	$(UTILDIR)/newAifBdiff.sh -a calibweightsflags.mod $(FC) -c $(BUGGY) \
          $(INC_PATHS) $(S)/CalibWeightsFlags.f90 $(FAFTER)
EngTbls.o: engtbls.mod
	$(UTILDIR)/mark_as_uptodate.sh -M $(MAKE) -t \
          -T EngTbls.o engtbls.mod 
engtbls.mod:
	$(UTILDIR)/newAifBdiff.sh -a engtbls.mod $(FC) -c $(DUSTIER) \
          $(INC_PATHS) $(S)/EngTbls.f90 $(FAFTER)
L1LogUtils.o: l1logutils.mod
	$(UTILDIR)/mark_as_uptodate.sh -M $(MAKE) -t \
          -T L1LogUtils.o l1logutils.mod
l1logutils.mod:
	$(UTILDIR)/newAifBdiff.sh -a l1logutils.mod $(FC) -c $(BUGGY) \
          $(INC_PATHS) $(S)/L1LogUtils.f90 $(FAFTER)
OutputL1B.o: outputl1b.mod
	$(UTILDIR)/mark_as_uptodate.sh -M $(MAKE) -t \
          -T OutputL1B.o outputl1b.mod 
outputl1b.mod:
	$(UTILDIR)/newAifBdiff.sh -a outputl1b.mod $(FC) -c $(DUSTY) $(INC_PATHS) $(S)/OutputL1B.f90 $(FAFTER)
OutputL1B_HDF4.o: outputl1b_hdf4.mod
	$(UTILDIR)/mark_as_uptodate.sh -M $(MAKE) -t \
          -T OutputL1B_HDF4.o outputl1b_hdf4.mod 
outputl1b_hdf4.mod:
	$(UTILDIR)/newAifBdiff.sh -a outputl1b_hdf4.mod $(FC) -c $(DUSTY) $(INC_PATHS) $(S)/OutputL1B_HDF4.f90 $(FAFTER)
parser_tables_l2cf.mod: Parser_Tables_L2CF.f90 Parser_L2CF.f9h
	$(UTILDIR)/newAifBdiff.sh -a parser_tables_l2cf.mod $(FC) -c $(DUSTY) $(INC_PATHS) $(srclib)/Parser_Tables_L2CF.f90 $(FAFTER)

Parser_L2CF.f9h $(srclib)/Parser_L2CF.f9h: $(UTILDIR)/lr/l2cf.grm $(INSTALLDIR)/lr
	$(INSTALLDIR)/lr \
          $(UTILDIR)/lr/l2cf.grm \
          $(srclib)/Parser_L2CF.f9h $(UTILDIR)/lr/l2cf.lls $(LRAFTER)
Parser_Tables_L2CF.o: Parser_Tables_L2CF.f90 Parser_L2CF.f9h
	$(UTILDIR)/mark_as_uptodate.sh -M $(MAKE) -t \
          -T Parser_Tables_L2CF.o parser_tables_l2cf.mod 
tree_checker.mod: tree_checker.f90 init_tables_module.mod
	$(UTILDIR)/newAifBdiff.sh -a tree_checker.mod $(FC) -c $(DUSTY) $(INC_PATHS) $(srclib)/tree_checker.f90 $(FAFTER)
tree_checker.o: tree_checker.mod tree_checker.f90
	$(UTILDIR)/mark_as_uptodate.sh -M $(MAKE) -t \
          -T tree_checker.o tree_checker.mod 

# l2
# This keeps switch_usage up-to-date with notes/switches
# (and similarly for options)
MLSL2.o: $(CONFDIR)/notes/switches $(CONFDIR)/notes/options $(S)/MLSL2.f90
	@rm -f MLSL2.f90 part_1 part_2 part_3 ; \
	sed -n '1,/=== (start of automatic usage lines) ===/ p' $(S)/MLSL2.f90 > part_1 ; \
	sed -n '/=== (end of automatic usage lines) ===/,$$ p' $(S)/MLSL2.f90 > part_3 ; \
	sed -n '/mlsl2 =========/,/mlsl2 ========/p' $(CONFDIR)/notes/switches > part_2 ; \
	sed '1 d; $$ d' part_2 > MLSL2.f90 ; \
	mv MLSL2.f90 part_2 ; \
	sed -n -f $(CONFDIR)/notes/switch_usage.sed part_2 > MLSL2.f90 ; \
	mv MLSL2.f90 part_2 ; \
	cat part_1 part_2 part_3 > MLSL2.f90 ; \
	rm -f part_1 part_2 part_3 ; \
	sed -n '1,/=== (start of automatic option lines) ===/ p' MLSL2.f90 > part_1 ; \
	sed -n '/=== (end of automatic option lines) ===/,$$ p' MLSL2.f90 > part_3 ; \
	sed -n '/mlsl2 =========/,/mlsl2 ========/p' $(CONFDIR)/notes/options > part_2 ; \
	sed '1 d; $$ d' part_2 > MLSL2.f90 ; \
	mv MLSL2.f90 part_2 ; \
	sed -n -f $(CONFDIR)/notes/switch_usage.sed part_2 > MLSL2.f90 ; \
	mv MLSL2.f90 part_2 ; \
	cat part_1 part_2 part_3 > MLSL2.f90
	$(FC) -c $(FOPTS) $(INC_PATHS) MLSL2.f90 $(FAFTER)
	rm -f MLSL2.f90 part_1 part_2 part_3

#tree_checker.mod: init_tables_module.mod
#	$(UTILDIR)/newAifBdiff.sh -a tree_checker.mod $(FC) -c $(DUSTY) $(INC_PATHS) $(srclib)/tree_checker.f90 $(FAFTER)

ifeq ($(short_name),l2)
# This uses the utility program init_gen to create two
# up-to-date prereqs for init_tables_module
$(S)/field_parm.f9h $(S)/field_add.f9h: $(S)/field_names.txt $(INSTALLDIR)/init_gen
	sort -f -u $(S)/field_names.txt | $(INSTALLDIR)/init_gen \
          -f"last_Spectroscopy_Field + 1" -l field_last -pF_ \
          -j field_names.txt \
          $(S)/field_parm.f9h $(S)/field_add.f9h field_indices

$(S)/lit_parm.f9h $(S)/lit_add.f9h: $(S)/lit_names.txt $(INSTALLDIR)/init_gen
	sort -f -u $(S)/lit_names.txt | $(INSTALLDIR)/init_gen \
          -f"last_Spectroscopy_Lit + 1" -l last_lit -pl_ \
          -j lit_names.txt \
          $(S)/lit_parm.f9h $(S)/lit_add.f9h lit_indices

init_tables_module.mod init_tables_module.o: $(CONFDIR)/lib/lit_add.f9h \
        $(CONFDIR)/lib/lit_parm.f9h
endif # short_name == l2

#	l3

l2interface.mod:
	$(UTILDIR)/newAifBdiff.sh -a l2interface.mod $(FC) -c $(DUSTY) $(INC_PATHS) $(S)/L2Interface.f90 $(FAFTER)
L2Interface.o: l2interface.mod
	$(UTILDIR)/mark_as_uptodate.sh -M $(MAKE) -t \
   -T L2Interface.o l2interface.mod 
l3dmdata.mod:
	$(UTILDIR)/newAifBdiff.sh -a l3dmdata.mod $(FC) -c $(DUSTY) $(INC_PATHS) $(S)/L3DMData.f90 $(FAFTER)
L3DMData.o: l3dmdata.mod
	$(UTILDIR)/mark_as_uptodate.sh -M $(MAKE) -t \
   -T L3DMData.o l3dmdata.mod 
#tree_checker.mod: init_tables_module.mod
#	$(UTILDIR)/newAifBdiff.sh -a tree_checker.mod $(FC) -c $(DUSTY) $(INC_PATHS) $(srclib)/tree_checker.f90 $(FAFTER)
getcf_m.mod: init_tables_module.mod tree_checker.mod
	$(UTILDIR)/newAifBdiff.sh -a getcf_m.mod $(FC) -c $(DUSTY) $(INC_PATHS) $(srclib)/getCF_m.f90 $(FAFTER)
getCF_m.o: getcf_m.mod
	$(UTILDIR)/mark_as_uptodate.sh -M $(MAKE) -t \
   -T getCF_m.o getcf_m.mod 

#  l3m
l3dzdata.mod:
	$(UTILDIR)/newAifBdiff.sh -a l3dzdata.mod $(FC) -c $(DUSTY) $(INC_PATHS) $(S)/L3DZData.f90 $(FAFTER)
L3DZData.o: l3dzdata.mod
	$(UTILDIR)/mark_as_uptodate.sh -M $(MAKE) -t \
   -T L3DZData.o l3dzdata.mod 

l3mzdata.mod:
	$(UTILDIR)/newAifBdiff.sh -a l3mzdata.mod $(FC) -c $(DUSTY) $(INC_PATHS) $(S)/L3MZData.f90 $(FAFTER)
L3MZData.o: l3mzdata.mod
	$(UTILDIR)/mark_as_uptodate.sh -M $(MAKE) -t \
   -T L3MZData.o l3mzdata.mod 

l3mmdata.mod:
	$(UTILDIR)/newAifBdiff.sh -a l3mmdata.mod $(FC) -c $(DUSTY) $(INC_PATHS) $(S)/L3MMData.f90 $(FAFTER)
L3MMData.o: l3mmdata.mod
	$(UTILDIR)/mark_as_uptodate.sh -M $(MAKE) -t \
   -T L3MMData.o l3mmdata.mod 

# lib

ifeq ($(short_name),lib)
# This uses the utility program init_gen to create two
# up-to-date prereqs for intrinsic
# Due to a strange error in level 2, we also create a sentinel file, 
# new_lit_names.txt
# whose existence will be detected later in l2 to trigger linking mlsl2 twice
# We must copy that sentienel file to each MLSCONFG
# because oherwise the .f9h will be already up-to-date and so no sentinel
# file will becreated if you change one MLSCONFG for another
$(S)/lit_parm.f9h $(S)/lit_add.f9h: $(S)/lit_names.txt $(INSTALLDIR)/init_gen
	sort -f -u $(S)/lit_names.txt | $(INSTALLDIR)/init_gen \
          -l last_intrinsic_lit -pl_ \
          -j lit_names.txt \
          $(S)/lit_parm.f9h $(S)/lit_add.f9h lit_indices
	echo "new lit_names" > new_lit_names.txt
	@for dir in $(EVERY_MLSCONFG); do \
          if [ -d ../$$dir ]; then \
            cp new_lit_names.txt ../$$dir 2>/dev/null; \
          fi; \
        done

$(S)/mol_parm.f9h $(S)/mol_add.f9h: $(TABLES_DIR)/sps_cross_ref_table.txt $(INSTALLDIR)/init_gen
	sed -n '/---\*/,$$ p' $(TABLES_DIR)/sps_cross_ref_table.txt | sed '1 d' | \
          sed '/H2O_R/ d' | awk ' {print $$2}' | sort -f -u| \
          $(INSTALLDIR)/init_gen \
          -f "last_rxx_molecule+1" -l last_molecule -pl_ \
          -j sps_cross_ref_table.txt \
          $(S)/mol_parm.f9h $(S)/mol_add.f9h lit_indices

machine.mod: $(MACH_DIR)/machine.f90
	rm -f machine.f90 ; \
	sed '/Start $(BAD_MATCH)/,/End $(BAD_MATCH)/ d' $(MACH_DIR)/machine.f90 \
     > machine.f90
	$(UTILDIR)/newAifBdiff.sh -a machine.mod $(FC) -c $(FOPTS) $(INC_PATHS) machine.f90 $(FAFTER)
	rm -f machine.f90
machine.o: machine.mod $(MACH_DIR)/machine.f90
	$(UTILDIR)/mark_as_uptodate.sh -M $(MAKE) -t \
   -T machine.o machine.mod 

ieee_arithmetic.o: ieee_arithmetic.mod ieee_arithmetic.f90
	$(UTILDIR)/mark_as_uptodate.sh -M $(MAKE) -t -T ieee_arithmetic.o ieee_arithmetic.mod
ieee_arithmetic.mod: ieee_arithmetic.f90
	$(UTILDIR)/newAifBdiff.sh -a ieee_arithmetic.mod $(FC) -c $(FOPTS) $(INC_PATHS) $(MACH_DIR)/ieee_arithmetic.f90 $(FAFTER)


intrinsic.o: $(S)/lit_parm.f9h $(S)/lit_add.f9h

dates_module.o: $(S)/DayMonthNames.f9h $(S)/LeapSecDates.f9h

$(S)/DayMonthNames.f9h: $(S)/DayMonthNames_${PREFERREDLANG}.f9h
	cp -a $(S)/DayMonthNames_${PREFERREDLANG}.f9h $(S)/DayMonthNames.f9h

NUMDATES := $(shell sed -n '2,$$ p' $(PGSTK)/../../database/common/TD/leapsec.dat | wc -l)

$(S)/LeapSecDates.f9h: $(PGSTK)/../../database/common/TD/leapsec.dat
# The first and last lines are merely Fortran boiler plate.
# The long sequence beginning with 'sed' is the real meat of the sandwich.
#
# Bugs and Limitations
# (1) If the format of the toolkit's leapsec.dat file ever changes
#     it will break spectacularly.
# (2) If the leapsec.dat file is merely checked, its log and
#     modification dates are updated and so make must rebuild LeapSecDates.f9h
#     We could ameliorate this shortcoming if we used the script
#     $(UTILDIR)/newAifBdiff.sh
	echo \
	  'character(len=*), dimension($(NUMDATES)), parameter :: LeapSecDates = (/ &' \
	  > $(S)/LeapSecDates.f9h
	sed -n '2,$$ p' $< | \
	  awk '{print $$1, $$2, $$3}' | sed 's/^/    \"/' | \
	  sed 's/$$/\", \&/' | \
	  sed '$$ s/,//' \
	  >> $(S)/LeapSecDates.f9h
	echo \
	  '    /) ' \
	  >> $(S)/LeapSecDates.f9h

endif # end short_name == lib

NOUNASS := $(shell echo ${DUSTY} | sed 's/\(--chk [aesx,]*\),u\([aesx,]*\)/\1\2/')

DUSTY_NO_OPT := $(shell echo ${DUSTY} | sed 's/-O[0-9]*//')
DONT_OPT := $(shell echo ${FOPTS} | sed 's/-O[0-9]*//; a-O0')

hdf.mod:
	$(UTILDIR)/newAifBdiff.sh -a hdf.mod $(FC) -c $(DUSTY) $(INC_PATHS) $(S)/Hdf.f90 $(FAFTER)
Hdf.o: hdf.mod
	$(UTILDIR)/mark_as_uptodate.sh -M $(MAKE) -t \
   -T Hdf.o hdf.mod 
string_table.mod:
	$(UTILDIR)/newAifBdiff.sh -a string_table.mod $(FC) -c $(BUGGY) $(INC_PATHS) $(S)/string_table.f90 $(FAFTER)
string_table.o: string_table.mod
	$(UTILDIR)/mark_as_uptodate.sh -M $(MAKE) -t \
   -T string_table.o string_table.mod 

# readGriddedUtils module doesn't play well with the lf95 undefined variable
# runtime checker because of endian swapping
readgriddedutils.mod:
	$(UTILDIR)/newAifBdiff.sh -a readgriddedutils.mod $(FC) -c $(NOUNASS) $(INC_PATHS) $(S)/readGriddedUtils.f90 $(FAFTER)
readGriddedUtils.o: readgriddedutils.mod
	$(UTILDIR)/mark_as_uptodate.sh -M $(MAKE) -t \
   -T readGriddedUtils.o readgriddedutils.mod 

else # end ifndef CASCADE
# i.e. CASCADE IS defined

#	l1

CalibWeightsFlags.o: $(S)/CalibWeightsFlags.f90
	$(FC) -c $(BUGGY) $(INC_PATHS) $(S)/CalibWeightsFlags.f90 $(FAFTER)
EngTbls.o: $(S)/EngTbls.f90
	$(FC) -c $(DUSTY) $(INC_PATHS) $(S)/EngTbls.f90 $(FAFTER)
L1LogUtils.o: $(S)/L1LogUtils.f90
	$(FC) -c $(BUGGY) $(INC_PATHS) $(S)/L1LogUtils.f90 $(FAFTER)
OutputL1B.o: $(S)/OutputL1B.f90
	$(FC) -c $(DUSTY) $(INC_PATHS) $(S)/OutputL1B.f90 $(FAFTER)
OutputL1B_HDF4.o: $(S)/OutputL1B_HDF4.f90
	$(FC) -c $(DUSTY) $(INC_PATHS) $(S)/OutputL1B_HDF4.f90 $(FAFTER)

# l2
# This keeps switch_usage up-to-date with notes/switches
# (and similarly for options)
MLSL2.o: $(CONFDIR)/notes/switches $(CONFDIR)/notes/options $(S)/MLSL2.f90
	@rm -f MLSL2.f90 part_1 part_2 part_3 ; \
	sed -n '1,/=== (start of automatic usage lines) ===/ p' $(S)/MLSL2.f90 > part_1 ; \
	sed -n '/=== (end of automatic usage lines) ===/,$$ p' $(S)/MLSL2.f90 > part_3 ; \
	sed -n '/mlsl2 =========/,/mlsl2 ========/p' $(CONFDIR)/notes/switches > part_2 ; \
	sed '1 d; $$ d' part_2 > MLSL2.f90 ; \
	mv MLSL2.f90 part_2 ; \
	sed -n -f $(CONFDIR)/notes/switch_usage.sed part_2 > MLSL2.f90 ; \
	mv MLSL2.f90 part_2 ; \
	cat part_1 part_2 part_3 > MLSL2.f90 ; \
	rm -f part_1 part_2 part_3 ; \
	sed -n '1,/=== (start of automatic option lines) ===/ p' MLSL2.f90 > part_1 ; \
	sed -n '/=== (end of automatic option lines) ===/,$$ p' MLSL2.f90 > part_3 ; \
	sed -n '/mlsl2 =========/,/mlsl2 ========/p' $(CONFDIR)/notes/options > part_2 ; \
	sed '1 d; $$ d' part_2 > MLSL2.f90 ; \
	mv MLSL2.f90 part_2 ; \
	sed -n -f $(CONFDIR)/notes/switch_usage.sed part_2 > MLSL2.f90 ; \
	mv MLSL2.f90 part_2 ; \
	cat part_1 part_2 part_3 > MLSL2.f90
	$(FC) -c $(FOPTS) $(INC_PATHS) MLSL2.f90 $(FAFTER)
	rm -f MLSL2.f90 part_1 part_2 part_3

ifeq ($(short_name),l2)
# This uses the utility program init_gen to create two
# up-to-date prereqs for init_tables_module
$(S)/field_parm.f9h $(S)/field_add.f9h: $(S)/field_names.txt $(INSTALLDIR)/init_gen
	sort -f -u $(S)/field_names.txt | $(INSTALLDIR)/init_gen \
          -f"last_Spectroscopy_Field + 1" -l field_last -pF_ \
          -j field_names.txt \
          $(S)/field_parm.f9h $(S)/field_add.f9h field_indices

$(S)/lit_parm.f9h $(S)/lit_add.f9h: $(S)/lit_names.txt $(INSTALLDIR)/init_gen
	sort -f -u $(S)/lit_names.txt | $(INSTALLDIR)/init_gen \
          -f"last_Spectroscopy_Lit + 1" -l last_lit -pl_ \
          -j lit_names.txt \
          $(S)/lit_parm.f9h $(S)/lit_add.f9h lit_indices

init_tables_module.o: $(CONFDIR)/lib/lit_add.f9h \
        $(CONFDIR)/lib/lit_parm.f9h
endif # end short_name==l2
#	l3

L2Interface.o: $(S)/L2Interface.f90
	$(FC) -c $(DUSTY) $(INC_PATHS) $(S)/L2Interface.f90 $(FAFTER)
L3DMData.o: $(S)/L3DMData.f90
	$(FC) -c $(DUSTY) $(INC_PATHS) $(S)/L3DMData.f90 $(FAFTER)
L3DMDiag.o: $(S)/L3DMDiag.f90
	$(FC) -c $(DUSTY) $(INC_PATHS) $(S)/L3DMDiag.f90 $(FAFTER)

#  l3m
L3DZData.o: $(S)/L3DZData.f90
	$(FC) -c $(DUSTY) $(INC_PATHS) $(S)/L3DZData.f90 $(FAFTER)

L3MZData.o: $(S)/L3MZData.f90
	$(FC) -c $(DUSTY) $(INC_PATHS) $(S)/L3MZData.f90 $(FAFTER)

L3MMData.o: $(S)/L3MMData.f90
	$(FC) -c $(DUSTY) $(INC_PATHS) $(S)/L3MMData.f90 $(FAFTER)

L3DZDiag.o: $(S)/L3DZDiag.f90
	$(FC) -c $(DUSTY) $(INC_PATHS) $(S)/L3DZDiag.f90 $(FAFTER)

L3MZDiag.o: $(S)/L3MZDiag.f90
	$(FC) -c $(DUSTY) $(INC_PATHS) $(S)/L3MZDiag.f90 $(FAFTER)

L3MMDiag.o: $(S)/L3MMDiag.f90
	$(FC) -c $(DUSTY) $(INC_PATHS) $(S)/L3MMDiag.f90 $(FAFTER)


# Custom build commands for some modules

ifeq ($(short_name),lib)
# This uses the utility program init_gen to create two
# up-to-date prereqs for intrinsic
$(S)/lit_parm.f9h $(S)/lit_add.f9h: $(S)/lit_names.txt $(INSTALLDIR)/init_gen
	sort -f -u $(S)/lit_names.txt | $(INSTALLDIR)/init_gen \
          -l last_intrinsic_lit -pl_ \
          -j lit_names.txt \
          $(S)/lit_parm.f9h $(S)/lit_add.f9h lit_indices
	echo "new lit_names" > new_lit_names.txt
	@for dir in $(EVERY_MLSCONFG); do \
          if [ -d ../$$dir ]; then \
            cp new_lit_names.txt ../$$dir 2>/dev/null; \
          fi; \
        done

$(S)/mol_parm.f9h $(S)/mol_add.f9h: $(TABLES_DIR)/sps_cross_ref_table.txt $(INSTALLDIR)/init_gen
	sed -n '/---\*/,$$ p' $(TABLES_DIR)/sps_cross_ref_table.txt | sed '1 d' | \
          sed '/H2O_R/ d' | awk ' {print $$2}' | sort -f -u| \
          $(INSTALLDIR)/init_gen \
          -f "last_rxx_molecule+1" -l last_molecule -pl_ \
          -j sps_cross_ref_table.txt \
          $(S)/mol_parm.f9h $(S)/mol_add.f9h lit_indices

machine.o: $(MACH_DIR)/machine.f90
	rm -f machine.f90 ; \
	sed '/Start $(BAD_MATCH)/,/End $(BAD_MATCH)/ d' $(MACH_DIR)/machine.f90 \
     > machine.f90
	$(FC) -c $(FOPTS) $(INC_PATHS) machine.f90 $(FAFTER)
	rm -f machine.f90

intrinsic.o: $(S)/lit_parm.f9h $(S)/lit_add.f9h

$(srclib)/Parser_L2CF.f9h: $(UTILDIR)/lr/l2cf.grm $(INSTALLDIR)/lr
	$(INSTALLDIR)/lr \
          $(UTILDIR)/lr/l2cf.grm \
          $(srclib)/Parser_L2CF.f9h $(UTILDIR)/lr/l2cf.lls $(LRAFTER)

endif # end short_name == lib

NOUNASS := $(shell echo ${FOPTS} | sed 's/\(--chk [aesx,]*\),u\([aesx,]*\)/\1\2/')

DUSTY_NO_OPT := $(shell echo ${DUSTY} | sed 's/-O[0-9]*//')
DONT_OPT := $(shell echo ${FOPTS} | sed 's/-O[0-9]*//; a-O0')

Hdf.o:
	$(FC) -c $(DUSTY) $(INC_PATHS) $(S)/Hdf.f90 $(FAFTER)
string_table.o:
	$(FC) -c $(BUGGY) $(INC_PATHS) $(S)/string_table.f90 $(FAFTER)

# readGriddedUtils module doesn't play well with the lf95 undefined variable
# runtime checker because of endian swapping
readGriddedUtils.o:
	$(FC) -c $(NOUNASS) $(INC_PATHS) $(S)/readGriddedUtils.f90 $(FAFTER)

endif  # end CASCADE
endif  # END not doc

# Some custom builds for the doc subdirectory
ifeq ($(short_name),doc)
# graphic for wvs-010
# Note this pattern for turning a .obj into LaTeX-useable files
# (the same pattern will be reused often)
# (1) an eps-type file
sparse.eps: sparse.obj

# (2) a pdf-type file
sparse.pdf: sparse.obj

wvs-010.dvi: wvs-010.tex sparse.eps

wvs-010.pdf: wvs-010.tex sparse.pdf

wvs-030-nearest.eps: wvs-030-nearest.obj

wvs-030-nearest.pdf: wvs-030-nearest.obj

wvs-030.dvi: wvs-030.tex wvs-030-nearest.eps

wvs-030.pdf: wvs-030.tex wvs-030-nearest.pdf

wvs-048-grid2.eps: wvs-048-grid2.obj

wvs-048-grid2.pdf: wvs-048-grid2.obj

wvs-048.dvi: wvs-048.tex wvs-048-grid2.eps

wvs-048.pdf: wvs-048.tex wvs-048-grid2.pdf

wvs-063-eq-circ.eps: wvs-063-eq-circ.obj

wvs-063-eq-circ.pdf: wvs-063-eq-circ.obj

wvs-063-grids_t.eps: wvs-063-grids_t.obj

wvs-063-grids_t.pdf: wvs-063-grids_t.obj

wvs-063-diag.eps: wvs-063-diag.obj

wvs-063-diag.pdf: wvs-063-diag.obj

wvs-063.dvi: wvs-063.tex wvs-063-eq-circ.eps wvs-063-grids_t.eps wvs-048-grid2.eps wvs-063-diag.eps

wvs-063.pdf: wvs-063.tex wvs-063-eq-circ.pdf wvs-063-grids_t.pdf wvs-048-grid2.pdf wvs-063-diag.pdf

wvs-074-scatter-1.eps: wvs-074-scatter.obj
	tgif -print -eps -page 1 -color wvs-074-scatter.obj
	mv wvs-074-scatter.eps wvs-074-scatter-1.eps

wvs-074-scatter-2.eps: wvs-074-scatter.obj
	tgif -print -eps -page 2 -color wvs-074-scatter.obj
	mv wvs-074-scatter.eps wvs-074-scatter-2.eps

wvs-074-scatter-1.pdf: wvs-074-scatter.obj
	tgif -print -pdf -page 1 -color wvs-074-scatter.obj
	mv wvs-074-scatter.pdf wvs-074-scatter-1.pdf

wvs-074-scatter-2.pdf: wvs-074-scatter.obj
	tgif -print -pdf -page 2 -color wvs-074-scatter.obj
	mv wvs-074-scatter.pdf wvs-074-scatter-2.pdf

wvs-074.dvi: wvs-074.tex wvs-074-scatter-1.eps wvs-074-scatter-2.eps

wvs-074.pdf: wvs-074.tex wvs-074-scatter-1.pdf wvs-074-scatter-2.pdf

wvs-083-grid.eps: wvs-083-grid.obj

wvs-083-grid.pdf: wvs-083-grid.obj

wvs-083.dvi: wvs-083.tex wvs-083-grid.eps

wvs-083.pdf: wvs-083.tex wvs-083-grid.pdf

wvs-089-pic.eps: wvs-089-pic.obj

wvs-089-pic.pdf: wvs-089-pic.obj

wvs-089.dvi: wvs-089.tex wvs-089-pic.eps

wvs-089.pdf: wvs-089.tex wvs-089-pic.pdf

wvs-095-eta.eps: wvs-095-eta.obj

wvs-095-eta.pdf: wvs-095-eta.obj

wvs-095.dvi: wvs-095.tex wvs-095-eta.eps
#	$(LATEX) wvs-095
#	$(LATEX) wvs-095

wvs-095.pdf: wvs-095.tex wvs-095-eta.pdf
#	pdf$(LATEX) wvs-095
#	pdf$(LATEX) wvs-095

wvs-111-MinLambda-not=0.eps: wvs-111-MinLambda-not=0.jpg

wvs-111-MinLambda=0.eps: wvs-111-MinLambda=0.jpg

wvs-111.dvi: wvs-111.tex wvs-111-MinLambda=0.eps wvs-111-MinLambda-not=0.eps

wvs-111.pdf: wvs-111.tex wvs-111-MinLambda=0.eps wvs-111-MinLambda-not=0.eps

wvs-126.dvi: wvs-126.tex wvs-126-1.eps

wvs-126.pdf: wvs-126.tex wvs-126-1.pdf

wvs-126-1.eps: wvs-126-1.obj

wvs-126-1.pdf: wvs-126-1.obj

wvs-128.dvi: wvs-128.tex wvs-128-QTM-2.eps wvs-128-QTM-3.eps

wvs-128.pdf: wvs-128.tex wvs-128-QTM-2.pdf wvs-128-QTM-3.pdf

wvs-128-QTM-2.eps: wvs-128-QTM-2.obj

wvs-128-QTM-2.pdf: wvs-128-QTM-2.obj

wvs-128-QTM-3.eps: wvs-128-QTM-3.obj

wvs-128-QTM-3.pdf: wvs-128-QTM-3.obj

wvs-133.dvi: wvs-133.tex wvs-133-fig.eps

wvs-133.pdf: wvs-133.tex wvs-133-fig.pdf

wvs-133-fig.eps: wvs-133-fig.obj

wvs-133-fig.pdf: wvs-133-fig.obj


wvs-134.dvi: wvs-134.tex wvs-134-1.eps

wvs-134.pdf: wvs-134.tex wvs-134-1.pdf

wvs-134-1.eps: wvs-134-1.obj

wvs-134-1.pdf: wvs-134-1.obj

wvs-136.dvi: wvs-136.tex wvs-136-lines.eps

wvs-136.pdf: wvs-136.tex wvs-136-lines.pdf

wvs-136-lines.eps: wvs-136-lines.obj

wvs-136-lines.pdf: wvs-136-lines.obj

wvs-137.dvi: wvs-137.tex wvs-137-1.eps wvs-137-2.eps wvs-137-3.eps

wvs-137.pdf: wvs-137.tex wvs-137-1.pdf wvs-137-2.pdf wvs-137-3.pdf

wvs-137-1.eps: wvs-137-1.obj

wvs-137-1.pdf: wvs-137-1.obj

wvs-137-2.eps: wvs-137-2.obj

wvs-137-2.pdf: wvs-137-2.obj

wvs-137-3.eps: wvs-137-3.obj

wvs-137-3.pdf: wvs-137-3.obj

wvs-138.dvi: wvs-138.tex wvs-138-1.eps

wvs-138.pdf: wvs-138.tex wvs-138-1.pdf

wvs-138-1.eps: wvs-138-1.obj

wvs-138-1.pdf: wvs-138-1.obj

wvs-139.dvi: wvs-139.tex wvs-139-1.eps

wvs-139.pdf: wvs-139.tex wvs-139-1.pdf

wvs-139-1.eps: wvs-139-1.obj

wvs-139-1.pdf: wvs-139-1.obj

wvs-141.dvi: wvs-141.tex wvs-141-1.eps

wvs-141.pdf: wvs-141.tex wvs-141-1.pdf

wvs-141-1.eps: wvs-141-1.obj

wvs-141-1.pdf: wvs-141-1.obj

wvs-142.dvi: wvs-142.tex wvs-142-1.eps

wvs-142.pdf: wvs-142.tex wvs-142-1.pdf

wvs-142-1.eps: wvs-142-1.obj

wvs-142-1.pdf: wvs-142-1.obj

wvs-143.dvi: wvs-143.tex wvs-143-1.eps wvs-143-2.eps

wvs-143.pdf: wvs-143.tex wvs-143-1.pdf wvs-143-2.pdf

wvs-143-1.eps: wvs-143-1.obj

wvs-143-1.pdf: wvs-143-1.obj

wvs-143-2.eps: wvs-143-2.obj

wvs-143-2.pdf: wvs-143-2.obj

wvs-144.dvi: wvs-144.tex Sparse.eps

wvs-144.pdf: wvs-144.tex Sparse.pdf

Sparse.eps: Sparse.obj

Sparse.pdf: Sparse.obj

wvs-145.dvi: wvs-145.tex wvs-145-1.eps wvs-145-2.eps

wvs-145.pdf: wvs-145.tex wvs-145-1.pdf wvs-145-2.pdf

wvs-145-1.eps: wvs-145-1.obj

wvs-145-1.pdf: wvs-145-1.obj

wvs-145-2.eps: wvs-145-2.obj

wvs-145-2.pdf: wvs-145-2.obj

wvs-146.dvi: wvs-146.tex wvs-146-1.eps wvs-146-2.eps wvs-146-3.eps

wvs-146.pdf: wvs-146.tex wvs-146-1.pdf wvs-146-2.pdf wvs-146-3.pdf

wvs-146-1.eps: wvs-146-1.obj

wvs-146-1.pdf: wvs-146-1.obj

wvs-146-2.eps: wvs-146-2.obj

wvs-146-2.pdf: wvs-146-2.obj

wvs-146-3.eps: wvs-146-3.obj

wvs-146-3.pdf: wvs-146-3.obj

wvs-147.dvi: wvs-147.tex \
             wvs-147-1.eps wvs-147-2.eps wvs-147-3.eps \
             wvs-147-4.eps wvs-147-5.eps wvs-147-6.eps \
             wvs-147-7.eps wvs-147-2-5.eps wvs-147-0.eps 

wvs-147.pdf: wvs-147.tex \
             wvs-147-1.pdf wvs-147-2.pdf wvs-147-3.pdf \
             wvs-147-4.pdf wvs-147-5.pdf wvs-147-6.pdf \
             wvs-147-7.pdf wvs-147-2-5.pdf wvs-147-0.pdf 

wvs-147-1.eps: wvs-147-1.obj

wvs-147-1.pdf: wvs-147-1.obj

wvs-147-2.eps: wvs-147-2.obj

wvs-147-2.pdf: wvs-147-2.obj

wvs-147-3.eps: wvs-147-3.obj

wvs-147-3.pdf: wvs-147-3.obj

wvs-147-0.eps: wvs-147-0.obj

wvs-147-0.pdf: wvs-147-0.obj

wvs-147-4.eps: wvs-147-4.obj

wvs-147-4.pdf: wvs-147-4.obj

wvs-147-5.eps: wvs-147-5.obj

wvs-147-5.pdf: wvs-147-5.obj

wvs-147-6.eps: wvs-147-6.obj

wvs-147-6.pdf: wvs-147-6.obj

wvs-147-7.eps: wvs-147-7.obj

wvs-147-7.pdf: wvs-147-7.obj

wvs-147-2-5.eps: wvs-147-2-5.obj

wvs-147-2-5.pdf: wvs-147-2-5.obj

wvs-149.dvi: wvs-149.tex wvs-149-1.txp wvs-149.bbl wvs-149-1.600pk

wvs-149.pdf: wvs-149.tex wvs-149-1.txp wvs-149.bbl wvs-149-1.600pk

wvs-149-1.txp: $(INSTALLDIR)/wvs-149
	$(INSTALLDIR)/wvs-149

wvs-149.aux: wvs-149.tex
	$(LATEX) wvs-149.tex

wvs-149.bbl: wvs-149.bib wvs-149.aux
	$(BIBTEX) wvs-149

wvs-149-1.600gf: wvs-149-1.txp
	mf '\mode=localfont; ' input wvs-149-1

wvs-149-1.600pk: wvs-149-1.600gf
	gftopk wvs-149-1.600gf

dnwt.ps: dnwt.dot
	dot -Tps dnwt.dot > dnwt.ps

dnwt.pdf: dnwt.ps
	ps2pdf dnwt.ps dnwt.pdf

wvs-150.dvi: wvs-150.tex dnwt.epsi

wvs-150.pdf: wvs-150.tex dnwt.pdf

wvs-151.dvi: wvs-151.tex wvs-151-1.eps wvs-151-2.eps wvs-151-3.eps

wvs-151.pdf: wvs-151.tex wvs-151-1.pdf wvs-151-2.pdf wvs-151-3.pdf

wvs-151-1.eps: wvs-151-1.obj

wvs-151-1.pdf: wvs-151-1.obj

wvs-151-2.eps: wvs-151-2.obj

wvs-151-2.pdf: wvs-151-2.obj

wvs-151-3.eps: wvs-151-3.obj

wvs-151-3.pdf: wvs-151-3.obj

dnwt.epsi: dnwt.ps
	ps2epsi dnwt.ps
        
wvs-154.dvi: wvs-154.tex wvs-154-1.eps

wvs-154.pdf: wvs-154.tex wvs-154-1.pdf

wvs-154-1.eps: wvs-154-1.obj

wvs-154-1.pdf: wvs-154-1.obj

wvs-155.dvi: wvs-155.tex Barycentric.eps wvs-151-QTM-1.eps

wvs-155.pdf: wvs-155.tex Barycentric.pdf wvs-151-QTM-1.pdf

Barycentric.eps: Barycentric.obj

Barycentric.pdf: Barycentric.obj

wvs-157.dvi: wvs-157.tex wvs-157-circ.eps wvs-146.dvi

wvs-157.pdf: wvs-157.tex wvs-157-circ.pdf wvs-146.pdf

wvs-157-circ.eps: wvs-157-circ.obj

wvs-157-circ.pdf: wvs-157-circ.obj

wvs-151-QTM-1.eps: wvs-151-QTM-1.obj

wvs-151-QTM-1.pdf: wvs-151-QTM-1.obj

wvs-160-1.eps: wvs-160-1.obj

wvs-160-1.pdf: wvs-160-1.obj

wvs-160-2.eps: wvs-160-2.obj

wvs-160-2.pdf: wvs-160-2.obj

wvs-160-facets.eps: wvs-160-facets.obj

wvs-160-facets.pdf: wvs-160-facets.obj

wvs-160-surf.eps: wvs-160-surf.obj

wvs-160-surf.pdf: wvs-160-surf.obj

wvs-160-cone.eps: wvs-160-cone.obj

wvs-160-cone.pdf: wvs-160-cone.obj

wvs-160-plane.eps: wvs-160-plane.obj

wvs-160-plane.pdf: wvs-160-plane.obj

wvs-160.dvi: wvs-160.tex wvs-160-1.eps wvs-160-2.eps wvs-160-facets.eps wvs-160-surf.eps wvs-160-cone.eps wvs-160-plane.eps

wvs-160.pdf: wvs-160.tex wvs-160-1.pdf wvs-160-2.pdf wvs-160-facets.pdf wvs-160-surf.pdf wvs-160-cone.pdf wvs-160-plane.pdf

# We discovered these by grep 'wvs-146' doc/*.tex
wvs-030.pdf wvs-131.pdf wvs-132.pdf wvs-157.pdf wvs-158.pdf wvs-159.pdf wvs-160.pdf: wvs-146.pdf

wvs-030.dvi wvs-131.dvi wvs-132.dvi wvs-157.dvi wvs-158.dvi wvs-159.dvi wvs-160.dvi: wvs-146.dvi

# We discovered these by grep 'externaldocument' doc/*.tex
wvs-030.pdf wvs-160.pdf: wvs-131.pdf

wvs-030.dvi wvs-160.dvi: wvs-131.dvi

wvs-102.pdf wvs-104.pdf : wvs-095.pdf wvs-100.pdf

wvs-102.dvi wvs-104.dvi : wvs-095.dvi wvs-100.dvi

wvs-160.pdf: wvs-146.pdf

wvs-160.dvi: wvs-146.dvi

wvs-128.aux: wvs-128.tex
	$(LATEX) wvs-128.tex

iy-005.aux: iy-005.tex
	$(LATEX) iy-005.tex

iy-005.dvi: iy-005.tex iy-005.bbl

iy-005.pdf: iy-005.tex iy-005.bbl

iy-005.bbl: yanovsky.bib iy-005.aux
	$(BIBTEX) iy-005

iy-006.aux: iy-006.tex
	$(LATEX) iy-006.tex

iy-006.dvi: iy-006.tex iy-006.bbl

iy-006.pdf: iy-006.tex iy-006.bbl

iy-006.bbl: yanovsky.bib iy-006.aux
	$(BIBTEX) iy-006


endif # end shortn_name == doc
# $Log$
# Revision 1.52  2020/06/10 21:01:51  pwagner
# Builds wvs-160 OK, now
#
# Revision 1.51  2020/04/21 23:14:08  pwagner
# Added more interdependencies among wvs in doc
#
# Revision 1.50  2020/04/21 00:07:27  pwagner
# Removed custom build for wvs-157; added more dependncies on wvs-146
#
# Revision 1.49  2020/01/15 00:23:38  pwagner
# Attempt to correct ties between wvs-157 and wvs-146
#
# Revision 1.48  2020/01/09 21:15:36  pwagner
# Can now build changed wvs-146.tex and new wvs-157.tex
#
# Revision 1.47  2019/12/23 18:09:34  pwagner
# Intel hates deferring size of parameter arrays
#
# Revision 1.46  2019/12/20 21:16:59  pwagner
# Builds .f9h files for lib/dates_module.f90
#
# Revision 1.45  2019/09/27 17:40:49  pwagner
# Builds wvs-155.tex with new wvs-151-QTM-1 dependencies (why name it that?)
#
# Revision 1.44  2019/09/24 18:04:39  pwagner
# Can now build wvs-155.tex
#
# Revision 1.43  2019/09/19 22:19:57  pwagner
# Use default build rules for .obj files
#
# Revision 1.42  2019/09/19 16:08:12  pwagner
# Added build commands for wvs-154
#
# Revision 1.41  2019/09/11 16:45:39  pwagner
# Augmented dependencies for wvs-063
#
# Revision 1.40  2019/09/10 17:37:47  pwagner
# Copes with new wvs-063-eq-circ.obj
#
# Revision 1.39  2019/06/18 23:37:22  pwagner
# Can now build iy-005.tex in doc
#
# Revision 1.38  2019/06/18 20:26:51  pwagner
# Added cmds to build iy-006.tex in doc
#
# Revision 1.37  2019/05/31 20:07:51  pwagner
# Needed actual build cmds for wvs-111
#
# Revision 1.36  2019/05/31 17:19:10  pwagner
# We can now build wvs-111 and wvs-139
#
# Revision 1.35  2018/11/30 22:31:43  pwagner
# Can now build wvs-147 and wvs-151
#
# Revision 1.34  2018/07/27 00:21:02  pwagner
# Now builds wvs-150.tex
#
# Revision 1.33  2018/07/11 16:46:07  pwagner
# Can now build wvs-149
#
# Revision 1.32  2018/06/04 23:04:05  pwagner
# Incorporates changes to build wvs-138
#
# Revision 1.31  2017/10/20 18:16:27  pwagner
# Added build cmds for wvs-145.tex
#
# Revision 1.30  2017/10/13 16:33:48  pwagner
# Added build cmds for wvs-143
#
# Revision 1.29  2017/09/25 23:20:15  pwagner
# Can now build wvs-144
#
# Revision 1.28  2017/09/22 17:58:52  pwagner
# Added wvs-142
#
# Revision 1.27  2017/09/06 18:59:05  pwagner
# Now builds wvs-141.tex successfully
#
# Revision 1.26  2017/03/28 23:58:10  pwagner
# Now builds  doc/wvs-137 properly
#
# Revision 1.25  2017/03/08 00:11:55  pwagner
# ncep_dao doesn't need custom, but readGriddedUtils does
#
# Revision 1.24  2016/11/04 23:26:58  pwagner
# Now builds  doc/wvs-136 properly
#
# Revision 1.23  2016/06/08 00:00:41  pwagner
# Prevent creating unwanted text files in lib, etc. with names of machine directories
#
# Revision 1.22  2016/05/31 23:41:16  pwagner
# Tried again to fix need to build mlsl2 twice; added graphic to wvs-030
#
# Revision 1.21  2016/04/06 22:15:21  pwagner
# Now builds  doc/wvs-134 properly
#
# Revision 1.20  2016/03/15 18:38:22  pwagner
# Added wvs-133
#
# Revision 1.19  2015/11/18 17:45:06  pwagner
# Now builds wvs-128.tex
#
# Revision 1.18  2015/09/18 23:38:34  pwagner
# Added custom builds for wvs-126.tex
#
# Revision 1.17  2015/04/23 17:51:38  whdaffer
# Added some comments to make the conditional includes a little easier
# to understand
#
# Revision 1.16  2014/06/16 20:27:21  pwagner
# Fixed bug in building lr
#
# Revision 1.15  2014/05/24 00:23:52  vsnyder
# Make Parser_L2CF.f9h in the right place
#
# Revision 1.14  2014/05/22 01:18:44  vsnyder
# Moved Parser_Tables_L2CF.f90 from lib to srclib
#
# Revision 1.13  2014/05/20 23:41:29  vsnyder
# New parser, simplified build
#
# Revision 1.12  2014/02/21 22:32:28  pwagner
# Uses Parser_Tables.wc to trigger two-phase build
#
# Revision 1.11  2014/01/29 21:03:41  pwagner
# Shorten make update
#
# Revision 1.10  2014/01/24 00:58:08  pwagner
# Removed certain amount of clumsiness
#
# Revision 1.9  2014/01/22 18:22:15  pwagner
# More changes, bug fixes
#
# Revision 1.8  2014/01/18 00:56:45  pwagner
# Redirect lr stdout
#
# Revision 1.7  2014/01/15 19:12:19  pwagner
# Compatible with new lr
#
# Revision 1.6  2013/12/03 00:03:29  pwagner
# Fixed bugs in building lr; added LRAFTER
#
# Revision 1.5  2013/10/28 23:12:17  pwagner
# May build parser.f9h
#
# Revision 1.4  2013/08/28 00:39:46  pwagner
# Lifted custom builds from MLSStrings.f90 MLSStringLists.f90
#
# Revision 1.3  2013/06/12 18:18:04  pwagner
# Changes to pass most FOPTS for Strings modules
#
# Revision 1.2  2012/04/20 00:43:26  pwagner
# Fixed NOUNASS so it works again with NAG
#
# Revision 1.1  2012/04/04 00:45:40  pwagner
# First commit
#
