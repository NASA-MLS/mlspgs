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
# Two useful utilities
$(INSTALLDIR)/init_gen: $(UTILDIR)/init_gen.f90
	$(UTILDIR)/build_f90_in_misc.sh -p init_gen \
	-d $(INSTALLDIR) -t $(TESTSDIR) -M $(MAKE) \
	-C $(MLSCFILE) \
   $(UTILDIR)/init_gen.f90
$(INSTALLDIR)/lr: $(UTILDIR)/lr/*.f90
	$(UTILDIR)/build_f90_in_misc.sh -d $(INSTALLDIR) -t $(TESTSDIR) \
	-c $(MLSCONFG) -p lr -M $(MAKE) -O short_name=lr_custom -m lib \
	-C $(MLSCFILE) $(UTILDIR)/lr/*.[Ff]9[0h]

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

endif
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
$(S)/lit_parm.f9h $(S)/lit_add.f9h: $(S)/lit_names.txt $(INSTALLDIR)/init_gen
	sort -f -u $(S)/lit_names.txt | $(INSTALLDIR)/init_gen \
          -l last_intrinsic_lit -pl_ \
          -j lit_names.txt \
          $(S)/lit_parm.f9h $(S)/lit_add.f9h lit_indices

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

$(srclib)/Parser_L2CF.f9h: $(UTILDIR)/lr/l2cf.grm $(INSTALLDIR)/lr
	$(INSTALLDIR)/lr \
          $(UTILDIR)/lr/l2cf.grm \
          $(srclib)/Parser_L2CF.f9h $(UTILDIR)/lr/l2cf.lls $(LRAFTER)

endif

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

# ncep_dao module doesn't play well with the lf95 undefined variable
# runtime checker because of endian swapping
ncep_dao.mod:
	$(UTILDIR)/newAifBdiff.sh -a ncep_dao.mod $(FC) -c $(NOUNASS) $(INC_PATHS) $(S)/ncep_dao.f90 $(FAFTER)
ncep_dao.o: ncep_dao.mod
	$(UTILDIR)/mark_as_uptodate.sh -M $(MAKE) -t \
   -T ncep_dao.o ncep_dao.mod 

# Intel's ifort compiler v11.x for 32-bit architecture has a string-handling bug
# unless treated very carefully
#MLSStringLists.o:
#	$(UTILDIR)/mark_as_uptodate.sh -M $(MAKE) -t -T MLSStringLists.o mlsstringlists.mod 
#mlsstringlists.mod: MLSStringLists.f90 mlscommon.mod mlsmessagemodule.mod \
#	mlssets.mod mlsstrings.mod
#	$(UTILDIR)/newAifBdiff.sh -a mlsstringlists.mod $(FC) -c $(DONT_OPT) $(PRE) $(INC_PATHS) $(S)/MLSStringLists.f90 $(FAFTER)
#MLSStrings.o: mlsstrings.mod MLSStrings.f90 
#	$(UTILDIR)/mark_as_uptodate.sh -M $(MAKE) -t -T MLSStrings.o mlsstrings.mod 
#mlsstrings.mod: MLSStrings.f90 mlscommon.mod ReadANumFromChars.f9h
#	$(UTILDIR)/newAifBdiff.sh -a mlsstrings.mod $(FC) -c $(DONT_OPT) $(PRE) $(INC_PATHS) $(S)/MLSStrings.f90 $(FAFTER)

else
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

endif
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

endif

NOUNASS := $(shell echo ${FOPTS} | sed 's/\(--chk [aesx,]*\),u\([aesx,]*\)/\1\2/')

DUSTY_NO_OPT := $(shell echo ${DUSTY} | sed 's/-O[0-9]*//')
DONT_OPT := $(shell echo ${FOPTS} | sed 's/-O[0-9]*//; a-O0')

Hdf.o:
	$(FC) -c $(DUSTY) $(INC_PATHS) $(S)/Hdf.f90 $(FAFTER)
string_table.o:
	$(FC) -c $(BUGGY) $(INC_PATHS) $(S)/string_table.f90 $(FAFTER)

# ncep_dao module doesn't play well with the lf95 undefined variable
# runtime checker because of endian swapping
ncep_dao.o:
	$(FC) -c $(NOUNASS) $(INC_PATHS) $(S)/ncep_dao.f90 $(FAFTER)

# Intel's ifort compiler v11.x for 32-bit architecture has a string-handling bug
# unless treated very carefully
#MLSStringLists.o:
#	$(FC) -c $(DONT_OPT) $(PRE) $(INC_PATHS) $(S)/MLSStringLists.f90 $(FAFTER)
#MLSStrings.o:
#	$(FC) -c $(DONT_OPT) $(PRE) $(INC_PATHS) $(S)/MLSStrings.f90 $(FAFTER)

endif
endif

# Some custom builds for the doc subdirectory
ifeq ($(short_name),doc)
# graphic for wvs-010
sparse.eps: sparse.obj
	tgif -print -eps -page 1 -color sparse.obj

# graphic for wvs-010
sparse.pdf: sparse.obj
	tgif -print -pdf -page 1 -color sparse.obj

wvs-010.dvi: wvs-010.tex sparse.eps

wvs-010.pdf: wvs-010.tex sparse.pdf

wvs-048-grid2.eps: wvs-048-grid2.obj
	tgif -print -eps -page 1 -color wvs-048-grid2.obj

wvs-048-grid2.pdf: wvs-048-grid2.obj
	tgif -print -pdf -page 1 -color wvs-048-grid2.obj

wvs-048.dvi: wvs-048.tex wvs-048-grid2.eps

wvs-048.pdf: wvs-048.tex wvs-048-grid2.pdf

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
	tgif -print -eps -page 1 -color wvs-083-grid.obj

wvs-083-grid.pdf: wvs-083-grid.obj
	tgif -print -pdf -page 1 -color wvs-083-grid.obj

wvs-083.dvi: wvs-083.tex wvs-083-grid.eps

wvs-083.pdf: wvs-083.tex wvs-083-grid.pdf

wvs-089-pic.eps: wvs-089-pic.obj
	tgif -print -eps -page 1 -color wvs-089-pic.obj

wvs-089-pic.pdf: wvs-089-pic.obj
	tgif -print -pdf -page 1 -color wvs-089-pic.obj

wvs-089.dvi: wvs-089.tex wvs-089-pic.eps

wvs-089.pdf: wvs-089.tex wvs-089-pic.pdf

wvs-095-eta.eps: wvs-095-eta.obj
	tgif -print -eps -page 1 -color wvs-095-eta.obj

wvs-095-eta.pdf: wvs-095-eta.obj
	tgif -print -pdf -page 1 -color wvs-095-eta.obj

wvs-095.dvi: wvs-095.tex wvs-095-eta.eps
#	latex wvs-095
#	latex wvs-095

wvs-095.pdf: wvs-095.tex wvs-095-eta.pdf
#	pdflatex wvs-095
#	pdflatex wvs-095
endif
# $Log$
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
