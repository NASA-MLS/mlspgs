# $Id$

L2OBJS = Construct.o ConstructQuantityTemplates.o \
       ConstructVectorTemplates.o dumper.o Fill.o get_3d_view.o \
       global_settings.o HGrid.o init_tables_module.o Join.o L2AUXData.o \
       L2GPData.o MLSL2.o obtainclimatology.o obtaindao.o obtain_mlscf.o \
       obtainncep.o open_init.o OutputAndClose.o OutputL2GP.o \
       QuantityTemplates.o ScanDivide.o tree_walker.o VGrid.o VectorsModule.o

all: mlsl2

Construct.o: $S/Construct.f90 ConstructQuantityTemplates.o \
          ConstructVectorTemplates.o dumper.o HGrid.o init_tables_module.o \
          $M/MLSMessageModule.o $M/MLSSignalNomenclature.o QuantityTemplates.o \
          $C/toggles_core.o $C/trace_m.o $C/tree.o $C/tree_types.o \
          VectorsModule.o VGrid.o
	$(FC) -c $(FOPTS) $I$C $I$S $I$M $S/Construct.f90

ConstructQuantityTemplates.o: $S/ConstructQuantityTemplates.f90 \
          $C/declaration_table.o HGrid.o init_tables_module.o $C/lexer_core.o \
          $M/MLSMessageModule.o $M/MLSSignalNomenclature.o $C/output_m.o \
          QuantityTemplates.o $C/string_table.o $C/toggles_core.o \
          $C/trace_m.o $C/tree.o VGrid.o
	$(FC) -c $(FOPTS) $I$C $I$S $I$M $S/ConstructQuantityTemplates.f90

ConstructVectorTemplates.o: $S/ConstructVectorTemplates.f90 \
          $M/Allocate_Deallocate.o init_tables_module.o $M/MLSMessageModule.o \
          QuantityTemplates.o $C/toggles_core.o $C/trace_m.o $C/tree.o \
          VectorsModule.o
	$(FC) -c $(FOPTS) $I$C $I$S $I$M $S/ConstructVectorTemplates.f90

dumper.o: $S/dumper.f90 $M/dump_0.o HGrid.o init_tables_module.o \
          $C/output_m.o QuantityTemplates.o $C/string_table.o VGrid.o
	$(FC) -c $(FOPTS) $I$C $I$S $I$M $S/dumper.f90

Fill.o: $S/Fill.f90 dumper.o init_tables_module.o $C/lexer_core.o \
          $C/output_m.o $C/toggles_core.o $C/trace_m.o $C/tree.o \
          $C/tree_types.o VectorsModule.o
	$(FC) -c $(FOPTS) $I$C $I$S $I$M $S/Fill.f90

get_3d_view.o: $S/get_3d_view.f90
	$(FC) -c $(FOPTS) $I$M $S/get_3d_view.f90

global_settings.o: $S/global_settings.f90 init_tables_module.o \
          $C/toggles_core.o $C/trace_m.o $C/tree.o
	$(FC) -c $(FOPTS) $I$C $I$S $I$M $S/global_settings.f90

HGrid.o: $S/HGrid.f90 $S/addItemToDatabase.f9h $C/expr_m.o \
          init_tables_module.o $C/lexer_core.o $C/output_m.o \
          $C/string_table.o $C/toggles_core.o $C/trace_m.o $C/tree.o
	$(FC) -c $(FOPTS) $I$C $I$S $I$M $S/HGrid.f90

init_tables_module.o: $S/init_tables_module.f90 $C/symbol_table.o \
          $C/symbol_types.o $C/tree.o $C/tree_types.o
	$(FC) -c $(FOPTS) $I$C $I$S $I$M $S/init_tables_module.f90

Join.o: $S/Join.f90 init_tables_module.o L2AUXData.o L2GPData.o \
          $C/lexer_core.o $M/MLSMessageModule.o $C/output_m.o \
          QuantityTemplates.o $C/tree.o $C/tree_types.o VectorsModule.o
	$(FC) -c $(FOPTS) $I$C $I$S $I$M $S/Join.f90

L2AUXData.o: $S/L2AUXData.f90 $S/addItemToDatabase.f9h \
          $M/Allocate_Deallocate.o $M/MLSMessageModule.o
	$(FC) -c $(FOPTS) $I$S $I$M $S/L2AUXData.f90

L2GPData.o: $S/L2GPData.f90 $S/addItemToDatabase.f9h $M/Allocate_Deallocate.o \
          $M/MLSMessageModule.o
	$(FC) -c $(FOPTS) $I$S $I$M $S/L2GPData.f90

MatrixModule.o: $S/MatrixModule.f90 $S/addItemToDatabase.f9h \
          $M/Allocate_Deallocate.o dump_0.o $M/MLSMessageModule.o \
          $C/output_m.o $C/string_table.o
	$(FC) -c $(FOPTS) $I$C $I$S $I$M $S/MatrixModule.f90

MLSL2.o: $S/MLSL2.f90 $C/declaration_table.o init_tables_module.o \
          $C/lexer_core.o machine.o open_init.o $C/parser.o $C/string_table.o \
          $C/toggles_core.o $C/tree.o $C/ tree_walker.o $C/units.o
	$(FC) -c $(FOPTS) $I$C $I$S $I$M $S/MLSL2.f90

obtainclimatology.o: $M/Allocate_Deallocate.o $S/obtainclimatology.f90 \
          $M/MLSMessageModule.o
	$(FC) -c $(FOPTS) $I$C $I$S $I$M $S/obtainclimatology.f90

obtaindao.o: $S/obtaindao.f90 $M/MLSMessageModule.o
	$(FC) -c $(FOPTS) $I$C $I$S $I$M $S/obtaindao.f90

obtain_mlscf.o: $S/obtain_mlscf.f90 $M/MLSMessageModule.o
	$(FC) -c $(FOPTS) $I$C $I$S $I$M $S/obtain_mlscf.f90

obtainncep.o: $S/obtainncep.f90 $M/MLSMessageModule.o
	$(FC) -c $(FOPTS) $I$C $I$S $I$M $S/obtainncep.f90

open_init.o: $S/open_init.f90 $M/MLSMessageModule.o $M/MLSSignalNomenclature.o \
          $C/string_table.o
	$(FC) -c $(FOPTS) $I$C $I$S $I$M $S/open_init.f90

OutputAndClose.o: $S/OutputAndClose.f90 $C/init_tables_module.o \
          L2AUXData.o L2GPData.o $C/lexer_core.o $M/MLSMessageModule.o \
          $C/output_m.o OutputL2GP.o $C/string_table.o $C/tree.o \
          $C/tree_types.o
	$(FC) -c $(FOPTS) $I$C $I$S $I$M $S/OutputAndClose.f90

OutputL2GP.o: $S/OutputL2GP.f90 $M/HDF-EOS.o $M/Allocate_Deallocate.o \
          L2GPData.o $M/MLSMessageModule.o $M/swapi.o
	$(FC) -c $(FOPTS) $I$C $I$S $I$M $S/OutputL2GP.f90

QuantityTemplates.o: $S/QuantityTemplates.f90 $S/addItemToDatabase.f9h \
          $M/Allocate_Deallocate.o $M/MLSMessageModule.o \
          $M/MLSSignalNomenclature.o
	$(FC) -c $(FOPTS) $I$C $I$S $I$M $S/QuantityTemplates.f90

ScanDivide.o: $S/ScanDivide.f90 dumper.o $C/expr_m.o init_tables_module.o \
          $M/MLSMessageModule.o $C/string_table.o $C/tree.o
	$(FC) -c $(FOPTS) $I$C $I$S $I$M $S/ScanDivide.f90

tree_walker.o: $S/tree_walker.f90 Construct.o $C/declaration_table.o dumper.o \
          global_settings.o Fill.o init_tables_module.o Join.o L2AUXData.o \
          L2GPData.o $C/lexer_core.o obtainclimatology.o \
          obtaindao.o obtainncep.o \
          open_init.o $C/output_m.o OutputAndClose.o QuantityTemplates.o \
          ScanDivide.o \
          $C/string_table.o $C/toggles_core.o $C/trace_m.o $C/tree.o \
          $C/tree_types.o VectorsModule.o
	$(FC) -c $(FOPTS) $I$C $I$S $I$M $S/tree_walker.f90

VectorsModule.o: $S/VectorsModule.f90 $S/addItemToDatabase.f9h \
          $M/Allocate_Deallocate.o dumper.o get_3d_view.o \
          $M/MLSMessageModule.o $M/MLSSignalNomenclature.o $C/output_m.o \
          QuantityTemplates.o
	$(FC) -c $(FOPTS) $I$C $I$S $I$M $S/VectorsModule.f90

VGrid.o: $S/VGrid.f90 $S/addItemToDatabase.f9h $C/expr_m.o \
          init_tables_module.o $C/lexer_core.o $C/output_m.o $C/string_table.o \
          $C/toggles_core.o $C/trace_m.o $C/tree.o
	$(FC) -c $(FOPTS) $I$C $I$S $I$M $S/VGrid.f90

# ========================================================================
# Temporary until the toolkit gets installed:
pgs_smf_generatestatusreport.o: $S/pgs_smf_generatestatusreport.f
	$(FC) -c $(FOPTS) $I$C $I$S $I$M $S/pgs_smf_generatestatusreport.f
PG =
#PG = pgs_smf_generatestatusreport.o
#PG = pgs_smf_generatestatusreport.o \
#     $(HOME)/mlspgs/SDPToolkit/src/IO/GEN/PGS_IO_Gen_OpenF.o \
#     $(HOME)/mlspgs/SDPToolkit/src/IO/GEN/PGS_IO_Gen_CloseF.o
#SDP = $M/SDPToolkit.o
#SDP = $M/
SDP =SDPToolkitSubstitute.o

GDapi.o: $L/hdfeos/src/GDapi.c
	cc -c -DNAGf90Fortran -I/van1/vsnyder/mlspgs/SDPToolkit/include/ \
           -I$L/hdf/linux/HDF4.1r3/include/ -I$L/hdf/linux/HDF4.1r3/include \
           -I$L/hdfeos/include/ $L/hdfeos/src/GDapi.c
GD=GDapi.o
# ========================================================================
#-L $L/lib/$P -lPGSTK \

mlsl2: $(L2OBJS) $(GD) $(PG) $(SDP) $C/libCF.a $M/libmls.a
	$(FC) -o mlsl2 $(L2OBJS) \
                      $(SDP) $(PG) $(GD) \
                      -L$C -lCF \
                      -L$M -lmls \
                      -L$(HOME)/mlspgs/SDPToolkit -lPGSTK \
		      -L$E/hdfeos/lib/$P -lhdfeos \
                      -L$(HOME)/mlspgs/l2/GCT -lGCT \
                      -L$L/hdf/$P/HDF4.1r3/lib -lmfhdf -ldf -ljpeg -lz

# ========================================================================

clean:
	rm -f *.o *.mod mlsl2

tar:
	rm -f parser.tgz
	tar hcvzf parser.tgz * --exclude \*.o --exclude \*.mod --exclude test

# ========================================================================

# $Log$
