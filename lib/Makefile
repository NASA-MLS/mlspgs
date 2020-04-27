#	Makefile
#
#  Lost? Type 'cd ../; make help -f MakeFC' from this directory
#
# Copyright 2005, by the California Institute of Technology. ALL
# RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
# commercial use must be negotiated with the Office of Technology Transfer
# at the California Institute of Technology.

# This software may be subject to U.S. export control laws. By accepting this
# software, the user agrees to comply with all applicable U.S. export laws and
# regulations. User has the responsibility to obtain export licenses, or other
# export authority as may be required before exporting such information to
# foreign countries or providing access to foreign persons.

# "$Id$"

short_name=lib

PROG=libmls.a

# We store all the common Makefile things in one shared file
SRCLIB=../srclib
-include $(SRCLIB)/Makefile.h

# This will set our MLS platform, compiler. etc.
# w/o complaining if it does nor exist yet
-include $(CONFDIR)/$(MLSCFILE)

#the subdirectory where testing programs were kept
#(outdated now)
testdir = tests

testsubdirs = $(testdir)/multiply $(testdir)/sparsify $(testdir)/test_cf_parser

OMIT_SRCS := $(OMIT) CCSDS_Time.f90 MLSHashes.f90

#Toolkit-free version requires *Substitute.f90
ifdef PGSTK
   OMIT_SRCS := $(OMIT_SRCS) PseudoToolkit.f90 SDPToolkitSubstitute.f90 \
   MLSMessageSubstitute.f90
endif

OMIT_PREFXD := $(shell ${REECHO} -prefix=-d $(OMIT_SRCS))

#These modules require custom builds, so omit build commands from dependencies
CUSTOM_BUILDS := Hdf.f90 string_table.f90 readGriddedUtils.f90
CUSTOM_PREFXD := $(shell ${REECHO} -prefix=-db $(CUSTOM_BUILDS))

ifdef HDFINC
   IHDFINC := -I $(HDFINC)
else
   IHDFINC := -I $(HDF)/../include
endif

ifdef HDF5INC
   IHDF5INC := -I $(HDF5INC)
   ADD2PATHS := $(ADD2PATHS) $(HDF5INC)
else
   IHDF5INC := -I $(HDF5)/../fortran/src -I $(HDF5)/../include
   ADD2PATHS := $(ADD2PATHS) $(HDF5)/../fortran/src $(HDF5)/../include
endif

ifdef HDFEOSINC
   IHDFEOSINC := -I $(HDFEOSINC)
else
   IHDFEOSINC := -I $(HDFEOS)/../../include
endif

ifdef HDFEOS5INC
   IHDFEOS5INC := $(HDFEOS5INC)/hdfeos5.inc
else
   IHDFEOS5INC := $(HDFEOS5)/../../include/hdfeos5.inc
endif

#----------------------- Build instructions

$(PROG): $(CONFDIR)/$(MLSCFILE) $(PLATFORMS)/$(MLSCONFG) \
         $(libblas_objs)/libmlspack.a hdfeos5.f9h
	@cd $(PLATFORMS)/$(MLSCONFG) && $(MAKE) $(PROG) \
      UNPREFXDPTHS="$(ADD2PATHS)" ; \
        return_status=`expr $$?` ; \
        if [ $$return_status != "0" ]; then \
		     exit 2 ; \
        elif [ "$(MARK_ALL_AS_UPTODATE)" != "no" ]; then \
		     $(MAKE) mark_all_as_uptodate \
           UNPREFXDPTHS="$(ADD2PATHS)" MARK_ALL_AS_UPTODATE= ;\
         fi 

hdfeos5.f9h: $(IHDFEOS5INC)
	@sed 's/^C/!/' $(IHDFEOS5INC) > hdfeos5.f9h

.PHONY: configure depends clean alltests tar tarz targ mostlyclean \
	all libmls.a ghostbuster update

tar:

configure:
	$(MLSBIN)/mlsconfigure -dc "$(CONFDIR)" -pc "$(PLATFORMS)" \
          -c "$(MLSF95)" -p "$(MLSPLAT)" -cf "$(MLSCFILE)" \
          $(RECFG)

$(CONFDIR)/$(MLSCFILE):
	$(MLSBIN)/mlsconfigure -dc "$(CONFDIR)" -pc "$(PLATFORMS)" \
          -c "$(MLSF95)" -p "$(MLSPLAT)" -cf "$(MLSCFILE)" \
          $(RECFG)

$(PLATFORMS)/$(MLSCONFG): $(CONFDIR)/$(MLSCFILE)
	$(MLSBIN)/mlsconfigure -dc "$(CONFDIR)" -pc "$(PLATFORMS)" \
          -c "$(MLSF95)" -p "$(MLSPLAT)" -cf "$(MLSCFILE)" -nc  \
          $(RECFG)

depends:
	if [ "$(PREPRO)" != "" ] ; then make $(PREPRO) ; fi
	$(MAKE) depends--mod
	if [ "$(PREPRO)" != "" ] ; then \
        temp=""; \
        for f in $(PREPRO); do temp="$$temp -e s/$$f/`basename $$f .f90`.F90/"; \
        done; \
        sed -i $$temp $(PLATFORMS)/$(MLSCONFG)/Makefile.dep; \
        rm -f $(PREPRO); \
        fi

depends--mod:
	$(MLSBIN)/makemakedep.sh -p $(PLATFORMS)/$(MLSCONFG)/Makefile.dep \
   $(OMIT_PREFXD) $(UP_OR_LOW) \
   $(CUSTOM_PREFXD) ./machines/$(MLSCONFG) $(srclib) $(libblas_sources)

ghostbuster:
	$(MLSBIN)/ghostbuster.sh -lib $(OMIT_PREFXD) \
           $(PLATFORMS)/$(MLSCONFG) ./ $(srclib) ./machines/$(MLSCONFG)
           
update:
	@cd .. ; \
   $(MAKE) update SUBDIRS="blas lib" tests_sub_dirs= JUSTSUBDIRS=yes

mostlyclean:
	@echo "Cleaning directory $(MLSCONFG)" ; \
	if [ "$(MLSCONFG)" != "" ] ; then  \
		echo "Deleting $(MLSCONFG) files" ; \
		rm -f $(PLATFORMS)/$(MLSCONFG)/* ; \
	fi ; \
	if [ -d "$(PLATFORMS)/$(MLSCONFG)" ] ; then  \
		echo "Deleting $(MLSCONFG) directory" ; \
		rmdir $(PLATFORMS)/$(MLSCONFG) ; \
	fi ; \
	rm -f Make.dep.*

clean:
	@echo "Cleaning objects from directory $(MLSCONFG)" ; \
	if [ "$(MLSCONFG)" != "" ] ; then  \
		echo "Deleting $(MLSCONFG) files" ; \
		rm -f $(PLATFORMS)/$(MLSCONFG)/*.[ao] ; \
		rm -f $(PLATFORMS)/$(MLSCONFG)/*.mod ; \
	fi; \

#	echo "Cleaning objects from $(testdir)"; \
#	$(MAKE) -f $(MakeFName) -C $(testdir) clean
	rm -f Make.dep.*

tar:
	rm -f ../lib.tar
	@echo "tar-ing lib"; \
	$(TAR) hcvf ../lib.tar * --exclude CVS\
	--exclude Make.dep\* --exclude lib/NAG.\* --exclude lib/LF95.\*

tarz: tar
	@echo "... then gzipping it"; \
	cd ..; rm -f lib.tgz; gzip lib.tar; mv lib.tar.gz lib.tgz

targ:
	rm -f ../lib.tgz
	$(TAR) hcvzf ../lib.tgz * --exclude \*.o --exclude \*.mod

#================= For testing purposes only =============
#alltests: $(testsubdirs)
alltests:
	$(MAKE) -f $(MakeFName) -C $(testdir)

$(testsubdirs):
	echo "Making $@"
	$(MAKE) -f $(MakeFName) -C $@
# $Log$
# Revision 2.3  2017/03/08 00:12:19  pwagner
# ncep_dao doesn't need custom, but readGriddedUtils does
#
# Revision 2.2  2016/03/02 00:56:12  pwagner
# Omits compiling PseudoToolkit.f90
#
# Revision 2.1  2015/08/10 23:20:33  pwagner
# After switching away from MakeFC
#
# Revision 2.68  2014/05/20 23:54:17  vsnyder
# parser.f90 no longer needs a custom build
#
# Revision 2.67  2014/02/21 22:32:00  pwagner
# Uses Parser_Tables.wc to trigger two-phase build
#
# Revision 2.66  2014/01/29 21:01:28  pwagner
# Correct various bugs; shorten make update
#
# Revision 2.65  2014/01/28 19:49:20  pwagner
# set UNPREFXDPTHS when building parser and Parser_Tables
#
# Revision 2.64  2014/01/24 00:58:47  pwagner
# Removed certain amount of clumsiness
#
# Revision 2.63  2014/01/22 18:22:44  pwagner
# More changes, bug fixes
#
# Revision 2.62  2014/01/18 00:55:40  pwagner
# Could this be the fix to building Parser_Tables we sought?
#
# Revision 2.61  2013/08/28 00:39:15  pwagner
# Lifted custom builds from MLSStrings.f90 MLSStringLists.f90
#
# Revision 2.60  2013/06/05 18:40:51  pwagner
# Remove vulnerability if temp already exists in lib
#
# Revision 2.59  2010/10/12 22:24:29  pwagner
# Uses Makefile.h from srclib
#
# Revision 2.58  2010/08/27 20:56:39  pwagner
# Uses optional HDF*INC Makefile variables
#
# Revision 2.57  2009/10/08 23:07:24  pwagner
# Added workarounds for buggy Intel v11 handling of strings
#
# Revision 2.56  2008/07/11 23:59:19  pwagner
# 1st changes to get sunstudio to link mlsl2
#
# Revision 2.55  2008/04/22 17:16:45  pwagner
# CASCADE no longer affects where Makefile.dep goes
#
# Revision 2.54  2007/05/31 18:31:03  pwagner
# Added Absoft to list of tested compilers
#
# Revision 2.53  2007/04/14 00:41:12  vsnyder
# Preprocess .F90 and .F9h files
#
# Revision 2.52  2007/04/13 20:15:18  pwagner
# Changes in line with Intel v9.1
#
# Revision 2.51  2007/01/12 00:23:02  pwagner
# Made way for MLSMessageSubstitute
#
# Revision 2.50  2005/06/23 22:10:11  pwagner
# Reworded Copyright statement
#
# Revision 2.49  2003/05/20 18:39:23  pwagner
# Uses DONTBUILDPREQS to stop some hunting
#
# Revision 2.48  2003/04/30 23:04:55  pwagner
# Now lets you make update in each subdirectory
#
# Revision 2.47  2003/02/24 21:16:01  pwagner
# Fixed bug where changes to .f9h were not triggering recompiles
#
# Revision 2.46  2002/09/13 22:54:17  pwagner
# Dependencies cross over into blas
#
# Revision 2.45  2002/08/23 20:33:05  pwagner
# Removed lib trickery L2GPData, MLSFiles, etc.
#
# Revision 2.44  2002/07/18 18:00:02  pwagner
# NOHUNTING default when .mod-mediated dependency
#
# Revision 2.43  2002/07/15 17:39:10  pwagner
# Custom builds for L2PC_m.f90 MLSFiles.f90 L2GPData.f90 so LF95 finishes in finite time
#
# Revision 2.42  2002/07/15 16:52:59  pwagner
# Exits make with status 2 if error in any level
#
# Revision 2.41  2002/07/11 22:12:11  pwagner
# Changes to support hdf5/hdfeos5
#
# Revision 2.40  2002/07/02 00:14:01  pwagner
# Skips marking uptodate after an error
#
# Revision 2.39  2002/06/21 23:24:09  pwagner
# More steps enroute to ending cascades of recompiles
#
# Revision 2.38  2002/06/21 00:15:28  pwagner
# Enroute to ending cascading recompiles
#
# Revision 2.37  2002/05/02 00:16:32  pwagner
# depends now creates Makefile.nodep, too
#
# Revision 2.36  2002/01/29 00:46:47  pwagner
# Builds he5lib when HDFVERSIONS defined
#
# Revision 2.35  2002/01/15 00:40:18  pwagner
# Now uses user-settable OMIT and internal OMIT_SRCS and OMIT_PREFIXD variables
#
# Revision 2.34  2001/11/15 00:09:24  pwagner
# Rebuilds on changes to mlspgs/blas
#
# Revision 2.33  2001/11/09 21:52:10  pwagner
# Dropped testdir from targets automatically cleaned
#
# Revision 2.32  2001/11/08 00:40:51  pwagner
# NOCLEANING avoids massive recompiles on changes to .configure
#
# Revision 2.31  2001/11/03 00:00:22  pwagner
# Added machines/mlsconfg to dependency search path
#
# Revision 2.30  2001/10/04 22:24:49  pwagner
# Limit suffix-based rules to (MLSCONFG)
#
# Revision 2.29  2001/09/26 17:24:59  pwagner
# Removed errors introduced with last changes; SIMPLIFYMAKE renamed NOHUNTING
#
# Revision 2.28  2001/09/25 20:27:15  pwagner
# New SIMPLIFYMAKE setting; other changes
#
# Revision 2.27  2001/08/15 21:01:05  pwagner
# Now excludes OMIT_SRCS from dependencies
#
# Revision 2.26  2001/07/03 17:18:21  pwagner
# Shield QuantityPVM from pvm-less lib
#
# Revision 2.25  2001/05/31 20:24:49  pwagner
# ghostbuster called with -lib option
#
# Revision 2.24  2001/05/14 23:58:43  pwagner
# default .configure name now set by environment
#
# Revision 2.23  2001/05/14 17:37:14  pwagner
# .configure now a variable file name
#
# Revision 2.22  2001/05/09 23:28:23  pwagner
# Optionally builds library w/o toolkit support
#
# Revision 2.21  2001/05/01 17:07:44  pwagner
# Now runs ghostbuster when updating
#
# Revision 2.20  2001/04/06 19:44:27  pwagner
# Fixed double-make problem
#
# Revision 2.19  2001/04/05 17:15:54  pwagner
# Added PROG to .PHONY list
#
# Revision 2.18  2001/04/02 17:38:50  pwagner
# make clean no longer causes bogus error messages in tests
#
# Revision 2.17  2001/03/12 22:46:57  pwagner
# Revised clean, tar; automatic symlink to Makefile
#
# Revision 2.16  2001/03/09 19:36:06  pwagner
# Added TAR=... option
#
# Revision 2.15  2001/03/06 18:17:57  pwagner
# clean more amitious--rms * from MLSCONFGl[13]/MakeFC
#
# Revision 2.14  2001/03/05 18:51:32  pwagner
# libmls.a again automatic target
#
# Revision 2.13  2001/02/09 18:38:20  pwagner
# Removed special treatment of Matrix_OBJS
#
# Revision 2.12  2001/02/08 23:16:46  pwagner
# Without PVM_OBJS.inc file
#
# Revision 2.11  2001/01/23 19:04:54  pwagner
# With optional PVM support
#
# Revision 2.10  2001/01/16 22:52:33  pwagner
# clean* and alltests go thru testdir
#
# Revision 2.9  2001/01/10 23:45:34  pwagner
# General housekeeping
#
# Revision 2.8  2000/11/22 21:48:20  pwagner
# Omitting MatrixModule_n from lib
#
# Revision 2.7  2000/11/21 00:40:05  pwagner
# machine.f90 out of lib/*.f90
#
# Revision 2.6  2000/11/17 22:09:30  pwagner
# conditional machine.f90 chmod
#
# Revision 2.5  2000/11/14 23:33:41  pwagner
# Excludes test*.f90
#
# Revision 2.4  2000/11/11 00:13:07  pwagner
# libmls.a now a phony target
#
# Revision 2.3  2000/11/02 23:24:28  pwagner
# clean_objs option now available
#
# Revision 2.2  2000/10/24 20:55:53  pwagner
# Added tar option
#
# Revision 2.1  2000/10/18 21:59:19  pwagner
# First commit
#

