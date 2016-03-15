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

short_name=l1

PROG=mlsl1log mlsl1g mlsl1t mlsl0sn

# We store all the common Makefile things in one shared file
SRCLIB=../srclib
-include $(SRCLIB)/Makefile.h

# This will set our MLS platform, compiler. etc.
# w/o complaining if it does nor exist yet
-include $(CONFDIR)/$(MLSCFILE)

PROG_prereqs := $(CONFDIR)/$(MLSCFILE) \
   $(PLATFORMS)/$(MLSCONFG) $(libmls_objs)/libmls.a \
   $(INSTALLDIR)/libutctotai.a

# Any sources to explicitly exclude from dependency list?  Include any
# source files that have a 'program' statement in them in this line. If
# you don't do this your first l1 link will fail with the opaque error
# message that there are multiple "'mains' defined"
OMIT_SRCS := $(OMIT) MLSL1.f90 test_streq.f90

OMIT_PREFXD := $(shell ${REECHO} -prefix=-d $(OMIT_SRCS))

#These are the source files for the drivers; they should be excluded from OBJS
DRIVERS=MLSL1G.f90 MLSL1log.f90 MLSL1T.f90 MLSL0SN.f90
DRIVERS_PREFXD := $(shell ${REECHO} -prefix=-do $(DRIVERS))

#These modules require custom builds, so omit build commands from dependencies
CUSTOM_BUILDS=OutputL1B.f90 OutputL1B_HDF4.f90 EngTbls.f90 \
  CalibWeightsFlags.f90 L1LogUtils.f90 
CUSTOM_PREFXD := $(shell ${REECHO} -prefix=-db $(CUSTOM_BUILDS))

ifdef HDF5INC
   ADD2PATHS := $(ADD2PATHS) $(HDF5INC)
else
   ADD2PATHS := $(ADD2PATHS) $(HDF5)/../fortran/src $(HDF5)/../include
endif

#----------------------- Build instructions

$(PROG):  $(PROG_prereqs)
	if [ ! -f Calibration.f9h ]; then \
          pwd; \
          ls Calibration*.f9h;\
          cp Calibration-normal.f9h Calibration.f9h; \
        fi; \
	#grep WinMAFs Calibration.f9h
	cat Calibration.f9h
	@cd $(PLATFORMS)/$(MLSCONFG) && $(MAKE) $(PROG) \
        UNPREFXDPTHS="$(ADD2PATHS)" ;\
        return_status=`expr $$?` ; \
        if [ $$return_status != "0" ]; then \
		     exit 2 ; \
        elif [ "$(MARK_ALL_AS_UPTODATE)" != "no" ]; then \
		     $(MAKE) mark_all_as_uptodate \
           UNPREFXDPTHS="$(ADD2PATHS)" MARK_ALL_AS_UPTODATE= ;\
         fi

$(PLATFORMS)/$(MLSCONFG)/mlsl0sn \
$(PLATFORMS)/$(MLSCONFG)/mlsl1log \
$(PLATFORMS)/$(MLSCONFG)/mlslg \
$(PLATFORMS)/$(MLSCONFG)/mlsl1t : $(PROG)
	@cd $(PLATFORMS)/$(MLSCONFG) && $(MAKE) $(PROG) \
        UNPREFXDPTHS="$(ADD2PATHS)" ;\
        return_status=`expr $$?` ; \
        if [ $$return_status != "0" ]; then \
		     exit 2 ; \
        elif [ "$(MARK_ALL_AS_UPTODATE)" != "no" ]; then \
		     $(MAKE) mark_all_as_uptodate \
           UNPREFXDPTHS="$(ADD2PATHS)" MARK_ALL_AS_UPTODATE= ;\
         fi

.PHONY: clean configure depends depends--mod ghostbuster install \
  internalupdate mlsl0sn mlsl1log mlsl1g mlsl1t mostlyclean \
  tar targ tarz update

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
          -c "$(MLSF95)" -p "$(MLSPLAT)" -cf "$(MLSCFILE)" -nc \
          $(RECFG)

depends:
	$(MLSBIN)/makemakedep.sh -p $(PLATFORMS)/$(MLSCONFG)/Makefile.dep \
    $(OMIT_PREFXD) $(DRIVERS_PREFXD) $(UP_OR_LOW) $(CUSTOM_PREFXD) \
	    $(libmls_sources) $(srclib)

ghostbuster:
	$(MLSBIN)/ghostbuster.sh $(OMIT_PREFXD) \
       $(PLATFORMS)/$(MLSCONFG) ./ $(srclib)

update:
	@cd .. ; \
   $(MAKE) update SUBDIRS="blas lib l1" tests_sub_dirs=

mostlyclean:
	@echo "Cleaning directory $(MLSCONFG)" ; \
	if [ "$(MLSCONFG)" != "" ] ; then  \
		echo "Deleting $(MLSCONFG) files" ; \
		rm -f $(PLATFORMS)/$(MLSCONFG)/* ; \
	fi ; \
	if [ -d "$(PLATFORMS)/$(MLSCONFG)" ] ; then  \
		echo "Deleting $(MLSCONFG) directory" ; \
		rmdir $(PLATFORMS)/$(MLSCONFG) ; \
	fi
	rm -f Make.dep.*

clean:
	@echo "Cleaning objects from directory $(MLSCONFG)" ; \
	if [ "$(MLSCONFG)" != "" ] ; then  \
		rm -f $(PLATFORMS)/$(MLSCONFG)/*.[ao] ; \
		rm -f $(PLATFORMS)/$(MLSCONFG)/*.mod ; \
		rm -f $(PLATFORMS)/$(MLSCONFG)/$(PROG) ; \
	fi
	rm -f Make.dep.*

internalupdate: 
	@echo "Updating the Makefile in ./"; \
    m_option="-machine" ; \
		$(MLSBIN)/mlsconfigure -dc "$(CONFDIR)" -pc "$(PLATFORMS)" \
         -cf "$(MLSCFILE)" -nc -nd $$m_option; \
      m_option=""

# install $(PROG) into $(INSTALLDIR)
install: $(PROG)
	@cd .. ; \
   $(MAKE) install LEVELS=l1

tar:
	rm -f ../$(PROG).tar
	@echo "tar-ing l1"; \
	rm -f ../$(PROG).tar
	$(TAR) hcvf ../$(PROG).tar * --exclude CVS\
	--exclude Make.dep\* --exclude NAG.\* --exclude LF95.\*

tarz: tar
	@echo "... then gzipping it"; \
	cd ..; rm -f $(PROG).tgz; gzip $(PROG).tar; mv $(PROG).tar.gz $(PROG).tgz

targ:
	rm -f ../$(PROG).tgz
	$(TAR) hcvzf ../$(PROG).tgz * --exclude \*.o --exclude \*.mod --exclude mlsl1

# $Log$
# Revision 2.10  2015/08/14 20:06:51  whdaffer
# A comment that I think explains the issue with new files with
# `program' statements failure
#
# Revision 2.9  2015/08/10 23:15:57  pwagner
# *** empty log message ***
#
# Revision 2.60  2015/04/23 23:04:37  pwagner
# Should fix problem with first checkout before Calibration.f9h initialized
#
# Revision 2.59  2015/04/23 17:46:27  whdaffer
# removed Makefile
#
# Revision 2.58  2010/10/12 22:24:29  pwagner
# Uses Makefile.h from srclib
#
# Revision 2.57  2010/08/27 20:56:39  pwagner
# Uses optional HDF*INC Makefile variables
#
# Revision 2.56  2008/07/11 23:59:19  pwagner
# 1st changes to get sunstudio to link mlsl2
#
# Revision 2.55  2008/04/22 17:20:21  pwagner
# CASCADE no longer affects where Makefile.dep goes
#
# Revision 2.54  2007/05/31 18:31:02  pwagner
# Added Absoft to list of tested compilers
#
# Revision 2.53  2007/04/13 20:15:18  pwagner
# Changes in line with Intel v9.1
#
# Revision 2.52  2006/04/01 01:11:58  pwagner
# Moved installation commands to ../MakeFC
#
# Revision 2.51  2006/03/28 21:54:33  pwagner
# Must turn off optimization to work around Lahey bug
#
# Revision 2.50  2006/03/27 19:18:30  pwagner
# Added new pge mlsl0sn to store switch settings
#
# Revision 2.49  2005/08/25 19:21:58  pwagner
# Removed MLSL1BOA as target
#
# Revision 2.48  2005/06/23 23:02:30  pwagner
# Reworded Copyright statement
#
# Revision 2.47  2004/03/30 00:21:50  pwagner
# Take care that lib is rebuilt after changes in machines directory
#
# Revision 2.46  2003/10/24 20:59:21  pwagner
# Added mlsl1boa utility as new target
#
# Revision 2.45  2003/06/03 20:38:24  pwagner
# Renamed mkdir_cascade mkpath
#
# Revision 2.44  2003/05/21 20:52:13  pwagner
# Added PROG to .PHONY list as appropriate
#
# Revision 2.43  2003/05/20 18:39:23  pwagner
# Uses DONTBUILDPREQS to stop some hunting
#
# Revision 2.42  2003/05/19 23:43:06  pwagner
# Defined MLSHOME used by install command
#
# Revision 2.41  2003/04/30 23:04:55  pwagner
# Now lets you make update in each subdirectory
#
# Revision 2.40  2003/03/06 19:45:58  pwagner
# Added libutctotai as explicit prerequisite
#
# Revision 2.39  2003/02/24 21:16:01  pwagner
# Fixed bug where changes to .f9h were not triggering recompiles
#
# Revision 2.38  2003/01/31 23:52:42  pwagner
# Changed for new l1 with 3 executables
#
# Revision 2.37  2002/11/18 19:26:16  pwagner
# Added Output_HDF4.f90 to CUSTOM_BUILDS
#
# Revision 2.36  2002/07/18 18:00:01  pwagner
# NOHUNTING default when .mod-mediated dependency
#
# Revision 2.35  2002/07/16 17:18:31  pwagner
# UNPREFXDPTHS adds hdf5 dirs for f95 to search
#
# Revision 2.34  2002/07/15 16:52:59  pwagner
# Exits make with status 2 if error in any level
#
# Revision 2.33  2002/07/02 00:14:01  pwagner
# Skips marking uptodate after an error
#
# Revision 2.32  2002/06/21 23:24:09  pwagner
# More steps enroute to ending cascades of recompiles
#
# Revision 2.31  2002/06/21 00:15:27  pwagner
# Enroute to ending cascading recompiles
#
# Revision 2.30  2002/06/19 19:39:46  pwagner
# Allows use of .mod and .f9h files from hdf5
#
# Revision 2.29  2002/06/18 18:53:35  pwagner
# Modified to correctly use he5lib modules
#
# Revision 2.28  2002/05/02 00:16:32  pwagner
# depends now creates Makefile.nodep, too
#
# Revision 2.27  2002/01/29 00:46:47  pwagner
# Builds he5lib when HDFVERSIONS defined
#
# Revision 2.26  2002/01/15 00:40:18  pwagner
# Now uses user-settable OMIT and internal OMIT_SRCS and OMIT_PREFIXD variables
#
# Revision 2.25  2001/11/08 00:40:51  pwagner
# NOCLEANING avoids massive recompiles on changes to .configure
#
# Revision 2.24  2001/10/04 22:25:36  pwagner
# Limit suffix-based rules to (MLSCONFG); add lib/lit_names.txt to dependency list
#
# Revision 2.23  2001/09/26 17:24:59  pwagner
# Removed errors introduced with last changes; SIMPLIFYMAKE renamed NOHUNTING
#
# Revision 2.22  2001/09/25 20:27:15  pwagner
# New SIMPLIFYMAKE setting; other changes
#
# Revision 2.21  2001/08/22 22:46:13  pwagner
# Made reecho.sh a prerequisite where appropriate
#
# Revision 2.20  2001/08/15 21:01:04  pwagner
# Now excludes OMIT_SRCS from dependencies
#
# Revision 2.19  2001/08/13 18:09:08  pwagner
# replaced shell echo -> shell {REECHO}
#
# Revision 2.18  2001/05/31 20:19:31  pwagner
# Added ghostbuster to .PHONY list
#
# Revision 2.17  2001/05/24 17:56:39  pwagner
# Expanded LIB_prereqs
#
# Revision 2.16  2001/05/17 20:02:59  pwagner
# Automatically rebuilds lib on change to machine.f90
#
# Revision 2.15  2001/05/14 23:58:42  pwagner
# default .configure name now set by environment
#
# Revision 2.14  2001/05/14 17:37:14  pwagner
# .configure now a variable file name
#
# Revision 2.13  2001/05/01 18:32:13  pwagner
# Added .configure to prereqs for lib
#
# Revision 2.12  2001/05/01 17:07:43  pwagner
# Now runs ghostbuster when updating
#
# Revision 2.11  2001/04/12 22:03:08  pwagner
# Changes to lib or fwdmdl trigger rebuilds even on suns
#
# Revision 2.10  2001/04/05 17:15:54  pwagner
# Added PROG to .PHONY list
#
# Revision 2.9  2001/03/13 19:07:59  pwagner
# clean now removes executables, too
#
# Revision 2.8  2001/03/12 22:46:57  pwagner
# Revised clean, tar; automatic symlink to Makefile
#
# Revision 2.7  2001/03/09 19:36:06  pwagner
# Added TAR=... option
#
# Revision 2.6  2001/03/06 18:17:57  pwagner
# clean more amitious--rms * from MLSCONFGl[13]/MakeFC
#
# Revision 2.5  2001/02/08 23:16:46  pwagner
# Without PVM_OBJS.inc file
#
# Revision 2.4  2001/01/23 23:27:38  pwagner
# With optional PVM support
#
# Revision 2.3  2001/01/10 23:45:33  pwagner
# General housekeeping
#
# Revision 2.2  2000/11/17 20:54:34  pwagner
# Using srclib_f90 as preq to lib
#
# Revision 2.1  2000/11/15 00:11:27  pwagner
# First commit
#
#


