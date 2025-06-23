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

short_name=l3

PROG=mlsl3d

# We store all the common Makefile things in one shared file
SRCLIB=../srclib
-include $(SRCLIB)/Makefile.h

# This will set our MLS platform, compiler. etc.
# w/o complaining if it does nor exist yet
-include $(CONFDIR)/$(MLSCFILE)

#Where else to look for USEd and INCLUDED files to add to dependency list
depends_where:=$(libmls_sources) $(srclib)
PROG_prereqs := $(CONFDIR)/$(MLSCFILE) $(PLATFORMS)/$(MLSCONFG) \
	  $(libmls_objs)/libmls.a $(REECHO) \
     $(INSTALLDIR)/libutctotai.a

#Any sources to explicitly exclude from dependency list?
#OMIT_SRCS := $(OMIT) -d synoptic_prog.f90
OMIT_SRCS := $(OMIT)

OMIT_PREFXD := $(shell ${REECHO} -prefix=-d $(OMIT_SRCS))

#These modules require custom builds, so omit build commands from dependencies
CUSTOM_BUILDS=L2Interface.f90 L3DMData.f90
CUSTOM_PREFXD := $(shell ${REECHO} -prefix=-db $(CUSTOM_BUILDS))
#----------------------- Build instructions

$(PROG): $(PROG_prereqs)
		@cd $(PLATFORMS)/$(MLSCONFG) && $(MAKE) $(PROG) \
      UNPREFXDPTHS="$(ADD2PATHS)" ; \
        return_status=`expr $$?` ; \
        if [ $$return_status != "0" ]; then \
		     exit 2 ; \
        elif [ "$(MARK_ALL_AS_UPTODATE)" != "no" ]; then \
		     $(MAKE) mark_all_as_uptodate \
           MARK_ALL_AS_UPTODATE= ;\
         fi

$(PLATFORMS)/$(MLSCONFG)/$(PROG): $(PROG)
		@cd $(PLATFORMS)/$(MLSCONFG) && $(MAKE) $(PROG) \
      UNPREFXDPTHS="$(ADD2PATHS)" ; \
        return_status=`expr $$?` ; \
        if [ $$return_status != "0" ]; then \
		     exit 2 ; \
        elif [ "$(MARK_ALL_AS_UPTODATE)" != "no" ]; then \
		     $(MAKE) mark_all_as_uptodate \
           MARK_ALL_AS_UPTODATE= ;\
         fi

.PHONY: clean configure depends depends--mod ghostbuster install \
  internalupdate mlsl3d mostlyclean \
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
   $(OMIT_PREFXD) $(UP_OR_LOW) \
   $(CUSTOM_PREFXD) $(depends_where)

ghostbuster:
	$(MLSBIN)/ghostbuster.sh -lib $(OMIT_PREFXD) \
         $(PLATFORMS)/$(MLSCONFG) ./ $(srclib)

update:
	@cd .. ; \
   $(MAKE) update SUBDIRS="blas lib l3" tests_sub_dirs=

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
	@echo "Updating the Makefile in $(SUBDIRS)"; \
    m_option="-machine" ; \
		$(MLSBIN)/mlsconfigure -dc "$(CONFDIR)" -pc "$(PLATFORMS)" \
         -cf "$(MLSCFILE)" -nc -nd $$m_option; \
      m_option=""

# install $(PROG) into $(INSTALLDIR)
install: $(PROG)
	@cd .. ; \
   $(MAKE) install LEVELS=l3

# install mlsl3d and mlsl3m into $(INSTALLDIR)
install-l3:
	@echo "Installing mlsl3d and mlsl3m into $(INSTALLDIR)"; \
	@cd .. && $(MAKE) install-l3

tar:
	rm -f ../$(PROG).tar
	@echo "tar-ing l3"; \
	$(TAR) hcvf ../$(PROG).tar * --exclude CVS\
	--exclude Make.dep\* --exclude NAG.\* --exclude LF95.\*

tarz: tar
	@echo "... then gzipping it"; \
	cd ..; rm -f $(PROG).tgz; gzip $(PROG).tar; mv $(PROG).tar.gz $(PROG).tgz

targ:
	rm -f ../$(PROG).tgz
	$(TAR) hcvzf ../$(PROG).tgz * --exclude \*.o \
   --exclude \*.mod --exclude $(PROG)

# $Log$
# Revision 1.56  2012/06/19 23:05:19  pwagner
# Uses parts from SRCLIB like the other Makefiles
#
# Revision 1.55  2010/08/27 20:56:39  pwagner
# Uses optional HDF*INC Makefile variables
#
# Revision 1.54  2008/07/11 23:59:19  pwagner
# 1st changes to get sunstudio to link mlsl2
#
# Revision 1.53  2008/04/22 17:20:37  pwagner
# CASCADE no longer affects where Makefile.dep goes
#
# Revision 1.52  2007/05/31 18:31:02  pwagner
# Added Absoft to list of tested compilers
#
# Revision 1.51  2007/04/13 20:15:18  pwagner
# Changes in line with Intel v9.1
#
# Revision 1.50  2006/04/01 01:14:02  pwagner
# Removed obsolete references to libhe5.a
#
# Revision 1.49  2005/06/23 23:02:30  pwagner
# Reworded Copyright statement
#
# Revision 1.48  2004/03/30 00:21:50  pwagner
# Take care that lib is rebuilt after changes in machines directory
#
# Revision 1.47  2003/06/23 22:51:16  pwagner
# Some fixes (I hope) to the prob of compiling in l3 twice from l3m
#
# Revision 1.46  2003/06/03 20:38:24  pwagner
# Renamed mkdir_cascade mkpath
#
# Revision 1.45  2003/05/21 20:52:13  pwagner
# Added PROG to .PHONY list as appropriate
#
# Revision 1.44  2003/05/20 18:39:23  pwagner
# Uses DONTBUILDPREQS to stop some hunting
#
# Revision 1.43  2003/05/19 23:43:06  pwagner
# Defined MLSHOME used by install command
#
# Revision 1.42  2003/04/30 23:04:55  pwagner
# Now lets you make update in each subdirectory
#
# Revision 1.41  2003/03/06 19:45:58  pwagner
# Added libutctotai as explicit prerequisite
#
# Revision 1.40  2003/02/24 21:16:01  pwagner
# Fixed bug where changes to .f9h were not triggering recompiles
#
# Revision 1.39  2002/08/22 23:33:33  pwagner
# Made ready for hdf5
#
# Revision 1.38  2002/07/18 18:00:02  pwagner
# NOHUNTING default when .mod-mediated dependency
#
# Revision 1.37  2002/07/15 16:52:59  pwagner
# Exits make with status 2 if error in any level
#
# Revision 1.36  2002/07/02 00:14:01  pwagner
# Skips marking uptodate after an error
#
# Revision 1.35  2002/06/21 23:24:09  pwagner
# More steps enroute to ending cascades of recompiles
#
# Revision 1.34  2002/06/21 00:15:28  pwagner
# Enroute to ending cascading recompiles
#
# Revision 1.33  2002/05/02 00:16:32  pwagner
# depends now creates Makefile.nodep, too
#
# Revision 1.32  2002/01/29 00:46:47  pwagner
# Builds he5lib when HDFVERSIONS defined
#
# Revision 1.31  2002/01/15 00:40:18  pwagner
# Now uses user-settable OMIT and internal OMIT_SRCS and OMIT_PREFIXD variables
#
# Revision 1.30  2001/11/08 00:40:51  pwagner
# NOCLEANING avoids massive recompiles on changes to .configure
#
# Revision 1.29  2001/10/04 22:25:36  pwagner
# Limit suffix-based rules to (MLSCONFG); add lib/lit_names.txt to dependency list
#
# Revision 1.28  2001/09/26 17:24:59  pwagner
# Removed errors introduced with last changes; SIMPLIFYMAKE renamed NOHUNTING
#
# Revision 1.27  2001/09/25 20:27:15  pwagner
# New SIMPLIFYMAKE setting; other changes
#
# Revision 1.26  2001/08/22 22:46:13  pwagner
# Made reecho.sh a prerequisite where appropriate
#
# Revision 1.25  2001/08/15 21:01:04  pwagner
# Now excludes OMIT_SRCS from dependencies
#
# Revision 1.24  2001/08/13 18:09:08  pwagner
# replaced shell echo -> shell {REECHO}
#
# Revision 1.23  2001/07/19 22:48:52  pwagner
# Builds libmlsl3.a; renamed prog to mlsl3d
#
# Revision 1.22  2001/05/31 20:19:31  pwagner
# Added ghostbuster to .PHONY list
#
# Revision 1.21  2001/05/24 17:56:39  pwagner
# Expanded LIB_prereqs
#
# Revision 1.20  2001/05/17 20:02:59  pwagner
# Automatically rebuilds lib on change to machine.f90
#
# Revision 1.19  2001/05/14 23:58:42  pwagner
# default .configure name now set by environment
#
# Revision 1.18  2001/05/14 17:37:14  pwagner
# .configure now a variable file name
#
# Revision 1.17  2001/05/01 18:32:13  pwagner
# Added .configure to prereqs for lib
#
# Revision 1.16  2001/05/01 17:07:44  pwagner
# Now runs ghostbuster when updating
#
# Revision 1.15  2001/04/12 22:03:08  pwagner
# Changes to lib or fwdmdl trigger rebuilds even on suns
#
# Revision 1.14  2001/04/05 17:15:54  pwagner
# Added PROG to .PHONY list
#
# Revision 1.13  2001/03/13 19:07:59  pwagner
# clean now removes executables, too
#
# Revision 1.12  2001/03/12 22:46:57  pwagner
# Revised clean, tar; automatic symlink to Makefile
#
# Revision 1.11  2001/03/09 19:36:06  pwagner
# Added TAR=... option
#
# Revision 1.10  2001/03/06 18:17:57  pwagner
# clean more amitious--rms * from MLSCONFGl[13]/MakeFC
#
# Revision 1.9  2001/03/01 18:20:46  pwagner
# Hide case statement
#
# Revision 1.8  2001/02/28 23:33:20  pwagner
# Build only if FFTW_ROOT defined
#
# Revision 1.7  2001/02/22 18:08:11  pwagner
# Includes L3DZData.f90
#
# Revision 1.6  2001/02/15 18:55:12  pwagner
# L3DZData.f90 not ready for compiling
#
# Revision 1.5  2001/02/08 23:16:46  pwagner
# Without PVM_OBJS.inc file
#
# Revision 1.4  2001/01/23 23:27:38  pwagner
# With optional PVM support
#
# Revision 1.3  2001/01/10 23:45:34  pwagner
# General housekeeping
#
# Revision 1.2  2000/11/17 20:56:31  pwagner
# Using srclib_f90 as preq to lib
#
# Revision 1.1  2000/11/15 18:02:58  pwagner
# First commit
#
#


