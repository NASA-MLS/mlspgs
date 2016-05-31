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

short_name=l2

PROG=mlsl2

# We store all the common Makefile things in one shared file
SRCLIB=../srclib
-include $(SRCLIB)/Makefile.h

# This will set our MLS platform, compiler. etc.
# w/o complaining if it does nor exist yet
-include $(CONFDIR)/$(MLSCFILE)

depends_where:=$(libmls_sources) $(SRCLIB) \
   $(libfwdmdl_sources) $(libcloud_sources)
PROG_prereqs := $(CONFDIR)/$(MLSCFILE) $(PLATFORMS)/$(MLSCONFG) \
	  $(libblas_objs)/libmlspack.a \
	  $(libmls_objs)/libmls.a \
	  $(libfwdmdl_objs)/libfwdmdl.a \
	  $(libcloud_objs)/libcloud.a $(REECHO) \
     $(INSTALLDIR)/libutctotai.a


#Any sources to explicitly exclude from dependency list?
OMIT_SRCS := $(OMIT) TestCholesky_m.f90 Test_Parse_Signals_m.f90

OMIT_PREFXD := $(shell ${REECHO} -prefix=-d $(OMIT_SRCS))

#These modules require custom builds, so omit build commands from dependencies
#CUSTOM_BUILDS=MLSL2.f90 init_tables_module.f90
CUSTOM_BUILDS=MLSL2.f90
CUSTOM_PREFXD := $(shell ${REECHO} -prefix=-db $(CUSTOM_BUILDS))

#----------------------- Build instructions

$(PROG): $(PROG_prereqs)
	$(MAKE) the_program
	if [ -f $(CONFDIR)/lib/$(MLSCONFG)/new_lit_names.txt ]; then \
	  rm -f $(CONFDIR)/l2/$(MLSCONFG)/mlsl2 \
	        $(CONFDIR)/lib/$(MLSCONFG)/new_lit_names.txt; \
	  $(MAKE) the_program; \
	fi

# For reasons not yet understood, executables built
# directly by make did not always
# run properly, sometimes generating errors during tree_checker
#mlsl2-delayed: $(PROG_prereqs)
#	$(MAKE) the_program
#	/bin/rm -f $(PLATFORMS)/$(MLSCONFG)/mlsl2

the_program: $(PROG_prereqs)
#$(PROG): $(PROG_prereqs)
	@cd $(PLATFORMS)/$(MLSCONFG) && $(MAKE) $(PROG) \
        UNPREFXDPTHS="$(ADD2PATHS)" MARK_ALL_AS_UPTODATE=no;\
        return_status=`expr $$?` ; \
        if [ $$return_status != "0" ]; then \
		     exit 2 ; \
        elif [ "$(MARK_ALL_AS_UPTODATE)" != "no" ]; then \
		     $(MAKE) mark_all_as_uptodate \
           UNPREFXDPTHS="$(ADD2PATHS)" MARK_ALL_AS_UPTODATE=yes ;\
        fi

.DELETE_ON_ERROR:

.PHONY: clean configure depends depends--mod ghostbuster install \
  internalupdate  mlsl2 mostlyclean \
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
	$(MLSBIN)/ghostbuster.sh $(OMIT_PREFXD) \
       $(PLATFORMS)/$(MLSCONFG) ./ $(SRCLIB)

update:
	@cd .. ; \
   $(MAKE) update SUBDIRS="blas lib fwdmdl cloudfwdm l2" tests_sub_dirs= JUSTSUBDIRS=yes

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

install:
	@cd .. ; \
   $(MAKE) install LEVELS=l2
tar:
	rm -f ../$(PROG).tar
	@echo "tar-ing l2"; \
	$(TAR) hcvf ../$(PROG).tar * --exclude CVS\
	--exclude Make.dep\* --exclude NAG.\* --exclude LF95.\*

tarz: tar
	@echo "... then gzipping it"; \
	cd ..; rm -f $(PROG).tgz; gzip $(PROG).tar; mv $(PROG).tar.gz $(PROG).tgz

targ:
	rm -f ../$(PROG).tgz
	$(TAR) hcvzf ../$(PROG).tgz * --exclude \*.o --exclude \*.mod --exclude mlsl2

# $Log$
# Revision 2.5  2016/05/31 23:42:03  pwagner
# Tried again to fix need to build mlsl2 twice
#
# Revision 2.4  2015/08/10 23:19:08  pwagner
# After switching away from MakeFC
#
# Revision 2.76  2014/01/29 21:00:27  pwagner
# Shorten make update
#
# Revision 2.75  2011/11/18 23:39:31  pwagner
# Removed double-link for mlsl2; corrected spelling error
#
# Revision 2.74  2011/09/02 18:35:36  pwagner
# Double-make workaround for tree_checker issue moved here
#
# Revision 2.73  2011/05/26 20:40:48  pwagner
# Small tweaks; do they fix anything?
#
# Revision 2.72  2010/10/12 22:24:29  pwagner
# Uses Makefile.h from srclib
#
# Revision 2.71  2010/08/27 20:56:39  pwagner
# Uses optional HDF*INC Makefile variables
#
# Revision 2.70  2008/07/11 23:59:19  pwagner
# 1st changes to get sunstudio to link mlsl2
#
# Revision 2.69  2008/04/22 17:18:03  pwagner
# CASCADE no longer affects where Makefile.dep goes
#
# Revision 2.68  2007/06/08 22:01:17  vsnyder
# Stuff for .F90 files
#
# Revision 2.67  2007/05/31 18:31:02  pwagner
# Added Absoft to list of tested compilers
#
# Revision 2.66  2007/04/13 20:15:18  pwagner
# Changes in line with Intel v9.1
#
# Revision 2.65  2006/04/01 01:14:02  pwagner
# Removed obsolete references to libhe5.a
#
# Revision 2.64  2005/06/23 23:02:30  pwagner
# Reworded Copyright statement
#
# Revision 2.63  2004/03/30 00:21:50  pwagner
# Take care that lib is rebuilt after changes in machines directory
#
# Revision 2.62  2003/08/01 20:01:41  pwagner
# (Re)Uses install target at root level to install
#
# Revision 2.61  2003/06/03 20:38:24  pwagner
# Renamed mkdir_cascade mkpath
#
# Revision 2.60  2003/05/21 20:52:13  pwagner
# Added PROG to .PHONY list as appropriate
#
# Revision 2.59  2003/05/20 18:39:23  pwagner
# Uses DONTBUILDPREQS to stop some hunting
#
# Revision 2.58  2003/05/19 23:43:06  pwagner
# Defined MLSHOME used by install command
#
# Revision 2.57  2003/04/30 23:04:55  pwagner
# Now lets you make update in each subdirectory
#
# Revision 2.56  2003/02/24 21:16:01  pwagner
# Fixed bug where changes to .f9h were not triggering recompiles
#
# Revision 2.55  2003/02/21 21:02:40  pwagner
# Eliminated TestCholesky_m.f90 Test_Parse_Signals_m.f90 from srcs
#
# Revision 2.54  2002/07/18 18:00:01  pwagner
# NOHUNTING default when .mod-mediated dependency
#
# Revision 2.53  2002/07/15 16:52:59  pwagner
# Exits make with status 2 if error in any level
#
# Revision 2.52  2002/07/12 20:34:54  pwagner
# init_tables_module does not require custom build
#
# Revision 2.51  2002/07/11 22:12:10  pwagner
# Changes to support hdf5/hdfeos5
#
# Revision 2.50  2002/07/02 00:14:01  pwagner
# Skips marking uptodate after an error
#
# Revision 2.49  2002/07/01 17:31:19  pwagner
# No longer compiles some sources twice
#
# Revision 2.48  2002/06/21 23:24:09  pwagner
# More steps enroute to ending cascades of recompiles
#
# Revision 2.47  2002/06/21 00:15:27  pwagner
# Enroute to ending cascading recompiles
#
# Revision 2.46  2002/05/02 00:16:32  pwagner
# depends now creates Makefile.nodep, too
#
# Revision 2.45  2002/04/29 17:43:29  pwagner
# Directs build of utctotai library in INSTALLDIR
#
# Revision 2.44  2002/01/31 00:34:33  pwagner
# Tries to get .mod files from he5lib
#
# Revision 2.43  2002/01/29 00:46:47  pwagner
# Builds he5lib when HDFVERSIONS defined
#
# Revision 2.42  2002/01/15 00:40:18  pwagner
# Now uses user-settable OMIT and internal OMIT_SRCS and OMIT_PREFIXD variables
#
# Revision 2.41  2001/11/15 00:10:13  pwagner
# Rebuilds on changes to mlspgs/blas
#
# Revision 2.40  2001/11/08 00:40:51  pwagner
# NOCLEANING avoids massive recompiles on changes to .configure
#
# Revision 2.39  2001/10/04 22:25:36  pwagner
# Limit suffix-based rules to (MLSCONFG); add lib/lit_names.txt to dependency list
#
# Revision 2.38  2001/09/26 17:24:59  pwagner
# Removed errors introduced with last changes; SIMPLIFYMAKE renamed NOHUNTING
#
# Revision 2.37  2001/09/25 20:27:15  pwagner
# New SIMPLIFYMAKE setting; other changes
#
# Revision 2.36  2001/08/22 22:46:13  pwagner
# Made reecho.sh a prerequisite where appropriate
#
# Revision 2.35  2001/08/15 21:01:04  pwagner
# Now excludes OMIT_SRCS from dependencies
#
# Revision 2.34  2001/08/13 18:09:08  pwagner
# replaced shell echo -> shell {REECHO}
#
# Revision 2.33  2001/07/25 23:19:58  pwagner
# Removed unnecessary lib...a from target prereqs
#
# Revision 2.32  2001/07/18 00:12:55  pwagner
# No longer makes libfwdmdl directly
#
# Revision 2.31  2001/05/31 23:27:51  pwagner
# Ended role of l2_NEEDS_fwdmdl; cloudfwdm prereq of l2
#
# Revision 2.30  2001/05/31 20:19:31  pwagner
# Added ghostbuster to .PHONY list
#
# Revision 2.29  2001/05/23 22:27:41  pwagner
# Expanded LIB_prereqs
#
# Revision 2.28  2001/05/17 20:02:59  pwagner
# Automatically rebuilds lib on change to machine.f90
#
# Revision 2.27  2001/05/14 23:58:42  pwagner
# default .configure name now set by environment
#
# Revision 2.26  2001/05/14 17:37:14  pwagner
# .configure now a variable file name
#
# Revision 2.25  2001/05/01 18:32:56  pwagner
# Added .configure to prereqs for lib, fwdmdl
#
# Revision 2.24  2001/05/01 17:07:44  pwagner
# Now runs ghostbuster when updating
#
# Revision 2.23  2001/04/12 17:17:27  pwagner
# Changes to lib or fwdmdl trigger rebuilds even on suns
#
# Revision 2.22  2001/04/06 19:45:17  pwagner
# Fixed prob of lib changes not triggering rebuild fwdmdl
#
# Revision 2.21  2001/04/05 17:15:54  pwagner
# Added PROG to .PHONY list
#
# Revision 2.20  2001/04/02 23:46:17  pwagner
# Now even LF95 should use write_meta
#
# Revision 2.19  2001/04/02 17:42:08  pwagner
# Builds mlsl2 only if PVM_ROOT defined
#
# Revision 2.18  2001/03/14 17:19:37  pwagner
# Omits write_meta.f90 if Lahey
#
# Revision 2.17  2001/03/13 19:07:59  pwagner
# clean now removes executables, too
#
# Revision 2.16  2001/03/12 22:46:57  pwagner
# Revised clean, tar; automatic symlink to Makefile
#
# Revision 2.15  2001/03/09 19:36:06  pwagner
# Added TAR=... option
#
# Revision 2.14  2001/03/08 00:53:21  pwagner
# l2_NEEDS_fwdmdl now set by default
#
# Revision 2.13  2001/03/07 23:27:24  pwagner
# Adds SnoopMLSL2 to Makefile.dep if PVM_ROOT defined
#
# Revision 2.12  2001/03/06 18:17:57  pwagner
# clean more amitious--rms * from MLSCONFGl[13]/MakeFC
#
# Revision 2.11  2001/03/05 21:48:50  pwagner
# fwdmdl searched for dependencies if l2_NEEDS_fwdmdl
#
# Revision 2.10  2001/02/12 17:43:13  pwagner
# RetrievalModule and ForwardModelInterface now included
#
# Revision 2.9  2001/02/08 23:16:46  pwagner
# Without PVM_OBJS.inc file
#
# Revision 2.8  2001/01/23 19:06:25  pwagner
# With optional PVM support
#
# Revision 2.7  2001/01/10 23:45:34  pwagner
# General housekeeping
#
# Revision 2.6  2000/11/18 00:11:06  pwagner
# No target with *.f90 to confuse sun
#
# Revision 2.5  2000/11/13 17:52:23  pwagner
# Forces rebuild of lib even from l2
#
# Revision 2.4  2000/11/03 18:16:43  pwagner
# Shares sources in srclib
#
# Revision 2.3  2000/11/02 23:26:09  pwagner
# clean_objs option now available
#
# Revision 2.2  2000/10/24 21:07:43  pwagner
# Added tar option
#
# Revision 2.1  2000/10/18 21:58:16  pwagner
# First commit
#


