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

short_name=cloudfwdm

PROG=libcloud.a

# We store all the common Makefile things in one shared file
SRCLIB=../srclib
-include $(SRCLIB)/Makefile.h

# This will set our MLS platform, compiler. etc.
# w/o complaining if it does nor exist yet
-include $(CONFDIR)/$(MLSCFILE)

#Where else to look for USEd and INCLUDED files to add to dependency list
depends_where:=$(libmls_sources) $(srclib) $(libfwdmdl_sources)

#Any sources to explicitly exclude from dependency list?
OMIT_SRCS := $(OMIT)

PROG_prereqs=$(PLATFORMS)/$(MLSCONFG) $(libfwdmdl_objs)/libfwdmdl.a

OMIT_PREFXD := $(shell ${REECHO} -prefix=-d $(OMIT_SRCS))

#These modules require custom builds, so omit build commands from dependencies
CUSTOM_BUILDS=
CUSTOM_PREFXD := $(shell ${REECHO} -prefix=-db $(CUSTOM_BUILDS))

ifdef HDF5INC
   ADD2PATHS := $(ADD2PATHS) $(HDF5INC)
else
   ADD2PATHS := $(ADD2PATHS) $(HDF5)/../fortran/src $(HDF5)/../include
endif

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

.PHONY: clean configure depends depends--mod ghostbuster mostlyclean \
  tar targ tarz update libcloud.a

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
   $(MAKE) update SUBDIRS="blas lib fwdmdl cloudfwdm" tests_sub_dirs=

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

tar:
	rm -f ../$(PROG).tar
	@echo "tar-ing $(PROG)"; \
	$(TAR) hcvf ../$(PROG).tar * --exclude CVS\
	--exclude Make.dep\* --exclude NAG.\* --exclude LF95.\*

tarz: tar
	@echo "... then gzipping it"; \
	cd ..; rm -f $(PROG).tgz; gzip $(PROG).tar; mv $(PROG).tar.gz $(PROG).tgz

targ:
	rm -f ../$(PROG).tgz
	$(TAR) hcvzf ../$(PROG).tgz * --exclude \*.o --exclude \*.mod --exclude mlsl2

# $Log$
# Revision 1.33  2010/10/12 22:25:03  pwagner
# Uses Makefile.h from srclib
#
# Revision 1.32  2010/08/27 20:56:39  pwagner
# Uses optional HDF*INC Makefile variables
#
# Revision 1.31  2008/07/11 23:59:19  pwagner
# 1st changes to get sunstudio to link mlsl2
#
# Revision 1.30  2008/04/22 17:20:01  pwagner
# CASCADE no longer affects where Makefile.dep goes
#
# Revision 1.29  2007/05/31 18:31:02  pwagner
# Added Absoft to list of tested compilers
#
# Revision 1.28  2007/04/13 20:15:18  pwagner
# Changes in line with Intel v9.1
#
# Revision 1.27  2005/06/23 23:02:30  pwagner
# Reworded Copyright statement
#
# Revision 1.26  2004/03/30 00:21:50  pwagner
# Take care that lib is rebuilt after changes in machines directory
#
# Revision 1.25  2003/05/21 20:52:13  pwagner
# Added PROG to .PHONY list as appropriate
#
# Revision 1.24  2003/05/20 18:39:21  pwagner
# Uses DONTBUILDPREQS to stop some hunting
#
# Revision 1.23  2003/04/30 23:04:55  pwagner
# Now lets you make update in each subdirectory
#
# Revision 1.22  2003/02/24 21:16:01  pwagner
# Fixed bug where changes to .f9h were not triggering recompiles
#
# Revision 1.21  2002/07/18 18:00:01  pwagner
# NOHUNTING default when .mod-mediated dependency
#
# Revision 1.20  2002/07/15 16:52:58  pwagner
# Exits make with status 2 if error in any level
#
# Revision 1.19  2002/07/11 22:12:10  pwagner
# Changes to support hdf5/hdfeos5
#
# Revision 1.18  2002/07/02 00:14:00  pwagner
# Skips marking uptodate after an error
#
# Revision 1.17  2002/06/21 23:24:08  pwagner
# More steps enroute to ending cascades of recompiles
#
# Revision 1.16  2002/06/21 00:15:27  pwagner
# Enroute to ending cascading recompiles
#
# Revision 1.15  2002/05/02 00:16:32  pwagner
# depends now creates Makefile.nodep, too
#
# Revision 1.14  2002/01/31 00:34:33  pwagner
# Tries to get .mod files from he5lib
#
# Revision 1.13  2002/01/15 00:40:18  pwagner
# Now uses user-settable OMIT and internal OMIT_SRCS and OMIT_PREFIXD variables
#
# Revision 1.12  2001/11/08 00:40:51  pwagner
# NOCLEANING avoids massive recompiles on changes to .configure
#
# Revision 1.11  2001/10/04 22:25:36  pwagner
# Limit suffix-based rules to (MLSCONFG); add lib/lit_names.txt to dependency list
#
# Revision 1.10  2001/09/26 17:24:58  pwagner
# Removed errors introduced with last changes; SIMPLIFYMAKE renamed NOHUNTING
#
# Revision 1.9  2001/09/25 20:27:15  pwagner
# New SIMPLIFYMAKE setting; other changes
#
# Revision 1.8  2001/09/07 22:23:52  pwagner
# Triggers build on changes in lib
#
# Revision 1.7  2001/08/22 22:46:13  pwagner
# Made reecho.sh a prerequisite where appropriate
#
# Revision 1.6  2001/08/15 21:01:04  pwagner
# Now excludes OMIT_SRCS from dependencies
#
# Revision 1.5  2001/08/13 18:09:08  pwagner
# replaced shell echo -> shell {REECHO}
#
# Revision 1.4  2001/08/07 23:38:50  pwagner
# Did away with -f77 option to dependency calculation
#
# Revision 1.3  2001/07/25 23:19:58  pwagner
# Removed unnecessary lib...a from target prereqs
#
# Revision 1.2  2001/05/31 20:24:49  pwagner
# ghostbuster called with -lib option
#
# Revision 1.1  2001/05/30 16:42:53  pwagner
# First commit
#
#


