#	Makefile
#
#  Lost? Type 'cd ../; make help -f MakeFC' from this directory
#
# Copyright (c) 2003, California Institute of Technology.  ALL RIGHTS RESERVED.
# U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

# "$Id$"

# ---------------------- The name of this file (is the right-hand-side below)
MakeFName = MakeFC

# where to store .configure
CONFDIR=..

ifdef MLSINSTALLDIR
   INSTALLDIR=$(MLSINSTALLDIR)
else
   INSTALLDIR=$(CONFDIR)/bin/$(MLSCONFG)
endif

# ---------------------- Paths to special directories
# This is where this Makefile is located
# Also beginning path for all relative paths defined below; e.g. CONFDIR
# This assumes therefore that you have summoned this via
#    cd $(MLSHOME)
#    make
# otherwise the following definition will fail
# You might try something like
#    make -f mlspgs/Makefile MLSHOME=$HOME/mlspgs
# but this has not been tested.
# Anyway, why would you want to
# make things more difficult?
MLSHOME := $(shell pwd)

# where configure and makemakedep.sh are stored
MLSBIN=../util

# where sources for libmls.a are stored
libmls_sources=../lib

# where objects and libmls.a are stored
libmls_objs=$(libmls_sources)/$(MLSCONFG)

# where sources (but not OBJS) shared with l3, etc. are stored
srclib=../srclib

# where sources (but not OBJS) shared with l3, etc. are stored
blaslib=../blas

# where sources (but not OBJS) used for hdf5 file versions are stored
he5lib=../he5lib

# where he5 objects and libhe5.a are stored
he5lib_objs=$(he5lib)/$(MLSCONFG)

# where the platforms directory is
PLATFORMS = .

# This is the default name of the configure file
# It may be overridden for multiple platforms, special options, etc.
ifdef HOSTMLSCFILE
   MLSCFILE=$(HOSTMLSCFILE)
else
   MLSCFILE=.configure
endif

# This will set our MLS platform, compiler. etc.
# w/o complaining if it does nor exist yet
-include $(CONFDIR)/$(MLSCFILE)

# All the implicit pattern-matching involving suffixes
# will be done in the uniquely-named sub-directories, not here--so
.SUFFIXES:

# This will filter out globbed file expansions that fail
# e.g., '*.f90 *.f' -> '*.f90" if there aren't any .f files
REECHO=$(MLSBIN)/reecho.sh

# This is in case tar doesn't summon GNU tar, you may still use it via
# something like
#    make -f MakeFC tar_sips TAR=/usr/local/gnu/bin/tar
TAR=tar

SHELL = /bin/sh

#--------------- Target definitions
# $PROG          -- automatic; built when target is not specified

#                           (functional)
# configure      -- runs mlsconfigure to (re)set compiler, platform, paths
#                  also (re)creates uniquely-named sub-sub-directories
#                  with custom Makefiles where *.o, *.mod, etc. built
# depends        -- runs makemakedep.sh to calculate dependencies
# mostlyclean    -- deletes uniquely-named sub-sub-directories
# clean          -- deletes only *.o, *.mod from sub-subs
# tar            -- creates a directory archive suitable for distribution
#

#----------------------------------------------------------
#wildcard seems to no longer work right (maybe some change to our suns?)
libmls_f90 := $(shell ${REECHO} $(libmls_sources)/*.f $(libmls_sources)/*.f90 \
  $(libmls_sources)/*.f9h $(srclib)/*.f9h $(libmls_sources)/*.c \
  $(libmls_sources)/lit_names.txt)
blaslib_f := $(shell ${REECHO} $(blaslib)/*.c $(blaslib)/*.f $(blaslib)/*.f90 \
  $(blaslib)/*.f9h $(srclib)/*.f9h)
he5lib_f := $(shell ${REECHO} $(he5lib)/*.c $(he5lib)/*.f $(he5lib)/*.f90)

PROG_prereqs := $(CONFDIR)/$(MLSCFILE) \
   $(PLATFORMS)/$(MLSCONFG) $(libmls_objs)/libmls.a \
   $(INSTALLDIR)/libutctotai.a

ifdef HDFVERSIONS
   PROG_prereqs := $(PROG_prereqs) $(he5lib_objs)/libhe5.a
   ADD2PATHS := $(ADD2PATHS) ../$(he5lib)/$(MLSCONFG) \
   $(HDF5)/../fortran/src $(HDF5)/../include
#   ADD2PATHS := $(ADD2PATHS) ../$(he5lib)/$(MLSCONFG)
endif

ifndef DONTBUILDPREQS
ifdef BLAS
   LIB_prereqs=$(libmls_f90) $(CONFDIR)/$(MLSCFILE) \
   $(libmls_sources)/machines/$(MLSCONFG)/machine.f90 \
   $(REECHO) $(libmls_objs)
else
   LIB_prereqs=$(libmls_f90) $(CONFDIR)/$(MLSCFILE) \
   $(libmls_sources)/machines/$(MLSCONFG)/machine.f90 \
   $(blaslib_f) $(REECHO) $(libmls_objs)
endif
HE5_LIB_prereqs=$(libmls_objs)/libmls.a $(he5lib_f)
endif

#Any sources to explicitly exclude from dependency list?
OMIT_SRCS := $(OMIT) MLSL1.f90

PROG=mlsl1log mlsl1g mlsl1t
# PROG=mlsl1

OMIT_PREFXD := $(shell ${REECHO} -prefix=-d $(OMIT_SRCS))

#These are the source files for the drivers; they should be excluded from OBJS
DRIVERS=MLSL1G.f90 MLSL1log.f90 MLSL1T.f90
DRIVERS_PREFXD := $(shell ${REECHO} -prefix=-do $(DRIVERS))

#CASCADE determines whether to create/use .mod dependencies:
#if CASCADE, use older .o dependencies, which entail cascades
#of massive recompilation
#These are configurations for which I have tested CASCADE-avoidance
CAN_AVOID="NAG.Linux LF95.Linux NAG.Sun Sun.Sun"
#CAN_AVOID=""
ifndef CASCADE
  WE_CAN := $(shell echo $(CAN_AVOID) | grep $(MLSF95).$(MLSPLAT))
ifndef WE_CAN
  CASCADE := yes
endif
endif

ifdef WE_CAN
NOHUNTING := yes
endif
#These modules require custom builds, so omit build commands from dependencies
CUSTOM_BUILDS=OutputL1B.f90 OutputL1B_HDF4.f90 EngTbls.f90
CUSTOM_PREFXD := $(shell ${REECHO} -prefix=-db $(CUSTOM_BUILDS))
ifeq ($(MLSF95),Abs)
  UP_OR_LOW := -mod UPPER
else
  UP_OR_LOW := -mod lower
endif

#MARK_ALL_AS_UPTODATE determines whether, on rebuilding any module,
#all modules, rebuilt or not, should be marked as up to date
#ifndef MARK_ALL_AS_UPTODATE
#  MARK_ALL_AS_UPTODATE := $(shell echo $(CAN_AVOID) | grep $(MLSF95).$(MLSPLAT))
#endif

#NOCLEANING determines value of RECFG passed to mlsconfigure script
#RECFG determines whether changing .configure triggers deleting
#existing sub-directory, and thus recompiling every source file
#Set NOCLEANING to "yes" if you want to preserve existing compiled objects
#even when changing .configure file
ifdef NOCLEANING
   RECFG = -nd
endif

#----------------------- Build instructions

$(PROG):  $(PROG_prereqs)
	@cd $(PLATFORMS)/$(MLSCONFG) && $(MAKE) $(PROG) \
        UNPREFXDPTHS="$(ADD2PATHS) $(HDF5)/../fortran/src $(HDF5)/../include" ;\
        return_status=`expr $$?` ; \
        if [ $$return_status != "0" ]; then \
		     exit 2 ; \
        elif [ "$(WE_CAN)" != "" -a "$(MARK_ALL_AS_UPTODATE)" != "no" ]; then \
		     $(MAKE) mark_all_as_uptodate \
           UNPREFXDPTHS="$(ADD2PATHS)" MARK_ALL_AS_UPTODATE= ;\
         fi

# $(PLATFORMS)/$(MLSCONFG)/$(PROG): $(PROG)
$(PLATFORMS)/$(MLSCONFG)/mlsl1log \
$(PLATFORMS)/$(MLSCONFG)/mlslg \
$(PLATFORMS)/$(MLSCONFG)/mlsl1t : $(PROG)
	@cd $(PLATFORMS)/$(MLSCONFG) && $(MAKE) $(PROG) \
        UNPREFXDPTHS="$(ADD2PATHS) $(HDF5)/../fortran/src $(HDF5)/../include" ;\
        return_status=`expr $$?` ; \
        if [ $$return_status != "0" ]; then \
		     exit 2 ; \
        elif [ "$(WE_CAN)" != "" -a "$(MARK_ALL_AS_UPTODATE)" != "no" ]; then \
		     $(MAKE) mark_all_as_uptodate \
           UNPREFXDPTHS="$(ADD2PATHS)" MARK_ALL_AS_UPTODATE= ;\
         fi

$(libmls_objs)/libmls.a: $(LIB_prereqs)
	@if [ "$(NOHUNTING)" != "" ] ; then \
	   rm -f $(libmls_objs)/libmls.a ; \
	fi ; \
	$(MAKE) -f $(MakeFName) libmls.a -C $(libmls_sources) \
     MARK_ALL_AS_UPTODATE=no DONTBUILDPREQS=TRUE

$(he5lib_objs)/libhe5.a: $(HE5_LIB_prereqs)
	@if [ "$(NOHUNTING)" != "" ] ; then \
	   rm -f $(he5lib_objs)/libhe5.a ; \
	fi ; \
	$(MAKE) -f $(MakeFName) libhe5.a -C $(he5lib) \
     MARK_ALL_AS_UPTODATE=no DONTBUILDPREQS=TRUE

$(INSTALLDIR)/libutctotai.a:
	$(MAKE) -f $(MakeFName) utctotai -C $(CONFDIR)

.PHONY: clean configure depends depends--mod ghostbuster install \
  internalsaymake internalupdate mlsl1log mlsl1g mlsl1t mostlyclean \
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
ifdef CASCADE
	$(MLSBIN)/makemakedep.sh -nodep -p Makefile.nodep \
      $(OMIT_PREFXD) $(DRIVERS_PREFXD)
	$(MLSBIN)/makemakedep.sh $(OMIT_PREFXD) $(DRIVERS_PREFXD) \
      $(libmls_sources) $(srclib) 
else
	$(MAKE) depends--mod
endif

depends--mod:
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

internalsaymake:
	@echo "Copying or linking MakeFC to Makefile in ./" ; \
	cp MakeFC Makefile

internalupdate: internalsaymake
	@echo "Updating the Makefile in ./"; \
    m_option="-machine" ; \
		$(MLSBIN)/mlsconfigure -dc "$(CONFDIR)" -pc "$(PLATFORMS)" \
         -cf "$(MLSCFILE)" -nc -nd $$m_option; \
      m_option=""

# install $(PROG) into $(INSTALLDIR)
install: $(PROG)
	@echo "Installing $(PROG) into $(INSTALLDIR)"; \
	if [ ! -d $(INSTALLDIR) ] ; then \
	   $(MLSBIN)/mkpath.sh $(INSTALLDIR) ; \
	fi ; \
	if [ -f $(PLATFORMS)/$(MLSCONFG)/mlsl1log ] ; then \
	   mv $(PLATFORMS)/$(MLSCONFG)/mlsl1log $(INSTALLDIR) ; \
	fi ; \
	if [ -f $(PLATFORMS)/$(MLSCONFG)/mlsl1t ] ; then \
	   mv $(PLATFORMS)/$(MLSCONFG)/mlsl1t $(INSTALLDIR) ; \
	fi ; \
	if [ -f $(PLATFORMS)/$(MLSCONFG)/mlsl1g ] ; then \
	   mv $(PLATFORMS)/$(MLSCONFG)/mlsl1g $(INSTALLDIR) ; \
      rm -f $(INSTALLDIR)/mlsl1.sh ; \
	   sed "s=mlsxxyyzz_1=mlsl1log=; s=mlsxxyyzz_2=mlsl1g=; \
        s=mlsxxyyzz_3=mlsl1t=; s=mlsbbiinn=$(INSTALLDIR)=; \
        s=mlshhoommee=$(MLSHOME)=; s=mlseexxttrraa==" \
        $(MLSBIN)/mlsl1prog.sh > $(INSTALLDIR)/mlsl1.sh ; \
      chmod a+x $(INSTALLDIR)/mlsl1.sh ; \
	fi

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


