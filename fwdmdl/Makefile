#	Makefile
#
# Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
# U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.
#
#special case:
#the 1st time through, before .configure exists yet,
#execute $(MLSBIN)/configure to create .configure

#i.e., if you call are building this mls
#by calling make with a different Makefile name,
#say: make -f MakeFC
#then rename the Makefiles in /l2 and /lib also MakeFC

# "$Id$"

MakeFName = MakeFC

# where to store .configure
CONFDIR=..

# where configure and makemakedep.sh are stored
MLSBIN=../util

# where sources for libmls.a are stored
libmls_sources=../lib

# where objects and libmls.a are stored
libmls_objs=$(libmls_sources)/$(MLSCONFG)

# where sources (but not OBJS) shared with l3, etc. are stored
srclib=../srclib

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
libmls_f90 := $(shell ${REECHO} $(libmls_sources)/*.f $(libmls_sources)/*.f90 $(libmls_sources)/*.c $(libmls_sources)/lit_names.txt)

ifdef BLAS
   LIB_prereqs=$(libmls_f90) $(CONFDIR)/$(MLSCFILE) \
   $(libmls_sources)/machines/$(MLSCONFG)/machine.f90 \
   $(REECHO)
else
   LIB_prereqs=$(libmls_f90) $(CONFDIR)/$(MLSCFILE) \
   $(libmls_sources)/machines/$(MLSCONFG)/machine.f90 \
   $(blaslib_f) $(REECHO)
endif

#Any sources to explicitly exclude from dependency list?
OMIT_SRCS := $(OMIT) 

#PROG=l2_test
PROG=libfwdmdl.a

OMIT_PREFXD := $(shell ${REECHO} -prefix=-d $(OMIT_SRCS))

#NOCLEANING determines value of RECFG passed to mlsconfigure script
#RECFG determines whether changing .configure triggers deleting
#existing sub-directory, and thus recompiling every source file
#Set NOCLEANING to "yes" if you want to preserve existing compiled objects
#even when changing .configure file
ifdef NOCLEANING
   RECFG = -nd
endif

#----------------------- Build instructions

$(PROG): $(CONFDIR)/$(MLSCFILE) \
   $(PLATFORMS)/$(MLSCONFG) $(libmls_objs)/libmls.a
	cd $(PLATFORMS)/$(MLSCONFG) && $(MAKE) $(PROG)

$(PLATFORMS)/$(MLSCONFG)/$(PROG): $(PROG)
	cd $(PLATFORMS)/$(MLSCONFG) && $(MAKE) $(PROG)

l2_test: $(PLATFORMS)/$(MLSCONFG)/$(PROG)
	cd $(PLATFORMS)/$(MLSCONFG) && $(MAKE) l2_test

$(libmls_objs)/libmls.a: $(LIB_prereqs) $(libmls_objs)
	@if [ "$(NOHUNTING)" != "" ] ; then \
	   rm -f $(libmls_objs)/libmls.a ; \
	fi ; \
	$(MAKE) -f $(MakeFName) libmls.a -C $(libmls_sources)

.PHONY: configure depends clean tar mostlyclean targ tarz libfwdmdl.a \
 l2_test ghostbuster

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
	$(MLSBIN)/makemakedep.sh -nodep $(OMIT_PREFXD) ; \
     mv Makefile.dep Makefile.nodep
	$(MLSBIN)/makemakedep.sh $(OMIT_PREFXD) $(libmls_sources) $(srclib)

ghostbuster:
	$(MLSBIN)/ghostbuster.sh -lib $(OMIT_PREFXD) \
           $(PLATFORMS)/$(MLSCONFG) ./ $(srclib)

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
		rm -f $(PLATFORMS)/$(MLSCONFG)/l2_test ; \
	fi
	rm -f Make.dep.*

tar:
	rm -f ../fwdmdl.tar
	@echo "tar-ing fwdmdl"; \
	$(TAR) hcvf ../fwdmdl.tar * --exclude CVS\
	--exclude Make.dep\* --exclude NAG.\* --exclude LF95.\*

tarz: tar
	@echo "... then gzipping it"; \
	cd ..; rm -f fwdmdl.tgz; gzip fwdmdl.tar; mv fwdmdl.tar.gz fwdmdl.tgz

targ:
	rm -f ../fwdmdl.tgz
	$(TAR) hcvzf ../fwdmdl.tgz * --exclude \*.o --exclude \*.mod --exclude l2_test

# $Log$
# Revision 2.6  2002/05/02 00:16:32  pwagner
# depends now creates Makefile.nodep, too
#
# Revision 2.5  2002/01/15 00:40:18  pwagner
# Now uses user-settable OMIT and internal OMIT_SRCS and OMIT_PREFIXD variables
#
# Revision 2.4  2001/11/08 00:40:51  pwagner
# NOCLEANING avoids massive recompiles on changes to .configure
#
# Revision 2.3  2001/10/04 22:25:36  pwagner
# Limit suffix-based rules to (MLSCONFG); add lib/lit_names.txt to dependency list
#
# Revision 2.2  2001/09/26 17:24:58  pwagner
# Removed errors introduced with last changes; SIMPLIFYMAKE renamed NOHUNTING
#
# Revision 2.1  2001/09/25 20:27:15  pwagner
# New SIMPLIFYMAKE setting; other changes
#
# Revision 2.0  2001/09/17 20:26:26  livesey
# New forward model
#
# Revision 1.23  2001/09/07 22:23:52  pwagner
# Triggers build on changes in lib
#
# Revision 1.22  2001/08/22 22:46:13  pwagner
# Made reecho.sh a prerequisite where appropriate
#
# Revision 1.21  2001/08/15 21:01:04  pwagner
# Now excludes OMIT_SRCS from dependencies
#
# Revision 1.20  2001/08/13 18:09:08  pwagner
# replaced shell echo -> shell {REECHO}
#
# Revision 1.19  2001/07/25 23:19:58  pwagner
# Removed unnecessary lib...a from target prereqs
#
# Revision 1.18  2001/05/31 20:24:49  pwagner
# ghostbuster called with -lib option
#
# Revision 1.17  2001/05/24 17:56:39  pwagner
# Expanded LIB_prereqs
#
# Revision 1.16  2001/05/17 20:02:59  pwagner
# Automatically rebuilds lib on change to machine.f90
#
# Revision 1.15  2001/05/14 23:58:42  pwagner
# default .configure name now set by environment
#
# Revision 1.14  2001/05/14 17:37:14  pwagner
# .configure now a variable file name
#
# Revision 1.13  2001/05/01 18:32:13  pwagner
# Added .configure to prereqs for lib
#
# Revision 1.12  2001/05/01 17:07:43  pwagner
# Now runs ghostbuster when updating
#
# Revision 1.11  2001/04/12 22:03:08  pwagner
# Changes to lib or fwdmdl trigger rebuilds even on suns
#
# Revision 1.10  2001/04/05 17:15:54  pwagner
# Added PROG to .PHONY list
#
# Revision 1.9  2001/03/13 19:07:59  pwagner
# clean now removes executables, too
#
# Revision 1.8  2001/03/12 22:46:57  pwagner
# Revised clean, tar; automatic symlink to Makefile
#
# Revision 1.7  2001/03/09 19:36:06  pwagner
# Added TAR=... option
#
# Revision 1.6  2001/03/06 17:36:29  pwagner
# libfwdmdl.a now the automatic target
#
# Revision 1.5  2001/02/01 00:47:02  pwagner
# Making libfwdmdl.a and l2_test
#
# Revision 1.4  2001/01/31 23:58:24  pwagner
# cConsistent with zvi's latest changes
#
# Revision 1.3  2001/01/10 23:45:33  pwagner
# General housekeeping
#
# Revision 1.2  2000/11/22 23:40:57  pwagner
# Fixed program name
#
# Revision 1.1  2000/11/22 23:34:06  pwagner
# First commit
#
# Revision 1.1  2000/11/22 22:51:32  pwagner
# First commit
#
#


