#	Makefile
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
#
#special case:
#the 1st time through, before .configure exists yet,
#execute $(MLSBIN)/configure to create .configure

#i.e., if you call are building this mls
#by calling make with a different Makefile name,
#say: make -f $(MakeFName)
#then rename the Makefiles in /l2 and /lib also $(MakeFName)

# "$Id$"

MakeFName = Makefile

# where to store .configure
CONFDIR=../..

# where configure and makemakedep.sh are stored
MLSBIN=../../util

# where sources for libmls.a are stored
libmls_sources=../../lib

# where objects and libmls.a are stored
libmls_objs=$(libmls_sources)/$(MLSCONFG)

# where sources (but not OBJS) shared with l3, etc. are stored
srclib=../../srclib

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

# This will filter out globbed file expansions that fail
# e.g., '*.f90 *.f' -> '*.f90" if there aren't any .f files
REECHO=$(MLSBIN)/reecho.sh
UNIQUE_NAME=$(MLSBIN)/unique_name.sh

# This is in case tar doesn't summon GNU tar, you may still use it via
# something like
#    make tar_sips TAR=/usr/local/gnu/bin/tar
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
libmls_f90 := $(shell ${REECHO} $(libmls_sources)/*.f $(libmls_sources)/*.f90 $(libmls_sources)/*.c)
blaslib_f := $(shell ${REECHO} $(blaslib)/*.c $(blaslib)/*.f $(blaslib)/*.f90)

ifdef BLAS
   LIB_prereqs=$(libmls_f90) $(CONFDIR)/$(MLSCFILE) \
   $(libmls_sources)/machines/$(MLSCONFG)/machine.f90
else
   LIB_prereqs=$(libmls_f90) $(CONFDIR)/$(MLSCFILE) \
   $(libmls_sources)/machines/$(MLSCONFG)/machine.f90 \
   $(blaslib_f)
endif
NEEDS_ITM="no"
PROG=test

ifdef HDF5INC
   ADD2PATHS := $(ADD2PATHS) $(HDF5INC)
else
   ADD2PATHS := $(ADD2PATHS) $(HDF5)/../fortran/src $(HDF5)/../include
endif

#----------------------- Build instructions

$(PROG): $(CONFDIR)/$(MLSCFILE) $(PLATFORMS)/$(MLSCONFG) $(libmls_objs)/libmls.a
		cd $(PLATFORMS)/$(MLSCONFG) && $(MAKE) $(PROG) \
        UNPREFXDPTHS="$(ADD2PATHS)" NEEDS_ITM=$(NEEDS_ITM)

$(libmls_objs)/libmls.a: $(LIB_prereqs)
	$(MAKE) -f $(MakeFName) libmls.a -C $(libmls_sources)

.PHONY: configure depends clean tar mostlyclean tarz targ test ghostbuster update

tar:

configure:
	$(MLSBIN)/mlsconfigure -dc "$(CONFDIR)" -pc "$(PLATFORMS)" \
          -c "$(MLSF95)" -p "$(MLSPLAT)" -cf "$(MLSCFILE)"

$(CONFDIR)/$(MLSCFILE):
	$(MLSBIN)/mlsconfigure -dc "$(CONFDIR)" -pc "$(PLATFORMS)" \
          -c "$(MLSF95)" -p "$(MLSPLAT)" -cf "$(MLSCFILE)"

$(PLATFORMS)/$(MLSCONFG): $(CONFDIR)/$(MLSCFILE)
	$(MLSBIN)/mlsconfigure -dc "$(CONFDIR)" -pc "$(PLATFORMS)" \
          -c "$(MLSF95)" -p "$(MLSPLAT)" -cf "$(MLSCFILE)" -nc

depends:
	temp_name_1=`${UNIQUE_NAME} 1` ;\
	$(MLSBIN)/makemakedep.sh -p $$temp_name_1 \
        $(libmls_sources) $(srclib) ; \
     mv $$temp_name_1 Makefile.dep

ghostbuster:
	$(MLSBIN)/ghostbuster.sh $(PLATFORMS)/$(MLSCONFG) ./ $(srclib)

update:
	@cd ../.. ; \
   $(MAKE) update SUBDIRS= tests_sub_dirs=tests/lib

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
	@echo "tar-ing l3"; \
	$(TAR) hcvf ../$(PROG).tar * --exclude CVS\
	--exclude Make.dep\* --exclude NAG.\* --exclude LF95.\*

tarz: tar
	@echo "... then gzipping it"; \
	cd ..; rm -f $(PROG).tgz; gzip $(PROG).tar; mv $(PROG).tar.gz $(PROG).tgz

targ:
	rm -f ../$(PROG).tgz
	$(TAR) hcvzf ../$(PROG).tgz * --exclude \*.o --exclude \*.mod --exclude test

# $Log$
# Revision 1.13  2010/08/27 20:56:39  pwagner
# Uses optional HDF*INC Makefile variables
#
# Revision 1.12  2009/10/05 16:44:49  pwagner
# Added hdfeos5 to extra paths so we can build aug_hdfeos5
#
# Revision 1.11  2008/08/20 18:40:06  pwagner
# Added NEEDS_ITM setting to control whether to link with init_tables_module
#
# Revision 1.10  2005/06/23 22:10:41  pwagner
# Reworded Copyright statement
#
# Revision 1.9  2003/04/30 23:04:55  pwagner
# Now lets you make update in each subdirectory
#
# Revision 1.8  2002/09/25 23:25:43  pwagner
# Uses unique names for temp files
#
# Revision 1.7  2002/07/24 20:32:48  pwagner
# Changes to support hdf5/hdfeos5
#
# Revision 1.6  2001/08/13 18:09:08  pwagner
# replaced shell echo -> shell {REECHO}
#
# Revision 1.5  2001/05/31 20:21:59  pwagner
# Added ghostbuster to .PHONY list
#
# Revision 1.4  2001/05/24 23:23:24  pwagner
# Fixed prob with LIB_prereqs being undefined
#
# Revision 1.3  2001/05/24 18:04:57  pwagner
# Expanded LIB_prereqs
#
# Revision 1.2  2001/05/18 23:30:22  pwagner
# Brought into harmony with latest non-tests versions
#
# Revision 1.1  2001/02/01 23:11:12  pwagner
# First commit
#
#


