#	Makefile.h
# An include file for all mainline Makefiles
# We assemble here all the dependencies, including hierarchical ones,
# to avoid having to repeat them when, e.g., building in l2 what
# we already had to say about when and how to build ib lib
#
# Copyright 2010, by the California Institute of Technology. ALL
# RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
# commercial use must be negotiated with the Office of Technology Transfer
# at the California Institute of Technology.

# This software may be subject to U.S. export control laws. By accepting this
# software, the user agrees to comply with all applicable U.S. export laws and
# regulations. User has the responsibility to obtain export licenses, or other
# export authority as may be required before exporting such information to
# foreign countries or providing access to foreign persons.

# "$Id$"

# ---------------------- The name of this file (is the right-hand-side below)
MakeFName = MakeFC

# ---------------------- Paths to special directories
# where to store .configure
CONFDIR=..

# where configure and makemakedep.sh are stored
MLSBIN=../util

# where sources for cfm are stored
cfm_sources=../cfm

# where objects and libcfm.a are stored
cfm_objs=$(cfm_sources)/$(MLSCONFG)

# where sources for l1 are stored
l1_sources=../l1

# where objects and libmls.a are stored
l1_objs=$(l1_sources)/$(MLSCONFG)

# where sources for l2 are stored
l2_sources=../l2

# where objects and libmls.a are stored
l2_objs=$(l2_sources)/$(MLSCONFG)

# where sources for libmls.a are stored
libmls_sources=../lib

# where objects and libmls.a are stored
libmls_objs=$(libmls_sources)/$(MLSCONFG)

# where sources (but not OBJS) shared with l3, etc. are stored
srclib=../srclib

# where sources (but not OBJS) shared with l3, etc. are stored
blaslib=../blas

# where blas objects and libmlspack.a are stored
libblas_objs=$(blaslib)/$(MLSCONFG)

# where sources for the forward model are stored
libfwdmdl_sources=../fwdmdl

# where forward model objects and libfwdmdl.a are stored
libfwdmdl_objs=$(libfwdmdl_sources)/$(MLSCONFG)

# where sources for the cloud forward model are stored
libcloud_sources=../cloudfwdm

# where cloud forward model objects and libfwdmdl.a are stored
libcloud_objs=$(libcloud_sources)/$(MLSCONFG)

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

ifdef MLSINSTALLDIR
   INSTALLDIR:=$(MLSINSTALLDIR)
else
ifndef INSTALLDIR
   INSTALLDIR:=$(CONFDIR)/bin/$(MLSCONFG)
endif
endif

ifdef HDF5INC
   ADD2PATHS := $(ADD2PATHS) $(HDF5INC)
else
   ADD2PATHS := $(ADD2PATHS) $(HDF5)/../fortran/src $(HDF5)/../include
endif

# Pattern-matching rules for preprocessing
.SUFFIXES: .f90 .F90

.F90.f90:
	cpp -D$(MLSF95) $< $@

ifeq ($(MLSF95),Abs)
  UP_OR_LOW := -mod UPPER
else
  UP_OR_LOW := -mod lower
endif

# This will filter out globbed file expansions that fail
# e.g., '*.f90 *.f' -> '*.f90" if there aren't any .f files
REECHO=$(MLSBIN)/reecho.sh
UNIQUE_NAME=$(MLSBIN)/unique_name.sh

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

# Files to get by preprocessing
export PREPRO=$(patsubst %.F90,%.f90,$(wildcard *.F90))

cfm_f90 := $(shell ${REECHO} $(cfm_sources)/*.f $(cfm_sources)/*.[Ff]9? $(cfm_sources)/*.c)
l1_f90 := $(shell ${REECHO} $(l1_sources)/*.f $(l1_sources)/*.[Ff]9? $(l1_sources)/*.c)
l2_f90 := $(shell ${REECHO} $(l2_sources)/*.f $(l2_sources)/*.[Ff]9? $(l2_sources)/*.c)
libmls_f90 := $(shell ${REECHO} $(libmls_sources)/*.f $(libmls_sources)/*.[Ff]9? \
  $(srclib)/*.f9h $(libmls_sources)/*.c \
  $(libmls_sources)/lit_names.txt)
blaslib_f := $(shell ${REECHO} $(blaslib)/*.c $(blaslib)/*.f $(blaslib)/*.[Ff]9?)
libfwdmdl_f90 := $(shell ${REECHO} $(libfwdmdl_sources)/*.c $(libfwdmdl_sources)/*.f $(libfwdmdl_sources)/*.[Ff]9?)
libcloud_f90 := $(shell ${REECHO} $(libcloud_sources)/*.f90 $(libcloud_sources)/*.f $(libcloud_sources)/*.c)

ifdef DONTBUILDPREQS
  BLAS_PREQS := $(blaslib_f)
  LIB_prereqs=$(libmls_f90)
  FM_LIB_prereqs=$(libfwdmdl_f90)
  CLD_LIB_prereqs=$(libcloud_f90)
  L1_LIB_prereqs=$(l1_f90)
  L2_LIB_prereqs=$(l2_f90)
  CFM_LIB_prereqs=$(cfm_f90)
else
  BLAS_PREQS := $(blaslib_f) $(libblas_objs)/Makefile
  LIB_prereqs=$(libmls_f90) $(CONFDIR)/$(MLSCFILE) \
     $(machine_f) \
     $(libblas_objs)/libmlspack.a $(REECHO) $(libmls_objs)/Makefile
  FM_LIB_prereqs=$(libfwdmdl_f90) $(libmls_objs)/libmls.a \
     $(libfwdmdl_objs)/Makefile
  CLD_LIB_prereqs=$(libcloud_f90) $(libmls_objs)/libmls.a \
     $(libfwdmdl_objs)/libfwdmdl.a \
     $(libcloud_objs)/Makefile
  L1_LIB_prereqs=$(l1_f90) $(libmls_objs)/libmls.a $(l1_objs)/Makefile
  L2_LIB_prereqs=$(l2_f90) $(libmls_objs)/libmls.a \
     $(libfwdmdl_objs)/libfwdmdl.a $(l2_objs)/Makefile
  CFM_LIB_prereqs=$(cfm_f90) $(libmls_objs)/libmls.a \
     $(libfwdmdl_objs)/libfwdmdl.a \
     $(libcloud_objs)/libcloud.a $(l2_objs)/mlsl2 $(cfm_objs)/Makefile
endif

UPTODATEMARKS=MARK_ALL_AS_UPTODATE=no DONTBUILDPREQS=TRUE

# A pseudo-target which we hope to build successfully
#all: $(PROG)
$(PROG): $(MakeFName)

$(libblas_objs)/libmlspack.a: $(BLAS_PREQS)
	$(MAKE) -f $(MakeFName) libmlspack.a -C $(blaslib) $(UPTODATEMARKS)

$(libmls_objs)/libmls.a: $(LIB_prereqs)
	echo Makefile.h $(MAKE) -f $(MakeFName) libmls.a -C $(libmls_sources) $(UPTODATEMARKS)
	$(MAKE) -f $(MakeFName) libmls.a -C $(libmls_sources) $(UPTODATEMARKS)

$(libfwdmdl_objs)/libfwdmdl.a: $(FM_LIB_prereqs)
	echo Makefile.h $(MAKE) -f $(MakeFName) libfwdmdl.a -C $(libfwdmdl_sources) $(UPTODATEMARKS)
	$(MAKE) -f $(MakeFName) libfwdmdl.a -C $(libfwdmdl_sources) $(UPTODATEMARKS)

$(libcloud_objs)/libcloud.a: $(CLD_LIB_prereqs)
	echo Makefile.h $(MAKE) -f $(MakeFName) libcloud.a -C $(libcloud_sources) $(UPTODATEMARKS)
	$(MAKE) -f $(MakeFName) libcloud.a -C $(libcloud_sources) $(UPTODATEMARKS)

$(l1_objs)/mlsl1: $(l1_LIB_prereqs)
	$(MAKE) -f $(MakeFName) mlsl1 -C $(l1_sources) $(UPTODATEMARKS)

$(l2_objs)/mlsl2: $(L2_LIB_prereqs)
	$(MAKE) -f $(MakeFName) mlsl2 -C $(l2_sources) $(UPTODATEMARKS)

$(cfm_objs)/libcfm_all.a: $(CFM_LIB_prereqs)
	echo Makefile.h $(MAKE) -f $(MakeFName) libcfm_all.a -C $(cfm_sources) $(UPTODATEMARKS)
	$(MAKE) -f $(MakeFName) libcfm_all.a -C $(cfm_sources) $(UPTODATEMARKS)

$(INSTALLDIR)/libutctotai.a:
	$(MAKE) -f $(MakeFName) utctotai -C $(CONFDIR) $(UPTODATEMARKS)

# $Log$
# Revision 1.1  2010/10/12 22:20:26  pwagner
# First commit
#
