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

short_name=idlcfm

PROG=server client

# We store all the common Makefile things in one shared file
SRCLIB=../srclib
-include $(SRCLIB)/Makefile.h

# This will set our MLS platform, compiler. etc.
# w/o complaining if it does nor exist yet
-include $(CONFDIR)/$(MLSCFILE)

#Where else to look for USEd and INCLUDED files to add to dependency list
depends_where:=$(libmls_sources) $(srclib) \
  $(libfwdmdl_sources) $(libcloud_sources) $(cfm_sources)                      
PROG_prereqs:=$(CONFDIR)/$(MLSCFILE) $(PLATFORMS)/$(MLSCONFG) \
  $(l2_objs)/mlsl2 \
  $(libmls_objs)/libmls.a \
  $(libfwdmdl_objs)/libfwdmdl.a \
  $(libcloud_objs)/libcloud.a \
  $(cfm_objs)/libcfm_all.a

OMIT_SRCS := $(OMIT) 

OMIT_PREFXD := $(shell ${REECHO} -prefix=-d $(OMIT_SRCS))

ifdef HDF5INC
   ADD2PATHS := $(ADD2PATHS) $(HDF5INC)
else
   ADD2PATHS := $(ADD2PATHS) $(HDF5)/../fortran/src $(HDF5)/../include
endif

#These are the source files for the drivers; they should be excluded from OBJS
DRIVERS=server.f90 client.f90
DRIVERS_PREFXD := $(shell ${REECHO} -prefix=-do $(DRIVERS))

#----------------------- Build instructions

all: $(PROG)

$(PROG): $(PROG_prereqs)
	cd $(PLATFORMS)/$(MLSCONFG) && $(MAKE) $(PROG) \
        UNPREFXDPTHS="$(ADD2PATHS)" \

$(PLATFORMS)/$(MLSCONFG)/server \
$(PLATFORMS)/$(MLSCONFG)/client: $(PROG_prereqs)
	cd $(PLATFORMS)/$(MLSCONFG) && $(MAKE) $(PROG) \
        UNPREFXDPTHS="$(ADD2PATHS)" \

.PHONY: configure depends clean tar mostlyclean tarz targ ghostbuster update

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
	$(MLSBIN)/makemakedep.sh  -p $$temp_name_1 \
     $(OMIT_PREFXD) $(DRIVERS_PREFXD) $(depends_where) ; \
     mv $$temp_name_1 Makefile.dep

ghostbuster:
	$(MLSBIN)/ghostbuster.sh -lib $(PLATFORMS)/$(MLSCONFG) ./ $(srclib)

update:
	@cd .. ; \
   $(MAKE) update SUBDIRS=idlcfm

install: $(PROG)
	@cp $(PLATFORMS)/$(MLSCONFG)/server \
     $(PLATFORMS)/$(MLSCONFG)/client $(INSTALLDIR); \
   echo "$(PROG) installed in $(INSTALLDIR)"

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
	@echo "tar-ing l2"; \
	$(TAR) hcvf ../$(PROG).tar * --exclude CVS\
	--exclude Make.dep\* --exclude NAG.\* --exclude LF95.\*

tarz: tar
	@echo "... then gzipping it"; \
	cd ..; rm -f $(PROG).tgz; gzip $(PROG).tar; mv $(PROG).tar.gz $(PROG).tgz

targ:
	rm -f ../$(PROG).tgz
	$(TAR) hcvzf ../$(PROG).tgz * --exclude \*.o --exclude \*.mod --exclude server

# $Log$
# Revision 1.3  2017/01/05 18:46:22  pwagner
# Changed method of disguising PseudoToolkit
#
# Revision 1.2  2016/03/03 18:56:29  pwagner
# Usual method of hiding files in lib fails from cfm; need long-term fix
#
# Revision 1.1  2015/08/10 23:21:50  pwagner
# After switching away from MakeFC
#
# Revision 1.1  2011/05/16 22:46:52  pwagner
# First commit
#
