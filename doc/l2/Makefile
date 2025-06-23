# "$Id$"
#	Makefile
#
# Copyright 2009, by the California Institute of Technology. ALL
# RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
# commercial use must be negotiated with the Office of Technology Transfer
# at the California Institute of Technology.

# This software may be subject to U.S. export control laws. By accepting this
# software, the user agrees to comply with all applicable U.S. export laws and
# regulations. User has the responsibility to obtain export licenses, or other
# export authority as may be required before exporting such information to
# foreign countries or providing access to foreign persons.

# ---------------------- The name of this file (is the right-hand-side below)
MakeFName = Makefile
# If the right-hand-side is not Makefile
# then to make mls with this file you must instead  of bare "make" use
#     make -f $(MakeFName)

# The name of this file must be used in each of $(SUBDIRS)
# i.e., if you are calling make with a different Makefile name,
#        say: make -f MakeFC
# then (re)name the Makefiles in /l2 /lib etc. also MakeFC

# where to store .configure; 
#CONFDIR=..
CONFDIR=${HOME}/mlspgs

UTILDIR=$(CONFDIR)/util

# reecho.sh can filter out globbed file expansions that fail
# e.g., '*.f90 *.f' -> '*.f90" if there aren't any .f files
# However we will be using it just to prefix any unprefixed preset paths
REECHO=$(UTILDIR)/reecho.sh

SUBDIRS := $(shell ${REECHO} -d cfm l2)

# Errors in the .tex sources require us to omit the following
OMIT_SRCS := wvs-010.tex wvs-074.tex wvs-051-neural_body.tex
OMIT_PREFXD := $(shell ${REECHO} -prefix=-d $(OMIT_SRCS))
OMIT_SRCS_REV := 'wvs-004\*.tex' 'wvs-051\*.tex' 'wvs-063\*.tex'
OMIT_PREFXD_REV := $(shell ${REECHO} -prefix=-d -escape -nlast $(OMIT_SRCS_REV))

# Set paths for finding source files

DVIDIR = ./dvi/
LOGDIR = ./log/
PDFDIR = ./pdf/

vpath %.tex ./
vpath %.dvi $(DVIDIR)
vpath %.l2cf $(CONFDIR)/l2/l2cf/lib
vpath %.log $(LOGDIR)
vpath %.pdf $(PDFDIR)

# Define some suffix rules

.SUFFIXES:
.SUFFIXES: .tex .dvi .pdf
.tex.dvi:
	echo "\end" | latex $<
	echo "\end" | latex $<
	mv $@ $(DVIDIR)
	mv $*.log $(LOGDIR)

.dvi.pdf:
	dvipdf $<
	mv $@ $(PDFDIR)

doc:

-include Makefile.dep

all: doc

doc: ${NEAROBJS} ${OBJS}
	date > doc; echo $(OBJS) >> doc

# ========================================================================

clean:
	rm -fr *.dvi *.pdf $(DVIDIR) $(PDFDIR)

depends:
	$(UTILDIR)/makemakedep.sh -s tex -o dvi -p Makefile.dep.1 \
     $(OMIT_PREFXD) $(OMIT_PREFXD_REV)
	sed 's/OBJS =/NEAROBJS =/' Makefile.dep.1 > Makefile.dep.2
	$(UTILDIR)/makemakedep.sh -s tex -o pdf -p Makefile.dep.1 \
     $(OMIT_PREFXD) $(OMIT_PREFXD_REV)
	cat Makefile.dep.1 Makefile.dep.2 > Makefile.dep
	if [ ! -d $(DVIDIR) ] ; then \
	  		mkdir $(DVIDIR) ;\
			fi ; \
	if [ ! -d $(PDFDIR) ] ; then \
	  		mkdir $(PDFDIR) ;\
			fi ; \
	if [ ! -d $(LOGDIR) ] ; then \
	  		mkdir $(LOGDIR) ;\
			fi
	 for dir in $(SUBDIRS); do \
	  $(MAKE) -C $$dir depends; \
	done

cfm:
	$(MAKE) -f $(MakeFName) -C $@

l2:
	$(MAKE) -f $(MakeFName) -C $@

tar:
	rm -f doc.tgz
	tar hcvzf docs.tgz * --exclude \*.dvi --exclude \*.pdf --exclude doc

update:
	$(MAKE) depends

.PHONY: all clean depends tar update l2 \

# ========================================================================
# $Log$
# Revision 1.2  2009/09/22 17:46:29  pwagner
# Can now ignore sources marked by revision number
#
# Revision 1.1  2009/09/18 19:47:31  pwagner
# First commit
#
