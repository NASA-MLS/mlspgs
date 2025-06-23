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


# This is the default name of the configure file
# It may be overridden for multiple platforms, special options, etc.
ifdef HOSTMLSCFILE
   MLSCFILE=$(HOSTMLSCFILE)
else
   MLSCFILE=.configure
endif

# where to store .configure; 
CONFDIR=../..

UTILDIR=$(CONFDIR)/util

# reecho.sh can filter out globbed file expansions that fail
# e.g., '*.f90 *.f' -> '*.f90" if there aren't any .f files
REECHO=$(UTILDIR)/reecho.sh

# Get any settings from our .configure
include $(CONFDIR)/$(MLSCFILE)

# Most source files match one of the following patterns
f9h_sources := $(shell ${REECHO} *.f9*)
tex_sources := $(shell ${REECHO} *.tex)

ifdef LATEXPATH
   BIBTEX=$(LATEXPATH)/bibtex
   LATEX=$(LATEXPATH)/latex
   DVIPDF=$(LATEXPATH)/dvipdfm
   PDFLATEX=$(LATEXPATH)/pdflatex
else
   BIBTEX=bibtex
   LATEX=latex
   DVIPDF=dvipdfm
   PDFLATEX=pdflatex
endif

# Set paths for finding source files

DVIDIR = ./dvi/
LOGDIR = ./log/
PDFDIR = ./pdf/
AUX=cfm.aux
BIB=cfmbib.bib
BBL=cfm.bbl

vpath %.f9h ./
vpath %.tex ./
vpath %.dvi $(DVIDIR)
vpath %.l2cf $(CONFDIR)/l2/l2cf/lib
vpath %.log $(LOGDIR)
vpath %.pdf $(PDFDIR)

# Define some suffix rules

.SUFFIXES:
.SUFFIXES: .tex .dvi .pdf
.tex.dvi:
	$(LATEX) $<
	$(LATEX) $<
	mv $@ $(DVIDIR)
	mv $*.log $(LOGDIR)

.dvi.pdf:
	# $(DVIPDF) $(DVIDIR)/$*.dvi
	$(PDFLATEX) $*.tex
	mv $@ $(PDFDIR)

# Bugs and limitations:
# (1) We are assuming latex and dvipdf are both available
# An alternative approach would skip building dvi files; e.g.
#   pdflatex cfm
#   bibtex cfm
#   pdflatex cfm
#   pdflatex cfm
# (2) Dependencies can be traced from .tex file using makemake.sh
# It cannot track .f9h nor .f90 files yet, which is why we use the 
# f9h_sources mechanism
# (3) Nor can it handle .bib files, and thus we resort to the BIB mechanism

cfm:

-include Makefile.dep

all: cfm

cfm: ${NEAROBJS} ${OBJS} ${BBL}
	date > cfm; echo $(OBJS) >> cfm; \
   touch ${BBL}; touch $(DVIDIR)/*; touch $(PDFDIR)/*

${AUX}: cfm.tex
	$(LATEX) cfm

${BBL}: ${BIB} ${AUX}
	$(BIBTEX) cfm

${NEAROBJS}: ${BBL} ${f9h_sources} ${tex_sources}
# ========================================================================

clean:
	rm -fr *.dvi *.pdf Makefile.dep $(DVIDIR) $(PDFDIR) ${AUX} ${BBL} cfm

depends:
	$(UTILDIR)/makemakedep.sh -s tex -o dvi -p Makefile.dep.1 -S cfm.tex \
     -i '\\input\{([^)]+)}' \
     $(OMIT_PREFXD) $(OMIT_PREFXD_REV)
	sed 's/OBJS =/NEAROBJS =/' Makefile.dep.1 > Makefile.dep.2
	$(UTILDIR)/makemakedep.sh -s tex -o pdf -p Makefile.dep.1 -S cfm.tex \
     -i '\\input\{([^)]+)}' \
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
tar:
	rm -f cfm.tgz
	tar hcvzf cfm.tgz * --exclude \*.dvi --exclude \*.pdf --exclude cfm

update: 
	$(MAKE) -f $(MakeFName) depends

.PHONY: all clean depends tar update cfm \

# ========================================================================
# $Log$
# Revision 1.6  2021/04/21 23:43:38  pwagner
# Added some robustness
#
# Revision 1.5  2020/04/08 21:48:17  pwagner
# Gets settings from .configure; corrected dependencies
#
# Revision 1.4  2019/06/27 21:08:56  pwagner
# make clean more thorough
#
# Revision 1.3  2019/06/24 17:49:53  pwagner
# Works more reliably like this
#
# Revision 1.2  2015/08/11 00:25:00  pwagner
# Fixed a few errors
#
# Revision 1.1  2015/08/11 00:19:55  pwagner
# After switching away from MakeFC
#
# Revision 1.4  2011/05/25 20:34:42  pwagner
# CONFDIR reverted to point to grandparent directory
#
# Revision 1.3  2010/06/30 22:07:46  pwagner
# make update removes sometimes troublesome .toc file
#
# Revision 1.2  2010/01/06 17:35:48  pwagner
# Does a better job of building only when it must
#
# Revision 1.1  2010/01/05 19:21:26  pwagner
# First commit
#
