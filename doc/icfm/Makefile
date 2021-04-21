# "$Id$"
#	Makefile
#
# Copyright 2011, by the California Institute of Technology. ALL
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

# where to store .configure; 
CONFDIR=../..

UTILDIR=$(CONFDIR)/util

# reecho.sh can filter out globbed file expansions that fail
# e.g., '*.f90 *.f' -> '*.f90" if there aren't any .f files
REECHO=$(UTILDIR)/reecho.sh

# Most source files match one of the following patterns
pro_sources := $(shell ${REECHO} *.pro)
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
AUX=icfm.aux
BIB=icfmbib.bib
BBL=icfm.bbl

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
#   pdflatex icfm
#   bibtex icfm
#   pdflatex icfm
#   pdflatex icfm
# (2) Dependencies can be traced from .tex file using makemake.sh
# It cannot track .f9h nor .f90 files yet, which is why we use the 
# f9h_sources mechanism
# (3) Nor can it handle .bib files, and thus we resort to the BIB mechanism

icfm:

-include Makefile.dep

all: icfm

icfm: ${NEAROBJS} ${OBJS} ${BBL}
	date > icfm; echo $(OBJS) >> icfm; \
   touch ${BBL}; touch $(DVIDIR)/*; touch $(PDFDIR)/*

${AUX}: icfm.tex
	latex icfm

${BBL}: ${BIB} ${AUX}
	bibtex icfm

${NEAROBJS}: ${BBL} ${pro_sources}
# ========================================================================

clean:
	rm -fr *.dvi *.pdf Makefile.dep $(DVIDIR) $(PDFDIR) ${AUX} ${BBL} icfm

depends:
	$(UTILDIR)/makemakedep.sh -s tex -o dvi -p Makefile.dep.1 -S icfm.tex \
     -i '\\input\{([^)]+)}' \
     $(OMIT_PREFXD) $(OMIT_PREFXD_REV)
	sed 's/OBJS =/NEAROBJS =/' Makefile.dep.1 > Makefile.dep.2
	$(UTILDIR)/makemakedep.sh -s tex -o pdf -p Makefile.dep.1 -S icfm.tex \
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
	rm -f icfm.tgz
	tar hcvzf icfm.tgz * --exclude \*.dvi --exclude \*.pdf --exclude icfm

update: 
	$(MAKE) -f $(MakeFName) depends

.PHONY: all clean depends tar update icfm \

# ========================================================================
# $Log$
# Revision 1.3  2019/06/27 21:09:08  pwagner
# make clean more thorough
#
# Revision 1.2  2019/06/24 17:45:48  pwagner
# Works more reliably like this
#
# Revision 1.1  2015/08/11 00:26:36  pwagner
# After switching away from MakeFC
#
# Revision 1.1  2011/05/31 22:10:47  pwagner
# First commit
#
