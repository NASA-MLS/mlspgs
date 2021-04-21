# "$Id$"
#	Makefile
#
#  Lost? Type 'make help -f MakeFC' from the directory where you see this file
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

# Help screens
# --------------- main help
# To see help on specific topics:
# Topic                               Command
# tex files in this directory         'make help--doc'
# tex files in cfm                    'make help--cfm'
# tex files in icfm                   'make help--icfm'
# tex files in l2cf                   'make help--l2cf'
# --------------- End main help
# --------------- doc help
# To build all the tex files in this directory type
# 'make' or 'make doc'
# to build just one, say wvs-022.pdf, type
# make wvs-022.pdf
# Look for the dvi and pdf products in the subdirectories named 'dvi' and 'pdf'
# --------------- End doc help
# --------------- ocfm help
# To build the documentation for the callable forward model, type
# 'make cfm'
# or else cd to the cfm subdirectory and just type 'make'
# From there look in the subdirectories named 'dvi' and 'pdf'
# --------------- End ocfm help
# --------------- icfm help
# To build the documentation for the idl-callable forward model, type
# 'make icfm'
# or else cd to the icfm subdirectory and just type 'make'
# From there look in the subdirectories named 'dvi' and 'pdf'
# --------------- End icfm help
# --------------- l2cf help
# To build the documentation for the l2cf, type
# 'make l2cf'
# or else cd to the l2cf subdirectory and just type 'make'
# From there look in the subdirectories named 'dvi' and 'pdf'
# --------------- End l2cf help
# ---------------------- The name of this file (is the right-hand-side below)
MakeFName = Makefile
# If the right-hand-side is not Makefile
# then to make mls with this file you must instead  of bare "make" use
#     make -f $(MakeFName)

# The name of this file must be used in each of $(SUBDIRS)
# i.e., if you are calling make with a different Makefile name,
#        say: make -f MakeFC
# then (re)name the Makefiles in /l2 /lib etc. also MakeFC

# This is the default name of the configure file
# It may be overridden for multiple platforms, special options, etc.
ifdef HOSTMLSCFILE
   MLSCFILE=$(HOSTMLSCFILE)
else
   MLSCFILE=.configure
endif

# where to store .configure; 
CONFDIR=..

UTILDIR=$(CONFDIR)/util

# Where non-standard build commands are kept
CUSTOM = $(CONFDIR)/customblds

# reecho.sh can filter out globbed file expansions that fail
# e.g., '*.f90 *.f' -> '*.f90" if there aren't any .f files
# However we will be using it just to prefix any unprefixed preset paths
REECHO=$(UTILDIR)/reecho.sh

# Get any settings from our .configure
include $(CONFDIR)/$(MLSCFILE)

short_name=doc

# l2 won't build yet (it lacks certain files, e.g. l2cf_Fill_BinMax.tex)
#SUBDIRS := $(shell ${REECHO} -d cfm)

# Errors in the .tex sources require us to omit the files named in
# tex.omitted
OMIT_SRCS := $(shell cat tex.omitted)
OMIT_PREFXD := $(shell ${REECHO} -prefix=-d $(OMIT_SRCS))
OMIT_SRCS_REV := 'wvs-004\*.tex' 'wvs-051\*.tex' 'wvs-063\*.tex'
OMIT_SRCS_SUB := 'iy-001-\*.tex' \
   'iy-002-\*.tex' 'iy-003-\*.tex' 'iy-004-\*.tex' 'iy-008-\*.tex' \
   'iy-009-\*.tex' 'iy-010-\*.tex'
OMIT_PREFXD_REV := $(shell ${REECHO} -prefix=-d -escape -nlast $(OMIT_SRCS_REV))
OMIT_PREFXD_SUB := $(shell ${REECHO} -prefix=-d -escape $(OMIT_SRCS_SUB))

# also where to install mls programs by default
ifdef MLSINSTALLDIR
   INSTALLDIR=$(MLSINSTALLDIR)
else
   INSTALLDIR=$(CONFDIR)/bin/$(MLSCONFG)
endif

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

TESTSDIR=$(CONFDIR)/tests

# Set paths for finding source files

# DOCDIR = $(CONFDIR)/doc
DOCDIR = $(shell pwd)
TEXDIR = $(DOCDIR)/tex
DVIDIR = $(DOCDIR)/dvi
LOGDIR = $(DOCDIR)/log
PDFDIR = $(DOCDIR)/pdf
MATH77DIR = $(DOCDIR)/Math77

vpath %.f $(MATH77DIR)
vpath %.tex $(DOCDIR)
vpath %.dvi $(DVIDIR)
vpath %.l2cf $(CONFDIR)/l2/l2cf/lib
vpath %.log $(LOGDIR)
vpath %.pdf $(PDFDIR)
vpath %.eps $(DOCDIR)
vpath %.jpg $(DOCDIR)

# Define some suffix rules

.SUFFIXES:
.SUFFIXES: .o .tex .dvi .pdf .eps .f .jpg .obj
.tex.dvi:
	echo "LATEXPATH $(LATEXPATH)"
	echo "LATEX $(LATEX)"
	echo "\def\dvidir{../dvi}" | cat - $< > $*.temp
	echo "\def\logdir{$(LOGDIR)}" | cat - $*.temp > $*.temp2
	echo "\def\pdfdir{../pdf}" | cat - $*.temp2 > $(TEXDIR)/$<
	echo "\end" | $(LATEX) $(TEXDIR)/$<; \
	echo "\end" | $(LATEX) $(TEXDIR)/$<
	mv $@ $(DVIDIR)

.dvi.pdf:
	# dvipdf $(DVIDIR)/$<
	if [ "$*" != "index" ]; then \
	  echo "\end" | $(LATEX) $(TEXDIR)/$*.tex; \
	  echo "\end" | $(LATEX) $(TEXDIR)/$*.tex; \
	  $(DVIPDF) $*.dvi; \
	  mv $*.dvi $(DVIDIR); \
	else \
	  $(PDFLATEX) $(TEXDIR)/$*.tex; \
	fi
	mv $@ $(PDFDIR)
	mv `${REECHO} $*.log $*.aux $*.out $*.nav` $(LOGDIR)
	rm -f $*.dvi $*.temp*

.f.o:
	$(FC) -c $(FOPTS) $(INC_PATHS) $< $(FAFTER)

.jpg.eps:
	convert $*.jpg $*.eps

.obj.eps:
	tgif -print -eps -page 1 -color $*.obj

.obj.pdf:
	tgif -print -pdf -page 1 -color $*.obj

all doc:

debug:
	echo "$(OMIT_SRCS)"
	echo "$(OMIT_PREFXD)"

-include Makefile.dep

all: doc cfm icfm

doc: ${NEAROBJS} ${OBJS} 
	date > doc; echo $(OBJS) >> doc

include $(CUSTOM)/$(MLSF95).$(MLSPLAT)

# ========================================================================

clean:
	rm -fr *.dvi *.pdf *.log *.aux *.bbl Makefile.dep $(DVIDIR) $(PDFDIR) $(LOGDIR) $(TEXDIR)

depends:
	$(UTILDIR)/makemakedep.sh -s tex -o dvi -p Makefile.dep.1 \
     $(OMIT_PREFXD) $(OMIT_PREFXD_REV) $(OMIT_PREFXD_SUB)
	sed 's/OBJS =/NEAROBJS =/' Makefile.dep.1 > Makefile.dep.2
	$(UTILDIR)/makemakedep.sh -s tex -o pdf -p Makefile.dep.1 \
     $(OMIT_PREFXD) $(OMIT_PREFXD_REV) $(OMIT_PREFXD_SUB)
	cat Makefile.dep.1 Makefile.dep.2 > Makefile.dep
	if [ ! -d $(TEXDIR) ] ; then \
	  		mkdir $(TEXDIR) ;\
			fi ; \
	if [ ! -d $(DVIDIR) ] ; then \
	  		mkdir $(DVIDIR) ;\
			fi ; \
	if [ ! -d $(PDFDIR) ] ; then \
	  		mkdir $(PDFDIR) ;\
			fi ; \
	if [ ! -d $(LOGDIR) ] ; then \
	  		mkdir $(LOGDIR) ;\
			fi
	 for dir in $(shell ${REECHO} -d cfm icfm l2cf); do \
	  $(MAKE) -f $(MakeFName) -C $$dir DOCDIR=$(DOCDIR)/$$dir depends; \
	done

help:
	@TOPIC=main; sed -n '/'$$TOPIC' help/,/End '$$TOPIC' help/ p' $(DOCDIR)/Makefile \
		| sed -n 's/^..//p' | sed '1 d; $$ d'

help--cfm:
	@TOPIC=ocfm; sed -n '/'$$TOPIC' help/,/End '$$TOPIC' help/ p' $(DOCDIR)/Makefile \
		| sed -n 's/^..//p' | sed '1 d; $$ d'

help--doc:
	@TOPIC=doc; sed -n '/'$$TOPIC' help/,/End '$$TOPIC' help/ p' $(DOCDIR)/Makefile \
		| sed -n 's/^..//p' | sed '1 d; $$ d'

help--icfm:
	@TOPIC=icfm; sed -n '/'$$TOPIC' help/,/End '$$TOPIC' help/ p' $(DOCDIR)/Makefile \
		| sed -n 's/^..//p' | sed '1 d; $$ d'

help--l2cf:
	@TOPIC=l2cf; sed -n '/'$$TOPIC' help/,/End '$$TOPIC' help/ p' $(DOCDIR)/Makefile \
		| sed -n 's/^..//p' | sed '1 d; $$ d'


cfm: $(DOCDIR)/cfm $(DOCDIR)/cfm/Makefile
	$(MAKE) -f $(MakeFName) -C $@ DOCDIR=$(DOCDIR)/cfm

icfm: $(DOCDIR)/icfm $(DOCDIR)/icfm/Makefile
	$(MAKE) -f $(MakeFName) -C $@ DOCDIR=$(DOCDIR)/icfm

l2:
	$(MAKE) -f $(MakeFName) -C $@  DOCDIR=$(DOCDIR)/l2

l2cf:
	$(MAKE) -f $(MakeFName) -C $@ DOCDIR=$(DOCDIR)/l2cf

tar:
	rm -f doc.tgz
	tar hcvzf docs.tgz * --exclude \*.dvi --exclude \*.pdf --exclude doc

update: 
	$(MAKE) -f $(MakeFName) depends

.PHONY: all clean depends internalupdate tar update l2 l2cf cfm icfm \

# ========================================================================
# $Log$
# Revision 1.12  2020/06/09 22:58:09  pwagner
# Attempt to fix problems caused by upgrading LaTeX and .nav files
#
# Revision 1.11  2020/04/21 00:08:57  pwagner
# externaldocument directed toward LOGDIR instead of cwd
#
# Revision 1.10  2019/09/19 22:19:05  pwagner
# Includes default build rules for .obj files
#
# Revision 1.9  2019/06/27 21:08:47  pwagner
# make clean more thorough
#
# Revision 1.8  2019/06/24 18:48:01  pwagner
# make clean is now more thorough
#
# Revision 1.7  2019/06/14 16:28:22  pwagner
# Seems to work better this way
#
# Revision 1.6  2019/05/31 20:05:13  pwagner
# May build eps filees from jpg (assuming convert is in path)
#
# Revision 1.5  2019/05/31 17:14:40  pwagner
# Fixed a problem preventing pdfs being built properly
#
# Revision 1.4  2018/09/20 22:53:21  pwagner
# Added new help--l2cf
#
# Revision 1.3  2018/08/07 20:36:48  pwagner
# Expanded list of .PHONY targets
#
# Revision 1.2  2018/07/11 16:45:18  pwagner
# Can now build l2cf, wvs-149
#
# Revision 1.1  2015/08/10 23:21:49  pwagner
# After switching away from MakeFC
#
# Revision 1.20  2014/11/21 21:20:00  pwagner
# OMIT_SRCS now got from tex.omitted
#
# Revision 1.19  2013/08/07 18:55:20  pwagner
# omits wvs-004.tex
#
# Revision 1.18  2012/09/12 00:38:15  pwagner
# Fixed ill-considered make updating subdirectories; added help; omoy wvs-111
#
# Revision 1.17  2012/04/20 22:39:18  pwagner
# Works around bugs in pdflatex href implementaion
#
# Revision 1.16  2012/04/04 16:31:52  pwagner
# Removed include of doc.h
#
# Revision 1.15  2012/04/04 00:47:01  pwagner
# Adds new builds for more wvs- sources
#
# Revision 1.14  2011/05/31 22:11:25  pwagner
# Added commands to build icfm
#
# Revision 1.13  2011/05/25 20:34:11  pwagner
# CONFDIR reverted to point to parent directory
#
# Revision 1.12  2010/12/07 00:22:21  pwagner
# Works correctly with current directory contents
#
# Revision 1.11  2010/07/14 20:54:39  pwagner
# Simple fix to problem redefining SUBDIRS
#
# Revision 1.10  2010/06/29 20:58:27  pwagner
# Must omit wvs-095.tex
#
# Revision 1.9  2010/06/03 23:32:47  pwagner
# Because wvs-083.tex needs a graphic
#
# Revision 1.8  2010/05/26 22:56:10  pwagner
# Should work better, even for index
#
# Revision 1.7  2010/02/26 02:01:49  pwagner
# Fixed bug in if statement
#
# Revision 1.6  2010/01/06 17:34:22  pwagner
# Should use MakeFC to make depends in SUBDIRS
#
# Revision 1.5  2010/01/05 19:24:29  pwagner
# Preserves Makefiles under subdirectories
#
# Revision 1.4  2009/12/16 19:33:44  pwagner
# Hope this fixes goldbrick
#
# Revision 1.3  2009/12/08 19:51:37  pwagner
# Working prior to adding cfm items
#
# Revision 1.2  2009/09/22 17:46:29  pwagner
# Can now ignore sources marked by revision number
#
# Revision 1.1  2009/09/18 19:47:31  pwagner
# First commit
#
