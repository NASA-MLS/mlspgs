#	Makefile
#
#  Lost? Type 'make help -f uptodate.make' from the directory where you see this file
#
# Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
# U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

# ---------------------- The name of this file (is the right-hand-side below)
MakeFName = uptodate.make
# All the implicit pattern-matching involving suffixes
# are irrelevant to us. We make our own rules.
.SUFFIXES:

SHELL = /bin/sh

MLSBIN=./util
REECHO=$(MLSBIN)/reecho.sh

# --------------- uptodate.make help
# Usage: make -f uptodate.make SOURCE=sourcefile TARGET=targetfile targetfile
# Purpose: keep targetfile up-to-date with sourcefile
# If target's modification date is newer or the same as sourcefile's--do nothing
# Otherwise, targetfile is presumed out-of-date--cp sourcefile over it

# --------------- End uptodate.make help
help:
	@sed -n '/'$(MakeFName)' help/,/End '$(MakeFName)' help/ p' $(MLSBIN)/uptodate.make \
		| sed -n 's/^..//p' | sed '1 d; $$ d'

# These are just dummy names
SOURCE=sourcefile
TARGET=targetfile

$(TARGET): $(SOURCE)
	cp $(SOURCE) $(TARGET)

.PHONY: help directory
#
# An obvious improvement would be
# a special target named 'directory'
# So that if given say
#  make -f uptodate.make SOURCEDIR=sourcedir TARGETDIR=targetdir directory
# we would form a list of all the files in sourcedir
# and check that all the ones in targetdir were up-to-date with them
# This would proceed probably along the following lines
#EVERY_FILE = $(shell ${REECHO} -dirn $(SOURCEDIR))
#directory:
#	for file in $(EVERY_FILE); do\
#	   $(MAKE) -f $(MakeFName) SOURCE=$(SOURCEDIR)/$$file \
#      TARGET=$(TARGETDIR)/$$file $(TARGETDIR)/$$file ;\
#	done

# $Log$

