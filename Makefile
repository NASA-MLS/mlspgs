# Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
# U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

# $Id$

SHELL=/bin/sh

# Note: From the point of view of "make", the IF-ELSE-FI block has to be all
# one rule.  This is because "make" sends each rule through a separate
# invocation of the shell.

all:
	@if [ ".$(SYS)" = "." ]; then \
          echo "Specify a compiler/platform directory using the SYS variable."; \
          echo "E.g. make SYS=NAG"; \
        else \
          cd $(HOME)/mlspgs/lib/$(SYS); make all; \
          cd $(HOME)/mlspgs/l2/$(SYS); make all; \
        fi

clean:
	@if [ ".$(SYS)" = "." ]; then \
          echo "Specify a compiler/platform directory using the SYS variable."; \
          echo "E.g. make SYS=NAG"; \
        else \
          cd $(HOME)/mlspgs/lib/$(SYS); make clean; \
          cd $(HOME)/mlspgs/l2/$(SYS); make clean; \
        fi

# $Log$
# Revision 1.1  2000/10/10 20:44:50  vsnyder
# Initial entry
#
