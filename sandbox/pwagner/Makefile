# Makefile
# testing a few ideas before incorporating them into mlspgs

# where to store .configure
CONFDIR=.

# This will set our MLS platform, compiler. etc.
# w/o complaining if it does nor exist yet
-include $(CONFDIR)/.configure

SHELL = /bin/sh

# where a.sh is stored
MLSBIN=.

# The program will be simply a file containing the value of a
PROG=value

.PHONY: configure depends clean tar clean_objs

$(PROG): $(CONFDIR)/.configure
	rm -f $(PROG)
	echo "$(a)" > $(PROG)
	@case "$(a)" in \
		"") \
		echo "a is blank" ; \
		;; \
		*) \
		echo "a is non-blank" ; \
		;; \
	esac 

configure:
	echo "Currently, a is $$a"
	$(MLSBIN)/a.sh
	echo "However, a is still $(a)"
	echo "and, buck-a is still $$a"

$(CONFDIR)/.configure:
	echo "Currently, a is $$a"
	$(MLSBIN)/a.sh
	echo "However, a is still $(a)"
	echo "and, buck-a is still $$a"

#
	
 

