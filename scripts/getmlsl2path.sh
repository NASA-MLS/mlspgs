#! /bin/bash
# $Id$
# This shell script prints the name of the desired directory for a given file type
# First work out where we are.
MACHINE=$(uname -n)
SUMS=0
SOUNDBARRIER=0
DESKTOP=1
case "$MACHINE" in
    headnode | sum* )
	SUMS=1
	DESKTOP=0 ;;
    mach* )
	SOUNDBARRIER=1
	DESKTOP=0 ;;
esac

# Now work out what we have been asked for.
while [ -n "$(echo $1)" ]; do
    if [ $1 == "--output" ]; then
	if [ $SOUNDBARRIER == 1 ]; then
	    echo "$HOME/v1.0.1"
	else
	    echo "/bigdata/livesey/v1.0.1"
	fi
    fi
    if [ $1 == "--l2pc" ]; then
	if [ $SUMS == 1 ]; then
	    echo "/scratch/livesey"
	fi
	if [ $SOUNDBARRIER == 1 ]; then
	    echo "$HOME/v1.0.1"
	fi
	if [ $DESKTOP == 1 ]; then
	    echo "/bigdata/livesey/v1.0.1"
	fi
    fi
    if [ $1 == "--truthl2gp" ]; then
	if [ $SOUNDBARRIER == 1 ]; then
	    echo "$HOME/v1.0.1"
	else
	    echo "/data/emls/l2gp/s5--t/1996"
	fi
    fi
    if [ $1 == "--corel2gp" ]; then
	if [ $SOUNDBARRIER == 1 ]; then
	    echo "$HOME/v1.0.1"
	else
	    echo "/bigdata/livesey/v1.0.1"
	fi
    fi
    if [ $1 == "--l1boa" ]; then
	if [ $SOUNDBARRIER == 1 ]; then
	    echo "$HOME/v1.0.1"
	else
	    echo "/data/emls/l1boa/s5--t/1996"
	fi
    fi
    if [ $1 == "--l1brad" ]; then
	if [ $SOUNDBARRIER == 1 ]; then
	    echo "$HOME/v1.0.1"
	else
	    echo "/bigdata/livesey/v1.0.1"
	fi
    fi
    if [ $1 == "--l2cal" ]; then
	if [ $SOUNDBARRIER == 1 ]; then
	    echo "$HOME/v1.0.1"
	else
	    echo "/data/emls/l2cal"
	fi
    fi
    if [ $1 == "--leapsec" ]; then
	echo "$HOME/emls/v1.0.1"
    fi
    shift
done

# $Log$
