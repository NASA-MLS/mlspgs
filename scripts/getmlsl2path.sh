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
    mach* | soundbarrier )
	SOUNDBARRIER=1
	DESKTOP=0 ;;
esac

VERSION=v1.0.1
SIMULATION=s5
YEAR=1996
FWMVERSION=''

# Now work out what we have been asked for.
while [ -n "$(echo $1)" ]; do
    # First the 'global' arguments
    if [ -z "${1##--version=*}" ]; then
	VERSION=${1#--version=}
    fi
    if [ -z "${1##--simulation=*}" ]; then
	SIMULATION=${1#--simulation=}
    fi
    if [ -z "${1##--year=*}" ]; then
	YEAR=${1#--year=}
    fi
    if [ -z "${1##--fwmVersion=*}" ]; then
	FWMVERSION=${1#--fwmVersion=}
    fi
    if [ $1 == "--dao" ]; then
	if [ $SOUNDBARRIER == 1 ]; then
	    echo "/work3/livesey"
	else
	    echo "/data/dao/tsyn3d_mis_p/geos4/1996"
	fi
    fi
    if [ $1 == "--output" ]; then
	if [ $SOUNDBARRIER == 1 ]; then
	    echo "$HOME/$VERSION"
	else
	    echo "/bigdata/livesey/$VERSION"
	fi
    fi
    if [ $1 == "--l2pc" ]; then
	if [ $SOUNDBARRIER == 1 ]; then
	    echo "/work3/livesey/"
	else
	    echo "/bigdata/livesey/$VERSION"
	fi
    fi
    if [ $1 == "--truthl2gp" ]; then
	if [ $SOUNDBARRIER == 1 ]; then
	    echo "$HOME/$VERSION"
	else
	    echo "/data/emls/l2gp/$SIMULATION--t/1996"
	fi
    fi
    if [ $1 == "--corel2gp" ]; then
	if [ $SOUNDBARRIER == 1 ]; then
	    echo "$HOME/$VERSION"
	else
	    echo "/bigdata/livesey/$VERSION"
	fi
    fi
    if [ $1 == "--l1boa" ]; then
	if [ $SOUNDBARRIER == 1 ]; then
	    echo "$HOME/$VERSION"
	else
	    echo "/data/emls/l1boa/$SIMULATION--t/$YEAR"
	fi
    fi
    if [ $1 == "--l1brad" ]; then
	if [ $SOUNDBARRIER == 1 ]; then
	    echo "$HOME/$VERSION"
	else
	    echo "/data/emls/l1brad/$SIMULATION/$FWMVERSION"
	fi
    fi
    if [ $1 == "--l2cal" ]; then
	if [ $SOUNDBARRIER == 1 ]; then
	    echo "$HOME/$VERSION"
	else
	    echo "/data/emls/l2cal"
	fi
    fi
    if [ $1 == "--leapsec" ]; then
	echo "$HOME/emls/$VERSION"
    fi
    shift
done

# $Log$
# Revision 1.11  2003/05/10 20:35:41  livesey
# Added DAO path.
#
# Revision 1.10  2003/05/09 23:17:10  livesey
# Bug fix.
#
# Revision 1.9  2003/05/08 22:57:53  livesey
# Tried to make it more generic
#
# Revision 1.8  2003/05/08 01:03:20  livesey
# Changed to s5--t
#
# Revision 1.7  2003/05/07 20:54:53  livesey
# Moved to v1.1.9
#
# Revision 1.6  2003/05/06 00:54:21  livesey
# Moved DAO data for the moment.
#
# Revision 1.5  2003/05/03 22:01:11  livesey
# Added DAO
#
# Revision 1.4  2003/02/19 21:58:01  livesey
# Moved l2pc files on soundbarrier
#
# Revision 1.3  2002/12/10 02:19:52  livesey
# New soundbarrier configuration
#
# Revision 1.2  2002/11/27 18:26:48  livesey
# Got rid of the distinction between sums and the rest of the scf.
# Probably don't need to worry about the l2pcs as much any more.
#
# Revision 1.1  2002/11/21 17:30:53  livesey
# First version
#
