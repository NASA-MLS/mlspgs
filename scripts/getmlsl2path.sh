#! /bin/bash
# $Id$
# This shell script prints the name of the desired directory for a given file type
# First work out where we are.
if [ -f ~/machine.id ]
then
  MACHINE=$(cat ~/machine.id)
else
  MACHINE="desktop"
fi
SUMS=0

VERSION=v1.0.1
SIMULATION=s5
YEAR=1996
FWMVERSION=''
MYUSER=$USER

# Now work out what we have been asked for.
while [ -n "$(echo $1)" ]; do
    # First the 'global' arguments
    if [ -z "${1##--user=*}" ]; then
	MYUSER=${1#--user=}
    fi
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
	case $MACHINE in
	    soundbarrier | scramjet )
		echo "/research1/$MYUSER";;
	    lightspeed )
		echo "$HOME/$VERSION";;
	    desktop )
		if [ $YEAR == "1996" ]; then
		    echo "/data/dao/tsyn3d_mis_p/geos4/$YEAR"
		else
		    echo "/data/dao/D4FAPMIS/$YEAR"
		fi
	esac
    fi
    if [ $1 == "--output" ]; then
	case $MACHINE in
	    soundbarrier | scramjet)
		echo "/research1/$MYUSER";;
	    lightspeed )
		echo "$HOME/$VERSION";;
	    desktop )
		echo "/bigdata/$MYUSER/$VERSION"
	esac
    fi
    if [ $1 == "--l2pc" ]; then
	case $MACHINE in
	    soundbarrier )
		echo "/research1/$MYUSER/";;
	    scramjet )
		echo "/state/partition1/livesey/";;
	    lightspeed )
		echo "/science/l2cal/";;
	    desktop )
		echo "/bigdata/$MYUSER/$VERSION";;
	esac
    fi
    if [ $1 == "--truthl2gp" ]; then
	case $MACHINE in
	    soundbarrier | scramjet)
		echo "/research1/$MYUSER/";;
	    lightspeed )
		echo "$HOME/$VERSION";;
	    desktop )
		echo "/testing/emls/l2gp/$SIMULATION--t/$YEAR";;
	esac
    fi
    if [ $1 == "--corel2gp" ]; then
	case $MACHINE in
	    soundbarrier | scramjet)
		echo "/research1/$MYUSER/";;
	    lightspeed )
		echo "$HOME/$VERSION";;
	    desktop )
	    echo "/bigdata/$MYUSER/$VERSION";;
	esac
    fi
    if [ $1 == "--l1boa" ]; then
	case $MACHINE in
	    soundbarrier | scramjet)
		echo "/research1/$MYUSER/";;
	    lightspeed )
		echo "$HOME/$VERSION";;
	    desktop )
		echo "/testing/emls/l1boa/$SIMULATION--t/$YEAR";;
	esac
    fi
    if [ $1 == "--l1brad" ]; then
	case $MACHINE in
	    soundbarrier | scramjet)
		echo "/research1/$MYUSER/";;
	    lightspeed )
		echo "$HOME/$VERSION";;
	    desktop )
		echo "/testing/emls/l1brad/$SIMULATION/$FWMVERSION";;
	esac
    fi
    if [ $1 == "--l2cal" ]; then
	case $MACHINE in
	    soundbarrier | scramjet)
		echo "/research1/$MYUSER/";;
	    lightspeed )
		echo "$HOME/$VERSION";;
	    desktop )
		echo "/data/emls/l2cal";;
	esac
    fi
    if [ $1 == "--leapsec" ]; then
	echo "$HOME/emls/$VERSION"
    fi
    if [ $1 == "--tmp" ]; then
	case $MACHINE in
	    soundbarrier | scramjet)
		echo "/research1/$MYUSER/";;
	    lightspeed )
		echo "$HOME/$VERSION";;
	    desktop )
		echo "/bigdata/$MYUSER/$VERSION";;
	esac
    fi
    if [ $1 == "--sidsrad" ]; then
	case $MACHINE in
	    soundbarrier | scramjet)
		echo "/research1/$MYUSER/";;
	    lightspeed )
		echo "$HOME/$VERSION";;
	    desktop )
		echo "/bigdata/$MYUSER/$VERSION";;
	esac
    fi
    shift
done

# $Log$
# Revision 1.31  2004/05/01 00:22:48  livesey
# Moved l2pcfiles on lightspeed
#
# Revision 1.30  2004/04/17 05:58:37  livesey
# Bug fix
#
# Revision 1.29  2004/04/17 05:50:03  livesey
# Moved the l2pcs on scramjet
#
# Revision 1.28  2004/04/08 17:15:17  livesey
# Bug fix
#
# Revision 1.27  2004/04/07 23:51:19  jonathan
# No longer requires ~/machine.id
#
# Revision 1.26  2004/03/10 22:36:41  livesey
# Added scramjet
#
# Revision 1.25  2004/02/09 19:23:03  livesey
# Made more lightspeed ready, gets machine name from ~/machine.id file.
#
# Revision 1.24  2004/02/09 18:47:05  livesey
# Typo!
#
# Revision 1.23  2004/02/09 18:37:58  livesey
# Added lightspeed.
#
# Revision 1.22  2004/01/22 21:15:45  livesey
# Bug fix from previous change
#
# Revision 1.21  2004/01/22 21:15:18  livesey
# Added ability to overwrite user.
#
# Revision 1.20  2003/11/27 01:28:22  livesey
# Bug fix
#
# Revision 1.19  2003/11/27 00:31:49  livesey
# Moved stuff for soundbarrier
#
# Revision 1.18  2003/11/14 21:25:45  livesey
# A little more (irrelevant as it turns out) intelligence for the dao
# data.
#
# Revision 1.17  2003/11/07 01:08:34  livesey
# More changes for soundbarrier
#
# Revision 1.16  2003/08/28 00:43:31  livesey
# Changed work3 to work1
#
# Revision 1.15  2003/05/22 17:41:27  livesey
# Added sids rad option
#
# Revision 1.14  2003/05/13 04:13:28  livesey
# Changed livesey's to $USERS and added --tmp option
#
# Revision 1.13  2003/05/10 22:23:44  livesey
# Tidied up year related issues
#
# Revision 1.12  2003/05/10 20:39:34  livesey
# Bug fix for dao
#
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
