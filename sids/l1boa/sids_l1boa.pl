#!/usr/local/bin/perl

# program:	sids_l1boa.pl

# description:	run the SIDS L1BOA program in the Toolkit environment

$cdir = `pwd`;
chomp($cdir);

print "Running sids_l1boa.pl ...\n\n";

print "using User Input File $cdir/sids_uif.txt ...\n\n";

print "executing SIDS L1BOA program ...\n\n";

print "Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.\n";

print "U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.\n\n";

#Log files suppressed

`PGS_PC_Shell.sh /user2/nakamura/L1BOA/sids_l1boa 1111 /user2/nakamura/L1BOA/PCF.sids 50`;

print "sids_l1boa.pl finished\n";

# end of sids_l1boa.pl script
