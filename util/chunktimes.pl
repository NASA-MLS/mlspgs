#!/usr/local/bin/perl
# chunktimes.pl
# Usage:
# chunktimes.pl mlsl2.log
# where mlsl2.log is the catenation of all the chunks' logs from a run
# in which each is mediated by a PCF and the switches include "time"
# so that each chunk outputs its own timing summary
# and phases are marked by special phase-setting Fill sections
#
# Output:
# A table of chunk timings similar to the following:
#chunk initptan updateptan inituth core coreplusr2 highcloud coreplusr3 coreplusr4 coreplusr5 (final)
#99 788.45 1288.13 2767.67 785.33 16979.92 2603.82 18267.70 11243.27 7146.03 68051.70
#276 656.02 1326.15 2349.71 783.88 30826.25 2608.30 14016.60 10608.34 6788.58 76286.80
#10 637.48 1267.51 2386.22 756.75 17297.47 2623.50 18777.74 11160.70 7106.07 68217.81
#100 604.73 1393.60 3215.92 788.31 30573.74 2614.37 13866.33 10129.87 6612.55 76165.96
#    .   .   .
#
# -----------------------------------------------------------------------
#      Options
# -nohead           skip writing header line to output
# -headonly         write only header line to output
# -head list        head table by comma-separated list of phases
#                     (instead of the default)
# -s2h              convert times from seconds to hours
# -h2s              convert times from hours to seconds
# Bugs and limitations:
# (1) Should be able to handle toolkitless runs
#      would need to change from $_[7] to $_[5] and
#      from $_[3] to $_[1]
# (2) list of phase names must be lower case
# (3) Wouldn't you like to be able to sort these by, e.g., (final) or chunk?
#      would need to master idea of perl references; see @results array
#      (got past this by letting chunktimes.sh call this perl script)
# Copyright (c) 2004, California Institute of Technology.  ALL RIGHTS RESERVED.
# U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

# "$Id$"
#
my $chunk;
my $headList;
my $headonly;
my $ishead;
my $lastChunk;
my $more_opts;
my $nohead;
my $numChunks;
my $time;
my $timeConvert;
my @results;
my $Sp;
my @Sps;
my %times;
#An array of possible phase names
#@Sps = split(',', "initptan,updateptan,inituth,core,coreplusr2,highcloud,coreplusr3,coreplusr4,coreplusr5,(final)");
#         * * * * * * * * * * * * * * * * * * * * * * * 
# Check command-line args for options -s etc.
$headList = "initptan,updateptan,inituth,core,coreplusr2,highcloud,coreplusr3,coreplusr4,coreplusr5,(final)";
$headonly = 0;
$nohead = 0;
$timeConvert = 1;
$more_opts = TRUE;
while ($more_opts) {
   if ($ARGV[0] =~ /^-nohead/) {
      $nohead = 1;
      shift;
   } elsif ($ARGV[0] =~ /^-headonly/) {
      $headonly = 1;
      shift;
   } elsif ($ARGV[0] =~ /^-head/) {
      $headList = $ARGV[1];
      shift;
      shift;
   } elsif ($ARGV[0] =~ /^-s2h/) {
      $timeConvert = 3600;
      shift;
   } elsif ($ARGV[0] =~ /^-h2s/) {
      $timeConvert = 1./3600;
      shift;
   } else {
      $more_opts = 0;
   }
}
@Sps = split(',', $headList);
for $Sp (@Sps) {
  $times{"$Sp"} = "$Sp";
  }
$times{chunk} = "chunk";
$results[0] = %times;
$lastChunk = 0;
# print $times{chunk}, "\n";
# print $times{initptan}, "\n";
# @keys = keys(%times);
# @values = values(%times);
# print @keys, "\n";
# print @values, "\n";
$ishead = TRUE;
while (<>) {
  if (/\(Slave/) {
    chomp;
    split;
    # print "After split: @_ \n";
    $chunk = $_[7];
    if ( $chunk eq "" ) {
      print $chunk;
    }
    elsif ( $chunk ne $lastChunk ) {
      push(@results, %times);
      # print "lastChunk: $lastChunk\n";
      # print "chunk: $chunk\n";
      # print values(%times), "\n";
      # for $Sp (values(%times)) {
      if (! ($ishead && $nohead)) {
        print "$times{chunk} ";
        for $Sp (@Sps) {
          if ( $ishead ) {
            print "$times{$Sp} ";
          } else {
            $time = $times{$Sp} / $timeConvert;
            # print "$time ";
            # printf("%.2f %s", $time, ' ');
            &PrintTime($time);
          }
        }
        print "\n";
      }
      if ($ishead && $headonly) {
        exit;
      }
      $lastChunk = $chunk;
      $times{chunk} = $chunk;
      $ishead = 0;
      if (/\(final\)/) {
        $times{"(final)"} = $_[3];
        # print "field: (final) set to $_[3] \n";
        }
    }
    $indx=0;
    tr/A-Z/a-z/;
    for $Sp (@Sps) {
      if (/\s$Sp\s/) {
        $times{$Sp} = $_[3];
        # print "field: $Sp set to $_[3] \n";
        # $indx++;
        }
    }
   }
}
#
# &PrintTime(time);
# --- print time nicely, with appropriate number of digits
#
sub PrintTime {
   local($ptime) = shift(@_);
   if ( $ptime < 0.0001 ) {
     printf("%.8f %s", $ptime, ' ');
   } elsif ( $ptime < 0.01 ) {
     printf("%.6f %s", $ptime, ' ');
   } elsif ( $ptime < 1. ) {
     printf("%.4f %s", $ptime, ' ');
   } elsif ( $ptime < 100. ) {
     printf("%.2f %s", $ptime, ' ');
#   } elsif ( $ptime < 10000. ) {
#     printf("%.0f %s", $ptime, ' ');
   } else {
     print "$ptime ";
   }
}
# $Log$
# Revision 1.1  2004/05/13 22:51:58  pwagner
# First commit
#
