#!/usr/bin/perl -w
# This is a program that generates an EOS MLS Level 1 fragment
# PCF file for the requested day.
# i.e., just the ML0ENG{1-6} and ML0SCI{1-2} and AURATTH and AUREPHMH files
# Author Ryan Fuller; May 11, 2004
# Author Vince Perun; October 25, 2006
# Hacker Anonymous Biped; Sep 7, 2009
use Cwd;
use Date::Manip;
use DBI;
use File::Path;

my($line);

print "Date to Process (YYYY-DDD): ";
chomp($line = <STDIN>);
unless($line =~ /[0-9]{4}[d-][0-9]{3}/) {
    die 'Day must be in the form "YYYY-DDD"';
}
my($dayString) = $line;
$dayString =~ /[d-]/;
my($doy) = $';
my($year, $month, $day) = &Date_NthDayOfYear($`, $');

my($inVersion) = "";
#print 'Version of input data (default="'. $inVersion . '"): ';
#chomp($line = <STDIN>);
#if($line ne "") {
#    $inVersion = $line;
#}

#my($inDir) = cwd . "/inputs/";
#print "Default input directory (default=$inDir): ";
#chomp($line = <STDIN>);
#if($line ne "") {
#    $inDir = $line;
#}

#my($outDir) = cwd . "/outputs_$dayString/";
#print "Default output directory (default=$outDir): ";
#chomp($line = <STDIN>);
#if($line ne "") {
#    $outDir = $line;
#}
my($outVersion) = "v02-10";
#print "Version of the resulting Level 1 Files (default=$outVersion): ";
#chomp($line = <STDIN>);
#if($line ne "") {
#    $outVersion = $line;
#}
my($cycle) = 1;
#print "Cycle of the resulting Level 1 Files (default=$cycle): ";
#chomp($line = <STDIN>);
#if($line ne "") {
#    $cycle = $line;
#}
#
#my($pcfFile) = "${dayString}_";
#if($inVersion ne "") {
#    $pcfFile .= "${inVersion}_";
#}
#$pcfFile .= "${outVersion}-c" . sprintf("%02d", $cycle) . ".pcf";
#print "Name of PCF file to be created (default=$pcfFile): ";
#chomp($line = <STDIN>);
#if($line ne "") {
#    $pcfFile = $line;
#}
#
#my($cfFile) = "l1-v210.cf";
#print "The configuration file to use (default=$cfFile): ";
#chomp($line = <STDIN>);
#if($line ne "") {
#    $cfFile = $line;
#}
#
#my($cfDir) = cwd;
#print "The location of the cf and pcf files (default=$cfDir): ";
#chomp($line = <STDIN>);
#if($line ne "") {
#    $cfDir = $line;
#}
#
#my($toolkit) = "SCF  TK5.2.11.A";
#print 'Toolkit Version String (default="' . $toolkit . '"): ';
#chomp($line = <STDIN>);
#if($line ne "") {
#    $toolkit = $line;
#}

# Can get into trouble when adding/subtracting a day over the boundary
# between DST/non-DST, e.g. subtracting 1 day from 2017/03/13 00:00:00
# returns 2017/03/11 23:00:00, because of the change from PST to PDT
# on 2013/03/12 02:00:00. Fix this by requiring all calculations to
# occur in UTC with a DSTFlag=0 (i.e. no Daylight Saving Times
# modifications.

$dateStr="$month/$day/$year 00:00:00 GMT";

my $date = new Date::Manip::Date;
my $delta=new Date::Manip::Delta;
$date->parse($dateStr);

# Calculate strings for today - 1
$delta->parse("-1 day");
print "delta is defined" if $delta;
print "date is defined" if $date;
my $d=$date->calc($delta);
my $pdsStartTime=$d->printf("%Y-%m-%d 22:00:00");
my $prevDay=$d->printf("%Y-%m-%d 00:00:00");


# Calculate strings for today + 1
$delta->parse("+1 day");
$d=$date->calc($delta);
my $pdsEndTime = $d->printf("%Y-%m-%d 02:00:00");
my $nextDay = $d->printf("%Y-%m-%d 23:59:59");

my $dayStart=$date->printf("%Y-%m-%d 00:00:00");
my $dayEnd=$date->printf("%Y-%m-%d 23:59:59");






# Open a connection to the database - Is global
# You will need import your database connection
# credentials here
$db = "";
$connectString = "";
($username, $password) = ("", "");
$dbh = DBI->connect($connectString, $username, $password);

# Now go through each APID and collect the available files over the 
# desired time range
my($i);
my(@eng1) = &GetPDSFiles("ML0ENG1", $inVersion, $pdsStartTime, $pdsEndTime);
my(@eng2) = &GetPDSFiles("ML0ENG2", $inVersion, $pdsStartTime, $pdsEndTime);
my(@eng3) = &GetPDSFiles("ML0ENG3", $inVersion, $pdsStartTime, $pdsEndTime);
my(@eng4) = &GetPDSFiles("ML0ENG4", $inVersion, $pdsStartTime, $pdsEndTime);
my(@eng5) = &GetPDSFiles("ML0ENG5", $inVersion, $pdsStartTime, $pdsEndTime);
my(@eng6) = &GetPDSFiles("ML0ENG6", $inVersion, $pdsStartTime, $pdsEndTime);
my(@sci1) = &GetPDSFiles("ML0SCI1", $inVersion, $pdsStartTime, $pdsEndTime);
my(@sci2) = &GetPDSFiles("ML0SCI2", $inVersion, $pdsStartTime, $pdsEndTime);
#my(@att) =  &GetAttFiles("AURATTH", $inVersion, $dayStart, $dayEnd);
my(@att) =  &GetAttFiles("AURATTH", $inVersion, $pdsStartTime, $pdsEndTime);
my(@eph) =  &GetEphFiles("AUREPHMH", $inVersion, $prevDay, $nextDay);

# Now get the leapsec/utcpole files
my($leapsecDir) = "~/database/common/TD";
#print "Location of file leapsec.dat (default=$leapsecDir, in toolkit): ";
#chomp($line = <STDIN>);
#if($line ne "") {
#    $leapsecDir = $line;
#}
my($utcpoleDir) = "~/database/common/CSC";
#print "Location of file utcpole.dat (default=$utcpoleDir, in toolkit): ";
#chomp($line = <STDIN>);
#if($line ne "") {
#    $utcpoleDir = $line;
#}


# Disconnect from the database and exit
$dbh->disconnect();

# Now we know we can make the PCF file, we need to go through the long and arduous 
# process of outputting it to a file.
$pcfFile = "PCF.fragment";
unless(open PCFOUT, "> $pcfFile") {
    die "cannot open $pcfFile: $!";
}

my($j, $k, @thisFileType, $pcfId, $pcfIdString, @lgIds, $shortName);
for($i = 0; $i < 10; $i++) {
    if($i == 0) {
	@thisFileType = @eng1;
	$pcfId = 21000;
	$shortName = "ML0ENG1";
	@lgIds = ("16899", "16902", "16911", "16912", "16937", "16939", "16957", 
		  "16956", "16989", "16988", "17006", "17005", "17072", "17071", 
		  "17093", "17092", "17115", "17114", "17128", "17127", "17145", 
		  "17144", "17163", "17162", "17216", "17215", "17200", "17201");
    } elsif($i == 1) {
	@thisFileType = @eng2;
	$pcfId = 21020;
	$shortName = "ML0ENG2";
	@lgIds = ("16895", "16896", "16916", "16915", "16933", "16932", "16964",
		  "16960", "16995", "16994", "17004", "17002", "17070", "17069",
		  "17095", "17094", "17111", "17110", "17130", "17129", "17149",
		  "17148", "17167", "17165", "17185", "17183", "17207", "17209");
    } elsif($i == 2) {
	@thisFileType = @eng3;
	$pcfId = 21040;
	$shortName = "ML0ENG3";
	@lgIds = ("16904", "16903", "16913", "16914", "16934", "16935", "16969",
		  "16968", "16986", "16987", "17020", "17019", "17076", "17075", 
		  "17087", "17086", "17119", "17117", "17132", "17131", "17155",
		  "17154", "17166", "17164", "17180", "17181", "17196", "17198");
    } elsif($i == 3) {
	@thisFileType = @eng4;
	$pcfId = 21060;
	$shortName = "ML0ENG4";
	@lgIds = ("16906", "16905", "16918", "16817", "16938", "16936", "16962", 
		  "16966", "16993", "16992", "17008", "17009", "17074", "17073",
		  "17099", "17098", "17107", "17106", "17135", "17136", "17147",
		  "17146", "17170", "17171", "17189", "17188", "17199", "17197");
    } elsif($i == 4) {
	@thisFileType = @eng5;
	$pcfId = 21080;
	$shortName = "ML0ENG5";
	@lgIds = ("16893", "16894", "16920", "16919", "16944", "16945", "16967",
		  "16965", "16984", "16985", "17012", "17011", "17077", "17078",
		  "17097", "17096", "17112", "17113", "17134", "17133", "17151",
		  "17150", "17169", "17168", "17179", "17178", "17203", "17202");
    } elsif($i == 5) {
	@thisFileType = @eng6;
	$pcfId = 21100;
	$shortName = "ML0ENG6";
	@lgIds = ("16900", "16901", "16922", "16921", "16943", "16942", "16959",
		  "16961", "16991", "16990", "17007", "17010", "17080", "17079",
		  "17091", "17089", "17109", "17108", "17138", "17137", "17194",
		  "17195", "17173", "17172", "17186", "17187", "17205", "17204");
    } elsif($i == 6) {
	@thisFileType = @sci1;
	$pcfId = 21120;
	$shortName = "ML0SCI1";
	@lgIds = ("16909", "16908", "16926", "16925", "16946", "16947", "16973", 
		  "16972", "16999", "16998", "17014", "17013", "17082", "17081", 
		  "17102", "17100", "17122", "17123", "17140", "17139", "17157", 
		  "17156", "17177", "17176", "17193", "17192", "17213", "17212");
    } elsif($i == 7) {
	@thisFileType = @sci2;
	$pcfId = 21140;
	$shortName = "ML0SCI2";
	@lgIds = ("16910", "16907", "16924", "16923", "16949", "16948", "16971", 
		  "16970", "16996", "16997", "17016", "17015", "17084", "17083", 
		  "17106", "17101", "17121", "17120", "17142", "17141", "17159", 
		  "17158", "17175", "17174", "17190", "17191", "17211", "17210");
    } elsif($i == 8) {
	@thisFileType = @eph;
	$pcfId = 10501;
	$shortName = "AUREPHMH";
	@lgIds = ("16975", "17283");
    } elsif($i == 9) {
	@thisFileType = @att;
	$pcfId = 10502;
	$shortName = "AURATTH";
	@lgIds = ("16981", "16982", "17003", "17085", "17104", "17124", 
		  "17284", "17284", "17288", "17285", "17287", "17286");
    }

    $k = 0;
    for($j = 0; $j <= $#thisFileType; $j++) {
	#$pcfIdString = sprintf("%d", $pcfId);
	if($thisFileType[$j]->{present} ne "No") {
	    if($i <= 7) {
		print PCFOUT "$pcfId|$thisFileType[$j]->{auxProductName}|" .
              	 	     "$thisFileType[$j]->{path}|||" . "|2\n";
			     #"LGID:$shortName:2:$lgIds[$k]|2\n";
		$k++;
		print PCFOUT "$pcfId|$thisFileType[$j]->{productName}|" .
		             "$thisFileType[$j]->{path}|||" . "|1\n";
			     #"LGID:$shortName:2:$lgIds[$k]|1\n";
		$k++;
		$pcfId++;
	    } else {
		$k = $#thisFileType - $j + 1;
		print PCFOUT "$pcfId|$thisFileType[$j]->{productName}|" .
		             "$thisFileType[$j]->{path}|||" . "|$k\n";
			     #"LGID:$shortName:2:$lgIds[$j]|$k\n";
	    }
	}
    }
}

# Now this part requires a bit of formatting
my($outVersionString) = $outVersion;
$outVersionString =~ s/[.]/-/g;
$outVersionString = sprintf("$outVersionString-c%02d_%04dd%03d", $cycle, $year, $doy);

close PCFOUT;



exit;


sub GetPDSFiles {
    
    my($shortName) = $_[0];
    my($version) = $_[1];
    my($startTime) = $_[2];
    my($endTime) = $_[3];

    my($i, $sql, $sth, $row, @rVal, $thisTime);

    # Initialize the hash reference
    for($i = 0; $i < 14; $i++) {
	$rVal[$i] = {   present => "No",
			productName => "",
			auxProductName => "",
			path => "",
			startTime => "",
			endTime => "",
		    };
    }

    $sql = 'SELECT * FROM scfFromSips WHERE startTime >= "' . $startTime . '" ' .
	   'AND endTime <= "' . $endTime . '" AND  shortName="' . $shortName . '" ';
    if($version ne "") {
	$sql .= 'AND version="' . $version . '" '; 
    }
    $sql .= 'ORDER BY addTime ASC';
    $sth = $dbh->prepare($sql);
    $sth->execute or die "Unable to execute query: $dbh->errstr\n";

    while($row = $sth->fetchrow_hashref) {

	# Figure out what index the result fits into
	$row->{startTime} =~ /\s[0-9]{2}:/;
	$thisTime = $&;
	$thisTime =~ s/[\s:]//g;
	if($thisTime eq "22" || $thisTime eq "23") {
	    # Now see what day
	    $row->{startTime} =~ /\s/;
	    $thisTime = $`;
	    $startTime =~ /\s/;
	    $i = $thisTime eq $` ? 0 : 12;
	} elsif($thisTime eq "00" || $thisTime eq "01") {
	    # Now see what day
	    $row->{endTime} =~ /\s/;
	    $thisTime = $`;
	    $endTime =~ /\s/;
	    $i = $thisTime eq $` ? 13 : 1;
	} elsif($thisTime eq "02" || $thisTime eq "03") {
	    $i = 2;
	} elsif($thisTime eq "04" || $thisTime eq "05") {
	    $i = 3;
	} elsif($thisTime eq "06" || $thisTime eq "07") {
	    $i = 4;
	} elsif($thisTime eq "08" || $thisTime eq "09") {
	    $i = 5;
	} elsif($thisTime eq "10" || $thisTime eq "11") {
	    $i = 6;
	} elsif($thisTime eq "12" || $thisTime eq "13") {
	    $i = 7;
	} elsif($thisTime eq "14" || $thisTime eq "15") {
	    $i = 8;
	} elsif($thisTime eq "16" || $thisTime eq "17") {
	    $i = 9;
	} elsif($thisTime eq "18" || $thisTime eq "19") {
	    $i = 10;
	} elsif($thisTime eq "20" || $thisTime eq "21") {
	    $i = 11;
	} else {
	    die "Uh Oh $thisTime $row->{startTime}\n";
	}
	$rVal[$i]{present} = "Yes";
	$rVal[$i]{productName} = $row->{productName};
	$rVal[$i]{auxProductName} = $row->{auxProductName};
	$rVal[$i]{path} = $row->{path};
	$rVal[$i]{startTime} = $row->{startTime};
	$rVal[$i]{endTime} = $row->{endTime};
    }

    # Figure out APID from Short Name
    my($apid);
    if($shortName eq "ML0ENG1") {
	$apid = "1732";
    } elsif($shortName eq "ML0ENG2") {
	$apid = "1734";
    } elsif($shortName eq "ML0ENG3") {
	$apid = "1736";
    } elsif($shortName eq "ML0ENG4") {
	$apid = "1738";
    } elsif($shortName eq "ML0ENG5") {
	$apid = "1740";
    } elsif($shortName eq "ML0ENG6") {
	$apid = "1742";
    } elsif($shortName eq "ML0SCI1") {
	$apid = "1744";
    } elsif($shortName eq "ML0SCI2") {
	$apid = "1746";
    } 

    # Now report the values, if the user doesn't like them then abort
    print "_______________________________________________________\n";
    print "Data for Short Name $shortName (aka APID $apid)\n";
    print "Start Time            EndTime\n";
    for($i = 0; $i <= $#rVal; $i++) {
	if($rVal[$i]->{present} =~ /Yes/i) {
	    print $rVal[$i]->{startTime} . "   " . $rVal[$i]->{endTime} . "\n";
	} else {
	    print "Missing\n";
	}
    }
#    print "Continue? (y/n): ";
#    chomp($line = <STDIN>);
#    unless($line =~ /y/i) {
#	die "PCF Generation has been cancelled\n";
#    }

    @rVal;        # The return value
}

sub GetAttFiles {
    
    my($shortName) = $_[0];
    my($version)   = $_[1];
    my($startTime) = $_[2];
    my($endTime)   = $_[3];

    my($i, $sql, $sth, $row, @rVal, $thisTime);

    # Initialize the hash reference
    for($i = 0; $i < 12; $i++) {
	$rVal[$i] = {   present => "No",
			productName => "",
			path => "",
			startTime => "",
			endTime => "",
		    };
    }

    $sql = 'SELECT * FROM scfFromSips WHERE startTime >= "' . $startTime . '" ' .
  	   'AND endTime <= "' . $endTime . '" AND  shortName="' . $shortName . '" ';
    if($version ne "") {
	$sql .= 'AND version="' . $version . '" ';
    }
    $sql .= 'ORDER BY addTime ASC';
    $sth = $dbh->prepare($sql);
    $sth->execute or die "Unable to execute query: $dbh->errstr\n";

    while($row = $sth->fetchrow_hashref) {

	# Figure out what index the result fits into
	$row->{startTime} =~ /\s[0-9]{2}:/;
	$thisTime = $&;
	$thisTime =~ s/[\s:]//g;
	
	if($thisTime eq "22" || $thisTime eq "23") {
	    # Now see what day
	    $row->{startTime} =~ /\s/;
	    $thisTime = $`;
	    $startTime =~ /\s/;
	    $i = $thisTime eq $` ? 0 : 12;
	} elsif($thisTime eq "00" || $thisTime eq "01") {
	    $i = 0;
	    # Now see what day
	    $row->{endTime} =~ /\s/;
	    $thisTime = $`;
	    $endTime =~ /\s/;
	    $i = $thisTime eq $` ? 13 : 1;
	} elsif($thisTime eq "02" || $thisTime eq "03") {
	    $i = 2;
	} elsif($thisTime eq "04" || $thisTime eq "05") {
	    $i = 3;
	} elsif($thisTime eq "06" || $thisTime eq "07") {
	    $i = 4;
	} elsif($thisTime eq "08" || $thisTime eq "09") {
	    $i = 5;
	} elsif($thisTime eq "10" || $thisTime eq "11") {
	    $i = 6;
	} elsif($thisTime eq "12" || $thisTime eq "13") {
	    $i = 7;
	} elsif($thisTime eq "14" || $thisTime eq "15") {
	    $i = 8;
	} elsif($thisTime eq "16" || $thisTime eq "17") {
	    $i = 9;
	} elsif($thisTime eq "18" || $thisTime eq "19") {
	    $i = 10;
	} elsif($thisTime eq "20" || $thisTime eq "21") {
	    $i = 11;
	} else {
	    die "Uh Oh $thisTime $row->{startTime}\n";
	}
	$rVal[$i]{present} = "Yes";
	$rVal[$i]{productName} = $row->{productName};
	$rVal[$i]{path} = $row->{path};
	$rVal[$i]{startTime} = $row->{startTime};
	$rVal[$i]{endTime} = $row->{endTime};
    }

    # Now report the values, if the user doesn't like them then abort
    print "_______________________________________________________\n";
    print "Attitude ($shortName) Data\n";
    print "Start Time            EndTime\n";
    for($i = 0; $i <= $#rVal; $i++) {
	if($rVal[$i]->{present} =~ /Yes/i) {
	    print $rVal[$i]->{startTime} . "   " . $rVal[$i]->{endTime} . "\n";
	} else {
	    print "Missing\n";
	}
    }
#    print "Continue? (y/n): ";
#    chomp($line = <STDIN>);
#    unless($line =~ /y/i) {
#	die "PCF Generation has been cancelled\n";
#    }

    @rVal;        # The return value
}

sub GetEphFiles {
    
    my($shortName) = $_[0];
    my($version)   = $_[1];
    my($startTime) = $_[2];
    my($endTime)   = $_[3];

    my($i, $sql, $sth, $row, @rVal, $thisStart, $thisEnd);

    # Initialize the hash reference
    for($i = 0; $i < 2; $i++) {
	$rVal[$i] = {   present => "No",
			productName => "",
			path => "",
			startTime => "",
			endTime => "",
		    };
    }

    $sql = 'SELECT * FROM scfFromSips WHERE startTime >= "' . $startTime . '" ' .
	   'AND endTime <= "' . $endTime . '" AND  shortName="' . $shortName . '" ';
    if($version ne "") {
	$sql .= 'AND version="' . $version . '" ';
    }
    $sql .= 'ORDER BY addTime ASC';
    $sth = $dbh->prepare($sql);
    $sth->execute or die "Unable to execute query: $dbh->errstr\n";

    while($row = $sth->fetchrow_hashref) {

	# Figure out what index the result fits into
	$row->{startTime} =~ /\s/;
	$thisTime = $`;
	$row->{startTime} =~ /\s/;
	$thisStart = $`;
	$row->{endTime} =~ /\s/;
	$thisEnd = $`;

	$startTime =~ /\s/;

	if($` eq $thisStart) {
	    $i = 0;
	} else {
	    $endTime =~ /\s/;
	    unless($` eq $thisEnd) {
		die "Bad Time $row->{startTime}\n";
	    }
	    $i = 1;
	}

	$rVal[$i]{present} = "Yes";
	$rVal[$i]{productName} = $row->{productName};
	$rVal[$i]{path} = $row->{path};
	$rVal[$i]{startTime} = $row->{startTime};
	$rVal[$i]{endTime} = $row->{endTime};
    }

    # Now report the values, if the user doesn't like them then abort
    print "_______________________________________________________\n";
    print "Ephemeris ($shortName) Data\n";
    print "Start Time            EndTime\n";
    for($i = 0; $i <= $#rVal; $i++) {
	if($rVal[$i]->{present} =~ /Yes/i) {
	    print $rVal[$i]->{startTime} . "   " . $rVal[$i]->{endTime} . "\n";
	} else {
	    print "Missing\n";
	}
    }
#    print "Continue? (y/n): ";
#    chomp($line = <STDIN>);
#    unless($line =~ /y/i) {
#	die "PCF Generation has been cancelled\n";
#    }

    @rVal;        # The return value
}
