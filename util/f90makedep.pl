#!/usr/bin/perl
#
# Usage: f90makedep.pl
#
# It has become a multipurpose tool, exploiting some of perl's
# text-handling prowess. In its main usage,
# generate a Makefile-ready dependency listing from the f90 sources 
# in the current directory. If there are any command line arguments,
# and   $ARG_DEPENDS   is set (which see)
# they will be interpreted as additional directories in which to
# find dependencies. In particular, if when a source file USEs a
# module, the file containing that USEd module will be printed among
# the list of dependencies if it is found in the current directory or
# else found among the directories named in the command line. If it is not
# found among any of these, it will be excluded from the dependencies.
# In addition, if  $SELF_DEPENDS  is set (which see) the source file
# will itself be included among the list of dependencies
#
# Main Usage:
# 	(1) f90makedep.pl [arg1 arg2 ..]
# Output:
# (printed to stdout)
# source1.o: source1.f90 source2.o .. $(d1)s1.o $(d1)s2.o .. $(d2)t1.o ..
# source2.o: source2.f90 source3.o .. $(d1)s1.o $(d1)s2.o .. $(d2)t1.o ..
# source3.o: source3.f90 .. $(d1)s1.o $(d1)s2.o .. $(d2)t1.o ..
# ...
# where the contents of the current directory (cwd), arg1, arg2 are
#    cwd         arg1       arg2
# source1.f90   s1.f90     t1.f90
# source2.f90   s2.f90     t2.f90
# source3.f90   s3.f90     t3.f90
#
# The theory behind the $(d1), $(d2) stuff is that you will define them
# with lines like 
# d1 = ../up_one/elsewhere/ 
# and similarly for d2 in the same Makefile where you plan to paste the output
#
# To run this, it will be necessary to modify the first line of this script
# to point to the actual location of Perl on your system.
#
# 	Variant (1)-mod: f90makedep.pl -mod case [arg1 arg2 ..]
#   where case is one of {U[PPER], l[ower]}
# Output:
# (printed to stdout)
# source1.mod Source1.o: Source1.f90 source2.mod .. $(d1)s1.mod \
#             $(d1)s2.mod .. $(d2)t1.mod ..
#     $(UTILDIR)/newAifBdiff.sh -a source1.mod \
#        $(FC) -c $(FOPTS) $(INC_PATHS) $(S)/Source1.f90
#   where we have shown the lower case. Note that the stem name preceding
# the .mod suffix is the module name, which may differ from the file name
# Also note the build command which comes below each target in this variant
# The build command is built out of parts $BUILD_PART(1 2)
#
# -----------------------------------------------------------------------
# Alternate Usages:
# 	(2) f90makedep.pl -s pattern1 -o pattern2
# Output:
# (printed to stdout) (assume that pattern1 is "tex" and pattern2 is "dvi"
# source1.dvi: source1.tex source2.tex .. $(d1)s1.tex $(d1)s2.tex .. $(d2)t1.tex ..
# source2.dvi: source2.tex .. $(d1)s1.tex $(d1)s2.tex .. $(d2)t1.tex ..
# ...
# (currently this alternative only works for m4 files;
#  a clever programmer could fix it to work for tex, too)
# (also, this usage ignores any subsequent command line args: thus
#   no additional directories)
#
# 	(3) f90makedep.pl -p very_long_line -d delimiter -t term_char
# Output: splits the very_long_line into 72-column long lines
# (printed to stdout)
#  
# 	(4) f90makedep.pl -f file_name -d delimiter -t term_char -c comment_char
#      reads each line of file_name
# Output: splits any of its long lines into 72-column long lines
# (except for comments)
# (printed to stdout)
#  
# -----------------------------------------------------------------------
#      Options
# -s pattern        suffix to be matched by source files
# -o pattern        suffix to be matched by object files
# -p text           text to be split among several lines
# -d char           delimiter between words of text
# -t char           line termination at split
# -f file_name      file containing text to be split
# -db file_name     don't add build commands below target file_name
#                    if "(1)-mod" variant
# -c char           comment character
# -mod case         use variant "(1)-mod" of dependency lines
#                    with module names all in case (UPPER, lower)
#
# P.A. Wagner (May 15 2002)
# Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
# U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.
# however ..
# based on f90mkmf from H. Pumphrey (check that this is true in fact)
#
# Based on Michael Wester <wester@math.unm.edu> February 16, 1995
# Cotopaxi (Consulting), Albuquerque, New Mexico
# Bugs and limitations
#(1) Does too much--should be split up into dependency maker and
#     separate swiss army knife
#(2) Still too specialized--works well only with NAG-like compilers
#(3) Build commands used with (1)-mod variant are grotesque
#    shouldn't they be configurable at least?
#
# "$Id$"

#         Settings which affect how the script operates
#         * * * * * * * * * * * * * * * * * * * * * * * 
# if TRUE, debug by printing extra stuff
$DEBUG = 0;

# if TRUE, don't split long lines containing a quote mark "'"
$DONT_SPLIT_QUOTES = TRUE;

# if TRUE, include source file in list of dependencies
$SELF_DEPENDS = TRUE;

# if variant (1)-mod, what case to use in module names
$var_1_mod = 0;
$mod_case = "lower";

# if TRUE, include source files in directories named in arg list
# in list of dependencies (if arguments not empty)
# (only for usage (1) currently)
$ARG_DEPENDS = TRUE;

# if TRUE, for source files in directories named in arg list
# prepend $(d1), $(d2), .. stuff; else leave bare
$dir1_PREPENDS = TRUE;

#         * * * * * * * * * * * * * * * * * * * * * * * 
# Check command-line args for options -s etc.
$more_opts = TRUE;
$source_ext = '';
$obj_ext = '';
$long_line = '';
$long_file = '';
$delimiter = ' ';
$term_char = "\\";
$comment_char = "#";
while ($more_opts) {
   if ($ARGV[0] =~ /^-db/) {
      $dont_build{$ARGV[1]} = 1;
      shift;
      shift;
   } elsif ($ARGV[0] =~ /^-s/) {
      $source_ext = $ARGV[1];
      shift;
      shift;
   } elsif ($ARGV[0] =~ /^-o/) {
      $obj_ext = $ARGV[1];
      shift;
      shift;
   } elsif ($ARGV[0] =~ /^-p/) {
      $long_line = $ARGV[1];
      shift;
      shift;
   } elsif ($ARGV[0] =~ /^-d/) {
      $delimiter = $ARGV[1];
      shift;
      shift;
   } elsif ($ARGV[0] =~ /^-f/) {
      $long_file = $ARGV[1];
      shift;
      shift;
   } elsif ($ARGV[0] =~ /^-c/) {
      $comment_char = $ARGV[1];
      shift;
      shift;
   } elsif ($ARGV[0] =~ /^-t/) {
      $term_char = $ARGV[1];
      shift;
      shift;
   } elsif ($ARGV[0] =~ /^-mod/) {
      $var_1_mod = TRUE;
      $mod_case = $ARGV[1];
      shift;
      shift;
   } else {
      $more_opts = 0;
   }
}
# Fool the script into writing to stdout instead of an actual file
open(MAKEFILE, ">-");
# Dependency listings
#
if ($source_ext ne '') {
#  note: the final argument in the subroutine call below works for m4 files only
  &MakeDependsGeneral("*.$source_ext", "$obj_ext", '\!include\(([^)]+)');
} elsif ($long_file ne '') {
  &ReadAndSplit("$long_file", "$delimiter", "$term_char", "$comment_char");
} elsif ($long_line ne '') {
  &SplitAndPrint("$long_line", "$delimiter", "$term_char");
} else {
  &MakeDependsf90($ARGV[1]);
  # &MakeDepends("*.f *.F", '^\s*include\s+["\']([^"\']+)["\']');
  # &MakeDepends("*.c",     '^\s*#\s*include\s+["\']([^"\']+)["\']');
}
#
# &PrintWords(current output column, extra tab?, word list);
# --- print space-separated words nicely
#
sub PrintWords {
   local($columns) = 78 - shift(@_);
   local($extratab) = shift(@_);
   local($wordlength);
   #
   print MAKEFILE @_[0];
   $columns -= length(shift(@_));
   foreach $word (@_) {
      $wordlength = length($word);
      if ($wordlength + 1 < $columns) {
         print MAKEFILE " $word";
         $columns -= $wordlength + 1;
         }
      else {
         #
         # Continue onto a new line
         #
         if ($extratab) {
            print MAKEFILE " \\\n\t\t$word";
            $columns = 62 - $wordlength;
            }
         else {
            print MAKEFILE " \\\n\t$word";
            $columns = 70 - $wordlength;
            }
         }
      }
   }

#
# &PrintElements(current output column, delimiter, line term, element list);
# --- print delimiter-separated words nicely, with line term at end of each line
#
sub PrintElements {
   local($columns) = 78 - shift(@_);
   local($delimiter) = shift(@_);
   local($line_term) = shift(@_);
   local($wordlength);
   #
   print MAKEFILE @_[0];
   $columns -= length(shift(@_));
   foreach $word (@_) {
      $wordlength = length($word);
      if ($wordlength + 1 < $columns) {
         print MAKEFILE "$delimiter" . "$word";
         $columns -= $wordlength + 1;
         }
      else {
         #
         # Continue onto a new line
         #
         print MAKEFILE "$delimiter" . "$line_term" . "\n\t$word";   
         $columns = 70 - $wordlength;     
         }                                
      }
   }

#
# &LanguageCompiler(compiler, sources); --- determine the correct language
#    compiler
#
sub LanguageCompiler {
   local($compiler) = &toLower(shift(@_));
   local(@srcs) = @_;
   #
   if (length($compiler) > 0) {
      CASE: {
         grep(/^$compiler$/, ("fc", "f77")) &&
            do { $compiler = "FC"; last CASE; };
         grep(/^$compiler$/, ("cc", "c"))   &&
            do { $compiler = "CC"; last CASE; };
         $compiler = "F90";
         }
      }
   else {
      CASE: {
         grep(/\.f90$/, @srcs)   && do { $compiler = "F90"; last CASE; };
         grep(/\.(f|F)$/, @srcs) && do { $compiler = "FC";  last CASE; };
         grep(/\.c$/, @srcs)     && do { $compiler = "CC";  last CASE; };
         $compiler = "???";
         }
      }
   $compiler;
   }

#
# &toLower(string); --- convert string into lower case
#
sub toLower {
   local($string) = @_[0];
   $string =~ tr/A-Z/a-z/;
   $string;
   }

#
# &toUpper(string); --- convert string into upper case
#
sub toUpper {
   local($string) = @_[0];
   $string =~ tr/a-z/A-Z/;
   $string;
   }

#
# &uniq(sorted word list); --- remove adjacent duplicate words
#
sub uniq {
   local(@words);
   foreach $word (@_) {
      if ($word ne $words[$#words]) {
         push(@words, $word);
         }
      }
   @words;
   }

#
# &MakeDepends(language pattern, include file sed pattern); --- dependency
#    maker
#
sub MakeDepends {
   local(@incs);
   local($lang) = @_[0];
   local($pattern) = @_[1];
   #
   foreach $file (<${lang}>) {
      open(FILE, $file) || warn "Cannot open $file: $!\n";
      while (<FILE>) {
         /$pattern/i && push(@incs, $1);
         }
      if (defined @incs) {
         $file =~ s/\.[^.]+$/.o/;
         print MAKEFILE "$file: ";
         &PrintWords(length($file) + 2, 0, @incs);
         print MAKEFILE "\n";
         undef @incs;
         }
      }
   }

#
# &MakeDependsGeneral(source language pattern, object extension,
#   include file sed pattern); --- dependency maker
#
sub MakeDependsGeneral {
   local(@incs);
   local($lang) = @_[0];
   local($obj_ext) = @_[1];
   local($pattern) = @_[2];
   #
   foreach $file (<${lang}>) {
      open(FILE, $file) || warn "Cannot open $file: $!\n";
      while (<FILE>) {
         /$pattern/i && push(@incs, $1);
         }
      if (defined @incs || $SELF_DEPENDS) {
         if ($SELF_DEPENDS) {
		     @incs = reverse(@incs);
           push(@incs, $file);
		     @incs = reverse(@incs);
         }
         $file =~ s/\.[^.]+$/.$obj_ext/;
         print MAKEFILE "$file: ";
         &PrintWords(length($file) + 2, 0, @incs);
         print MAKEFILE "\n";
         undef @incs;
         }
      }
   }

#
# &ReadAndSplit(file to print, delimiter, term char, comment char
#
sub ReadAndSplit {
   local(@incs);
   local($line);
   local($file_name) = @_[0];
   local($delimiter) = @_[1];
   local($term_char) = @_[2];
   local($comment_char) = @_[3];
   #
   open(LONGFILE, "$file_name");
   while (<LONGFILE>) {
      $line = $_;
      if ($line =~ /^\s*$comment_char/) {
         print MAKEFILE $line;
      } elsif ($DONT_SPLIT_QUOTES && $line =~ /'/) {
         print MAKEFILE $line;
      } else {
        chomp($line);
        @incs = split("$delimiter","$line");
        &PrintElements(0, "$delimiter", "$term_char", @incs); 
        print MAKEFILE "\n";                                 
      }
   }
   close LONGFILE;
}

#
# &SplitAndPrint(line to print, delimiter, term char at end of each line
#
sub SplitAndPrint {
   local(@incs);
   local($line) = @_[0];
   local($delimiter) = @_[1];
   local($term_char) = @_[2];
   #
   @incs = split("$delimiter","$line");
   &PrintElements(0, "$delimiter", "$term_char", @incs);      
   print MAKEFILE "\n";                                      
   }

#
# &MakeDependsf90(f90 compiler); --- FORTRAN 90 dependency maker
#
sub MakeDependsf90 {
   local($compiler) = &toLower(@_[0]);
   local(@dependencies);
#  To store object file name keyed by module name
   local(%filename);
#  To module file name keyed by object name
   local(%modulename_by_file);
#  To module file name keyed by module name
   local(%modulename_by_module);
   local(@incs);
   local(@modules);
   local($objfile);
   local($arg_number);
   local($prependant);
   local($key);
   local($value);
   local($BUILD_PART1) = '$(UTILDIR)/newAifBdiff.sh -a ';
   local($BUILD_PART2) = '$(FC) -c $(FOPTS) $(INC_PATHS) $(S)/';
   #
   # Associate each module with the name of the file that contains it
   #
   foreach $file (<*.f90>) {
      open(FILE, $file) || warn "Cannot open $file: $!\n";
      while (<FILE>) {
         # This bit of filtering is to exclude statements that begin 
         # module procedure ...
         if ( /^\s*module\s+[pP][rR][oO][cC][eE][dD][uU][rR][eE]/ ) {    
            value = $key;                                                
         } else {                                                        
           /^\s*module\s+([^\s!]+)/i &&
            ($filename{&toLower($1)} = $file) =~ s/\.f90$/.o/;
         }
       }
   }
   if ($DEBUG) {
      print "object file name keyed by module name (in home directory)\n";
      while (($key,$value) = each %filename) {     
                   print "$key=$value\n";          
           }                                       
      print "\n";
   }

   #
   # Do the same in every directory named in arg list
   #
   if ($ARG_DEPENDS) {
   #
   # chop($cwd = \`pwd\`);
   # The above line, although found among man pages, won't pass muster
   #
   # The following is a poor way to fill $cwd with the current working directory
   # Ask a more perl-savvy programmer how to do this w/o creating temp file
     system "pwd > temp";
     open (TEMPFILE, temp);
     while (<TEMPFILE>) {
       chop;
       $cwd = $_;
       }
     unlink temp;
 #   Now loop over list of directories which are command-line arguments
    $arg_number = 0;
    while ($#ARGV >= $[) {
   	die "Can't cd to $ARGV[$[]:$!\n"
        	unless chdir $ARGV[$[] ;
        shift @ARGV;
        $arg_number++;
        if ($dir1_PREPENDS) {
           $prependant = '$(dir' . "$arg_number" . ')';
        } else {        
            $prependant = '';
       }
       foreach $file (<*.f90>) {
          open(FILE, $file) || warn "Cannot open $file: $!\n";
          ($objfile = "$prependant" . "$file") =~ s/\.f90$/.o/;
          $key = '';
          while (<FILE>) {
             # This bit of filtering is to exclude statements that begin 
             # module procedure ...
             if ( /^\s*module\s+[pP][rR][oO][cC][eE][dD][uU][rR][eE]/ ) {
                value = $key;
             } else {
               /^\s*module\s+([^\s!]+)/i &&
                ($key = &toLower($1));
    #            ($filename{&toLower($1)} = "$prependant" . "$file") =~ s/\.f90$/.o/;
             }
           }
          # Now add to our hash, but only if this is a new key  
          # W a r n i n g   W a r n i n g
          # This assumes the following behavior
          # if %hash defined for $key, then $hash{$key} != '' 
          # otherwise $hash{$key} = '' 
          # if instead $hash{$key} = '' for a defined $key
          # or else $hash{$key} dumps core for an undefined $key
          # then we will need to revisit and revise this
          $value = "$filename{$key}.suffix";
          # print "file ($file), key ($key), value ($value) \n";
          if ($value =~ /^\.suf/) {
            $filename{$key} = $objfile;
          }
       }
       chdir $cwd;
      }
    }

   # Find module name given object file name
   # We can't simply use the following:
   # %modulename_by_file = reverse %filename;
   # for two reasons
   # (1) filenames from outside directories come with "$(dirn)" prepended
   # (2) Depending on $mod_case we may want modulename lower or UPPER
   # So the 1st thing we do is acquire key, value pair from filename
   while (($key, $value) = each %filename) {
         # Next, value  => $(prependant)file
         if ( $value =~ /\(dir/) {
           $_ = $value;
            /\((.*?)\)(.*)/ && ($prependant=$1) && ($file=$2);
           if ($mod_case =~ /^U/) {
               $modulename_by_file{$file} = &toUpper($key);
               $modulename_by_module{$key} = &toUpper($key);
           } else {
               $modulename_by_file{$file} = &toLower($key);
               $modulename_by_module{$key} = &toLower($key);
           }
           # Finally, modulename_by_file  = $(prependant)module
           $prependant = '$(' . $prependant . ')';
           $modulename_by_file{$file} = $prependant .  $modulename_by_file{$file};
           $modulename_by_module{$key} = $prependant .  $modulename_by_module{$key};
         } else {
           if ($mod_case =~ /^U/) {
               $modulename_by_file{$value} = &toUpper($key);
               $modulename_by_module{$key} = &toUpper($key);
           } else {
               $modulename_by_file{$value} = &toLower($key);
               $modulename_by_module{$key} = &toLower($key);
           }
         }
       }  

   if ($DEBUG) {
      print "object file name keyed by module name \n";
      while (($key,$value) = each %filename) {     
                   print "$key=$value\n";          
           }                                       
      print "\n";
      print "module file name keyed by object name \n";
      while (($key,$value) = each %modulename_by_file) {   
                   print "$key=$value\n";          
           }                                       
      print "\n";
      print "module file name keyed by module name \n";
      while (($key,$value) = each %modulename_by_module) {  
                   print "$key=$value\n";          
           }                                       
      print "\n";
   }

   #
   # Print the dependencies of each file that has one or more include's or
   # references one or more modules
   #
   foreach $file (<*.f90>) {
      open(FILE, $file);
      while (<FILE>) {
         /^\s*include\s+["\']([^"\']+)["\']/i && push(@incs, $1);
         /^\s*use\s+([^\s,!]+)/i && push(@modules, &toLower($1));
         }
      if (defined @incs || defined @modules) {
         ($objfile = $file) =~ s/\.f90$/.o/;
         undef @dependencies;
         if ( $var_1_mod ) {
#           print MAKEFILE "$modulename_by_file{$objfile}.mod $objfile: ";
           print MAKEFILE "$objfile: ";
           if ("$modulename_by_file{$objfile}.mod" !~ /^\./) {
              print MAKEFILE "$modulename_by_file{$objfile}.mod ";
           }
           foreach $module (@modules) {
              $value = "$modulename_by_module{$module}.mod";
              if ($value !~ /^\./) {
                push(@dependencies, $value);
              }
           }
         } else {
           print MAKEFILE "$objfile: ";
           foreach $module (@modules) {
              push(@dependencies, $filename{$module});
           }
         }
         @dependencies = &uniq(sort(@dependencies));
         if ($SELF_DEPENDS) {
#         if ($SELF_DEPENDS && ($var_1_mod != TRUE)) {
#		     @dependencies = reverse(push(reverse(@dependencies), $file));
		     @dependencies = reverse(@dependencies);
                    push(@dependencies, $file);
		     @dependencies = reverse(@dependencies);
         }
         &PrintWords(length($objfile) + 2, 0,
                     @dependencies, &uniq(sort(@incs)));
         print MAKEFILE "\n";
         if ( $var_1_mod && ("$modulename_by_file{$objfile}.mod" !~ /^\./)) {
           print MAKEFILE "$modulename_by_file{$objfile}.mod: \n";
#           print MAKEFILE "$modulename_by_file{$objfile}.mod: $file\n";
           if ($dont_build{$file} != 1) {
             print MAKEFILE "\t";
             print MAKEFILE $BUILD_PART1;
             print MAKEFILE "$modulename_by_file{$objfile}.mod ";
             print MAKEFILE $BUILD_PART2;
             print MAKEFILE $file;
             print MAKEFILE "\n";
             print MAKEFILE "$objfile: \n";
             print MAKEFILE "\t";
             print MAKEFILE $BUILD_PART1;
             print MAKEFILE "$modulename_by_file{$objfile}.mod ";
             print MAKEFILE $BUILD_PART2;
             print MAKEFILE $file;
             print MAKEFILE "\n";
           }
         }
         undef @incs;
         undef @modules;
         }
      }
   }
# $Log$
# Revision 1.3  2002/04/09 19:39:36  pwagner
# Alternate usages 2-4 (depending on options) added
#
# Revision 1.2  2000/11/02 23:22:38  pwagner
# Dependencies may cross directories
#
# Revision 1.1  2000/10/27 23:21:26  pwagner
# First commit
#
