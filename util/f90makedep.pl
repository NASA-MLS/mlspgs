#!/usr/bin/perl
#
# Usage: f90makedep.pl
#
# Generate a Makefile-ready dependency listing from the f90 sources 
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
# Usage:
# 	f90makedep.pl [arg1 arg2 ..]
# Output:
# (printed to stdout)
# source1.o: source1.f90 source2.f90 .. $(d1)s1.f90 $(d1)s2.f90 .. $(d2)t1.f90 ..
# source2.o: source2.f90 source3.f90 .. $(d1)s1.f90 $(d1)s2.f90 .. $(d2)t1.f90 ..
# source3.o: source3.f90 .. $(d1)s1.f90 $(d1)s2.f90 .. $(d2)t1.f90 ..
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
# P.A. Wagner (October 31 2000)
# (see mlsconfigure.sh for copyright statement)
#
# Based on f90mkmf from H. Pumphrey (check that this is true in fact)
#
# Based on Michael Wester <wester@math.unm.edu> February 16, 1995
# Cotopaxi (Consulting), Albuquerque, New Mexico
#
# "$Id$"

#         Settings which affect how the script operates
#         * * * * * * * * * * * * * * * * * * * * * * * 
# if TRUE, include source file in list of dependencies
$SELF_DEPENDS = TRUE;

# if TRUE, include source files in directories named in arg list
# in list of dependencies (if arguments not empty)
$ARG_DEPENDS = TRUE;

# if TRUE, for source files in directories named in arg list
# prepend $(d1), $(d2), .. stuff; else leave bare
$dir1_PREPENDS = TRUE;

#         * * * * * * * * * * * * * * * * * * * * * * * 
# Fool the script into writing to stdout instead of an actual file
open(MAKEFILE, ">-");
# Dependency listings
#
&MakeDependsf90($ARGV[1]);
# &MakeDepends("*.f *.F", '^\s*include\s+["\']([^"\']+)["\']');
# &MakeDepends("*.c",     '^\s*#\s*include\s+["\']([^"\']+)["\']');

#
# &PrintWords(current output column, extra tab?, word list); --- print words
#    nicely
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
# &MakeDependsf90(f90 compiler); --- FORTRAN 90 dependency maker
#
sub MakeDependsf90 {
   local($compiler) = &toLower(@_[0]);
   local(@dependencies);
   local(%filename);
   local(@incs);
   local(@modules);
   local($objfile);
   local($arg_number);
   local($prependant);
   #
   # Associate each module with the name of the file that contains it
   #
   foreach $file (<*.f90>) {
      open(FILE, $file) || warn "Cannot open $file: $!\n";
      while (<FILE>) {
         /^\s*module\s+([^\s!]+)/i &&
            ($filename{&toLower($1)} = $file) =~ s/\.f90$/.o/;
         }
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
           while (<FILE>) {
             /^\s*module\s+([^\s!]+)/i &&
                ($filename{&toLower($1)} = "$prependant" . "$file") =~ s/\.f90$/.o/;
             }
          }
       chdir $cwd;
      }
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
         print MAKEFILE "$objfile: ";
         undef @dependencies;
         foreach $module (@modules) {
            push(@dependencies, $filename{$module});
            }
         @dependencies = &uniq(sort(@dependencies));
         if ($SELF_DEPENDS) {
#		@dependencies = reverse(push(reverse(@dependencies), $file));
		@dependencies = reverse(@dependencies);
               push(@dependencies, $file);
		@dependencies = reverse(@dependencies);
         }
         &PrintWords(length($objfile) + 2, 0,
                     @dependencies, &uniq(sort(@incs)));
         print MAKEFILE "\n";
         undef @incs;
         undef @modules;
         }
      }
   }
# $Log$
# Revision 1.1  2000/10/27 23:21:26  pwagner
# First commit
#
