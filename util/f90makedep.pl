#!/usr/bin/perl
#
# Usage: f90makedep.pl
#
# Generate a Makefile-ready dependency listing from the f90 sources 
# in the current directory. 
# To run this, it will be necessary to modify the first line of this script
# to point to the actual location of Perl on your system.
#
# Based on f90mkmf from H. Pumphrey
#
# Based on Michael Wester <wester@math.unm.edu> February 16, 1995
# Cotopaxi (Consulting), Albuquerque, New Mexico
#
# "$Id$"

# if TRUE, include source file in list of dependencies
$SELF_DEPENDS = TRUE;

# Fool the script into writing to stdout instead of an actual file
open(MAKEFILE, ">-");
# print "Why do I quit prematurely?";
# exit;
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
# Revision 1.2  2000/03/11 00:41:59  vsnyder
# Inserted CVS id and log
#
