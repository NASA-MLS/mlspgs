# Usage: f90GhostFiles.pl 'b_dir' 's_dir1' 's_dir2' ..
#
# Generate a list of files with .o and .mod suffixes in
# binary-directory b_dir for which corresponding files with .f90
# suffixes are not found in any of the source-directories s_dirn. 
# The most likely reason being that the original source file
# has been moved or deleted, yet the persistent .o or .mod
# file still "haunts" the binary-directory.
#
# Output:
# (printed to stdout)
# Ghost1.o Ghost2.o .. Ghostn.o ghost1.mod ghost2.mod .. ghostn.mod
#
# P.A. Wagner (April 25 2001)
# Copyright 2005, by the California Institute of Technology. ALL
# RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
# commercial use must be negotiated with the Office of Technology Transfer
# at the California Institute of Technology.

# This software may be subject to U.S. export control laws. By accepting this
# software, the user agrees to comply with all applicable U.S. export laws and
# regulations. User has the responsibility to obtain export licenses, or other
# export authority as may be required before exporting such information to
# foreign countries or providing access to foreign persons.
#
# "$Id$"

#         Settings which affect how the script operates
#         * * * * * * * * * * * * * * * * * * * * * * * 
# if TRUE, include source files with .c suffixes
$CC_SOURCES = TRUE;

# if TRUE, include source files with .f suffixes
$F77_SOURCES = TRUE;

# if TRUE, format output suitable for inclusion in a Makefile
$MAKEFILE_FORMAT = FALSE;

$debug = 0;				# Print more if 1

#         * * * * * * * * * * * * * * * * * * * * * * * 
# Fool the script into writing to stdout instead of an actual file
open(MAKEFILE, ">-");

print STDERR "Binary directory: $ARGV[0] \n" unless !($debug);
print STDERR "Source directory: $ARGV[1] \n" unless !($debug);

# GhostFiles
#
&GhostFiles($ARGV[1], $ARGV[2]);

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
# GhostFiles; --- lists dead .o and .mod files
#
sub GhostFiles {
   local($source_dir) = $ARGV[$[ + 1];
   local($binary_dir) = $ARGV[$[];
   local($current_source_dir_num);
   local($extra_source_dir);
   local(@f90_sources);
   local($num_source_dirs);
   local(@objects);
   local($recasedF90Source);
   local(@source_dirs);
   local(@sources);
   local($modulename);
   local(@modules);
   local($objfile);
   local($cwd);
   local($arg_number);
   local($num_sources);
   local($is_ghost);
   local($got_module);

print STDERR "Num args: $#ARGV \n" unless !($debug);
print STDERR "Source directory: $source_dir \n" unless !($debug);
print STDERR "Binary directory: $binary_dir \n" unless !($debug);

   if ($#ARGV > 1) {
      $extra_source_dir = $ARGV[$[ + 2];
   } else {
      $extra_source_dir = ' ';
   }
print STDERR "Extra source directory: $extra_source_dir \n" unless !($debug);
   $num_source_dirs = $#ARGV;
#   $source_dirs[0] = $source_dir;
#   $source_dirs[1] = $extra_source_dir;
   $current_source_dir_num = 0;
   
   while ($current_source_dir_num < $num_source_dirs) {
#  source directories
      $source_dirs[$current_source_dir_num] = $ARGV[$current_source_dir_num + 1];;
      $current_source_dir_num++;
   }     # end while loop over source_dirs
   
   #
   chomp($cwd = `pwd`);
   $num_sources = 0;

   $num_source_dirs = $#ARGV;
   $current_source_dir_num = 0;
   
   while ($current_source_dir_num < $num_source_dirs) {
#  source directories
   	die "Can't cd to source_dir:$!\n"
        	unless chdir $source_dirs[$current_source_dir_num] ;
   # Associate each module with the name of the file that contains it
   #
   @f90_sources = <*.[fF]90>;
   foreach $file (@f90_sources) {
      $sources[$num_sources] = $file;
      ($recasedF90Source = $file) =~ s/\.F90$/.f90/;
      ($objfile = $recasedF90Source) =~ s/\.f90$/.o/;
      $objects[$num_sources] = $objfile;
    $got_module = "false";
     open(FILE, $file) || warn "Cannot open $file: $!\n";
print STDERR "file: $file \n" unless !($debug);
      while (<FILE>) {
         if ($got_module ne "true") {
#print STDERR "line: $_ \n" unless !($debug);
          if (/^\s*module\s+(\w+)/i) {
            /^\s*module\s+(\w+)/i &&
            ($modulename = &toLower($1) . ".f90") =~ s/\.f90$/.mod/;
             $got_module = "true";
#print STDERR "line where got module was: $_ \n" unless !($debug);
#print STDERR "module was: $modulename \n" unless !($debug);
             }
            }
         }
      $modules[$num_sources] = $modulename;
      $num_sources++;
#      exit;
      }

print STDERR "num_sources: $num_sources \n \n" unless !($debug);
print STDERR "sources: @sources \n \n" unless !($debug);
print STDERR "objects: @objects \n \n" unless !($debug);
print STDERR "modules: @modules \n \n" unless !($debug);

   # In case any .c or .f files
   #
   if ($CC_SOURCES) {
    foreach $file (<*.c>) {
      $sources[$num_sources] = $file;
      ($objfile = $file) =~ s/\.c$/.o/;
      $objects[$num_sources] = $objfile;
# Of course, they will produce no .mod files when compiled
      $modules[$num_sources] = ' ';
      $num_sources++;
      }
     }

   if ($F77_SOURCES) {
    foreach $file (<*.f>) {
      $sources[$num_sources] = $file;
      ($objfile = $file) =~ s/\.f$/.o/;
      $objects[$num_sources] = $objfile;
# Of course, they will produce no .mod files when compiled
      $modules[$num_sources] = ' ';
      $num_sources++;
      }
     }

       chdir $cwd;

   $current_source_dir_num++;
   }     # end while loop over source_dirs

   	die "Can't cd to binary_dir:$!\n"
        	unless chdir $binary_dir ;

   # Check for ghosts among the .o suffixes
   #
   undef @ghosts;
   foreach $file (<*.o>) {
    $arg_number = 0;
    $is_ghost = "true";
      while ($is_ghost ne "false" && $arg_number < $num_sources) {
         if ($file ne $objects[$arg_number]) {
         } else {
          $is_ghost = "false";
            }
 #        print STDERR "file, object[i]: $file, $objects[$arg_number], $is_ghost  \n" unless !($debug);
        $arg_number++;
         }
    if ($is_ghost eq "true") {
            push(@ghosts, $file);
         }
      }

   # Check for ghosts among the .mod suffixes
   #
   foreach $file (<*.mod>) {
    $arg_number = 0;
    $is_ghost = "true";
      while ($is_ghost ne "false" && $arg_number < $num_sources) {
#       We have decided to store only lower-case module names in modules, so ..
#         if ($file ne $modules[$arg_number]) {
         if (&toLower($file) ne $modules[$arg_number]) {
         } else {
          $is_ghost = "false";
            }
 #        print STDERR "file, modules[i]: $file, $modules[$arg_number], $is_ghost  \n" unless !($debug);
        $arg_number++;
         }
    if ($is_ghost eq "true") {
            push(@ghosts, $file);
         }
      }

   if ($MAKEFILE_FORMAT) {
    &PrintWords(1, 0,
                     &uniq(sort(@ghosts)));
     print MAKEFILE "\n";
   } else {
     print MAKEFILE &uniq(sort(@ghosts));
   }
  }
# $Log$
# Revision 1.7  2007/03/06 21:30:10  pwagner
# Deleted outmoded lines concerning 'actua;l location of perl'
#
# Revision 1.6  2005/06/23 22:22:46  pwagner
# Reworded Copyright statement
#
# Revision 1.5  2004/11/03 19:09:33  pwagner
# perl scripts now get launched via perl rather than as stand-alone executables
#
# Revision 1.4  2004/03/18 17:58:09  pwagner
# Ended foolish use of temp file to learn cwd
#
# Revision 1.3  2003/09/24 19:18:42  pwagner
# Restored to proper functioning despite new perl
#
# Revision 1.2  2001/11/26 23:46:45  pwagner
# Now works better with Absoft--test lower-case against module names
#
# Revision 1.1  2001/05/01 17:03:31  pwagner
# First commit
#
