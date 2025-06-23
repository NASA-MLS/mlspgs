/* Copyright 2005, by the California Institute of Technology. ALL
   RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
   commercial use must be negotiated with the Office of Technology Transfer
   at the California Institute of Technology.

   This software may be subject to U.S. export control laws. By accepting this
   software, the user agrees to comply with all applicable U.S. export laws and
   regulations. User has the responsibility to obtain export licenses, or other
   export authority as may be required before exporting such information to
   foreign countries or providing access to foreign persons.
*/

/* We can't usefully access stat() or readdir() directly from Fortran using
   interoperability features because their results are structs for which the
   order and sizes of members, and the presence or absence of padding, are not
   standardized.
*/

#include <dirent.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <stdio.h>

char Id[] = "$Id$";

/* Get the file type using stat() */

int FileType ( const char *PathName, int *Error )
/* Returns S_IFSOCK   0140000   socket
           S_IFLNK    0120000   symbolic link
           S_IFREG    0100000   regular file
           S_IFBLK    0060000   block device
           S_IFDIR    0040000   directory
           S_IFCHR    0020000   character device
           S_IFIFO    0010000   FIFO
shifted right 12 bits (i.e., eliminating the four low-order octal zeros) if
PathName is one of these, else -1 (usually because PathName doesn't exist).
Because we use stat() instead of lstat(), S_IFLNK is never returned.
If the result value is -1, Error is set equal to errno, else it's zero.
*/
{
  struct stat sb;
  if ( stat ( PathName, &sb ) < 0 )
  { *Error = errno;
    return -1;
  }
  *Error = 0;
  return ( ( sb.st_mode & S_IFMT ) >> 12 );
}

/* Read the next item from a directory opened by opendir().  If the
   returned file type is a link (DT_LNK == S_IFLNK >> 12), you can use
   FileType to inquire the type of file the link refers to. */

int ReadDir ( DIR *DirHandle, char **FileName, int *FileType )
{
  struct dirent *dp;
  dp = readdir ( DirHandle );
  if ( dp == 0 )
  { FileName = 0;
    return -1;
  };
  if ( dp->d_reclen <= 0 )
  { FileName = 0;
    return -1;
  };
  if ( dp->d_name == 0 )
  { FileName = 0;
    return -1;
  };
  *FileName = dp->d_name;
  *FileType = dp->d_type;
  return dp->d_reclen;
}

/* $Log$
 * Revision 2.1  2016/10/21 22:53:17  vsnyder
 * Initial commit
 *

*/
