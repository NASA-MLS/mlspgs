! Copyright (c) 2000, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!===============================================================================
module HDFEOS               ! F90 interface to HDF-EOS.
!===============================================================================

  implicit none
  public


  !------------------- RCS Ident Info -----------------------
  character(len=130), private :: Id = &
    & "$Id$"
  character(len=*), parameter, private :: ModuleName= &
    & "$RCSfile$"
  !----------------------------------------------------------

  ! Now define f90 interfaces for some HDF-EOS.

  interface

  ! grid interfaces
   integer function GDATTACH ( GDFID, GRIDNAME )
      integer, intent(in) :: GDFID
      character (len=*), intent(in) :: GRIDNAME
    end function GDATTACH

    integer function GDCLOSE ( FILE_ID )
      integer, intent(in) :: FILE_ID
    end function GDCLOSE

    integer function GDDETACH ( GRIDID )
      integer, intent(in) :: GRIDID
    end function GDDETACH

    integer function GDDIMINFO (GRIDID,DIMNAME)
       integer,intent(in)::GRIDID
       character(len=*),intent(IN)::DIMNAME
    end function GDDIMINFO

    integer function GDFLDINFO (GRIDID,FIELDNAME,RANK,DIMS,NUMBERTYPE,DIMLIST)
       integer,intent(in)::GRIDID
       character(len=*),intent(IN)::FIELDNAME
       integer,intent(out),dimension(*):: DIMS
       integer,intent(out)             :: RANK, NUMBERTYPE
       character(len=*),intent(OUT)::DIMLIST
    end function GDFLDINFO

    integer function GDINQATTRS (GRIDID,ATTRLIST,STRINGBUFFERSIZE)
       integer,intent(in)::GRIDID
       character(len=*),intent(out)::ATTRLIST
       integer,intent(out)::STRINGBUFFERSIZE
    end function GDINQATTRS
    
    integer function GDINQDIMS (GRIDID,DIMNAME,DIMS)
       integer,intent(in)::GRIDID
       character(len=*),intent(out)::DIMNAME
       integer,intent(out),dimension(*)::DIMS
    end function GDINQDIMS
    
    integer function GDINQFLDS (GRIDID, FIELDLIST, RANK, NUMBERTYPE)
       integer,intent(in)::GRIDID
       character(len=*),intent(out)::FIELDLIST
       integer,intent(out),dimension(*)::RANK, NUMBERTYPE
    end function GDINQFLDS
    
    integer function GDINQGRID (FILENAME, GRIDLIST, STRINGBUFFERSIZE)
       character (len=*), intent(in) :: FILENAME
       character(len=*),intent(out)::GRIDLIST
       integer,intent(out)::STRINGBUFFERSIZE
    end function GDINQGRID
    
    integer function GDNENTRIES (GRIDID, ENTRYCODE, STRINGBUFFERSIZE)
       integer,intent(in)::GRIDID
       integer,intent(in)::ENTRYCODE
       integer,intent(out)::STRINGBUFFERSIZE
    end function GDNENTRIES
    
    integer function GDOPEN ( FILENAME, ACCESS_MODE )
      character (len=*), intent(in) :: FILENAME
      integer, intent(in) :: ACCESS_MODE
    end function GDOPEN

  ! swath interfaces
    integer function SWATTACH ( SWFID, SWATHNAME )
      integer, intent(in) :: SWFID
      character (len=*), intent(in) :: SWATHNAME
    end function SWATTACH

    integer function SWCLOSE ( FILE_ID )
      integer, intent(in) :: FILE_ID
    end function SWCLOSE

    integer function SWCREATE ( SWFID, SWATHNAME )
      integer, intent(in) :: SWFID
      character (len=*), intent(in) :: SWATHNAME
    end function SWCREATE

    integer function SWDEFDFLD ( SWATHID, FIELDNAME, DIMLIST, NUMBERTYPE, &
      & MERGE )
      integer, intent(in) :: SWATHID
      character (len=*), intent(in) :: FIELDNAME
      character (len=*), intent(in) :: DIMLIST
      integer, intent(in) :: NUMBERTYPE
      integer, intent(in) :: MERGE
    end function SWDEFDFLD

    integer function SWDEFDIM ( SWATHID, FIELDNAME, DIM )
      integer, intent(in)  :: SWATHID
      character (len=*), intent(in) :: FIELDNAME
      integer, intent(in) :: DIM
    end function SWDEFDIM 

    integer function SWDEFGFLD ( SWATHID, FIELDNAME, DIMLIST, NUMBERTYPE, &
      & MERGE )
      integer, intent(in) :: SWATHID
      character (len=*), intent(in) :: FIELDNAME
      character (len=*), intent(in) :: DIMLIST
      integer, intent(in) :: NUMBERTYPE
      integer, intent(in) :: MERGE
    end function SWDEFGFLD

    integer function SWDETACH ( SWID )
      integer, intent(in) :: SWID
    end function SWDETACH

    integer function SWOPEN ( FILENAME, ACCESS_MODE )
      character (len=*), intent(in) :: FILENAME
      integer, intent(in) :: ACCESS_MODE
    end function SWOPEN

    integer function SWINQSWATH (FILENAME,SWATHLIST,STRBUFSIZE)
      character (len=*), intent(in) :: FILENAME
      character (len=*), intent(out) :: SWATHLIST
      integer, intent(out):: STRBUFSIZE
    end function SWINQSWATH

    integer function SWINQDIMS (SWATHID,DIMNAME,DIMS)
       integer,intent(in)::SWATHID
       character(len=*),intent(out)::DIMNAME
       integer,intent(out),dimension(*)::DIMS
    end function SWINQDIMS
    
    integer function SWINQDFLDS (SWATHID, FIELDLIST, RANK, NUMBERTYPE)
       integer,intent(in)::SWATHID
       character(len=*),intent(out)::FIELDLIST
       integer,intent(out),dimension(*)::RANK, NUMBERTYPE
    end function SWINQDFLDS
    
    integer function SWDIMINFO (SWATHID,DIMNAME)
       integer,intent(in)::SWATHID
       character(len=*),intent(IN)::DIMNAME
    end function SWDIMINFO


  end interface

!PAGE_BREAK
!---------------------------------
! Access types for gdnentries
!---------------------------------

   INTEGER, PARAMETER :: HDFE_NENTDIM            =     0
   INTEGER, PARAMETER :: HDFE_NENTDFLD           =     4

!====================
end module HDFEOS
!====================

!
! $Log$
! Revision 2.10  2001/03/15 00:36:03  pwagner
! Added gdfldinfo
!
! Revision 2.9  2001/03/14 00:33:42  pwagner
! Added gdinqattrs
!
! Revision 2.8  2001/03/03 00:45:25  pwagner
! New grid interfaces--except for gdrdfld
!
! Revision 2.7  2001/02/24 01:01:26  pwagner
! Started adding gd routine interfaces
!
! Revision 2.6  2001/02/23 17:27:20  pwagner
! Added Access types for gdnentries
!
! Revision 2.5  2001/02/02 21:39:56  pwagner
! Added swinqdflds
!
! Revision 2.4  2000/09/17 20:49:33  pumphrey
! Removed references to SWRDFLD as these are handled by swapi.f90
!
! Revision 2.3  2000/09/14 11:16:00  pumphrey
! HCP added SWRDFLD_REAL etc as ways to call SWRDFLD. Is there a better way?
!
! Revision 2.2  2000/09/13 16:56:22  pumphrey
! HCP added  swdiminfo and swinqdims -- more to do, though.
!
! Revision 2.1  2000/09/13 16:28:49  pumphrey
! HCP added SWINQSWATH
!
! Revision 2.0  2000/09/05 17:41:06  dcuddy
! Change revision to 2.0
!
! Revision 1.5  2000/09/02 01:58:30  vsnyder
! Cosmetic changes
!
