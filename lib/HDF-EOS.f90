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
    
    integer function SWDIMINFO (SWATHID,DIMNAME)
       integer,intent(in)::SWATHID
       character(len=*),intent(IN)::DIMNAME
    end function SWDIMINFO


  end interface

!====================
end module HDFEOS
!====================

!
! $Log$
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
