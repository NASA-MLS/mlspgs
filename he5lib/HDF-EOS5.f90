! Copyright (c) 2000, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!===========================================================================
module HDFEOS5               ! F90 interface to HDF-EOS.
!===========================================================================

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
    integer function HE5_SWATTACH ( SWFID, SWATHNAME )
      integer, intent(in) :: SWFID
      character (len=*), intent(in) :: SWATHNAME
    end function HE5_SWATTACH

    integer function HE5_SWCLOSE ( FILE_ID )
      integer, intent(in) :: FILE_ID
    end function HE5_SWCLOSE
    integer function HE5_SWCREATE ( SWFID, SWATHNAME )
      integer, intent(in) :: SWFID
      character (len=*), intent(in) :: SWATHNAME
    end function HE5_SWCREATE

    integer function HE5_SWDEFDFLD ( SWATHID, FIELDNAME, DIMLIST, MAXDIMLIST,&
         NUMBERTYPE,MERGE )
      integer, intent(in) :: SWATHID
      character (len=*), intent(in) :: FIELDNAME
      character (len=*), intent(in) :: DIMLIST
      character (len=*), intent(in) :: MAXDIMLIST
      integer, intent(in) :: NUMBERTYPE
      integer, intent(in) :: MERGE
    end function HE5_SWDEFDFLD

    integer function HE5_SWDEFDIM ( SWATHID, FIELDNAME, DIM )
      integer, intent(in)  :: SWATHID
      character (len=*), intent(in) :: FIELDNAME
      integer, intent(in) :: DIM
    end function HE5_SWDEFDIM 

    integer function HE5_SWDEFGFLD ( SWATHID, FIELDNAME, DIMLIST, MAXDIMLIST,&
         NUMBERTYPE,  MERGE )
      integer, intent(in) :: SWATHID
      character (len=*), intent(in) :: FIELDNAME
      character (len=*), intent(in) :: DIMLIST
      character (len=*), intent(in) :: MAXDIMLIST
      integer, intent(in) :: NUMBERTYPE
      integer, intent(in) :: MERGE
    end function HE5_SWDEFGFLD

    integer function HE5_SWDETACH ( SWID )
      integer, intent(in) :: SWID
    end function HE5_SWDETACH

    integer function HE5_SWOPEN ( FILENAME, ACCESS_MODE )
      character (len=*), intent(in) :: FILENAME
      integer, intent(in) :: ACCESS_MODE
    end function HE5_SWOPEN

    integer function HE5_SWINQSWATH (FILENAME,SWATHLIST,STRBUFSIZE)
      character (len=*), intent(in) :: FILENAME
      character (len=*), intent(out) :: SWATHLIST
      integer, intent(out):: STRBUFSIZE
    end function HE5_SWINQSWATH

    integer function HE5_SWINQDIMS (SWATHID,DIMNAME,DIMS)
       integer,intent(in)::SWATHID
       character(len=*),intent(out)::DIMNAME
       integer,intent(out),dimension(*)::DIMS
    
    end function HE5_SWINQDIMS
    
    integer function HE5_SWDIMINFO (SWATHID,DIMNAME)
       integer,intent(in)::SWATHID
       character(len=*),intent(IN)::DIMNAME
    end function HE5_SWDIMINFO

  end interface



!====================
end module HDFEOS5
!====================
