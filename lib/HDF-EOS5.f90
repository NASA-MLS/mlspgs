! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

!===========================================================================
module HDFEOS5               ! F90 interface to HDF-EOS5.
!===========================================================================
  implicit none
  public

  !------------------- RCS Ident Info -----------------------
  character(len=130), private :: Id = &
    & "$Id$"
  character(len=*), parameter, private :: ModuleName= &
    & "$RCSfile$"
  private :: not_used_here 
  !----------------------------------------------------------

  !The following is included only because hte original file
  ! $(HDFEOS5)/../../include/hdfeos5.inc isn't f95-compatible
  ! as of Toolkit version 5.2.8
  include 'hdfeos5.f9h'

  ! Now define f90 interfaces for some HDF-EOS.
  ! Warning: as you are calling C directly, make sure that your 
  ! array args are like this.
  !integer function HE5_SWGONKULATE(ARRAY)
  !  integer,intent(in),dimension(*)::ARRAY
  !end function HE5_SWGONKULATE!
  ! NOT like this.
  !integer function HE5_SWGONKULATE(ARRAY)
  !  integer,intent(in),dimension(:)::ARRAY
  !end function HE5_SWGONKULATE!  ^
  !                               |
  !                        OY! NO! NOT like this! No (:)s already!

  interface
    integer function HE5_GDCLOSE ( FILE_ID )
      integer, intent(in) :: FILE_ID
    end function HE5_GDCLOSE

    integer function HE5_GDOPEN ( FILENAME, ACCESS_MODE )
      character (len=*), intent(in) :: FILENAME
      integer, intent(in) :: ACCESS_MODE
    end function HE5_GDOPEN

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

    integer function HE5_SWDEFCHUNK(SWATHID,CHUNKRANK,CHUNKDIMS)
      ! use hdf5_params ! Note: all integers in HDF-EOS5 fortran interface are
      ! C long. Hsize_t is not known about.
      integer, intent(in) :: SWATHID,CHUNKRANK
      integer, intent(in),dimension(*) :: CHUNKDIMS
    end function HE5_SWDEFCHUNK

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
contains 
  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module HDFEOS5
!====================

! $Log$
! Revision 2.2  2002/10/08 00:09:09  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.1  2002/07/11 22:18:26  pwagner
! First commit in this directory--welcome old friendshe5*.f90
!
! Revision 1.4  2002/05/28 23:11:23  pwagner
! Changed to comply with hdf5.1.4.3/hdfeos5.1.2
!
! Revision 1.3  2002/01/29 00:49:26  pwagner
! Added he5_gd(open)(close) functions
!
