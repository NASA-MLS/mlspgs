! Copyright (c) 2004, California Institute of Technology.  ALL RIGHTS RESERVED.
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

  !The following is included only because the original file
  ! $(HDFEOS5)/../../include/hdfeos5.inc isn't f95-compatible
  ! as of Toolkit version 5.2.8
  include 'hdfeos5.f9h'

  ! Now define f90 interfaces for some HDF-EOS5.
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
    integer function HE5_EHINQGLATTS ( FILE_ID, ATTRLIST, LISTSIZE )
      integer, intent(in)              :: FILE_ID
      character(len=*), intent(out)    :: attrlist
      integer, intent(out)             :: listsize
    end function HE5_EHINQGLATTS

    integer function HE5_GDCLOSE ( FILE_ID )
      integer, intent(in) :: FILE_ID
    end function HE5_GDCLOSE

    integer function HE5_GDOPEN ( FILENAME, ACCESS_MODE )
      character (len=*), intent(in) :: FILENAME
      integer, intent(in) :: ACCESS_MODE
    end function HE5_GDOPEN

    integer function HE5_GDATTACH ( GDFID, GRIDNAME )
      integer, intent(in) :: GDFID
      character (len=*), intent(in) :: GRIDNAME
    end function HE5_GDATTACH

    integer function HE5_GDCREATE ( GDFID, GRIDNAME,&
      & XDIMSIZE,YDIMSIZE,UPLEFT,LOWRIGHT )
      integer, intent(in) :: GDFID
      character (len=*), intent(in) :: GRIDNAME
       integer,intent(in)::XDIMSIZE
       integer,intent(in)::YDIMSIZE
       double precision,dimension(2),intent(in)::UPLEFT
       double precision,dimension(2),intent(in)::LOWRIGHT
    end function HE5_GDCREATE

    integer function HE5_GDDETACH ( GDID )
      integer, intent(in) :: GDID
    end function HE5_GDDETACH

    integer function HE5_GDINQGRID (FILENAME,GRIDLIST,STRBUFSIZE)
      character (len=*), intent(in) :: FILENAME
      character (len=*), intent(out) :: GRIDLIST
      integer, intent(out):: STRBUFSIZE
    end function HE5_GDINQGRID

    integer function HE5_GDINQDIMS (GRIDID,DIMNAME,DIMS)
       integer,intent(in)::GRIDID
       character(len=*),intent(out)::DIMNAME
       integer,intent(out),dimension(*)::DIMS
    
    end function HE5_GDINQDIMS
    
    integer function HE5_GDDIMINFO (GRIDID,DIMNAME)
       integer,intent(in)::GRIDID
       character(len=*),intent(IN)::DIMNAME
    end function HE5_GDDIMINFO

    integer function HE5_GDGRIDINFO (GRIDID,XDIMSIZE,YDIMSIZE,UPLEFT,LOWRIGHT)
       integer,intent(in)::GRIDID
       integer,intent(out)::XDIMSIZE
       integer,intent(out)::YDIMSIZE
       double precision,dimension(2),intent(out)::UPLEFT
       double precision,dimension(2),intent(out)::LOWRIGHT
    end function HE5_GDGRIDINFO

    integer function HE5_GDFLDINFO (GRIDID,FIELDNAME,&
     & RANK,DIMS,NUMBERTYPE,DIMLIST, MAXDIMLIST)
       integer,intent(in)::GRIDID
       character(len=*),intent(IN)::FIELDNAME
       integer,intent(out),dimension(*):: DIMS
       integer,intent(out)             :: RANK, NUMBERTYPE
       character(len=*),intent(OUT)::DIMLIST
       character(len=*),intent(OUT)::MAXDIMLIST
    end function HE5_GDFLDINFO

    integer function HE5_GDNENTRIES (GRIDID, ENTRYCODE, STRINGBUFFERSIZE)
       integer,intent(in)::GRIDID
       integer,intent(in)::ENTRYCODE
       integer,intent(out)::STRINGBUFFERSIZE
    end function HE5_GDNENTRIES
    
    integer function HE5_GDPROJINFO (GRIDID,&
      & PROJCODE,ZONECODE,SPHERECODE,PROJPARM)
       integer,intent(in)::GRIDID
       integer,intent(out)::PROJCODE
       integer,intent(out)::ZONECODE
       integer,intent(out)::SPHERECODE
       double precision,dimension(*),intent(out)::PROJPARM
    end function HE5_GDPROJINFO

    integer function HE5_GDINQFLDS (GRIDID, FIELDLIST, RANK, NUMBERTYPE)
       integer,intent(in)::GRIDID
       character(len=*),intent(out)::FIELDLIST
       integer,intent(out),dimension(*)::RANK, NUMBERTYPE
    end function HE5_GDINQFLDS
    
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

    integer function HE5_SWFLDINFO (SWATHID,FIELDNAME,&
     & RANK,DIMS,NUMBERTYPE,DIMLIST, MAXDIMLIST)
       integer,intent(in)::SWATHID
       character(len=*),intent(IN)::FIELDNAME
       integer,intent(out),dimension(*):: DIMS, NUMBERTYPE
       integer,intent(out)             :: RANK
       character(len=*),intent(OUT)::DIMLIST
       character(len=*),intent(OUT)::MAXDIMLIST
    end function HE5_SWFLDINFO

    integer function HE5_SWINQDFLDS (SWATHID, FIELDLIST, RANK, NUMBERTYPE)
       integer,intent(in)::SWATHID
       character(len=*),intent(out)::FIELDLIST
       integer,intent(out),dimension(*)::RANK, NUMBERTYPE
    end function HE5_SWINQDFLDS
    
    integer function HE5_SWSETALIAS (SWATHID,FIELDNAME,ALIASLIST)
       integer,intent(in)::SWATHID
       character(len=*),intent(IN)::FIELDNAME
       character(len=*),intent(IN)::ALIASLIST
    end function HE5_SWSETALIAS

  end interface

!====================
contains 
  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module HDFEOS5
!====================

! $Log$
! Revision 2.9  2004/02/26 21:58:23  pwagner
! Added HE5_EHINQGLATTS
!
! Revision 2.8  2003/09/03 22:38:42  pwagner
! Added HE5_SWFLDINFO
!
! Revision 2.7  2003/06/03 20:39:39  pwagner
! Added some grid apis
!
! Revision 2.6  2003/04/21 19:30:07  pwagner
! Added interface to HE5_SWINQDFLDS
!
! Revision 2.5  2003/04/11 23:32:23  pwagner
! Moved he5_swsetfill he5_ehwrglatt interfaces
!
! Revision 2.4  2003/03/07 00:34:18  pwagner
! Added HE5_SWSETFILL
!
! Revision 2.3  2002/12/11 22:21:47  pwagner
! Added HE5_SWSETALIAS
!
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
