! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

!===========================================================================
module HDFEOS5               ! F90 interface to HDF-EOS5.
!===========================================================================
  implicit none
  ! use hdf5, only: hsize_t, hssize_t, size_t ! Someday should be import
  public

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

  !The following is included only because the original file
  ! $(HDFEOS5)/../../include/hdfeos5.inc isn't f95-compatible
  ! as of Toolkit version 5.2.8
  include 'hdfeos5.f9h'

  ! To switch to/from hdfeos5.1.6(+) swap comments between next 2 lines
  integer, parameter         :: MLS_CHARTYPE = HE5T_NATIVE_SCHAR
  ! integer, parameter         :: MLS_CHARTYPE = HE5T_CHARSTRING

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

  ! ------------------------------------------------------------------------
  ! Bugs and limitations: 
  ! (1) Due to misguided Fortran standard, we are not permitted module-level
  ! variables in the interfaces below. Fortran 2008 will permit
  ! the import command in place of use, but until then we tediously add
  ! multiple instances of use ..
  ! (2) (don't) use hdf5_params !
  ! some integers in HDF-EOS5 fortran interface are
  ! C long (see, e.g. SWapi.c). Hsize_t as declared in hdf5 is ignored.
  ! Right now we're use-ing size_t from the hdf5 library as an ad hoc type.
  ! A better plan would be to declare a module variable, either here 
  ! or in MLSKinds, corresponding to the hdfeos wrappers, calling
  ! it hdfeos5size_t
  ! integer, parameter, public         :: hdfeos5size_t = SIZE_T
  ! ------------------------------------------------------------------------
  
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
      use hdf5, only: hsize_t, hssize_t, size_t
      integer, intent(in) :: GDFID
      character (len=*), intent(in) :: GRIDNAME
       integer(kind=size_t),intent(in)::XDIMSIZE
       integer(kind=size_t),intent(in)::YDIMSIZE
       double precision,dimension(2),intent(in)::UPLEFT
       double precision,dimension(2),intent(in)::LOWRIGHT
    end function HE5_GDCREATE

    integer function HE5_GDDETACH ( GDID )
      integer, intent(in) :: GDID
    end function HE5_GDDETACH

    integer function HE5_GDINQGRID (FILENAME,GRIDLIST,STRBUFSIZE)
      use hdf5, only: hsize_t, hssize_t, size_t
      character (len=*), intent(in) :: FILENAME
      character (len=*), intent(out) :: GRIDLIST
      integer(kind=size_t), intent(out):: STRBUFSIZE
    end function HE5_GDINQGRID

    integer function HE5_GDINQDIMS (GRIDID,DIMNAME,DIMS)
      use hdf5, only: hsize_t, hssize_t, size_t
       integer,intent(in)::GRIDID
       character(len=*),intent(out)::DIMNAME
       integer(kind=size_t),intent(out),dimension(*)::DIMS
    
    end function HE5_GDINQDIMS
    
    integer function HE5_GDDIMINFO (GRIDID,DIMNAME)
       integer,intent(in)::GRIDID
       character(len=*),intent(IN)::DIMNAME
    end function HE5_GDDIMINFO

    integer function HE5_GDGRIDINFO (GRIDID,XDIMSIZE,YDIMSIZE,UPLEFT,LOWRIGHT)
      use hdf5, only: hsize_t, hssize_t, size_t
       integer,intent(in)::GRIDID
       integer(kind=size_t),intent(out)::XDIMSIZE
       integer(kind=size_t),intent(out)::YDIMSIZE
       double precision,dimension(2),intent(out)::UPLEFT
       double precision,dimension(2),intent(out)::LOWRIGHT
    end function HE5_GDGRIDINFO

    integer function HE5_GDFLDINFO (GRIDID,FIELDNAME,&
     & RANK,DIMS,NUMBERTYPE,DIMLIST, MAXDIMLIST)
      use hdf5, only: hsize_t, hssize_t, size_t
       integer,intent(in)::GRIDID
       character(len=*),intent(IN)::FIELDNAME
       integer(kind=size_t),intent(out),dimension(*):: DIMS
       integer,intent(out)             :: RANK, NUMBERTYPE
       character(len=*),intent(OUT)::DIMLIST
       character(len=*),intent(OUT)::MAXDIMLIST
    end function HE5_GDFLDINFO

    integer function HE5_GDNENTRIES (GRIDID, ENTRYCODE, STRINGBUFFERSIZE)
      use hdf5, only: hsize_t, hssize_t, size_t
       integer,intent(in)::GRIDID
       integer,intent(in)::ENTRYCODE
       integer(kind=size_t),intent(out)::STRINGBUFFERSIZE
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
      use hdf5, only: hsize_t, hssize_t, size_t
      integer, intent(in) :: SWATHID,CHUNKRANK
      integer(kind=size_t), intent(in),dimension(*) :: CHUNKDIMS
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
      use hdf5, only: hsize_t, hssize_t, size_t
      integer, intent(in)  :: SWATHID
      character (len=*), intent(in) :: FIELDNAME
      integer(kind=size_t), intent(in) :: DIM
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
      use hdf5, only: hsize_t, hssize_t, size_t
      character (len=*), intent(in) :: FILENAME
      character (len=*), intent(out) :: SWATHLIST
      integer(kind=size_t), intent(out):: STRBUFSIZE
    end function HE5_SWINQSWATH

    integer function HE5_SWINQDIMS (SWATHID,DIMNAME,DIMS)
      use hdf5, only: hsize_t, hssize_t, size_t
       integer,intent(in)::SWATHID
       character(len=*),intent(out)::DIMNAME
       integer(kind=size_t),intent(out),dimension(*)::DIMS
    end function HE5_SWINQDIMS
    
    integer function HE5_SWDIMINFO (SWATHID,DIMNAME)
       integer,intent(in)::SWATHID
       character(len=*),intent(IN)::DIMNAME
    end function HE5_SWDIMINFO

    integer function HE5_SWFLDINFO (SWATHID,FIELDNAME,&
     & RANK,DIMS,NUMBERTYPE,DIMLIST, MAXDIMLIST)
      use hdf5, only: hsize_t, hssize_t, size_t
       integer,intent(in)::SWATHID
       character(len=*),intent(IN)::FIELDNAME
       integer(kind=size_t),intent(out),dimension(*):: DIMS
       integer,intent(out),dimension(*):: NUMBERTYPE
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
!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module HDFEOS5
!====================

! $Log$
! Revision 2.13  2009/09/29 23:33:12  pwagner
! Changes needed by 64-bit build
!
! Revision 2.12  2009/06/23 18:25:42  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.11  2005/06/22 17:25:48  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 2.10  2004/03/24 23:52:11  pwagner
! Added MLS_CHARTYPE
!
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
