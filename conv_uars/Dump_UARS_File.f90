! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

program Dump_UARS_file

  use Rad_File_Contents, only: Lvl1_Hdr_t, Limb_Hdr_t, Limb_Stat_t, Limb_Rad_t, &
                             & Limb_OA_t, Pad_T
  use Swap_OA_Rec_m, only: Swap_Lvl1_Hdr_Rec
  use UARS_Dumps, only: Dump_Limb_Hdr, Dump_Lvl1_Hdr, Dump_Limb_OA, &
    & Dump_Limb_Rad, Dump_Limb_Stat

  implicit NONE

  type(lvl1_hdr_t) :: lvl1_hdr
  type(limb_hdr_t) :: limb_hdr
  type(limb_stat_t) :: limb_stat
  type(limb_rad_t) :: limb_rad
  type(limb_oa_t) :: limb_oa
  type(pad_t) :: pad

  character(255) :: IoMSG, Line
  integer :: Details, MAFs, N, Recl, Stat, Total
  character :: YN ! for yes or no input

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: ModuleName= &
       "$RCSfile$"
  !---------------------------------------------------------------------------

  inquire ( iolength = recl ) limb_hdr, limb_stat, limb_rad, limb_oa, pad

  call get_command_argument ( 1, line )
  do
    if ( line == '' ) then
      write ( *, '(a)', advance='no' ) 'Enter UARS file name to dump: '
      read ( *, '(a)', end=9 ) line
    end if
    open ( 10, file=trim(line), access='direct', status='old', &
         & form='unformatted', recl=recl, iostat=stat, iomsg=iomsg )
    if ( stat /= 0 ) then
      print '(3a,i0/2a)', 'Unable to open "', trim(line), '" IOSTAT = ', stat, &
                      & 'IOMSG = ', trim(iomsg)
      line = ''
      cycle
    end if

    read ( 10, rec=1, iostat=stat, iomsg=iomsg ) lvl1_hdr
    if ( stat /= 0 ) then
      print '(a,i0,a,i0,2a)', 'Unable to read record ', 1, &
        & ' IOSTAT = ', stat, ' IOMSG = ', trim(iomsg)
      cycle
    end if
    call swap_lvl1_hdr_rec ( lvl1_hdr )
    if ( lvl1_hdr%rec_type == 'H' ) then
      mafs = lvl1_hdr%nummmaf
      print '(2(a,i0))', 'Number of MAFS ', mafs, ' Total number of records ', mafs+4
      exit
    end if
    print '(a)', 'First record appears not to be a UARS file header record'
  end do

  do

    write ( *, '(a)', advance='no' ) 'Enter record number to dump: '
    read ( *, *, iostat=stat, iomsg=iomsg, end=8 ) n
    if ( stat /= 0 ) then
      print '(a/a,i0,2a)', 'Improper record number, probably not numeric.', &
        & 'IOSTAT = ', stat, ' IOMSG = ', trim(iomsg)
      cycle
    end if
    read ( 10, rec=n, iostat=stat, iomsg=iomsg ) lvl1_hdr
    if ( stat /= 0 ) then
      print '(a,i0,a/a,i0,2a)', 'Unable to read record ', n, &
        & ', probably out of range ', &
        & 'IOSTAT = ', stat, ' IOMSG = ', trim(iomsg)
      cycle
    end if
    select case ( lvl1_hdr%rec_type )
    case ( 'H' ) ! Header record
      call dump_lvl1_hdr ( lvl1_hdr, .true. )
    case ( 'D' ) ! Data record
      print '(a)', 'Data record'
      read ( 10, rec=n, iostat=stat, iomsg=iomsg ) &
        & limb_hdr, limb_stat, limb_rad, limb_oa
      call dump_limb_hdr ( limb_hdr, .true. )
      write ( *, '(a)', advance='no' ) 'Limb stat details (0 = nothing): '
      read ( *, *, end=8 ) details
      call dump_limb_stat ( limb_stat, .true., details )
      write ( *, '(a)', advance='no' ) 'Dump 90x2x32 radiances (y/n)? '
      read ( *, *, end=8 ) yn
      if ( yn == 'Y' .or. yn == 'y' ) call dump_limb_rad ( limb_rad, .true. )
      write ( *, '(a)', advance='no' ) 'Dump OA (y/n)? '
      read ( *, *, end=8 ) yn
      if ( yn == 'Y' .or. yn == 'y' ) call dump_limb_oa ( limb_oa, .true. )
    case ( 'T' ) ! Trailer record
      call dump_lvl1_hdr ( lvl1_hdr, .true. )
    case default ! What is it?
    end select

  end do

8 print *
9 continue

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

end program Dump_UARS_file

! $Log$
