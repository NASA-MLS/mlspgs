module AntennaPatterns_m

  ! Read the antenna patterns file.

  use Allocate_Deallocate, only: Allocate_Test, Deallocate_Test
  use MLSCommon, only: R8
  use MLSMessageModule, only: MLSMessage, MLSMSG_Allocate, MLSMSG_DeAllocate, &
    & MLSMSG_Error
  use MLSSignals_m, only: MaxSigLen, Signals

  ! More USEs below in each procedure, if they're only used therein.

  private
  ! Public procedures:
  public :: Open_Antenna_Patterns_File, Read_Antenna_Patterns_File
  public :: Close_Antenna_Patterns_File
! public :: Close_Antenna_Patterns_File, Destroy_Antenna_Patterns_Database
  public :: Dump_Antenna_Patterns_Database

  type, public :: AntennaPattern_T
    real(r8) :: Lambda
    real(r8), dimension(:), pointer :: Aaap => NULL()
    real(r8), dimension(:), pointer :: D1aap => NULL()
    real(r8), dimension(:), pointer :: D2aap => NULL()
    character(len=MaxSigLen), pointer, dimension(:) :: Signals => NULL()
  end type AntennaPattern_T

  ! The filter shape database:
  type(AntennaPattern_T), dimension(:), pointer, public :: AntennaPatterns => NULL()

  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
    & "$Id$"
  character (len=len(idParm)) :: Id = idParm
  character (len=*), parameter, private :: ModuleName = &
    & "$RCSfile$"
  !---------------------------------------------------------------------------

contains

  ! ------------------------------------  Open_Antenna_Patterns_File  -----
  subroutine Open_Antenna_Patterns_File ( Filename, Lun )

    character(len=*), intent(in) :: Filename ! Name of the filter shape file
    integer, intent(out) :: Lun              ! Logical unit number to read it

    logical :: Exist, Opened
    integer :: Status

    do lun = 20, 99
      inquire ( unit=lun, exist=exist, opened=opened )
      if ( exist .and. .not. opened ) exit
    end do
    if ( opened .or. .not. exist ) call MLSMessage ( MLSMSG_Error, moduleName, &
      & "No logical unit numbers available" )
    open ( unit=lun, file=filename, status='old', form='formatted', &
      & access='sequential', iostat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
      & "Unable to open filter shapes file " // Filename )
  end subroutine Open_Antenna_Patterns_File

  ! ------------------------------------  Read_Antenna_Patterns_File  -----
  subroutine Read_Antenna_Patterns_File ( Lun, Spec_Indices )
    use Machine, only: IO_Error
    use Parse_Signal_m, only: Parse_Signal
    use Toggles, only: Gen, Levels, Switches, Toggle
    use Trace_M, only: Trace_begin, Trace_end
    use Units, only: Pi

    integer, intent(in) :: Lun               ! Logical unit number to read it
    integer, intent(in) :: Spec_Indices(:)   ! Needed by Parse_Signal, q.v.

    real(r8), parameter :: Pi2 = 2.0_r8 * Pi

    integer :: DataBaseSize                  ! How many antenna patterns?
    integer :: HowManyPoints(size(signals))  ! for each pattern
    integer :: HowManySignals(size(signals)) ! for each pattern
    integer :: I, J, K, L, N                 ! Loop inductors, subscripts
    real(r8) :: Lambda
    real(r8) :: LambdaX2Pi                   ! 2 * Pi * Lambda
    real(r8) :: Q                            ! Factor used to scale derivatives
    character(len=MaxSigLen) :: SigName      ! Signal Name
    integer :: Status                        ! From read or allocate
    integer, pointer, dimension(:) :: Signal_Indices => NULL()   ! From Parse_Signal, q.v.
    real(r8) :: V(6)                         ! To read a line from the file

    if ( toggle(gen) ) call trace_begin ( "Read_Antenna_Patterns_File" )

!   if ( associated(AntennaPatterns) ) call destroy_antenna_patterns_database


    ! First, read through the file and count how much stuff is there.
    read ( lun, '(a)', end=98, err=99, iostat=status ) sigName
    dataBaseSize = 0
outer1: do
      dataBaseSize = dataBaseSize + 1
      if ( dataBaseSize > size(howManySignals) ) &
        & call MLSMessage ( MLSMSG_Error, moduleName, &
          & "More patterns in the file than signals in the database" )
      howManySignals(dataBaseSize) = 0
      do ! Count how many signals there are
        sigName = adjustl(sigName)
        if ( verify(sigName(1:1), '0123456789.+-') == 0 ) exit ! a number
        call parse_signal ( sigName, signal_indices, spec_indices )
        if ( .not. associated(signal_indices) ) &
          call MLSMessage ( MLSMSG_Error, moduleName, &
            & trim(sigName) // " is not a valid signal." )
        howManySignals(dataBaseSize) = howManySignals(dataBaseSize) + 1
        read ( lun, '(a)', end=98, err=99, iostat=status ) sigName
      end do
      read ( sigName, *, err=99, iostat=status ) lambda
      howManyPoints(dataBaseSize) = 0
      do ! Count how many points in the pattern
        read ( lun, '(a)', err=99, iostat=status ) sigName
        if ( status < 0 ) exit outer1
        sigName = adjustl(sigName)
        if ( verify(sigName(1:1), '0123456789.+-') /= 0 ) exit ! not a number
        read ( sigName, *, err=99, iostat=status ) v ! Check the numbers
        howManyPoints(dataBaseSize) = howManyPoints(dataBaseSize) + 1
      end do
    end do outer1

    ! Now we know how big the database is, and each part of it.
    allocate ( antennaPatterns(dataBaseSize), stat=status )
    if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
      & MLSMSG_Allocate // "AntennaPatterns" )
    rewind ( lun )
    do i = 1, dataBaseSize
      call Allocate_Test ( antennaPatterns(i)%signals, howManySignals(i), &
        & "AntennaPatterns(?)%Signals", moduleName )
      call Allocate_Test ( antennaPatterns(i)%aaap, 6*howManyPoints(i), &
        & "AntennaPatterns(?)%Aaap", moduleName )
      call Allocate_Test ( antennaPatterns(i)%d1aap, 6*howManyPoints(i), &
        & "AntennaPatterns(?)%D1aap", moduleName )
      call Allocate_Test ( antennaPatterns(i)%d2aap, 6*howManyPoints(i), &
        & "AntennaPatterns(?)%D2aap", moduleName )
      do j = 1, howManySignals(i)
        read ( lun, '(a)', err=99, iostat=status ) antennaPatterns(i)%signals(j)
      end do ! j
      read ( lun, *, err=99, iostat=status ) lambda
      antennaPatterns(i)%lambda = lambda
      lambdaX2Pi = pi2 * lambda
      do j = 1, howManyPoints(i)
        read ( lun, *, err=99, iostat=status ) v
        k = 2 * j - 1
        l = 2 * howManyPoints(i) + k
        n = 4 * howManyPoints(i) + k
        antennaPatterns(i)%aaap(k:k+1) = v(1:2)
        antennaPatterns(i)%aaap(l:l+1) = v(3:4)
        antennaPatterns(i)%aaap(n:n+1) = v(5:6)

        ! First derivative field:     i*Q * F(S), i = Sqrt(-1)

        q = (k-1) * lambdaX2Pi
        antennaPatterns(i)%d1aap(k)    = -v(2) * q
        antennaPatterns(i)%d1aap(k+1)  =  v(1) * q
        antennaPatterns(i)%d1aap(l)    = -v(4) * q
        antennaPatterns(i)%d1aap(l+1)  =  v(3) * q
        antennaPatterns(i)%d1aap(n)    = -v(6) * q
        antennaPatterns(i)%d1aap(n+1)  =  v(5) * q

        ! Second derivative field:    (i*Q)**2 * F(S), i = Sqrt(-1)

        q = -q * q
        antennaPatterns(i)%d2aap(k:k+1)  =  v(1:2) * q
        antennaPatterns(i)%d2aap(l:l+1)  =  v(3:4) * q
        antennaPatterns(i)%d2aap(n:n+1)  =  v(5:6) * q
      end do ! j
    end do ! i

    if ( toggle(gen) ) then
      if ( levels(gen) > 0 .or. index(switches,'A') /= 0 ) &
        & call dump_Antenna_Patterns_database
      call trace_end ( "Read_Antenna_Patterns_File" )
    end if

    return
  98 call MLSMessage ( MLSMSG_Error, moduleName, "Unexpected end-of-file" )
  99 call io_error ( "While reading the filter shape file", status )
     call MLSMessage ( MLSMSG_Error, moduleName, "Input error" )
  end subroutine Read_Antenna_Patterns_File

  ! -----------------------------------  Close_Antenna_Patterns_File  -----
  subroutine Close_Antenna_Patterns_File ( Lun )
    close ( lun )
  end subroutine Close_Antenna_Patterns_File

  ! -----------------------------  Destroy_Antenna_Patterns_Database  -----
! subroutine Destroy_Antenna_Patterns_Database
!   integer :: I, Status
!   do i = 1, size(AntennaPatterns)
!     call deallocate_test ( AntennaPatterns(i)%signals, &
!       & "AntennaPatterns(?)%Signals", moduleName )
!     call deallocate_test ( AntennaPatterns(i)%aaap, &
!       & "AntennaPatterns(?)%aaap", moduleName )
!     call deallocate_test ( AntennaPatterns(i)%d1aap, &
!       & "AntennaPatterns(?)%D1aap", moduleName )
!     call deallocate_test ( AntennaPatterns(i)%d2aap, &
!       & "AntennaPatterns(?)%D2aap", moduleName )
!   end do ! i
!   deallocate ( AntennaPatterns, stat=status )
!   if ( status /= 0 ) call MLSMessage ( MLSMSG_Error, moduleName, &
!     MLSMSG_DeAllocate // "AntennaPatterns" )
! end subroutine Destroy_Antenna_Patterns_Database

  ! --------------------------------  Dump_Antenna_Patterns_Database  -----
  subroutine Dump_Antenna_Patterns_Database
    use Dump_0, only: Dump
  use Output_m, only: Blanks, Output

    integer :: I, J                ! Subscripts, loop inductors
    call output ( 'Filter Shapes: SIZE = ' )
    call output ( size(AntennaPatterns), advance='yes' )
    do i = 1, size(AntennaPatterns)
      call output ( i, 4 )
      call output ( ':    Signal = ' )
      call output ( trim(antennaPatterns(i)%signals(1)), advance='yes' )
      do j = 2, size(antennaPatterns(i)%signals)
        call blanks ( 18 )
        call output ( trim(antennaPatterns(i)%signals(i)), advance='yes' )
      end do
      call output ( ' Lambda = ' )
      call output ( antennaPatterns(i)%lambda, advance='yes' )
      call dump ( antennaPatterns(i)%aaap, name='Aaap' )
      call dump ( antennaPatterns(i)%aaap, name='D1aap' )
      call dump ( antennaPatterns(i)%aaap, name='D2aap' )
    end do ! i
  end subroutine Dump_Antenna_Patterns_Database

end module AntennaPatterns_m

! $Log$
! Revision 1.1  2001/03/30 02:37:44  vsnyder
! Initial Commit
!
