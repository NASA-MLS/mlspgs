module DUMP_K_RECORDS_M
  use FORMAT_KEY_M, only: FORMAT_KEY
  use L2PC_FILE_PARAMETERS, only: MSVD => max_no_sv_derivatives
  use L2PC_FILE_STRUCTURES, only: DUMP_K_STAR, L2PC_HEADER_ONE, &
                                  L2PC_HEADER_TWO, L2PC_HEADER_TRI
  use L2PCDIM, only: MNP => max_no_phi
  use MLSCommon, only: I4, R4
  use STRINGS, only: SQZSTR
  implicit NONE
  private
  public :: DUMP_K_RECORDS

!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
!---------------------------------------------------------------------------

contains

!
  Subroutine DUMP_K_RECORDS ( jch, time_i, band, scale_factor, header1, &
 &           header2, header3, no_phi_vec, k_star_all, cycle, Ier )

    integer(i4), intent(in) :: JCH
    integer(i4), intent(in) :: TIME_I
    integer(i4), intent(in) :: BAND
    real(r4), intent(in) :: SCALE_FACTOR(*)
    type(l2pc_header_one), intent(in) :: HEADER1
    type(l2pc_header_two), intent(in) :: HEADER2
    type(l2pc_header_tri), intent(in) :: HEADER3 ! not used
    integer(i4), intent(in) :: NO_PHI_VEC(*)
    real(r4), intent(in) :: K_STAR_ALL(msvd,mnp,*)
    integer(i4), intent(inout) :: CYCLE
    integer(i4), intent(out) :: IER
!
    Character(len=8) :: A8
    character(len=2) :: ABAND
    Character(len=40) :: AKEY
    type(dump_k_star) :: DUMP_K
    Integer(i4), save :: EL1
    character(len=*), parameter :: FDUM = 'k_star_dummy.dat'
    Integer(i4) :: IZ, J
    Integer(i4), save :: JB1
    Integer(i4) :: JP, K
    Character(len=40) :: Key
    Integer(i4) :: N
    Integer(i4), save  :: NDR
    Integer(i4) :: NF, NO_C
    Integer(i4), save :: N1= -1, PB = -2, PJCH = 500
    Integer(i4), save :: PL1
    Integer(i4) :: SV_I
    Integer(i4), save :: TI = -1
    Real(r4) :: Q, R, S
!
    Ier = 0
    if (band /= pb) pjch = jch + 5
    if (ti /= time_i) pjch = jch + 5
!
    if (jch < pjch) then
      Ndr = 0
      pb = band
      pjch = jch
      ti = time_i
      cycle = -1
      inquire (iolength=j) dump_k
      Close (97, iostat=Ier)
      Open (97, file=Fdum, status='UNKNOWN', access='DIRECT', recl=j, &
   &        form='UNFORMATTED', iostat=Ier)
      if (Ier /= 0) Return
    endif
!
! Write the records by pointing and description
! Load in the state vector
!
    aband = '_b'
    write (aband(2:2),'(i1)') band
!
    k = (band - 1) * header1%no_channels_per_band
    dump_k%first_channel_number = k + 1
!
    dump_k%derivatives(1:header1%no_pointings) = 0.0
!
! *** FIRST Derivative records handling:
!
! Set k_star key:
!
    Key(1:) = ' '
    Key = 'D_10N_A_NA_THE_C_SVNAMEXXEL_PHI+NF_B'         ! 1st deriv.
    Key(01:08) = header1%avail_keys(time_i)
!
    if (n1 < 0) then
      el1 = 0
      pl1 = 0
      n1 = index(Key,'SVNAMEXX')
      if (n1 > 0) then
        el1 = n1 + 8
        pl1 = index(Key,'PHI+') + 3
      endif
      jb1 = index(Key,'_B')
    endif
!
    Key(jb1:jb1+1) = aband
!
    do k = 1, header1%no_sv_components
!
! Write the type(k_star) record
!
! Check if user expects this derivative is non-zero
!
      a8(1:8) =' '
      dump_k%sv_name1 = a8
      dump_k%sv_name2 = a8
!
      dump_k%log_psv1 = -5.0
      dump_k%log_psv2 = -5.0
!
      if (header1%sv_rtrvl_by_band(k,band)) then
!
        a8 = header1%sv_components(k)(1:8)
        Key(n1:n1+7) = a8(1:8)                   ! Update the key
!
! Record the alphanumeric SV name, get the first element index
!
        dump_k%sv_name1 = a8
        sv_i = header1%sv_component_first_elmnt_index(k) - 1
        iz = max(1,no_phi_vec(sv_i+1))
!
! Run over all coefficients within the SV element
!
        no_c = header1%no_elmnts_per_sv_component(k)
        do n = 1, no_c
!
! Set the element number with the key
!
          write (Key(el1:el1+1),'(i2.2)') n
!
          sv_i = sv_i + 1
          dump_k%log_psv1 = header2%tri_basis_vert_grid(sv_i)
!
! Loop over all applicable Phi's
!
          jp = (iz + 1) / 2
          do nf = 1, iz
!
! Set the Phi index number with the key
!
            j = nf - jp
            Key(pl1:pl1)='+'
            if (j < 0) Key(pl1:pl1) = '-'
!
!           write (Key(pl1+1:pl1+2),'(i2.2)') iabs(j)
!
            AKey(1:)=' '
            write (Akey,'(i2.2)') iabs(j)
            Key(pl1+1:pl1+2) = AKey(1:2)
!
            s = 1.0
            AKey(1:)=' '
            AKey = Key
            Call SqzStr(AKey)
            j = index(AKey,'PHI')-4
            if (AKey(j-1:j) == '_W' .or. AKey(j-1:j) == '_N' .or. &
   &           AKey(j-1:j) == '_V') s = 0.0
!
! Create derivative, scale the derivatives
!
            do j = 1, header1%no_pointings
              r = scale_factor(j) * s + 1.0 - s
              q = k_star_all(sv_i,nf,j) * r
              dump_k%derivatives(j) = q
            end do
!
            AKey(1:)=' '
            AKey = Key
            Call format_key(AKey,Ier)
            if (Ier /= 0) Return
            dump_k%Key = AKey
!
            Ndr = Ndr + 1
            write (97,rec=Ndr) dump_k
!
          end do
!
        end do                     ! SV coefficients (n)
!
      endif
!
    end do                         ! SV cycle (k)
!
    if (cycle < 1) cycle = Ndr
!
    Return
  End Subroutine DUMP_K_RECORDS
end module DUMP_K_RECORDS_M

! $Log$
