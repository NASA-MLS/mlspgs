Module Bill_GasAbsorption

  use MLSCommon, only: r8, rp, ip
  use L2PC_PFA_STRUCTURES, only: SLABS_STRUCT
  use SpectroscopyCatalog_m, only: CATALOG_T, LINES
  use CREATE_BETA_M, only: CREATE_BETA
  Implicit NONE
  private
  PUBLIC :: get_beta_bill

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter, private :: IdParm = &
  "$Id$"
  character (len=len(idParm)) :: Id = IdParm
  character (len=*), parameter, private :: ModuleName= &
  "$RCSfile$"
!---------------------------------------------------------------------------
contains

  SUBROUTINE get_beta_bill (T, PB, F, RH, VMR, ABSC, NS, Catalog )

!==============================================================
!      CALCULATE CLEAR-SKY ABSORPTION COEFFICIENT AT F AND T
!      USING BILL'S NEW SPECTROSCOPY DATA. 
!      LATEST UPDATE: J.JIANG, NOVEMBER 10, 2001
!==============================================================

    !-----------------
    ! INPUTS
    !-----------------
    INTEGER :: NS

    REAL(r8) :: F                            ! FREQUENCY IN GHz
    REAL(r8) :: FF                           ! FREQUENCY IN MHz
    REAL(r8) :: T                            ! TEMPERATURE (K)
    REAL(r8) :: P                            ! DRY AIR PARTIAL PRESSURE (hPa)
    REAL(r8) :: PB                           ! TOTAL AIR PRESSURE (hPa)
!    REAL(r8) :: VP                           ! VAPOR PARTIAL PRESSURE (hPa)
    REAL(r8) :: VMR(NS)                      ! MINOR SPECIES 1-O3
    REAL(r8) :: VMR_H2O                      ! H2O VOLUME MIXING RATIO
    REAL(r8) :: VMR_O2                       ! O2 VOLUME MIXING RATIO

    Type(Catalog_T), INTENT(IN) :: Catalog(:)
    Type (slabs_struct), DIMENSION(:), POINTER :: gl_slabs

    !------------------
    ! OUTPUTS
    !------------------
    REAL(r8), INTENT(out) :: ABSC            ! ABSORPTION COEFFICIENT (1/m)

    !------------------
    ! LOCAL VARIABLES
    !-----------------

    REAL(r8) :: B                            ! BETA (1/m/ppv)
    REAL(r8) :: RH                           ! Relative Humidity

    Integer(ip) :: n_sps, n_path, i, j, k, m, nl, no_of_lines, Spectag, status
    REAL(rp) :: bb, vp, v0, vm, tm, tp, bp, bm
    REAL(rp), allocatable, dimension(:) :: LineWidth
!-----------------------------------------------------------------------------

    B=0._r8
    FF = F*1000._r8
    n_sps = Size(Catalog)
    allocate ( gl_slabs (n_sps), stat=status )
    DO i = 1, n_sps
      Spectag = Catalog(i)%Spec_Tag
      no_of_lines =  size(Catalog(i)%Lines)
      gl_slabs(i)%no_lines =  no_of_lines
      Allocate(LineWidth(no_of_lines))
      do k = 1, no_of_lines
        m = Catalog(i)%Lines(k)
        LineWidth(k) = Lines(m)%W
      end do

      CALL create_beta(Spectag, Catalog(i)%continuum, PB, T, &
        &  FF, no_of_lines, LineWidth, gl_slabs(i)%v0s, gl_slabs(i)%x1,  &
        &  gl_slabs(i)%y, gl_slabs(i)%yi, gl_slabs(i)%slabs1,  bb,       &
        &  gl_slabs(i)%dslabs1_dv0, DBETA_DW=v0,DBETA_DN=vp,DBETA_DV=vm )

      B = B + bb

      DEAllocate(LineWidth)
      DEAllocate(gl_slabs)
    ENDDO

    !  spec_tags(l_h2o)          / 00018003 /
    !  spec_tags(l_o2)           / 00032001 /
    !  spec_tags(l_o_18_o)       / 00034001 /
    !  spec_tags(l_h2o_18)       / 00020003 /
    !  spec_tags(l_o3)           / 00048004 /

  END SUBROUTINE get_beta_bill 
 
End Module Bill_GasAbsorption

! $Log$




