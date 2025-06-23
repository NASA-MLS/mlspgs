! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module SLABS_SW_M

  ! Single-Line Absorption Software

  use MLSKinds, only: R8, RP
  use SpectroscopyCatalog_m, only: Catalog_T
  use Constants, only: SqrtPi

  implicit NONE

  private

!--------------------------------------------------  SLABS_STRUCT  -----
! This structure contains the "slabs preps arrays."  These are the
! frequency-independent terms in the cross section.

  type, public :: SLABS_STATE ! State stuff in Slabs_Struct
    real(r8) :: v0
    real(r8) :: v0s
    real(r8) :: x1
    real(r8) :: y
    real(r8) :: yi
    real(r8) :: slabs1
    real(r8) :: dslabs1_dv0 ! / slabs1
    logical :: polarized
  end type

  type, public :: SLABS_DERIV ! Derivative stuff in Slabs_Struct
    ! Contribution of dx1_dv0 and dy_dv0 in d Beta / d Nu0 cancel.
    ! See Slabs_DSpectral.
!   real(r8) :: dx1_dv0
!   real(r8) :: dy_dv0
    real(r8) :: dv0s_dT    ! not * 1 / v0s
    real(r8) :: dx1_dT     ! / x1
    real(r8) :: dy_dT      ! / y
    real(r8) :: dyi_dT     ! / yi
    real(r8) :: dslabs1_dT ! / slabs1
  end type

  type, public :: SLABS_STRUCT
    type(catalog_t), pointer :: Catalog ! everything else is same size
    !                                     as catalog%lines
    logical :: UseYi ! Are any d(:)%yi > 0?
    type(slabs_state), dimension(:), allocatable :: S ! State stuff
    type(slabs_deriv), dimension(:), allocatable :: D ! Derivative stuff
  end type SLABS_STRUCT

  ! Routines to manipulate slabs structs

  interface DUMP
    module procedure Dump_Slabs_Struct, Dump_Slabs_Struct_2D
  end interface

  public :: AllocateOneSlabs, AllocateSlabs
  public :: DeallocateAllSlabs, DeallocateOneSlabs
  public :: DestroyCompleteSlabs
  public :: Dump, Dump_Slabs_Struct, Dump_Slabs_Struct_2D

  ! Routines to compute Betas and their derivatives, and to get ready to do so:
  public :: Get_GL_Slabs_Arrays
  public :: Slabs, Slabs_dSpectral, Slabs_dAll, Slabs_dT, Slabs_Lines
  public :: Slabs_Lines_dAll, Slabs_lines_dSpectral, Slabs_Lines_dT
  public :: Slabs_Lines_dT_path
  public :: Slabs_Prep, Slabs_Prep_Struct, Slabs_Prep_Struct_Offset, Slabs_Prep_dT
  public :: Slabswint, Slabswint_dAll, Slabswint_dT
  public :: Slabswint_Lines, Slabswint_Lines_dAll, Slabswint_Lines_dSpectral
  public :: Slabswint_Lines_dT
  public :: Voigt_Lorentz, DVoigt_Spectral, DVoigt_Spectral_Lines

  real(rp), parameter :: OneOvSPi = 1.0_rp / sqrtPi  ! 1.0/Sqrt(Pi)

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: ModuleName= &
    & "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------
contains
!---------------------------------------------------------------------------

  ! -------------------------------------------  AllocateOneSlabs  -----
  subroutine AllocateOneSlabs ( Slabs, Catalog, InName, TempDer )
    ! Allocates the items in a slabs structure
    use Allocate_Deallocate, only: Test_Allocate
    use, intrinsic :: ISO_C_Binding, only: C_Intptr_t, C_Loc
    type (slabs_struct), intent(inout), target :: slabs ! Slabs to allocate
    type (catalog_t), target, intent(in) :: Catalog
    character(len=*), intent(in) :: InName      ! Who wants it
    logical, intent(in), optional :: TempDer    ! "Allocate temperature
                                                !  derivative fields"

    ! Local variables
    integer(c_intptr_t) :: Addr         ! For tracing
    logical :: MyDer
    integer :: NL, Stat

    ! Executable code
    myDer = .false.
    if ( present(tempDer) ) myDer = tempDer
    if ( associated(catalog%lines) ) then
      nl = size(catalog%lines)
    else
      nl = 0
    end if

    slabs%catalog => catalog
    allocate ( slabs%s(nl), stat=stat )
    addr = 0
    if ( stat == 0 .and. nl > 0 ) addr = transfer(c_loc(slabs%s(1)), addr)
    call test_allocate ( stat, inName, "Slabs%S", (/ 1 /), (/ nl /), &
      & storage_size(slabs%s) / 8, address=addr )
    if ( myDer ) then
      allocate ( slabs%d(nl), stat=stat )
      if ( stat == 0 .and. nl > 0 ) addr = transfer(c_loc(slabs%d(1)), addr)
      call test_allocate ( stat, inName, "Slabs%D", (/ 1 /), (/ nl /), &
      & storage_size(slabs%d) / 8, address=addr )
    end if
    if ( nl /= 0 ) then
      slabs%s%v0s = 0.0_r8
      slabs%s%x1 = 0.0_r8
      slabs%s%y = 0.0_r8
      slabs%s%yi = 0.0_r8
      slabs%s%slabs1 = 0.0_r8
      slabs%s%dslabs1_dv0 = 0.0_r8
      slabs%s%polarized = .false.
      if ( myDer ) then
        slabs%d%dv0s_dT = 0.0_r8
        slabs%d%dx1_dT = 0.0_r8
        slabs%d%dy_dT = 0.0_r8
        slabs%d%dyi_dT = 0.0_r8
        slabs%d%dslabs1_dT = 0.0_r8
      end if
    end if
  end subroutine AllocateOneSlabs

  ! --------------------------------------------  AllocateSlabs  ----------
  subroutine AllocateSlabs ( Slabs, No_Ele, Catalog, Caller, TempDer )
  ! Allocate an array of slabs structures, and then the items in each one

    use Allocate_Deallocate, only: Test_Allocate
    use, intrinsic :: ISO_C_Binding, only: C_Intptr_t, C_Loc

    type (slabs_struct), dimension(:,:), allocatable, target :: Slabs
    integer, intent(in) :: No_Ele
    type (catalog_t), dimension(:), target, intent(in) :: Catalog
    character(len=*), intent(in) :: Caller
    logical, intent(in), optional :: TempDer    ! "Allocate temperature
                                                !  derivative fields"

    integer(c_intptr_t) :: Addr         ! For tracing
    integer :: I, J

    allocate ( slabs(no_ele, size(catalog)), stat=i )
    addr = 0
    if ( i == 0 ) then
      if ( size(slabs) > 0 ) addr = transfer(c_loc(slabs(1,1)), addr)
    end if
    call test_allocate ( i, caller, 'Slabs', (/1,1/), (/no_ele,size(catalog)/), &
      & storage_size(slabs) / 8 )

    do i = 1, size(catalog)
      do j = 1, no_ele
        call AllocateOneSlabs ( slabs(j,i), catalog(i), caller, TempDer )
      end do
    end do

  end subroutine AllocateSlabs

  ! ------------------------------------------ DeallocateAllSlabs ---------
  subroutine DeallocateAllSlabs ( Slabs, inName )
    ! Allocates the items in a slabs
    type (slabs_struct), intent(inout), dimension(:,:) :: Slabs
    character(len=*), intent(in) :: InName

    integer :: I
    integer :: J
    ! Executable code
    do i = 1, size(slabs,2)
      do j = 1, size(slabs,1)
        call DeallocateOneSlabs ( slabs(j,i), inName )
      end do
    end do
  end subroutine DeallocateAllSlabs

  ! ------------------------------------------ DeallocateOneSlabs ---------
  subroutine DeallocateOneSlabs ( slabs, inName )
    ! DeAllocates the items in a slabs structure
    use Allocate_Deallocate, only: Test_DeAllocate
    use, intrinsic :: ISO_C_Binding, only: C_Intptr_t, C_Loc
    type (slabs_struct), intent(inout), target :: slabs ! Slabs to deallocate
    character (len=*), intent(in) :: inName ! ModuleName of caller

    integer(c_intptr_t) :: Addr         ! For tracing
    integer :: S, Stat

    ! Executable code
    if ( allocated(slabs%s) ) then
      s = size(slabs%s) * storage_size(slabs%s) / 8
      addr = 0
      if ( s > 0 ) addr = transfer(c_loc(slabs%s(1)), addr)
      deallocate ( slabs%s, stat=stat )
      call test_deallocate ( stat, inName, "Slabs%S", s, address=addr )
    end if
    if ( allocated(slabs%d) ) then
      s = size(slabs%d) * storage_size(slabs%d) / 8
      addr = 0
      if ( s > 0 ) addr = transfer(c_loc(slabs%d(1)), addr)
      deallocate ( slabs%d, stat=stat )
      call test_deallocate ( stat, inName, "Slabs%D", s, address=addr )
    end if
  end subroutine DeallocateOneSlabs

  ! ------------------------------------------- DestroyCompleteSlabs -----
  subroutine DestroyCompleteSlabs ( Slabs )
    ! Destroys all the components of a slabs
    use Allocate_Deallocate, only: Test_Deallocate
    use, intrinsic :: ISO_C_Binding, only: C_Intptr_t, C_Loc
    type (slabs_struct), dimension(:,:), allocatable, target :: Slabs

    integer(c_intptr_t) :: Addr         ! For tracing
    integer :: I, S
    ! Executable code
    call deallocateAllSlabs ( slabs, moduleName )
    s = size(slabs) * storage_size(slabs) / 8
    addr = 0
    if ( s > 0 ) addr = transfer(c_loc(slabs(1,1)), addr)
    deallocate ( slabs, stat=i )
    call test_deallocate ( i, moduleName, 'slabs', s, address=addr )
  end subroutine DestroyCompleteSlabs

  ! ------------------------------------------  Dump_Slabs_Struct  -----
  subroutine Dump_Slabs_Struct ( The_Slabs_Struct, Name )

    use Dump_0, only: Dump
    use Intrinsic, only: Lit_indices
    use Output_m, only: NewLine, Output
    use String_Table, only: Display_String

    type(slabs_struct), intent(in) :: The_Slabs_Struct
    character(len=*), intent(in), optional :: Name

    integer :: NL

    call output ( 'Slabs_Struct ' )
    if ( present(name) ) call output ( trim(name) )
    nl = size(the_slabs_struct%catalog%lines)
    if ( nl == 0 ) then
      call output ( ' is empty', advance='yes' )
    else
      call output ( '', advance='yes' )
      if ( the_slabs_struct%catalog%species_name /= 0 ) then
        call output ( 'Species ' )
        call display_string ( the_slabs_struct%catalog%species_name )
      end if
      call output ( 'Molecule ' )
      call display_string ( lit_indices(the_slabs_struct%catalog%molecule) )
      if ( the_slabs_struct%useYi ) call output ( ', use Yi' )
      call newLine
      call dump ( the_slabs_struct%s(:nl)%v0s, name='v0s' )
      call dump ( the_slabs_struct%s(:nl)%x1, name='x1' )
      call dump ( the_slabs_struct%s(:nl)%y, name='y' )
      call dump ( the_slabs_struct%s(:nl)%yi, name='yi' )
      call dump ( the_slabs_struct%s(:nl)%slabs1, name='slabs1' )
      call dump ( the_slabs_struct%s(:nl)%dslabs1_dv0, name='dslabs1_dv0' )
      if ( allocated (the_slabs_struct%d) ) then
        call dump ( the_slabs_struct%d(:nl)%dv0s_dT, name='dv0s_dT' )
        call dump ( the_slabs_struct%d(:nl)%dx1_dT, name='dx1_dT' )
        call dump ( the_slabs_struct%d(:nl)%dy_dT, name='dy_dT' )
        call dump ( the_slabs_struct%d(:nl)%dyi_dT, name='dyi_dT' )
        call dump ( the_slabs_struct%d(:nl)%dslabs1_dT, name='dslabs1_dT' )
      end if
    end if

  end subroutine Dump_Slabs_Struct

  ! ---------------------------------------  Dump_Slabs_Struct_2D  -----
  subroutine Dump_Slabs_Struct_2D ( The_Slabs_Struct, Name )

    use Output_m, only: Output

    type(slabs_struct), intent(in) :: The_Slabs_Struct(:,:)
    character(len=*), intent(in), optional :: Name

    integer :: I, J

    call output ( 'Slabs Struct' )
    if ( present(name) ) call output ( ' ' // trim(name) )
    call output ( ', SIZE = ' )
    call output ( size(the_slabs_struct,1) )
    call output ( ' X ' )
    call output ( size(the_slabs_struct,2), advance='yes' )
    do j = 1, size(the_slabs_struct,2)
      do i = 1, size(the_slabs_struct,1)
        call output ( 'Item ' )
        call output ( i )
        call output ( ', ' )
        call output ( j, advance='yes' )
        call dump ( the_slabs_struct(i,j) )
      end do
    end do

  end subroutine Dump_Slabs_Struct_2D

  ! --------------------------------------------  dVoigt_spectral  -----
  elemental subroutine dVoigt_spectral ( dNu, Nu0, x1, yi, y, w, t, tanh1, slabs1, SwI, &
                         &  dslabs1_dNu0, dSwI_dw, dSwI_dn, dSwI_dNu0 )

! Compute the Voigt function and its first derivatives with respect
! to spectral parameters: w, n & Nu0

! NOTE: Before calling this routine, the user needs to call slabs_prep()
!       routine to compute dslabs1_dNu0

! NOTE: In here and in all other routines in this module, 
!       tanh1 = tanh(h * nu / ( 2.0 * k * T ) )

    use Voigt_m, only: Simple_Voigt

    real(r8), intent(in) :: dnu, nu0
    real(rp), intent(in) :: x1, yi, y, w, t, tanh1, slabs1             
    real(rp), intent(in), optional :: dslabs1_dNu0 ! 1 / slabs1 * d(slabs1)/dNu0

    real(rp), intent(out) :: SwI, dSwI_dw,dSwI_dn
    real(rp), intent(out) :: dSwI_dNu0  ! 1 / Slabs1 * d Slabs1 / d Nu0 ???

    real(rp) :: x, u, v, du_dx, du_dy, dv_dx, dv_dy, q, b, g, z, r    

    real(rp) :: dx_dv0, du_dv0, dv_dv0, vvw, slabs2

    x = x1 * dNu                                                          
    call simple_voigt ( x, y, u, v )  

!  Van Vleck - Wieskopf (VVW) line shape with Voigt

    q = 1.0_rp + dNu / Nu0

    b = x1 * (2.0_r8 * Nu0 + dNu)
    g = b * b + y * y
    z = (y - b * yi) / g
    r = z * OneOvSPi + yi * v
    vvw = (u + r) * q
    slabs2 = slabs1 * tanh1
    SwI = slabs2 * vvw

    du_dx = 2.0_rp * (y * v - x * u)
    du_dy = 2.0_rp * (y * u + x * v - OneOvSPi)

    dv_dx = -du_dy         ! Cauchy-Riemann equation
    dv_dy =  du_dx         ! Cauchy-Riemann equation

! Compute the derivative of SwI w.r.t. w

    dSwI_dw = q * slabs2* (y/w) * (du_dy + yi*du_dx + &
                                &   OneOvSPi*((1.0_rp-2.0_rp*z*y)/g))

! Compute the derivative of SwI w.r.t. n

    dSwI_dn = q * slabs2 * y * Log(3.0d2/t) * (du_dy + yi * dv_dy)

! Finaly, compute the derivative of SwI w.r.t. Nu0

! ***** Analytically *****

!    dq_dv0 = -(Nu0+dNu)/(Nu0*Nu0)
    dx_dv0 = -x1
    du_dv0 = du_dx * dx_dv0
    dv_dv0 = dv_dx * dx_dv0
!    db_dv0 = x1
!    dg_dv0 = 2.0_rp * b*db_dv0
!    dz_dv0 = (-yi*db_dv0-z*dg_dv0)/g
!    dr_dv0 = dz_dv0*OneOvSPi+yi*dv_dv0
!    dvvw_dv0 = (du_dv0+dr_dv0)*q + dq_dv0*(u+r)
!    dvvw_dv0 = (du_dv0+dr_dv0)*q
    dSwI_dNu0 = slabs2*q*(du_dv0 + yi * dv_dv0) - swi / Nu0
    if (present(dslabs1_dNu0)) dSwI_dNu0 = dSwI_dNu0 + swi * dslabs1_dNu0

  end subroutine dVoigt_spectral

  ! --------------------------------------  dVoigt_spectral_Lines  -----
  elemental subroutine dVoigt_spectral_Lines ( dNu, Slabs, t, tanh1, &
                                  &  SwI, dSwI_dw, dSwI_dn, dSwI_dNu0 )

! Compute the sums of the Voigt function and its first derivatives with respect
! to spectral parameters: w, n & Nu0 for all lines in Slabs.

! NOTE: Before calling this routine, the user needs to call slabs_prep()
!       routine to compute dslabs1_dNu0

! NOTE: In here and in all other routines in this module, 
!       tanh1 = tanh( h * nu / ( 2.0 * k * T ) )

    use SpectroscopyCatalog_m, only: Lines
    use Voigt_m, only: Simple_Voigt

    real(r8), intent(in) :: dnu
    type(slabs_struct), intent(in) :: Slabs
    real(rp), intent(in) :: T
    real(rp), intent(in) :: Tanh1

    real(rp), intent(out) :: SwI, dSwI_dw, dSwI_dn, dSwI_dNu0               

    integer :: L  ! Line index

    real(r8) :: Nu0
    real(rp) :: x1, yi, y, w, slabs1             
    real(rp) :: dslabs1_dNu0                        

    real(rp) :: x, u, v, SwI1, du_dx, du_dy, dv_dx, dv_dy, q, b, g, z, r    

    real(rp) :: dx_dv0, du_dv0, dv_dv0, vvw, slabs2

    SwI = 0.0_rp
    dSwI_dw = 0.0_rp
    dSwI_dn = 0.0_rp
    dSwI_dNu0 = 0.0_rp     
    do l = 1, size(slabs%s)

      nu0 = slabs%s(l)%v0s
      x1 = slabs%s(l)%x1
      yi = slabs%s(l)%yi
      y = slabs%s(l)%y
      w = lines(slabs%catalog%lines(l))%w
      slabs1 = slabs%s(l)%slabs1
      dslabs1_dNu0 = slabs%s(l)%dslabs1_dv0

      x = x1 * dNu                                                          
      call simple_voigt ( x, y, u, v )  

  !  Van Vleck - Wieskopf (VVW) line shape with Voigt

      q = 1.0_rp + dNu / Nu0

      b = x1 * (2.0_r8 * Nu0 + dNu)
      g = b * b + y * y
      z = (y - b * yi) / g
      r = z * OneOvSPi + yi * v
      vvw = (u + r) * q
      slabs2 = slabs1 * tanh1
      SwI1 = slabs2 * vvw

      du_dx = 2.0_rp * (y * v - x * u)
      du_dy = 2.0_rp * (y * u + x * v - OneOvSPi)

      dv_dx = -du_dy         ! Cauchy-Riemann equation
      dv_dy =  du_dx         ! Cauchy-Riemann equation

  ! Compute the derivative of SwI w.r.t. w

      dSwI_dw = dSwI_dw + q * slabs2* (y/w) * (du_dy + yi*du_dx + &
                                  &   OneOvSPi*((1.0_rp-2.0_rp*z*y)/g))

  ! Compute the derivative of SwI w.r.t. n

      dSwI_dn = dSwI_dn + q * slabs2 * y * Log(3.0d2/t) * (du_dy + yi * dv_dy)

  ! Finaly, compute the derivative of SwI w.r.t. Nu0

  ! ***** Analytically *****

  !    dq_dv0 = -(Nu0+dNu)/(Nu0*Nu0)
      dx_dv0 = -x1
      du_dv0 = du_dx * dx_dv0
      dv_dv0 = dv_dx * dx_dv0
  !    db_dv0 = x1
  !    dg_dv0 = 2.0_rp * b*db_dv0
  !    dz_dv0 = (-yi*db_dv0-z*dg_dv0)/g
  !    dr_dv0 = dz_dv0*OneOvSPi+yi*dv_dv0
  !    dvvw_dv0 = (du_dv0+dr_dv0)*q + dq_dv0*(u+r)
  !    dvvw_dv0 = (du_dv0+dr_dv0)*q
  !    dSwI_dNu0 = dslabs1_dNu0*vvw + slabs2*dvvw_dv0
      dSwI_dNu0 = dSwI_dNu0 + swi1 * (dslabs1_dNu0/slabs1 &
                     - 1.0_r8/Nu0) + slabs2*q*(du_dv0 + yi * dv_dv0)

      SwI = SwI + SwI1

    end do

  end subroutine dVoigt_spectral_Lines

  ! ------------------------------------------------------  Slabs  -----
  real(rp) elemental function Slabs ( Nu, v0, v0s, x1, tanh1, slabs1, y )

    use Voigt_m, only: Real_Simple_Voigt

    real(r8), intent(in) :: Nu    ! Frequency
    real(r8), intent(in) :: v0    ! Line center frequency
    real(r8), intent(in) :: v0s   ! Pressure-shifted line center frequency
    real(rp), intent(in) :: x1
    real(rp), intent(in) :: tanh1 ! tanh( h nu / (2 k T) )
    real(rp), intent(in) :: slabs1
    real(rp), intent(in) :: y

! If the molecular transition and temperature have not changed but
! frequency has enter here.

    real(rp) :: u

    call real_simple_voigt ( x1*real(nu-v0s,rp), y, u )

!  Van Vleck - Wieskopf line shape with Voigt, Added Mar/2/91, Bill
!  Modified code to include interference: June/3/1992 (Bill + Zvi)
!  Modified code to correct a sign error (introduced in last change)
!  (Bill + Zvi, July/7/92)

!{ Let $V(a,y)$ be the Voigt function ({\tt u} above), $\delta = \nu-\nu_{0_s}$,
!  $\sigma = \nu + \nu_{0_s}$, $a = x_1 \delta, b = x_1 \sigma$, and
!  $D = \frac1{b^2 + y^2}$.
!  Then {\tt Slabs = } $ S_1 \frac{\nu}{\nu_0}
!   \tanh\left(\frac{h \nu}{2 k T}\right)
!    \left( V(a,y) + V(b,y) \right)$.  $b$ is always huge, so we approximate
!  $V(b,y)$ with one term of an asymptotic expansion, \emph{viz.}
!  $V(b,y) \sim \frac{y D}{\sqrt{\pi}}$.

    Slabs = slabs1 * real(nu / v0, rp) * tanh1 * &
      & (u + OneOvSPi*y/((x1*(nu+v0s))**2 + y*y))

  end function Slabs

  ! ---------------------------------------------  Slabs_dSpectral  -----
  elemental subroutine Slabs_dSpectral ( Nu, Nu0, Nu0s, x1, yi, y, w, t, tanh1, &
                         & VelCor, slabs1, &
                         & dslabs1_dNu0, SwI, dSwI_dw, dSwI_dn, dSwI_dNu0 )

  ! Compute Single Line Absorption and its first derivatives with respect
  ! to spectral parameters w, n & Nu0

  ! NOTE: slabs_prep() must compute Nu0s, x1, yi, y, w, slabs1 and dslabs1_dNu0
  !       in advance

    use Voigt_m, only: Simple_Voigt

    real(r8), intent(in) :: Nu, Nu0  ! Measurement frequency, line center
    real(r8), intent(in) :: Nu0s     ! Pressure shifted line center
    real(rp), intent(in) :: x1, yi, y, w ! Spectral parameters
    real(rp), intent(in) :: t, tanh1 ! temperature, tanh(h*nu/(2*k*T))
    real(rp), intent(in) :: VelCor   ! Doppler velocity correction
    real(rp), intent(in) :: Slabs1   ! frequency-independent absorption terms
    real(rp), intent(in) :: dSlabs1_dNu0 ! d Slabs1 / d Nu0 * 1 / slabs1

    real(rp), intent(out) :: SwI, dSwI_dw, dSwI_dn, dSwI_dNu0 ! Absorption, derivs

    real(rp) :: x !, y  ! Re(z), Im(z)
    real(rp) :: u, v    ! Re(w), Im(w)
    real(rp) :: du, dv  ! Re(dw/dz), Im(dw/dz)
    real(rp) :: d, denom, dH1_dNu0, dx_dNu0, dy_dNu0, g, h1, h2, numer, r
    real(rp) :: s2, sigmaX1, y_g_dHdy, y2

    !{ The absorption is
    !
    ! \begin{equation}\begin{split}
    !  \beta =\,& G H 
    ! \text{ where }
    !  G = S_1 \tanh \left( \frac{h \nu}{2 k T} \right) \frac{\nu}{\nu_0}\,,
    !  H = (H_1 + H_2)\,,\\
    !  H_1 =\,& V(x_1 \delta, y) + y_i L(x_1 \delta, y)\,,
    !  H_2 = V(x_1 \sigma, y) - y_i L(x_1 \sigma, y)\,,\\
    !  \delta =\,& \nu - {\nu_0}_s\,, \sigma = \nu + {\nu_0}_s\,,
    !  w(z) = V(\Re z, \Im z) + i L(\Re z, \Im z) =
    !  \exp(-z^2) \, \text{erfc}(-iz)\,.\\
    !  H_2 \sim\,& \frac{n}{\sqrt\pi d} \text{ where }
    !  n = y - x_1 \text{ and } d = x_1^2 \sigma^2 + y^2\,.
    ! \end{split}\end{equation}
    ! $S_1$ is given by the {\tt slabs1} variable.  The approximation for
    ! $H_2$ arises from a one-term asymptotic approximation for $w(z)$, which
    ! is acceptable because $\sigma$ is known to be large. See {\tt
    ! slabs\_prep} below for definitions of $x_1$ and $y$.

    x = x1 * ( nu - nu0s )
    call simple_voigt ( x, y, u, v, du, dv ) ! get w(z) and dw/dz

    !{ We could use a one-term asymptotic approximation for $L(x_1 \delta, y)$
    !  and $\frac{\partial L}{\partial y}$ because $y_i$ is known to be small,
    !  and therefore we don't need much precision for $L$, but $L$ is needed
    !  for calculation of $\frac{\partial V}{\partial y}$ so using {\tt
    !  Simple\_Voigt} is in fact a saving.

    r = 1.0 / nu0
    sigmaX1 = x1 * ( nu + nu0s )
    s2 = sigmaX1**2
    y2 = y**2
    d = 1.0 / ( s2 + y2 )
    denom = oneOvSpi * d
    numer = y - sigmaX1 * yi
    h1 = u + yi * v
    h2 = numer * denom

    g = slabs1 * tanh1 * nu * r
    SwI = g * ( h1 + h2 )

    !{ The derivatives w.r.t $w$ and $n$ are very similar:
    ! \begin{equation}
    !  \frac{\partial \beta}{\partial w} =
    !  \frac{\partial \beta}{\partial y} \frac{\partial y}{\partial w} =
    !  \frac{y}{w} \frac{\partial \beta}{\partial y} =
    !  \frac{y}{w}\, G\, \frac{\partial H}{\partial y} \text{ and }
    !  \frac{\partial \beta}{\partial n} =
    !  \frac{\partial \beta}{\partial y} \frac{\partial y}{\partial n} =
    !  y\, \frac{300}T\, G\, \frac{\partial H}{\partial y}
    ! \end{equation}
    !%
    ! where
    !%
    ! \begin{equation}
    !  \frac{\partial H}{\partial y} =
    !    \frac{\partial u}{\partial y} + y_i\frac{\partial v}{\partial y} +
    !    \frac1{\sqrt\pi d} \frac{\partial n}{\partial y} -
    !    \frac{n}{\sqrt\pi d^2} \frac{\partial d}{\partial y}=
    !    \frac{\partial u}{\partial y} + y_i\frac{\partial v}{\partial y}
    !    + \frac1d \left(\frac1{\sqrt\pi} - 2 y H_2 \right)\,,
    ! \end{equation}
    ! 
    ! and, by applying the Cauchy-Riemann equations
    ! 
    ! \begin{equation}
    !  \frac{\partial u}{\partial y} = - \frac{\partial v}{\partial x} =
    !  - \Im\, \frac{\partial z}{\partial x} \frac{\text{d} w}{\text{d} z} =
    !  - \Im\, \frac{\text{d} w}{\text{d} z} \text{ and }
    !  \frac{\partial v}{\partial y} = \frac{\partial u}{\partial x} =
    !    \Re\, \frac{\partial z}{\partial x} \frac{\text{d} w}{\text{d} z} =
    !    \Re\, \frac{\text{d} w}{\text{d} z}\,.
    ! \end{equation}

    y_g_dHdy = y * g * ( - dv + yi * du + denom - 2.0 * d * y * h2 )
    dSwI_dw = y_g_dHdy / w
    dSwI_dn = y_g_dHdy * 300.0 / T

    !{\begin{equation}\begin{split}
    !  \frac{\partial \beta}{\partial \nu_0} =\,&
    !   \beta \left( \frac1{S_1} \frac{\partial S_1}{\partial \nu_0} -
    !                \frac1{\nu_0} \right) + G \frac{\partial H}{\partial \nu_0} =
    !  \frac\beta{S_1}\frac{\partial S_1}{\partial \nu_0} - \frac{G H_1}{\nu_0} -
    !  \frac{G H_2}{\nu_0} + G \frac{\partial H_1}{\partial \nu_0} +
    !                        G \frac{\partial H_2}{\partial \nu_0}\,.\\
    !  \frac{\partial H_1}{\partial \nu_0} =\,& \Re \frac{\partial w}{\partial \nu_0} +
    !  y_i \Im \frac{\partial w}{\partial \nu_0}\,.\,\,
    !  \text{ Using the chain rule and } z = x_1 ( \nu - {\nu_0}_s ) + i y\,,\\
    !  \frac{\partial w}{\partial \nu_0} =\,&
    !   \frac{\text{d} w}{\text{d} z} \frac{\partial z}{\partial \nu_0}
    !  = (a + i b ) \left( \frac{\partial x}{\partial \nu_0} +
    !                       i \frac{\partial y}{\partial \nu_0} \right)
    !  = a \frac{\partial x}{\partial \nu_0} - b \frac{\partial y}{\partial \nu_0}
    !  + i \left[ b \frac{\partial x}{\partial \nu_0} + a \frac{\partial y}{\partial \nu_0}
    !      \right] \text{, and}\\
    ! \frac{\partial H_1}{\partial \nu_0}
    !  =\,& a \frac{\partial x}{\partial \nu_0} - b \frac{\partial y}{\partial \nu_0}
    !  + y_i \left[ b \frac{\partial x}{\partial \nu_0} + a \frac{\partial y}{\partial \nu_0}
    !      \right]\,,
    !  \text{ where}\\
    ! a =\,& \Re \frac{\text{d} w}{\text{d} z}\,,
    ! b =    \Im \frac{\text{d} w}{\text{d} z}\,,
    ! \frac{\partial y}{\partial \nu_0} = -\frac{y}{\nu_0}\,,
    ! \frac{\partial x_1}{\partial \nu_0} = -\frac{x_1}{\nu_0}\,,
    ! \frac{\partial {\nu_0}_s}{\partial \nu_0} = v_c \text{ and } \\
    ! \frac{\partial x}{\partial \nu_0}
    !  =\,& \frac{\partial x_1}{\partial \nu_0} (\nu - {\nu_0}_s)
    !    - x_1 \frac{\partial {\nu_0}_s}{\partial \nu_0} =
    !    -\frac{x}{\nu_0} - x_1 v_c\,.
    ! \end{split}\end{equation}

    !{\begin{equation}\begin{split}
    !  \frac{\partial H_2}{\partial \nu_0} =\,&
    !   \frac1{\sqrt\pi d}\frac{\partial n}{\partial \nu_0} -
    !   \frac{n}{\sqrt\pi d^2}\frac{\partial d}{\partial \nu_0}\,.\\
    !  \frac{\partial n}{\partial \nu_0} =\,& \frac{\partial y}{\partial \nu_0} -
    !   \frac{\partial x_1}{\partial \nu_0} \sigma y_i -
    !   x_1 \frac{\partial \sigma}{\partial \nu_0} y_i -
    !   x_1 \sigma \frac{\partial y_i}{\partial \nu_0} =
    !   -\frac{y}{\nu_0} + \frac{x_1}{\nu_0} \sigma y_i
    !    - x_1 \frac{\partial {\nu_0}_s}{\partial \nu_0} y_i - 0 \\
    !   =\,& -\frac{n}{\nu_0} - x_1 v_c y_i\,.\\
    !  \frac{\partial d}{\partial \nu_0} =\,&
    !   2 x_1 \frac{\partial x_1}{\partial \nu_0} \sigma^2
    !   + 2 x_1^2 \sigma \frac{\partial \sigma}{\partial \nu_0}
    !   + 2 y \frac{\partial y}{\partial \nu_0} =
    !  -2\frac{x_1^2 \sigma^2}{\nu_0} - 2 x_1^2 \sigma v_c
    !   -2 \frac{y^2}{\nu_0}\\
    !  =\,& -2\frac{d}{\nu_0} - 2 x_1^2 \sigma v_c\,. \text{ Thus}\\
    !  \frac{\partial H_2}{\partial \nu_0} =\,&
    !   -\frac{n}{\sqrt\pi d \nu_0} -\frac{x_1 v_c y_i}{\sqrt\pi d}
    !   +\frac{n}{\sqrt\pi d^2}\left(\frac{2 d}{\nu_0} + 2 x_1^2 \sigma v_c \right)\\
    !  =\,& \frac{H_2}{\nu_0} - \frac{x_1 v_c y_i}{\sqrt\pi d}
    !   +2 H_2 \frac{x_1^2 \sigma v_c}d\,.\\
    ! \end{split}\end{equation}
    !
    ! Noticing that the $\frac{G H_2}{\nu_0}$ terms arising from
    ! $\frac{\partial G}{\partial \nu_0} H_2$ and
    ! $G \frac{\partial H_2}{\partial \nu_0}$ cancel, we finally have
    !
    ! \begin{equation}
    ! \frac{\partial \beta}{\partial \nu_0} =
    !  \beta \frac1{S_1} \frac{\partial S_1}{\partial \nu_0}
    !   - G \left[\frac{H_1}{\nu_0} - \frac{\partial H_1}{\partial \nu_0}
    !    + \frac{x_1 v_c}d
    !     \left( \frac{y_i}{\sqrt\pi}- 2 H_2 x_1 \sigma \right) \right]\,.
    ! \end{equation}

    dx_dNu0 = - r * x - x1 * velcor
    dy_dNu0 = - r * y
    dH1_dNu0 = du * dx_dNu0 - dv * dy_dNu0 + yi * ( dv * dx_dNu0 + du * dy_dNu0 )
    dSwI_dNu0 = SwI * dslabs1_dNu0 - & ! remember dslabs1_dNu0 is divided by slabs1
      & g * ( r * h1 - dH1_dNu0 + &
      &       d * x1 * velCor * ( yi * oneOvSpi - 2.0 * h2 * sigmaX1 ) )

  end subroutine Slabs_dSpectral

  ! ---------------------------------------  Slabs_lines_dSpectral  -----
  elemental subroutine Slabs_lines_dSpectral ( Nu, Slabs, T, tanh1, VelCor, &
                                & NoPolarized, SwI, dSwI_dw, dSwI_dn, dSwI_dNu0 )

  ! Compute Single Line Absorption and its first derivatives with respect
  ! to spectral parameters w, n & Nu0, without interference.

  ! NOTE: slabs_prep() must compute Nu0s, x1, yi, y, w, slabs1 and dslabs1_dNu0
  !       in advance

    use SpectroscopyCatalog_m, only: Lines
    use Voigt_m, only: Simple_Voigt

    real(r8), intent(in) :: Nu       ! Measurement frequency
    type(slabs_struct), intent(in) :: Slabs ! Frequency-independent stuff, including
      ! Nu0s     ! Pressure shifted line center
      ! x1, yi, y, w ! Spectral parameters
      ! Slabs1   ! frequency-independent absorption terms
      ! dSlabs1_dNu0 ! d Slabs1 / d Nu0 * 1 / slabs1
    real(rp), intent(in) :: T, tanh1 ! temperature, tanh(h*nu/(2*k*T))
    real(rp), intent(in) :: VelCor   ! Doppler velocity correction
    logical, intent(in) :: NoPolarized ! Don't evaluate for polarized lines

    real(rp), intent(out) :: SwI, dSwI_dw, dSwI_dn, dSwI_dNu0 ! Absorption, derivs

    integer :: L         ! Line index
    real(r8) :: Nu0s     ! Pressure shifted line center from Slabs structure
    real(rp) :: x1, y    ! Spectral parameters from Slabs structure
    real(rp) :: x !, y   ! Re(z), Im(z)
    real(rp) :: u, v     ! Re(w), Im(w)
    real(rp) :: du, dv   ! Re(dw/dz), Im(dw/dz)
    real(rp) :: d, denom, dx_dNu0, dy_dNu0, g, h1, h2, r
    real(rp) :: s2, sigmaX1, Swi_up, y_g_dHdy, y2

    ! See Slabs_dSpectral for TeXnicalities

    SwI = 0.0_rp
    dSwI_dw = 0.0_rp
    dSwI_dn = 0.0_rp
    dSwI_dNu0 = 0.0_rp

    do l = 1, size(slabs%s)

      if ( noPolarized .and. slabs%s(l)%polarized ) cycle

      nu0s = slabs%s(l)%v0s
      x1 = slabs%s(l)%x1
      y = slabs%s(l)%y

      x = x1 * ( nu - nu0s )
      call simple_voigt ( x, y, u, v, du, dv ) ! get w(z) and dw/dz

      r = 1.0 / slabs%s(l)%v0
      sigmaX1 = x1 * ( nu + nu0s )
      s2 = sigmaX1**2
      y2 = y**2
      d = 1.0 / ( s2 + y2 )
      denom = oneOvSpi * d
      g = slabs%s(l)%slabs1 * tanh1 * nu * r

      dx_dNu0 = - r * x - x1 * velcor
      dy_dNu0 = - r * y

      h1 = u
      h2 = y * denom

      SwI_up = g * ( h1 + h2 )

      y_g_dHdy = y * g * ( - dv + denom - 2.0 * d * y * h2 )
      dSwI_dw = dSwI_dw + y_g_dHdy / lines(slabs%catalog%lines(l))%w
      dSwI_dn = dSwI_dn + y_g_dHdy * 300.0 / T

      dSwI_dNu0 = dSwI_dNu0 + &
        & SwI_up * slabs%s(l)%dslabs1_dv0 - & ! remember dslabs1_dv0 is divided by slabs1
        & g * ( r * h1 - ( du * dx_dNu0 - dv * dy_dNu0 ) - &
        &       2.0 * d * x1 * velCor * h2 * sigmaX1 )

      SwI = SwI + SwI_up

    end do

  end subroutine Slabs_Lines_dSpectral

  ! ----------------------------------------------------  Slabs_dT  -----
  elemental subroutine Slabs_dT ( Nu, v0, v0s, x1, tanh1, slabs1, y, &
    &                           dv0s_dT, dx1_dT, dtanh_dT, dslabs1_dT, dy_dT, &
    &                   Slabs, dSlabs_dT )

  ! Compute single-line absorption and its derivative w.r.t. temperature.

!   use Voigt_m, only: D_Real_Simple_Voigt
    use Voigt_m, only: Simple_Voigt

    real(r8), intent(in) :: Nu, v0, v0s, dv0s_dT
    real(rp), intent(in) :: x1
    real(rp), intent(in) :: tanh1      ! tanh(h nu / 2 k T), 
    real(rp), intent(in) :: slabs1, y
    real(rp), intent(in) :: dx1_dT     ! 1/x1 dx1 / dT
    real(rp), intent(in) :: dtanh_dT   ! 1/tanh(...) d tanh(h nu / 2 k T) / dT
      ! = h nu / ( 2 k t^2 ) (tanh(...) - 1/tanh(...) )
    real(rp), intent(in) :: dslabs1_dT ! 1/slabs1 dslabs1 / dT
    real(rp), intent(in) :: dy_dT      ! 1/y dy / dT
    real(rp), intent(out) :: Slabs, dSlabs_dT

    real(rp) :: C       ! Terms common to the two parts of dSlabs_dT
    real(rp) :: D       ! 1 / (SigmaX1**2 + y**2)
    real(rp) :: Delta   ! Nu-v0s
    real(rp) :: Du      ! du/dT
    real(rp) :: Da      ! d(x1*delta)
    real(rp) :: Sa, Sb  ! parts of Slabs
    real(rp) :: Sigma   ! Nu+v0s
    real(rp) :: SigmaX1 ! sigma * x1
    real(rp) :: U       ! Voigt
    real(rp) :: Y2      ! y**2

    delta = Nu - v0s
    da = x1 * ( delta * dx1_dT - dv0s_dT ) ! remember, dx1_dT = 1/x1 dx1 / dT
    call simple_voigt ( x1*delta, y, u, du=du, dx=da, dy=y*dy_dT )
!   call D_Real_Simple_Voigt ( x1*delta, y, da, y*dy_dT, u, du )

!{ Let $V(a,y)$ be the Voigt function ({\tt u} above), $\delta = \nu-\nu_{0_s}$,
!  $\sigma = \nu + \nu_{0_s}$, $a = x_1 \delta, b = x_1 \sigma$, and
!  $D = \frac1{b^2 + y^2}$.
!  Then {\tt Slabs = } $ S_1 \frac{\nu}{\nu_0}
!   \tanh\left(\frac{h \nu}{2 k T}\right)
!    \left( V(a,y) + V(b,y) \right)$.
!  This doesn't include the interference terms $y_i ( L(a,y) - L(b,y) )$.
!  $b$ is always huge, so we approximate
!  $V(b,y)$ with one term of an asymptotic expansion, \emph{viz.}
!  $V(b,y) \sim \frac{y D}{\sqrt{\pi}}$.

    sigma = nu + v0s
    sigmaX1 = sigma * x1
    y2 = y*y
    d = 1.0_rp / ( sigmaX1**2 + y2 )
    sa = slabs1 * real(nu / v0,rp) * tanh1
    sb = sa * OneOvSPi * y * d
!   sb = slabs1 * real(nu / v0,rp) * tanh1
!   sa = sb * u
!   sb = sb * OneOvSPi * y * d
    Slabs = sa * u + sb
!   Slabs = sa + sb

!{ The Fadeeva function $w(z)$, where $z = a + i y$, can be written as $V(a,y) +
!  i L(a,y)$.  $V(a,y)$ is frequently called the Voigt function ({\tt u}
!  above).  All we want for $S$ is its real part, so we don't need $L(a,y)$ for
!  $S$.  For $\frac{\partial S}{\partial T}$ we want the real part of the
!  derivative; this requires the real part of the derivative of $w(z)$, not just
!  Voigt.
!
!  Write $S = S_a + S_b$ where
!  $S_a = S_1 \frac{\nu}{\nu_0} \tanh\left(\frac{h \nu}{2 k T}\right) V(a,y)$ and
!  $S_b = S_1 \frac{\nu}{\nu_0} \tanh\left(\frac{h \nu}{2 k T}\right) \frac{y D}{\sqrt{\pi}}$.\\
!  Then
!  $\frac{\partial S}{\partial T} = \frac{\partial S_a}{\partial T} +
!   \frac{\partial S_b}{\partial T}$, where\\
!  $\frac1{S_a}\frac{\partial S_a}{\partial T} =
!   \frac1{S_1}\frac{\partial S_1}{\partial T} +
!   \frac1{\tanh\left(\frac{h \nu}{2 k T}\right)}\frac{\partial}{\partial T} \tanh\left(\frac{h \nu}{2 k T}\right) +
!   \frac1{V(a,y)} \Re \frac{\partial w(z)}{\partial T}$ and\\
!  $\frac1{S_b}\frac{\partial S_b}{\partial T} = 
!   \frac1{S_1}\frac{\partial S_1}{\partial T} +
!   \frac1{\tanh\left(\frac{h \nu}{2 k T}\right)}\frac{\partial}{\partial T} \tanh\left(\frac{h \nu}{2 k T}\right) +
!    \frac1y \frac{\partial y}{\partial T} - 
!    2 D \left( x_1 \sigma \left( x_1 \frac{\partial \nu_{0_s}}{\partial T} +
!     \sigma \frac{\partial x_1}{\partial T} \right) +
!     y \frac{\partial y}{\partial T} \right)$.\\
!    Notice that the first two terms of
!    $\frac1{S_a}\frac{\partial S_a}{\partial T}$ and
!    $\frac1{S_b}\frac{\partial S_b}{\partial T}$ are the same.\\
!  For $\Re \frac{\partial w(z)}{\partial T}$ we need
!  $\frac{\partial a}{\partial T} = -x_1 \frac{\partial \nu_{0_s}}{\partial T} +
!   \delta \frac{\partial x_1}{\partial T}$ and $\frac{\partial y}{\partial T}$.

    c = dSlabs1_dT + dtanh_dT
    dSlabs_dT = sa * ( c * u + du ) + &
!   dSlabs_dT = sa * ( c + du / u ) + &
      &         sb * ( c + dy_dT - 2.0_rp * D * &
      &           ( sigmaX1 * ( x1 * dv0s_dT + sigmaX1 * dx1_dT ) + y2 * dy_dT ) )

  end subroutine Slabs_dT

  ! --------------------------------------------------  Slabs_dAll  -----
! elemental &
  subroutine Slabs_dAll ( Nu, v0, v0s, x1, y, w, T, tanh1, slabs1, &
    &                     dv0s_dT, dx1_dT, dtanh_dT, dslabs1_dT, dy_dT, &
    &                     dslabs1_dNu0, velCor, &
    &                     Slabs, dSlabs_dT, dSlabs_dw, dSlabs_dn, dSlabs_dNu0 )

  ! Compute single-line absorption and its derivatives w.r.t. temperature
  ! and spectroscopy parameters w, n and nu_0.  No interference.

    use Voigt_m, only: Simple_Voigt

    real(r8), intent(in) :: Nu, v0, v0s, dv0s_dT
    real(rp), intent(in) :: x1, y, w   ! Spectral parameters
    real(rp), intent(in) :: T          ! Temperature, K
    real(rp), intent(in) :: tanh1      ! tanh(h nu / 2 k T), 
    real(rp), intent(in) :: slabs1     ! Frequency-independent strength
    real(rp), intent(in) :: dx1_dT     ! 1/x1 dx1 / dT
    real(rp), intent(in) :: dtanh_dT   ! 1/tanh(...) d tanh(h nu / 2 k T) / dT
      ! = h nu / ( 2 k t^2 ) (tanh(...) - 1/tanh(...) )
    real(rp), intent(in) :: dslabs1_dT ! 1/slabs1 dslabs1 / dT
    real(rp), intent(in) :: dy_dT      ! 1/y dy / dT
    real(rp), intent(in) :: dslabs1_dNu0 ! 1/slabs1 dslabs1 / dNu0
    real(rp), intent(in) :: velCor     ! Doppler correction = d v0s / d v0
    real(rp), intent(out) :: Slabs, dSlabs_dT, dSlabs_dw, dSlabs_dn, dSlabs_dNu0

    real(rp) :: A, B    ! Re(dw/dz), Im(dw/dz)
    real(rp) :: C       ! Terms common to the two parts of dSlabs_dT
    real(rp) :: D       ! 1 / (SigmaX1**2 + y**2)
    real(rp) :: Delta   ! Nu - v0s
    real(rp) :: Denom   ! D / Sqrt(pi)
    real(rp) :: Du      ! du/dT 
    real(rp) :: Dx      ! d(x1*delta)
    real(rp) :: Dx_dNu0 ! dx / dNu0
    real(rp) :: Dy_dNu0 ! dy / dNu0
    real(rp) :: G       ! Slabs1 * tanh1 * nu / v0
    real(rp) :: H2      ! y * Denom
    real(rp) :: Q       ! nu/v0
    real(rp) :: R       ! 1/v0
    real(rp) :: Sa, Sb  ! parts of Slabs
    real(rp) :: SigmaX1 ! ( Nu + v0s ) * x1
    real(rp) :: S2      ! ( sigma * x1 ) ** 2
    real(rp) :: U, V    ! Voigt, Lorentz
    real(rp) :: X       ! X1 * Delta
    real(rp) :: Y2      ! y**2
    real(rp) :: Y_G_dHdy ! Y * G * (du/dy + dH2/dy)

! See Slabs_dT and Slabs_dSpectral for TeXnicalities.

    ! Absorption

    delta = Nu - v0s
    x = x1 * delta
    dx = x * dx1_dT - x1 * dv0s_dT ! remember, dx1_dT = 1/x1 dx1 / dT
    call simple_voigt ( x, y, u, v, a, b )
    du = a * dx - b * y * dy_dT

    r = 1 / v0
    q = nu * r
    sigmaX1 = (nu + v0s) * x1
    s2 = sigmaX1**2
    y2 = y*y
    d = 1.0_rp / ( s2 + y2 )
    denom = OneOvSPi * d
    g = slabs1 * q * tanh1
    h2 = y * denom
    sa = g * u
    sb = g * h2
    Slabs = sa + sb

    ! Temperature derivative

    c = dSlabs1_dT + dtanh_dT
    dSlabs_dT = g * ( c * u + du ) + &
!   dSlabs_dT = sa * ( c + du / u ) + &
      &         sb * ( c + dy_dT - 2.0_rp * d * &
      &           ( sigmaX1 * ( x1 * dv0s_dT + sigmaX1 * dx1_dT ) + y2 * dy_dT ) )

    ! Spectroscopy derivatives

    y_g_dHdy = y * g * ( - b + denom - 2.0 * d * y * h2 )
    dSlabs_dw = y_g_dHdy / w
    dSlabs_dn = y_g_dHdy * 300.0 / T

    dx_dNu0 = - r * x - x1 * velcor
    dy_dNu0 = - r * y
    dSlabs_dNu0 = Slabs * dslabs1_dNu0 - & ! remember dslabs1_dNu0 is divided by slabs1
      & g * ( r * u - a * dx_dNu0 + b * dy_dNu0 - &
      &       2.0 * d * x1 * velCor * h2 * sigmaX1 )

  end subroutine Slabs_dAll

  ! ------------------------------------------------  Slabs_Lines  -----
  elemental function Slabs_Lines ( Nu, Slabs, tanh1, NoPolarized ) result ( Beta )

    use Voigt_m, only: Real_Simple_Voigt

    real(r8), intent(in) :: Nu    ! Frequency
    type(slabs_struct), intent(in) :: Slabs ! Frequency-independent stuff
    real(rp), intent(in) :: tanh1 ! tanh( h nu / (2 k T) )
    ! "Don't do line(L) if slabs%catalog%polarized(L)"
    logical, intent(in) :: NoPolarized

    real(rp) :: Beta ! Output

! If the molecular transition and temperature have not changed but
! frequency has, enter here to calculate sum of Beta for all lines in SLABS.

    integer :: L      ! Line Index
    real(rp) :: U     ! Voigt = real part of Fadeeva
    real(r8) :: V0S   ! Pressure-shifted line center
    real(r8) :: X1, Y ! Doppler width, ratio Pressure to Doppler widths

!{ Let $V(a,y)$ be the Voigt function ({\tt u} above), $\delta = \nu-\nu_{0_s}$,
!  $\sigma = \nu + \nu_{0_s}$, $a = x_1 \delta$, and
!  $D = \frac1{\sigma^2 x_1^2 + y^2}$.
!  Then {\tt Slabs = } $ S_1 \frac{\nu}{\nu_0}
!   \tanh\left(\frac{h \nu}{2 k T}\right)
!    \left( V(a,y) + \frac{y D}{\sqrt{\pi}} \right)$.

    beta = 0.0_rp
    if ( .not. noPolarized ) then
      do l = 1, size(slabs%s)
        v0s = slabs%s(l)%v0s
        x1 = slabs%s(l)%x1
        y = slabs%s(l)%y
        call real_simple_voigt ( x1*real(nu-v0s,rp), y, u )

        beta = beta + slabs%s(l)%slabs1 * &
          &           real(nu / slabs%s(l)%v0, rp) * tanh1 * &
          & (u + OneOvSPi*y/((x1*(nu+v0s))**2 + y*y))

      end do
    else
      do l = 1, size(slabs%s)
        if ( slabs%s(l)%polarized ) cycle
        v0s = slabs%s(l)%v0s
        x1 = slabs%s(l)%x1
        y = slabs%s(l)%y
        call real_simple_voigt ( x1*real(nu-v0s,rp), y, u )

        beta = beta + slabs%s(l)%slabs1 * &
          &           real(nu / slabs%s(l)%v0, rp) * tanh1 * &
          & (u + OneOvSPi*y/((x1*(nu+v0s))**2 + y*y))

      end do
    end if

  end function Slabs_Lines

  ! ---------------------------------------------  Slabs_Lines_dT  -----
!  elemental &
  pure &
  subroutine Slabs_Lines_dT ( Nu, Slabs, Tanh1, dTanh_dT, &
    &                         Beta, dBeta_dT, NoPolarized )

  ! Compute single-line absorption and its derivative w.r.t. temperature
  ! for all lines in the Slabs structure.

    use Voigt_m, only: D_Real_Simple_Voigt
!   use Voigt_m, only: Simple_Voigt

    real(r8), intent(in) :: Nu    ! Frequency
    type(slabs_struct), intent(in) :: Slabs ! Frequency-independent stuff
    real(rp), intent(in) :: Tanh1 ! tanh( h nu / (2 k T) )
    real(rp), intent(in) :: dTanh_dT ! -h nu / (2 k T^2) 1/tanh(...) dTanh(...)/dT

    real(rp), intent(out) :: Beta, dBeta_dT

    ! "Don't do line(L) if slabs%catalog%polarized(L)"
    logical, intent(in) :: NoPolarized

    real(rp) :: C       ! Terms common to the two parts of dSlabs_dT
    real(rp) :: D       ! 1 / (SigmaX1**2 + y**2)
    real(rp) :: Delta   ! Nu-v0s
    real(rp) :: Du      ! du/dT
    real(rp) :: Da      ! d(x1*delta)
    real(rp) :: Dv0s_dT ! dv0s / dT
    real(rp) :: Dx1_dT  ! 1/x1 dx1/dT
    real(rp) :: Dy_dT   ! 1/y dy/dT
    integer :: L        ! Line index
    real(rp) :: Sa, Sb  ! parts of Slabs
    real(rp) :: SigmaX1 ! (nu + v0s) * x1
    real(rp) :: U       ! Voigt
    real(r8) :: V0S     ! Pressure-shifted line center
    real(rp) :: X1, Y   ! Doppler width, ratio Pressure to Doppler widths
    real(rp) :: Y2      ! y**2

! See Slabs_dT for TeXnicalities

    Beta = 0.0_rp
    dBeta_dT = 0.0_rp

    do l = 1, size(slabs%s) ! == size(slabs%d) == size(slabs%catalog%lines)

      if ( noPolarized .and. slabs%s(l)%polarized ) cycle

      v0s = slabs%s(l)%v0s
      x1 = slabs%s(l)%x1
      y = slabs%s(l)%y
      sa = slabs%s(l)%slabs1 * (nu / slabs%s(l)%v0) * tanh1
      dv0s_dT = slabs%d(l)%dv0s_dT
      dx1_dT = slabs%d(l)%dx1_dT
      dy_dT = slabs%d(l)%dy_dT
      c = slabs%d(l)%dSlabs1_dT + dtanh_dT
      delta = Nu - v0s
      da = x1 * ( delta * dx1_dT - dv0s_dT )
      call D_Real_Simple_Voigt ( x1*delta, y, da, y*dy_dT, u, du )
      sigmaX1 = (nu + v0s) * x1
      y2 = y * y
      d = 1.0_rp / ( sigmaX1**2 + y2 )
      sb = sa * OneOvSPi * y * d
      beta = beta + sa * u + sb

      dBeta_dT = dBeta_dT + sa * ( u * c + du ) &
        &                 + sb * ( c + dy_dT - 2.0_rp * d * &            
        &                          ( sigmaX1 * ( x1 * dv0s_dT + &        
        &                            sigmaX1 * dx1_dT ) + y2 * dy_dT ) )

    end do

  end subroutine Slabs_Lines_dT

  ! ----------------------------------------  Slabs_Lines_dT_path  -----
!  elemental &
  pure &
  subroutine Slabs_Lines_dT_path ( Nu, Slabs, Path_Inds, Tanh1, dTanh_dT, &
    &                              Ratio, Beta, dBeta_dT, NoPolarized )

  ! Compute single-line absorption and its derivative w.r.t. temperature
  ! for all lines in the Slabs structure.

    use Voigt_m, only: D_Real_Simple_Voigt
!   use Voigt_m, only: Simple_Voigt

    real(r8), intent(in) :: Nu          ! Frequency
    type(slabs_struct), intent(in) :: Slabs(:) ! Frequency-independent stuff
    integer, intent(in) :: Path_Inds(:) ! Indices for Slabs
    real(rp), intent(in) :: Tanh1(:)    ! tanh( h nu / (2 k T) )
    real(rp), intent(in) :: dTanh_dT(:) ! -h nu / (2 k T^2) 1/tanh(...) dTanh(...)/dT
    real(rp), intent(in) :: Ratio       ! Isotope ratio

    real(rp), intent(inout) :: Beta(:), dBeta_dT(:)

    ! "Don't do line(L) if slabs%catalog%polarized(L)"
    logical, intent(in) :: NoPolarized
    real :: B, dB, th, dTh
    real(rp) :: C       ! Terms common to the two parts of dSlabs_dT
    real(rp) :: D       ! 1 / (SigmaX1**2 + y**2)
    real(rp) :: Delta   ! Nu-v0s
    real(rp) :: Du      ! du/dT
    real(rp) :: Da      ! d(x1*delta)
    real(rp) :: Dv0s_dT ! dv0s / dT
    real(rp) :: Dx1_dT  ! 1/x1 dx1/dT
    real(rp) :: Dy_dT   ! 1/y dy/dT
    integer :: J, K     ! Path indices
    integer :: L        ! Line index
    real(rp) :: Sa, Sb  ! parts of Slabs
    real(rp) :: SigmaX1 ! (nu + v0s) * x1
    real(rp) :: U       ! Voigt
    real(r8) :: V0S     ! Pressure-shifted line center
    real(rp) :: X1, Y   ! Doppler width, ratio Pressure to Doppler widths
    real(rp) :: Y2      ! y**2

! See Slabs_dT for TeXnicalities
    do j = 1, size(path_inds)
      k = path_inds(j)
      b = 0.0  ! Beta for one path slot
      db = 0.0 ! dBeta for one path slot
      th = tanh1(j)
      dTh = dTanh_dT(j)
      do l = 1, size(slabs(k)%s) ! == size(slabs(k)%d) == size(slabs(k)%catalog%lines)

        if ( noPolarized .and. slabs(k)%s(l)%polarized ) cycle

        v0s = slabs(k)%s(l)%v0s
        x1 = slabs(k)%s(l)%x1
        y = slabs(k)%s(l)%y
        sa = slabs(k)%s(l)%slabs1 * (nu / slabs(k)%s(l)%v0) * th
        dv0s_dT = slabs(k)%d(l)%dv0s_dT
        dx1_dT = slabs(k)%d(l)%dx1_dT
        dy_dT = slabs(k)%d(l)%dy_dT
        c = slabs(k)%d(l)%dSlabs1_dT + dTh
        delta = Nu - v0s
        da = x1 * ( delta * dx1_dT - dv0s_dT )
        call D_Real_Simple_Voigt ( x1*delta, y, da, y*dy_dT, u, du )
        sigmaX1 = (nu + v0s) * x1
        y2 = y * y
        d = 1.0_rp / ( sigmaX1**2 + y2 )
        sb = sa * OneOvSPi * y * d
        b = b + sa * u + sb

        db = db + sa * ( u * c + du ) &
          &                 + sb * ( c + dy_dT - 2.0_rp * d * &            
          &                          ( sigmaX1 * ( x1 * dv0s_dT + &        
          &                            sigmaX1 * dx1_dT ) + y2 * dy_dT ) )

      end do
      beta(j) = beta(j) + ratio*b
      dBeta_dT(j) = dBeta_dT(j) + ratio*db
    end do

  end subroutine Slabs_Lines_dT_path

  ! ---------------------------------------------  Slabs_Lines_dAll  -----
  elemental &
  subroutine Slabs_Lines_dAll ( Nu, Slabs, T, tanh1, &
    &                     dtanh_dT, velCor, NoPolarized, &
    &                     Beta, dBeta_dT, dBeta_dw, dBeta_dn, dBeta_dNu0 )

  ! Compute single-line absorption and its derivatives w.r.t. temperature
  ! and spectroscopy parameters w, n and nu_0.  No interference.

    use SpectroscopyCatalog_m, only: Lines
    use Voigt_m, only: Simple_Voigt

    real(r8), intent(in) :: Nu         ! Measurement frequency
    type(slabs_struct), intent(in) :: Slabs ! Frequency-independent stuff
    real(rp), intent(in) :: T          ! Temperature, K
    real(rp), intent(in) :: tanh1      ! tanh(h nu / 2 k T), 
    real(rp), intent(in) :: dtanh_dT   ! 1/tanh(...) d tanh(h nu / 2 k T) / dT
      ! = h nu / ( 2 k t^2 ) (tanh(...) - 1/tanh(...) )
    real(rp), intent(in) :: velCor     ! Doppler correction = d v0s / d v0
    logical, intent(in) :: NoPolarized ! Don't evaluate for polarized lines
    real(rp), intent(out) :: Beta, dBeta_dT, dBeta_dw, dBeta_dn, dBeta_dNu0

    ! Stuff from Slabs structure
    real(rp) :: dslabs1_dNu0 ! 1/slabs1 dslabs1 / dNu0
    real(rp) :: dv0s_dT    ! dv0s / dT ( not 1/V0s dv0s / dT !)
    real(rp) :: dx1_dT     ! 1/x1 dx1 / dT
    real(rp) :: dy_dT      ! 1/y dy / dT
    real(r8) :: v0s        ! Pressure-shifted line center
    real(rp) :: x1, y, w   ! Spectral parameters

    integer :: L        ! Line index
    real(rp) :: A, B    ! Re(dw/dz), Im(dw/dz)
    real(rp) :: Beta_up ! Beta update
    real(rp) :: C       ! Terms common to the two parts of dSlabs_dT
    real(rp) :: D       ! 1 / (SigmaX1**2 + y**2)
    real(rp) :: Delta   ! Nu - v0s
    real(rp) :: Denom   ! D / Sqrt(pi)
    real(rp) :: Du      ! du/dT 
    real(rp) :: Dx      ! d(x1*delta)
    real(rp) :: Dx_dNu0 ! dx / dNu0
    real(rp) :: Dy_dNu0 ! dy / dNu0
    real(rp) :: G       ! Slabs1 * tanh1 * nu / v0
    real(rp) :: H2      ! y * Denom
    real(rp) :: Q       ! nu/v0
    real(rp) :: R       ! 1/v0
    real(rp) :: Sa, Sd  ! parts of Slabs
    real(rp) :: SigmaX1 ! ( Nu + v0s ) * x1
    real(rp) :: S2      ! ( sigma * x1 ) ** 2
    real(rp) :: U, V    ! Voigt, Lorentz
    real(rp) :: X       ! X1 * Delta
    real(rp) :: Y2      ! y**2
    real(rp) :: Y_dYdT  ! y * dy/dT
    real(rp) :: Y_G_dHdy ! Y * G * (du/dy + dH2/dy)

! See Slabs_dT and Slabs_dSpectral for TeXnicalities.

    Beta = 0.0
    dBeta_dT = 0.0
    dBeta_dw = 0.0
    dBeta_dn = 0.0
    dBeta_dNu0 = 0.0

    do l = 1, size(slabs%s)

      if ( noPolarized .and. slabs%s(l)%polarized ) cycle

      dslabs1_dNu0 = slabs%s(l)%dslabs1_dv0
      dv0s_dT = slabs%d(l)%dv0s_dT
      dx1_dT = slabs%d(l)%dx1_dT
      dy_dT = slabs%d(l)%dy_dT
      v0s = slabs%s(l)%v0s
      w = lines(slabs%catalog%lines(l))%w
      x1 = slabs%s(l)%x1
      y = slabs%s(l)%y

      ! Absorption

      delta = Nu - v0s
      x = x1 * delta
      dx = x * dx1_dT - x1 * dv0s_dT ! remember, dx1_dT = 1/x1 dx1 / dT
      call simple_voigt ( x, y, u, v, a, b )
      y_dydT = y * dy_dT
      du = a * dx - b * y_dydT

      r = 1.0 / slabs%s(l)%v0
      q = nu * r
      sigmaX1 = ( nu + v0s ) * x1
      s2 = sigmaX1**2
      y2 = y*y
      d = 1.0_rp / ( s2 + y2 )
      denom = OneOvSPi * d
      g = slabs%s(l)%slabs1 * q * tanh1
      c = slabs%d(l)%dSlabs1_dT + dtanh_dT
      h2 = y * denom
      sa = g * u
      sd = g * h2

      beta_up = sa + sd

      ! Temperature derivative

      dbeta_dT = dbeta_dT &
!       &      + sa * ( c + du / u ) &
        &      + g * ( u * c + du ) &
        &      + sd * ( c + dy_dT &
        &      - 2.0_rp * d * ( sigmaX1 * ( x1 * dv0s_dT + sigmaX1 * dx1_dT ) + y2 * dy_dT )  )

      ! Spectroscopy derivatives

      y_g_dHdy = y * g * ( - b + denom - 2.0 * d * y * h2 )
      dBeta_dw = dBeta_dw + y_g_dHdy / w
      dBeta_dn = dBeta_dn + y_g_dHdy * 300.0 / T

      dx_dNu0 = - r * x - x1 * velcor
      dy_dNu0 = - r * y
      dBeta_dNu0 = dBeta_dNu0 + &
        & beta_up * dslabs1_dNu0 - & ! remember dslabs1_dNu0 is divided by slabs1
        &   g * ( r * u - a * dx_dNu0 + b * dy_dNu0 - &
        &         2.0 * d * x1 * velCor * h2 * sigmaX1 )

      beta = beta + beta_up

    end do

  end subroutine Slabs_Lines_dAll

  ! --------------------------------------------------  Slabswint  -----
  real(rp) elemental function Slabswint ( Nu, v0, v0s, x1, tanh1, slabs1, y, yi )

    use Voigt_m, only: Real_Simple_Voigt

    real(r8), intent(in) :: Nu    ! Frequency
    real(r8), intent(in) :: v0    ! Line center frequency
    real(r8), intent(in) :: v0s   ! Pressure-shifted line center frequency
    real(rp), intent(in) :: x1
    real(rp), intent(in) :: tanh1 ! tanh( h nu / (2 k T) )
    real(rp), intent(in) :: slabs1
    real(rp), intent(in) :: y
    real(rp), intent(in) :: yi

! If the molecular transition and temperature have not changed but
! frequency has enter here.

    real(rp) :: a, sigmaX1, u, y2

    a = x1 * real(nu-v0s,rp)
    call real_simple_voigt ( a, y, u )

!  Van Vleck - Wieskopf line shape with Voigt, Added Mar/2/91, Bill
!  Modified code to include interference: June/3/1992 (Bill + Zvi)
!  Modified code to correct a sign error (introduced in last change)
!  (Bill + Zvi, July/7/92)

!{ Let $V(a,y)$ be the Voigt function ({\tt u} above), $\delta = \nu-\nu_{0_s}$,
!  $\sigma = \nu + \nu_{0_s}$, $a = x_1 \delta$, $b = x_1 \sigma$,
!  $D_1 = \frac1{\sqrt\pi} \frac1{a^2 + y^2}$ and
!  $D_2 = \frac1{\sqrt\pi} \frac1{b^2 + y^2}$.
!  Then
!   {\tt Slabswint = } $ S_1 \frac{\nu}{\nu_0}
!   \tanh\left(\frac{h \nu}{2 k T}\right)
!    \left( V(a,y) + y_i a D_1
!     + (y - b y_i) D_2 \right)$.
!  The term $y_i a D_1$ is a one-term asymptotic approximation to
!  $L(a, y)$, where $L$ is the imaginary part of Fadeeva.  This is
!  acceptable for all $a$ and $y$ because $y_i$ is known (so far) to
!  be small relative to 1.
!  The term $y D_2$ is a one-term asymptotic approximation to $V(b,
!  y)$, while the term $b D_2$ is a is a one-term asymptotic
!  approximation to $L(b, y)$.  These are both acceptable for all
!  arguments because we know that $b$ is large.

    sigmaX1 = x1 * (nu + v0s)
    y2 = y*y
    Slabswint = slabs1 * real(Nu / v0, rp) * tanh1 * &
      & (u + OneOvSPi*((y - sigmaX1*yi)/(sigmaX1*sigmaX1 + y2) + yi*a/(a*a+y2)))

  end function Slabswint

  ! -----------------------------------------------  Slabswint_dT  -----
  elemental subroutine Slabswint_dT ( Nu, v0, v0s, x1, tanh1, slabs1, y, yi, &
    &                              dv0s_dT, dx1_dT, dtanh_dT, dslabs1_dT, &
    &                              dy_dT, dyi_dT, &
    &                              Slabswint, dSlabs_dT )

    use Voigt_m, only: D_Simple_Voigt

    real(r8), intent(in) :: Nu, v0, v0s, dv0s_dT
    real(rp), intent(in) :: x1
    real(rp), intent(in) :: tanh1      ! tanh(h nu / 2 k T), 
    real(rp), intent(in) :: slabs1, y, yi
    real(rp), intent(in) :: dx1_dT     ! 1/x1 dx1 / dT
    real(rp), intent(in) :: dtanh_dT   ! 1/tanh(...) d tanh(h nu / 2 k T) / dT
      ! = h nu / ( 2 k t^2 ) (tanh(...) - 1/tanh(...) )
    real(rp), intent(in) :: dslabs1_dT ! 1/slabs1 dslabs1 / dT
    real(rp), intent(in) :: dy_dT      ! 1/y dy / dT
    real(rp), intent(in) :: dyi_dT     ! 1/yi dyi / dT

    real(rp), intent(out) :: Slabswint, dSlabs_dT

! If the molecular transition and temperature have not changed but
! frequency has enter here.

    real(rp) :: A               ! x1 * delta
    real(rp) :: C               ! Terms common to the parts of dSlabs_dT
    real(rp) :: D               ! 1 / (a**2 + y**2)
    real(rp) :: DD              ! 1/D dD/dT, 
    real(rp) :: Da              ! da / dT
    real(rp) :: Delta           ! Nu-v0s
    real(rp) :: DU, DV          ! du/dT, dv/dT
    real(rp) :: Sa, Sb, Sc, Sd  ! parts of SlabsWint
    real(rp) :: Sigma           ! Nu+v0s
    real(rp) :: SigmaX1         ! sigma * x1
    real(rp) :: U, V            ! w(z) = u + i v
    real(rp) :: Y2              ! Y**2

    delta = nu - v0s
    a = x1 * delta
    da = x1 * ( delta * dx1_dT - dv0s_dT ) ! remember, dx1_dT = 1/x1 dx1 / dT
    call D_Simple_Voigt ( a, y, u, v, du, dv, da, y*dy_dT )

!  Van Vleck - Wieskopf line shape with Voigt, Added Mar/2/91, Bill
!  Modified code to include interference: June/3/1992 (Bill + Zvi)
!  Modified code to correct a sign error (introduced in last change)
!  (Bill + Zvi, July/7/92)

!{ Let $V(a,y)$ be the Voigt function ({\tt u} above),
!  $L(a,y)$ be the imaginary part of the Fadeeva function ({\tt v} above),
!  $\delta = \nu-\nu_{0_s}$, $\sigma = \nu + \nu_{0_s}$, $a = x_1 \delta$,
!  and $D = \frac1{\sigma^2 x_1^2 + y^2}$.
!  Then {\tt Slabswint = } $ S_1 \frac{\nu}{\nu_0}
!  \tanh\left(\frac{h \nu}{2 k T}\right)
!   \left( V(a,y) + y_i L(a,y)
!     + \frac{(y - \sigma x_1 y_i) D}{\sqrt{\pi}} \right)$.
!  The last term arises from a one-term asymptotic expansion of
!  $V(x_1\sigma,y) - y_i L(x_1\sigma,y)$.  This is acceptable because
!  $\sigma$ is known to be large.  We usually use a one-term asymptotic
!  expansion of $L(a,y)$ as well, because $y_i$ is known to be small.  In
!  this case, however, $L(a,y)$ is needed to compute
!  $\frac{\partial V(a,y)}{\partial T}$, so using the general method for
!  $L(a,y)$ actually represents a saving.
!
!  Write $S = S_a + S_b - S_c + S_d$ where
!  $S_a = S_1 \frac{\nu}{\nu_0} \tanh\left(\frac{h \nu}{2 k T}\right) V(a,y)$,
!  $S_b = S_1 \frac{\nu}{\nu_0} \tanh\left(\frac{h \nu}{2 k T}\right) y_i L(a,y)$,
!  $S_c = S_1 \frac{\nu}{\nu_0} \tanh\left(\frac{h \nu}{2 k T}\right) \frac{\sigma x_1 y_i D}{\sqrt{\pi}}$, and
!  $S_d = S_1 \frac{\nu}{\nu_0} \tanh\left(\frac{h \nu}{2 k T}\right) \frac{y D}{\sqrt{\pi}}$.\\

    sigma = nu + v0s
    sigmaX1 = sigma * x1
    y2 = y * y
    d = 1.0_rp / ( sigmaX1**2 + y2 )
    c = slabs1 * real(nu / v0,rp) * tanh1
    sa = c * u
!   sb = c * yi * v
    sb = c * yi
    c = c * d * OneOvSPi
    sc = c * sigmaX1 * yi
    sd = c * y
    Slabswint = sa + sb * v - sc + sd

!{ Then
!  $\frac{\partial S}{\partial T} = \frac{\partial S_a}{\partial T} +
!   \frac{\partial S_b}{\partial T} - \frac{\partial S_c}{\partial T} +
!   \frac{\partial S_d}{\partial T}$, where\\
!  $\frac1{S_a}\frac{\partial S_a}{\partial T} =
!   \frac1{S_1}\frac{\partial S_1}{\partial T} +
!   \frac1{\tanh\left(\frac{h \nu}{2 k T}\right)}\frac{\partial}{\partial T} \tanh\left(\frac{h \nu}{2 k T}\right) +
!   \frac1{V(a,y)} \Re \frac{\partial w(z)}{\partial T}$,\\
!  $\frac1{S_b}\frac{\partial S_b}{\partial T} = 
!   \frac1{S_1}\frac{\partial S_1}{\partial T} +
!   \frac1{\tanh\left(\frac{h \nu}{2 k T}\right)}\frac{\partial}{\partial T} \tanh\left(\frac{h \nu}{2 k T}\right) +
!   \frac1{L(a,y)} \Im \frac{\partial w(z)}{\partial T}$,\\
!  $\frac1{S_c}\frac{\partial S_c}{\partial T} = 
!   \frac1{S_1}\frac{\partial S_1}{\partial T} +
!   \frac1{\tanh\left(\frac{h \nu}{2 k T}\right)}\frac{\partial}{\partial T} \tanh\left(\frac{h \nu}{2 k T}\right) +
!   \frac1{D}\frac{\partial D}{\partial T} +
!   \frac1{y_i}\frac{\partial y_i}{\partial T} +
!   \frac1{x_1}\frac{\partial x_1}{\partial T} +
!   \frac1{\sigma}\frac{\partial \nu_{0_s}}{\partial T}$, and\\
!  $\frac1{S_d}\frac{\partial S_d}{\partial T} = 
!   \frac1{S_1}\frac{\partial S_1}{\partial T} +
!   \frac1{\tanh\left(\frac{h \nu}{2 k T}\right)}\frac{\partial}{\partial T} \tanh\left(\frac{h \nu}{2 k T}\right) +
!   \frac1{D}\frac{\partial D}{\partial T} +
!   \frac1y \frac{\partial y}{\partial T}$, where\\
!  $\frac1{D}\frac{\partial D}{\partial T} = -2 D \left( x_1 \sigma \left( x_1 \frac{\partial \nu_{0_s}}{\partial T} +
!     \sigma \frac{\partial x_1}{\partial T} \right) +
!     y \frac{\partial y}{\partial T} \right)$.\\

!{ Notice that the first two terms of
!    $\frac1{S_a}\frac{\partial S_a}{\partial T}$,
!    $\frac1{S_b}\frac{\partial S_b}{\partial T}$,
!    $\frac1{S_c}\frac{\partial S_c}{\partial T}$ and
!    $\frac1{S_d}\frac{\partial S_d}{\partial T}$ are the same, and
!  the third terms of
!    $\frac1{S_c}\frac{\partial S_c}{\partial T}$ and
!    $\frac1{S_d}\frac{\partial S_d}{\partial T}$ are the same.\\
!  For $\frac{\partial w(z)}{\partial T}$ we need
!  $\frac{\partial a}{\partial T} = -x_1 \frac{\partial \nu_{0_s}}{\partial T} +
!   \delta \frac{\partial x_1}{\partial T}$ and $\frac{\partial y}{\partial T}$.

!{ We don't actually calculate things this way, because if $\delta \equiv 0$
!  then $L(a,y) \equiv 0$.

    c = dSlabs1_dT + dtanh_dT
    dd = -2.0_rp * d * ( sigmaX1 * ( x1 * dv0S_dT + sigmaX1 * dx1_dT ) + y2 * dy_dT )
    dSlabS_dT = c  * ( c * u + du ) & ! sa * ( c + du / u )
      &       + sb * ( c * v + dv ) & ! sb * v * ( c + dv / v )
      &       - sc * ( c + dd + dyi_dT + dx1_dT + dv0S_dT / sigma ) &
      &       + sd * ( c + dd + dy_dT )

  end subroutine Slabswint_dT

  ! ----------------------------------------------  Slabswint_dAll  -----
  elemental &
  subroutine Slabswint_dAll ( Nu, v0, v0s, x1, yi, y, w, T, tanh1, slabs1, &
    &                         dv0s_dT, dx1_dT, dtanh_dT, dslabs1_dT, dy_dT, dyi_dT, &
    &                         dslabs1_dNu0, velCor, &
    &              Slabswint, dSlabs_dT, dSlabs_dw, dSlabs_dn, dSlabs_dNu0 )

    use Voigt_m, only: Simple_Voigt

    real(r8), intent(in) :: Nu, v0, v0s, dv0s_dT
    real(rp), intent(in) :: x1, yi, y, w ! Spectral parameters
    real(rp), intent(in) :: T          ! Temperature, K
    real(rp), intent(in) :: tanh1      ! tanh(h nu / 2 k T), 
    real(rp), intent(in) :: slabs1     ! Frequency-independent strength
    real(rp), intent(in) :: dx1_dT     ! 1/x1 dx1 / dT
    real(rp), intent(in) :: dtanh_dT   ! 1/tanh(...) d tanh(h nu / 2 k T) / dT
      ! = h nu / ( 2 k t^2 ) (tanh(...) - 1/tanh(...) )
    real(rp), intent(in) :: dslabs1_dT ! 1/slabs1 dslabs1 / dT
    real(rp), intent(in) :: dy_dT      ! 1/y dy / dT
    real(rp), intent(in) :: dyi_dT     ! 1/yi dyi / dT
    real(rp), intent(in) :: dslabs1_dNu0 ! 1/slabs1 dslabs1 / dNu0
    real(rp), intent(in) :: velCor     ! Doppler correction = d v0s / d v0
    real(rp), intent(out) :: Slabswint, dSlabs_dT, dSlabs_dw, dSlabs_dn, dSlabs_dNu0

    real(rp) :: A, B    ! Re(dw/dz), Im(dw/dz)
    real(rp) :: C       ! Terms common to the two parts of dSlabs_dT
    real(rp) :: D       ! 1 / (SigmaX1**2 + y**2)
    real(rp) :: DD      ! 1/D dD/dT, 
    real(rp) :: Delta   ! Nu - v0s
    real(rp) :: Denom   ! D / Sqrt(pi)
    real(rp) :: DH1_dNu0 ! dH1 / dNu0
    real(rp) :: Du, Dv  ! du/dT, du/dT
    real(rp) :: Dx      ! d(x1*delta)
    real(rp) :: Dx_dNu0 ! dx / dNu0
    real(rp) :: Dy_dNu0 ! dy / dNu0
    real(rp) :: G       ! Slabs1 * tanh1 * nu / v0
    real(rp) :: H2      ! y * Denom
    real(rp) :: Q       ! nu/v0
    real(rp) :: R       ! 1/v0
    real(rp) :: Sa, Sb, Sc, Sd  ! parts of Slabswint
    real(rp) :: Sigma   ! Nu + v0s
    real(rp) :: SigmaX1 ! Sigma * x1
    real(rp) :: S2      ! ( sigma * x1 ) ** 2
    real(rp) :: U, V    ! Voigt, Lorentz
    real(rp) :: X       ! X1 * Delta
    real(rp) :: Y2      ! y**2
    real(rp) :: Y_dYdT  ! y * dy/dT
    real(rp) :: Y_G_dHdy ! Y * G * (du/dy + dH2/dy)

! See Slabs_dT and Slabs_dSpectral for TeXnicalities.

    ! Absorption

    delta = Nu - v0s
    x = x1 * delta
    dx = x * dx1_dT - x1 * dv0s_dT ! remember, dx1_dT = 1/x1 dx1 / dT
    call simple_voigt ( x, y, u, v, a, b )
    y_dydT = y * dy_dT
    du = a * dx - b * y_dydT
    dv = a * y_dydT + b * dx

    r = 1 / v0
    q = nu * r
    sigma = nu + v0s
    sigmaX1 = sigma * x1
    s2 = sigmaX1**2
    y2 = y*y
    d = 1.0_rp / ( s2 + y2 )
    denom = OneOvSPi * d
    g = slabs1 * q * tanh1
    c = g * denom
    h2 = y * denom
    sa = g * u
    sb = g * yi
!   sb = g * yi * v
    sc = c * yi * sigmaX1
    sd = g * h2
    Slabswint = sa + sb * v - sc + sd
!   Slabswint = sa + sb - sc + sd

    ! Temperature derivative

    c = dSlabs1_dT + dtanh_dT
    dd = -2.0_rp * d * ( sigmaX1 * ( x1 * dv0S_dT + sigmaX1 * dx1_dT ) + y2 * dy_dT )
    dSlabS_dT = g * ( u * c + du ) &
!   dSlabS_dT = sa * ( c + du / du ) &
      &       + sb * ( c * v + dv ) &
!     &       + sb * ( c + dv / v ) &
      &       - sc * ( c + dd + dyi_dT + dx1_dT + dv0S_dT / sigma ) &
      &       + sd * ( c + dd + dy_dT )

    ! Spectroscopy derivatives

    y_g_dHdy = y * g * ( - b + yi * dx + denom - 2.0 * d * y * h2 )
    dSlabs_dw = y_g_dHdy / w
    dSlabs_dn = y_g_dHdy * 300.0 / T

    dx_dNu0 = - r * x - x1 * velcor
    dy_dNu0 = - r * y
    dH1_dNu0 = a * dx_dNu0 - b * dy_dNu0 + yi * ( b * dx_dNu0 + a * dy_dNu0 )
    dSlabs_dNu0 = Slabswint * dslabs1_dNu0 - & ! remember dslabs1_dNu0 is divided by slabs1
      & g * ( r * u - dH1_dNu0 + &
      &       d * x1 * velCor * ( yi * oneOvSpi - 2.0 * h2 * sigmaX1 ) )

  end subroutine Slabswint_dAll

  ! --------------------------------------------  Slabswint_Lines  -----
  elemental function Slabswint_Lines ( Nu, Slabs, tanh1, NoPolarized ) result ( Beta )

    use Voigt_m, only: Real_Simple_Voigt

    real(r8), intent(in) :: Nu    ! Frequency
    type(slabs_struct), intent(in) :: Slabs ! Frequency-independent stuff
    real(rp), intent(in) :: tanh1 ! tanh( h nu / (2 k T) )
    ! "Don't do line(L) if slabs%catalog%polarized(L)"
    logical, intent(in) :: NoPolarized

    real(rp) :: Beta ! Output

! If the molecular transition and temperature have not changed but
! frequency has, enter here to calculate sum of Beta for all lines in SLABS.

    real(rp) :: A        ! First argument for real_simple_voigt
    integer :: L         ! Line Index
    real(rp) :: SigmaX1  ! (nu + nu0s) * x1
    real(rp) :: U        ! Voigt = real part of Fadeeva
    real(r8) :: V0S      ! Pressure-shifted line center
    real(rp) :: X1, Y    ! Doppler width, ratio Pressure to Doppler widths
    real(rp) :: Yi       ! Interference coefficient
    real(rp) :: Y2       ! y**2

!{ Let $V(a,y)$ be the Voigt function ({\tt u} above), $\delta = \nu-\nu_{0_s}$,
!  $\sigma = \nu + \nu_{0_s}$, $a = x_1 \delta$, $b = x_1 \sigma$,
!  $D_1 = \frac1{\sqrt\pi} \frac1{a^2 + y^2}$ and
!  $D_2 = \frac1{\sqrt\pi} \frac1{b^2 + y^2}$.
!  For $y_i \leq 10^{-6}$,
!   {\tt Slabswint = } $ S_1 \frac{\nu}{\nu_0}
!   \tanh\left(\frac{h \nu}{2 k T}\right)
!   \tanh\left(\frac{h \nu}{2 k T}\right)
!    \left( V(a,y) + y D_2 \right)$, and for
!  $y_i > 10^{-6}$, 
!   {\tt Slabswint = } $ S_1 \frac{\nu}{\nu_0}
!   \tanh\left(\frac{h \nu}{2 k T}\right)
!    \left( V(a,y) + y_i a D_1
!     + (y - b y_i) D_2 \right)$.
!  The term $y_i a D_1$ is a one-term asymptotic approximation to
!  $L(a, y)$, where $L$ is the imaginary part of Fadeeva.  This is
!  acceptable for all $a$ and $y$ because $y_i$ is known (so far) to
!  be small relative to 1.
!  The term $y D_2$ is a one-term asymptotic approximation to $V(b,
!  y)$, while the term $b D_2$ is a is a one-term asymptotic
!  approximation to $L(b, y)$.  These are both acceptable for all
!  arguments because we know that $b$ is large.

    beta = 0.0_rp
    do l = 1, size(slabs%s)
      if ( noPolarized .and. slabs%s(l)%polarized ) cycle
      v0s = slabs%s(l)%v0s
      x1 = slabs%s(l)%x1
      y = slabs%s(l)%y
      yi = slabs%s(l)%yi
      a = x1 * real(nu-v0s,rp)
      call real_simple_voigt ( a, y, u )

      sigmaX1 = x1 * (nu + v0s)
      y2 = y*y
      if ( abs(yi) > 1.0e-6_rp ) then ! Include interference effect
        beta = beta + slabs%s(l)%slabs1 * &
          &           real(Nu / slabs%s(l)%v0, rp) * tanh1 * &
          & (u + OneOvSPi*((y - sigmaX1*yi)/(sigmaX1*sigmaX1 + y2) + yi*a/(a*a+y2)))
      else
        beta = beta + slabs%s(l)%slabs1 * &
          &           real(Nu / slabs%s(l)%v0, rp) * tanh1 * &
          & (u + OneOvSPi*(y/(sigmaX1*sigmaX1 + y2)))
      end if
    end do

  end function Slabswint_Lines

  ! -----------------------------------------  Slabswint_Lines_dT  -----
  elemental &
  subroutine Slabswint_Lines_dT ( Nu, Slabs, tanh1, dTanh_dT, &
    &                             Beta, dBeta_dT, NoPolarized )

  ! Compute single-line absorption and its derivative w.r.t. temperature,
  ! with interference, for all lines in the Slabs structure.

    use Voigt_m, only: D_Simple_Voigt

    real(r8), intent(in) :: Nu    ! Frequency
    type(slabs_struct), intent(in) :: Slabs ! Frequency-independent stuff
    real(rp), intent(in) :: Tanh1 ! tanh( h nu / (2 k T) )
    real(rp), intent(in) :: dTanh_dT ! -h nu / (2 k T^2) 1/tanh(...) dTanh(...)/dT

    real(rp), intent(out) :: Beta, dBeta_dT

    ! "Don't do line(L) if slabs%catalog%polarized(L)"
    logical, intent(in) :: NoPolarized

! If the molecular transition and temperature have not changed but
! frequency has enter here.

    real(rp) :: A               ! x1 * delta
    real(rp) :: C1, C2, C3      ! Common terms
    real(rp) :: D               ! 1 / (SigmaX1**2 + y**2)
    real(rp) :: DD              ! 1/D dD/dT
    real(rp) :: Da              ! da / dT
    real(rp) :: Delta           ! Nu-v0s
    real(rp) :: DU, DV          ! du/dT, dv/dT
    real(rp) :: Dv0s_dT         ! dv0s / dT
    real(rp) :: Dx1_dT          ! 1/x1 dx1/dT
    real(rp) :: Dy_dT           ! 1/y dy/dT
    real(rp) :: Dyi_dT          ! 1/yi dyi/dT
    integer :: L                ! Line index
    real(rp) :: Sa, Sb, Sc, Sd  ! parts of SlabsWint
    real(rp) :: Sigma           ! Nu+v0s
    real(rp) :: SigmaX1         ! sigma * x1
    real(rp) :: U, V            ! w(z) = u + i v
    real(r8) :: V0S             ! Pressure-shifted line center
    real(rp) :: X1, Y           ! Doppler width, ratio Pressure to Doppler widths
    real(rp) :: Yi              ! Interference coefficient
    real(rp) :: Y2              ! Y**2

    ! See Slabswint_dT for TeXnicalities.

    beta = 0.0_rp
    dBeta_dT = 0.0_rp
    do l = 1, size(slabs%s)

      if ( noPolarized .and. slabs%s(l)%polarized ) cycle

      v0s = slabs%s(l)%v0s
      x1 = slabs%s(l)%x1
      y = slabs%s(l)%y
      yi = slabs%s(l)%yi
      dv0s_dT = slabs%d(l)%dv0s_dT
      dx1_dT = slabs%d(l)%dx1_dT
      dy_dT = slabs%d(l)%dy_dT
      dyi_dT = slabs%d(l)%dyi_dT
      delta = nu - v0s
      a = x1 * delta
      da = x1 * ( delta * dx1_dT - dv0s_dT )
      call D_Simple_Voigt ( a, y, u, v, du, dv, da, y*dy_dT )

      sigma = nu + v0s
      sigmaX1 = sigma * x1
      y2 = y * y
      d = 1.0_rp / ( sigmaX1**2 + y2 )
      c1 = slabs%s(l)%slabs1 * real(nu / slabs%s(l)%v0,rp) * tanh1
      c2 = c1 * d * OneOvSPi
      sa = c1 * u
      sd = c2 * y
      c3 = slabs%d(l)%dSlabs1_dT + dtanh_dT
      dd = -2.0_rp * d * ( sigmaX1 * ( x1 * dv0s_dT + sigmaX1 * dx1_dT ) + y2 * dy_dT )
      beta = beta + sa + sd
!     dBeta_dT = dBeta_dT + sa * ( c3 + du / u ) + sd * ( c3 + dd + dy_dT )
      dBeta_dT = dBeta_dT + c1 * ( c3 * u + du ) + sd * ( c3 + dd + dy_dT )
      if ( abs(yi) > 1.0e-6_rp ) then
        sb = c1 * yi
        sc = c2 * sigmaX1 * yi
        beta = beta + sb * v - sc

        dBeta_dT = dBeta_dT &
          &      + sb * ( c3 * v + dv ) & ! sb * v * ( c3 + dv / v )
          &      - sc * ( c3 + dd + dyi_dT + dx1_dT + dv0s_dT / sigma ) ! &
      end if

    end do

  end subroutine Slabswint_Lines_dT

  ! ---------------------------------------  Slabswint_lines_dSpectral  -----
  elemental subroutine Slabswint_lines_dSpectral ( Nu, Slabs, T, tanh1, VelCor, &
                                & NoPolarized, SwI, dSwI_dw, dSwI_dn, dSwI_dNu0 )

  ! Compute Single Line Absorption and its first derivatives with respect
  ! to spectral parameters w, n & Nu0

  ! NOTE: slabs_prep() must compute Nu0s, x1, yi, y, w, slabs1 and dslabs1_dNu0
  !       in advance

    use SpectroscopyCatalog_m, only: Lines
    use Voigt_m, only: Simple_Voigt

    real(r8), intent(in) :: Nu       ! Measurement frequency
    type(slabs_struct), intent(in) :: Slabs ! Frequency-independent stuff, including
      ! Nu0s     ! Pressure shifted line center
      ! x1, yi, y, w ! Spectral parameters
      ! Slabs1   ! frequency-independent absorption terms
      ! dSlabs1_dNu0 ! d Slabs1 / d Nu0 * 1 / slabs1
    real(rp), intent(in) :: T, tanh1 ! temperature, tanh(h*nu/(2*k*T))
    real(rp), intent(in) :: VelCor   ! Doppler velocity correction
    logical, intent(in) :: NoPolarized ! Don't evaluate for polarized lines

    real(rp), intent(out) :: SwI, dSwI_dw, dSwI_dn, dSwI_dNu0 ! Absorption, derivs

    integer :: L         ! Line index
    real(r8) :: Nu0s     ! Pressure shifted line center from Slabs structure
    real(rp) :: x1, yi, y ! Spectral parameters from Slabs structure
    real(rp) :: x !, y   ! Re(z), Im(z)
    real(rp) :: u, v     ! Re(w), Im(w)
    real(rp) :: du, dv   ! Re(dw/dz), Im(dw/dz)
    real(rp) :: d, denom, dH1_dNu0, dx_dNu0, dy_dNu0, g, h1, h2, r
    real(rp) :: s2, sigmaX1, Swi_up, y_g_dHdy, y2

    ! See Slabs_dSpectral for TeXnicalities

    SwI = 0.0_rp
    dSwI_dw = 0.0_rp
    dSwI_dn = 0.0_rp
    dSwI_dNu0 = 0.0_rp

    do l = 1, size(slabs%s)

      if ( noPolarized .and. slabs%s(l)%polarized ) cycle

      nu0s = slabs%s(l)%v0s
      x1 = slabs%s(l)%x1
      yi = slabs%s(l)%yi
      y = slabs%s(l)%y

      x = x1 * ( nu - nu0s )
      call simple_voigt ( x, y, u, v, du, dv ) ! get w(z) and dw/dz

      r = 1.0 / slabs%s(l)%v0
      sigmaX1 = x1 * ( nu + nu0s )
      s2 = sigmaX1**2
      y2 = y**2
      d = 1.0 / ( s2 + y2 )
      denom = oneOvSpi * d
      g = slabs%s(l)%slabs1 * tanh1 * nu * r

      dx_dNu0 = - r * x - x1 * velcor
      dy_dNu0 = - r * y

      if ( yi > 1.0e-6 ) then

        h1 = u + yi * v
        h2 = ( y - sigmaX1 * yi ) * denom

        SwI_up = g * ( h1 + h2 )

        y_g_dHdy = y * g * ( - dv + yi * du + denom - 2.0 * d * y * h2 )
        dSwI_dw = dSwI_dw + y_g_dHdy / lines(slabs%catalog%lines(l))%w
        dSwI_dn = dSwI_dn + y_g_dHdy * 300.0 / T

        dH1_dNu0 = du * dx_dNu0 - dv * dy_dNu0 + yi * ( dv * dx_dNu0 + du * dy_dNu0 )
        dSwI_dNu0 = dSwI_dNu0 + &
          & SwI_up * slabs%s(l)%dslabs1_dv0 - & ! remember dslabs1_dv0 is divided by slabs1
          & g * ( r * h1 - dH1_dNu0 + &
          &       d * x1 * velCor * ( yi * oneOvSpi - 2.0 * h2 * sigmaX1 ) )

      else

        h1 = u
        h2 = y * denom

        SwI_up = g * ( h1 + h2 )

        y_g_dHdy = y * g * ( - dv + denom - 2.0 * d * y * h2 )
        dSwI_dw = dSwI_dw + y_g_dHdy / lines(slabs%catalog%lines(l))%w
        dSwI_dn = dSwI_dn + y_g_dHdy * 300.0 / T

        dSwI_dNu0 = dSwI_dNu0 + &
          & SwI_up * slabs%s(l)%dslabs1_dv0 - & ! remember dslabs1_dv0 is divided by slabs1
          & g * ( r * h1 - ( du * dx_dNu0 - dv * dy_dNu0 ) - &
          &       2.0 * d * x1 * velCor * h2 * sigmaX1 )

      end if

      SwI = SwI + SwI_up

    end do

  end subroutine Slabswint_Lines_dSpectral

  ! ----------------------------------------  Slabswint_Lines_dAll  -----
  elemental &
  subroutine Slabswint_Lines_dAll ( Nu, Slabs, T, tanh1, &
    &                         dtanh_dT, velCor, NoPolarized, &
    &                   Beta, dBeta_dT, dBeta_dw, dBeta_dn, dBeta_dNu0 )

    use SpectroscopyCatalog_m, only: Lines
    use Voigt_m, only: Simple_Voigt

    real(r8), intent(in) :: Nu         ! Measurement frequency
    type(slabs_struct), intent(in) :: Slabs ! Frequency-independent stuff
    real(rp), intent(in) :: T          ! Temperature, K
    real(rp), intent(in) :: tanh1      ! tanh(h nu / 2 k T), 
    real(rp), intent(in) :: dtanh_dT   ! 1/tanh(...) d tanh(h nu / 2 k T) / dT
      ! = h nu / ( 2 k t^2 ) (tanh(...) - 1/tanh(...) )
    real(rp), intent(in) :: velCor     ! Doppler correction = d v0s / d v0
    logical, intent(in) :: NoPolarized ! Don't evaluate for polarized lines
    real(rp), intent(out) :: Beta, dBeta_dT, dBeta_dw, dBeta_dn, dBeta_dNu0

    ! Stuff from Slabs structure
    real(rp) :: dslabs1_dNu0 ! 1/slabs1 dslabs1 / dNu0
    real(rp) :: dv0s_dT    ! dv0s / dT ( not 1/V0s dv0s / dT !)
    real(rp) :: dx1_dT     ! 1/x1 dx1 / dT
    real(rp) :: dy_dT      ! 1/y dy / dT
    real(rp) :: dyi_dT     ! 1/yi dyi / dT
    real(r8) :: v0s        ! Pressure-shifted line center
    real(rp) :: x1, yi, y, w ! Spectral parameters

    ! Other local variables
    integer :: L        ! Line index
    real(rp) :: A, B    ! Re(dw/dz), Im(dw/dz)
    real(rp) :: Beta_up ! Beta update
    real(rp) :: C2, C3  ! Common terms
    real(rp) :: D       ! 1 / (SigmaX1**2 + y**2)
    real(rp) :: DD      ! 1/D dD/dT, 
    real(rp) :: Delta   ! Nu - v0s
    real(rp) :: Denom   ! D / Sqrt(pi)
    real(rp) :: DH1_dNu0 ! dH1 / dNu0
    real(rp) :: Du, Dv  ! du/dT, du/dT
    real(rp) :: Dx      ! d(x1*delta)
    real(rp) :: Dx_dNu0 ! dx / dNu0
    real(rp) :: Dy_dNu0 ! dy / dNu0
    real(rp) :: G       ! Slabs1 * tanh1 * nu / v0
    real(rp) :: H2      ! y * Denom
    real(rp) :: Q       ! nu/v0
    real(rp) :: R       ! 1/v0
    real(rp) :: Sa, Sb, Sc, Sd  ! parts of Slabswint
    real(rp) :: Sigma   ! Nu + v0s
    real(rp) :: SigmaX1 ! Sigma * x1
    real(rp) :: S2      ! ( sigma * x1 ) ** 2
    real(rp) :: U, V    ! Voigt, Lorentz
    real(rp) :: X       ! X1 * Delta
    real(rp) :: Y2      ! y**2
    real(rp) :: Y_dYdT  ! y * dy/dT
    real(rp) :: Y_G_dHdy ! Y * G * (du/dy + dH2/dy)

! See Slabs_dT and Slabs_dSpectral for TeXnicalities.

    Beta = 0.0
    dBeta_dT = 0.0
    dBeta_dw = 0.0
    dBeta_dn = 0.0
    dBeta_dNu0 = 0.0

    do l = 1, size(slabs%s)

      if ( noPolarized .and. slabs%s(l)%polarized ) cycle

      dslabs1_dNu0 = slabs%s(l)%dslabs1_dv0
      dv0s_dT = slabs%d(l)%dv0s_dT
      dx1_dT = slabs%d(l)%dx1_dT
      dy_dT = slabs%d(l)%dy_dT
      dyi_dT = slabs%d(l)%dyi_dT
      v0s = slabs%s(l)%v0s
      w = lines(slabs%catalog%lines(l))%w
      x1 = slabs%s(l)%x1
      yi = slabs%s(l)%yi
      y = slabs%s(l)%y

      ! Absorption

      delta = Nu - v0s
      x = x1 * delta
      dx = x * dx1_dT - x1 * dv0s_dT ! remember, dx1_dT = 1/x1 dx1 / dT
      call simple_voigt ( x, y, u, v, a, b )
      y_dydT = y * dy_dT
      du = a * dx - b * y_dydT
      dv = a * y_dydT + b * dx

      r = 1.0 / slabs%s(l)%v0
      q = nu * r
      sigma = nu + v0s
      sigmaX1 = sigma * x1
      s2 = sigmaX1**2
      y2 = y*y
      d = 1.0_rp / ( s2 + y2 )
      denom = OneOvSPi * d
      g = slabs%s(l)%slabs1 * q * tanh1
      c2 = g * denom
      c3 = slabs%d(l)%dSlabs1_dT + dtanh_dT
      dd = -2.0_rp * d * ( sigmaX1 * ( x1 * dv0s_dT + sigmaX1 * dx1_dT ) + y2 * dy_dT )
      h2 = y * denom
      sa = g * u
      sd = g * h2

      if ( abs(yi) > 1.0e-6_rp ) then

        sb = g * yi
!       sb = g * yi * v
        sc = c2 * yi * sigmaX1
        beta_up = sa + sb * v - sc + sd
!       beta_up = sa + sb - sc + sd

        ! Temperature derivative

        dBeta_dT = dBeta_dT &
          &      + g * ( c3 * u + du ) &
!         &      + sa * ( c3 + du / u ) &
          &      + sb * ( c3 * v + dv ) &
!         &      + sb * ( c3 + dv / v ) &
          &      - sc * ( c3 + dd + dyi_dT + dx1_dT + dv0S_dT / sigma ) &
          &      + sd * ( c3 + dd + dy_dT )

        ! Spectroscopy derivatives

        y_g_dHdy = y * g * ( - b + yi * dx + denom - 2.0 * d * y * h2 )
        dBeta_dw = dBeta_dw + y_g_dHdy / w
        dBeta_dn = dBeta_dn + y_g_dHdy * 300.0 / T

        dx_dNu0 = - r * x - x1 * velcor
        dy_dNu0 = - r * y
        dH1_dNu0 = a * dx_dNu0 - b * dy_dNu0 + yi * ( b * dx_dNu0 + a * dy_dNu0 )
        dBeta_dNu0 = dBeta_dNu0 &
          & + beta_up * dslabs1_dNu0 - & ! remember dslabs1_dNu0 is divided by slabs1
          &   g * ( r * u - dH1_dNu0 + &
          &         d * x1 * velCor * ( yi * oneOvSpi - 2.0 * h2 * sigmaX1 ) )

      else

        beta_up = sa + sd

        ! Temperature derivative

        dbeta_dT = dbeta_dT &
          &      + g * ( c3 * u + du ) &
!         &      + sa * ( c3 + du / u ) &
          &      + sd * ( c3 + dd + dy_dT )

        ! Spectroscopy derivatives

        y_g_dHdy = y * g * ( - b + denom - 2.0 * d * y * h2 )
        dBeta_dw = dBeta_dw + y_g_dHdy / w
        dBeta_dn = dBeta_dn + y_g_dHdy * 300.0 / T

        dx_dNu0 = - r * x - x1 * velcor
        dy_dNu0 = - r * y
        dBeta_dNu0 = dBeta_dNu0 + &
          & beta_up * dslabs1_dNu0 - & ! remember dslabs1_dNu0 is divided by slabs1
          &   g * ( r * u - a * dx_dNu0 + b * dy_dNu0 - &
          &         2.0 * d * x1 * velCor * h2 * sigmaX1 )

      end if

      beta = beta + beta_up

    end do

  end subroutine Slabswint_Lines_dAll

  ! ----------------------------------------------  Voigt_Lorentz  -----

  elemental subroutine Voigt_Lorentz ( dNu,  Nu0,  x1,  yi,  y,  w,  t,  tanh1, slabs1,  &
                         &   VL, dslabs1_dNu0,  dVL_dw,  dVL_dn,  dVL_dNu0 )

! Compute the Voigt/Lorentz function and its first derivatives with respect
! to spectral parameters: w, n & Nu0

! NOTE: Before calling this routine, the user needs to call slabs_prep()
!       routine to compute dslabs1_dNu0

    use Voigt_m, only: Simple_Voigt

    real(r8), intent(in) :: dNu, nu0
    real(rp), intent(in) :: x1, yi, y, w, t, tanh1, slabs1
    real(rp), intent(in) :: dslabs1_dNu0 ! 1 / Slabs1 * d Slabs1 / d Nu0 

    real(rp), intent(out) :: VL, dVL_dw, dVL_dn, dVL_dNu0

    real(rp) :: xj, zj, q, y2, u, v, up1, up2, dn1, dn2, dup1, &
     &          dup2, ddn1, ddn2, dy_dw, dy_dn, dSum_dw, dSum_dn, &
     &          slabs2

    q = 1.0_rp + dNu / Nu0

    y2 = y * y
    xj = x1 * dNu
    zj = x1 * (2.0_r8 * Nu0 + dNu)
    dn1 = zj * zj + y2
    up1 = y - zj * yi

! Van Vleck - Wieskopf (VVW) line shape with Voigt

    call simple_voigt ( xj, y, u, v )
    dup1 = up1 * OneOvSPi / dn1 + yi * v
    ddn2 = u + dup1
    slabs2 = slabs1 * tanh1
    VL = slabs2 * ddn2 * q            ! This is the Voigt + VVW correction

    dn2 = xj * xj + y2
    up2 = y + yi * xj

    dy_dw = y / w
    dup1 = dy_dw
    ddn1 = 2.0 * y * dy_dw
    dup2 = dy_dw
    ddn2 = 2.0 * y * dy_dw

    dSum_dw = (dn1*dup1-up1*ddn1)/(dn1*dn1) + &
   &          (dn2*dup2-up2*ddn2)/(dn2*dn2)

    dVL_dw = OneOvSPi * slabs2 * q * dSum_dw

    dy_dn = y * Log(300.0/t)
    dup1 = dy_dn
    ddn1 = 2.0 * y * dy_dn
    dup2 = dy_dn
    ddn2 = 2.0 * y * dy_dn
    dSum_dn = (dn1*dup1-up1*ddn1)/(dn1*dn1) + &
   &          (dn2*dup2-up2*ddn2)/(dn2*dn2)

    dVL_dn = OneOvSPi * slabs2 * q * dSum_dn

!    dup2 =  yi * x1               !  x1 = -dxj_dNu0
!    ddn2 = -2.0 * xj * x1         !  x1 = -dxj_dNu0
!    dSum_dNu0 = (dn2*dup2-up2*ddn2)/(dn2*dn2)
!    dq_dNu0 = -(dNu+Nu0)/(Nu0*Nu0)

!    dVL_dNu0 = OneOvSPi * (dslabs1_dNu0 * q * Sum          + &
!              &         slabs2 * dq_dNu0 * Sum + &
!              &         slabs2 * q * dSum_dNu0)

!    dVL_dNu0 = OneOvSPi * slabs2 * q * dSum_dNu0
    dVL_dNu0 = VL * (dslabs1_dNu0 - 1.0_rp / Nu0) &
           & + OneOvSPi * slabs2 * q * ( &
           &  (-yi*x1*(xj*xj+y2)+2.0_rp*(y+yi*xj)*xj*x1)/(xj*xj+y2)**2 &
           & + (yi*x1*  dn1     -2.0_rp*   up1   *x1*zj)/dn1**2)

  end subroutine Voigt_Lorentz

  ! -------------------------------------------------  Slabs_prep  -----
  pure subroutine Slabs_prep ( t, m, v0, el, w, ps, p, n, ns, i, q, delta, gamma, &
                      &   n1, n2, velCor, useYi, &
                      &   v0s, x1, y, yi, slabs1, dslabs1 )

! This function computes a single line type absorption coefficient
! WITH INTERFERENCE ! using predominantly data from the Pickett catalogue.

! ** UPDATED: Jul/3/97  To include Hugh Pumphrey's Pressure Shift effects
! >>2004-03-18 WV Snyder Use 1 - exp(-v0/300/Boltzmhz) in denominator instead
!                        of 1 - exp(-v0s/300/Boltzmhz).
! >>2005-06-16 WV Snyder Use nu0s/velCor instead of nu0s in Beta_v.

    use Physics, only: H_OVER_K, k, SpeedOfLight
    use Constants, only: Ln10, Sqrtln2, SqrtPi

! inputs:

    real(rp), intent(in) :: T        ! Temperature K
    real(r8), intent(in) :: M        ! Molecular mass amu
    real(r8), intent(in) :: V0       ! Line center frequency MHz
    real(r8), intent(in) :: El       ! Lower state energy cm-1
    real(r8), intent(in) :: W        ! Collision broadening parameter
                                     ! MHz/mbar at 300 K
    real(r8), intent(in) :: Ps       ! Pressure shift parameter in MHz/mbar
    real(rp), intent(in) :: P        ! Pressure mbar
    real(r8), intent(in) :: N        ! Temperature power dependence of w
    real(r8), intent(in) :: Ns       ! Temperature power dependence of ps
    real(r8), intent(in) :: I        ! Integrated spectral intensity
                                     ! Log(nm**2 MHz) at 300 K
    real(r8), intent(in) :: Q(3)     ! Logarithm of the partition function
                                     ! At 300 , 225 , and 150 K
    real(r8), intent(in) :: Delta    ! Delta interference coefficient at 300K 1/mb
    real(r8), intent(in) :: Gamma    ! Gamma               "
    real(r8), intent(in) :: N1       ! Temperature dependency of delta
    real(r8), intent(in) :: N2       ! Temperature dependency of gamma
    real(r8), intent(in) :: VelCor   ! Doppler velocity correction term
    logical, intent(in) :: UseYi     ! delta + gamma > 0.0

! outputs:

    real(r8), intent(out) :: V0s     ! Pressure shifted line position
    real(r8), intent(out) :: X1      ! Sqrt(Ln(2))/Doppler half width MHz
    real(r8), intent(out) :: Y       ! Sqrt(Ln(2))*collision width /
                                     !             doppler width
    real(r8), intent(out) :: Yi      ! Interference contribution
    real(r8), intent(out) :: Slabs1  ! Frequency independent piece of slabs
    real(r8), intent(out) :: Dslabs1 ! Derivative of slabs1 w.r.t. v0 / slabs1

!  The above outputs along with frequency offset are used with routine
!  SLABSWINT to compute a Single Line ABSorption in 1/Km units. With unit
!  mixing ratio.

! Internal constants:

!  i2abs    - converts intensity into absorption
!  dc       - sqrt(amu/K) used to calculate doppler width
!  boltzcm  - boltzmann constant cm-1/K = k / (hc) 1e6/100
!  boltzmhz - boltzmann constant MHz/K = k / h
!  sqrtln2  - sqrt(ln(2))

    real(rp), parameter :: I2abs = sqrtln2 / ( sqrtPi * 1.0e13 * k )
!   real(rp), parameter :: I2abs = 3.402155052e9_rp ! using above constants
!   real(rp), parameter :: I2abs = 3.402136078e9_rp ! Zvi's original value
    real(rp), parameter :: Dc = 3.58116514e-7_rp ! sqrt(1000 k ln 4 avogadro) / c
!   real(rp), parameter :: Dc = 3.58117369e-7_rp ! Zvi's original value
    real(rp), parameter :: BoltzMHz = 1.0_rp / H_over_k
    real(rp), parameter :: Boltzcm = boltzMHz / SpeedOfLight * 1.0e6 / 100.0
    real(rp), parameter :: Oned300 = 1.0_rp/300.0_rp

    real(rp), parameter :: LT2 = 2.35218251811136_rp      ! Log10(225)
    real(rp), parameter :: LT3 = 2.47712125471966_rp      ! Log10(300)
    real(rp), parameter :: LN300 = lt3 * ln10

    real(rp), parameter :: Tl1 = 0.176091259055681_rp     ! Log10(225/150)
    real(rp), parameter :: Tl2 = 0.124938736608300_rp     ! Log10(300/225)

! Internal data:

    real(rp) :: betae, betav, expd, expn, log_T, onedt
    real(rp) :: Q_Log_a, Q_Log_b, t3t, Wd, z1, z2

! The action begins here

    onedt = 1.0_rp / t
    log_T = log(t)
    t3t = ln300 - log_T ! log(300/T)

!{ $y_i = p \left( \delta \left( \frac{300}T \right)^{n_1} +
!                  \gamma \left( \frac{300}T \right)^{n_2} \right)$.

    yi = 0.0
    if ( useYi ) yi = p * (delta*exp(n1*t3t) + gamma*exp(n2*t3t))

!{ $w_d = \nu_0 d_c \sqrt{\frac{T}M}$.

    Wd = v0 * dc * Sqrt(t/m)

!{ $x_1 = \frac{\sqrt{\ln 2}}{w_d}$.
!  $\frac{\partial x_1}{\partial \nu_0} = -\frac{x_1}{\nu_0}$.

    x1 = real(sqrtln2,rp) / Wd

!{ $y = x_1 w p \left( \frac{300}T \right) ^n$.
!  $\frac{\partial y}{\partial \nu_0} = -\frac{y}{\nu_0}$.

    y = x1 * w * p * exp(n*t3t)

!{ $\nu_{0_s} = \nu_0 + p_s p \left( \frac{300}T \right) ^{n_s}$.
!  $\frac{\partial \nu_{0_s}}{\partial \nu_0} = 0$.

    v0s = v0 + ps * p * exp(ns*t3t) ! VelCor is multiplied below

!{ $\beta_e = 10^{-4} \frac{h}k c\, e_l$.
!  $\beta_v = \frac{h}k \frac{\nu_{0_s}}{v_c}$.

    betae = el / boltzcm
    betav = v0s / boltzmhz

!{ Now $\nu_{0_s} = v_c \left[ \nu_0 + p_s p \left( \frac{300}T \right) ^{n_s} \right]$.
!  $\frac{\partial \nu_{0_s}}{\partial \nu_0} = v_c$.

    v0s = velcor * v0s

    if ( t < 225.0_rp ) then
      q_log_b = (q(2)-q(3)) / tl1
      q_log_a = q(2) - q(1) - lt2 * q_log_b
    else
      q_log_b = (q(1)-q(2)) / tl2
      q_log_a = -lt3 * q_log_b
    end if
    q_log_b = 1.0 + q_log_b ! This is more useful below

    expd = EXP(-v0*(oned300/boltzmhz)) ! H
    expn = EXP(-betav*onedt) ! G
    z1 = 1.0 + expn          ! 1 + G
    z2 = 1.0 - expd          ! 1 - H

!{ $S = \frac{I_2 \, p \, 10^{i - a - b \log_{10} T
!   + \frac{\beta_e}{\ln 10}\left(\frac1{300}-\frac1T\right)}
!      (1+e^{-\frac{\beta_v}T})}
!   {T w_d \left(1 - e^{-\frac{h \nu_0}{300 k}}\right)} =
!  I_2 \, p \, e^{(i-a) \ln 10 -(b+1) \log T +
!    \beta_e \left( \frac1{300} - \frac1T \right)}
!  \frac{(1+G)} {w_d (1-H)}$, where
!  $G = e^{-\frac{\beta_v}T}$, $H = e^{-\frac{h \nu_0}{300 k}}$ and
!  $S =$ {\tt slabs1}.

    slabs1 = i2abs * p * &
      & exp((i - q_log_a)*ln10 - q_log_b * log_T + betae * (oned300 - onedt)) &
      & * z1 / (Wd * z2)

!{ $\frac1S\,\frac{\partial S}{\partial \nu_0} =
!   - \left[ \left( \frac{G}{(1+G)T}
!   + \frac{H}{300(1-H)} \right) \frac{h}k + \frac1{\nu_0} \right]$.

    dslabs1 = - ( ( expn / (t * z1) + expd / (300.0_rp * z2)) / &
      & boltzmhz + 1.0 / v0 )

  end subroutine Slabs_prep

  ! ------------------------------------------  Slabs_prep_struct  -----
  elemental subroutine Slabs_prep_struct ( T, P, Catalog, VelCor, Derivs, Slabs )
  ! Fill all the fields of the Slabs structure

    use SpectroscopyCatalog_m, only: Lines

    ! inputs:

    real(rp), intent(in) :: T        ! Temperature K
    real(rp), intent(in) :: P        ! Pressure
    type(catalog_t), intent(in) :: Catalog ! The spectroscopy
    real(rp), intent(in) :: VelCor   ! Doppler velocity correction term, 
                                     ! 1 - losVel / C
    logical, intent(in) :: Derivs    ! "Setup for derivative calculations"

    ! output:

    type(slabs_struct), intent(inout) :: Slabs ! inout so as not to clobber
                                     ! pointer associations

    integer :: I ! A loop index
    integer :: L ! Index in the Lines array

    slabs%useYi = .false.
    !ocl independent
    !ocl temp(l)
    do i = 1, size(catalog%lines)
      l = catalog%lines(i)
      slabs%useYi = slabs%useYi .or. lines(l)%useYi
      slabs%s(i)%v0 = lines(l)%v0
      if ( associated(slabs%catalog%polarized) ) &
        & slabs%s(i)%polarized = slabs%catalog%polarized(i)
      if ( derivs ) then
        call slabs_prep_dT ( t, catalog%mass, &
          & lines(l)%v0, lines(l)%el, lines(l)%w, lines(l)%ps, p, &
          & lines(l)%n, lines(l)%ns, lines(l)%str, catalog%QLOG(1:3), &
          & lines(l)%delta, lines(l)%gamma, lines(l)%n1, lines(l)%n2, &
          & velCor, lines(l)%useYi, &
          & slabs%s(i)%v0s, slabs%s(i)%x1, slabs%s(i)%y, &
          & slabs%s(i)%yi, slabs%s(i)%slabs1, &
          & slabs%s(i)%dslabs1_dv0, &
          & slabs%d(i)%dv0s_dT, slabs%d(i)%dx1_dT, &
          & slabs%d(i)%dy_dT, slabs%d(i)%dyi_dT, &
          & slabs%d(i)%dslabs1_dT )
      else
        call slabs_prep ( t, catalog%mass, &
          & lines(l)%v0, lines(l)%el, lines(l)%w, lines(l)%ps, p, &
          & lines(l)%n, lines(l)%ns, lines(l)%str, catalog%QLOG(1:3), &
          & lines(l)%delta, lines(l)%gamma, lines(l)%n1, lines(l)%n2, &
          & velCor, lines(l)%useYi, &
          & slabs%s(i)%v0s, slabs%s(i)%x1, slabs%s(i)%y, &
          & slabs%s(i)%yi, slabs%s(i)%slabs1, &
          & slabs%s(i)%dslabs1_dv0 )
      end if
    end do ! i = 1, size(catalog%lines)

  end subroutine Slabs_prep_struct

  ! -----------------------------------  Slabs_prep_struct_offset  -----
  elemental subroutine Slabs_prep_struct_offset ( T, P, Catalog, VelCor, &
    & Derivs, Slabs, DV0, DW, DN )
  ! Fill all the fields of the Slabs structure

    use SpectroscopyCatalog_m, only: Lines

    ! inputs:

    real(rp), intent(in) :: T        ! Temperature K
    real(rp), intent(in) :: P        ! Pressure
    type(catalog_t), intent(in) :: Catalog ! The spectroscopy
    real(rp), intent(in) :: VelCor   ! Doppler velocity correction term, 
                                     ! 1 - losVel / C
    logical, intent(in) :: Derivs    ! "Setup for derivative calculations"
    real(r8), intent(in) :: DV0, DW, DN  ! Offsets from catalog%v0, %w, %n
    ! output:

    type(slabs_struct), intent(inout) :: Slabs ! inout so as not to clobber
                                     ! pointer associations

    integer :: I ! A loop index
    integer :: L ! Index in the Lines array

    slabs%useYi = .false.
    !ocl independent
    !ocl temp(l)
    do i = 1, size(catalog%lines)
      l = catalog%lines(i)
      slabs%useYi = slabs%useYi .or. lines(l)%useYi
      slabs%s(i)%v0 = lines(l)%v0
      if ( derivs ) then
        call slabs_prep_dT ( t, catalog%mass, &
          & lines(l)%v0+dv0, lines(l)%el, lines(l)%w+dw, lines(l)%ps, p, &
          & lines(l)%n+dn, lines(l)%ns, lines(l)%str, catalog%QLOG(1:3), &
          & lines(l)%delta, lines(l)%gamma, lines(l)%n1, lines(l)%n2, &
          & velCor, lines(l)%useYi, &
          & slabs%s(i)%v0s, slabs%s(i)%x1, slabs%s(i)%y, &
          & slabs%s(i)%yi, slabs%s(i)%slabs1, &
          & slabs%s(i)%dslabs1_dv0, &
          & slabs%d(i)%dv0s_dT, slabs%d(i)%dx1_dT, &
          & slabs%d(i)%dy_dT, slabs%d(i)%dyi_dT, &
          & slabs%d(i)%dslabs1_dT )
      else
        call slabs_prep ( t, catalog%mass, &
          & lines(l)%v0+dv0, lines(l)%el, lines(l)%w+dw, lines(l)%ps, p, &
          & lines(l)%n+dn, lines(l)%ns, lines(l)%str, catalog%QLOG(1:3), &
          & lines(l)%delta, lines(l)%gamma, lines(l)%n1, lines(l)%n2, &
          & velCor, lines(l)%useYi, &
          & slabs%s(i)%v0s, slabs%s(i)%x1, slabs%s(i)%y, &
          & slabs%s(i)%yi, slabs%s(i)%slabs1, &
          & slabs%s(i)%dslabs1_dv0 )
      end if
    end do ! i = 1, size(catalog%lines)

  end subroutine Slabs_prep_struct_offset

  ! ----------------------------------------------  Slabs_prep_DT  -----
  pure subroutine Slabs_prep_dT ( t, m, v0, el, w, ps, p, n, ns, i, q, delta, gamma, &
                         &   n1, n2, velCor, useYi, &
                         &   v0s, x1, y, yi, slabs1, dslabs1_dv0, &
                         &   dv0s_dT, dx1_dT, dy_dT, dyi_dT, dslabs1_dT )

! This function computes a single line type absorption coefficient
! WITH INTERFERENCE ! using predominantly data from the Pickett catalogue.
! Compute the derivatives with respect to temperature, too.

! ** UPDATED: Jul/3/97  To include Hugh Pumphrey's Pressure Shift effects

    use Physics, only: H_OVER_K, k, SpeedOfLight
    use Constants, only: Ln10, Sqrtln2, SqrtPi

! inputs:

    real(rp), intent(in) :: T        ! Temperature K
    real(r8), intent(in) :: M        ! Molecular mass amu
    real(r8), intent(in) :: V0       ! Line center frequency MHz
    real(r8), intent(in) :: El       ! Lower state energy cm-1
    real(r8), intent(in) :: W        ! Collision broadening parameter
                                     ! MHz/mbar at 300 K
    real(r8), intent(in) :: Ps       ! Pressure shift parameter in MHz/mbar
    real(rp), intent(in) :: P        ! Pressure mbar
    real(r8), intent(in) :: N        ! Temperature power dependence of w
    real(r8), intent(in) :: Ns       ! Temperature power dependence of ps
    real(r8), intent(in) :: I        ! Integrated spectral intensity
                                     ! Log(nm**2 MHz) at 300 K
    real(r8), intent(in) :: Q(3)     ! Logarithm of the partition function
                                     ! At 300 , 225 , and 150 K
    real(r8), intent(in) :: Delta    ! Delta interference coefficient at 300K 1/mb
    real(r8), intent(in) :: Gamma    ! Gamma               "
    real(r8), intent(in) :: N1       ! Temperature dependency of delta
    real(r8), intent(in) :: N2       ! Temperature dependency of gamma
    real(rp), intent(in) :: VelCor   ! Doppler velocity correction term
    logical, intent(in) :: UseYi     ! delta + gamma > 0

! outputs:

    real(r8), intent(out) :: V0s     ! Pressure shifted line position
    real(r8), intent(out) :: X1      ! Sqrt(Ln(2))/Doppler half width MHz
    real(r8), intent(out) :: Y       ! Sqrt(Ln(2))*collision width /
                                     !             doppler width
    real(r8), intent(out) :: Yi      ! Interference contribution
    real(r8), intent(out) :: Slabs1  ! Frequency independent piece of slabs
    real(r8), intent(out) :: Dslabs1_dv0 ! Derivative of slabs1 w.r.t. v0 / slabs1

!  The above outputs along with frequency offset are used with routine
!  SLABSWINT to compute a Single Line ABSorption in 1/Km units. With unit
!  mixing ratio.

! Derivatives with respect to temperature:

    real(r8), intent(out) :: dv0s_dT    ! dv0s/dT
    real(r8), intent(out) :: dx1_dT     ! 1/x1 dx1/dT
    real(r8), intent(out) :: dy_dT      ! 1/y dy/dT
    real(r8), intent(out) :: dyi_dT     ! 1/yi dyi/dT
    real(r8), intent(out) :: dslabs1_dT ! 1/slabs1 dslabs1/dT

! Internal constants:

!  i2abs    - converts intensity into absorption
!  dc       - sqrt(amu/K) used to calculate doppler width
!  boltzcm  - boltzmann constant cm-1/K
!  boltzmhz - boltzmann constant MHz/K = k/h
!  sqrtln2  - sqrt(ln(2))

    real(rp), parameter :: I2abs = sqrtln2 / ( sqrtPi * 1.0e13 * k )
!   real(rp), parameter :: I2abs = 3.402155052e9_rp ! using above constants
!   real(rp), parameter :: I2abs = 3.402136078e9_rp ! Zvi's original value
    real(rp), parameter :: Dc = 3.58116514e-7_rp ! sqrt(1000 k ln 4 avogadro) / c
!   real(rp), parameter :: Dc = 3.58117369e-7_rp ! Zvi's original value
    real(rp), parameter :: BoltzMHz = 1.0_rp / H_over_k
    real(rp), parameter :: Boltzcm = boltzMHz / SpeedOfLight * 1.0e6 / 100.0
    real(rp), parameter :: Oned300 = 1.0_rp/300.0_rp

    real(rp), parameter :: LT2 = 2.35218251811136_rp      ! Log10(225)
    real(rp), parameter :: LT3 = 2.47712125471966_rp      ! Log10(300)
    real(rp), parameter :: LN300 = lt3 * ln10

    real(rp), parameter :: Tl1 = 0.176091259055681_rp     ! Log10(225/150)
    real(rp), parameter :: Tl2 = 0.124938736608300_rp     ! Log10(300/225)

! Internal data:

    real(rp) :: Betae, Betav, dBetav_dT, onedt, expn, expd
    real(rp) :: Q_Log_a, Q_Log_b ! Q_Log = a + b * log10(T)
    real(rp) :: Log_T         ! log(t)
    real(rp) :: t3t           ! log(300/t) = log(300) - log(t)
    real(rp) :: Wd, DWd_dT    ! dWd_dT is actually -dWd_dT/Wd
    real(rp) :: Z1, Z2        ! Temporaries

! The action begins here

    onedt = 1.0_rp / t
    log_t = log(t)
    t3t = ln300 - log_t  ! log(300/T)

!{ $y_i = p \left( \delta \left( \frac{300}T \right)^{n_1} +
!                  \gamma \left( \frac{300}T \right)^{n_2} \right)$.
!  $\frac{\partial y_i}{\partial T} = -\frac{p}T \left(
!                    n_1 \delta \left( \frac{300}T \right)^{n_1} +
!                    n_2 \gamma \left( \frac{300}T \right)^{n_2} \right)$.

    if ( useYi ) then
      z1 = delta*exp(n1*t3t)
      z2 = gamma*exp(n2*t3t)
      yi = ( z1 + z2 )
      dyi_dT = -onedt * ( n1 * z1 + n2 * z2 ) / yi ! 1/yi dyi/dT
      yi = p * yi
    else ! yi == 0.0
      ! If yi == 0.0, dyi_dT will necessarily be zero.  The 1/yi cancels
      ! a yi in a numerator where dyi_dT is used, so 0.0 is the correct
      ! result.  We don't need a fancy l'Hospital argument to justify it.
      yi = 0.0
      dyi_dT = 0.0
    end if

!{ $w_d = \nu_0 d_c \sqrt{\frac{T}M}$.  The $\nu_0$ term should
!  really be $\nu$.  We approximate $\nu$ by $\nu_0$ so that we can use
!  this routine outside the frequency loop.  Thus
!  $-\frac1{w_d}\frac{\partial w_d}{\partial T} = 
!   - \frac1{2 T}$.
!  $-\frac1{w_d}\frac{\partial w_d}{\partial T}$ is what's actually
!  useful later.

    Wd = v0 * dc * sqrt(t/m)
    dWd_dT = - 0.5 * onedt ! Actually -dWd_dT/Wd

!{ $x_1 = \frac{\sqrt{\ln 2}}{w_d}$.
!  $\frac1{x_1}\frac{\partial x_1}{\partial T} =
!   -\frac1{w_d} \frac{\partial w_d}{\partial T}
!   = -\frac1{2T}$.
!  We don't calculate $x$ = $x_1 ( \nu - \nu_{0_s} )$ here because it
!  depends on frequency.  Here's $\frac1x \frac{\partial x}{\partial T} =
!  \frac1{x_1}\frac{\partial x_1}{\partial T} - \frac1{\nu - \nu_{0_s}}
!  \frac{\partial \nu_{0_s}}{\partial T}$ anyway, for reference.
!  $\frac{\partial x_1}{\partial \nu_0} = -\frac{x_1}{\nu_0}$.

    x1 = real(sqrtln2,rp) / Wd
    dx1_dT = dWd_dT ! 1/x1 dx1/dT

!{ $y = x_1 w p \left( \frac{300}T \right)^n$.
!  $\frac1y \frac{\partial y}{\partial T} =
!    \left( \frac1{x_1} \frac{\partial x_1}{\partial T} - \frac{n}T \right)
!    = -\frac1{2T}(1+2n)$.
!  $\frac{\partial y}{\partial \nu_0} = -\frac{y}{\nu_0}$.

    y = x1 * w * p * exp(n*t3t)
    dy_dT = ( dx1_dT - n * onedt ) ! 1/y dy/dT

    if ( t < 225.0_rp ) then
      q_log_b = (q(2)-q(3)) / tl1
      q_log_a = q(2) - q(1) - lt2 * q_log_b
    else
      q_log_b = (q(1)-q(2)) / tl2
      q_log_a = -lt3 * q_log_b
    end if
    q_log_b = -q_log_b - 1.0 ! This is what's interesting later

!{ $\nu_{0_s} = \nu_0 + p_s p \left( \frac{300}T \right)^{n_s}$.
!  $\frac{\partial \nu_{0_s}}{\partial T} = \frac{-n_s}T ( \nu_{0_s} - \nu_0 )$.
!  $\frac{\partial \nu_{0_s}}{\partial \nu_0} = 0$.

    if ( ps /= 0.0_r8 ) then
      v0s = ps * p * exp(ns*t3t)
      dv0s_dT = -ns * v0s * onedt
      v0s = v0 + v0s
    else
      v0s = v0
      dv0s_dT = 0.0
    end if

!{ $\beta_e = 10^{-4} \frac{h}k c\, e_l$. $\beta_v = \frac{h}k \nu_{0_s}$.
!  $\frac{\partial \beta_v}{\partial T} =
!    \frac{h}k \frac{\partial \nu_{0_s}}{\partial T}$.

    betae = el / boltzcm
    betav = v0s / boltzmhz ! should not be velocity corrected
    dBetav_dT = dv0s_dT / boltzmhz

!{ Now, $\nu_{0_s} = v_c \left[ \nu_0 + p_s p \left( \frac{300}T \right)^{n_s} \right]$.
!  $\frac{\partial \nu_{0_s}}{\partial T} = \frac{-n_s}T ( \nu_{0_s} - v_c \nu_0 )$.
!  $\frac{\partial \nu_{0_s}}{\partial \nu_0} = v_c$.

    v0s = velcor * v0s
    dv0s_dT = velcor * dv0s_dT

!{ Write {\tt slabs1} $= \frac{I_2 \, p \, 10^{i - a - b \log_{10} T
!   + \frac{\beta_e}{\ln 10}\left(\frac1{300}-\frac1T\right)}
!      (1+e^{-\frac{\beta_v}T})}
!   {T w_d \left(1 - e^{\frac{\beta_v}{300}}\right)}$
!  as $S = f \frac{T^{-b-1} e^{-\frac{\beta_e}T} (1+G)}
!                                    {w_d}$, where
!  $f = \frac{I_2 \, p \, e^{(i-a) \ln 10 + \frac{\beta_e}{300}}}{1-H}$
!  is independent of T,
!  $H = e^{-\frac{h \nu_0}{300 k}}$ and $G = e^{-\frac{\beta_v}T}$.  Then
!  $\frac1S \frac{\partial S}{\partial T} =
!    \frac1T \left( -b -1 -G_1 \frac{\partial \beta_v}{\partial T} +
!     \frac1T \left[ \beta_e + G_1 \beta_v \right ] \right )
!    - \frac1{w_d}\frac{\partial w_d}{\partial T}$, where
!  $G_1 = \frac{G}{1+G}$.

    expd = EXP(-v0*(oned300/boltzmhz)) ! H
    expn = EXP(-betav*onedt)   ! G
    z1 = 1.0 + expn            ! 1 + G
    z2 = 1.0 - expd            ! 1 - H

    ! This is rearranged to reduce the number of references to "exp".
    slabs1 = i2abs * p * z1 / ( wd * z2 ) * &
      & exp((i-q_log_a)*ln10 + betae*(oned300 -onedt) + q_log_b * log_t )

    z1 = expn / z1             ! G1
    ! 1/slabs1 dslabs1/dT:
    dslabs1_dT = onedt * ( q_log_b - z1 * dBetav_dT + &
      & onedt * ( betae + z1 * betav ) ) + &
      & dWd_dT ! Remember dWd_dT is really -dWd_dT/Wd

!{ $\frac1S \,\frac{\partial S}{\partial \nu_0} =
!   -\left[ \left( \frac{G_1}T
!    + \frac{H_1}{300} \right) \frac{h}k + \frac1{\nu_0} \right] $
!   where $G_1 = \frac{G}{1+G}$ and $H_1 = \frac{H}{1-H}$.

    dslabs1_dv0 = - ( (z1 * onedt + expd / z2 * oned300) / boltzmhz + &
      & 1.0 / v0 )

  end subroutine Slabs_prep_dT

  ! ----------------------------------------  Get_GL_Slabs_Arrays  -----
  !ocl disjoint
  pure &
  subroutine Get_GL_Slabs_Arrays ( P_path, T_path, Vel_Rel, GL_Slabs, &
    & Do_1D, LineCenter, LineCenter_ix, LineWidth, LineWidth_ix, &
    & LineWidth_TDep, LineWidth_TDep_ix, T_der_flags )

    use Molecules, only: IsExtinction
    use SpectroscopyCatalog_m, only: Lines

    real(rp), intent(in) :: p_path(:) ! Pressure in hPa or mbar
    real(rp), intent(in) :: t_path(:)

    real(rp), intent(in) :: Vel_Rel   ! Vel_Z / speedOfLight

    ! GL_Slabs needs to have been created by AllocateSlabs
    type (slabs_struct), intent(inout) :: GL_Slabs(:,:)

    logical, intent(in) :: Do_1D

    ! Line parameter offsets from catalog (path x molecule) -- intent(in):
    real(rp), optional, intent(in) :: LineCenter(:,:), LineWidth(:,:), &
      & LineWidth_TDep(:,:)
    ! Which molecules are to be offset, and where are they in Line....  Zero
    ! means no offset, otherwise, second subscript for Line... above.
    integer, optional, pointer :: LineCenter_ix(:), LineWidth_ix(:), &
      & LineWidth_TDep_ix(:)
    logical, intent(in), optional :: t_der_flags(:) ! do derivatives if present

!  ----------------
!  Local variables:
!  ----------------

    type(Catalog_T), pointer :: Catalog
    integer :: i, j, k, n, n_sps, nl, no_ele
    logical :: DoCenter, DoWidth, DoWidth_TDep, Offset, Temp_Der
    real(r8) :: DV0, DW, DN  !  Offsets to center, width, width_TDep

    real(rp) :: vel_z_correction

! Begin code:

    doCenter = .false.; doWidth = .false.; doWidth_TDep = .false.
    if ( present(lineCenter) .and. present(lineCenter_ix) ) &
      & doCenter = size(lineCenter) > 0 .and. associated(lineCenter_ix)
    if ( present(lineWidth) .and. present(lineWidth_ix) ) &
      & doWidth = size(lineWidth) > 0 .and. associated(lineWidth_ix)
    if ( present(lineWidth_TDep) .and. present(lineWidth_TDep_ix) ) &
      & doWidth_TDep = size(lineWidth_TDep) > 0 .and. associated(lineWidth_TDep_ix)
    offset = doCenter .or. doWidth .or. doWidth_TDep

    no_ele = size(p_path)
    n_sps = size(gl_slabs,2)
    n = no_ele
    if ( Do_1D ) n = n / 2

    ! opposite sign convention here from ATBD
    Vel_z_correction = 1.0_rp - vel_rel

    do i = 1, n_sps

      catalog => gl_slabs(1,i)%catalog ! gl_slabs(:,i)%catalog are all the same
      gl_slabs(:,i)%useYi = any(lines(catalog%lines)%useYi)

      nl = Size(catalog%Lines)
      if ( nl == 0 ) cycle

      if ( isExtinction(catalog%molecule) ) cycle

      do j = 1, n
        temp_der = present(t_der_flags)
        if ( temp_der ) temp_der = t_der_flags(j)
        if ( offset ) then
          dv0 = 0.0; dw = 0.0; dn = 0.0
          if ( doCenter ) then
            if ( lineCenter_ix(i) /= 0 ) dv0 = lineCenter(j,lineCenter_ix(i))
          end if
          if ( doWidth ) then
            if ( lineWidth_ix(i) /= 0 ) dw = lineWidth(j,lineWidth_ix(i))
          end if
          if ( doWidth_TDep ) then
            if ( lineWidth_TDep_ix(i) /= 0 ) &
              dn = lineWidth_TDep(j,lineWidth_TDep_ix(i))
          end if
          call slabs_prep_struct_offset ( t_path(j), p_path(j), catalog, &
            & Vel_z_correction, temp_der, gl_slabs(j,i), dv0, dw, dn )
        else
          call slabs_prep_struct ( t_path(j), p_path(j), catalog, &
            &                      Vel_z_correction, temp_der, gl_slabs(j,i) )
        end if
      end do ! j = 1, n

      if ( Do_1D ) then
        ! fill in grid points on other side with above value
        !ocl temp(k)
        do j = no_ele, no_ele/2+1, -1
          k = no_ele - j + 1
          gl_slabs(j,i)%s%v0s         = gl_slabs(k,i)%s%v0s
          gl_slabs(j,i)%s%x1          = gl_slabs(k,i)%s%x1
          gl_slabs(j,i)%s%y           = gl_slabs(k,i)%s%y
          gl_slabs(j,i)%s%yi          = gl_slabs(k,i)%s%yi 
          gl_slabs(j,i)%s%slabs1      = gl_slabs(k,i)%s%slabs1 
          gl_slabs(j,i)%s%dslabs1_dv0 = gl_slabs(k,i)%s%dslabs1_dv0
        end do ! j = no_ele, no_ele/2+1, -1

        if ( present(t_der_flags) ) then
          !ocl temp(k)
          do j = no_ele, no_ele/2+1, -1
            if ( t_der_flags(j) ) then ! do derivative stuff
              k = no_ele - j + 1
              gl_slabs(j,i)%d%dv0s_dT    = gl_slabs(k,i)%d%dv0s_dT
              gl_slabs(j,i)%d%dx1_dT     = gl_slabs(k,i)%d%dx1_dT
              gl_slabs(j,i)%d%dy_dT      = gl_slabs(k,i)%d%dy_dT
              gl_slabs(j,i)%d%dyi_dT     = gl_slabs(k,i)%d%dyi_dT
              gl_slabs(j,i)%d%dslabs1_dT = gl_slabs(k,i)%d%dslabs1_dT
            end if
          end do ! j = no_ele, no_ele/2+1, -1
        end if
      end if

    end do              ! On i = 1, n_sps

  end subroutine Get_GL_Slabs_Arrays

!=====================================================================

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module SLABS_SW_M

! $Log$
! Revision 2.67  2018/04/19 01:59:16  vsnyder
! Compute address of allocatable/deallocatable for tracking.  Remove USE
! statements for unused names.
!
! Revision 2.66  2016/10/24 22:16:38  vsnyder
! Make Slabs allocatable instead of a pointer
!
! Revision 2.65  2015/03/28 02:11:31  vsnyder
! Added stuff to trace allocate/deallocate addresses
!
! Revision 2.64  2014/09/05 21:27:29  vsnyder
! More complete and accurate allocate/deallocate size tracking
!
! Revision 2.63  2013/05/09 01:02:48  vsnyder
! Add useYi to dump
!
! Revision 2.62  2011/11/11 00:42:06  vsnyder
! Use IsExtinction array from Molecules module
!
! Revision 2.61  2009/06/23 18:26:11  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.60  2009/05/13 20:03:02  vsnyder
! Get constants from Constants, kinds from MLSKinds
!
! Revision 2.59  2008/10/03 16:30:48  livesey
! Added EXTINCTIONV2
!
! Revision 2.58  2008/05/20 00:23:50  vsnyder
! Receive Vel/C instead of Vel
!
! Revision 2.57  2008/02/29 01:57:37  vsnyder
! Use MLSKinds instead of MLSCommon
!
! Revision 2.56  2007/12/04 23:40:12  vsnyder
! Don't look for polarized if it's not allocated
!
! Revision 2.55  2007/05/23 22:41:49  vsnyder
! Change line data from struct of arrays to array of structs to improve
! cache locality.
!
! Revision 2.54  2006/12/04 21:17:28  vsnyder
! Reorganize FullForwardModel to use automatic arrays instead of allocating
! pointer arrays.  Requires testing for zero size instead of testing for
! associated in several subsidiary procedures.
!
! Revision 2.53  2006/09/01 00:59:45  vsnyder
! "Catalog" argument of AllocateSlabs needs TARGET attribute so that
! slabs(i)%catalog does not become undefined when AllocateOneSlabs returns
!
! Revision 2.52  2006/07/29 03:02:11  vsnyder
! Send callers module name from AllocateSlabs to AllocateOneSlabs
!
! Revision 2.51  2006/05/05 22:20:58  vsnyder
! Remove duplicate USE statements
!
! Revision 2.50  2006/03/25 00:27:46  vsnyder
! Avoid dividing by zeroes that weren't avoided in the previous revision
!
! Revision 2.49  2006/01/26 03:05:51  vsnyder
! Avoid dividing by zero
!
! Revision 2.48  2005/09/17 00:48:09  vsnyder
! Don't look at an array that might not be there, plus some cannonball polishing
!
! Revision 2.47  2005/09/03 01:21:33  vsnyder
! Spectral parameter offsets stuff
!
! Revision 2.46  2005/08/03 18:02:31  vsnyder
! Some spectroscopy derivative stuff, finish v0s modification
!
! Revision 2.45  2005/07/06 02:17:21  vsnyder
! Revisions for spectral parameter derivatives
!
! Revision 2.44  2005/06/09 02:34:16  vsnyder
! Move stuff from l2pc_pfa_structures to slabs_sw_m
!
! Revision 2.43  2005/03/29 01:58:17  vsnyder
! Make stuff pure
!
! Revision 2.42  2004/12/28 00:26:40  vsnyder
! Remove unreferenced declaration
!
! Revision 2.41  2004/12/13 20:55:36  vsnyder
! Make Slabs_Prep and Slabs_Prep_dT public.  Add Slabs_Prep_Struct.  Polish
! some TeXnicalities.  Revise Slabswint_dT and Slabswint_Lines_dT not to
! compute interference if |yi| < 1.0e-6.  Added UseYi argument to Slabs_Prep
! and Slabs_Prep_dT.  Use Slabs_Prep_Struct from Get_GL_Slabs_Arrays.
!
! Revision 2.40  2004/09/23 20:08:47  vsnyder
! Finish correcting divide by zero
!
! Revision 2.39  2004/09/16 22:16:21  vsnyder
! Avoid dividing by zero in Slabswint_dT also
!
! Revision 2.38  2004/09/16 20:24:23  vsnyder
! Avoid dividing by zero in Slabswing_Lines_dT
!
! Revision 2.37  2004/09/01 01:14:48  vsnyder
! Correct 'not_used_here' routine
!
! Revision 2.36  2004/08/05 20:59:32  vsnyder
! Don't do any calculations for gl_slabs with no lines
!
! Revision 2.35  2004/05/11 02:52:43  vsnyder
! Remove USE for Pi, which isn't referenced
!
! Revision 2.34  2004/04/24 02:26:54  vsnyder
! Move Voigt stuff to its own module
!
! Revision 2.33  2004/04/20 00:48:06  vsnyder
! Only use Taylor really close to the origin
!
! Revision 2.32  2004/04/19 21:02:19  vsnyder
! Use Taylor instead of CDrayson near the origin
!
! Revision 2.31  2004/04/17 00:37:00  vsnyder
! Analytic temperature derivatives
!
! Revision 2.30  2004/04/06 23:40:21  vsnyder
! Do slabs_prep where derivatives not requested in get_gl_slabs_arrays
! instead of doing nothing.
!
! Revision 2.29  2004/04/02 00:59:24  vsnyder
! Get catalog from slabs structure
!
! Revision 2.28  2004/03/30 02:25:08  vsnyder
! Comment out call to slabs_prep_dt until n1==0 etc are worked out
!
! Revision 2.27  2004/03/27 03:35:27  vsnyder
! Add pointer to catalog in slabs_struct.  Use it so as not to need to drag
! line centers and line widths around.  Write slabs_lines and slabswint_lines
! to get sum of beta over all lines; put slabs_struct instead of its components
! in the calling sequence.
!
! Revision 2.26  2004/03/20 03:17:44  vsnyder
! Steps along the way toward analytic temperature derivatives
!
! Revision 2.24  2003/07/09 22:46:24  vsnyder
! Futzing
!
! Revision 2.23  2003/07/08 00:09:18  vsnyder
! Inlined several functions
!
! Revision 2.22  2003/07/04 02:49:03  vsnyder
! Simplify interface to Get_GL_Slabs_Arrays
!
! Revision 2.21  2003/06/18 14:45:00  bill
! added subsetting feature for T-ders
!
! Revision 2.20  2003/06/13 21:28:20  bill
! fixed/improved some bugs with line shape derivative computations
!
! Revision 2.19  2003/05/19 19:58:07  vsnyder
! Remove USEs for unreferenced symbols, remove unused local variables
!
! Revision 2.18  2003/05/16 23:53:05  livesey
! Now uses molecule indices rather than spectags
!
! Revision 2.17  2003/05/09 19:25:31  vsnyder
! Expect T+DT instead of T and DT separately in Get_GL_Slabs_Arrays
!
! Revision 2.16  2003/05/05 23:00:26  livesey
! Merged in feb03 newfwm branch
!
! Revision 2.15.2.2  2003/02/27 00:57:20  vsnyder
! Cosmetic changes, get rid of declared but unused variables
!
! Revision 2.15.2.1  2003/02/13 17:29:26  bill
! abs coeff obeys detailed balance
!
! Revision 2.15  2003/01/16 19:41:42  jonathan
! tested version: in 1D case, compute only each element along the LOS path before tangent point, and fill otherside accordingly
!
! Revision 2.14  2003/01/16 19:08:43  jonathan
! testing
!
! Revision 2.13  2003/01/16 18:50:20  jonathan
! For 1D FWM compute first element along the LOS path then fill other grid points with value of the first grid point
!
! Revision 2.12  2003/01/16 18:04:12  jonathan
! add Do_1D option to get_gl_slabs_arrays
!
! Revision 2.11  2003/01/10 21:55:26  vsnyder
! Move SpeedOfLight from Geometry to Units
!
! Revision 2.10  2002/12/20 20:22:59  vsnyder
! Cosmetic changes
!
! Revision 2.9  2002/12/03 00:34:23  vsnyder
! Test optional argument presence before using them
!
! Revision 2.8  2002/10/08 17:08:06  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 2.7  2002/10/02 21:06:03  vsnyder
! Get SpeedOfLight from Geometry module
!
! Revision 2.6  2002/09/12 23:00:04  vsnyder
! Cosmetic changes, move USEs from module scope to procedure scope
!
! Revision 2.5  2002/08/05 17:51:15  jonathan
! debug
!
! Revision 2.4  2001/12/14 23:43:44  zvi
! Modification for Grouping concept
!
! Revision 2.3  2001/11/30 01:18:11  zvi
! Correcting a minor bug
!
! Revision 2.1  2001/10/17 22:01:00  zvi
! Eliminate computation of: ns
!
! Revision 2.0  2001/09/17 20:26:27  livesey
! New forward model
!
! Revision 1.4.2.4  2001/09/12 21:38:54  zvi
! Added CVS stuff
!
! Revision 1.4.2.3  2001/09/12 00:05:44  livesey
! Corrected sign of velocity correction
!
! Revision 1.4.2.2  2001/09/10 10:02:32  zvi
! Cleanup..comp_path_entities_m.f90
!
! Revision 1.1  2001/01/31 18:12:06  Z.Shippony
! Initial conversion to Fortran 90
