! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Min_Zeta_m

  implicit NONE
  private
  public :: Get_Min_Zeta

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

contains

  subroutine Get_Min_Zeta ( P_Basis, H_Ref, T_Ref, Zeta_T, P_Grid, Tan_Pt, H_t, &
    &                       Min_Zeta, Min_Phi, Min_Index )

    use GLNP, only: NGP1
    use MLSKinds, only: IP, RP

    real(rp), intent(in) :: P_Basis(:) ! Temperature basis phis
    real(rp), intent(in) :: H_Ref(:)   ! Tangent row of Metrics H_Ref
    real(rp), intent(in) :: T_Ref(:)   ! Tangent row of Metrics T_Ref
    real(rp), intent(in) :: Zeta_T     ! Zeta at the tangent
    real(rp), intent(in) :: P_Grid(:)  ! Phi on the fine path
    integer, intent(in) :: Tan_Pt      ! Tangent index in P_Grid
    real(rp), intent(in) :: H_T        ! H at the tangent

    real(rp), intent(out) :: Min_Zeta  ! Minimum zeta
    real(rp), intent(out) :: Min_Phi   ! Phi where minimum zeta occurs
    integer(ip), intent(out) :: Min_Index  ! Min_Phi is between
    !            P_Grid(min_index) and P_Grid(min_index+1) if min_index > 0,
    !            else there is no local zeta minimum other than within one
    !            point of the tangent point.

    !{ Get the minimum zeta along the line of sight.
    !
    ! Start with $\Delta\zeta = 14.8 \Delta H/T$.  Take a first-order expansion
    ! around a reference point $\phi_i$, with derivatives along constant-$\zeta$
    ! surfaces, giving
    !
    ! \begin{equation*}
    ! \Delta\zeta = \zeta - \zeta_i = 
    !  14.8 \frac{\Delta H + \frac{\partial H}{\partial \phi}{\Delta \phi}}
    !            {T_i + \frac{\partial T}{\partial \phi}{\Delta \phi}}
    ! \end{equation*}
    !
    ! Substitute $\hat\phi = \phi-\phi_t$, $\tilde\phi = \phi_t-\phi_i$,
    ! $\Delta \phi = \phi - \phi_i = \hat\phi+\tilde\phi$,
    ! $\Delta H = H - H_t = H_t(\sec \hat\phi - 1)$,
    !
    ! \begin{equation*}
    ! \frac1{H_t} \frac{\partial H}{\partial\phi} \approx G =
    ! \frac1{H_t} \frac{H_{i+1}-H_i}{\phi_{i+1}-\phi_i} \text{ and }
    ! \frac{\partial T}{\partial\phi} \approx \Gamma =
    ! \frac{T_{i+1}-T_i}{\phi_{i+1}-\phi_i} \text{ giving}
    ! \end{equation*}
    !
    ! \begin{equation*}
    ! \Delta \zeta \approx \zeta - \zeta_i = 14.8 H_t
    !  \frac{\sec \hat\phi - 1 + G(\hat\phi+\tilde\phi)}
    !       {T_i + \Gamma(\hat\phi+\tilde\phi)} \text{ from which}
    ! \end{equation*}
    !
    ! \begin{equation*}
    ! \frac{\text{d}\zeta}{\text{d}\phi} \approx
    ! \frac{14.8 H_t}{\left(T_i + \Gamma(\hat\phi+\tilde\phi)\right)^2}
    ! \left[ \left( \sec \hat\phi \tan \hat\phi + G \right)
    !             (T_i + \Gamma(\hat\phi+\tilde\phi)) -
    !        \Gamma\left(\sec \hat\phi - 1 + G(\hat\phi+\tilde\phi)\right)
    ! \right]\,.
    ! \end{equation*}
    !
    ! Write the second factor in terms of sines and cosines, take only
    ! the numerator and set it to zero, giving
    !
    ! \begin{equation*}
    ! (\sin\hat\phi + G \cos^2 \hat\phi)
    !  \left(T_i + \Gamma(\hat\phi+\tilde\phi)\right) -
    ! \Gamma \left( \cos \hat\phi -\cos^2 \hat\phi
    !  (1-G(\hat\phi+\tilde\phi))\right) = 0\,.
    ! \end{equation*}
    !
    ! Expand to second order in $\hat\phi$ giving
    !
    ! \begin{equation*}
    ! \left(\frac12 \Gamma - G T_i\right) \hat\phi^2 + 
    ! ( T_i + \Gamma \tilde\phi ) \hat\phi + G T_i = 0\,.
    ! \end{equation*}
    !
    ! Since we expect $|\frac12\Gamma-G T_i| << |T_i + \Gamma \tilde\phi|$,
    ! write the solution as
    !
    ! \begin{equation*}
    ! \hat\phi \approx \frac{-2 G T_i}
    !  {T_i + \Gamma \tilde \phi \pm
    !   \sqrt{(T_i+\Gamma \tilde \phi)^2 - 2 G T_i(\Gamma-2 G T_i)}}
    ! \end{equation*}
    !
    ! and take the positive sign.
    !
    ! See wvs-045.

    real(rp) :: D     ! tgp**2 - gt2*(gamma-gt2)
    real(rp) :: D_P_Basis ! P_basis(i+1) - P_Basis(i)
    real(rp) :: G ! Difference approximation to d(log H)/d(phi) for constant zeta
    real(rp) :: Gamma ! Difference approximation to dT/d(phi) for constant zeta
    real(rp) :: GT2   ! 2 * G * T_ref(i)
    integer :: I      ! Subscript and loop inductor
    real(rp) :: P     ! Phi_T - P_Basis(i), \tilde\phi in the TeXnicalities,
                      ! then \Delta\phi = \phi - \phi_i = \hat\phi + \tilde\phi.
    real(rp) :: Phi_T ! Phi at the tangent
    real(rp) :: Test_P_Min ! Candidate min_phi
    real(rp) :: Test_Z_Min ! Candidate min_zeta
    real(rp) :: TGP   ! T_ref(i) + Gamma * P

    logical, parameter :: Throw = .true. ! Throw out min zeta closer to the
    !                                       tangent than the first GL point

    min_index = -1 ! Assume there are no local zeta minima
    min_zeta = huge(min_zeta)

    ! Find the interval of P_Basis containing Phi_T
    phi_t = p_grid(tan_pt)
    do i = 0, size(p_basis)-1
      if ( phi_t <= p_basis(i+1) ) exit
    end do
    if ( i < 1 .or. i >= size(p_basis) ) return ! Phi_T not within range

    ! P_Basis(i) < Phi_T <= P_Basis(i+1) here.
    ! Look between either p_basis(i-1) and p_basis(i+1) or between
    ! p_basis(i) and p_basis(i+2) depending on whether phi_t is closer to
    ! p_basis(i) or p_basis(i+1), but don't look outside p_basis.
    if ( phi_t - p_basis(i) > p_basis(i+1) - phi_t ) i = i + 1
    do i = max(i-1,1), min(i,size(p_basis)-1)
      d_p_basis = p_basis(i+1) - p_basis(i)
      g = ( h_ref(i+1) - h_ref(i) ) / ( h_t * d_p_basis )
      gamma = ( t_ref(i+1) - t_ref(i) ) / d_p_basis
      p = phi_t - p_basis(i)

      !{ $\frac{\text{d}\zeta}{\text{d}\phi} = 0$ occurs for
      ! $\hat\phi = \phi-\phi_t \approx
      !  \frac{-2 G T_i}{T_i + \Gamma \hat\phi +
      !                  \sqrt{(T_i + \Gamma \tilde\phi)^2 - 
      !                        2 G T_i ( \Gamma - 2 G T_i)}}$
      ! where $\tilde\phi = \phi_t - \phi_i$ = {\tt P}.
      gt2 = g * t_ref(i) * 2.0
      tgp = t_ref(i) + gamma*p
      d = tgp**2 - gt2*(gamma-gt2)
      if ( d < 0.0 ) cycle ! Complex solutions

      test_p_min = -gt2 / ( tgp + sqrt(d) )
      p = test_p_min + p ! \Delta\phi = \phi - \phi_i = \hat\phi + \tilde\phi
      if ( p < 0.0 .or. p > d_p_basis ) cycle ! Outside range

      !{ {\tt test\_z\_min} = 
      ! $\Delta \zeta \approx \zeta - \zeta_i = 
      !  14.8 H_t \frac{\sec \hat\phi - 1 + G(\hat\phi+\tilde\phi)}
      !                {T_i + \Gamma(\hat\phi+\tilde\phi)}$
      test_z_min = 14.8 * h_t * (1.0/cos(test_p_min)-1.0 + G * p ) / &
        &                     (t_ref(i) + gamma*p )

      if ( test_z_min > 0.0 ) cycle ! zeta not < zeta_t, so no local minimum

      if ( test_z_min < min_zeta ) then
        min_zeta = test_z_min
        min_phi = test_p_min
        min_index = i
      end if

    end do

    ! Now add phi_t to get min_phi and zeta_t to get the min_zeta, and
    ! find where min_phi is in p_grid.
    if ( min_index > 0 ) then
      min_phi = min_phi + phi_t
      min_zeta = min_zeta + zeta_t
      if ( min_phi < phi_t ) then
        if ( min_phi >= p_grid(tan_pt-1) .and. throw ) then
          ! Too close to tangent point
          min_index = 0
          return
        end if
        min_index = tan_pt - 1
        do
          min_index = min_index - 1
          if ( min_index <= 0 ) return ! off the end of the path
          if ( p_grid(min_index) <= min_phi ) exit
        end do
      else
        if ( min_phi <= p_grid(tan_pt+ngp1+1) .and. throw ) then
          ! Too close to tangent point.  Remember that the tangent is duplicated
          min_index = 0
          return
        end if
        min_index = tan_pt + ngp1 + 1
        do
          min_index = min_index + 1
          if ( min_index > size(p_grid)-1 ) then ! off the end of the path
            min_index = 0
            return
          end if
          if ( min_phi <= p_grid(min_index+1) ) exit
        end do
      end if
    end if

  end subroutine Get_Min_Zeta

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Min_Zeta_m

! $Log$
! Revision 2.4  2013/05/21 23:54:17  vsnyder
! NG fine-grid points between coarse grid tangent points
!
! Revision 2.3  2009/06/23 18:26:10  pwagner
! Prevent Intel from optimizing ident string away
!
! Revision 2.2  2007/06/08 22:05:48  vsnyder
! More work on min zeta
!
! Revision 2.1  2006/12/21 01:33:31  vsnyder
! Initial commit
!
