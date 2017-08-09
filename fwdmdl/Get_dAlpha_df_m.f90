! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Get_dAlpha_df_m
  
  implicit NONE
  private
  public :: Get_dAlpha_df, Get_d2Alpha_df2

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------
contains
!---------------------------------------------------------------------------

!-------------------------------------------------  Get_dAlpha_df  -----

  subroutine Get_dAlpha_df ( Sps_path, Beta_Path, dBeta_df, Grids_f, dAlpha_df )

    use LOAD_SPS_DATA_M, ONLY: GRIDS_T
    use MLSKinds, only: RP

    real(rp), intent(in) :: Sps_path(:,:)   ! Path mixing ratios.  Path X Sps.
    real(rp), intent(in) :: Beta_Path(:,:)  ! Path betas.  Path X Sps.
    real(rp), intent(in) :: dBeta_df(:,:)   ! Path beta derivatives w.r.t.
                                            ! mixing ratio.  Path X #Sps for
                                            ! which it exists.  The second
                                            ! subscripts here are
                                            ! grids_f%where_dBeta_df where
                                            ! those are nonzero, else the
                                            ! array is not referenced.
    type(grids_t), intent(in) :: Grids_f    ! Only lin_log and mol components
                                            ! per sps, not state-vector values.
    real(rp), intent(out) :: dAlpha_df(:,:) ! Path X Sps.

    !{ Compute $\frac{\partial \alpha(s)}{\partial f^k(s)}$ or
    !  $\hat{f}^k(s) \frac{\partial \alpha(s)}{\partial \hat{f}^k(s)}$.
    !  Ultimately, we want $\frac{\partial \alpha(s)}{\partial f^k_{lm}}
    !  = \frac{\partial \alpha(s)}{\partial f^k(s)}
    !  \frac{\partial f^k(s)}{\partial f^k_{lm}}$ or
    !  $\frac{\partial \alpha(s)}{\partial \hat{f}^k(s)}
    !  \frac{\partial \hat{f}^k(s)}{\partial f^k_{lm}}$.  There are four cases,
    !  depending upon whether linear interpolation is used from the solution
    !  grid to the path,
    !  $f^k(s) = \sum_{lm} \eta^k_{lm}(s) f^k_{lm}$, or logarithmic
    !  interpolation is used,
    !  $\hat{f}^k(s) = \exp\left( \sum_{lm} \eta^k_{lm}(s) \ln f^k_{lm} \right)$,
    !  and whether $\beta$ depends upon $f$.
    !  \begin{equation*}
    !  \begin{array}{l|l|l|l}
    !    \alpha(s) & \frac{\partial\alpha(s)}{\partial f^k(s)} \text{ or }
    !                \frac{\partial\alpha(s)}{\partial \hat{f}^k(s)}
    !              & \frac{\partial f^k(s)}{\partial f^k_{lm}} \text{ or }
    !                \frac{\partial \hat{f}^k(s)}{\partial f^k_{lm}}
    !              & \text{Computed here} \\[5pt]
    !    \hline &&&\\[-7pt]
    !    \sum_k f^k(s) \beta^k(s) & \beta^k(s) & \eta^k_{lm}(s)
    ! &    \frac{\partial\alpha(s)}{\partial f^k(s)}
    !\\  \sum_k f^k(s) \beta^k(s,f^k(s)) &
    !      \beta^k(s,f^k(s)) + f^k(s)
    !       \frac{\partial\beta^k(s,f^k(s))}{\partial f^k(s)}
    ! &    \eta^k_{lm}(s) &
    !      \frac{\partial\alpha(s)}{\partial f^k(s)}
    !\\  \sum_k \hat{f}^k(s) \beta^k(s) & \beta^k(s)
    ! &    \hat{f}^k(s) \frac{\eta^k_{lm}(s)}{f^k_{lm}}
    ! &    \hat{f}^k(s) \frac{\partial\alpha(s)}{\partial \hat{f}^k(s)}
    !\\  \sum_k \hat{f}^k(s) \beta^k(s,\hat{f}^k(s))
    ! &    \beta^k(s,\hat{f}^k(s)) + \hat{f}^k(s)
    !       \frac{\partial\beta^k(s,\hat{f}^k(s))}{\partial \hat{f}^k(s)}
    ! &    \hat{f}^k(s) \frac{\eta^k_{lm}(s)}{f^k_{lm}} &
    !      \hat{f}^k(s) \frac{\partial\alpha(s)}{\partial \hat{f}^k(s)}
    !\\
    !  \end{array}
    !  \end{equation*}
    !  All that's left to be done to what's computed here when we want
    !  $\frac{\partial \alpha(s)}{\partial f^k_{lm}}$ or
    !  $\frac{\partial \alpha(s)}{\partial \hat{f}^k_{lm}}$ is to multiply by
    !  $\eta^k_{lm}(s)$ or $\frac{\eta^k_{lm}(s)}{f^k_{lm}}$.  The reason
    !  for this separation is that all of what is computed here is eventually
    !  necessary to compute the desired derivatives, but in the several
    !  integrals $\int_{\zeta_j}^{\zeta_{j-1}}
    !     \frac{\partial \alpha(s)}{\partial f^k_{lm}} \,\text{d} s$
    !  that {\tt dRad_tran_df} needs to compute, $\eta^k_{lm}(s)$ is only
    !  nonzero for a few values of $s$.  For TScat computations,
    !  $\frac{\partial \omega_{ij}(s)}{\partial f^k_{lm}}$ is needed, and
    !  this depends upon $\frac{\partial \alpha(s)}{\partial f^k_{lm}}$.
    !  See wvs-095 and wvs-102.

    integer :: I_dBeta_df ! Second subscript of dBeta_df
    integer :: Sps_I      ! Subscript for molecules

    i_dBeta_df = 0

    do sps_i = 1, ubound(Grids_f%mol,1)

      i_dBeta_df = grids_f%where_dBeta_df(sps_i) ! Which column of dBeta_df?
      if ( i_dBeta_df /= 0 ) then
        ! beta depends upon f
        dAlpha_df(:,sps_i) = beta_path(:,sps_i) + &
                           & sps_path(:,sps_i) * dBeta_df(:,i_dBeta_df)
      else ! beta does not depend upon f
        dAlpha_df(:,sps_i) = beta_path(:,sps_i)
      end if
      if ( grids_f%lin_log(sps_i) ) &
        & dAlpha_df(:,sps_i) = sps_path(:,sps_i) * dAlpha_df(:,sps_i)

    end do ! sps_i

  end subroutine Get_dAlpha_df


!-----------------------------------------------  Get_d2Alpha_df2  -----

  subroutine Get_d2Alpha_df2 ( Sps_path, Beta_Path, dBeta_df, Grids_f, &
                             & d2Alpha_df2 )

    use LOAD_SPS_DATA_M, ONLY: GRIDS_T
    use MLSKinds, only: RP

    real(rp), intent(in) :: Sps_path(:,:)     ! Path mixing ratios.  Path X Sps.
    real(rp), intent(in) :: Beta_Path(:,:)    ! Path betas.  Path X Sps.
    real(rp), intent(in) :: dBeta_df(:,:)     ! Path beta derivatives w.r.t.
                                              ! mixing ratio.  Path X #Sps for
                                              ! which it exists.  The second
                                              ! subscripts here are
                                              ! grids_f%where_dBeta_df where
                                              ! those are nonzero, else the
                                              ! array is not referenced.
    type(grids_t), intent(in) :: Grids_f      ! For lin_log and mol components.
    real(rp), intent(out) :: d2Alpha_df2(:,:) ! Path X Sps.

    !{ Compute almost everything we need for
    !  $\frac{\partial^2 \alpha(s)}
    !        {\partial f^k_{lm}\partial f^k_{\tilde l \tilde m}}
    !  = \frac{\partial^2 \alpha(s)}{\partial f^k(s)^2}
    !  \frac{\partial f^k(s)}{\partial f^k_{lm}}
    !  \frac{\partial f^k(s)}{\partial f^k_{\tilde l \tilde m}}$ or
    !  $\frac{\partial^2 \alpha(s)}{\partial \hat{f}^k(s)^2}
    !  \frac{\partial \hat{f}^k(s)}{\partial f^k_{lm}}
    !  \frac{\partial f^k(s)}{\partial f^k_{\tilde l \tilde m}}$.
    !  $\frac{\partial^2 \alpha(s)}
    !        {\partial f^k_{lm}\partial f^{\tilde k}_{\tilde l \tilde m}} =0$
    !  for $k \neq \tilde k$.
    !  There are four cases, depending upon whether linear interpolation is
    !  used from the solution grid to the path,
    !  $f^k(s) = \sum_{lm} \eta^k_{lm}(s) f^k_{lm}$, or logarithmic
    !  interpolation is used,
    !  $\hat{f}^k(s) = \exp\left( \sum_{lm} \eta^k_{lm}(s) \ln f^k_{lm} \right)$,
    !  and whether $\beta$ depends upon $f$.  The terms involving
    !  $\frac{\partial^2 \beta}{\partial \hat f(s)^2}$ that are shown in wvs-102
    !  do not appear here because (so far) we have no species for which that
    !  is nonzero.
    !  \begin{equation*}
    !  \begin{array}{l|l|l}
    !  \alpha(s)
    ! & \text{What's computed here}
    ! & \text{Multiply by this to get}
    !   \frac{\partial^2 \alpha(s)}
    !        {\partial f^k_{lm}\partial f^k_{\tilde l \tilde m}}
    ! \\[10pt]
    !  \hline
    ! &
    !\\[-5pt]
    !  \sum_k f(s)\, \beta(s) 
    ! & 0
    ! & \text{nothing}
    !\\
    !  \sum_k f(s)\, \beta(s,f(s))
    ! & 2 \frac{\partial \beta}{\partial f(s)}
    ! & \eta^k_{lm}(s) \eta^k_{\tilde l \tilde m}(s)
    !\\
    !  \sum_k \hat{f}(s)\, \beta(s)
    ! &\beta \hat{f}(s)
    ! & \frac{\eta^k_{lm}(s)}{f^k_{lm} f^k_{\tilde l \tilde m}}
    !   (\eta^k_{\tilde l \tilde m}(s) - d_{lm, \tilde{l} \tilde{m}})
    !\\
    !  \sum_k \hat{f}(s)\, \beta(s,\hat{f}(s))
    ! &\hat{f}(s) \left( \beta +
    !       3 \hat{f}(s)\frac{\partial \beta}
    !                        {\partial \hat{f}(s)} \right)
    ! & \frac{\eta^k_{lm}(s)}{f^k_{lm} f^k_{\tilde l \tilde m}}
    !   (\eta^k_{\tilde l \tilde m}(s) - d_{lm, \tilde{l} \tilde{m}})
    !  \end{array}
    ! \end{equation*}
    ! \begin{eqnarray*}
    ! \mbox{where } \ \ d_{lm, \tilde{l} \tilde{m}} \ = \ 
    ! \begin{cases}  
    ! 1, \ \ \ l = \tilde{l} \ \mbox{and} \ m = \tilde{m} \\ 
    ! 0, \ \ \ \mbox{otherwise,} \end{cases} \ \ \ \ \mbox{is the Dirac delta function.}
    ! \end{eqnarray*}
    !  All that's left to be done to what's computed here when we want
    !  $\frac{\partial^2 \alpha(s)}
    !        {\partial f^k_{lm}\partial f^k_{\tilde l \tilde m}}$ or
    !  $\frac{\partial^2 \alpha(s)}
    !        {\partial \hat{f}^k_{lm}\partial \hat{f}^k_{\tilde l \tilde m}}$
    !  is to multiply by
    !  $\eta^k_{lm}(s)\eta^k_{\tilde l \tilde m}(s)$, and then divide by
    !  $f^k_{lm} f^k_{\tilde l \tilde m}$ in the last two cases.
    !  The reason for this separation is that all of what is computed here is
    !  eventually necessary to compute the desired derivatives, but in the
    !  several integrals $\int_{\zeta_j}^{\zeta_{j-1}}
    !     \frac{\partial^2 \alpha(s)}
    !          {\partial f^k_{lm}\partial f^k_{\tilde l \tilde m}} \,\text{d} s$
    !  that {\tt d2Rad_tran_df2} needs to compute, $\eta^k_{lm}(s)$ is only
    !  nonzero for a few values of $s$.

    integer :: I_dBeta_df ! Second subscript of dBeta_df
    integer :: Sps_I      ! Subscript for molecules

    i_dBeta_df = 0

    do sps_i = 1, ubound(Grids_f%mol,1)

      i_dBeta_df = grids_f%where_dBeta_df(sps_i) ! Which column of dBeta_df?
      select case ( merge(1,0,i_dBeta_df /= 0) + &
                  & merge(2,0,grids_f%lin_log(sps_i)) )
      case ( 0 ) ! f linear, beta doesn't depend upon f
        d2Alpha_df2(:,sps_i) = 0.0
      case ( 1 ) ! f linear, beta depends upon f
        d2Alpha_df2(:,sps_i) = 2.0 * dBeta_df(:,i_dBeta_df)
      case ( 2 ) ! f logarithmic, beta doesn't depend upon f
        d2Alpha_df2(:,sps_i) = beta_path(:,sps_i) * sps_path(:,sps_i)
      case ( 3 ) ! f logarithmic, beta depends upon f
        d2Alpha_df2(:,sps_i) = sps_path(:,sps_i) * ( beta_path(:,sps_i) + &
          & 3.0 * sps_path(:,sps_i) * dBeta_df(:,i_dBeta_df) )
      end select

    end do ! sps_i

  end subroutine Get_d2Alpha_df2

!----------------------------------------------------------------------
  logical function not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), save :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, not_used_here ! .mod files sometimes change if PRINT is added
  end function not_used_here

end module Get_dAlpha_df_m

! $Log$
! Revision 2.8  2017/08/09 20:21:06  vsnyder
! Spiff a comment
!
! Revision 2.7  2011/07/08 21:09:40  yanovsky
! Include subroutine Get_d2Alpha_df2 in a list of public subroutines
!
! Revision 2.6  2011/06/16 17:06:40  yanovsky
! Tex updates for second derivatives of alpha in Get_d2Alpha_df2 subroutine
!
! Revision 2.5  2011/03/17 00:00:25  vsnyder
! Simplify Get_d2Alpha_df2, more fiddling TeXnicalities
!
! Revision 2.4  2011/03/16 00:39:00  vsnyder
! Repair TeXnicalities, add second derivative
!
! Revision 2.3  2011/03/11 03:08:37  vsnyder
! Correct TeXnicalities
!
! Revision 2.2  2011/02/12 02:59:06  vsnyder
! Get column of dBeta_df from Grids_f, TeXnicalities
!
! Revision 2.1  2011/02/05 01:16:56  vsnyder
! Initial commit
!
