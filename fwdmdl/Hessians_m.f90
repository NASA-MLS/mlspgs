! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Hessians_m

  implicit NONE
  private

  public :: d2Rad_Tran_dF2,  Get_d2Alpha_df2
  public :: Convolve_Other_Second_Deriv

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------
contains
!---------------------------------------------------------------------------

!-----------------------------------------------  Get_d2Alpha_df2  -----

  subroutine Get_d2Alpha_df2 ( Sps_path, Beta_Path, dBeta_df, Grids_f, &
                             & d2Alpha_df2 )

    use Load_SPS_Data_M, only: Grids_T
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

!------------------------------------------------  D2Rad_Tran_df2  -----
! This is the radiative transfer second derivative wrt mixing ratio model

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!                                                              !!!!!
!!!!! This has been revised to use a sparse representation of      !!!!!
!!!!! d_delta_df, but IT HAS NOT BEEN TESTED!                      !!!!!
!!!!!                                                              !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine D2Rad_Tran_df2 ( gl_inds, del_zeta, Grids_f, eta_fzp, do_gl,    &
                            & del_s, ref_cor, ds_dz_gw, inc_rad_path,        &
                            & d2Alpha_df2_c, d2Alpha_df2_f, i_start, tan_pt, &
                            & i_stop, d_delta_df, d2_delta_df2, d2rad_df2 )

    use Load_SPS_Data_m, only: Grids_t
    use MLSKinds, only: RP
    use Sparse_m, only: Sparse_t

! Inputs

    integer, intent(in) :: gl_inds(:)        ! Gauss-Legendre grid indices
    real(rp), intent(in) :: del_zeta(:)      ! path -log(P) differences on the
      !              main grid. This is for the whole coarse path, not just
      !              the part up to the black-out
    type (Grids_T), intent(in) :: Grids_f    ! All the coordinates
    class(sparse_t), intent(in) :: Eta_FZP(:) ! Interpolating coefficients
                                             ! from state vector to combined
                                             ! coarse & fine path for each sps
    logical, intent(in) :: do_gl(:)          ! A logical indicating where to
      !                                        do gl integrations
    real(rp), intent(in) :: del_s(:)        ! unrefracted path length.
    real(rp), intent(in) :: ref_cor(:)      ! refracted to unrefracted path
      !                                       length ratios.
    real(rp), intent(in) :: ds_dz_gw(:)     ! path length wrt zeta derivative *
      !              gw on the entire grid. Only the gl_inds part is used.
    real(rp), intent(in) :: inc_rad_path(:)  ! incremental radiance along the
                                             ! path.  t_script * tau.
    real(rp), intent(in) :: d2Alpha_df2_c(:,:) ! On the coarse path
    real(rp), intent(in) :: d2Alpha_df2_f(:,:) ! On the GL path
    integer, intent(in) :: I_start           ! path_start_index + 1
    integer, intent(in) :: tan_pt            ! Tangent point index
    integer, intent(in) :: I_stop            ! path stop index
    class (sparse_t), intent(in) :: d_delta_df(:) ! derivative of delta wrt
      !              mixing ratio state vector element.  One per SPS

! Outputs

    real(rp), intent(inout) :: d2_delta_df2(:,:,:) ! path x sve x sve.  Second 
      !               derivative of delta wrt mixing ratio state vector element.
    real(rp), intent(out) :: d2rad_df2(:,:)    ! second derivative of radiances wrt
                                               ! mixing ratio state vector element. (K)

! Internals

    real(rp) :: d_delta_df_q(size(del_s)), d_delta_df_r(size(del_s))
    integer :: First_q, First_r      ! Row number of first nonzero
    integer :: i_begin
    integer :: q, r                  ! state vector indices
    integer :: sparse_q_col, sparse_r_col
    integer :: sps_i, sps_j          ! species indices

! Begin code

    call get_all_d2_delta_df2( tan_pt, gl_inds, del_zeta, Grids_f, eta_fzp,    &
                             & do_gl, del_s, ref_cor, ds_dz_gw, d2Alpha_df2_c, &
                             & d2Alpha_df2_f, d2_delta_df2 )

    do sps_i = 1, ubound(Grids_f%l_z,1)

      do sps_j = 1, ubound(Grids_f%l_z,1)
    
        ! do q = 1,  Grids_f%l_v(sps_i)
        do q = Grids_f%l_v(sps_i-1)+1, Grids_f%l_v(sps_i)

          sparse_q_col = q - Grids_f%l_v(sps_i-1)
          if ( d_delta_df(sps_i)%cols(sparse_q_col) == 0 ) then
            d2rad_df2(q,:) = 0.0
            d2rad_df2(:,q) = 0.0
            cycle
          end if

          call d_delta_df(sps_i)%get_col ( sparse_q_col, d_delta_df_q )
          do first_q = 1, size(d_delta_df_q)
            if ( d_delta_df_q(first_q) /= 0 ) exit
          end do

          do r = Grids_f%l_v(sps_j-1)+1, Grids_f%l_v(sps_j)

            sparse_r_col = r - Grids_f%l_v(sps_j-1)
            call d_delta_df(sps_j)%get_col ( sparse_r_col, d_delta_df_r )
            do first_r = 1, size(d_delta_df_r)
              if ( d_delta_df_r(first_r) /= 0 ) exit
            end do

            i_begin = max(i_start,min(first_q,first_r,i_stop))

            call d2scrt_dx2 ( tan_pt, d_delta_df_q, d_delta_df_r, d2_delta_df2(:,q,r), &
                            & inc_rad_path, i_begin, i_stop, d2rad_df2(q,r) )

            call d_delta_df(sps_j)%clear_col ( sparse_r_col, d_delta_df_r )
          end do ! r
          call d_delta_df(sps_j)%clear_col ( sparse_q_col, d_delta_df_q )

        end do ! q

      end do ! sps_j

    end do ! sps_i

  end subroutine D2Rad_Tran_df2

!-------------------------------------------- Get_All_d2_Delta_df2 -----

  subroutine Get_All_d2_Delta_df2 ( tan_pt_c, gl_inds, del_zeta, Grids_f,  &
                              & eta_fzp, do_gl, del_s, ref_cor, ds_dz_gw,  &
                              & d2Alpha_df2_c, d2Alpha_df2_f, d2_delta_df2 )

    use Load_SPS_Data_M, only: Grids_T
    use MLSKinds, only: RP
    use Rad_Tran_m, only: Get_Do_Calc_Indexed, Get_Inds
    use Sparse_m, only: Sparse_t

! Inputs

    integer, intent(in) :: Tan_Pt_C          ! Index of tangent point in coarse
      !                                        path
    integer, intent(in) :: gl_inds(:)        ! Gauss-Legendre grid indices
    real(rp), intent(in) :: del_zeta(:)      ! path -log(P) differences on the
      !              main grid.  This is for the whole coarse path, not just
      !              the part up to the black-out
    type (Grids_T), intent(in) :: Grids_f    ! All the coordinates
    class (sparse_t), intent(in) :: Eta_FZP(:) ! Interpolating coefficients
                                             ! from state vector coordinates to
                                             ! combined coarse & fine path for
                                             ! each sps.
    logical, intent(in) :: do_gl(:)          ! A logical indicating where to
      !                                        do gl integrations
    real(rp), intent(in) :: ref_cor(:)       ! refracted to unrefracted path
      !                                        length ratios.
    real(rp), intent(in) :: del_s(:)         ! unrefracted path length.
    real(rp), intent(in) :: ds_dz_gw(:)      ! path length wrt zeta derivative *
      !              gw on the entire grid.  Only the gl_inds part is used.
    real(rp), intent(in) :: d2Alpha_df2_c(:,:)  ! On the coarse path
    real(rp), intent(in) :: d2Alpha_df2_f(:,:)  ! On the GL path

! Outputs

    real(rp), intent(inout) :: d2_delta_df2(:,:,:) ! path x sve x sve.  Second Derivative
      !              of delta wrt mixing ratio state vector elements. (K)
      !              Initially set to zero by caller.

! Internals

    integer :: n_inds_q, n_inds_r
    integer :: no_to_gl_q, no_to_gl_r
    integer :: sps_i, sps_j          ! species indices
    integer :: sps_n
    integer :: q, r                  ! state vector indices: sv_i, sv_j
    integer :: diracDelta            !   =1 if q=r;  =0 otherwise
    integer, target, dimension(1:size(del_s)) ::  all_inds_B_q,  all_inds_B_r
    integer, target, dimension(1:size(del_s)) :: more_inds_B_q, more_inds_B_r
    integer, pointer :: all_inds_q(:), all_inds_r(:)  ! all_inds => part of all_inds_B;
                                     ! Indices on GL grid for stuff
                                     ! used to make GL corrections
    integer, dimension(n_eta_rows(eta_fzp,.true.)) :: inds_q, inds_r ! Indices
                                     ! on coarse path where do_calc.
    integer, pointer :: more_inds_q(:), more_inds_r(:) ! more_inds => part of more_inds_B;
                                     ! Indices on the coarse path where GL
                                     ! corrections get applied.
    integer :: Sparse_q, Sparse_r    ! Column indices within Eta_FZP(sps_i)
                                     ! and Eat_FZP(sps_j)

    real(rp), dimension(n_eta_rows(eta_fzp,.true.)) :: Eta_FZP_q, Eta_FZP_r
    logical, dimension(n_eta_rows(eta_fzp,.true.)) :: NZ_FZP_q, NZ_FZP_r
    real(rp) :: singularity(1:size(del_s)) ! integrand on left edge of coarse
                                     ! grid panel -- singular at tangent pt.
    logical :: do_calc_q(1:size(del_s)) ! Flags on coarse path where do_calc_c
                                     ! or (do_gl and any corresponding
                                     ! nz_fzp_q).
    logical :: do_calc_r(1:size(del_s)) ! Flags on coarse path where do_calc_c
                                     ! or (do_gl and any corresponding
                                     ! nz_fzp_r).

! Begin code

    eta_fzp_q = 0
    eta_fzp_r = 0
    nz_fzp_q = .false.
    nz_fzp_r = .false.

    sps_n = ubound(Grids_f%l_z,1)

    do sps_i = 1, sps_n

      do sps_j = 1, sps_n

        do q = Grids_f%l_v(sps_i-1)+1, Grids_f%l_v(sps_i)

          if ( .not. Grids_f%deriv_flags(q) ) cycle

          sparse_q = q - Grids_f%l_v(sps_i-1)
          if ( eta_fzp(sps_i)%cols(sparse_q) /= 0 ) &
            & call eta_fzp(sps_i)%get_col ( sparse_q, eta_fzp_q, nz_fzp_q )

          ! find where the non zeros are along the path (for q)

          call get_do_calc_indexed ( size(do_gl), tan_pt_c, nz_fzp_q, &
            & gl_inds, do_gl, do_calc_q, n_inds_q, inds_q )

          if ( n_inds_q == 0 ) cycle

          do r = Grids_f%l_v(sps_j-1)+1, Grids_f%l_v(sps_j)

            sparse_r = r - Grids_f%l_v(sps_j-1)
            if ( eta_fzp(sps_j)%cols(sparse_r) == 0 ) cycle
            call eta_fzp(sps_j)%get_col ( sparse_r, eta_fzp_r, nz_fzp_r )

            ! Skip the masked derivatives, according to the l2cf inputs

            no_to_gl_q = count(do_gl(inds_q))

            all_inds_q =>  all_inds_B_q(1:no_to_gl_q)
            more_inds_q => more_inds_B_q(1:no_to_gl_q)

            ! see if anything needs to be gl-d (for q)
            if ( no_to_gl_q > 0 ) &
              & call get_inds ( do_gl, do_calc_q, more_inds_q, all_inds_q )

            ! find where the non zeros are along the path (for r)

            call get_do_calc_indexed ( size(do_gl), tan_pt_c, nz_fzp_r, &
              & gl_inds, do_gl, do_calc_r, n_inds_r, inds_r )

            if ( n_inds_r == 0 ) cycle

            no_to_gl_r = count(do_gl(inds_r))

            all_inds_r =>  all_inds_B_r(1:no_to_gl_r)
            more_inds_r => more_inds_B_r(1:no_to_gl_r)

            if ( no_to_gl_r > 0 ) &
            ! see if anything needs to be gl-d (for r)
              & call get_inds ( do_gl, do_calc_r, more_inds_r, all_inds_r )

     ! IGOR - Might not need to enter this subroutine for NON log basis sps, 
     ! since d2_delta_df2 for linear basis is zero.

            ! For molecules in logarithmic basis, calculate d2_delta_df2:

            if ( grids_f%lin_log(sps_i) ) then

              if( sps_i == sps_j ) then    ! otherwise, d2_delta_df2 = 0

                ! For same species, the following quantities should be the same
                ! for different sve:
                !   inds, all_inds, more_inds, sps.  Thus, only q quantities
                ! are passed.

                if ( q == r ) then
                  diracDelta = 1.0
                else
                  diracDelta = 0.0
                end if

                  if ( eta_fzp(sps_i)%cols(sparse_q) /= 0 .and. &
                   & eta_fzp(sps_j)%cols(sparse_r) /= 0 ) then ! Only nonzero columns
                  
                    call get_d2_delta_df2( diracDelta, inds_q, gl_inds,        &
                      & all_inds_q, more_inds_q, eta_fzp_q, eta_fzp_r,         &
                      & d2Alpha_df2_c(:,sps_i), d2Alpha_df2_f(:,sps_i), del_s, &
                      & del_zeta, ds_dz_gw, grids_f%values(q),                 &
                      & grids_f%values(r), singularity, d2_delta_df2(:,q,r),   &
                      & ref_cor )

                  end if
              end if   ! sps_i == sps_j

            end if   ! lin_log

            if ( eta_fzp(sps_j)%cols(sparse_r) /= 0 ) &
              & call eta_fzp(sps_i)%clear_col ( sparse_r, eta_fzp_r, nz_fzp_r )

          end do   ! r

          
          if ( eta_fzp(sps_j)%cols(sparse_r) /= 0 ) &
            & call eta_fzp(sps_j)%clear_col ( sparse_q, eta_fzp_q, nz_fzp_q )

        end do   ! q

      end do   ! sps_j

    end do   ! sps_i

  end subroutine Get_All_d2_Delta_df2

! .............................................  Get_d2_delta_df2  .....
  subroutine Get_d2_delta_df2 ( diracDelta, Inds, GL_Inds, All_inds, &
    & More_inds, eta_fzp_q, eta_fzp_r, d2Alpha_df2_path_c, d2Alpha_df2_path_f, &
    & Del_s, Del_Zeta, ds_dz_gw, Grids_v_q, Grids_v_r, Singularity,            &
    & d2_delta_df2, Ref_cor )

    ! Get d2_delta_df2 for the case of lin_log species for which beta
    ! does not depend upon mixing ratio.

    use GLNP, only: NG, NGP1
    use MLSKinds, only: RP

    integer, intent(in) :: diracDelta   !   =1 if q=r;  =0 otherwise
    integer, intent(in) :: Inds(:)      ! Indices on coarse path needing calc
    integer, intent(in) :: GL_Inds(:)   ! Gauss-Legendre grid indices
    integer, intent(in) :: All_inds(:)  ! Indices on GL grid for stuff
                                        ! used to make GL corrections
    integer, intent(in) :: More_inds(:) ! Indices on the coarse path where
                                        ! GL corrections get applied.
    real(rp), intent(in) :: eta_fzp_q(*), eta_fzp_r(*)  ! representation basis function.
    real(rp), intent(in) :: d2Alpha_df2_path_c(*)  ! d2Alpha_df2 on coarse grid.
    real(rp), intent(in) :: d2Alpha_df2_path_f(*)  ! d2Alpha_df2 on GL grid.
    real(rp), intent(in) :: Del_s(:)    ! unrefracted path length.
    real(rp), intent(in) :: Del_zeta(:) ! path -log(P) differences on the
                       !  main grid.  This is for the whole coarse path, not just
                       !  the part up to the black-out
    real(rp), intent(in) :: ds_dz_gw(:) ! ds/dh * dh/dz * GL weights
    real(rp), intent(in) :: Grids_v_q, Grids_v_r     ! Grids_f%values(sv_i),  
                                                     ! Grids_f%values(sv_j)
    real(rp), intent(out) :: singularity(:) ! integrand on left edge of coarse
                               ! grid panel -- singular at tangent pt.
                               ! Actually just work space we don't want
                               ! to allocate on every invocation.
    real(rp), intent(inout) :: d2_delta_df2(:) ! Second Derivative of delta w.r.t.
                               ! Sps_Path.  intent(inout) so the unreferenced
                               ! elements do not become undefined.
    real(rp), intent(in), optional :: ref_cor(:)  ! refracted to unrefracted path
                                                  !  length ratios.

    integer :: AA, GA, I, II, III

    do i = 1, size(inds)

      ii = inds(i)
      iii = ii*ngp1 - ng

      singularity(ii) = d2Alpha_df2_path_c(ii) * eta_fzp_q(iii) * (eta_fzp_r(iii) - diracDelta)
      d2_delta_df2(ii) = singularity(ii) * del_s(ii)

    end do ! i


    do i = 1, size(all_inds)

      aa = all_inds(i)
      ga = gl_inds(aa)
      ii = more_inds(i)

      d2_delta_df2(ii) = d2_delta_df2(ii) + &
        & del_zeta(ii) * &
        & sum( (eta_fzp_q(ga:ga+ng-1) * (eta_fzp_r(ga:ga+ng-1) - diracDelta) * &
             & d2Alpha_df2_path_f(aa:aa+ng-1) - singularity(ii)) &
             & * ds_dz_gw(ga:ga+ng-1) )

    end do

    ! Refraction correction
    if( present(ref_cor) ) d2_delta_df2(inds) = ref_cor(inds) * d2_delta_df2(inds)

    d2_delta_df2(inds) = d2_delta_df2(inds) * exp(-grids_v_q) * exp(-grids_v_r)

  end subroutine Get_d2_delta_df2

! ......................................  Get_d2_delta_df2_linlog  .....
  subroutine Get_d2_delta_df2_linlog ( diracDelta, Inds, GL_Inds, All_inds, &
    & More_inds, eta_fzp_q, eta_fzp_r, Sps_path, Beta_path_c, Beta_path_f,  &
    & Del_s, Del_Zeta, ds_dz_gw, Ref_cor, Grids_v_q, Grids_v_r, &
    & Singularity, d2_delta_df2 )

    ! Get d2_delta_df2 for the case of lin_log species for which beta
    ! does not depend upon mixing ratio.

    use GLNP, only: NG, NGP1
    use MLSKinds, only: RP

    integer, intent(in) :: diracDelta   !   =1 if q=r;  =0 otherwise
    integer, intent(in) :: Inds(:)   ! Indices on coarse path needing calc
    integer, intent(in) :: GL_Inds(:)   ! Gauss-Legendre grid indices
    integer, intent(in) :: All_inds(:)  ! Indices on GL grid for stuff
                                        ! used to make GL corrections
    integer, intent(in) :: More_inds(:) ! Indices on the coarse path where
                                        ! GL corrections get applied.
    real(rp), intent(in) :: eta_fzp_q(*), eta_fzp_r(*)  ! representation basis function.
    real(rp), intent(in) :: Sps_path(:) ! exp(Path mixing ratios)
    real(rp), intent(in) :: Beta_path_c(*)  ! cross section on coarse grid.
    real(rp), intent(in) :: Beta_path_f(*)  ! cross section on GL grid.
    real(rp), intent(in) :: Del_s(:)    ! unrefracted path length.
    real(rp), intent(in) :: Del_zeta(:) ! path -log(P) differences on the
      !              main grid.  This is for the whole coarse path, not just
      !              the part up to the black-out
    real(rp), intent(in) :: ds_dz_gw(:) ! ds/dh * dh/dz * GL weights
    real(rp), intent(in) :: ref_cor(:)  ! refracted to unrefracted path
                                        !  length ratios.
    real(rp), intent(in) :: Grids_v_q, Grids_v_r     ! Grids_f%values(sv_i)
    real(rp), intent(out) :: singularity(:) ! integrand on left edge of coarse
                               ! grid panel -- singular at tangent pt.
                               ! Actually just work space we don't want
                               ! to allocate on every invocation.
    real(rp), intent(inout) :: d2_delta_df2(:) ! Second Derivative of delta w.r.t.
                               ! Sps_Path.  intent(inout) so the unreferenced
                               ! elements do not become undefined.

    integer :: AA, GA, I, II, III

    do i = 1, size(inds)

      ii = inds(i)
      iii = ii*ngp1 - ng

      singularity(ii) = eta_fzp_q(iii) * (eta_fzp_r(iii) - diracDelta) * &
                      & sps_path(iii) * beta_path_c(ii)
      d2_delta_df2(ii) = singularity(ii) * del_s(ii)

    end do ! i


    do i = 1, size(all_inds)

      aa = all_inds(i)
      ga = gl_inds(aa)
      ii = more_inds(i)

      d2_delta_df2(ii) = d2_delta_df2(ii) + &
        & del_zeta(ii) * &
        & sum( (eta_fzp_q(ga:ga+ng-1) * (eta_fzp_r(ga:ga+ng-1) - diracDelta) * &
             & sps_path(ga:ga+ng-1) * beta_path_f(aa:aa+ng-1) - singularity(ii)) &
             & * ds_dz_gw(ga:ga+ng-1) )

    end do

    ! Refraction correction
    d2_delta_df2(inds) = ref_cor(inds) * d2_delta_df2(inds) * exp(-grids_v_q) * exp(-grids_v_r)

  end subroutine Get_d2_delta_df2_linlog

!----------------------------------------------------  D2SCRT_DX2  -----
! Compute the scalarized condensed radiative transfer second derivatives,
! without derivatives of the differential Temperatures w.r.t. species.

  subroutine D2SCRT_DX2 ( TAN_PT, D_DELTA_DX_q, D_DELTA_DX_r, D2_DELTA_DX2_qr, &
                        &  INC_RAD_PATH, I_START, I_END, D2RAD_DX2 )

    use MLSKinds, only: IP, RP

! Linear basis: \\
!{ $\frac{ \partial^2 I(x) }{ \partial f_q^k \partial f_{\tilde{q}}^{\tilde{k}} } = 
! \sum_{i=1}^{2N} \triangle B_i W_{i, q}^{k} W_{i, \tilde{q}}^{\tilde{k}} \, \tau_i$
! Logarithmic basis: \\

! Inputs


    integer, intent(in) :: Tan_pt               ! Tangent point index in inc_rad_path
    real(rp), intent(in) :: d_delta_dx_q(:)     ! path opacity derivatives wrt x_q
    real(rp), intent(in) :: d_delta_dx_r(:)     ! path opacity derivatives wrt x_r
    real(rp), intent(in) :: d2_delta_dx2_qr(:)  ! path opacity derivatives wrt x_q and x_r
    real(rp), intent(in) :: inc_rad_path(:)     ! incremental radiance along the
                                                ! path.  t_script * tau.

    integer(ip), intent(in) :: i_start    ! where non-zeros on the path begin
    integer(ip), intent(in) :: i_end      ! where non-zeros on path end

! Output

    real(rp), intent(out) :: d2rad_dx2      ! radiance second derivative wrt x

! internals

    integer(ip) :: i, n_path
    real(rp) :: wq
    real(rp) :: wr
    real(rp) :: dwqr

    d2rad_dx2 = 0.0_rp
    wq = 0.0_rp
    wr = 0.0_rp
    dwqr = 0.0_rp
    n_path = size(inc_rad_path)

!{ $W_{i,q} = \sum_{j=2}^i \frac{\partial \Delta \delta_{j \rightarrow j - 1}}
!                           {\partial x_q}$,
!  where the derivatives are given by d_delta_dx_q and d_delta_dx_r.

    do i = max(2,i_start), min(i_end,tan_pt)
      wq   = wq   + d_delta_dx_q(i)
      wr   = wr   + d_delta_dx_r(i)
      dwqr = dwqr + d2_delta_dx2_qr(i)
      d2rad_dx2 = d2rad_dx2 + inc_rad_path(i) * (wq * wr - dwqr)
    end do

    if ( i_end <= tan_pt ) return

! Tangent point or bounce off the earth's surface

    if ( i_start <= tan_pt ) d2rad_dx2 = d2rad_dx2 + inc_rad_path(tan_pt+1) * (wq * wr - dwqr)
    ! IGOR: Why wq,wr,dwqr are not recalculated here???

! Same as above, but don't add in the zero-emission layer at the 
! tangent point or earth's surface

    do i = max(i_start,tan_pt+2) , min(n_path,i_end)
      wq   = wq   + d_delta_dx_q(i-1)
      wr   = wr   + d_delta_dx_r(i-1)
      dwqr = dwqr + d2_delta_dx2_qr(i-1)
      d2rad_dx2 = d2rad_dx2 + inc_rad_path(i) * (wq * wr - dwqr)
    end do

  end subroutine D2SCRT_DX2

  ! -------------------------------------  Convolve_Other_Second_Deriv -----
  subroutine Convolve_Other_Second_Deriv ( Convolve_Support, MAF, Channel, &
             & SbRatio, Update, Radiance, Qtys, Grids_f, &
             & MIF_Times, DeadTime, d2I_df2, Hessian, RowFlags )

    use Convolve_All_m, only: LoadMatrixValue
    use Dump_0, only: Dump
    use ForwardModelConfig, only: QtyStuff_T
    use Fov_Convolve_m, only: Convolve_Support_T
    use HessianModule_0, only: Dump
    use HessianModule_1, only: Hessian_t
    use HighOutput, only: OutputNamedValue
    use Load_sps_data_m, only: Grids_T, Dump
    use MLSFillValues, only: isNaN
    use MLSKinds, only: R8, RP, RV
    use MatrixModule_1, only: FindBlock
    use MLSStringLists, only: switchDetail
    use MLSMessageModule, only: MLSMessage, MLSMsg_Warning
    use Output_m, only: ResumeOutput, SuspendOutput
    use Toggles, only: Switches
    use VectorsModule, only: VectorValue_T

    ! Inputs
    type(convolve_support_t), intent(in) :: Convolve_Support
    integer, intent(in) :: MAF
    integer, intent(in) :: Channel
    real(r8), intent(in) :: SbRatio    ! Sideband ratio
    logical, intent(in) :: Update      ! "add to Jacobian, don't overwrite"
    type (VectorValue_T), intent(in) :: Radiance ! Only to get some indices
    type(QtyStuff_T), intent(in) :: Qtys(:)
    type (Grids_T), intent(in) :: Grids_f
    real(rv), pointer :: MIF_Times(:,:) ! Disassociated if no scan average, q.v.
    real(rv), pointer :: DeadTime(:,:)  ! Disassociated if no scan average, q.v.
    real(rp), intent(in) :: d2I_df2(:,:,:) ! mixing ratio derivatives or any
    !                                    parameter for which a simple
    !                                    convolution will suffice

    ! Outputs
    type (Hessian_t), intent(inout) :: Hessian
    logical, intent(inout) :: rowFlags(:) ! Flag to calling code

    ! Local variables
    integer :: Row, Col1, Col2
    integer :: err
    integer :: JF_I, JF_J
    integer :: KI, KJ
    integer :: NFZ_I, NFZ_J ! # frequency elements times # z elements in species
    integer :: NoChans
    integer :: SPS_I, SPS_J ! species indices
    integer :: q, r         ! indices of final elements in state vectors
    real(r8) :: d2rad_df2_out( size(convolve_support%del_chi_out), &
                             & size(d2i_df2,dim=2), size(d2i_df2,dim=3))
    real(rv), pointer :: MIF_Times_for_MAF(:)

    nullify ( MIF_Times_for_MAF )
    if ( associated(MIF_Times) ) MIF_Times_for_MAF => MIF_Times(:,maf)
    noChans = Radiance%template%noChans

    ! do the convolution
    if ( any(isNaN(d2i_df2)) ) then
      call dump( d2i_df2, 'd2i_df2', options='-H' )
      call MLSMessage( MLSMSG_Warning, ModuleName // 'Convolve_Other_Second_Deriv', &
        & 'NaNs found in d2i_df2' )
    endif

    call fov_convolve_3d ( convolve_support, d2i_df2, MIF_Times_for_MAF, &
                         & DeadTime, grids_f%deriv_flags, d2rad_df2_out )

    ! load second derivatives into Hessian
    ! First, find index location in Hessian and set the derivative flags

    row = FindBlock( Hessian%row, radiance%index, maf )
    rowFlags(row) = .TRUE.

    if ( switchDetail( switches, 'hnan', options='-fc' ) > -1 ) then
      call outputNamedValue( 'MAF', MAF )
      call outputNamedValue( 'Channel', Channel )
      call outputNamedValue( 'row', row )
      call outputNamedValue( 'noChans', noChans )
      call outputNamedValue( 'size(qtys)', size(qtys) )
      if ( .not. update) then
        call dump( Grids_f, details=2 )
      endif
    endif

    do sps_i = 1, size(qtys)

      if ( .not. qtys(sps_i)%foundInFirst ) cycle

      q = grids_f%l_v(sps_i-1)
      nfz_i = (Grids_f%l_f(sps_i) - Grids_f%l_f(sps_i-1)) * &
            & (Grids_f%l_z(sps_i) - Grids_f%l_z(sps_i-1))

      do jf_i = Grids_f%windowStart(sps_i), Grids_f%windowfinish(sps_i)

        col1 = FindBlock ( Hessian%col, qtys(sps_i)%qty%index, jf_i)

        do ki = 1, nfz_i

          q = q + 1

          do sps_j = 1, size(qtys)

            if ( .not. qtys(sps_j)%foundInFirst ) cycle

            r = grids_f%l_v(sps_j-1)
            nfz_j = (Grids_f%l_f(sps_j) - Grids_f%l_f(sps_j-1)) * &
                  & (Grids_f%l_z(sps_j) - Grids_f%l_z(sps_j-1))

            do jf_j = Grids_f%windowStart(sps_j), Grids_f%windowfinish(sps_j)

              col2 = FindBlock ( Hessian%col, qtys(sps_j)%qty%index, jf_j)
        
              call getFullBlock_Hessian ( hessian, row, col1, col2, 'atmospheric' )

              do kj = 1, nfz_j
            
                r = r + 1

                ! load derivatives for this (zeta & phi) if needed:

                if ( Grids_f%deriv_flags(q) .and. Grids_f%deriv_flags(r) ) then

                  call loadMatrixValue ( d2rad_df2_out(:,q,r), &
                  & hessian%block(row,col1,col2)%values(channel::noChans,ki,kj), &
                  & sbRatio, update, ERR )

                  select case ( err ) 
                  case (1)
                    call MLSMessage( MLSMSG_Warning, ModuleName, &
                      & 'd2rad_df2_out bigger than hessian slice' )
                  case (2)
                    call MLSMessage( MLSMSG_Warning, ModuleName, &
                      & 'd2rad_df2_out smaller than hessian slice' )
                  case (4)
                    call MLSMessage( MLSMSG_Warning, ModuleName, &
                      & 'NaNs found in d2rad_df2_out' )
                  case (5)
                    call MLSMessage( MLSMSG_Warning, ModuleName, &
                      & 'd2rad_df2_out bigger than hessian slice; NaNs, too' )
                  case (6)
                    call MLSMessage( MLSMSG_Warning, ModuleName, &
                      & 'd2rad_df2_out smaller than hessian slice; NaNs, too' )
                  case default
                  end select
                endif
              end do
              if ( switchDetail( switches, 'hnan', options='-fc' ) > -1) then
                call suspendOutput
                call dump( hessian%block(row,col1,col2), details=-3 )
                call resumeOutput
              endif
   
            end do
        
          end do

        end do

      end do

    end do

  end subroutine Convolve_Other_Second_Deriv

  ! ---------------------------------------  GetFullBlock_Hessian  -----
  subroutine GetFullBlock_Hessian ( Hessian, Row, Col1, Col2, What )
    use HessianModule_0, only: H_ABSENT, H_FULL, RH
    use HessianModule_1, only: CREATEBLOCK, HESSIAN_T
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    type (Hessian_t), intent(inout) :: Hessian
    integer, intent(in) :: Row, Col1, Col2
    character(len=*), intent(in) :: What

    select case ( Hessian%block(row,col1,col2)%kind )
      case ( h_absent )
        ! If profiling shows we're spending too much time filling the block,
        ! we could send in a signal and only fill for the channels for which
        ! we are not computing results.
        call CreateBlock ( Hessian, row, col1, col2, h_full, inittuples=0, fill=0.0_rh )
      case ( h_full )
      case default
        call MLSMessage ( MLSMSG_Error, ModuleName, &
        & 'Wrong matrix block type for ' // what // ' second derivative matrix' )
    end select

  end subroutine GetFullBlock_Hessian

  ! --------------------------------------------  FOV_Convolve_3d  -----
  subroutine FOV_Convolve_3d ( Convolve_Support, d2I_df2, MIF_Times, DeadTime, &
    & deriv_flag, d2Rad_df2_out )

    use FOV_Convolve_m, only: Convolve_Support_t, FOV_Convolve_1d
    use MLSKinds, only: Rp, Rv

    ! Inputs

    type(convolve_support_t), intent(in) :: Convolve_Support
    real(rp), intent(in) :: d2i_df2(:,:,:) ! mixing ratio SECOND derivatives or any
    !                                 parameter where a simple convolution
    !                                 will suffice
    real(rv), pointer :: MIF_Times(:)   ! Disassociated if no scan average, q.v.
    real(rv), pointer :: DeadTime(:,:)  ! Disassociated if no scan average, q.v.
    logical, optional, intent(in) :: deriv_flag(:) ! Flag to indicate which of
    !                                 d2i_df2 are to be calculated.  Deafult true.

    ! outputs

    real(rp), intent(out) :: d2rad_df2_out(:,:,:) ! output radiance
    !                                 derivatives for input di_df.

    integer :: I, J, N_Coeffs

    ! nominally the mixing ratio SECOND derivatives but can be used for any
    ! quantity requiring a simple convolution.

    n_coeffs = size(d2i_df2,dim=2)

    do i = 1, n_coeffs
      do j = 1, n_coeffs
        if ( present(deriv_flag) ) then
          if ( (.not. deriv_flag(i)) .or. (.not. deriv_flag(j)) ) then
            d2rad_df2_out(:,i,j) = 0.0
            cycle
          end if
        end if

        call fov_convolve_1d ( convolve_support, d2i_df2(:,i,j), MIF_Times, DeadTime, d2rad_df2_out(:,i,j) )

      end do
    end do

  end subroutine FOV_Convolve_3d

!----------------------------------------------------  N_Eta_Rows  -----
  pure integer function N_Eta_Rows ( Eta_FZP, Need ) result ( N )
    use Sparse_m, only: Sparse_t
    class(sparse_t), intent(in) :: Eta_FZP(:) ! Interpolating coefficients
                                           ! from state vector to combined
                                           ! coarse & fine path for each sps
    logical, intent(in) :: Need            ! Need to do the computation
    n = 0
    if ( need ) then ! need space for stuff; Some stuff is needed only for TScat
      n = maxval(eta_fzp%nRows)
    end if
  end function N_Eta_Rows

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

end module Hessians_m

! $Log$
! Revision 2.1  2018/05/14 23:33:48  vsnyder
! Initial commit
!
