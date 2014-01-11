! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Compute_Model_Plane_m

  ! This module contains procedures to compute the plane in which the
  ! models evaluate the radiative-transfer equation, where that
  ! plane intersects grids of the atmospheric state, and interpolation
  ! coefficients between lines-of-sight and those grids.

  implicit none
  private

  public :: Compute_Model_Plane


  !---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: ModuleName= &
    "$RCSfile$"
  private :: not_used_here 
  !---------------------------------------------------------------------------

contains ! ============= Public Procedures =============================

  ! ----------------------------------------  Compute_Model_Plane  -----
  subroutine Compute_Model_Plane ( FwdModelExtra, Config, MAF, Normal, &
    & InOrbitPlane )

    ! Compute the normal to the model plane if it's not the orbit plane

    use Constants, only: Deg2Rad
    use ForwardModelConfig, only: ForwardModelConfig_t
    use ForwardModelVectorTools, only: GetQuantityForForwardModel
    use Geometry, only: To_Cart
    use Init_Tables_Module, only: L_Azimuth, L_ScVelECR
    use MLSKinds, only: RV
    use MLSNumerics, only: Cross ! Cross product of two 3-vectors
    use QuantityTemplates, only: RT
    use VectorsModule, only: Vector_t, VectorValue_t

    type(vector_t), intent(in) :: FwdModelExtra
    type(forwardModelConfig_t), intent(in) :: Config
    integer, intent(in) :: MAF
    real(rv), intent(out) :: Normal(3)   ! Normal vector to model plane, XYZ
      ! where the X axis pierces the prime meridian, the Y axis pierces
      ! the +90 degree meridian, and the Z axis pierces the north pole.
    logical, intent(out) :: InOrbitPlane ! Model plane is orbit plane

    real(rv), parameter :: Tol = 0.001 ! If |sin(azimuth)| < tol, the
      ! model plane is in the orbit plane.

    real(rv) :: Azimuth ! AzimuthQuantity%values(1,1)
    type(vectorValue_t), pointer :: AzimuthQuantity  ! angle by which the
      ! model plane is rotated counterclockwise from the spacecraft
      ! velocity direction about the spacecraft position vector
    real(rv) :: Caz, Saz  ! Cos(azimuth), Sin(azimuth)
    real(rv) :: P(3)      ! Normal to orbit plane
    real(rt) :: S(3)      ! Vector to the SC
    real(rt) :: V(3)      ! Unit normal parallel to SC velocity
    type(vectorValue_t), pointer :: Velocity ! Use associate construct ???

    ! The model plane is assumed to be the orbit plane if either there
    ! is no azimuth quantity in the Extra vector, or |sin(azimuth)| < tol.

    inOrbitPlane = .true.
    azimuthQuantity => GetQuantityForForwardModel ( fwdModelExtra, &
                    & noError=.true., quantityType=l_azimuth, config=config )
    if ( .not. associated(azimuthQuantity) ) return
    azimuth = azimuthQuantity%values(1,1)
    saz = sin(azimuth*deg2rad)
    if ( abs(saz) < tol ) return

    caz = cos(azimuth*deg2rad)    
    inOrbitPlane = .false.

    !{ Given $\mathbf{S}$, the spacecraft position and $\mathbf{V}$, the
    !  spacecraft velocity, the vector $\mathbf{P} = \mathbf{S} \times
    !  \mathbf{V}$ is normal to the orbit plane.  The vector $\mathbf{T}$
    !  that is in the orbit plane and perpendicular to $\mathbf{S}$ and
    !  $\mathbf{P}$ is given by $\mathbf{T} = \mathbf{S} \times \mathbf{P}
    !  = ( \mathbf{S} \cdot \mathbf{V}) \mathbf{S} - ( \mathbf{S} \cdot
    !  \mathbf{S} ) \mathbf{V} = ( \mathbf{S} \cdot \mathbf{V}) \mathbf{S}
    !  -  \mathbf{V}$ (because $|\mathbf{S}| = 1$). ($\mathbf{T}$ would be
    !  $-\mathbf{V}$ if the orbit were circular.) We want a vector in the
    !  $\mathbf{P} - \mathbf{T}$ plane, at an angle $\alpha$ ({\tt
    !  Azimuth}) from $\mathbf{P}$.  This vector $\mathbf{N}$ ({\tt
    !  Normal}) is given by $\mathbf{N} = \sin(\alpha) \mathbf{T} +
    !  \cos(\alpha) \mathbf{P}$.  If $\alpha = 0$, $\mathbf{N} =
    !  \mathbf{P}$, the normal to the orbit plane.  If $\alpha =
    !  \frac\pi2$, $\mathbf{N} = \mathbf{T} \approx -\mathbf{V}$, i.e.,
    !  $\alpha$ is an anti-clockwise rotation, from the spacecraft velocity
    !  vector, about the spacecraft position vector.

    ! Get geolocation information from velocity quantity.
    ! This will generate an error message if there is no velocity quantity.

!a    associate ( &
    velocity => &
              & GetQuantityForForwardModel ( fwdModelExtra, noError=.false., &
              & quantityType=l_scVelECR, config=config ) ! )
      call to_cart ( (/ velocity%template%geodLat(config%referenceMIF,MAF), &
                     &  velocity%template%lon(config%referenceMIF,MAF), 0.0_rt /), s )
      v = velocity%value3(1:3,config%referenceMIF,MAF)   ! V
!a    end associate
    s = s / sqrt(dot_product(s,s)) ! Unit S
    v = v / sqrt(dot_product(v,v)) ! Unit V
    p = cross(s, v)                ! right-handed unit normal to the orbit plane
    normal = dot_product(s,v) - v  ! in orbit plane, the vector T above
    ! Rotate Normal about S to desired azimuth.  If Azimuth is zero, this
    ! just gives back the orbit plane normal.
    normal =  normal * saz + p * caz
  end subroutine Compute_Model_Plane

!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Compute_Model_Plane_m

! $Log$
! Revision 2.5  2014/01/11 02:23:28  vsnyder
! Get the MIF from Config everywhere instead of an undefined local variable
!
! Revision 2.4  2014/01/11 02:22:16  vsnyder
! Get the MIF from Config instead of an undefined local variable
!
! Revision 2.3  2013/08/16 02:35:34  vsnyder
! Get all geolocation from SCVelECR quantity
!
! Revision 2.2  2013/07/18 01:11:32  vsnyder
! Replace scVel with scVelECR
!
! Revision 2.1  2013/07/15 16:34:54  vsnyder
! Initial commit
!
