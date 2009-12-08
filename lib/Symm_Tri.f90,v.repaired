! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Symm_Tri

  ! Several routines for dealing with symmetric tridiagonal systems, or ones
  ! that have symmetric outliers in the corners.  See wvs-086.

  use MLSKinds, only: R4, R8

  implicit NONE

!  private

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

! Factor_Symm_Tri ( DIAG, H, D, U )
! Factor_Symm_Tri_Plus_1 ( DIAG, H, D, R, U, C )
! Solve_Factored_Symm_Tri ( H, D, U, Z [, B] )
! Solve_Factored_Symm_Tri_Plus_1 ( H, D, R, U, C, Z [, B] )
! Solve_Symm_Tri ( DIAG, H, Z [, B] )
! Solve_Symm_Tri_Plus_1 ( DIAG, H, Z [, B] )

  public :: Factor_Symm_Tri, Factor_Symm_Tri_Plus_1
  public :: Solve_Factored_Symm_Tri, Solve_Factored_Symm_Tri_Plus_1
  public :: Solve_Symm_Tri, Solve_Symm_Tri_Plus_1

  interface Factor_Symm_Tri
    module procedure Factor_Symm_Tri_r4, Factor_Symm_Tri_r8
    module procedure Factor_Symm_Tri_inPlace_r4, Factor_Symm_Tri_inPlace_r8
    ! Plus 1
    module procedure Factor_Symm_Tri_Plus_1_r4, Factor_Symm_Tri_Plus_1_r8
    ! Factor in place:
    module procedure Factor_Symm_Tri_Plus_1X_r4, Factor_Symm_Tri_Plus_1X_r8
  end interface

  interface Factor_Symm_Tri_Plus_1
    module procedure Factor_Symm_Tri_Plus_1_r4, Factor_Symm_Tri_Plus_1_r8
    ! Factor in place:
    module procedure Factor_Symm_Tri_Plus_1X_r4, Factor_Symm_Tri_Plus_1X_r8
  end interface

  interface Solve_Factored_Symm_Tri
    module procedure Solve_Factored_Symm_Tri_r4, Solve_Factored_Symm_Tri_r8
    module procedure Solve_Factored_Symm_Tri_a_r4, Solve_Factored_Symm_Tri_a_r8
    module procedure Solve_Factored_Symm_Tri_X_r4, Solve_Factored_Symm_Tri_X_r8
    module procedure Solve_Factored_Symm_Tri_Xa_r4, Solve_Factored_Symm_Tri_Xa_r8
    ! Plus 1
    module procedure Solve_Factored_Symm_Tri_Plus_1s
    module procedure Solve_Factored_Symm_Tri_Plus_1d
    module procedure Solve_Fact_Symm_Tri_Plus_1_a_s
    module procedure Solve_Fact_Symm_Tri_Plus_1_a_d
    module procedure Solve_Factored_Symm_Tri_Plus1Xs
    module procedure Solve_Factored_Symm_Tri_Plus1Xd
    module procedure Solve_Fact_Symm_Tri_Plus_1X_as
    module procedure Solve_Fact_Symm_Tri_Plus_1X_ad
  end interface

  interface Solve_Factored_Symm_Tri_Plus_1
    module procedure Solve_Factored_Symm_Tri_Plus_1s
    module procedure Solve_Factored_Symm_Tri_Plus_1d
    module procedure Solve_Fact_Symm_Tri_Plus_1_a_s
    module procedure Solve_Fact_Symm_Tri_Plus_1_a_d
    module procedure Solve_Factored_Symm_Tri_Plus1Xs
    module procedure Solve_Factored_Symm_Tri_Plus1Xd
    module procedure Solve_Fact_Symm_Tri_Plus_1X_as
    module procedure Solve_Fact_Symm_Tri_Plus_1X_ad
  end interface

  interface Solve_Symm_Tri
    module procedure Solve_Symm_Tri_r4, Solve_Symm_Tri_r8
    module procedure Solve_Symm_Tri_a_r4, Solve_Symm_Tri_a_r8
    module procedure Solve_Symm_Tri_X_r4, Solve_Symm_Tri_X_r8
    module procedure Solve_Symm_Tri_X_a_r4, Solve_Symm_Tri_X_a_r8
  end interface

  interface Solve_Symm_Tri_Plus_1
    module procedure Solve_Symm_Tri_Plus_1_r4, Solve_Symm_Tri_Plus_1_r8
    module procedure Solve_Symm_Tri_Plus_1_a_r4, Solve_Symm_Tri_Plus_1_a_r8
    module procedure Solve_Symm_Tri_Plus_1X_r4, Solve_Symm_Tri_Plus_1X_r8
  end interface

contains

! ----------------------------------------------  Factor_Symm_Tri  -----
  subroutine Factor_Symm_Tri_r4 ( DIAG, H, D, U )
    ! Factor a symmetric tridiagonal matrix

    ! There are at least two ways to LU factor a symmetric tridiagonal matrix.
    ! The one with ones on the diagonal of U (the Crout factorization) is
    ! simple to compute:

    ! The subdiagonal of L is H, the subdiagonal of the matrix.
    ! The diagonal elements of L follow the simple recurrence
    ! \beta_1 = d_1
    ! \beta_i = d_i - h_i^2 / \beta_{i-1}

    ! The diagonal elements of U are 1
    ! The i'th superdiagonal element, where i indexes the row, is
    ! h_i / \beta_i.

    ! The one with ones on the diagonal of L (the Doolittle factorization) is
    ! also simple to compute, but I had to choose only one of them, and I
    ! chose Crout.

    integer, parameter :: RK = R4

    include "Factor_Symm_Tri.f9h"

  end subroutine Factor_Symm_Tri_r4

  subroutine Factor_Symm_Tri_r8 ( DIAG, H, D, U )
    ! Factor a symmetric tridiagonal matrix

    ! There are at least two ways to LU factor a symmetric tridiagonal matrix.
    ! The one with ones on the diagonal of U (the Crout factorization) is
    ! simple to compute:

    ! The subdiagonal of L is H, the subdiagonal of the matrix.
    ! The diagonal elements of L follow the simple recurrence
    ! \beta_1 = d_1
    ! \beta_i = d_i - h_i^2 / \beta_{i-1}

    ! The diagonal elements of U are 1
    ! The i'th superdiagonal element, where i indexes the row, is
    ! h_i / \beta_i.

    ! The one with ones on the diagonal of L (the Doolittle factorization) is
    ! also simple to compute, but I had to choose only one of them, and I
    ! chose Crout.

    integer, parameter :: RK = R8

    include "Factor_Symm_Tri.f9h"

  end subroutine Factor_Symm_Tri_r8

! --------------------------------------  Factor_Symm_Tri_inPlace  -----
  subroutine Factor_Symm_Tri_inPlace_r4 ( D, H, U )
    ! Factor a symmetric tridiagonal matrix in place
    integer, parameter :: RK = R4

    include "Factor_Symm_Tri_inPlace.f9h"

  end subroutine Factor_Symm_Tri_inPlace_r4

  subroutine Factor_Symm_Tri_inPlace_r8 ( D, H, U )
    ! Factor a symmetric tridiagonal matrix in place
    integer, parameter :: RK = R8

    include "Factor_Symm_Tri_inPlace.f9h"

  end subroutine Factor_Symm_Tri_inPlace_r8

! ---------------------------------------  Factor_Symm_Tri_Plus_1  -----
  subroutine Factor_Symm_Tri_Plus_1_r4 ( DIAG, H, D, R, U, C )
    ! Factor a symmetric tridiagonal matrix augmented by the same
    ! nonzero element in the top right and bottom left corners (so it's
    ! still symmetric).

    integer, parameter :: RK = R4

    include "Factor_Symm_Tri_Plus_1.f9h"

  end subroutine Factor_Symm_Tri_Plus_1_r4

  subroutine Factor_Symm_Tri_Plus_1_r8 ( DIAG, H, D, R, U, C )
    ! Factor a symmetric tridiagonal matrix augmented by the same
    ! nonzero element in the top right and bottom left corners (so it's
    ! still symmetric).

    integer, parameter :: RK = R8

    include "Factor_Symm_Tri_Plus_1.f9h"

  end subroutine Factor_Symm_Tri_Plus_1_r8

! -------------------------------------  Factor_Symm_Tri_Plus_1X  -----
  subroutine Factor_Symm_Tri_Plus_1X_r4 ( D, H, R, U, C )
    ! Factor a symmetric tridiagonal matrix augmented by the same
    ! nonzero element in the top right and bottom left corners (so it's
    ! still symmetric), in place.

    integer, parameter :: RK = R4

    include "Factor_Symm_Tri_Plus_1X.f9h"

  end subroutine Factor_Symm_Tri_Plus_1X_r4

  subroutine Factor_Symm_Tri_Plus_1X_r8 ( D, H, R, U, C )
    ! Factor a symmetric tridiagonal matrix augmented by the same
    ! nonzero element in the top right and bottom left corners (so it's
    ! still symmetric), in place.

    integer, parameter :: RK = R8

    include "Factor_Symm_Tri_Plus_1X.f9h"

  end subroutine Factor_Symm_Tri_Plus_1X_r8

! --------------------------------------  Solve_Factored_Symm_Tri  -----
  subroutine Solve_Factored_Symm_Tri_r4 ( H, D, U, Z, B )

    ! Solve the system A b = L U b = z, where A has been factored as L U,
    ! and its factors are represented by H, D, and U.

    integer, parameter :: RK = R4

    include "Solve_Factored_Symm_Tri.f9h"

  end subroutine Solve_Factored_Symm_Tri_r4

  subroutine Solve_Factored_Symm_Tri_r8 ( H, D, U, Z, B )

    ! Solve the system A b = L U b = z, where A has been factored as L U,
    ! and its factors are represented by H, D, and U.

    integer, parameter :: RK = R8

    include "Solve_Factored_Symm_Tri.f9h"

  end subroutine Solve_Factored_Symm_Tri_r8

! --------------------------------------  Solve_Factored_Symm_Tri  -----
  subroutine Solve_Factored_Symm_Tri_a_r4 ( H, D, U, Z, B )

    ! Solve the system A b = L U b = z, where A has been factored as L U,
    ! and its factors are represented by H, D, and U.

    integer, parameter :: RK = R4

    include "Solve_Factored_Symm_Tri_a.f9h"

  end subroutine Solve_Factored_Symm_Tri_a_r4

  subroutine Solve_Factored_Symm_Tri_a_r8 ( H, D, U, Z, B )

    ! Solve the system A b = L U b = z, where A has been factored as L U,
    ! and its factors are represented by H, D, and U.

    integer, parameter :: RK = R8

    include "Solve_Factored_Symm_Tri_a.f9h"

  end subroutine Solve_Factored_Symm_Tri_a_r8


! --------------------------------------  Solve_Factored_Symm_Tri_X  ---
  subroutine Solve_Factored_Symm_Tri_X_r4 ( H, D, U, Z )

    ! Solve the system A b = L U b = z, where A has been factored as L U,
    ! and its factors are represented by H, D, and U.

    integer, parameter :: RK = R4

    include "Solve_Factored_Symm_Tri_X.f9h"

  end subroutine Solve_Factored_Symm_Tri_X_r4

  subroutine Solve_Factored_Symm_Tri_X_r8 ( H, D, U, Z )

    ! Solve the system A b = L U b = z, where A has been factored as L U,
    ! and its factors are represented by H, D, and U.

    integer, parameter :: RK = R8

    include "Solve_Factored_Symm_Tri_X.f9h"

  end subroutine Solve_Factored_Symm_Tri_X_r8

! --------------------------------------  Solve_Factored_Symm_Tri_Xa  ---
  subroutine Solve_Factored_Symm_Tri_Xa_r4 ( H, D, U, Z )

    ! Solve the system A b = L U b = z, where A has been factored as L U,
    ! and its factors are represented by H, D, and U.  Solve with
    ! multiple RHS vectors.  The solution is returned in Z.

    integer, parameter :: RK = R4

    include "Solve_Factored_Symm_Tri_Xa.f9h"

  end subroutine Solve_Factored_Symm_Tri_Xa_r4

  subroutine Solve_Factored_Symm_Tri_Xa_r8 ( H, D, U, Z )

    ! Solve the system A b = L U b = z, where A has been factored as L U,
    ! and its factors are represented by H, D, and U.  Solve with
    ! multiple RHS vectors.  The solution is returned in Z.

    integer, parameter :: RK = R8

    include "Solve_Factored_Symm_Tri_Xa.f9h"

  end subroutine Solve_Factored_Symm_Tri_Xa_r8

! -------------------------------  Solve_Factored_Symm_Tri_Plus_1  -----
  subroutine Solve_Factored_Symm_Tri_Plus_1s ( H, D, R, U, C, Z, B )

    ! Solve the system A b = L U b = z, where A has been factored as L U,
    ! and its factors are represented by H, D, R, U, and C.

    integer, parameter :: RK = kind(1.0e0)

    include "Solve_Factored_Symm_Tri_Plus_1.f9h"

  end subroutine Solve_Factored_Symm_Tri_Plus_1s

  subroutine Solve_Factored_Symm_Tri_Plus_1d ( H, D, R, U, C, Z, B )

    ! Solve the system A b = L U b = z, where A has been factored as L U,
    ! and its factors are represented by H, D, R, U, and C.

    integer, parameter :: RK = kind(1.0d0)

    include "Solve_Factored_Symm_Tri_Plus_1.f9h"

  end subroutine Solve_Factored_Symm_Tri_Plus_1d

! ---------------------------------  Solve_Fact_Symm_Tri_Plus_1_a  -----
  subroutine Solve_Fact_Symm_Tri_Plus_1_a_s ( H, D, R, U, C, Z, B )

    ! Solve the system A b = L U b = z, where A has been factored as L U,
    ! and its factors are represented by H, D, R, U, and C.  Solve with
    ! multiple RHS vectors.

    integer, parameter :: RK = kind(1.0e0)

    include "Solve_Fact_Symm_Tri_Plus_1_a.f9h"

  end subroutine Solve_Fact_Symm_Tri_Plus_1_a_s

  subroutine Solve_Fact_Symm_Tri_Plus_1_a_d ( H, D, R, U, C, Z, B )

    ! Solve the system A b = L U b = z, where A has been factored as L U,
    ! and its factors are represented by H, D, R, U, and C.  Solve with
    ! multiple RHS vectors.

    integer, parameter :: RK = kind(1.0d0)

    include "Solve_Fact_Symm_Tri_Plus_1_a.f9h"

  end subroutine Solve_Fact_Symm_Tri_Plus_1_a_d

! -------------------------------  Solve_Factored_Symm_Tri_Plus_1X  ----
  subroutine Solve_Factored_Symm_Tri_Plus1Xs ( H, D, R, U, C, Z )

    ! Solve the system A b = L U b = z, where A has been factored as L U,
    ! and its factors are represented by H, D, R, U, and C.  The solution
    ! is returned in Z.

    integer, parameter :: RK = kind(1.0e0)

    include "Solve_Factored_Symm_Tri_Plus1X.f9h"

  end subroutine Solve_Factored_Symm_Tri_Plus1Xs

  subroutine Solve_Factored_Symm_Tri_Plus1Xd ( H, D, R, U, C, Z )

    ! Solve the system A b = L U b = z, where A has been factored as L U,
    ! and its factors are represented by H, D, R, U, and C.  The solution
    ! is returned in Z.

    integer, parameter :: RK = kind(1.0d0)

    include "Solve_Factored_Symm_Tri_Plus1X.f9h"

  end subroutine Solve_Factored_Symm_Tri_Plus1Xd

! ---------------------------------  Solve_Fact_Symm_Tri_Plus_1X_a  ----
  subroutine Solve_Fact_Symm_Tri_Plus_1X_as ( H, D, R, U, C, Z )

    ! Solve the system A b = L U b = z, where A has been factored as L U,
    ! and its factors are represented by H, D, R, U, and C.  The solution
    ! is returned in Z.  Solve with multiple RHS vectors.

    integer, parameter :: RK = kind(1.0e0)

    include "Solve_Fact_Symm_Tri_Plus_1X_a.f9h"

  end subroutine Solve_Fact_Symm_Tri_Plus_1X_as

  subroutine Solve_Fact_Symm_Tri_Plus_1X_ad ( H, D, R, U, C, Z )

    ! Solve the system A b = L U b = z, where A has been factored as L U,
    ! and its factors are represented by H, D, R, U, and C.  The solution
    ! is returned in Z.  Solve with multiple RHS vectors.

    integer, parameter :: RK = kind(1.0d0)

    include "Solve_Fact_Symm_Tri_Plus_1X_a.f9h"

  end subroutine Solve_Fact_Symm_Tri_Plus_1X_ad

! -----------------------------------------------  Solve_Symm_Tri  -----
  subroutine Solve_Symm_Tri_R4 ( DIAG, H, Z, B )

    ! Solve A b = z all at once -- don't return the factors.
    ! DIAG(1:n) is the diagonal of the matrix
    ! H(1:n-1) is the super- and subdiagonal, H(n) is the corner
    ! Z(1:n) is the RHS
    ! B(1:n) is output

    integer, parameter :: RK = R4

    include "Solve_Symm_Tri.f9h"

  end subroutine Solve_Symm_Tri_r4

  subroutine Solve_Symm_Tri_R8 ( DIAG, H, Z, B )

    ! Solve A b = z all at once -- don't return the factors.
    ! DIAG(1:n) is the diagonal of the matrix
    ! H(1:n-1) is the super- and subdiagonal, H(n) is the corner
    ! Z(1:n) is the RHS
    ! B(1:n) is output

    integer, parameter :: RK = R8

    include "Solve_Symm_Tri.f9h"

  end subroutine Solve_Symm_Tri_r8

! ---------------------------------------------  Solve_Symm_Tri_a  -----
  subroutine Solve_Symm_Tri_a_r4 ( DIAG, H, Z, B )

    ! Solve A b = z all at once -- don't return the factors.
    ! DIAG(1:n) is the diagonal of the matrix
    ! H(1:n-1) is the super- and subdiagonal, H(n) is the corner
    ! Z(1:n) is the RHS
    ! B(1:n) is output

    integer, parameter :: RK = R4

    include "Solve_Symm_Tri_a.f9h"

  end subroutine Solve_Symm_Tri_a_r4

  subroutine Solve_Symm_Tri_a_r8 ( DIAG, H, Z, B )

    ! Solve A b = z all at once -- don't return the factors.
    ! DIAG(1:n) is the diagonal of the matrix
    ! H(1:n-1) is the super- and subdiagonal, H(n) is the corner
    ! Z(1:n) is the RHS
    ! B(1:n) is output

    integer, parameter :: RK = R8

    include "Solve_Symm_Tri_a.f9h"

  end subroutine Solve_Symm_Tri_a_r8

! ---------------------------------------------  Solve_Symm_Tri_X  -----
  subroutine Solve_Symm_Tri_X_r4 ( DIAG, H, Z )

    ! Solve A b = z all at once -- don't return the factors.
    !   The solution b is returned in Z.
    ! DIAG(1:n) is the diagonal of the matrix
    ! H(1:n-1) is the super- and subdiagonal, H(n) is the corner
    ! Z(1:n) is the RHS on input and the solution on output

    integer, parameter :: RK = R4

    include "Solve_Symm_Tri_X.f9h"

  end subroutine Solve_Symm_Tri_X_r4

  subroutine Solve_Symm_Tri_X_r8 ( DIAG, H, Z )

    ! Solve A b = z all at once -- don't return the factors.
    !   The solution b is returned in Z.
    ! DIAG(1:n) is the diagonal of the matrix
    ! H(1:n-1) is the super- and subdiagonal, H(n) is the corner
    ! Z(1:n) is the RHS on input and the solution on output

    integer, parameter :: RK = R8

    include "Solve_Symm_Tri_X.f9h"

  end subroutine Solve_Symm_Tri_X_r8

! -------------------------------------------  Solve_Symm_Tri_X_a  -----
  subroutine Solve_Symm_Tri_X_a_r4 ( DIAG, H, Z )

    ! Solve A b = z all at once -- don't return the factors.
    !   The solution b is returned in Z.
    ! DIAG(1:n) is the diagonal of the matrix
    ! H(1:n-1) is the super- and subdiagonal, H(n) is the corner
    ! Z(1:n) is the RHS on input and the solution on output

    integer, parameter :: RK = R4

    include "Solve_Symm_Tri_X_a.f9h"

  end subroutine Solve_Symm_Tri_X_a_r4

  subroutine Solve_Symm_Tri_X_a_r8 ( DIAG, H, Z )

    ! Solve A b = z all at once -- don't return the factors.
    !   The solution b is returned in Z.
    ! DIAG(1:n) is the diagonal of the matrix
    ! H(1:n-1) is the super- and subdiagonal, H(n) is the corner
    ! Z(1:n) is the RHS on input and the solution on output

    integer, parameter :: RK = R8

    include "Solve_Symm_Tri_X_a.f9h"

  end subroutine Solve_Symm_Tri_X_a_r8

! -----------------------------------------------  Solve_Symm_Tri_Plus_1  -----
  subroutine Solve_Symm_Tri_Plus_1_r4 ( DIAG, H, Z, B )

    ! Solve A b = z all at once -- don't return the factors.
    ! DIAG(1:n) is the diagonal of the matrix
    ! H(1:n-1) is the super- and subdiagonal, H(n) is the corner
    ! Z(1:n) is the RHS
    ! B(1:n) is output

    integer, parameter :: RK = R4

    include "Solve_Symm_Tri_Plus_1.f9h"

  end subroutine Solve_Symm_Tri_Plus_1_r4

  subroutine Solve_Symm_Tri_Plus_1_r8 ( DIAG, H, Z, B )

    ! Solve A b = z all at once -- don't return the factors.
    ! DIAG(1:n) is the diagonal of the matrix
    ! H(1:n-1) is the super- and subdiagonal, H(n) is the corner
    ! Z(1:n) is the RHS
    ! B(1:n) is output

    integer, parameter :: RK = R8

    include "Solve_Symm_Tri_Plus_1.f9h"

  end subroutine Solve_Symm_Tri_Plus_1_r8

! -----------------------------------------------  Solve_Symm_Tri_Plus_1_a  -----
  subroutine Solve_Symm_Tri_Plus_1_a_r4 ( DIAG, H, Z, B )

    ! Solve A b = z all at once -- don't return the factors.
    ! DIAG(1:n) is the diagonal of the matrix
    ! H(1:n-1) is the super- and subdiagonal, H(n) is the corner
    ! Z(1:n) is the RHS
    ! B(1:n) is output

    integer, parameter :: RK = R4

    include "Solve_Symm_Tri_Plus_1_a.f9h"

  end subroutine Solve_Symm_Tri_Plus_1_a_r4

  subroutine Solve_Symm_Tri_Plus_1_a_r8 ( DIAG, H, Z, B )

    ! Solve A b = z all at once -- don't return the factors.
    ! DIAG(1:n) is the diagonal of the matrix
    ! H(1:n-1) is the super- and subdiagonal, H(n) is the corner
    ! Z(1:n) is the RHS
    ! B(1:n) is output

    integer, parameter :: RK = R8

    include "Solve_Symm_Tri_Plus_1_a.f9h"

  end subroutine Solve_Symm_Tri_Plus_1_a_r8

! -----------------------------------------------  Solve_Symm_Tri_Plus_1X  -----
  subroutine Solve_Symm_Tri_Plus_1X_r4 ( DIAG, H, Z )

    ! Solve A b = z all at once -- don't return the factors.
    ! DIAG(1:n) is the diagonal of the matrix
    ! H(1:n-1) is the super- and subdiagonal, H(n) is the corner
    ! Z(1:n) is the RHS
    ! B(1:n) is output

    integer, parameter :: RK = R4

    include "Solve_Symm_Tri_Plus_1X.f9h"

  end subroutine Solve_Symm_Tri_Plus_1X_r4

  subroutine Solve_Symm_Tri_Plus_1X_r8 ( DIAG, H, Z )

    ! Solve A b = z all at once -- don't return the factors.
    ! DIAG(1:n) is the diagonal of the matrix
    ! H(1:n-1) is the super- and subdiagonal, H(n) is the corner
    ! Z(1:n) is the RHS
    ! B(1:n) is output

    integer, parameter :: RK = R8

    include "Solve_Symm_Tri_Plus_1X.f9h"

  end subroutine Solve_Symm_Tri_Plus_1X_r8

! -----------------------------------------------  Solve_Symm_Tri_Plus_1X_a  -----
  subroutine Solve_Symm_Tri_Plus_1X_a_r4 ( DIAG, H, Z )

    ! Solve A b = z all at once -- don't return the factors.
    ! DIAG(1:n) is the diagonal of the matrix
    ! H(1:n-1) is the super- and subdiagonal, H(n) is the corner
    ! Z(1:n) is the RHS
    ! B(1:n) is output

    integer, parameter :: RK = R4

    include "Solve_Symm_Tri_Plus_1X_a.f9h"

  end subroutine Solve_Symm_Tri_Plus_1X_a_r4

  subroutine Solve_Symm_Tri_Plus_1X_a_r8 ( DIAG, H, Z )

    ! Solve A b = z all at once -- don't return the factors.
    ! DIAG(1:n) is the diagonal of the matrix
    ! H(1:n-1) is the super- and subdiagonal, H(n) is the corner
    ! Z(1:n) is the RHS
    ! B(1:n) is output

    integer, parameter :: RK = R8

    include "Solve_Symm_Tri_Plus_1X_a.f9h"

  end subroutine Solve_Symm_Tri_Plus_1X_a_r8

  logical function not_used_here()
!---------------------------- RCS Ident Info -------------------------------
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), save :: Id = idParm
!---------------------------------------------------------------------------
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, id ! .mod files sometimes change if PRINT is added
  end function not_used_here

end module Symm_Tri

! $Log$
! Revision 2.3  2009/12/08 21:26:28  vsnyder
! Combined some generics
!
! Revision 2.2  2009/12/02 01:17:36  vsnyder
! Add in-place solvers
!
! Revision 2.1  2009/11/18 00:07:46  vsnyder
! Initial commit
!
