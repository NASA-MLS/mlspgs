module Gemm_M

  implicit NONE
  public

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

  interface GEMM
    subroutine DGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, &
      &                BETA, C, LDC )
      character ::        TRANSA, TRANSB
      integer ::          M, N, K, LDA, LDB, LDC
      double precision :: ALPHA, BETA
      double precision :: A( LDA, * ), B( LDB, * ), C( LDC, * )
    end subroutine DGEMM
    subroutine SGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, &
      &                BETA, C, LDC )
      character ::        TRANSA, TRANSB
      integer ::          M, N, K, LDA, LDB, LDC
      real             :: ALPHA, BETA
      real             :: A( LDA, * ), B( LDB, * ), C( LDC, * )
    end subroutine SGEMM
  end interface

end module Gemm_M

! $Log$
