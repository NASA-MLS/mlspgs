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
    module procedure DGEMM, SGEMM
  end interface

contains
  subroutine DGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, &
    &                BETA, C, LDC )
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    integer, parameter :: RK = kind(0.0d0)
    include "gemm.f9h"
  end subroutine DGEMM

  subroutine SGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, &
    &                BETA, C, LDC )
    use MLSMessageModule, only: MLSMessage, MLSMSG_Error
    integer, parameter :: RK = kind(0.0e0)
    include "gemm.f9h"
  end subroutine SGEMM
end module Gemm_M

!$Log: !
