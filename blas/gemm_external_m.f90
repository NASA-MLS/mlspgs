! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

module Gemm_M

  implicit NONE
  public

!---------------------------- RCS Module Info ------------------------------
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
  private :: not_used_here 
!---------------------------------------------------------------------------

! === (start of toc) ===
! gemm        product of two matrices
! === (end of toc) ===

! === (start of api) ===
! gemm ( char TRANSA, char TRANSB, int M, int N, int K, val ALPHA, 
!                  val A(LDA,*), int LDA, val B(LDB,*), int LDB, 
!                  val BETA, val C(LDC,*), int LDC )
!      where val can be  one of the types:
!      { real, double precision }
! === (end of api) ===
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

contains 
!--------------------------- end bloc --------------------------------------
  logical function not_used_here()
  character (len=*), parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)) :: Id = idParm
    not_used_here = (id(1:1) == ModuleName(1:1))
    print *, Id ! .mod files sometimes change if PRINT is added
  end function not_used_here
!---------------------------------------------------------------------------

end module Gemm_M

! $Log$
! Revision 1.4  2005/06/22 19:45:19  pwagner
! Reworded Copyright statement, moved rcs id
!
! Revision 1.3  2002/10/10 23:48:57  pwagner
! Added idents to survive zealous Lahey optimizer
!
! Revision 1.2  2001/11/13 23:51:20  vsnyder
! Fix CVS's Log comment
!
