! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module Gemm_M

  implicit NONE
  public

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
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
  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module Gemm_M

! $Log$
! Revision 1.2  2001/11/13 23:51:20  vsnyder
! Fix CVS's Log comment
!
