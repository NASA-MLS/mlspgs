! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module Gemv_M

  implicit NONE
  public

!---------------------------- RCS Ident Info -------------------------------
  character (len=*), private, parameter :: IdParm = &
       "$Id$"
  character (len=len(idParm)), private :: Id = idParm
  character (len=*), private, parameter :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------

  interface GEMV
    subroutine DGEMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX, &
        &                BETA, Y, INCY )
      double precision   ALPHA, BETA
      integer            INCX, INCY, LDA, M, N
      character          TRANS
      double precision   A( LDA, * ), X( * ), Y( * )
    end subroutine DGEMV
    subroutine SGEMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX, &
        &                BETA, Y, INCY )
      real               ALPHA, BETA
      integer            INCX, INCY, LDA, M, N
      character          TRANS
      real               A( LDA, * ), X( * ), Y( * )
    end subroutine SGEMV
  end interface

end module Gemv_M

!$Log$
