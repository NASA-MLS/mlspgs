! Copyright (c) 2002, California Institute of Technology.  ALL RIGHTS RESERVED.
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
  private :: not_used_here 
!---------------------------------------------------------------------------

! === (start of toc) ===
! gemv        product of a matrix and a vector
! === (end of toc) ===

! === (start of api) ===
! gemv ( char TRANS, int M, int N, 
!                  val1 ALPHA, val1 A(LDA,*), int LDA, 
!                  val2 X(*), int incx, val2 BETA, val Y(*), int incy )
!      where val1, val2 can be any of the types:
!      { real, double precision }
! === (end of api) ===
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
    module procedure  SDGEMV, DSGEMV
  end interface

  contains
    subroutine SDGEMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX, &
        &                BETA, Y, INCY )
    ! Matrix-vector product of mixed type
      real               ALPHA, BETA
      integer            INCX, INCY, LDA, M, N
      character          TRANS
      real               A( LDA, * )
      double precision   X( * ), Y( * )
      ! double precision, dimension(M, N) :: a_double
      
      if ( m <= 0 .or. n <= 0 ) return
      ! allocate(a_double(m,n))
      ! a_double = a(1:m, 1:n)
      call dgemv(TRANS, m, n, alpha*1.d0, dble(a(1:m,1:n)), m, x, incx, &
        & beta*1.d0, y, incy)
      ! call dgemv(TRANS, m, n, alpha*1.d0, a_double, m, x, incx, &
      ! deallocate(a_double)
    end subroutine SDGEMV
  
    subroutine DSGEMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX, &
        &                BETA, Y, INCY )
    ! Matrix-vector product of mixed type
      double precision               ALPHA, BETA
      integer            INCX, INCY, LDA, M, N
      character          TRANS
      double precision               A( LDA, * )
      real   X( * ), Y( * )
      double precision, dimension(M) :: Y_DOUBLE
      ! double precision, dimension(N) :: X_DOUBLE
      
      if ( m <= 0 .or. n <= 0 ) return
      ! allocate(y_double(m), x_double(n))
      ! x_double = x(1:(n-1)*abs(incx)+1:incx)
      y_double = y(1:(m-1)*abs(incy)+1:incy)
      call dgemv(TRANS, m, n, alpha, a, m, dble(x(1:(n-1)*abs(incx)+1:incx)), incx, &
        & beta, y_double, 1)
      ! call dgemv(TRANS, m, n, alpha, a, m, x_double, 1, &
      y(1:(m-1)*abs(incy)+1:incy) = y_double
      ! deallocate(y_double, x_double)
    end subroutine DSGEMV
  
  logical function not_used_here()
    not_used_here = (id(1:1) == ModuleName(1:1))
  end function not_used_here

end module Gemv_M

!$Log$
!Revision 1.4  2002/09/18 22:45:34  pwagner
!Replaced internal arrays a-double, x_double with intrinsic dble function
!
!Revision 1.3  2002/09/13 22:52:13  pwagner
!Replaced pointers with automatic arrays
!
!Revision 1.2  2002/09/13 18:06:10  pwagner
!Added mixed type gemv: sdgemv and dsgemv when matrix and vector precisions differ
!
!Revision 1.1  2001/11/14 00:23:28  vsnyder
!Initial commit
!
