! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

! Since this is NOT a module, we can link it or ATLAS or LAPACK
! without causing a cascade.



  subroutine DGEMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX, &
        &                BETA, Y, INCY )
    integer, parameter :: RKM = kind(0.0d0), RKV = kind(0.0d0)
    include "gemv.f9h"
  end subroutine DGEMV

  subroutine DSGEMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX, &
        &                BETA, Y, INCY )
    integer, parameter :: RKM = kind(0.0d0), RKV = kind(0.0e0)
    include "gemv.f9h"
  end subroutine DSGEMV

! $Log$
! Revision 1.3  2002/10/01 20:59:17  bwknosp
! Added Id and RCS info
!
! Revision 1.2  2002/09/13 22:53:12  pwagner
! Exploits matmul flexibility in accepting mixed types
!
! Revision 1.1  2001/11/14 00:23:28  vsnyder
! Initial commit
!
