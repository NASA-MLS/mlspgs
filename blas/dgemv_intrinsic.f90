! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

! Since this is NOT a module, we can link it or ATLAS or LAPACK
! without causing a cascade.

! $RCSfile$

  subroutine DGEMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX, &
        &                BETA, Y, INCY )
    integer, parameter :: RK = kind(0.0d0)
    include "gemv.f9h"
  end subroutine DGEMV

! $Log$
