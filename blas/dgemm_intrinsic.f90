! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

! This is NOT a module so that we can link it in place of ATLAS or LAPACK
! without causing a cascade.

! $RCSfile$

  subroutine DGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, &
    &                BETA, C, LDC )
    integer, parameter :: RK = kind(0.0d0)
    include "gemm.f9h"
  end subroutine DGEMM

! $Log$
! Revision 1.1  2001/11/12 18:59:20  vsnyder
! Initial commit
!
