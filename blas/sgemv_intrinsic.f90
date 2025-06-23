! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

! Since this is NOT a module, we can link it or ATLAS or LAPACK
! without causing a cascade.

  subroutine SGEMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX, &
        &                BETA, Y, INCY )
    integer, parameter :: RKM = kind(0.0e0), RKV = kind(0.0e0)
    include "gemv.f9h"
  end subroutine SGEMV

  subroutine SDGEMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX, &
        &                BETA, Y, INCY )
    integer, parameter :: RKM = kind(0.0e0), RKV = kind(0.0d0)
    include "gemv.f9h"
  end subroutine SDGEMV

! $Log$
! Revision 1.4  2002/10/01 22:40:16  pwagner
! Removed RCS Ident Block; already included in .h file
!
! Revision 1.3  2002/10/01 20:59:21  bwknosp
! Added Id and RCS info
!
! Revision 1.2  2002/09/13 22:53:12  pwagner
! Exploits matmul flexibility in accepting mixed types
!
! Revision 1.1  2001/11/14 00:23:28  vsnyder
! Initial commit
!
