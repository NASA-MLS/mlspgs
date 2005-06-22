! Copyright 2005, by the California Institute of Technology. ALL
! RIGHTS RESERVED. United States Government Sponsorship acknowledged. Any
! commercial use must be negotiated with the Office of Technology Transfer
! at the California Institute of Technology.

! This software may be subject to U.S. export control laws. By accepting this
! software, the user agrees to comply with all applicable U.S. export laws and
! regulations. User has the responsibility to obtain export licenses, or other
! export authority as may be required before exporting such information to
! foreign countries or providing access to foreign persons.

! This is NOT a module so that we can link it in place of ATLAS or LAPACK
! without causing a cascade.


  subroutine DGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, &
    &                BETA, C, LDC )
    integer, parameter :: RK = kind(0.0d0)
    include "gemm.f9h"
  end subroutine DGEMM

! $Log$
! Revision 1.4  2002/10/01 22:40:15  pwagner
! Removed RCS Ident Block; already included in .h file
!
! Revision 1.3  2002/10/01 20:59:11  bwknosp
! Added Id and RCS info
!
! Revision 1.2  2001/11/13 21:16:02  vsnyder
! Get rid of dependence on MLSMessage
!
! Revision 1.1  2001/11/12 18:59:20  vsnyder
! Initial commit
!
