      SUBROUTINE DCHOL (A, NDA, N, B, S, TOL, IERR)
c>> 1996-03-30 DCHOL Krogh  Added external statement.
c>> 1994-10-20 DCHOL Krogh  Changes to use M77CON
c>> 1991-06-26 DCHOL Krogh  Initial MATH77 version.
c
c     C. L. Lawson, JPL, 1970 April 30
c     Extensively modified by F. T. Krogh, JPL, 1991, September 23.
c
c Solution of positive definite system using Cholesky decomposition.
c
c Suppose the overdetermined system to be solved in the sense of least
c squares is   C*X = D
c
c Suppose normal equations are formed and stored as follows.
c
c                    A = (C**T)*C
c                    B = (C**T)*D
c                    S = (D**T)*D
c
c Given B and the upper triangle of A this subroutine computes the
c solution vector X, storing it in place of B.
c
c On input S should either = (D**T)D, or be zero.  If the former, then
c on return it is the euclidean norm of (CX - D).  Else the zero value
c is returned.
c
c The upper triangle of A is replaced by the upper triangular Cholesky
c matrix, R, satisfying  (R**T)*R = A.
c
c The value of TOL should be 10**(-k-1) where k is the minimum number of
c significant digits in A and B.   If TOL is < the relative precision in
c the computer's arithmetic, it is effectively replaced by that relative
c precision.  If some R(I,I)**2 is < TOL * (corresponding diagonal of
c A), then IERR is set to the value of I for the algebraically smallest
c value of R(I,I).  If R(I,I)**2 is .le. zero, then everything in the
c I-th row of R is set to zero, as is the I-th component of the
c solution, and IERR is replaced by -IERR to flag that this occured.
c
c--D replaces "?": ?CHOL
c
c ********************* Variable definitions ***************************
c
c A      Input matrix, only upper triangle is used.  Replace by factor.
c B      Input right hand side, output solution.
c D1MACH Library routine to get characteristics of floating point
c  numbers.  D1MACH(4) = smallest number, x, such that 1.0 + x .ne. 1.0.
c G      Temporary for accumulating sums.
c GMIN   Smallest value seen for G when working on a diagonal entry.
c I      Index used to access matrix entries.
c IERR   Used to return status, see above.
c J      Index used to access matrix entries.
c K      Index used to access matrix entries.
c N      Order of matrix.
c NDA    Declared first dimension for A.
c S      Formal argument, see above.
c SM     Temporary variable used for accumulating a sum.
c SM1    Temporary variable used for accumulating a sum.
c TOL    Input value for the tolerance, see above.
c TOLI   Internal value for the tolerance.
c TOLM   Machine tolerance, lower limit for TOLI.
c TSQ    TOLI**2.
c ZERO   Parameter equal to zero.
c
c ******************** Variable declarations ***************************
c
      external         D1MACH
      integer          NDA, N, IERR
      double precision A(NDA, N), B(N), S, TOL
      integer          I, J, K
      double precision ZERO, SM, SM1, G, TSQ, GMIN, TOLI, TOLM, D1MACH
      parameter (ZERO = 0.D0)
      save TOLM
      data TOLM /ZERO/
c
c ******************** Start of executable code ************************
c
      if (N .le. 0) return
      if (TOLM .eq. ZERO) TOLM = D1MACH(4)
      TOLI = max(TOL, TOLM)
      TSQ = TOLI**2
      SM1 = ZERO
      GMIN = A(1, 1)
      IERR = 0
      DO 50 I = 1, N
         DO 30 J = I, N
            SM = ZERO
            DO 10 K = 1, I - 1
               SM = SM + A(K, I) * A(K, J)
   10       continue
            G = A(I, J) - SM
            if ( J .ne. I) then
               A(I, J) = G / A(I, I)
            else
               if ( G .lt. TSQ * abs(A(J, J)) ) then
                  if (G .le. GMIN) then
                     GMIN = G
                     IERR = J
                     if (G .le. ZERO) IERR = -J
                  end if
                  if (G .le. ZERO) then
                     do 20 K = J, N
                        A(I, K) = ZERO
   20                continue
                     B(I) = ZERO
                     go to 50
                  end if
               end if
               A(J, J) = sqrt(G)
            end if
   30    continue
c
c                        Solve next row of first lower triangular system
         SM = ZERO
         DO 40 K = 1, I - 1
            SM = SM + A(K, I) * B(K)
   40    continue
         B(I) = (B(I) - SM) / A(I, I)
         SM1 = SM1 + B(I) * B(I)
   50 continue
c                        Get the solution norm
      S = sqrt(max(ZERO, S - SM1))
c
c                        Solve the second lower triangular system.
      DO 80  I = N, 1, -1
         if (A(I, I) .ne. ZERO) then
            SM = ZERO
            DO 70 J = I+1, N
               SM = SM + A(I, J) * B(J)
   70       continue
            B(I) = (B(I) - SM) / A(I, I)
         end if
   80 continue
      end
