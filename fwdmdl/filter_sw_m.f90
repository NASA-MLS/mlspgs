!
module FILTER_SW_M
  use UNITS, only: filter_unit
  use MLSCommon, only: I4, R8
  use STRINGS, only: LEFTJ, SQZSTR, STRUPR, STRLWR
  implicit NONE

!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
       "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
       "$RCSfile$"
!---------------------------------------------------------------------------
contains

!--------------------------------------------------------------------

SUBROUTINE load_filter(ich,ier,fdir,prim,lf,ael,aer,apc,asc,fa,x1,xl, &
                       xr,xn,nda,nsa,jcos)

!  This routine loads the filter parameters into memory. It has to be called
!  whenever the channel is changed,to initialize all the needed parameters.
!  (Before calling the Filter function). Once the parameters are loaded,the
!  user can now compute the filter's response by calling the Filter function.

!   Input: Fdir - String.  Name of the directory containing: *.ftn  files
!          Ich  - Integer. Channel number
!         prim  - String. Must be 'p' for Primary,or 'i' for Image

!  Output: Ier - Error status flag. Any non-zero value means trouble.

INTEGER(i4), PARAMETER :: mdeg = 10
INTEGER(i4), PARAMETER :: maxc = 130

INTEGER(i4), INTENT(IN)  :: ich
Integer(i4), INTENT(OUT) :: nda,nsa,jcos,lf,ier

Real(r8), INTENT(OUT) :: fa,x1,xl,xr,xn
Real(r8), INTENT(OUT) :: aer(*),ael(*),asc(*),apc(*)
!
CHARACTER (LEN=*), INTENT(IN) :: fdir
CHARACTER (LEN=*), INTENT(IN) :: prim

Real(r8) :: area
Real(r8) :: xlhs(2),xrhs(2),cer(3),cel(3)
Real(r8) :: scoeff(maxc),pcoeff(mdeg)
INTEGER(i4) :: j,k,ns,ne,ndeg,io

CHARACTER (LEN=80) :: fln
CHARACTER (LEN=2) :: api

NAMELIST/in/area,xlhs,xrhs,ne,cel,cer,ndeg,pcoeff,ns,scoeff

! Begin code:

  ier = 1
  jcos = 1

  api = 'xx'
  fln = prim
  CALL leftj(fln)
  CALL strlwr(fln)
  IF(fln(1:1) == 'p') THEN
    api = '_p'
  ELSE IF(fln(1:1) == 'i') THEN
    api = '_i'
  ELSE
    PRINT *,'** Error: Ileagal option ',prim,'. Must be: p or: i'
    RETURN
  END IF

!  Load the filter namelist in

  fln(1:)=' '
  lf = LEN_TRIM(fdir)
  WRITE(fln,'(A,''FLTR'',I2.2,A,''.FTN'')') fdir(1:lf),ich,api
  k = LEN_TRIM(fln)
  CALL strlwr(fln)

  OPEN(filter_unit,FILE=fln(1:k),STATUS='OLD',IOSTAT=io)
  IF(io /= 0) THEN
    PRINT *,'** Couldn''T OPEN FILE: ',fln(1:k)
    CALL errmsg(' ',io)
    GOTO 99
  END IF

  READ(UNIT=filter_unit,nml=in,IOSTAT=io)
  CLOSE(filter_unit,IOSTAT=j)
  IF(io /= 0) THEN
    PRINT *,'** Error while reading Namelist in file: ',fln(1:k)
    CALL errmsg(' ',io)
    GOTO 99
  END IF

  ier = 0
  nsa = ns
  nda = ndeg

  asc(1:nsa) = scoeff(1:nsa)
  apc(1:nda) = pcoeff(1:nda)

  ael(1:3) = cel(1:3)
  aer(1:3) = cer(1:3)

  fa = area
  x1 = xlhs(1)
  xl = xlhs(2)
  xr = xrhs(1)
  xn = xrhs(2)

99 CLOSE(filter_unit,IOSTAT=j)

  RETURN
END SUBROUTINE load_filter

!--------------------------------------------------------------------

SUBROUTINE get_filter_param(ich,a1,an,sf,ier,fdir,prim,lf)

!  This routine acts exactly like Load_Filter but in addition it returns
!  the left-hand & right-hand limits of the filter range,RELATIVE TO THE
!  FILTER'S CENTER FREQUENCY !!

!   Input: Ich - Integer,Channel number
!          Fdir - String.  Name of the directory containing: *.ftn  files
!          prim - String. Must be 'p' for Primary,or 'i' for Image

!  Output: a1 - Left limit
!          an - Right limit
!          sf - The area under the filter
!          lf - Lenght of string: fdir
!          Ier - Error status flag. Any non-zero value means trouble.

REAL(r8), INTENT(OUT)         :: a1, an, sf
INTEGER(i4), INTENT(IN)       :: ich
INTEGER(i4), INTENT(OUT)      :: ier, lf

CHARACTER (LEN=*), INTENT(IN) :: fdir
CHARACTER (LEN=*), INTENT(IN) :: prim

INTEGER(i4), PARAMETER :: mdeg = 10
INTEGER(i4), PARAMETER :: maxc = 130

Integer(i4) :: nda,nsa,jcos
Real(r8)    :: fa,x1,xl,xr,xn,ael(3),aer(3),asc(maxc),apc(mdeg)

  Call load_filter(ich,ier,fdir,prim,lf,ael,aer,apc,asc,fa,x1,xl, &
                   xr,xn,nda,nsa,jcos)
  IF(ier /= 0) RETURN

  a1 = x1
  an = xn
  sf = fa

  RETURN
END SUBROUTINE get_filter_param

!---------------------------------------------------------------------

REAL(r8) FUNCTION filter(f,ich,a1,an,sf,ier,fdir,prim,lf)

!  This function computes the filter's response.

!   Input: f - A double precision frequency,RELATIVE TO THE CENTRAL
!              FREQUENCY FOR THAT CHANNEL !!!
!          Ich - Optional. Integer,Channel number
!          Fdir - Optional. String. Name of the directory containing:
!                 *.ftn  files
!          prim - Optional. String. Must be 'p' for Primary,or 'i' for Image

!  Output: filter - Double precision Filter's response for that frequency
!          a1 - Left limit
!          an - Right limit
!          sf - The area under the filter
!          lf - Lenght of string: fdir
!          Ier - Error status flag. Any non-zero value means trouble.

REAL(r8), INTENT(IN) :: f
REAL(r8), INTENT(OUT), Optional :: a1, an, sf

Integer(i4), Intent(IN), optional :: ich
Integer(i4), Intent(OUT), optional :: ier, lf

CHARACTER (LEN=*), INTENT(IN), optional :: fdir
CHARACTER (LEN=*), INTENT(IN), optional :: prim

INTEGER(i4), PARAMETER :: mdeg = 10
INTEGER(i4), PARAMETER :: maxc = 130

Real(r8), SAVE :: fa,x1,xl,xr,xn
Real(r8), SAVE :: ael(3),aer(3)
Real(r8), SAVE :: asc(maxc),apc(mdeg)

Integer(i4), SAVE :: nda,nsa,jcos

Integer(i4), SAVE :: fcheck = -1

Real(r8) :: q,z,zs,zp,vf

! Begin code:

  filter = 0.0_r8
  if(PRESENT(ich) .and. PRESENT(fdir) .and. PRESENT(prim)) then
    fcheck = ich
    Call load_filter(ich,ier,fdir,prim,lf,ael,aer,apc,asc, &
                     fa,x1,xl,xr,xn,nda,nsa,jcos)
    if(ier /= 0) Stop
    a1 = x1
    an = xn
    sf = fa
  else if(fcheck < 1) then
    Print *,'** Error in function: Filter **'
    Print *,'   Subroutine LOAD_FILTER was not called for channel: ',ich
    Print *,'   BEFORE calling the Filter function ... Program ABORT'
    STOP
  endif

  IF(f < xl) THEN                 ! Left hand side exponential
    z = f - x1
    q = Pol(z,ael,3)
    vf = EXP(q)
  ELSE IF(f > xr) THEN            ! Right hand side exponential
    z = f - xn
    q = Pol(z,aer,3)
    vf = EXP(q)
  ELSE                             ! Middle,polynomial + Cosine series
    zs = 0.0D0
    zp = Pol(f,apc,nda)
    IF(jcos > 0) zs = CosSum(f,xl,xr,asc,nsa)
    vf = zp + zs
  END IF
!
  filter = vf

  RETURN
END FUNCTION filter

!--------------------------------------------------------------------

REAL(r8) FUNCTION Pol(x,vec,n)

!  Evaluates a polynomial,given the coefficients,degree and x


INTEGER(i4), INTENT(IN) :: n
REAL(r8), INTENT(IN)  :: x, vec(*)

REAL(r8) :: s
INTEGER(i4) :: i

! Begin code:

  s = 0.0_r8
  DO i = n, 1, -1
    s = x * s + vec(i)
  END DO

  Pol = s

  RETURN
END FUNCTION Pol

!---------------------------------------------------------------------

REAL(r8) FUNCTION CosSum(f,xl,xr,scoeff,n)

!  Evaluates a cosine series expansion for given f

INTEGER(i4), INTENT(IN) :: n
REAL(r8), INTENT(IN) :: f, xl, xr, scoeff(*)

REAL(r8), Parameter :: Pi = 3.1415926535897932385D0     ! Pi

INTEGER(i4) :: j
REAL(r8) :: df, fj, sum

! Begin code:

  sum = 0.5_r8 * scoeff(1)
  df = Pi * (f - xl) / (xr - xl)

  DO j = 2, n
    fj = (j - 1) * df
    sum = sum + scoeff(j) * COS(fj)
  END DO

  CosSum = sum

  RETURN
END FUNCTION CosSum

end module FILTER_SW_M
! $Log$
! Revision 1.2  2000/07/06 00:11:44  zvi
!  This is the Freeze version of Jun/24/2000
!
! Revision 1.1 2000/06/09 00:08:14  Z.Shippony
! Initial conversion to Fortran 90
