module PTG_FRQ_LOAD_M
  use MLSCommon, only: I4, R8
  use L2_TEST_STRUCTURES_M
  use D_HUNT_M, only: HUNT
  implicit NONE
!---------------------------- RCS Ident Info -------------------------------
  CHARACTER (LEN=256) :: Id = &
    "$Id$"
  CHARACTER (LEN=*), PARAMETER :: ModuleName= &
    "$RCSfile$"
!---------------------------------------------------------------------------
contains
!---------------------------------------------------------------------------
! ===========================================     ptg_frq_load =====
! This subprogram loads the pointing frequency gridding by band
!
  SUBROUTINE ptg_frq_load(FMC, FMI, Ier)
!---------------------------------------------------------------------------

Type(fwd_mdl_config), INTENT(IN OUT) :: FMC
Type(fwd_mdl_info), INTENT(IN OUT) :: FMI
Integer(i4), INTENT(OUT) :: ier
!
! ---- Local variables ---------
!
Integer(i4) :: i, k, kk, jp, io, l

Real(r8) :: q, r, dummy(200)
!
Character (LEN=80) :: Fnd, Line
!
! Load the pointing vs. frequencies database for the given band
! (needed for frequency averaging)
!
  Ier = 0
  Fnd(1:) = ' '
  Fnd = FMC%B
!
  k = size(FMI%tan_press)
  ALLOCATE(FMI%no_ptg_frq(k),STAT=ier)
  IF(ier /= 0) ALLOCATE(FMI%ptg_frq_grid(k),STAT=ier)
  IF(ier /= 0) THEN
    PRINT *,'** ALLOCATION error for no_ptg_frq or ptg_frq_grid ..'
    PRINT *,'   Allocation STAT =',ier
    Return
  ENDIF

  kk = -1
  FMI%no_ptg_frq(1:k) = 0
!
  CLOSE(32,iostat=i)
  OPEN(32,file=Fnd,action='READ',status='OLD',iostat=io)
  if(io /= 0) goto 10
!
! First entry in the file is the 'Band' frequency. All the rest are
! relative to this (center) frequency for this band
!
  Read(32,*,iostat=io) q
  if(io /= 0) goto 10
!
  DO
!
    Read(32,*,iostat=io) r, jp
    if(io > 0) goto 10
    if(io /= 0) EXIT

    k = -1
    Call Hunt(r,FMI%tan_press,size(FMI%tan_press),k,i)
    IF(ABS(r-FMI%tan_press(i)) < ABS(r-FMI%tan_press(k))) k = i
!
    DEALLOCATE(FMI%ptg_frq_grid(k)%values,STAT=i)

    if(FMC%do_frqavg) then
      FMI%no_ptg_frq(k) = jp
      ALLOCATE(FMI%ptg_frq_grid(k)%values(jp),STAT=i)
    else
      FMI%no_ptg_frq(k) = 1
      ALLOCATE(FMI%ptg_frq_grid(k)%values(2),STAT=i)
    endif

    IF(i /= 0) THEN
      ier = i
      PRINT *,'** Error: ALLOCATION error for ptg_frq_grid ..'
      PRINT *,'   tan_hts index:',k,' STAT =',ier
      do l = 1, k
        DEALLOCATE(FMI%ptg_frq_grid(l)%values,STAT=i)
      end do
      goto 99
    ENDIF

    Read(32,*,iostat=io) (dummy(i),i=1,jp)
    if(io /= 0) goto 10
    if(kk < 0) kk = k

    if(FMC%Zfrq > 0.0) dummy(1) = FMC%Zfrq - q
!
! Add 'band' frequency to ptg_frq_grid to convert to absolute grid
!
    jp = FMI%no_ptg_frq(k)
    FMI%ptg_frq_grid(k)%values(1:jp) = dummy(1:jp) + q
!
  END DO
!
  if(kk > 1) then
    jp = FMI%no_ptg_frq(kk)
    do k = 1, kk-1
      DEALLOCATE(FMI%ptg_frq_grid(k)%values,STAT=i)
      ALLOCATE(FMI%ptg_frq_grid(k)%values(jp),STAT=i)
      IF(i /= 0) THEN
        ier = i
        PRINT *,'** Error: ALLOCATION error for ptg_frq_grid ..'
        PRINT *,'   tan_hts index:',k,' STAT =',ier
        do l = 1, k
          DEALLOCATE(FMI%ptg_frq_grid(l)%values,STAT=i)
        end do
        goto 99
      ENDIF
      FMI%no_ptg_frq(k) = jp
      FMI%ptg_frq_grid(k)%values(1:jp) = &
     &           FMI%ptg_frq_grid(kk)%values(1:jp)
    end do
  endif
!
  do k = 1, size(FMI%tan_press)
    if(FMI%no_ptg_frq(k) < 1) then
      i = k
      do while(i > 1 .and. FMI%no_ptg_frq(i) < 1)
        i = i - 1
      end do
      jp = FMI%no_ptg_frq(i)
      if(jp < 1) CYCLE
      FMI%no_ptg_frq(k) = jp
      DEALLOCATE(FMI%ptg_frq_grid(k)%values,STAT=l)
      ALLOCATE(FMI%ptg_frq_grid(k)%values(jp),STAT=l)
      FMI%ptg_frq_grid(k)%values(1:jp) = &
     &           FMI%ptg_frq_grid(i)%values(1:jp)
    endif
  end do
!
 10 CLOSE(32,iostat=i)
    if(io > 0) then
      ier = io
      goto 99
    else
      io = 0
    endif
!
 99  CLOSE(32,iostat=i)
!
     if(io /= 0) then
       Ier = iabs(io)
       Call ErrMsg(Line,io)
     endif

     Return

  END SUBROUTINE PTG_FRQ_LOAD

end module PTG_FRQ_LOAD_M
! $Log$
! Revision 1.5  2001/03/20 11:03:16  zvi
! Fixing code for "real" data run, increase dim. etc.
!
! Revision 1.3  2001/03/09 00:51:28  zvi
! Fixed ier not set in ptg_frq_load
!
! Revision 1.2  2001/03/09 00:40:32  zvi
! Correcting an error in HUNT routine
!
! Revision 1.1  2001/03/06 09:28:28  zvi
! *** empty log message ***
!
! Revision 1.1 2000/06/09 00:08:13  Z.Shippony
! Initial conversion to Fortran 90
