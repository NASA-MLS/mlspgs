!=================================
PROGRAM wait_test ! tests subroutine wait
!=================================

   USE MLSCommon , ONLY: r8
   USE time_m , ONLY: wait, time_now, WAIT_LOOP_LIMITS, &
   & RETRY, INIT_RETRY, TRY_AGAIN, RETRY_SUCCESS, TOO_MANY_RETRIES   
   
   IMPLICIT NONE

!------------------- RCS Ident Info -----------------------
   CHARACTER(LEN=130) :: Id = &                                                    
   "$Id$"
   CHARACTER (LEN=*), PARAMETER :: ModuleName= "$RCSfile$"
!----------------------------------------------------------

! Brief description of program
! This program tests the wait subroutine.

! To use this, copy it into
! mlspgs/tests/lib
! then enter "make depends" followed by "make"


! Then run it, entering how long you want to wait

! Variables

   real(r8) :: theWait, t0, t1
   integer :: ErrTyp, max_retries, shall_i
   real :: delay, maxTime

 	print *, 'Loop limits in wait routine (or enter 0 to use default)'  
	read(*, '(i10)') ErrTyp   
	IF(ErrTyp > 0) THEN  
      WAIT_LOOP_LIMITS = ErrTyp      
	ENDIF                        
	
	DO

! Prompt for input

 	  print *, 'How long to wait'
		read(*, '(f10.0)') theWait
			IF(theWait <= 0.d0) THEN
   		print *, 'GAME OVER'
			exit
		ENDIF
	
!	Process theString
    call time_now(t0)
    print *, 'Waiting ...'
    call wait(theWait, ErrTyp)
    print *, 'Wait over!'

    call time_now(t1)                         
    if(ErrTyp==0) then                        
      print *, 'You wanted to wait', theWait  
      print *, 'You actually waited', t1-t0   
    else                                      
      print *, '(Error in wait routine)'      
      print *, 'You wanted to wait', theWait  
      print *, 'You actually waited', t1-t0   
    endif                                     

	ENDDO
	DO
! Prompt for input

 	  print *, 'delay, max retries, max time for retries'
		read(*,*) delay, max_retries, maxTime
			IF(delay <= 0.d0) THEN
   		print *, 'GAME OVER'
			exit
		ENDIF
    call init_retry(SUCCESSFUL_RESULT=0)
    do
       call home(ErrTyp)
       shall_i = retry(ErrTyp, delay=delay, max_retries=max_retries)
       if ( shall_i /= try_again) exit
    enddo
    if ( shall_i /= RETRY_SUCCESS ) then
       call exception_handler(shall_i)
    endif
	ENDDO
contains
  subroutine home ( result )
    integer, intent(out) :: result
    character(len=*), parameter :: file_name='temp.file'
    open(1, file=file_name, iostat=result, status='old', form='formatted')
    if ( result /= 0 ) then
      print *, 'Sorry, did not find ', file_name
    else    
      print *, 'Whoo-hoo-hoo, found ', file_name
    endif
    close(1)
  end subroutine home

  subroutine exception_handler ( result )
    integer, intent(in) :: result
    select case (result)
    case ( TRY_AGAIN )
      print *, 'Died despite being prompted to try again'
    case ( RETRY_SUCCESS )
      print *, 'Died despite succeeding'
    case ( TOO_MANY_RETRIES )
      print *, 'Died after too many retries'
    case default
      print *, 'Died with unknown result: ', result
    end select
  end subroutine exception_handler

!==================
END PROGRAM wait_test
!==================

!# $Log$
!#
