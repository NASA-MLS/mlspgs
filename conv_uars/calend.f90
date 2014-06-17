subroutine calend(yyyy, ddd, mm, dd)

!====CALEND WHEN GIVEN A VALID YEAR, YYYY, AND DAY OF THE YEAR, DDD,
!        RETURNS THE MONTH, MM, AND DAY OF THE MONTH, DD.
!        SEE ACM ALGORITHM 398, TABLELESS DATE CONVERSION, BY
!        DICK STONE, CACM 13(10):621.

 integer, intent(in)   :: yyyy
 integer, intent(in)   :: ddd
 integer, intent(out)  :: mm
 integer, intent(out)  :: dd

 integer :: t

 t = 0
 if(modulo(yyyy, 4) == 0) t = 1

!------THE FOLLOWING STATEMENT IS NECESSARY IF YYYY IS < 1900 OR > 2100.
 if(modulo(yyyy, 400) /= 0 .and. modulo(yyyy, 100) == 0) t = 0

 dd = ddd
 if(ddd > 59+t) dd = dd + 2 - t
 mm = ((dd+91)*100)/3055
 dd = (dd+91) - (mm*3055)/100
 mm = mm - 2
!------MM WILL BE CORRECT IFF DDD IS CORRECT FOR YYYY.
 if(mm >= 1 .and. mm <= 12) return
 write(unit=*,fmt="(a,i11,a)")  &
 "$$CALEND: DAY OF THE YEAR INPUT =",ddd," IS OUT OF RANGE."
 stop
 end subroutine calend
