program test_rl

use MLSStrings

implicit none
logical :: eof
character(len=80)::filename,inline
integer::linecount
print*,"Enter a file name"
read(unit=*,fmt="(a)")filename
open(unit=1,file=filename,status="old",action="read")

lineloop:do linecount=1,1000
    call ReadCompleteLineWithoutComments(1,inline,eof=eof)
    print*,"----------------L-i-n-e-",linecount,"-------E-O-F-=",eof
    print*,inline
    print*,"----------E-n-d-L-i-n-e-",linecount,"--------------"
    if (eof) then
       exit lineloop
    end if

end do lineloop

print*,"Finished"

end program test_rl
