program test_dates
use dates_module ! Module to handle date conversions
integer::eudtf,ds,eudtf2!,eudtf0
integer,dimension(3)::foo
character(len=30)::caldate
character(len=8):: ccsdsb
foo=(/1,2,3/)

print*,"Enter a date in ccsds A format (yyyy-mm-dd)"
read(unit=*,fmt="(a)")caldate
eudtf=cal2eudtf(caldate,perm=(/3,2,1/))
print*,"EUDTF is ",eudtf
ds=eudtf2daysince(eudtf,1991254)
print*,"UARS day:",ds
ccsdsb=ccsdsa2b(caldate)
print*,"B format is ",ccsdsb
caldate=ccsdsb2a(ccsdsb)

print*,"converts back as",caldate

print*,"Enter a eudtf date"
read*,eudtf
caldate=eudtf2cal(eudtf,perm=foo)
print*,"Cal date is -->",caldate,"<--"
eudtf=cal2eudtf(caldate)
ds=eudtf2daysince(eudtf,1991254)
print*,"UARS day:",ds
print*,"Converts back as",eudtf
print*,"Enter a calendar date"
read(unit=*,fmt="(a)")caldate
eudtf=cal2eudtf(caldate)
print*,"EUDTF=",eudtf
ds=eudtf2daysince(eudtf,1991254)
print*,"UARS day:",ds

caldate=eudtf2cal(eudtf,sep="@")
print*,"Converts back as---",caldate,"---"

print*,"Enter a series of dates  in cal format (e.g. 1 Jun 1993)"

do 
    if (eudtf == 2000001) then
       exit
    endif
    read(unit=*,fmt="(a)")caldate
    eudtf=cal2eudtf(caldate)
    ds=eudtf2daysince(eudtf,1991254)
    eudtf2=daysince2eudtf(ds,1991254)
    print*,caldate(1:14)," UARS day:",ds,"EUDTF is ",eudtf,eudtf2,&
         " err",eudtf2-eudtf
    
end do


print*,"Enter a series of dates in CCSDS  format (1990-001 ends)"

do 
    if (eudtf == 1900001) then
       exit
    endif
    read*,caldate
    eudtf=ccsds2eudtf(caldate)
    ds=ccsds2tai(caldate)
    eudtf2=daysince2eudtf(ds,1993001)

    print*,caldate,ds,eudtf,eudtf2,eudtf-eudtf2
    
enddo


!print*,cal2eudtf(caldate)! why does this cause a recursive io reference?


end program test_dates
