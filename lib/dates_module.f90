! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.

module dates_module

  ! Converts dates between various formats. 
  ! Cal is calendar date as in "25 Jan 1998" or "25 1 1998"
  ! Separators can be any non-alphanumerc character. Leading 0s no problem
  ! If you want the three items in an order that is not day/month/year then
  ! you can supply the optional second argument perm. This must be a 3-element
  ! integer array with elements 1,2,3 in any order
  ! EUDTF is Extended UDTF -- Y2K compliant variant of UDTF
  ! EUDTF is an integer of form yyyyddd with yyyy==year and ddd==day of year
  ! (I invented this, it isn;t a standard)

  ! Since starting this, I discover that the standard text format for EOS
  ! will be CCSDS format. This comes in two sorts: 
  ! Form A:  yyyy-mm-dd   where mm= month, dd=day in month
  ! Form B:  yyyy-DDD     where DDD=day of year. No spaces. 
  ! These are clearly special cases of calendar format, and eudtf format
  ! but we provide convenience calls to change from form A to form B

  ! The SDP Toolkit has some sort of date format that fits in a double
  ! precision number. It is called the TAI93 format and is the number 
  ! of seconds since midnight, 1 Jan 1993. We provide a format called
  ! TAIDATE , which is the integer number of days since a particular 
  ! midnight. This can be used to get UARS day or TAI93 time to an accuracy
  ! of 1 Day.


  use MLSStrings
  use MLSMessageModule

  implicit none
  private

  !Here are the provided functions 
  public:: eudtf2cal,cal2eudtf,lastday,ccsdsa2b,ccsdsb2a,eudtf2daysince
  public:: daysince2eudtf,ccsds2tai,ccsds2eudtf,days_in_year
  private::isleap
  
  character(len=40),private,parameter::moduleNameIn=&
       "$RCSFile: dates_module.f90,v $"
  !Here are some useful definitions of the properties of months
  character(len=3),private,dimension(12),parameter::monthnames=&
       (/ "Jan","Feb","Mar","Apr","May","Jun","Jul","Aug",&
       "Sep","Oct","Nov","Dec" /)
  character(len=3),private,dimension(12),parameter::capmonthnames=&
       (/ "JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG",&
       "SEP","OCT","NOV","DEC" /)
  integer,private,parameter,dimension(12)::nonleap_days_in_month=&
       (/31,28,31,30,31,30,31,31,30,31,30,31 /)

contains
  function ccsds2eudtf(ccsds) result (eudtf)
    !Converts CCSDS dates to eudtf
    !----args----!
    character(len=*),intent(in)::ccsds
    !---function--result---!
    integer::eudtf
    !---local vars-----!
    character(len=30)::ccsdsb,ccsdsi
    integer::year,day
    !----Executable----!
    ccsdsi=adjustl(ccsds)
    ! Check if it looks like CCSDSA format and, if so,  bludgeon into B format
    if (ccsdsi(8:8)=="-") then
       ccsdsb=ccsdsa2b(ccsdsi)
    else
       ccsdsb=ccsdsi
    endif
    ! Get year and day out of B format
    read(unit=ccsdsb(1:4),fmt="(i4)")year
    read(unit=ccsdsb(6:8),fmt="(i3)")day
    eudtf=year*1000+day
  end function ccsds2eudtf

  function ccsds2tai(ccsds) result (tai)
    character(len=*),intent(in)::ccsds
    !---function--result---!
    integer::tai
    !----local -----!
    integer:: eudtf
    eudtf=ccsds2eudtf(ccsds)
    tai=eudtf2daysince(eudtf,1993001)
  end function ccsds2tai


  function lastday(imonth) result (day)
    integer,intent(in)::imonth
    integer::day
    if(imonth < 1 .or. imonth >12) then
       call MLSMessage( MLSMSG_Warning,moduleNameIn,&
       "in function lastday: month is out of range")
       day=31
    else
       day=nonleap_days_in_month(imonth)
    endif
  end function lastday
  ! Private function to test whether a year is a leap year or not.
  function isleap(year) result(leap)
    integer,intent(in)::year
    logical::leap
    if (modulo(year,100)==0 ) then
       if(modulo(year,400)==0 ) then
          leap=.true.
       else
          leap=.false.
       endif
    else if (modulo(year,4)==0 ) then
       leap=.true.
    else
       leap=.false.
    endif
  end function isleap


  function days_in_year(year) result(days) 
    !Private function returns no. of days in a year
    integer,intent(in)::year
    integer::days
    if ( isleap(year)) then
       days=366
    else
       days=365
    endif

  end function days_in_year

  function eudtf2daysince(eudtf,eudtf0) result(daysince)
    ! converts eudtf to days since a fixed day ( also given as EUDTF format)
    ! set eudtf0 to 11 Sept 1991 to get UARS days
    ! Set eudtf0 to 1 Jan 1993 for TAI day
    !-----args------------!
    integer,intent(in) :: eudtf,eudtf0
    !----result-----------!
    integer::daysince
    !------local vars ----!
    integer::year,year0,day,day0,daysinyear,direction
    !----------Executable----!

    year0=eudtf0/1000
    year=eudtf/1000
    day0=modulo(eudtf0,1000)
    day=modulo(eudtf,1000)
    daysince=0

    !    print*,"Start",year0,year,day0,day,daysince
    if (year < year0) then 
       
       direction=year
       year=year0
       year0=direction
       
       direction=day
       day=day0
       day0=direction

       direction=-1
    else
       direction=1
    endif

    yrsloop:do
        daysinyear=days_in_year(year0)
        if (year0 == year) then 
           daysince=daysince+day-day0
           exit yrsloop
        endif
        daysince=daysince+daysinyear-day0+1
        day0=1
        year0=year0+1
        !           print*,"In loop",year0,year,day0,day,daysince

    end do yrsloop
    !        print*,"Final",year0,year,day0,day,daysince

    daysince=daysince*direction
  end function eudtf2daysince

  function daysince2eudtf(daysince,eudtf0) result (eudtf)
    ! Converts days since a given eudtf date to an eudtf date
    !----args-----!
    integer,intent(in)::daysince,eudtf0
    !---result---!
    integer::eudtf
    !---locals---!
    integer::year,days,year0,day0,daysinyear
    !---Executable-!
    year0=eudtf0/1000
    day0=modulo(eudtf0,1000)
    year=year0
    days=daysince+day0 ! days is now days since beginning of year0
    daysinyear=days_in_year(year)
    yrsloop:do 
!        print*,"In yearsloop: days=",days," year=",year
        if (days > days_in_year(year)) then 
           days=days-days_in_year(year)
           year=year+1
        else if (days <= 0)then
           year=year-1
           days=days+days_in_year(year)
        else
           exit yrsloop 
        endif 
    enddo yrsloop
    eudtf=1000*year+days
  end function daysince2eudtf

  function eudtf2cal(eudtf,perm,num,sep) result (cal)
    ! Public function to convert eudtf (yyyyddd) to cal date
    ! Returns a character(len=11) 

    !---------args-----------
    integer,intent(in) :: eudtf
    integer,dimension(:),intent(in),optional::perm
    logical,intent(in),optional::num
    character(len=*),intent(in),optional::sep
    !-------- Result---------
    character(len=11)::cal
    !------Local vars--------
    integer:: year,dayofyear,month,dayofmonth,daysinyear,j,cumul_days,i
    !    logical:: isleap ! is this necessary???? No, it is an error
    integer,dimension(12)::days_in_month
    character(len=20),dimension(3)::tmpstring
    character(len=1)::sep_char
    character(len=5)::str1,str2
    integer,dimension(3)::order
    logical::num_month
    !--------Executable bit-------------
    !Check for optional args!
    if(present(num)) then
       num_month=num
    else
       num_month=.false.
    endif
    if(present(sep)) then
       sep_char=sep(1:1)
    else
       sep_char=" "
    endif


    ! This is necessary for consistency if perm=3,1,2 or 2,3,1
    ! just doing order=perm is the Wrong Thing here.
    if(present(perm)) then
       do i=1,3
           do j=1,3
               if(perm(j)==i) then
                  order(i)=j
               endif
           enddo
       enddo
    else
       order=(/ 1,2,3 /)
    endif
    year=eudtf/1000
    dayofyear=modulo(eudtf,1000)
    if(year <1) then ! Trap bad year
       call MLSMessage(MLSMSG_Warning,moduleNameIn,&
            "Module dates_module,function eudtf2cal: year <1" // &
            "I Can not do BC dates. Why on earth do you want one?") 
       cal="01 Jan 0001"
       return
    endif
    days_in_month=nonleap_days_in_month
    if (isleap(year)) then
       daysinyear=366
       days_in_month(2)=29 !Sodding February
    else
       daysinyear=365
    endif
    if (dayofyear < 1 .or. dayofyear > daysinyear) then
       write(str1,fmt="(i5)")dayofyear
       dayofyear=max(1,dayofyear)
       dayofyear=min(dayofyear,daysinyear)
       write(str2,fmt="(i5)")dayofyear
       call MLSMessage( MLSMSG_Warning,moduleNameIn,&
       " in function eudtf2cal: day "//str1//" is out of range."//&
       "Setting it to "//str2 )
    endif
    ! Year and day are now in range
    ! Work out which month this day is in
    cumul_days=0
    do j=1,12
        cumul_days=cumul_days+days_in_month(j)
        if (dayofyear<=cumul_days) then
           cumul_days=cumul_days-days_in_month(j)
           month=j
           exit
        endif
    enddo
    dayofmonth=dayofyear-cumul_days
    ! put the three pieces of the date into 3 elements of tmpstring
    ! in order indicated by order (i.e. by the perm arg if present)
    write(unit=tmpstring(order(1)),fmt="(i2.2)")dayofmonth
    if(num_month)then
       write(unit=tmpstring(order(2)),fmt="(i2.2)")month
    else
       tmpstring(order(2))=monthnames(month)
    endif
    write(unit=tmpstring(order(3)),fmt=*)year
    ! Stick the three bits into a string to be returned
    cal=trim(adjustl(tmpstring(1)))//sep_char//&
         trim(adjustl(tmpstring(2)))//sep_char//&
         trim(adjustl(tmpstring(3)))
    return
  end function eudtf2cal

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

  !Public function to convert cal dates (dd mmm yyyy) to eudtf (yyyyddd)
  function cal2eudtf(caldate,perm) result (eudtf)
    !------- Arguments ---------!
    character(len=*),intent(in)::caldate
    integer,dimension(:),intent(in),optional::perm
    !------- function result --------!
    integer::eudtf
    !---------Locals---------------!
    character(len=len(caldate))::dpdate ! auto-length string ???OK???
    integer::year,dayofyear,month,dayofmonth,i,j,iw
    integer::daysinyear
    integer,dimension(3)::order
    character(len=5)::errstr
    character(len=3)::monthstring
    character(len=50)::fmtstring,instring
    character(len=20),dimension(3)::tmpstring
    character(len=1)::tchar
    integer,dimension(12)::days_in_month
    !--------Executable statements--------!
    ! Get year, month and day out of string
    ! To be robust, we want the spacing to not matter and to allow
    ! things other than blank to be separators
    ! We also need to be able to handle DMY, MDY and YMD dates
    if (present(perm)) then
       order=perm
    else
       order=(/1,2,3/)
    endif
    do j=1,3
        tmpstring(j)=" "
    enddo
    dpdate=depunctuate(caldate) !get rid of punctuation
    instring=" "//adjustl(dpdate)
    iw=0
    l1:do j=1,len_trim(instring)-1
        tchar=instring(j+1:j+1)
        if(tchar /= " " .and. instring(j:j) == " ") then
           iw=iw+1
        endif
        if (iw > 3) then
           call MLSMessage( MLSMSG_Warning,moduleNameIn,&
           "in fn cal2eudtf: Warning: date"//caldate//" contains >3 words")
           exit l1
        endif
        tmpstring(order(iw))(j:j)=tchar
    enddo l1
    do j=1,3
        tmpstring(j)=adjustl(tmpstring(j))
    enddo
    fmtstring=" "
    read(unit=tmpstring(1),fmt=*)dayofmonth
    read(unit=tmpstring(2),fmt="(a)")monthstring
    read(unit=tmpstring(3),fmt=*)year
    !    print*,"dayofmonth=",dayofmonth,"Month= -->",monthstring,"<--","yr=",year
    if(year <1) then ! Trap bad year
       call MLSMessage( MLSMSG_Warning,moduleNameIn,&
            "Module dates_module,function cal2eudtf: year <1" // &
            "I Can not do BC dates. Why on earth do you want one?" ) 
       year=0000
       return
    endif

    days_in_month=nonleap_days_in_month
    if (isleap(year)) then
       daysinyear=366
       days_in_month(2)=29 !Sodding February
    else
       daysinyear=365
    endif

    !Get number of this month
    j=verify(monthstring," 0123456789")
    if(j==0) then !month was provided as number
       read(unit=monthstring,fmt=*)month
    else! Month provided as letters
       monthstring=capitalize(monthstring)
       month=0
       mcloop:do i=1,12
           if (monthstring==capmonthnames(i)) then
              month=i
              exit mcloop
           endif
       enddo mcloop
       if(month==0) then 
          call MLSMessage(MLSMSG_Warning,moduleNameIn,&
               "in function cal2eudtf: Cannot interpret month "//monthstring)
       endif
    endif
    if (month < 1 .or. month > 12) then
       write(unit=errstr,fmt="(i5)")month
        call MLSMessage(MLSMSG_Warning,moduleNameIn,&
             "in function cal2eudtf Month "//errstr//" Not valid. "//&
             "Setting it to 1 (==Jan)")
        month=1
    endif


    dayofyear=sum(days_in_month(1:month-1)) + dayofmonth
    eudtf=year*1000 + dayofyear

  end function cal2eudtf

  function ccsdsa2b(a) result(b)
    !-----Arg--------!
    character(len=*),intent(in)::a
    !-----Function result-----!
    character(len=8)::b
    !------locals--------!
    integer::eudtf
    eudtf=cal2eudtf(a,perm=(/3,2,1/))
    write(unit=b(1:4),fmt="(i4.4)")eudtf/1000
    write(unit=b(6:8),fmt="(i3.3)")modulo(eudtf,1000)
    b(5:5)="-"
  end function ccsdsa2b

  function ccsdsb2a(b) result(a)
    !-----Arg--------!
    character(len=*),intent(in)::b
    !-----Function result-----!
    character(len=10)::a
    !------locals--------!
    integer::eudtf,year
    character(len=20)::btr

    btr=adjustl(b)
    read(unit=btr(1:4),fmt="(i4)")year
    read(unit=btr(6:8),fmt="(i3)")eudtf
    eudtf=eudtf+1000*year
    a=eudtf2cal(eudtf,perm=(/3,2,1/),sep="-",num=.true.)

  end function ccsdsb2a

end module dates_module
