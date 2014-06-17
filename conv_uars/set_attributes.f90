subroutine set_attributes (sdid, yrdoy, start_time, end_time)

USE PCFHdr, ONLY: h5_writeglobalattr, GlobalAttributes, FillTAI93Attribute

implicit none

integer, intent (in) :: sdid, yrdoy, start_time(2), end_time(2)

integer :: day, doy, month, year, yr_range(2), doy_range(2)
integer :: hrs, mins, secs

character (len=17) :: start_utc, end_utc
CHARACTER (LEN=*), PARAMETER :: &
     DateFmt = "(I4, '-', I3.3, 'T', I2.2, ':', I2.2, ':', I2.2)"
CHARACTER (LEN=*), PARAMETER :: inst_name = 'MLS UARS'

! year and doy of file:

year = yrdoy / 1000 + 1900
doy = mod (yrdoy, 1000)

! convert DOY to Month, Day:

call calend (year, doy, month, day)
print *, 'year, doy, month, day: ', year, doy, month, day

! convert start millisecs to HMS:

call ms_to_hms (start_time(2), hrs, mins, secs)

yr_range(1) = start_time(1) / 1000 + 1900
doy_range(1) = mod (start_time(1), 1000)
write (start_utc, fmt=DateFmt) yr_range(1), doy_range(1), hrs, mins, secs

! convert stop millisecs to HMS:

call ms_to_hms (end_time(2), hrs, mins, secs)

yr_range(2) = end_time(1) / 1000 + 1900
doy_range(2) = mod (end_time(1), 1000)
write (end_utc, fmt=DateFmt) yr_range(2), doy_range(2), hrs, mins, secs

! Store global attributes

GlobalAttributes%StartUTC = Start_UTC
GlobalAttributes%EndUTC = End_UTC
GlobalAttributes%ProcessLevel = '1'
GlobalAttributes%InstrumentName = 'MLS UARS'
GlobalAttributes%PGEVersion = 'Test'
GlobalAttributes%GranuleYear =  year
GlobalAttributes%GranuleMonth = month
GlobalAttributes%GranuleDay = day

CALL FillTAI93Attribute

! Write global attributes:

call h5_writeglobalattr (sdid, skip_if_already_there=.false.)

end subroutine set_attributes
