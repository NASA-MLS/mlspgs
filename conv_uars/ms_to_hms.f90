subroutine ms_to_hms (ms, hrs, mins, secs)

! convert millisecs to hours, minutes and seconds

integer, intent (in) :: ms
integer, intent (out) :: hrs, mins, secs

hrs = ms / 3600000
mins = mod (ms, 3600000) / 60000
secs = mod (mod (ms, 3600000), 60000) / 1000

end subroutine ms_to_hms
