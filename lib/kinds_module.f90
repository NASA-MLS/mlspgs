
module kinds_module
! This module contains kind parameters for the usual real and integer types 
! It is used by a number of routines supplied by HCP. 
implicit none

integer,public,parameter:: i1=selected_int_kind(2)
integer,public,parameter:: i2=selected_int_kind(4)
integer,public,parameter:: i4=selected_int_kind(7)
integer,public,parameter:: r4=selected_real_kind(5)
integer,public,parameter:: r8=selected_real_kind(13)

end module kinds_module
