! Copyright (c) 1999, California Institute of Technology.  ALL RIGHTS RESERVED.
! U.S. Government Sponsorship under NASA Contract NAS7-1407 is acknowledged.


! --------------------------------------------------
 subroutine read_dao (fname, vname, data_array)
! --------------------------------------------------
!------------------- RCS Ident Info -----------------------
CHARACTER(LEN=130) :: Id = &                                                    
"$Id$"
!----------------------------------------------------------

! Brief description of program
! This subroutine reads the DAO correlative file and returns
! the data_array to the caller

use GEOSAVG_EOS_1
implicit none

! Arguments

      character*(*), intent(IN) :: fname, vname

! - - - local declarations - - -

  integer*4 sd_id, sds_id, status
  integer*4 sds_index
  integer*4 start(0:3), edges(0:3), stride(0:3)
  real*4 data_array(0:XDIM, 0:YDIM, 0:ZDIM)
  integer i,j,k
  integer file_id, gd_id
  integer, external :: sdstart, gdopen, gdattach, gdrdfld, gddetach, &
                       & gdclose

! - - - begin - - -
  sd_id = sdstart(fname, DFACC_RDONLY)
  file_id = gdopen(fname, DFACC_RDONLY)

  if (file_id < 0) then
    write (*, *) "Could not open ",fname
    stop '(-1)'
  end if

  gd_id = gdattach(file_id,"EOSGRID")

  if (gd_id < 0) then
    write (*,*) "Could not open ",fname
    stop '(-1)'
  end if

    start(0) = 0
    start(1) = 0
    start(2) = 0
    start(3) = 0

    stride(0) = 1
    stride(1) = 1
    stride(2) = 1
    stride(3) = 1

    edges(0) = 1
    edges(1) = ZDIM
    edges(2) = YDIM
    edges(3) = XDIM


!   In this subroutine, we read the entire field.  By manipulating the start 
!   and edges arrays, it is possible to read a subset of the entire array.  
!   For example, to read a 3D section defined by x=100,224 y=50,149 
!   z=15,16 you would set the start and edges arrays to the following:

!   start(0) = 0    time start location
!   start(1) = 15   z-dim start location
!   start(2) = 50   y-dim start location
!   start(3) = 100  x-dim start location

!   edges(0) = 1    time length
!   edges(1) = 2    z-dim length
!   edges(2) = 100  y-dim length
!   edges(3) = 125  x-dim length



  status = gdrdfld(gd_id, vname, start, stride, edges, data_array)

  write (*,*) "Read status=",status
  status = gddetach(gd_id)
  status = gdclose(file_id)

 
return
end
