module TOGGLES

  implicit NONE
  private

  integer, parameter, public :: MAX_TOGGLE = 255
  integer, public :: TOGGLE(0:MAX_TOGGLE) = 0 ! Indexed by IACHAR(...)

end module TOGGLES
