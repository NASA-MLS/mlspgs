!INTEGER(kind=2) FUNCTION SwapShort (i2)
integer*2 FUNCTION SwapShort (i2)

  !INTEGER(kind=2) :: i2
  integer*2 :: i2
  CHARACTER(len=1) :: cbuf(2)

  cbuf = TRANSFER (i2, cbuf)
  SwapShort = TRANSFER (cbuf(2:1:-1), i2)

END FUNCTION SwapShort
