; this function makes a scalar out of a vector
; it will take the first element if input
; has more than one element
  FUNCTION scalarize,vector
  scalar = vector
  RETURN,scalar(0)
  END
