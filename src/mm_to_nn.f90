FUNCTION mm_to_nn( mm, S ) RESULT(idx)
  IMPLICIT NONE
  INTEGER :: idx
  INTEGER :: mm
  INTEGER :: S
  IF(mm > S/2) THEN 
    idx = mm - S
  ELSE
    idx = mm
  ENDIF
END FUNCTION 

