SUBROUTINE timestamp( messages )
  IMPLICIT NONE 
  CHARACTER(*) :: messages
  !
  INTEGER :: values(8)
  INTEGER :: year, month, day, hour, minute, second
  CHARACTER(9) :: get_month_name

  CALL date_and_time(VALUES=values)
  year   = values(1)
  month  = values(2)
  day    = values(3)
  hour   = values(5)
  minute = values(6)
  second = values(7)

  WRITE(*,*)
  WRITE(*,fmt=333) trim(messages), ' ', hour, ':', minute, ':', second,&
                   ' , ', day, ' ', trim(get_month_name(month)), ' ', year

  333 FORMAT(1x,A,4(A,I2),3A,I4)

END SUBROUTINE 
