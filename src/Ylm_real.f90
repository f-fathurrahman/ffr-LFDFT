FUNCTION Ylm_real( l, m, R ) RESULT(ylm)

  USE m_constants, ONLY : PI
  IMPLICIT NONE 
  REAL(8) :: ylm
  INTEGER :: l, m
  REAL(8) :: R(3)
  REAL(8) :: cost, sint, phi
  REAL(8) :: Rmod
  REAL(8), PARAMETER :: SMALL = 1.0d-9

  Rmod = sqrt( R(1)**2 + R(2)**2 + R(3)**2 )
  IF( Rmod < SMALL ) THEN
    cost = 0.d0
  ELSE
    cost = R(3)/Rmod
  ENDIF 
  !
  ! beware the arc tan, it is defined modulo pi
  !
  IF( R(1) > SMALL ) THEN
    phi = atan( R(2)/R(1) )
  ELSEIF( R(1) < -SMALL ) THEN
    phi = atan( R(2)/R(1) ) + PI
  ELSE
    phi = sign( PI/2.d0, R(2) )
  ENDIF
  sint = sqrt( max(0d0, 1.d0 - cost**2) )


  ylm = 0.d0
  IF( l == 0 ) THEN
    ylm = 0.5d0*sqrt(1.d0/PI)
    RETURN

  ELSEIF( l == 1 ) THEN
    ! py
    IF( m == -1 ) THEN
      ylm = 0.5d0*sqrt(3.d0/PI)*sint*sin(phi)
      RETURN
    ! pz
    ELSEIF( m == 0 ) THEN
      ylm = 0.5d0*sqrt(3.d0/pi)*cost
      RETURN 
    ! px
    ELSEIF( m == 1 ) THEN
      ylm = 0.5d0*sqrt(3.d0/pi)*sint*cos(phi)
      RETURN 
    ENDIF 

  ELSEIF( l == 2 ) THEN 
    ! dxy
    IF( m == -2 ) THEN 
      ylm = sqrt(15.d0/16.d0/PI) * sint**2 * sin(2.d0*phi)
      RETURN 
    ! dyz
    ELSEIF( m == -1 ) THEN 
      ylm = sqrt(15.d0/(4.d0*PI))*cost*sint*sin(phi)
      RETURN 
    ! dz2
    ELSEIF( m == 0 ) THEN 
      ylm = 0.25d0*sqrt(5.d0/PI)*( 3.d0*cost**2 - 1.d0 )
      RETURN 
    ! dxz
    ELSEIF( m == 1 ) THEN 
      ylm = sqrt(15.d0/4.d0/PI)*cost*sint*cos(phi)
      RETURN 
    ! dx2-y2
    ELSEIF( m == 2 ) THEN 
      ylm = 0.5d0*sqrt(15.d0/4.d0/PI) * sint**2 * cos(2.d0*phi)
      RETURN 
    ENDIF

  ELSEIF( l == 3 ) THEN

    IF( m == -3 ) THEN 
      ylm = 0.25d0*sqrt(35.d0/2.d0/PI) * sint**3 * sin(3.d0*phi)
      RETURN 

    ELSEIF( m == -2 ) THEN 
      ylm = 0.25d0*sqrt(105.d0/PI)* sint**2 *cost * sin(2.d0*phi)
      RETURN 

    ELSEIF( m == -1 ) THEN 
      ylm = 0.25d0*sqrt(21.d0/2.d0/PI)*sint*( 5d0*cost**2 - 1.d0 )*sin(phi)
      RETURN 

    ELSEIF( m == 0 ) THEN 
      ylm = 0.25d0*sqrt(7.d0/PI)*( 5.d0*cost**3 - 3.d0*cost )
      RETURN 

    ELSEIF( m == 1 ) THEN 
      ylm = 0.25*sqrt(21.d0/2.d0/PI)*sint*( 5.d0*cost**2 - 1.d0 )*cos(phi)
      RETURN 

    ELSEIF( m == 2 ) THEN 
      ylm = 0.25d0*sqrt(105.d0/PI) * sint**2 * cost * cos(2.d0*phi)
      RETURN 

    ELSEIF( m == 3 ) THEN 
      ylm = 0.25d0*sqrt(35.d0/2.d0/PI) * sint**3 * cos(3.d0*phi)
      RETURN 
    ENDIF

  ENDIF 

END FUNCTION 


