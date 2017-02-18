FUNCTION mathieu( Vm, L, x )
  USE m_constants, ONLY : PI
  IMPLICIT NONE
  !
  REAL(8) :: mathieu
  REAL(8) :: Vm, L, x

  mathieu = Vm + Vm*cos( 2*PI*x/L )
END FUNCTION

!
FUNCTION dmathieu( Vm, L, x )
  USE m_constants, ONLY : PI
  IMPLICIT NONE
  !
  REAL(8) :: dmathieu
  REAL(8) :: Vm, L, x

  dmathieu = -2.d0*PI*Vm*sin(2*PI*x/L)/L
END FUNCTION

FUNCTION d2mathieu( Vm, L,x )
  USE m_constants, ONLY : PI
  IMPLICIT NONE
  !
  REAL(8) :: d2mathieu
  REAL(8) :: Vm, L, x

 d2mathieu = -4.d0*pi**2 * Vm*cos(2.d0*pi*x/L)/L**2
END FUNCTION


!------------------------------------------------------------------------------
PROGRAM t_fun
!------------------------------------------------------------------------------
  USE m_constants, ONLY : PI
  USE m_LF1d
  IMPLICIT NONE
  !
  INTEGER :: N = 201
  REAL(8) :: L = 5.d0
  !
  TYPE(LF1d_t) :: LF
  REAL(8), ALLOCATABLE :: coefs(:), dcoefs(:), d2coefs(:)
  !
  INTEGER :: ii, jj
  ! Functions
  REAL(8) :: mathieu, dmathieu, d2mathieu
  !
  REAL(8) :: Vm
  INTEGER :: NPTS_PLOT=201
  REAL(8) :: xx, yy

  Vm = 1.5

  ! Initialize the basis functions
  CALL init_LF1d_p( LF, N, 0.d0, L )

  ! Calculate coefficients of expansion for the function
  ALLOCATE( coefs(N) )
  ALLOCATE( dcoefs(N) )
  ALLOCATE( d2coefs(N) )

  DO ii=1,N
    coefs(ii) = sqrt(L/N)* mathieu( Vm, L, LF%grid(ii) ) 
    WRITE(11,*) LF%grid(ii), sqrt(N/L)*coefs(ii)
  ENDDO

  dcoefs = matmul( LF%D1jl, coefs )
  d2coefs = matmul( LF%D2jl, coefs )

  ! Write to file for plotting
  DO ii=1,NPTS_PLOT
    !
    !xx = (ii-1)*L/NPTS_PLOT
    xx = LF%grid(ii)  ! only works for NPTS_PLOT = N
    yy = 0.d0
    DO jj=1,N
      yy = yy + eval_LF1d_p( LF, jj, xx )*coefs(jj)
    ENDDO
    WRITE(12,*) xx
    WRITE(13,*) yy
    WRITE(14,*) mathieu( Vm, L, xx )
    !
    ! First derivative
    !
    yy = 0.d0
    DO jj=1,N
      !yy = yy + eval_LF1d_p_d1( LF, jj, xx )*coefs(jj)
      yy = yy + eval_LF1d_p( LF, jj, xx )*dcoefs(jj)
    ENDDO
    WRITE(15,*) yy
    WRITE(16,*) dmathieu( Vm, L, xx )
    !
    ! Second derivative
    !
    yy = 0.d0
    DO jj=1,N
      yy = yy + eval_LF1d_p( LF, jj, xx )*d2coefs(jj)
    ENDDO
    WRITE(17,*) yy
    WRITE(18,*) d2mathieu( Vm, L, xx )
  ENDDO

  ! Free memory
  CALL dealloc_LF1d( LF )
  DEALLOCATE( coefs, dcoefs, d2coefs )


END PROGRAM

