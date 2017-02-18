! efefer 10 January 2016

MODULE gbl_poisson
  USE m_LF1d
  IMPLICIT NONE
  !
  INTEGER :: N
  REAL(8) :: L
  TYPE(LF1d_t) :: LF
  !
END MODULE


SUBROUTINE plot_rho_phi( NPLOT, A, B )
  USE m_constants, ONLY : PI
  IMPLICIT NONE
  !
  INTEGER :: NPLOT
  REAL(8) :: A, B
  REAL(8) :: x, rho, phi
  !
  INTEGER :: ip
  !
  REAL(8) :: eval_rho, eval_phi


  DO ip = 1, NPLOT
    x = A + (ip-1)*(B-A)/(NPLOT-1)
    rho = eval_rho( x )
    phi = -4.d0*PI*eval_phi( x )
    WRITE(111,*) x, rho
    WRITE(112,*) x, phi
  ENDDO


END SUBROUTINE


FUNCTION eval_rho( x ) RESULT( f )
  USE m_constants, ONLY : PI
  USE gbl_poisson, ONLY : L
  IMPLICIT NONE
  REAL(8) :: x, f
  REAL(8) :: t1, t2

  t1 = sin(2.d0*PI*x/L)**2 + cos(2.d0*PI*x/L)
  t2 = exp( -cos(2*PI*x/L) )
  !
  f = 4.d0*PI**2 * t1 * t2 / L**2

END FUNCTION


FUNCTION eval_phi( x ) RESULT(f)
  USE m_constants, ONLY : PI
  USE gbl_poisson, ONLY : L
  IMPLICIT NONE
  REAL(8) :: x, f

  f = exp( -cos( 2.d0*PI*x/L ) )
END FUNCTION


!------------------------------------------------------------------------------
PROGRAM t_poisson
!------------------------------------------------------------------------------
  USE m_constants, ONLY : PI
  USE gbl_poisson
  IMPLICIT NONE
  !
  REAL(8), ALLOCATABLE :: rho(:), phi(:)
  INTEGER :: ip
  ! For LAPACK
  REAL(8), ALLOCATABLE :: work(:), tmpD2jl(:,:)
  INTEGER, ALLOCATABLE :: ipiv(:)
  INTEGER :: lwork, info
  !
  REAL(8) :: eval_rho, eval_phi
  
  N = 35
  L = 16.d0
  
  ! Initialize the basis functions
  CALL init_LF1d_p( LF, N, 0.d0, L )

  !CALL info_LF1d( LF, .TRUE. )

  ALLOCATE( rho(N), phi(N) )
  DO ip = 1, N
    rho(ip) = eval_rho( LF%grid(ip) )
    WRITE(113,*) LF%grid(ip), rho(ip) 
  ENDDO

  CALL plot_rho_phi( 300, 0, L )

  ALLOCATE( tmpD2jl(N,N) )
  tmpD2jl(:,:) = LF%D2jl(:,:)
  ALLOCATE( ipiv(N) )
  lwork = N
  ALLOCATE( work(lwork) )

  CALL dsytrf( 'L', N, tmpD2jl,N, ipiv, work, lwork, info )
  !WRITE(*,*) 'info from dsytrf = ', info
  !WRITE(*,*) 'work(1) from dsytrf = ', work(1)

  phi(:) = -4.d0*PI*rho(:)
  CALL dsytrs( 'L',N, 1, tmpD2jl,N, ipiv, phi,N, info ) 

  DO ip = 1, N
    WRITE(114,*) LF%grid(ip), phi(ip)
  ENDDO

  DEALLOCATE( rho, phi )
  DEALLOCATE( tmpD2jl, ipiv, work )
END PROGRAM

