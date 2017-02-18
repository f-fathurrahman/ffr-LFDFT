! ffr: 3 May 2015
! experimental implementation of periodic Lagrange basis functions
! in one dimension

MODULE m_LF1d

USE m_constants, ONLY: PI
IMPLICIT NONE


!------------------------------------------------------------------------------
TYPE LF1d_t
  INTEGER :: N
  REAL(8) :: L
  REAL(8) :: A, B
  REAL(8) :: h
  REAL(8), ALLOCATABLE :: grid(:)
  REAL(8), ALLOCATABLE :: D1jl(:,:), D2jl(:,:)
  ! Transformation matrix, for the moment only applicable to periodic LF
  COMPLEX(8), ALLOCATABLE :: Tal(:,:)
END TYPE
!------------------------------------------------------------------------------



  CONTAINS


SUBROUTINE init_transform_matrix_p( LF )
  IMPLICIT NONE
  !
  TYPE(LF1d_t) :: LF
  !
  INTEGER :: alpha, l
  INTEGER :: kl
  ! N must be an odd number

  kl = -(LF%N-1)/2
  !kl = 0
  DO l = 1, LF%N
    DO alpha = 1, LF%N
      !WRITE(*,*) 'kl = ', kl
      LF%Tal(alpha,l) = exp( 2.d0*PI*(0.d0,1.d0)*kl*LF%grid(alpha)/LF%L ) / sqrt(dble(LF%N))
    ENDDO
    kl = kl + 1
  ENDDO
END SUBROUTINE




! Initialize Langrange function with cluster boundary condition
! The LFs are defined on [A,B] with B > A
!------------------------------------------------------------------------------
SUBROUTINE init_LF1d_c( LF, N, A, B )
!------------------------------------------------------------------------------
  IMPLICIT NONE
  !
  TYPE(LF1d_t) :: LF
  INTEGER :: N
  REAL(8) :: A, B
  !
  INTEGER :: ii
  !
  LF%N = N
  LF%A = A
  LF%B = B
  !
  ! Check if B > A
  !
  IF( B < A ) THEN
    WRITE(*,*) 'Error in init_LF1d_c: B should be larger that A', A, B
    STOP
  ENDIF
  LF%L = B - A
  LF%h = (B - A)/(N + 1) ! Note that this is different from the periodic LF
  ALLOCATE( LF%grid(N) )
  ! Initializatio
  DO ii=1,N
    LF%grid(ii) = A + ii*(B-A)/(N+1)
  ENDDO
  !
  ! FIXME: Currently only D2jl matrix is initialized
  ALLOCATE( LF%D2jl(N,N) )
  CALL init_deriv_matrix_c( LF%D2jl, LF )
  !
  WRITE(*,'(1x,A,F18.10)') 'Allocated: cluster 1d LBF grid, h = ', LF%h

END SUBROUTINE


!-----------------------------------------
SUBROUTINE init_deriv_matrix_c( D2jl, LF )
!-----------------------------------------
  IMPLICIT NONE
  !
  TYPE(LF1d_t) :: LF
  REAL(8) :: D2jl(:,:)
  !
  REAL(8) :: L, t1, t2, pre
  INTEGER :: N, nnm, nnp
  INTEGER :: ll, jj

  L = LF%L
  N = LF%N
  !
  ! Diagonal
  !
  pre = -PI**2/(2*L**2)
  DO ll = 1, N
    t1 = ( 2*(N+1)**2 + 1 )/3.d0
    t2 = sin( PI*ll/(N+1) )**2
    D2jl(ll,ll) = pre*(t1 - 1.0/t2)
  ENDDO
  !
  ! Off-diagonal
  !
  DO ll = 1, N
    DO jj = ll+1, N
      nnm = ll - jj
      nnp = ll + jj
      pre = -PI**2/(2*L**2)*(-1)**nnm
      t1 = sin( PI*nnm/2.d0/(N+1) )**2
      t2 = sin( PI*nnp/2.d0/(N+1) )**2
      !
      D2jl(ll,jj) = pre*(1.d0/t1 - 1.d0/t2)
      D2jl(jj,ll) = pre*(1.d0/t1 - 1.d0/t2)
    ENDDO
  ENDDO
END SUBROUTINE


! Evaluate the value of L(x) for arbitrary x according to the analytic form
! of the basis function.
! TODO: Rename to eval_lbd1dp
!------------------------------------------------------------------------------
FUNCTION eval_LF1d_p(LF, ibf, x) RESULT(ff)
!------------------------------------------------------------------------------
  IMPLICIT NONE
  ! arguments
  TYPE(LF1d_t) :: LF
  INTEGER :: ibf
  REAL(8) :: x, ff
  !
  REAL(8) :: pre1, L
  INTEGER :: ii, N
  !
  N = LF%N
  L = LF%L
  pre1 = 1.d0/sqrt(N*L)
  ff = 0.d0
  DO ii=1,LF%N
    ff = ff + cos( PI*( 2*ii - N - 1 )*( x - LF%grid(ibf) )/ L )
  ENDDO
  ff = ff*pre1
END FUNCTION


FUNCTION eval_LF1d_c( LF, ibf, x ) RESULT(ff)
  IMPLICIT NONE
  !
  REAL(8) :: ff
  TYPE(LF1d_t) :: LF
  INTEGER :: ibf
  REAL(8) :: x
  !
  REAL(8) :: ki, pre, L, A
  INTEGER :: N, ii

  N = LF%N
  L = LF%L
  A = LF%A
  pre = 2.d0/sqrt( (N+1)*L )
  !pre = 2.d0/(N+1)
  ff = 0.d0
  DO ii = 1, N
    ki = PI*ii/L
    ff = ff + sin( ki*(LF%grid(ibf) - A) ) * sin( ki*(x - A) )
  ENDDO
  ff = ff*pre
END FUNCTION


!------------------------------------
SUBROUTINE init_LF1d_sinc( LF, N, h )
!------------------------------------
  IMPLICIT NONE
  !
  TYPE(LF1d_t) :: LF
  INTEGER :: N
  REAL(8) :: h
  INTEGER :: ip

  LF%N = N
  LF%h = h
  ALLOCATE( LF%grid(N) )
  LF%A = -(N-1)/2.d0*h
  LF%B =  (N-1)/2.d0*h
  DO ip = 1, N   ! BEWARE: integer division
    LF%grid(ip) = LF%A + (ip-1)*h
  ENDDO
  ! L is not important in this case?
  LF%L = LF%B - LF%A

  ! FIXME: Currently only D2jl matrix is initialized
  ALLOCATE( LF%D2jl(N,N) )
  CALL init_deriv_matrix_sinc( LF )
  !
  WRITE(*,'(1x,A,2F18.10)') 'Allocated: cluster 1d LF sinc, h, L = ', LF%h, LF%L
END SUBROUTINE

SUBROUTINE init_deriv_matrix_sinc( LF )
  IMPLICIT NONE
  !
  TYPE(LF1d_t) :: LF
  INTEGER :: i, j
  
  ! Diagonal part
  DO i = 1, LF%N
    LF%D2jl(i,i) = -pi**2/3.d0/LF%h**2
  ENDDO
  ! Off-diagonal
  DO j=1,LF%N
    DO i=j+1,LF%N
      LF%D2jl(i,j) = -2.d0*(-1)**(i-j)/(LF%grid(i)-LF%grid(j))**2
      LF%D2jl(j,i) = LF%D2jl(i,j)
    ENDDO
  ENDDO
END SUBROUTINE


FUNCTION eval_LF1d_sinc( LF, ibf, x ) RESULT(ff)
  IMPLICIT NONE
  !
  TYPE(LF1d_t) :: LF
  REAL(8) :: ff, x
  INTEGER :: ibf
  !
  REAL(8) :: h
  
  h = LF%h
  ff = sin( PI*(x - LF%grid(ibf))/h ) / (PI*(x-LF%grid(ibf))) * h /sqrt(h)
END FUNCTION


SUBROUTINE info_LF1d( LF, l_print_grid )
  IMPLICIT NONE
  !
  TYPE(LF1d_t) :: LF
  LOGICAL :: l_print_grid
  !
  INTEGER :: ip
  !
  WRITE(*,*)
  WRITE(*,*) 'LF1d info:'
  WRITE(*,'(1x,A,I10)') 'N = ', LF%N
  WRITE(*,'(1x,A,F16.5)') 'L = ', LF%L
  WRITE(*,'(1x,A,F16.5)') 'h = ', LF%h
  IF( l_print_grid ) THEN
    WRITE(*,*) 'Grid points:'
    DO ip=1, LF%N
      WRITE(*,'(1x,F18.10)') LF%grid(ip)
    ENDDO
  ENDIF
END SUBROUTINE


!------------------------------------------------------------------------------
SUBROUTINE dealloc_LF1d( LF )
!------------------------------------------------------------------------------
  IMPLICIT NONE
  !
  TYPE(LF1d_t) :: LF

  DEALLOCATE( LF%grid, LF%D2jl )
  IF( allocated(LF%D1jl) ) DEALLOCATE( LF%D1jl )
END SUBROUTINE


END MODULE


