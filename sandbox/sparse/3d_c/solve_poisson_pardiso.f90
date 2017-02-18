INCLUDE 'mkl_pardiso.f90'

SUBROUTINE solve_poisson_pardiso( vecin, vecout )
  USE mkl_pardiso
  USE gbl_poisson, ONLY : Nx
  USE gbl_laplacian
  IMPLICIT NONE
  !
  REAL(8) :: vecin(Nx**3), vecout(Nx**3)  ! SPECIAL CASE for Nx=Ny=Nz
  ! Internal solver memory POINTER 
  TYPE(MKL_PARDISO_HANDLE), ALLOCATABLE  :: pt(:)
  ! All other variables
  INTEGER :: maxfct, mnum, mtype, phase, n, nrhs, error, msglvl, nnz
  INTEGER :: error1
  INTEGER, ALLOCATABLE :: iparm(:)
  INTEGER :: i, idum(1)
  REAL(8) :: ddum(1)

  n = Nx**3  ! SPECIAL CASE: Nx = Ny = Nz
  nrhs = 1
  nnz = size( Laplacian%values )
  maxfct = 1
  mnum = 1

  !
  ! Set up PARDISO control parameter
  !
  ALLOCATE( iparm ( 64 ) )

  DO i = 1, 64
    iparm(i) = 0
  ENDDO 

  iparm(1) = 1 ! no solver default
  iparm(2) = 2 ! fill-in reordering from METIS
  iparm(4) = 0 ! no iterative-direct algorithm
  iparm(5) = 0 ! no user fill-in reducing permutation
  iparm(6) = 0 ! =0 solution on the first n compoments of x
  iparm(8) = 9 ! numbers of iterative refinement steps
  iparm(10) = 13 ! perturbe the pivot elements with 1E-13
  iparm(11) = 1 ! use nonsymmetric permutation and scaling MPS
  iparm(13) = 0 ! maximum weighted matching algorithm is switched-off (default for symmetric). Try iparm(13) = 1 in case of inappropriate accuracy
  iparm(14) = 0 ! Output: number of perturbed pivots
  iparm(18) = -1 ! Output: number of nonzeros in the factor LU
  iparm(19) = -1 ! Output: Mflops for LU factorization
  iparm(20) = 0 ! Output: Numbers of CG Iterations

  error  = 0 ! initialize error flag
  msglvl = 1 ! print statistical information
  mtype  = -2 ! symmetric, indefinite
  
  !.. Initiliaze the internal solver memory pointer. This is only
  ! necessary for the FIRST call of the PARDISO solver.
  
  ALLOCATE ( pt ( 64 ) )
  do i = 1, 64
     pt( i )%DUMMY =  0 
  end do
  
  !.. Reordering and Symbolic Factorization, This step also allocates
  ! all memory that is necessary for the factorization
  
  phase = 11 ! only reordering and symbolic factorization
  
  CALL pardiso (pt, maxfct, mnum, mtype, phase, n, &
                Laplacian%values, Laplacian%rowIdx, Laplacian%column, &
                idum, nrhs, iparm, msglvl, ddum, ddum, error)
      
  WRITE(*,*) 'Reordering completed ... '
  IF (error /= 0) THEN
     WRITE(*,*) 'The following ERROR was detected: ', error
     GOTO 1000
  END IF
  WRITE(*,*) 'Number of nonzeros in factors = ',iparm(18)
  WRITE(*,*) 'Number of factorization MFLOPS = ',iparm(19)
  
  !.. Factorization.
  phase = 22 ! only factorization
  CALL pardiso (pt, maxfct, mnum, mtype, phase, n, &
                Laplacian%values, Laplacian%rowIdx, Laplacian%column, &
                idum, nrhs, iparm, msglvl, ddum, ddum, error)
  WRITE(*,*) 'Factorization completed ... '
  IF (error /= 0) THEN
     WRITE(*,*) 'The following ERROR was detected: ', error
     GOTO 1000
  ENDIF
  
  !.. Back substitution and iterative refinement
  iparm(8) = 2 ! max numbers of iterative refinement steps
  phase = 33 ! only factorization
  CALL pardiso (pt, maxfct, mnum, mtype, phase, n, &
                Laplacian%values, Laplacian%rowIdx, Laplacian%column, &
                idum, nrhs, iparm, msglvl, vecin, vecout, error)
  WRITE(*,*) 'Solve completed ... '
  IF (error /= 0) THEN
     WRITE(*,*) 'The following ERROR was detected: ', error
     GOTO 1000
  ENDIF
        
  1000 CONTINUE
  !.. Termination and release of memory
  phase = -1 ! release internal memory
  CALL pardiso (pt, maxfct, mnum, mtype, phase, n, ddum, idum, idum, &
  idum, nrhs, iparm, msglvl, ddum, ddum, error1)
  
  IF ( ALLOCATED( iparm ) )   DEALLOCATE( iparm )
  
  IF (error1 /= 0) THEN
     WRITE(*,*) 'The following ERROR on release stage was detected: ', error1
     STOP 1
  ENDIF
  
  !IF ( error /= 0 ) STOP 1
  !STOP 0
END SUBROUTINE
