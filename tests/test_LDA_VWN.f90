! A simple test program to check whether my LDA_VWN is
! giving the same answer as the Julia version (which in turn is
! converted from the original Octave code by T.A. Arias
PROGRAM test_LDA_VWN

  IMPLICIT NONE
  REAL(8), ALLOCATABLE :: Rhoe(:)
  REAL(8), ALLOCATABLE :: epsxc(:), depsxc(:)
  INTEGER :: Npoints

  Npoints = 5
  ALLOCATE( Rhoe(Npoints) )
  ALLOCATE( epsxc(Npoints) )
  ALLOCATE( depsxc(Npoints) )

  Rhoe(:) = 1.5d0

  CALL excVWN( Npoints, Rhoe, epsxc )
  CALL excpVWN( Npoints, Rhoe, depsxc )

  WRITE(*,*) epsxc
  WRITE(*,*) depsxc

  DEALLOCATE( Rhoe )
  DEALLOCATE( epsxc )
  DEALLOCATE( depsxc )

END PROGRAM 

