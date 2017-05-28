PROGRAM test_eggbox

  USE m_PsPot, ONLY : PsPot_Dir
  USE m_LF3d, ONLY : xyz2lin => LF3d_xyz2lin, &
                     lingrid => LF3d_lingrid, &
                     Npoints => LF3d_Npoints, &
                     dVol => LF3d_dVol, &
                     LL => LF3d_LL
  USE m_hamiltonian, ONLY : V_ps_loc
  USE m_atoms, ONLY : atpos => AtomicCoords
  USE m_constants, ONLY : ANG2BOHR, PI
  IMPLICIT NONE 
  INTEGER :: Narg
  INTEGER :: NN(3)
  REAL(8) :: AA(3), BB(3)
  CHARACTER(64) :: arg_N
  INTEGER :: ip, N_in, ix, iy, iz
  REAL(8) :: integ, scal
  REAL(8), ALLOCATABLE :: psi(:)

  Narg = iargc()
  IF( Narg /= 1 ) THEN 
    WRITE(*,*) 'ERROR: exactly one arguments must be given: N'
    STOP 
  ENDIF 

  CALL getarg( 1, arg_N )
  READ(arg_N, *) N_in

  CALL init_atoms_xyz('ATOM.xyz')
  ! so that coord given in xyz file is in bohr
  atpos(:,:) = atpos(:,:)/ANG2BOHR  

  ! Override PsPot_Dir
  PsPot_Dir = '../../HGH/'
  CALL init_PsPot()

  !
  NN = (/ N_in, N_in, N_in /)
  AA = (/ 0.d0, 0.d0, 0.d0 /)
  BB = (/ 16.d0, 16.d0, 16.d0 /)
  CALL init_LF3d_p( NN, AA, BB )

  CALL init_strfact_shifted()

  CALL alloc_hamiltonian()

  CALL init_V_ps_loc_G()

  ALLOCATE( psi(Npoints) )

  ! constant
  DO ip = 1,Npoints
    psi(ip) = sin( lingrid(1,ip)*2.d0*PI/LL(1) ) * &
              cos( lingrid(2,ip)*2.d0*PI/LL(2) ) * &
              sin( lingrid(3,ip)*2.d0*PI/LL(3) )
  ENDDO 

  integ = sum( V_ps_loc(:)*psi(:)**2 ) * dVol
  WRITE(*,'(1x,A,2F18.10)') 'pos, integ = ', atpos(2,1), integ

  ix = 1
  iz = 1
  WRITE(*,*) 'ix iz = ', ix, iz
  scal = maxval( abs(V_ps_loc) )/maxval( abs(psi) )
  WRITE(*,*) 'scal = ', scal
  DO iy = 1,NN(2)
    ip = xyz2lin(ix,iy,iz)
    WRITE(N_in,'(3F22.12)') lingrid(2,ip), V_ps_loc(ip), psi(ip)*scal*10.d0
  ENDDO 


  DEALLOCATE( psi )

  CALL dealloc_hamiltonian()
  CALL dealloc_LF3d()
  CALL dealloc_PsPot()
  CALL dealloc_atoms()

END PROGRAM

