PROGRAM test_V_NL

  USE m_PsPot, ONLY : PsPot_Dir
  USE m_LF3d, ONLY : Npoints => LF3d_Npoints, &
                     lingrid => LF3d_lingrid, &
                     LL => LF3d_LL
  USE m_atoms, ONLY : atpos => AtomicCoords
  USE m_Ps_HGH, ONLY : hgh_eval_proj_R
  USE m_PsPot, ONLY : Ps => Ps_HGH_Params
  IMPLICIT NONE 
  INTEGER :: Narg
  INTEGER :: NN(3)
  REAL(8) :: AA(3), BB(3)
  CHARACTER(64) :: filexyz, arg_N
  INTEGER :: ip, N_in, ia, l, m, iprj, isp
  INTEGER :: iargc  ! pgf90 
  REAL(8) :: Ylm_real
  REAL(8) :: dr_vec(3)
  REAL(8) :: dr
  REAL(8), ALLOCATABLE :: betaNL(:)
  INTEGER :: Np_beta

  Narg = iargc()
  IF( Narg /= 2 ) THEN 
    WRITE(*,*) 'ERROR: exactly two arguments must be given:'
    WRITE(*,*) '       N and path to structure file'
    STOP 
  ENDIF 

  CALL getarg( 1, arg_N )
  READ(arg_N, *) N_in

  CALL getarg( 2, filexyz )

  CALL init_atoms_xyz(filexyz)

  ! Override PsPot_Dir
  PsPot_Dir = '../HGH/'
  CALL init_PsPot()

  !
  NN = (/ N_in, N_in, N_in /)
  AA = (/ 0.d0, 0.d0, 0.d0 /)
  BB = (/ 16.d0, 16.d0, 16.d0 /)
  CALL init_LF3d_p( NN, AA, BB )

  CALL info_atoms()
  CALL info_PsPot()
  CALL info_LF3d()


  CALL dealloc_atoms()
  CALL dealloc_PsPot()
  CALL dealloc_LF3d()
END PROGRAM 

