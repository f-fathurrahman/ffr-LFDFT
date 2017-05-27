MODULE m_input_vars

  IMPLICIT NONE 
  
  ! CONTROL
  CHARACTER(128) :: pseudo_dir
  REAL(8) :: etot_conv_thr
  INTEGER, PARAMETER :: IU = 100
  NAMELIST /CONTROL/ pseudo_dir, etot_conv_thr

  ! SYSTEM
  REAL(8) :: A, B, C
  INTEGER :: nr1, nr2, nr3
  INTEGER :: nat, ntyp
  INTEGER :: ibrav
  NAMELIST /SYSTEM/ A, B, C, nr1, nr2, nr3, nat, ntyp, ibrav

  ! ELECTRONS
  CHARACTER(56) :: KS_Solve
  CHARACTER(56) :: cg_beta
  INTEGER :: electron_maxstep
  REAL(8) :: mixing_beta
  CHARACTER(56) :: diagonalization
  NAMELIST /ELECTRONS/ KS_Solve, cg_beta, electron_maxstep, mixing_beta, diagonalization

  ! ATOMIC_SPECIES
  CHARACTER(5), ALLOCATABLE :: species(:)
  REAL(8), ALLOCATABLE :: masses(:)
  CHARACTER(128), ALLOCATABLE :: pp_name(:)

  ! ATOMIC_POSITIONS
  CHARACTER(5), ALLOCATABLE :: in_atmsymb(:)
  REAL(8), ALLOCATABLE :: in_pos(:,:)


END MODULE 

