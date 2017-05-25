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

  ! ATOMIC_SPECIES
  CHARACTER(5), ALLOCATABLE :: species(:)
  REAL(8), ALLOCATABLE :: Masses(:)
  CHARACTER(128), ALLOCATABLE :: pp_name(:)

  ! ATOMIC_POSITIONS
  CHARACTER(5), ALLOCATABLE :: in_atmsymb(:)
  REAL(8), ALLOCATABLE :: in_pos(:,:)

END MODULE 

