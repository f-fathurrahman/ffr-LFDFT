MODULE m_options

  IMPLICIT NONE 

  ! Options for controlling how beta is calculated
  INTEGER :: CG_BETA = 2
  ! 1 => Fletcher-Reeves
  ! 2 => Polak-Ribiere
  ! 3 => Hestenes-Stiefel
  ! 4 => Dai-Yuan

  ! whether free nabla2 after constructing Laplacian matrix or not
  ! set to .TRUE. to minimize memory usage
  LOGICAL :: FREE_NABLA2 = .TRUE.

  REAL(8) :: ETHR_EVALS = 1.d-3
  REAL(8) :: ETHR_EVALS_LAST = 1.0d-10
  INTEGER :: IALG_DIAG = 1
  ! 1 => Davidson v1 (from PWSCF)
  ! 2 => Davidson v2
  ! 3 => LOBPCG

  ! type of mixing to use for the potential
  integer :: mixtype
  ! mixing type description
  character(256) :: mixdescr
  ! adaptive mixing parameter
  REAL(8) :: beta0
  REAL(8) :: betamax
  ! subspace dimension for Broyden mixing
  INTEGER :: mixsdb
  ! Broyden mixing parameters alpha and w0
  REAL(8) :: broydpm(2)


END MODULE 

