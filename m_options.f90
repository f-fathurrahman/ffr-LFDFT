MODULE m_options

  IMPLICIT NONE 

  INTEGER :: I_KS_SOLVE = 1
  ! 1 => KS_solve_Emin_pcg
  ! 2 => SCF

  ! Options for controlling how beta is calculated
  INTEGER :: I_CG_BETA = 2
  ! 1 => Fletcher-Reeves
  ! 2 => Polak-Ribiere
  ! 3 => Hestenes-Stiefel
  ! 4 => Dai-Yuan

  ! whether free nabla2 after constructing Laplacian matrix or not
  ! set to .TRUE. to minimize memory usage
  LOGICAL :: FREE_NABLA2 = .TRUE.

  REAL(8) :: ETHR_EVALS = 1.d-3
  REAL(8) :: ETHR_EVALS_LAST = 1.0d-10

  INTEGER :: I_ALG_DIAG = 2
  ! 1 => Davidson v1 (from PWSCF)
  ! 2 => Davidson v2
  ! 3 => LOBPCG

  INTEGER :: Emin_NiterMax = 500
  INTEGER :: SCF_NiterMax = 100

  REAL(8) :: Emin_ETOT_CONV_THR = 1.d-6
  REAL(8) :: SCF_ETOT_CONV_THR = 1.d-6

  ! mixing beta
  REAL(8) :: SCF_betamix = 0.2d0
  ! type of mixing to use for the potential
  INTEGER :: mixtype
  ! mixing type description
  CHARACTER(256) :: mixdescr
  ! adaptive mixing parameter
  REAL(8) :: beta0
  REAL(8) :: betamax
  ! subspace dimension for Broyden mixing
  INTEGER :: mixsdb
  ! Broyden mixing parameters alpha and w0
  REAL(8) :: broydpm(2)


END MODULE 

