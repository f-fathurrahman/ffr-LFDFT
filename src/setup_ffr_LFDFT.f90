!!>
!!> \section{Subroutine \texttt{setup\_ffr\_LFDFT()}}
!!>
!!> This subroutine prepares various tasks before actually solving the Kohn-Sham
!!> equation.
!!>
!!>
SUBROUTINE setup_ffr_LFDFT()

  USE m_options, ONLY : FREE_NABLA2, I_POISSON_SOLVE
  USE m_input_vars, ONLY : assume_isolated
  USE m_LF3d, ONLY : Npoints => LF3d_Npoints
  USE m_states, ONLY : Nstates, &
                       evals => KS_evals, &
                       evecs => KS_evecs 
  !
  INTEGER :: Narg   ! number of argument
  INTEGER :: iargc  ! needed for several compilers
  CHARACTER(64) :: filein

!!> \begin{itemize}
!!>
!!> \item
!!> We first check the number of argument(s) give to the program using built-in
!!> function \texttt{iargc()} and save the result to variable \texttt{Narg}.
!!> Currently, we only support one argument, i.e. path to input file.
!!> The program will stop and display error message if \texttt{Narg}
!!> is not equal to one.
  Narg = iargc()
  IF( Narg /= 1 ) THEN
    WRITE(*,*)
    WRITE(*,*) 'ERROR: exactly one argument must be given: input file path'
    STOP
  ENDIF

!!> \item
!!> We get the actual argument using built-in subroutine \texttt{getarg()}.
  CALL getarg( 1, filein )

!!> \item
!!> We read the input file using subroutine \texttt{read\_input()}
  CALL read_input( filein )

!!> \item
!!> Various options, such as convergence criteria, choice of algorithms, etc which are
!!> given in the input file, will be converted to internal variables (mostly defined in
!!> module \texttt{m\_options}).
  CALL setup_options()

!!> \item
!!> The following subroutine will initialize global variables related to basis function
!!> and grid points, molecular or crystalline structures, and pseudopotentials.
  CALL setup_from_input()

!!> \item
!!> The following calls will output information about molecular or crystalline
!!> structures, pseudopotentials, and basis function and grid points.
  CALL info_atoms()
  CALL info_PsPot()
  CALL info_LF3d()

!!> \item
!!> This subroutine initialize nonlocal pseudopotential projectors. It
!!> must be called after variables from \texttt{m\_LF3d} are initialized as they are
!!> defined on grid points.
!!>
!!> TODO/FIXME: To be consistent, this call should be made in pseudopotential setup.
!!> (probably via subroutine \texttt{setup\_from\_input()})
  CALL init_betaNL()

!!> \item
!!> This call will determined number of occupied and unoccupied states.
!!> It will also initialize occupation numbers.
!!>
!!> Note: Memories for eigenvectors and eigenvalues of KS equations are not
!!> allocated here.
  CALL init_states()

!!> \item
!!> Allocate KS eigenvectors and eigenvalues.
!!>
!!> This step should be done in some wrapper subroutine, however, for current use-case
!!> this is sufficient.
  ALLOCATE( evecs(Npoints,Nstates), evals(Nstates) )

!!> \item
!!> This call will initialize and calculate structure factor $S_{f}(\mathbf{G})$.
!!> This is required for periodic LF.
!!>
!!> TODO/FIXME: This call should be made in \texttt{setup\_from\_input} or anther
!!> subroutine.
  CALL init_strfact_shifted()

!!> \item
!!> Ewald energy is calculated here.
!!>
!!> FIXME: Ewald energy calculation should be called everytime atomic positions
!!> are updated, for example in geometry optimization or molecular dynamics.
  IF( assume_isolated == 'sinc' ) THEN
    CALL calc_E_NN()
  ELSE
    !CALL calc_Ewald_qe()
    CALL calc_Ewald()
  ENDIF

!!> \item
!!> The following call will initialize various global variables (arrays) needed to
!!> define Hamiltonian. It is mainly used for storing potential terms.
  CALL alloc_hamiltonian()

!!> \item
!!> Allocate local pseudopotential. For periodic sinc LF the potential is constructed
!!> directly on real space grid. For periodic LF, the potential is first constructed on
!!> reciprocal space and then transformed to real space grid via inverse FFT.
  IF( assume_isolated == 'sinc' ) THEN
    CALL init_V_ps_loc()
  ELSE
    CALL init_V_ps_loc_G()
  ENDIF

!!> \item
!!> Laplacian matrix $\nabla^2$ is initialized here.
  CALL init_nabla2_sparse()

!!> \item
!!> Here we construct ILU0 preconditioner based on kinetic matrix.
  CALL init_ilu0_prec()

!!> \item
!!> If the option \texttt{FREE\_NABLA2} is \texttt{.TRUE.}, then memory for storing
!!> Laplacian matrix is freed immediately. In this way, application of kinetic operator
!!> will be done by a matrix-free algorithm instead of using sparse-matrix multiplication.
  IF( FREE_NABLA2 ) THEN
    CALL dealloc_nabla2_sparse()
  ENDIF

!!> \item
!!> This is a special setup for Poisson solver in the case of non-periodic
!!> system. Currrently there are two different methods: (1) interpolating scaling
!!> function (ISF) method and (2) Direct Algorithm for Gravitation and Electrostatic (DAGE).
  IF( I_POISSON_SOLVE == 1 ) THEN
    CALL init_Poisson_solve_ISF()
  ELSEIF( I_POISSON_SOLVE == 2 ) THEN
    CALL init_Poisson_solve_DAGE()
  ENDIF
END SUBROUTINE
!!> \end{itemize}
