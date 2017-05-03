SUBROUTINE init_ilu0_prec()

  USE m_ilu0_prec
  USE m_nabla2_sparse, ONLY : nzval => nabla2_nzval, &
                              rowval => nabla2_rowval, &
                              colptr => nabla2_colptr
  USE m_LF3d, ONLY : Npoints => LF3d_Npoints, &
                     NN => LF3d_NN

  IMPLICIT NONE 
  INTEGER :: Nx, Ny, Nz
  INTEGER :: ierr

  Nx = NN(1)
  Ny = NN(2)
  Nz = NN(3)

  WRITE(*,*)
  WRITE(*,*) 'ILU0 preconditioner is based on kinetic matrix'

  ALLOCATE( alu_ilu0(Npoints*(Nx+Ny+Nz-2)) )
  ALLOCATE( jlu_ilu0(Npoints*(Nx+Ny+Nz-2)) )
  ALLOCATE( ju_ilu0(Npoints) )
  ALLOCATE( iw_ilu0(Npoints) )
  !
  CALL ilu0( Npoints, -0.5d0*nzval, rowval, colptr, alu_ilu0, jlu_ilu0, ju_ilu0, iw_ilu0, ierr )

END SUBROUTINE 

