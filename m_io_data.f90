! various subroutines to write binary data

MODULE m_io_data

  IMPLICIT NONE 

  INTEGER, PARAMETER :: IU_EVECS=101
  INTEGER, PARAMETER :: IU_EVALS=102
  INTEGER, PARAMETER :: IU_RHO=103
  INTEGER, PARAMETER :: IU_KSPOT=104

END MODULE 

! FIXME: only adapted for periodic grid !
SUBROUTINE write_data3d_xsf( dat, filexsf )
  USE m_io_data
  USE m_LF3d, ONLY : Npoints => LF3d_Npoints, &
                     NN => LF3d_NN, &
                     AA => LF3d_AA, BB => LF3d_BB, &
                     origin => LF3d_GRID_SHIFT
  USE m_atoms, ONLY : atpos => AtomicCoords, &
                      SpeciesSymbols, atm2species, Natoms
  USE m_constants, ONLY : ANG2BOHR
  !
  IMPLICIT NONE 
  !
  CHARACTER(*) :: filexsf
  REAL(8) :: dat(Npoints)
  INTEGER, PARAMETER :: unitxsf = 123
  REAL(8) :: LatVecs(3,3)

  LatVecs(:,:) = 0.d0
  LatVecs(1,1) = BB(1) - AA(1)
  LatVecs(2,2) = BB(2) - AA(2)
  LatVecs(3,3) = BB(3) - AA(3)
  WRITE(*,*) 'LatVecs = ', LatVecs(1,1)
  WRITE(*,*) 'LatVecs = ', LatVecs(2,2)
  WRITE(*,*) 'LatVecs = ', LatVecs(3,3)

  ! conversion to angstrom is done in xsf_* subroutines

  OPEN(unit=unitxsf, file=filexsf, form='formatted')
  CALL xsf_struct( LatVecs, Natoms, atpos, SpeciesSymbols, atm2species, unitxsf )
  CALL xsf_fast_datagrid_3d( dat, NN(1), NN(2), NN(3), NN(1), NN(2), NN(3), &
            origin, LatVecs, unitxsf)
  CLOSE(unitxsf)
END SUBROUTINE 


! FIXME This should be more delicate than this !
!
! Npoints and Nstates that are read from file should be compared
! to input data
!
SUBROUTINE read_KS_evecs(filname)
  USE m_io_data
  USE m_LF3d, ONLY : Npoints => LF3d_Npoints
  USE m_states, ONLY : KS_evecs, Nstates
  IMPLICIT NONE 
  CHARACTER(*) :: filname
  INTEGER :: in_Npoints, in_Nstates

  OPEN(unit=IU_EVECS, file=filname , action='read', form='unformatted')
  
  READ(IU_EVECS) in_Npoints
  READ(IU_EVECS) in_Nstates
  WRITE(*,*) 'in_Npoints, in_Npoints = ', in_Npoints, in_Nstates
  
  IF( allocated(KS_evecs) ) THEN 
    IF( size(KS_evecs,1) /= in_Npoints .OR. size(KS_evecs,2) /= in_Nstates ) THEN 
      WRITE(*,*)
      WRITE(*,*) 'ERROR reading KS_evecs: size not match'
      WRITE(*,'(1x,A,2I8)') 'Read size: ', in_Npoints, in_Nstates
      WRITE(*,*) 'Allocated size: ', Npoints, Nstates
      STOP 
    ENDIF 
  ELSE 
    Npoints = in_Npoints
    Nstates = in_Nstates
    ALLOCATE( KS_evecs(Npoints,Nstates) )
  ENDIF 

  READ(IU_EVECS) KS_evecs
  CLOSE(IU_EVECS)
END SUBROUTINE 


SUBROUTINE write_KS_evecs(filname)
  USE m_io_data
  USE m_LF3d, ONLY : Npoints => LF3d_Npoints
  USE m_states, ONLY : KS_evecs, Nstates
  IMPLICIT NONE 
  CHARACTER(*) :: filname
  INTEGER, PARAMETER :: IU_WFC = 55

  OPEN(unit=IU_EVECS, file=filname , action='write', form='unformatted')
  WRITE(IU_EVECS) Npoints
  WRITE(IU_EVECS) Nstates
  WRITE(IU_EVECS) KS_evecs
  CLOSE(IU_EVECS)
END SUBROUTINE 


