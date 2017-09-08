! various subroutines to write binary data

MODULE m_io_data

  IMPLICIT NONE 

CONTAINS 

! FIXME: only adapted for periodic grid !
SUBROUTINE write_data3d_xsf( dat, filexsf )
  USE m_LF3d, ONLY : Npoints => LF3d_Npoints, &
                     NN => LF3d_NN, &
                     grid_x => LF3d_grid_x, &
                     grid_y => LF3d_grid_y, &
                     grid_z => LF3d_grid_z, &
                     AA => LF3d_AA, BB => LF3d_BB
  USE m_atoms, ONLY : atpos => AtomicCoords, &
                      SpeciesSymbols, atm2species, Natoms
  !
  IMPLICIT NONE 
  !
  CHARACTER(*) :: filexsf
  REAL(8) :: dat(Npoints)
  INTEGER, PARAMETER :: unitxsf = 123
  REAL(8) :: LatVecs(3,3)
  REAL(8) :: origin(3)

  LatVecs(:,:) = 0.d0
  LatVecs(1,1) = BB(1) - AA(1)
  LatVecs(2,2) = BB(2) - AA(2)
  LatVecs(3,3) = BB(3) - AA(3)

  origin(1) = 0.5d0*( grid_x(2) - grid_x(1) )
  origin(2) = 0.5d0*( grid_y(2) - grid_y(1) )
  origin(3) = 0.5d0*( grid_z(2) - grid_z(1) )

  WRITE(*,*) 'filexsf = ', filexsf
  OPEN( unit=unitxsf, file=filexsf )
  CALL xsf_struct( LatVecs, Natoms, atpos, SpeciesSymbols, atm2species, unitxsf )
  CALL xsf_fast_datagrid_3d( dat, NN(1), NN(2), NN(3), NN(1), NN(2), NN(3), &
            origin, LatVecs, unitxsf)
  CLOSE(unitxsf)
END SUBROUTINE 



SUBROUTINE write_KS_evecs(filname)
  USE m_states, ONLY : KS_evecs
  IMPLICIT NONE 
  CHARACTER(*) :: filname
  INTEGER, PARAMETER :: IU = 55

  OPEN( unit=55, file=filname , action='write', form='unformatted' )
  WRITE(IU) KS_evecs
  CLOSE(IU)
END SUBROUTINE 


END MODULE 

