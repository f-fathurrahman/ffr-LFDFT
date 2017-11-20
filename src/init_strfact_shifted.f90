! Calculate structure factor
! FIXME: needed anyway even for sinc ??
SUBROUTINE init_strfact_shifted()

  USE m_atoms, ONLY : Na => Natoms, &
                      Xpos => AtomicCoords, &
                      Nspecies, &
                      atm2species, &
                      strf => StructureFactor
  USE m_LF3d, ONLY : Ng => LF3d_Npoints, &
                     Gv => LF3d_Gv, &
                     grid_x => LF3d_grid_x, &
                     grid_y => LF3d_grid_y, &
                     grid_z => LF3d_grid_z
  IMPLICIT NONE
  !
  INTEGER :: ia, isp, ig
  REAL(8) :: GX, shiftx, shifty, shiftz

  WRITE(*,*)
  WRITE(*,*) 'Calculating structure factor: shifted to Fourier grid'

  shiftx = 0.5d0*( grid_x(2) - grid_x(1) )
  shifty = 0.5d0*( grid_y(2) - grid_y(1) )
  shiftz = 0.5d0*( grid_z(2) - grid_z(1) )

  ALLOCATE( strf(Ng,Nspecies) )

  strf(:,:) = cmplx(0.d0,0.d0,kind=8)
  DO ia = 1,Na
    isp = atm2species(ia)
    DO ig = 1,Ng
      GX = (Xpos(1,ia)-shiftx)*Gv(1,ig) + (Xpos(2,ia)-shifty)*Gv(2,ig) + (Xpos(3,ia)-shiftz)*Gv(3,ig)
      strf(ig,isp) = strf(ig,isp) + cmplx( cos(GX), -sin(GX), kind=8 )
    ENDDO
  ENDDO

  CALL flush(6)

END SUBROUTINE
