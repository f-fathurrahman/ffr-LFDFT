! ffr
! This is intended to work for periodic case

SUBROUTINE test_write_xsf()
  USE m_LF1d, ONLY : eval_LF1d_p
  USE m_globals, ONLY : LF, Npoints, rho, ATOMS
  IMPLICIT NONE
  REAL(8), ALLOCATABLE :: dat(:,:)
  REAL(8), ALLOCATABLE :: gridx(:), gridy(:), gridz(:)
  INTEGER :: Nx, Ny, Nz
  INTEGER :: ix, iy, iz, ip, jp
  REAL(8) :: x0, y0, z0, Lx, Ly, Lz
  REAL(8), ALLOCATABLE :: rho_interp(:)
  REAL(8) :: dVol, xx, yy, zz
  REAL(8) :: latvec(3,3)

  Nx = LF%LFx%N
  Ny = LF%LFy%N
  Nz = LF%LFz%N

  ALLOCATE( gridx(Nx) )
  ALLOCATE( gridy(Ny) )
  ALLOCATE( gridz(Nz) )

  ! Assumption: Npoints = Nx * Ny * Nz
  ALLOCATE( dat(3,Npoints) )

  x0 = LF%LFx%A
  y0 = LF%LFy%A
  z0 = LF%LFz%A
  !
  Lx = LF%LFx%L
  Ly = LF%LFy%L
  Lz = LF%LFz%L
  !
  latvec(:,:) = 0.d0
  latvec(1,1) = Lx
  latvec(2,2) = Ly
  latvec(3,3) = Lz
  !
  ip = 0
  !WRITE(*,*) 'Nx, Ny, Nz:', Nx, Ny, Nz
  !WRITE(*,*) 'Lx, Ly, Lz:', Lx, Ly, Lz
  !
  ! FIXME It does not necessarily like this
  DO iz = 1, Nz
  DO iy = 1, Ny
  DO ix = 1, Nx
    ip = ip + 1
    dat(1,ip) = x0 + (ix-1)*Lx/Nx  ! end points are not included
    dat(2,ip) = y0 + (iy-1)*Ly/Ny
    dat(3,ip) = z0 + (iz-1)*Lz/Nz
    !WRITE(*,'(1x,I8,3F18.10)') ip, dat(1:3,ip)
  ENDDO
  ENDDO
  ENDDO

  ALLOCATE( rho_interp(Npoints) )
  dVol = LF%LFx%h * LF%LFy%h * LF%LFz%h
  ! BEWARE: The following loop is VERY time consuming, especially for large data
!  WRITE(*,*) 'Interpolating data ...'
!  DO ip = 1, Npoints ! FIXME it does not have to be Npoints
!    xx = dat(1,ip)
!    yy = dat(2,ip)
!    zz = dat(3,ip)
!    rho_interp(ip) = 0.d0
!    WRITE(*,*) 'Point:', ip
!    DO jp = 1, Npoints
!      ix = LF%lin2xyz(1,jp)
!      iy = LF%lin2xyz(2,jp)
!      iz = LF%lin2xyz(3,jp)
!      ! FIXME we take it to be periodic case
!      rho_interp(ip) = rho_interp(ip) + rho(jp)*sqrt(dVol)* eval_LF1d_p(LF%LFx,ix,xx) * &
!          eval_LF1d_p(LF%LFy,iy,yy) * eval_LF1d_p(LF%LFz,iz,zz)
!    ENDDO
!    !WRITE(*,'(1x,I8,F18.10)') ip, rho_interp(ip)
!  ENDDO
!  WRITE(*,*) 'Done interpolating data'

  OPEN(unit=133,file='rho.xsf',form='formatted')  ! FIXME unit=133 is HARDCODED
  CALL xsf_struct (1.d0, latvec, ATOMS%Natoms, ATOMS%positions, &
     ATOMS%atmSymb, 133)
  ! This is appropriate for FFT grid (periodic data)
  !CALL xsf_fast_datagrid_3d(rho_interp, LF%LFx%N, LF%LFy%N, LF%LFz%N, &
  !   LF%LFx%N, LF%LFy%N, LF%LFz%N, latvec, 1.d0, 133)
  CALL xsf_datagrid_3d(rho, Nx, Ny, Nz, 1.d0, 1.d0, 1.d0, LF%lingrid(:,1), &
     latvec(1,:), latvec(2,:), latvec(3,:), 1.d0, 133)
  CLOSE(133)

  DEALLOCATE( rho_interp )
  DEALLOCATE( dat )
END SUBROUTINE


