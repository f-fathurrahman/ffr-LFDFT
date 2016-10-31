!
! Copyright (C) 2003 Tone Kokalj
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! This file holds XSF (=Xcrysden Structure File) utilities.
! Routines written by Tone Kokalj on Mon Jan 27 18:51:17 CET 2003
!
! minor modification by ffr

! -------------------------------------------------------------------
!   this routine writes the crystal structure in XSF format
! -------------------------------------------------------------------
SUBROUTINE xsf_struct (alat, at, nat, tau, atm, ounit)
  USE m_constants, ONLY : ANG2BOHR
  IMPLICIT NONE
  INTEGER, PARAMETER :: DP = 8
  INTEGER          :: nat, ounit
  CHARACTER(len=*) :: atm(*)  ! atomic symbol
  real(DP)    :: alat, tau (3, nat), at (3, 3)
  ! --
  INTEGER          :: i, j, n
  real(DP)    :: at1 (3, 3)
  ! convert lattice vectors to ANGSTROM units ...
  DO i=1,3
     DO j=1,3
        at1(j,i) = at(j,i)*alat/ANG2BOHR
     ENDDO
  ENDDO

  WRITE(ounit,*) 'CRYSTAL'
  WRITE(ounit,*) 'PRIMVEC'
  WRITE(ounit,'(2(3F15.9/),3f15.9)') at1
  WRITE(ounit,*) 'PRIMCOORD'
  WRITE(ounit,*) nat, 1

  DO n=1,nat
     write(*,*) n, trim(atm(n))
     ! positions are in Angstroms
     WRITE(ounit,'(a3,3x,3f15.9)') atm(n), &
          tau(1,n)*alat/ANG2BOHR, &
          tau(2,n)*alat/ANG2BOHR, &
          tau(3,n)*alat/ANG2BOHR
  ENDDO
  RETURN
END SUBROUTINE xsf_struct



! -------------------------------------------------------------------
!   this routine writes the 3D scalar field (i.e. uniform mesh of points)
!   in XSF format using the FFT mesh (i.e. fast write)
! -------------------------------------------------------------------
SUBROUTINE xsf_fast_datagrid_3d &
     (rho, nr1, nr2, nr3, nr1x, nr2x, nr3x, at, alat, ounit)
  USE m_constants, ONLY : ANG2BOHR
  IMPLICIT NONE
  INTEGER, PARAMETER :: DP=8
  INTEGER       :: nr1x, nr2x, nr3x, nr1, nr2, nr3, ounit
  real(DP) :: alat, at (3, 3), rho(nr1x,nr2x,nr3x)
  ! --
  INTEGER       :: i1, i2, i3, ix, iy, iz, count, i, &
       ind_x(10), ind_y(10),ind_z(10)

  ! XSF scalar-field header
  WRITE(ounit,'(a)') 'BEGIN_BLOCK_DATAGRID_3D'
  WRITE(ounit,'(a)') '3D_PWSCF'
  WRITE(ounit,'(a)') 'DATAGRID_3D_UNKNOWN'

  ! number of points in each direction
  WRITE(ounit,*) nr1+1, nr2+1, nr3+1
  ! origin
  WRITE(ounit,'(3f10.6)') 0.0d0, 0.0d0, 0.0d0
  ! 1st spanning (=lattice) vector
  WRITE(ounit,'(3f12.6)') (alat*at(i,1)/ANG2BOHR,i=1,3) ! in ANGSTROMS
  ! 2nd spanning (=lattice) vector
  WRITE(ounit,'(3f12.6)') (alat*at(i,2)/ANG2BOHR,i=1,3)
  ! 3rd spanning (=lattice) vector
  WRITE(ounit,'(3f12.6)') (alat*at(i,3)/ANG2BOHR,i=1,3)

  count=0
  DO i3=0,nr3
     !iz = mod(i3,nr3)
     iz = mod(i3,nr3) + 1

     DO i2=0,nr2
        !iy = mod(i2,nr2)
        iy = mod(i2,nr2) + 1

        DO i1=0,nr1
           !ix = mod(i1,nr1)
           ix = mod(i1,nr1) + 1

           !ii = (1+ix) + iy*nr1x + iz*nr1x*nr2x
           IF (count<6) THEN
              count = count + 1
              !ind(count) = ii
           ELSE
              WRITE(ounit,'(6e13.5)') &
                   (rho(ind_x(i),ind_y(i),ind_z(i)),i=1,6)
              count=1
              !ind(count) = ii
           ENDIF
           ind_x(count) = ix
           ind_y(count) = iy
           ind_z(count) = iz
        ENDDO
     ENDDO
  ENDDO
  WRITE(ounit,'(6e13.5:)') (rho(ind_x(i),ind_y(i),ind_z(i)),i=1,count)
  WRITE(ounit,'(a)') 'END_DATAGRID_3D'
  WRITE(ounit,'(a)') 'END_BLOCK_DATAGRID_3D'
  RETURN
END SUBROUTINE xsf_fast_datagrid_3d




SUBROUTINE xsf_datagrid_2d (rho, nx, ny, m1, m2, x0, e1, e2, alat, ounit)
  USE m_constants, ONLY : ANG2BOHR
  IMPLICIT NONE
  INTEGER, PARAMETER :: DP=8
  INTEGER       :: nx, ny, ounit
  real(DP) :: m1, m2, alat, x0(3), e1(3), e2(3), rho(2, nx, ny)
  ! --
  INTEGER       :: ix, iy, count, i, ind_x(10), ind_y(10)

  ! XSF scalar-field header
  WRITE(ounit,'(a)') 'BEGIN_BLOCK_DATAGRID_2D'
  WRITE(ounit,'(a)') '2D_PWSCF'
  WRITE(ounit,'(a)') 'DATAGRID_2D_UNKNOWN'

  ! number of points in each direction
  WRITE(ounit,*) nx, ny
  ! origin
  WRITE(ounit,'(3f10.6)') (alat*x0(i)/ANG2BOHR,i=1,3) ! in ANGSTROMS
  ! 1st spanning (=lattice) vector
  WRITE(ounit,'(3f10.6)') (alat*e1(i)*m1/ANG2BOHR,i=1,3) ! in ANGSTROMS
  ! 2nd spanning (=lattice) vector
  WRITE(ounit,'(3f10.6)') (alat*e2(i)*m2/ANG2BOHR,i=1,3) ! in ANGSTROMS

  count=0
  DO iy=1,ny
     DO ix=1,nx
        IF (count < 6) THEN
           count = count + 1
        ELSE
           WRITE(ounit,'(6e13.5)') (rho(1,ind_x(i),ind_y(i)),i=1,6)
           count=1
        ENDIF
        ind_x(count) = ix
        ind_y(count) = iy
     ENDDO
  ENDDO

  WRITE(ounit,'(6e13.5:)') (rho(1,ind_x(i),ind_y(i)),i=1,count)
  WRITE(ounit,'(a)') 'END_DATAGRID_2D'
  WRITE(ounit,'(a)') 'END_BLOCK_DATAGRID_2D'
  RETURN
END SUBROUTINE xsf_datagrid_2d



SUBROUTINE xsf_datagrid_3d &
     (rho, nx, ny, nz, m1, m2, m3, x0, e1, e2, e3, alat, ounit)
  USE m_constants, ONLY : ANG2BOHR
  IMPLICIT NONE
  INTEGER, PARAMETER :: DP=8
  INTEGER       :: nx, ny, nz, ounit
  real(DP) :: m1, m2, m3, alat, x0(3), e1(3),e2(3),e3(3), rho(nx, ny, nz)
  ! --
  INTEGER       :: ix, iy, iz, count, i, ind_x(10), ind_y(10), ind_z(10)

  ! XSF scalar-field header
  WRITE(ounit,'(a)') 'BEGIN_BLOCK_DATAGRID_3D'
  WRITE(ounit,'(a)') '3D_PWSCF'
  WRITE(ounit,'(a)') 'DATAGRID_3D_UNKNOWN'

  ! number of points in each direction
  WRITE(ounit,*) nx, ny, nz
  ! origin
  WRITE(ounit,'(3f10.6)') (alat*x0(i)/ANG2BOHR,i=1,3) ! in ANGSTROMS
  ! 1st spanning (=lattice) vector
  WRITE(ounit,'(3f10.6)') (alat*e1(i)*m1/ANG2BOHR,i=1,3) ! in ANGSTROMS
  ! 2nd spanning (=lattice) vector
  WRITE(ounit,'(3f10.6)') (alat*e2(i)*m2/ANG2BOHR,i=1,3) ! in ANGSTROMS
  ! 3rd spanning (=lattice) vector
  WRITE(ounit,'(3f10.6)') (alat*e3(i)*m3/ANG2BOHR,i=1,3)

  count=0
  DO iz=1,nz
     DO iy=1,ny
        DO ix=1,nx
           IF (count<6) THEN
              count = count + 1
           ELSE
              WRITE(ounit,'(6e13.5)') (rho(ind_x(i),ind_y(i),ind_z(i)),i=1,6)
              count=1
           ENDIF
           ind_x(count) = ix
           ind_y(count) = iy
           ind_z(count) = iz
        ENDDO
     ENDDO
  ENDDO

  WRITE(ounit,'(6e13.5:)') (rho(ind_x(i),ind_y(i),ind_z(i)),i=1,count)
  WRITE(ounit,'(a)') 'END_DATAGRID_3D'
  WRITE(ounit,'(a)') 'END_BLOCK_DATAGRID_3D'
  RETURN
END SUBROUTINE xsf_datagrid_3d

