!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! modified by Fadjar Fathurrahman for ffr-LFDFT (2017)
!
!-----------------------------------------------------------------------
SUBROUTINE calc_Ewald_qe()
!-----------------------------------------------------------------------
  !
  ! Calculates Ewald energy with both G- and R-space terms.
  ! Determines optimal alpha. Should hopefully work for any structure.
  !
  USE m_constants, ONLY : tpi => TWOPI
  USE m_atoms, ONLY : nat => Natoms, &
                      ntyp => Nspecies, &
                      ityp => atm2species, &
                      zv => AtomicValences, &
                      tau => AtomicCoords, &
                      strf => StructureFactor
  USE m_LF3d, ONLY : gg => LF3d_G2, &
                     ngm => LF3d_Npoints, &
                     LL => LF3d_LL, &
                     NN => LF3d_NN
  USE m_energies, ONLY : E_nn
  IMPLICIT NONE
  INTEGER, PARAMETER :: DP=8
  !
  !   first the dummy variables
  !
  INTEGER :: gstart
  ! input: number of atoms in the unit cell
  ! input: number of different types of atoms
  ! input: the type of each atom
  ! input: number of plane waves for G sum
  ! input: first non-zero G vector

  LOGICAL :: gamma_only

  REAL(DP) :: at(3,3), bg(3,3), omega, alat, gcutm
  ! input: the direct lattice vectors
  ! input: the reciprocal lattice vectors
  ! input: the volume of the unit cell
  ! input: lattice parameter
  ! input: cut-off of g vectors
  !
  !    here the local variables
  !
  INTEGER, PARAMETER :: mxr = 50
  ! the maximum number of R vectors included in r
  INTEGER :: ng, nr, na, nb, nt, nrm
  ! counter over reciprocal G vectors
  ! counter over direct vectors
  ! counter on atoms
  ! counter on atoms
  ! counter on atomic types
  ! number of R vectors included in r sum

  REAL(DP) :: charge, ewaldg, ewaldr, dtau (3), alpha, &
       r (3, mxr), r2 (mxr), rmax, rr, upperbound, fact
  ! total ionic charge in the cell
  ! length in reciprocal space
  ! ewald energy computed in reciprocal space
  ! ewald energy computed in real space
  ! the difference tau_s - tau_s'
  ! alpha term in ewald sum
  ! input of the rgen routine ( not used here )
  ! the square modulus of R_j-tau_s-tau_s'
  ! the maximum radius to consider real space sum
  ! buffer variable
  ! used to optimize alpha
  COMPLEX(DP) :: rhon
!  REAL(DP), EXTERNAL :: qe_erfc

  ! setup at and bg
  at(:,:) = 0.d0
  at(1,1) = LL(1)
  at(2,2) = LL(2)
  at(3,3) = LL(3)

  omega = LL(1)*LL(2)*LL(3)

!  WRITE(*,*) 'TPI = ', TPI
!  WRITE(*,*)

  bg(:,:) = 0.d0
  bg(1,1) = TPI/LL(1)
  bg(2,2) = TPI/LL(2)
  bg(3,3) = TPI/LL(3)

  gcutm = maxval( NN )*TPI  !! ???? FIXME

  alat = 1.d0
  gamma_only = .FALSE.
  gstart = 2

  charge = 0.d0
  DO na = 1, nat
     charge = charge + zv( ityp(na) )
  ENDDO
  alpha = 2.9d0
100 alpha = alpha - 0.1d0
  !
  ! choose alpha in order to have convergence in the sum over G
  ! upperbound is a safe upper bound for the error in the sum over G
  !
  IF( alpha <= 0.d0) THEN 
    WRITE(*,*) 'ERROR in calculating Ewald energy:'
    WRITE(*,*) 'optimal alpha not found'
    STOP 
  ENDIF 
  !
  ! beware of unit of gcutm
  upperbound = 2.d0*charge**2*sqrt(2.d0*alpha/tpi) * erfc(sqrt(gcutm/4.d0/alpha))  
  IF(upperbound > 1.0d-7) GOTO 100
  !
  ! G-space sum here.
  ! Determine if this processor contains G=0 and set the constant term
  !
  IF(gstart==2) THEN 
    ewaldg = - charge**2 / alpha / 4.0d0
  ELSE 
    ewaldg = 0.0d0
  ENDIF 

  ! gamma_only should be .FALSE. for our case
  IF(gamma_only) THEN 
    fact = 2.d0
  ELSE
    fact = 1.d0
  ENDIF 

  DO ng = gstart, ngm
    rhon = (0.d0, 0.d0)
    DO nt = 1, ntyp
      rhon = rhon + zv(nt)*CONJG(strf(ng, nt))
    ENDDO
    ewaldg = ewaldg + fact*abs(rhon) **2 * exp( -gg(ng)/alpha/4.d0 )/ gg(ng)
  ENDDO
  ewaldg = 2.d0 * tpi / omega * ewaldg
  !
  !  Here add the other constant term
  !
  IF (gstart==2) THEN 
    DO na = 1, nat
      ewaldg = ewaldg - zv(ityp(na))**2 * sqrt(8.d0/tpi*alpha)
    ENDDO
  ENDIF 
  !
  ! R-space sum here (only for the processor that contains G=0)
  !
  ewaldr = 0.d0
  IF( gstart==2 ) THEN 
    rmax = 4.d0 / sqrt(alpha) / alat
    !
    ! with this choice terms up to ZiZj*erfc(4) are counted (erfc(4)=2x10^-8
    !
    DO na = 1, nat
      DO nb = 1, nat
        dtau(:) = tau(:,na) - tau(:,nb)
        !
        ! generates nearest-neighbors shells
        !
        CALL rgen( dtau, rmax, mxr, at, bg, r, r2, nrm )
        !
        ! and sum to the real space part
        !
        DO nr = 1, nrm
          rr = sqrt (r2 (nr) ) * alat
          ewaldr = ewaldr + zv (ityp (na) ) * zv (ityp (nb) ) * erfc( sqrt (alpha) * rr) / rr
        ENDDO 
      ENDDO 
    ENDDO 
  ENDIF 
  
  E_nn = 0.5d0*(ewaldg + ewaldr)

END SUBROUTINE 

