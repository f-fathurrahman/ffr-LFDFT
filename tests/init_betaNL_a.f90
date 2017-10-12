SUBROUTINE init_betaNL_a()

  USE m_LF3d, ONLY : LL => LF3d_LL, &
                     dVol => LF3d_dVol, &
                     LF3d_TYPE, LF3d_PERIODIC
  USE m_PsPot, ONLY : betaNL, NbetaNL, &
                      Ps => Ps_HGH_Params
  USE m_atoms, ONLY : atpos => AtomicCoords, Natoms, atm2species
  USE m_Ps_HGH, ONLY : hgh_eval_proj_R
  USE m_grid_atom_cube
  IMPLICIT NONE
  INTEGER :: ia, isp, l, m, iprj
  INTEGER :: Np_beta, ip, ibeta
  REAL(8) :: dr_vec(3)
  REAL(8) :: dr
  REAL(8) :: Ylm_real
  REAL(8) :: nrm

  ALLOCATE( betaNL(Npoints_a,NbetaNL) )

  WRITE(*,*)
  WRITE(*,*) 'Initializing betaNL functions'


  ! loop structure must be the same as in init_PsPot
  ibeta = 0
  DO ia = 1,Natoms
    isp = atm2species(ia)
    DO l = 0,Ps(isp)%lmax
      DO iprj = 1,Ps(isp)%Nproj_l(l)
        DO m = -l,l
          ibeta = ibeta + 1
          Np_beta = 0
          DO ip = 1,Npoints_a
            !
            IF( LF3d_TYPE == LF3d_PERIODIC ) THEN 
              ! FIXME This will only works for atpos within the cell ??
              CALL calc_dr_periodic_1pnt( LL, atpos(:,ia), grid_a(:,ip), dr_vec )
              dr = sqrt( dr_vec(1)**2 + dr_vec(2)**2 + dr_vec(3)**2 )
            ELSE 
              ! Note that calc_dr_1pnt already took the square root
              CALL calc_dr_1pnt( atpos(:,ia), grid_a(:,ip), dr )
            ENDIF 
            !
            IF( dr <= Ps(isp)%rcut_NL(l) ) THEN
              Np_beta = Np_beta + 1
              betaNL(ip,ibeta) = hgh_eval_proj_R( Ps(isp), l, iprj, dr ) * Ylm_real( l, m, dr_vec )
            ENDIF
          ENDDO
          nrm = sum(betaNL(:,ibeta)**2)*dVol_a
          WRITE(*,'(1x,A,I5,I8,F18.10)') 'ibeta, Np_beta, integ = ', ibeta, Np_beta, nrm
        ENDDO ! iprj
      ENDDO ! m
    ENDDO ! l
  ENDDO

  flush(6)

END SUBROUTINE
