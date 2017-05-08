SUBROUTINE init_betaNL()

  USE m_LF3d, ONLY : Npoints => LF3d_Npoints, &
                     lingrid => LF3d_lingrid, &
                     LL => LF3d_LL
  USE m_PsPot, ONLY : betaNL, NbetaNL, &
                      Ps => Ps_HGH_Params
  USE m_atoms, ONLY : atpos => AtomicCoords, Natoms, atm2species
  USE m_Ps_HGH, ONLY : hgh_eval_proj_R
  IMPLICIT NONE 
  INTEGER :: ia, isp, l, m, iprj
  INTEGER :: Np_beta, ip, ibeta
  REAL(8) :: dr_vec(3)
  REAL(8) :: dr
  REAL(8) :: Ylm_real
  

  ALLOCATE( betaNL(Npoints,NbetaNL) )

  ibeta = 0
  DO ia = 1,Natoms
    isp = atm2species(ia)
    DO l = 0,Ps(isp)%lmax
      DO iprj = 1,Ps(isp)%Nproj_l(l)
        DO m = -l,l
          ibeta = ibeta + 1
          Np_beta = 0
          DO ip = 1,Npoints
            CALL calc_dr_periodic_1pnt( LL, atpos(:,ia), lingrid(:,ip), dr_vec )
            dr = sqrt( dr_vec(1)**2 + dr_vec(2)**2 + dr_vec(3)**2 )
            IF( dr <= Ps(isp)%rcut_NL(l) ) THEN 
              Np_beta = Np_beta + 1
              betaNL(ip,ibeta) = hgh_eval_proj_R( Ps(isp), l, iprj, dr ) * Ylm_real( l, m, dr_vec )
            ENDIF 
          ENDDO 
          WRITE(*,'(1x,A,I5,I8)') 'ibeta, Np_beta = ', ibeta, Np_beta
        ENDDO ! m
      ENDDO ! iprj
    ENDDO ! l
  ENDDO 

  WRITE(*,*) 'sum(betaNL) = ', sum(betaNL)

END SUBROUTINE 

