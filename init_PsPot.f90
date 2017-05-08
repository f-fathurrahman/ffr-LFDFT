SUBROUTINE init_PsPot()

  USE m_atoms, ONLY : SpeciesSymbols, Nspecies, AtomicValences, &
                      atm2species, Natoms
  USE m_PsPot
  USE m_Ps_HGH

  IMPLICIT NONE 
  INTEGER :: isp, NprojTot, iprj, ii, l, ia

  ! Use default
  ALLOCATE( PsPot_FilePath(Nspecies) )
  ALLOCATE( Ps_HGH_Params(Nspecies) )

  WRITE(*,*) 'Pass here ...'

  iprj = 0
  DO isp = 1,Nspecies

    PsPot_FilePath(isp) = trim(PsPot_Dir) // trim(SpeciesSymbols(isp)) // '.hgh'
    CALL init_Ps_HGH_Params( Ps_HGH_Params(isp), PsPot_FilePath(isp) )
    
    AtomicValences(isp) = Ps_HGH_Params(isp)%zval

  ENDDO 

  NprojTotMax = 0
  DO isp=1,Nspecies
    NprojTot = sum( Ps_HGH_Params(isp)%Nproj_l )
    WRITE(*,*) Ps_HGH_Params(isp)%Nproj_l
    WRITE(*,*) 'isp, NprojTot = ', isp, NprojTot
    IF( NprojTotMax < NprojTot ) NprojTotMax = NprojTot
  ENDDO 
  WRITE(*,*) 'NprojTotMax = ', NprojTotMax

  ALLOCATE( w_NL(Nspecies,NprojTotMax) )
  w_NL(:,:) = 0.d0
  DO isp=1,Nspecies
    DO l = 0,Ps_HGH_Params(isp)%lmax
      DO ii = 1,Ps_HGH_Params(isp)%Nproj_l(l)
        iprj = iprj + 1
        w_NL(isp,iprj) = Ps_HGH_Params(isp)%h(l,ii,ii)
      ENDDO 
    ENDDO 
  ENDDO 
 
  NbetaNL = 0
  DO ia = 1,Natoms
    isp = atm2species(ia)
    NbetaNL = NbetaNL + Ps_HGH_Params(isp)%Nproj_l(0)*1 + &
                        Ps_HGH_Params(isp)%Nproj_l(1)*3 + &
                        Ps_HGH_Params(isp)%Nproj_l(2)*5 + &
                        Ps_HGH_Params(isp)%Nproj_l(3)*7
  ENDDO 
  WRITE(*,*) 'NbetaNL = ', NbetaNL

END SUBROUTINE 

