SUBROUTINE init_PsPot()

  USE m_atoms, ONLY : SpeciesSymbols, Nspecies, AtomicValences, &
                      atm2species, Natoms
  USE m_PsPot, ONLY : Ps_HGH_Params, PsPot_FilePath, PsPot_Dir, &
                      NbetaNL, prj2beta
  USE m_Ps_HGH, ONLY : init_Ps_HGH_Params

  IMPLICIT NONE 
  INTEGER :: isp, iprj, l, ia, m

  ! Use default
  ALLOCATE( PsPot_FilePath(Nspecies) )
  ALLOCATE( Ps_HGH_Params(Nspecies) )

  DO isp = 1,Nspecies
    ! initialize HGH pseudopotentials for each atomic species
    PsPot_FilePath(isp) = trim(PsPot_Dir) // trim(SpeciesSymbols(isp)) // '.hgh'
    CALL init_Ps_HGH_Params( Ps_HGH_Params(isp), PsPot_FilePath(isp) )
    ! Set atomic valences
    AtomicValences(isp) = Ps_HGH_Params(isp)%zval
  ENDDO 

  ALLOCATE( prj2beta(1:3,1:Natoms,0:3,-3:3) )
  prj2beta(:,:,:,:) = -1
  NbetaNL = 0
  DO ia = 1,Natoms
    isp = atm2species(ia)
    DO l = 0,Ps_HGH_Params(isp)%lmax
      DO iprj = 1,Ps_HGH_Params(isp)%Nproj_l(l)
        DO m = -l,l
          NbetaNL = NbetaNL + 1
          prj2beta(iprj,ia,l,m) = NbetaNL
        ENDDO ! m
      ENDDO ! iprj
    ENDDO ! l
  ENDDO 

END SUBROUTINE 

