SUBROUTINE init_PsPot()

  USE m_atoms, ONLY : SpeciesSymbols, Nspecies, AtomicValences
  USE m_PsPot
  USE m_Ps_HGH

  IMPLICIT NONE 
  INTEGER :: isp

  ! Use default
  ALLOCATE( PsPot_FilePath(Nspecies) )
  ALLOCATE( Ps_HGH_Params(Nspecies) )

  DO isp = 1,Nspecies

    PsPot_FilePath(isp) = trim(PsPot_Dir) // trim(SpeciesSymbols(isp)) // '.hgh'
    !WRITE(*,*) isp, PsPot_FilePath(isp)

    CALL init_Ps_HGH_Params( Ps_HGH_Params(isp), PsPot_FilePath(isp) )
    
    AtomicValences(isp) = Ps_HGH_Params(isp)%zval

  ENDDO 

END SUBROUTINE 

