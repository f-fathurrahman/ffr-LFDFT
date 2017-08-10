SUBROUTINE calc_E_NN()

  USE m_atoms, ONLY : atpos => AtomicCoords, &
                      atm2species, Natoms, &
                      Zv => AtomicValences
  USE m_energies, ONLY : E_NN
  IMPLICIT NONE 
  INTEGER :: ia, ja, isp, jsp
  REAL(8) :: r, dx2, dy2, dz2

  E_NN = 0.d0
  DO ia = 1,Natoms
   isp = atm2species(ia)
   DO ja = 1,Natoms
     IF( ia /= ja ) THEN 
       jsp = atm2species(ja)
       dx2 = ( atpos(1,ia) - atpos(1,ja) )**2
       dy2 = ( atpos(2,ia) - atpos(2,ja) )**2
       dz2 = ( atpos(3,ia) - atpos(3,ja) )**2
       r = sqrt( dx2 + dy2 + dz2 )
       E_NN = E_NN + Zv(isp)*Zv(jsp)/r
      ENDIF 
   ENDDO 
  ENDDO
  E_NN = 0.5d0*E_NN
  
END SUBROUTINE 

