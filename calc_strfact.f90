! Calculate structure factor
SUBROUTINE calc_strfact( Na, Xpos, Nspecies, atm2species, Ng, Gv, strf )
  
  IMPLICIT NONE
  INTEGER :: Na
  REAL(8) :: Xpos(3,Na)
  INTEGER :: Nspecies
  INTEGER :: atm2species(Na) ! mapping from idx atm to idx species
  INTEGER :: Ng
  REAL(8) :: Gv(3,Ng)
  COMPLEX(8) :: strf(Ng,Nspecies)
  !
  INTEGER :: ia, isp, ig
  REAL(8) :: GX

  strf(:,:) = cmplx(0.d0,0.d0,kind=8)
  DO ia = 1,Na
    isp = atm2species(ia)
    DO ig = 1,Ng
      GX = Xpos(1,ia)*Gv(1,ig) + Xpos(2,ia)*Gv(2,ig) + Xpos(3,ia)*Gv(3,ig)
      strf(ig,isp) = strf(ig,isp) + cmplx( cos(GX), -sin(GX), kind=8 )
    ENDDO 
  ENDDO 

END SUBROUTINE 

