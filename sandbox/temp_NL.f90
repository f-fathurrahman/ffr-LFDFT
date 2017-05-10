
  !NprojTotMax = 0
  !DO isp=1,Nspecies
  !  NprojTot = sum( Ps_HGH_Params(isp)%Nproj_l )
  !  WRITE(*,*) Ps_HGH_Params(isp)%Nproj_l
  !  WRITE(*,*) 'isp, NprojTot = ', isp, NprojTot
  !  IF( NprojTotMax < NprojTot ) NprojTotMax = NprojTot
  !ENDDO 
!  iprj = 0
!  DO isp=1,Nspecies
!    DO l = 0,Ps_HGH_Params(isp)%lmax
!      DO ii = 1,Ps_HGH_Params(isp)%Nproj_l(l)
!        iprj = iprj + 1
!      ENDDO 
!    ENDDO 
!  ENDDO
!  NprojTotMax = iprj
!  WRITE(*,*) 'NprojTotMax = ', NprojTotMax

!  ALLOCATE( w_NL(Nspecies,NprojTotMax) )
!  w_NL(:,:) = 0.d0
!  iprj = 0
!  DO isp=1,Nspecies
!    DO l = 0,Ps_HGH_Params(isp)%lmax
!      DO ii = 1,Ps_HGH_Params(isp)%Nproj_l(l)
!        iprj = iprj + 1
!        w_NL(isp,iprj) = Ps_HGH_Params(isp)%h(l,ii,ii)
!      ENDDO 
!    ENDDO 
!  ENDDO 
