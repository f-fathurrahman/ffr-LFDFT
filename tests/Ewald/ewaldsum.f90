PROGRAM ewaldsum

! Ewaldsum fortran 90 code
!   4/23/01  NAWH
!   Input :   T1x  T1y   T1z    (Cartesian components of lattice vectors)
!         :   T2x  T2y   T2z
!         :   T3x  T3y   T3z
!         :   nions             (Number of ions in unit cell)
!         :   q(i), tau(1:3,i), i=1,nion  (charge and fractional position
!                                             of ions)
!             Note:  ion positions = tau(1,i)*T1 + tau(2,i)*T2 + tau(3,i)*T3
!         :   eps    (error tolerance)

  IMPLICIT NONE

  REAL(8) :: pi,eps, t1(3), t2(3), t3(3), g1(3),g2(3),g3(3),volcry, arg,x,gexp
  REAL(8) :: eta, totalcharge, g1m, g2m, g3m, t1m, t2m, t3m,gcut, tmax ,ebsl,seta
  REAL(8) :: tpi,glast2, con, con2 , cccc , ewald , v(3), w(3), rmag2 , prod
  INTEGER  :: nion, i, j, k, ng, nt , mmm1, mmm2, mmm3 , a , b
  REAL(8), ALLOCATABLE :: q(:), tau(:,:)

  Write(6,*) ' Enter T1x, T1y, T1z in bohr units'
  Read(5,*) t1(1),t1(2),t1(3)
  Write(6,*) ' Enter T2x, T2y, T2z in bohr units'
  Read(5,*) t2(1),t2(2),t2(3)
  Write(6,*) ' Enter T3x, T3y, T3z in bohr units'
  Read(5,*) t3(1),t3(2),t3(3)
  
  pi = acos(-1.d0)
  volcry   = t1(1)*(t2(2)*t3(3)-t2(3)*t3(2)) +  &
             t1(2)*(t2(3)*t3(1)-t2(1)*t3(3)) +  &
             t1(3)*(t2(1)*t3(2)-t2(2)*t3(1))

  g1(1) = 2 * pi * (t2(2)*t3(3)-t2(3)*t3(2))/volcry
  g1(2) = 2 * pi * (t2(3)*t3(1)-t2(1)*t3(3))/volcry
  g1(3) = 2 * pi * (t2(1)*t3(2)-t2(2)*t3(1))/volcry
  g2(1) = 2 * pi * (t3(2)*t1(3)-t3(3)*t1(2))/volcry
  g2(2) = 2 * pi * (t3(3)*t1(1)-t3(1)*t1(3))/volcry
  g2(3) = 2 * pi * (t3(1)*t1(2)-t3(2)*t1(1))/volcry
  g3(1) = 2 * pi * (t1(2)*t2(3)-t1(3)*t2(2))/volcry
  g3(2) = 2 * pi * (t1(3)*t2(1)-t1(1)*t2(3))/volcry
  g3(3) = 2 * pi * (t1(1)*t2(2)-t1(2)*t2(1))/volcry
  
  volcry = abs(volcry)

  t1m = SQRT(DOT_PRODUCT(t1,t1))
  t2m = SQRT(DOT_PRODUCT(t2,t2))
  t3m = SQRT(DOT_PRODUCT(t3,t3))
  g1m = SQRT(DOT_PRODUCT(g1,g1))
  g2m = SQRT(DOT_PRODUCT(g2,g2))
  g3m = SQRT(DOT_PRODUCT(g3,g3))

  Write(6,*) 'Input total number of ions in unit cell'
  Read(5,*) nion

  Allocate(q(nion),tau(3,nion))

  Write(6,*) 'For each ion, input q and tau (in fractional coordinates of T)'

  do i=1,nion
    read(5,*) q(i),tau(1,i),tau(2,i),tau(3,i)
  enddo

  Write(6,*) ' Input Gcut for reciprocal lattice sum and error tolerance'
  Read(5,*) gcut , ebsl

   tpi = 2.d0*pi
   con = volcry/(4.d0*pi)
   con2 = (4.d0*pi)/volcry
   glast2 = gcut*gcut
   gexp = -log(ebsl)
   eta = glast2/gexp

   Write(6,*) ' eta value for this calculation' , eta
   cccc = sqrt(eta/pi)

   x = 0.d0
  DO i = 1,nion 
    x = x + q(i)**2
  ENDDO

  totalcharge = sum(q) 
  WRITE(6,*) ' Total charge = ', totalcharge
  ewald = -cccc*x - 4.d0*pi*(totalcharge**2)/(volcry*eta)

   tmax = sqrt(2.d0*gexp/eta)
   seta = sqrt(eta)/2.d0

      mmm1=tmax/t1m+1.5d0
      mmm2=tmax/t2m+1.5d0  
      mmm3=tmax/t3m+1.5d0  
      

  Write (6,*) ' lattice summation indices -- ', mmm1,mmm2,mmm3
   do a = 1,nion
      do b = 1,nion
         v(:) = (tau(1,a)-tau(1,b))*t1(:) + (tau(2,a)-tau(2,b))*t2(:) &
              + (tau(3,a)-tau(3,b))*t3(:)
         prod=q(a)*q(b)
         do i = -mmm1, mmm1
            do j = -mmm2, mmm2
               do k = -mmm3, mmm3
                    if ((a.ne.b).or.((abs(i)+abs(j)+abs(k)).ne.0)) then
                     w(:) = v(:) + i*t1 + j*t2 + k*t3
                     rmag2 = sqrt(DOT_PRODUCT(w,w))
                     arg=rmag2*seta 
                       ewald = ewald + prod*erfc(arg)/rmag2
                    endif
               enddo
             enddo
          enddo
       enddo
    enddo

      mmm1=gcut/g1m+1.5
      mmm2=gcut/g2m+1.5
      mmm3=gcut/g3m+1.5
      
  Write(6,*) ' Reciprocal lattice summation indices --', mmm1,mmm2,mmm3
   do i = -mmm1, mmm1
      do j = -mmm2, mmm2
         do k = -mmm3, mmm3
              if ((abs(i)+abs(j)+abs(k)).ne.0) then
               w(:) = i*g1(:) + j*g2(:) + k*g3(:)
               rmag2=DOT_PRODUCT(w,w)
               x=con2*exp(-rmag2/eta)/rmag2
               do a = 1,nion
                  do b = 1,nion
                     v(:) = tau(:,a)-tau(:,b)
                     prod = q(a)*q(b)
                     arg=tpi*(i*v(1)+j*v(2)+k*v(3))
                     ewald=ewald + x*prod*cos(arg)
                  enddo
               enddo
              endif
          enddo
       enddo
   ENDDO

 Write(6,'(1x,A,F18.10)') 'Ewald energy in Ry', ewald
 WRITE(*,'(1x,A,F18.10)') 'Ewald energy in Ha', 2.d0*ewald

END PROGRAM 
