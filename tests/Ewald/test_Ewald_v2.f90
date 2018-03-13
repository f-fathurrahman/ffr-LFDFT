MODULE ewald_sum

IMPLICIT NONE 

REAL(8) :: L_vec1(3), L_vec2(3), L_vec3(3)
REAL(8) :: g1(3), g2(3), g3(3)
REAL(8) :: volume, eta
INTEGER :: lr(3), lg(3)
  
CONTAINS 

SUBROUTINE ewald_setup( gcut, error )

  IMPLICIT NONE
  REAL(8) :: pi,eps, arg,x,gexp
  REAL(8) :: gcut, error

  pi = acos(-1.d0)
  
  volume   = L_vec1(1)*(L_vec2(2)*L_vec3(3)-L_vec2(3)*L_vec3(2)) +  &
             L_vec1(2)*(L_vec2(3)*L_vec3(1)-L_vec2(1)*L_vec3(3)) +  &
             L_vec1(3)*(L_vec2(1)*L_vec3(2)-L_vec2(2)*L_vec3(1))

  g1(1) = 2 * pi * (L_vec2(2)*L_vec3(3)-L_vec2(3)*L_vec3(2))/volume
  g1(2) = 2 * pi * (L_vec2(3)*L_vec3(1)-L_vec2(1)*L_vec3(3))/volume
  g1(3) = 2 * pi * (L_vec2(1)*L_vec3(2)-L_vec2(2)*L_vec3(1))/volume
  g2(1) = 2 * pi * (L_vec3(2)*L_vec1(3)-L_vec3(3)*L_vec1(2))/volume
  g2(2) = 2 * pi * (L_vec3(3)*L_vec1(1)-L_vec3(1)*L_vec1(3))/volume
  g2(3) = 2 * pi * (L_vec3(1)*L_vec1(2)-L_vec3(2)*L_vec1(1))/volume
  g3(1) = 2 * pi * (L_vec1(2)*L_vec2(3)-L_vec1(3)*L_vec2(2))/volume
  g3(2) = 2 * pi * (L_vec1(3)*L_vec2(1)-L_vec1(1)*L_vec2(3))/volume
  g3(3) = 2 * pi * (L_vec1(1)*L_vec2(2)-L_vec1(2)*L_vec2(1))/volume
  
  volume = abs(volume)

  eta = -gcut**2/log(error)

  lr(1) = sqrt(-2*log(error)/eta)/SQRT(DOT_PRODUCT(L_vec1,L_vec1)) + 1.5d0
  lr(2) = sqrt(-2*log(error)/eta)/SQRT(DOT_PRODUCT(L_vec2,L_vec2)) + 1.5d0
  lr(3) = sqrt(-2*log(error)/eta)/SQRT(DOT_PRODUCT(L_vec3,L_vec3)) + 1.5d0
      
  lg(1) = gcut/SQRT(DOT_PRODUCT(g1,g1)) + 1.5d0
  lg(2) = gcut/SQRT(DOT_PRODUCT(g2,g2)) + 1.5d0
  lg(3) = gcut/SQRT(DOT_PRODUCT(g3,g3)) + 1.5d0

END SUBROUTINE 


SUBROUTINE ewaldsum( nion, atpos, q ,ewald )

  IMPLICIT NONE
  
  REAL(8) :: pi,eps, arg,x,gexp
  REAL(8) :: gcut ,error
  REAL(8) :: ewald , v(3), w(3), r2 , prod
  INTEGER :: nion, i, j, k, ng, nt , l1, l2, l3 , a , b
  REAL(8) :: q(nion), atpos(3,nion)

  pi = acos(-1.d0)

  ewald = -sqrt(eta/pi)*sum(q**2)-4*pi*(sum(q)**2)/(volume*eta)

  DO a = 1,nion
    DO b = 1,nion
      v(:) = atpos(:,a) - atpos(:,b)
      prod = q(a)*q(b)
      DO i = -lr(1), lr(1)
      DO j = -lr(2), lr(2)
      DO k = -lr(3), lr(3)
        IF ( (a /= b) .or. ( ( abs(i) + abs(j) + abs(k) ) /= 0) ) THEN 
          w(:) = v(:) + i*L_vec1 + j*L_vec2 + k*L_vec3
          r2 = sqrt(DOT_PRODUCT(w,w))
          arg = r2*sqrt(eta)/2.d0 
          ewald = ewald + prod*erfc(arg)/r2
        ENDIF 
      ENDDO 
      ENDDO 
      ENDDO 
    ENDDO 
  ENDDO

  DO i = -lg(1), lg(1)
  DO j = -lg(2), lg(2)
  DO k = -lg(3), lg(3)
    IF( ( abs(i) + abs(j) + abs(k) ) /= 0 ) THEN 
      w(:) = i*g1(:) + j*g2(:) + k*g3(:)
      r2 = dot_product(w,w)
      x = 4.d0*pi/volume*exp(-r2/eta)/r2
      DO a = 1,nion
        DO b = 1,nion
          prod = q(a)*q(b)
          v(:) = atpos(:,a) - atpos(:,b)
          arg = v(1)*w(1) + v(2)*w(2) + v(3)*w(3)                     
          ewald = ewald + x*prod*cos(arg)
        ENDDO 
      ENDDO 
    ENDIF 
  ENDDO 
  ENDDO 
  ENDDO 

END SUBROUTINE 

END MODULE EWALD_SUM

PROGRAM test_ewald
  USE ewald_sum
  IMPLICIT NONE 
  REAL(8) :: gcut ,error,h,su,ewald,p(3),pp,suo,sut,dx
  INTEGER :: nion, i, i1,i2,i3,n,x
  REAL(8), ALLOCATABLE :: q(:), atpos(:,:)
  CHARACTER(128) :: filename
  REAL(8) :: alat

  CALL getarg(1,filename)
  WRITE(*,*) 'Reading from input file: ', trim(filename)

  OPEN(1,file=filename)
  
  READ(1,*) alat
  READ(1,*) L_vec1(1), L_vec1(2), L_vec1(3)
  READ(1,*) L_vec2(1), L_vec2(2), L_vec2(3)
  READ(1,*) L_vec3(1), L_vec3(2), L_vec3(3)
  READ(1,*) nion

  ALLOCATE( q(nion), atpos(3,nion) )
  
  DO i=1,nion
    READ(1,*) q(i), atpos(1,i), atpos(2,i), atpos(3,i)
  ENDDO
  CLOSE(1)

  WRITE(*,*)
  WRITE(*,*) 'Read data:'
  WRITE(*,*)
  WRITE(*,*) 'Lattice vectors:'
  WRITE(*,*)
  WRITE(*,*) 'alat = ', alat
  WRITE(*,*) L_vec1(1), L_vec1(2), L_vec1(3)
  WRITE(*,*) L_vec2(1), L_vec2(2), L_vec2(3)
  WRITE(*,*) L_vec3(1), L_vec3(2), L_vec3(3)
  WRITE(*,*)
  WRITE(*,*) 'Ion data:'
  DO i=1,Nion
    WRITE(*,'(F5.2,3F18.10)') q(i), atpos(1,i), atpos(2,i), atpos(3,i)
  ENDDO

  gcut = 2.d0
  error = 0.000000001d0
  
  CALL ewald_setup( gcut, error )
  CALL ewaldsum( nion, atpos, q, ewald )
  
  WRITE(*,*)
  WRITE(*,'(1x,A,F18.10)') 'Calculated Ewald energy = ', ewald/alat

!  DO x=0,4999
!    dx = 0.001d0*x + 0.001d0
!    WRITE(101,*) dx, 1/dx, erf(dx)/dx, (1.d0-erf(dx))/dx
!  ENDDO

END PROGRAM 

