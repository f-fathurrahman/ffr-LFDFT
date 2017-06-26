
! Copyright (C) 2011 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

SUBROUTINE mixbroyden(iscl,n,msd,alpha,w0,nu,mu,f,df,u,a,d)
  IMPLICIT NONE 
  ! arguments
  INTEGER, INTENT(in) :: iscl,n,msd
  REAL(8), INTENT(in) :: alpha,w0
  REAL(8), INTENT(inout) :: nu(n),mu(n,2)
  REAL(8), INTENT(inout) :: f(n,2),df(n,msd)
  REAL(8), INTENT(inout) :: u(n,msd)
  REAL(8), INTENT(inout) :: a(msd,msd)
  REAL(8), INTENT(out) :: d
  ! local variables
  INTEGER :: jc,kp,kc
  INTEGER :: k,l,m,info
  REAL(8) :: t1
  ! automatic arrays
  INTEGER :: ipiv(msd)
  REAL(8) c(msd),beta(msd,msd),gamma(msd)
  REAL(8) work(msd)
  ! external functions
  REAL(8) :: dnrm2
  EXTERNAL dnrm2

  WRITE(*,*) 'Broyden mixing (ELK)'

  IF( n < 1 ) THEN 
    WRITE(*,*)
    WRITE(*,'("Error(mixbroyden): n < 1 : ",I8)') n
    WRITE(*,*)
    STOP 
  ENDIF 

  IF (msd < 2) THEN 
    WRITE(*,*)
    WRITE(*,'("Error(mixbroyden): msd < 2 : ",I8)') msd
    WRITE(*,*)
    STOP
  ENDIF 
  
  ! initialise mixer
  IF( iscl <= 0 ) THEN 
    CALL dcopy(n,nu,1,mu(:,1),1)
    CALL dcopy(n,nu,1,mu(:,2),1)
    f(:,1) = 0.d0
    df(:,1) = 0.d0
    u(:,1) = 0.d0
    a(:,:) = 0.d0
    d = 1.d0
    RETURN 
  ENDIF 

  ! current subspace dimension
  m = min(iscl+1,msd)
  ! current index modulo m
  jc = mod(iscl,m) + 1
  ! previous index modulo 2
  kp = mod(iscl-1,2) + 1
  ! current index modulo 2
  kc = mod(iscl,2) + 1
  f(:,kc) = nu(:) - mu(:,kp)
  d = sum(f(:,kc)**2)
  d = sqrt(d/dble(n))
  df(:,jc)=f(:,kc)-f(:,kp)
  t1 = dnrm2(n,df(:,jc),1)
  !
  IF( t1 > 1.d-8 ) t1=1.d0/t1
  
  CALL dscal(n,t1,df(:,jc),1)
  
  u(:,jc) = alpha*df(:,jc) + t1*(mu(:,kp)-mu(:,kc))
  
  DO k=1,m
    c(k) = dot_product(df(:,k),f(:,kc))
  ENDDO 

  DO k=1,m
    a(k,jc) = dot_product(df(:,jc),df(:,k))
    a(jc,k) = a(k,jc)
  ENDDO 

  beta(:,:) = a(:,:)
  DO k=1,m
    beta(k,k) = beta(k,k) + w0**2
  ENDDO 

  ! invert beta
  CALL dgetrf(m,m,beta,msd,ipiv,info)
  IF( info == 0 ) CALL dgetri(m,beta,msd,ipiv,work,m,info)
  
  IF (info /= 0) THEN 
    WRITE(*,*)
    WRITE(*,'("Error(mixbroyden): could not invert matrix")')
    WRITE(*,*)
    STOP 
  ENDIF 
  
  DO l=1,m
    gamma(l)=0.d0
    DO k=1,m
      gamma(l)=gamma(l)+c(k)*beta(k,l)
    ENDDO 
  ENDDO 

  nu(:) = mu(:,kp) + alpha*f(:,kc)
  DO l=1,m
    CALL daxpy(n,-gamma(l),u(:,l),1,nu,1)
  ENDDO 
  call dcopy(n,nu,1,mu(:,kc),1)
  RETURN 
END SUBROUTINE 

