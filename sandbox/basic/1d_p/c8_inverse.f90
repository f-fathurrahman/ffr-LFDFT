
!------------------------------
subroutine c8_inverse(dim1,A)
!------------------------------
  implicit none
  ! Local variables
  INTEGER :: dim1,dim2
  COMPLEX(8), INTENT(INOUT) :: A(dim1,dim1)
  ! Local variables
  INTEGER :: lwork
  COMPLEX(8), ALLOCATABLE :: work(:)
  INTEGER, ALLOCATABLE :: ipiv(:)
  INTEGER :: lda
  INTEGER :: info
 
  ! Factorize A
  lda = dim1
  ALLOCATE( ipiv(dim1) )
  CALL zgetrf(dim1,dim2,A,lda,ipiv,info)
  IF(info /= 0) THEN
    WRITE(*,*) 'Error in c8_inverse: zgetrf returned info = ', info
    STOP
  ENDIF

  ! Calculate inverse of A
  lwork = 2*dim1
  ALLOCATE(work(lwork))
  CALL zgetri(dim1,A,lda,ipiv,work,lwork,info)
  IF(info /= 0) THEN
    WRITE(*,*) 'Error in c8_inverse: zgetri returned info = ', info
    STOP
  ENDIF

  ! Free memory
  DEALLOCATE(ipiv)
  DEALLOCATE(work)
END SUBROUTINE


