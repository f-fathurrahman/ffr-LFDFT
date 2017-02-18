SUBROUTINE write_grad( Nrow, Ncol, grad, iter )
  IMPLICIT NONE
  !
  INTEGER :: Nrow, Ncol, iter
  REAL(8) :: grad(Nrow,Ncol)
  !
  INTEGER :: ios, iu, ic, ir
  CHARACTER(32) :: filnam, str1

  IF( iter < 10 ) THEN
    WRITE(str1,'(I1)') iter
  ELSEIF( iter < 100 ) THEN
    WRITE(str1,'(I2)') iter
  ENDIF
  filnam = 'grad_'//str1
  !WRITE(*,*) 'filnam = ', filnam
  iu = 100 + iter
  OPEN( unit=iu, file=filnam, action='write', iostat=ios )
  DO ir = 1, Nrow
    DO ic = 1, Ncol
      WRITE(iu,'(F18.10)', advance='no') grad(ir,ic)
    ENDDO
    WRITE(iu,*)
  ENDDO
  CLOSE(iu)
END SUBROUTINE

