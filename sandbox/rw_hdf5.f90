! efefer, 4 January 2016

SUBROUTINE write_matrix_hdf5( mat, Nrow, Ncol, filename, dsetname )
  USE hdf5
  IMPLICIT NONE
  !
  INTEGER :: Nrow, Ncol
  REAL(8) :: mat(Nrow,Ncol)
  CHARACTER(*) :: filename, dsetname
  !
  INTEGER(HID_T) :: file_id, dset_id, dspace_id
  INTEGER :: error
  INTEGER(HSIZE_T) :: data_dims(2)
  INTEGER :: rank

  CALL h5open_f( error )

  CALL h5fcreate_f( filename, H5F_ACC_TRUNC_F, file_id, error )

  !
  rank = 2
  data_dims = (/ Nrow, Ncol /)
  CALL h5screate_simple_f( rank, data_dims, dspace_id, error )
  
  CALL h5dcreate_f( file_id, dsetname, H5T_NATIVE_DOUBLE, dspace_id, dset_id, error )
  
  ! FIXME Need this?
  !CALL h5dopen_f( file_id, dsetname, dset_id, error )

  ! write data
  CALL h5dwrite_f( dset_id, H5T_NATIVE_DOUBLE, mat, data_dims, error )

  CALL h5dclose_f( dset_id, error )  ! close the dataset
  CALL h5sclose_f( dspace_id, error ) ! close the space
  CALL h5fclose_f( file_id, error )  ! close the file
  CALL h5close_f( error )  ! close the INTERFACE

END SUBROUTINE


!PROGRAM test
!  IMPLICIT NONE
!  !
!  INTEGER :: Nrow, Ncol
!  REAL(8), ALLOCATABLE :: matrix(:,:)
!  !
!  INTEGER :: ii, jj
!
!  Nrow = 4
!  Ncol = 4
!
!  ALLOCATE( matrix(Nrow,Ncol) )
!
!  DO jj=1,Ncol
!    DO ii=jj,Nrow
!      matrix(ii,jj) = 1.5d0*(ii+jj)/Nrow/Ncol
!      matrix(jj,ii) = matrix(ii,jj)
!    ENDDO
!  ENDDO
!
!  DO ii=1,Nrow
!    WRITE(*,*)
!    DO jj=1,Ncol
!      WRITE(*,'(1x,F18.10)', advance='NO') matrix(ii,jj)
!    ENDDO
!  ENDDO
!  WRITE(*,*)
!
!  CALL write_matrix_hdf5( matrix, Nrow, Ncol, 'matrix.h5', 'test-matrix' )
!
!  DEALLOCATE( matrix )
!END PROGRAM
!

