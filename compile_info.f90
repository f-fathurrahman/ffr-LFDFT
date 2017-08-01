!
! FIXME: Search info about predefined preprocessor macros

SUBROUTINE compile_info()
!PROGRAM compile_info
  IMPLICIT NONE 

#ifdef __PGI
  WRITE(*,*)
  WRITE(*,*) 'Compiled using PGI Fortran compiler'
#endif

#ifdef __INTEL_COMPILER
  WRITE(*,*)
  WRITE(*,*) 'Compiled using Intel Fortran compiler'
#endif

#ifdef __GFORTRAN__
  WRITE(*,*)
  WRITE(*,*) 'Compiled using GNU Fortran compiler'
#endif

#ifdef __G95__
  WRITE(*,*)
  WRITE(*,*) 'Compiled using G95 Fortran compiler'
#endif

#ifdef __SUNPRO_F95
  WRITE(*,*)
  WRITE(*,*) 'Compiled using Sun F95 Compiler'
#endif

END SUBROUTINE 
!END PROGRAM 

