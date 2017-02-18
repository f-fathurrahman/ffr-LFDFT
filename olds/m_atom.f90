! ffr, modified from Octopus


MODULE m_atom

  IMPLICIT NONE

  INTEGER, PARAMETER :: VALCONF_STRING_LENGTH=80

  TYPE valconf_t
    INTEGER           :: z
    CHARACTER(len=3)  :: symbol
    INTEGER           :: typ    !< 0 for the most normal valence configuration, 1 for semicore.
    INTEGER           :: p        !< number of orbitals.
    INTEGER           :: n(12)     !< n quantum number
    INTEGER           :: l(12)     !< l quantum number
    REAL(8)           :: occ(12,2) !< occupations of each level
  END TYPE valconf_t


CONTAINS


  ! ---------------------------------------------------------
  !> Subroutines to write and read valence configurations.
  SUBROUTINE valconf_null(c)
    TYPE(valconf_t), INTENT(OUT) :: c

    c%z = 0
    c%symbol = ''
    c%typ = 0
    c%p = 0
    c%n = 0
    c%l = 0
    c%occ = 0.d0

  END SUBROUTINE valconf_null


  ! ---------------------------------------------------------
  SUBROUTINE read_valconf(s, c)
    CHARACTER(len=VALCONF_STRING_LENGTH), intent(in) :: s
    TYPE(valconf_t), intent(out) :: c

    INTEGER :: j
    CHARACTER(len=1) :: lvalues(1:6)

    READ(s,'(i2,1x,a2,i1,1x,i1,1x,6(i1,a1,f6.3,1x))') c%z, c%symbol, c%typ, c%p,&
         (c%n(j),lvalues(j),c%occ(j,1),j=1,c%p)
    DO j = 1, c%p
       SELECT CASE(lvalues(j))
       CASE('s'); c%l(j) = 0
       CASE('p'); c%l(j) = 1
       CASE('d'); c%l(j) = 2
       CASE('f'); c%l(j) = 3
       CASE default
         WRITE(*,*) 'ERROR: in read_valconf'
         STOP
       END SELECT
    ENDDO

  END SUBROUTINE

END MODULE

