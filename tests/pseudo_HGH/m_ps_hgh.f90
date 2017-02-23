MODULE m_ps_hgh

  USE m_atom
  IMPLICIT NONE

  !> The following data type contains:
  !!   (a) the pseudopotential parameters, as read from a *.hgh file,
  !!   (b) auxiliary intermediate functions, to store stuff before passing it to the "ps" variable.
  TYPE ps_hgh_t
    !< HGH parameters.
    CHARACTER(len=5) :: atom_name
    INTEGER          :: z_val
    REAL(8)          :: rlocal
    REAL(8)          :: rc(0:3)
    REAL(8)          :: c(1:4)
    REAL(8)          :: h(0:3, 1:3, 1:3)
    REAL(8)          :: k(0:3, 1:3, 1:3)

    TYPE(valconf_t)  :: conf
    INTEGER          :: l_max     !< Maximum l for the Kleinman-Bylander component.

    REAL(8), ALLOCATABLE :: vlocal(:) !< Local potential
    REAL(8), ALLOCATABLE :: kb(:,:,:) !< KB projectors
    REAL(8), ALLOCATABLE :: kbr(:)    !< KB radii
    REAL(8), ALLOCATABLE :: rphi(:,:), eigen(:)

    !> Logarithmic grid parameters
    !type(logrid_t) :: g
  END TYPE ps_hgh_t


CONTAINS


  FUNCTION hgh_eval_proj( psp, l, i, r ) RESULT( pil )
    TYPE(ps_hgh_t) :: psp
    INTEGER :: l, i  ! l can be zero, i varies from 1 to 3
    REAL(8) :: r
    REAL(8) :: pil
    !
    REAL(8) :: num, denum
    !
    num = sqrt(2.d0)*r**(l + 2*(i-1)) * exp( -(r/psp%rc(l))**2 )
    denum = (psp%rc(l))**( l + (4*i-1)/2 ) * gamma( dble(l) + (4*i-1)/2.d0 )
    pil = num/denum
  END FUNCTION


  FUNCTION hgh_eval_Vloc( psp, r ) RESULT(Vloc)
    TYPE(ps_hgh_t) :: psp
    REAL(8) :: r
    !
    INTEGER :: i
    REAL(8) :: rrloc ! r/rlocal
    REAL(8) :: term1, Vloc

    rrloc = r/psp%rlocal
    term1 = psp%c(1)
    DO i = 2, 4
      term1 = term1 + psp%c(i)*rrloc**(2*(i-1))
    ENDDO
    Vloc = -psp%z_val/r * erf( rrloc/sqrt(2.d0) ) + exp(-0.5d0*rrloc**2)*term1
  END FUNCTION


  ! ---------------------------------------------------------
  SUBROUTINE hgh_init(psp, filename)
    TYPE(ps_hgh_t), INTENT(INOUT) :: psp
    CHARACTER(len=*), INTENT(IN)  :: filename

    INTEGER :: iunit, i
    LOGICAL :: found

    INQUIRE(file=filename, exist=found)
    IF(.NOT.found) THEN
      WRITE(*,*) 'ERROR: Pseudopotential file: ', trim(filename), ' is not found.'
      STOP
    ENDIF

    iunit = 303  ! just saying
    OPEN(unit=iunit,file=filename, action='read', form='formatted', status='old')
    i = load_params(iunit, psp)
    IF(i /= 0) THEN
      WRITE(*,*) 'ERROR: Reading psfile: ', trim(filename), '.'
      STOP
    ENDIF
    CLOSE(iunit)

    ! Finds out psp%l_max. The most special cases are H, He, Li_sc and Be_sc, where psp%l_max = -1.
    psp%l_max = 0
    DO WHILE(psp%rc(psp%l_max) > 0.01d0)
      psp%l_max = psp%l_max + 1
      IF(psp%l_max > 3) EXIT
    ENDDO
    psp%l_max = psp%l_max - 1

  END SUBROUTINE


  ! ---------------------------------------------------------
  FUNCTION load_params(iunit, params)
    INTEGER, INTENT(IN)  :: iunit ! where to read from
    TYPE(ps_hgh_t), INTENT(out) :: params ! should INOUT instead?
    INTEGER :: load_params ! 0 if success, 1 otherwise.

    INTEGER :: i, ios, j, k
    CHARACTER(len=VALCONF_STRING_LENGTH) :: line

    ! Set initially everything to zero.
    params%c(1:4) = 0.d0
    params%rlocal = 0.d0
    params%rc = 0.d0
    params%h = 0.d0
    params%k = 0.d0

    ! get valence configuration
    READ(iunit,'(a)') line
    CALL read_valconf(line, params%conf)

    ! Reads the file in a hopefully smart way
    ios = 1
    j = 5
    READ(iunit,'(a)') line
    DO WHILE((ios /= 0) .AND. (j > 0))
      j = j - 1
      READ(line, *, iostat=ios) params%atom_name, params%z_val, params%rlocal, params%c(1:j)
    ENDDO
    IF (j<1) READ(line, *, iostat=ios) params%atom_name, params%z_val, params%rlocal
    IF ( ios /= 0 ) then
      load_params = 1
      RETURN
    ENDIF

    READ(iunit,'(a)', iostat = ios) line
    IF(ios /= 0) THEN
      load_params = 0
      RETURN
    ENDIF
    ios = 1
    j = 4
    DO WHILE((ios /= 0) .AND. (j > 0))
      j = j - 1
      READ(line, *, iostat=ios) params%rc(0), (params%h(0, i, i), i = 1, j)
    ENDDO
    IF(j < 0) THEN
      load_params = 2
      RETURN
    ENDIF

    kloop: DO k = 1, 3
      READ(iunit, '(a)', iostat=ios) line
      IF(ios /= 0) EXIT kloop
      ios = 1
      j = 4
      DO WHILE((ios /= 0) .AND. (j > 0))
        j = j - 1
        READ(line, *, iostat=ios) params%rc(k), (params%h(k, i, i), i = 1, j)
      ENDDO
      IF(params%rc(k) == 0.d0) EXIT kloop
      READ(iunit, '(a)') line
      ios = 1
      j = 4
      DO WHILE((ios /= 0) .AND. (j>0))
        j = j - 1
        READ(line, *, iostat=ios) (params%k(k, i, i), i = 1, 3)
      ENDDO
    ENDDO kloop

    ! Fill in the rest of the parameter matrices...
    ! Fill in the rest of the parameter matrices...
    params%h(0, 1, 2) = -0.5d0      * sqrt(3.d0/5.d0)       * params%h(0, 2, 2)
    params%h(0, 1, 3) =  0.5d0      * sqrt(5.d0/21.d0)      * params%h(0, 3, 3)
    params%h(0, 2, 3) = -0.5d0      * sqrt(100.d0/63.d0)    * params%h(0, 3, 3)
    params%h(1, 1, 2) = -0.5d0      * sqrt(5.d0/7.d0)       * params%h(1, 2, 2)
    params%h(1, 1, 3) =  1.d0/6.d0  * sqrt(35.d0/11.d0)     * params%h(1, 3, 3)
    params%h(1, 2, 3) = -1.d0/6.d0  * (14.d0 / sqrt(11.d0)) * params%h(1, 3, 3)
    params%h(2, 1, 2) = -0.5d0      * sqrt(7.d0/9.d0)       * params%h(2, 2, 2)
    params%h(2, 1, 3) =  0.5d0      * sqrt(63.d0/143.d0)    * params%h(2, 3, 3)
    params%h(2, 2, 3) = -0.5d0      * (18.d0/sqrt(143.0d0)) * params%h(2, 3, 3)

    params%k(0, 1, 2) = -0.5d0      * sqrt(3.d0/5.d0)       * params%k(0, 2, 2)
    params%k(0, 1, 3) =  0.5d0      * sqrt(5.d0/21.d0)      * params%k(0, 3, 3)
    params%k(0, 2, 3) = -0.5d0      * sqrt(100.d0/63.d0)    * params%k(0, 3, 3)
    params%k(1, 1, 2) = -0.5d0      * sqrt(5.d0/7.d0)       * params%k(1, 2, 2)
    params%k(1, 1, 3) =  1.d0/6.d0  * sqrt(35.d0/11.d0)     * params%k(1, 3, 3)
    params%k(1, 2, 3) = -1.d0/6.d0  * (14.d0/sqrt(11.d0))   * params%k(1, 3, 3)
    params%k(2, 1, 2) = -0.5d0      * sqrt(7.d0/9.d0)       * params%k(2, 2, 2)
    params%k(2, 1, 3) =  0.5d0      * sqrt(63.0/143.0)      * params%k(2, 3, 3)
    params%k(2, 2, 3) = -0.5d0      * (18.d0/sqrt(143.d0))  * params%k(2, 3, 3)

    ! Parameters are symmetric.
    DO k = 0, 3
      DO i = 1, 3
        DO j = i + 1, 3
          params%h(k, j, i) = params%h(k, i, j)
          params%k(k, j, i) = params%k(k, i, j)
        ENDDO
      ENDDO
    ENDDO

    load_params = 0
  END FUNCTION load_params


  SUBROUTINE hgh_info( ps )
    IMPLICIT NONE
    !
    TYPE(ps_hgh_t) :: ps
    INTEGER :: i, j, l

    WRITE(*,*) 'atom_name = ', trim(ps%atom_name)
    WRITE(*,'(1x,A,I5)') 'z_val = ', ps%z_val
    WRITE(*,'(1x,A,I5)') 'l_max = ', ps%l_max
    WRITE(*,'(1x,A,F18.10)') 'rlocal = ', ps%rlocal
    !
    IF( ps%l_max > 0 ) THEN
      WRITE(*,*) 'rc (for non-local potential) = '
      DO i = 0, ps%l_max
        WRITE(*,'(6x,I5,F18.10)') i, ps%rc(i)
      ENDDO
    ENDIF
    !
    WRITE(*,*) 'c (for local potential ) = '
    DO i = 1, 4
      WRITE(*,'(6x,I5,F18.10)') i, ps%c(i)
    ENDDO
    !
    DO l = 0, ps%l_max
      WRITE(*,*)
      WRITE(*,*) 'Matrix h for l = ', l
      DO i = 1, 3
        WRITE(*,*)
        DO j = 1, 3
          WRITE(*,'(F18.10)',advance='no') ps%h(l,i,j)
        ENDDO
      ENDDO
      WRITE(*,*)
    ENDDO
  END SUBROUTINE

END MODULE

