! adapted from Octopus distribution
! 2 Feb 2016

MODULE m_ps_upf

  USE m_atom
  IMPLICIT NONE


  TYPE ps_upf_t
    TYPE(valconf_t)    :: conf

    INTEGER :: kb_nc
    INTEGER :: l_local
    REAL(8) :: local_radius
    REAL(8), ALLOCATABLE :: kb_radius(:)

    !> The contents of the file:
    !!Header
    INTEGER      :: version  !< UPF file version number
    CHARACTER(3) :: symbol   !< Element label
    CHARACTER(2) :: pstype     !< Pseudo type (NC or US). In octopus only NC is implemented.
    REAL(8)      :: z_val    !< z valence
    INTEGER      :: l_max    !< maximum angular momentum component
    INTEGER      :: n_proj   !< number of projectors
    INTEGER      :: n_wfs    !< number of wavefunctions
    INTEGER, ALLOCATABLE :: n(:)
    INTEGER, ALLOCATABLE :: l(:)
    REAL(8), ALLOCATABLE :: occ(:)

    !>Radial mesh
    INTEGER :: np
    REAL(8), ALLOCATABLE :: r(:)
    REAL(8), ALLOCATABLE :: drdi(:)

    !>nlcc
    LOGICAL        :: nlcc
    REAL(8), ALLOCATABLE :: core_density(:)

    !>KB projectors
    REAL(8), ALLOCATABLE :: v_local(:)
    INTEGER, ALLOCATABLE :: proj_l(:)
    INTEGER, ALLOCATABLE :: proj_np(:)
    REAL(8), ALLOCATABLE :: proj(:,:)
    REAL(8), ALLOCATABLE :: e(:)

    !>Wavefunctions
    REAL(8), ALLOCATABLE :: wfs(:,:)

    !>Valence charge
    REAL(8), ALLOCATABLE :: rho(:)
    
    !>Extra information
    REAL(8), ALLOCATABLE :: proj_j(:)

  END TYPE ps_upf_t


CONTAINS

  SUBROUTINE ps_upf_init(ps_upf, filename)
    TYPE(ps_upf_t),   INTENT(INOUT) :: ps_upf
    CHARACTER(len=*), INTENT(IN)    :: filename

    CHARACTER(len=256) :: filename2
    INTEGER :: iunit, l
    LOGICAL :: found
    LOGICAL, ALLOCATABLE :: found_l(:)    

    ! Find out the file and read it.
    filename2 = trim(filename) // '.UPF'
    INQUIRE(file=filename2, exist=found)
    IF(.NOT. found) THEN
      WRITE(*,*) 'Error: pseudopotential file is not found'
      STOP
    ENDIF

    WRITE(*,*) 'Reading pseudopotential from file:', trim(filename2)
    
    iunit = 333
    OPEN(unit=iunit, file=filename2, action='read', form='formatted', status='old')
    
    CALL ps_upf_file_read(iunit, ps_upf)
    
    CLOSE(iunit)

    !Build valence configuration
!    call valconf_null(ps_upf%conf)
!    ps_upf%conf%symbol = ps_upf%symbol
!    ps_upf%conf%type = 0
!    ps_upf%conf%p = ps_upf%n_wfs
!    ps_upf%conf%n(1:ps_upf%n_wfs) = ps_upf%n(1:ps_upf%n_wfs)
!    ps_upf%conf%l(1:ps_upf%n_wfs) = ps_upf%l(1:ps_upf%n_wfs)
!    ps_upf%conf%occ(1:ps_upf%n_wfs,1) = ps_upf%occ(1:ps_upf%n_wfs)
!
!    ps_upf%l_max = maxval(ps_upf%l)
!
!    !Check if the local component is one of the angular momentum channels
!    SAFE_ALLOCATE(found_l(0:ps_upf%l_max))
!    found_l = .true.
!    do l = 0, ps_upf%l_max
!      if (any(ps_upf%proj_l == l)) then
!        found_l(l) = .false.
!      end if
!    end do
!    if (count(found_l) /= 1) then
!      ps_upf%l_local = -1
!    else
!      do l = 0, ps_upf%l_max
!        if (found_l(l)) then
!          ps_upf%l_local = l
!          exit
!        end if
!      end do
!    end if
!
!    ! Define the KB-projector cut-off radii
!    call ps_upf_cutoff_radii(ps_upf)
!
!    ! check norm of rphi
!    call ps_upf_check_rphi(ps_upf)

  END SUBROUTINE

  
  ! ---------------------------------------------------------
  SUBROUTINE ps_upf_end(ps_upf)
    TYPE(ps_upf_t), INTENT(INOUT) :: ps_upf

    DEALLOCATE(ps_upf%kb_radius)

    DEALLOCATE(ps_upf%n)
    DEALLOCATE(ps_upf%l)
    DEALLOCATE(ps_upf%occ)
    DEALLOCATE(ps_upf%core_density)
    DEALLOCATE(ps_upf%r)
    DEALLOCATE(ps_upf%drdi)
    DEALLOCATE(ps_upf%v_local)
    DEALLOCATE(ps_upf%proj_l)
    DEALLOCATE(ps_upf%proj_np)
    DEALLOCATE(ps_upf%proj)
    DEALLOCATE(ps_upf%e)
    DEALLOCATE(ps_upf%wfs)
    DEALLOCATE(ps_upf%rho)
    DEALLOCATE(ps_upf%proj_j)

  END SUBROUTINE ps_upf_end



  SUBROUTINE ps_upf_file_read(iunit, ps_upf)
    INTEGER,        INTENT(IN)    :: iunit
    TYPE(ps_upf_t), intent(inout) :: ps_upf

    INTEGER :: ip, np, i, ir, idummy, ii, jj, n_dij
    CHARACTER(len=2) :: nl
    CHARACTER(len=80) :: dummy
    LOGICAL :: ok
    REAL(8) :: fdummy

    ps_upf%kb_nc = 1

    !Header info
    CALL init_tag(iunit, 'PP_HEADER', .true.)

    READ(iunit,*) ps_upf%version, dummy  ! n        "Version Number"
    READ(iunit,*) ps_upf%symbol, dummy   ! psd      "Element"
    READ(iunit,*) ps_upf%pstype, dummy     ! US|NC|PAW    "Ultrasoft|Norm conserving|Projector-augmented"
    IF(ps_upf%pstype /= "NC") THEN
      WRITE(*,*) 'Error: only norm-conserving PP is supported'
    ENDIF
    READ(iunit,*) ps_upf%nlcc, dummy   ! nlcc     "Nonlinear Core Correction"
    READ(iunit,*) dummy                ! dft      "Exch-Corr"
    READ(iunit,*) ps_upf%z_val, dummy  ! zp       "Z valence"
    READ(iunit,*) dummy                ! etotps   "Total Energy"
    READ(iunit,*) dummy                ! ecutwfc, ecutrho     "Suggested Cutoff for wfc and rho"
    READ(iunit,*) dummy                ! lmax     "Max angular momentum component", THIS IS NOT THE LMAX WE NEED
    READ(iunit,*) np, dummy            ! mesh     "Number of points in mesh"
    READ(iunit,*) ps_upf%n_wfs, ps_upf%n_proj, dummy !  natwfc, nbeta   "Number of wavefunctions, projectors"
    READ(iunit,*) dummy                ! "Wavefunctions   nl   l   occ"
    ALLOCATE(ps_upf%n(1:ps_upf%n_wfs))
    ALLOCATE(ps_upf%l(1:ps_upf%n_wfs))
    ALLOCATE(ps_upf%occ(1:ps_upf%n_wfs))

    ! els(1)      lchi(1)      oc(1)
    !  ...
    ! els(natwfc) lchi(natwfc) oc(natwfc)
    DO i = 1, ps_upf%n_wfs
      READ(iunit,*) nl, ps_upf%l(i), ps_upf%occ(i)
      ps_upf%n(i) = iachar(nl(1:1)) - iachar('0')
      ! some pseudopotentials do not have a number, but just a letter,
      ! so we assume the level is 1a
      IF(ps_upf%n(i) < 1 .OR. ps_upf%n(i) > 9)  ps_upf%n(i) = 1
    ENDDO

    CALL check_end_tag(iunit, 'PP_HEADER')

    !Mesh info
    CALL init_tag(iunit, 'PP_MESH', .true.)

    CALL init_tag(iunit, 'PP_R', .false.)
    READ(iunit,*) fdummy
    IF(fdummy /= 0.d0) THEN
      ps_upf%np = np + 1
      ip = 2
    ELSE
      ps_upf%np = np
      ip = 1
    ENDIF
    ALLOCATE(ps_upf%r(1:ps_upf%np))
    ps_upf%r(1) = 0.d0
    ps_upf%r(ip) = fdummy
    CALL init_tag(iunit, 'PP_R', .true.)
    READ(iunit,*) (ps_upf%r(ir), ir = ip, ps_upf%np)
    CALL check_end_tag(iunit, 'PP_R')

    CALL init_tag(iunit, 'PP_RAB', .false.)
    ALLOCATE(ps_upf%drdi(1:ps_upf%np))
    READ(iunit,*) (ps_upf%drdi(ir), ir = ip, ps_upf%np)
    CALL check_end_tag(iunit, 'PP_RAB')
    IF(ip == 2) ps_upf%drdi(1) = 0.d0
    CALL check_end_tag(iunit, 'PP_MESH')

    !Non-linear core-corrections
    IF(ps_upf%nlcc) THEN
      CALL init_tag(iunit, 'PP_NLCC', .true.)
      ALLOCATE(ps_upf%core_density(1:ps_upf%np))
      READ(iunit,*) (ps_upf%core_density(ir), ir = ip, ps_upf%np)
      CALL check_end_tag(iunit, 'PP_NLCC')
      IF(ip == 2) ps_upf%core_density(1) = linear_extrapolate(ps_upf%r(1), ps_upf%r(2), &
           ps_upf%r(3), ps_upf%core_density(2), ps_upf%core_density(3))
    ELSE
      !ps_upf%core_density
      ! do nothing
    ENDIF

    !Local component
    CALL init_tag(iunit, 'PP_LOCAL', .true.)
    ALLOCATE(ps_upf%v_local(1:ps_upf%np))
    READ(iunit,*) (ps_upf%v_local(ir), ir = ip, ps_upf%np)
    IF(ip == 2) ps_upf%v_local(1) = linear_extrapolate(ps_upf%r(1), ps_upf%r(2), &
           ps_upf%r(3), ps_upf%v_local(2), ps_upf%v_local(3))
    CALL check_end_tag(iunit, 'PP_LOCAL')

    !Non-local components
    call init_tag(iunit, 'PP_NONLOCAL', .true.)

    ALLOCATE(ps_upf%proj(1:ps_upf%np, 1:ps_upf%n_proj))
    ALLOCATE(ps_upf%proj_l(1:ps_upf%n_proj))
    ALLOCATE(ps_upf%proj_np(1:ps_upf%n_proj))
    ps_upf%proj = 0.d0
    DO i = 1, ps_upf%n_proj
      CALL init_tag(iunit, "PP_BETA", .false.)
      READ(iunit,*) idummy, ps_upf%proj_l(i), dummy
      READ(iunit,*) ps_upf%proj_np(i)
      READ(iunit,*) (ps_upf%proj(ir, i), ir = ip, ps_upf%proj_np(i)+ip-1)
      IF(ip == 2) ps_upf%proj(1, i) = 0.d0 !linear_extrapolate(ps_upf%r(1), ps_upf%r(2), &
      !ps_upf%r(3), ps_upf%proj(2, i), ps_upf%proj(3, i))
      CALL check_end_tag(iunit, 'PP_BETA', ok = ok)
      IF(.NOT. ok) THEN
        ! This is a 'new' UPF file generated by Quantum Espresso that
        ! is not compatible with the Quantum Espresso
        ! specifications. We have to skip one line (one line was read
        ! by check_end_tag already).
        READ(iunit,*) dummy
        CALL check_end_tag(iunit, 'PP_BETA')
      ENDIF
    ENDDO
    
    CALL init_tag(iunit, 'PP_DIJ', .false.)
    ALLOCATE(ps_upf%e(1:ps_upf%n_proj))
    ps_upf%e = M_ZERO
    READ(iunit,*) n_dij, dummy
    DO i = 1, n_dij
      READ(iunit,*) ii, jj, ps_upf%e(ii)
      IF(ii /= jj) THEN
        WRITE(*,*) 'Error while reading pseudopotential data'
        STOP
      ENDIF
    ENDDO
    CALL check_end_tag(iunit, 'PP_DIJ')

    CALL check_end_tag(iunit, 'PP_NONLOCAL')

    !Pseudo wavefunctions
    CALL init_tag(iunit, 'PP_PSWFC', .true.)
    ALLOCATE(ps_upf%wfs(1:ps_upf%np, 1:ps_upf%n_wfs))
    DO i = 1, ps_upf%n_wfs
      READ(iunit,*) dummy
      READ(iunit,*) (ps_upf%wfs(ir, i), ir = ip, ps_upf%np)
      IF (ip == 2) ps_upf%wfs(1, i) = 0.d0 !linear_extrapolate(ps_upf%r(1), ps_upf%r(2), &
      !ps_upf%r(3), ps_upf%wfs(2, i), ps_upf%wfs(3, i))
    ENDDO
    CALL check_end_tag(iunit, "PP_PSWFC")

    !Valence charge
    CALL init_tag(iunit, "PP_RHOATOM", .true.)
    ALLOCATE(ps_upf%rho(1:ps_upf%np))
    READ(iunit,*) (ps_upf%rho(ir), ir = ip, ps_upf%np)
      IF(ip == 2) ps_upf%rho(1) = 0.d0 !linear_extrapolate(ps_upf%r(1), ps_upf%r(2), &
      !ps_upf%r(3), ps_upf%rho(2), ps_upf%rho(3))
    CALL check_end_tag(iunit, "PP_RHOATOM")

    !Extra information
    IF(tag_isdef(iunit, "PP_ADDINFO")) THEN
      ps_upf%kb_nc = 2
      CALL init_tag(iunit, "PP_ADDINFO", .true.)
      DO i = 1, ps_upf%n_wfs
        READ(iunit,*) dummy
      ENDDO
      ALLOCATE(ps_upf%proj_j(1:ps_upf%n_proj))
      DO i = 1, ps_upf%n_proj
        READ(iunit,*) dummy, ps_upf%proj_j(i)
      ENDDO
      READ(iunit,*) dummy
      CALL check_end_tag(iunit, "PP_ADDINFO")
    ELSE 
      !nullify(ps_upf%proj_j)
      ! do nothing
    ENDIF

  END SUBROUTINE ps_upf_file_read




  SUBROUTINE init_tag(iunit, string, go_back)
    INTEGER,          INTENT(IN) :: iunit
    CHARACTER(len=*), INTENT(IN) :: string
    LOGICAL,          INTENT(IN) :: go_back

    INTEGER :: ios
    CHARACTER(len=80) :: string2

    IF(go_back) REWIND(iunit)
    DO
      READ(iunit, *, iostat=ios) string2
      IF(ios > 0) THEN
        WRITE(*,*) 'Error in subroutine init_tag'
        STOP
      ELSEIF (iostat == -1) then
        WRITE(*,*) 'No ', trim(string), ' tag found.'
        WRITE(*,*) 'Please check that this is a valid UPF file.'
        WRITE(*,*) '(Note that version 2.0 or any later version of the UPF file format'
        WRITE(*,*) ' are not yet supported).'
        STOP
      ENDIF
      IF(string_matches('<'//string//'>', string2) ) EXIT
    ENDDO

  END SUBROUTINE init_tag



  SUBROUTINE check_end_tag(iunit, string, ok)
    INTEGER,           intent(in)  :: iunit
    character(len=*),  intent(in)  :: string
    logical, optional, intent(out) :: ok

    integer :: ios
    character(len=80) :: string2
      
    READ(iunit, '(a)', iostat = ios) string2
    IF((.not. string_matches('</'//string//'>', string2)) .or. ios /= 0) THEN
      IF(present(ok)) THEN
        ok = .false.
      ELSE
        WRITE(*,*) 'Could not find closing tag </', string, '>'
      ENDIF
    ELSE
      IF(present(ok)) ok = .true.
    ENDIF

  END SUBROUTINE check_end_tag



  LOGICAL FUNCTION string_matches(string1, string2)
    CHARACTER(len=*), INTENT(IN) :: string1, string2
      
    INTEGER :: l1, l2, l
    
    string_matches = .false.

    l1 = len_trim(string1)
    l2 = len_trim(string2)
    DO l = 1, (l2 - l1 + 1)
      IF(string1(1:l1) == string2(l:(l+l1-1))) THEN
        string_matches = .true.
        EXIT
      ENDIF
    ENDDO
      
  END FUNCTION string_matches



  REAL(8) FUNCTION linear_extrapolate(x0, x1, x2, y1, y2) RESULT(y0)
    REAL(8), INTENT(IN) :: x0, x1, x2, y1, y2
    
    REAL(8) :: mm

    mm = (y2 - y1)/(x2 - x1)
    y0 = y1 + mm*(x0 - x1)
  END FUNCTION linear_extrapolate



  LOGICAL FUNCTION tag_isdef(iunit, string)
    integer,          intent(in) :: iunit
    character(len=*), intent(in) :: string
      
    integer :: ios
    character(len=80) :: string2
    
    rewind(iunit)
    do
      read(iunit, *, iostat = ios) string2
      IF(ios > 0) THEN
        WRITE(*,*) 'Error in tag_isdef'
        STOP
      ELSEIF (ios < 0) THEN
        tag_isdef = .false.
        EXIT
      ENDIF
      if (string_matches("<"//string//">", string2) ) then
        tag_isdef = .true.
        exit
      end if
    end DO

  end function tag_isdef


END MODULE

