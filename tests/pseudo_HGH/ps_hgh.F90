!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
!!
!! This program is free software; you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation; either version 2, or (at your option)
!! any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with this program; if not, write to the Free Software
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!
!! $Id: ps_hgh.F90 14362 2015-06-24 17:02:22Z xavier $

module ps_hgh_m
  !< For information about the Hartwinger-Goedecker-Hutter pseudopotentials, take a look at:
  !!  (1) S. Goedecker, M. Teter and J. Hutter, Phys. Rev. B 54, 1703 (1996).
  !!  (2) C. Hartwinger, S. Goedecker and J. Hutter, Phys. Rev. B 58, 3641 (1998).
  use atomic_m
  use logrid_m
  use m_constants, only : PI

  implicit none

  !> The following data type contains:
  !!   (a) the pseudopotential parameters, as read from a *.hgh file,
  !!   (b) auxiliary intermediate functions, to store stuff before passing it to the "ps" variable.
  type hgh_t
    !< HGH parameters.
    character(len=5) :: atom_name
    integer          :: z_val
    REAL(8)            :: rlocal
    REAL(8)            :: rc(0:3)
    REAL(8)            :: c(1:4)
    REAL(8)            :: h(0:3, 1:3, 1:3)
    REAL(8)            :: k(0:3, 1:3, 1:3)

    type(valconf_t)  :: conf
    integer          :: l_max     !< Maximum l for the Kleinman-Bylander component.

    REAL(8), pointer   :: vlocal(:) !< Local potential
    REAL(8), pointer   :: kb(:,:,:) !< KB projectors
    REAL(8), pointer   :: kbr(:)    !< KB radii
    REAL(8), pointer   :: rphi(:,:), eigen(:)

    !> Logarithmic grid parameters
    type(logrid_t) :: g
  end type hgh_t

  REAL(8), parameter :: eps = (1.0e-8)

  interface vlocalr
    module procedure vlocalr_scalar, vlocalr_vector
  end interface vlocalr

  interface projectorr
    module procedure projectorr_scalar, projectorr_vector
  end interface projectorr

contains

  ! ---------------------------------------------------------
  subroutine hgh_init(psp, filename)
    type(hgh_t),    intent(inout) :: psp
    character(len=*), intent(in)  :: filename

    integer :: iunit, i

    write(*,*) 'Reading pseudopotential from file:', trim(filename)

    iunit = 333
    open(iunit, file=trim(filename), action='read', form='formatted', status='old')
    i = load_params(iunit, psp)
    if(i /= 0) then
      write(*,*) 'Error reading hgh file'
      stop
    end if
    close(iunit)

    ! Finds out psp%l_max. The most special cases are H, He, Li_sc and Be_sc, where psp%l_max = -1.
    psp%l_max = 0
    do while(psp%rc(psp%l_max) > (0.01))
      psp%l_max = psp%l_max + 1
      if(psp%l_max > 3) exit
    end do
    psp%l_max = psp%l_max - 1

    ! Initializes the logarithmic grid. Parameters are hard-coded.
    call logrid_init(psp%g, LOGRID_PSF, 3.0d-2, 4.0d-4, 431)

    ! Allocation of stuff.
    ALLOCATE(psp%vlocal(1:psp%g%nrval))
    psp%vlocal = 0.d0
    if(psp%l_max >= 0) then
      ALLOCATE(psp%kbr(0:psp%l_max))
      ALLOCATE(psp%kb(1:psp%g%nrval, 0:psp%l_max, 1:3))
      psp%kbr = 0.d0
      psp%kb = 0.d0
    end if
    ALLOCATE(psp%rphi(1:psp%g%nrval, 1:psp%conf%p))
    ALLOCATE(psp%eigen(1:psp%conf%p))
    psp%rphi = 0.d0
    psp%eigen = 0.d0

  end subroutine hgh_init


  ! ---------------------------------------------------------
  subroutine hgh_end(psp)
    type(hgh_t), intent(inout) :: psp

    if(psp%l_max >= 0) then
      DEALLOCATE(psp%kbr)
      DEALLOCATE(psp%kb)
    end if
    DEALLOCATE(psp%vlocal)
    DEALLOCATE(psp%rphi)
    DEALLOCATE(psp%eigen)
    call logrid_end(psp%g)

  end subroutine hgh_end


  ! ---------------------------------------------------------
  subroutine hgh_process(psp)
    type(hgh_t), intent(inout) :: psp

    integer :: l, i
    REAL(8), pointer :: ptr(:)

    ! Fixes the local potential
    ptr => vlocalr(psp%g%rofi, psp)
    psp%vlocal(1:psp%g%nrval) = ptr(1:psp%g%nrval)
    DEALLOCATE(ptr)

    ! And the projectors
    do l = 0, psp%l_max
      do i = 1, 3
        ptr => projectorr(psp%g%rofi, psp, i, l)
        psp%kb(1:psp%g%nrval, l, i) = ptr(1:psp%g%nrval)
        DEALLOCATE(ptr)
      end do
    end do

    ! get the pseudoatomic eigenfunctions (WARNING: This is not correctly done yet: "some" wavefunctions
    ! are obtained, but not the real ones!!!
    !call solve_schroedinger(psp, ierr)
    !if(ierr /= 0) then ! If the wavefunctions could not be found, we set its number to zero.
    !  write(message(1),'(a)') 'The algorithm that calculates atomic wavefunctions could not'
    !  write(message(2),'(a)') 'do its job. The program will continue, but expect poor'
    !  write(message(3),'(a)') 'convergence properties.'
    !  call messages_warning(3)
    !  psp%conf%p = 0
    !end if

    ! Define the KB-projector cut-off radii
    call get_cutoff_radii(psp)

  end subroutine hgh_process


  ! ---------------------------------------------------------
  function load_params(unit, params)
    integer,     intent(in)  :: unit        ! where to read from
    type(hgh_t), intent(out) :: params      ! obvious
    integer                  :: load_params ! 0 if success,
    ! 1 otherwise.

    integer :: i, iostat, j, k
    character(len=VALCONF_STRING_LENGTH) :: line

    ! Set initially everything to zero.
    params%c(1:4) = 0.d0
    params%rlocal = 0.d0
    params%rc = 0.d0
    params%h = 0.d0
    params%k = 0.d0

    ! get valence configuration
    read(unit,'(a)') line
    call read_valconf(line, params%conf)

    ! Reads the file in a hopefully smart way
    iostat = 1
    j = 5
    read(unit,'(a)') line
    do while((iostat /= 0) .and. (j > 0))
      j = j - 1
      read(line, *, iostat=iostat) params%atom_name, params%z_val, params%rlocal, params%c(1:j)
    end do
    if(j<1) read(line, *, iostat=iostat) params%atom_name, params%z_val, params%rlocal
    if( iostat /= 0 ) then
      load_params = 1
      return
    end if

    read(unit,'(a)', iostat = iostat) line
    if(iostat /= 0) then
      load_params = 0
      return
    end if
    iostat = 1
    j = 4
    do while((iostat /= 0) .and. (j > 0))
      j = j - 1
      read(line, *, iostat=iostat) params%rc(0), (params%h(0, i, i), i = 1, j)
    end do
    if(j < 0) then
      load_params = 2
      return
    end if

    kloop: do k = 1, 3
      read(unit, '(a)', iostat = iostat) line
      if(iostat /= 0) exit kloop
      iostat = 1
      j = 4
      do while((iostat /= 0) .and. (j > 0))
        j = j - 1
        read(line, *, iostat = iostat) params%rc(k), (params%h(k, i, i), i = 1, j)
      end do
      if(params%rc(k) == 0.d0) exit kloop
      read(unit, '(a)') line
      iostat = 1
      j = 4
      do while((iostat /= 0) .and. (j>0))
        j = j - 1
        read(line, *, iostat = iostat) (params%k(k, i, i), i = 1, 3)
      end do
    end do kloop

    ! Fill in the rest of the parameter matrices...
    ! Fill in the rest of the parameter matrices...
    params%h(0, 1, 2) = -0.5d0     * sqrt(3.d0/5.d0)      * params%h(0, 2, 2)
    params%h(0, 1, 3) =  0.5d0     * sqrt(5.d0/21.d0)     * params%h(0, 3, 3)
    params%h(0, 2, 3) = -0.5d0     * sqrt(100.d0/63.d0)   * params%h(0, 3, 3)
    params%h(1, 1, 2) = -0.5d0     * sqrt(5.d0/7.d0)      * params%h(1, 2, 2)
    params%h(1, 1, 3) =  1.d0/6.d0 * sqrt(35.d0/11.d0)    * params%h(1, 3, 3)
    params%h(1, 2, 3) = -1.d0/6.d0 * (14.d0/sqrt(11.d0))  * params%h(1, 3, 3)
    params%h(2, 1, 2) = -0.5d0     * sqrt(7.d0/9.d0)      * params%h(2, 2, 2)
    params%h(2, 1, 3) =  0.5d0     * sqrt(63.d0/143.d0)   * params%h(2, 3, 3)
    params%h(2, 2, 3) = -0.5d0     * (18.d0/sqrt(143.d0)) * params%h(2, 3, 3)

    params%k(0, 1, 2) = -0.5d0     * sqrt(3.d0/5.d0)      * params%k(0, 2, 2)
    params%k(0, 1, 3) =  0.5d0     * sqrt(5.d0/21.d0)     * params%k(0, 3, 3)
    params%k(0, 2, 3) = -0.5d0     * sqrt(100.d0/63.d0)   * params%k(0, 3, 3)
    params%k(1, 1, 2) = -0.5d0     * sqrt(5.d0/7.d0)      * params%k(1, 2, 2)
    params%k(1, 1, 3) =  1.d0/6.d0 * sqrt(35.0/11.d0)     * params%k(1, 3, 3)
    params%k(1, 2, 3) = -1.d0/6.d0 * (14.d0/sqrt(11.d0))  * params%k(1, 3, 3)
    params%k(2, 1, 2) = -0.5d0     * sqrt(7.d0/9.d0)      * params%k(2, 2, 2)
    params%k(2, 1, 3) =  0.5d0     * sqrt(63.d0/143.d0)   * params%k(2, 3, 3)
    params%k(2, 2, 3) = -0.5d0     * (18.d0/sqrt(143.d0)) * params%k(2, 3, 3)


    ! Parameters are symmetric.
    do k = 0, 3
      do i = 1, 3
        do j = i + 1, 3
          params%h(k, j, i) = params%h(k, i, j)
          params%k(k, j, i) = params%k(k, i, j)
        end do
      end do
    end do

    load_params = 0
  end function load_params


  ! ---------------------------------------------------------
  subroutine get_cutoff_radii(psp)
    type(hgh_t), intent(inout)     :: psp

    integer  :: ir, l, i
    REAL(8) :: dincv, tmp
    REAL(8), parameter :: threshold = 1.0d-4

    do l = 0, psp%l_max
      tmp = 0.d0
      do i = 1, 3
        do ir = psp%g%nrval, 2, -1
          dincv = abs(psp%kb(ir, l, i))
          if(dincv > threshold) exit
        end do
        tmp = psp%g%rofi(ir + 1)
        psp%kbr(l) = max(tmp, psp%kbr(l))
      end do
    end do

  end subroutine get_cutoff_radii


  ! ---------------------------------------------------------
  ! Local pseudopotential, both in real and reciprocal space.
  function vlocalr_scalar(r, p)
    type(hgh_t), intent(in)     :: p
    REAL(8), intent(in)           :: r
    REAL(8)                       :: vlocalr_scalar

    REAL(8) :: r1, r2, r4, r6

    r1 = r/p%rlocal
    r2 = r1**2
    r4 = r2**2
    r6 = r4*r2

    if(r < 1.0d-7) then
      vlocalr_scalar = - (2.d0 * p%z_val)/(sqrt(2.d0*PI)*p%rlocal) + p%c(1)
      return
    end if

    ! using erf from intrinsic function
    vlocalr_scalar = - (p%z_val/r)*erf(r1/sqrt(2.d0))   &
      + exp( -0.5d0*r2 ) *    &
      ( p%c(1) + p%c(2)*r2 + p%c(3)*r4 + p%c(4)*r6 )

  end function vlocalr_scalar


  ! ---------------------------------------------------------
  function vlocalr_vector(r, p)
    type(hgh_t), intent(in)      :: p
    REAL(8), intent(in)            :: r(:)
    REAL(8), pointer               :: vlocalr_vector(:)

    integer :: i


    ALLOCATE(vlocalr_vector(1:size(r)))
    do i = 1, size(r)
      vlocalr_vector(i) = vlocalr_scalar(r(i), p)
    end do

  end function vlocalr_vector


  ! ---------------------------------------------------------
  function vlocalg(g, p)
    type(hgh_t), intent(in)     :: p
    REAL(8), intent(in)           :: g
    REAL(8)                       :: vlocalg

    REAL(8) :: g1, g2, g4, g6

    g1 = g*p%rlocal
    g2 = g1*g1
    g4 = g2*g2
    g6 = g4*g2

    vlocalg = -(4.d0*PI*p%z_val/g**2) * exp( -g2/2.d0) + &
      sqrt(8.d0*PI**3) * p%rlocal**3 * exp( -g2/2.d0) * &
      ( p%c(1) + p%c(2)*(3.d0 - g2) + p%c(3)*(15.d0 - 10.d0*g2 + g4) + &
      p%c(4)*(105.d0 - 105.d0*g2 + 21.d0*g4 - g6) )

  end function vlocalg


  ! ---------------------------------------------------------
  function projectorr_scalar(r, p, i, l)
    type(hgh_t), intent(in)     :: p
    REAL(8), intent(in)           :: r
    integer, intent(in)         :: i, l
    REAL(8)                       :: projectorr_scalar

    REAL(8) :: x, y, rr

    x = l + real(4*i-1, kind=8)/2.d0
    !y = loct_gamma(x)
    y = gamma(x)  ! use intrinsic function ?
    x = sqrt(y)
    if(l==0 .and. i==1) then
      rr = 1.d0
    else
      rr = r ** (l + 2*(i-1))
    end if

    projectorr_scalar = sqrt(2.d0) * rr * exp(-r**2/(2.d0*p%rc(l)**2)) / &
      (  p%rc(l)**(l + real(4*i-1, kind=8)/2.d0) * x )

  end function projectorr_scalar


  ! ---------------------------------------------------------
  function projectorr_vector(r, p, i, l)
    type(hgh_t), intent(in)     :: p
    REAL(8), intent(in)           :: r(:)
    integer, intent(in)         :: i, l
    REAL(8), pointer              :: projectorr_vector(:)

    integer :: j

    ALLOCATE(projectorr_vector(1:size(r)))
    do j=1, size(r)
      projectorr_vector(j) = projectorr_scalar(r(j), p, i, l)
    end do

  end function projectorr_vector


  ! ---------------------------------------------------------
  function projectorg(g, p, i, l)
    type(hgh_t), intent(in)     :: p
    REAL(8), intent(in)           :: g
    integer, intent(in)         :: i, l
    REAL(8)                       :: projectorg

    REAL(8) :: pif, ex

    pif = PI**(5.d0/4.d0)

    ex = exp( 0.5d0*(g*p%rc(l))**2 )

    projectorg = 0.d0

    select case(l)
    case(0)
      select case(i)
      case(1)
        projectorg = ( 4.d0*sqrt(2.d0*p%rc(0)**3)*pif ) / ex
      case(2)
        projectorg = ( sqrt((8.0)*2*p%rc(0)**3/(15.0))*pif * &
          (3.d0 - (g*p%rc(0))**2) ) / ex
      case(3)
        projectorg = ( 16.d0*sqrt(2.d0*p%rc(0)**3/105.d0) * pif * &
          (15.d0 - 10.d0*g**2*p%rc(0)**2 + g**4*p%rc(0)**2) ) / (3.d0*ex)
      end select

    case(1)
      select case(i)
      case(1)
        projectorg = ( (8.0)*sqrt(p%rc(1)**5/3.d0)*pif*g ) / ex
      case(2)
        projectorg = ( (16.0)*sqrt(p%rc(1)**5/(105.0))* pif * g * &
          ( 5.d0 - (g*p%rc(1))**2 ) ) / ex
      case(3)
        projectorg = ( (32.0)*sqrt(p%rc(1)**5/(1155.0))* pif * g * &
          ( (35.0) - (14.0)*g**2*p%rc(1)**2 + (g*p%rc(1))**4 ) ) / &
          (3.d0*ex)
      end select

    case(2)
      select case(i)
      case(1)
        projectorg = ( (8.0) * sqrt(2.d0*p%rc(2)**7/(15.0)) * pif * g**2 ) / ex
      case(2)
        projectorg = ( (16.0) * sqrt(2.d0*p%rc(2)**7/(105.0)) * pif * g**2 * &
          ((7.0) - g**2*p%rc(2)**2) ) / (3.d0*ex)
      case(3)
        projectorg = 0.d0 ! ??
      end select

    case(3)
      ! This should be checked. Probably will not be needed in an near future...
    end select

  end function projectorg

  ! ---------------------------------------------------------
!  subroutine hgh_debug(psp, dir)
!    type(hgh_t), intent(in) :: psp
!    character(len=*), intent(in) :: dir
!
!    integer :: hgh_unit, loc_unit, dat_unit, kbp_unit, wav_unit, i, l, k
!    character(len=256) :: dirname
!
!    ! Open files.
!    dirname = trim(dir)//'/hgh.'//trim(psp%atom_name)
!    call io_mkdir(trim(dir))
!    hgh_unit = io_open(trim(dirname)//'/hgh', action='write')
!    loc_unit = io_open(trim(dirname)//'/local', action='write')
!    dat_unit = io_open(trim(dirname)//'/info', action='write')
!    kbp_unit = io_open(trim(dirname)//'/nonlocal', action='write')
!    wav_unit = io_open(trim(dirname)//'/wave', action='write')
!
!    ! Writes down the input file, to be checked against !SHARE_OCTOPUS/pseudopotentials/HGH/ATOM_NAME.hgh
!    write(hgh_unit,'(a5,i6,5f12.6)') psp%atom_name, psp%z_val, psp%rlocal, psp%c(1:4)
!    write(hgh_unit,'(  11x,4f12.6)') psp%rc(0), (psp%h(0,i,i), i = 1, 3)
!    do k = 1, 3
!      write(hgh_unit,'(  11x,4f12.6)') psp%rc(k), (psp%h(k, i, i), i = 1, 3)
!      write(hgh_unit,'(  23x,4f12.6)')            (psp%k(k, i, i), i = 1, 3)
!    end do
!
!    ! Writes down some info.
!    write(dat_unit,'(a,i3)')        'lmax  = ', psp%l_max
!    if(psp%l_max >= 0) then
!      write(dat_unit,'(a,4f14.6)') 'kbr   = ', psp%kbr
!    end if
!    write(dat_unit,'(a,5f14.6)')    'eigen = ', psp%eigen
!    write(dat_unit,'(a,5f14.6)')    'occ   = ', psp%conf%occ(1:psp%conf%p, 1)
!    ! Writes down local part.
!    do i = 1, psp%g%nrval
!      write(loc_unit, *) psp%g%rofi(i), psp%vlocal(i)
!    end do
!
!    ! Writes down nonlocal part.
!    if(psp%l_max >=0) then
!      do i = 1, psp%g%nrval
!        write(kbp_unit, '(10es14.4)') psp%g%rofi(i), ( (psp%kb(i, l, k) ,k = 1, 3),l = 0, !psp%l_max)
!      end do
!    end if
!
!    ! And the pseudo-wavefunctions.
!    do i = 1, psp%g%nrval
!      write(wav_unit, *) psp%g%rofi(i), (psp%rphi(i, l), l = 1, psp%conf%p)
!    end do
!
!    ! Closes files and exits
!    call io_close(hgh_unit)
!    call io_close(loc_unit)
!    call io_close(wav_unit)
!    call io_close(dat_unit)
!    call io_close(kbp_unit)
!
!    POP_SUB(hgh_debug)
!  end subroutine hgh_debug

end module ps_hgh_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
