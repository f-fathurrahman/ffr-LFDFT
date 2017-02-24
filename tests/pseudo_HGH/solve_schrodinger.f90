! ---------------------------------------------------------
subroutine solve_schroedinger(psp, ierr)
  type(hgh_t), intent(inout) :: psp
  integer,     intent(out)   :: ierr

  integer :: iter, ir, l, nnode, nprin, i, j, irr, n, k
  REAL(8) :: vtot, a2b4, diff, nonl
  REAL(8), allocatable :: prev(:, :), rho(:, :), ve(:, :)
  REAL(8), parameter :: tol = (1.0e-4)
  REAL(8) :: e, z, dr, rmax
  REAL(8), allocatable :: s(:), hato(:), g(:), y(:)

  ierr = 0

  ! Allocations...
  ALLOCATE(   s(1:psp%g%nrval))
  ALLOCATE(hato(1:psp%g%nrval))
  ALLOCATE(   g(1:psp%g%nrval))
  ALLOCATE(   y(1:psp%g%nrval))
  ALLOCATE(prev(1:psp%g%nrval, 1:1))
  ALLOCATE( rho(1:psp%g%nrval, 1:1))
  ALLOCATE(  ve(1:psp%g%nrval, 1:1))
  hato = 0.d0
  g = 0.d0
  y = 0.d0
  rho = 0.d0
  ve = 0.d0

  ! These numerical parameters have to be fixed for egofv to work.
  s(2:psp%g%nrval) = psp%g%drdi(2:psp%g%nrval)*psp%g%drdi(2:psp%g%nrval)
  s(1) = s(2)
  a2b4 = (1.d0/4.d0)*psp%g%a**2

  ! Let us be a bit informative.
  write(*,*) 'Calculating atomic pseudo-eigenfunctions for species ' // psp%atom_name // '....'

  ! "Double" self-consistent loop: nonlocal and XC parts have to be calculated
  ! self-consistently.
  diff = 1d5
  iter = 0
  self_consistent: do while( diff > tol )
    prev = rho
    iter = iter + 1
    do n = 1, psp%conf%p
      l = psp%conf%l(n)
      do ir = 2, psp%g%nrval
        vtot = 2*psp%vlocal(ir) + ve(ir, 1) + dble(l*(l + 1))/(psp%g%rofi(ir)**2)
        nonl = 0.d0
        if(iter>2 .and. psp%l_max >=0 .and. psp%rphi(ir, n) > (1.0e-7)) then
          do i = 1, 3
            do j = 1, 3
              do irr = 2, psp%g%nrval
                nonl = nonl + psp%h(l, i, j)*psp%kb(ir, l, i)* &
                  psp%g%drdi(irr)*psp%g%rofi(irr)*psp%rphi(irr, n)*psp%kb(irr,l,j)
              end do
            end do
          end do
          nonl = 2*nonl/psp%rphi(ir, n)*psp%g%rofi(ir)
        end if
        vtot = vtot + nonl
        hato(ir) = vtot*s(ir) + a2b4
      end do
      hato(1) = hato(2)
      ! We will assume there is only the possibility of two equal l numbers...
      nnode = 1
      do k = 1, n - 1
        if(psp%conf%l(k)==psp%conf%l(n)) nnode = 2
      end do
      nprin = l + 1
      if(iter == 1) then
        e = -((psp%z_val/dble(nprin))**2)
        z = psp%z_val
      else
        e = psp%eigen(n)
        z = psp%z_val
      end if
      dr = -(1.0e5)
      rmax = psp%g%rofi(psp%g%nrval)
      call egofv(hato, s, psp%g%nrval, e, g, y, l, z, &
        real(psp%g%a, 8), real(psp%g%b, 8), rmax, nprin, nnode, dr, ierr)
      if(ierr /= 0) exit self_consistent
      psp%eigen(n) = e

      psp%rphi(2:psp%g%nrval, n) = g(2:psp%g%nrval) * sqrt(psp%g%drdi(2:psp%g%nrval))
      psp%rphi(1, n) = psp%rphi(2, n)
    end do
    rho = 0.d0
    do n = 1, psp%conf%p
      rho(1:psp%g%nrval, 1) = rho(1:psp%g%nrval, 1) + psp%conf%occ(n,1)*psp%rphi(1:psp%g%nrval, n)**2
    end do
    if(iter>1) rho = 0.5d0*(rho + prev)
    diff = sqrt(sum(psp%g%drdi(2:psp%g%nrval)*(rho(2:psp%g%nrval, 1)-prev(2:psp%g%nrval, 1))**2))
    call atomhxc('LDA', psp%g, 1, rho, ve)

  end do self_consistent

  if(ierr  ==  0) then
    !  checking normalization of the calculated wavefunctions
    !do l = 0, psp%l_max_occ
    do n = 1, psp%conf%p
      e = sqrt(sum(psp%g%drdi(2:psp%g%nrval)*psp%rphi(2:psp%g%nrval, n)**2))
      e = abs(e - 1.d0)
      if (e > (1.0e-5)) then
        write(*, '(a,i2,a)') "Eigenstate for n = ", n , ' is not normalized'
        write(*, '(a, f12.6,a)') '(abs(1-norm) = ', e, ')'
      end if
    end do

    ! Output in Ha and not in stupid Rydbergs.
    psp%eigen = psp%eigen / 2.d0
  end if

end subroutine solve_schroedinger
