! Regression against frozen original Python-driver / Fortran-EM500 outputs.
! This verification program itself requires only Fortran and LAPACK.
program check_reference
  use fb_precision, only: dp
  use fb_quadrature, only: gauss_legendre
  use fb_deuteron, only: deuteron_result, solve_deuteron
  use em500_force, only: initialize_em500, sd_kernel
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  type(deuteron_result) :: result
  integer :: unit, ios, count, i, n
  real(dp) :: pmax, energy, ps, pd, kinetic, potential
  character(len=512) :: line
  logical :: checked_wavefunction

  call check_quadrature()
  call check_potential()
  call open_reference('app/tests/reference/reference_bound_states.dat', unit)
  read(unit,'(a)',iostat=ios) line
  call require_read(ios)
  read(unit,'(a)',iostat=ios) line
  call require_read(ios)
  read(unit,*,iostat=ios) count
  call require_read(ios)
  if (count /= 8) error stop 'Unexpected number of reference bound-state cases'
  checked_wavefunction = .false.
  write(*,'(a)') '# N pmax_MeV energy_MeV probability_S probability_D kinetic_MeV potential_MeV residual_MeV'
  do i = 1, count
    read(unit,*,iostat=ios) n, pmax, energy, ps, pd, kinetic, potential
    call require_read(ios)
    call solve_deuteron(n, pmax, result)
    call assert_close('bound-state energy', result%energy_MeV, energy, 1.0e-8_dp)
    call assert_close('S probability', result%probability_s, ps, 1.0e-9_dp)
    call assert_close('D probability', result%probability_d, pd, 1.0e-9_dp)
    call assert_close('kinetic expectation', result%kinetic_MeV, kinetic, 1.0e-7_dp)
    call assert_close('potential expectation', result%potential_MeV, potential, 1.0e-7_dp)
    call assert_close('normalization', result%normalization, 1.0_dp, 1.0e-12_dp)
    call assert_close('T+V=E', result%kinetic_MeV+result%potential_MeV, result%energy_MeV, 1.0e-8_dp)
    call assert_close('eigenvector residual', result%residual_MeV, 0.0_dp, 1.0e-8_dp)
    write(*,'(i5,1x,f9.3,7(1x,es23.15))') n, pmax, result%energy_MeV, &
      result%probability_s, result%probability_d, result%kinetic_MeV, &
      result%potential_MeV, result%residual_MeV
    if (n == 100 .and. abs(pmax-1000.0_dp) < 1.0e-12_dp) then
      call check_wavefunction(result)
      checked_wavefunction = .true.
    end if
  end do
  close(unit)
  if (.not. checked_wavefunction) error stop 'N100/pmax1000 wavefunction test was not run'
  write(*,'(a)') 'PASS: quadrature, potential direction/reciprocity, 8 bound states, signed S/D wavefunction.'
contains
  subroutine open_reference(path, reference_unit)
    character(len=*), intent(in) :: path
    integer, intent(out) :: reference_unit
    integer :: status
    open(newunit=reference_unit,file=path,status='old',action='read',iostat=status)
    if (status /= 0) then
      write(*,'(a)') 'Cannot open '//path//'; run this program from the code directory.'
      error stop 'Missing reference data'
    end if
  end subroutine open_reference

  subroutine require_read(status)
    integer, intent(in) :: status
    if (status /= 0) error stop 'Reference file is truncated or malformed'
  end subroutine require_read

  subroutine assert_close(label, actual, expected, tolerance)
    character(len=*), intent(in) :: label
    real(dp), intent(in) :: actual, expected, tolerance
    if (.not. ieee_is_finite(actual) .or. .not. ieee_is_finite(expected)) &
      error stop 'Nonfinite value in reference comparison'
    if (abs(actual-expected) > tolerance) then
      write(*,'(a)') 'FAIL: '//label
      write(*,'(a,es24.16)') 'actual:    ', actual
      write(*,'(a,es24.16)') 'reference: ', expected
      write(*,'(a,es24.16)') 'tolerance: ', tolerance
      error stop 'Reference comparison failed'
    end if
  end subroutine assert_close

  subroutine check_quadrature()
    integer, parameter :: orders(5) = [4,5,20,96,1000]
    integer :: icase, nq, k
    real(dp), allocatable :: x(:), w(:)
    real(dp) :: value, exact, max_error
    max_error = 0.0_dp
    do icase = 1, size(orders)
      nq = orders(icase)
      allocate(x(nq),w(nq))
      call gauss_legendre(nq, 0.0_dp, 1000.0_dp, x, w)
      if (any(x(2:nq) <= x(1:nq-1))) error stop 'Quadrature nodes are not ordered'
      do k = 0, min(2*nq-1,20)
        value = sum(w*(x/1000.0_dp)**k)/1000.0_dp
        exact = 1.0_dp/real(k+1,dp)
        call assert_close('analytic Gauss moment', value, exact, 2.0e-13_dp)
        max_error = max(max_error,abs(value-exact))
      end do
      deallocate(x,w)
    end do
    write(*,'(a,es12.4)') 'PASS: analytic Gauss moments, max error = ', max_error
  end subroutine check_quadrature

  subroutine check_potential()
    integer :: reference_unit, status, ncase, icase, row, col, k
    real(dp) :: pout, pin, ref(2,2), actual(2,2), reverse(2,2), tolerance, max_error
    character(len=512) :: header
    call initialize_em500()
    call open_reference('app/tests/reference/reference_potential.dat', reference_unit)
    do k = 1, 3
      read(reference_unit,'(a)',iostat=status) header
      call require_read(status)
    end do
    read(reference_unit,*,iostat=status) ncase
    call require_read(status)
    if (ncase /= 12) error stop 'Unexpected number of reference potential cases'
    max_error = 0.0_dp
    do icase = 1, ncase
      ! Explicit subscripts are important: text is row-major SS,SD,DS,DD.
      read(reference_unit,*,iostat=status) pout, pin, ref(1,1), ref(1,2), ref(2,1), ref(2,2)
      call require_read(status)
      call sd_kernel(pout,pin,actual)
      call sd_kernel(pin,pout,reverse)
      do col = 1, 2
        do row = 1, 2
          tolerance = 1.0e-18_dp+1.0e-12_dp*abs(ref(row,col))
          call assert_close('potential element (outgoing row / incoming column)', &
            actual(row,col), ref(row,col), tolerance)
          call assert_close('potential reciprocity with exchanged momenta', &
            actual(row,col), reverse(col,row), tolerance)
          max_error = max(max_error,abs(actual(row,col)-ref(row,col)))
        end do
      end do
    end do
    close(reference_unit)
    write(*,'(a,es12.4)') 'PASS: 12 ordered potential samples, max absolute error = ', max_error
  end subroutine check_potential

  subroutine check_wavefunction(state)
    type(deuteron_result), intent(in) :: state
    integer :: reference_unit, status, j
    real(dp) :: p(100), w(100), s(100), d(100), density_s, density_d, norm, distance
    character(len=512) :: header
    call open_reference('app/tests/reference/momentum_wavefunction_N100_pmax1000.csv', reference_unit)
    read(reference_unit,'(a)',iostat=status) header
    call require_read(status)
    do j = 1, 100
      ! List-directed Fortran input accepts the comma-separated numeric rows.
      read(reference_unit,*,iostat=status) p(j), w(j), s(j), d(j), density_s, density_d
      call require_read(status)
      call assert_close('reference momentum node', state%p(j), p(j), 1.0e-10_dp)
      call assert_close('reference quadrature weight', state%w(j), w(j), 1.0e-10_dp)
      call assert_close('reference S density consistency', density_s, p(j)**2*s(j)**2, 1.0e-13_dp)
      call assert_close('reference D density consistency', density_d, p(j)**2*d(j)**2, 1.0e-13_dp)
    end do
    close(reference_unit)
    norm = sum(w*p**2*(s**2+d**2))
    call assert_close('reference wavefunction normalization', norm, 1.0_dp, 1.0e-10_dp)
    if (s(1) <= 0.0_dp .or. state%psi_s(1) <= 0.0_dp) error stop 'Unexpected ground-state global phase'
    ! Compare signed S and D together; independently flipping D cannot pass.
    distance = sqrt(sum(w*p**2*((state%psi_s-s)**2+(state%psi_d-d)**2)))
    call assert_close('signed S/D weighted wavefunction distance', distance, 0.0_dp, 1.0e-8_dp)
    write(*,'(a,es12.4)') 'PASS: N100 signed S/D wavefunction distance = ', distance
  end subroutine check_wavefunction
end program check_reference
