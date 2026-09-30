! Coupled 3S1-3D1 deuteron, in the momentum convention of the original
! EM500 interface: T=p**2/m and integral dp p**2 V; there is no extra 2/pi.
module fb_deuteron
  use fb_precision, only: dp
  use fb_quadrature, only: gauss_legendre
  use em500_force, only: initialize_em500, sd_kernel, mass_parameter_MeV
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: deuteron_result, solve_deuteron

  type :: deuteron_result
    integer :: n = 0
    real(dp) :: energy_MeV = 0.0_dp, binding_MeV = 0.0_dp
    real(dp) :: probability_s = 0.0_dp, probability_d = 0.0_dp
    real(dp) :: normalization = 0.0_dp
    real(dp) :: kinetic_MeV = 0.0_dp, potential_MeV = 0.0_dp, residual_MeV = 0.0_dp
    real(dp) :: potential_asymmetry = 0.0_dp, potential_relative_asymmetry = 0.0_dp
    real(dp) :: hamiltonian_asymmetry_MeV = 0.0_dp
    real(dp), allocatable :: p(:), w(:), psi_s(:), psi_d(:)
    real(dp), allocatable :: weighted_eigenvector(:)
    real(dp), allocatable :: potential_kernel(:,:), hamiltonian(:,:)
  end type deuteron_result

  ! LP64 BLAS/LAPACK: default Fortran INTEGER is used intentionally.
  interface
    subroutine dsyev(jobz, uplo, n, a, lda, eigenvalues, work, lwork, info)
      import dp
      character(len=1), intent(in) :: jobz, uplo
      integer, intent(in) :: n, lda, lwork
      real(dp), intent(inout) :: a(lda,*)
      real(dp), intent(out) :: eigenvalues(*), work(*)
      integer, intent(out) :: info
    end subroutine dsyev
  end interface
contains
  subroutine solve_deuteron(n, pmax, result)
    integer, intent(in) :: n
    real(dp), intent(in) :: pmax
    type(deuteron_result), intent(out) :: result
    integer :: m, i, j, a, b, row, col, info, lwork
    real(dp) :: block(2,2), potential_scale, work_query(1)
    real(dp), allocatable :: factor(:), kinetic(:), eigenvalues(:), eigenvectors(:,:), work(:), hy(:)

    if (n < 4 .or. n > 1000) error stop 'solve_deuteron: require 4 <= n <= 1000'
    if (.not. ieee_is_finite(pmax)) error stop 'solve_deuteron: pmax must be finite'
    if (pmax <= 0.0_dp .or. pmax > 1000.0_dp) &
      error stop 'solve_deuteron: require 0 < pmax <= 1000 MeV'
    if (.not. ieee_is_finite(mass_parameter_MeV) .or. mass_parameter_MeV <= 0.0_dp) &
      error stop 'solve_deuteron: invalid twice-reduced-mass parameter'

    m = 2*n
    result%n = n
    allocate(result%p(n), result%w(n), result%psi_s(n), result%psi_d(n))
    allocate(result%weighted_eigenvector(m), result%potential_kernel(m,m), result%hamiltonian(m,m))
    allocate(factor(m), kinetic(m), eigenvalues(m), eigenvectors(m,m), hy(m))
    call gauss_legendre(n, 0.0_dp, pmax, result%p, result%w)
    factor(1:n) = result%p*sqrt(result%w)
    factor(n+1:m) = factor(1:n)
    kinetic(1:n) = result%p**2/mass_parameter_MeV
    kinetic(n+1:m) = kinetic(1:n)

    call initialize_em500()
    ! Evaluate both momentum orders independently: do not impose symmetry
    ! before the reciprocity checks below have validated the actual force.
    ! Index ordering is all S nodes followed by all D nodes.
    do j = 1, n
      do i = 1, n
        call sd_kernel(result%p(i), result%p(j), block)
        if (.not. all(ieee_is_finite(block))) error stop 'solve_deuteron: nonfinite EM500 matrix element'
        do b = 1, 2
          col = (b-1)*n+j
          do a = 1, 2
            row = (a-1)*n+i
            result%potential_kernel(row,col) = block(a,b)
            result%hamiltonian(row,col) = factor(row)*block(a,b)*factor(col)
          end do
        end do
      end do
    end do
    do i = 1, m
      result%hamiltonian(i,i) = result%hamiltonian(i,i)+kinetic(i)
    end do
    potential_scale = max(tiny(1.0_dp), maxval(abs(result%potential_kernel)))
    result%potential_asymmetry = maxval(abs(result%potential_kernel-transpose(result%potential_kernel)))
    result%potential_relative_asymmetry = result%potential_asymmetry/potential_scale
    result%hamiltonian_asymmetry_MeV = &
      maxval(sum(abs(result%hamiltonian-transpose(result%hamiltonian)),dim=2))
    if (result%potential_relative_asymmetry > 1.0e-9_dp) &
      error stop 'solve_deuteron: EM500 reciprocity check failed (relative tolerance 1e-9)'
    if (result%hamiltonian_asymmetry_MeV > 1.0e-8_dp) &
      error stop 'solve_deuteron: weighted Hamiltonian Hermiticity check failed (1e-8 MeV)'
    ! Remove only verified roundoff asymmetry. Keep the original potential
    ! matrix in the result for independent inspection.
    result%hamiltonian = 0.5_dp*(result%hamiltonian+transpose(result%hamiltonian))
    eigenvectors = result%hamiltonian
    call dsyev('V', 'U', m, eigenvectors, m, eigenvalues, work_query, -1, info)
    if (info /= 0) error stop 'solve_deuteron: LAPACK DSYEV workspace query failed'
    if (.not. ieee_is_finite(work_query(1))) error stop 'solve_deuteron: invalid LAPACK workspace size'
    if (work_query(1) < 1.0_dp .or. work_query(1) > real(huge(lwork)-1,dp)) &
      error stop 'solve_deuteron: LAPACK workspace size outside INTEGER range'
    lwork = max(3*m-1, ceiling(work_query(1)))
    allocate(work(lwork))
    call dsyev('V', 'U', m, eigenvectors, m, eigenvalues, work, lwork, info)
    if (info < 0) error stop 'solve_deuteron: invalid argument to LAPACK DSYEV'
    if (info > 0) error stop 'solve_deuteron: LAPACK DSYEV failed to converge'
    if (.not. all(ieee_is_finite(eigenvalues))) error stop 'solve_deuteron: nonfinite eigenvalues'
    result%energy_MeV = eigenvalues(1)
    if (result%energy_MeV >= 0.0_dp) error stop 'solve_deuteron: no negative-energy bound state at this grid'
    result%binding_MeV = -result%energy_MeV
    result%weighted_eigenvector = eigenvectors(:,1)
    if (result%weighted_eigenvector(1) < 0.0_dp) &
      result%weighted_eigenvector = -result%weighted_eigenvector
    result%psi_s = result%weighted_eigenvector(1:n)/factor(1:n)
    result%psi_d = result%weighted_eigenvector(n+1:m)/factor(n+1:m)
    result%probability_s = sum(result%weighted_eigenvector(1:n)**2)
    result%probability_d = sum(result%weighted_eigenvector(n+1:m)**2)
    result%normalization = result%probability_s+result%probability_d
    result%kinetic_MeV = sum(kinetic*result%weighted_eigenvector**2)
    hy = matmul(result%hamiltonian,result%weighted_eigenvector)
    result%residual_MeV = sqrt(sum((hy-result%energy_MeV*result%weighted_eigenvector)**2))
    result%potential_MeV = dot_product(result%weighted_eigenvector, &
      hy-kinetic*result%weighted_eigenvector)
    if (.not. all(ieee_is_finite(result%weighted_eigenvector))) &
      error stop 'solve_deuteron: nonfinite eigenvector'
  end subroutine solve_deuteron
end module fb_deuteron
