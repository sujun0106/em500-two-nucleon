! Gauss-Legendre quadrature on a finite interval; no external dependencies.
module fb_quadrature
  use fb_precision, only: dp
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: gauss_legendre
contains
  subroutine gauss_legendre(n, a, b, x, w)
    integer, intent(in) :: n
    real(dp), intent(in) :: a, b
    real(dp), intent(out) :: x(n), w(n)
    integer :: i, iteration
    real(dp) :: z, znew, pn, pnm1, derivative, halfwidth, pi
    logical :: converged

    if (n < 1) error stop 'gauss_legendre: n must be positive'
    if (.not. ieee_is_finite(a) .or. .not. ieee_is_finite(b)) &
      error stop 'gauss_legendre: finite endpoints required'
    if (b <= a) error stop 'gauss_legendre: require a < b'
    halfwidth = 0.5_dp*(b-a)
    pi = acos(-1.0_dp)
    do i = 1, (n+1)/2
      z = cos(pi*(real(i,dp)-0.25_dp)/(real(n,dp)+0.5_dp))
      converged = .false.
      do iteration = 1, 100
        call legendre_value(n, z, pn, pnm1)
        derivative = real(n,dp)*(z*pn-pnm1)/(z*z-1.0_dp)
        znew = z-pn/derivative
        if (abs(znew-z) <= 4.0_dp*epsilon(z)) then
          z = znew
          converged = .true.
          exit
        end if
        z = znew
      end do
      if (.not. converged) error stop 'gauss_legendre: root iteration did not converge'
      ! For an odd order the middle node is exactly the origin.
      if (i == n+1-i) z = 0.0_dp
      call legendre_value(n, z, pn, pnm1)
      derivative = real(n,dp)*(z*pn-pnm1)/(z*z-1.0_dp)
      x(i) = a+halfwidth*(1.0_dp-z)
      x(n+1-i) = a+halfwidth*(1.0_dp+z)
      w(i) = 2.0_dp*halfwidth/((1.0_dp-z*z)*derivative*derivative)
      w(n+1-i) = w(i)
    end do
    if (.not. all(ieee_is_finite(x)) .or. .not. all(ieee_is_finite(w))) &
      error stop 'gauss_legendre: nonfinite quadrature result'
    if (any(w <= 0.0_dp)) error stop 'gauss_legendre: nonpositive weight'
    if (any(x <= a) .or. any(x >= b)) error stop 'gauss_legendre: node outside open interval'
  end subroutine gauss_legendre

  subroutine legendre_value(n, z, pn, pnm1)
    integer, intent(in) :: n
    real(dp), intent(in) :: z
    real(dp), intent(out) :: pn, pnm1
    integer :: k
    real(dp) :: pold, pnew
    pn = z
    pnm1 = 1.0_dp
    do k = 2, n
      pold = pnm1
      pnm1 = pn
      pnew = ((2.0_dp*real(k,dp)-1.0_dp)*z*pn-real(k-1,dp)*pold)/real(k,dp)
      pn = pnew
    end do
  end subroutine legendre_value
end module fb_quadrature
