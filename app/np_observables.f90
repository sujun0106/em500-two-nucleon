! Unpolarized np observables from the original EM500 partial-wave S matrices.
! Angles refer to the outgoing neutron in the CM frame; phi=0.
! q: MeV/c; spin amplitudes: fm; d(sigma)/dOmega: mb/sr; integrated sigma: mb.
! The EM500 radial-rephased basis must be changed to the plain momentum basis:
! S_plain(Lout,Lin) = i**Lout * S_radial(Lout,Lin) * conjg(i**Lin).
module fb_np_observables
  use fb_precision, only: dp
  use em500_parameters, only: hbarc_MeV_fm
  use fb_np_scattering, only: np_channel_result
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: spin_amplitudes, differential_cross_section
  public :: total_elastic_cross_section, optical_total_cross_section
  real(dp), parameter :: pi = acos(-1.0_dp)
  complex(dp), parameter :: imaginary_unit = (0.0_dp,1.0_dp)
contains
  subroutine validate_waves(waves)
    type(np_channel_result), intent(in) :: waves(:)
    integer :: w, previous, a, l, spin, j, nc
    do w = 1, size(waves)
      spin = waves(w)%spin
      j = waves(w)%j
      nc = waves(w)%nc
      if ((spin /= 0 .and. spin /= 1) .or. j < 0 .or. nc < 1 .or. nc > 2) &
        error stop 'np observables: invalid spin, J or channel count'
      do a = 1, nc
        l = waves(w)%ls(a)
        if (l < 0 .or. j < abs(l-spin) .or. j > l+spin) &
          error stop 'np observables: L and spin cannot couple to J'
      end do
      if (nc == 2) then
        if (waves(w)%ls(2) <= waves(w)%ls(1)) &
          error stop 'np observables: L channels must be sorted and distinct'
        if (modulo(waves(w)%ls(2)-waves(w)%ls(1),2) /= 0) &
          error stop 'np observables: each block must conserve parity'
      end if
      if (.not. all(ieee_is_finite(real(waves(w)%s(1:nc,1:nc),dp))) .or. &
          .not. all(ieee_is_finite(aimag(waves(w)%s(1:nc,1:nc))))) &
        error stop 'np observables: nonfinite S matrix'
      do previous = 1, w-1
        if (spin == waves(previous)%spin .and. j == waves(previous)%j .and. &
            modulo(waves(w)%ls(1),2) == modulo(waves(previous)%ls(1),2)) &
          error stop 'np observables: duplicate spin/J/parity block'
      end do
    end do
  end subroutine validate_waves

  function wave_number(q, hbarc) result(k)
    real(dp), intent(in) :: q
    real(dp), intent(in), optional :: hbarc
    real(dp) :: k, hc
    hc = hbarc_MeV_fm
    if (present(hbarc)) hc = hbarc
    if (.not. ieee_is_finite(q) .or. .not. ieee_is_finite(hc)) &
      error stop 'np observables: q and hbarc must be finite'
    if (q <= 0.0_dp .or. hc <= 0.0_dp) &
      error stop 'np observables: q and hbarc must be positive'
    k = q/hc
  end function wave_number

  pure function log_factorial(n) result(value)
    integer, intent(in) :: n
    real(dp) :: value
    value = log_gamma(real(n+1,dp))
  end function log_factorial

  pure function clebsch(j1,m1,j2,m2,j,m) result(value)
    ! Integer angular momenta; Condon--Shortley convention, Racah sum.
    ! For the present spin=0,1 amplitudes there are at most three summands.
    integer, intent(in) :: j1,m1,j2,m2,j,m
    real(dp) :: value, logpref, term, correction, adjusted, updated
    integer :: low, high, z
    value = 0.0_dp
    if (m1+m2 /= m .or. abs(m1) > j1 .or. abs(m2) > j2 .or. abs(m) > j) return
    if (j < abs(j1-j2) .or. j > j1+j2) return
    logpref = 0.5_dp*(log(real(2*j+1,dp)) + &
      log_factorial(j+j1-j2)+log_factorial(j-j1+j2)+log_factorial(j1+j2-j) &
      -log_factorial(j1+j2+j+1)+log_factorial(j+m)+log_factorial(j-m) &
      +log_factorial(j1-m1)+log_factorial(j1+m1)+log_factorial(j2-m2)+log_factorial(j2+m2))
    low = max(0,j2-j-m1,j1+m2-j)
    high = min(j1+j2-j,j1-m1,j2+m2)
    correction = 0.0_dp
    do z = low, high
      term = real((-1)**z,dp)*exp(logpref-log_factorial(z)-log_factorial(j1+j2-j-z) &
        -log_factorial(j1-m1-z)-log_factorial(j2+m2-z) &
        -log_factorial(j-j2+m1+z)-log_factorial(j-j1-m2+z))
      adjusted = term-correction
      updated = value+adjusted
      correction = (updated-value)-adjusted
      value = updated
    end do
  end function clebsch

  pure function spherical_harmonic_phi0(l,m,cosine) result(ylm)
    ! Real Y_lm(theta,0), with associated Legendre Condon--Shortley phase.
    integer, intent(in) :: l,m
    real(dp), intent(in) :: cosine(:)
    real(dp) :: ylm(size(cosine)), pmm(size(cosine)), plm(size(cosine))
    real(dp) :: previous(size(cosine)), next_value(size(cosine)), norm
    integer :: order, a, degree
    order = abs(m)
    ylm = 0.0_dp
    if (order > l) return
    pmm = 1.0_dp
    do a = 1, order
      pmm = -real(2*a-1,dp)*sqrt(max(0.0_dp,1.0_dp-cosine*cosine))*pmm
    end do
    if (l == order) then
      plm = pmm
    else
      previous = pmm
      plm = real(2*order+1,dp)*cosine*pmm
      do degree = order+2, l
        next_value = (real(2*degree-1,dp)*cosine*plm-real(degree+order-1,dp)*previous) &
          /real(degree-order,dp)
        previous = plm
        plm = next_value
      end do
    end if
    norm = sqrt(real(2*l+1,dp)/(4.0_dp*pi)* &
      exp(log_factorial(l-order)-log_factorial(l+order)))
    if (m < 0) norm = norm*real((-1)**order,dp)
    ylm = norm*plm
  end function spherical_harmonic_phi0

  subroutine spin_amplitudes(waves,q,theta,fs,ft,hbarc,convention)
    ! Fixed-axis spin amplitudes, not helicity amplitudes.
    ! ft(:,mout+2,min+2) uses spin magnetic order -1,0,1; fs is the singlet.
    ! Default convention is the original EM500 radial-rephased basis.
    ! For artificial/reference S blocks one may specify 'plain_momentum'.
    type(np_channel_result), intent(in) :: waves(:)
    real(dp), intent(in) :: q,theta(:)
    complex(dp), intent(out) :: fs(size(theta)), ft(size(theta),3,3)
    real(dp), intent(in), optional :: hbarc
    character(len=*), intent(in), optional :: convention
    real(dp) :: k, cosine(size(theta)), c_in, c_out
    complex(dp) :: delta_s, coefficient, prefactor
    integer :: w, a, b, lout, lin, spin, j, mi, mo, ml
    logical :: rephase
    call validate_waves(waves)
    k = wave_number(q,hbarc)
    if (.not. all(ieee_is_finite(theta))) &
      error stop 'np observables: finite CM angles required'
    if (any(theta < 0.0_dp) .or. any(theta > 180.0_dp)) &
      error stop 'np observables: CM angles must be between 0 and 180 degrees'
    rephase = .true.
    if (present(convention)) then
      select case(trim(convention))
      case('radial_rephased')
        rephase = .true.
      case('plain_momentum')
        rephase = .false.
      case default
        error stop 'np observables: unknown orbital phase convention'
      end select
    end if
    cosine = cos(theta*pi/180.0_dp)
    fs = (0.0_dp,0.0_dp)
    ft = (0.0_dp,0.0_dp)
    prefactor = sqrt(4.0_dp*pi)/(2.0_dp*imaginary_unit*k)
    do w = 1, size(waves)
      spin = waves(w)%spin
      j = waves(w)%j
      do a = 1, waves(w)%nc
        lout = waves(w)%ls(a)
        do b = 1, waves(w)%nc
          lin = waves(w)%ls(b)
          delta_s = waves(w)%s(a,b)
          if (rephase) delta_s = imaginary_unit**lout*delta_s*conjg(imaginary_unit**lin)
          if (a == b) delta_s = delta_s-1.0_dp
          coefficient = prefactor*sqrt(real(2*lin+1,dp))*delta_s
          if (abs(coefficient) == 0.0_dp) cycle
          do mi = -spin, spin
            c_in = clebsch(lin,0,spin,mi,j,mi)
            do mo = -spin, spin
              ml = mi-mo
              c_out = clebsch(lout,ml,spin,mo,j,mi)
              if (abs(c_in*c_out) == 0.0_dp) cycle
              if (spin == 0) then
                fs = fs+coefficient*c_in*c_out*spherical_harmonic_phi0(lout,ml,cosine)
              else
                ft(:,mo+2,mi+2) = ft(:,mo+2,mi+2) &
                  +coefficient*c_in*c_out*spherical_harmonic_phi0(lout,ml,cosine)
              end if
            end do
          end do
        end do
      end do
    end do
  end subroutine spin_amplitudes

  function differential_cross_section(waves,q,theta,hbarc,convention) result(sigma)
    ! Spin average: 1/4 singlet and 3/4 triplet; 1 fm^2 = 10 mb.
    type(np_channel_result), intent(in) :: waves(:)
    real(dp), intent(in) :: q,theta(:)
    real(dp), intent(in), optional :: hbarc
    character(len=*), intent(in), optional :: convention
    real(dp) :: sigma(size(theta))
    complex(dp) :: fs(size(theta)), ft(size(theta),3,3)
    integer :: a,b
    call spin_amplitudes(waves,q,theta,fs,ft,hbarc,convention)
    sigma = abs(fs)**2
    do b = 1,3
      do a = 1,3
        sigma = sigma+abs(ft(:,a,b))**2
      end do
    end do
    sigma = 2.5_dp*sigma
  end function differential_cross_section

  function total_elastic_cross_section(waves,q,hbarc) result(sigma)
    ! Angularly integrated elastic cross section directly from the S blocks.
    ! Isospin labels enforce antisymmetry; they are not a degeneracy factor.
    type(np_channel_result), intent(in) :: waves(:)
    real(dp), intent(in) :: q
    real(dp), intent(in), optional :: hbarc
    real(dp) :: sigma, k
    complex(dp) :: delta_s(2,2)
    integer :: w, a, nc
    call validate_waves(waves)
    k = wave_number(q,hbarc)
    sigma = 0.0_dp
    do w = 1, size(waves)
      nc = waves(w)%nc
      delta_s = waves(w)%s
      do a = 1,nc
        delta_s(a,a) = delta_s(a,a)-1.0_dp
      end do
      sigma = sigma+real(2*waves(w)%j+1,dp)*sum(abs(delta_s(1:nc,1:nc))**2)
    end do
    sigma = sigma*10.0_dp*pi/(4.0_dp*k*k)
  end function total_elastic_cross_section

  function optical_total_cross_section(waves,q,hbarc,convention) result(sigma)
    ! Optical theorem. Equals elastic only for unitary S, as in a real NN
    ! Hamiltonian without pion-production or other inelastic channels.
    type(np_channel_result), intent(in) :: waves(:)
    real(dp), intent(in) :: q
    real(dp), intent(in), optional :: hbarc
    character(len=*), intent(in), optional :: convention
    real(dp) :: sigma, k
    complex(dp) :: fs(1), ft(1,3,3)
    integer :: a
    k = wave_number(q,hbarc)
    call spin_amplitudes(waves,q,[0.0_dp],fs,ft,hbarc,convention)
    sigma = aimag(fs(1))
    do a = 1,3
      sigma = sigma+aimag(ft(1,a,a))
    end do
    sigma = sigma*10.0_dp*pi/k
  end function optical_total_cross_section
end module fb_np_observables
