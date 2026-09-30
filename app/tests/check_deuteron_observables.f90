! Independent analytic Hankel-transform tests: Gaussian S and S+D states.
! The mixed-state quadrupole catches an isolated D-wave sign error.
program check_deuteron_observables
  use fb_precision, only: dp
  use fb_quadrature, only: gauss_legendre
  use fb_deuteron, only: deuteron_result
  use fb_deuteron_observables, only: deuteron_observables_result, calculate_deuteron_observables
  use em500_force, only: hbarc_MeV_fm
  implicit none
  type(deuteron_result) :: toy
  type(deuteron_observables_result) :: obs
  real(dp), parameter :: pi=acos(-1.0_dp),sigma=80.0_dp
  real(dp) :: a,exact_u,exact_w,ns,nd,max_s,max_d,rr,q_expected,density_origin
  integer :: k
  toy%n=240
  toy%binding_MeV=2.224575_dp
  allocate(toy%p(toy%n),toy%w(toy%n),toy%psi_s(toy%n),toy%psi_d(toy%n))
  call gauss_legendre(toy%n,0.0_dp,1000.0_dp,toy%p,toy%w)
  toy%psi_s=2/(pi**0.25_dp*sigma**1.5_dp)*exp(-toy%p**2/(2*sigma**2))
  toy%psi_d=0
  a=hbarc_MeV_fm/sigma
  call calculate_deuteron_observables(toy,2400,24.0_dp,obs)
  call near(obs%coordinate_norm,1.0_dp,1.e-11_dp,'S Gaussian Parseval')
  call near(obs%point_radius_fm,sqrt(1.5_dp)*a/2,1.e-11_dp,'S Gaussian point radius')
  call near(obs%quadrupole_IA_fm2,0.0_dp,1.e-14_dp,'S Gaussian quadrupole')
  density_origin=8/(pi**1.5_dp*a**3)
  call near(obs%density_per_nucleon(0),density_origin,1.e-12_dp,'density Jacobian/origin')
  toy%psi_s=toy%psi_s/sqrt(2.0_dp)
  toy%psi_d=4/(sqrt(15.0_dp)*pi**0.25_dp*sigma**3.5_dp)*toy%p**2* &
    exp(-toy%p**2/(2*sigma**2))/sqrt(2.0_dp)
  call calculate_deuteron_observables(toy,2400,24.0_dp,obs)
  call near(obs%coordinate_norm,1.0_dp,1.e-11_dp,'S+D Gaussian Parseval')
  call near(obs%coordinate_d_probability,0.5_dp,1.e-11_dp,'S+D Gaussian D norm')
  call near(obs%point_radius_fm,sqrt(2.5_dp)*a/2,1.e-11_dp,'S+D Gaussian point radius')
  q_expected=(sqrt(120.0_dp)-7)*a*a/80
  call near(obs%quadrupole_IA_fm2,q_expected,1.e-11_dp,'S+D Gaussian quadrupole/relative phase')
  ns=2/(pi**0.25_dp*a**1.5_dp)/sqrt(2.0_dp)
  nd=4/(sqrt(15.0_dp)*pi**0.25_dp*a**3.5_dp)/sqrt(2.0_dp)
  max_s=0; max_d=0
  do k=0,obs%n_radial_intervals
    rr=obs%r(k)
    exact_u=ns*rr*exp(-rr**2/(2*a*a))
    exact_w=nd*rr**3*exp(-rr**2/(2*a*a))
    max_s=max(max_s,abs(obs%u(k)-exact_u))
    max_d=max(max_d,abs(obs%w(k)-exact_w))
  end do
  call near(max_s,0.0_dp,1.e-11_dp,'analytic S Hankel transform')
  call near(max_d,0.0_dp,1.e-11_dp,'analytic D Hankel transform/sign')
  print '(a)', 'PASS: independent analytic Gaussian S and S+D transforms, norm, radius, density and quadrupole.'
contains
  subroutine near(actual,expected,tolerance,label)
    real(dp),intent(in) :: actual,expected,tolerance
    character(len=*),intent(in) :: label
    if (abs(actual-expected)>tolerance) then
      print *,trim(label),actual,expected,tolerance
      error stop 'deuteron observables analytic check failed'
    end if
  end subroutine near
end program check_deuteron_observables
