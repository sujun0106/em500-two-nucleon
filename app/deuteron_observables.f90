! Coordinate-space, one-body/point-particle observables of the deuteron.
! This is a postprocessor of the UNCHANGED EM500 momentum-space solution.
! r is the neutron-proton separation; R=r/2 is a nucleon's distance from
! the equal-mass centre of mass. Neither density is an electromagnetic density.
module fb_deuteron_observables
  use fb_precision, only: dp
  use fb_deuteron, only: deuteron_result
  use em500_force, only: mass_parameter_MeV, hbarc_MeV_fm
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: deuteron_observables_result, calculate_deuteron_observables, write_deuteron_observables
  real(dp), parameter :: pi=acos(-1.0_dp)
  ! Measured one-nucleon properties used ONLY by one-body estimates;
  ! these constants do not modify the fixed EM500 interaction.
  real(dp), parameter :: proton_mu_muN=2.79284734463_dp, neutron_mu_muN=-1.91304276_dp
  real(dp), parameter :: proton_charge_radius_fm=0.84075_dp, neutron_charge_radius_squared_fm2=-0.1155_dp
  real(dp), parameter :: proton_mass_MeV=938.27208943_dp

  type :: deuteron_observables_result
    integer :: n_radial_intervals=0
    real(dp) :: rmax_fm=0, dr_fm=0, binding_MeV=0, gamma_fm_inverse=0
    real(dp) :: coordinate_norm=0, momentum_norm=0, norm_difference=0
    real(dp) :: coordinate_d_probability=0, momentum_d_probability=0
    real(dp) :: relative_radius_fm=0, point_radius_fm=0, quadrupole_IA_fm2=0
    real(dp) :: magnetic_moment_IA_muN=0, charge_radius_estimate_fm=0, darwin_foldy_fm2=0
    real(dp) :: asymptotic_s_fm_mhalf=0, asymptotic_d_to_s=0
    real(dp) :: tail_fit_rmin_fm=20, tail_fit_rmax_fm=25
    real(dp) :: tail_fit_s_relative_spread=0, tail_fit_d_relative_spread=0
    logical :: tail_fit_available=.false.
    real(dp), allocatable :: r(:), u(:), w(:), radial_probability(:)
    real(dp), allocatable :: R_cm(:), density_per_nucleon(:), density_total_nucleons(:)
  end type deuteron_observables_result
contains
  pure real(dp) function spherical_j0(x) result(y)
    real(dp), intent(in) :: x
    real(dp) :: z
    if (abs(x)<0.1_dp) then
      z=x*x
      y=1-z/6+z*z/120-z*z*z/5040+z**4/362880
    else
      y=sin(x)/x
    end if
  end function spherical_j0

  pure real(dp) function spherical_j2(x) result(y)
    real(dp), intent(in) :: x
    real(dp) :: z
    if (abs(x)<0.1_dp) then
      z=x*x
      y=z/15-z*z/210+z*z*z/7560-z**4/498960
    else
      y=(3/x**2-1)*sin(x)/x-3*cos(x)/x**2
    end if
  end function spherical_j2

  pure real(dp) function simpson(y,h) result(integral)
    real(dp), intent(in) :: y(0:),h
    integer :: n
    n=ubound(y,1)
    integral=h/3*(y(0)+y(n)+4*sum(y(1:n-1:2))+2*sum(y(2:n-2:2)))
  end function simpson

  subroutine calculate_deuteron_observables(momentum,n_radial_intervals,rmax_fm,obs)
    type(deuteron_result), intent(in) :: momentum
    integer, intent(in) :: n_radial_intervals
    real(dp), intent(in) :: rmax_fm
    type(deuteron_observables_result), intent(out) :: obs
    real(dp) :: factor,x,derivative_u_at_zero,as_fit,ad_fit,shape_s,shape_d,den_s,den_d
    real(dp) :: as_min,as_max,ad_min,ad_max,rr,z,mean_r2
    integer :: i,k,n,ntail
    if (n_radial_intervals<100 .or. n_radial_intervals>1000000 .or. modulo(n_radial_intervals,2)/=0) &
      error stop 'Deuteron observables: radial interval count must be even and in [100,1000000]'
    if (.not. ieee_is_finite(rmax_fm)) error stop 'Deuteron observables: nonfinite rmax'
    if (rmax_fm<=0 .or. rmax_fm>1000) error stop 'Deuteron observables: require 0 < rmax <= 1000 fm'
    if (momentum%n<4 .or. .not. allocated(momentum%p)) error stop 'Deuteron observables: missing momentum state'
    if (momentum%binding_MeV<=0) error stop 'Deuteron observables: binding energy must be positive'
    if (.not. all(ieee_is_finite(momentum%psi_s)) .or. .not. all(ieee_is_finite(momentum%psi_d))) &
      error stop 'Deuteron observables: nonfinite wavefunction'
    n=n_radial_intervals
    obs%n_radial_intervals=n
    obs%rmax_fm=rmax_fm
    obs%dr_fm=rmax_fm/real(n,dp)
    obs%binding_MeV=momentum%binding_MeV
    obs%gamma_fm_inverse=sqrt(mass_parameter_MeV*momentum%binding_MeV)/hbarc_MeV_fm
    allocate(obs%r(0:n),obs%u(0:n),obs%w(0:n),obs%radial_probability(0:n))
    allocate(obs%R_cm(0:n),obs%density_per_nucleon(0:n),obs%density_total_nucleons(0:n))
    factor=sqrt(2/pi)/hbarc_MeV_fm**1.5_dp
    obs%u=0; obs%w=0
    ! EM500's psi_rad is related to the plain momentum angular basis by
    ! psi_plain,L=i**L psi_rad,L (also used by np_observables.f90).
    ! Fourier e**(i p.r/hbar) supplies a further i**L, hence (-1)**L.
    ! L=0 and L=2 are BOTH EVEN: both Hankel integrals have a PLUS sign.
    ! Adding an isolated minus sign to the D wave would mix conventions.
    do k=0,n
      obs%r(k)=k*obs%dr_fm
      do i=1,momentum%n
        x=momentum%p(i)*obs%r(k)/hbarc_MeV_fm
        obs%u(k)=obs%u(k)+momentum%w(i)*momentum%p(i)**2*momentum%psi_s(i)*spherical_j0(x)
        obs%w(k)=obs%w(k)+momentum%w(i)*momentum%p(i)**2*momentum%psi_d(i)*spherical_j2(x)
      end do
      obs%u(k)=factor*obs%r(k)*obs%u(k)
      obs%w(k)=factor*obs%r(k)*obs%w(k)
    end do
    obs%radial_probability=obs%u**2+obs%w**2
    obs%R_cm=obs%r/2
    ! Per nucleon: 4*pi*R^2*rho(R)*dR = [u(r)^2+w(r)^2]*dr, r=2R.
    ! At R=0 use the analytic limiting S-wave slope, avoiding 0/0.
    derivative_u_at_zero=factor*sum(momentum%w*momentum%p**2*momentum%psi_s)
    obs%density_per_nucleon(0)=2/pi*derivative_u_at_zero**2
    obs%density_per_nucleon(1:n)=2/pi*obs%radial_probability(1:n)/obs%r(1:n)**2
    obs%density_total_nucleons=2*obs%density_per_nucleon
    obs%momentum_norm=sum(momentum%w*momentum%p**2*(momentum%psi_s**2+momentum%psi_d**2))
    obs%momentum_d_probability=sum(momentum%w*momentum%p**2*momentum%psi_d**2)
    obs%coordinate_norm=simpson(obs%radial_probability,obs%dr_fm)
    obs%norm_difference=obs%coordinate_norm-obs%momentum_norm
    obs%coordinate_d_probability=simpson(obs%w**2,obs%dr_fm)
    ! Do NOT renormalize the finite coordinate grid: its norm loss is diagnostic.
    mean_r2=simpson(obs%r**2*obs%radial_probability,obs%dr_fm)
    obs%relative_radius_fm=sqrt(mean_r2)
    obs%point_radius_fm=obs%relative_radius_fm/2
    obs%quadrupole_IA_fm2=simpson(obs%r**2*obs%w*(sqrt(8.0_dp)*obs%u-obs%w),obs%dr_fm)/20
    obs%magnetic_moment_IA_muN=proton_mu_muN+neutron_mu_muN- &
      1.5_dp*(proton_mu_muN+neutron_mu_muN-0.5_dp)*obs%momentum_d_probability
    obs%darwin_foldy_fm2=0.75_dp*(hbarc_MeV_fm/proton_mass_MeV)**2
    obs%charge_radius_estimate_fm=sqrt(obs%point_radius_fm**2+proton_charge_radius_fm**2+ &
      neutron_charge_radius_squared_fm2+obs%darwin_foldy_fm2)
    ! Tail amplitudes: fixed-gamma linear least squares over 20 <= r <= 25 fm.
    ! Pointwise coefficient ranges expose residual nuclear tails/quadrature noise.
    as_fit=0; ad_fit=0; den_s=0; den_d=0; ntail=0
    as_min=huge(1.0_dp); as_max=-huge(1.0_dp)
    ad_min=huge(1.0_dp); ad_max=-huge(1.0_dp)
    if (rmax_fm>=obs%tail_fit_rmax_fm) then
      do k=0,n
        rr=obs%r(k)
        if (rr<obs%tail_fit_rmin_fm .or. rr>obs%tail_fit_rmax_fm) cycle
        z=obs%gamma_fm_inverse*rr
        shape_s=exp(-z)
        shape_d=shape_s*(1+3/z+3/z**2)
        as_fit=as_fit+obs%u(k)*shape_s; den_s=den_s+shape_s**2
        ad_fit=ad_fit+obs%w(k)*shape_d; den_d=den_d+shape_d**2
        as_min=min(as_min,obs%u(k)/shape_s); as_max=max(as_max,obs%u(k)/shape_s)
        ad_min=min(ad_min,obs%w(k)/shape_d); ad_max=max(ad_max,obs%w(k)/shape_d)
        ntail=ntail+1
      end do
      if (ntail>=5 .and. den_s>tiny(1.0_dp) .and. den_d>tiny(1.0_dp)) then
        as_fit=as_fit/den_s; ad_fit=ad_fit/den_d
        if (abs(as_fit)>tiny(1.0_dp)) then
          obs%tail_fit_available=.true.
          obs%asymptotic_s_fm_mhalf=as_fit
          obs%asymptotic_d_to_s=ad_fit/as_fit
          obs%tail_fit_s_relative_spread=(as_max-as_min)/abs(as_fit)
          obs%tail_fit_d_relative_spread=(ad_max-ad_min)/max(abs(ad_fit),tiny(1.0_dp))
        end if
      end if
    end if
    if (.not. all(ieee_is_finite(obs%radial_probability))) error stop 'Deuteron observables: nonfinite transform'
  end subroutine calculate_deuteron_observables

  subroutine write_deuteron_observables(obs,prefix)
    type(deuteron_observables_result), intent(in) :: obs
    character(len=*), intent(in), optional :: prefix
    character(len=:), allocatable :: base
    integer :: unit,k
    base='deuteron'
    if (present(prefix)) base=prefix
    open(newunit=unit,file=base//'_observables.dat',status='replace',action='write')
    write(unit,'(a)') '# Fixed EM500; point-particle/equal-mass one-body observables. No charge-density folding or two-body currents.'
    write(unit,'(a)') '# charge_radius_estimate: finite nucleon sizes + Darwin-Foldy only; omits spin-orbit/boost/two-body terms.'
    write(unit,'(a)') '# Coordinate integrals are NOT renormalized; norm_difference diagnoses transform/tail errors.'
    write(unit,'(a)') '# D-state probability is representation dependent, not an experimental observable.'
    write(unit,'(a,i0)') 'n_radial_intervals = ',obs%n_radial_intervals
    write(unit,'(a,es25.16e3)') 'rmax_relative_fm = ',obs%rmax_fm
    write(unit,'(a,es25.16e3)') 'dr_relative_fm = ',obs%dr_fm
    write(unit,'(a,es25.16e3)') 'binding_MeV = ',obs%binding_MeV
    write(unit,'(a,es25.16e3)') 'gamma_fm_inverse = ',obs%gamma_fm_inverse
    write(unit,'(a,es25.16e3)') 'momentum_norm = ',obs%momentum_norm
    write(unit,'(a,es25.16e3)') 'coordinate_norm = ',obs%coordinate_norm
    write(unit,'(a,es25.16e3)') 'norm_difference = ',obs%norm_difference
    write(unit,'(a,es25.16e3)') 'momentum_d_probability = ',obs%momentum_d_probability
    write(unit,'(a,es25.16e3)') 'coordinate_d_probability = ',obs%coordinate_d_probability
    write(unit,'(a,es25.16e3)') 'relative_radius_fm = ',obs%relative_radius_fm
    write(unit,'(a,es25.16e3)') 'point_radius_fm = ',obs%point_radius_fm
    write(unit,'(a,es25.16e3)') 'quadrupole_IA_fm2 = ',obs%quadrupole_IA_fm2
    write(unit,'(a,es25.16e3)') 'magnetic_moment_IA_muN = ',obs%magnetic_moment_IA_muN
    write(unit,'(a,es25.16e3)') 'charge_radius_estimate_fm = ',obs%charge_radius_estimate_fm
    write(unit,'(a,es25.16e3)') 'proton_mu_muN = ',proton_mu_muN
    write(unit,'(a,es25.16e3)') 'neutron_mu_muN = ',neutron_mu_muN
    write(unit,'(a,es25.16e3)') 'proton_charge_radius_fm = ',proton_charge_radius_fm
    write(unit,'(a,es25.16e3)') 'neutron_charge_radius_squared_fm2 = ',neutron_charge_radius_squared_fm2
    write(unit,'(a,es25.16e3)') 'proton_mass_for_DF_MeV = ',proton_mass_MeV
    write(unit,'(a,es25.16e3)') 'darwin_foldy_fm2 = ',obs%darwin_foldy_fm2
    write(unit,'(a,l1)') 'tail_fit_available = ',obs%tail_fit_available
    write(unit,'(a,es25.16e3)') 'asymptotic_s_fm_mhalf = ',obs%asymptotic_s_fm_mhalf
    write(unit,'(a,es25.16e3)') 'asymptotic_d_to_s = ',obs%asymptotic_d_to_s
    write(unit,'(a,es25.16e3)') 'tail_fit_rmin_fm = ',obs%tail_fit_rmin_fm
    write(unit,'(a,es25.16e3)') 'tail_fit_rmax_fm = ',obs%tail_fit_rmax_fm
    write(unit,'(a,es25.16e3)') 'tail_fit_s_relative_spread = ',obs%tail_fit_s_relative_spread
    write(unit,'(a,es25.16e3)') 'tail_fit_d_relative_spread = ',obs%tail_fit_d_relative_spread
    close(unit)
    open(newunit=unit,file=base//'_coordinate_wavefunction.csv',status='replace',action='write')
    write(unit,'(a)') 'r_relative_fm,u_fm_mhalf,w_fm_mhalf,radial_probability_fm_inverse,'// &
      'R_cm_fm,rho_per_nucleon_fm_m3,rho_total_nucleons_fm_m3'
    do k=0,obs%n_radial_intervals
      write(unit,'(*(es25.16e3,:,","))') obs%r(k),obs%u(k),obs%w(k),obs%radial_probability(k), &
        obs%R_cm(k),obs%density_per_nucleon(k),obs%density_total_nucleons(k)
    end do
    close(unit)
  end subroutine write_deuteron_observables
end module fb_deuteron_observables
