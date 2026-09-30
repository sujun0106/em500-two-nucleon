! pp observables from Coulomb-modified nuclear S matrices in the EM500 radial basis.
! The input S_C is relative to Coulomb waves: pure Coulomb therefore has S_C=I.
! Point-Coulomb plus strong interaction only; eta is supplied by the solver.
! q: MeV/c, theta: CM degrees, amplitudes: fm, cross section: mb/sr.
module fb_pp_observables
  use fb_precision, only: dp
  use em500_parameters, only: hbarc_MeV_fm
  use fb_np_scattering, only: np_channel_result
  use fb_np_observables, only: spin_amplitudes
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: pp_spin_amplitudes, pp_differential_cross_section, coulomb_phases
  real(dp), parameter :: pi=acos(-1.0_dp)
  complex(dp), parameter :: iu=(0.0_dp,1.0_dp)
contains
  subroutine coulomb_phases(eta,lmax,sigma)
    ! sigma_l=arg Gamma(l+1+i eta), continued using the Gamma recurrence.
    ! Lanczos g=7, 9 coefficients; real part of Gamma argument is exactly 1.
    ! The overall branch is immaterial for exp[i(sigma_l+sigma_l')].
    real(dp), intent(in) :: eta
    integer, intent(in) :: lmax
    real(dp), intent(out) :: sigma(0:lmax)
    real(dp), parameter :: c(0:8)=[0.99999999999980993_dp,676.5203681218851_dp, &
      -1259.1392167224028_dp,771.32342877765313_dp,-176.61502916214059_dp, &
      12.507343278686905_dp,-0.13857109526572012_dp,9.9843695780195716e-6_dp, &
      1.5056327351493116e-7_dp]
    complex(dp) :: z,x,t,loggamma
    integer :: a
    if (.not. ieee_is_finite(eta) .or. lmax < 0) &
      error stop 'Coulomb phases require finite eta and nonnegative lmax'
    z=iu*eta
    x=cmplx(c(0),0.0_dp,dp)
    do a=1,8
      x=x+c(a)/(z+real(a,dp))
    end do
    t=z+7.5_dp
    loggamma=0.5_dp*log(2.0_dp*pi)+(z+0.5_dp)*log(t)-t+log(x)
    sigma(0)=aimag(loggamma)
    if (eta == 0.0_dp) sigma(0)=0.0_dp
    do a=1,lmax
      sigma(a)=sigma(a-1)+atan2(eta,real(a,dp))
    end do
  end subroutine coulomb_phases

  subroutine pp_spin_amplitudes(waves,q,eta,theta,fs,ft)
    ! Fixed-axis coupled-spin basis: ft(:,mout+2,min+2), m=-1,0,1.
    ! Only Pauli-allowed T=1 waves are accepted: L+S is even.
    ! Identical-particle exchange doubles the allowed nuclear partial waves.
    ! Coulomb singlet is f_C(theta)+f_C(pi-theta), triplet is their difference.
    type(np_channel_result), intent(in) :: waves(:)
    real(dp), intent(in) :: q,eta,theta(:)
    complex(dp), intent(out) :: fs(size(theta)),ft(size(theta),3,3)
    type(np_channel_result) :: effective(size(waves))
    real(dp), allocatable :: sigma(:)
    real(dp) :: k,s2(size(theta)),c2(size(theta))
    complex(dp) :: delta_s,fc(size(theta)),exchange(size(theta))
    integer :: w,a,b,lmax
    if (.not. ieee_is_finite(q) .or. q <= 0.0_dp) &
      error stop 'pp observables require finite positive CM momentum'
    if (.not. ieee_is_finite(eta) .or. eta < 0.0_dp) &
      error stop 'pp observables require finite nonnegative eta'
    if (.not. all(ieee_is_finite(theta))) error stop 'pp angles must be finite'
    ! The physical Coulomb cross section diverges at both exchange endpoints.
    ! Reject them explicitly rather than outputting a fictitious finite value.
    if (any(theta <= 0.0_dp) .or. any(theta >= 180.0_dp)) &
      error stop 'pp CM angles must be strictly inside (0,180) degrees'
    lmax=0
    do w=1,size(waves)
      if (waves(w)%isospin /= 1) error stop 'pp observables require T=1 waves'
      if (waves(w)%nc < 1 .or. waves(w)%nc > 2) error stop 'pp invalid channel count'
      if (any(waves(w)%ls(1:waves(w)%nc) < 0)) error stop 'pp negative orbital momentum'
      do a=1,waves(w)%nc
        if (modulo(waves(w)%ls(a)+waves(w)%spin,2) /= 0) &
          error stop 'pp wave violates Pauli antisymmetry'
        lmax=max(lmax,waves(w)%ls(a))
      end do
    end do
    allocate(sigma(0:lmax))
    call coulomb_phases(eta,lmax,sigma)
    effective=waves
    do w=1,size(waves)
      do b=1,waves(w)%nc
        do a=1,waves(w)%nc
          delta_s=waves(w)%s(a,b)
          if (a == b) delta_s=delta_s-1.0_dp
          effective(w)%s(a,b)=exp(iu*(sigma(waves(w)%ls(a))+sigma(waves(w)%ls(b))))*delta_s
          if (a == b) effective(w)%s(a,b)=effective(w)%s(a,b)+1.0_dp
        end do
      end do
    end do
    ! Existing routine also checks spin/J coupling, duplicate blocks and finite S,
    ! and applies the original EM500 radial-to-plain-momentum orbital phases.
    call spin_amplitudes(effective,q,theta,fs,ft)
    fs=2.0_dp*fs
    ft=2.0_dp*ft
    if (eta == 0.0_dp) return
    k=q/hbarc_MeV_fm
    s2=sin(theta*pi/360.0_dp)**2
    ! Use complementary angles to preserve symmetry near an exchange endpoint.
    c2=sin((180.0_dp-theta)*pi/360.0_dp)**2
    fc=-eta/(2.0_dp*k*s2)*exp(iu*(2.0_dp*sigma(0)-eta*log(s2)))
    exchange=-eta/(2.0_dp*k*c2)*exp(iu*(2.0_dp*sigma(0)-eta*log(c2)))
    fs=fs+fc+exchange
    do a=1,3
      ft(:,a,a)=ft(:,a,a)+fc-exchange
    end do
  end subroutine pp_spin_amplitudes

  function pp_differential_cross_section(waves,q,eta,theta) result(dsdo)
    ! No extra 1/2 in a point differential cross section. A finite integrated
    ! identical-particle event cross section uses one hemisphere (or 1/2 of 4pi).
    ! Point Coulomb has no finite integral over a hemisphere touching theta=0.
    type(np_channel_result), intent(in) :: waves(:)
    real(dp), intent(in) :: q,eta,theta(:)
    real(dp) :: dsdo(size(theta))
    complex(dp) :: fs(size(theta)),ft(size(theta),3,3)
    integer :: a,b
    call pp_spin_amplitudes(waves,q,eta,theta,fs,ft)
    dsdo=abs(fs)**2
    do b=1,3
      do a=1,3
        dsdo=dsdo+abs(ft(:,a,b))**2
      end do
    end do
    dsdo=2.5_dp*dsdo
  end function pp_differential_cross_section
end module fb_pp_observables
