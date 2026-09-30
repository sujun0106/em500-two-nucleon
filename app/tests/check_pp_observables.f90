! Physical pp normalization, exchange, Coulomb and spin tests independent of EM500.
program check_pp_observables
  use fb_precision, only: dp
  use em500_parameters, only: hbarc_MeV_fm
  use fb_np_scattering, only: np_channel_result
  use fb_np_observables, only: spin_amplitudes, differential_cross_section
  use fb_pp_observables
  use fb_quadrature, only: gauss_legendre
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  implicit none
  real(dp), parameter :: pi=acos(-1.0_dp)
  complex(dp), parameter :: iu=(0.0_dp,1.0_dp)
  integer :: checks=0
  character(len=64) :: argument
  if (command_argument_count() > 0) then
    call get_command_argument(1,argument)
    call invalid_case(trim(argument))
    error stop 'Invalid-input test returned normally'
  end if
  call gamma_test()
  call mott_test()
  call strong_s_wave_test()
  call coulomb_interference_test()
  call coulomb_p_wave_phase_test()
  call exchange_tensor_test()
  call coupled_test()
  write(*,'(a,i0,a)') 'PASS pp observables: ',checks,' independent checks'
contains
  subroutine require(condition,message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) then
      write(*,'(a)') 'FAIL '//message
      error stop 1
    end if
    checks=checks+1
  end subroutine require

  subroutine identity(wave,spin,j,ls)
    type(np_channel_result), intent(out) :: wave
    integer, intent(in) :: spin,j,ls(:)
    integer :: a
    wave%spin=spin;wave%j=j;wave%ls=0;wave%nc=size(ls);wave%isospin=1
    wave%ls(1:size(ls))=ls
    wave%s=(0.0_dp,0.0_dp)
    do a=1,size(ls)
      wave%s(a,a)=1.0_dp
    end do
  end subroutine identity

  subroutine gamma_test()
    real(dp) :: sigma(0:32),negative(0:32)
    integer :: l
    call coulomb_phases(1.0_dp,32,sigma)
    ! High precision Gamma definition, independently evaluated, frozen here.
    call require(abs(sigma(0)+0.30164032046753319788753165779689654_dp) < 3.e-14_dp, &
      'sigma0 at eta=1 from independently evaluated Gamma(1+i)')
    call require(abs(sigma(1)-0.48375784292991511172812918802297918_dp) < 3.e-14_dp, &
      'sigma1 at eta=1 from independently evaluated Gamma(2+i)')
    call coulomb_phases(-1.0_dp,32,negative)
    call require(maxval(abs(sigma+negative)) < 3.e-14_dp,'Gamma conjugation for real eta')
    call coulomb_phases(0.0_dp,32,sigma)
    call require(maxval(abs(sigma)) == 0.0_dp,'zero-charge Coulomb phases exactly zero')
  end subroutine gamma_test

  subroutine mott_test()
    type(np_channel_result) :: waves(0)
    real(dp) :: theta(7),s2(7),c2(7),observed(7),expected(7),q,eta,k
    complex(dp) :: fs(1),ft(1,3,3)
    theta=[0.5_dp,13.0_dp,45.0_dp,90.0_dp,135.0_dp,167.0_dp,179.5_dp]
    q=60.0_dp;eta=0.23_dp;k=q/hbarc_MeV_fm
    s2=sin(theta*pi/360.0_dp)**2;c2=cos(theta*pi/360.0_dp)**2
    expected=10.0_dp*(eta/(2.0_dp*k))**2 &
      *(1.0_dp/s2**2+1.0_dp/c2**2-cos(eta*log(s2/c2))/(s2*c2))
    observed=pp_differential_cross_section(waves,q,eta,theta)
    call require(maxval(abs(observed/expected-1.0_dp)) < 1.e-12_dp, &
      'pure Coulomb equals analytic unpolarized spin-1/2 Mott cross section')
    call require(maxval(abs(observed/observed(7:1:-1)-1.0_dp)) < 2.e-13_dp, &
      'pure Coulomb exchange theta <-> pi-theta symmetry')
    call pp_spin_amplitudes(waves,q,eta,[90.0_dp],fs,ft)
    call require(maxval(abs(ft)) < 1.e-14_dp,'pure Coulomb triplet amplitude vanishes at 90 degrees')
    call require(abs(observed(4)-10.0_dp*(eta/k)**2) < 2.e-13_dp, &
      '90-degree Mott normalization contains no spurious factor 1/2')
  end subroutine mott_test

  subroutine strong_s_wave_test()
    type(np_channel_result) :: wave(1)
    real(dp) :: phase,q,k,expected,observed(3),x(32),w(32),theta(32),integrated
    phase=0.43_dp;q=100.0_dp;k=q/hbarc_MeV_fm
    call identity(wave(1),0,0,[0])
    wave(1)%s(1,1)=exp(2.0_dp*iu*phase)
    expected=10.0_dp*sin(phase)**2/k**2
    observed=pp_differential_cross_section(wave,q,0.0_dp,[1.0_dp,63.0_dp,179.0_dp])
    call require(maxval(abs(observed/expected-1.0_dp)) < 2.e-14_dp, &
      'strong-only singlet S wave isotropy and pp spin/exchange normalization')
    call gauss_legendre(32,0.0_dp,1.0_dp,x,w)
    theta=acos(x)*180.0_dp/pi
    integrated=2.0_dp*pi*sum(w*pp_differential_cross_section(wave,q,0.0_dp,theta))
    call require(abs(integrated/(2.0_dp*pi*expected)-1.0_dp) < 2.e-14_dp, &
      'identical-particle event cross section integrates one hemisphere')
  end subroutine strong_s_wave_test

  subroutine coulomb_interference_test()
    type(np_channel_result) :: wave(1)
    real(dp) :: phase,q,k,theta(3),s2,c2,expected(3),observed(3),sigma0
    complex(dp) :: a,b,nuclear
    integer :: i
    theta=[30.0_dp,72.0_dp,134.0_dp];phase=0.43_dp;q=100.0_dp;k=q/hbarc_MeV_fm
    sigma0=-0.30164032046753319788753165779689654_dp
    call identity(wave(1),0,0,[0]);wave(1)%s(1,1)=exp(2.0_dp*iu*phase)
    nuclear=2.0_dp*exp(iu*(2.0_dp*sigma0+phase))*sin(phase)/k
    do i=1,3
      s2=sin(theta(i)*pi/360.0_dp)**2;c2=cos(theta(i)*pi/360.0_dp)**2
      a=-exp(iu*(2.0_dp*sigma0-log(s2)))/(2.0_dp*k*s2)
      b=-exp(iu*(2.0_dp*sigma0-log(c2)))/(2.0_dp*k*c2)
      expected(i)=2.5_dp*(abs(nuclear+a+b)**2+3.0_dp*abs(a-b)**2)
    end do
    observed=pp_differential_cross_section(wave,q,1.0_dp,theta)
    call require(maxval(abs(observed/expected-1.0_dp)) < 3.e-14_dp, &
      'analytic singlet S-wave nuclear-Coulomb interference and phase sign')
  end subroutine coulomb_interference_test


  subroutine coulomb_p_wave_phase_test()
    ! Independent central l=1 limit requires sigma_1 in the nuclear part but
    ! sigma_0 in the analytic Coulomb amplitude: a relative phase-sensitive test.
    type(np_channel_result) :: waves(3)
    real(dp) :: phase,q,k,theta(3),s2,c2,expected(3),observed(3),sigma0,sigma1
    complex(dp) :: a,b,nuclear
    integer :: i
    theta=[37.0_dp,72.0_dp,123.0_dp];phase=0.25_dp;q=100.0_dp;k=q/hbarc_MeV_fm
    sigma0=-0.30164032046753319788753165779689654_dp
    sigma1=0.48375784292991511172812918802297918_dp
    call identity(waves(1),1,0,[1]);call identity(waves(2),1,1,[1]);call identity(waves(3),1,2,[1,3])
    do i=1,3
      waves(i)%s(1,1)=exp(2.0_dp*iu*phase)
    end do
    do i=1,3
      s2=sin(theta(i)*pi/360.0_dp)**2;c2=cos(theta(i)*pi/360.0_dp)**2
      a=-exp(iu*(2.0_dp*sigma0-log(s2)))/(2.0_dp*k*s2)
      b=-exp(iu*(2.0_dp*sigma0-log(c2)))/(2.0_dp*k*c2)
      nuclear=6.0_dp*cos(theta(i)*pi/180.0_dp)*exp(iu*(2.0_dp*sigma1+phase))*sin(phase)/k
      expected(i)=2.5_dp*(abs(a+b)**2+3.0_dp*abs(nuclear+a-b)**2)
    end do
    observed=pp_differential_cross_section(waves,q,1.0_dp,theta)
    call require(maxval(abs(observed/expected-1.0_dp)) < 1.e-13_dp, &
      'central P-wave Coulomb interference preserves relative sigma1-sigma0 phase')
  end subroutine coulomb_p_wave_phase_test

  subroutine tensor_matrix(qx,qz,matrix)
    real(dp), intent(in) :: qx,qz
    real(dp), intent(out) :: matrix(4,4)
    real(dp) :: sigq(2,2)
    integer :: a,b,c,d
    sigq(1,:)=[qz,qx];sigq(2,:)=[qx,-qz]
    do b=1,2
      do a=1,2
        do d=1,2
          do c=1,2
            matrix(2*(a-1)+b,2*(c-1)+d)=sigq(a,c)*sigq(b,d)
          end do
        end do
      end do
    end do
  end subroutine tensor_matrix

  subroutine exchange_tensor_test()
    ! Start from the direct Pauli tensor T=(sigma1.Q)(sigma2.Q), Q=pout-pin.
    ! Compare explicitly antisymmetrized momentum-space matrix elements with
    ! Pauli-allowed LSJ expansion. This checks exchange in a spin-dependent force.
    type(np_channel_result) :: waves(3)
    real(dp) :: theta(6),angle,singlet(4,1),triplet(4,3),direct(4,4),exchange(4,4)
    real(dp) :: expected_s(1,1),expected_t(3,3),error
    complex(dp) :: fs(6),ft(6,3,3)
    integer :: i
    call identity(waves(1),0,0,[0]);call identity(waves(2),1,0,[1]);call identity(waves(3),1,1,[1])
    waves(1)%s(1,1)=1.0_dp-iu*pi*(-8.0_dp*pi)
    waves(2)%s(1,1)=1.0_dp-iu*pi*(8.0_dp*pi)
    waves(3)%s(1,1)=1.0_dp-iu*pi*(-16.0_dp*pi/3.0_dp)
    theta=[1.0_dp,13.0_dp,45.0_dp,90.0_dp,143.0_dp,179.0_dp]
    ! q=hbarc makes the observable routine's wave number equal to 1/fm.
    call pp_spin_amplitudes(waves,hbarc_MeV_fm,0.0_dp,theta,fs,ft)
    singlet(:,1)=[0.0_dp,1.0_dp/sqrt(2.0_dp),-1.0_dp/sqrt(2.0_dp),0.0_dp]
    triplet=0.0_dp;triplet(4,1)=1.0_dp;triplet(2:3,2)=1.0_dp/sqrt(2.0_dp);triplet(1,3)=1.0_dp
    error=0.0_dp
    do i=1,size(theta)
      angle=theta(i)*pi/180.0_dp
      call tensor_matrix(sin(angle),cos(angle)-1.0_dp,direct)
      call tensor_matrix(-sin(angle),-cos(angle)-1.0_dp,exchange)
      expected_s=-2.0_dp*pi*pi*matmul(matmul(transpose(singlet),direct+exchange),singlet)
      expected_t=-2.0_dp*pi*pi*matmul(matmul(transpose(triplet),direct-exchange),triplet)
      error=max(error,abs(fs(i)-expected_s(1,1)),maxval(abs(ft(i,:,:)-expected_t)))
    end do
    call require(error < 3.e-12_dp,'tensor force equals explicit Pauli-space fermion exchange')
    write(*,'(a,es13.5)') 'Independent tensor-exchange amplitude error: ',error
  end subroutine exchange_tensor_test

  subroutine coupled_test()
    type(np_channel_result) :: wave(1)
    real(dp) :: angle,rotation(2,2),theta(7),sigma(7)
    complex(dp) :: diag(2,2),fs(7),ft(7,3,3),ns(7),nt(7,3,3)
    call identity(wave(1),1,2,[1,3])
    angle=0.23_dp
    rotation(1,:)=[cos(angle),-sin(angle)];rotation(2,:)=[sin(angle),cos(angle)]
    diag=0.0_dp;diag(1,1)=exp(0.74_dp*iu);diag(2,2)=exp(-0.36_dp*iu)
    wave(1)%s=matmul(matmul(rotation,diag),transpose(rotation))
    theta=[1.0_dp,13.0_dp,45.0_dp,90.0_dp,135.0_dp,167.0_dp,179.0_dp]
    sigma=pp_differential_cross_section(wave,200.0_dp,0.14_dp,theta)
    call require(maxval(abs(sigma/sigma(7:1:-1)-1.0_dp)) < 4.e-13_dp, &
      'coupled P-F nuclear plus Coulomb unpolarized exchange symmetry')
    call pp_spin_amplitudes(wave,200.0_dp,0.0_dp,theta,fs,ft)
    call spin_amplitudes(wave,200.0_dp,theta,ns,nt)
    call require(max(maxval(abs(fs-2.0_dp*ns)),maxval(abs(ft-2.0_dp*nt))) < 1.e-14_dp, &
      'zero-charge allowed nuclear amplitudes double for pp')
  end subroutine coupled_test

  subroutine invalid_case(name)
    character(len=*), intent(in) :: name
    type(np_channel_result) :: wave(1)
    real(dp) :: result(1),eta,q,theta(1)
    call identity(wave(1),0,0,[0])
    eta=0.1_dp;q=100.0_dp;theta=30.0_dp
    select case(name)
    case('invalid_zero_angle');theta=0.0_dp
    case('invalid_pi_angle');theta=180.0_dp
    case('invalid_eta');eta=-0.1_dp
    case('invalid_q');q=0.0_dp
    case('invalid_isospin');wave(1)%isospin=0
    case('invalid_pauli');wave(1)%ls(1)=1;wave(1)%j=1
    case('invalid_s');wave(1)%s(1,1)=cmplx(ieee_value(1.0_dp,ieee_quiet_nan),0.0_dp,dp)
    case('invalid_nan_angle');theta=ieee_value(1.0_dp,ieee_quiet_nan)
    case default;error stop 'Unknown invalid-input case'
    end select
    result=pp_differential_cross_section(wave,q,eta,theta)
    write(*,*) result
  end subroutine invalid_case
end program check_pp_observables
