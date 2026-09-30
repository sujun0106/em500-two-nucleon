! Independent angular-physics checks: no EM500 evaluations required.
! Run normally for the positive suite; invalid_* arguments must exit nonzero.
program check_np_observables
  use fb_precision, only: dp
  use fb_np_scattering, only: np_channel_result
  use fb_np_observables
  use fb_quadrature, only: gauss_legendre
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  implicit none
  real(dp), parameter :: pi = acos(-1.0_dp)
  complex(dp), parameter :: iu = (0.0_dp,1.0_dp)
  character(len=64) :: argument
  integer :: checks
  checks = 0
  if (command_argument_count() > 0) then
    call get_command_argument(1,argument)
    call invalid_case(trim(argument))
    error stop 'Invalid-input test unexpectedly returned normally'
  end if
  call s_wave_test()
  call central_test()
  call coupled_integral_test()
  call tensor_born_test()
  call absorption_test()
  call zero_scattering_test()
  write(*,'(a,i0,a)') 'PASS np observables: ',checks,' independent checks'
contains
  subroutine require(condition,message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) then
      write(*,'(a)') 'FAIL '//message
      error stop 1
    end if
    checks = checks+1
  end subroutine require

  subroutine initialize_wave(wave,spin,j,ls)
    type(np_channel_result), intent(out) :: wave
    integer, intent(in) :: spin,j,ls(:)
    integer :: a
    wave%spin = spin
    wave%j = j
    wave%nc = size(ls)
    wave%ls = 0
    wave%ls(1:size(ls)) = ls
    wave%isospin = modulo(1-ls(1)-spin,2)
    wave%s = (0.0_dp,0.0_dp)
    do a = 1,size(ls)
      wave%s(a,a) = (1.0_dp,0.0_dp)
    end do
  end subroutine initialize_wave

  subroutine phase_wave(wave,spin,j,ls,phases,rotation)
    type(np_channel_result), intent(out) :: wave
    integer, intent(in) :: spin,j,ls(:)
    real(dp), intent(in) :: phases(:)
    real(dp), intent(in), optional :: rotation
    complex(dp) :: diag(2,2)
    real(dp) :: r(2,2), angle
    integer :: a
    call initialize_wave(wave,spin,j,ls)
    diag = (0.0_dp,0.0_dp)
    do a = 1,size(ls)
      diag(a,a) = exp(2.0_dp*iu*phases(a))
    end do
    angle = 0.0_dp
    if (present(rotation)) angle = rotation
    r(1,:) = [cos(angle),-sin(angle)]
    r(2,:) = [sin(angle),cos(angle)]
    if (size(ls) == 2) then
      wave%s = matmul(matmul(r,diag),transpose(r))
    else
      wave%s = diag
    end if
  end subroutine phase_wave

  subroutine s_wave_test()
    type(np_channel_result) :: waves(2)
    real(dp) :: q,ds,dt,expected,observed(5),total,optical
    q = 1.2_dp
    ds = 0.3_dp
    dt = 0.7_dp
    call phase_wave(waves(1),0,0,[0],[ds])
    call phase_wave(waves(2),1,1,[0,2],[dt,0.0_dp])
    expected = 10.0_dp*(sin(ds)**2+3.0_dp*sin(dt)**2)/(4.0_dp*q*q)
    observed = differential_cross_section(waves,q,[0.0_dp,31.0_dp,90.0_dp,159.0_dp,180.0_dp], &
      hbarc=1.0_dp,convention='plain_momentum')
    total = total_elastic_cross_section(waves,q,1.0_dp)
    optical = optical_total_cross_section(waves,q,1.0_dp,'plain_momentum')
    call require(maxval(abs(observed-expected)) < 2.0e-14_dp*expected,'S-wave isotropy, spin weights and units')
    call require(abs(total-4.0_dp*pi*expected) < 1.0e-12_dp,'S-wave analytic elastic total')
    call require(abs(optical-total) < 1.0e-12_dp,'S-wave optical theorem')
    write(*,'(a,es13.5)') 'S-wave maximum differential error: ',maxval(abs(observed-expected))
  end subroutine s_wave_test

  function central_phase(l) result(phase)
    integer, intent(in) :: l
    real(dp) :: phase
    real(dp), parameter :: values(0:6) = [0.5_dp,-0.25_dp,0.11_dp,0.07_dp,-0.03_dp,0.01_dp,0.005_dp]
    phase = 0.0_dp
    if (l >= 0 .and. l <= 6) phase = values(l)
  end function central_phase

  subroutine central_test()
    type(np_channel_result) :: waves(23)
    real(dp) :: theta(6),q,x,previous,current,next_value,error
    complex(dp) :: fs(6),ft(6,3,3),scalar(6),expected
    integer :: count,j,l,a,mo,mi
    count = 0
    do j = 0,7
      count = count+1
      call phase_wave(waves(count),0,j,[j],[central_phase(j)])
      count = count+1
      if (j == 0) then
        call phase_wave(waves(count),1,j,[1],[central_phase(1)])
      else
        call phase_wave(waves(count),1,j,[j-1,j+1],[central_phase(j-1),central_phase(j+1)])
        count = count+1
        call phase_wave(waves(count),1,j,[j],[central_phase(j)])
      end if
    end do
    theta = [0.0_dp,19.0_dp,67.0_dp,90.0_dp,143.0_dp,180.0_dp]
    q = 1.4_dp
    scalar = (0.0_dp,0.0_dp)
    do a = 1,size(theta)
      x = cos(theta(a)*pi/180.0_dp)
      previous = 1.0_dp
      current = x
      do l = 0,6
        if (l == 0) then
          next_value = 1.0_dp
        else if (l == 1) then
          next_value = x
        else
          next_value = (real(2*l-1,dp)*x*current-real(l-1,dp)*previous)/real(l,dp)
          previous = current
          current = next_value
        end if
        scalar(a) = scalar(a)+real(2*l+1,dp)*(exp(2.0_dp*iu*central_phase(l))-1.0_dp) &
          /(2.0_dp*iu*q)*next_value
      end do
    end do
    call spin_amplitudes(waves,q,theta,fs,ft,1.0_dp,'plain_momentum')
    error = maxval(abs(fs-scalar))
    do mi = 1,3
      do mo = 1,3
        do a = 1,size(theta)
          expected = (0.0_dp,0.0_dp)
          if (mo == mi) expected = scalar(a)
          error = max(error,abs(ft(a,mo,mi)-expected))
        end do
      end do
    end do
    call require(error < 2.0e-13_dp,'spin-independent central amplitude equals scalar Legendre series')
    write(*,'(a,es13.5)') 'Central Legendre-series maximum amplitude error: ',error
  end subroutine central_test

  subroutine coupled_integral_test()
    type(np_channel_result) :: waves(7), radial(7)
    real(dp) :: x(32),w(32),theta(32),angular,direct,optical,rephase_error
    complex(dp) :: fs(32),ft(32,3,3),fsr(32),ftr(32,3,3)
    integer :: b,a,i
    call phase_wave(waves(1),0,0,[0],[0.23_dp])
    call phase_wave(waves(2),0,1,[1],[-0.19_dp])
    call phase_wave(waves(3),1,0,[1],[0.47_dp])
    call phase_wave(waves(4),1,1,[1],[-0.21_dp])
    call phase_wave(waves(5),1,1,[0,2],[0.66_dp,-0.17_dp],0.29_dp)
    call phase_wave(waves(6),1,2,[1,3],[0.37_dp,0.11_dp],-0.18_dp)
    call phase_wave(waves(7),1,3,[2,4],[-0.22_dp,0.01_dp],0.1_dp)
    call gauss_legendre(32,-1.0_dp,1.0_dp,x,w)
    theta = acos(x)*180.0_dp/pi
    angular = 2.0_dp*pi*sum(w*differential_cross_section(waves,200.0_dp,theta,convention='plain_momentum'))
    direct = total_elastic_cross_section(waves,200.0_dp)
    optical = optical_total_cross_section(waves,200.0_dp,convention='plain_momentum')
    call require(abs(angular-direct) < 1.0e-11_dp,'coupled waves: angular integral equals direct elastic total')
    call require(abs(optical-direct) < 1.0e-11_dp,'coupled waves: optical theorem')
    radial = waves
    do i = 1,size(waves)
      do b = 1,waves(i)%nc
        do a = 1,waves(i)%nc
          radial(i)%s(a,b) = conjg(iu**waves(i)%ls(a))*waves(i)%s(a,b)*iu**waves(i)%ls(b)
        end do
      end do
    end do
    call spin_amplitudes(waves,200.0_dp,theta,fs,ft,convention='plain_momentum')
    call spin_amplitudes(radial,200.0_dp,theta,fsr,ftr)
    rephase_error = max(maxval(abs(fs-fsr)),maxval(abs(ft-ftr)))
    call require(rephase_error < 1.0e-13_dp,'orbital convention conversion preserves physical amplitudes')
    write(*,'(a,3es17.8)') 'Coupled angular/direct/optical totals (mb): ',angular,direct,optical
  end subroutine coupled_integral_test

  subroutine tensor_born_test()
    ! Independent Pauli-operator benchmark: T=(sigma1.Q)(sigma2.Q).
    ! Q=pout-pin, |pout|=|pin|=1, mass=1. These LSJ Born coefficients
    ! are in the EM500 radial basis. Their S-D signs matter at every angle.
    type(np_channel_result) :: waves(5)
    real(dp) :: tensor(2,2,5),theta(6),angle,sigq(2,2),direct_t(4,4)
    real(dp) :: singlet(4,1),triplet(4,3),direct_s(1,1),direct_triplet(3,3),error,wrong_error
    complex(dp) :: fs(6),ft(6,3,3),wrong_s(1),wrong_t(1,3,3)
    integer :: i,a,b,c,d,row,column
    call initialize_wave(waves(1),0,0,[0])
    call initialize_wave(waves(2),0,1,[1])
    call initialize_wave(waves(3),1,0,[1])
    call initialize_wave(waves(4),1,1,[1])
    call initialize_wave(waves(5),1,1,[0,2])
    tensor = 0.0_dp
    tensor(1,1,1) = -8.0_dp*pi
    tensor(1,1,2) = 8.0_dp*pi/3.0_dp
    tensor(1,1,3) = 8.0_dp*pi
    tensor(1,1,4) = -16.0_dp*pi/3.0_dp
    tensor(1,1,5) = 8.0_dp*pi/3.0_dp
    tensor(1,2,5) = -sqrt(2.0_dp)*8.0_dp*pi/3.0_dp
    tensor(2,1,5) = tensor(1,2,5)
    do i = 1,size(waves)
      waves(i)%s = waves(i)%s-iu*pi*tensor(:,:,i)
    end do
    theta = [0.0_dp,13.0_dp,45.0_dp,90.0_dp,143.0_dp,180.0_dp]
    call spin_amplitudes(waves,1.0_dp,theta,fs,ft,hbarc=1.0_dp)
    singlet(:,1) = [0.0_dp,1.0_dp/sqrt(2.0_dp),-1.0_dp/sqrt(2.0_dp),0.0_dp]
    triplet = 0.0_dp
    triplet(4,1) = 1.0_dp
    triplet(2:3,2) = 1.0_dp/sqrt(2.0_dp)
    triplet(1,3) = 1.0_dp
    error = 0.0_dp
    do i = 1,size(theta)
      angle = theta(i)*pi/180.0_dp
      sigq(1,:) = [cos(angle)-1.0_dp,sin(angle)]
      sigq(2,:) = [sin(angle),1.0_dp-cos(angle)]
      do b = 1,2
        do a = 1,2
          row = 2*(a-1)+b
          do d = 1,2
            do c = 1,2
              column = 2*(c-1)+d
              direct_t(row,column) = sigq(a,c)*sigq(b,d)
            end do
          end do
        end do
      end do
      direct_s = -2.0_dp*pi*pi*matmul(matmul(transpose(singlet),direct_t),singlet)
      direct_triplet = -2.0_dp*pi*pi*matmul(matmul(transpose(triplet),direct_t),triplet)
      error = max(error,abs(fs(i)-direct_s(1,1)),maxval(abs(ft(i,:,:)-direct_triplet)))
    end do
    call require(error < 2.0e-12_dp,'tensor Born amplitude equals direct Pauli operator in spin space')
    call require(max(abs(fs(1)),maxval(abs(ft(1,:,:)))) < 2.0e-13_dp,'tensor Born forward amplitude vanishes')
    call spin_amplitudes(waves,1.0_dp,[0.0_dp],wrong_s,wrong_t,1.0_dp,'plain_momentum')
    wrong_error = maxval(abs(wrong_t))
    call require(wrong_error > 100.0_dp,'wrong convention produces the expected large false forward amplitude')
    write(*,'(a,es13.5)') 'Pauli-tensor independent amplitude error: ',error
    write(*,'(a,es13.5)') 'Counterexample with omitted rephasing, false forward amplitude: ',wrong_error
  end subroutine tensor_born_test

  subroutine absorption_test()
    type(np_channel_result) :: wave(1)
    real(dp) :: elastic,total
    call initialize_wave(wave(1),0,0,[0])
    wave(1)%s = (0.0_dp,0.0_dp)
    elastic = total_elastic_cross_section(wave,1.0_dp,1.0_dp)
    total = optical_total_cross_section(wave,1.0_dp,1.0_dp,'plain_momentum')
    call require(abs(elastic-10.0_dp*pi/4.0_dp) < 1.0e-13_dp,'absorptive S=0 elastic normalization')
    call require(abs(total-2.0_dp*elastic) < 1.0e-13_dp,'absorptive S=0 total is twice elastic')
  end subroutine absorption_test

  subroutine zero_scattering_test()
    type(np_channel_result) :: waves(0),identity(1)
    real(dp) :: zero(3)
    zero = differential_cross_section(waves,1.0_dp,[0.0_dp,90.0_dp,180.0_dp])
    call require(maxval(abs(zero)) == 0.0_dp,'empty wave expansion has zero cross section')
    call initialize_wave(identity(1),1,1,[0,2])
    zero = differential_cross_section(identity,1.0_dp,[0.0_dp,90.0_dp,180.0_dp])
    call require(maxval(abs(zero)) == 0.0_dp,'identity S has zero cross section')
  end subroutine zero_scattering_test

  subroutine invalid_case(name)
    character(len=*), intent(in) :: name
    type(np_channel_result) :: wave(1),pair(2)
    real(dp) :: result(1)
    call phase_wave(wave(1),0,0,[0],[0.3_dp])
    select case(name)
    case('invalid_duplicate')
      pair(1) = wave(1)
      pair(2) = wave(1)
      result = differential_cross_section(pair,1.0_dp,[30.0_dp])
    case('invalid_coupling')
      wave(1)%spin = 1
      result = differential_cross_section(wave,1.0_dp,[30.0_dp])
    case('invalid_angle')
      result = differential_cross_section(wave,1.0_dp,[181.0_dp])
    case('invalid_q')
      result = differential_cross_section(wave,0.0_dp,[30.0_dp])
    case('invalid_convention')
      result = differential_cross_section(wave,1.0_dp,[30.0_dp],convention='unspecified')
    case('invalid_s')
      wave(1)%s(1,1) = cmplx(ieee_value(1.0_dp,ieee_quiet_nan),0.0_dp,dp)
      result = differential_cross_section(wave,1.0_dp,[30.0_dp])
    case default
      error stop 'Unknown invalid-input test case'
    end select
    write(*,*) result
  end subroutine invalid_case
end program check_np_observables
