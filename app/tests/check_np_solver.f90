program check_np_solver
  use fb_precision, only: dp
  use fb_quadrature, only: gauss_legendre
  use fb_np_scattering, only: np_channel_result, solve_k_matrix, neutron_lab_momentum, solve_np_channel
  use em500_force, only: initialize_em500, np_kernel, sd_kernel, orbital_channels, mass_parameter_MeV
  use idaho_chiral_potential, only: chp
  implicit none
  real(dp) :: nodes(70),weights(70),grid(71),g(71),v(71,71),q,beta,mass,strength,integral,exact,rho,error
  real(dp) :: qs(4),nodes24(24),weights24(24),v25(25,25),pi
  real(dp) :: rotation(2,2),coupling(2,2),coupled_v(142,142),eigen(2),expected(2,2),angle
  real(dp) :: block(2,2),reverse(2,2),raw(6),old_sd(2,2),mp,mn,sinv,kallen,err_kin
  complex(dp) :: exact_s
  integer :: i,j,a,b,nc,ls(2),jj,spin,iso
  logical :: coupled
  character(len=40) :: mode
  type(np_channel_result) :: r
  pi=acos(-1.0_dp)
  call get_command_argument(1,mode)
  if (len_trim(mode) > 0) then
    select case(trim(mode))
    case('reject_pauli')
      call orbital_channels(0,0,0,.false.,nc,ls)
    case('reject_j0_triplet_uncoupled')
      call orbital_channels(1,0,1,.false.,nc,ls)
    case('reject_pole_node')
      call gauss_legendre(70,0.0_dp,6.0_dp,nodes,weights)
      v=0.0_dp
      call solve_k_matrix(nodes,weights,v,nodes(3),1.0_dp,6.0_dp,r)
    case('reject_nonreciprocal')
      call gauss_legendre(70,0.0_dp,6.0_dp,nodes,weights)
      v=0.0_dp; v(1,2)=0.1_dp
      call solve_k_matrix(nodes,weights,v,0.2_dp,1.0_dp,6.0_dp,r)
    case default
      error stop 'Unknown test mode'
    end select
    error stop 'Expected input rejection did not happen'
  end if
  call gauss_legendre(70,0.0_dp,6.0_dp,nodes,weights)
  qs=[0.2_dp,0.9_dp,2.0_dp,5.5_dp]
  mass=3.0_dp; beta=1.0_dp; strength=0.7_dp
  do i=1,size(qs)
    q=qs(i)
    grid(1:70)=nodes; grid(71)=q
    g=1.0_dp/(grid**2+beta**2)
    do j=1,71
      v(:,j)=-strength*g*g(j)
    end do
    call solve_k_matrix(nodes,weights,v,q,mass,6.0_dp,r)
    integral=yamaguchi_integral(q,beta,6.0_dp,mass)
    exact=-strength/(q*q+beta*beta)**2/(1.0_dp+strength*integral)
    error=abs(r%k(1,1)-exact)/max(1.0_dp,abs(exact))
    rho=pi*mass*q/2.0_dp
    exact_s=cmplx(1.0_dp,-rho*exact,dp)/cmplx(1.0_dp,rho*exact,dp)
    if (error > 2.0e-10_dp .or. abs(r%s(1,1)-exact_s) > 1.0e-10_dp) error stop 'Yamaguchi analytic mismatch'
    if (r%residual_relative > 1.0e-12_dp .or. r%unitarity_error > 1.0e-12_dp) error stop 'Yamaguchi diagnostics'
    write(*,'(A,F6.2,3ES16.7)') 'Yamaguchi q,K-error,residual,unitarity: ',q,error,r%residual_relative,r%unitarity_error
  end do
  call gauss_legendre(24,0.0_dp,7.0_dp,nodes24,weights24)
  q=4.1_dp; mass=2.0_dp; strength=-0.04_dp
  v25=strength
  call solve_k_matrix(nodes24,weights24,v25,q,mass,7.0_dp,r)
  integral=mass*(-7.0_dp+q/2.0_dp*log((7.0_dp+q)/(7.0_dp-q)))
  exact=strength/(1.0_dp-strength*integral)
  if (abs(r%k(1,1)-exact) > 1.0e-12_dp) error stop 'Constant potential finite cutoff logarithm'
  write(*,'(A,ES16.7)') 'Constant potential exact error: ',abs(r%k(1,1)-exact)
  q=0.7_dp; beta=1.1_dp; mass=2.0_dp; angle=0.23_dp
  rotation(1,:)=[cos(angle),-sin(angle)]; rotation(2,:)=[sin(angle),cos(angle)]
  do b=1,2
    do a=1,2
      coupling(a,b)=0.2_dp*rotation(a,1)*rotation(b,1)+0.9_dp*rotation(a,2)*rotation(b,2)
    end do
  end do
  grid(1:70)=nodes; grid(71)=q
  g=1.0_dp/(grid**2+beta**2)
  do b=1,2
    do a=1,2
      do j=1,71
        coupled_v((a-1)*71+1:a*71,(b-1)*71+j)=-coupling(a,b)*g*g(j)
      end do
    end do
  end do
  call solve_k_matrix(nodes,weights,coupled_v,q,mass,6.0_dp,r)
  integral=yamaguchi_integral(q,beta,6.0_dp,mass)
  eigen=-[0.2_dp,0.9_dp]/(1.0_dp+[0.2_dp,0.9_dp]*integral)/(q*q+beta*beta)**2
  do b=1,2
    do a=1,2
      expected(a,b)=rotation(a,1)*eigen(1)*rotation(b,1)+rotation(a,2)*eigen(2)*rotation(b,2)
    end do
  end do
  error=maxval(abs(r%k-expected))
  if (error > 1.0e-11_dp) error stop 'Rotated coupled Yamaguchi mismatch'
  if (maxval(abs(r%phase_deg-atan(-pi*mass*q/2.0_dp*eigen)*180.0_dp/pi)) > 1.0e-9_dp) &
    error stop 'Rotated coupled eigenphase mismatch'
  if (abs(r%mixing_deg-angle*180.0_dp/pi) > 1.0e-9_dp) error stop 'Mixing convention mismatch'
  write(*,'(A,2ES16.7)') 'Coupled analytic K and mixing errors: ',error,abs(r%mixing_deg-angle*180.0_dp/pi)
  ! Verify every allowed channel family, including the upstream J=0 triplet exception.
  call initialize_em500()
  do jj=0,8
    do i=1,3
      if (jj == 0 .and. i == 2) cycle
      spin=merge(0,1,i==1)
      coupled=i==3
      if (coupled) then
        iso=modulo(jj+1,2)
      else
        iso=modulo(jj+spin+1,2)
      end if
      call orbital_channels(spin,jj,iso,coupled,nc,ls)
      call np_kernel(73.0_dp,241.0_dp,spin,jj,iso,coupled,block)
      call np_kernel(241.0_dp,73.0_dp,spin,jj,iso,coupled,reverse)
      call chp(73.0_dp,241.0_dp,coupled,spin,jj,iso,0,raw)
      expected=0.0_dp
      if (coupled .and. jj > 0) then
        expected(1,1)=raw(4); expected(1,2)=raw(6)
        expected(2,1)=raw(5); expected(2,2)=raw(3)
      else if (coupled) then
        expected(1,1)=raw(3)
      else
        expected(1,1)=raw(spin+1)
      end if
      if (any(block /= expected)) error stop 'Raw EM500 mapping mismatch'
      error=maxval(abs(block-transpose(reverse)))/max(maxval(abs(block)),tiny(1.0_dp))
      if (error > 1.0e-8_dp) error stop 'EM500 momentum reciprocity mismatch'
    end do
  end do
  call sd_kernel(73.0_dp,241.0_dp,old_sd)
  call np_kernel(73.0_dp,241.0_dp,1,1,0,.true.,block)
  if (any(old_sd /= block)) error stop 'Deuteron wrapper changed'
  write(*,'(A)') 'EM500 mapping/reciprocity: 26 allowed J<=8 channel families passed'
  mp=938.272_dp; mn=939.5653_dp
  q=neutron_lab_momentum(100.0_dp,'relativistic')
  sinv=(mp+mn)**2+2.0_dp*mp*100.0_dp
  kallen=(sinv-(mp+mn)**2)*(sinv-(mp-mn)**2)
  err_kin=abs(q*q-kallen/(4.0_dp*sinv))
  if (err_kin > 1.0e-8_dp) error stop 'Invariant two-body kinematics mismatch'
  q=neutron_lab_momentum(100.0_dp)
  if (abs(q*q-mass_parameter_MeV*mp/(mp+mn)*100.0_dp) > 1.0e-8_dp) error stop 'Nonrelativistic kinematics mismatch'
  write(*,'(A,ES16.7)') 'Invariant two-body q^2 error (MeV^2): ',err_kin
  call solve_np_channel(1,0,1,.true.,32,1000.0_dp,q,r)
  if (r%nc /= 1 .or. r%ls(1) /= 1) error stop '3P0 solved channel mapping'
  if (r%residual_relative > 1.0e-10_dp .or. r%unitarity_error > 1.0e-11_dp) error stop '3P0 diagnostics'
  call solve_np_channel(1,2,1,.true.,32,1000.0_dp,q,r)
  if (r%nc /= 2 .or. any(r%ls /= [1,3])) error stop '3P2-F2 solved channel mapping'
  if (r%residual_relative > 1.0e-10_dp .or. r%unitarity_error > 1.0e-11_dp) error stop '3P2-F2 diagnostics'
  write(*,'(A)') 'PASS: analytic PV, finite cutoff, coupled phases, raw mapping, kinematics, EM500 channels'
contains
  function yamaguchi_integral(q,beta,cutoff,mass) result(value)
    real(dp), intent(in) :: q,beta,cutoff,mass
    real(dp) :: value,a,c,jpv,i1,i2
    a=q*q/(q*q+beta*beta)**2
    c=-beta*beta/(q*q+beta*beta)
    jpv=log((cutoff+q)/(cutoff-q))/(2.0_dp*q)
    i1=atan(cutoff/beta)/beta
    i2=cutoff/(2.0_dp*beta*beta*(cutoff*cutoff+beta*beta))+atan(cutoff/beta)/(2.0_dp*beta**3)
    value=mass*(a*jpv+a*i1+c*i2)
  end function yamaguchi_integral
end program check_np_solver
