program check_pp_solver
  use fb_precision, only: dp
  use fb_quadrature, only: gauss_legendre
  use fb_np_scattering, only: np_channel_result,solve_k_matrix
  use em500_force, only: initialize_em500,pp_kernel,orbital_channels
  use idaho_chiral_potential, only: chp
  use fb_pp_scattering, only: solve_pp_channel,proton_lab_momentum,coulomb_functions,spherical_bessel_j
  implicit none
  type(np_channel_result) :: result,reference,almost_neutral
  real(dp) :: q,x,eta,raw(6),block(2,2),expected(2,2),reverse(2,2),error,maximum
  real(dp) :: f(0:32),g(0:32),fp(0:32),gp(0:32),nodes(48),weights(48),grid(49),v(98,98)
  real(dp),parameter :: xs(3)=[0.3_dp,3.0_dp,30.0_dp]
  real(dp),parameter :: energies(3)=[3.037_dp,98.0_dp,286.7_dp]
  integer :: i,j,a,b,c,jj,spin,nc,ls(2),ie
  logical :: coupled
  ! Wrapper = raw upstream pp(Tz=-1), with no use of np charge parameters.
  call initialize_em500()
  do jj=0,8
    do c=1,3
      spin=merge(0,1,c==1); coupled=c==3
      if(c==1 .and. modulo(jj,2)/=0) cycle
      if(c==2 .and. (jj==0 .or. modulo(jj,2)/=1)) cycle
      if(c==3 .and. modulo(jj,2)/=0) cycle
      call pp_kernel(73.0_dp,241.0_dp,spin,jj,coupled,block)
      call pp_kernel(241.0_dp,73.0_dp,spin,jj,coupled,reverse)
      call chp(73.0_dp,241.0_dp,coupled,spin,jj,1,-1,raw)
      expected=0.0_dp
      if(coupled .and. jj>0) then
        expected(1,1)=raw(4); expected(1,2)=raw(6); expected(2,1)=raw(5); expected(2,2)=raw(3)
      else if(coupled) then
        expected(1,1)=raw(3)
      else
        expected(1,1)=raw(spin+1)
      end if
      if(any(block/=expected)) error stop 'pp wrapper raw Tz mapping mismatch'
      error=maxval(abs(block-transpose(reverse)))/max(maxval(abs(block)),tiny(1.0_dp))
      if(error>1.0e-8_dp) error stop 'pp potential reciprocity mismatch'
    end do
  end do
  write(*,'(A)') 'pp wrapper: charge-dependent T=1,Tz=-1 mapping/reciprocity passed'
  ! Free Riccati F0=x*j_l versus independently implemented Miller recurrence.
  do i=1,3
    x=xs(i)
    call coulomb_functions(x,0.0_dp,32,f,g,fp,gp)
    do j=0,32
      error=abs(f(j)-x*spherical_bessel_j(j,x))/max(abs(f(j)),1.0e-100_dp)
      if(error>2.0e-11_dp) error stop 'Free F versus spherical Bessel mismatch'
    end do
  end do
  write(*,'(A)') 'Free Coulomb functions versus Miller spherical Bessel j_L: passed'
  ! Pure point Coulomb: matched nuclear S_C should be I. Finite pmax=1000
  ! leaves a measured ~1.2e-5 maximum S residual; this is NOT exact zero.
  maximum=0.0_dp
  do ie=1,3
    q=proton_lab_momentum(energies(ie))
    do jj=0,4,2
      call solve_pp_channel(0,jj,.false.,160,1000.0_dp,q,16.0_dp,192,result,include_nuclear=.false.)
      error=abs(result%s(1,1)-1.0_dp)
      maximum=max(maximum,error)
      if(error>2.0e-5_dp) error stop 'Pure Coulomb VP residual exceeds finite-cutoff tolerance'
      if(result%unitarity_error>1.0e-12_dp) error stop 'Pure Coulomb unitarity failed'
    end do
  end do
  write(*,'(A,ES16.7)') 'Pure Coulomb max |S_C-I|, pmax=1000 MeV (tail error): ',maximum
  ! At alpha=0 the pp result must equal the independent existing PV solver
  ! driven directly with pp force and pp reduced mass, including coupling.
  q=proton_lab_momentum(98.0_dp)
  call gauss_legendre(48,0.0_dp,1000.0_dp,nodes,weights)
  grid(1:48)=nodes; grid(49)=q
  do j=1,49
    do i=1,49
      call pp_kernel(grid(i),grid(j),1,2,.true.,block)
      do b=1,2
        do a=1,2
          v((a-1)*49+i,(b-1)*49+j)=block(a,b)
        end do
      end do
    end do
  end do
  call solve_k_matrix(nodes,weights,v,q,938.272_dp,1000.0_dp,reference)
  call solve_pp_channel(1,2,.true.,48,1000.0_dp,q,16.0_dp,192,result,alpha_value=0.0_dp)
  error=maxval(abs(reference%s-result%s))
  if(error>1.0e-12_dp) error stop 'Neutral pp direct-LS bridge failed'
  ! Small nonzero alpha enters and tests the actual matching path, unlike
  ! alpha=0 which deliberately bypasses Coulomb matching.
  call solve_pp_channel(1,2,.true.,48,1000.0_dp,q,16.0_dp,192,almost_neutral,alpha_value=1.0e-11_dp)
  if(maxval(abs(almost_neutral%s-result%s))>1.0e-8_dp) error stop 'VP alpha->0 coupled limit failed'
  write(*,'(A,2ES16.7)') 'Alpha=0 bridge and alpha->0 matching errors: ',error,maxval(abs(almost_neutral%s-result%s))
  if(almost_neutral%unitarity_error>1.0e-12_dp) error stop 'Coupled Coulomb matching unitarity failed'
  write(*,'(A)') 'PASS: pp Tz, Pauli channels, Bessel normalization, pure Coulomb, neutral bridge, coupled matching'
end program check_pp_solver
