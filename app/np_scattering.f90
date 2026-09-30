! Real principal-value np Lippmann--Schwinger solver, NN only.
! Momenta/masses in MeV, V/K/T in MeV^-2; radial measure dk*k^2 (NO 2/pi).
! This is the Fortran migration of fewbody_research/scattering.py.
module fb_np_scattering
  use fb_precision, only: dp
  use fb_quadrature, only: gauss_legendre
  use em500_force, only: initialize_em500, np_kernel, orbital_channels, mass_parameter_MeV
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: np_channel_result, solve_np_channel, solve_k_matrix, neutron_lab_momentum
  ! Same masses as the frozen EM500 configuration, not new fitted parameters.
  real(dp), parameter :: mp_MeV=938.272_dp, mn_MeV=939.5653_dp
  type :: np_channel_result
    integer :: spin=-1,j=-1,isospin=-1,nc=0,ls(2)=-1
    real(dp) :: q=0.0_dp,k(2,2)=0.0_dp
    complex(dp) :: t(2,2)=(0.0_dp,0.0_dp),s(2,2)=(0.0_dp,0.0_dp)
    real(dp) :: phase_deg(2)=0.0_dp,mixing_deg=0.0_dp
    real(dp) :: residual_relative=0.0_dp,reciprocity_error_relative=0.0_dp,unitarity_error=0.0_dp
  end type np_channel_result
  interface
    subroutine dgesv(n,nrhs,a,lda,ipiv,b,ldb,info)
      import dp
      integer, intent(in) :: n,nrhs,lda,ldb
      integer, intent(out) :: ipiv(*),info
      real(dp), intent(inout) :: a(lda,*),b(ldb,*)
    end subroutine dgesv
    subroutine dsyev(jobz,uplo,n,a,lda,w,work,lwork,info)
      import dp
      character, intent(in) :: jobz,uplo
      integer, intent(in) :: n,lda,lwork
      integer, intent(out) :: info
      real(dp), intent(inout) :: a(lda,*)
      real(dp), intent(out) :: w(*),work(*)
    end subroutine dsyev
  end interface
contains
  function neutron_lab_momentum(energy,kinematics) result(q)
    real(dp), intent(in) :: energy
    character(len=*), intent(in), optional :: kinematics
    real(dp) :: q,q2
    character(len=32) :: mode
    if (.not. ieee_is_finite(energy) .or. energy <= 0.0_dp) &
      error stop 'Require positive finite incident neutron laboratory energy'
    mode='nonrelativistic'
    if (present(kinematics)) mode=kinematics
    select case(trim(mode))
    case('nonrelativistic')
      q2=mass_parameter_MeV*mp_MeV/(mp_MeV+mn_MeV)*energy
    case('relativistic')
      ! Exact two-body q only: the LS spectral parameter remains q^2/(2mu).
      q2=mp_MeV**2*energy*(energy+2.0_dp*mn_MeV)/((mp_MeV+mn_MeV)**2+2.0_dp*mp_MeV*energy)
    case default
      error stop 'kinematics must be nonrelativistic or relativistic'
    end select
    q=sqrt(q2)
    if (.not. ieee_is_finite(q)) error stop 'Nonfinite incident momentum'
  end function neutron_lab_momentum

  subroutine solve_np_channel(spin,j,isospin,coupled,n_grid,pmax,q,result)
    integer, intent(in) :: spin,j,isospin,n_grid
    logical, intent(in) :: coupled
    real(dp), intent(in) :: pmax,q
    type(np_channel_result), intent(out) :: result
    real(dp), allocatable :: nodes(:),weights(:),grid(:),v(:,:)
    real(dp) :: block(2,2)
    integer :: nc,ls(2),n1,i,h,a,b
    call orbital_channels(spin,j,isospin,coupled,nc,ls)
    if (n_grid < 4 .or. n_grid > 1000) error stop 'Scattering requires n_grid in [4,1000]'
    if (.not. ieee_is_finite(pmax) .or. pmax <= 0.0_dp .or. pmax > 1000.0_dp) &
      error stop 'EM500 momentum cutoff must be finite in (0,1000] MeV'
    if (.not. ieee_is_finite(q) .or. q <= 0.0_dp .or. q >= pmax) &
      error stop 'Scattering requires finite 0 < q < pmax'
    n1=n_grid+1
    allocate(nodes(n_grid),weights(n_grid),grid(n1),v(nc*n1,nc*n1))
    call gauss_legendre(n_grid,0.0_dp,pmax,nodes,weights)
    grid(1:n_grid)=nodes
    grid(n1)=q
    call initialize_em500()
    ! Evaluate both momentum orders independently: do not assume reciprocity
    ! or symmetrize V before checking the upstream mapping and numerical values.
    do h=1,n1
      do i=1,n1
        call np_kernel(grid(i),grid(h),spin,j,isospin,coupled,block)
        do b=1,nc
          do a=1,nc
            v((a-1)*n1+i,(b-1)*n1+h)=block(a,b)
          end do
        end do
      end do
    end do
    call solve_k_matrix(nodes,weights,v,q,mass_parameter_MeV,pmax,result)
    result%spin=spin; result%j=j; result%isospin=isospin; result%ls=ls
  end subroutine solve_np_channel

  subroutine solve_k_matrix(nodes,weights,v,q,mass,cutoff,result)
    ! Each channel is ordered [k_1,...,k_N,q]. Principal-value subtraction:
    ! J = atanh(q/cutoff)/q
    ! d_i = m*w_i*k_i^2/(q^2-k_i^2)
    ! d_q = m*q^2*(J - sum_i w_i/(q^2-k_i^2))
    ! Solve (I-V*diag(d))*K(:,q)=V(:,q), including the finite-cutoff J.
    real(dp), intent(in) :: nodes(:),weights(:),v(:,:),q,mass,cutoff
    type(np_channel_result), intent(out) :: result
    real(dp), parameter :: reciprocity_tolerance=1.0e-8_dp
    real(dp), allocatable :: denom(:),d(:),a(:,:),original_a(:,:),rhs(:,:),solution(:,:)
    integer, allocatable :: pivots(:)
    real(dp) :: jpv,raw_k(2,2),scale,rho,rotation(2,2),eigenvalues(2),work(8),tempcol(2),temp
    real(dp) :: identity(2,2),pi
    complex(dp) :: spectral_t,closure(2,2)
    integer :: n,n1,nc,dim,i,h,c,b,info,onshell(2)
    result=np_channel_result()
    n=size(nodes); n1=n+1
    if (n < 1 .or. size(weights) /= n) error stop 'PV quadrature requires matching nonempty vectors'
    if (.not. all(ieee_is_finite(nodes)) .or. .not. all(ieee_is_finite(weights))) &
      error stop 'Nonfinite PV quadrature'
    if (.not. ieee_is_finite(q) .or. .not. ieee_is_finite(mass) .or. .not. ieee_is_finite(cutoff)) &
      error stop 'Nonfinite PV momentum, mass, or cutoff'
    if (mass <= 0.0_dp .or. q <= 0.0_dp .or. q >= cutoff) error stop 'PV requires mass>0 and 0<q<cutoff'
    if (any(weights <= 0.0_dp) .or. any(nodes <= 0.0_dp) .or. any(nodes >= cutoff)) &
      error stop 'PV requires positive weights and nodes inside (0,cutoff)'
    if (any(nodes(2:n) <= nodes(1:n-1))) error stop 'PV nodes must increase strictly'
    if (any(abs(nodes-q) <= 1.0e-11_dp*max(nodes,q))) &
      error stop 'On-shell momentum coincides with a quadrature node; change n_grid'
    dim=size(v,1)
    if (size(v,2) /= dim .or. (dim /= n1 .and. dim /= 2*n1)) &
      error stop 'PV potential must have one or two augmented channel blocks'
    if (.not. all(ieee_is_finite(v))) error stop 'Nonfinite PV potential'
    nc=dim/n1
    scale=max(maxval(abs(v)),tiny(1.0_dp))
    if (maxval(abs(v-transpose(v)))/scale > reciprocity_tolerance) &
      error stop 'Potential violates reciprocity V_ab(p,k)=V_ba(k,p)'
    allocate(denom(n),d(dim),a(dim,dim),original_a(dim,dim),rhs(dim,nc),solution(dim,nc),pivots(dim))
    denom=q*q-nodes*nodes
    jpv=atanh(q/cutoff)/q
    do c=1,nc
      d((c-1)*n1+1:c*n1-1)=mass*weights*nodes*nodes/denom
      d(c*n1)=mass*q*q*(jpv-sum(weights/denom))
      onshell(c)=c*n1
      rhs(:,c)=v(:,onshell(c))
    end do
    do h=1,dim
      a(:,h)=-v(:,h)*d(h)
      a(h,h)=a(h,h)+1.0_dp
    end do
    original_a=a
    solution=rhs
    call dgesv(dim,nc,a,dim,pivots,solution,dim,info)
    if (info < 0) error stop 'DGESV invalid argument in PV solver'
    if (info > 0) error stop 'PV K system is singular (possible K-matrix pole)'
    if (.not. all(ieee_is_finite(solution))) error stop 'Nonfinite PV K solution'
    result%residual_relative=sqrt(sum((matmul(original_a,solution)-rhs)**2))/max(sqrt(sum(rhs**2)),tiny(1.0_dp))
    raw_k=0.0_dp
    do c=1,nc
      raw_k(c,1:nc)=solution(onshell(c),:)
    end do
    scale=max(maxval(abs(raw_k(1:nc,1:nc))),tiny(1.0_dp))
    result%reciprocity_error_relative=maxval(abs(raw_k-transpose(raw_k)))/scale
    if (result%reciprocity_error_relative > reciprocity_tolerance) &
      error stop 'On-shell K fails reciprocity; refine mesh/check conditioning'
    ! Remove roundoff antisymmetry only after the diagnostic has passed.
    result%k=0.5_dp*(raw_k+transpose(raw_k))
    pi=acos(-1.0_dp)
    rho=pi*mass*q/2.0_dp
    identity=0.0_dp
    do c=1,nc
      identity(c,c)=1.0_dp
    end do
    rotation=identity
    eigenvalues=0.0_dp
    if (nc == 1) then
      eigenvalues(1)=result%k(1,1)
    else
      rotation=result%k
      call dsyev('V','U',2,rotation,2,eigenvalues,work,size(work),info)
      if (info /= 0) error stop 'DSYEV failed for on-shell K'
      ! Same convention as Python: first eigenvector dominated by lower L,
      ! and R=[[cos(eps),-sin(eps)],[sin(eps),cos(eps)]]. No energy tracking.
      if (abs(rotation(1,2)) > abs(rotation(1,1))) then
        tempcol=rotation(:,1); rotation(:,1)=rotation(:,2); rotation(:,2)=tempcol
        temp=eigenvalues(1); eigenvalues(1)=eigenvalues(2); eigenvalues(2)=temp
      end if
      if (rotation(1,1) < 0.0_dp) rotation(:,1)=-rotation(:,1)
      if (rotation(1,1)*rotation(2,2)-rotation(1,2)*rotation(2,1) < 0.0_dp) rotation(:,2)=-rotation(:,2)
      result%mixing_deg=atan2(rotation(2,1),rotation(1,1))*180.0_dp/pi
    end if
    ! T=K/(I+i*rho*K), S=I-2i*rho*T. A real orthogonal K spectrum
    ! avoids introducing a separate complex linear-algebra dependency.
    result%t=(0.0_dp,0.0_dp)
    do c=1,nc
      spectral_t=eigenvalues(c)/cmplx(1.0_dp,rho*eigenvalues(c),kind=dp)
      do b=1,nc
        do i=1,nc
          result%t(i,b)=result%t(i,b)+rotation(i,c)*spectral_t*rotation(b,c)
        end do
      end do
      result%phase_deg(c)=atan(-rho*eigenvalues(c))*180.0_dp/pi
    end do
    result%s=cmplx(identity,0.0_dp,kind=dp)-cmplx(0.0_dp,2.0_dp*rho,kind=dp)*result%t
    closure=matmul(conjg(transpose(result%s)),result%s)-cmplx(identity,0.0_dp,kind=dp)
    result%unitarity_error=maxval(abs(closure(1:nc,1:nc)))
    result%nc=nc; result%q=q
  end subroutine solve_k_matrix
end module fb_np_scattering
