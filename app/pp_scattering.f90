! pp scattering: real momentum-space LS + Vincent--Phatak Coulomb matching.
! V is in MeV^-2, radial measure dp*p^2; all range inputs are fm.
module fb_pp_scattering
  use fb_precision, only: dp
  use em500_force, only: initialize_em500, pp_kernel, orbital_channels, hbarc_MeV_fm
  use fb_np_scattering, only: np_channel_result, solve_k_matrix
  use fb_quadrature, only: gauss_legendre
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: solve_pp_channel, proton_lab_momentum, pp_coulomb_eta, pp_alpha_effective
  public :: coulomb_functions, spherical_bessel_j
  real(dp), parameter :: mp_MeV=938.272_dp, alpha_default=1.0_dp/137.0359895_dp
  real(dp), parameter :: pi=acos(-1.0_dp)
  interface
    subroutine coul90(x,eta,xlmin,lrange,fc,gc,fcp,gcp,kfn,ifail)
      import dp
      real(dp) :: x,eta,xlmin,fc(0:*),gc(0:*),fcp(0:*),gcp(0:*)
      integer :: lrange,kfn,ifail
    end subroutine coul90
    subroutine dgesv(n,nrhs,a,lda,ipiv,b,ldb,info)
      import dp
      integer,intent(in) :: n,nrhs,lda,ldb
      integer,intent(out) :: ipiv(*),info
      real(dp),intent(inout) :: a(lda,*),b(ldb,*)
    end subroutine dgesv
    subroutine dsyev(jobz,uplo,n,a,lda,w,work,lwork,info)
      import dp
      character,intent(in) :: jobz,uplo
      integer,intent(in) :: n,lda,lwork
      integer,intent(out) :: info
      real(dp),intent(inout) :: a(lda,*)
      real(dp),intent(out) :: w(*),work(*)
    end subroutine dsyev
  end interface
contains
  function proton_lab_momentum(energy) result(q)
    real(dp),intent(in) :: energy
    real(dp) :: q
    if (.not.ieee_is_finite(energy) .or. energy<=0.0_dp) error stop 'pp requires positive laboratory energy'
    ! Equal masses: exact relativistic and nonrelativistic CM momenta coincide.
    q=sqrt(mp_MeV*energy/2.0_dp)
  end function proton_lab_momentum

  function pp_alpha_effective(q,alpha_value,relativistic_coulomb) result(alpha)
    real(dp),intent(in) :: q
    real(dp),intent(in),optional :: alpha_value
    logical,intent(in),optional :: relativistic_coulomb
    real(dp) :: alpha
    logical :: relativistic
    alpha=alpha_default
    if (present(alpha_value)) alpha=alpha_value
    if (.not.ieee_is_finite(alpha) .or. alpha<0.0_dp) error stop 'pp alpha must be finite and nonnegative'
    relativistic=.true.
    if (present(relativistic_coulomb)) relativistic=relativistic_coulomb
    if (relativistic) alpha=alpha*(mp_MeV**2+2.0_dp*q*q)/(mp_MeV*sqrt(mp_MeV**2+q*q))
  end function pp_alpha_effective

  function pp_coulomb_eta(q,alpha_value,relativistic_coulomb) result(eta)
    real(dp),intent(in) :: q
    real(dp),intent(in),optional :: alpha_value
    logical,intent(in),optional :: relativistic_coulomb
    real(dp) :: eta
    if (.not.ieee_is_finite(q) .or. q<=0.0_dp) error stop 'pp eta requires positive finite q'
    eta=mp_MeV*pp_alpha_effective(q,alpha_value,relativistic_coulomb)/(2.0_dp*q)
  end function pp_coulomb_eta

  subroutine coulomb_functions(x,eta,lmax,f,g,fp,gp)
    real(dp),intent(in) :: x,eta
    integer,intent(in) :: lmax
    real(dp),intent(out) :: f(0:lmax),g(0:lmax),fp(0:lmax),gp(0:lmax)
    real(dp) :: xx,ee,lmin
    integer :: ifail,kfn,nl
    if (x<=1.0e-7_dp .or. lmax<0) error stop 'Invalid Coulomb function arguments'
    xx=x; ee=eta; lmin=0.0_dp; nl=lmax; kfn=0; ifail=0
    call coul90(xx,ee,lmin,nl,f,g,fp,gp,kfn,ifail)
    if (ifail/=0) error stop 'Barnett COUL90 failed'
    if (.not.all(ieee_is_finite(f)) .or. .not.all(ieee_is_finite(g))) error stop 'Nonfinite Coulomb function'
    if (maxval(abs(fp*g-f*gp-1.0_dp))>1.0e-7_dp) error stop 'Coulomb Wronskian check failed'
  end subroutine coulomb_functions

  pure function spherical_bessel_j(l,x) result(value)
    integer,intent(in) :: l
    real(dp),intent(in) :: x
    real(dp) :: value,j0,j1,a,b,c,kept,scale
    integer :: m,n
    if (abs(x)<1.0e-7_dp) then
      value=1.0_dp
      do n=1,l
        value=value*x/real(2*n+1,dp)
      end do
      value=value*(1.0_dp-x*x/(2.0_dp*real(2*l+3,dp)))
      return
    end if
    j0=sin(x)/x
    j1=(j0-cos(x))/x
    if (l==0) then
      value=j0; return
    else if (l==1) then
      value=j1; return
    end if
    if (x>real(l+10,dp)) then
      a=j0; b=j1
      do n=1,l-1
        c=real(2*n+1,dp)*b/x-a; a=b; b=c
      end do
      value=b
    else
      m=max(l+30,int(x)+30)
      a=0.0_dp; b=1.0_dp; kept=0.0_dp
      do n=m,1,-1
        c=real(2*n+1,dp)*b/x-a
        a=b; b=c
        if (n-1==l) kept=b
        if (abs(b)>1.0e100_dp) then
          a=a*1.0e-100_dp; b=b*1.0e-100_dp; kept=kept*1.0e-100_dp
        end if
      end do
      if (abs(j0)>=abs(j1)) then
        scale=j0/b
      else
        scale=j1/a
      end if
      value=kept*scale
    end if
  end function spherical_bessel_j

  subroutine solve_pp_channel(spin,j,coupled,n_grid,pmax,q,radius_fm,n_radial,result, &
      include_nuclear,alpha_value,relativistic_coulomb)
    integer,intent(in) :: spin,j,n_grid,n_radial
    logical,intent(in) :: coupled
    real(dp),intent(in) :: pmax,q,radius_fm
    type(np_channel_result),intent(out) :: result
    logical,intent(in),optional :: include_nuclear,relativistic_coulomb
    real(dp),intent(in),optional :: alpha_value
    real(dp),allocatable :: nodes(:),weights(:),grid(:),v(:,:),r(:),wr(:),basis(:,:),weighted(:,:),vc(:,:)
    real(dp) :: block(2,2),alpha,eta
    logical :: nuclear
    integer :: nc,ls(2),n1,i,h,a,b,t,l
    call orbital_channels(spin,j,1,coupled,nc,ls)
    if (n_grid<4 .or. n_grid>1000 .or. n_radial<8 .or. n_radial>4000) error stop 'Invalid pp quadrature sizes'
    if (pmax<=q .or. pmax>1000.0_dp .or. q<=0.0_dp) error stop 'Require 0<q<pmax<=1000 MeV'
    if (radius_fm<1.0_dp .or. radius_fm>100.0_dp) error stop 'Require 1<=pp matching radius<=100 fm'
    n1=n_grid+1
    allocate(nodes(n_grid),weights(n_grid),grid(n1),v(nc*n1,nc*n1))
    call gauss_legendre(n_grid,0.0_dp,pmax,nodes,weights)
    grid(1:n_grid)=nodes; grid(n1)=q
    v=0.0_dp
    nuclear=.true.
    if(present(include_nuclear)) nuclear=include_nuclear
    if(nuclear) then
      call initialize_em500()
      do h=1,n1
        do i=1,n1
          call pp_kernel(grid(i),grid(h),spin,j,coupled,block)
          do b=1,nc
            do a=1,nc
              v((a-1)*n1+i,(b-1)*n1+h)=block(a,b)
            end do
          end do
        end do
      end do
    end if
    alpha=pp_alpha_effective(q,alpha_value,relativistic_coulomb)
    if(alpha>0.0_dp) then
      allocate(r(n_radial),wr(n_radial),basis(n1,n_radial),weighted(n1,n_radial),vc(n1,n1))
      call gauss_legendre(n_radial,0.0_dp,radius_fm,r,wr)
      do a=1,nc
        l=ls(a)
        do t=1,n_radial
          do i=1,n1
            basis(i,t)=spherical_bessel_j(l,grid(i)*r(t)/hbarc_MeV_fm)
          end do
          weighted(:,t)=basis(:,t)*wr(t)*r(t)*2.0_dp*alpha/(pi*hbarc_MeV_fm**2)
        end do
        vc=matmul(weighted,transpose(basis))
        v((a-1)*n1+1:a*n1,(a-1)*n1+1:a*n1)=v((a-1)*n1+1:a*n1,(a-1)*n1+1:a*n1)+vc
      end do
    end if
    call solve_k_matrix(nodes,weights,v,q,mp_MeV,pmax,result)
    result%spin=spin; result%j=j; result%isospin=1; result%ls=ls
    if(alpha>0.0_dp) then
      eta=mp_MeV*alpha/(2.0_dp*q)
      call match_coulomb(result,eta,radius_fm)
    end if
  end subroutine solve_pp_channel

  subroutine match_coulomb(result,eta,radius_fm)
    type(np_channel_result),intent(inout) :: result
    real(dp),intent(in) :: eta,radius_fm
    real(dp),allocatable :: f(:),g(:),fp(:),gp(:),f0(:),g0(:),fp0(:),gp0(:)
    real(dp) :: u(2,2),up(2,2),acoef(2,2),bcoef(2,2),left(2,2),right(2,2)
    real(dp) :: kd(2,2),kc(2,2),rho,x,eye(2,2),scale,rotation(2,2),eigen(2),work(8),tempcol(2),temp
    complex(dp) :: stemp,closure(2,2)
    integer :: nc,lmax,a,b,l,info,ipiv(2)
    nc=result%nc; lmax=maxval(result%ls(1:nc))
    allocate(f(0:lmax),g(0:lmax),fp(0:lmax),gp(0:lmax),f0(0:lmax),g0(0:lmax),fp0(0:lmax),gp0(0:lmax))
    x=result%q*radius_fm/hbarc_MeV_fm
    call coulomb_functions(x,eta,lmax,f,g,fp,gp)
    call coulomb_functions(x,0.0_dp,lmax,f0,g0,fp0,gp0)
    rho=pi*mp_MeV*result%q/2.0_dp
    kd=-rho*result%k
    eye=0.0_dp
    do a=1,nc
      eye(a,a)=1.0_dp
    end do
    u=0.0_dp; up=0.0_dp; acoef=0.0_dp; bcoef=0.0_dp
    do a=1,nc
      l=result%ls(a)
      do b=1,nc
        u(a,b)=f0(l)*eye(a,b)+g0(l)*kd(a,b)
        up(a,b)=fp0(l)*eye(a,b)+gp0(l)*kd(a,b)
        ! W=F'G-FG'=1: U=F*A+G*B; K_C=B*A^{-1}.
        acoef(a,b)=g(l)*up(a,b)-gp(l)*u(a,b)
        bcoef(a,b)=fp(l)*u(a,b)-f(l)*up(a,b)
      end do
    end do
    left=transpose(acoef); right=transpose(bcoef)
    call dgesv(nc,nc,left,2,ipiv,right,2,info)
    if(info/=0) error stop 'Singular VP Coulomb matching matrix'
    kc=transpose(right)
    scale=max(maxval(abs(kc(1:nc,1:nc))),1.0e-12_dp)
    if(maxval(abs(kc-transpose(kc)))/scale>1.0e-6_dp) error stop 'VP Coulomb K violates reciprocity'
    kc=0.5_dp*(kc+transpose(kc))
    ! Stored K/T/phase/S ALL refer to the Coulomb-distorted nuclear part.
    result%k=-kc/rho
    rotation=eye; eigen=0.0_dp
    if(nc==1) then
      eigen(1)=kc(1,1)
    else
      rotation=kc
      call dsyev('V','U',2,rotation,2,eigen,work,size(work),info)
      if(info/=0) error stop 'VP K eigendecomposition failed'
      if(abs(rotation(1,2))>abs(rotation(1,1))) then
        tempcol=rotation(:,1); rotation(:,1)=rotation(:,2); rotation(:,2)=tempcol
        temp=eigen(1); eigen(1)=eigen(2); eigen(2)=temp
      end if
      if(rotation(1,1)<0.0_dp) rotation(:,1)=-rotation(:,1)
      if(rotation(1,1)*rotation(2,2)-rotation(1,2)*rotation(2,1)<0.0_dp) rotation(:,2)=-rotation(:,2)
      result%mixing_deg=atan2(rotation(2,1),rotation(1,1))*180.0_dp/pi
    end if
    result%s=(0.0_dp,0.0_dp)
    do l=1,nc
      stemp=cmplx(1.0_dp,eigen(l),dp)/cmplx(1.0_dp,-eigen(l),dp)
      result%phase_deg(l)=atan(eigen(l))*180.0_dp/pi
      do b=1,nc
        do a=1,nc
          result%s(a,b)=result%s(a,b)+rotation(a,l)*stemp*rotation(b,l)
        end do
      end do
    end do
    result%t=(cmplx(eye,0.0_dp,dp)-result%s)/cmplx(0.0_dp,2.0_dp*rho,dp)
    closure=matmul(conjg(transpose(result%s)),result%s)-cmplx(eye,0.0_dp,dp)
    result%unitarity_error=maxval(abs(closure(1:nc,1:nc)))
  end subroutine match_coulomb
end module fb_pp_scattering
