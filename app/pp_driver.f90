! Production pp output driver. Nuclear S matrices are relative to Coulomb waves.
module fb_pp_driver
  use fb_precision, only: dp
  use fb_np_scattering, only: np_channel_result
  use fb_pp_scattering, only: solve_pp_channel, proton_lab_momentum, pp_coulomb_eta
  use fb_pp_observables, only: pp_spin_amplitudes, pp_differential_cross_section
  use, intrinsic :: iso_fortran_env, only: output_unit, compiler_version, compiler_options
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: calculate_pp
contains
  subroutine calculate_pp(energy,kinematics,n_grid,pmax,jmax,theta_min,theta_max,n_angles,write_amplitudes,radius,n_radial)
    real(dp), intent(in) :: energy,pmax,theta_min,theta_max,radius
    character(len=*), intent(in) :: kinematics
    integer, intent(in) :: n_grid,jmax,n_angles,n_radial
    logical, intent(in) :: write_amplitudes
    type(np_channel_result), allocatable :: waves(:)
    integer, allocatable :: ends(:)
    real(dp), allocatable :: theta(:),dsdo(:),reflected(:)
    complex(dp), allocatable :: fs(:),ft(:,:,:)
    real(dp) :: q,eta,symmetry,maxres,maxunit,maxrec,check_values(5)
    integer :: j,n,count,w,a,b,u,i,mi,mo
    logical :: relc
    relc=trim(kinematics)=='relativistic'
    q=proton_lab_momentum(energy)
    eta=pp_coulomb_eta(q,relativistic_coulomb=relc)
    allocate(waves(2*(jmax+1)),ends(0:jmax))
    count=0
    do j=0,jmax
      if (mod(j,2)==0) then
        count=count+1
        call solve_pp_channel(0,j,.false.,n_grid,pmax,q,radius,n_radial,waves(count),relativistic_coulomb=relc)
        count=count+1
        call solve_pp_channel(1,j,.true.,n_grid,pmax,q,radius,n_radial,waves(count),relativistic_coulomb=relc)
      else
        count=count+1
        call solve_pp_channel(1,j,.false.,n_grid,pmax,q,radius,n_radial,waves(count),relativistic_coulomb=relc)
      end if
      ends(j)=count
      write(*,'(a,i0,a,i0)') 'Completed pp J=',j,'; blocks=',count
      flush(output_unit)
    end do
    maxres=maxval(waves(:count)%residual_relative)
    maxunit=maxval(waves(:count)%unitarity_error)
    maxrec=maxval(waves(:count)%reciprocity_error_relative)
    if(maxres>1e-8_dp.or.maxunit>1e-8_dp.or.maxrec>1e-8_dp) error stop 'pp residual/unitarity check failed'
    allocate(theta(n_angles),dsdo(n_angles),reflected(n_angles),fs(n_angles),ft(n_angles,3,3))
    do i=1,n_angles
      theta(i)=theta_min+(theta_max-theta_min)*real(i-1,dp)/real(max(1,n_angles-1),dp)
    end do
    dsdo=pp_differential_cross_section(waves(:count),q,eta,theta)
    reflected=pp_differential_cross_section(waves(:count),q,eta,180.0_dp-theta)
    if (.not.all(ieee_is_finite(dsdo)).or.any(dsdo<0)) error stop 'Invalid pp cross section'
    symmetry=maxval(abs(dsdo-reflected)/max(dsdo,1e-30_dp))
    if (symmetry>1e-8_dp) error stop 'pp identical-particle symmetry failed'
    open(newunit=u,file='pp_angular_distribution.csv',status='new')
    write(u,'(a)') 'theta_cm_deg,dsigma_domega_mb_sr'
    do i=1,n_angles
      write(u,'(*(es25.16e3,:,","))')theta(i),dsdo(i)
    end do
    close(u)
    open(newunit=u,file='pp_nuclear_s.csv',status='new')
    write(u,'(a)') 'spin,J,isospin,L_out,L_in,S_coulomb_nuclear_real,S_coulomb_nuclear_imag'
    do w=1,count
      do b=1,waves(w)%nc
        do a=1,waves(w)%nc
          write(u,'(5(i0,","),*(es25.16e3,:,","))') waves(w)%spin,waves(w)%j,1, &
            waves(w)%ls(a),waves(w)%ls(b),real(waves(w)%s(a,b),dp),aimag(waves(w)%s(a,b))
        end do
      end do
    end do
    close(u)
    open(newunit=u,file='pp_partial_waves.csv',status='new')
    write(u,'(a)') 'spin,J,isospin,n_channels,L1,L2,residual_relative,reciprocity_error_relative,unitarity_error'
    do w=1,count
      write(u,'(6(i0,","),*(es25.16e3,:,","))') waves(w)%spin,waves(w)%j,1,waves(w)%nc, &
        waves(w)%ls,waves(w)%residual_relative,waves(w)%reciprocity_error_relative,waves(w)%unitarity_error
    end do
    close(u)
    open(newunit=u,file='pp_j_convergence.csv',status='new')
    write(u,'(a)') 'J_max,dsigma_10_mb_sr,dsigma_30_mb_sr,dsigma_60_mb_sr,dsigma_90_mb_sr,dsigma_150_mb_sr'
    do j=0,jmax
      check_values=pp_differential_cross_section(waves(:ends(j)),q,eta,[10.0_dp,30.0_dp,60.0_dp,90.0_dp,150.0_dp])
      write(u,'(i0,",",*(es25.16e3,:,","))')j,check_values
    end do
    close(u)
    if(write_amplitudes) then
      call pp_spin_amplitudes(waves(:count),q,eta,theta,fs,ft)
      open(newunit=u,file='pp_spin_amplitudes.csv',status='new')
      write(u,'(a)') 'theta_cm_deg,spin,m_out,m_in,f_real_fm,f_imag_fm'
      do i=1,n_angles
        write(u,'(es25.16e3,",",3(i0,","),*(es25.16e3,:,","))')theta(i),0,0,0,real(fs(i),dp),aimag(fs(i))
        do mi=-1,1
          do mo=-1,1
            write(u,'(es25.16e3,",",3(i0,","),*(es25.16e3,:,","))')theta(i),1,mo,mi, &
              real(ft(i,mo+2,mi+2),dp),aimag(ft(i,mo+2,mi+2))
          end do
        end do
      end do
      close(u)
    end if
    open(newunit=u,file='pp_summary.txt',status='new')
    write(u,'(a)') 'EM500 Fortran pp scattering v0.4.0'
    write(u,'(a)') 'calculation = pp_scattering'
    write(u,'(a)') 'system = proton incident on stationary proton'
    write(u,'(a)') 'interaction = EM2003 Idaho N3LO500 pp + point Coulomb'
    write(u,'(a)') 'coulomb_method = Vincent-Phatak finite-radius matching'
    write(u,'(a)') 'kinematics = '//trim(kinematics)
    call scalar(u,'energy_lab_MeV',energy);call scalar(u,'q_MeV',q);call scalar(u,'eta',eta)
    call scalar(u,'pmax_MeV',pmax);call scalar(u,'coulomb_radius_fm',radius)
    write(u,'(a,i0)') 'n_grid = ',n_grid
    write(u,'(a,i0)') 'J_max = ',jmax
    write(u,'(a,i0)') 'n_coulomb_points = ',n_radial
    write(u,'(a,i0)') 'number_of_blocks = ',count
    call scalar(u,'maximum_reciprocity_error_relative',maxrec)
    if(relc) then
      write(u,'(a)') 'coulomb_parameter_convention = relativistic effective alpha(q); static nonrelativistic LS propagator'
    else
      write(u,'(a)') 'coulomb_parameter_convention = static alpha'
    end if
    call scalar(u,'maximum_residual_relative',maxres);call scalar(u,'maximum_unitarity_error',maxunit)
    call scalar(u,'theta_reflection_relative_error',symmetry)
    write(u,'(a)') 'angle_definition = outgoing proton, center-of-mass frame'
    write(u,'(a)') 'identical_particles = singlet symmetric and triplet antisymmetric; no 1/2 in point differential cross section'
    write(u,'(a)') 'full_solid_angle_integral = divergent for point Coulomb; no finite optical-total closure asserted'
    write(u,'(a)') 'nuclear_s_basis = radial_rephased, Coulomb-relative; observables add Coulomb phases and exchange'
    write(u,'(a)') 'not_included = additional two-photon, vacuum-polarization, magnetic-moment electromagnetic corrections; pion-production channels'
    write(u,'(a)') 'Compiler = '//compiler_version()
    write(u,'(a)') 'Compiler options = '//compiler_options()
    close(u)
    write(*,'(a,f18.9)') 'Proton lab energy MeV = ',energy
    write(*,'(a,es14.6)') 'Angular reflection error = ',symmetry
    write(*,'(a)') 'Outputs: pp_summary.txt, pp_angular_distribution.csv, pp_nuclear_s.csv'
  end subroutine calculate_pp
  subroutine scalar(u,key,value)
    integer,intent(in)::u
    character(len=*),intent(in)::key
    real(dp),intent(in)::value
    write(u,'(a," = ",es25.16e3)')key,value
  end subroutine scalar
end module fb_pp_driver
