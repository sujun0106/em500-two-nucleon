! Shared command-line entry for EM500 deuteron and np elastic scattering.
program fewbody
  use fb_precision, only: dp
  use fb_pp_driver, only: calculate_pp
  use fb_deuteron, only: deuteron_result, solve_deuteron
  use fb_deuteron_observables, only: deuteron_observables_result, calculate_deuteron_observables, write_deuteron_observables
  use em500_force, only: mass_parameter_MeV, hbarc_MeV_fm
  use fb_quadrature, only: gauss_legendre
  use fb_np_scattering, only: np_channel_result, solve_np_channel, neutron_lab_momentum
  use fb_np_observables, only: spin_amplitudes, differential_cross_section, &
    total_elastic_cross_section, optical_total_cross_section
  use, intrinsic :: iso_fortran_env, only: error_unit, output_unit, int64, compiler_version, compiler_options
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  integer :: n_grid, ios, u, jmax, n_angles, n_coulomb_points, coordinate_n_intervals
  real(dp) :: pmax_MeV, energy_lab_MeV, theta_min_deg, theta_max_deg, coulomb_radius_fm, coordinate_rmax_fm
  logical :: write_matrices, write_amplitudes, write_coordinate
  character(len=1024) :: input_file
  character(len=64) :: calculation, system, interaction, kinematics
  character(len=8) :: started_date
  character(len=10) :: started_time
  character(len=5) :: started_zone
  character(len=512) :: message
  namelist /calculation_input/ calculation, system, interaction, n_grid, pmax_MeV, write_matrices, &
    energy_lab_MeV, kinematics, jmax, theta_min_deg, theta_max_deg, n_angles, write_amplitudes, &
    coulomb_radius_fm, n_coulomb_points, write_coordinate, coordinate_rmax_fm, coordinate_n_intervals
  n_grid=100; pmax_MeV=1000.0_dp; write_matrices=.false.
  calculation='deuteron_bound_state'; system='np'; interaction='EM500'
  energy_lab_MeV=-1.0_dp; kinematics='nonrelativistic'; jmax=8
  theta_min_deg=0.0_dp; theta_max_deg=180.0_dp; n_angles=181; write_amplitudes=.true.
  coulomb_radius_fm=15.0_dp; n_coulomb_points=192
  write_coordinate=.false.; coordinate_rmax_fm=40.0_dp; coordinate_n_intervals=4000
  input_file='input.nml'
  if (command_argument_count()>1) error stop 'Usage: ./app/fewbody.x [input.nml]'
  if (command_argument_count()==1) call get_command_argument(1,input_file)
  call require_clean_outputs()
  call snapshot_input(trim(input_file))
  call date_and_time(date=started_date,time=started_time,zone=started_zone)
  call write_status('running')
  open(newunit=u,file='run_input.nml',status='old',action='read',iostat=ios,iomsg=message)
  if (ios/=0) error stop 'Cannot read input snapshot'
  read(u,nml=calculation_input,iostat=ios,iomsg=message)
  close(u)
  if (ios/=0) then
    write(error_unit,'(a)') 'Invalid calculation_input namelist: '//trim(message)
    error stop 1
  end if
  if(trim(system)/='np'.and.trim(system)/='pp') error stop 'Supported systems: np, pp'
  if(trim(calculation)/='pp_scattering'.and.trim(system)/='np') error stop 'deuteron/np calculation requires system=np'
  if (trim(interaction)/='EM500') error stop 'Only interaction=EM500 is supported'
  if (n_grid<4 .or. n_grid>1000) error stop 'Require 4 <= n_grid <= 1000'
  if (.not. ieee_is_finite(pmax_MeV)) error stop 'pmax_MeV must be finite'
  if (pmax_MeV<=0.0_dp .or. pmax_MeV>1000.0_dp) error stop 'Require 0 < pmax_MeV <= 1000'
  select case (trim(calculation))
  case ('deuteron_bound_state')
    if (energy_lab_MeV/=-1.0_dp) error stop 'Deuteron has no incident energy; omit energy_lab_MeV'
    if(write_coordinate) then
      if(.not.ieee_is_finite(coordinate_rmax_fm).or.coordinate_rmax_fm<=0) error stop 'Invalid coordinate radius'
      if(coordinate_n_intervals<2.or.mod(coordinate_n_intervals,2)/=0) error stop 'Coordinate intervals must be positive and even'
    end if
    call open_output(u,'deuteron_input.nml'); write(u,nml=calculation_input); call close_output(u)
    call calculate_deuteron()
  case ('np_scattering')
    if (.not. ieee_is_finite(energy_lab_MeV)) error stop 'Incident energy must be finite'
    if (energy_lab_MeV<=0.0_dp) error stop 'np_scattering requires energy_lab_MeV > 0'
    if (trim(kinematics)/='nonrelativistic' .and. trim(kinematics)/='relativistic') &
      error stop 'kinematics must be nonrelativistic or relativistic'
    if (jmax<0 .or. jmax>39) error stop 'Require 0 <= jmax <= 39'
    if (n_angles<1 .or. n_angles>100001) error stop 'Require 1 <= n_angles <= 100001'
    if (.not. all(ieee_is_finite([theta_min_deg,theta_max_deg]))) error stop 'Angles must be finite'
    if (theta_min_deg<0.0_dp .or. theta_max_deg>180.0_dp .or. theta_min_deg>theta_max_deg) &
      error stop 'Require 0 <= theta_min_deg <= theta_max_deg <= 180'
    if (n_angles==1 .and. theta_min_deg/=theta_max_deg) error stop 'One angle requires equal min/max'
    if (write_matrices) error stop 'write_matrices is for deuteron; np on-shell matrices are always exported'
    call open_output(u,'np_input.nml'); write(u,nml=calculation_input); call close_output(u)
    call calculate_np()
  case ('pp_scattering')
    if(trim(system)/='pp') error stop 'pp_scattering requires system=pp'
    if(.not.ieee_is_finite(energy_lab_MeV).or.energy_lab_MeV<=0) error stop 'pp requires positive finite energy'
    if(trim(kinematics)/='relativistic'.and.trim(kinematics)/='nonrelativistic') error stop 'Invalid pp kinematics'
    if(jmax<0.or.jmax>39) error stop 'Require 0 <= jmax <= 39'
    if(n_angles<1.or.n_angles>100001) error stop 'Require 1 <= n_angles <= 100001'
    if(.not.all(ieee_is_finite([theta_min_deg,theta_max_deg,coulomb_radius_fm]))) error stop 'Nonfinite pp setting'
    if(theta_min_deg<=0.or.theta_max_deg>=180.or.theta_min_deg>theta_max_deg) error stop 'pp Coulomb requires 0 < angles < 180'
    if(n_angles==1.and.theta_min_deg/=theta_max_deg) error stop 'One angle requires equal bounds'
    if(coulomb_radius_fm<=0.or.n_coulomb_points<8) error stop 'Invalid Coulomb radius/radial quadrature'
    if(write_matrices) error stop 'write_matrices applies only to deuteron'
    call open_output(u,'pp_input.nml');write(u,nml=calculation_input);call close_output(u)
    call calculate_pp(energy_lab_MeV,kinematics,n_grid,pmax_MeV,jmax,theta_min_deg,theta_max_deg,n_angles, &
      write_amplitudes,coulomb_radius_fm,n_coulomb_points)
  case default
    error stop 'Supported calculations: deuteron_bound_state, np_scattering, pp_scattering'
  end select
  call write_status('completed')
contains
  subroutine calculate_deuteron()
    type(deuteron_result) :: r
    type(deuteron_observables_result) :: obs
    integer :: i,j,a,b,row,col
  call solve_deuteron(n_grid,pmax_MeV,r)
  if(write_coordinate) then
    call calculate_deuteron_observables(r,coordinate_n_intervals,coordinate_rmax_fm,obs)
    call write_deuteron_observables(obs)
  end if

  call open_output(u,'deuteron_summary.txt')
  write(u,'(a)') 'EM500 Fortran deuteron v0.5.0'
  write(u,'(a)') 'calculation = '//trim(calculation)
  write(u,'(a)') 'system = '//trim(system)
  write(u,'(a)') 'interaction = '//trim(interaction)
  write(u,'(a)') 'Interaction description = bare EM2003 Idaho N3LO500, fixed NN-only parameters'
  write(u,'(a,i0)') 'n_grid_per_channel = ',n_grid
  write(u,'(a,i0)') 'hamiltonian_dimension = ',2*n_grid
  call scalar(u,'pmax_MeV',pmax_MeV)
  call scalar(u,'mass_parameter_MeV',mass_parameter_MeV)
  call scalar(u,'hbarc_MeV_fm',hbarc_MeV_fm)
  call scalar(u,'energy_MeV',r%energy_MeV)
  call scalar(u,'binding_MeV',r%binding_MeV)
  call scalar(u,'probability_S',r%probability_s)
  call scalar(u,'probability_D',r%probability_d)
  call scalar(u,'normalization',r%normalization)
  call scalar(u,'kinetic_MeV',r%kinetic_MeV)
  call scalar(u,'potential_MeV',r%potential_MeV)
  call scalar(u,'residual_MeV',r%residual_MeV)
  call scalar(u,'potential_asymmetry_MeV_inverse2',r%potential_asymmetry)
  call scalar(u,'potential_relative_asymmetry',r%potential_relative_asymmetry)
  call scalar(u,'hamiltonian_asymmetry_MeV',r%hamiltonian_asymmetry_MeV)
  write(u,'(a)') 'measure = integral dp p^2; no additional 2/pi'
  write(u,'(a)') 'phase = upstream radial_rephased; first S component positive'
  write(u,'(a)') 'wavefunction = psi_L(p); reduced u_L(p)=p*psi_L(p) also exported'
  write(u,'(a)') 'D probability is representation-dependent, not an experimental observable.'
  write(u,'(a)') 'Deuteron binding energy was used in interaction calibration.'
  write(u,'(a)') 'Compiler = '//compiler_version()
  write(u,'(a)') 'Compiler options = '//compiler_options()
  call close_output(u)

  call open_output(u,'deuteron_wavefunction.csv')
  write(u,'(a)') 'p_MeV,weight_MeV,psi_S_MeV^-1.5,psi_D_MeV^-1.5,'// &
      'u_S_MeV^-0.5,u_D_MeV^-0.5,prob_density_S_per_MeV,prob_density_D_per_MeV'
  do i=1,n_grid
    write(u,'(*(es25.16e3,:,","))') r%p(i),r%w(i),r%psi_s(i),r%psi_d(i), &
        r%p(i)*r%psi_s(i),r%p(i)*r%psi_d(i), &
        (r%p(i)*r%psi_s(i))**2,(r%p(i)*r%psi_d(i))**2
  end do
  call close_output(u)
  if (write_matrices) then
    call open_output(u,'deuteron_matrices.csv')
    write(u,'(a)') 'row,column,L_out,L_in,p_out_MeV,p_in_MeV,V_MeV^-2,H_MeV'
    do a=1,2
      do i=1,n_grid
        row=(a-1)*n_grid+i
        do b=1,2
          do j=1,n_grid
            col=(b-1)*n_grid+j
            write(u,'(4(i0,","),3(es25.16e3,","),es25.16e3)') &
              row,col,2*(a-1),2*(b-1),r%p(i),r%p(j), &
              r%potential_kernel(row,col),r%hamiltonian(row,col)
          end do
        end do
      end do
    end do
    call close_output(u)
  end if
  write(*,'(a,f19.12)') 'Energy E_d (MeV)    = ',r%energy_MeV
  write(*,'(a,f19.12)') 'Binding B_d (MeV)   = ',r%binding_MeV
  write(*,'(a,f19.12)') 'D probability      = ',r%probability_d
  write(*,'(a,es14.6)') 'Eigen residual MeV = ',r%residual_MeV
  write(*,'(a)') 'Outputs: deuteron_summary.txt, deuteron_wavefunction.csv, deuteron_input.nml'
  if (write_matrices) write(*,'(a)') 'Additional output: deuteron_matrices.csv'
  end subroutine calculate_deuteron
  subroutine calculate_np()
    type(np_channel_result), allocatable :: waves(:)
    real(dp), allocatable :: theta(:), dsdo(:), gx(:), gw(:), gt(:), gd(:)
    complex(dp), allocatable :: fs(:), ft(:,:,:)
    real(dp) :: q, pi, elastic, optical, angular, err_ang, err_opt, maxres, maxrec, maxunit
    integer :: jj, b, count, i, ao, ai, m, mo, mi, unit, ng, index_end
    pi=acos(-1.0_dp)
    q=neutron_lab_momentum(energy_lab_MeV,trim(kinematics))
    if (q>=pmax_MeV) error stop 'On-shell q must lie below the momentum integration cutoff'
    allocate(waves(3*jmax+2))
    count=0
    do jj=0,jmax
      count=count+1
      call solve_np_channel(0,jj,1-mod(jj,2),.false.,n_grid,pmax_MeV,q,waves(count))
      count=count+1
      call solve_np_channel(1,jj,1-mod(jj,2),.true.,n_grid,pmax_MeV,q,waves(count))
      if (jj>0) then
        count=count+1
        call solve_np_channel(1,jj,mod(jj,2),.false.,n_grid,pmax_MeV,q,waves(count))
      end if
      write(*,'(a,i0,a,i0)') 'Completed np J=',jj,'; blocks=',count
      flush(output_unit)
    end do
    maxres=maxval(waves%residual_relative); maxrec=maxval(waves%reciprocity_error_relative)
    maxunit=maxval(waves%unitarity_error)
    if (maxres>1.0e-8_dp .or. maxunit>1.0e-10_dp) error stop 'np numerical residual/unitarity check failed'
    allocate(theta(n_angles),dsdo(n_angles),fs(n_angles),ft(n_angles,3,3))
    do i=1,n_angles
      theta(i)=theta_min_deg+(theta_max_deg-theta_min_deg)*real(i-1,dp)/real(max(1,n_angles-1),dp)
    end do
    dsdo=differential_cross_section(waves,q,theta)
    elastic=total_elastic_cross_section(waves,q)
    optical=optical_total_cross_section(waves,q)
    ! The independent angular integral covers all angles regardless of output window.
    ng=max(64,2*jmax+8); allocate(gx(ng),gw(ng),gt(ng),gd(ng))
    call gauss_legendre(ng,-1.0_dp,1.0_dp,gx,gw)
    gt=acos(gx)*180.0_dp/pi
    gd=differential_cross_section(waves,q,gt)
    angular=2.0_dp*pi*sum(gw*gd)
    err_ang=abs(angular-elastic)/max(abs(elastic),tiny(1.0_dp))
    err_opt=abs(optical-elastic)/max(abs(elastic),tiny(1.0_dp))
    if (.not. all(ieee_is_finite(dsdo))) error stop 'Nonfinite np angular distribution'
    if (err_ang>1.0e-9_dp .or. err_opt>1.0e-9_dp) error stop 'np angular/optical consistency check failed'
    call open_output(unit,'np_partial_waves.csv')
    write(unit,'(a)') 'spin,J,isospin,n_channels,L1,L2,eigenphase1_deg,eigenphase2_deg,mixing_deg,'// &
      'residual_relative,reciprocity_error_relative,unitarity_error'
    do b=1,count
      write(unit,'(6(i0,","),*(es25.16e3,:,","))') waves(b)%spin,waves(b)%j,waves(b)%isospin, &
        waves(b)%nc,waves(b)%ls(1),waves(b)%ls(2),waves(b)%phase_deg,waves(b)%mixing_deg, &
        waves(b)%residual_relative,waves(b)%reciprocity_error_relative,waves(b)%unitarity_error
    end do
    call close_output(unit)
    call open_output(unit,'np_on_shell.csv')
    write(unit,'(a)') 'spin,J,isospin,L_out,L_in,K_MeV^-2,T_real_MeV^-2,T_imag_MeV^-2,S_real,S_imag'
    do b=1,count
      do ao=1,waves(b)%nc
        do ai=1,waves(b)%nc
          write(unit,'(5(i0,","),*(es25.16e3,:,","))') waves(b)%spin,waves(b)%j,waves(b)%isospin, &
            waves(b)%ls(ao),waves(b)%ls(ai),waves(b)%k(ao,ai), &
            real(waves(b)%t(ao,ai),dp),aimag(waves(b)%t(ao,ai)), &
            real(waves(b)%s(ao,ai),dp),aimag(waves(b)%s(ao,ai))
        end do
      end do
    end do
    call close_output(unit)
    call open_output(unit,'np_angular_distribution.csv')
    write(unit,'(a)') 'theta_cm_deg,dsigma_domega_mb_sr'
    do i=1,n_angles
      write(unit,'(*(es25.16e3,:,","))') theta(i),dsdo(i)
    end do
    call close_output(unit)
    if (write_amplitudes) then
      call spin_amplitudes(waves,q,theta,fs,ft)
      call open_output(unit,'np_spin_amplitudes.csv')
      write(unit,'(a)') 'theta_cm_deg,spin,m_out,m_in,f_real_fm,f_imag_fm'
      do i=1,n_angles
        write(unit,'(es25.16e3,",",3(i0,","),*(es25.16e3,:,","))') theta(i),0,0,0,real(fs(i),dp),aimag(fs(i))
        do mo=-1,1
          do mi=-1,1
            write(unit,'(es25.16e3,",",3(i0,","),*(es25.16e3,:,","))') theta(i),1,mo,mi, &
              real(ft(i,mo+2,mi+2),dp),aimag(ft(i,mo+2,mi+2))
          end do
        end do
      end do
      call close_output(unit)
    end if
    call open_output(unit,'np_j_convergence.csv')
    write(unit,'(a)') 'J_max,total_elastic_mb,optical_total_mb,dsigma_0_mb_sr,dsigma_90_mb_sr,dsigma_180_mb_sr'
    do jj=0,jmax
      index_end=3*jj+2
      write(unit,'(i0,",",*(es25.16e3,:,","))') jj, &
        total_elastic_cross_section(waves(:index_end),q),optical_total_cross_section(waves(:index_end),q), &
        differential_cross_section(waves(:index_end),q,[0.0_dp,90.0_dp,180.0_dp])
    end do
    call close_output(unit)
    call open_output(unit,'np_summary.txt')
    write(unit,'(a)') 'EM500 Fortran np elastic scattering v0.5.0'
    write(unit,'(a)') 'calculation = np_scattering'
    write(unit,'(a)') 'system = neutron incident on stationary proton'
    write(unit,'(a)') 'interaction = bare EM2003 Idaho N3LO500; NN-only; no inelastic channels'
    write(unit,'(a)') 'kinematics = '//trim(kinematics)
    call scalar(unit,'energy_lab_MeV',energy_lab_MeV)
    call scalar(unit,'q_MeV',q)
    call scalar(unit,'LS_energy_parameter_MeV',q*q/mass_parameter_MeV)
    call scalar(unit,'mass_parameter_MeV',mass_parameter_MeV)
    call scalar(unit,'hbarc_MeV_fm',hbarc_MeV_fm)
    call scalar(unit,'pmax_MeV',pmax_MeV)
    write(unit,'(a,i0)') 'n_grid = ',n_grid
    write(unit,'(a,i0)') 'J_max = ',jmax
    write(unit,'(a,i0)') 'number_of_blocks = ',count
    call scalar(unit,'total_elastic_mb',elastic)
    call scalar(unit,'optical_total_mb',optical)
    call scalar(unit,'angular_integral_mb',angular)
    call scalar(unit,'angular_integral_relative_error',err_ang)
    call scalar(unit,'optical_relative_error',err_opt)
    call scalar(unit,'maximum_residual_relative',maxres)
    call scalar(unit,'maximum_reciprocity_error_relative',maxrec)
    call scalar(unit,'maximum_unitarity_error',maxunit)
    write(unit,'(a)') 'angle_definition = outgoing neutron, center-of-mass frame'
    write(unit,'(a)') 'on_shell_basis = radial_rephased; amplitudes convert by U_L=i^L to plain momentum'
    write(unit,'(a)') 'phase_definition = principal eigenphases modulo 180 deg; channel-1 dominated first'
    write(unit,'(a)') 'single-channel L2=-1; phase2=0 and mixing=0 are unused placeholders'
    write(unit,'(a)') 'relativistic kinematics changes q only, not the LS interaction or dynamics'
    write(unit,'(a)') 'unitarity/optical checks do not establish grid/J convergence or physical model accuracy'
    write(unit,'(a)') 'Compiler = '//compiler_version()
    write(unit,'(a)') 'Compiler options = '//compiler_options()
    call close_output(unit)
    write(*,'(a,f18.9)') 'Neutron lab energy MeV = ',energy_lab_MeV
    write(*,'(a,f18.9)') 'Total elastic mb      = ',elastic
    write(*,'(a,es14.6)') 'Angular closure error = ',err_ang
    write(*,'(a)') 'Outputs: np_summary.txt, np_angular_distribution.csv, np_on_shell.csv, np_partial_waves.csv'
  end subroutine calculate_np
  ! A direct invocation must never silently overwrite an unarchived run.
  subroutine require_clean_outputs()
    character(len=32), parameter :: names(*) = [character(len=32) :: &
      'run_input.nml','run_status.txt','deuteron_input.nml','deuteron_summary.txt', &
      'deuteron_wavefunction.csv','deuteron_matrices.csv','deuteron_status.txt', &
      'deuteron_observables.dat','deuteron_coordinate_wavefunction.csv', &
      'np_input.nml','np_summary.txt','np_partial_waves.csv','np_on_shell.csv', &
      'np_angular_distribution.csv','np_spin_amplitudes.csv','np_j_convergence.csv', &
      'pp_input.nml','pp_summary.txt','pp_partial_waves.csv','pp_nuclear_s.csv', &
      'pp_angular_distribution.csv','pp_spin_amplitudes.csv','pp_j_convergence.csv']
    integer :: k
    logical :: exists
    do k=1,size(names)
      inquire(file=trim(names(k)),exist=exists)
      if (exists) then
        write(error_unit,'(a)') 'Existing run file: '//trim(names(k))
        write(error_unit,'(a)') 'Archive or move the existing run files before starting a new calculation.'
        error stop 1
      end if
    end do
  end subroutine require_clean_outputs

  subroutine snapshot_input(path)
    character(len=*), intent(in) :: path
    integer :: source, target, stat, chunk
    integer(int64) :: byte_count, position
    character(len=8192) :: buffer
    character(len=512) :: msg
    open(newunit=source,file=path,status='old',action='read',access='stream', &
        form='unformatted',iostat=stat,iomsg=msg)
    if (stat/=0) then
      write(error_unit,'(a)') 'Cannot read input: '//path//' : '//trim(msg)
      error stop 1
    end if
    inquire(unit=source,size=byte_count,iostat=stat,iomsg=msg)
    if (stat/=0 .or. byte_count<0_int64) error stop 'Cannot determine input file size'
    open(newunit=target,file='run_input.nml',status='new',action='write',access='stream', &
        form='unformatted',iostat=stat,iomsg=msg)
    if (stat/=0) then
      write(error_unit,'(a)') 'Cannot create run_input.nml: '//trim(msg)
      error stop 1
    end if
    position=1_int64
    do while (position<=byte_count)
      chunk=int(min(int(len(buffer),int64),byte_count-position+1_int64))
      read(source,pos=position,iostat=stat,iomsg=msg) buffer(:chunk)
      if (stat/=0) then
        write(error_unit,'(a)') 'Cannot copy input: '//trim(msg)
        error stop 1
      end if
      write(target) buffer(:chunk)
      position=position+int(chunk,int64)
    end do
    close(source)
    call close_output(target)
  end subroutine snapshot_input

  subroutine write_status(state)
    character(len=*), intent(in) :: state
    integer :: unit, stat
    character(len=8) :: finished_date
    character(len=10) :: finished_time
    character(len=5) :: finished_zone
    character(len=512) :: msg
    if (state=='running') then
      call open_output(unit,'run_status.txt')
    else
      ! Replace only the status record created by this invocation.
      open(newunit=unit,file='run_status.txt',status='replace',action='write',iostat=stat,iomsg=msg)
      if (stat/=0) then
        write(error_unit,'(a)') 'Cannot update run_status.txt: '//trim(msg)
        error stop 1
      end if
    end if
    write(unit,'(a)') 'version=0.5.0'
    write(unit,'(a)') 'calculation='//trim(calculation)
    write(unit,'(a)') 'system='//trim(system)
    write(unit,'(a)') 'interaction='//trim(interaction)
    write(unit,'(a)') 'started_date='//started_date
    write(unit,'(a)') 'started_time='//started_time
    write(unit,'(a)') 'started_zone='//started_zone
    if (state=='completed') then
      call date_and_time(date=finished_date,time=finished_time,zone=finished_zone)
      write(unit,'(a)') 'finished_date='//finished_date
      write(unit,'(a)') 'finished_time='//finished_time
      write(unit,'(a)') 'finished_zone='//finished_zone
    end if
    write(unit,'(a)') 'status='//state
    call close_output(unit)
  end subroutine write_status

  subroutine open_output(unit,path)
    integer, intent(out) :: unit
    character(len=*), intent(in) :: path
    integer :: stat
    character(len=512) :: msg
    open(newunit=unit,file=path,status='new',action='write',iostat=stat,iomsg=msg)
    if (stat/=0) then
      write(error_unit,'(a)') 'Cannot write '//trim(path)//' : '//trim(msg)
      error stop 1
    end if
  end subroutine open_output

  subroutine close_output(unit)
    integer, intent(in) :: unit
    integer :: stat
    character(len=512) :: msg
    close(unit,iostat=stat,iomsg=msg)
    if (stat/=0) then
      write(error_unit,'(a)') 'Cannot finish writing output: '//trim(msg)
      error stop 1
    end if
  end subroutine close_output

  subroutine scalar(unit,key,value)
    integer, intent(in) :: unit
    character(len=*), intent(in) :: key
    real(dp), intent(in) :: value
    write(unit,'(a," = ",es25.16e3)') key,value
  end subroutine scalar
end program fewbody
