! Direct Fortran interface to the unchanged Idaho kernel: no C/Python bridge.
module em500_force
  use fb_precision, only: dp
  use em500_parameters, only: parameters, mass_parameter_MeV, hbarc_MeV_fm
  use idaho_chiral_potential, only: chp
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  public :: initialize_em500, sd_kernel, np_kernel, pp_kernel, orbital_channels, mass_parameter_MeV, hbarc_MeV_fm
  logical, save :: initialized = .false.
contains
  subroutine initialize_em500()
    interface
      subroutine chp_preset_Idaho_N3LO(p)
        import dp
        real(dp), intent(in) :: p(0:34)
      end subroutine chp_preset_Idaho_N3LO
    end interface
    if (storage_size(1.0_dp) /= 64) error stop 'EM500 requires 64-bit real arithmetic'
    if (.not. initialized) call chp_preset_Idaho_N3LO(parameters)
    initialized = .true.
  end subroutine initialize_em500

  subroutine orbital_channels(spin,j,isospin,coupled,nc,ls)
    integer, intent(in) :: spin,j,isospin
    logical, intent(in) :: coupled
    integer, intent(out) :: nc,ls(2)
    if (spin < 0 .or. spin > 1 .or. isospin < 0 .or. isospin > 1 .or. j < 0 .or. j > 39) &
      error stop 'Invalid EM500 np spin/J/isospin: require S,T=0,1 and J=0..39'
    ls=-1
    nc=1
    if (coupled) then
      if (spin /= 1) error stop 'Coupled EM500 channel requires spin=1'
      if (j == 0) then
        ! Upstream places the single 3P0 channel in the coupled/raw ++ slot.
        ls(1)=1
      else
        nc=2
        ls=[j-1,j+1]
      end if
    else
      if (spin == 1 .and. j == 0) error stop 'EM500 3P0 requires coupled=.true.'
      ls(1)=j
    end if
    if (any(modulo(ls(1:nc)+spin+isospin,2) /= 1)) &
      error stop 'EM500 channel violates nucleon antisymmetry: L+S+T must be odd'
  end subroutine orbital_channels

  subroutine np_kernel(pout,pin,spin,j,isospin,coupled,v)
    real(dp), intent(in) :: pout,pin
    integer, intent(in) :: spin,j,isospin
    logical, intent(in) :: coupled
    real(dp), intent(out) :: v(2,2)
    real(dp) :: raw(6)
    integer :: nc,ls(2)
    if (.not. initialized) error stop 'Call initialize_em500 before np_kernel'
    call orbital_channels(spin,j,isospin,coupled,nc,ls)
    if (.not. ieee_is_finite(pout) .or. .not. ieee_is_finite(pin)) &
      error stop 'Nonfinite momentum in EM500'
    if (min(pout,pin) < 0.0_dp .or. pout+pin == 0.0_dp) &
      error stop 'Unsupported EM500 momentum: negative or both zero'
    raw=0.0_dp
    ! Rows outgoing L; columns incoming L. Upstream argument comments are
    ! reversed: executable assignments establish pout=final, pin=initial.
    ! Preserve the original radial_rephased convention; observables handle it.
    call chp(pout,pin,coupled,spin,j,isospin,0,raw)
    if (.not. all(ieee_is_finite(raw))) error stop 'Nonfinite EM500 kernel'
    v=0.0_dp
    if (coupled .and. j > 0) then
      v(1,1)=raw(4); v(1,2)=raw(6)
      v(2,1)=raw(5); v(2,2)=raw(3)
    else if (coupled) then
      v(1,1)=raw(3)
    else if (spin == 0) then
      v(1,1)=raw(1)
    else
      v(1,1)=raw(2)
    end if
  end subroutine np_kernel


  subroutine pp_kernel(pout,pin,spin,j,coupled,v)
    ! Charge-dependent EM500 strong pp interaction, T=1,Tz=-1 upstream.
    ! Long-range Coulomb is supplied separately by fb_pp_scattering.
    real(dp), intent(in) :: pout,pin
    integer, intent(in) :: spin,j
    logical, intent(in) :: coupled
    real(dp), intent(out) :: v(2,2)
    real(dp) :: raw(6)
    integer :: nc,ls(2)
    if (.not. initialized) error stop 'Call initialize_em500 before pp_kernel'
    call orbital_channels(spin,j,1,coupled,nc,ls)
    if (.not. ieee_is_finite(pout) .or. .not. ieee_is_finite(pin)) &
      error stop 'Nonfinite momentum in EM500 pp'
    if (min(pout,pin) < 0.0_dp .or. pout+pin == 0.0_dp) &
      error stop 'Unsupported EM500 pp momentum'
    call chp(pout,pin,coupled,spin,j,1,-1,raw)
    if (.not. all(ieee_is_finite(raw))) error stop 'Nonfinite EM500 pp kernel'
    v=0.0_dp
    if (coupled .and. j > 0) then
      v(1,1)=raw(4); v(1,2)=raw(6)
      v(2,1)=raw(5); v(2,2)=raw(3)
    else if (coupled) then
      v(1,1)=raw(3)
    else if (spin == 0) then
      v(1,1)=raw(1)
    else
      v(1,1)=raw(2)
    end if
  end subroutine pp_kernel

  subroutine sd_kernel(pout,pin,v)
    real(dp), intent(in) :: pout,pin
    real(dp), intent(out) :: v(2,2)
    ! Preserve the public deuteron interface, with exactly the same mapping.
    call np_kernel(pout,pin,1,1,0,.true.,v)
  end subroutine sd_kernel
end module em500_force
