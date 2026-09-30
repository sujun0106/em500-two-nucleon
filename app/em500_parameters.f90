! Fixed EM2003 Idaho N3LO500 parameters; not refitted in this project.
! Exact names, values and source hashes: docs/来源/em500.json and docs/来源/upstream.json.
module em500_parameters
  use fb_precision, only: dp
  implicit none
  private
  public :: parameters, mass_parameter_MeV, hbarc_MeV_fm
  real(dp), parameter :: mass_parameter_MeV = 938.9182046406256_dp
  real(dp), parameter :: hbarc_MeV_fm = 197.327_dp
  real(dp), parameter :: parameters(35) = [ &
    -0.145286_dp, & ! 1 Ct_1S0pp_idaho_n3lo
    -0.147167_dp, & ! 2 Ct_1S0np_idaho_n3lo
    -0.146285_dp, & ! 3 Ct_1S0nn_idaho_n3lo
    -0.118972496_dp, & ! 4 Ct_3S1_idaho_n3lo
    2.38_dp, & ! 5 C_1S0_idaho_n3lo
    1.487_dp, & ! 6 C_3P0_idaho_n3lo
    0.656_dp, & ! 7 C_1P1_idaho_n3lo
    -0.63_dp, & ! 8 C_3P1_idaho_n3lo
    0.76_dp, & ! 9 C_3S1_idaho_n3lo
    0.826_dp, & ! 10 C_3S1_3D1_idaho_n3lo
    -0.538_dp, & ! 11 C_3P2_idaho_n3lo
    -2.545_dp, & ! 12 Dh_1S0_idaho_n3lo
    -16.0_dp, & ! 13 D_1S0_idaho_n3lo
    0.245_dp, & ! 14 D_3P0_idaho_n3lo
    5.25_dp, & ! 15 D_1P1_idaho_n3lo
    2.35_dp, & ! 16 D_3P1_idaho_n3lo
    7.000000000000001_dp, & ! 17 Dh_3S1_idaho_n3lo
    6.55_dp, & ! 18 D_3S1_idaho_n3lo
    -2.8_dp, & ! 19 D_3D1_idaho_n3lo
    2.25_dp, & ! 20 Dh_3S1_3D1_idaho_n3lo
    6.61_dp, & ! 21 D_3S1_3D1_idaho_n3lo
    -1.77_dp, & ! 22 D_1D2_idaho_n3lo
    -1.46_dp, & ! 23 D_3D2_idaho_n3lo
    2.295_dp, & ! 24 D_3P2_idaho_n3lo
    -0.465_dp, & ! 25 D_3P2_3F2_idaho_n3lo
    5.66_dp, & ! 26 D_3D3_idaho_n3lo
    1.29_dp, & ! 27 gA_idaho_n3lo
    -0.81_dp, & ! 28 c1_idaho_n3lo
    2.8_dp, & ! 29 c2_idaho_n3lo
    -3.2_dp, & ! 30 c3_idaho_n3lo
    5.4_dp, & ! 31 c4_idaho_n3lo
    3.06_dp, & ! 32 d1_plus_d2_idaho_n3lo
    -3.27_dp, & ! 33 d3_idaho_n3lo
    0.45_dp, & ! 34 d5_idaho_n3lo
    -5.65_dp ] ! 35 d14_minus_d15_idaho_n3lo
end module em500_parameters
