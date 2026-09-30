! Common arithmetic kind; matches the original EM500 REAL(8) on supported compilers.
module fb_precision
  use, intrinsic :: iso_fortran_env, only: real64
  implicit none
  private
  integer, parameter, public :: dp = real64
end module fb_precision
