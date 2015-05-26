module precision
  !! Provides kind attributes for setting machine-compiler-independent precision for real numbers. Here, we define "single precision" to mean 32-bit precision, "double precision" to mean 64-bit precision, and "quadruple precision" to mean 128-bit precision. This ensures that precision is preserved regardless of the compiler or computer architecture.

  use ISO_FORTRAN_ENV
  
  implicit none
  
  integer, parameter :: sp = REAL32 !selected_real_kind(6, 37)
  !! Single precision: 32-bit real
  integer, parameter :: dp = REAL64 !selected_real_kind(15, 307)
  !! Double precision: 64-bit real
  integer, parameter :: qp = REAL128 !selected_real_kind(33, 4931)
  !! Quadruple precision: 128-bit real

  integer, parameter :: kd = dp
  !! Internal kind used within this library. All real kinds in this library are defined with this parameter.
  
contains 
  
end module precision
  
