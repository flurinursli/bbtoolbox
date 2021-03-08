MODULE m_precisions

! Purpose:
!   to provide a set of constants defining compiler-independent precision for integer, real and complex numbers. Precision at compile-
!   time is determined by preprocessor flag "DOUBLE_PREC".
!
! Revisions:
!     Date                    Description of change
!     ====                    =====================
!   02/09/20                  original version
!

USE, INTRINSIC :: iso_fortran_env
USE, INTRINSIC :: iso_c_binding

IMPLICIT none

PUBLIC

! set FORTRAN precision
INTEGER, PARAMETER :: r32 = real32
INTEGER, PARAMETER :: r64 = real64
INTEGER, PARAMETER :: i32 = int32
INTEGER, PARAMETER :: i64 = int64

! set precision for C/C++ interfaces
INTEGER, PARAMETER :: c_i32 = c_int
INTEGER, PARAMETER :: c_i64 = c_long
INTEGER, PARAMETER :: c_c32 = c_float_complex
INTEGER, PARAMETER :: c_c64 = c_double_complex
INTEGER, PARAMETER :: c_r64 = c_double
INTEGER, PARAMETER :: c_r32 = c_float

! set default precisions at compile-time
! #ifdef DOUBLE_PREC
! INTEGER, PARAMETER :: r__   = r64
! INTEGER, PARAMETER :: c_r__ = c_r64
! INTEGER, PARAMETER :: c_c__ = c_c64
! #else
! INTEGER, PARAMETER :: r__   = r32
! INTEGER, PARAMETER :: c_r__ = c_r32
! INTEGER, PARAMETER :: c_c__ = c_c32
! #endif

END MODULE m_precisions
