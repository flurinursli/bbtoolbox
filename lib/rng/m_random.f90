MODULE m_random

  ! Purpose:
  !   To draw random numbers, independent of hardware and compiler, according to the following distributions:
  !      1) uniform, in the interval [a, b) with a < b
  !      2) normal, with mean "mu" and standard deviation "sigma" (> 0)
  !
  !   The calling sequence is as follows:
  !     a) subroutine "setup_rng" (called only once) to set distribution, seed number and starting position along the stream.
  !     b) subroutine "rng" (called as many times as necessary)
  !
  !   Parallelization at module level is based on OPENMP and is enabled when "rng" is called from a serial region. In this case, the
  !   random number calculation is distributed internally amongst threads by means of the "block-splitting" method.
  !   Similarly, parallelization at the calling program can be obtained by using MPI, OPENMP, etc. as long as each process/thread
  !   initialise the random number generator indicating the amount of variates to be skipped (i.e. already drawn). For OPENMP, this
  !   implies placing both calls to "setup_rng" and "rng" inside parallel regions.
  !   Both approaches work out as expected as long as nested OPENMP parallel regions are disabled (default).
  !
  !   Precision is determined at compile-time ("-DDOUBLE_PREC" flag)
  !
  !   Routines rely on the TRNG random number generator library (version 4.0 or higher)
  !
  ! Revisions:
  !     Date                    Description of change
  !     ====                    =====================
  !   02/09/20                  original version
  !

  USE                :: omp_lib
  USE, NON_INTRINSIC :: m_precisions

  IMPLICIT none

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  PRIVATE

  PUBLIC :: setup_rng, rng, rng_error

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  ! set precisions of input/output arguments at compile-time
#ifdef DOUBLE_PREC
  INTEGER, PARAMETER :: r__   = r64
  INTEGER, PARAMETER :: c_r__ = c_r64
#else
  INTEGER, PARAMETER :: r__   = r32
  INTEGER, PARAMETER :: c_r__ = c_r32
#endif

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  ! interface to TRNG C++ routines
  INTERFACE

    SUBROUTINE setup_trng(seed, skip) BIND(c, name="setup_trng")
      USE, INTRINSIC     :: iso_c_binding
      USE, NON_INTRINSIC :: m_precisions
      INTEGER(c_i32), VALUE :: seed
      INTEGER(c_i64), VALUE :: skip
    END SUBROUTINE setup_trng

    SUBROUTINE uni01(a, b, npts, rand) BIND(c, name="uni01")
      USE, INTRINSIC     :: iso_c_binding
      USE, NON_INTRINSIC :: m_precisions
      IMPLICIT none
#ifdef DOUBLE_PREC
      INTEGER, PARAMETER :: c_r__ = c_r64
#else
      INTEGER, PARAMETER :: c_r__ = c_r32
#endif
      REAL(c_r__),                     VALUE       :: a, b
      INTEGER(c_i32),                  VALUE       :: npts
      REAL(c_r__),    DIMENSION(npts), INTENT(OUT) :: rand
    END SUBROUTINE uni01

    SUBROUTINE uni(a, b, npts, rand) BIND(c, name="uni")
      USE, INTRINSIC     :: iso_c_binding
      USE, NON_INTRINSIC :: m_precisions
      IMPLICIT none
#ifdef DOUBLE_PREC
      INTEGER, PARAMETER :: c_r__ = c_r64
#else
      INTEGER, PARAMETER :: c_r__ = c_r32
#endif
      REAL(c_r__),                     VALUE       :: a, b
      INTEGER(c_i32),                  VALUE       :: npts
      REAL(c_r__),    DIMENSION(npts), INTENT(OUT) :: rand
    END SUBROUTINE uni

    SUBROUTINE norm(mean, std, npts, rand) BIND(c, name="norm")
      USE, INTRINSIC     :: iso_c_binding
      USE, NON_INTRINSIC :: m_precisions
      IMPLICIT none
#ifdef DOUBLE_PREC
      INTEGER, PARAMETER :: c_r__ = c_r64
#else
      INTEGER, PARAMETER :: c_r__ = c_r32
#endif
      REAL(c_r__),                     VALUE       :: mean, std
      INTEGER(c_i32),                  VALUE       :: npts
      REAL(c_r__),    DIMENSION(npts), INTENT(OUT) :: rand
    END SUBROUTINE norm

  END INTERFACE

  INTERFACE rng
    MODULE PROCEDURE rng_1d, rng_2d, rng_3d
  END INTERFACE rng

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  ! first and second input parameter for distributions
  REAL(c_r__)    :: c_par1, c_par2

  ! seed number
  INTEGER(c_i32) :: c_seed

  ! random variates to be skipped
  INTEGER(c_i64) :: c_skip

  !$omp threadprivate(c_skip)

  ! generic pointer to distribution routines
  PROCEDURE(uni), POINTER :: fun

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  CONTAINS

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE setup_rng(ok, dist, p1, p2, seed, skip)

      ! Purpose:
      !   To set up the random number generator by providing the desired distribution and its parameters, the seed number and,
      !   optionally, the amount of random numbers already extracted. The latter parameter is useful for parallel execution relying
      !   on the "block splitting" technique (see TRNG manual).
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   02/09/20                  original version
      !

      INTEGER(i32),               INTENT(OUT) :: ok                 !< error flag
      CHARACTER(len=*),           INTENT(IN)  :: dist               !< desired distribution
      REAL(r__),                  INTENT(IN)  :: p1, p2             !< distribution parameters
      INTEGER(i32),               INTENT(IN)  :: seed               !< seed number
      INTEGER(i32),     OPTIONAL, INTENT(IN)  :: skip               !< amount of random numbers already consumed (e.g. by a thread)

      !-----------------------------------------------------------------------------------------------------------------------------

      ok = 0

      c_skip = 0

      IF (PRESENT(skip)) c_skip = skip

      c_seed = seed

      c_par1 = p1
      c_par2 = p2

      SELECT CASE(dist)
        CASE('uniform')
          IF ( (p1 .eq. 0._r__) .and. (p2 .eq. 1._r__) ) THEN
            fun => uni01
          ELSE
            fun => uni
          ENDIF
          ! for uniform distributions, condition "a < b" must hold
          IF (p1 .ge. p2) ok = 2
        CASE('normal')
          fun => norm
          ! for normal distributions, std.dev. must be positive
          IF (p2 .le. 0._r__) ok = 3
        CASE DEFAULT
          ok = 1
      END SELECT

      ! setup trng generator
      CALL setup_trng(c_seed, c_skip)

    END SUBROUTINE setup_rng

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    FUNCTION rng_error(ierr) RESULT(msg)

      ! Purpose:
      !   to translate an error code into a text for all routines in current module.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   02/09/20                  original version
      !

      INTEGER(i32),             INTENT(IN) :: ierr
      CHARACTER(:), ALLOCATABLE            :: msg

      !-----------------------------------------------------------------------------------------------------------------------------

      SELECT CASE(ierr)
        CASE(1)
          msg = 'wrong value for argument "dist" in setup_rng: unknown distribution type'

        CASE(2)
          msg = 'wrong value for argument "p1" and "p2" in setup_rng: uniform distributions require "p1 < p2"'

        CASE(3)
          msg = 'wrong value for argument "p2" in setup_rng: normal distributions must have positive std.dev. ("p2 > 0")'

      END SELECT

    END FUNCTION rng_error

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE rng_1d(x, nt)

      ! Purpose:
      !   to draw a sequence of random numbers when "x" is a vector. If input is not contiguous in memory, it is copied in&out. "nt"
      !   can be used to set the number of threads for parallel execution.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   02/09/20                  original version
      !

      REAL(r__), DIMENSION(:), CONTIGUOUS, TARGET, INTENT(OUT) :: x

#include "rng_incl.f90"

    END SUBROUTINE rng_1d

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE rng_2d(x, nt)

      ! Purpose:
      !   to draw a sequence of random numbers when "x" is a 2d matrix. If input is not contiguous in memory, it is copied in&out.
      !   "nt" can be used to set the number of threads for parallel execution.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   02/09/20                  original version
      !

      REAL(r__), DIMENSION(:,:), CONTIGUOUS, TARGET, INTENT(OUT) :: x

#include "rng_incl.f90"

    END SUBROUTINE rng_2d

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE rng_3d(x, nt)

      ! Purpose:
      !   to draw a sequence of random numbers when "x" is a 3d matrix. If input is not contiguous in memory, it is copied in&out.
      !   "nt" can be used to set the number of threads for parallel execution.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   02/09/20                  original version
      !

      REAL(r__), DIMENSION(:,:,:), CONTIGUOUS, TARGET, INTENT(OUT) :: x

#include "rng_incl.f90"

    END SUBROUTINE rng_3d

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *


END MODULE m_random
