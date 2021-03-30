#ifdef DOUBLE_PREC
MODULE m_interpolation_r64
#else
MODULE m_interpolation_r32
#endif

  ! Purpose:
  !   To provide routines to interpolate functions of one, two and three variables defined along monotonic vectors. Supported
  !   interpolation methods are linear and nearest neighbour. At nodes outside the input grid, values can either be set to zero, to
  !   the value of the closest node or linearly extrapolated.
  !
  !   The calling sequence is as follows:
  !     a) subroutine "setup_interpolation" (called only once), to set interpolation and extrapolation methods
  !     b) subroutine "interpolate" (called as many times as necessary) to interpolate data
  !
  !   If interpolation occurs on a set of points (i.e. vectors), routines are parallelized at module level based on OPENMP. Routines
  !   interpolating at a single point are thread-safe and can be parallelized in the calling program.
  !
  !   Precision is determined at compile time by adding (or not) the "-DDOUBLE_PREC" flag.
  !
  ! Revisions:
  !     Date                    Description of change
  !     ====                    =====================
  !   02/09/20                  original version
  !

  USE, NON_INTRINSIC :: m_precisions

  IMPLICIT none

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  PRIVATE

  PUBLIC :: setup_interpolation, interpolate, interpolation_error

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  INTERFACE interpolate
    MODULE PROCEDURE interp_1d_scalar, interp_2d_scalar, interp_3d_scalar, interp_1d_vector, interp_2d_vector, interp_3d_vector
  END INTERFACE interpolate

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  PROCEDURE(linear), POINTER :: fun => NULL(), extrafun => NULL()

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  ! set precisions of input/output arguments at compile-time
#ifdef DOUBLE_PREC
  INTEGER, PARAMETER :: r__ = r64
#else
  INTEGER, PARAMETER :: r__ = r32
#endif

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  CONTAINS

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE setup_interpolation(method, extrapolation, ok)

      ! Purpose:
      !   To initialise the interpolation module by defining interpolation ("method") and extrapolation ("extrapolation") methods.
      !   For the former, valid arguments are "linear" and (nearest) "neighbour". Extrapolation can be "linear", "zero" or "copy".
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   02/09/20                  original version
      !

      CHARACTER(*), INTENT(IN)  :: method, extrapolation
      INTEGER(i32), INTENT(OUT) :: ok

      !-----------------------------------------------------------------------------------------------------------------------------

      ok = 0

      SELECT CASE(method)
        CASE('linear')
          fun => linear
        CASE('neighbour')
          fun => neighbour
        CASE DEFAULT
          ok = 1
      END SELECT

      SELECT CASE(extrapolation)
        CASE('linear')
          extrafun => extralin
        CASE('zero')
          extrafun => extrazero
        CASE('copy')
          extrafun => extracopy
        !CASE('none')
        CASE DEFAULT
          ok = 2
      END SELECT

    END SUBROUTINE

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE interp_core(x, y, xo, yo, ix)

      ! Purpose:
      !   To interpolate function "y" defined at "x" at single point "xo". Integer "ix" represents the index of "x" such that
      !   "x(ix) <= xo < x(ix + 1)". If "ix=0" or "ix=size(x)", "xo" is outside input range "x" and extrapolation will occur. On input,
      !   "ix" is used as initial position (hint) for bisection algorithm, while on output it is updated to the actual index value.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   02/09/20                  original version
      !

      REAL(r__),    DIMENSION(:),       INTENT(IN)    :: x
      REAL(r__),    DIMENSION(SIZE(x)), INTENT(IN)    :: y
      REAL(r__),                        INTENT(IN)    :: xo
      REAL(r__),                        INTENT(OUT)   :: yo
      INTEGER(i32),                     INTENT(INOUT) :: ix
      INTEGER(i32)                                    :: n
      REAL(r__),    DIMENSION(2)                      :: xl, yl

      !-----------------------------------------------------------------------------------------------------------------------------

      ix = hunt(x, xo, ix)

      n = SIZE(x)

      IF ( (ix .eq. 0) .or. (ix .eq. n) ) THEN
        yo = extrafun(x, y, xo, ix)          !< does not work with PGI when nested
        ! yo = extralin(x, y, xo, ix)
      ELSE
        yo = fun(x, y, xo, ix)               !< does not work with PGI when nested
        ! yo = linear(x, y, xo, ix)
      ENDIF

    END SUBROUTINE interp_core

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    PURE FUNCTION linear(x, y, xo, i) RESULT(yo)

      ! Purpose:
      !   To perform linear interpolation.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   02/09/20                  original version
      !

      REAL(r__),    DIMENSION(:),       INTENT(IN) :: x
      REAL(r__),    DIMENSION(SIZE(x)), INTENT(IN) :: y
      REAL(r__),                        INTENT(IN) :: xo
      INTEGER(i32),                     INTENT(IN) :: i
      REAL(r__)                                    :: yo
      REAL(r__),    DIMENSION(2)                   :: xl, yl

      !-----------------------------------------------------------------------------------------------------------------------------

      xl = x(i:i + 1)
      yl = y(i:i + 1)

      yo = yl(1) + (xo - xl(1)) * (yl(2) - yl(1)) / (xl(2) - xl(1))

    END FUNCTION linear

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    PURE FUNCTION neighbour(x, y, xo, i) RESULT(yo)

      ! Purpose:
      !   To perform nearest neighbour interpolation.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   02/09/20                  original version
      !

      REAL(r__),    DIMENSION(:),       INTENT(IN) :: x
      REAL(r__),    DIMENSION(SIZE(x)), INTENT(IN) :: y
      REAL(r__),                        INTENT(IN) :: xo
      INTEGER(i32),                     INTENT(IN) :: i
      REAL(r__)                                    :: yo

      !-----------------------------------------------------------------------------------------------------------------------------

      IF ( ABS(xo - x(i)) .le. ABS(xo - x(i + 1)) ) THEN
        yo = y(i)
      ELSE
        yo = y(i + 1)
      ENDIF

    END FUNCTION neighbour

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    PURE FUNCTION extrazero(x, y, xo, i) RESULT(yo)

      ! Purpose:
      !   To extrapolate by setting value to zero.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   02/09/20                  original version
      !

      REAL(r__),    DIMENSION(:),       INTENT(IN) :: x
      REAL(r__),    DIMENSION(SIZE(x)), INTENT(IN) :: y
      REAL(r__),                        INTENT(IN) :: xo
      INTEGER(i32),                     INTENT(IN) :: i
      REAL(r__)                                    :: yo

      !-----------------------------------------------------------------------------------------------------------------------------

      yo = 0._r__

    END FUNCTION extrazero

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    PURE FUNCTION extracopy(x, y, xo, i) RESULT(yo)

      ! Purpose:
      !   To extrapolate by copying closest value
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   02/09/20                  original version
      !

      REAL(r__),    DIMENSION(:),       INTENT(IN) :: x
      REAL(r__),    DIMENSION(SIZE(x)), INTENT(IN) :: y
      REAL(r__),                        INTENT(IN) :: xo
      INTEGER(i32),                     INTENT(IN) :: i
      REAL(r__)                                    :: yo
      REAL(r__),    DIMENSION(2)                   :: xl, yl

      !-----------------------------------------------------------------------------------------------------------------------------

      IF (i .eq. 0) THEN
        yo = y(1)
      ELSE
        yo = y(SIZE(x))
      ENDIF

    END FUNCTION extracopy

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    PURE FUNCTION extralin(x, y, xo, i) RESULT(yo)

      ! Purpose:
      !   To apply linear extrapolation.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   02/09/20                  original version
      !

      REAL(r__),    DIMENSION(:),       INTENT(IN) :: x
      REAL(r__),    DIMENSION(SIZE(x)), INTENT(IN) :: y
      REAL(r__),                        INTENT(IN) :: xo
      INTEGER(i32),                     INTENT(IN) :: i
      REAL(r__)                                    :: yo

      !-----------------------------------------------------------------------------------------------------------------------------

      IF (i .eq. 0) THEN
        yo = linear(x, y, xo, 1)
      ELSEIF (i .eq. SIZE(x)) THEN
        yo = linear(x, y, xo, SIZE(x) - 1)
      ENDIF

    END FUNCTION extralin

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE interp_1d_scalar(x, y, xo, yo)

      ! Purpose:
      !   To interpolate a function of one variable at a single point.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   02/09/20                  original version
      !

      REAL(r__),    DIMENSION(:), INTENT(IN)    :: x, y
      REAL(r__),                  INTENT(IN)    :: xo
      REAL(r__),                  INTENT(OUT)   :: yo
      INTEGER(i32)                              :: pos

      !-----------------------------------------------------------------------------------------------------------------------------

      pos = 0

      CALL interp_core(x, y, xo, yo, pos)

    END SUBROUTINE interp_1d_scalar

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE interp_1d_vector(x, y, xo, yo)

      ! Purpose:
      !   To interpolate a function of one variable at a set of points.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   02/09/20                  original version
      !

      REAL(r__),    DIMENSION(:),        INTENT(IN)    :: x, y
      REAL(r__),    DIMENSION(:),        INTENT(IN)    :: xo
      REAL(r__),    DIMENSION(SIZE(xo)), INTENT(OUT)   :: yo
      INTEGER(i32)                                     :: i, i0

      !-----------------------------------------------------------------------------------------------------------------------------

      i0 = 0

      !$omp parallel do default(shared) firstprivate(i0) private(i)
      DO i = 1, SIZE(xo)
        CALL interp_core(x, y, xo(i), yo(i), i0)
      ENDDO
      !$omp end parallel do

    END SUBROUTINE interp_1d_vector

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE interp_2d_scalar(x, y, f, xo, yo, fo)

      ! Purpose:
      !   To interpolate a function of two variables at a single point.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   02/09/20                  original version
      !

      REAL(r__),    DIMENSION(:),               INTENT(IN)    :: x, y
      REAL(r__),    DIMENSION(SIZE(x),SIZE(y)), INTENT(IN)    :: f
      REAL(r__),                                INTENT(IN)    :: xo, yo
      REAL(r__),                                INTENT(OUT)   :: fo
      INTEGER(i32)                                            :: n, j0, j1, i0, tmp
      REAL(r__),    DIMENSION(2)                              :: fy

      !-----------------------------------------------------------------------------------------------------------------------------

      n = SIZE(y)

      j0 = hunt(y, yo, 0)

      IF (j0 .eq. 0) THEN
        j0 = 1
      ELSEIF (j0 .eq. n) THEN
        j0 = n - 1
      ENDIF

      j1 = j0 + 1

      i0  = 0
      tmp = 0

      CALL interp_core(x, f(:, j0), xo, fy(1), i0)
      CALL interp_core(x, f(:, j1), xo, fy(2), i0)

      CALL interp_core(y(j0:j1), fy, yo, fo, tmp)

    END SUBROUTINE interp_2d_scalar

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE interp_2d_vector(x, y, f, xo, yo, fo)

      ! Purpose:
      !   To interpolate a function of two variables at a set of points.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   02/09/20                  original version
      !

      REAL(r__),    DIMENSION(:),                 INTENT(IN)  :: x, y
      REAL(r__),    DIMENSION(SIZE(x),SIZE(y)),   INTENT(IN)  :: f
      REAL(r__),    DIMENSION(:),                 INTENT(IN)  :: xo, yo
      REAL(r__),    DIMENSION(SIZE(xo),SIZE(yo)), INTENT(OUT) :: fo
      INTEGER(i32)                                            :: i, j, n
      INTEGER(i32)                                            :: i0, j0, j1, tmp
      REAL(r__),    DIMENSION(2)                              :: fy

      !-----------------------------------------------------------------------------------------------------------------------------

      n = SIZE(y)

      j0 = 0

      !$omp parallel do default(shared) firstprivate(j0) private(i, j, j1, i0, tmp, fy)
      DO j = 1, SIZE(yo)

        j0 = hunt(y, yo(j), j0)

        IF (j0 .eq. 0) THEN
          j0 = 1
        ELSEIF (j0 .eq. n) THEN
          j0 = n - 1
        ENDIF

        j1 = j0 + 1

        i0  = 0
        tmp = 0

        DO i = 1, SIZE(xo)
          CALL interp_core(x, f(:, j0), xo(i), fy(1), i0)
          CALL interp_core(x, f(:, j1), xo(i), fy(2), i0)
          CALL interp_core(y(j0:j1), fy, yo(j), fo(i, j), tmp)
        ENDDO

      ENDDO
      !$omp end parallel do

    END SUBROUTINE interp_2d_vector

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE interp_3d_scalar(x, y, z, f, xo, yo, zo, fo)

      ! Purpose:
      !   To interpolate a function of three variables at a single point.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   02/09/20                  original version
      !

      REAL(r__),    DIMENSION(:),                       INTENT(IN)  :: x, y, z
      REAL(r__),    DIMENSION(SIZE(x),SIZE(y),SIZE(z)), INTENT(IN)  :: f
      REAL(r__),                                        INTENT(IN)  :: xo, yo, zo
      REAL(r__),                                        INTENT(OUT) :: fo
      INTEGER(i32)                                                  :: m, n
      INTEGER(i32)                                                  :: tmp, i0, j0, j1, k0, k1
      REAL(r__),    DIMENSION(2)                                    :: fy, fz

      !-----------------------------------------------------------------------------------------------------------------------------

      m = SIZE(y)
      n = SIZE(z)

      k0 = hunt(z, zo, 0)

      IF (k0 .eq. 0) THEN
        k0 = 1
      ELSEIF (k0 .eq. n) THEN
        k0 = n - 1
      ENDIF

      k1 = k0 + 1

      n = SIZE(y)

      j0 = hunt(y, yo, 0)

      IF (j0 .eq. 0) THEN
        j0 = 1
      ELSEIF (j0 .eq. n) THEN
        j0 = n - 1
      ENDIF

      j1 = j0 + 1

      i0  = 0
      tmp = 0

      CALL interp_core(x, f(:, j0, k0), xo, fy(1), i0)
      CALL interp_core(x, f(:, j1, k0), xo, fy(2), i0)

      CALL interp_core(y(j0:j1), fy, yo, fz(1), tmp)

      CALL interp_core(x, f(:, j0, k1), xo, fy(1), i0)
      CALL interp_core(x, f(:, j1, k1), xo, fy(2), i0)

      CALL interp_core(y(j0:j1), fy, yo, fz(2), tmp)

      CALL interp_core(z(k0:k1), fz, zo, fo, tmp)

    END SUBROUTINE interp_3d_scalar

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE interp_3d_vector(x, y, z, f, xo, yo, zo, fo)

      ! Purpose:
      !   To interpolate a function of three variables at a set of points.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   02/09/20                  original version
      !

      REAL(r__),    DIMENSION(:),                          INTENT(IN)  :: x, y, z
      REAL(r__),    DIMENSION(SIZE(x),SIZE(y),SIZE(z)),    INTENT(IN)  :: f
      REAL(r__),    DIMENSION(:),                          INTENT(IN)  :: xo, yo, zo
      REAL(r__),    DIMENSION(SIZE(xo),SIZE(yo),SIZE(zo)), INTENT(OUT) :: fo
      INTEGER(i32)                                                     :: i, j, k, n, m
      INTEGER(i32)                                                     :: i0, j0, j1, k0, k1, tmp
      REAL(r__),    DIMENSION(2)                                       :: fy, fz

      !-----------------------------------------------------------------------------------------------------------------------------

      m = SIZE(y)
      n = SIZE(z)

      k0 = 0

      !$omp parallel do default(shared) firstprivate(k0) private(i, j, k, k1, j0, j1, i0, tmp, fy, fz)
      DO k = 1, SIZE(zo)

        k0 = hunt(z, zo(k), k0)

        IF (k0 .eq. 0) THEN
          k0 = 1
        ELSEIF (k0 .eq. n) THEN
          k0 = n - 1
        ENDIF

        k1 = k0 + 1

        j0 = 0

        DO j = 1, SIZE(yo)

          j0 = hunt(y, yo(j), j0)

          IF (j0 .eq. 0) THEN
            j0 = 1
          ELSEIF (j0 .eq. m) THEN
            j0 = m - 1
          ENDIF

          j1 = j0 + 1

          i0  = 0
          tmp = 0

          DO i = 1, SIZE(xo)

            CALL interp_core(x, f(:, j0, k0), xo(i), fy(1), i0)
            CALL interp_core(x, f(:, j1, k0), xo(i), fy(2), i0)

            CALL interp_core(y(j0:j1), fy, yo(j), fz(1), tmp)

            CALL interp_core(x, f(:, j0, k1), xo(i), fy(1), i0)
            CALL interp_core(x, f(:, j1, k1), xo(i), fy(2), i0)

            CALL interp_core(y(j0:j1), fy, yo(j), fz(2), tmp)

            CALL interp_core(z(k0:k1), fz, zo(k), fo(i, j, k), tmp)

          ENDDO

        ENDDO

      ENDDO
      !$omp end parallel do


    END SUBROUTINE interp_3d_vector

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    FUNCTION hunt(xx, x, hint) RESULT(jlo)

      ! Purpose:
      !   To compute index "jlo" such that "xx(jlo) <= x < x(jlo + 1)". If result is zero or size(xx), "x" is outside input range.
      !   If "hint" is not zero, it is used as a starting point for search and execution time is sensibly reduced if "x" is in its
      !   neighbourhood.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   02/09/20                  original version, based on Press et al. (1996)
      !

      REAL(r__),    DIMENSION(:), INTENT(IN) :: xx
      REAL(r__),                  INTENT(IN) :: x
      INTEGER(i32),               INTENT(IN) :: hint
      INTEGER(i32)                           :: n, inc, jhi, jm, jlo
      LOGICAL                                :: ascnd

      !-----------------------------------------------------------------------------------------------------------------------------

      n = SIZE(xx)

      ascnd = (xx(n) .ge. xx(1))

      jlo = hint

      IF (jlo .le. 0 .or. jlo .gt. n) THEN

        jlo = 0
        jhi = n + 1

      ELSE

        inc = 1

        ! hunt up
        IF (x .ge. xx(jlo) .eqv. ascnd) THEN

          DO
            jhi = jlo + inc

            IF (jhi .gt. n) THEN
              jhi = n + 1
              EXIT
            ELSE
              IF (x .lt. xx(jhi) .eqv. ascnd) EXIT
              jlo = jhi
              inc = inc + inc
            ENDIF

          ENDDO

        ! hunt down
        ELSE

          jhi = jlo

          DO
            jlo = jhi - inc

            IF (jlo .lt. 1) THEN
              jlo = 0
              EXIT
            ELSE
              IF (x .ge. xx(jlo) .eqv. ascnd) EXIT
              jhi = jlo
              inc = inc + inc
            ENDIF

          ENDDO

        ENDIF

      ENDIF

      ! bisection
      DO

        IF (jhi - jlo .le. 1) THEN
          IF (x .eq. xx(n)) jlo = n - 1
          IF (x .eq. xx(1)) jlo = 1
          EXIT
        ELSE
          jm = (jhi + jlo) / 2

          IF (x .ge. xx(jm) .eqv. ascnd) THEN
            jlo = jm
          ELSE
            jhi = jm
          ENDIF
        ENDIF

      ENDDO

    END FUNCTION hunt

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    FUNCTION interpolation_error(ierr) RESULT(msg)

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
          msg = 'wrong argument in setup_interpolation: unknown interpolation method ("method")'

        CASE(2)
          msg = 'wrong argument in setup_interpolation: unknown extrapolation method ("extrapolation")'

      END SELECT

    END FUNCTION interpolation_error

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

#ifdef DOUBLE_PREC
END MODULE m_interpolation_r64
#else
END MODULE m_interpolation_r32
#endif
