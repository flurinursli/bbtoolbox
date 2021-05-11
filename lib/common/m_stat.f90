MODULE m_stat

  USE                :: omp_lib
  USE, NON_INTRINSIC :: m_precisions

  IMPLICIT none

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  PRIVATE

  PUBLIC :: variance, mean, parallel_variance

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  ! set precisions of input/output arguments at compile-time
#ifdef DOUBLE_PREC
  INTEGER, PARAMETER :: r__ = r64
#else
  INTEGER, PARAMETER :: r__ = r32
#endif

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  INTERFACE variance
    MODULE PROCEDURE variance_1d, variance_2d, variance_3d
  END INTERFACE

  INTERFACE mean
    MODULE PROCEDURE mean_1d, mean_2d, mean_3d
  END INTERFACE

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  CONTAINS

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    REAL(r__) FUNCTION variance_1d(r)

      ! Purpose:
      !   To compute variance according to the compensated-summation version of the two-pass algorithm. Calculations in double-precision.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   04/05/20                  original version
      !

      REAL(r__),   DIMENSION(:), INTENT(IN) :: r                      !< vector of data
      INTEGER(i32)                          :: i, n
      REAL(r64)                             :: mu, s1, s2, x, v

      !-----------------------------------------------------------------------------------------------------------------------------

      n = SIZE(r)

      mu = mean(r)

      s1 = 0._r64
      s2 = 0._r64

      DO i = 1, n
        x = REAL(r(i), r64) - mu
        s1 = s1 + x
        s2 = s2 + x**2
      ENDDO

      s1 = (s1**2) / REAL(n, r64)

      v = (s2 - s1) / REAL(n - 1, r64)

      variance_1d = REAL(v, r__)

    END FUNCTION variance_1d

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    REAL(r__) FUNCTION variance_2d(r)

      ! Purpose:
      !   To compute variance according to the compensated-summation version of the two-pass algorithm. Calculations in double-precision.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   04/05/20                  original version
      !

      REAL(r__),   DIMENSION(:,:), INTENT(IN) :: r                    !< array of data
      INTEGER(i32)                            :: i, j
      INTEGER(i32)                            :: nx, ny
      REAL(r64)                               :: mu, s1, s2, x, v

      !-----------------------------------------------------------------------------------------------------------------------------

      nx = SIZE(r, 1)
      ny = SIZE(r, 2)

      mu = mean(r)

      s1 = 0._r64
      s2 = 0._r64

      DO j = 1, ny
        DO i = 1, nx
          x = REAL(r(i, j), r64) - mu
          s1 = s1 + x
          s2 = s2 + x**2
        ENDDO
      ENDDO

      s1 = (s1**2) / REAL(SIZE(r), r64)

      v = (s2 - s1) / REAL(SIZE(r) - 1, r64)

      variance_2d = REAL(v, r__)

    END FUNCTION variance_2d

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    REAL(r__) FUNCTION variance_3d(r)

      ! Purpose:
      !   To compute variance according to the compensated-summation version of the two-pass algorithm. Calculations in double-precision.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   04/05/20                  original version
      !

      REAL(r__),   DIMENSION(:,:,:), INTENT(IN) :: r                    !< array of data
      INTEGER(i32)                              :: i, j, k
      INTEGER(i32)                              :: nx, ny, nz
      REAL(r64)                                 :: mu, s1, s2, x, v

      !-----------------------------------------------------------------------------------------------------------------------------

      nx = SIZE(r, 1)
      ny = SIZE(r, 2)
      nz = SIZE(r, 3)

      mu = mean(r)

      s1 = 0._r64
      s2 = 0._r64

      DO k = 1, nz
        DO j = 1, ny
          DO i = 1, nx
            x = REAL(r(i, j, k), r64) - mu
            s1 = s1 + x
            s2 = s2 + x**2
          ENDDO
        ENDDO
      ENDDO

      s1 = (s1**2) / REAL(SIZE(r), r64)

      v = (s2 - s1) / REAL(SIZE(r) - 1, r64)

      variance_3d = REAL(v, r__)

    END FUNCTION variance_3d

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    REAL(r__) FUNCTION mean_1d(r)

      ! Purpose:
      !   To compute average. Calculations in double-precision to reduce risk of cancellation.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   04/05/20                  original version
      !

      REAL(r__),  DIMENSION(:), INTENT(IN) :: r               !< vector of data
      INTEGER(i32)                          :: i
      INTEGER(i32)                          :: n
      REAL(r64)                            :: v, c

      !-----------------------------------------------------------------------------------------------------------------------------

      n = SIZE(r)

      v = 0._r64
      c = 1._r64

      DO i = 1, n
        v = v + (REAL(r(i), r64) - v) / c
        c = c + 1._r64
      ENDDO

      ! return mean at desired precision
      mean_1d = REAL(v, r__)

    END FUNCTION mean_1d

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    REAL(r__) FUNCTION mean_2d(r)

      ! Purpose:
      !   To compute average. Calculations in double-precision to reduce risk of cancellation.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   04/05/20                  original version
      !

      REAL(r__),    DIMENSION(:,:), INTENT(IN) :: r                       !< array of data
      INTEGER(i32)                             :: i, j
      INTEGER(i32)                             :: nx, ny
      REAL(r64)                                :: v, c

      !-----------------------------------------------------------------------------------------------------------------------------

      nx = SIZE(r, 1)
      ny = SIZE(r, 2)

      v = 0._r64
      c = 1._r64

      DO j = 1, ny
        DO i = 1, nx
          v = v + (REAL(r(i, j), r64) - v) / c
          c = c + 1._r64
        ENDDO
      ENDDO

      ! return mean with desired precision
      mean_2d = REAL(v, r__)

    END FUNCTION mean_2d

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    REAL(r__) FUNCTION mean_3d(r)

      ! Purpose:
      !   To compute average. Calculations in double-precision to reduce risk of cancellation.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   04/05/20                  original version
      !

      REAL(r__),  DIMENSION(:,:,:), INTENT(IN) :: r                       !< array of data
      INTEGER(i32)                              :: i, j, k
      INTEGER(i32)                              :: nx, ny, nz
      REAL(r64)                                :: v, c

      !-----------------------------------------------------------------------------------------------------------------------------

      nx = SIZE(r, 1)
      ny = SIZE(r, 2)
      nz = SIZE(r, 3)

      v = 0._r64
      c = 1._r64

      DO k = 1, nz
        DO j = 1, ny
          DO i = 1, nx
            v = v + (REAL(r(i, j, k), r64) - v) / c
            c = c + 1._r64
          ENDDO
        ENDDO
      ENDDO

      ! return mean with desired precision
      mean_3d = REAL(v, r__)

    END FUNCTION mean_3d

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE parallel_variance(varset, avgset, nset, var, avg)

      ! Purpose:
      !   To compute variance and mean of a sequence of variance and mean values, each computed from a set of "nset" elements, based
      !   on the algorithm of Chan et al. (1979). Internal calculations in double-precision.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   04/05/20                  original version
      !

      REAL(r__),   DIMENSION(:), INTENT(IN)  :: varset                             !< set of variance values
      REAL(r__),   DIMENSION(:), INTENT(IN)  :: avgset                             !< set of mean values
      INTEGER(i32), DIMENSION(:), INTENT(IN)  :: nset                               !< number of points for each set
      REAL(r__),                 INTENT(OUT) :: var                                !< resulting variance
      REAL(r__),                 INTENT(OUT) :: avg                                !< resulting average
      INTEGER(i32)                            :: i
      REAL(r64)                              :: sig, mu, n, delta, m1, m2, m

      !-----------------------------------------------------------------------------------------------------------------------------

      ! variance, mean and number of points of first set
      sig = REAL(varset(1), r64)
      mu  = REAL(avgset(1), r64)
      n   = REAL(nset(1), r64)

      ! loop over set of variance/mean values
      DO i = 2, SIZE(nset)

        delta = REAL(avgset(i), r64) - mu

        m1 = sig * (n - 1._r64)
        m2 = REAL(varset(i), r64) * (nset(i) - 1._r64)

        m = m1 + m2 + delta**2 * n * REAL(nset(i), r64) / (n + REAL(nset(i), r64))

        ! resulting mean
        mu = (mu * n + REAL(avgset(i), r64) * REAL(nset(i), r64)) / (n + REAL(nset(i), r64))

        ! resulting number of points
        n = n + REAL(nset(i), r64)

        ! resulting variance
        sig = m / (n - 1._r64)

      ENDDO

      ! return with desired precision
      var = REAL(sig, r__)
      avg = REAL(mu, r__)

    END SUBROUTINE parallel_variance

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *


END MODULE m_stat
