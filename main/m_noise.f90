MODULE m_noise

  USE, INTRINSIC     :: ieee_arithmetic, ONLY: ieee_is_nan
  USE, NON_INTRINSIC :: m_precisions
  USE, NON_INTRINSIC :: m_logfile
  USE, NON_INTRINSIC :: m_fft_real                        !< note: transform are always in r64
  USE, NON_INTRINSIC :: m_random
  USE, NON_INTRINSIC :: m_toolbox, ONLY: input, timeseries

  IMPLICIT none

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  PRIVATE

  PUBLIC :: generate_noise, generate_noise_cdra

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  REAL(r32), PARAMETER :: TWOPI = 3.14159265358979323846_r64 * 2._r64

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  CONTAINS

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE generate_noise(ok, rank)

      INTEGER(i32),                             INTENT(OUT) :: ok
      INTEGER(i32),                             INTENT(IN)  :: rank
      COMPLEX(r32)                                          :: sij
      COMPLEX(r32), ALLOCATABLE, DIMENSION(:)               :: z
      INTEGER(i32)                                          :: i, l, p, comp, fr, rcvr, n, npts, row, col, lu, error
      LOGICAL                                               :: is_col_odd, is_row_odd, is_singular
      REAL(r32)                                             :: freq, df, eta, lambda, phi, sigma, mean
      REAL(r32),                 DIMENSION(2)               :: rnorm
      REAL(r32),    ALLOCATABLE, DIMENSION(:,:)             :: fc_a, fc_b, sii, c, v

      !-----------------------------------------------------------------------------------------------------------------------------

      ok = 0
      error = 0

      npts = SIZE(timeseries%sp%time)

      CALL make_fftw_plan([npts])

      df = 1._r32 / (npts * timeseries%sp%dt)

      npts = npts / 2 + 1

      ALLOCATE(z(npts))
      ALLOCATE(fc_a(npts, SIZE(input%receiver)), fc_b(npts, SIZE(input%receiver)), sii(npts, SIZE(input%receiver)))

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! ------------------------------------------------- initial time-series  -----------------------------------------------------

      ! first time-series is fundamentally white noise with rms amplitude ~ 1
      ! CALL setup_rng(ok, 'uniform', -SQRT(3._r32), SQRT(3._r32), input%coda%seed)
      CALL setup_rng(ok, 'normal', 0., 1., input%coda%seed)

      ! components for same receiver are not inter-correlated
      CALL rng(timeseries%sp%x(:, 1))
      CALL rng(timeseries%sp%y(:, 1))
      CALL rng(timeseries%sp%z(:, 1))

      CALL setup_rng(ok, 'normal', 0._r32, 1._r32, input%coda%seed)

      DO comp = 1, 3

#ifdef ERROR_TRAP
        IF (error .ne. 0) CYCLE
#endif

        ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --
        ! ---------------------------- conditionally simulate ground motion at rcvr-th receiver  -----------------------------------
        DO rcvr = 2, SIZE(input%receiver)

#ifdef ERROR_TRAP
          IF (error .ne. 0) CYCLE
#endif

          IF (comp .eq. 1) THEN
            CALL fft(timeseries%sp%x(:, rcvr-1), z)
          ELSEIF (comp .eq. 2) THEN
            CALL fft(timeseries%sp%y(:, rcvr-1), z)
          ELSE
            CALL fft(timeseries%sp%z(:, rcvr-1), z)
          ENDIF

          ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
          ! ------------------------------------------ add PSD of previous time-series  --------------------------------------------
          DO fr = 1, npts
            fc_a(fr, rcvr-1) = REAL(z(fr))
            fc_b(fr, rcvr-1) = AIMAG(z(fr))
            sii(fr, rcvr-1)  = ABS(z(fr))**2
          ENDDO

          sii(:, rcvr) = sii(:, 1)         !< assign PSD of first time-series to rcvr (used in [rcvr*2-1 rcvr*2] columns/rows)

          n = 2*rcvr                       !< size of covariance matrix

          ALLOCATE(c(n, n), v(n, n))

          ! copy zero frequency: all time-series are fully correlated at this frequency
          fc_a(1, rcvr) = fc_a(1, 1)
          fc_b(1, rcvr) = fc_b(1, 1)

          ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
          ! ---------------------------------- find fourier coefficients at fr-th frequency  ---------------------------------------
          DO fr = 2, npts

            freq = (fr - 1) * df

            ! initialise covariance matrix and its inverse
            DO col = 1, n
              DO row = 1, n
                c(row, col) = 0._r32
                v(row, col) = 0._r32
              ENDDO
            ENDDO

            ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
            ! ----------------------------------------------- build covariance matrix  ---------------------------------------------
            DO col = 1, n

              l = CEILING(col / 2._r32)            !< receiver index for columns

              is_col_odd = MOD(col, 2) .ne. 0

              DO row = 1, n

                p = CEILING(row / 2._r32)           !< receiver index for rows

                ! inter-station distance (must be in meter)
                eta = HYPOT(input%receiver(p)%x - input%receiver(l)%x, input%receiver(p)%y - input%receiver(l)%y)

                ! cross spectral density
                sij = CMPLX(lagged_coherency(freq, eta) * SQRT(sii(fr,l) * sii(fr,p)), 0._r32)

                is_row_odd = MOD(row, 2) .ne. 0

                ! odd columns
                IF (is_col_odd) THEN

                  IF (row .gt. col + 1) THEN          !< fill entries below diagonal (first element below diagonal is always 0)
                    IF (is_row_odd) THEN
                      c(row, col) = REAL(sij)       !< K
                    ELSE
                      c(row, col) = -AIMAG(sij)      !< -Q
                    ENDIF
                  ENDIF

                ! even columns
                ELSE

                  IF (row .gt. col) THEN          !< fill entries below diagonal
                    IF (is_row_odd) THEN
                      c(row, col) = AIMAG(sij)   !< Q
                    ELSE
                      c(row, col) = REAL(sij)    !< K
                    ENDIF
                  ENDIF

                ENDIF

                ! diagonal elements
                IF (row .eq. col) c(row, col) = sii(fr, p)   !< G

              ENDDO
            ENDDO

            ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
            ! ----------------------------------------------- invert covariance matrix  --------------------------------------------
            CALL matinv(ok, c, v, is_singular)

            CALL rng(rnorm)        !< draw two numbers from normal distribution

            ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
            ! --------------------------------------------- compute fourier coefficients  ------------------------------------------
            IF (is_singular) THEN

              ! when matrix is singular (i.e. coherency almost 1), fourier coefficients are identical to those of previous timeseries
              fc_a(fr, rcvr) = fc_a(fr, rcvr-1)
              fc_b(fr, rcvr) = fc_b(fr, rcvr-1)

            ELSE

              phi    = 0._r32
              lambda = 0._r32

              l = 2*rcvr - 1

              DO i = 1, rcvr-1
                phi    = phi    + ( v(2*i - 1, l)     * fc_a(fr, i) -  &
                                    v(2*i - 1, l + 1) * fc_b(fr, i) )
                lambda = lambda + ( v(2*i - 1, l + 1) * fc_a(fr, i) +  &
                                    v(2*i - 1, l)     * fc_b(fr, i) )
              ENDDO

              sigma = SQRT(1._r32 / v(l, l))

#ifdef ERROR_TRAP
              IF (ieee_is_nan(sigma)) error = 1
#endif

              mean = -phi / v(l, l)

              fc_a(fr, rcvr) = mean + rnorm(1)*sigma

              mean = -lambda / v(l, l)

              fc_b(fr, rcvr) = mean + rnorm(2)*sigma

            ENDIF

          ENDDO       !< end loop over frequency

          ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
          ! ----------------------------------------- build spectrum, return timeseries  -------------------------------------------

          DO fr = 1, npts
            z(fr) = CMPLX(fc_a(fr, rcvr), fc_b(fr, rcvr))
          ENDDO

          IF (comp .eq. 1) THEN
            CALL ifft(timeseries%sp%x(:, rcvr), z)
          ELSEIF (comp .eq. 2) THEN
            CALL ifft(timeseries%sp%y(:, rcvr), z)
          ELSEIF (comp .eq. 3) THEN
            CALL ifft(timeseries%sp%z(:, rcvr), z)
          ENDIF

          DEALLOCATE(c, v)

        ENDDO     !< end loop over receivers

      ENDDO  !< end loop over components

      npts = SIZE(timeseries%sp%time)

      CALL destroy_fftw_plan([npts])

!#ifdef DEBUG
      OPEN(newunit = lu, file = 'noise_dbg.txt', status = 'unknown', form = 'formatted', access = 'sequential', action = 'write', &
           iostat = ok)
      DO i = 1, SIZE(timeseries%sp%time)
        WRITE(lu, *) timeseries%sp%time(i), timeseries%sp%x(i, :)
      ENDDO
      CLOSE(lu)
!#endif

#ifdef ERROR_TRAP
      IF (error .ne. 0) THEN
        IF (rank .eq. 0) CALL report_error('generate_noise - ERROR: matrix is singular')
        ok = error
      ENDIF
#endif

    END SUBROUTINE generate_noise

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE matinv(ok, a, inv, is_singular)

      ! Purpose:
      !   to find the inverse "inv" of symmetric matrix "a", i.e. find "inv" such that "a*inv = id" where "id" is the identity matrix.
      !   "is_singular" is set to true if "a" is singular.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   08/03/21                  original version
      !

      INTEGER(i32),                                              INTENT(OUT)   :: ok
      REAL(r32),                 DIMENSION(:,:),                 INTENT(INOUT) :: a
      REAL(r32),                 DIMENSION(SIZE(a,1),SIZE(a,2)), INTENT(OUT)   :: inv
      LOGICAL,                                                   INTENT(OUT)   :: is_singular
      INTEGER(i32)                                                             :: i, j, n, info
      INTEGER(i32)                                                             :: lwork
      INTEGER(i32),              DIMENSION(SIZE(a,1))                          :: ipiv, iwork
      REAL(r32)                                                                :: rcond
      REAL(r32),                 DIMENSION(SIZE(a,1))                          :: ferr, berr
      REAL(r32),                 DIMENSION(SIZE(a,1),SIZE(a,2))                :: af, id
      REAL(r32),    ALLOCATABLE, DIMENSION(:)                                  :: work

      !-----------------------------------------------------------------------------------------------------------------------------

      is_singular = .false.

      n = SIZE(a, 1)

      ALLOCATE(work(1))

      lwork = -1

      ! CALL ssysv('L', n, n, a, n, ipiv, inv, n, work, lwork, info)
      CALL ssysvx('N', 'L', n, n, a, n, af, n, ipiv, id, n, inv, n, rcond, ferr, berr, work, lwork, iwork, info)

      lwork = NINT(work(1))

      DEALLOCATE(work)
      ALLOCATE(work(lwork))

      DO j = 1, n
        DO i = 1, n
          inv(i, j) = 0._r32
          id(i, j)  = 0._r32
        ENDDO
      ENDDO

      DO i = 1, n
        inv(i, i) = 1._r32
        id(i, i)  = 1._r32
      ENDDO

      ! CALL ssysv('L', n, n, a, n, ipiv, inv, n, work, lwork, info)
      CALL ssysvx('N', 'L', n, n, a, n, af, n, ipiv, id, n, inv, n, rcond, ferr, berr, work, lwork, iwork, info)

      IF (info .gt. 0) is_singular = .true.

      ok = info

    END SUBROUTINE matinv

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    REAL(r32) FUNCTION lagged_coherency(freq, eta)

      ! Purpose:
      !   to compute (lagged) coherency at a given frequency "freq" and separation distance "eta" (in meters) according to the model
      !   of Luco and Wang.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   08/03/21                  original version
      !

      REAL(r32), INTENT(IN) :: freq, eta

      !-----------------------------------------------------------------------------------------------------------------------------

      lagged_coherency = EXP(-(input%coda%alpha * TWOPI * freq * eta)**2)

    END FUNCTION

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE generate_noise_cdra(ok, rank)

      INTEGER(i32),                             INTENT(OUT) :: ok
      INTEGER(i32),                             INTENT(IN)  :: rank
      COMPLEX(r32)                                          :: sij
      COMPLEX(r32), ALLOCATABLE, DIMENSION(:)               :: z
      INTEGER(i32)                                          :: i, l, p, comp, fr, rcvr, n, npts, row, col, lu, error
      LOGICAL                                               :: is_col_odd, is_row_odd, is_singular
      REAL(r32)                                             :: freq, df, eta, lambda, phi, sigma, mean
      REAL(r32),    ALLOCATABLE, DIMENSION(:)               :: rnorm
      REAL(r32),    ALLOCATABLE, DIMENSION(:,:)             :: fc_a, fc_b, sii, c

      !-----------------------------------------------------------------------------------------------------------------------------

      ok = 0
      error = 0

      npts = SIZE(timeseries%sp%time)

      CALL make_fftw_plan([npts])

      df = 1._r32 / (npts * timeseries%sp%dt)

      npts = npts / 2 + 1

      n = SIZE(input%receiver)

      ALLOCATE(z(npts))
      ALLOCATE(fc_a(npts, n), fc_b(npts, n), sii(npts, n))
      ALLOCATE(rnorm(n), c(n, n))

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! ------------------------------------------------- initial time-series  -----------------------------------------------------

      CALL setup_rng(ok, 'normal', 0._r32, 1._r32, input%coda%seed)

      DO comp = 1, 3

#ifdef ERROR_TRAP
        IF (error .ne. 0) CYCLE
#endif

        fc_a(1, :) = 0._r32
        fc_b(1, :) = 0._r32

        ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --
        ! ---------------------------------- find fourier coefficients at fr-th frequency  -----------------------------------------
        DO fr = 2, npts

          freq = (fr - 1) * df

          ! initialise covariance matrix and its inverse
          DO col = 1, n
            DO row = 1, n
              c(row, col) = 0._r32
            ENDDO
          ENDDO

          ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
          ! ----------------------------------------------- build covariance matrix  ---------------------------------------------
          DO col = 1, n
            DO row = col, n

              ! inter-station distance (must be in meter)
              eta = HYPOT(input%receiver(row)%x - input%receiver(col)%x, input%receiver(row)%y - input%receiver(col)%y)

              ! cross spectral density
              c(row, col) =  2._r32 * lagged_coherency(freq, eta)

            ENDDO
          ENDDO

          ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
          ! ----------------------------------------------- invert covariance matrix  --------------------------------------------
          CALL cholsk(ok, c)

          CALL rng(rnorm)        !< draw two numbers from normal distribution

          fc_a(fr, :) = MATMUL(c, rnorm)

          CALL rng(rnorm)

          fc_b(fr, :) = MATMUL(c, rnorm)

        ENDDO       !< end loop over frequency

        ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
        ! ----------------------------------------- build spectrum, return timeseries  -------------------------------------------

        DO rcvr = 1, n

          DO fr = 1, npts
            z(fr) = CMPLX(fc_a(fr, rcvr), fc_b(fr, rcvr)) * (2*npts) * timeseries%sp%dt
          ENDDO

          IF (comp .eq. 1) THEN
            CALL ifft(timeseries%sp%x(:, rcvr), z)
          ELSEIF (comp .eq. 2) THEN
            CALL ifft(timeseries%sp%y(:, rcvr), z)
          ELSEIF (comp .eq. 3) THEN
            CALL ifft(timeseries%sp%z(:, rcvr), z)
          ENDIF

        ENDDO     !< end loop over receivers

      ENDDO  !< end loop over components

      npts = SIZE(timeseries%sp%time)

      CALL destroy_fftw_plan([npts])

!#ifdef DEBUG
      OPEN(newunit = lu, file = 'noise_dbg.txt', status = 'unknown', form = 'formatted', access = 'sequential', action = 'write', &
           iostat = ok)
      DO i = 1, SIZE(timeseries%sp%time)
        WRITE(lu, *) timeseries%sp%time(i), timeseries%sp%x(i, :)
      ENDDO
      CLOSE(lu)
!#endif

#ifdef ERROR_TRAP
      IF (error .ne. 0) THEN
        IF (rank .eq. 0) CALL report_error('generate_noise - ERROR: matrix is singular')
        ok = error
      ENDIF
#endif

    END SUBROUTINE generate_noise_cdra

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE cholsk(ok, a)

      INTEGER(i32),                      INTENT(OUT)   :: ok
      REAL(r32),    DIMENSION(:,:),      INTENT(INOUT) :: a
      INTEGER(i32)                                     :: n, row, col

      !-----------------------------------------------------------------------------------------------------------------------------

      n = SIZE(a, 1)

      CALL spotrf('L', n, a, n, ok)

      DO col = 1, n
        DO row = 1, col - 1
          a(row, col) = 0._r32
        ENDDO
      ENDDO

    END SUBROUTINE cholsk

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *


END MODULE m_noise
