MODULE m_filter

  ! description:
  !   collection of subprograms to filter time-series according to iir or fir filters.
  !   available iir filters are butterworth, chebychev type i and ii and elliptic.
  !   available fir filters are [...]
  !   accepted filter types are lowpass, highpass, bandpass and bandstop. note that bandpass and bandstop filters have an order twice
  !   higher than specified in input.
  !
  !   the general CALL sequence is:
  !     a) setiir / setfir
  !     b) iir / fir
  !
  !   filter coefficients and OPTIONAL initial conditions for zero-phase filtering are stored in global arrays private to the module.
  !   modules required by "filter" are listed below this header.
  !
  ! author:
  !   w. imperatori, walter.imperatori@sed.ethz.ch
  !
  ! references:
  !   parks & burrus, "digital filter design", 1987
  !   gnu octave, www.gnu.org/software/octave/
  !
  ! version:
  !   0.1 -> initial version, only iir filters implemented
  !

  !$ USE             :: omp_lib
  USE, NON_INTRINSIC :: m_precisions
  USE, NON_INTRINSIC :: m_llsq_r64        !< double precision
  USE, NON_INTRINSIC :: m_fft_real        !< double precision

  IMPLICIT none

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  PRIVATE

  PUBLIC :: make_iir_plan, destroy_iir_plan, iir, filter_error, impz, freqz, get_iir_coefficients, set_iir_coefficients

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  ! set precisions of input/output arguments at compile-time
#ifdef DOUBLE_PREC
  INTEGER, PARAMETER :: r__ = r64
#else
  INTEGER, PARAMETER :: r__ = r32
#endif

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  REAL(r64),    PARAMETER :: PI = 3.14159265358979323846_r64
  REAL(r64),    PARAMETER :: BIG = HUGE(0._r64)
  COMPLEX(r64), PARAMETER :: EMPTY = CMPLX(BIG, BIG, r64)           !< constant for empty complex vector


  ! transfer function coefficients of filter ("b" is numerator, "a" is denominator)
  REAL(r64), ALLOCATABLE, DIMENSION(:) :: iira, iirb

  ! initial conditions for zero-phase filtering
  REAL(r64), ALLOCATABLE, DIMENSION(:) :: zi

  !$omp threadprivate (iira, iirb, zi)

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

  CONTAINS

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE make_iir_plan(ok, design, dt, fc, type, order, ripple, zphase, fir)

      ! Purpose:
      !   To set up the coefficients of an IIR filter, including initial conditions to avoid transients at startup for zero-phase
      !   filtering. If requested, the filter impulse response is computed and stored to disk ("filter_impulse_response.txt"). Filter
      !   coefficients and (optional) initial conditions are computed and stored in private arrays accessible by routine "iir" for
      !   the actual filtering.
      !
      !   Specify filter designs with "design". Accepted values are "butter" (Butterworth), "cheby1" or "cheby2" (Chebychev type i
      !   or ii) and "ellip" (Elliptic).
      !
      !   "dt" is the time-step of the time-series to be filtered.
      !
      !   "fc" is a one/two-elements vector specifying corner frequencies (Hz) in the range (0, Nyquist).
      !
      !   Possible filter types ("type") are "low" (low-pass), "high" (high-pass), "pass" (band-pass) and "stop" (bandstop).
      !
      !   "order" expresses the number of poles and must be larger than zero. This value affects the minimum length of the sequence
      !   to be filtered if zero-phase filtering is selected.
      !
      !   "ripple" is a one/two-elements vector ("ripple") specifies the desired ripple for Chebychev and Elliptic filters. Units
      !   are assumed in dB if larger than 1 and linear otherwise.
      !
      !   For zero-phase (two-pass) filters, set "zphase" to true.
      !
      !   The impulse response of the filter is saved to disk if "fir" is set to true.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   02/09/20                  original version
      !

      INTEGER(i32),                                      INTENT(OUT) :: ok
      CHARACTER(*),                                      INTENT(IN)  :: design            !< filter design
      REAL(r__),                               OPTIONAL, INTENT(IN)  :: dt                !< time-step
      REAL(r__),                 DIMENSION(:), OPTIONAL, INTENT(IN)  :: fc                !< corner frequency
      CHARACTER(*),                            OPTIONAL, INTENT(IN)  :: type              !< filter type
      INTEGER(i32),                            OPTIONAL, INTENT(IN)  :: order             !< filter order
      REAL(r__),                 DIMENSION(:), OPTIONAL, INTENT(IN)  :: ripple            !< ripple for chebychev/elliptic filters
      LOGICAL,                                 OPTIONAL, INTENT(IN)  :: zphase            !< single- or double-pass (zerophase)
      LOGICAL,                                 OPTIONAL, INTENT(IN)  :: fir               !< output filter impulse response
      INTEGER(i32)                                                   :: i, n
      INTEGER(i32)                                                   :: fd, ft
      REAL(r64)                                                      :: fs
      REAL(r64)                                                      :: gain
      REAL(r64),    ALLOCATABLE, DIMENSION(:)                        :: rip
      REAL(r64),    ALLOCATABLE, DIMENSION(:)                        :: w
      REAL(r__),    ALLOCATABLE, DIMENSION(:)                        :: freq, amp, phase
      COMPLEX(r64), ALLOCATABLE, DIMENSION(:)                        :: zero, pole

      !-----------------------------------------------------------------------------------------------------------------------------

      ok = 0

      ! raise error if not all mandatory arguments are present
      IF (.not.PRESENT(dt))    ok = 1
      IF (.not.PRESENT(type))  ok = 2
      IF (.not.PRESENT(fc))    ok = 3
      IF (.not.PRESENT(order)) ok = 4

      IF (ok .ne. 0) RETURN

      ! check filter design
      fd = 0

      IF (design .eq. 'butter') fd = 1
      IF (design .eq. 'cheby1') fd = 2
      IF (design .eq. 'cheby2') fd = 3
      IF (design .eq. 'ellip')  fd = 4

      ! error: filter type not supported
      IF (fd .eq. 0) ok = 5

      IF (ok .ne. 0) RETURN

      ! check ripple PARAMETER
      IF ( (fd .ge. 2) .and. (fd .le. 4) .and. .not.PRESENT(ripple) ) ok = 6

      IF (PRESENT(ripple)) THEN
        IF ( (fd .ge. 2) .and. (fd .le. 3) .and. (SIZE(ripple) .ne. 1) ) ok = 7
        IF ( (fd .eq. 4) .and. (SIZE(ripple) .ne. 2) ) ok = 8
      ENDIF

      IF (ok .ne. 0) RETURN

      ! filter order must be a positive number
      IF (order .le. 0) ok = 9

      ! check filter type
      ft = 0

      IF (type .eq. 'low')  ft = 1
      IF (type .eq. 'high') ft = 2
      IF (type .eq. 'pass') ft = 3
      IF (type .eq. 'stop') ft = 4

      IF (ft .eq. 0) ok = 10

      IF (ok .ne. 0) RETURN

      n = SIZE(fc)

      ALLOCATE(w(n))

      ! samping rate
      fs = 1._r64 / dt

      ! normalised frequency
      w = fc / (fs * 0.5_r64)

      ! check size of corner frequency vector
      IF ( (ft .ge. 3) .and. (n .ne. 2) ) ok = 11    ! text = 'for bandpass/bandstop, two corner frequencies are needed'
      IF ( (ft .le. 2) .and. (n .ne. 1) ) ok = 12    ! text = 'for lowpass/highpass, one corner frequency is needed'

      ! check that corner frequency is above zero and smaller than nyquist frequency
      IF (ANY(fc .le. 0._r__) .or. ANY(w .gt. 1._r64)) ok = 13

      IF (ok .ne. 0) RETURN

      IF (fd .eq. 4) THEN
        ALLOCATE(rip(2))
      ELSE
        ALLOCATE(rip(1))
      ENDIF

      rip(:) = 0._r64

      IF (PRESENT(ripple)) rip = ripple

      ! make sure that two-components ripple is either in linear units or decibels
      IF (SIZE(rip) .eq. 2) THEN
        IF ( ((rip(1) .lt. 1._r64) .and. (rip(2) .gt. 1._r64)) .or. ((rip(1) .gt. 1._r64) .and. (rip(2) .lt. 1._r64)) ) ok = 14
      ENDIF

      IF (ok .ne. 0) RETURN

      ! "ripple" is assumed to be in linear units if smaller than unit, in decibel otherwise.
      IF (rip(1) .lt. 1._r__) THEN

        ! transform ripple from percent to decibels (eq. 7.35 and 7.52 of parks & burrus)
        ! e.g. for cheby1, 0.05 (a 5% ripple, oscillating between 1 and 0.95) is equivalent to 0.44553 db
        ! e.g. for cheby2, 0.05 (a 5% ripple, oscillating between 0 and 0.05) is equivalent to 26.0206 db
        IF (fd .eq. 2) rip(1) = -20._r64 * log10(1._r64 - rip(1))
        IF (fd .eq. 4) rip(1) = -20._r64 * log10(1._r64 - rip(1))
        IF (fd .eq. 3) rip(1) = -20._r64 * log10(rip(1))
        IF (fd .eq. 4) rip(2) = -20._r64 * log10(rip(2))

      ENDIF

      ! analog, prewarped frequency
      w = TAN(PI * w / 2._r64)

      ! prepare arrays for zeros and poles large enough for bandpass/bandstop doubling
      ALLOCATE(zero(order * 2), pole(order * 2))

      zero = EMPTY
      pole = EMPTY

      ! compute poles, zeroes and gain
      IF (fd .eq. 1) CALL butter(order, zero, pole, gain)
      IF (fd .eq. 2) CALL cheby1(order, rip(1), zero, pole, gain)
      IF (fd .eq. 3) CALL cheby2(order, rip(1), zero, pole, gain)
      IF (fd .eq. 4) CALL ellip(order, rip, zero, pole, gain)

      ! CALL report('zeroes', zero)
      ! CALL report('poles', pole)
      ! print*, 'gain ', REAL(gain, r32)

      ! transform to filter of desired cutoff frequency
      CALL xtrans(zero, pole, gain, w, ft, ok)

      IF (ok .ne. 0) THEN
        ok = 15
        RETURN
      ENDIF


      ! print*, ''
      ! print*, 'after xtrans'
      ! CALL report('zeroes', zero)
      ! CALL report('poles', pole)
      ! print*, 'gain ', REAL(gain, r32)

      ! analog to digital conversion
      CALL bilinear(zero, pole, gain)

      ! print*, ''
      ! print*, 'after bilinear'
      ! CALL report('zeroes', zero)
      ! CALL report('poles', pole)
      ! print*, 'gain ', REAL(gain, r32)

      n = COUNT(pole .ne. EMPTY) + 1

      !$omp parallel default(shared)
      ALLOCATE(iira(n), iirb(n))

      ! now convert to filter coefficients
      iirb = gain * poly(zero(1:n-1))             !< numerator
      iira = poly(pole(1:n-1))                    !< denominator
      !$omp end parallel

      ! clean up
      DEALLOCATE(pole, zero)
      DEALLOCATE(w, rip)

      ! compute initial conditions
      IF (PRESENT(zphase)) THEN
        IF (zphase .eqv. .true.) THEN

          ok = 0

          !$omp parallel default(shared) reduction(max: ok)
          ALLOCATE(zi(n - 1))

          zi = initial_conditions(iirb, iira, ok)
          !$omp end parallel

          IF (ok .ne. 0) RETURN

        ENDIF
      ENDIF

      ! write filter impulse response to disk if requested
      IF (PRESENT(fir)) THEN
        IF (fir .eqv. .true.) THEN

          ! evaluate filter impulse response
          CALL freqz(freq, amp, phase, dt = dt)

          OPEN(1, file = 'filter_impulse_response.txt', status = 'unknown', form = 'formatted', iostat = ok)

          IF (ok .ne. 0) THEN
            ok = 17
            RETURN
          ENDIF

          DO i = 1, SIZE(freq)
            WRITE(1, *) freq(i), amp(i), phase(i)
          ENDDO

          CLOSE(1)

          ! clean up
          DEALLOCATE(freq, amp, phase)

        ENDIF
      ENDIF

    END SUBROUTINE make_iir_plan

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE destroy_iir_plan()

      ! Purpose:
      !   To deallocate filter coefficients computed in "make_iir_plan"
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   02/09/20                  original version
      !

      !-----------------------------------------------------------------------------------------------------------------------------

      ! always check allocation state to avoid run-time problems if one calls "destroy_iir_plan" before "make_iir_plan"
      !$omp parallel default(shared)
      IF (ALLOCATED(iira)) DEALLOCATE(iira)
      IF (ALLOCATED(iirb)) DEALLOCATE(iirb)
      IF (ALLOCATED(zi))   DEALLOCATE(zi)
      !$omp end parallel

    END SUBROUTINE destroy_iir_plan

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    FUNCTION iir(fun, ok)

      ! Purpose:
      !   to filter a sequence "fun" with an IIR filter whose coefficients and (optional) initial conditions (to avoid transients at
      !   startup) were computed in "make_iir_plan". This routine dynamically allocates a vector for better versatility. Calculations
      !   are always in double-precision.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   02/09/20                  original version
      !

      REAL(r__),                 DIMENSION(:),        INTENT(IN)  :: fun                   !< input time-series
      INTEGER(i32),                                   INTENT(OUT) :: ok                    !< error flag
      INTEGER(i32)                                                :: i, l, n, nt
      INTEGER(i32)                                                :: ierr
      REAL(r__),                 DIMENSION(SIZE(fun))             :: iir
      REAL(r64),    ALLOCATABLE, DIMENSION(:)                     :: x

      !-----------------------------------------------------------------------------------------------------------------------------

      ierr = 0

      ! raise error if filter coefficients are not available
      IF (.not.ALLOCATED(iirb) .or. .not.ALLOCATED(iira)) ierr = 18

      ! length of time-series
      n = SIZE(fun)

      ! initialise for one-pass
      l = 1

      ! for zero-phase (two-pass) filtering, set "l" to length of filter coefficients
      IF (ALLOCATED(zi)) l = SIZE(iirb)

      ! set length of transient
      nt = 3 * (l - 1)

      ! raise error if input sequence isn't long enough (relevant for zero-phase filtering)
      IF (n .le. nt) ierr = 19

      ! return immediately if input is not correct
      IF (ierr .ne. 0) THEN
        ok = MAX(ok, ierr)
        RETURN
      ENDIF

      ALLOCATE(x(n + 2*nt))

      ! copy time-series into double precision array, properly concatenated for zero-phase filtering. First and last loop are skipped
      ! (because "nt = 0") for one-pass only filtering.
      DO i = 1, nt
        x(i) = 2._r64 * REAL(fun(1), r64) - REAL(fun(nt + 2 - i), r64)
      ENDDO

      DO i = 1, n
        x(nt + i) = REAL(fun(i), r64)
      ENDDO

      DO i = 1, nt
        x(n + nt + i) = 2._r64 * REAL(fun(n), r64) - REAL(fun(n - i), r64)
      ENDDO

      ! two-pass filtering
      IF (l .ne. 1) THEN

        CALL xfilter(x, iirb, iira, zi)                !< first pass
        x = fliplr(x)                                  !< reverse time-series

        CALL xfilter(x, iirb, iira, zi)                !< second pass
        x = fliplr(x)                                  !< reverse once more

      ! single-pass filtering
      ELSE

        CALL xfilter(x, iirb, iira)

      ENDIF

      ! return central part
      DO i = 1, n
        iir(i) = REAL(x(nt + i), r__)
      ENDDO

      ! DEALLOCATE(x)

    END FUNCTION iir

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    PURE SUBROUTINE xfilter(x, fb, fa, z)

      ! Purpose:
      !   to filter a sequence "x" according to filter coefficients "fb" and "fa" (the latter must be empty for FIR filters) and
      !   optional initial conditions "z" to avoid transients at startup.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   02/09/20                  original version
      !

      REAL(r64),   DIMENSION(:),                INTENT(INOUT) :: x                !< time-series
      REAL(r64),   DIMENSION(:),                INTENT(IN)    :: fb               !< filter coefficients (numerator)
      REAL(r64),   DIMENSION(:),      OPTIONAL, INTENT(IN)    :: fa               !< filter coefficients (denominator)
      REAL(r64),   DIMENSION(:),      OPTIONAL, INTENT(IN)    :: z                !< initial conditions of system
      INTEGER(i32)                                            :: i, j, n
      REAL(r64),   DIMENSION(SIZE(x))                         :: y
      REAL(r64),   DIMENSION(SIZE(fb))                        :: b
      REAL(r64),   DIMENSION(SIZE(fb))                        :: a
      REAL(r64),   DIMENSION(SIZE(fb))                        :: zi

      !-----------------------------------------------------------------------------------------------------------------------------

      ! number of coefficients
      n = SIZE(fb)

      a(:)  = 0._r64
      zi(:) = 0._r64

      ! normalisation
      IF (PRESENT(fa)) THEN

        DO i = 1, n
          b(i) = fb(i) / fa(1)
        ENDDO

        DO i = 1, SIZE(fa)
          a(i) = fa(i) / fa(1)
        ENDDO

      ELSE

        DO i = 1, n
          b(i) = fb(i)
        ENDDO

        a(1) = 1._r64

      ENDIF

      ! assign existing initial conditions
      IF (PRESENT(z)) THEN
        DO i = 1, n - 1
          zi(i) = z(i) * x(1)
        ENDDO
      ENDIF

      DO j = 1, SIZE(x)

        y(j) = b(1) * x(j) + zi(1)

        DO i = 2, n
          zi(i - 1) = b(i) * x(j) + zi(i) - a(i) * y(j)
        ENDDO

        x(j) = y(j)

      ENDDO

      !z = zi(1:SIZE(z)-1)

    END SUBROUTINE xfilter

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    FUNCTION fliplr(x) RESULT(y)

      ! Purpose:
      !   to return a vector in reversed order, i.e. from its last to its first element.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   02/09/20                  original version
      !

      REAL(r64),   DIMENSION(:),      INTENT(INOUT) :: x
      INTEGER(i32)                                  :: i, n
      REAL(r64),   DIMENSION(SIZE(x))               :: y

      !-----------------------------------------------------------------------------------------------------------------------------

      n = SIZE(x)

      DO i = 1, n
        y(i) = x(n - i + 1)
      ENDDO

    END FUNCTION fliplr

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    ! SUBROUTINE report(message, x)
    !
    !   ! print filter coefficients in single precision to screen
    !
    !   CHARACTER(len=*),               INTENT(IN) :: message
    !   COMPLEX(r64),     DIMENSION(:), INTENT(IN) :: x
    !   INTEGER(i32)                               :: i
    !
    !   !-----------------------------------------------------------------------------------------------------------------------------
    !
    !   print*, message
    !
    !   DO i = 1, COUNT(x .ne. empty)
    !     WRITE(stdout, *), REAL(x(i), r32), REAL(AIMAG(x(i)), r32)
    !   ENDDO
    !
    ! END SUBROUTINE report

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE butter(order, zero, pole, gain)

      ! Purpose:
      !   to compute poles&zeroes of a prototype Butterworth filter.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   02/09/20                  original version
      !

      INTEGER(i32),               INTENT(IN)  :: order
      COMPLEX(r64), DIMENSION(:), INTENT(OUT) :: zero, pole
      REAL(r64),                  INTENT(OUT) :: gain
      INTEGER(i32)                            :: i
      REAL(r64)                               :: arg, n

      !-----------------------------------------------------------------------------------------------------------------------------

      n = order * 2._r64

      DO i = 1, order
        arg     = (2._r64 * (i - 1) + 1) * PI / n
        pole(i) = CMPLX(-SIN(arg), COS(arg), r64)
        zero(i) = EMPTY
      ENDDO

      IF (MOD(order, 2) .ne. 0) pole((order + 1) / 2) = -1._r64

      gain = 1._r64

    END SUBROUTINE butter

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE cheby1(order, rip, zero, pole, gain)

      ! Purpose:
      !   to compute poles&zeroes of a prototype Chebychev filter, type i).
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   02/09/20                  original version
      !

      INTEGER(i32),               INTENT(IN)  :: order
      REAL(r64),                  INTENT(IN)  :: rip
      COMPLEX(r64), DIMENSION(:), INTENT(OUT) :: zero, pole
      REAL(r64),                  INTENT(OUT) :: gain
      INTEGER(i32)                            :: i
      REAL(r64)                               :: arg, n
      REAL(r64)                               :: x, y, eps, v0

      !-----------------------------------------------------------------------------------------------------------------------------

      n = order * 2._r64

      eps = SQRT(10._r64**(0.1_r64 * rip) - 1._r64)

      v0 = asinh(1._r64 / eps) / REAL(order, r64)

      x = -sinh(v0)
      y = -cosh(v0)

      DO i = 1, order
        arg     = (2._r64 * (i - 1) + 1) * PI / n
        pole(i) = CMPLX(x * SIN(arg), y * COS(arg), r64)
        zero(i) = EMPTY
      ENDDO

      gain = PRODUCT(-pole(1:order))

      ! correct gain
      IF (MOD(order, 2) .eq. 0) gain = gain / 10._r64**(0.05_r64 * rip)

    END SUBROUTINE cheby1

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE cheby2(order, rip, zero, pole, gain)

      ! Purpose:
      !   to compute poles&zeroes of a prototype Chebychev filter, type ii).
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   02/09/20                  original version
      !

      INTEGER(i32),               INTENT(IN)  :: order
      REAL(r64),                  INTENT(IN)  :: rip
      COMPLEX(r64), DIMENSION(:), INTENT(OUT) :: zero, pole
      REAL(r64),                  INTENT(OUT) :: gain
      INTEGER(i32)                            :: i, k, l
      REAL(r64)                               :: arg, n
      REAL(r64)                               :: x, y, eps, v0
      REAL(r64)                               :: re, im, d

      !-----------------------------------------------------------------------------------------------------------------------------

      n = order * 2._r64

      eps = 1._r64 / SQRT(10._r64**(0.1_r64 * rip) - 1._r64)

      v0 = asinh(1._r64 / eps) / REAL(order, r64)

      x = -sinh(v0)
      y = -cosh(v0)

      k = 0

      DO i = 1, order

        l       = 2 * (i - 1) + 1
        arg     = REAL(l, r64) * PI / n
        im      = COS(arg)

        ! skip zeroes at infinity (e.g. when arg = PI/2, leading to COS(arg) = 0)
        IF (l .ne. order) THEN
          k       = k + 1
          zero(k) = CMPLX(0._r64, 1._r64 / COS(arg), r64)
        ENDIF

        re      = x * SIN(arg)
        im      = y * im
        d       = re**2 + im**2
        pole(i) = CMPLX(re / d, im / d, r64)

      ENDDO

      gain = PRODUCT(-pole(1:order)) / PRODUCT(-zero(1:k))

    END SUBROUTINE cheby2

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE ellip(order, rip, zero, pole, gain)

      ! Purpose:
      !   to compute poles&zeroes of a prototype Elliptic filter. "rip" must contain the passband and stopband ripple in dB.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   02/09/20                  original version
      !

      INTEGER(i32),                         INTENT(IN)  :: order
      REAL(r64),    DIMENSION(2),           INTENT(IN)  :: rip
      COMPLEX(r64), DIMENSION(:),           INTENT(OUT) :: zero, pole
      REAL(r64),                            INTENT(OUT) :: gain
      COMPLEX(r64)                                      :: u, v0
      COMPLEX(r64),               PARAMETER             :: j = CMPLX(0._r64, 1._r64, r64)
      INTEGER(i32)                                      :: i, i0, r
      REAL(r64)                                         :: eps
      REAL(r64)                                         :: k, k1
      REAL(r64),                  PARAMETER             :: tol = EPSILON(1._r64)

      !-----------------------------------------------------------------------------------------------------------------------------

      eps = SQRT(10._r64**(0.1_r64 * rip(1)) - 1._r64)                   !< eq. 7.90

      k1  = eps / SQRT(10._r64**(0.1_r64 * rip(2)) - 1._r64)             !< eq. 7.91

      ! compute modulus "k"
      k = ellip2k(k1, order, tol)

      r = MOD(order, 2)

      i0 = 1

      IF (r .eq. 1) i0 = 2

      ! compute "v0" term
      u  = j / eps
      v0 = -j * arcsn(u, k1, tol)                                        !< eq. 7.81
      v0 = v0 / REAL(order, r64)

      DO i = i0, order - 1, 2

        u = CMPLX(REAL(i, r64), 0._r64, r64)
        u = u / REAL(order, r64)

        ! pair of conjugate zeros
        zero(i - r)     = j / (k * sn(u, k, tol))
        zero(i - r + 1) = CONJG(zero(i - r))

        u = u + j * v0

        ! conjugate pair of poles
        pole(i - r)     = sn(u, k, tol) * j
        pole(i - r + 1) = CONJG(pole(i - r))

      ENDDO

      ! add pole for i = 0 (order odd)
      IF (r .eq. 1) THEN

        u           = j * v0
        pole(order) = j * sn(u, k, tol)

      ENDIF

      gain = 10**(-rip(1) * 0.05_r64)
      gain = gain**(1 - r)
      gain = ABS(gain * PRODUCT(pole(1:order)) / PRODUCT(zero(1:order - r)))

    END SUBROUTINE ellip

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE landen(k, tol, v)

      ! Purpose:
      !   to compute a Landen descending transformation to evalute Elliptic integrals
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   02/09/20                  original version
      !

      REAL(r64),                              INTENT(IN)    :: k
      REAL(r64),                              INTENT(IN)    :: tol
      REAL(r64),   ALLOCATABLE, DIMENSION(:), INTENT(INOUT) :: v
      INTEGER(i32)                                          :: i, n                    !< counters
      REAL(r64)                                             :: x
      REAL(r64),   ALLOCATABLE, DIMENSION(:)                :: temp

      !-----------------------------------------------------------------------------------------------------------------------------

      ALLOCATE(v(1))

      IF ( (k .eq. 0._r64) .or. (k .eq. 1._r64) ) THEN

        v = k

      ELSE

        x = k

        n = 0

        DO WHILE (x .gt. tol)

          x = x**2 / (1._r64 + SQRT(1._r64 - x**2))**2

          ALLOCATE(temp(n + 1))

          DO i = 1, n
            temp(i) = v(i)
          ENDDO

          CALL move_alloc(temp, v)

          v(n + 1) = x

          n = n + 1

        ENDDO

      ENDIF

    END SUBROUTINE landen

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE ceifk(k, tol, ek, ekc)

      ! Purpose:
      !   to compute the complete elliptic integral of the first kind for elliptic module "k" and complementary elliptic module.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   02/09/20                  original version
      !

      REAL(r64),                                     INTENT(IN)  :: k
      REAL(r64),                                     INTENT(IN)  :: tol
      REAL(r64),                                     INTENT(OUT) :: ek                !< solution for modulus "k"
      REAL(r64),                                     INTENT(OUT) :: ekc               !< solution for complementary modulus
      REAL(r64)                                                  :: kc
      REAL(r64),                           PARAMETER             :: kmin = 1.e-06_r64
      REAL(r64),                           PARAMETER             :: kmax = SQRT(1._r64 - kmin**2)
      REAL(r64)                                                  :: l
      REAL(r64), ALLOCATABLE, DIMENSION(:)                       :: v

      !-----------------------------------------------------------------------------------------------------------------------------

      ! complementary module
      kc = SQRT(1._r64 - k**2)

      ! exactly one
      IF (k .eq. 1._r64) THEN

        ek = BIG

      ! near one
      ELSEIF (k .gt. kmax) THEN

        l = -LOG(kc / 4._r64)

        ek = l + (l - 1) * kc**2 / 4._r64

      ELSE

        CALL landen(k, tol, v)

        ek = PRODUCT(1._r64 + v) * PI / 2._r64

        DEALLOCATE(v)

      ENDIF

      ! exactly zero
      IF (k .eq. 0._r64) THEN

        ekc = BIG

      ! near one
      ELSEIF (k .gt. kmin) THEN

        l = -LOG(k / 4._r64)

        ekc = l + (l - 1) * k**2 / 4._r64

      ELSE

        CALL landen(kc, tol, v)

        ekc = PRODUCT(1._r64 + v) * PI / 2._r64

        DEALLOCATE(v)

      ENDIF

    END SUBROUTINE ceifk

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    FUNCTION ellip2k(k1, order, tol) RESULT(k)

      ! Purpose:
      !   to return the module k.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   02/09/20                  original version
      !

      REAL(r64),    INTENT(IN) :: k1
      INTEGER(i32), INTENT(IN) :: order
      REAL(r64),    INTENT(IN) :: tol
      COMPLEX(r64)             :: u
      INTEGER(i32)             :: i, l
      REAL(r64)                :: k, kc, k1c
      REAL(r64)                :: ek1, ek1c
      REAL(r64)                :: q, s, n, d
      REAL(r64),    PARAMETER  :: kmin = 1.e-06_r64

      !-----------------------------------------------------------------------------------------------------------------------------

      l = order / 2

      ! nome expansion for small "k1"
      IF (k1 .lt. kmin) THEN

        CALL ceifk(k1, tol, ek1, ek1c)

        q = EXP(-PI * ek1c / ek1)                !< nome
        q = q**(1._r64 / REAL(order, r64))

        n = 0._r64
        d = 0._r64

        ! expansion based on 7 terms
        DO i = 1, 7
          n = n + q**(i**2 + i)
          d = d + q**(i**2)
        ENDDO

        k = 4._r64 * SQRT(q) * ((1._r64 + n) / (1._r64 + 2._r64 * d))**2

      ! exact solution
      ELSE

        k1c = SQRT(1._r64 - k1**2)                           !< complementary module of "k1"

        s = 1._r64

        DO i = 1, l
          u = CMPLX(2._r64 * i - 1._r64, 0._r64, r64)         !< argument of FUNCTION "sn"
          u = u / REAL(order, r64)
          s = s * REAL(sn(u, k1c, tol), r64)
        ENDDO

        kc = k1c**order * s**4                               !< complementary module of "k"

        k = SQRT(1._r64 - kc**2)

      ENDIF

    END FUNCTION ellip2k

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    FUNCTION sn(u, k, tol)

      ! Purpose:
      !   to compute the Jacobi elliptic function sn.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   02/09/20                  original version
      !

      COMPLEX(r64),                           INTENT(IN) :: u
      REAL(r64),                              INTENT(IN) :: k
      REAL(r64),                              INTENT(IN) :: tol
      INTEGER(i32)                                       :: i
      COMPLEX(r64)                                       :: sn
      REAL(r64),    ALLOCATABLE, DIMENSION(:)            :: v

      !-----------------------------------------------------------------------------------------------------------------------------

      CALL landen(k, tol, v)

      sn = SIN(u * PI / 2._r64)

      ! reverse loop order to perform an ascending landen transformation
      DO i = SIZE(v), 1, -1
        sn = (1._r64 + v(i)) * sn / (1._r64 + v(i) * sn**2)
      ENDDO

    END FUNCTION sn

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    FUNCTION arcsn(w, k, tol)

      ! Purpose:
      !   to compute the Jacobi elliptic function arcsn.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   02/09/20                  original version
      !

      COMPLEX(r64),                           INTENT(IN) :: w
      REAL(r64),                              INTENT(IN) :: k
      REAL(r64),                              INTENT(IN) :: tol
      COMPLEX(r64)                                       :: arcsn
      COMPLEX(r64)                                       :: x, u
      INTEGER(i32)                                       :: i
      REAL(r64)                                          :: v1
      REAL(r64)                                          :: si, sr
      REAL(r64)                                          :: ratio, ek, ekc
      REAL(r64),    ALLOCATABLE, DIMENSION(:)            :: v

      !-----------------------------------------------------------------------------------------------------------------------------

      CALL landen(k, tol, v)

      x = w

      DO i = 1, SIZE(v)

        IF (i .eq. 1) THEN
          v1 = k
        ELSE
          v1 = v(i - 1)
        ENDIF

        x = x / (1._r64 + SQRT(1._r64 - x**2 * v1**2)) * 2._r64 / (1._r64 + v(i))

      ENDDO

      u = 2._r64 / PI * acos(x)

      IF (REAL(x, r64) .eq. 1._r64) u = CMPLX(0._r64, 0._r64)

      CALL ceifk(k, tol, ek, ekc)

      ratio = ekc / ek

      sr = MOD(REAL(u, r64), 4._r64)

      IF (ABS(sr) .gt. 2._r64) sr = sr - 4._r64 * SIGN(1._r64, sr)

      si = MOD(AIMAG(u), ratio * 2._r64)

      IF (ABS(si) .gt. ratio) si = si - ratio * 2._r64 * SIGN(1._r64, si)

      arcsn = 1._r64 - CMPLX(sr, si)

    END FUNCTION arcsn

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE xtrans(zero, pole, gain, w, ft, ok)

      ! Purpose:
      !   to compute the analog low-pass to digital low-pass, high-pass, bandpass and bandstop transform. In general filters have
      !   same number of poles and zeroes, except Elliptic filters for which the number of poles may be larger than the number of
      !   zeroes. Note that the effective number of poles and zeroes is increased for bandpass and bandstop filters.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   02/09/20                  original version
      !

      COMPLEX(r64), DIMENSION(:), INTENT(INOUT) :: zero, pole            !< poles and zeroes
      REAL(r64),                  INTENT(INOUT) :: gain                  !< gain
      REAL(r64),    DIMENSION(:), INTENT(IN)    :: w                     !< normalised corner frequency
      INTEGER(i32),               INTENT(IN)    :: ft                    !< filter type (1=lp, 2=hp, 3=bp, 4=bs)
      INTEGER(i32),               INTENT(OUT)   :: ok                    !< error flag
      COMPLEX(r64)                              :: b, x
      INTEGER(i32)                              :: i, k
      INTEGER(i32)                              :: np, nz
      LOGICAL                                   :: zempty, pempty
      REAL(r64)                                 :: fl, fh

      !-----------------------------------------------------------------------------------------------------------------------------

      ok = 0

      ! number of poles and zeroes
      np = COUNT(pole .ne. EMPTY)
      nz = COUNT(zero .ne. EMPTY)

      ! determine if vector "zero" and "pole" are completely empty
      zempty = nz .eq. 0
      pempty = np .eq. 0

      ! we need to have at least as many poles as zeroes
      IF (nz > np) THEN
        ok = 15
        RETURN
      ENDIF

      ! deal with low-pass and high-pass filters
      IF (SIZE(w) .eq. 1) THEN

        ! low-pass filter
        IF (ft .eq. 1) THEN

          gain = gain * (1._r64 / w(1))**(nz - np)

          DO i = 1, np
            pole(i) = w(1) * pole(i)
          ENDDO
          DO i = 1, nz
            zero(i) = w(1) * zero(i)
          ENDDO

        ! high-pass filter
        ELSEIF (ft .eq. 2) THEN

          IF (zempty) THEN
            gain = gain * REAL(1._r64 / PRODUCT(-pole(1:np)), r64)
          ELSEIF (pempty) THEN
            gain = gain * REAL(PRODUCT(-zero(1:nz)), r64)
          ELSE
            gain = gain * REAL(PRODUCT(-zero(1:nz)) / PRODUCT(-pole(1:np)), r64)
          ENDIF

          DO i = 1, np
            pole(i) = w(1) / pole(i)
          ENDDO

          IF (zempty) THEN

            DO i = 1, np
              zero(i) = CMPLX(0._r64, 0._r64, r64)
            ENDDO

          ELSE

            DO i = 1, nz
              zero(i) = w(1) / zero(i)
            ENDDO

            ! add zeroes such that nz equals np
            IF (np .gt. nz) THEN
              DO i = nz + 1, np
                zero(i) = CMPLX(0._r64, 0._r64, r64)
              ENDDO
            ENDIF

          ENDIF

        ENDIF

      ! deal with bandpass and bandstop filters
      ELSE

        fl = w(1)     !< lower corner frequency
        fh = w(2)     !< higher corner frequency

        ! pass-band filter
        IF (ft .eq. 3) THEN

          gain = gain * (1._r64 / (fh - fl))**(nz - np)

          DO i = 1, np
            b            = pole(i) * (fh - fl) / 2._r64
            x            = SQRT(b**2 - fh*fl)
            pole(i)      = b + x
            pole(np + i) = b - x
          ENDDO

          IF (zempty) THEN

            DO i = 1, np
              zero(i) = CMPLX(0._r64, 0._r64, r64)
            ENDDO

          ELSE

            DO i = 1, nz
              b            = zero(i) * (fh - fl) / 2._r64
              x            = SQRT(b**2 - fh*fl)
              zero(i)      = b + x
              zero(nz + i) = b - x
            ENDDO

            ! add zeroes such that nz equals np
            IF (np .gt. nz) THEN
              DO i = 1, np - nz
                zero(i + 2*nz) = CMPLX(0._r64, 0._r64, r64)
              ENDDO
            ENDIF

          ENDIF

        ! stop-band filter
        ELSEIF (ft .eq. 4) THEN

          IF (zempty) THEN
            gain = gain * REAL(1._r64 / PRODUCT(-pole(1:np)), r64)
          ELSEIF (pempty) THEN
            gain = gain * REAL(PRODUCT(-zero(1:nz)), r64)
          ELSE
            gain = gain * REAL(PRODUCT(-zero(1:nz)) / PRODUCT(-pole(1:np)), r64)
          ENDIF

          DO i = 1, np
            b            = (fh - fl) / 2._r64 / pole(i)
            x            = SQRT(b**2 - fh*fl)
            pole(i)      = b + x
            pole(np + i) = b - x
          ENDDO

          IF (zempty) THEN

            DO i = 1, 2*np
              k = (-1)**(1 + MOD(i, 2))
              zero(i) = CMPLX(0._r64, SQRT(fh*fl), r64) * (-k)
            ENDDO

          ELSE

            DO i = 1, nz
              b            = (fh - fl) / 2._r64 / zero(i)
              x            = SQRT(b**2 - fh*fl)
              zero(i)      = b + x
              zero(nz + i) = b - x
            ENDDO

            IF (np .gt. nz) THEN

              DO i = 1, 2*(np - nz)
                k = (-1)**(1 + MOD(i, 2))
                zero(i + 2*nz) = CMPLX(0._r64, SQRT(fh*fl), r64) * (-k)
              ENDDO

            ENDIF

          ENDIF

        ENDIF

      ENDIF

    END SUBROUTINE xtrans

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE bilinear(zero, pole, gain)

      ! Purpose:
      !   to compute the bilinear transformation.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   02/09/20                  original version
      !

      COMPLEX(r64), DIMENSION(:), INTENT(INOUT) :: zero, pole            !< poles and zeroes
      REAL(r64),                  INTENT(INOUT) :: gain                  !< gain
      INTEGER(i32)                              :: i
      INTEGER(i32)                              :: np, nz

      !-----------------------------------------------------------------------------------------------------------------------------

      ! number of poles and zeroes
      np = COUNT(pole .ne. EMPTY)
      nz = COUNT(zero .ne. EMPTY)

      gain = gain * REAL(PRODUCT(1._r64 - zero(1:nz)) / PRODUCT(1._r64 - pole(1:np)), r64)

      DO i = 1, np
        pole(i) = (1._r64 + pole(i)) / (1._r64 - pole(i))
      ENDDO

      ! "zero" vector is empty
      IF (nz .eq. 0) THEN

        DO i = 1, np
          zero(i) = CMPLX(-1._r64, 0._r64, r64)
        ENDDO

      ELSE

        DO i = 1, nz
          zero(i) = (1._r64 + zero(i)) / (1._r64 - zero(i))
        ENDDO

        DO i = nz + 1, np
          zero(i) = CMPLX(-1._r64, 0._r64, r64)
        ENDDO

      ENDIF

    END SUBROUTINE bilinear

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    FUNCTION poly(x) RESULT(v)

      ! Purpose:
      !   to convert poles and zeroes to filter coefficients.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   02/09/20                  original version
      !

      COMPLEX(r64), DIMENSION(:),        INTENT(IN) :: x               !< poles or zeroes
      INTEGER(i32)                                  :: i, j, n
      COMPLEX(r64), DIMENSION(SIZE(x)+1)            :: y, c
      REAL(r64),    DIMENSION(SIZE(x)+1)            :: v

      !-----------------------------------------------------------------------------------------------------------------------------

      n = SIZE(x)

      DO i = 1, n + 1
        y(i) = CMPLX(0._r64, 0._r64, r64)
      ENDDO

      y(1) = CMPLX(1._r64, 0._r64, r64)

      DO j = 1, n
        DO i = 1, j
          c(i + 1) = y(i + 1) - x(j) * y(i)
        ENDDO
        DO i = 2, j + 1
          y(i) = c(i)          !< plug result back into y
        ENDDO
      ENDDO

!       DO j = 1, n
!         y(2:(j + 1)) = y(2:(j + 1)) - x(j) * y(1:j)
!       ENDDO

      v = REAL(y, r64)

    END FUNCTION poly

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE freqz(freq, amp, phase, dt, npts)

      ! Purpose:
      !   to compute the impulse response of a filter in terms of modulus ("amp") and unwrapped phase ("phase") at desired frequency
      !   points "freq" defined by parameters "dt" (time-step) and "npts" (number of points). In particular a frequency vector is
      !   defined as follows:
      !
      !   a) both "dt" and "npts" arguments are given => df = 1/(npts*dt),   fmax = 0.5/dt
      !   b) only argument "dt" is given              => df = 1/(1024*dt),   fmax = 0.5/dt
      !   c) only argument "npts" is given            => df = 1/(2*PI*npts), fmax = PI
      !   d) nor "npts" nor "dt" are given            => df = 1/(2*PI*1024), fmax = PI
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   02/09/20                  original version
      !

      REAL(r__),    ALLOCATABLE, DIMENSION(:),                     INTENT(OUT) :: freq
      REAL(r__),    ALLOCATABLE, DIMENSION(:),                     INTENT(OUT) :: amp
      REAL(r__),    ALLOCATABLE, DIMENSION(:),                     INTENT(OUT) :: phase
      REAL(r__),                                         OPTIONAL, INTENT(IN)  :: dt
      INTEGER(i32),                                      OPTIONAL, INTENT(IN)  :: npts
      INTEGER(i32)                                                             :: i, l                  !< counters
      INTEGER(i32)                                                             :: nfft                  !< number of points fft
      REAL(r__)                                                                :: sf                    !< sampling frequency
      REAL(r64),    ALLOCATABLE, DIMENSION(:)                                  :: num, den
      COMPLEX(r64)                                                             :: h
      COMPLEX(r64), ALLOCATABLE, DIMENSION(:)                                  :: u, v

      !-----------------------------------------------------------------------------------------------------------------------------

      ! set number of points in fft according to user-defined or default value
      IF (PRESENT(npts)) THEN
        nfft = npts
      ELSE
        nfft = 1024
      ENDIF

      ! prepare fft plan
      CALL make_fftw_plan([nfft])

      ! default sampling frequency
      sf = 2._r__ * PI

      ! update according to time-step
      IF (PRESENT(dt)) sf = 1._r__ / dt

      ALLOCATE(freq(nfft/2 + 1))

      ! scale sampling frequency by fft points
      sf = sf / REAL(nfft, r__)

      ! build frequency vector
      DO i = 1, nfft/2 + 1
        freq(i) = (i - 1) * sf
      ENDDO

      ! now prepare arrays of filter coefficient for fft
      ALLOCATE(num(nfft), den(nfft))
      ALLOCATE(u(nfft/2 + 1), v(nfft/2 + 1))

      num(:) = 0._r64
      den(:) = 0._r64

      IF (ALLOCATED(iira)) THEN

        DO i = 1, SIZE(iira)
          den(i) = CMPLX(iira(i), 0._r64, r64)
        ENDDO

      ! if "iira" coefficients are not available, we are dealing with a fir filter
      ELSE

        den(1) = CMPLX(1._r64, 0._r64, r64)

      ENDIF

      DO i = 1, SIZE(iirb)
        num(i) = CMPLX(iirb(i), 0._r64, r64)
      ENDDO

      ! forward transforms
      CALL fft(num, u)
      CALL fft(den, v)

      ALLOCATE(amp(nfft/2 + 1), phase(nfft/2 + 1))

      DO i = 1, nfft/2 + 1

        h = u(i) / v(i)                                 !< filter impulse response

        ! modulus
        amp(i) = ABS(h)

        ! phase
        phase(i) = ATAN2(AIMAG(h), REAL(h))

        ! note that, for zero-phase filtering, magnitude is squared and phase is zero
        ! IF (ALLOCATED(zi)) THEN
        !   amp(i)   = amp(i)**2
        !   phase(i) = 0._r64
        ! ENDIF

      ENDDO

      ! unwrap phase...
      phase = unwrap(phase)

      ! ... and convert to degrees
      phase = rad2deg(phase)

      CALL destroy_fftw_plan([nfft])

      DEALLOCATE(num, den, u, v)

    END SUBROUTINE freqz

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    PURE FUNCTION unwrap(x)

      ! Purpose:
      !   to unwrap the input phase "x".
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   02/09/20                  original version
      !

      REAL(r__),   DIMENSION(:), INTENT(IN) :: x
      INTEGER(i32)                          :: i, n
      REAL(r__)                             :: twopi
      REAL(r__),   DIMENSION(SIZE(x))       :: unwrap
      REAL(r__),   DIMENSION(SIZE(x))       :: d, peaks

      !-----------------------------------------------------------------------------------------------------------------------------

      twopi = 2._r__ * PI

      n = SIZE(x)

      ! first order difference
      d(1) = 0._r__

      DO i = 2, n
        d(i) = x(i - 1) - x(i)
      ENDDO

      peaks(:) = 0._r__

      ! compute corrections (deltas) to unwrap phase
      WHERE(d > PI)
        peaks = nint(ABS(d) / twopi) * twopi
      ELSEWHERE(d < -PI)
        peaks = -nint(ABS(d) / twopi) * twopi
      ENDWHERE

      ! integral of corrections (deltas => steps)
      unwrap = cumsum(peaks)

      ! sum corrections to input phase
      DO i = 1, n
        unwrap(i) = x(i) + unwrap(i)
      ENDDO

    END FUNCTION unwrap

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    PURE FUNCTION cumsum(x) RESULT(v)

      ! Purpose:
      !   to compute the cumulative sum of a vector.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   02/09/20                  original version
      !

      REAL(r__),   DIMENSION(:),      INTENT(IN) :: x
      REAL(r__),   DIMENSION(SIZE(x))            :: v
      INTEGER(i32)                               :: i, n

      !-----------------------------------------------------------------------------------------------------------------------------

      n = SIZE(x)

      v(1) = x(1)

      DO i = 2, n
        v(i) = v(i - 1) + x(i)
      ENDDO

    END FUNCTION cumsum

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    ELEMENTAL FUNCTION rad2deg(x) RESULT(y)

      ! Purpose:
      !   to convert radiants into degrees.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   02/09/20                  original version
      !

      REAL(r__), INTENT(IN) :: x
      REAL(r__)             :: y

      !-----------------------------------------------------------------------------------------------------------------------------

      y = x * 180._r__ / PI

    END FUNCTION rad2deg

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    FUNCTION initial_conditions(b, a, ok) RESULT(v)

      ! Purpose:
      !   to compute initial conditions based on filter coefficients to avoid transient at startup for zero-phase filtering
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   02/09/20                  original version
      !

      REAL(r64),                 DIMENSION(:),                   INTENT(IN)  :: b                    !< filter coefficients (numerator)
      REAL(r64),                 DIMENSION(:),                   INTENT(IN)  :: a                    !< filter coefficients (denominator)
      INTEGER(i32),                                              INTENT(OUT) :: ok
      INTEGER(i32)                                                           :: i, j, n
      REAL(r64)                                                              :: rsq
      REAL(r64),                 DIMENSION(SIZE(b)-1)                        :: v
      REAL(r64),                 DIMENSION(SIZE(b)-1, SIZE(b)-1)             :: m, id

      !-----------------------------------------------------------------------------------------------------------------------------

      ok = 0

      m(:,:) = 0._r64

      n = SIZE(b) - 1

      DO i = 1, n
        m(i, 1) = -a(i + 1)
      ENDDO

      DO j = 2, n
        m(j - 1, j) = 1._r64
      ENDDO

      id(:,:) = 0._r64

      DO j = 1, n
        id(j, j) = 1._r64
      ENDDO

      DO j = 1, n
        DO i = 1, n
          m(i, j) = id(i, j) - m(i, j)
        ENDDO
      ENDDO

      ! input/output vector for svd
      DO i = 1, n
        v(i) = b(i + 1) - b(1) * a(i + 1)
      ENDDO

      ! solve linear least-square problem
      CALL llsq_solver(m, v, ok)

      IF (ok .ne. 0) ok = 16

    END FUNCTION initial_conditions

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE impz(x)

      ! Purpose:
      !   to evaluate the impulse response of a IIR filter by filtering an impulse whose size is given by "x".
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   02/09/20                  original version
      !

      REAL(r__),   DIMENSION(:),      INTENT(OUT) :: x                    !< impulse response
      INTEGER(i32)                                :: l, n                 !< counters
      REAL(r__),   DIMENSION(SIZE(x))             :: delta
      INTEGER(i32)                                :: ok

      !-----------------------------------------------------------------------------------------------------------------------------

      n = SIZE(x)

      l = n / 2 + 1

      ! define impulse
      delta    = 0._r__
      delta(l) = 1._r__

      ! filter delta to compute impulse response
      x = iir(delta, ok)

    END SUBROUTINE impz

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE get_iir_coefficients(a, b, z)

      ! Purpose:
      !   to return the filter coefficients created when calling "make_iir_plan".
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   02/09/20                  original version
      !

      REAL(r64), ALLOCATABLE, DIMENSION(:),           INTENT(OUT) :: a, b
      REAL(r64), ALLOCATABLE, DIMENSION(:), OPTIONAL, INTENT(OUT) :: z

      !-----------------------------------------------------------------------------------------------------------------------------

      a = iira
      b = iirb

      IF (PRESENT(z) .and. ALLOCATED(zi)) z = zi

    END SUBROUTINE get_iir_coefficients

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE set_iir_coefficients(a, b, z)

      ! Purpose:
      !   to assign filter coefficients (created by calling "make_iir_plan").
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   02/09/20                  original version
      !

      REAL(r64), DIMENSION(:),           INTENT(IN) :: a, b
      REAL(r64), DIMENSION(:), OPTIONAL, INTENT(IN) :: z

      !-----------------------------------------------------------------------------------------------------------------------------

      !$omp parallel default(shared)

      iira = a
      iirb = b

      IF (PRESENT(z)) zi = z

! print*, omp_get_thread_num(), size(iira), size(iirb), ALLOCATED(zi)

      !$omp end parallel

    END SUBROUTINE set_iir_coefficients

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *


    FUNCTION filter_error(ierr) RESULT(msg)

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
          msg = 'missing argument in make_iir_plan: time-step ("dt") not specified'

        CASE(2)
          msg = 'missing argument in make_iir_plan: filter type ("type") not specified'

        CASE(3)
          msg = 'missing argument in make_iir_plan: corner frequency vector ("fc") not specified'

        CASE(4)
          msg = 'missing argument in make_iir_plan: filter order ("order") not specified'

        CASE(5)
          msg = 'wrong argument in make_iir_plan: unknown filter design ("design")'

        CASE(6)
          msg = 'missing argument in make_iir_plan: ripple PARAMETER ("ripple") not specified for Chebychev/Elliptic filters'

        CASE(7)
          msg = 'wrong argument size in make_iir_plan: size("ripple") = 1 for Chebychev filters'

        CASE(8)
          msg = 'wrong argument size in make_iir_plan: size("ripple") = 2 for Elliptic filters'

        CASE(9)
          msg = 'wrong argument size in make_iir_plan: filter order ("order") must be greater than 0'

        CASE(10)
          msg = 'wrong argument in make_iir_plan: unknown filter type ("type")'

        CASE(11)
          msg = 'wrong argument size in make_iir_plan: size("fc") = 2 for bandpass/bandstop filters'

        CASE(12)
          msg = 'wrong argument size in make_iir_plan: size("fc") = 1 for lowpass/highpass filters'

        CASE(13)
          msg = 'wrong argument in make_iir_plan: corner frequency ("fc") must be in the interval (0, nyquist frequency)'

        CASE(14)
          msg = 'wrong argument in make_iir_plan: ripple vector ("ripple") must be either in linear units or dB'

        CASE(15)
          msg = 'error in xtrans: found more zeroes than poles'

        CASE(16)
          msg = 'error in initial_conditions: llq_solver failed'

        CASE(17)
          msg = 'error in make_iir_plan: could not open file "filter_impulse_response.txt"'

        CASE(18)
          msg = 'error in iir: filter coefficients not available. Call make_iir_plan first'

        CASE(19)
          msg = 'error in iir: input sequence too short for current filter order'

      END SELECT

    END FUNCTION filter_error

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

END MODULE m_filter
