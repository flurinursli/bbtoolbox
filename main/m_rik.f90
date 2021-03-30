MODULE m_rik

  USE, NON_INTRINSIC :: m_precisions
  USE, NON_INTRINSIC :: m_interpolation
  USE, NON_INTRINSIC :: m_logfile
  USE, NON_INTRINSIC :: m_random
  USE, NON_INTRINSIC :: m_source
  USE, NON_INTRINSIC :: m_strings
  USE, NON_INTRINSIC :: m_toolbox, ONLY: input
! #ifdef MPI
!   USE :: mpi
! #endif

  IMPLICIT none

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  PRIVATE

  PUBLIC :: rik
  !PUBLIC ::

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  INTEGER(i32), PARAMETER :: PDFNU = 200, PDFNV = 100

  INTEGER(i32)                            :: subtot
  REAL(i32),    ALLOCATABLE, DIMENSION(:) :: su, sv, sradius

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  CONTAINS

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE rik(ok, pl)

      INTEGER(i32),                                      INTENT(OUT) :: ok
      INTEGER(i32),                                      INTENT(IN)  :: pl
      INTEGER(i32)                                                   :: nu, nv, submin, submax, subtot, i, slevel, n
      INTEGER(i32), ALLOCATABLE, DIMENSION(:)                        :: nsubs
      REAL(r32)                                                      :: dh, umin, vmin, umax, vmax, du, dv
      REAL(r32),                 DIMENSION(PDFNU, PDFNV)             :: cpdf, slip, tslip

      !-----------------------------------------------------------------------------------------------------------------------------

      ok = 0

      ! average grid-step in input slip model
      dh = (plane(pl)%u(3) - plane(pl)%u(2)) + (plane(pl)%v(3) - plane(pl)%v(2))
      dh = dh / 2

      nu = SIZE(plane(pl)%u)
      nv = SIZE(plane(pl)%v)

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! --------------------------------------------------- slip-based cpdf  -------------------------------------------------------

      ! define strong motion area: in our case we consider the whole input fault plane excluding only edges where slip was set to 0
      umin = plane(pl)%u(2)
      vmin = plane(pl)%v(2)
      umax = plane(pl)%u(nu - 1)
      vmax = plane(pl)%v(nv - 1)

      CALL fillpdf(pl, umin, umax, vmin, vmax, cpdf)                  !< slip-based cpdf based on strong-motion area

print*, umin, umax, vmin, vmax, minval(cpdf), maxval(cpdf)

      ! limit minimum subsource size to average mesh step
      submin = 2
      submax = NINT((vmax - vmin) / (dutr + dvtr) / 4)

      ALLOCATE(nsubs(submax))

      nsubs(:) = 0

      subtot = 0

      DO i = submin, submax
        nsubs(i) = NINT((2*i - 1) * (umax - umin) / (vmax - vmin))            !< subsources at i-th level
        subtot = subtot + nsubs(i)                                            !< total subsources
      ENDDO

      ALLOCATE(su(subtot), sv(subtot), sradius(subtot))

      du = plane(pl)%length / (PDFNU - 1)
      dv = plane(pl)%width / (PDFNV - 1)

      DO i = 1, PDFNV
        v(i) = (i - 1) * dv + plane(pl)%v(1)
      ENDDO

      DO i = 1, PDFNU
        u(i) = (i - 1) * du + plane(pl)%u(1)
      ENDDO

      CALL setup_interpolation('linear', 'zero', ok)

      CALL interpolate(plane(pl)%u, plane(pl)%v, plane(pl)%tslip, u, v, tslip)

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! ------------------------------------------------ generate slip model  ------------------------------------------------------

      CALL setup_rng(ok, 'uniform', 0._r32, 1._r32, input%source%seed)         !< initialise random number generator

      DO     !< cycle until RIK slip model correlates well with input slip model

        n = 0

        DO slevel = submin, submax
          DO i = 1, nsubs(slevel)         !< cycle over all subsources for current level

            n = n + 1

            sradius(n) = (vmax - vmin) / (2 * slevel)          !< subsource radius

            DO

              CALL rng(x)

              pos = MINLOC(ABS(cpdf - x))

              su(n) = (pos(1) - 1) * du + umin      !< subsource position
              sv(n) = (pos(2) - 1) * dv + vmin

              ! accept subsource if located inside strong motion area. This check is necessary because "x" can be very close to 0
              IF ( (su(n) .gt. umin) .and. (su(n) .lt. umax) .and. (sv(n) .gt. vmin) .and. (sv(n) .lt. vmax) ) EXIT

            ENDDO

          ENDDO
        ENDDO

        slip = total_slip(u, v)

        mu(1) = mean(slip)
        mu(2) = mean(tslip)

        std(1) = sigma(slip)
        std(2) = sigma(tslip)

        DO j = 1, PDFNV
          DO i = 1, PDFNU
            cross = cross + (slip(i, j) - mu(1)) * (tslip(i, j) - mu(2))
          ENDDO
        ENDDO

        cross = cross / PRODUCT(std) / (PDFNU * PDFNV)     !< zero-normalized cross-correlation

print*, cross, mu, std

        IF (cross .ge. input%source%correlation) EXIT

      ENDDO


      lc = input%source%l0 * (umax - umin)          !< rescale "l0" after strong motion area



    END SUBROUTINE rik

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE fillpdf(pl, umin, umax, vmin, vmax, cpdf)

      ! Purpose:
      !   to compute the cumulative and probability density function (cdf, pdf) for the "pl"-th fault plane, whose strong motion
      !   area covers the area delimited by "umin", "umax", "vmin", "vmax" (on-fault coordinates). The fault plane is sampled with
      !   "PDFNU * PDFNV" cells ("pdf" and "cpdf" refer to cell center).
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   08/03/21                  original version
      !

      INTEGER(i32),                        INTENT(IN)  :: pl
      REAL(r32),                           INTENT(IN)  :: umin, umax, vmin, vmax
      REAL(r32),    DIMENSION(:,:),        INTENT(OUT) :: cpdf
      INTEGER(i32)                                     :: i, j
      REAL(r32)                                        :: du, dv, cumul
      REAL(r32),    DIMENSION(PDFNU)                   :: u
      REAL(r32),    DIMENSION(PDFNV)                   :: v
      REAL(r32),    DIMENSION(PDFNU,PDFNV)             :: pdf

      !-----------------------------------------------------------------------------------------------------------------------------

      ! create regular grid covering
      du = plane(pl)%length / (PDFNU - 1)
      dv = plane(pl)%width / (PDFNV - 1)

      DO i = 1, PDFNV
        v(i) = (i - 1) * dv + plane(pl)%v(1)
      ENDDO

      DO i = 1, PDFNU
        u(i) = (i - 1) * du + plane(pl)%u(1)
      ENDDO

      CALL setup_interpolation('linear', 'zero', ok)

      CALL interpolate(plane(pl)%u, plane(pl)%v, plane(pl)%tslip, u, v, pdf)

      cumul = 0._r32

      ! compute "cpdf" setting values to zero if outside strong motion area
      DO j = 1, PDFNV
        DO i = 1, PDFNU
          IF ( ((u(i) .lt. umin) .and. (u(i) .gt. umax)) .or. ((v(j) .lt. vmin) .and. (v(j) .gt. vmax)) ) pdf(i, j) = 0._r32
          cumul = cumul + pdf(i, j)
          cpdf(i, j) = cumul
        ENDDO
      ENDDO

      ! normalize to 1
      DO j = 1, PDFNV
        DO i = 1, PDFNU
          cpdf(i, j) = cpdf(i, j) / cumul
        ENDDO
      ENDDO

    END SUBROUTINE fillpdf

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    ELEMENTAL REAL(r32) FUNCTION total_slip(u, v)

      ! Purpose:
      !   to compute slip at fault position "u" "v" due to contribution of all subsources.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   08/03/21                  original version
      !

      REAL(r32),   INTENT(IN) :: u, v
      INTEGER(i32)            :: i
      REAL(r32)               :: x

      !-----------------------------------------------------------------------------------------------------------------------------

      total_slip = 0._r32

      DO i = 1, subtot

        x = sradius(i)**2 - (u - su(i))**2 - (v - sv(i))**2

        IF (x .gt. 0._r32) total_slip = total_slip + SQRT(x)

      ENDDO

    END FUNCTION total_slip

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *



    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *


END MODULE m_rik
