MODULE m_rik

  USE                :: omp_lib
  USE, NON_INTRINSIC :: m_precisions
  USE, NON_INTRINSIC :: m_interpolation_r32
  USE, NON_INTRINSIC :: m_logfile
  USE, NON_INTRINSIC :: m_random
  USE, NON_INTRINSIC :: m_source
  USE, NON_INTRINSIC :: m_stat
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
  REAL(r32),    PARAMETER :: PI = 3.14159265358979323846_r64
  REAL(r32),    PARAMETER :: DEG_TO_RAD = PI / 180._r32
  REAL(r32),    PARAMETER :: TWOPI = 2._r64 * PI


  INTEGER(i32)                            :: subtot
  REAL(i32),    ALLOCATABLE, DIMENSION(:) :: su, sv, sradius

  TYPE :: grd
    REAL(r32), ALLOCATABLE, DIMENSION(:) :: slip, rupture, rise
  END TYPE grd

  TYPE(grd), ALLOCATABLE, DIMENSION(:,:) :: node

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  CONTAINS

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE rik(ok, pl, vel, iter)

      INTEGER(i32),                                     INTENT(OUT) :: ok
      INTEGER(i32),                                     INTENT(IN)  :: pl, vel, iter
      CHARACTER(:), ALLOCATABLE                                     :: fo
      INTEGER(i32)                                                  :: nu, nv, submin, submax, i, j, slevel, n, lu, seed
      INTEGER(i32),              DIMENSION(2)                       :: pos
      INTEGER(i32), ALLOCATABLE, DIMENSION(:)                       :: nsubs
      LOGICAL                                                       :: is_estimated
      REAL(r32)                                                     :: dh, umin, vmin, umax, vmax, du, dv, cross, lc
      REAL(r32),                 DIMENSION(1)                       :: x
      REAL(r32),                 DIMENSION(2)                       :: mu, std
      REAL(r32),                 DIMENSION(PDFNU)                   :: u
      REAL(r32),                 DIMENSION(PDFNV)                   :: v
      REAL(r32),                 DIMENSION(PDFNU,PDFNV)             :: cpdf, slip, tslip

      !-----------------------------------------------------------------------------------------------------------------------------

      ok = 0

      ! set "seed" such that random numbers depend on fault plane number and iteration
      seed = input%source%seed + (iter - 1) * SIZE(plane) + pl

      IF (ALLOCATED(node)) DEALLOCATE(node)

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

      ! limit minimum subsource size to average triangular mesh step
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

      IF (input%advanced%verbose .eq. 2) THEN
        CALL update_log(num2char('<RIK sources>', justify='c', width=30) +   &
                        num2char('Min radius', width=15, justify='r') + '|' + &
                        num2char('Max radius', width=15, justify='r') + '|' + &
                        num2char('Total', width=15, justify='r') + '|')

        CALL update_log(num2char('', width=30)  +  &
                        num2char((vmax-vmin)/submax/2., width=15, justify='r', notation='f', precision=1) + '|' + &
                        num2char((vmax-vmin)/submin/2., width=15, justify='r', notation='f', precision=1) + '|' + &
                        num2char(subtot, width=15, justify='r') + '|', blankline=.false.)
      ENDIF

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

      CALL setup_rng(ok, 'uniform', 0._r32, 1._r32, seed)         !< initialise random number generator

      DO     !< cycle until RIK slip model correlates well with input slip model

        n = 0

        DO slevel = submin, submax
          DO i = 1, nsubs(slevel)         !< cycle over all subsources for current level

            n = n + 1

            sradius(n) = (vmax - vmin) / (2 * slevel)          !< subsource radius

            DO

              CALL rng(x, 1)

              pos = MINLOC(ABS(cpdf - x(1)))

              su(n) = (pos(1) - 1) * du + umin        !< subsource position
              sv(n) = (pos(2) - 1) * dv + vmin

              ! accept subsource if located inside strong motion area. This check is necessary because "x" can be very close to 0
              IF ( (su(n) .gt. umin) .and. (su(n) .lt. umax) .and. (sv(n) .gt. vmin) .and. (sv(n) .lt. vmax) ) EXIT

            ENDDO

          ENDDO
        ENDDO

        !$omp parallel do default(shared) private(i, j)
        DO j = 1, PDFNV
          DO i = 1, PDFNU
            slip(i, j) = total_slip(u(i), v(j))
          ENDDO
        ENDDO
        !$omp end parallel do

        mu(1) = mean(slip)
        mu(2) = mean(tslip)

        std(1) = SQRT(variance(slip))
        std(2) = SQRT(variance(tslip))

        cross = 0._r32

        DO j = 1, PDFNV
          DO i = 1, PDFNU
            cross = cross + (slip(i, j) - mu(1)) * (tslip(i, j) - mu(2))
          ENDDO
        ENDDO

        cross = cross / PRODUCT(std) / (PDFNU * PDFNV)     !< zero-normalized cross-correlation

        IF (input%advanced%verbose .eq. 2) THEN
          CALL update_log(num2char('<slip correlation>', justify='c', width=30) +    &
                          num2char(cross, width=15, notation='f', precision=3, justify='r') + '|', blankline=.false.)
        ENDIF

        IF (cross .ge. input%source%correlation) EXIT

      ENDDO

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! -------------------------------------------- add rupture time if missing ---------------------------------------------------

      IF (plane(pl)%rupture(2, 2) .lt. 0._r32) THEN
        is_estimated = .true.
      ELSE
        is_estimated = .false.
      ENDIF

      IF (is_estimated) CALL missing_rupture(ok, pl, vel)

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! ------------------------------------------- generate rupture model on grid -------------------------------------------------

      lc = input%source%l0 * (umax - umin)          !< rescale "l0" after strong motion area

      ALLOCATE(node(nugr, nvgr))

      CALL rupture_on_grid(pl, vel, seed, lc)

      DEALLOCATE(su, sv, sradius)

      IF (is_estimated) THEN
        DO i = 1, SIZE(plane)
          plane(i)%rupture(:,:) = -999._r32       !< reset to not-available
        ENDDO
      ENDIF

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
      INTEGER(i32)                                     :: i, j, ok
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

    REAL(r32) FUNCTION total_slip(u, v)

      ! Purpose:
      !   to compute slip at fault position "u" "v" as contribution of all subsources.
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

        x = MAX(0._r32, x)

        total_slip = total_slip + SQRT(x)

        !IF (x .gt. 0._r32) total_slip = total_slip + SQRT(x)

      ENDDO

    END FUNCTION total_slip

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE rupture_on_grid(pl, vel, seed, lc)

      INTEGER(i32),                                  INTENT(IN) :: pl, vel, seed
      REAL(r32),                                     INTENT(IN) :: lc
      INTEGER(i32)                                              :: iv, iu, nsubs, n, i, skip, ok
      REAL(r32)                                                 :: x, v, du, u, rise, rupture, phi, r, subvr, sd, z0
      REAL(r32),                 DIMENSION(3)                   :: lvec, hvec
      REAL(r32),                 DIMENSION(subtot)              :: sub2node, pu, pv, subrupt
      REAL(r32),                 DIMENSION(2,subtot)            :: z
      REAL(r32),    ALLOCATABLE, DIMENSION(:)                   :: depth, vs, vsgrad

      !-----------------------------------------------------------------------------------------------------------------------------

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! -------------------------------------------- generate nucleation points  ---------------------------------------------------

      CALL setup_interpolation('linear', 'zero', ok)

      CALL setup_rng(ok, 'uniform', 0._r32, 1._r32, seed)

      CALL rng(z)      !< this is parallelized if openmp flag is set

      !$omp parallel do default(shared) private(phi, r, i)
      DO i = 1, subtot

        phi = z(1, i) * TWOPI
        r   = SQRT(z(2, i)) * sradius(i)          !< not sure why we have SQRT here

        pu(i) = r * COS(phi) + su(i)
        pv(i) = r * SIN(phi) + sv(i)

        pu(i) = MAX(pu(i), umingr)
        pu(i) = MIN(pu(i), umaxgr)
        pv(i) = MAX(pv(i), vmingr)
        pv(i) = MIN(pv(i), vmaxgr)

        ! compute rupture time from hypocenter to random point inside subsource
        CALL interpolate(plane(pl)%u, plane(pl)%v, plane(pl)%rupture, pu(i), pv(i), subrupt(i))

      ENDDO
      !$omp end parallel do

      IF (input%advanced%verbose .eq. 2) THEN
        CALL update_log(num2char('<min/max rupt-to-point>', justify='c', width=30) +   &
                        num2char(num2char(MINVAL(subrupt), notation='f', width=6, precision=2) + ', ' +  &
                        num2char(MAXVAL(subrupt), notation='f', width=6, precision=2), width=15, justify='r') + '|',   &
                        blankline=.false.)
      ENDIF

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! ---------------------------------------- slip, rupture and rise-time at nodes ----------------------------------------------

      sd = SIN(plane(pl)%dip * DEG_TO_RAD)
      z0 = plane(pl)%z

      ASSOCIATE(model => input%velocity(vel))
        depth  = model%depth
        vs     = model%vs
        vsgrad = model%vsgrad
      END ASSOCIATE

      ! cycle over nodes
      !$omp parallel do default(shared) private(iv, v, du, iu, u, x, sub2node, nsubs, n, i, subvr, rise, rupture)
      DO iv = 1, nvgr

        v = vmingr + (iv - 1) * dvtr

        du = (1 - MOD(iv, 2)) * dutr / 2._r32          !< "du" is 0 for odd "iv", dutr/2 for even "iv"

        DO iu = 1, nugr + MOD(iv, 2) - 1               !< "i" tops at "nugr" for odd "j", at "nugr-1" for even "j"

          u = umingr + (iu - 1) * dutr + du

          ! compute difference between radius and node-to-source distance (positive only when node is inside source)
          DO i = 1, subtot
            x = sradius(i)**2 - (u - su(i))**2 - (v - sv(i))**2
            sub2node(i) = MAX(0._r32, x)                                 !< "0" if outside, "> 0" if inside
          ENDDO

          nsubs = subtot - COUNT(sub2node .eq. 0._r32)          !< find number of subsources affecting current node

          ! no subsources act on current node: set slip to 0 and interpolate rupture time...
          IF (nsubs .eq. 0) THEN

            CALL interpolate(plane(pl)%u, plane(pl)%v, plane(pl)%rupture, u, v, rupture)

            ALLOCATE(node(iu,iv)%slip(1), node(iu,iv)%rupture(1), node(iu,iv)%rise(1))

            node(iu, iv)%slip = 0._r32
            node(iu, iv)%rupture = rupture
            node(iu, iv)%rise = 1._r32

          !... otherwise evaluate and store the contribution of each subsource
          ELSE

            ALLOCATE(node(iu,iv)%slip(nsubs), node(iu,iv)%rupture(nsubs), node(iu,iv)%rise(nsubs))

            n = 0

            ! define slip, rupture time and rise-time on current node
            DO i = 1, subtot

              IF (sub2node(i) .eq. 0._r32) CYCLE       !< skip if current subsource is outside

              n = n + 1

              subvr = meanvr(sv(i), sradius(i), z0, sd, depth, vs, vsgrad) * input%source%vrfact

              IF (2*sradius(i) .ge. lc) THEN

                CALL interpolate(plane(pl)%u, plane(pl)%v, plane(pl)%rupture, u, v, rupture)

                rise = input%source%aparam * lc / subvr

              ELSE

                ! add rupture time from random point to current node
                rupture = subrupt(i) + SQRT((u - pu(n))**2 + (v - pv(n))**2) / subvr

                rise = input%source%aparam * 2 * sradius(i) / subvr

              ENDIF

              node(iu, iv)%slip(n) = SQRT(sub2node(i))      !< note that slip here is simply proportional to final slip
              node(iu, iv)%rupture(n) = rupture
              node(iu, iv)%rise(n) = rise

            ENDDO

          ENDIF

        ENDDO
      ENDDO
      !$omp end parallel do

      IF (input%advanced%verbose .eq. 2) THEN

        lvec(:) = HUGE(0._r32)
        hvec(:) = -lvec(:)

        !$omp parallel do default(shared) private(iu, iv) reduction(min:lvec) reduction(max:hvec)
        DO iv = 1, nvgr
          DO iu = 1, nugr + MOD(iv, 2) - 1
            lvec = [MIN(lvec(1), MINVAL(node(iu, iv)%slip)), MIN(lvec(2), MINVAL(node(iu, iv)%rupture)),    &
                    MIN(lvec(3), MINVAL(node(iu, iv)%rise))]
            hvec = [MAX(hvec(1), MAXVAL(node(iu, iv)%slip)), MAX(hvec(2), MAXVAL(node(iu, iv)%rupture)),    &
                    MAX(hvec(3), MAXVAL(node(iu, iv)%rise))]
          ENDDO
        ENDDO
        !$omp end parallel do

        CALL update_log(num2char('<Min/Max at nodes>', justify='c', width=30) + num2char('Slip', width=15, justify='r') + '|' + &
                        num2char('Rupture', width=15, justify='r') + '|' + num2char('Rise', width=15, justify='r') + '|')

        CALL update_log(num2char('', width=30)  +  &
                        num2char(num2char(lvec(1), notation='f', width=6, precision=1) + ', ' +   &
                                num2char(hvec(1), notation='f', width=6, precision=1), width=15, justify='r') + '|' + &
                        num2char(num2char(lvec(2), notation='f', width=6, precision=1) + ', ' +   &
                                num2char(hvec(2), notation='f', width=6, precision=1), width=15, justify='r') + '|' + &
                        num2char(num2char(lvec(3), notation='f', width=6, precision=2) + ', ' +   &
                                num2char(hvec(3), notation='f', width=6, precision=2), width=15, justify='r') + '|',  &
                        blankline=.false.)

      ENDIF

    END SUBROUTINE rupture_on_grid

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    REAL(r32) FUNCTION meanvr(v, radius, z0, sd, depth, vs, vsgrad)

      ! Purpose:
      !   to compute the average rupture time inside a circular subsource of radius "radius" and centered at "v" (downdip) on a
      !   fault plane having minimum depth "z0", sin(dip) "sd" and embedded in a velocity model whose properties are given by "depth"
      !   "vs" and "vsgrad".
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   08/03/21                  original version
      !

      REAL(r32),                 INTENT(IN) :: v, radius, z0, sd
      REAL(r32),   DIMENSION(:), INTENT(IN) :: depth, vs, vsgrad
      INTEGER(i32)                          :: i
      INTEGER(i32),              PARAMETER  :: npts = 11
      REAL(r32)                             :: dz, z, tau

      !-----------------------------------------------------------------------------------------------------------------------------

      tau = 0._r32

      dz = 2._r32 * radius / npts          !< subdivide circular subsource in "npts"-1 segments

      DO i = 1, npts

        z = v - radius + (i - 0.5_r32) * dz         !< on-fault coordinates (segment midpoint)
        z = z0 + z * sd                             !< absolute coordinates

        z = MAX(0._r32, z)          !< avoid small negative values: is this needed??

        tau = tau + dz / vinterp(depth, vs, vsgrad, z)       !< add time needed by rupture front to cover distance "dz"

      ENDDO

      meanvr = 2._r32 * radius / tau           !< return average rupture speed (assuming vr=vs)

    END FUNCTION meanvr

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    REAL(r32) FUNCTION vinterp(ztop, phys, gradient, z)

      ! Purpose:
      !   to return a physical property (e.g. velocity, density, etc) at depth "z" given a model where property "phys", gradient
      !   "gradient" and depth to layer "ztop" are defined for each layer.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   08/03/21                  original version
      !

      REAL(r32),   DIMENSION(:), INTENT(IN) :: ztop, phys, gradient
      REAL(r32),                 INTENT(IN) :: z
      INTEGER(i32)                          :: i

      !-----------------------------------------------------------------------------------------------------------------------------

      DO i = SIZE(ztop), 1, -1
        IF (z .ge. ztop(i)) THEN
          vinterp = phys(i) + gradient(i)*(z - ztop(i))
          EXIT
        ENDIF
      ENDDO

    END FUNCTION vinterp

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *



END MODULE m_rik
