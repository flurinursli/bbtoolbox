MODULE m_isochron

  USE                :: omp_lib
  USE, NON_INTRINSIC :: m_precisions
  USE, NON_INTRINSIC :: m_interpolation_r32
  USE, NON_INTRINSIC :: m_compgeo
  USE, NON_INTRINSIC :: m_stat
  USE, NON_INTRINSIC :: m_source
  USE, NON_INTRINSIC :: m_rik
  USE, NON_INTRINSIC :: m_roughness
  USE, NON_INTRINSIC :: m_toolbox
  USE, NON_INTRINSIC :: m_strings
  USE, NON_INTRINSIC :: m_logfile

  IMPLICIT none

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  PRIVATE

  PUBLIC :: solve_isochron_integral, node2disk

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  REAL(r32), PARAMETER :: PI = 3.14159265358979323846_r64
  REAL(r32), PARAMETER :: DEG_TO_RAD = PI / 180._r32
  REAL(r32), PARAMETER :: BIG = HUGE(0._r32)

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  PROCEDURE(rik_at_nodes), POINTER :: nodefun

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  CONTAINS

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE solve_isochron_integral(ok, pl, vel, iter)

      INTEGER(i32),                           INTENT(OUT) :: ok
      INTEGER(i32),                           INTENT(IN)  :: pl, vel, iter
      INTEGER(i32)                                        :: ref, totnutr, i, j, icr, seed
      INTEGER(i32),              DIMENSION(3)             :: iuc, ivc
      REAL(r32),                 DIMENSION(3)             :: slip, dslip, sslip, u, v, rake
      REAL(r32),    ALLOCATABLE, DIMENSION(:)             :: depth, rho, vs, vp, vsgrad, vpgrad

      !-----------------------------------------------------------------------------------------------------------------------------

      depth  = input%velocity(vel)%depth
      rho    = input%velocity(vel)%rho
      vs     = input%velocity(vel)%vs
      vp     = input%velocity(vel)%vp
      vsgrad = input%velocity(vel)%vsgrad
      vpgrad = input%velocity(vel)%vpgrad

      ! set "seed" such that random numbers depend on fault plane number and iteration
      seed = input%source%seed + (iter - 1) * SIZE(plane) + pl

      DO ref = 1, SIZE(nvtr)         !< loop over mesh refinements

        ALLOCATE(nodes(nugr(ref), nvgr(ref)))

        !CALL nodefun(ref, pl, vel, seed)               !< define slip, rupture time and rise time on mesh nodes

        totnutr = nvtr(ref) * (2 * nutr(ref) - 1)  !< total triangles in a row

        !$omp parallel do default(shared) private(i, j, iuc, ivc, slip)
        DO j = 1, nvtr(ref)
          DO i = 1, totnutr

            CALL cornr(j, i, iuc, ivc)        !< corner indices for current triangle

            slip(1) = SUM(nodes(iuc(1), ivc(1))%slip)
            slip(2) = SUM(nodes(iuc(2), ivc(2))%slip)
            slip(3) = SUM(nodes(iuc(3), ivc(3))%slip)

            IF (ALL(slip .eq. 0._r32)) CYCLE          !< jump to next triangle if current has zero slip

            CALL cornr2uv(iuc, ivc, ref, u, v)        !< on-fault corners coordinates

            ! CALL interpolate(plane(pl)%u, plane(pl)%v, plane(pl)%sslip, u, v, sslip)
            ! CALL interpolate(plane(pl)%u, plane(pl)%v, plane(pl)%dslip, u, v, dslip)
            !
            ! rake(1) = ATAN2(-dslip(1), sslip(1))      !< rake in radians
            ! rake(2) = ATAN2(-dslip(2), sslip(2))
            ! rake(3) = ATAN2(-dslip(3), sslip(3))

          ENDDO
        ENDDO
        !$omp end parallel do

        CALL dealloc_nodes()

      ENDDO        !< end loop over mesh refinements


    END SUBROUTINE solve_isochron_integral

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE node2disk(ok, pl, vel, iter)

      ! Purpose:
      !   to write to disk the rupture parameters defined in a mesh array.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   08/03/21                  original version
      !

      INTEGER(i32),                          INTENT(OUT) :: ok
      INTEGER(i32),                          INTENT(IN)  :: pl, vel, iter
      CHARACTER(:), ALLOCATABLE                          :: fo
      INTEGER(i32)                                       :: i, j, ref, lu, icr, totnutr, seed
      INTEGER(i32),             DIMENSION(3)             :: iuc, ivc
      REAL(r32)                                          :: rho, beta, mu, m0, moment, strike, dip
      REAL(r32),                DIMENSION(3)             :: u, v, w, x, y, z, slip, rise, rupture, rake, nrl

      !-----------------------------------------------------------------------------------------------------------------------------

      CALL setup_interpolation('linear', 'zero', ok)

      fo = 'node_' + num2char(pl) + '_' + num2char(vel) + '_' + num2char(iter) + '.bin'

      OPEN(newunit = lu, file = fo, status = 'replace', form = 'unformatted', access = 'stream', action = 'write', IOSTAT = ok)

      IF (ok .ne. 0) THEN
        CALL report_error('Error while opening file' + TRIM(fo))
        RETURN
      ENDIF

      WRITE(lu, POS=1) SUM([nvtr * (2 * nutr - 1)])         !< total number of triangles (all refinements)

      WRITE(lu) SIZE(nvtr)

      moment = 0._r32

      nodefun => rik_at_nodes

      IF (input%source%is_point) nodefun => ptrsrc_at_nodes

      ! set "seed" such that random numbers depend on fault plane number and iteration
      seed = input%source%seed + (iter - 1) * SIZE(plane) + pl

      DO ref = 1, SIZE(nvtr)         !< loop over mesh refinements

        ALLOCATE(nodes(nugr(ref), nvgr(ref)))

        CALL nodefun(ref, pl, vel, seed)                !< define slip, rupture time and rise time on mesh nodes

        totnutr = 2 * nutr(ref) - 1                     !< total triangles in a row

        m0 = 0._r32

        WRITE(lu) nvtr(ref), nutr(ref)

        CALL fault_roughness(ok, ref, pl, iter)

        !$omp parallel do ordered default(shared) private(i, j, iuc, ivc, slip, rise, rake, nrl, rupture, icr, u, v, w, x, y, z)  &
        !$omp& private(beta, rho, mu) reduction(+:m0) collapse(2)
        DO j = 1, nvtr(ref)
          DO i = 1, totnutr

            CALL cornr(j, i, iuc, ivc)        !< corner indices for current triangle

            DO icr = 1, 3
              slip(icr) = SUM(nodes(iuc(icr), ivc(icr))%slip)
              rise(icr) = mean(nodes(iuc(icr), ivc(icr))%rise)
              rupture(icr) = MINVAL(nodes(iuc(icr), ivc(icr))%rupture)
            ENDDO

            CALL cornr2uv(iuc, ivc, ref, u, v)        !< on-fault coordinates

            DO icr = 1, 3
              w(icr) = roughness(iuc(icr), ivc(icr))
            ENDDO

            CALL uvw2xyz(pl, u, v, w, x, y, z)            !< get cartesian coordinates

            DO icr = 1, 3
              z(icr) = MAX(MIN_DEPTH, z(icr))            !< make sure we never breach the free-surface (clip roughness)
            ENDDO

            !$omp ordered
            WRITE(lu) x, y, z, slip, rise, rupture         !< values at each corner
            !$omp end ordered

            CALL normal2tri(x, y, z, nrl)

            ! always have normal pointing upward (i.e. negative)
            IF (MOD(i + (j-1)*totnutr, 2) == 0) nrl(:) = -nrl(:)

            strike = plane(pl)%strike * DEG_TO_RAD
            dip    = plane(pl)%dip    * DEG_TO_RAD

            CALL interpolate(plane(pl)%u, plane(pl)%v, plane(pl)%rake, u, v, rake)

            CALL perturbed_mechanism(strike, dip, rake, nrl)            !< all values in radians

            !$omp ordered
            WRITE(lu) strike, dip, rake, w     !< strike and dip are scalar
            !$omp end ordered

            beta = 0._r32
            rho  = 0._r32

            DO icr = 1, 3
              beta = beta + vinterp(input%velocity(vel)%depth, input%velocity(vel)%vs, input%velocity(vel)%vsgrad, z(icr))
              rho  = rho  + vinterp(input%velocity(vel)%depth, input%velocity(vel)%rho, input%velocity(vel)%rhograd, z(icr))
            ENDDO

            rho = rho / 3._r32
            beta = beta / 3._r32

            mu = rho * beta**2
            m0 = m0 + mu * mean(slip)

          ENDDO
        ENDDO
        !$omp end parallel do

        CALL dealloc_nodes()

        moment = moment + m0 * dutr(ref) * dvtr(ref) * 0.5_r32

      ENDDO        !< end loop over mesh refinements

      WRITE(lu) plane(pl)%targetm0 / moment

      CLOSE(lu, iostat = ok)

      IF (ok .ne. 0) THEN
        CALL report_error('Error while closing file ' + fo)
        RETURN
      ENDIF

    END SUBROUTINE node2disk

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE perturbed_mechanism(strike, dip, rake, nrl)

      ! Purpose:
      !   to return strike dip and rake of a triangle (described by its normal "nrl") whose initial orientation has been changed and
      !   now described by its normal "nrl". In input, "strike" and "dip" are two scalars and "rake" a vector representing the
      !   initial values at each corner, while in output they refer to the new orientation.
      !   The algorithm first computes the normal and slip vectors for the unperturbed case (as given by strike, dip and rake). These
      !   are used to define the pressure axis P and then (by plugging in the new normal) the new rake angle.
      !   WARNING: input values are expected to follow the N, E, D coordinates system, with dip <= 90. Internally the E, N, U system
      !            is assumed. "nrl" must be normalized to 1.
      !   Output values are in the range [0, pi] or [-pi, 0] (slip and rake) or [0, pi] dip.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   08/03/21                  original version
      !

      REAL(r32),                 INTENT(INOUT) :: strike, dip
      REAL(r32),   DIMENSION(3), INTENT(INOUT) :: rake
      REAL(r32),   DIMENSION(3), INTENT(IN)    :: nrl
      INTEGER(i32)                             :: i
      REAL(r32)                                :: p, cs, ss, cd, sd, cr, sr, arg, hn2
      REAL(r32),   DIMENSION(3)                :: unrl, ul, l, h

      !-----------------------------------------------------------------------------------------------------------------------------

      ! strike vector (perturbed case)
      h(1) = -nrl(1)   !< east
      h(2) = nrl(2)    !< north
      h(3) = 0._r32

      hn2 = NORM2(h)

      cs = COS(strike)
      ss = SIN(strike)
      cd = COS(dip)
      sd = SIN(dip)

      ! normal vector (unperturbed case)
      unrl(1) = cs * sd     !< east
      unrl(2) = -ss * sd    !< north
      unrl(3) = cd          !< up

      DO i = 1, 3

        cr = COS(rake(i))
        sr = SIN(rake(i))

        ! slip vector (unperturbed case)
        ul(1) = ss * cr - cs * cd * sr    !< east
        ul(2) = cs * cr + ss * cd * sr    !< north
        ul(3) = sd * sr                   !< up

        ! slip vector for triangle with normal "nrl", where P axis is based on unperturbed geometry
        ! sqrt(2)*P = unrl - ul
        ! l = nrl - sqrt(2)*P = nrl - (unrl - ul)
        ! note that we reverted (1,2) indices for input normal as it assumes N-E-D reference system
        l(1) = nrl(2) - (unrl(1) - ul(1))
        l(2) = nrl(1) - (unrl(2) - ul(2))
        l(3) = -nrl(3) - (unrl(3) - ul(3))

        p = NORM2(l) * hn2

        arg = (h(1)*l(1) + h(2)*l(2)) / p

        arg = SIGN(MIN(1._r32, ABS(arg)), arg)          !< prevent roundoff errors, enforcing |arg| <= 1

        ! sign of rake angle depends on vertical component of slip vector
        rake(i) = SIGN(ACOS(arg), l(3))       !< values in range [0 pi] or [-pi 0]

      ENDDO

      ! remember that ATAN2(y,x) is in the interval [0,pi] when y>=0, in the interval (-pi, 0) otherwise

      ! perturbed strike in the range [0 pi] or [-pi 0]
      strike = ATAN2(-dble(nrl(1)), dble(nrl(2)))

      ! IF (strike .lt. 0._r32) strike = 2._r32 * PI + strike

      ! perturbed dip in the range [0 pi]
      dip = ATAN2(HYPOT(dble(nrl(1)), dble(nrl(2))), -dble(nrl(3)))    !< revert sign vertical because input reference system is N-E-D

    END SUBROUTINE perturbed_mechanism

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *



    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *



END MODULE m_isochron
