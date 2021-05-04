MODULE m_isochron

  USE                :: omp_lib
  USE, NON_INTRINSIC :: m_precisions
  USE, NON_INTRINSIC :: m_interpolation_r32
  USE, NON_INTRINSIC :: m_compgeo
  USE, NON_INTRINSIC :: m_source
  USE, NON_INTRINSIC :: m_rik
  USE, NON_INTRINSIC :: m_toolbox
  USE, NON_INTRINSIC :: m_strings
  USE, NON_INTRINSIC :: m_logfile

  IMPLICIT none

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  PRIVATE

  PUBLIC :: solve_isochron_integral, node2disk

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

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

            CALL interpolate(plane(pl)%u, plane(pl)%v, plane(pl)%sslip, u, v, sslip)
            CALL interpolate(plane(pl)%u, plane(pl)%v, plane(pl)%dslip, u, v, dslip)

            rake(1) = ATAN2(-dslip(1), sslip(1))      !< rake in radians
            rake(2) = ATAN2(-dslip(2), sslip(2))
            rake(3) = ATAN2(-dslip(3), sslip(3))

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
      REAL(r32)                                          :: slip, rise, rupture, rho, beta, mu, m0, moment
      REAL(r32),                DIMENSION(3)             :: u, v, x, y, z

      !-----------------------------------------------------------------------------------------------------------------------------

      fo = 'node_' + num2char(pl) + '_' + num2char(vel) + '_' + num2char(iter) + '.bin'

      OPEN(newunit = lu, file = fo, status = 'replace', form = 'unformatted', access = 'stream', action = 'write', IOSTAT = ok)

      IF (ok .ne. 0) THEN
        CALL report_error('Error while opening file' + TRIM(fo))
        RETURN
      ENDIF

      WRITE(lu, POS=1) SUM([nvtr * (2 * nutr - 1)])         !< total number of triangles (all refinements)

      moment = 0._r32

      nodefun => rik_at_nodes

      IF (input%source%is_point) nodefun => ptrsrc_at_nodes

      ! set "seed" such that random numbers depend on fault plane number and iteration
      seed = input%source%seed + (iter - 1) * SIZE(plane) + pl

      DO ref = 1, SIZE(nvtr)         !< loop over mesh refinements

        ALLOCATE(nodes(nugr(ref), nvgr(ref)))

        CALL nodefun(ref, pl, vel, seed)               !< define slip, rupture time and rise time on mesh nodes

        totnutr = 2 * nutr(ref) - 1                     !< total triangles in a row

        m0 = 0._r32

        !$omp parallel do ordered default(shared) private(i, j, iuc, ivc, slip, rise, rupture, icr, u, v, x, y, z, beta, rho, mu)  &
        !$omp& reduction(+:m0) collapse(2)
        DO j = 1, nvtr(ref)
          DO i = 1, totnutr

            CALL cornr(j, i, iuc, ivc)        !< corner indices for current triangle

            slip    = 0._r32
            rise    = 0._r32
            rupture = 0._r32

            DO icr = 1, 3
              slip    = slip + SUM(nodes(iuc(icr), ivc(icr))%slip)
              rise    = rise + SUM(nodes(iuc(icr), ivc(icr))%rise)
              rupture = rupture + MINVAL(nodes(iuc(icr), ivc(icr))%rupture)
            ENDDO

            slip    = slip / 3._r32
            rise    = rise / 3._r32
            rupture = rupture / 3._r32

            CALL cornr2uv(iuc, ivc, ref, u, v)        !< on-fault coordinates

            CALL uv2xyz(pl, u, v, x, y, z)            !< cartesian coordinates

            !!$omp critical
            !$omp ordered
            WRITE(lu) x, y, z, slip, rise, rupture
            !$omp end ordered
            !!$omp end critical

            beta = 0._r32
            rho  = 0._r32

            DO icr = 1, 3
              beta = beta + vinterp(input%velocity(vel)%depth, input%velocity(vel)%vs, input%velocity(vel)%vsgrad, z(icr))
              rho  = rho  + vinterp(input%velocity(vel)%depth, input%velocity(vel)%rho, input%velocity(vel)%rhograd, z(icr))
            ENDDO

            rho = rho / 3._r32
            beta = beta / 3._r32

            mu = rho * beta**2
            m0 = m0 + mu * slip

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



    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *



    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *



END MODULE m_isochron
