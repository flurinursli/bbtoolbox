MODULE m_isochron

  USE, NON_INTRINSIC :: m_precisions
  USE, NON_INTRINSIC :: m_source
  USE, NON_INTRINSIC :: m_rik
  USE, NON_INTRINSIC :: m_toolbox

  IMPLICIT none

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  PRIVATE

  PUBLIC :: solve_isochron_integral

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --


  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  CONTAINS

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE solve_isochron_integral(ok, pl, vel)

      INTEGER(i32),                           INTENT(OUT) :: ok
      INTEGER(i32),                           INTENT(IN)  :: pl, vel
      INTEGER(i32)                                        :: ref, totnutr, i, j
      INTEGER(i32),              DIMENSION(3)             :: iuc, ivc
      REAL(r32),                 DIMENSION(3)             :: slip, dslip, sslip, u, v, rake
      REAL(r32),    ALLOCATABLE, DIMENSION(:)             :: depth, rho, vs, vp, vsgrad, vpgrad
      TYPE(grd),    ALLOCATABLE, DIMENSION(:,:)           :: nodes

      !-----------------------------------------------------------------------------------------------------------------------------

      depth  = input%velocity(vel)%depth
      rho    = input%velocity(vel)%rho
      vs     = input%velocity(vel)%vs
      vp     = input%velocity(vel)%vp
      vsgrad = input%velocity(vel)%vsgrad
      vpgrad = input%velocity(vel)%vpgrad

      DO ref = 1, SIZE(nvtr)         !< loop over mesh refinements

        ALLOCATE(nodes(nugr(ref), nvgr(ref)))

        !CALL nodefun(ref, nodes)               !< define slip, rupture time and rise time on mesh nodes

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

        DEALLOCATE(nodes)

      ENDDO        !< end loop over mesh refinements


    END SUBROUTINE solve_isochron_integral

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

END MODULE m_isochron
