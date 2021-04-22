MODULE m_eikonal

  USE, NON_INTRINSIC :: m_precisions

  IMPLICIT none

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  PRIVATE

  PUBLIC :: fast_marching

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  ! set precisions of input/output arguments at compile-time
#ifdef DOUBLE_PREC
  INTEGER, PARAMETER :: r__ = r64
#else
  INTEGER, PARAMETER :: r__ = r32
#endif

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  INTERFACE

    SUBROUTINE fun_c(gridsize, gridstep, nsrc, src, values, velocity, arrivals) BIND(c, name="fmm_c")
      USE, INTRINSIC     :: iso_c_binding
      USE, NON_INTRINSIC :: m_precisions
#ifdef DOUBLE_PREC
      INTEGER, PARAMETER :: c_r__ = c_r64
#else
      INTEGER, PARAMETER :: c_r__ = c_r32
#endif
      INTEGER(c_int), DIMENSION(2), INTENT(IN) :: gridsize
      REAL(c_r__),    DIMENSION(2), INTENT(IN) :: gridstep
      INTEGER(c_int),               VALUE      :: nsrc
      TYPE(c_ptr),                  VALUE      :: src, values, velocity, arrivals
    END SUBROUTINE fun_c

  END INTERFACE

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  CONTAINS

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE fast_marching(ok, source, values, gridstep, velocity, traveltime)

      ! Purpose:
      !   To solve the eikonal equation for a 2D medium discretized according to "gridstep" and having variable "velocity", where
      !   indices of a single point-source are given in "source". First dimension represent rows, while second is columns.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   08/03/21                  original version
      !

      INTEGER(i32),                                     INTENT(OUT) :: ok
      INTEGER(i32), DIMENSION(:,:), CONTIGUOUS, TARGET, INTENT(IN)  :: source
      REAL(r__),    DIMENSION(:),   CONTIGUOUS, TARGET, INTENT(IN)  :: values
      REAL(r__),    DIMENSION(2),                       INTENT(IN)  :: gridstep
      REAL(r__),    DIMENSION(:,:), CONTIGUOUS, TARGET, INTENT(IN)  :: velocity       !< get copied if not contiguous
      REAL(r__),    DIMENSION(:,:), CONTIGUOUS, TARGET, INTENT(OUT) :: traveltime
      INTEGER(i32)                                                  :: i, nsrc
      INTEGER(i32), DIMENSION(2)                                    :: gridsize

      !-----------------------------------------------------------------------------------------------------------------------------

      ok = 0

      gridsize = [SIZE(velocity, 1), SIZE(velocity, 2)]

      DO i = 1, SIZE(source, 2)
        IF ( (source(1, i) .lt. 1) .or. (source(1, i) .gt. gridsize(1)) .or. (source(2, i) .lt. 1) .or.     &
             (source(2, i) .gt. gridsize(2)) ) THEN
          ok = 1
          RETURN
        ENDIF
      ENDDO

      nsrc = SIZE(source, 2)

      CALL fun_c(gridsize, gridstep, nsrc, C_LOC(source), C_LOC(values), C_LOC(velocity), C_LOC(traveltime))

    END SUBROUTINE fast_marching

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *


    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

END MODULE m_eikonal
