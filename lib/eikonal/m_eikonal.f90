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

    SUBROUTINE fun_c(gridsize, gridstep, src, velocity, arrivals) BIND(c, name="fmm_c")
      USE, INTRINSIC     :: iso_c_binding
      USE, NON_INTRINSIC :: m_precisions
#ifdef DOUBLE_PREC
      INTEGER, PARAMETER :: c_r__ = c_r64
#else
      INTEGER, PARAMETER :: c_r__ = c_r32
#endif
      INTEGER(c_int), DIMENSION(2), INTENT(IN) :: gridsize
      REAL(c_r__),    DIMENSION(2), INTENT(IN) :: gridstep
      INTEGER(c_int), DIMENSION(2), INTENT(IN) :: src
      TYPE(c_ptr),                  VALUE      :: velocity, arrivals
    END SUBROUTINE fun_c

  END INTERFACE

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  CONTAINS

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE fast_marching(ok, source, gridstep, velocity, traveltime)

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
      INTEGER(i32), DIMENSION(2),                       INTENT(IN)  :: source
      REAL(r__),    DIMENSION(2),                       INTENT(IN)  :: gridstep
      REAL(r__),    DIMENSION(:,:), CONTIGUOUS, TARGET, INTENT(IN)  :: velocity       !< get copied if not contiguous
      REAL(r__),    DIMENSION(:,:), CONTIGUOUS, TARGET, INTENT(OUT) :: traveltime
      INTEGER(i32), DIMENSION(2)                                    :: gridsize

      !-----------------------------------------------------------------------------------------------------------------------------

      ok = 0

      gridsize = [SIZE(velocity, 1), SIZE(velocity, 2)]

      IF ( (source(1) .lt. 1) .or. (source(1) .gt. gridsize(1)) .or. (source(2) .lt. 1) .or. (source(2) .gt. gridsize(2)) ) THEN
        ok = 1
        RETURN
      ENDIF

      CALL fun_c(gridsize, gridstep, source, C_LOC(velocity), C_LOC(traveltime))

    END SUBROUTINE fast_marching

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *


    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

END MODULE m_eikonal
