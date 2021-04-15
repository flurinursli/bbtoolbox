MODULE m_compgeo

  USE, NON_INTRINSIC :: m_precisions

  IMPLICIT none

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  PRIVATE

  PUBLIC :: is_triangle_intersection, rotate

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  ! set precisions of input/output arguments at compile-time
#ifdef DOUBLE_PREC
  INTEGER, PARAMETER :: r__   = r64
  ! INTEGER, PARAMETER :: c_r__ = c_r64
#else
  INTEGER, PARAMETER :: r__   = r32
  ! INTEGER, PARAMETER :: c_r__ = c_r32
#endif

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  INTERFACE

    FUNCTION tritri(v0, v1, v2, u0, u1, u2) BIND(c, name="NoDivTriTriIsect")
      USE, INTRINSIC     :: iso_c_binding
      USE, NON_INTRINSIC :: m_precisions
#ifdef DOUBLE_PREC
      INTEGER, PARAMETER :: c_r__ = c_r64
#else
      INTEGER, PARAMETER :: c_r__ = c_r32
#endif
      REAL(c_r__),    DIMENSION(3), INTENT(IN) :: v0, v1, v2, u0, u1, u2
      LOGICAL(c_bool)                          :: tritri
    END FUNCTION tritri

  END INTERFACE

  INTERFACE rotate
    MODULE PROCEDURE rotate_scalar, rotate_vec
  END INTERFACE rotate

  CONTAINS

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    LOGICAL FUNCTION is_triangle_intersection(v0, v1, v2, u0, u1, u2)

      REAL(r__), DIMENSION(3), INTENT(IN) :: v0, v1, v2
      REAL(r__), DIMENSION(3), INTENT(IN) :: u0, u1, u2

      !-----------------------------------------------------------------------------------------------------------------------------

      is_triangle_intersection = tritri(v0, v1, v2, u0, u1, u2)

    END FUNCTION is_triangle_intersection

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE rotate_scalar(x, y, z, ax, ay, az)

      ! Purpose:
      !   to perform the 3-2-1 intrinsic rotation of a point whose initial position is given by "x", "y" and "z", according to angles
      !   "ax", "ay" and "az" respect to the x-, y- and z-axis. Angles must be given in radians. Axis and angles follow the righ-hand
      !   convention (x is north, y is east, z points downwards).
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   08/03/21                  original version
      !

      REAL(r__),                INTENT(INOUT) :: x, y, z
      REAL(r__),                INTENT(IN)    :: ax, ay, az
      REAL(r__), DIMENSION(1)                 :: xv, yv, zv

      !-----------------------------------------------------------------------------------------------------------------------------

      xv = [x]
      yv = [y]
      zv = [z]

      CALL rotate(xv, yv, zv, ax, ay, az)

      x = xv(1)
      y = yv(1)
      z = zv(1)

    END SUBROUTINE rotate_scalar

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE rotate_vec(x, y, z, ax, ay, az)

      ! Purpose:
      !   to perform the 3-2-1 intrinsic rotation of a point whose initial position is given by "x", "y" and "z", according to angles
      !   "ax", "ay" and "az" respect to the x-, y- and z-axis. Angles must be given in radians. Axis and angles follow the righ-hand
      !   convention (x is north, y is east, z points downwards).
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   08/03/21                  original version
      !

      REAL(r__),   DIMENSION(:),  INTENT(INOUT) :: x, y, z
      REAL(r__),                  INTENT(IN)    :: ax, ay, az
      INTEGER(i32)                              :: i
      REAL(r__),   DIMENSION(3)                 :: v
      REAL(r__),   DIMENSION(3,3)               :: mx, my, mz

      !-----------------------------------------------------------------------------------------------------------------------------

      mx(:,1) = [1._r32, 0._r32, 0._r32]
      mx(:,2) = [0._r32, COS(ax), SIN(ax)]
      mx(:,3) = [0._r32, -SIN(ax), COS(ax)]

      my(:,1) = [COS(ay), 0._r32, -SIN(ay)]
      my(:,2) = [0._r32, 1._r32, 0._r32]
      my(:,3) = [SIN(ay), 0._r32, COS(ay)]

      mz(:,1) = [COS(az), SIN(az), 0._r32]
      mz(:,2) = [-SIN(az), COS(az), 0._r32]
      mz(:,3) = [0._r32, 0._r32, 1._r32]

      DO i = 1, SIZE(x)

        v = [x(i), y(i), z(i)]

        v = MATMUL(mx, v)
        v = MATMUL(my, v)
        v = MATMUL(mz, v)

        x(i) = v(1)
        y(i) = v(2)
        z(i) = v(3)

      ENDDO

    END SUBROUTINE rotate_vec

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *



    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *



    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *


END MODULE m_compgeo
