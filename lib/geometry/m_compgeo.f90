MODULE m_compgeo

  USE                :: omp_lib
  USE, NON_INTRINSIC :: m_precisions

  IMPLICIT none

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  PRIVATE

  PUBLIC :: is_triangle_intersection, rotate, trinterpolate, normal2tri

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

  INTERFACE rotate
    MODULE PROCEDURE rotate_scalar, rotate_vector
  END INTERFACE rotate

  INTERFACE cplus

    FUNCTION tritri(v0, v1, v2, u0, u1, u2) BIND(c, name="NoDivTriTriIsect")
      USE, INTRINSIC     :: iso_c_binding
      USE, NON_INTRINSIC :: m_precisions
! #ifdef DOUBLE_PREC
!       INTEGER, PARAMETER :: c_r__ = c_r64
! #else
!       INTEGER, PARAMETER :: c_r__ = c_r32
! #endif
      REAL(c_r64),    DIMENSION(3), INTENT(IN) :: v0, v1, v2, u0, u1, u2
      LOGICAL(c_bool)                          :: tritri
    END FUNCTION tritri

  END INTERFACE cplus

  CONTAINS

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    LOGICAL FUNCTION is_triangle_intersection(v0, v1, v2, u0, u1, u2)

      ! Purpose:
      !   to test whether two triangles, defined by vertices "v0", "v1", "v2", "u0", "u1" and "u2", intersect each other.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   08/03/21                  original version
      !

      REAL(r__), DIMENSION(3), INTENT(IN) :: v0, v1, v2
      REAL(r__), DIMENSION(3), INTENT(IN) :: u0, u1, u2

      !-----------------------------------------------------------------------------------------------------------------------------

      ! "tritri" works in always in double precision
#ifdef DOUBLE_PREC
      is_triangle_intersection = tritri(v0, v1, v2, u0, u1, u2)
#else
      is_triangle_intersection = tritri(REAL(v0, r64), REAL(v1, r64), REAL(v2, r64), REAL(u0, r64), REAL(u1, r64), REAL(u2, r64))
#endif

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

      REAL(r__),                  INTENT(INOUT) :: x, y, z
      REAL(r__),                  INTENT(IN)    :: ax, ay, az
      REAL(r__),   DIMENSION(3)                 :: v
      REAL(r__),   DIMENSION(3,3)               :: mx, my, mz

      !-----------------------------------------------------------------------------------------------------------------------------

      mx(:,1) = [1._r__, 0._r__, 0._r__]
      mx(:,2) = [0._r__, COS(ax), SIN(ax)]
      mx(:,3) = [0._r__, -SIN(ax), COS(ax)]

      my(:,1) = [COS(ay), 0._r__, -SIN(ay)]
      my(:,2) = [0._r__, 1._r__, 0._r__]
      my(:,3) = [SIN(ay), 0._r__, COS(ay)]

      mz(:,1) = [COS(az), SIN(az), 0._r__]
      mz(:,2) = [-SIN(az), COS(az), 0._r__]
      mz(:,3) = [0._r__, 0._r__, 1._r__]

      v = [x, y, z]

      v = MATMUL(mx, v)
      v = MATMUL(my, v)
      v = MATMUL(mz, v)

      x = v(1)
      y = v(2)
      z = v(3)

    END SUBROUTINE rotate_scalar

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE rotate_vector(x, y, z, ax, ay, az)

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

      mx(:,1) = [1._r__, 0._r__, 0._r__]
      mx(:,2) = [0._r__, COS(ax), SIN(ax)]
      mx(:,3) = [0._r__, -SIN(ax), COS(ax)]

      my(:,1) = [COS(ay), 0._r__, -SIN(ay)]
      my(:,2) = [0._r__, 1._r__, 0._r__]
      my(:,3) = [SIN(ay), 0._r__, COS(ay)]

      mz(:,1) = [COS(az), SIN(az), 0._r__]
      mz(:,2) = [-SIN(az), COS(az), 0._r__]
      mz(:,3) = [0._r__, 0._r__, 1._r__]

      DO i = 1, SIZE(x)

        v = [x(i), y(i), z(i)]

        v = MATMUL(mx, v)
        v = MATMUL(my, v)
        v = MATMUL(mz, v)

        x(i) = v(1)
        y(i) = v(2)
        z(i) = v(3)

      ENDDO

    END SUBROUTINE rotate_vector

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE tribar(p, a, b, c, bar)

      ! Purpose:
      !   to return the baricentric coordinates "bar" of point "p" inside a triangle determined by vertices "a/b/c".
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   08/03/21                  original version
      !

      REAL(r__),   DIMENSION(3), INTENT(IN)  :: p, a, b, c
      REAL(r__),   DIMENSION(3), INTENT(OUT) :: bar
      INTEGER(i32)                           :: i
      REAL(r__)                              :: d00, d01, d11, d20, d21, den, u, v, w
      REAL(r__),   DIMENSION(3)              :: v0, v1, v2

      !-----------------------------------------------------------------------------------------------------------------------------

      DO i = 1, 3
        v0(i) = b(i) - a(i)
        v1(i) = c(i) - a(i)
        v2(i) = p(i) - a(i)
      ENDDO

      d00 = v0(1)*v0(1) + v0(2)*v0(2) + v0(3)*v0(3)
      d01 = v0(1)*v1(1) + v0(2)*v1(2) + v0(3)*v1(3)
      d11 = v1(1)*v1(1) + v1(2)*v1(2) + v1(3)*v1(3)
      d20 = v2(1)*v0(1) + v2(2)*v0(2) + v2(3)*v0(3)
      d21 = v2(1)*v1(1) + v2(2)*v1(2) + v2(3)*v1(3)

      den = d00 * d11 - d01 * d01

      v = (d11 * d20 - d01 * d21) / den
      w = (d00 * d21 - d01 * d20) / den
      u = 1._r__ - v - w

      bar = [u, v, w]

    END SUBROUTINE tribar

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    REAL(r__) FUNCTION trinterpolate(p, a, b, c, f)

      ! Purpose:
      !   to interpolate a field "f" defined at the vertices "a/b/c" of a triangle at point "p" inside the triangle itself.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   08/03/21                  original version
      !

      REAL(r__), DIMENSION(3), INTENT(IN) :: p, a, b, c, f
      REAL(r__), DIMENSION(3)             :: u

      !-----------------------------------------------------------------------------------------------------------------------------

      CALL tribar(p, a, b, c, u)

      trinterpolate = DOT_PRODUCT(u, f)

    END FUNCTION trinterpolate

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE normal2tri(x, y, z, n)

      ! Purpose:
      !   to return the (unit) normal to a triangle whose vertices coordinates are "x(1),y(1),z(1)", "x(2),y(2),z(2)" and
      !   "x(3),y(3),z(3)"
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   08/03/21                  original version
      !

      REAL(r__),   DIMENSION(3), INTENT(IN)  :: x, y, z
      REAL(r__),   DIMENSION(3), INTENT(OUT) :: n
      INTEGER(i32)                           :: i
      REAL(r__)                              :: v

      !-----------------------------------------------------------------------------------------------------------------------------

      n(1) = ( (y(2)-y(1)) * (z(3)-z(1)) - (y(3)-y(1)) * (z(2)-z(1)) )
      n(2) = ( (z(2)-z(1)) * (x(3)-x(1)) - (x(2)-x(1)) * (z(3)-z(1)) )
      n(3) = ( (x(2)-x(1)) * (y(3)-y(1)) - (x(3)-x(1)) * (y(2)-y(1)) )

      v = NORM2(n)

      DO i = 1, 3
        n(i) = n(i) / v
      ENDDO

    END SUBROUTINE normal2tri

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    ! SUBROUTINE normal2tri(xo, yo, zo, no)
    !
    !   ! Purpose:
    !   !   to return the (unit) normal to a triangle whose vertices coordinates are "x(1),y(1),z(1)", "x(2),y(2),z(2)" and
    !   !   "x(3),y(3),z(3)"
    !   !
    !   ! Revisions:
    !   !     Date                    Description of change
    !   !     ====                    =====================
    !   !   08/03/21                  original version
    !   !
    !
    !   REAL(r32),   DIMENSION(3), INTENT(IN)  :: xo, yo, zo
    !   REAL(r32),   DIMENSION(3), INTENT(OUT) :: no
    !   INTEGER(i32)                           :: i
    !   REAL(r64)                              :: v
    !   REAL(r64),   DIMENSION(3)              :: x, y, z, n
    !
    !   !-----------------------------------------------------------------------------------------------------------------------------
    !
    !   x = xo
    !   y = yo
    !   z = zo
    !
    !   n(1) = ( (y(2)-y(1)) * (z(3)-z(1)) - (y(3)-y(1)) * (z(2)-z(1)) )
    !   n(2) = ( (z(2)-z(1)) * (x(3)-x(1)) - (x(2)-x(1)) * (z(3)-z(1)) )
    !   n(3) = ( (x(2)-x(1)) * (y(3)-y(1)) - (x(3)-x(1)) * (y(2)-y(1)) )
    !
    !   v = NORM2(n)
    !
    !   DO i = 1, 3
    !     n(i) = n(i) / v
    !   ENDDO
    !
    !   no = n
    !
    ! END SUBROUTINE normal2tri


    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *


END MODULE m_compgeo
