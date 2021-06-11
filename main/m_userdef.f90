MODULE m_userdef

  USE, NON_INTRINSIC :: m_precisions

  IMPLICIT none

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  PRIVATE

  PUBLIC :: read_usr, write_usr

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  CONTAINS

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE read_usr(ok, lu, time, x, y, z)

      ! Purpose:
      !   to read an arbitrary ASCII file, connected as unit "lu", in order to provide three-component long-period timeseries stored
      !   in an arbitrary format. A vector representing time samples ("time") must be returned as well. On exit, "ok" is not zero if
      !   an error occurred.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   08/03/21                  original version
      !

      INTEGER(i32),                            INTENT(OUT) :: ok
      INTEGER(i32),                            INTENT(IN)  :: lu
      REAL(r32),    ALLOCATABLE, DIMENSION(:), INTENT(OUT) :: time, x, y, z

      !-----------------------------------------------------------------------------------------------------------------------------

      ok = 0

    END SUBROUTINE read_usr

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE write_usr(ok, lu, time, x, y, z)

      ! Purpose:
      !   to write an ASCII file, connected as unit "lu", having arbitrary format. Three-component long-period timeseries ("x", "y"
      !   and "z") and a time vector ("time") can be written to the file. On exit, "ok" is not zero if an error occurred while
      !   writing.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   08/03/21                  original version
      !

      INTEGER(i32),               INTENT(OUT) :: ok
      INTEGER(i32),               INTENT(IN)  :: lu
      REAL(r32),    DIMENSION(:), INTENT(IN)  :: time, x, y, z

      !-----------------------------------------------------------------------------------------------------------------------------

      ok = 0

    END SUBROUTINE write_usr

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *


END MODULE m_userdef
