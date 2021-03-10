MODULE m_source

  USE, NON_INTRINSIC :: m_precisions
  USE, NON_INTRINSIC :: m_strings
  USE, NON_INTRINSIC :: m_parser
  USE, NON_INTRINSIC :: m_logfile
  USE, NON_INTRINSIC :: m_roughness
  USE, NON_INTRINSIC :: m_toolbox, ONLY: input, missing_arg

  IMPLICIT none

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  PRIVATE

  PUBLIC :: hypocenter, plane, setup_source

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  REAL(r32), PARAMETER :: PTSRC_FACTOR = 1._r32 / 10

  TYPE :: hyp
    INTEGER(i32) :: plane
    REAL(r32)    :: u, v
  END TYPE hyp

  TYPE :: src
    INTEGER(i32)                             :: nu, nv                !< points along strike and dip
    REAL(r32)                                :: length, width
    REAL(r32)                                :: xref, yref, zref      !< location of reference point
    REAL(r32)                                :: strike, dip
    REAL(r32)                                :: targetm0
    REAL(r32),   ALLOCATABLE, DIMENSION(:,:) :: slip, rise_time, rupture_time

  END TYPE src

  TYPE(hyp)                            :: hypocenter
  TYPE(src), ALLOCATABLE, DIMENSION(:) :: plane

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  CONTAINS

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE setup_source(ok, rank)

      ! Purpose:
      !   to setup the source parameters, either point-like or extended (from external file).
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   08/03/21                  original version
      !

      INTEGER(i32), INTENT(OUT) :: ok
      INTEGER(i32), INTENT(IN)  :: rank

      !-----------------------------------------------------------------------------------------------------------------------------

      IF (input%source%is_point) THEN



      ELSE
        IF (rank .eq. 0) CALL read_fsp_file(ok)
      ENDIF

    END SUBROUTINE setup_source

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE read_fsp_file(ok)

      INTEGER(i32), INTENT(OUT) :: ok
      CHARACTER(256)            :: buffer
      INTEGER(i32)              :: lu, i, n, ntw, nsbfs, ios
      REAL(r32)                 :: strike, dip, rake, m0, du, dv, length, width
      REAL(r32),  ALLOCATABLE, DIMENSION(:) :: numeric

      !-----------------------------------------------------------------------------------------------------------------------------

      ok = 0

      OPEN(newunit = lu, file = TRIM(input%source%file), status = 'old', form = 'formatted', access = 'sequential',  &
           action = 'read', iostat = ok)

      IF (ok .ne. 0) THEN
        CALL report_error('Error while opening file' + TRIM(input%source%file))
        RETURN
      ENDIF

      CALL parse(ok, n, lu, 'Nsg', ['=', ' '], '% Invs')    !< number of planes
      CALL missing_arg(ok, is_empty(n), 'Number of fault segments not found')

      IF (ok .ne. 0) RETURN

      ALLOCATE(plane(n))

      CALL parse(ok, ntw, lu, 'Ntw', ['=', ' '], '% Invs')     !< number of time-windows
      CALL parse(ok, du, lu, 'Dx', ['=', ' '], '% Invs')       !< subfault size along strike
      CALL parse(ok, dv, lu, 'Dz', ['=', ' '], '% Invs')       !< subfault size down-dip

      CALL parse(ok, hypocenter%u, lu, '(X)', ['=', ' '], '%')         !< hypocenter position along-strike
      CALL parse(ok, hypocenter%v, lu, '(Z)', ['=', ' '], '%')         !< hypocenter position down-dip
      CALL parse(ok, hypocenter%plane, lu, 'SEG', ['#', ' '], '%')     !< hypocenter position along-strike

      CALL parse(ok, length, lu, 'LEN', ['=',' '], '% Size')            !< total fault length
      CALL parse(ok, width, lu, 'WID', ['=',' '], '% Size')             !< fault width
      CALL parse(ok, m0, lu, 'Mo', ['=', ' '], '% Size')                !< total moment

      plane(:)%length   = length
      plane(:)%width    = width
      plane(:)%targetm0 = m0

      roughness%length = length

      CALL parse(ok, strike, lu, 'STRK', ['=', ' '], '% Mech')
      CALL parse(ok, dip, lu, 'DIP', ['=', ' '], '% Mech')
      CALL parse(ok, rake, lu, 'RAKE', ['=', ' '], '% Mech')

      plane(:)%strike = strike
      plane(:)%dip    = dip

      ! cycle over planes to override some parameters (e.g. strike, length) and assign others (slip, z)
      DO i = 1, SIZE(plane)

        CALL parse(ok, strike, lu, 'STRIKE', ['=', ' '], '% SEGMENT', i)
        CALL parse(ok, dip, lu, 'DIP', ['=', ' '], '% SEGMENT', i)

        IF (.not.is_empty(strike)) plane(i)%strike = strike
        IF (.not.is_empty(dip))    plane(i)%dip    = dip

        CALL parse(ok, length, lu, 'LEN', ['=', ' '], '%', nfield = i + 2)
        CALL parse(ok, width, lu, 'WID', ['=', ' '], '%', nfield = i + 1)

        IF (.not.is_empty(length)) plane(i)%length = length
        IF (.not.is_empty(width))  plane(i)%width  = width

        ! determine number of subfaults
        plane(i)%nu = NINT(plane(i)%length / du)
        plane(i)%nv = NINT(plane(i)%width / dv)

        ! read explicitly number of subfaults and raise error for mismatch
        CALL parse(ok, nsbfs, lu, 'Nsbfs', ['=', ' '], '%', nfield = i)

        IF (nsbfs .ne. plane(i)%nu*plane(i)%nv) THEN
          CALL report_error('Number of subfaults mismatching for fault plane ' + num2char(i))
          ok = 1
          RETURN
        ENDIF

        ! add a couple of extra points to the number of subfaults
        plane(i)%nu = plane(i)%nu + 2
        plane(i)%nv = plane(i)%nv + 2


! down here
        DO
          READ(lu, '(A256)', iostat = ios) buffer

          IF (ios .eq. -1) EXIT                       !< EOF reached: quit loop

          CALL parse_columns(buffer, numeric)

          IF (ALLOCATED(numeric)) print*, numeric
        ENDDO
        REWIND(lu)

      ENDDO


    END SUBROUTINE read_fsp_file

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE broadcast()


    END SUBROUTINE broadcast

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *


END MODULE m_source
