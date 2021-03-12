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
  REAL(r32), PARAMETER :: PI = 3.14159265358979323846_r64
  REAL(r32), PARAMETER :: DEG_TO_RAD = PI / 180

  TYPE :: hyp
    INTEGER(i32) :: plane        !< plane number where hypocenter is located
    REAL(r32)    :: u, v         !< on-fault coordinates
  END TYPE hyp

  TYPE :: src
    INTEGER(i32)                             :: nu, nv                             !< points along strike and dip (segment-specific)
    REAL(r32)                                :: length, width, strike, dip         !< segment-specific
    REAL(r32)                                :: targetm0                           !< total moment (sum over all segments)
    REAL(r32)                                :: x, y, z                            !< position of u=v=0
    REAL(r32),   ALLOCATABLE, DIMENSION(:,:) :: rise, rupture, sslip, dslip
    REAL(r32),   ALLOCATABLE, DIMENSION(:,:) :: u, v                               !< on-fault coordinates
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

      ! Purpose:
      !   to read the source properties specified in a FSP file.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   08/03/21                  original version
      !

      INTEGER(i32),                            INTENT(OUT) :: ok
      CHARACTER(256)                                       :: buffer
      INTEGER(i32)                                         :: lu, i, j, k, n, ntw, nsbfs, ios, index, islip, irake, irupt
      INTEGER(i32)                                         :: irise, iplane, ilon, ilat, iz
      REAL(r32)                                            :: strike, dip, rake, m0, du, dv, length, width, ew, ns
      REAL(r32),     ALLOCATABLE, DIMENSION(:)             :: numeric

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

      CALL parse(ok, strike, lu, 'STRK', ['=', ' '], '% Mech')
      CALL parse(ok, dip, lu, 'DIP', ['=', ' '], '% Mech')
      CALL parse(ok, rake, lu, 'RAKE', ['=', ' '], '% Mech')

      rake = rake * DEG_TO_RAD

      ! assign some properties to all segments (required if we have one segment only)
      plane(:)%strike = strike
      plane(:)%dip    = dip
      plane(:)%length = length
      plane(:)%width  = width
      plane(:)%targetm0 = m0

      ! cycle over planes to override some parameters (e.g. strike, length) and assign others (slip, z)
      DO k = 1, SIZE(plane)

        CALL parse(ok, strike, lu, 'STRIKE', ['=', ' '], '% SEGMENT', k)
        CALL parse(ok, dip, lu, 'DIP', ['=', ' '], '% SEGMENT', k)

        ! overwrite if we have segment-specific values
        IF (.not.is_empty(strike)) plane(k)%strike = strike
        IF (.not.is_empty(dip))    plane(k)%dip    = dip

        CALL parse(ok, length, lu, 'LEN', ['=', ' '], '%', nfield = k + 2)
        CALL parse(ok, width, lu, 'WID', ['=', ' '], '%', nfield = k + 1)

        ! ...as above
        IF (.not.is_empty(length)) plane(k)%length = length
        IF (.not.is_empty(width))  plane(k)%width  = width

        ! determine number of subfaults
        plane(k)%nu = NINT(plane(k)%length / du)
        plane(k)%nv = NINT(plane(k)%width / dv)

        ! read explicitly number of subfaults and raise error for mismatch
        CALL parse(ok, nsbfs, lu, 'Nsbfs', ['=', ' '], '%', nfield = k)

        IF (nsbfs .ne. plane(k)%nu*plane(k)%nv) THEN
          CALL report_error('Number of subfaults mismatching for fault plane ' + num2char(k))
          ok = 1
          RETURN
        ENDIF

        ! add "edge" points: here slip will be set to zero
        plane(k)%nu = plane(k)%nu + 2
        plane(k)%nv = plane(k)%nv + 2

        ALLOCATE(plane(k)%rise(plane(k)%nu, plane(k)%nv), plane(k)%rupture(plane(k)%nu, plane(k)%nv))
        ALLOCATE(plane(k)%sslip(plane(k)%nu, plane(k)%nv), plane(k)%dslip(plane(k)%nu, plane(k)%nv))

        plane(k)%slip(:) = 0._r32

        iplane = 0

        ! read and store physical quantities associated to current subfault
        DO
          READ(lu, '(A256)', iostat = ios) buffer

          IF (ios .eq. -1) EXIT                       !< EOF reached: quit loop

          IF (INDEX(buffer, 'Y==NS') .ne. 0) THEN

            iplane = iplane + 1           !< need to count plane number because "parse" rewind file

            IF (iplane .eq. k) THEN

              ilon  = parse_split_index(buffer, 'LON', ' ') - 1
              ilat  = parse_split_index(buffer, 'LAT', ' ') - 1
              iz    = parse_split_index(buffer, 'Z', ' ') - 1
              islip = parse_split_index(buffer, 'SLIP', ' ') - 1
              irake = parse_split_index(buffer, 'RAKE', ' ') - 1
              irupt = parse_split_index(buffer, 'RUPT', ' ') - 1
              irise = parse_split_index(buffer, 'RISE', ' ') - 1

              READ(lu, '(A256)', iostat = ios) buffer

              DO j = 2, plane(k)%nv - 1
                DO i = 2, plane(k)%nu - 1

                  READ(lu, '(A256)', iostat = ios) buffer

                  CALL parse_split(buffer, numeric)

                  ! get position (top-center) of first subfault: in our system, this is were u=v=0
                  ! move from geographical to utm coordinates (y is "east-west" as lon, x is "north-south" as lat)
                  IF (i .eq. 2 .and. j .eq. 2) THEN
                    plane(k)%z = numeric(iz)
                    CALL geo2utm(numeric(ilon), numeric(ilat), input%origin%lon, input%origin%lat, plane(i)%y, plane(i)%x)
                  ENDIF

                  IF (irake .gt. 0) rake = numeric(irake) * DEG_TO_RAD

                  ! get along-strike and down-dip components of slip
                  plane(k)%sslip(i, j) = numeric(islip) * COS(rake)
                  plane(k)%dslip(i, j) = numeric(islip) * -SIN(rake)          !< in our system slip is positive down-dip

                  IF (irupt .gt. 0) plane(k)%rupture(i, j) = numeric(irupt)
                  IF (irise .gt. 0) plane(k)%rise(i, j) = numeric(irise)

                ENDDO
              ENDDO

              EXIT

            ENDIF

          ENDIF

        ENDDO

        REWIND(lu)

        ALLOCATE(plane(k)%u(plane(k)%nu), plane(k)%v(plane(k)%nv))

        ! fill arrays with on-fault coordinates U-V: origin is set at upper left edge
        DO i = 2, plane(k)%nu - 1
          !plane(k)%u(i) =
        ENDDO

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
