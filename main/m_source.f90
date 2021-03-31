MODULE m_source

  USE, NON_INTRINSIC :: m_precisions
  USE, NON_INTRINSIC :: m_strings
  USE, NON_INTRINSIC :: m_parser
  USE, NON_INTRINSIC :: m_logfile
  !USE, NON_INTRINSIC :: m_roughness
  USE, NON_INTRINSIC :: m_toolbox, ONLY: input, geo2utm, missing_arg
#ifdef MPI
  USE                :: mpi
#endif

  IMPLICIT none

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  PRIVATE

  PUBLIC :: hypocenter, plane
  PUBLIC :: setup_source, meshing, dutr, dvtr

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  INTEGER(i32)         :: nutr, nvtr                        !< triangles along u, v
  REAL(r32)            :: dutr, dvtr, umingr, vmingr        !< triangle base/height, minimum u,v values for triangular mesh

  REAL(r32), PARAMETER :: PTSRC_FACTOR = 1._r32 / 10._r32
  REAL(r32), PARAMETER :: PI = 3.14159265358979323846_r64
  REAL(r32), PARAMETER :: DEG_TO_RAD = PI / 180._r32

  TYPE :: hyp
    INTEGER(i32) :: plane        !< plane number where hypocenter is located
    REAL(r32)    :: u, v         !< on-fault coordinates
  END TYPE hyp

  TYPE :: seg
    REAL(r32), ALLOCATABLE, DIMENSION(:) :: mrf
  END TYPE seg

  TYPE :: src
    REAL(r32)                                :: length, width, strike, dip         !< segment-specific
    REAL(r32)                                :: targetm0                           !< total moment (sum over all segments)
    REAL(r32)                                :: x, y, z                            !< absolute position of u=v=0
    REAL(r32),   ALLOCATABLE, DIMENSION(:)   :: u, v                               !< on-fault coordinates
    REAL(r32),   ALLOCATABLE, DIMENSION(:,:) :: rise, rupture, sslip, dslip, tslip
    TYPE(seg),   ALLOCATABLE, DIMENSION(:,:) :: src
  END TYPE src

  TYPE(hyp)                            :: hypocenter
  TYPE(src), ALLOCATABLE, DIMENSION(:) :: plane

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  CONTAINS

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE setup_source(ok, rank, ntasks)

      ! Purpose:
      !   to setup the source parameters, either point-like or extended (from external file).
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   08/03/21                  original version
      !

      INTEGER(i32), INTENT(OUT) :: ok
      INTEGER(i32), INTENT(IN)  :: rank, ntasks
      INTEGER(i32)              :: i, ierr
      !REAL(r32)                 :: du, dt, beta

      !-----------------------------------------------------------------------------------------------------------------------------

      IF (input%source%is_point) THEN

        ! for point-sources, set a minimum amount of common parameters... other mesh-specific parameters will be set later on in
        ! subroutine "meshing"

        ALLOCATE(plane(1))

        plane(1)%targetm0 = input%source%m0
        plane(1)%strike   = input%source%strike
        plane(1)%dip      = input%source%dip

        ALLOCATE(plane(1)%u(3), plane(1)%v(3))
        ALLOCATE(plane(1)%rise(3,3), plane(1)%sslip(3,3), plane(1)%dslip(3,3), plane(1)%rupture(3,3))

        plane(1)%sslip(:,:) = 0._r32
        plane(1)%dslip(:,:) = 0._r32

        ! slip in this context does not matter, as we will rescale output according to input "m0"
        plane(1)%sslip(2, 2) = COS(input%source%rake * DEG_TO_RAD)
        plane(1)%dslip(2, 2) = -SIN(input%source%rake * DEG_TO_RAD)          !< in our system slip is positive down-dip

        plane(1)%rise(:,:) = input%source%freq

        ! this is where u=v=0 (i.e. middle of fault)
        plane(1)%x = input%source%x
        plane(1)%y = input%source%y
        plane(1)%z = input%source%z

      ELSE

        IF (rank .eq. 0) THEN
          IF (INDEX(TRIM(input%source%file), '.fsp') .ne. 0) THEN
            CALL read_fsp_file(ok)
          ELSEIF (INDEX(TRIM(input%source%file), '.srf') .ne. 0) THEN
            CALL read_srf_file(ok)
          ENDIF
        ENDIF

#ifdef MPI
        CALL mpi_bcast(ok, 1, mpi_int, 0, mpi_comm_world, ierr)
#endif

        IF (ok .ne. 0) RETURN

#ifdef MPI
        CALL broadcast()
#endif

      ENDIF

      IF (rank .eq. ntasks - 1) CALL echo()

    END SUBROUTINE setup_source

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE echo()

      ! Purpose:
      !   to report source parameters to standard output.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   08/03/21                  original version
      !

      USE, INTRINSIC     :: iso_fortran_env, only: compiler_version

      INTEGER(i32) :: i

      !---------------------------------------------------------------------------------------------------------------------------------

      IF (input%source%is_point) THEN

      ELSE

        CALL update_log(num2char('Number of fault segments', width=30, fill='.') + num2char(SIZE(plane), width=15, justify='r') + &
                        '|')

        IF (plane(1)%targetm0 .lt. 0._r32) THEN
          CALL update_log(num2char('Total moment (Nm)', width=30, fill='.') + num2char('NA', width=15, justify='r') + '|', &
                          blankline=.false.)
        ELSE
          CALL update_log(num2char('Total moment (Nm)', width=30, fill='.') +   &
                          num2char(plane(1)%targetm0, width=15, notation='s', precision=3, justify='r') + '|', blankline=.false.)
        ENDIF

        CALL update_log(num2char('Fault segment', width=30, fill='.') +   &
                        num2char('Length', width=15, justify='r') + '|' + &
                        num2char('Width', width=15, justify='r') + '|' + &
                        num2char('Strike', width=15, justify='r') + '|' + &
                        num2char('Dip', width=15, justify='r') + '|' + &
                        num2char('Max Slip', width=15, justify='r') + '|')

        DO i = 1, SIZE(plane)

          CALL update_log(num2char(i, width=30, justify='c')  +  &
                          num2char(plane(i)%length, width=15, notation='f', precision=2, justify='r') + '|' + &
                          num2char(plane(i)%width, width=15, notation='f', precision=2, justify='r') + '|' + &
                          num2char(plane(i)%strike, width=15, notation='f', precision=2, justify='r') + '|' + &
                          num2char(plane(i)%dip, width=15, notation='f', precision=2, justify='r') + '|' +  &
                          num2char(MAXVAL(euclid(plane(i)%sslip, plane(i)%dslip)), width=15, notation='f', precision=3,   &
                                   justify='r') + '|', blankline=.false.)
        ENDDO

      ENDIF

    END SUBROUTINE echo

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE read_fsp_file(ok)

      ! Purpose:
      !   to read the source properties specified in a FSP file (see http://equake-rc.info/SRCMOD/fileformats/fsp/ for a description).
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   08/03/21                  original version
      !

      INTEGER(i32),                            INTENT(OUT) :: ok
      CHARACTER(256)                                       :: buffer
      INTEGER(i32)                                         :: lu, i, j, k, n, nu, nv, ntw, nsbfs, ios, index, islip, irake, irupt
      INTEGER(i32)                                         :: irise, iplane, ilon, ilat, iz
      REAL(r32)                                            :: strike, dip, rake, m0, du, dv, length, width, ew, ns
      REAL(r32),     ALLOCATABLE, DIMENSION(:)             :: numeric

      !-----------------------------------------------------------------------------------------------------------------------------

      ok = 0

      OPEN(newunit = lu, file = TRIM(input%source%file), status = 'old', form = 'formatted', access = 'sequential',  &
           action = 'read', IOSTAT = ok)

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

      du = du * 1.E+03_r32        !km to m
      dv = dv * 1.E+03_r32

      ! try to read hypocenter position (on-fault u-v coordinates)
      IF (n .eq. 1) THEN
        CALL parse(ok, hypocenter%u, lu, 'HypX', ['=', ' '], '%')
        CALL parse(ok, hypocenter%v, lu, 'HypZ', ['=', ' '], '%')
        hypocenter%plane = 1
      ELSE
        CALL parse(ok, hypocenter%u, lu, '(X)', ['=', ' '], '%')
        CALL parse(ok, hypocenter%v, lu, '(Z)', ['=', ' '], '%')
        CALL parse(ok, hypocenter%plane, lu, 'SEG', ['#', ' '], '%')
      ENDIF

      CALL parse(ok, length, lu, 'LEN', ['=',' '], '% Size')            !< total fault length
      CALL parse(ok, width, lu, 'WID', ['=',' '], '% Size')             !< fault width
      CALL parse(ok, m0, lu, 'Mo', ['=', ' '], '% Size')                !< total moment

      CALL parse(ok, strike, lu, 'STRK', ['=', ' '], '% Mech')
      CALL parse(ok, dip, lu, 'DIP', ['=', ' '], '% Mech')
      CALL parse(ok, rake, lu, 'RAKE', ['=', ' '], '% Mech')

      rake = rake * DEG_TO_RAD

      IF (is_empty(m0)) m0 = -1._r32

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

        plane(k)%length = plane(k)%length * 1.E+03_r32   !< km to m
        plane(k)%width = plane(k)%width * 1.E+03_r32

        ! determine number of subfaults
        nu = NINT(plane(k)%length / du)
        nv = NINT(plane(k)%width / dv)

        ! read explicitly number of subfaults and raise error for mismatch
        CALL parse(ok, nsbfs, lu, 'Nsbfs', ['=', ' '], '%', nfield = k)

        IF (nsbfs .ne. nu*nv) THEN
          CALL report_error('Number of subfaults mismatching for fault plane ' + num2char(k))
          ok = 1
          RETURN
        ENDIF

        ! add "edge" points: here slip will be set to zero
        nu = nu + 2
        nv = nv + 2

        ALLOCATE(plane(k)%rise(nu, nv), plane(k)%rupture(nu, nv))
        ALLOCATE(plane(k)%sslip(nu, nv), plane(k)%dslip(nu, nv), plane(k)%tslip(nu, nv))

        DO j = 1, nv
          DO i = 1, nu
            plane(k)%sslip(i, j) = 0._r32
            plane(k)%dslip(i, j) = 0._r32
            plane(k)%tslip(i, j) = 0._r32
          ENDDO
        ENDDO

        iplane = 0

        ! read and store physical quantities associated to current subfault
        DO
          READ(lu, '(A256)', IOSTAT = ios) buffer

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

              READ(lu, '(A256)', IOSTAT = ios) buffer

              DO j = 2, nv - 1
                DO i = 2, nu - 1

                  READ(lu, '(A256)', IOSTAT = ios) buffer

                  CALL parse_split(buffer, numeric)

                  ! get position (top-center) of first subfault: in our system, this is were u=v=0
                  ! move from geographical to utm coordinates (y is "east-west" as lon, x is "north-south" as lat)
                  IF (i .eq. 2 .and. j .eq. 2) THEN
                    plane(k)%z = numeric(iz)
                    ! CALL geo2utm(numeric(ilon), numeric(ilat), input%origin%lon, input%origin%lat, plane(i)%y, plane(i)%x)
                  ENDIF

                  IF (irake .gt. 0) rake = numeric(irake) * DEG_TO_RAD

                  ! get along-strike and down-dip components of slip
                  plane(k)%sslip(i, j) = numeric(islip) * COS(rake)
                  plane(k)%dslip(i, j) = -numeric(islip) * SIN(rake)          !< in our system slip is positive down-dip
                  plane(k)%tslip(i, j) = numeric(islip)                       !< total slip

                  IF (irupt .gt. 0) plane(k)%rupture(i, j) = numeric(irupt)
                  IF (irise .gt. 0) plane(k)%rise(i, j) = numeric(irise)

                ENDDO
              ENDDO

              EXIT

            ENDIF

          ENDIF

        ENDDO

        REWIND(lu)

        ! at this point rise-time may be in seconds or empty: convert to "freq" or assign constant input value
        IF (irise .gt. 0) THEN
          plane(k)%rise(:,:) = 2._r32 * PI / plane(k)%rise(:,:)         !< assume "fc=1/rise"
        ELSE
          plane(k)%rise(:,:) = input%source%freq
        ENDIF

        ! rupture time may be empty... if so, set first element to negative value
        plane(k)%rupture(1, 1) = -999._r32

        ALLOCATE(plane(k)%u(nu), plane(k)%v(nv))

        ! in FSP geometry, reference point (u=v=0) is (first) top cell

        ! set on-fault u-coordinate: first and last point are edge points, second point correspond to reference point u=0
        plane(k)%u(1) = -du / 2._r32
        DO i = 2, nu - 1
          plane(k)%u(i) = (i - 2) * du
        ENDDO
        plane(k)%u(nu) = plane(k)%u(nu - 1) + du / 2._r32

        ! first and last point are edge points, first point correspond to reference point v=0
        plane(k)%v(1) = 0._r32
        DO i = 2, nv - 1
          plane(k)%v(i) = (i - 2) * dv + dv / 2._r32
        ENDDO
        plane(k)%v(nv) = plane(k)%v(nv - 1) + dv / 2._r32

        ! make sure uppermost down-dip "edge" point is always slightly (1m) below free-surface
        !plane(k)%z = MAX(1._r32, plane(k)%z + MAX(0._r32, SIN(plane(k)%dip * DEG_TO_RAD) * dv / 2._r32 - plane(k)%z))
        plane(k)%z = MAX(1._r32, plane(k)%z * 1.E+03_r32)

      ENDDO

      CLOSE(lu, IOSTAT = ok)

      IF (ok .ne. 0) THEN
        CALL report_error('Error while closing file ' + TRIM(input%source%file))
        RETURN
      ENDIF

    END SUBROUTINE read_fsp_file

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE read_srf_file(ok)

      ! Purpose:
      !   to read the source properties specified in a SRF file v2.0 (see http://scec.usc.edu/scecpedia/Standard_Rupture_Format for
      !   a description).
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   08/03/21                  original version
      !

      INTEGER(i32),                            INTENT(OUT) :: ok
      CHARACTER(256)                                       :: buffer
      INTEGER(i32)                                         :: n, lu, pl, lino, nt1, nt2, nt3, i, j, nu, nv, points
      REAL(r32)                                            :: version, lon, lat, z, strike, dip, area, tinit, dt, vs, density
      REAL(r32)                                            :: rake, slip1, slip2, slip3, shyp, dhyp

      !-----------------------------------------------------------------------------------------------------------------------------

      ok = 0

      OPEN(newunit = lu, file = TRIM(input%source%file), status = 'old', form = 'formatted', access = 'sequential',  &
           action = 'read', IOSTAT = ok)

      IF (ok .ne. 0) THEN
        CALL report_error('Error while opening file' + TRIM(input%source%file))
        RETURN
      ENDIF

      READ(lu, *, IOSTAT = ok) version

      IF (NINT(version) .ne. 2) THEN
        CALL report_error('read_srf_file - ERROR: file version is not 2')
        ok = 1
        RETURN
      ENDIF

      CALL parse(ok, n, lu, 'PLANE', [' ', ' '], com = '#')    !< number of planes

      IF (is_empty(n)) THEN
        CALL report_error('read_srf_file - ERROR: number of fault segments not found')
        ok = 1
        RETURN
      ENDIF

      IF (n .le. 0) THEN
        CALL report_error('read_srf_file - ERROR: number of fault segments cannot be 0 or less')
        ok = 1
        RETURN
      ENDIF

      ALLOCATE(plane(n))

      ! skip lines until we step onto the header
      DO
        READ(lu, *, IOSTAT = ok) buffer
        IF ( (INDEX(buffer, '#') .eq. 0) .and. (INDEX(buffer, 'PLANE') .gt. 0) ) EXIT
      ENDDO

      ! read header to find number of subfaults along strike and dip directions
      DO pl = 1, n

        READ(lu, *, IOSTAT = ok) lon, lat, nu, nv, plane(pl)%length, plane(pl)%width
        READ(lu, *, IOSTAT = ok) plane(pl)%strike, plane(pl)%dip, plane(pl)%z, shyp, dhyp

        nu = nu + 2          !< add edge points
        nv = nv + 2

        ALLOCATE(plane(pl)%rupture(nu, nv), plane(pl)%sslip(nu, nv), plane(pl)%dslip(nu, nv), plane(pl)%src(nu, nv))
        ALLOCATE(plane(pl)%rise(nu, nv), plane(pl)%u(nu), plane(pl)%v(nv), plane(pl)%tslip(nu, nv))

        DO j = 1, nv
          DO i = 1, nu
            plane(pl)%sslip(i, j) = 0._r32
            plane(pl)%dslip(i, j) = 0._r32
            plane(pl)%tslip(i, j) = 0._r32
          ENDDO
        ENDDO

      ENDDO

      ! read properties for each plane
      DO pl = 1, n

        nv = SIZE(plane(pl)%v)
        nu = SIZE(plane(pl)%u)

        plane(pl)%targetm0 = 0._r32

        READ(lu, *, IOSTAT = ok) buffer, points

        IF (points .ne. (nu-2)*(nv-2)) THEN
          CALL report_error('read_srf_file - ERROR: number of subfaults in header does not correspond to subfaults in data')
          ok = 1
          RETURN
        ENDIF

        DO j = 2, nv - 1
          DO i = 2, nu - 1

            READ(lu, *, IOSTAT = ok) lon, lat, z, strike, dip, area, tinit, dt, vs, density
            READ(lu, *, IOSTAT = ok) rake, slip1, nt1, slip2, nt2, slip3, nt3

            ! get position (top-center) of first subfault: in our system, this is were u=v=0
            ! move from geographical to utm coordinates (y is "east-west" as lon, x is "north-south" as lat)
            IF (i .eq. 2 .and. j .eq. 2) THEN
              !CALL geo2utm(lon, lat, input%origin%lon, input%origin%lat, plane(pl)%y, plane(pl)%x)
            ENDIF

            plane(pl)%targetm0 = plane(pl)%targetm0 + area * slip1 * density * vs**2        !< scalar moment (dyne*cm)

            rake = rake * DEG_TO_RAD              !< convert deg to rad
            slip1 = slip1 * 1.0e-02_r32           !< convert cm to m

            ! get along-strike and down-dip components of slip
            plane(pl)%sslip(i, j) = slip1 * COS(rake)
            plane(pl)%dslip(i, j) = -slip1 * SIN(rake)          !< in our system slip is positive down-dip
            plane(pl)%tslip(i, j) = slip1                      !< total slip

            plane(pl)%rupture(i, j) = tinit

            ALLOCATE(plane(pl)%src(i, j)%mrf(nt1 + 1))

            ! mrf values are arranged in columns
            DO lino = 1, nt1 / 6
              READ(lu, *, IOSTAT = ok) plane(pl)%src(i, j)%mrf((lino-1)*6 + 1: (lino-1)*6 + 6)
            ENDDO

            ! read remaining values
            IF (nt1 .ne. (lino-2)*6 + 6) READ(lu, *, IOSTAT = ok) plane(pl)%src(i, j)%mrf((lino-1)*6 + 1: nt1)

            plane(pl)%src(i, j)%mrf(nt1 + 1) = 0._r32     !< append a zero

            ! scale mrf such that its integral is one, convert from cm/s to m/s
            plane(pl)%src(i, j)%mrf = plane(pl)%src(i, j)%mrf * 1.0e-02_r32 / slip1

          ENDDO
        ENDDO

        plane(pl)%targetm0 = plane(pl)%targetm0 * 1.0e-07          !< convert from dyne*cm to N*m

        plane(pl)%targetm0 = MAX(plane(pl)%targetm0, -1._r32)

        ! make sure uppermost down-dip "edge" point is always slightly (1m) below free-surface
        !plane(k)%z = MAX(1._r32, plane(k)%z + MAX(0._r32, SIN(plane(k)%dip * DEG_TO_RAD) * dv / 2._r32 - plane(k)%z))
        plane(pl)%z = MAX(1._r32, plane(pl)%z)

      ENDDO

      CLOSE(lu, IOSTAT = ok)

      IF (ok .ne. 0) THEN
        CALL report_error('read_srf_file - ERROR: could not close file ' + TRIM(input%source%file))
        RETURN
      ENDIF

    END SUBROUTINE read_srf_file

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE broadcast()

      ! Purpose:
      !   to broadcast all source parameters to all processes. Make use of some manual packing/unpacking to minimize mpi calls.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   08/03/21                  original version
      !

#ifdef MPI

      INTEGER(i32)                            :: ierr, i, n
      INTEGER(i32),              DIMENSION(2) :: np
      INTEGER(i32), ALLOCATABLE, DIMENSION(:) :: intg
      REAL(r32),    ALLOCATABLE, DIMENSION(:) :: float

      !-----------------------------------------------------------------------------------------------------------------------------

      intg = [1, hypocenter%plane]

      IF (ALLOCATED(plane)) intg = [SIZE(plane), hypocenter%plane]

      CALL mpi_bcast(intg, 2, mpi_int, 0, mpi_comm_world, ierr)

      n = intg(1)

      hypocenter%plane = intg(2)

      float = [hypocenter%u, hypocenter%v]

      CALL mpi_bcast(float, 2, mpi_real, 0, mpi_comm_world, ierr)

      hypocenter%u = float(1)
      hypocenter%v = float(2)

      IF (.not.ALLOCATED(plane)) ALLOCATE(plane(n))

      float = [plane(:)%length, plane(:)%width, plane(:)%strike, plane(:)%dip, plane(:)%targetm0, plane(:)%x, plane(:)%y,   &
               plane(:)%z]

      CALL mpi_bcast(float, SIZE(float), mpi_real, 0, mpi_comm_world, ierr)

      plane(:)%length = float(1:n)
      plane(:)%width  = float(n+1:2*n)
      plane(:)%strike = float(2*n+1:3*n)
      plane(:)%dip    = float(3*n+1:4*n)
      plane(:)%targetm0 = float(4*n+1:5*n)
      plane(:)%x = float(5*n+1:6*n)
      plane(:)%y = float(6*n+1:7*n)
      plane(:)%z = float(7*n+1:8*n)

      DO i = 1, SIZE(plane)

        IF (ALLOCATED(plane(i)%u)) np = [SIZE(plane(i)%u), SIZE(plane(i)%v)]

        CALL mpi_bcast(np, 2, mpi_int, 0, mpi_comm_world, ierr)

        IF (.not.ALLOCATED(plane(i)%u)) THEN
          ALLOCATE(plane(i)%u(np(1)), plane(i)%v(np(2)), plane(i)%sslip(np(1),np(2)), plane(i)%dslip(np(1),np(2)))
          ALLOCATE(plane(i)%rise(np(1),np(2)), plane(i)%rupture(np(1),np(2)))
        ENDIF

        float = [plane(i)%u, plane(i)%v]

        CALL mpi_bcast(float, SIZE(float), mpi_real, 0, mpi_comm_world, ierr)

        plane(i)%u = float(1:np(1))
        plane(i)%v = float(np(1)+1:np(1)+np(2))

        CALL mpi_bcast(plane(i)%sslip, PRODUCT(np), mpi_real, 0, mpi_comm_world, ierr)
        CALL mpi_bcast(plane(i)%dslip, PRODUCT(np), mpi_real, 0, mpi_comm_world, ierr)
        CALL mpi_bcast(plane(i)%rise, PRODUCT(np), mpi_real, 0, mpi_comm_world, ierr)
        CALL mpi_bcast(plane(i)%rupture, PRODUCT(np), mpi_real, 0, mpi_comm_world, ierr)

      ENDDO

#endif

    END SUBROUTINE broadcast

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    REAL(r32) ELEMENTAL FUNCTION euclid(x, y)

      ! Purpose:
      !   to compute euclidean distance without under/overflow.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   08/03/21                  original version
      !

      REAL(r32), INTENT(IN) :: x, y

      !-----------------------------------------------------------------------------------------------------------------------------

      euclid = HYPOT(x, y)

    END FUNCTION euclid

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    REAL(r32) ELEMENTAL FUNCTION brune(t, freq)

      ! Purpose:
      !   to compute the brune source time function (equal to 1 for "t" -> +inf) as defined in Eq.6 of Beresnev&Atkinson, 1997. Note
      !   that "freq=2*pi*fc", where "fc" is the corner frequency.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   08/03/21                  original version
      !

      REAL(r32), INTENT(IN) :: t, freq
      REAL(r32)             :: x

      !-----------------------------------------------------------------------------------------------------------------------------

      IF (t .le. 0._r32) THEN
        brune = 0._r32
      ELSE
        x = freq * t
        brune = 1._r32 - EXP(-x) * (1._r32 + x)
      ENDIF

    END FUNCTION brune

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    REAL(r32) ELEMENTAL FUNCTION dbrune(t, freq)

      ! Purpose:
      !   to compute the brune moment rate time function ("dbrune"-> 0 for "t" -> +inf) as defined in Eq.8 of Beresnev&Atkinson,
      !   1997. Note that "freq=2*pi*fc", where "fc" is the corner frequency.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   08/03/21                  original version
      !

      REAL(r32), INTENT(IN) :: t, freq
      REAL(r32)             :: x

      !-----------------------------------------------------------------------------------------------------------------------------

      IF (t .le. 0._r32) THEN
        dbrune = 0._r32
      ELSE
        x = freq * t
        dbrune = freq * x * EXP(-x)
      ENDIF

    END FUNCTION dbrune

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE meshing(pl, vel)

      ! Purpose:
      !   to discretize the "pl"-th fault segment for the "vel"-th velocity model using a triangular mesh. Triangle size is
      !   determined by the number of points per minimum wavelength. Inside the mesh, triangles are organized as follows:
      !    _ _ _ _ _ _ _ _ _ _ _
      !    \  1 / \  3 / \  5 /
      !     \  / 2 \  / 4 \  /         -> j-th (odd) row
      !      \/__ _ \/_ _ _\/
      !     / \  2 / \  4 / \
      !    / 1 \  / 3 \  / 5 \         -> (j+1)-th (even) row
      !   /_ _ _\/ _ _ \/_ _ _\
      !
      ! Corners are numbered as follows, depending whether the triangle points downwards or upwards:
      !
      !   1_ _ _ 3      2
      !    \    /      / \
      !     \  /      /   \
      !      \/      /_ _ _\
      !      2      1       3
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   08/03/21                  original version
      !

      INTEGER(i32), INTENT(IN) :: pl, vel
      INTEGER(i32)             :: nu, nv, nugr, nvgr
      REAL(r32)                :: du, dv, dt, beta, umaxgr, vmaxgr

      !-----------------------------------------------------------------------------------------------------------------------------

      IF (input%source%is_point) THEN

        ASSOCIATE(model => input%velocity(vel))
          beta = vinterp(model%depth, model%vs, model%vsgrad, plane(pl)%z)           !< s-wave velocity at source depth
        END ASSOCIATE

        ! first condition for point-sources: lambda_min >> d, where d is diagonal of a square fault
        ! diagonal = l*sqrt(2) = 2*sqrt(2)*du (l = 2*du)
        du = PTSRC_FACTOR * beta / input%coda%fmax / SQRT(2._r32) / 2._r32

        plane(pl)%u = [-du, 0._r32, du]
        plane(pl)%v = [-dv, 0._r32, dv]

        ! second condition for point-source: max(rupture time) << tr, where tr is the shortes period observed (i.e. 1/fmax)
        dt = PTSRC_FACTOR / input%coda%fmax

        plane(pl)%rupture(:, 1) = dt * [1._r32, 1._r32/SQRT(2._r32), 1._r32]                   !< proceed along-strike (along u)
        plane(pl)%rupture(:, 2) = dt * [1._r32/SQRT(2._r32), 0._r32, 1._r32/SQRT(2._r32)]
        plane(pl)%rupture(:, 3) = dt * [1._r32, 1._r32/SQRT(2._r32), 1._r32]

        ! make sure uppermost down-dip "edge" point is always slightly (1m) below free-surface
        plane(pl)%z = MAX(1._r32, plane(pl)%z + MAX(0._r32, SIN(plane(pl)%dip * DEG_TO_RAD) * du - plane(pl)%z))

      ELSE

        du = mesh_spacing(pl, vel)

      ENDIF

      nu = SIZE(plane(pl)%u)
      nv = SIZE(plane(pl)%v)

      ! determine grid extension
      umingr = plane(pl)%u(1)
      umaxgr = plane(pl)%u(nu)
      vmingr = plane(pl)%v(1)
      vmaxgr = plane(pl)%v(nv)

      ! find roughly how many triangles with desired size we can cram onto the fault plane
      nutr = (umaxgr - umingr) / du  + 0.5_r32
      nvtr = (vmaxgr - vmingr) / du  + 0.5_r32

      ! gets actual triangle base (du) and heigth (dv)
      dutr = (umaxgr - umingr) / nutr
      dvtr = (vmaxgr - vmingr) / nvtr

      ! total number of triangles along a row
      ! nutr = 2*nutr - 1

      ! compute maximum number of vertexes. note: nugr is (nutr + 1) for odd rows and (nutr) for
      ! even rows
      nugr = nutr + 1
      nvgr = nvtr + 1

      IF (input%advanced%verbose .eq. 2) THEN

        CALL update_log(num2char('<grid extension>', justify='c', width=30) + num2char('Umin', width=15, justify='r') + '|' +  &
                        num2char('Umax', width=15, justify='r') + '|' + num2char('Vmin', width=15, justify='r') + '|' + &
                        num2char('Vmax', width=15, justify='r') + '|', blankline=.false.)

        CALL update_log(num2char('', width=30, justify='c')  +  &
                        num2char(umingr, width=15, notation='f', precision=2, justify='r') + '|' + &
                        num2char(umaxgr, width=15, notation='f', precision=2, justify='r') + '|' + &
                        num2char(vmingr, width=15, notation='f', precision=2, justify='r') + '|' + &
                        num2char(vmaxgr, width=15, notation='f', precision=2, justify='r') + '|', blankline=.false.)

        CALL update_log(num2char('<triangles>', justify='c', width=30) + num2char('Width', width=15, justify='r') + '|' + &
                        num2char('Heigth', width=15, justify='r') + '|' + num2char('Along u', width=15, justify='r') + '|' + &
                        num2char('Along v', width=15, justify='r') + '|' + num2char('Totals', width=15, justify='r') + '|',  &
                        blankline=.false.)

        CALL update_log(num2char('', width=30)  +  &
                        num2char(dutr, width=15, notation='f', precision=2, justify='r') + '|' + &
                        num2char(dvtr, width=15, notation='f', precision=2, justify='r') + '|' + &
                        num2char(nutr, width=15, justify='r') + '|' + &
                        num2char(nvtr, width=15, justify='r') + '|' +  &
                        num2char(nvtr*(2*nutr-1), width=15, justify='r') + '|', blankline=.false.)

      ENDIF

    END SUBROUTINE meshing

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    REAL(r32) FUNCTION mesh_spacing(pl, vel)

      ! Purpose:
      !   to compute mesh spacing (i.e. triangle side) for "pl"-th fault plane and "vel"-th velocity model such that the requirement
      !   of having "input%advanced%pmw" points per minimum wavelength is fulfilled.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   08/03/21                  original version
      !

      INTEGER(i32), INTENT(IN) :: pl, vel
      INTEGER(i32)             :: layer
      REAL(r32)                :: ztop, zmax, vmin

      !-----------------------------------------------------------------------------------------------------------------------------

      ztop = plane(pl)%z
      zmax = ztop + plane(pl)%width * SIN(plane(pl)%dip * DEG_TO_RAD)

      vmin = HUGE(0._r32)

      ! loop over layers to find minimum shear wave speed on fault plane
      ASSOCIATE(model => input%velocity(vel))

        DO layer = 1, SIZE(model%depth) - 1
          vmin = MIN(vmin, model%vsgrad(layer) * (ztop - model%depth(layer)) + model%vs(layer))    !< top
          vmin = MIN(vmin, model%vsgrad(layer) * (MIN(model%depth(layer + 1), zmax) - model%depth(layer)) + model%vs(layer))  !< bottom
          ztop = model%depth(layer + 1)
          IF (zmax .lt. ztop) EXIT
        ENDDO

        IF (zmax .ge. ztop) THEN
          vmin = MIN(vmin, model%vsgrad(layer) * (ztop - model%depth(layer)) + model%vs(layer))    !< top
          vmin = MIN(vmin, model%vsgrad(layer) * (zmax - model%depth(layer)) + model%vs(layer))    !< bottom
        ENDIF

      END ASSOCIATE

      mesh_spacing = vmin / input%coda%fmax / input%advanced%pmw

    END FUNCTION mesh_spacing

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    REAL(r32) FUNCTION time_stepping(pl, vel, du)

      ! Purpose:
      !   to compute integration time-step for "pl"-th fault plane, "vel"-th velocity model and grid-step "du" such that the
      !   requirement of having at least "input%advanced%avecuts" cuts per triangle is fulfilled.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   08/03/21                  original version
      !

      INTEGER(i32), INTENT(IN) :: pl, vel
      REAL(r32),    INTENT(IN) :: du
      INTEGER(i32)             :: layer
      REAL(r32)                :: ztop, zmax, vmax

      !-----------------------------------------------------------------------------------------------------------------------------

      ztop = plane(pl)%z + plane(pl)%v(1) * SIN(plane(pl)%dip * DEG_TO_RAD)     !< for point-source, "z" (u=v=0) is center, not top
      zmax = ztop + plane(pl)%width * SIN(plane(pl)%dip * DEG_TO_RAD)

      vmax = 0._r32

      ! loop over layers to find maximum shear wave speed on fault plane
      ASSOCIATE(model => input%velocity(vel))

        DO layer = 1, SIZE(model%depth) - 1
          vmax = MAX(vmax, model%vsgrad(layer) * (ztop - model%depth(layer)) + model%vs(layer))    !< top
          vmax = MAX(vmax, model%vsgrad(layer) * (MIN(model%depth(layer + 1), zmax) - model%depth(layer)) + model%vs(layer))  !< bottom
          ztop = model%depth(layer + 1)
          IF (zmax .lt. ztop) EXIT
        ENDDO

        IF (zmax .ge. ztop) THEN
          vmax = MAX(vmax, model%vsgrad(layer) * (ztop - model%depth(layer)) + model%vs(layer))    !< top
          vmax = MAX(vmax, model%vsgrad(layer) * (zmax - model%depth(layer)) + model%vs(layer))    !< bottom
        ENDIF

      END ASSOCIATE

      time_stepping = du / vmax / input%source%vrfact / input%advanced%avecuts

    END FUNCTION time_stepping

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    REAL(r32) FUNCTION vinterp(ztop, phys, gradient, z)

      ! Purpose:
      !   to return a physical property (e.g. velocity, density, etc) at depth "z" given a model where property "phys", gradient
      !   "gradient" and depth to layer "ztop" are defined for each layer.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   08/03/21                  original version
      !

      REAL(r32),   DIMENSION(:), INTENT(IN) :: ztop, phys, gradient
      REAL(r32),                 INTENT(IN) :: z
      INTEGER(i32)                          :: i

      !-----------------------------------------------------------------------------------------------------------------------------

      DO i = SIZE(ztop), 1, -1
        IF (z .ge. ztop(i)) THEN
          vinterp = phys(i) + gradient(i)*(z - ztop(i))
          EXIT
        ENDIF
      ENDDO

    END FUNCTION vinterp

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE cornr(irow, itri, iuc, ivc)

      ! Purpose:
      !   to return corner indices "iuc" (along u) and "ivc" (along v) for the "itri"-th triangle in the "irow"-th row. Corners are
      !   numbered sich that corner 1 has the least u value, corner 3 has the greatest and corner 2 is the intermediate.
      !
      ! Corners are numbered as follows, depending whether the triangle points downwards ("itri" odd) or upwards ("itri" even) and
      ! the "irow"-th row:
      !
      ! e.g. "irow"=1, "itri"=1 and 2
      !
      !    1_ _ _2        2                       1_ _ _1          1
      !    \    /        / \                      \    /          / \
      !     \  /        /   \       IUC            \  /          /   \      IVC
      !      \/        /_ _ _\                      \/          /_ _ _\
      !      1        1       2                     2          2       2
      !
      ! e.g. "irow"=2, "itri"=1 and 2
      !
      !       1       1_ _ _2                    1_ _ _1          1
      !      / \      \    /                     \    /          / \
      !     /   \      \  /         IUC           \  /          /   \      IVC
      !    /_ _ _\      \/                         \/          /_ _ _\
      !   1       2     2                          2          2       2
      !
      ! Corners indices for "iuc" and "ivc" are as follows, depending whether the triangle points downwards or upwards:
      !
      !   1_ _ _ 3      2
      !    \    /      / \
      !     \  /      /   \
      !      \/      /_ _ _\
      !      2      1       3
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   08/03/21                  original version
      !

      INTEGER(i32),               INTENT(IN)  :: irow, itri
      INTEGER(i32), DIMENSION(3), INTENT(OUT) :: iuc, ivc
      INTEGER(i32)                            :: i1, i2

      !-----------------------------------------------------------------------------------------------------------------------------

      i1 = (itri + 1) / 2
      i2 = itri / 2

      IF (MOD(irow, 2) .ne. 0) THEN               !< odd row
         IF (MOD(itri, 2) .ne. 0) THEN            !< odd triangle
            iuc(1) = i1
            ivc(1) = irow
            iuc(2) = i1
            ivc(2) = irow + 1
            iuc(3) = i1 + 1
            ivc(3) = irow
         ELSE
            iuc(1) = i2
            ivc(1) = irow + 1
            iuc(2) = i2 + 1
            ivc(2) = irow
            iuc(3) = i2 + 1
            ivc(3) = irow + 1
         ENDIF
      ELSE                                        !< even row
         IF (MOD(itri, 2) .ne. 0) THEN            !< odd triangle
            iuc(1) = i1
            ivc(1) = irow + 1
            iuc(2) = i1
            ivc(2) = irow
            iuc(3) = i1 + 1
            ivc(3) = irow + 1
         ELSE
            iuc(1) = i2
            ivc(1) = irow
            iuc(2) = i2 + 1
            ivc(2) = irow + 1
            iuc(3) = i2 + 1
            ivc(3) = irow
         ENDIF
      ENDIF

    END SUBROUTINE cornr

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE cornr2uv(icr, iuc, ivc, u, v)

      ! Purpose:
      !   to return position of "icr" corner in terms of on-fault coordinates "u" "v". Arrays "iuc" and "ivc" must be computed by
      !   subroutine "cornr"
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   08/03/21                  original version
      !

      INTEGER(i32),               INTENT(IN)  :: icr
      INTEGER(i32), DIMENSION(3), INTENT(IN)  :: iuc, ivc
      REAL(r32),                  INTENT(OUT) :: u, v
      INTEGER(i32)                            :: iu, iv
      REAL(r32)                               :: du

      !-----------------------------------------------------------------------------------------------------------------------------

      iu = iuc(icr)
      iv = ivc(icr)

      du = 0._r32

      IF (MOD(iv, 2) .eq. 0) du = dutr / 2._r32

      u = umingr + (iu - 1) * dutr + du
      v = vmingr + (iv - 1) * dvtr

    END SUBROUTINE cornr2uv

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *


    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *


    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *


END MODULE m_source
