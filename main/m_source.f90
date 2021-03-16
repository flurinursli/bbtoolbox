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

  PUBLIC :: hypocenter, plane, setup_source

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  REAL(r32), PARAMETER :: PTSRC_FACTOR = 1._r32 / 10._r32
  REAL(r32), PARAMETER :: PI = 3.14159265358979323846_r64
  REAL(r32), PARAMETER :: DEG_TO_RAD = PI / 180._r32

  TYPE :: hyp
    INTEGER(i32) :: plane        !< plane number where hypocenter is located
    REAL(r32)    :: u, v         !< on-fault coordinates
  END TYPE hyp

  TYPE :: src
    REAL(r32)                                :: length, width, strike, dip         !< segment-specific
    REAL(r32)                                :: targetm0                           !< total moment (sum over all segments)
    REAL(r32)                                :: x, y, z                            !< absolute position of u=v=0
    REAL(r32),   ALLOCATABLE, DIMENSION(:)   :: u, v                               !< on-fault coordinates
    REAL(r32),   ALLOCATABLE, DIMENSION(:,:) :: rise, rupture, sslip, dslip
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

        ALLOCATE(plane(1))

        ! first condition for point-sources: lambda_min >> l, where l is the linear extension (diagonal) of a square fault
        ! du = PTSRC_FACTOR * beta / input%coda%fmax / SQRT(2._r32) / 2._r32

        ! second condition for point-source: max(rupture time) << tr, where tr is the shortes period observed (i.e. 1/fmax)
        ! dt = PTSRC_FACTOR / input%coda%fmax

        ! define all quantities on a small fault whose length=width=2*du
        ! plane(1)%u = [-du, 0._r32, du]
        ! plane(1)%v = [-du, 0._r32, du]

        plane(1)%targetm0 = input%source%m0
        plane(1)%strike   = input%source%strike
        plane(1)%dip      = input%source%dip

        ALLOCATE(plane(1)%u(3), plane(1)%v(3))
        ALLOCATE(plane(1)%rise(3,3), plane(1)%sslip(3,3), plane(1)%dslip(3,3), plane(1)%rupture(3,3))

        ! define rupture times assuming rupture initiate in the middle of the small fault
        ! plane(1)%rupture(:, 1) = [dt, dt/SQRT(2._r32), dt]                   !< proceed along-strike (along u)
        ! plane(1)%rupture(:, 2) = [dt/SQRT(2._r32), 0._r32, dt/SQRT(2._r32)]
        ! plane(1)%rupture(:, 3) = [dt, dt/SQRT(2._r32), dt]

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

        ! actual fault size is connected to velocity model, which may be receiver-dependent... therefore these will be set later on

        ! make sure uppermost down-dip "edge" point is always slightly (1m) below free-surface
        ! plane(1)%z = MAX(1._r32, plane(1)%z + MAX(0._r32, SIN(plane(1)%dip * DEG_TO_RAD) * du - plane(1)%z))

      ELSE
        IF (rank .eq. 0) CALL read_fsp_file(ok)
#ifdef MPI
        CALL mpi_bcast(ok, 1, mpi_int, 0, mpi_comm_world, ierr)
        IF (ok .ne. 0) CALL mpi_abort(mpi_comm_world, ok, ierr)
        CALL broadcast()
#endif
      ENDIF

      IF (rank .eq. ntasks - 1) CALL echo_source()



    END SUBROUTINE setup_source

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE echo_source()

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

        IF (is_empty(plane(1)%targetm0)) THEN
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

    END SUBROUTINE echo_source

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
      INTEGER(i32)                                         :: lu, i, j, k, n, nu, nv, ntw, nsbfs, ios, index, islip, irake, irupt
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
        ALLOCATE(plane(k)%sslip(nu, nv), plane(k)%dslip(nu, nv))

        DO j = 1, nv
          DO i = 1, nu
            plane(k)%sslip(i, j) = 0._r32
            plane(k)%dslip(i, j) = 0._r32
          ENDDO
        ENDDO

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

              DO j = 2, nv - 1
                DO i = 2, nu - 1

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
                  plane(k)%dslip(i, j) = -numeric(islip) * SIN(rake)          !< in our system slip is positive down-dip

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
        plane(i)%rupture(1, 1) = -999._r32

        ALLOCATE(plane(k)%u(nu), plane(k)%v(nv))

        ! set on-fault u-coordinate: first and last point are edge points, second point correspond to reference point u=0
        plane(k)%u(1) = -du / 2._r32
        DO i = 2, nu - 1
          plane(k)%u(i) = (i - 2) * du
        ENDDO
        plane(k)%u(nu) = plane(k)%u(nu - 1) + du / 2._r32

        ! first and last point are edge points, second point correspond to reference point v=0
        plane(k)%v(1) = -dv / 2._r32
        DO i = 2, nv - 1
          plane(k)%v(i) = (i - 2) * dv
        ENDDO
        plane(k)%v(nv) = plane(k)%v(nv - 1) + dv / 2._r32

        ! make sure uppermost down-dip "edge" point is always slightly (1m) below free-surface
        plane(k)%z = MAX(1._r32, plane(k)%z + MAX(0._r32, SIN(plane(k)%dip * DEG_TO_RAD) * dv / 2._r32 - plane(k)%z))

      ENDDO

      CLOSE(lu, iostat = ok)

      IF (ok .ne. 0) THEN
        CALL report_error('Error while closing file ' + TRIM(input%source%file))
        RETURN
      ENDIF

    END SUBROUTINE read_fsp_file

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



END MODULE m_source
