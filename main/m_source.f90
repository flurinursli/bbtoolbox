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
    REAL(r32)                                :: x, y, z                            !< position of u=v=0
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
      INTEGER(i32)              :: ierr

      !-----------------------------------------------------------------------------------------------------------------------------

      IF (input%source%is_point) THEN



      ELSE
        IF (rank .eq. 0) CALL read_fsp_file(ok)
#ifdef MPI
        CALL mpi_bcast(ok, 1, mpi_int, 0, mpi_comm_world, ierr)
        IF (ok .ne. 0) CALL mpi_abort(mpi_comm_world, ok, ierr)
        CALL broadcast()
#endif

        !IF (rank .eq. ntasks - 1) CALL echo_source()

      ENDIF

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

      INTEGER(i32) :: i, j

      !---------------------------------------------------------------------------------------------------------------------------------

      IF (input%source%is_point) THEN

      ELSE




      CALL update_log(num2char('Input folder LF seismograms', width=30, fill='.') +  &
                      num2char(TRIM(input%input%folder), width=79, justify = 'c') + '|')
      CALL update_log(num2char('Format LF seismograms', width=30, fill='.') +  &
                      num2char(TRIM(input%input%format), width=79, justify='c') + '|', blankline = .false.)
      CALL update_log(num2char('Variable LF seismograms', width=30, fill='.') +  &
                      num2char(TRIM(input%input%variable), width=79, justify='c') + '|', blankline = .false.)

      CALL update_log(num2char('Output folder HF seismograms', width=30, fill='.') +  &
                      num2char(TRIM(input%output%folder), width=79, justify = 'c') + '|')
      CALL update_log(num2char('Variable HF seismograms', width=30, fill='.') +  &
                      num2char(TRIM(input%output%variable), width=79, justify='c') + '|', blankline = .false.)

      CALL update_log(num2char('Coda parameters', width=30, fill='.') + num2char('Fmax', width=15, justify='r') + '|' +  &
                      num2char('Matching fr.', width=15, justify='r') + '|' + num2char('Bandwidth', width=15, justify='r') + '|' + &
                      num2char('Seed', width=15, justify='r') + '|' + num2char('Samples', width=15, justify='r') + '|')
      CALL update_log(num2char('', width=30) + &
                      num2char(input%coda%fmax, notation='f', precision=3, width=15, justify='r')      + '|' +  &
                      num2char(input%coda%matching, notation='f', precision=3, width=15, justify='r')  + '|' +  &
                      num2char(input%coda%bandwidth, notation='f', precision=3, width=15, justify='r') + '|' +  &
                      num2char(input%coda%seed, width=15, justify='r')                                 + '|' +  &
                      num2char(input%coda%samples, width=15, justify='r')                              + '|', blankline = .false.)

      IF (input%coda%add_coherency) THEN
        CALL update_log(num2char('Coherency parameters', width=30, fill='.') +     &
                        num2char('Model', width=15, justify='r')       + '|' +     &
                        num2char('Alpha', width=15, justify='r')       + '|' +     &
                        num2char('Threshold', width=15, justify='r')   + '|')
        CALL update_log(num2char('', width=30) +  &
                        num2char(input%coda%model, width=15, justify='r') + '|' + &
                        num2char(input%coda%alpha, notation='s', precision=3, width=15, justify='r') + '|' + &
                        num2char(input%coda%threshold, notation='f', precision=3, width=15, justify='r') + '|', blankline = .false.)
      ENDIF

      IF (input%source%is_point) THEN
        IF (input%origin%is_geo) THEN
          CALL update_log(num2char('Source parameters', width=30, fill='.') +     &
                          num2char('lon', width=15, justify='r') + '|' +  &
                          num2char('lat', width=15, justify='r') + '|' +  &
                          num2char('z', width=15, justify='r') + '|' +  &
                          num2char('Strike', width=15, justify='r') + '|' +  &
                          num2char('Dip', width=15, justify='r')    + '|')
          CALL update_log(num2char('', width=30) + &
                          num2char(input%source%lon, notation='f', width=15, precision=3, justify='r') + '|' + &
                          num2char(input%source%lat, notation='f', width=15, precision=3, justify='r') + '|' + &
                          num2char(input%source%z, notation='f', width=15, precision=3, justify='r') + '|' + &
                          num2char(input%source%strike, notation='f', width=15, precision=1, justify='r') + '|' + &
                          num2char(input%source%dip, notation='f', width=15, precision=1, justify='r')    + '|', blankline=.false.)
        ELSE
          CALL update_log(num2char('Source parameters', width=30, fill='.') +     &
                          num2char('x', width=15, justify='r') + '|' +  &
                          num2char('y', width=15, justify='r') + '|' +  &
                          num2char('z', width=15, justify='r') + '|' +  &
                          num2char('Strike', width=15, justify='r') + '|' +  &
                          num2char('Dip', width=15, justify='r')    + '|')
          CALL update_log(num2char('', width=30) + &
                          num2char(input%source%x, notation='f', width=15, precision=3, justify='r') + '|' + &
                          num2char(input%source%y, notation='f', width=15, precision=3, justify='r') + '|' + &
                          num2char(input%source%z, notation='f', width=15, precision=3, justify='r') + '|' + &
                          num2char(input%source%strike, notation='f', width=15, precision=1, justify='r') + '|' + &
                          num2char(input%source%dip, notation='f', width=15, precision=1, justify='r')    + '|', blankline=.false.)
        ENDIF
        CALL update_log(num2char('(continued)', width=30, fill='.', justify='c') +  &
                        num2char('Rake', width=15, justify='r')   + '|' + &
                        num2char('Moment', width=15, justify='r')     + '|' +  &
                        num2char('Type', width=15, justify='r') + '|' +  &
                        num2char('Frequency', width=15, justify='r') + '|', blankline = .false.)
        CALL update_log(num2char('', width=30) +  &
                        num2char(input%source%rake, notation='f', width=15, precision=1, justify='r')   + '|' + &
                        num2char(input%source%m0, notation='s', width=15, precision=3, justify='r')     + '|' + &
                        num2char(input%source%type, width=15, justify='r') + '|' + &
                        num2char(input%source%freq, notation='f', width=15, precision=3, justify='r') + '|', blankline = .false.)
      ELSE
        CALL update_log(num2char('Rupture file', width=30, fill='.')   +   &
                        num2char(TRIM(input%source%file), width=79, justify='c') + '|')
        CALL update_log(num2char('Rupture parameters', width=30, fill='.')  + &
                        num2char('Type', width=15, justify='r') + '|' + &
                        num2char('Frequency', width=15, justify='r') + '|' + &
                        num2char('Roughness', width=15, justify='r') + '|', blankline = .false.)
        IF (input%source%add_roughness) THEN
          CALL update_log(num2char('', width=30) +  &
                          num2char(input%source%type, width=15, justify='r') + '|' + &
                          num2char(input%source%freq, width=15, notation='f', precision=3, justify='r') + '|' + &
                          num2char(input%source%roughness, width=15, notation='f', precision=1, justify='r') + '|',  &
                          blankline = .false.)
        ELSE
          CALL update_log(num2char('', width=30) +  &
                          num2char(input%source%type, width=15, justify='r') + '|' + &
                          num2char(input%source%freq, width=15, notation='f', precision=3, justify='r') + '|' + &
                          num2char('None', width=15, justify='r') + '|', blankline = .false.)
        ENDIF
        IF (input%source%add_rik) THEN
          CALL update_log(num2char('RIK parameters', width=30, fill='.') +  &
                          num2char('Correlation', width=15, justify='r') + '|' + &
                          num2char('L0', width=15, justify='r')      + '|' + &
                          num2char('Aparam', width=15, justify='r')  + '|' + &
                          num2char('Vrfact', width=15, justify='r')  + '|' + &
                          num2char('Seed', width=15, justify='r')    + '|' + &
                          num2char('Samples', width=15, justify='r') + '|')
          CALL update_log(num2char('', width=30) +  &
                          num2char(input%source%correlation, width=15, notation='f', precision=2, justify='r') + '|' + &
                          num2char(input%source%l0, width=15, notation='f', precision=1, justify='r')          + '|' + &
                          num2char(input%source%aparam, width=15, notation='f', precision=1, justify='r')      + '|' + &
                          num2char(input%source%vrfact, width=15, notation='f', precision=2, justify='r')      + '|' + &
                          num2char(input%source%seed, width=15, justify='r')                                   + '|' + &
                          num2char(input%source%samples, width=15, justify='r') + '|', blankline = .false.)
        ENDIF

      ENDIF

      DO i = 1, SIZE(input%velocity)
        CALL update_log(num2char('Velocity Model ' + num2char(i), width=30, fill='.') +   &
                        num2char('Vp', width=15, justify='r') + '|' + &
                        num2char('Vs', width=15, justify='r') + '|' + &
                        num2char('Rho', width=15, justify='r') + '|' + &
                        num2char('Depth', width=15, justify='r') + '|')
        DO j = 1, SIZE(input%velocity(i)%vp)
          CALL update_log(num2char('', width=30)  +  &
                          num2char(input%velocity(i)%vp(j), width=15, notation='f', precision=3, justify='r') + '|' + &
                          num2char(input%velocity(i)%vs(j), width=15, notation='f', precision=3, justify='r') + '|' + &
                          num2char(input%velocity(i)%rho(j), width=15, notation='f', precision=3, justify='r') + '|' + &
                          num2char(input%velocity(i)%depth(j), width=15, notation='f', precision=3, justify='r') + '|',  &
                          blankline = .false.)
        ENDDO
      ENDDO

      DO i = 1, SIZE(input%attenuation)
        CALL update_log(num2char('Attenuation Model ' + num2char(i), width=30, fill='.') +   &
                        num2char('Low Freq', width=15, justify='r') + '|' + &
                        num2char('High Freq', width=15, justify='r') + '|' + &
                        num2char('Gpp', width=15, justify='r') + '|' + &
                        num2char('Gps', width=15, justify='r') + '|' + &
                        num2char('Gss', width=15, justify='r') + '|' + &
                        num2char('b', width=15, justify='r') + '|')
        DO j = 1, SIZE(input%attenuation(i)%gpp)
          CALL update_log(num2char('', width=30)  +  &
                          num2char(input%attenuation(i)%lcut(j), width=15, notation='f', precision=1, justify='r') + '|' + &
                          num2char(input%attenuation(i)%hcut(j), width=15, notation='f', precision=1, justify='r') + '|' + &
                          num2char(input%attenuation(i)%gpp(j), width=15, notation='f', precision=3, justify='r') + '|' + &
                          num2char(input%attenuation(i)%gps(j), width=15, notation='f', precision=3, justify='r') + '|' + &
                          num2char(input%attenuation(i)%gss(j), width=15, notation='f', precision=3, justify='r') + '|' + &
                          num2char(input%attenuation(i)%b(j), width=15, notation='f', precision=3, justify='r') + '|',  &
                          blankline = .false.)
        ENDDO
      ENDDO

      IF (input%origin%is_geo) THEN
        CALL update_log(num2char('Receivers List', width=30, fill='.') + &
                        num2char('lon', width=15, justify='r') + '|' +  &
                        num2char('lat', width=15, justify='r') + '|' +  &
                        num2char('z', width=15, justify='r') + '|' +  &
                        num2char('File', width=15, justify='r') + '|' + &
                        num2char('Vel Model', width=15, justify='r') + '|' + &
                        num2char('Att Model', width=15, justify='r') + '|')
      ELSE
        CALL update_log(num2char('Receivers List', width=30, fill='.') + &
                        num2char('x', width=15, justify='r') + '|' +  &
                        num2char('y', width=15, justify='r') + '|' +  &
                        num2char('z', width=15, justify='r') + '|' +  &
                        num2char('File', width=15, justify='r') + '|' + &
                        num2char('Vel Model', width=15, justify='r') + '|' + &
                        num2char('Att Model', width=15, justify='r') + '|')
      ENDIF

      IF (input%origin%is_geo) THEN
        DO i = 1, SIZE(input%receiver)
          CALL update_log(num2char('', width=30)  +   &
                          num2char(input%receiver(i)%lon, width=15, notation='f', precision=3, justify='r') + '|' + &
                          num2char(input%receiver(i)%lat, width=15, notation='f', precision=3, justify='r') + '|' + &
                          num2char(input%receiver(i)%z, width=15, notation='f', precision=3, justify='r') + '|' + &
                          num2char(TRIM(input%receiver(i)%file), width=15, justify='r') + '|' + &
                          num2char(input%receiver(i)%velocity, width=15, justify='r')   + '|' + &
                          num2char(input%receiver(i)%attenuation, width=15, justify='r') + '|', blankline = .false.)
        ENDDO
      ELSE
        DO i = 1, SIZE(input%receiver)
          CALL update_log(num2char('', width=30)  +   &
                          num2char(input%receiver(i)%x, width=15, notation='f', precision=3, justify='r') + '|' + &
                          num2char(input%receiver(i)%y, width=15, notation='f', precision=3, justify='r') + '|' + &
                          num2char(input%receiver(i)%z, width=15, notation='f', precision=3, justify='r') + '|' + &
                          num2char(TRIM(input%receiver(i)%file), width=15, justify='r') + '|' + &
                          num2char(input%receiver(i)%velocity, width=15, justify='r')   + '|' + &
                          num2char(input%receiver(i)%attenuation, width=15, justify='r') + '|', blankline = .false.)
        ENDDO
      ENDIF

      CALL update_log(num2char('Advanced settings', width=30, fill='.') +   &
                      num2char('Pmw', width=15, justify = 'r') + '|' +   &
                      num2char('Avg cuts', width=15, justify = 'r')           + '|' +   &
                      num2char('Vr factor', width=15, justify = 'r')          + '|' +   &
                      num2char('Sheets', width=15, justify = 'r')             + '|' +   &
                      num2char('Waves', width=15, justify = 'r')              + '|' +   &
                      num2char('Verbose', width=15, justify = 'r')            + '|')

      CALL update_log(num2char('', width=30) + &
                      num2char(input%advanced%pmw, width=15, justify='r')     + '|' + &
                      num2char(input%advanced%avecuts, width=15, justify='r') + '|' + &
                      num2char(input%advanced%vrfact, width=15, notation='f', precision=2, justify='r') + '|' +  &
                      num2char(input%advanced%sheets, width=15, justify='r')  + '|' + &
                      num2char(input%advanced%waves, width=15, justify='r') + '|' + &
                      num2char(input%advanced%verbose, width=15, justify='r') + '|', blankline=.false.)

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

        ALLOCATE(plane(k)%u(nu), plane(k)%v(nv))

        ! set on-fault u-coordinate: first and last point are edge points, second point correspond to reference point u=0
        plane(k)%u(1) = -du / 2._r32
        DO i = 2, nu - 1
          plane(k)%u(i) = (i - 2) * du
        ENDDO
        plane(k)%u(nu) = plane(k)%u(nu - 1) + du / 2._r32

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


END MODULE m_source
