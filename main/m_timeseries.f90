MODULE m_timeseries

  USE, NON_INTRINSIC :: m_precisions
  USE, NON_INTRINSIC :: m_logfile
  USE, NON_INTRINSIC :: m_fft_real
  USE, NON_INTRINSIC :: m_interpolation_r32
  USE, NON_INTRINSIC :: m_strings
  USE, NON_INTRINSIC :: m_toolbox, ONLY: input
  USE, NON_INTRINSIC :: m_userdef, ONLY: parse_usr => read_usr, write_usr
#ifdef MPI
  USE :: mpi
#endif

  IMPLICIT none

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  PRIVATE

  PUBLIC :: timeseries, amplification
  PUBLIC :: read_lp, seis2disk, differentiate, stitch, load_amplification, amplify

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  INTEGER(i32), PARAMETER :: SW4_HEADER_LINE = 13
  REAL(r32),    PARAMETER :: PI = 3.14159265358979323846_r64

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  CHARACTER(2) :: m_cp

  INTEGER(i32) :: m_iter, m_rcvr

  TYPE :: cmp
    REAL(r32)                                :: dt
    REAL(r32), ALLOCATABLE, DIMENSION(:)     :: time
    REAL(r32), ALLOCATABLE, DIMENSION(:,:,:) :: xyz
  END TYPE cmp

  TYPE :: tsr
    TYPE(cmp) :: lp, sp, bb, cd
  END TYPE tsr

  TYPE(tsr) :: timeseries            !< timeseries%lp%x(npts, nrec), timeseries%lp%time(npts), timeseries%lp%dt

  TYPE :: amp
    REAL(r32), ALLOCATABLE, DIMENSION(:) :: frequency, value
  END TYPE amp

  TYPE(amp), ALLOCATABLE, DIMENSION(:) :: amplification

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  PROCEDURE(parse_sw4), POINTER :: psubr => NULL()
  PROCEDURE(write_sw4), POINTER :: wsubr => NULL()
  PROCEDURE(integrate), POINTER :: pfun => NULL()

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  CONTAINS

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE read_lp(ok, rank, ntasks)

      ! Purpose:
      !   to read long-period timeseries. On exit, "ok" is not zero if an error occurred.
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

      IF (rank .eq. 0) CALL disk2seis(ok)

#ifdef MPI
      CALL mpi_bcast(ok, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
#endif

      IF (ok .ne. 0) RETURN

#ifdef MPI
      CALL broadcast()
#endif

      timeseries%lp%dt = timeseries%lp%time(2) - timeseries%lp%time(1)

      ! IF (rank .eq. ntasks - 1) CALL echo()

    END SUBROUTINE read_lp

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE load_amplification(ok, rank, ntasks)

      ! Purpose:
      !   to read (optional) amplification curves stored in plain ASCII files. On exit, "ok" is not zero if an error occurred.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   01/07/21                  original version
      !

      INTEGER(i32), INTENT(OUT) :: ok
      INTEGER(i32), INTENT(IN)  :: rank, ntasks
      INTEGER(i32)              :: ierr, rcvr, n
      LOGICAL                   :: is_allocated

      !-----------------------------------------------------------------------------------------------------------------------------

      IF (rank .eq. 0) CALL parse_amplification(ok)

#ifdef MPI
      CALL mpi_bcast(ok, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
#endif

      IF (ok .ne. 0) RETURN

#ifdef MPI

      is_allocated = ALLOCATED(amplification)

      CALL mpi_bcast(is_allocated, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

      IF (.not.is_allocated) RETURN       !< do nothing if amplification curves are not available

      DO rcvr = 1, SIZE(input%receiver)

        IF (rank .eq. 0) n = SIZE(amplification(rcvr)%frequency)

        CALL mpi_bcast(n, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

        IF (rank .ne. 0) ALLOCATE(amplification(rcvr)%frequency(n), amplification(rcvr)%value(n))

        CALL mpi_bcast(amplification(rcvr)%frequency, n, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
        CALL mpi_bcast(amplification(rcvr)%value, n, MPI_REAL, 0, MPI_COMM_WORLD, ierr)

      ENDDO

#endif

    END SUBROUTINE load_amplification

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE parse_amplification(ok)

      ! Purpose:
      !   to read (optional) amplification curves stored in plain ASCII files. On exit, "ok" is not zero if an error occurred.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   01/07/21                  original version
      !

      INTEGER(i32),             INTENT(OUT) :: ok
      CHARACTER(256)                        :: msg
      CHARACTER(:), ALLOCATABLE             :: fo
      INTEGER(i32)                          :: lu, rcvr, n, lino, ierr, idum
      REAL(r32)                             :: rdum

      !-----------------------------------------------------------------------------------------------------------------------------

      ok = 0

      fo = TRIM(input%input%amplification)

      IF (fo .eq. 'none') RETURN

      ALLOCATE(amplification(SIZE(input%receiver)))

      DO rcvr = 1, SIZE(input%receiver)

        fo = TRIM(input%input%amplification) + '/' + TRIM(input%receiver(rcvr)%file) + '.txt'

        OPEN(NEWUNIT = lu, FILE = fo, STATUS = 'old', FORM = 'formatted', ACCESS = 'sequential', action = 'read', IOSTAT = ok,   &
             IOMSG = msg)

        IF (ok .ne. 0) THEN
          CALL report_error('parse_amplification - ERROR: ' + TRIM(msg))
          RETURN
        ENDIF

        n = 0

        DO
          READ(lu, *, IOSTAT = ok, IOMSG = msg)
          IF (ok .ne. 0) THEN
            IF (IS_IOSTAT_END(ok)) THEN
              EXIT
            ELSE
              CALL report_error('parse_amplification - ERROR: ' + TRIM(msg))
              RETURN
            ENDIF
          ENDIF
          n = n + 1
        ENDDO

        REWIND(lu, IOSTAT = ok, IOMSG = msg)

        IF (ok .ne. 0) THEN
          CALL report_error('parse_amplification - ERROR: ' + TRIM(msg))
          RETURN
        ENDIF

        IF (n .le. 0) THEN
          CALL report_error('parse_amplification - ERROR: file ' + TRIM(fo) + ' is empty')
          ok = 1
          RETURN
        ENDIF

        ALLOCATE(amplification(rcvr)%frequency(n), amplification(rcvr)%value(n))

        DO lino = 1, n
          READ(lu, *, IOSTAT = ok, IOMSG = msg) amplification(rcvr)%frequency(lino), amplification(rcvr)%value(lino),       &
                                                rdum, rdum, idum
          IF (ok .ne. 0) THEN
            CALL report_error('parse_amplification - ERROR: ' + TRIM(msg))
            RETURN
          ENDIF
        ENDDO

        CLOSE(lu, IOSTAT = ok, IOMSG = msg)

        IF (ok .ne. 0) THEN
          CALL report_error('parse_amplification - ERROR: ' + TRIM(msg))
          RETURN
        ENDIF


amplification(rcvr)%value(:) = 2.


      ENDDO

    END SUBROUTINE parse_amplification

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE disk2seis(ok)

      ! Purpose:
      !   to load long-period timeseries stored to disk into memory. On exit, "ok" is not zero if an error occurred.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   08/03/21                  original version
      !

      INTEGER(i32),                            INTENT(OUT) :: ok
      CHARACTER(256)                                       :: msg
      CHARACTER(:),  ALLOCATABLE                           :: fo, extension
      INTEGER(i32)                                         :: lu, rcvr
      REAL(r32)                                            :: dt
      REAL(r32),     ALLOCATABLE, DIMENSION(:)             :: time, x, y, z

      !-----------------------------------------------------------------------------------------------------------------------------

      fo = TRIM(input%input%format)

      IF ( (fo .eq. 'sw4') .or. (fo .eq. 'SW4') ) THEN
        extension = '.txt'
        psubr => parse_sw4

      ELSEIF ( (fo .eq. 'txt') .or. (fo .eq. 'TXT') ) THEN
        extension = '.txt'
        psubr => parse_plaintxt

      ELSEIF ( (fo .eq. 'usr') .or. (fo .eq. 'USR') ) THEN
        extension = '.txt'
        psubr => parse_usr

      ELSE
        CALL report_error('disk2seis - ERROR: ' + fo + ' is an invalid format descriptor for long-period timeseries')
        ok = 1

      ENDIF

      ! set processing function such that input timeseries are velocities
      IF (TRIM(input%input%variable) .eq. 'displacement') pfun => differentiate
      IF (TRIM(input%input%variable) .eq. 'acceleration') pfun => integrate


      ! IF (input%input%variable .ne. input%output%variable) THEN
      !
      !   IF (TRIM(input%input%variable) .eq. 'displacement') pfun => differentiate
      !
      !   IF (TRIM(input%input%variable) .eq. 'velocity') THEN
      !     IF (TRIM(input%output%variable) .eq. 'displacement') THEN
      !       pfun => integrate
      !     ELSE
      !       pfun => differentiate
      !     ENDIF
      !   ENDIF
      !
      !   IF (TRIM(input%input%variable) .eq. 'acceleration') pfun => integrate
      !
      ! ENDIF

      DO rcvr = 1, SIZE(input%receiver)

        fo = TRIM(input%input%folder) + '/' + TRIM(input%receiver(rcvr)%file) + extension

        OPEN(NEWUNIT = lu, FILE = fo, STATUS = 'old', FORM = 'formatted', ACCESS = 'sequential', ACTION = 'read', IOSTAT = ok,  &
             IOMSG = msg)

        IF (ok .ne. 0) THEN
          CALL report_error('disk2seis - ERROR: ' + TRIM(msg))
          RETURN
        ENDIF

        CALL psubr(ok, lu, time, x, y, z)

        IF (ok .ne. 0) THEN
          CALL report_error('disk2seis - ERROR: could not read file ' + fo)
          RETURN
        ENDIF

        CLOSE(lu, IOSTAT = ok, IOMSG = msg)

        IF (ok .ne. 0) THEN
          CALL report_error('disk2seis - ERROR: ' + TRIM(msg))
          RETURN
        ENDIF

        dt = time(2) - time(1)

        IF (ASSOCIATED(pfun)) THEN
          x = pfun(x, dt)
          y = pfun(y, dt)
          z = pfun(z, dt)

          ! IF (TRIM(input%input%variable) .ne. 'velocity') THEN
          !   IF ( (TRIM(input%output%variable) .eq. 'acceleration') .and. (TRIM(input%output%variable) .eq. 'displacement') ) THEN
          !     x = pfun(x, dt)
          !     y = pfun(y, dt)
          !     z = pfun(z, dt)
          !   ENDIF
          ! ENDIF

        ENDIF

        ASSOCIATE(n => SIZE(input%receiver))

          ! we make sure that a) all timeseries have same number of points
          !                   b) all timeseries have same time vector

          IF (rcvr .eq. 1) THEN
            ALLOCATE(timeseries%lp%xyz(SIZE(x), 3, n))
          ELSE
            IF (SIZE(timeseries%lp%xyz, 1) .ne. SIZE(x)) THEN
              CALL report_error('disk2seis - ERROR: all long-period timeseries must have same number of points')
              ok = 1
              RETURN
            ENDIF
            IF (ANY(timeseries%lp%time .ne. time)) THEN
              CALL report_error('disk2seis - ERROR: all long-period timeseries must have identical time vector')
              ok = 1
              RETURN
            ENDIF
          ENDIF

          timeseries%lp%xyz(:, 1, rcvr) = x
          timeseries%lp%xyz(:, 2, rcvr) = y
          timeseries%lp%xyz(:, 3, rcvr) = z

          timeseries%lp%time = time

        END ASSOCIATE

        DEALLOCATE(time, x, y, z)

      ENDDO

      pfun => NULL()

    END SUBROUTINE disk2seis

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE broadcast()

      ! Purpose:
      !   to broadcast all timeseries to all processes.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   08/03/21                  original version
      !

      INTEGER(i32)               :: rcvr, ierr
      INTEGER(i32), DIMENSION(3) :: n

      !-----------------------------------------------------------------------------------------------------------------------------

#ifdef MPI

      n = SHAPE(timeseries%lp%xyz)

      CALL mpi_bcast(n, 3, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

      IF (.not.ALLOCATED(timeseries%lp%xyz)) THEN
        ALLOCATE(timeseries%lp%xyz(n(1),n(2),n(3)), timeseries%lp%time(n(1)))
      ENDIF

      CALL mpi_bcast(timeseries%lp%xyz, n(1)*n(2)*n(3), MPI_REAL, 0, MPI_COMM_WORLD, ierr)
      CALL mpi_bcast(timeseries%lp%time, n(1), MPI_REAL, 0, MPI_COMM_WORLD, ierr)

#endif

    END SUBROUTINE broadcast

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE parse_sw4(ok, lu, time, x, y, z)

      ! Purpose:
      !   to read ASCII sw4/wpp output files and return associated time-series. On exit, "ok" is not zero if an error occurred.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   08/03/21                  original version
      !

      INTEGER(i32),                             INTENT(OUT) :: ok
      INTEGER(i32),                             INTENT(IN)  :: lu
      REAL(r32),     ALLOCATABLE, DIMENSION(:), INTENT(OUT) :: time, x, y, z
      CHARACTER(256)                                        :: msg
      INTEGER(i32)                                          :: ios, n, lino

      !-----------------------------------------------------------------------------------------------------------------------------

      n = 0

      DO
        READ(lu, *, IOSTAT = ok, IOMSG = msg)
        IF (ok .ne. 0) THEN
          IF (IS_IOSTAT_END(ok)) THEN
            EXIT
          ELSE
            CALL report_error('parse_sw4 - ERROR: ' + TRIM(msg))
            RETURN
          ENDIF
        ENDIF
        n = n + 1
      ENDDO

      REWIND(1, IOSTAT = ok, IOMSG = msg)

      IF (ok .ne. 0) THEN
        CALL report_error('parse_sw4 - ERROR: ' + TRIM(msg))
        RETURN
      ENDIF

      n = n - SW4_HEADER_LINE

      IF (n .le. 0) THEN
        CALL report_error('parse_sw4 - ERROR: file is empty or contains no data')
        ok = 1
        RETURN
      ENDIF

      ALLOCATE(time(n), x(n), y(n), z(n))

      DO lino = 1, SW4_HEADER_LINE
        READ(lu, *, IOSTAT = ok, IOMSG = msg)
        IF (ok .ne. 0) THEN
          CALL report_error('parse_sw4 - ERROR: ' + TRIM(msg))
          RETURN
        ENDIF
      ENDDO

      DO lino = 1, n
        READ(lu, *, IOSTAT = ok, IOMSG = msg) time(lino), x(lino), y(lino), z(lino)
        IF (ok .ne. 0) THEN
          CALL report_error('parse_sw4 - ERROR: ' + TRIM(msg))
          RETURN
        ENDIF
      ENDDO

    END SUBROUTINE parse_sw4

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE parse_plaintxt(ok, lu, time, x, y, z)

      ! Purpose:
      !   to read a plain ASCII file and return associated time-series. The file is supposed to have four columns and an arbitrary
      !   number of rows, where the first column is the time samples and all other columns are the x, y and z component of motion,
      !   respectively. On exit, "ok" is not zero if an error occurred.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   08/03/21                  original version
      !

      INTEGER(i32),                               INTENT(OUT) :: ok
      INTEGER(i32),                               INTENT(IN)  :: lu
      REAL(r32),     ALLOCATABLE, DIMENSION(:),   INTENT(OUT) :: time, x, y, z
      CHARACTER(256)                                          :: msg
      INTEGER(i32)                                            :: ios, n, lino

      !-----------------------------------------------------------------------------------------------------------------------------

      n = 0

      DO
        READ(lu, IOSTAT = ok, IOMSG = msg)
        IF (ok .ne. 0) THEN
          IF (IS_IOSTAT_END(ok)) THEN
            EXIT
          ELSE
            CALL report_error('parse_plaintxt - ERROR: ' + TRIM(msg))
            RETURN
          ENDIF
        ENDIF
        n = n + 1
      ENDDO

      REWIND(1, IOSTAT = ok, IOMSG = msg)

      IF (ok .ne. 0) THEN
        CALL report_error('parse_plaintxt - ERROR: ' + TRIM(msg))
        RETURN
      ENDIF

      IF (n .le. 0) THEN
        CALL report_error('parse_plaintxt - ERROR: file is empty')
        ok = 1
        RETURN
      ENDIF

      ALLOCATE(time(n), x(n), y(n), z(n))

      DO lino = 1, n
        READ(lu, IOSTAT = ok, IOMSG = msg) time(lino), x(lino), y(lino), z(lino)
        IF (ok .ne. 0) THEN
          CALL report_error('parse_plaintxt - ERROR: ' + TRIM(msg))
          RETURN
        ENDIF
      ENDDO

    END SUBROUTINE parse_plaintxt

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    FUNCTION differentiate(x, dt) RESULT(y)

      ! Purpose:
      !   to differentiate timeseries "x" having time-step "dt".
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   08/03/21                  original version
      !

      REAL(r32),   DIMENSION(:),      INTENT(IN) :: x
      REAL(r32),                      INTENT(IN) :: dt
      REAL(r32),   DIMENSION(SIZE(x))            :: y
      INTEGER(i32)                               :: n, i

      !-----------------------------------------------------------------------------------------------------------------------------

      n = SIZE(x)

      y(1) = (x(2) - x(1)) / dt           !< forward difference

      DO i = 2, n - 1
        y(i) = 0.5_r32 * (x(i + 1) - x(i - 1)) / dt       !< central difference
      ENDDO

      y(n) = (x(n) - x(n - 1)) / dt     !< backward difference

    END FUNCTION differentiate

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    FUNCTION integrate(x, dt) RESULT(y)

      ! Purpose:
      !   to compute the cumulative integral of timeseries "x" having time-step "dt".
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   08/03/21                  original version
      !

      REAL(r32),   DIMENSION(:),      INTENT(IN) :: x
      REAL(r32),                      INTENT(IN) :: dt
      INTEGER(i32)                               :: n, i
      REAL(r32)                                  :: step
      REAL(r32),   DIMENSION(SIZE(x))            :: y

      !-----------------------------------------------------------------------------------------------------------------------------

      n = SIZE(x)

      y(1) = 0._r32

      DO i = 2, n
        y(i) = y(i - 1) + x(i) * dt
      ENDDO

    END FUNCTION integrate

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE seis2disk(ok, iter, seis)

      ! Purpose:
      !   to store hybrid timeseries to disk. On exit, "ok" is not zero if an error occurred.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   08/03/21                  original version
      !

      INTEGER(i32),                           INTENT(OUT) :: ok
      INTEGER(i32),                           INTENT(IN)  :: iter
      CHARACTER(2),                           INTENT(IN)  :: seis
      CHARACTER(:), ALLOCATABLE                           :: fo, extension
      INTEGER(i32)                                        :: lu, rcvr, ic
      REAL(r32)                                           :: dt
      REAL(r32),    ALLOCATABLE, DIMENSION(:)             :: time, x, y, z

      !-----------------------------------------------------------------------------------------------------------------------------

      fo = TRIM(input%output%format)

      IF ( (fo .eq. 'sw4') .or. (fo .eq. 'SW4') ) THEN
        extension = '.txt'
        wsubr => write_sw4

      ELSEIF ( (fo .eq. 'txt') .or. (fo .eq. 'TXT') ) THEN
        extension = '.txt'
        wsubr => write_plaintxt

      ELSEIF ( (fo .eq. 'paz') .or. (fo .eq. 'PAZ') ) THEN
        extension = '.txt'
        wsubr => write_bafu

      ELSEIF ( (fo .eq. 'usr') .or. (fo .eq. 'USR') ) THEN
        extension = '.txt'
        wsubr => write_usr

      ELSE
        CALL report_error('seis2disk - ERROR: ' + fo + ' is an invalid format descriptor for output files')
        ok = 1

      ENDIF

      IF ( (seis .ne. 'lp') .and. (seis .ne. 'sp') .and. (seis .ne. 'bb') .and. (seis .ne. 'cd') ) THEN
        CALL report_error('seis2disk - ERROR: ' + seis + ' is an invalid descriptor for timeseries')
        ok = 1
      ENDIF

      ! set processing function such that output timeseries are transformed from velocity as necessary
      IF (TRIM(input%output%variable) .eq. 'displacement') pfun => integrate
      IF (TRIM(input%output%variable) .eq. 'acceleration') pfun => differentiate

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! ------------------------------------------ BAFU format (one component per file) --------------------------------------------

      IF ( (fo .eq. 'paz') .or. (fo .eq. 'PAZ') ) THEN

        DO rcvr = 1, SIZE(input%receiver)

          DO ic = 1, 3

            IF (ic .eq. 1) THEN
              m_cp = 'NS'
            ELSEIF (ic .eq. 2) THEN
              m_cp = 'EW'
            ELSE
              m_cp = 'UD'
            ENDIF

            fo = TRIM(input%output%folder) + '/' + uppercase(input%output%variable(1:3)) + '_' + num2char(iter) + '_' +    &
                 TRIM(input%receiver(rcvr)%file) + '.' + m_cp + '.' + seis + extension

            OPEN(NEWUNIT = lu, FILE = fo, STATUS = 'replace', FORM = 'formatted', access= 'sequential', action= 'write', IOSTAT= ok)

            IF (ok .ne. 0) THEN
              CALL report_error('seis2disk - ERROR: opening file ' + TRIM(fo) + ' failed with code ' + num2char(ok))
              RETURN
            ENDIF

            IF (seis .eq. 'lp') THEN

              ASSOCIATE(x => timeseries%lp%xyz(:, ic, rcvr), time => timeseries%lp%time, dt => timeseries%lp%dt)

                IF (ASSOCIATED(pfun)) x = pfun(x, dt) * 100._r32       !< move from m to cm

                m_iter = iter         !< pass values to module aliases
                m_rcvr = rcvr

                CALL wsubr(ok, lu, time, x, x, x)

              END ASSOCIATE

            ELSEIF (seis .eq. 'sp') THEN

              ASSOCIATE(x => timeseries%sp%xyz(:, ic, rcvr), time => timeseries%sp%time, dt => timeseries%sp%dt)

                IF (ASSOCIATED(pfun)) x = pfun(x, dt) * 100._r32        !< move from m to cm

                m_iter = iter         !< pass values to module aliases
                m_rcvr = rcvr

                CALL wsubr(ok, lu, time, x, x, x)

              END ASSOCIATE

            ELSEIF (seis .eq. 'bb') THEN

              ASSOCIATE(x => timeseries%bb%xyz(:, ic, rcvr), time => timeseries%bb%time, dt => timeseries%bb%dt)

                IF (ASSOCIATED(pfun)) x = pfun(x, dt) * 100._r32        !< move from m to cm

                m_iter = iter         !< pass values to module aliases
                m_rcvr = rcvr

                CALL wsubr(ok, lu, time, x, x, x)

              END ASSOCIATE

            ELSEIF (seis .eq. 'cd') THEN

              ASSOCIATE(x => timeseries%cd%xyz(:, ic, rcvr), time => timeseries%cd%time, dt => timeseries%cd%dt)

                IF (ASSOCIATED(pfun)) x = pfun(x, dt) * 100._r32        !< move from m to cm

                m_iter = iter         !< pass values to module aliases
                m_rcvr = rcvr

                CALL wsubr(ok, lu, time, x, x, x)

              END ASSOCIATE

            ENDIF

            IF (ok .ne. 0) THEN
              CALL report_error('seis2disk - ERROR: could not write to file ' + fo)
              RETURN
            ENDIF

            CLOSE(lu, IOSTAT = ok)

            IF (ok .ne. 0) THEN
              CALL report_error('seis2disk - ERROR: closing file ' + fo + ' returned error code ' + num2char(ok))
              RETURN
            ENDIF

          ENDDO      !< end loop over components

        ENDDO    !< end loop over receivers

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! -------------------------------------- all other formats (three component per file) ----------------------------------------

      ELSE

        DO rcvr = 1, SIZE(input%receiver)

          fo = TRIM(input%output%folder) + '/' + TRIM(input%receiver(rcvr)%file) + '_' + seis + '_' + num2char(iter) + extension

          OPEN(NEWUNIT = lu, FILE = fo, STATUS = 'replace', FORM = 'formatted', ACCESS = 'sequential', action = 'write', IOSTAT= ok)

          IF (ok .ne. 0) THEN
            CALL report_error('seis2disk - ERROR: opening file ' + TRIM(fo) + ' failed with code ' + num2char(ok))
            RETURN
          ENDIF

          IF (seis .eq. 'lp') THEN

            ASSOCIATE(x => timeseries%lp%xyz(:, 1, rcvr), y => timeseries%lp%xyz(:, 2, rcvr), z => timeseries%lp%xyz(:, 3, rcvr),  &
              time => timeseries%lp%time, dt => timeseries%lp%dt)

              IF (ASSOCIATED(pfun)) THEN
                x = pfun(x, dt)
                y = pfun(y, dt)
                z = pfun(z, dt)
              ENDIF

              m_iter = iter         !< pass values to module aliases
              m_rcvr = rcvr

              CALL wsubr(ok, lu, time, x, y, z)

            END ASSOCIATE

          ELSEIF (seis .eq. 'sp') THEN

            ASSOCIATE(x => timeseries%sp%xyz(:, 1, rcvr), y => timeseries%sp%xyz(:, 2, rcvr), z => timeseries%sp%xyz(:, 3, rcvr),  &
              time => timeseries%sp%time, dt => timeseries%sp%dt)

              IF (ASSOCIATED(pfun)) THEN
                x = pfun(x, dt)
                y = pfun(y, dt)
                z = pfun(z, dt)
              ENDIF

              m_iter = iter         !< pass values to module aliases
              m_rcvr = rcvr

              CALL wsubr(ok, lu, time, x, y, z)

            END ASSOCIATE

          ELSEIF (seis .eq. 'bb') THEN

            ASSOCIATE(x => timeseries%bb%xyz(:, 1, rcvr), y => timeseries%bb%xyz(:, 2, rcvr), z => timeseries%bb%xyz(:, 3, rcvr),  &
                      time => timeseries%bb%time, dt => timeseries%bb%dt)

              IF (ASSOCIATED(pfun)) THEN
                x = pfun(x, dt)
                y = pfun(y, dt)
                z = pfun(z, dt)
              ENDIF

              m_iter = iter         !< pass values to module aliases
              m_rcvr = rcvr

              CALL wsubr(ok, lu, time, x, y, z)

            END ASSOCIATE

          ELSEIF (seis .eq. 'cd') THEN

            ASSOCIATE(x => timeseries%cd%xyz(:, 1, rcvr), y => timeseries%cd%xyz(:, 2, rcvr), z => timeseries%cd%xyz(:, 3, rcvr),  &
              time => timeseries%cd%time, dt => timeseries%cd%dt)

              IF (ASSOCIATED(pfun)) THEN
                x = pfun(x, dt)
                y = pfun(y, dt)
                z = pfun(z, dt)
              ENDIF

              m_iter = iter         !< pass values to module aliases
              m_rcvr = rcvr

              CALL wsubr(ok, lu, time, x, y, z)

            END ASSOCIATE

          ENDIF

          IF (ok .ne. 0) THEN
            CALL report_error('seis2disk - ERROR: could not write to file ' + fo)
            RETURN
          ENDIF

          CLOSE(lu, IOSTAT = ok)

          IF (ok .ne. 0) THEN
            CALL report_error('seis2disk - ERROR: closing file ' + fo + ' returned error code ' + num2char(ok))
            RETURN
          ENDIF

        ENDDO

      ENDIF

      pfun => NULL()

    END SUBROUTINE seis2disk

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE write_sw4(ok, lu, time, x, y, z)

      ! Purpose:
      !   to write ASCII sw4-style output files where timeseries are stored. On exit, "ok" is not zero if an error occurred.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   08/03/21                  original version
      !

      INTEGER(i32),                           INTENT(OUT) :: ok
      INTEGER(i32),                           INTENT(IN)  :: lu
      REAL(r32),                DIMENSION(:), INTENT(IN)  :: time, x, y, z
      CHARACTER(5)                                        :: zone
      CHARACTER(8)                                        :: date
      CHARACTER(10)                                       :: hms
      CHARACTER(:), ALLOCATABLE                           :: timestamp, lonlat, xy, msg
      INTEGER(i32)                                        :: lino

      !-----------------------------------------------------------------------------------------------------------------------------

      CALL DATE_AND_TIME(date, hms, zone)

      timestamp = date(7:8) + '/' + date(5:6) + '/' + date(1:4) + ' ' + hms(1:2) + ':' + hms(3:4) + ':' + hms(5:10) + ' (UTC ' +  &
                  zone(1:3) + ':' + zone(4:5) + ')'
      lonlat    = num2char(input%receiver(m_rcvr)%lon, notation='f', width=9,  precision=3) +   &
                  num2char(input%receiver(m_rcvr)%lat, notation='f', width=9,  precision=3)
      xy        = num2char(input%receiver(m_rcvr)%x,   notation='f', width=15, precision=1) +   &
                  num2char(input%receiver(m_rcvr)%y,   notation='f', width=15, precision=1)

      WRITE(lu, '(A)', IOSTAT = ok, IOMSG = msg) '# Author: BBTool'
      IF (ok .ne. 0) THEN
        CALL report_error('write_sw4 - ERROR: ' + msg)
        RETURN
      ENDIF

      WRITE(lu, '(A)', IOSTAT = ok, IOMSG = msg) '# Scenario: ' + num2char(m_iter)
      IF (ok .ne. 0) THEN
        CALL report_error('write_sw4 - ERROR: ' + msg)
        RETURN
      ENDIF

      WRITE(lu, '(A)', IOSTAT = ok, IOMSG = msg) '# Date: ' + timestamp
      IF (ok .ne. 0) THEN
        CALL report_error('write_sw4 - ERROR: ' + msg)
        RETURN
      ENDIF

      WRITE(lu, '(A)', IOSTAT=ok, IOMSG=msg) '# Bandwidth (Hz): ' + num2char(input%coda%fmax, notation='f', width=6, precision=1)
      IF (ok .ne. 0) THEN
        CALL report_error('write_sw4 - ERROR: ' + msg)
        RETURN
      ENDIF

      WRITE(lu, '(A)', IOSTAT = ok, IOMSG = msg) '# Station: ' + TRIM(input%receiver(m_rcvr)%file)
      IF (ok .ne. 0) THEN
        CALL report_error('write_sw4 - ERROR: ' + msg)
        RETURN
      ENDIF

      WRITE(lu, '(A)', IOSTAT = ok, IOMSG = msg) '# Location (WGS84 longitude, latitude) (deg): ' + lonlat
      IF (ok .ne. 0) THEN
        CALL report_error('write_sw4 - ERROR: ' + msg)
        RETURN
      ENDIF

      WRITE(lu, '(A)', IOSTAT = ok, IOMSG = msg) '# Location (northing, easting) (m): ' + xy
      IF (ok .ne. 0) THEN
        CALL report_error('write_sw4 - ERROR: ' + msg)
        RETURN
      ENDIF

      WRITE(lu, '(A)', IOSTAT = ok, IOMSG = msg) '# Ground motion: ' + TRIM(input%output%variable)
      IF (ok .ne. 0) THEN
        CALL report_error('write_sw4 - ERROR: ' + msg)
        RETURN
      ENDIF

      WRITE(lu, '(A)', IOSTAT = ok, IOMSG = msg) '# nColumns: 4'
      IF (ok .ne. 0) THEN
        CALL report_error('write_sw4 - ERROR: ' + msg)
        RETURN
      ENDIF

      WRITE(lu, '(A)', IOSTAT = ok, IOMSG = msg) '# Column 1: Time (s)'
      IF (ok .ne. 0) THEN
        CALL report_error('write_sw4 - ERROR: ' + msg)
        RETURN
      ENDIF

      WRITE(lu, '(A)', IOSTAT = ok, IOMSG = msg) '# Column 2: NS (x)'
      IF (ok .ne. 0) THEN
        CALL report_error('write_sw4 - ERROR: ' + msg)
        RETURN
      ENDIF

      WRITE(lu, '(A)', IOSTAT = ok, IOMSG = msg) '# Column 3: EW (y)'
      IF (ok .ne. 0) THEN
        CALL report_error('write_sw4 - ERROR: ' + msg)
        RETURN
      ENDIF

      WRITE(lu, '(A)', IOSTAT = ok, IOMSG = msg) '# Column 4: UD (z)'
      IF (ok .ne. 0) THEN
        CALL report_error('write_sw4 - ERROR: ' + msg)
        RETURN
      ENDIF

      DO lino = 1, SIZE(time)
#ifdef DOUBLE_PREC
        WRITE(lu, '(4(ES23.13))', IOSTAT = ok, IOMSG = msg) time(lino), x(lino), y(lino), z(lino)
#else
        WRITE(lu, '(4(ES15.5))', IOSTAT = ok, IOMSG = msg) time(lino), x(lino), y(lino), z(lino)
#endif
        IF (ok .ne. 0) THEN
          CALL report_error('write_sw4 - ERROR: ' + msg)
          RETURN
        ENDIF
      ENDDO

    END SUBROUTINE write_sw4

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE write_plaintxt(ok, lu, time, x, y, z)

      ! Purpose:
      !   to write a plain ASCII file where time-series are stored. The file is supposed to have four columns and an arbitrary
      !   number of rows, where the first column is the time samples and all other columns are the x, y and z component of motion,
      !   respectively. On exit, "ok" is not zero if an error occurred.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   08/03/21                  original version
      !

      INTEGER(i32),                           INTENT(OUT) :: ok
      INTEGER(i32),                           INTENT(IN)  :: lu
      REAL(r32),                DIMENSION(:), INTENT(IN)  :: time, x, y, z
      CHARACTER(:), ALLOCATABLE                           :: msg
      INTEGER(i32)                                        :: lino

      !-----------------------------------------------------------------------------------------------------------------------------

      DO lino = 1, SIZE(time)
#ifdef DOUBLE_PREC
        WRITE(lu, '(4(ES23.13))', IOSTAT = ok, IOMSG = msg) time(lino), x(lino), y(lino), z(lino)
#else
        WRITE(lu, '(4(ES15.5))', IOSTAT = ok, IOMSG = msg) time(lino), x(lino), y(lino), z(lino)
#endif
        IF (ok .ne. 0) THEN
          CALL report_error('write_plaintxt - ERROR: ' + msg)
          RETURN
        ENDIF
      ENDDO

    END SUBROUTINE write_plaintxt

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE write_bafu(ok, lu, time, x, y, z)

      ! Purpose:
      !   to write ASCII output files where timeseries are stored based on the format specifications provided by F.Panzera for the
      !   BAFU project. On exit, "ok" is not zero if an error occurred. Note that in this case "x", "y" and "z" all refer to the
      !   same vector.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   02/07/21                  original version
      !

      INTEGER(i32),                           INTENT(OUT) :: ok
      INTEGER(i32),                           INTENT(IN)  :: lu
      REAL(r32),                DIMENSION(:), INTENT(IN)  :: time, x, y, z
      CHARACTER(5)                                        :: zone
      CHARACTER(8)                                        :: date
      CHARACTER(10)                                       :: hms
      CHARACTER(:), ALLOCATABLE                           :: timestamp, msg
      INTEGER(i32)                                        :: lino

      !-----------------------------------------------------------------------------------------------------------------------------

      CALL DATE_AND_TIME(date, hms, zone)

      timestamp = date(7:8) + '/' + date(5:6) + '/' + date(1:4) + ' ' + hms(1:2) + ':' + hms(3:4) + ':' + hms(5:10) + ' (UTC ' +  &
      zone(1:3) + ':' + zone(4:5) + ')'

      timestamp = date(1:4) + date(5:6) + date(7:8) + hms(1:2) + hms(3:4) + hms(5:6)

      WRITE(lu, '(A)', IOSTAT = ok, IOMSG = msg) 'K-NET and KiK-net strong motion database'
      IF (ok .ne. 0) THEN
        CALL report_error('write_bafu - ERROR: ' + msg)
        RETURN
      ENDIF

      WRITE(lu, '(A)', IOSTAT = ok, IOMSG = msg) 'KiK-net_ID: ' + timestamp
      IF (ok .ne. 0) THEN
        CALL report_error('write_bafu - ERROR: ' + msg)
        RETURN
      ENDIF

      WRITE(lu, '(A)', IOSTAT = ok, IOMSG = msg) 'Event_time_HH:MM:SS: ' + hms(1:2) + hms(3:4) + hms(5:6)
      IF (ok .ne. 0) THEN
        CALL report_error('write_bafu - ERROR: ' + msg)
        RETURN
      ENDIF

      WRITE(lu, '(A)', IOSTAT = ok, IOMSG = msg) 'Station_code: ' + TRIM(input%receiver(m_rcvr)%file)
      IF (ok .ne. 0) THEN
        CALL report_error('write_bafu - ERROR: ' + msg)
        RETURN
      ENDIF

      WRITE(lu, '(A)', IOSTAT = ok, IOMSG = msg) 'Stream: ' + m_cp
      IF (ok .ne. 0) THEN
        CALL report_error('write_bafu - ERROR: ' + msg)
        RETURN
      ENDIF

      WRITE(lu, '(A)', IOSTAT = ok, IOMSG = msg) 'NPTS: ' + num2char(SIZE(time))
      IF (ok .ne. 0) THEN
        CALL report_error('write_bafu - ERROR: ' + msg)
        RETURN
      ENDIF

      WRITE(lu, '(A)', IOSTAT = ok, IOMSG = msg) 'Sample_rate_(s): ' + num2char(time(2) - time(1))
      IF (ok .ne. 0) THEN
        CALL report_error('write_bafu - ERROR: ' + msg)
        RETURN
      ENDIF

      msg = 'cm'

      IF (input%output%variable .eq. 'velocity') THEN
        msg = msg + '/s'
      ELSEIF (input%output%variable .eq. 'acceleration') THEN
        msg = msg + '/s^2'
      ENDIF

      WRITE(lu, '(A)', IOSTAT = ok, IOMSG = msg) 'Units: ' + msg
      IF (ok .ne. 0) THEN
        CALL report_error('write_bafu - ERROR: ' + msg)
        RETURN
      ENDIF

      WRITE(lu, '(A)', IOSTAT = ok, IOMSG = msg) 'Filter_type: BUTTERWORTH'
      IF (ok .ne. 0) THEN
        CALL report_error('write_bafu - ERROR: ' + msg)
        RETURN
      ENDIF

      WRITE(lu, '(A)', IOSTAT = ok, IOMSG = msg) 'Filter order: 2'
      IF (ok .ne. 0) THEN
        CALL report_error('write_bafu - ERROR: ' + msg)
        RETURN
      ENDIF

      WRITE(lu, '(A)', IOSTAT = ok, IOMSG = msg) 'Low_cut_frequency_(Hz): 0'
      IF (ok .ne. 0) THEN
        CALL report_error('write_bafu - ERROR: ' + msg)
        RETURN
      ENDIF

      WRITE(lu, '(A)', IOSTAT = ok, IOMSG = msg) 'High_cut_frequency_(Hz): ' + num2char(input%coda%fmax, precision=1, notation='f')
      IF (ok .ne. 0) THEN
        CALL report_error('write_bafu - ERROR: ' + msg)
        RETURN
      ENDIF

      WRITE(lu, '(A)', IOSTAT = ok, IOMSG = msg) 'Trigger: Normal_triggered'
      IF (ok .ne. 0) THEN
        CALL report_error('write_bafu - ERROR: ' + msg)
        RETURN
      ENDIF

      DO lino = 1, SIZE(time)
#ifdef DOUBLE_PREC
        WRITE(lu, '(4(ES23.13))', IOSTAT = ok, IOMSG = msg) x(lino)
#else
        WRITE(lu, '(4(ES15.5))', IOSTAT = ok, IOMSG = msg) x(lino)
#endif
        IF (ok .ne. 0) THEN
          CALL report_error('write_bafu - ERROR: ' + msg)
          RETURN
        ENDIF
      ENDDO

    END SUBROUTINE write_bafu

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE stitch(ok)

      ! Purpose:
      !   to merge long- and short-period timeseries for receiver "recvr" together. On exit, "ok" is not zero if an error occurred.
      !   The subroutine follows the procedure outlined in Mai&Beroza (2003) based on a couple of matched filters (Hann windows)
      !   whose characteristics are controlled by "input%coda%matching" and "input%coda%bandwidth".
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   08/03/21                  original version
      !

      INTEGER(i32),                           INTENT(OUT) :: ok
      CHARACTER(:), ALLOCATABLE                           :: fo
      INTEGER(i32)                                        :: i, ic, rcvr, npts, l, lu, i0, n
      COMPLEX(r32), ALLOCATABLE, DIMENSION(:)             :: lspectrum, sspectrum
      REAL(r32)                                           :: df
      REAL(r32),    ALLOCATABLE, DIMENSION(:)             :: lpres, freq, win

      !-----------------------------------------------------------------------------------------------------------------------------

      CALL setup_interpolation('linear', 'zero', ok)

      IF (ok .ne. 0) RETURN

      ASSOCIATE(lp => timeseries%lp%xyz, sp => timeseries%sp%xyz, ptime => timeseries%lp%time, stime => timeseries%sp%time,   &
                bb => timeseries%bb%xyz)

        npts = SIZE(stime, 1)

        CALL make_fftw_plan([npts])

        ALLOCATE(lpres(npts))
        ALLOCATE(lspectrum(npts/2 + 1), sspectrum(npts/2 + 1), freq(npts/2 + 1), win(npts/2 + 1))

        df = 1._r32 / (timeseries%sp%dt * npts)

        ! frequency vector
        DO i = 1, npts/2 + 1
          freq(i) = (i - 1) * df
        ENDDO

        ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --
        ! ----------------------------------------------- build matched filter -----------------------------------------------------

        n = NINT(2._r32 * input%coda%bandwidth / df) + 1
        l = ((n - 1) / 2) + 1                                 !< size (in points) of rising part of window

        i0 = NINT(input%coda%matching / df) + 1

        DO i = 1, i0 - (l / 2) - 1
          win(i) = 0._r32
        ENDDO

        DO i = i0 - (l / 2), i0 - (l / 2) + (l - 1)
          win(i) = SIN(PI * (i - i0 + (l / 2)) / (n - 1))**2
        ENDDO

        DO i = i0 - (l / 2) + (l - 1) + 1, npts / 2 + 1
          win(i) = 1._r32
        ENDDO

        DO rcvr = 1, SIZE(input%receiver)
          DO ic = 1, 3

            CALL interpolate(ptime, lp(:, ic, rcvr), stime, lpres)            !< resample long-period timeseries

            CALL fft(lpres, lspectrum)
            CALL fft(sp(:, ic, rcvr), sspectrum)

            DO i = 1, npts/2 + 1
              sspectrum(i) = lspectrum(i) * (1._r32 - win(i)) + sspectrum(i) * win(i)      !< spectrum hybrid timeseries
            ENDDO

            CALL ifft(bb(:, ic, rcvr), sspectrum)

          ENDDO
        ENDDO

        CALL destroy_fftw_plan([npts])

      END ASSOCIATE

#ifdef DEBUG

      fo = 'matched.txt'

      OPEN(NEWUNIT = lu, FILE = fo, STATUS = 'replace', FORM = 'formatted', ACCESS = 'sequential', action = 'write', IOSTAT= ok)

      IF (ok .ne. 0) THEN
        CALL report_error('Error while opening file' + TRIM(fo))
        RETURN
      ENDIF

      DO i = 1, SIZE(freq)
        WRITE(lu, *) freq(i), win(i)
      ENDDO

      CLOSE(lu, IOSTAT = ok)

      IF (ok .ne. 0) THEN
        CALL report_error('Error while closing file ' + fo)
        RETURN
      ENDIF

#endif

    END SUBROUTINE stitch

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE amplify(ok)

      ! Purpose:
      !   to apply user-defined amplification curves (if available) to every broadband timeseries.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   08/03/21                  original version
      !

      INTEGER(i32),                           INTENT(OUT) :: ok
      COMPLEX(r32), ALLOCATABLE, DIMENSION(:)             :: spectrum
      INTEGER(i32)                                        :: npts, i, rcvr, ic
      REAL(r32)                                           :: df
      REAL(r32),   ALLOCATABLE, DIMENSION(:)              :: freq, x

      !-----------------------------------------------------------------------------------------------------------------------------

      ok = 0

      IF (TRIM(input%input%amplification) .eq. 'none') RETURN             !< do nothing if amplification curves are not available

      CALL setup_interpolation('linear', 'zero', ok)

      IF (ok .ne. 0) RETURN

      npts = SIZE(timeseries%bb%time)

      ALLOCATE(freq(npts/2 + 1), x(npts/2 + 1), spectrum(npts/2 + 1))

      CALL make_fftw_plan([npts])

      df = 1._r32 / (timeseries%bb%dt * npts)

      ! frequency vector
      DO i = 1, npts/2 + 1
        freq(i) = (i - 1) * df
      ENDDO

      DO rcvr = 1, SIZE(input%receiver)

        CALL interpolate(amplification(rcvr)%frequency, amplification(rcvr)%value, freq, x)

        DO ic = 1, 3

          CALL fft(timeseries%bb%xyz(:, ic, rcvr), spectrum)

          DO i = 1, npts/2 + 1
            spectrum(i) = spectrum(i) * x(i)
          ENDDO

          CALL ifft(timeseries%bb%xyz(:, ic, rcvr), spectrum)

        ENDDO

      ENDDO

      CALL destroy_fftw_plan([npts])

    END SUBROUTINE amplify

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

END MODULE m_timeseries
