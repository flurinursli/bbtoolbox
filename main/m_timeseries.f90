MODULE m_timeseries

  USE, NON_INTRINSIC :: m_precisions
  USE, NON_INTRINSIC :: m_logfile
  USE, NON_INTRINSIC :: m_strings
  USE, NON_INTRINSIC :: m_toolbox, ONLY: input
  USE, NON_INTRINSIC :: m_userdef, ONLY: parse_usr => read_usr
#ifdef MPI
  USE :: mpi
#endif

  IMPLICIT none

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  PRIVATE

  PUBLIC :: timeseries
  PUBLIC :: read_lp

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  INTEGER(i32), PARAMETER :: SW4_HEADER_LINE = 13

  TYPE :: cmp
    REAL(r32)                              :: dt
    REAL(r32), ALLOCATABLE, DIMENSION(:)   :: time              
    REAL(r32), ALLOCATABLE, DIMENSION(:,:) :: x, y, z
  END TYPE cmp

  TYPE :: tsr
    TYPE(cmp) :: lp, sp, bb
  END TYPE tsr

  TYPE(tsr) :: timeseries            !< timeseries%lp%x(npts, nrec), timeseries%lp%time(npts), timeseries%lp%dt

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  PROCEDURE(parse_sw4), POINTER :: psubr
  PROCEDURE(integrate), POINTER :: pfun

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

    SUBROUTINE disk2seis(ok)

      ! Purpose:
      !   to load long-period timeseries stored to disk into memory. On exit, "ok" is not zero if an error occurred.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   08/03/21                  original version
      !

      INTEGER(i32),                             INTENT(OUT) :: ok
      CHARACTER(:), ALLOCATABLE                             :: fo, extension
      INTEGER(i32)                                          :: lu, rcvr
      INTEGER(i32),              DIMENSION(0:2)             :: request
      REAL(r32)                                             :: dt
      REAL(r32),    ALLOCATABLE, DIMENSION(:)               :: time, x, y, z

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

      ! set processing function
      IF (input%input%variable .ne. input%output%variable) THEN

        IF (TRIM(input%input%variable) .eq. 'displacement') pfun => differentiate

        IF (TRIM(input%input%variable) .eq. 'velocity') THEN
          IF (TRIM(input%output%variable) .eq. 'displacement') THEN
            pfun => integrate
          ELSE
            pfun => differentiate
          ENDIF
        ENDIF

        IF (TRIM(input%input%variable) .eq. 'acceleration') pfun => integrate

      ENDIF

      DO rcvr = 1, SIZE(input%receiver)

        fo = TRIM(input%input%folder) + '/' + TRIM(input%receiver(rcvr)%file) + extension

        OPEN(newunit = lu, file = fo, status = 'old', form = 'formatted', access = 'sequential', action = 'read', iostat = ok)

        IF (ok .ne. 0) THEN
          CALL report_error('disk2seis - ERROR: opening file ' + TRIM(fo) + ' failed with code ' + num2char(ok))
          RETURN
        ENDIF

        CALL psubr(ok, lu, time, x, y, z)

        IF (ok .ne. 0) THEN
          CALL report_error('disk2seis - ERROR: could not read file ' + fo)
          RETURN
        ENDIF

        CLOSE(lu, IOSTAT = ok)

        IF (ok .ne. 0) THEN
          CALL report_error('disk2seis - ERROR: closing file ' + fo + ' returned error code ' + num2char(ok))
          RETURN
        ENDIF

        dt = time(2) - time(1)

        IF (ASSOCIATED(pfun)) THEN
          x = pfun(x, dt)
          y = pfun(y, dt)
          z = pfun(z, dt)

          IF (TRIM(input%input%variable) .ne. 'velocity') THEN
            IF ( (TRIM(input%output%variable) .eq. 'acceleration') .and. (TRIM(input%output%variable) .eq. 'displacement') ) THEN
              x = pfun(x, dt)
              y = pfun(y, dt)
              z = pfun(z, dt)
            ENDIF
          ENDIF

        ENDIF

        ASSOCIATE(n => SIZE(input%receiver))

          ! we make sure that a) all timeseries have same number of points
          !                   b) all timeseries have same time vector

          IF (rcvr .eq. 1) THEN
            ALLOCATE(timeseries%lp%x(SIZE(x), n), timeseries%lp%y(SIZE(x), n), timeseries%lp%z(SIZE(x), n))
          ELSE
            IF (SIZE(timeseries%lp%x, 1) .ne. SIZE(x)) THEN
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

          timeseries%lp%x(:, rcvr) = x
          timeseries%lp%y(:, rcvr) = y
          timeseries%lp%z(:, rcvr) = z

          timeseries%lp%time = time

        END ASSOCIATE

        DEALLOCATE(time, x, y, z)

      ENDDO

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
      INTEGER(i32), DIMENSION(2) :: n

      !-----------------------------------------------------------------------------------------------------------------------------

#ifdef MPI

      n = SHAPE(timeseries%lp%x)

      CALL mpi_bcast(n, 2, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

      IF (.not.ALLOCATED(timeseries%lp%x)) THEN
        ALLOCATE(timeseries%lp%x(n(1),n(2)), timeseries%lp%y(n(1),n(2)), timeseries%lp%z(n(1),n(2)), timeseries%lp%time(n(1)))
      ENDIF

      CALL mpi_bcast(timeseries%lp%x, n(1)*n(2), MPI_REAL, 0, MPI_COMM_WORLD, ierr)
      CALL mpi_bcast(timeseries%lp%y, n(1)*n(2), MPI_REAL, 0, MPI_COMM_WORLD, ierr)
      CALL mpi_bcast(timeseries%lp%z, n(1)*n(2), MPI_REAL, 0, MPI_COMM_WORLD, ierr)

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

      INTEGER(i32),                              INTENT(OUT) :: ok
      INTEGER(i32),                              INTENT(IN)  :: lu
      REAL(r32),    ALLOCATABLE, DIMENSION(:),   INTENT(OUT) :: time, x, y, z
      INTEGER(i32)                                           :: ios, n, lino

      !-----------------------------------------------------------------------------------------------------------------------------

      n   = 0

      DO
        READ(lu, IOSTAT = ok)
        IF (ok .ne. 0) THEN
          IF (IS_IOSTAT_END(ok)) THEN
            EXIT
          ELSE
            RETURN
          ENDIF
        ENDIF
        n = n + 1
      ENDDO

      REWIND(1, IOSTAT = ok)

      IF (ok .ne. 0) RETURN

      n = n - SW4_HEADER_LINE

      IF (n .le. 0) THEN
        CALL report_error('parse_sw4 - ERROR: file is empty or contains no data')
        ok = 1
        RETURN
      ENDIF

      ALLOCATE(time(n), x(n), y(n), z(n))

      DO lino = 1, SW4_HEADER_LINE
        READ(lu, IOSTAT = ok)
        IF (ok .ne. 0) RETURN
      ENDDO

      DO lino = 1, n
        READ(lu, IOSTAT = ok) time(lino), x(lino), y(lino), z(lino)
        IF (ok .ne. 0) RETURN
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

      INTEGER(i32),                              INTENT(OUT) :: ok
      INTEGER(i32),                              INTENT(IN)  :: lu
      REAL(r32),    ALLOCATABLE, DIMENSION(:),   INTENT(OUT) :: time, x, y, z
      INTEGER(i32)                                           :: ios, n, lino

      !-----------------------------------------------------------------------------------------------------------------------------

      n   = 0

      DO
        READ(lu, IOSTAT = ok)
        IF (ok .ne. 0) THEN
          IF (IS_IOSTAT_END(ok)) THEN
            EXIT
          ELSE
            RETURN
          ENDIF
        ENDIF
        n = n + 1
      ENDDO

      REWIND(1, IOSTAT = ok)

      IF (ok .ne. 0) RETURN

      IF (n .le. 0) THEN
        CALL report_error('parse_plaintxt - ERROR: file is empty')
        ok = 1
        RETURN
      ENDIF

      ALLOCATE(time(n), x(n), y(n), z(n))

      DO lino = 1, n
        READ(lu, IOSTAT = ok) time(lino), x(lino), y(lino), z(lino)
        IF (ok .ne. 0) RETURN
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
      !   to compute the cumulative integration of timeseries "x" having time-step "dt" based on the trapezoidal rule.
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

      step = dt * 0.5_r32

      DO i = 2, n
        y(i) = y(i - 1) + step * (x(i) - x(i - 1))
      ENDDO

    END FUNCTION integrate

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *



    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

END MODULE m_timeseries
