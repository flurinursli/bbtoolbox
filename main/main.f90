PROGRAM main

  ! Purpose:
  !   to invert coda waves and return best fitting model parameters based on RTT. Sampling of parameters space relies on the NA by
  !   Sambridge. The following inversion setups are possible:
  !
  !     a) inversion of each single coda recorded by a receiver
  !     b) single inversion of all coda recorded by a receiver
  !     c) single inversion of coda generated by the same event and recorded by a group of receivers
  !
  !   The code read recordings in miniseed ASCII format, filter and downsample them, extract envelope and then run the inversion. It
  !   is assumed that recordings first time sample corresponds to earthquake origin time, the latter being set to zero (P/S-wave
  !   arrival times must be based on this convention).
  !
  ! Revisions:
  !     Date                    Description of change
  !     ====                    =====================
  !   18/12/20                  original version
  !
! mpirun -np 1 -x OMP_NUM_THREADS=8 --bind-to socket ./debug.exe deby.par


#ifdef MPI
  USE, NON_INTRINSIC :: mpi
#endif
  USE, NON_INTRINSIC :: m_precisions
  USE, NON_INTRINSIC :: m_toolbox
  USE, NON_INTRINSIC :: m_source
  USE, NON_INTRINSIC :: m_isochron
  USE, NON_INTRINSIC :: m_logfile
  USE, NON_INTRINSIC :: m_strings
  USE, NON_INTRINSIC :: m_noise
  USE, NON_INTRINSIC :: m_rik
  USE, NON_INTRINSIC :: m_timeseries

  IMPLICIT none

  INTEGER(i32)              :: ierr, rank, ntasks, ok, npts, nrecs, i, i0, i1, iter, pl, vel, band, rec
  REAL(r64),   DIMENSION(2) :: tictoc

  !---------------------------------------------------------------------------------------------------------------------------------

  rank   = 0
  ntasks = 1

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --
  ! ----------------------------------------------------- initialization -----------------------------------------------------------

#ifdef MPI
  CALL mpi_init(ierr)
  CALL mpi_comm_rank(mpi_comm_world, rank, ierr)
  CALL mpi_comm_size(mpi_comm_world, ntasks, ierr)
#endif

  CALL watch_start(tictoc(1))

  ! IF (world_rank .eq. 0) CALL set_log_module(ok, screen = .true.)
  CALL set_log_module(ok, screen = .true., errclr = 'red')      !< errors are printed in red

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --
  ! --------------------------------- read main input file and stop execution if error is raised -----------------------------------

  CALL read_input_file(ok, rank, ntasks)

#ifdef MPI
  IF (ok .ne. 0) CALL mpi_abort(mpi_comm_world, ok, ierr)
#else
  STOP
#endif

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --
  ! ------------------------------------------------- setup the source model  ------------------------------------------------------

  CALL setup_source(ok, rank, ntasks)

#ifdef MPI
  IF (ok .ne. 0) CALL mpi_abort(mpi_comm_world, ok, ierr)
#else
  STOP
#endif

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --
  ! ------------------------------------------------- load input time-series  ------------------------------------------------------

  ! CALL read_lp(ok, rank, ntasks)

#ifdef MPI
  IF (ok .ne. 0) CALL mpi_abort(mpi_comm_world, ok, ierr)
#else
  STOP
#endif

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --
  ! ------------------------------------------------------- time stuff -------------------------------------------------------------

  timeseries%sp%dt = 0.5_r32 / 16 !(input%coda%fmax * 2)
  timeseries%cd%dt = 0.5_r32 / 16 !(input%coda%fmax * 2) 

  npts = SIZE(timeseries%lp%time)

  !npts = NINT(timeseries%lp%time(npts) / timeseries%sp%dt) + 1
  npts = NINT(40._r32 / timeseries%sp%dt) + 1

  nrecs = SIZE(input%receiver)

  ALLOCATE(timeseries%sp%xyz(npts, 3, nrecs))
  ALLOCATE(timeseries%cd%xyz(npts, 3, nrecs))

  ALLOCATE(timeseries%sp%time(npts), timeseries%cd%time(npts))

  DO i = 1, npts
    timeseries%sp%time(i) = (i - 1) * timeseries%sp%dt
    timeseries%cd%time(i) = (i - 1) * timeseries%cd%dt
  ENDDO

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --
  ! ---------------------------------------------------- start simulations ---------------------------------------------------------

  ! distribute simulations amongst MPI processes
  CALL split_task(input%coda%samples, ntasks, rank, i0, i1)

  DO iter = i0, i1            !< every iteration is characterized by different source properties (rik, roughness)

    ! reset timeseries
    timeseries%sp%xyz(:,:,:) = 0._r32

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    ! --------------------------------------------- generate coda with random phase ------------------------------------------------

    CALL generate_noise(ok, iter, rank)          !< coda are in the frequency range [0, fmax] in timeseries%cd

    DO vel = 1, SIZE(input%velocity)

      nrecs = COUNT(input%receiver(:)%velocity .eq. vel)

      IF (nrecs .eq. 0) CYCLE       !< skip velocity model if no receivers are associated to it

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! --------------------------------------------- sum contribution over planes -------------------------------------------------

      DO pl = 1, SIZE(plane)

        ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --
        ! --------------------------------------- sum contribution over frequency bands --------------------------------------------

        DO band = 1, SIZE(input%attenuation(1)%gss)

          IF (input%attenuation(1)%lcut(band) .ge. input%coda%fmax) CYCLE         !< skip band if higher than fmax

          CALL update_log(num2char('Simulation for', width=30, fill='.') +         &
                          num2char('Iteration #', width=15, justify='r')   + '|' +   &
                          num2char('Segment #', width=15, justify='r')  + '|' +   &
                          num2char('Velocity #', width=15, justify='r') + '|' +   &
                          num2char('Band #', width=15, justify='r')     + '|')

          CALL update_log(num2char('', width=30) + num2char(iter, width=15, justify='r') + '|' +  &
                                                   num2char(pl, width=15, justify='r')   + '|' +  &
                                                   num2char(vel, width=15, justify='r')  + '|' +  &
                                                   num2char(band, width=15, justify='r') + '|', blankline=.false.)

          ! create mesh for current plane, velocity model and frequency band
          CALL meshing(pl, vel, band)

          ! add rik rupture model
          IF (input%source%add_rik) CALL rik(ok, pl, vel, iter)

#ifdef DEBUG
          CALL node2disk(ok, pl, vel, iter)
#endif

          DO rec = 1, SIZE(input%receiver)
            CALL solve_isochron_integral(ok, rec, band, pl, vel, iter)
          ENDDO

          ! maxband = band        !< highest active frequency band

        ENDDO     !< end loop over frequency bands

        ! CALL correct4impz(ok, pl, maxband)         !< correct amplitude due to filter bank

      ENDDO    !< end loop over fault planes

    ENDDO   !< end loop over velocity models

    CALL correct4impz(ok, iter)

    !CALL add_coda(irec)
    !CALL stitch(irec)

    CALL seis2disk(ok, iter, 'sp')

  ENDDO  !< end loop over iterations

  CALL watch_stop(tictoc(1))

  tictoc = tictoc / 60._r32

  IF (rank .eq. 0) CALL update_log('Program completed in' + num2char(tictoc(1), notation='f', width=10, precision=1) + ' minutes')

  CALL mpi_finalize(ierr)

END PROGRAM main

! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
!===================================================================================================================================
! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
