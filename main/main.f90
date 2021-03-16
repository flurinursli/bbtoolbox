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
  USE, NON_INTRINSIC :: m_logfile
  USE, NON_INTRINSIC :: m_strings

  IMPLICIT none

  INTEGER(i32)              :: ierr, rank, ntasks, ok, i0, i1, iter, io, ip
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
  CALL set_log_module(ok, screen = .true.)

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --
  ! --------------------------------- read main input file and stop execution if error is raised -----------------------------------

  IF (rank .eq. 0) CALL read_input_file(ok)

#ifdef MPI
  CALL mpi_bcast(ok, 1, mpi_int, 0, mpi_comm_world, ierr)

  IF (ok .ne. 0) CALL mpi_abort(mpi_comm_world, ok, ierr)
#endif

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --
  ! ------------------------------------------- broadcast relevant input parameters  -----------------------------------------------

#ifdef MPI
  CALL broadcast()
#endif

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --
  ! ------------------------------------------------- echo input parameters  -------------------------------------------------------

  IF (rank .eq. ntasks - 1) CALL echo_input()

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --
  ! ------------------------------------------------- setup the source model  ------------------------------------------------------

  CALL setup_source(ok, rank, ntasks)

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --
  ! ------------------------------------------------------- time stuff -------------------------------------------------------------

  dt = 0.25_r32 / input%coda%fmax             !< fnyq is twice fmax


  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --
  ! --------------------------------------------- generate coda with random phase --------------------------------------------------

  CALL generate_noise(ok)

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --
  ! ---------------------------------------------------- start simulations ---------------------------------------------------------

  CALL split_task(input%source%samples, ntasks, rank, i0, i1)

  DO iter = i0, i1            !< every iteration is characterized by different source properties

    DO ip = 1, SIZE(plane)

      !CALL setup_mesh(ip)
      !CALL roughness(ip, iter)

      DO io = 1, SIZE(input%receiver)
        !CALL quake(io, ip)      !< compute direct wave contributions and stack envelopes (one component only)
      ENDDO

    ENDDO

    DO io = 1, SIZE(input%receiver)
      !CALL add_coda(io)
      !CALL stitch(io)
      !CALL save2disk(io, iter)
    ENDDO

  ENDDO

  CALL watch_stop(tictoc(1))

  tictoc = tictoc / 60._r32

  IF (rank .eq. 0) CALL update_log('Program completed in' + num2char(tictoc(1), notation='f', width=10, precision=1) + ' minutes')

END PROGRAM main

! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
!===================================================================================================================================
! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
