MODULE m_isochron

  USE, NON_INTRINSIC :: m_precisions
  USE, NON_INTRINSIC :: m_interpolation_r32
  USE, NON_INTRINSIC :: m_compgeo
  USE, NON_INTRINSIC :: m_filter
  USE, NON_INTRINSIC :: m_llsq_r32
  USE, NON_INTRINSIC :: m_stat
  USE, NON_INTRINSIC :: m_source
  USE, NON_INTRINSIC :: m_rik
  USE, NON_INTRINSIC :: m_rtt
  USE, NON_INTRINSIC :: m_roughness
  USE, NON_INTRINSIC :: m_toolbox
  USE, NON_INTRINSIC :: m_strings
  USE, NON_INTRINSIC :: m_timeseries
  USE, NON_INTRINSIC :: m_logfile
  USE, NON_INTRINSIC :: m_wkbj
#ifdef MPI
  USE :: mpi
#endif

  IMPLICIT none

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  PRIVATE

  PUBLIC :: solve_isochron_integral, node2disk, correct4impz

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

#ifdef MPI
  INTEGER(i32), PARAMETER :: COMM = MPI_COMM_SELF
#else
  INTEGER(i32), PARAMETER :: COMM = 0
#endif

  REAL(r32), PARAMETER :: PI = 3.14159265358979323846_r64
  REAL(r32), PARAMETER :: DEG_TO_RAD = PI / 180._r32
  REAL(r32), PARAMETER :: BIG = HUGE(0._r32)
  REAL(r32), PARAMETER :: DP = 0.002_r32, DPSM = 0.07_r32         !< free-surface smoothing parameters
  REAL(r32), PARAMETER :: DC = 0.5_r32                            !< resolution coda tablau (traveltime)

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  INTEGER(i32),              DIMENSION(2)       :: minsheets, maxsheets
  REAL(r32),    ALLOCATABLE, DIMENSION(:)       :: shooting                 !< depth of shooting points

  TYPE :: co
    INTEGER(i32), ALLOCATABLE, DIMENSION(:,:) :: displs, counts
    REAL(r32),    ALLOCATABLE, DIMENSION(:)   :: direct
    REAL(r32),    ALLOCATABLE, DIMENSION(:,:) :: lotrvt, uptrvt
    REAL(r32),    ALLOCATABLE, DIMENSION(:,:) :: envelope
  END TYPE

  TYPE(co) :: coda

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  PROCEDURE(rik_at_nodes), POINTER :: nodefun

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  CONTAINS

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE solve_isochron_integral(ok, rec, band, pl, vel, iter)

      ! Purpose:
      !   to evaluate the isochron integral for plane "pl" embedded in velocity model "vel" and iteration "iter". Integration occurs
      !   on a plane whose initial geometry and mechanism are randomized based on fractal roughness. Contributions of direct P- and
      !   S-waves are calculated for each triangle making up the composite (refined) mesh.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   08/03/21                  original version
      !

      USE, NON_INTRINSIC :: m_fft_real

      INTEGER(i32),                             INTENT(OUT) :: ok
      INTEGER(i32),                             INTENT(IN)  :: rec, band, pl, vel, iter
      COMPLEX(r32)                                          :: fxpfzs, fzpfxs, g31, g21, g23, ga, gb
      COMPLEX(r32),              DIMENSION(3,3)             :: fp, fs, g
      COMPLEX(r32), ALLOCATABLE, DIMENSION(:)               :: fxp, fzp, tur, tui, tv
      INTEGER(i32)                                          :: n, i, j, ref, icr, totnutr, seed, src, sheet, wtp, it1, it2, ic
      INTEGER(i32)                                          :: it, cut, model, npts, noslip
      INTEGER(i32),              DIMENSION(3)               :: iuc, ivc, abv, blw, shot
      ! INTEGER(i32), ALLOCATABLE, DIMENSION(:,:)             :: ncuts, active
      REAL(r32)                                             :: m0, moment, strike, dip, urec, vrec, wrec, taumin, taumax, sd, cd, dt
      REAL(r32)                                             :: srfvel, c, scale
      REAL(r32)                                             :: avepath, avetrvt, scattering, attenuation, t, a0
      REAL(r32),                 DIMENSION(2)               :: ncuts, active
      REAL(r32),                 DIMENSION(3)               :: u, v, w, x, y, z, slip, rupture, rake, nrl, rho, alpha, beta, mu
      REAL(r32),                 DIMENSION(3)               :: p, q, path, trvt, repi, tau, sr, cr, velloc
      REAL(r32),    ALLOCATABLE, DIMENSION(:)               :: fsp, mrf, envelope, stack, time
      REAL(r32),    ALLOCATABLE, DIMENSION(:,:)             :: rtri, itri, rseis, iseis, seis
      REAL(r64)                                             :: tictoc
      REAL(r64),                 DIMENSION(5)               :: timer
      REAL(r64),                 DIMENSION(10)              :: otimer

      !-----------------------------------------------------------------------------------------------------------------------------

      ok = 0

      ! don't do anything if velocity model is not the right one for current receiver
      IF (input%receiver(rec)%velocity .ne. vel) RETURN         !< consider only receivers for current velocity model "vel"

      ASSOCIATE(model => input%attenuation(1), fmax => input%coda%fmax, dt => timeseries%sp%dt)

        IF (input%advanced%verbose .eq. 2) THEN
          CALL update_log(num2char('Integration for recv', width=30, fill='.') + num2char(rec, width=15, justify='r') + '|')
          CALL update_log(num2char('<freq band>', width=30, justify='c') +   &
                          num2char(num2char(model%lcut(band), notation='f', width=6, precision=1) + ', ' +   &
                                   num2char(MIN(model%hcut(band), fmax), notation='f', width=6, precision=1),           &
                          width=15, justify='r') + '|')
        ENDIF

      END ASSOCIATE

      model = input%receiver(rec)%attenuation      !< attenuation model

#ifdef PERF
      timer(:)  = 0._r64
      otimer(:) = 0._r64
#endif

      CALL setup_interpolation('linear', 'zero', ok)

      IF (ok .ne. 0) RETURN

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! ---------------------------------------- lookup table free-surface coefficients --------------------------------------------

      CALL lookup_fs(ok, rec, fsp, fxp, fzp)

      IF (ok .ne. 0) RETURN

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! ---------------------------------------- allocate arrays to handle seismograms ---------------------------------------------

      dt = timeseries%sp%dt !* (input%coda%fmax / MIN(input%attenuation(1)%hcut(band), fmax))

      IF (input%source%is_point) dt = MAXVAL(plane(pl)%rupture) / (input%advanced%avecuts * 4)

      npts = NINT(timeseries%sp%time(SIZE(timeseries%sp%time)) / dt) + 2

      ALLOCATE(seis(npts, 3), time(npts))

      seis(:, :) = 0._r32

      DO it = 1, npts
        time(it) = (it - 1) * dt
      ENDDO

      IF (input%advanced%verbose .eq. 2) THEN
        CALL update_log(num2char('<time step>', width=30, justify='c') +    &
                        num2char(dt, notation='f', width=15, precision=4, justify='r') + '|')
      ENDIF

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! --------------------------------------- solve integral for each mesh refinement --------------------------------------------

      moment = 0._r32      !< keep track of actual moment, this will be used to rescale synthetics to target moment

      ! select function to compute rupture at grid nodes, depending whether we are dealing with extended- or point-sources
      nodefun => rik_at_nodes

      IF (input%source%is_point) nodefun => ptrsrc_at_nodes

      ! set "seed" such that random numbers depend on fault plane number and iteration
      seed = input%source%seed + (iter - 1) * SIZE(plane) + pl

      ALLOCATE(rseis(npts, 3), iseis(npts, 3), rtri(npts, 3), itri(npts, 3))
      ALLOCATE(mrf(npts), tur(npts/2+1), tui(npts/2+1), tv(npts/2+1))

      ! loop over mesh refinements
      DO ref = 1, SIZE(nvtr)

        IF (input%advanced%verbose .eq. 2) THEN
          CALL update_log(num2char('<mesh refinement>', width=30, justify='c') + num2char(ref, width=15, justify='r') + '|')
        ENDIF

        ! these stack contributions from each mesh refinement
        DO ic = 1, 3
          DO it = 1, npts
            rseis(it, ic) = 0._r32
            iseis(it, ic) = 0._r32
          ENDDO
        ENDDO

        noslip    = 0       !< total triangles having zero slip
        ncuts(:)  = 0       !< total isochron cuts for each wave type
        active(:) = 0       !< total triangles that were cut

        ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --
        ! ------------------------------------------ define RIK model on mesh nodes ------------------------------------------------

#ifdef PERF
        CALL watch_start(tictoc, COMM)
#endif

        ALLOCATE(nodes(nugr(ref), nvgr(ref)))

        CALL nodefun(ref, pl, vel, seed)                !< define slip, rupture time and rise time on mesh nodes

#ifdef PERF
        CALL watch_stop(tictoc, COMM)
        timer(1) = timer(1) + tictoc
#endif

        totnutr = 2 * nutr(ref) - 1                     !< total triangles in a row

        m0 = 0._r32

#ifdef PERF
        CALL watch_start(tictoc, COMM)
#endif

        ! generate fault plane roughness
        CALL fault_roughness(ok, ref, pl, iter)

#ifdef PERF
        CALL watch_stop(tictoc, COMM)
        timer(2) = timer(2) + tictoc
#endif

#ifdef PERF
        CALL watch_start(tictoc, COMM)
#endif

        ! shoot rays between receivers and a set of sources spanning (vertically) current mesh
        CALL rayshooting(ok, ref, pl, vel)

#ifdef PERF
        CALL watch_stop(tictoc, COMM)
        timer(3) = timer(3) + tictoc
#endif

#ifdef PERF
        CALL watch_start(tictoc, COMM)
#endif

        ! fill lookup table for envelopes
        CALL lookup_coda(ok, ref, rec, band, pl, vel, time)

#ifdef PERF
        CALL watch_stop(tictoc, COMM)
        timer(3) = timer(3) + tictoc
#endif


        ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --
        ! ------------------------------------------- allocate memory for FFT/IFFT -------------------------------------------------

        CALL make_fftw_plan([npts])

        !$omp parallel do default(shared) private(i, j, iuc, ivc, icr, slip, u, v, w, x, y, z, nrl, strike, dip, rake, sd, cd)   &
        !$omp private(sr, cr, tictoc, rupture, alpha, beta, rho, mu, ic, src, shot, mrf, urec, vrec, wrec, repi)     &
        !$omp private(rtri, itri, wtp, abv, blw, velloc, srfvel, scattering, sheet, p, path, trvt, q, tau, taumin, taumax, it1)   &
        !$omp private(it2, g, attenuation, cut, c, it, tur, tui, tv, a0)   &
        !$omp reduction(+: m0, rseis, iseis, otimer, ncuts, active, noslip)
        DO j = 1, nvtr(ref)
          DO i = 1, totnutr

            ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
            ! ---------------------------------------- check corners have non-zero slip --------------------------------------------

            CALL cornr(j, i, iuc, ivc)        !< corner indices for current triangle

            DO icr = 1, 3
              slip(icr) = SUM(nodes(iuc(icr), ivc(icr))%slip)
            ENDDO

            noslip = noslip + 1

            IF (ALL(slip .eq. 0._r32)) CYCLE          !< jump to next triangle if current has zero slip

            noslip = noslip - 1

            ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
            ! --------------------------------------------- get perturbed mechanism ------------------------------------------------

#ifdef PERF
            CALL watch_start(tictoc, COMM)
#endif

            CALL cornr2uv(iuc, ivc, ref, u, v)        !< on-fault coordinates

            DO icr = 1, 3
              w(icr) = roughness(iuc(icr), ivc(icr))      !< add roughness (always zero for point-sources)
            ENDDO

            CALL uvw2xyz(pl, u, v, w, x, y, z)            !< get cartesian coordinates

            DO icr = 1, 3
              z(icr) = MAX(MIN_DEPTH, z(icr))             !< make sure we never breach the free-surface (clip roughness)
            ENDDO

            CALL normal2tri(x, y, z, nrl)

            ! always have normal pointing upward (i.e. negative z direction)
            IF (MOD(i + (j-1)*totnutr, 2) == 0) nrl(:) = -nrl(:)

            strike = plane(pl)%strike * DEG_TO_RAD
            dip    = plane(pl)%dip    * DEG_TO_RAD

            CALL interpolate(plane(pl)%u, plane(pl)%v, plane(pl)%rake, u, v, rake)

            CALL perturbed_mechanism(strike, dip, rake, nrl)            !< get new strike, dip, rake (all values in radians)

            sd = SIN(dip)
            cd = COS(dip)

            DO icr = 1, 3
              sr(icr) = SIN(rake(icr))
              cr(icr) = COS(rake(icr))
            ENDDO

#ifdef PERF
            CALL watch_stop(tictoc, COMM)
            otimer(1) = otimer(1) + tictoc
#endif

            ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
            ! ----------------------------------------------- corner parameters ----------------------------------------------------

#ifdef PERF
            CALL watch_start(tictoc, COMM)
#endif

            DO icr = 1, 3
              rupture(icr) = MINVAL(nodes(iuc(icr), ivc(icr))%rupture)
              alpha(icr)   = vinterp(input%velocity(vel)%depth, input%velocity(vel)%vp, input%velocity(vel)%vpgrad, z(icr))
              beta(icr)    = vinterp(input%velocity(vel)%depth, input%velocity(vel)%vs, input%velocity(vel)%vsgrad, z(icr))
              rho(icr)     = vinterp(input%velocity(vel)%depth, input%velocity(vel)%rho, input%velocity(vel)%rhograd, z(icr))
              mu(icr)      = rho(icr) * beta(icr)**2         !< rigidity
            ENDDO

#ifdef PERF
            CALL watch_stop(tictoc, COMM)
            otimer(2) = otimer(2) + tictoc
#endif

            ! move from "x-y" to "u-t" coordinates for corners
            u = x
            v = y
            w = 0._r32
            CALL rotate(u, v, w, 0._r32, 0._r32, -strike)

            ! find shooting points above/below
            DO icr = 1, 3
              DO src = 1, SIZE(shooting) - 1
                IF ( (z(icr) .ge. shooting(src)) .and. (z(icr) .le. shooting(src + 1)) ) shot(icr) = src
              ENDDO
            ENDDO

! print*, 'shot ', shot, ' - ', shooting, ' - ', z

            ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
            ! ---------------------------------------- moment rate function of triangle --------------------------------------------

#ifdef PERF
            CALL watch_start(tictoc, COMM)
#endif

            IF (input%source%is_point) THEN
              mrf = dbrune(time, input%source%freq)
            ELSE
              CALL mrf_rik(iuc, ivc, dt, rupture, mrf)     !< moment rate function (averaged over corners)
            ENDIF

            CALL normalize(mrf, dt)

            CALL fft(mrf, tv)                              !< ... and its spectrum

#ifdef PERF
            CALL watch_stop(tictoc, COMM)
            otimer(3) = otimer(3) + tictoc
#endif

            ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
            ! -------------------------------------- triangle contribution to ground motion ----------------------------------------

            urec = input%receiver(rec)%x
            vrec = input%receiver(rec)%y
            wrec = 0._r32

            CALL rotate(urec, vrec, wrec, 0._r32, 0._r32, -strike)       !< move to "u-t" coordinates for receiver

            DO icr = 1, 3
              repi(icr) = HYPOT(urec - u(icr), vrec - v(icr)) / 1000._r32        !< epicentral distance (corner-receiver)
            ENDDO

            ! initialize arrays with triangle contribution
            DO ic = 1, 3
              DO it = 1, npts
                rtri(it, ic) = 0._r32
                itri(it, ic) = 0._r32
              ENDDO
            ENDDO

            ! loop over wave types
            DO wtp = 1, 2

              ! abv(:) = 1       !< reset sheet-related index
              ! blw(:) = 1

              IF (wtp .eq. 1) THEN
                velloc(:)  = alpha(:) / 1000._r32
                srfvel     = input%velocity(input%receiver(rec)%velocity)%vp(1) / 1000._r32
                scattering = input%attenuation(model)%gpp(band) + input%attenuation(model)%gps(band)             !< P
              ELSE
                velloc(:)  = beta(:) / 1000._r32
                srfvel     = input%velocity(input%receiver(rec)%velocity)%vs(1) / 1000._r32
                scattering = input%attenuation(model)%gss(band) + input%attenuation(model)%gps(band) / 6._r32    !< S
              ENDIF

              DO sheet = 1, maxsheets(wtp)

#ifdef PERF
                CALL watch_start(tictoc, COMM)
#endif

                CALL raypar_at_corners(sheet, shooting, repi, z, shot, wtp, p, path, trvt, q)

#ifdef PERF
                CALL watch_stop(tictoc, COMM)
                otimer(4) = otimer(4) + tictoc
#endif

                ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
                ! ------------------------------------------ time integration limits -----------------------------------------------


! print*, 'r ', rec, wtp, sheet, ' x ', mean(repi), ' - ', path, ' - ', trvt


                IF (ANY(trvt .eq. 0._r32)) CYCLE       !< all corners must belong to same sheet

                tau(:) = trvt(:) + rupture(:)          !< sum travel-time and rupture time

                taumax = MAXVAL(tau)
                taumin = MINVAL(tau)

                it1 = NINT(taumin / dt) + 1

                it2 = it1 + (taumax - taumin) / dt - 1      !< when (taumax-taumin) < dt, it2 < it1

                it2 = MIN(npts, it2)                        !< limit to max number of time points


! print*, 't ', rec, wtp, sheet, ' - ', taumin, taumax, ' - ', it1, it2


                ! jump to next sheet if no isochron is spanned (i.e. no cuts)
                IF (it2 .lt. it1) CYCLE

                active(wtp) = active(wtp) + 1

                ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
                ! ------------------------------------------ evaluate integral kernel ----------------------------------------------

#ifdef PERF
                CALL watch_start(tictoc, COMM)
#endif

                CALL integral_kernel(wtp, urec, vrec, u, v, velloc, srfvel, sd, cd, sr, cr, p, fsp, fxp, fzp, g)

                ! complete integral kernel
                DO icr = 1, 3
                  g(1, icr) = g(1, icr) * slip(icr) * mu(icr) * q(icr)
                  g(2, icr) = g(2, icr) * slip(icr) * mu(icr) * q(icr)
                  g(3, icr) = g(3, icr) * slip(icr) * mu(icr) * q(icr)
                ENDDO

#ifdef PERF
                CALL watch_stop(tictoc, COMM)
                otimer(5) = otimer(5) + tictoc
#endif

                ! term for direct wave amplitude attenuation
                attenuation = EXP(-(mean(path) * scattering + input%attenuation(model)%b(band) * mean(trvt)) / 2._r32)

#ifdef PERF
                CALL watch_start(tictoc, COMM)
#endif

                ! "rtri" and "itri" contain direct wave contribution for current triangle and receiver
                CALL time_integration(ref, it1, it2, dt, attenuation, tau, g, c, cut, rtri, itri)

                ncuts(wtp) = ncuts(wtp) + cut

! print*, 'c ', rec, wtp, sheet, ' - ', cut, ' - ', ncuts(wtp), ' - ', active(wtp)


#ifdef PERF
                CALL watch_stop(tictoc, COMM)
                otimer(6) = otimer(6) + tictoc
#endif

                ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
                ! ------------------------------------------ add coda contribution -----------------------------------------------

                ! (averaged) amplitude of direct wave for isotropic radiation (0.63) at free-surface (2)
                ! a0 = (2._r32 * ABS(mean(q)) * mean(mu) * mean(slip) * c) * 0.63_r32 * 2._r32
                !
                ! avepath = mean(path)
                ! avetrvt = mean(trvt)

                ! a0 = a0 / direct(, sheet, wtp)

                !   DO it = 1, npts
                !     rtri(it, 1) = rtri(it, 1) + envelope(it) * timeseries%cd%x(it, rec)
                !     rtri(it, 2) = rtri(it, 2) + envelope(it) * timeseries%cd%y(it, rec)
                !     rtri(it, 3) = rtri(it, 3) + envelope(it) * timeseries%cd%z(it, rec)
                !   ENDDO

              ENDDO     !< end loop over sheets

            ENDDO    !< end loop over wave types

            ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --
            ! --------------------------------------------- convolution with MRF -------------------------------------------------

#ifdef PERF
            CALL watch_start(tictoc, COMM)
#endif

            DO ic = 1, 3

              CALL fft(rtri(:, ic), tur)
              CALL fft(itri(:, ic), tui)

              ! multiply spectra
              DO it = 1, npts/2 + 1
                tur(it) = tur(it) * tv(it)
                tui(it) = tui(it) * tv(it)
              ENDDO

              CALL ifft(rtri(:, ic), tur)
              CALL ifft(itri(:, ic), tui)

            ENDDO

#ifdef PERF
            CALL watch_stop(tictoc, COMM)
            otimer(9) = otimer(9) + tictoc
#endif

            ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --
            ! -------------------------------------------- rotate & stack --------------------------------------------------------

            ! move back from "u/t" (fault-parallel/fault-normal) to "x/y" coordinates
            CALL rotate(rtri(:, 1), rtri(:, 2), rtri(:, 3), 0._r32, 0._r32, strike)
            CALL rotate(itri(:, 1), itri(:, 2), itri(:, 3), 0._r32, 0._r32, strike)

            DO ic = 1, 3
              DO it = 1, npts
                rseis(it, ic) = rseis(it, ic) + rtri(it, ic)
                iseis(it, ic) = iseis(it, ic) + itri(it, ic)
              ENDDO
            ENDDO

            m0 = m0 + mean(mu * slip)        !< average moment contribution (area is included later)

          ENDDO
        ENDDO
        !$omp end parallel do

        CALL dealloc_nodes()

        moment = moment + m0 * dutr(ref) * dvtr(ref) * 0.5_r32      !< sum is over mesh refinements

        CALL destroy_fftw_plan([npts])

        ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --
        ! ----------------------------------------------- Hilbert transform --------------------------------------------------------

#ifdef PERF
        CALL watch_start(tictoc, COMM)
#endif

        DO ic = 1, 3
          CALL hilbert(iseis(:, ic))
          DO it = 1, npts
            seis(it, ic) = seis(it, ic) + (rseis(it, ic) - iseis(it, ic)) * dt         !< dt is from convolution
          ENDDO
        ENDDO


#ifdef PERF
        CALL watch_stop(tictoc, COMM)
        timer(4) = timer(4) + tictoc
#endif

        IF (input%advanced%verbose .eq. 2) THEN

          CALL update_log(num2char('<cuts>', justify='c', width=30) +           &
                          num2char('P', width=15, justify='r') + '|' +  num2char('S', width=15, justify='r') + '|' +  &
                          num2char('Skipped', width=15, justify='r') + '|')

          CALL update_log(num2char('', width=30)  +  &
                          num2char(ncuts(1)/active(1), notation='f', width=15, precision=1, justify='r') + '|'  +   &
                          num2char(ncuts(2)/active(2), notation='f', width=15, precision=1, justify='r') + '|'  +   &
                          num2char(num2char(noslip) + ' (' +  &
                                   num2char(100*noslip/(totnutr*nvtr(ref))) + '%)', width=15, justify='r') + '|', blankline=.false.)

        ENDIF

print*, 'area ', dutr(ref) * dvtr(ref) * 0.5_r32

      ENDDO        !< end loop over mesh refinements

      scale = plane(pl)%targetm0 / moment          !< scaling factor to scale to desired moment

print*, 'scale ', scale

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! ------------------------------------- filter, differentiate, scale & stack -------------------------------------------------

#ifdef PERF
      CALL watch_start(tictoc, COMM)
#endif

      ALLOCATE(stack(SIZE(timeseries%sp%time)))

      ASSOCIATE(model => input%attenuation(1), fmax => input%coda%fmax)

        CALL make_iir_plan(ok, 'butter', dt, [model%lcut(band), MIN(model%hcut(band), fmax)], 'pass', 2, zphase = .true.)

        IF (ok .ne. 0) THEN
          CALL report_error(filter_error(ok))
          RETURN
        ENDIF

        DO ic = 1, 3

          seis(:, ic) = iir(seis(:, ic), ok)                       !< filter

          seis(:, ic) = differentiate(seis(:, ic), dt)             !< move from displacement to velocity

          CALL interpolate(time, seis(:, ic), timeseries%sp%time, stack)

          IF (ic .eq. 1) THEN
            timeseries%sp%x(:, rec) = timeseries%sp%x(:, rec) + stack(:) * scale
          ELSEIF (ic .eq. 2) THEN
            timeseries%sp%y(:, rec) = timeseries%sp%y(:, rec) + stack(:) * scale
          ELSE
            timeseries%sp%z(:, rec) = timeseries%sp%z(:, rec) + stack(:) * scale
          ENDIF

        ENDDO

        CALL destroy_iir_plan()

      END ASSOCIATE

#ifdef PERF
      CALL watch_stop(tictoc, COMM)
      timer(5) = timer(5) + tictoc
#endif


#ifdef PERF
      print*, 'Time RIK nodes: ', real(timer(1))
      print*, 'Time roughness: ', real(timer(2))
      print*, 'Time shooting: ', real(timer(3))

      print*, 'Time mechanism: ', real(otimer(1))
      print*, 'Time corners: ', real(otimer(2))
      print*, 'Time MRF + FFT: ', real(otimer(3))
      print*, 'Time ray corners: ', real(otimer(4))
      print*, 'Time kernel: ', real(otimer(5))
      ! print*, 'Time envelope: ', real(timer(9))
      print*, 'Time intg: ', real(otimer(6))
      print*, 'Time iir: ', real(otimer(8))
      print*, 'Time conv MRF: ', real(otimer(9))

      print*, 'Time Hilbert: ', real(timer(4))
      print*, 'Time Filter: ', real(timer(5))
#endif

    END SUBROUTINE solve_isochron_integral

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE time_integration(ref, it1, it2, dt, attenuation, tau, g, c, ncuts, rtri, itri)

      INTEGER(i32),                 INTENT(IN)  :: ref, it1, it2
      REAL(r32),                    INTENT(IN)  :: dt, attenuation
      REAL(r32),    DIMENSION(3),   INTENT(IN)  :: tau
      COMPLEX(r32), DIMENSION(3,3), INTENT(IN)  :: g
      REAL(r32),                    INTENT(OUT) :: c
      INTEGER(i32),                 INTENT(OUT) :: ncuts
      REAL(r32),    DIMENSION(:,:), INTENT(OUT) :: rtri, itri
      COMPLEX(r32)                              :: g31, g21, g23, ga, gb
      INTEGER(i32)                              :: ic, it, icut
      REAL(r32)                                 :: tau31, tau21, tau23, du, dv, p13, p12, p32, dl, t
      REAL(r32),    DIMENSION(3)                :: dtau

      !-----------------------------------------------------------------------------------------------------------------------------

      ncuts = 0

      ! differences in tau between corners
      tau31 = tau(3) - tau(1)
      tau21 = tau(2) - tau(1)
      tau23 = tau(2) - tau(3)

      ! apparent velocity of the tau curves
      c = HYPOT(tau31 / dutr(ref), (tau(2) - 0.5_r32 * (tau(1) + tau(3))) / dvtr(ref))
      c = 1._r32 / c

      ! loop over components of motion
      DO ic = 1, 3

        g31 = g(ic, 3) - g(ic, 1)
        g21 = g(ic, 2) - g(ic, 1)
        g23 = g(ic, 2) - g(ic, 3)

        ! loop over time (cuts)
        DO it = it1, it2

          t = (it - 1) * dt

          dtau(1) = t - tau(1)
          dtau(2) = t - tau(2)
          dtau(3) = t - tau(3)

          ! determine which corner is cut. Set "icut=0" when isochron is outside triangle
          IF (SIGN(1._r32, dtau(1)) .eq. SIGN(1._r32, dtau(2))) THEN
            IF (SIGN(1._r32, dtau(2)) .eq. SIGN(1._r32, dtau(3))) THEN
              icut = 0
            ELSE
              icut = 3
            ENDIF
          ELSE
            IF (SIGN(1._r32, dtau(1)) .eq. SIGN(1._r32, dtau(3))) THEN
              icut = 2
            ELSE
              icut = 1
            ENDIF
          ENDIF

          IF ( (dtau(1) .eq. 0._r32) .and. (dtau(2) .eq. 0._r32) ) icut = 3
          IF ( (dtau(1) .eq. 0._r32) .and. (dtau(3) .eq. 0._r32) ) icut = 2
          IF ( (dtau(2) .eq. 0._r32) .and. (dtau(3) .eq. 0._r32) ) icut = 1

          ! keep track of cuts
          IF ( (ic .eq. 1) .and. (icut .ne. 0) ) ncuts = ncuts + 1

          ! choose the correct integration formulae depending on which corner of the triangle is cut. ga and gb are the
          ! values of the integrand at the boundaries of the triangle where intersected by the integration contour.
          ! dl is the length of the contour within the triangle.

          SELECT CASE (icut)
          CASE(0)               !< isochron cuts a corner
            ga = 0._r32
            gb = 0._r32
            du = 0._r32
            dv = 0._r32
          CASE(1)
            p13 = dtau(1) / tau31
            p12 = dtau(1) / tau21
            ga  = g(ic, 1) + p13 * g31
            gb  = g(ic, 1) + p12 * g21
            dv  = dvtr(ref) * p12
            du  = dutr(ref) * (p13 - 0.5_r32 * p12)
          CASE(2)
            p12 = dtau(1) / tau21
            p32 = dtau(3) / tau23
            ga  = g(ic, 3) + p32 * g23
            gb  = g(ic, 1) + p12 * g21
            dv  = dvtr(ref) * (p12 - p32)
            du  = dutr(ref) * (1._r32 - 0.5_r32 * p12 - 0.5_r32 * p32)
          CASE(3)
            p13 = dtau(1) / tau31
            p32 = dtau(3) / tau23
            ga  = g(ic, 1) + p13 * g31
            gb  = g(ic, 3) + p32 * g23
            dv  = dvtr(ref) * p32
            du  = dutr(ref) * (1._r32 - p13 - 0.5_r32 * p32)
          END SELECT

          dl = HYPOT(du, dv)

          ! result of integration (real and imaginary parts)
          rtri(it, ic) = rtri(it, ic) + c * REAL(ga + gb, r32) * dl * 0.5_r32 * attenuation
          itri(it, ic) = itri(it, ic) + c * AIMAG(ga + gb)     * dl * 0.5_r32 * attenuation

        ENDDO      !< end loop over time (cuts)

      ENDDO      !< end loop over components

    END SUBROUTINE time_integration

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE integral_kernel(wtp, urec, vrec, u, v, velloc, srfvel, sd, cd, sr, cr, pn, fsp, fxp, fzp, g)

      INTEGER(i32),                 INTENT(IN)  :: wtp
      REAL(r32),                    INTENT(IN)  :: urec, vrec
      REAL(r32),    DIMENSION(3),   INTENT(IN)  :: u, v, velloc
      REAL(r32),                    INTENT(IN)  :: srfvel, sd, cd
      REAL(r32),    DIMENSION(3),   INTENT(IN)  :: sr, cr, pn
      REAL(r32),    DIMENSION(:),   INTENT(IN)  :: fsp
      COMPLEX(r32), DIMENSION(:),   INTENT(IN)  :: fxp, fzp
      COMPLEX(r32), DIMENSION(3,3), INTENT(OUT) :: g
      COMPLEX(r32)                              :: fxpfzs, fzpfxs
      COMPLEX(r32), DIMENSION(3,3)              :: fp, fs
      INTEGER(i32)                              :: icr, ic
      REAL(r32)                                 :: dx, dy, dr, cpsi, spsi, sloloc, p, signp, sthf, cthf, stho, ctho
      REAL(r32)                                 :: rn, rv, ru, bn, cn, bu, cu, bv, cv
      REAL(r32),    DIMENSION(3)                :: ir, itheta, iphi
      REAL(r32),    DIMENSION(3,3)              :: rvec, thvec, phvec

      !-----------------------------------------------------------------------------------------------------------------------------

      DO icr = 1, 3

        dx = urec - u(icr)
        dy = vrec - v(icr)

        dr = HYPOT(dx, dy)

        ! cos and sin of angle psi
        cpsi = dx / dr
        spsi = dy / dr

        sloloc = 1._r32 / velloc(icr)        !< slowness at corner

        p = pn(icr)

        ! "p" is "-p" for downgoing rays
        IF (p .gt. sloloc) p = -(2._r32 * sloloc - p)

        ! get sign of ray (+ if upgoing, - if downgoing)
        signp = SIGN(1._r32, p)

        p = ABS(p)         !< from here onwards we need only its absolute value

        ! sine of angle eta at source
        sthf = velloc(icr) * p

#ifdef ERROR_TRAP
        IF (sthf .gt. 1.01_r32) CALL report_error('solve_isochron_integral - ERROR: sin(eta) sensibly larger than 1')
#endif

        sthf = MIN(sthf, 1._r32)     !< handle roundoff errors

        ! cosine of angle eta at source
        cthf = signp * SQRT(1._r32 - sthf**2)

        ! dot products
        ! the fault normal points in the -y. direction for a 90 degree dip, and a positive u component of slip
        ! correspond to a left lateral strike slip motion
        rn = sthf * spsi * sd + cthf * cd          !< r.n (direction r and normal to plane)
        rv = cthf * sd - sthf * spsi * cd          !< r.v (directions r and dip-slip)
        ru = - sthf * cpsi                         !< r.u (directions r and along-strike)

        ! radiation pattern-dependent kernel stuff
        IF (wtp .eq. 1) THEN
          ir(icr) = 2._r32 * rn * (ru * cr(icr) - rv * sr(icr))
        ELSE
          ! dot products
          bn = (cthf * spsi * sd - sthf * cd)
          cn = (cpsi * sd)
          bu = -cthf * cpsi
          cu = spsi
          bv = -cthf * spsi * cd - sthf * sd
          cv = -cpsi * cd

          itheta(icr) = (rn * bu + bn * ru) * cr(icr) - (rn * bv + bn * rv) * sr(icr)
          iphi(icr)   = (rn * cu + cn * ru) * cr(icr) - (rn * cv + cn * rv) * sr(icr)
        ENDIF

        stho = p * srfvel

#ifdef ERROR_TRAP
        IF (stho .gt. 1.01_r32) CALL report_error('solve_isochron_integral - ERROR: sin(eta) sensibly larger than 1')
#endif

        stho = MIN(stho, 1._r32)          !< handle roundoff errors
        ctho = SQRT(1._r32 - stho**2)

        ! determine stuff that depends on the component of motion. rvec, thvec, and phvec are the components of the r,
        ! theta, and phi (unit vectors at the observer)
        rvec(1, icr)  = -stho * cpsi
        rvec(2, icr)  = -stho * spsi
        rvec(3, icr)  =  ctho
        thvec(1, icr) = -ctho * cpsi
        thvec(2, icr) = -ctho * spsi
        thvec(3, icr) = -stho
        phvec(1, icr) = spsi
        phvec(2, icr) = -cpsi
        phvec(3, icr) = 0._r32

        ! interpolate free-surface amplification coefficients at current ray parameter
        CALL interpolate(fsp, fxp, p, fxpfzs)
        CALL interpolate(fsp, fzp, p, fzpfxs)

        ! free-surface amplification coefficients for incident p (fp) and sv (fs) waves for each corner and component
        ! of motion (x,y,z)
        fp(1, icr) = fxpfzs
        fp(2, icr) = fxpfzs
        fp(3, icr) = fzpfxs
        fs(1, icr) = fzpfxs
        fs(2, icr) = fzpfxs
        fs(3, icr) = fxpfzs

        ! store integral kernel in "g"
        IF (wtp .eq. 1) THEN
          DO ic = 1, 3
            g(ic, icr) = ir(icr) * rvec(ic, icr) * fp(ic, icr)          !< P
          ENDDO
        ELSE
          DO ic = 1, 3
            g(ic, icr) = itheta(icr) * thvec(ic, icr) * fs(ic, icr) + iphi(icr) * phvec(ic, icr) * 2._r32    !< SV + SH
          ENDDO
        ENDIF

      ENDDO

    END SUBROUTINE integral_kernel

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE raypar_at_corners(sheet, shooting, repi, z, shot, wtp, p, path, trvt, q)

      INTEGER(i32),               INTENT(IN)    :: sheet
      REAL(r32),    DIMENSION(:), INTENT(IN)    :: shooting
      REAL(r32),    DIMENSION(3), INTENT(IN)    :: repi, z
      INTEGER(i32), DIMENSION(3), INTENT(IN)    :: shot
      INTEGER(i32),               INTENT(IN)    :: wtp
      REAL(r32),    DIMENSION(3), INTENT(OUT)   :: p, path, trvt, q
      INTEGER(i32)                              :: icr
      REAL(r32),    DIMENSION(2)                :: zshots, po, ro, to, qo

      !-----------------------------------------------------------------------------------------------------------------------------

      trvt(:) = 0._r32

      DO icr = 1, 3

        CALL interpray(sheet, repi(icr), shot(icr), wtp, po(1), ro(1), to(1), qo(1))        !< shot-point above
        CALL interpray(sheet, repi(icr), shot(icr) + 1, wtp, po(2), ro(2), to(2), qo(2))    !< shot-point below

        ! interpolate values at depth "z(icr)" only if corner is fully inside sheet
        IF (ALL(to .ne. 0._r32)) THEN

          zshots = [shooting(shot(icr)), shooting(shot(icr) + 1)]      !< shooting points depth vector

          CALL interpolate(zshots, po, z(icr), p(icr))
          CALL interpolate(zshots, ro, z(icr), path(icr))
          CALL interpolate(zshots, to, z(icr), trvt(icr))
          CALL interpolate(zshots, qo, z(icr), q(icr))

        ENDIF

      ENDDO

    END SUBROUTINE raypar_at_corners

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE lookup_fs(ok, rec, fsp, fxp, fzp)

      ! Purpose:
      !   to compute a lookup table for free-surface coefficients "fxp" and "fzp" at "fsp" points.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   08/03/21                  original version
      !

      INTEGER(i32),                            INTENT(OUT) :: ok
      INTEGER(i32),                            INTENT(IN)  :: rec
      REAL(r32),    ALLOCATABLE, DIMENSION(:), INTENT(OUT) :: fsp
      COMPLEX(r32), ALLOCATABLE, DIMENSION(:), INTENT(OUT) :: fxp, fzp
      INTEGER(i32)                                         :: n
      REAL(r32)                                            :: alpha, beta

      !-----------------------------------------------------------------------------------------------------------------------------

      alpha = input%velocity(input%receiver(rec)%velocity)%vp(1) / 1000._r32
      beta  = input%velocity(input%receiver(rec)%velocity)%vs(1) / 1000._r32

      n = (1._r32 / beta) / DP + 1      !< coefficients are computed at "n" slowness points

      ALLOCATE(fxp(n), fzp(n), fsp(n))

      ! determine coefficients
      CALL setfs(ok, alpha, beta, 1, fsp, fxp, fzp)

    END SUBROUTINE lookup_fs

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE hilbert(x)

      ! Purpose:
      !   to compute the Hilbert transform of timeseries "x".
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   08/03/21                  original version
      !

      USE, NON_INTRINSIC :: m_fft_cmplx

      REAL(r32),    DIMENSION(:),      INTENT(INOUT) :: x
      COMPLEX(r32), DIMENSION(SIZE(x))               :: analytic, spectrum
      INTEGER(i32)                                   :: i, n
      REAL(r32),    DIMENSION(SIZE(x))               :: h

      !-----------------------------------------------------------------------------------------------------------------------------

      n = SIZE(x)

      CALL make_fftw_plan([n])

      ! "h" is weighting function to compute analytic signal
      h(:) = 0._r32
      h(1) = 1._r32

      IF (MOD(n, 2) .eq. 0) THEN
        h(n/2 + 1) = 1._r32
        h(2:n/2)   = 2._r32
      ELSE
        h(2:(n + 1)/2) = 2._r32
      ENDIF

      CALL fft(x, spectrum)

      DO i = 1, n
        spectrum(i) = spectrum(i) * h(i)
      ENDDO

      CALL ifft(analytic, spectrum)

      DO i = 1, n
        x(i) = AIMAG(analytic(i))
      ENDDO

      CALL destroy_fftw_plan([n])

    END SUBROUTINE hilbert

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE setfs(ok, vp, vs, irort, fsp, fxp, fzp)

      INTEGER(i32),               INTENT(OUT) :: ok
      REAL(r32),                  INTENT(IN)  :: vp, vs
      INTEGER(i32),               INTENT(IN)  :: irort
      REAL(r32),    DIMENSION(:), INTENT(OUT) :: fsp
      COMPLEX(r32), DIMENSION(:), INTENT(OUT) :: fxp, fzp
      COMPLEX(r32)                            :: cosi, cosj, cicj, q, d, pupd, c, cs, cc
      INTEGER(i32)                            :: i, n
      REAL(r32)                               :: p, ps, g, gs, sini, sinj

      !-----------------------------------------------------------------------------------------------------------------------------

      ok = 0

      n = SIZE(fxp)

      ! compute free-surface amplification coefficients
      DO i = 1, n

        p  = (i - 1) * DP
        ps = p**2
        g  = 1._r32 / vs**2 - 2._r32 * ps
        gs = g**2

        sini = vp * p
        sinj = vs * p
        cosi = SQRT(CMPLX(1._r32 - sini**2, 0._r32))
        cosj = SQRT(CMPLX(1._r32 - sinj**2, 0._r32))

        IF (AIMAG(cosi) .lt. 0._r32) THEN
          cosi = CMPLX(REAL(cosi, r32), -AIMAG(cosi))
        ENDIF
        IF (AIMAG(cosj) .lt. 0._r32) THEN
          cosj = CMPLX(REAL(cosj), -AIMAG(cosj))
        ENDIF

        cicj = cosi * cosj

        q = 4._r32 * ps * cicj / vp / vs

        ! d is the Rayleigh denominator
        d = gs + q
        c = 4._r32 / vp / vs * g / d

        cs = c * sini * sinj
        cc = c * cicj

        pupd = (q - gs) / d

        fxp(i) = 1._r32 + pupd + cc
        fzp(i) = 1._r32 - pupd + cs

        fsp(i) = p

        ! check that for p < 1/vp, all coefficients must be real
        IF (p .le. 1._r32/vp) THEN
          IF ( (AIMAG(fxp(i)) .ne. 0._r32) .or. (AIMAG(fzp(i)) .ne. 0._r32) ) THEN
            CALL report_error('setfs - ERROR: coefficient fxp and/or fzp not purely real for p < 1/vp')
            ok = 1
            RETURN
          ENDIF
        ENDIF

      ENDDO

      ! smooth coefficients
      CALL rnmnc(ok, fxp, irort)

      IF (ok .ne. 0) RETURN

      CALL rnmnc(ok, fzp, irort)

      IF (ok .ne. 0) RETURN

    END SUBROUTINE setfs

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE rnmnc(ok, a, irort)

      INTEGER(i32),                            INTENT(OUT)   :: ok
      COMPLEX(r32),              DIMENSION(:), INTENT(INOUT) :: a
      INTEGER(i32),                            INTENT(IN)    :: irort
      COMPLEX(r32)                                           :: adif, sum
      COMPLEX(r32), ALLOCATABLE, DIMENSION(:)                :: ansr
      INTEGER(i32)                                           :: i, j, na, nrmn, ntot, nend, nmid, nshift
      REAL(r32)                                              :: wsum
      REAL(r32),    ALLOCATABLE, DIMENSION(:)                :: w

      ! ----------------------------------------------------------------------------------------------------------------------------

      ok = 0

      na = SIZE(a)

      nrmn = DPSM / DP + 2          !< points in the running mean window

      ntot = na + 2 * (nrmn - 1)

      IF (nrmn .le. 0) THEN
        CALL report_error('rnmc - ERROR: running mean window too short')
        ok = 1
        RETURN
      ENDIF

      IF (nrmn .ge. na) THEN
        CALL report_error('rnmc - ERROR: running mean window too long')
        ok = 1
        RETURN
      ENDIF

      IF ( (irort .lt. 0) .and. (irort .gt. 1) ) THEN
        CALL report_error('rnmc - ERROR: unsupported window')
        ok = 1
        RETURN
      ENDIF

      ALLOCATE(ansr(ntot), w(nrmn))

      ! load a into answer, shifted
      DO i = 1, na
        ansr(nrmn - 1 + i) = a(i)
      ENDDO

      ! determine end values by linear interpolation
      adif = a(2) - a(1)
      DO  i = 1, nrmn - 1
        ansr(i) = a(1) + adif * (i - nrmn)
      ENDDO

      nend = na + nrmn - 1
      adif = a(na) - a(na - 1)
      DO i = na + nrmn, ntot
        ansr(i) = a(na) + adif * (i - nend)
      ENDDO

      ! running mean
      IF (irort .eq. 0) THEN

        ! rectangular weighting
        DO i = 1, nend
          sum = 0._r32
          DO j = 1, nrmn
            sum = sum + ansr(i + j - 1)
          ENDDO
          ansr(i) = sum / nrmn
        ENDDO

      ELSE

        ! triangular weighting.  figure out weights
        nmid = (nrmn + 1) / 2
        DO i = 1, nmid
          w(i) = 2._r32 * i - 1
        ENDDO

        IF (nmid .lt. nrmn) THEN
          DO i = nmid + 1, nrmn
            w(i) = w(nrmn - i + 1)
          ENDDO
        ENDIF

        ! normalize weights to 1
        wsum = 0._r32
        DO i = 1, nrmn
          wsum = wsum + w(i)
        ENDDO
        DO i = 1, nrmn
          w(i) = w(i) / wsum
        ENDDO

        ! running mean triangular weighting
        DO i = 1, nend
          sum = 0._r32
          DO j = 1, nrmn
            sum = sum + ansr(i + j - 1) * w(j)
          ENDDO
          ansr(i) = sum
        ENDDO

      ENDIF

      ! shift to centered windows.  exact for nrmn odd
      nshift = (nrmn - 1) / 2
      IF (nshift .gt. 0) THEN
        DO i = 1, na
          ansr(i) = ansr(i + nshift)
        ENDDO
      ENDIF

      ! return only first na elements
      DO i = 1, na
        a(i) = ansr(i)
      ENDDO

    END SUBROUTINE rnmnc

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE node2disk(ok, pl, vel, iter)

      ! Purpose:
      !   to write to disk the rupture parameters defined for plane "pl" embedded in velocity model "vel" and iteration "iter".
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   08/03/21                  original version
      !

      INTEGER(i32),                          INTENT(OUT) :: ok
      INTEGER(i32),                          INTENT(IN)  :: pl, vel, iter
      CHARACTER(:), ALLOCATABLE                          :: fo
      INTEGER(i32)                                       :: i, j, ref, lu, icr, totnutr, seed
      INTEGER(i32),             DIMENSION(3)             :: iuc, ivc
      REAL(r32)                                          :: m0, moment, strike, dip, rho, beta
      REAL(r32),                DIMENSION(3)             :: u, v, w, x, y, z, slip, rise, rupture, rake, nrl, mu

      !-----------------------------------------------------------------------------------------------------------------------------

      CALL setup_interpolation('linear', 'zero', ok)

      IF (ok .ne. 0) RETURN

      fo = 'node_' + num2char(pl) + '_' + num2char(vel) + '_' + num2char(iter) + '.bin'

      OPEN(newunit = lu, file = fo, status = 'replace', form = 'unformatted', access = 'stream', action = 'write', IOSTAT = ok)

      IF (ok .ne. 0) THEN
        CALL report_error('Error while opening file' + TRIM(fo))
        RETURN
      ENDIF

      WRITE(lu, POS=1) SUM([nvtr * (2 * nutr - 1)])         !< total number of triangles (all refinements)

      WRITE(lu) SIZE(nvtr)

      moment = 0._r32

      nodefun => rik_at_nodes

      IF (input%source%is_point) nodefun => ptrsrc_at_nodes

      ! set "seed" such that random numbers depend on fault plane number and iteration
      seed = input%source%seed + (iter - 1) * SIZE(plane) + pl

      DO ref = 1, SIZE(nvtr)         !< loop over mesh refinements

        ALLOCATE(nodes(nugr(ref), nvgr(ref)))

        CALL nodefun(ref, pl, vel, seed)                !< define slip, rupture time and rise time on mesh nodes

        totnutr = 2 * nutr(ref) - 1                     !< total triangles in a row

        m0 = 0._r32

        WRITE(lu) nvtr(ref), nutr(ref)

        CALL fault_roughness(ok, ref, pl, iter)

        ! CALL rayshooting(ok, ref, pl, vel)
        ! !stop

        !$omp parallel do ordered default(shared) private(i, j, iuc, ivc, slip, rise, rake, nrl, rupture, icr, u, v, w, x, y, z)  &
        !$omp& private(beta, rho, mu) reduction(+:m0) collapse(2)
        DO j = 1, nvtr(ref)
          DO i = 1, totnutr

            CALL cornr(j, i, iuc, ivc)        !< corner indices for current triangle

            DO icr = 1, 3
              slip(icr) = SUM(nodes(iuc(icr), ivc(icr))%slip)
              rise(icr) = mean(nodes(iuc(icr), ivc(icr))%rise)
              rupture(icr) = MINVAL(nodes(iuc(icr), ivc(icr))%rupture)
            ENDDO

            CALL cornr2uv(iuc, ivc, ref, u, v)        !< on-fault coordinates

            DO icr = 1, 3
              w(icr) = roughness(iuc(icr), ivc(icr))       !< this is always zero for point-sources
            ENDDO

            CALL uvw2xyz(pl, u, v, w, x, y, z)            !< get cartesian coordinates

            DO icr = 1, 3
              z(icr) = MAX(MIN_DEPTH, z(icr))             !< make sure we never breach the free-surface (clip roughness)
            ENDDO

            !$omp ordered
            WRITE(lu) x, y, z, slip, rise, rupture         !< values at each corner
            !$omp end ordered

            CALL normal2tri(x, y, z, nrl)

            ! always have normal pointing upward (i.e. negative)
            IF (MOD(i + (j-1)*totnutr, 2) == 0) nrl(:) = -nrl(:)

            strike = plane(pl)%strike * DEG_TO_RAD
            dip    = plane(pl)%dip    * DEG_TO_RAD

            CALL interpolate(plane(pl)%u, plane(pl)%v, plane(pl)%rake, u, v, rake)

            CALL perturbed_mechanism(strike, dip, rake, nrl)            !< all values in radians

            !$omp ordered
            WRITE(lu) strike, dip, rake, w     !< strike and dip are scalar
            !$omp end ordered

            DO icr = 1, 3
              beta    = vinterp(input%velocity(vel)%depth, input%velocity(vel)%vs, input%velocity(vel)%vsgrad, z(icr))
              rho     = vinterp(input%velocity(vel)%depth, input%velocity(vel)%rho, input%velocity(vel)%rhograd, z(icr))
              mu(icr) = rho * beta**2
            ENDDO

            m0 = m0 + mean(mu * slip)

          ENDDO
        ENDDO
        !$omp end parallel do

        CALL dealloc_nodes()

        moment = moment + m0 * dutr(ref) * dvtr(ref) * 0.5_r32

      ENDDO        !< end loop over mesh refinements

      WRITE(lu) plane(pl)%targetm0 / moment

      CLOSE(lu, iostat = ok)

      IF (ok .ne. 0) THEN
        CALL report_error('Error while closing file ' + fo)
        RETURN
      ENDIF

    END SUBROUTINE node2disk

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE perturbed_mechanism(strike, dip, rake, nrl)

      ! Purpose:
      !   to return strike dip and rake of a triangle (described by its normal "nrl") whose initial orientation has been changed and
      !   now described by its normal "nrl". In input, "strike" and "dip" are two scalars and "rake" a vector representing the
      !   initial values at each corner, while in output they refer to the new orientation.
      !   The algorithm first computes the normal and slip vectors for the unperturbed case (as given by strike, dip and rake). These
      !   are used to define the pressure axis P and then (by plugging in the new normal) the new rake angle.
      !   WARNING: input values are expected to follow the N, E, D coordinates system, with dip <= 90. Internally the E, N, U system
      !            is assumed. "nrl" must be normalized to 1.
      !   Output values are in the range [0, pi] or [-pi, 0] (slip and rake) or [0, pi] dip.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   08/03/21                  original version
      !

      REAL(r32),                 INTENT(INOUT) :: strike, dip
      REAL(r32),   DIMENSION(3), INTENT(INOUT) :: rake
      REAL(r32),   DIMENSION(3), INTENT(IN)    :: nrl
      INTEGER(i32)                             :: i
      REAL(r32)                                :: p, cs, ss, cd, sd, cr, sr, arg, hn2
      REAL(r32),   DIMENSION(3)                :: unrl, ul, l, h

      !-----------------------------------------------------------------------------------------------------------------------------

      ! strike vector (perturbed case)
      h(1) = -nrl(1)   !< east
      h(2) = nrl(2)    !< north
      h(3) = 0._r32

      hn2 = NORM2(h)

      cs = COS(strike)
      ss = SIN(strike)
      cd = COS(dip)
      sd = SIN(dip)

      ! normal vector (unperturbed case)
      unrl(1) = cs * sd     !< east
      unrl(2) = -ss * sd    !< north
      unrl(3) = cd          !< up

      DO i = 1, 3

        cr = COS(rake(i))
        sr = SIN(rake(i))

        ! slip vector (unperturbed case)
        ul(1) = ss * cr - cs * cd * sr    !< east
        ul(2) = cs * cr + ss * cd * sr    !< north
        ul(3) = sd * sr                   !< up

        ! slip vector for triangle with normal "nrl", where P axis is based on unperturbed geometry
        ! sqrt(2)*P = unrl - ul
        ! l = nrl - sqrt(2)*P = nrl - (unrl - ul)
        ! note that we reverted (1,2) indices for input normal as it assumes N-E-D reference system
        l(1) = nrl(2) - (unrl(1) - ul(1))
        l(2) = nrl(1) - (unrl(2) - ul(2))
        l(3) = -nrl(3) - (unrl(3) - ul(3))

        p = NORM2(l) * hn2

        arg = (h(1)*l(1) + h(2)*l(2)) / p

        arg = SIGN(MIN(1._r32, ABS(arg)), arg)          !< prevent roundoff errors, enforcing |arg| <= 1

        ! sign of rake angle depends on vertical component of slip vector
        rake(i) = SIGN(ACOS(arg), l(3))       !< values in range [0 pi] or [-pi 0]

      ENDDO

      ! remember that ATAN2(y,x) is in the interval [0,pi] when y>=0, in the interval (-pi, 0) otherwise

      ! perturbed strike in the range [0 pi] or [-pi 0]
      strike = ATAN2(-REAL(nrl(1), r64), REAL(nrl(2), r64))

      ! IF (strike .lt. 0._r32) strike = 2._r32 * PI + strike

      ! perturbed dip in the range [0 pi]
      dip = ATAN2(HYPOT(REAL(nrl(1), r64), REAL(nrl(2), r64)), -REAL(nrl(3), r64))    !< revert sign vertical because input reference system is N-E-D

    END SUBROUTINE perturbed_mechanism

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE rayshooting(ok, ref, pl, vel)

      ! Purpose:
      !   to compute a set of ray-related quantities (ray parameter, horizontal distance, traveltime, ray amplitude) for a set of
      !   point-like sources covering fault plane "pl" embedded in velocity model "vel". These quantities refer to direct P- and S-
      !   waves.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   08/03/21                  original version
      !

      INTEGER(i32),                           INTENT(OUT) :: ok
      INTEGER(i32),                           INTENT(IN) :: ref, pl, vel
      CHARACTER(:), ALLOCATABLE                          :: fo
      INTEGER(i32)                                       :: iface, i, j, wavetype, rec, sheets, totnutr, icr, layers, lu
      INTEGER(i32),              DIMENSION(3)            :: iuc, ivc
      REAL(r32)                                          :: zmin, zmax, delta, vtop, vbottom, vs, dz, zo, alpha, beta, rho, rmax
      REAL(r32),                              PARAMETER  :: GAP = 10._r32
      REAL(r32),                 DIMENSION(3)            :: u, v, w, x, y, z
      REAL(r32),    ALLOCATABLE, DIMENSION(:)            :: distance
      REAL(r32),    ALLOCATABLE, DIMENSION(:,:)          :: velocity
      REAL(r64),                 DIMENSION(1)            :: tictoc


      !-----------------------------------------------------------------------------------------------------------------------------

      ok = 0

      IF (ALLOCATED(shooting)) DEALLOCATE(shooting)
      IF (ALLOCATED(wkbj)) DEALLOCATE(wkbj)

      delta = MAXVAL(ABS(roughness)) * COS(plane(pl)%dip * DEG_TO_RAD)

      w = [0._r32, 0._r32, 0._r32]

      CALL cornr(1, 1, iuc, ivc)                !< corner indices for first triangle first row
      CALL cornr2uv(iuc, ivc, ref, u, v)        !< on-fault coordinates
      CALL uvw2xyz(pl, u, v, w, x, y, z)        !< get cartesian coordinates

      zmin = MINVAL(z)

      CALL cornr(nvtr(ref), 1, iuc, ivc)        !< corner indices for first triangle last row
      CALL cornr2uv(iuc, ivc, ref, u, v)        !< on-fault coordinates
      CALL uvw2xyz(pl, u, v, w, x, y, z)        !< get cartesian coordinates

      zmax = MAXVAL(z)

      ASSOCIATE(model => input%velocity(vel))

        vtop    = vinterp(model%depth, model%vs, model%vsgrad, zmin)        !< vs at upper edge
        vbottom = vinterp(model%depth, model%vs, model%vsgrad, zmax)        !< vs at lower edge

        vs = MIN(vtop, vbottom)     !< min shear wave speed over current mesh

        dz = vs / input%coda%fmax / input%advanced%pmw        !< desired spacing between shooting points
        ! dz = vs / input%coda%fmax

        ! add max deviation from planarity
        zmin = MAX(MIN_DEPTH, zmin - delta)                   !< don't go too close free-surface (i.e. above MIN_DEPTH)
        zmax = zmax + delta

        ALLOCATE(distance(SIZE(model%depth)))

        ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --
        ! ---------------------------------------------- define shooting points ----------------------------------------------------

        ! handle shallowest shooting point
        distance = ABS(zmin - model%depth)
        iface    = MINLOC(distance, DIM=1)

        ! shooting points cannot coincide with velocity discontinuities
        IF (distance(iface) .lt. GAP) zmin = MAX(MIN_DEPTH, model%depth(iface) - GAP*2)   !< move it just above interface if too close

        shooting = [zmin]

        zo = MIN(zmin + dz, zmax)            !< move to following source depth (can be directly "zmax" if very close to "zmin")

        DO

          distance = ABS(zo - model%depth)
          iface    = MINLOC(distance, DIM=1)

          IF (distance(iface) .lt. GAP) zo = model%depth(iface) + GAP*2        !< move source just below interface

          shooting = [shooting, zo]        !< append shooting point

          IF (zo .ge. zmax) EXIT           !< we have covered depth interval as needed

          zo = zo + dz

        ENDDO

      END ASSOCIATE

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! -------------------------------------------------- max s/r distance --------------------------------------------------------

      totnutr = 2 * nutr(ref) - 1                     !< total triangles in a row

      rmax = 0._r32

      !$omp parallel do default(shared) private(i, j, iuc, ivc, u, v, w, x, y, z, rec) reduction(max: rmax)
      DO j = 1, nvtr(ref)
        DO i = 1, totnutr

          CALL cornr(j, i, iuc, ivc)                 !< corner indices for current triangle
          CALL cornr2uv(iuc, ivc, ref, u, v)         !< on-fault coordinates

          DO icr = 1, 3
            w(icr) = roughness(iuc(icr), ivc(icr))
          ENDDO

          CALL uvw2xyz(pl, u, v, w, x, y, z)            !< get cartesian coordinates

          ! compute max source-receiver distance
          DO icr = 1, 3
            DO rec = 1, SIZE(input%receiver)
              rmax = MAX(rmax, HYPOT(x(icr) - input%receiver(rec)%x, y(icr) - input%receiver(rec)%y))
            ENDDO
          ENDDO

        ENDDO
      ENDDO
      !$omp end parallel do

      IF (input%advanced%verbose .eq. 2) THEN

        zmin = BIG
        zmax = -BIG

        DO i = 2, SIZE(shooting)
          zo   = shooting(i) - shooting(i - 1)
          zmin = MIN(zmin, zo)
          zmax = MAX(zmax, zo)
        ENDDO

        CALL update_log(num2char('<ray shooting>', justify='c', width=30) + num2char('Points', width=15, justify='r') + '|' +  &
                        num2char('Rmax', width=15, justify='r') + '|' +  &
                        num2char('Zmin', width=15, justify='r') + '|' + num2char('Zmax', width=15, justify='r') + '|' +  &
                        num2char('Avg Dz', width=15, justify='r') + '|')

        CALL update_log(num2char('', width=30, justify='c')  +  &
                        num2char(SIZE(shooting), width=15, justify='r') + '|' + &
                        num2char(rmax, width=15, notation='f', precision=1, justify='r') + '|' + &
                        num2char(MINVAL(shooting), width=15, notation='f', precision=2, justify='r') + '|' + &
                        num2char(MAXVAL(shooting), width=15, notation='f', precision=2, justify='r') + '|' +  &
                        ! num2char(num2char(zmin, notation='f', width=6, precision=2) + ', ' +   &
                        !          num2char(zmax, notation='f', width=6, precision=2), width=15, justify='r') + '|',blankline=.false.)
                        num2char(mean([zmin, zmax]), notation='f', width=15, precision=2, justify='r') + '|', blankline=.false.)

      ENDIF

      rmax = 1.2 * rmax / 1000._r32            !< slightly increase maximum distance

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! --------------------------------------------- prepare velocity model -------------------------------------------------------

      ASSOCIATE(model => input%velocity(vel))

        layers = SIZE(model%depth)

        ALLOCATE(velocity(4, 2*layers))

        DO i = 1, layers

          IF (i .gt. 1) THEN

            zo = model%depth(i) - GAP

            alpha = vinterp(model%depth, model%vp, model%vpgrad, zo)
            beta  = vinterp(model%depth, model%vs, model%vsgrad, zo)
            rho   = vinterp(model%depth, model%rho, model%rhograd, zo)

            velocity(:, (i - 1)*2) = [zo, alpha, beta, rho] / 1000._r32         !< slightly above interface

          ENDIF

          zo = model%depth(i) + GAP

          IF (i .eq. 1) zo = 0._r32

          alpha = vinterp(model%depth, model%vp, model%vpgrad, zo)
          beta  = vinterp(model%depth, model%vs, model%vsgrad, zo)
          rho   = vinterp(model%depth, model%rho, model%rhograd, zo)

          velocity(:, (i - 1)*2 + 1) = [zo, alpha, beta, rho]  / 1000._r32      !< slightly below interface

        ENDDO

        zo = 100000._r32

        alpha = vinterp(model%depth, model%vp, model%vpgrad, zo)
        beta  = vinterp(model%depth, model%vs, model%vsgrad, zo)
        rho   = vinterp(model%depth, model%rho, model%rhograd, zo)

        velocity(:, 2*layers) = [zo, alpha, beta, rho] / 1000._r32

      END ASSOCIATE

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! ----------------------------------------------------- shoot rays -----------------------------------------------------------

#ifdef DEBUG
      fo = 'spred_' + num2char(ref) + '_' + num2char(pl) + '_' + num2char(vel) + '.txt'

      OPEN(newunit = lu, file = TRIM(fo), status = 'unknown', form = 'formatted', access = 'sequential', action = 'write',  &
           iostat = ok)

      IF (ok .ne. 0) THEN
        CALL report_error('Error while opening file' + TRIM(fo))
        RETURN
      ENDIF
#endif

      ALLOCATE(wkbj(SIZE(shooting), 2))

#ifdef PERF
      CALL watch_start(tictoc(1), COMM)
#endif

      DO wavetype = 1, 2

        maxsheets(wavetype) = 0
        minsheets(wavetype) = HUGE(0)

        DO i = 1, SIZE(shooting)

          CALL spred(ok, velocity, rmax, shooting(i) / 1000._r32, i, wavetype, lu, sheets)

          maxsheets(wavetype) = MAX(maxsheets(wavetype), sheets)     !< max number of sheets for each wave type
          minsheets(wavetype) = MIN(minsheets(wavetype), sheets)

        ENDDO

      ENDDO

#ifdef PERF
      CALL watch_stop(tictoc(1), COMM)
#endif

#ifdef DEBUG
      CLOSE(lu, iostat = ok)

      IF (ok .ne. 0) THEN
        CALL report_error('Error while closing file ' + fo)
        RETURN
      ENDIF
#endif

      IF (input%advanced%verbose .eq. 2) THEN

        CALL update_log(num2char('<sheets (min, max)>', justify='c', width=30) +       &
                        num2char('P-wave', width=15, justify='r') + '|' + num2char('S-wave', width=15, justify='r') + '|')

        CALL update_log(num2char('', width=30, justify='c')  +  &
                        num2char(num2char(minsheets(1), width=2) + ', ' +   &
                                 num2char(maxsheets(1), width=2), width=15, justify='r') + '|' + &
                        num2char(num2char(minsheets(2), width=2) + ', ' +   &
                                 num2char(maxsheets(2), width=2), width=15, justify='r') + '|',  blankline=.false.)

      ENDIF

! #ifdef PERF
!       IF (input%advanced%verbose .eq. 2) THEN
!         CALL update_log(num2char('<<elapsed time>>', justify='c', width=30) + num2char(tictoc(1), width=15, notation='f',   &
!                         precision=3, justify='r') + '|', blankline=.false.)
!       ENDIF
! #endif

    END SUBROUTINE rayshooting

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE correct4impz(ok, iter)

      USE, NON_INTRINSIC :: m_fft_real

      INTEGER(i32),                           INTENT(OUT) :: ok
      INTEGER(i32),                           INTENT(IN)  :: iter
      CHARACTER(:), ALLOCATABLE                           :: fo
      COMPLEX(r32), ALLOCATABLE, DIMENSION(:)             :: spectrum
      INTEGER(i32)                                        :: band, rec, i, lu, npts
      REAL(r32),    ALLOCATABLE, DIMENSION(:)             :: freq, amp, phase, respz, blob

      !-----------------------------------------------------------------------------------------------------------------------------

      ASSOCIATE(model => input%attenuation(1), fmax => input%coda%fmax, dt => timeseries%sp%dt)

        ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --
        ! ----------------------------------------- overall filter response --------------------------------------------------------

        DO band = 1, SIZE(model%lcut)

          IF (model%lcut(band) .ge. fmax) EXIT

          CALL make_iir_plan(ok, 'butter', dt, [model%lcut(band), MIN(model%hcut(band), fmax)], 'pass', 2, zphase = .true.)

! print*, 'correct4impz ', ok, [model%lcut(1), MIN(model%hcut(band), fmax)], dt

          CALL freqz(freq, amp, phase, dt, SIZE(timeseries%sp%time))

          IF (band .eq. 1) THEN
            respz = amp**2
          ELSE
            respz = respz + amp**2
          ENDIF

          CALL destroy_iir_plan()

        ENDDO

        band = band - 1

        ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --
        ! ----------------------------------------- compensate filter response -----------------------------------------------------

        IF (band .gt. 1) THEN

          ! we want the composite filter reponse to look identical to that of a filter between lowest fcut and highest fcut(/fmax)
          CALL make_iir_plan(ok, 'butter', dt, [model%lcut(1), MIN(model%hcut(band), fmax)], 'pass', 2, zphase = .true.)

! print*, 'correct4impz ', ok, [model%lcut(1), MIN(model%hcut(band), fmax)], dt

          CALL freqz(freq, amp, phase, dt, SIZE(timeseries%sp%time))

          CALL destroy_iir_plan()

          npts = SIZE(timeseries%sp%time)

          blob = respz

          ! correct filter response for non-negligibile amplitude levels
          DO i = 1, SIZE(respz)
            IF (respz(i) .gt. 1.E-20) respz(i) = amp(i)**2 / respz(i)    !< **2 is because we used two-pass filters
          ENDDO

#ifdef DEBUG

          fo = 'respz_' + num2char(iter) + '.txt'

          OPEN(newunit = lu, file = fo, status = 'replace', form = 'formatted', access = 'sequential', action = 'write', iostat= ok)

          IF (ok .ne. 0) THEN
            CALL report_error('Error while opening file' + TRIM(fo))
            RETURN
          ENDIF

          DO i = 1, SIZE(respz)
            WRITE(lu, *) freq(i), blob(i), respz(i), amp(i)**2       !< blob*respz should be equal to amp**2
          ENDDO

          CLOSE(lu, iostat = ok)

          IF (ok .ne. 0) THEN
            CALL report_error('Error while closing file ' + fo)
            RETURN
          ENDIF

#endif

          CALL make_fftw_plan([npts])

          ALLOCATE(spectrum(npts/2 + 1))

          ! apply spectral correction function
          DO rec = 1, SIZE(input%receiver)

            CALL fft(timeseries%sp%x(:, rec), spectrum)
            spectrum(:) = spectrum(:) * respz(:)
            CALL ifft(timeseries%sp%x(:, rec), spectrum)

            CALL fft(timeseries%sp%y(:, rec), spectrum)
            spectrum(:) = spectrum(:) * respz(:)
            CALL ifft(timeseries%sp%y(:, rec), spectrum)

            CALL fft(timeseries%sp%z(:, rec), spectrum)
            spectrum(:) = spectrum(:) * respz(:)
            CALL ifft(timeseries%sp%z(:, rec), spectrum)

          ENDDO

          CALL destroy_fftw_plan([npts])

        ENDIF

      END ASSOCIATE


    END SUBROUTINE correct4impz

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE lookup_coda(ok, ref, rec, band, pl, vel, time)

      INTEGER(i32),                             INTENT(OUT) :: ok
      INTEGER(i32),                             INTENT(IN)  :: ref, rec, band, vel, pl
      REAL(r32),                 DIMENSION(:),  INTENT(IN)  :: time
      CHARACTER(:), ALLOCATABLE                             :: fo
      INTEGER(i32)                                          :: totnutr, i, j, icr, src, wave, sheet, c, lu, n
      INTEGER(i32),              DIMENSION(3)               :: iuc, ivc, shot
      REAL(r32)                                             :: avg, weight, r, tp, ts, vs, wp, ws, gsp, t, vratio
      REAL(r32),                 DIMENSION(2,2)             :: pbounds, tbounds
      REAL(r32),                 DIMENSION(3)               :: p, path, trvt, q, u, v, w, x, y, z, repi

      !-----------------------------------------------------------------------------------------------------------------------------

      ok = 0

      totnutr = 2 * nutr(ref) - 1                     !< total triangles in a row

      pbounds(1,:) = BIG
      pbounds(2,:) = -BIG

      tbounds(1,:) = BIG
      tbounds(2,:) = -BIG

      vratio = 0._r32
      n = 0

!     !$omp parallel do default(shared) private(i, j, iuc, ivc, u, v, w, icr, src, x, y, z, rec, shot, repi, p, path, trvt)  &
!     !$omp private(q, weight, c, avg) reduction(max: pbounds(2), tbounds(2)) reduction(min: pbounds(1), tbounds(1)) &
!     !$omp reduction(sum: vratio, n)
      DO j = 1, nvtr(ref)
        DO i = 1, totnutr

          CALL cornr(j, i, iuc, ivc)                 !< corner indices for current triangle
          CALL cornr2uv(iuc, ivc, ref, u, v)         !< on-fault coordinates

          DO icr = 1, 3
            w(icr) = roughness(iuc(icr), ivc(icr))
          ENDDO

          CALL uvw2xyz(pl, u, v, w, x, y, z)            !< get cartesian coordinates

          DO icr = 1, 3
            z(icr) = MAX(MIN_DEPTH, z(icr))             !< make sure we never breach the free-surface (clip roughness)
          ENDDO

          ! find shooting points above/below
          DO icr = 1, 3
            DO src = 1, SIZE(shooting) - 1
              IF ( (z(icr) .ge. shooting(src)) .and. (z(icr) .le. shooting(src + 1)) ) shot(icr) = src
            ENDDO
          ENDDO

          DO icr = 1, 3
            repi(icr) = HYPOT(input%receiver(rec)%x - x(icr), input%receiver(rec)%y - y(icr)) / 1000._r32        !< epicentral distance (corner-receiver)
          ENDDO

          DO sheet = 1, MAXVAL(maxsheets)

            DO wave = 1, 2

              CALL raypar_at_corners(sheet, shooting, repi, z, shot, wave, p(:, wave), path(:, wave), trvt(:, wave), q(:, wave))

              IF (ALL(trvt(:, wave) .ne. 0._r32)) THEN
                pbounds(1, wave) = MIN(pbounds(1, wave), mean(path(:, wave)))
                pbounds(2, wave) = MAX(pbounds(2, wave), mean(path(:, wave)))
                tbounds(1, wave) = MIN(tbounds(1, wave), mean(trvt(:, wave)))
                tbounds(2, wave) = MAX(tbounds(2, wave), mean(trvt(:, wave)))
              ENDIF

            ENDDO

            IF (ALL(trvt) .ne. 0._r32) THEN
              vratio = vratio + mean(path(:, 1)) / mean(trvt(:, 1)) * mean(trvt(:, 2)) / mean(path(:, 2))   !< alpha/beta
              n      = n + 1
            ENDDO

          ENDDO

        ENDDO
      ENDDO  !< end loop over triangles
!     !$omp end parallel do

      vratio = vratio / n         !< average alpha/beta ratio

      coda%lopath = pbounds(1, :)
      coda%uppath = pbounds(2, :)

      coda%lotrvt = tbounds(1, :)
      coda%uptrvt = tbounds(2, :)

      IF (input%advanced%verbose .eq. 2) THEN
        CALL update_log(num2char('<coda tableau>', justify='c', width=30) +    &
                        num2char('Min trvt', width=15, justify='r') + '|' + num2char('Max trvt', width=15, justify='r') + '|' + &
                        num2char('Min path', width=15, justify='r') + '|' + num2char('Max path', width=15, justify='r') + '|' + &
                        num2char('<vp/vs>', width=15, justify='r')  + '|')

        CALL update_log(num2char('', width=30) +  &
                        num2char(coda%lotrvt(sheet, 1), width=15, justify='r', notation='f', precision=2) + '|' +    &
                        num2char(coda%uptrvt(sheet, 1), width=15, justify='r', notation='f', precision=2) + '|' +    &
                        num2char(coda%counts(sheet, 1), width=15, justify='r') + '|', blankline=.false.)





      IF (input%advanced%verbose .eq. 2) THEN
        CALL update_log(num2char('<coda P>', justify='c', width=30) + num2char('Sheet', width=15, justify='r') + '|' + &
                        num2char('Min trvt', width=15, justify='r') + '|' + num2char('Max trvt', width=15, justify='r') + '|' + &
                        num2char('Envelopes', width=15, justify='r') + '|')

        DO sheet = 1, maxsheets(1)
          CALL update_log(num2char('', width=30) + num2char(sheet, width=15, justify='r')   + '|' +    &
                          num2char(coda%lotrvt(sheet, 1), width=15, justify='r', notation='f', precision=2) + '|' +    &
                          num2char(coda%uptrvt(sheet, 1), width=15, justify='r', notation='f', precision=2) + '|' +    &
                          num2char(coda%counts(sheet, 1), width=15, justify='r') + '|', blankline=.false.)
        ENDDO

        CALL update_log(num2char('<coda S>', justify='c', width=30) + num2char('Sheet', width=15, justify='r') + '|' + &
                        num2char('Min trvt', width=15, justify='r') + '|' + num2char('Max trvt', width=15, justify='r') + '|' + &
                        num2char('Envelopes', width=15, justify='r') + '|')

        DO sheet = 1, maxsheets(2)
          CALL update_log(num2char('', width=30) + num2char(sheet, width=15, justify='r')   + '|' +    &
                          num2char(coda%lotrvt(sheet, 2), width=15, justify='r', notation='f', precision=2) + '|' +    &
                          num2char(coda%uptrvt(sheet, 2), width=15, justify='r', notation='f', precision=2) + '|' +    &
                          num2char(coda%counts(sheet, 2), width=15, justify='r') + '|', blankline=.false.)
        ENDDO

      ENDIF


print*, 'a0 ', 1./a0
print*, 'b0 ', b0
print*, SUM(coda%counts)
!

      IF (ALLOCATED(coda%envelope)) DEALLOCATE(coda%envelope, coda%direct)

      ALLOCATE(coda%envelope(SIZE(time), SUM(coda%counts)), coda%direct(SUM(coda%counts)))

      ASSOCIATE(model => input%attenuation(input%receiver(rec)%attenuation))

        gsp = model%gps(band) / 6._r32

!       !$omp parallel do default(shared) private(wave, sheet, i, r, t, tp, ts, vs, wp, ws, d0) reduction(max: ok) collapse(3)
        DO wave = 1, 2
          DO sheet = 1, MAXVAL(maxsheets)
            DO i = 1, coda%counts(sheet, wave)

              t = coda%lotrvt(sheet, wave) + (i - 1) * DC
              r = (t - b0(sheet, wave)) / a0(sheet, wave)

              IF (wave .eq. 1) THEN
                tp = t
                ts = a0(sheet, 2) * r + b0(sheet, 2)
              ELSE
                ts = t
                tp = a0(sheet, 1) * r + b0(sheet, 1)
              ENDIF

              vs = r / ts

              IF (wave .eq. 1) THEN
                wp = 1._r32
                ws = 0._r32
              ELSE
                wp = 0._r32
                ws = 1._r32
              ENDIF

              c = coda%displs(sheet, wave) + i

              ! CALL rtt(time, tp, ts, model%gpp(band), model%gps(band), gsp, model%gss(band), model%b(band), vs, wp, ws,   &
              !           0.25_r32, d0, coda%envelope(:, c), ok)

              coda%direct(c) = d0(wave)

print*, i, sheet, wave, ' x ', c, ' -- ', tp, ts, ' - ', t, NINT( (t-coda%lotrvt(sheet, wave)) / DC ) + 1

            ENDDO
          ENDDO
        ENDDO
!        !$omp end parallel do

      END ASSOCIATE

print*, 'ok coda ', ok

    END SUBROUTINE lookup_coda

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *


END MODULE m_isochron
