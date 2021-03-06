MODULE m_isochron

  ! Purpose:
  !   this module is the core of the broadband tool and contains a collection of subprograms to carry out the isochron integration
  !   and related helper functions
  !
  !   solve_isochron_integral -> does the actual integration
  !   integration_step        -> set the integration step based on desired accuracy ("cuts")
  !   hilbert                 -> compute the Hilbert transform
  !
  !
  !   lookup_fs, setfs, rnmc  -> produce a lookup table of (complex-valued) P-SV free-surface amplification coefficients
  !
  !   rayshooting, raycorner  -> compute ray parameters (spreading, traveltime, etc.) for each triangle
  !
  !   perturbed_mechanism     -> return strike, dip and rake for each geometrically perturbed (rough) triangle
  !
  !   node2disk, write_srf    -> write rupture model to disk
  !
  !   correct4impz            -> correct for filter response at several frequency bands (e.g. 1-2, 2-4, 4-8, etc)
  !
  !
  ! Revisions:
  !     Date                    Description of change
  !     ====                    =====================
  !   08/03/21                  original version
  !

  USE                :: omp_lib
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

  PUBLIC :: solve_isochron_integral, node2disk, correct4impz, write_srf

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

#ifdef MPI
  INTEGER(i32), PARAMETER :: COMM = MPI_COMM_SELF
#else
  INTEGER(i32), PARAMETER :: COMM = 0
#endif

  REAL(r32), PARAMETER :: PI = 3.14159265358979323846_r64
  REAL(r32), PARAMETER :: RAD_TO_DEG = 180._r32 / PI
  REAL(r32), PARAMETER :: BIG = HUGE(0._r32)
  REAL(r32), PARAMETER :: DP = 0.0002_r32, DPSM = 0.007_r32            !< free-surface smoothing parameters
  REAL(r32), PARAMETER :: DCP = 1._r32, DCT = 0.5_r32                  !< resolution coda tableau (path in km, traveltime in s)
  REAL(r32), PARAMETER :: ISQRT3 = 1._r32 / SQRT(3._r32)
  REAL(r32), PARAMETER :: GSS_MIN = 0.0001_r32                         !< lowest non-zero scattering coefficient

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  INTEGER(i32),              DIMENSION(2) :: minsheets, maxsheets
  REAL(r32),    ALLOCATABLE, DIMENSION(:) :: shooting                 !< depth of shooting points

  TYPE :: co
    REAL(r32),              DIMENSION(2)     :: lotrvt, uptrvt, lopath, uppath
    REAL(r32), ALLOCATABLE, DIMENSION(:,:)   :: pdirect, sdirect
    REAL(r32), ALLOCATABLE, DIMENSION(:,:,:) :: penvelope, senvelope
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

      CHARACTER(:), ALLOCATABLE                             :: fo

      COMPLEX(r32)                                          :: fxpfzs, fzpfxs, g31, g21, g23, ga, gb
      COMPLEX(r32),              DIMENSION(3,3)             :: fp, fs, g
      COMPLEX(r32), ALLOCATABLE, DIMENSION(:)               :: fxp, fzp, tur, tui, tv

      INTEGER(i32)                                          :: i, j, ref, icr, totnutr, seed, src, sheet, wtp, it1, it2, ic, shift
      INTEGER(i32)                                          :: it, icut, model, npts, noslip, ipc, itc, lu

      INTEGER(i32),              DIMENSION(2)               :: weight
      INTEGER(i32),              DIMENSION(3)               :: iuc, ivc, shot

      REAL(r32)                                             :: m0, moment, strike, dip, urec, vrec, wrec, taumin, taumax, sd, cd
      REAL(r32)                                             :: srfvel, scale
      REAL(r32)                                             :: avepath, avetrvt, scattering, attenuation, a0, area

      REAL(r32)                                             :: dx, dy, dr, cpsi, spsi, sloloc, pn, signp, sthf, cthf
      REAL(r32)                                             :: rn, rv, ru, bn, cn, bu, cu, bv, cv, stho, ctho

      REAL(r32)                                             :: tau31, tau21, tau23, c, t, du, dv, p13, p12, p32, dl, dt

      REAL(r32),                 DIMENSION(2)               :: ncuts, tricross
      REAL(r32),                 DIMENSION(3)               :: u, v, w, x, y, z, slip, rupture, rake, nrl, rho, alpha, beta, mu
      REAL(r32),                 DIMENSION(3)               :: p, q, path, trvt, repi, tau, dtau, sr, cr, velloc

      REAL(r32),                 DIMENSION(3)               :: ir, itheta, iphi
      REAL(r32),                 DIMENSION(3,3)             :: rvec, thvec, phvec

      REAL(r32),    ALLOCATABLE, DIMENSION(:)               :: fsp, mrf, envelope, stack, time, cseis, ctri
      REAL(r32),    ALLOCATABLE, DIMENSION(:,:)             :: rtri, itri, rseis, iseis, seis, nseis

      REAL(r64)                                             :: tictoc
      REAL(r64),                 DIMENSION(6)               :: stopwatch
      REAL(r64),                 DIMENSION(10)              :: ostopwatch

      !-----------------------------------------------------------------------------------------------------------------------------

      ok = 0

      ! don't do anything if velocity model is not the right one for current receiver
      IF (input%receiver(rec)%velocity .ne. vel) RETURN

#ifdef PERF
      stopwatch(:)  = 0._r64
      ostopwatch(:) = 0._r64
#endif

      weight(:) = 1

      SELECT CASE(input%advanced%waves)
        CASE(1)
          weight(2) = 0                 !< neglect S-waves contribution
        CASE(2)
          weight(1) = 0                 !< neglect P-waves contribution
      END SELECT

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! ---------------------------------------- lookup table free-surface coefficients --------------------------------------------

      ASSOCIATE(model => input%attenuation(1), fmax => input%coda%fmax)

        IF (input%advanced%verbose .eq. 2) THEN
          CALL update_log(num2char('Integration for recv', width=30, fill='.') + num2char(rec, width=15, justify='r') + '|')

          CALL update_log(num2char('<freq band>', width=30, justify='c') +                                               &
                          num2char(num2char(model%lcut(band),            notation='f', width=6, precision=1) + ', ' +    &
                                   num2char(MIN(model%hcut(band), fmax), notation='f', width=6, precision=1),            &
                          width=15, justify='r') + '|')
        ENDIF

      END ASSOCIATE

      model = input%receiver(rec)%attenuation      !< attenuation model

      CALL setup_interpolation('linear', 'zero', ok)

      IF (ok .ne. 0) RETURN

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! ---------------------------------------- lookup table free-surface coefficients --------------------------------------------

      CALL lookup_fs(ok, rec, fsp, fxp, fzp)

      IF (ok .ne. 0) RETURN

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! ---------------------------------------- allocate arrays to handle seismograms ---------------------------------------------

      ALLOCATE(seis(SIZE(timeseries%sp%time), 3))

      seis(:, :) = 0._r32

#ifdef DEBUG
      ALLOCATE(cseis(SIZE(timeseries%sp%time)))
      cseis(:) = 0._r32
#endif

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! --------------------------------------- solve integral for each mesh refinement --------------------------------------------

      moment = 0._r32      !< keep track of actual moment, this will be used to rescale synthetics to target moment

      ! select function to compute rupture at grid nodes, depending whether we are dealing with extended- or point-sources
      nodefun => rik_at_nodes

      IF (input%source%is_point) nodefun => ptrsrc_at_nodes

      ! set "seed" such that random numbers depend on fault plane number and iteration
      seed = input%source%seed + (iter - 1) * SIZE(plane) + pl

      ! loop over mesh refinements
      DO ref = 1, SIZE(nvtr)

        IF (input%advanced%verbose .eq. 2) THEN
          CALL update_log(num2char('<mesh refinement>', width=30, justify='c') + num2char(ref, width=15, justify='r') + '|')
        ENDIF

        ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --
        ! ------------------------------------------- set integration time-step ----------------------------------------------------

        CALL integration_step(ok, ref, rec, pl, vel, iter, dt)

        IF (input%advanced%verbose .eq. 2) THEN
          CALL update_log(num2char('<time step>', width=30, justify='c') +    &
                          num2char(dt, notation='f', width=15, precision=4, justify='r') + '|')
        ENDIF

        npts = NINT(timeseries%sp%time(SIZE(timeseries%sp%time)) / dt) + 2

        ALLOCATE(time(npts), mrf(npts), tur(npts/2+1), tui(npts/2+1), tv(npts/2+1))
        ALLOCATE(rseis(npts, 3), iseis(npts, 3), rtri(npts, 3), itri(npts, 3))

#ifdef DEBUG
        ALLOCATE(ctri(npts))
#endif

        DO it = 1, npts
          time(it) = (it - 1) * dt
        ENDDO

        ! these stack contributions from each mesh refinement
        DO ic = 1, 3
          DO it = 1, npts
            rseis(it, ic) = 0._r32
            iseis(it, ic) = 0._r32
          ENDDO
        ENDDO

        ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --
        ! --------------------------------------------------- interpolate noise ----------------------------------------------------

        ALLOCATE(nseis(npts, 3))

        ! no need to filter before since integration "dt" is always smaller than noise (or output) "dt"
        DO ic = 1, 3
          CALL interpolate(timeseries%cd%time, timeseries%cd%xyz(:, ic, rec), time, nseis(:, ic))
        ENDDO

        noslip      = 0       !< total triangles having zero slip
        ncuts(:)    = 0       !< total isochron cuts for each wave type
        tricross(:) = 0       !< total triangles that were cut

        ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --
        ! ------------------------------------------ define RIK model on mesh nodes ------------------------------------------------

#ifdef PERF
        CALL watch_start(tictoc, COMM)
#endif

        ALLOCATE(nodes(nugr(ref), nvgr(ref)))

        CALL nodefun(ref, pl, vel, seed)                !< define slip, rupture time and rise time on mesh nodes

#ifdef PERF
        CALL watch_stop(tictoc, COMM)
        stopwatch(1) = stopwatch(1) + tictoc
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
        stopwatch(2) = stopwatch(2) + tictoc
#endif

#ifdef PERF
        CALL watch_start(tictoc, COMM)
#endif

        ! shoot rays between receivers and a set of sources spanning (vertically) current mesh
        CALL rayshooting(ok, ref, pl, vel)

#ifdef PERF
        CALL watch_stop(tictoc, COMM)
        stopwatch(3) = stopwatch(3) + tictoc
#endif

#ifdef PERF
        CALL watch_start(tictoc, COMM)
#endif

        ! fill lookup table for envelopes
        CALL lookup_coda(ok, ref, rec, band, pl, vel, time)

#ifdef PERF
        CALL watch_stop(tictoc, COMM)
        stopwatch(4) = stopwatch(4) + tictoc
#endif

#ifdef DEBUG
        ctri(:) = 0._r32
#endif

        area = dutr(ref) * dvtr(ref) * 0.5_r32

        ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --
        ! ------------------------------------------- allocate memory for FFT/IFFT -------------------------------------------------

        CALL make_fftw_plan([npts])

        !$omp parallel do default(shared) private(i, j, iuc, ivc, icr, slip, u, v, w, x, y, z, nrl, strike, dip, rake, sd, cd)   &
        !$omp private(sr, cr, tictoc, rupture, alpha, beta, rho, mu, ic, src, shot, mrf, urec, vrec, wrec, repi)     &
        !$omp private(rtri, itri, wtp, velloc, srfvel, scattering, sheet, p, path, trvt, q, tau, taumin, taumax, it1)   &
        !$omp private(it2, attenuation, it, tur, tui, tv, a0, ipc, itc, shift)   &
        !$omp private(fxpfzs, fzpfxs, fp, fs, dx, dy, dr, cpsi, spsi, sloloc, pn, signp, sthf, cthf, rn, rv, ru, bn, cn, bu, cu)  &
        !$omp private(bv, cv, stho, ctho, ir, itheta, iphi, rvec, thvec, phvec, g)     &
        !$omp private(icut, g31, g21, g23, ga, gb, tau31, tau21, tau23, dtau, c, t, du, dv, p13, p12, p32, dl)  &
#ifdef DEBUG
        !$omp reduction(+: ctri) &
#endif
        !$omp reduction(+: m0, rseis, iseis, ostopwatch, ncuts, tricross, noslip) collapse(2)
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

            ! this makes triangle degenerate (all "z" = MIN_DEPTH) if all "z" are above MIN_DEPTH
            ! DO icr = 1, 3
            !   z(icr) = MAX(MIN_DEPTH, z(icr))             !< make sure we never breach the free-surface (clip roughness)
            ! ENDDO

            IF (ANY(z .lt. MIN_DEPTH)) THEN
              w(:) = 0._r32                          !< simply set roughness to zero and recompute x, y, z
              CALL uvw2xyz(pl, u, v, w, x, y, z)
            ENDIF

            CALL normal2tri(x, y, z, nrl)

            ! always have normal pointing upward (i.e. negative z direction)
            IF (MOD(i + (j-1)*totnutr, 2) == 0) nrl(:) = -nrl(:)

            strike = plane(pl)%strike !* DEG_TO_RAD
            dip    = plane(pl)%dip    !* DEG_TO_RAD

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
            ostopwatch(1) = ostopwatch(1) + tictoc
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
            ostopwatch(2) = ostopwatch(2) + tictoc
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
            ostopwatch(3) = ostopwatch(3) + tictoc
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

                CALL raycorner(sheet, shooting, repi, z, shot, wtp, p, path, trvt, q)

#ifdef PERF
                CALL watch_stop(tictoc, COMM)
                ostopwatch(4) = ostopwatch(4) + tictoc
#endif

                ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
                ! ------------------------------------------ time integration limits -----------------------------------------------

                IF (ANY(trvt .eq. 0._r32)) CYCLE       !< all corners must belong to same sheet

                tau(:) = trvt(:) + rupture(:)          !< sum travel-time and rupture time

                taumax = MAXVAL(tau)
                taumin = MINVAL(tau)

                it1 = NINT(taumin / dt) + 1

                it2 = it1 + (taumax - taumin) / dt - 1      !< when (taumax-taumin) < dt, it2 < it1

                it2 = MIN(npts, it2)                        !< limit to max number of time points

                ! jump to next sheet if no isochron is spanned (i.e. no cuts)
                IF (it2 .lt. it1) CYCLE

                tricross(wtp) = tricross(wtp) + 1

                ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
                ! ------------------------------------------ evaluate integral kernel ----------------------------------------------

#ifdef PERF
                CALL watch_start(tictoc, COMM)
#endif

                DO icr = 1, 3

                  dx = urec - u(icr)
                  dy = vrec - v(icr)

                  dr = HYPOT(dx, dy)

                  ! cos and sin of angle psi
                  cpsi = dx / dr
                  spsi = dy / dr

                  sloloc = 1._r32 / velloc(icr)        !< slowness at corner

                  pn = p(icr)

                  ! "p" is "-p" for downgoing rays
                  IF (pn .gt. sloloc) pn = -(2._r32 * sloloc - pn)

                  ! get sign of ray (+ if upgoing, - if downgoing)
                  signp = SIGN(1._r32, pn)

                  pn = ABS(pn)         !< from here onwards we need only its absolute value

                  ! sine of angle eta at source
                  sthf = velloc(icr) * pn

#ifdef ERROR_TRAP
                  IF (sthf .gt. 1.01_r32) CALL report_error('solve_isochron_integral - ERROR: sin(eta) sensibly larger than 1')
#endif

                  sthf = MIN(sthf, 1._r32)     !< handle roundoff errors

                  ! cosine of angle eta at source
                  cthf = signp * SQRT(1._r32 - sthf**2)

                  ! dot products
                  ! the fault normal points in the -y direction for a 90 degree dip, and a positive u component of slip
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

                  stho = pn * srfvel

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

                  ! interpolate free-surface amplification coefficients at ray parameter "pn"
                  CALL interpolate(fsp, fxp, pn, fxpfzs)
                  CALL interpolate(fsp, fzp, pn, fzpfxs)

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

                  ! complete integral kernel
                  DO ic = 1, 3
                    g(ic, icr) = g(ic, icr) * slip(icr) * mu(icr) * q(icr)
                  ENDDO

                ENDDO

#ifdef PERF
                CALL watch_stop(tictoc, COMM)
                ostopwatch(5) = ostopwatch(5) + tictoc
#endif

                ! term for direct wave amplitude attenuation
                attenuation = EXP(-(mean(path) * scattering + input%attenuation(model)%b(band) * mean(trvt)) / 2._r32)

                ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
                ! ---------------------------------------------- time integration --------------------------------------------------

#ifdef PERF
                CALL watch_start(tictoc, COMM)
#endif

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
                    IF ( (ic .eq. 1) .and. (icut .ne. 0) ) ncuts(wtp) = ncuts(wtp) + 1

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
                    rtri(it, ic) = rtri(it, ic) + c * REAL(ga + gb, r32) * dl * 0.5_r32 * attenuation * weight(wtp)
                    itri(it, ic) = itri(it, ic) + c * AIMAG(ga + gb)     * dl * 0.5_r32 * attenuation * weight(wtp)

                  ENDDO      !< end loop over time (cuts)

                ENDDO      !< end loop over components

#ifdef PERF
                CALL watch_stop(tictoc, COMM)
                ostopwatch(6) = ostopwatch(6) + tictoc
#endif

                ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
                ! ------------------------------------------ add coda contribution -------------------------------------------------

                ipc = NINT( (mean(path) - coda%lopath(wtp)) / DCP ) + 1         !< index closest path value
                itc = NINT( (mean(trvt) - coda%lotrvt(wtp)) / DCT ) + 1         !< index closest traveltime value

                shift = NINT(mean(rupture) / dt) + 1 - 1         !< additional time-shift due to rupture (in samples)

                ! here below coda contributions due to direct P- and S-waves are treated separately. In any case, envelopes are
                ! scaled such that amplitude of scattered direct wave is given by an isotropic source (0.44 for P, 0.60 for S - see
                ! Boore&Boatwright 1984) at free-surface (~2, see book of Sato&Fehler). SQRT(3) is there because isotropic radiation
                ! is distributed evenly amongst the three directions of motion. Product "envelope * noise" is independent of the
                ! time-shift due to rupture.

!                 IF (wtp .eq. 1) THEN
!
!                   a0 = (2 * ABS(mean(q)) * mean(mu) * mean(slip)) * area * attenuation * 2._r32 * 0.33_r32 * ISQRT3 / 2.
!                   a0 = a0 / maxsheets(wtp)                            !< each sheet contributes equally
!                   a0 = a0 / coda%pdirect(itc, ipc) * weight(wtp)
!
!                   DO it = 1, npts - shift
!                     rtri(it + shift, 1) = rtri(it + shift, 1) + coda%penvelope(it, itc, ipc) * a0 * nseis(it, 1)
!                     rtri(it + shift, 2) = rtri(it + shift, 2) + coda%penvelope(it, itc, ipc) * a0 * nseis(it, 2)
!                     rtri(it + shift, 3) = rtri(it + shift, 3) + coda%penvelope(it, itc, ipc) * a0 * nseis(it, 3)
!                   ENDDO
!
! #ifdef DEBUG
!                   DO it = 1, npts - shift
!                     ctri(it + shift) = ctri(it + shift) + coda%penvelope(it, itc, ipc) * a0
!                   ENDDO
! #endif
!
!                 ELSE      !< contribution of both SH- and SV-waves

                IF (wtp .eq. 2) THEN

                  a0 = (ABS(mean(q)) * mean(mu) * mean(slip)) * area * attenuation * 2._r32 * 0.6_r32 * ISQRT3 / 2._r32
                  a0 = a0 / maxsheets(wtp)                            !< each sheet contributes equally
                  a0 = a0 / coda%sdirect(itc, ipc) * weight(wtp)

                  DO it = 1, npts - shift
                    rtri(it + shift, 1) = rtri(it + shift, 1) + coda%senvelope(it, itc, ipc) * a0 * nseis(it, 1)
                    rtri(it + shift, 2) = rtri(it + shift, 2) + coda%senvelope(it, itc, ipc) * a0 * nseis(it, 2)
                    rtri(it + shift, 3) = rtri(it + shift, 3) + coda%senvelope(it, itc, ipc) * a0 * nseis(it, 3)
                  ENDDO

#ifdef DEBUG
                  DO it = 1, npts - shift
                    ctri(it + shift) = ctri(it + shift) + coda%senvelope(it, itc, ipc) * a0
                  ENDDO
#endif

                ENDIF

              ENDDO     !< end loop over sheets

            ENDDO    !< end loop over wave types

            ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
            ! --------------------------------------------- convolution with MRF ---------------------------------------------------

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
            ostopwatch(9) = ostopwatch(9) + tictoc
#endif

            ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
            ! -------------------------------------------- rotate & stack ----------------------------------------------------------

            ! move back from "u/t" (fault-parallel/fault-normal) to "x/y" coordinates
            CALL rotate(rtri(:, 1), rtri(:, 2), rtri(:, 3), 0._r32, 0._r32, strike)
            CALL rotate(itri(:, 1), itri(:, 2), itri(:, 3), 0._r32, 0._r32, strike)

            DO ic = 1, 3
              DO it = 1, npts
                rseis(it, ic) = rseis(it, ic) + rtri(it, ic)
                iseis(it, ic) = iseis(it, ic) + itri(it, ic)
              ENDDO
            ENDDO

            m0 = m0 + mean(mu) * mean(slip)      !< average moment contribution (area is included later)

          ENDDO
        ENDDO
        !$omp end parallel do

        CALL dealloc_nodes()

        moment = moment + m0 * area      !< sum is over mesh refinements

        CALL destroy_fftw_plan([npts])

        ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --
        ! ------------------------------------ Hilbert transform, filter, differentiate & stack ------------------------------------

#ifdef PERF
        CALL watch_start(tictoc, COMM)
#endif

        ASSOCIATE(model => input%attenuation(1), fmax => input%coda%fmax)

          IF (band .eq. 1) THEN
            CALL make_iir_plan(ok, 'butter', dt, [0.25_r32, MIN(model%hcut(band), fmax)], 'pass', 2, zphase = .true.)
          ELSE
            CALL make_iir_plan(ok, 'butter', dt, [model%lcut(band), MIN(model%hcut(band), fmax)], 'pass', 2, zphase = .true.)
          ENDIF

          IF (ok .ne. 0) THEN
            CALL report_error(filter_error(ok))
            RETURN
          ENDIF

          ALLOCATE(stack(SIZE(timeseries%sp%time)))

          DO ic = 1, 3

            CALL hilbert(iseis(:, ic))

            DO it = 1, npts
              rseis(it, ic) = (rseis(it, ic) + iseis(it, ic)) * dt               !< "dt" is from convolution
            ENDDO

            rseis(:, ic) = iir(rseis(:, ic), ok)                                 !< filter

            IF (ok .ne. 0) THEN
              CALL report_error(filter_error(ok))
              RETURN
            ENDIF

            rseis(:, ic) = differentiate(rseis(:, ic), dt)                       !< move from displacement to velocity

            CALL interpolate(time, rseis(:, ic), timeseries%sp%time, stack)      !< resampling (most likely downsampling)

            DO it = 1, SIZE(timeseries%sp%time)
              seis(it, ic) = seis(it, ic) + stack(it)
            ENDDO

          ENDDO

          CALL destroy_iir_plan()

        END ASSOCIATE

#ifdef PERF
        CALL watch_stop(tictoc, COMM)
        stopwatch(5) = stopwatch(5) + tictoc
#endif

        IF (input%advanced%verbose .eq. 2) THEN

          CALL update_log(num2char('<cuts>', justify='c', width=30) +           &
                          num2char('P', width=15, justify='r') + '|' +  num2char('S', width=15, justify='r') + '|' +  &
                          num2char('Skipped', width=15, justify='r') + '|')

          CALL update_log(num2char('', width=30)  +  &
                          num2char(ncuts(1)/tricross(1), notation='f', width=15, precision=1, justify='r') + '|'  +   &
                          num2char(ncuts(2)/tricross(2), notation='f', width=15, precision=1, justify='r') + '|'  +   &
                          num2char(num2char(noslip) + ' (' +  &
                                   num2char(100*noslip/(totnutr*nvtr(ref))) + '%)', width=15, justify='r') + '|', blankline=.false.)

        ENDIF

#ifdef DEBUG

        CALL interpolate(time, ctri(:), timeseries%sp%time, stack)      !< resampling

        DO it = 1, SIZE(timeseries%sp%time)
          cseis(it) = cseis(it) + stack(it)
        ENDDO

        DEALLOCATE(ctri)

#endif

        DEALLOCATE(time, mrf, tur, tui, tv, stack)
        DEALLOCATE(rseis, iseis, rtri, itri, nseis)

      ENDDO        !< end loop over mesh refinements

      scale = plane(pl)%targetm0 / moment          !< scaling factor to scale to desired moment
      scale = scale * 1.E-15_r32                   !< take into account the fact that units for "q" were km, g/cm^3, km/s

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! ---------------------------------------------- scale to desired moment -----------------------------------------------------

#ifdef PERF
      CALL watch_start(tictoc, COMM)
#endif

      DO ic = 1, 3
        timeseries%sp%xyz(:, ic, rec) = timeseries%sp%xyz(:, ic, rec) + seis(:, ic) * scale
      ENDDO

#ifdef PERF
      CALL watch_stop(tictoc, COMM)
      stopwatch(6) = stopwatch(6) + tictoc
#endif


#ifdef PERF
      print*, 'Time RIK nodes: ', real(stopwatch(1))
      print*, 'Time roughness: ', real(stopwatch(2))
      print*, 'Time shooting: ', real(stopwatch(3))
      print*, 'Time coda: ', real(stopwatch(4))

      print*, 'Time mechanism: ', real(ostopwatch(1))
      print*, 'Time corners: ', real(ostopwatch(2))
      print*, 'Time MRF + FFT: ', real(ostopwatch(3))
      print*, 'Time ray corners: ', real(ostopwatch(4))
      print*, 'Time kernel: ', real(ostopwatch(5))
      ! print*, 'Time envelope: ', real(stopwatch(9))
      print*, 'Time intg: ', real(ostopwatch(6))
      print*, 'Time iir: ', real(ostopwatch(8))
      print*, 'Time conv MRF: ', real(ostopwatch(9))

      print*, 'Time Hilbert: ', real(stopwatch(5))
      print*, 'Time Filter: ', real(stopwatch(6))
#endif

#ifdef DEBUG

      fo = 'envelope_rec' + num2char(rec) + '_band' + num2char(band) + '_plane' + num2char(pl) + '_vel' + num2char(vel) + '.txt'

      OPEN(newunit = lu, file = fo, status = 'replace', form = 'formatted', access = 'sequential', action = 'write', IOSTAT = ok)

      DO it = 1, SIZE(timeseries%sp%time)
        WRITE(lu, *) timeseries%sp%time(it), cseis(it) * scale
      ENDDO

      CLOSE(lu)

#endif

    END SUBROUTINE solve_isochron_integral

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE raycorner(sheet, shooting, repi, z, shot, wtp, p, path, trvt, q)

      ! Purpose:
      !   to compute ray parameter "p", travelled distance "path", traveltime "trvt" and spreading factor "q" for a corner at depth
      !   "z" and referred to wave type "wtp" and sheet "sheet", by interpolating values pre-computed at discrete depth levels.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   08/03/21                  original version
      !

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

    END SUBROUTINE raycorner

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

      ! Purpose:
      !   to compute a set of complex-valued free-surface amplification coefficients for P- and SV-waves in the radial ("fxp") and
      !   vertical ("fzp") direction for a set of ray parameters "fsp". These coefficients are smoothed with a window determined by
      !   "irort".
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   08/03/21                  original version
      !

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

      ! Purpose:
      !   to smooth complex free-surface amplification coefficients "a" based on a window controlled by "irort". Triangular windows
      !   (i.e. "irort = 1") seems to return best results.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   08/03/21                  original version
      !

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

    SUBROUTINE node2disk(ok, band, pl, vel, iter)

      ! Purpose:
      !   to write to disk the rupture parameters defined for frequency babd "band", plane "pl" embedded in velocity model "vel" and
      !   iteration "iter".
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   08/03/21                  original version
      !

      INTEGER(i32),                          INTENT(OUT) :: ok
      INTEGER(i32),                          INTENT(IN)  :: band, pl, vel, iter
      CHARACTER(:), ALLOCATABLE                          :: fo
      INTEGER(i32)                                       :: i, j, ref, lu, icr, totnutr, seed
      INTEGER(i32),             DIMENSION(3)             :: iuc, ivc
      REAL(r32)                                          :: m0, moment, strike, dip, rho, beta
      REAL(r32),                DIMENSION(3)             :: u, v, w, x, y, z, slip, rise, rupture, rake, nrl, mu

      !-----------------------------------------------------------------------------------------------------------------------------

      CALL setup_interpolation('linear', 'zero', ok)

      IF (ok .ne. 0) RETURN

      fo = 'node_band' + num2char(band) + '_pl' + num2char(pl) + '_vel' + num2char(vel) + '_iter' + num2char(iter) + '.bin'

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

            ! this makes triangle degenerate (all "z" = MIN_DEPTH) if all "z" are above MIN_DEPTH
            ! DO icr = 1, 3
            !   z(icr) = MAX(MIN_DEPTH, z(icr))             !< make sure we never breach the free-surface (clip roughness)
            ! ENDDO

            IF (ANY(z .lt. MIN_DEPTH)) THEN
              w(:) = 0._r32                          !< simply set roughness to zero and recompute x, y, z
              CALL uvw2xyz(pl, u, v, w, x, y, z)
            ENDIF

            !$omp ordered
            WRITE(lu) x, y, z, slip, rise, rupture         !< values at each corner
            !$omp end ordered

            CALL normal2tri(x, y, z, nrl)

            ! always have normal pointing upward (i.e. negative)
            IF (MOD(i + (j-1)*totnutr, 2) == 0) nrl(:) = -nrl(:)

            strike = plane(pl)%strike !* DEG_TO_RAD
            dip    = plane(pl)%dip    !* DEG_TO_RAD

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

            m0 = m0 + mean(mu) * mean(slip)

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

    SUBROUTINE write_srf(ok, pl, vel, iter)

      ! Purpose:
      !   to write finite-fault rupture properties in a SRF file v2.0 (see http://scec.usc.edu/scecpedia/Standard_Rupture_Format for
      !   a description).
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   08/03/21                  original version
      !

      INTEGER(i32),                           INTENT(OUT) :: ok
      INTEGER(i32),                           INTENT(IN)  :: pl, vel, iter
      CHARACTER(:), ALLOCATABLE                           :: fo
      INTEGER(i32)                                        :: lu, seed, npts, it, ref, totnutr, i, j, icr, offset, nt
      INTEGER(i32),              DIMENSION(3)             :: iuc, ivc
      REAL(r32)                                           :: m0, dt, area, strike, dip
      REAL(r32),                 DIMENSION(3)             :: slip, u, v, w, x, y, z, nrl, rake, rupture, beta, rho, lon, lat
      REAL(r32),    ALLOCATABLE, DIMENSION(:)             :: time, mrf

#ifdef DEBUG
      CHARACTER(3)                                        :: version_
      INTEGER(i32)                                        :: offset_, totnutr_, nt_, nt2_, nt3_
      REAL(r32)                                           :: strike_, dip_, area_, dt_, beta_, lon_, lat_, z_, rupture_, rho_
      REAL(r32)                                           :: rake_, slip_, slip2_, slip3_
      REAL(r32),    ALLOCATABLE, DIMENSION(:)             :: mrf_
#endif

      !-----------------------------------------------------------------------------------------------------------------------------

      fo = 'srf_' + 'vel_' + num2char(vel) + '_' + 'iter_' + num2char(iter) + '.bin'

      OPEN(newunit = lu, file = fo, status = 'replace', form = 'unformatted', access = 'stream', action = 'write', IOSTAT = ok)

      IF (ok .ne. 0) THEN
        CALL report_error('Error while opening file' + TRIM(fo))
        RETURN
      ENDIF

      IF (pl .eq. 1) THEN
        WRITE(lu, POS=1) '2.0'                           !< 3 bytes
        ! WRITE(lu)        'PLANE', SIZE(plane)          !< 5 + 4 bytes
      ENDIF

      ! WRITE(lu, POS=3 + 9 + 44*(pl - 1) + 1) lon, lat, nstk, ndip, len, width, stk, dip, dtop, shyp, dhyp       !< 11 * 4 bytes

      INQUIRE(lu, POS=offset)        !< get position of where to write next item

#ifdef DEBUG
      offset_ = offset
#endif

      WRITE(lu, POS=offset) 'POINTS', SUM(nvtr * (2*nutr - 1))         !< total subfaults for current plane (all refinements)

      m0 = 0._r32

      ! select function to compute rupture at grid nodes, depending whether we are dealing with extended- or point-sources
      nodefun => rik_at_nodes

      IF (input%source%is_point) nodefun => ptrsrc_at_nodes

      ! set "seed" such that random numbers depend on fault plane number and iteration
      seed = input%source%seed + (iter - 1) * SIZE(plane) + pl

      ! loop over mesh refinements
      DO ref = 1, SIZE(nvtr)

        dt = timeseries%sp%dt

        npts = SIZE(timeseries%sp%time)

        ALLOCATE(time(npts), mrf(npts))

        DO it = 1, npts
          time(it) = (it - 1) * dt
        ENDDO

        ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --
        ! ------------------------------------------ define RIK model on mesh nodes ------------------------------------------------

        ALLOCATE(nodes(nugr(ref), nvgr(ref)))

        CALL nodefun(ref, pl, vel, seed)                !< define slip, rupture time and rise time on mesh nodes

        totnutr = 2 * nutr(ref) - 1                     !< total triangles in a row

        ! generate fault plane roughness
        CALL fault_roughness(ok, ref, pl, iter)

        area = dutr(ref) * dvtr(ref) * 0.5_r32 * 1.E+04_r32          !< triangle area in cm^2

        INQUIRE(lu, POS=offset)

!        !$omp parallel do ordered default(shared) private(i, j, iuc, ivc, slip, u, v, w, x, y, z, nrl, strike, dip, rake, icr)   &
!        !$omp private(rupture, beta, rho, mrf, lon, lat, nt) reduction(+: m0) collapse(2)
        DO j = 1, nvtr(ref)
          DO i = 1, totnutr

            ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
            ! ---------------------------------------- check corners have non-zero slip --------------------------------------------

            CALL cornr(j, i, iuc, ivc)        !< corner indices for current triangle

            DO icr = 1, 3
              slip(icr) = SUM(nodes(iuc(icr), ivc(icr))%slip) * 100._r32          !< slip in cm
            ENDDO

            CALL cornr2uv(iuc, ivc, ref, u, v)        !< on-fault coordinates

            DO icr = 1, 3
              w(icr) = roughness(iuc(icr), ivc(icr))      !< add roughness (always zero for point-sources)
            ENDDO

            CALL uvw2xyz(pl, u, v, w, x, y, z)            !< get cartesian coordinates (m)

            ! this makes triangle degenerate (all "z" = MIN_DEPTH) if all "z" are above MIN_DEPTH
            ! DO icr = 1, 3
            !   z(icr) = MAX(MIN_DEPTH, z(icr))             !< make sure we never breach the free-surface (clip roughness)
            ! ENDDO

            IF (ANY(z .lt. MIN_DEPTH)) THEN
              w(:) = 0._r32                          !< simply set roughness to zero and recompute x, y, z
              CALL uvw2xyz(pl, u, v, w, x, y, z)
            ENDIF

            CALL normal2tri(x, y, z, nrl)

            ! always have normal pointing upward (i.e. negative z direction)
            IF (MOD(i + (j-1)*totnutr, 2) == 0) nrl(:) = -nrl(:)

            strike = plane(pl)%strike !* DEG_TO_RAD
            dip    = plane(pl)%dip    !* DEG_TO_RAD

            CALL interpolate(plane(pl)%u, plane(pl)%v, plane(pl)%rake, u, v, rake)

            CALL perturbed_mechanism(strike, dip, rake, nrl)            !< get new strike, dip, rake (all values in radians)

            strike = strike * RAD_TO_DEG
            dip    = dip * RAD_TO_DEG

            ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
            ! ----------------------------------------------- corner parameters ----------------------------------------------------

            DO icr = 1, 3
              rupture(icr) = MINVAL(nodes(iuc(icr), ivc(icr))%rupture)
              beta(icr)    = vinterp(input%velocity(vel)%depth, input%velocity(vel)%vs, input%velocity(vel)%vsgrad, z(icr))
              rho(icr)     = vinterp(input%velocity(vel)%depth, input%velocity(vel)%rho, input%velocity(vel)%rhograd, z(icr))
              rake(icr)    = rake(icr) * RAD_TO_DEG
            ENDDO

            ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
            ! ---------------------------------------- moment rate function of triangle --------------------------------------------

            nt = npts

            IF (input%source%is_point) THEN
              mrf = dbrune(time, input%source%freq)
            ELSE
              CALL mrf_rik(iuc, ivc, dt, rupture, mrf, nt)     !< moment rate function (averaged over corners)
            ENDIF

            CALL normalize(mrf, dt)

            DO icr = 1, 3
              CALL utm2geo(y(icr), x(icr), input%origin%lon, input%origin%lat, lon(icr), lat(icr))     !< from UTM to geo coordinates
            ENDDO

            ! all values are referred to triangle (subfault) center. Note that only one thread at a time will first find the next
            ! position in the file, then append a first point block and then a second one. The following thread will be asked to
            ! find the next position in the file and then write. This implies that subfaults will be written in any order.
!            !$omp ordered
            INQUIRE(lu, POS=offset)

            WRITE(lu, POS=offset) mean(lon), mean(lat), mean(z)*1.E-03_r32, strike, dip, area, mean(rupture), dt,    &
                                  mean(beta)*1.E+02_r32, mean(rho)*1.E-02_r32, mean(rake), mean(slip)*1.E+02_r32,    &
                                  nt, 0._r32, 0, 0._r32, 0

            DO it = 1, nt
              WRITE(lu) mrf(it) * 1.E+02_r32             !< append MRF values
            ENDDO
!            !$omp end ordered

          ENDDO
        ENDDO
!        !$omp end parallel do

        CALL dealloc_nodes()

        DEALLOCATE(time, mrf)

      ENDDO        !< end loop over mesh refinements

      CLOSE(lu, iostat = ok)

      IF (ok .ne. 0) THEN
        CALL report_error('Error while closing file ' + fo)
        RETURN
      ENDIF


#ifdef DEBUG

      OPEN(newunit = lu, file = fo, status = 'old', form = 'unformatted', access = 'stream', action = 'read', IOSTAT = ok)

      IF (ok .ne. 0) THEN
        CALL report_error('Error while opening file' + TRIM(fo))
        RETURN
      ENDIF

      IF (pl .eq. 1) THEN

        READ(lu, POS=1) version_                         !< 3 bytes

        IF (version_ .ne. '2.0') THEN
          CALL report_error('"version" mismatch: ' + version_)
          ok = 1
          RETURN
        ENDIF

      ENDIF

      READ(lu, POS=offset_ + 6) totnutr_         !< total subfaults for current plane (all refinements)
      totnutr = SUM(nvtr * (2*nutr - 1))

      IF (assert([totnutr_], totnutr)) THEN
        CALL report_error('total number of subfaults ("points") mismatch: ' + num2char(totnutr_) + ', ' + num2char(totnutr))
        ok = 1
        RETURN
      ENDIF

      ! select function to compute rupture at grid nodes, depending whether we are dealing with extended- or point-sources
      ! nodefun => rik_at_nodes
      !
      ! IF (input%source%is_point) nodefun => ptrsrc_at_nodes

      ! set "seed" such that random numbers depend on fault plane number and iteration
      seed = input%source%seed + (iter - 1) * SIZE(plane) + pl

      ! loop over mesh refinements
      DO ref = 1, SIZE(nvtr)

        dt = timeseries%sp%dt

        npts = SIZE(timeseries%sp%time)

        ALLOCATE(time(npts), mrf(npts), mrf_(npts))

        DO it = 1, npts
          time(it) = (it - 1) * dt
        ENDDO

        ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --
        ! ------------------------------------------ define RIK model on mesh nodes ------------------------------------------------

        ALLOCATE(nodes(nugr(ref), nvgr(ref)))

        CALL nodefun(ref, pl, vel, seed)                !< define slip, rupture time and rise time on mesh nodes

        totnutr = 2 * nutr(ref) - 1                     !< total triangles in a row

        ! generate fault plane roughness
        CALL fault_roughness(ok, ref, pl, iter)

        area = dutr(ref) * dvtr(ref) * 0.5_r32 * 1.E+04_r32          !< triangle area in cm^2

        DO j = 1, nvtr(ref)
          DO i = 1, totnutr

            ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
            ! ---------------------------------------- check corners have non-zero slip --------------------------------------------

            CALL cornr(j, i, iuc, ivc)        !< corner indices for current triangle

            DO icr = 1, 3
              slip(icr) = SUM(nodes(iuc(icr), ivc(icr))%slip) * 100._r32          !< slip in cm
            ENDDO

            CALL cornr2uv(iuc, ivc, ref, u, v)        !< on-fault coordinates

            DO icr = 1, 3
              w(icr) = roughness(iuc(icr), ivc(icr))      !< add roughness (always zero for point-sources)
            ENDDO

            CALL uvw2xyz(pl, u, v, w, x, y, z)            !< get cartesian coordinates (m)

            ! this makes triangle degenerate (all "z" = MIN_DEPTH) if all "z" are above MIN_DEPTH
            ! DO icr = 1, 3
            !   z(icr) = MAX(MIN_DEPTH, z(icr))             !< make sure we never breach the free-surface (clip roughness)
            ! ENDDO

            IF (ANY(z .lt. MIN_DEPTH)) THEN
              w(:) = 0._r32                          !< simply set roughness to zero and recompute x, y, z
              CALL uvw2xyz(pl, u, v, w, x, y, z)
            ENDIF

            CALL normal2tri(x, y, z, nrl)

            ! always have normal pointing upward (i.e. negative z direction)
            IF (MOD(i + (j-1)*totnutr, 2) == 0) nrl(:) = -nrl(:)

            strike = plane(pl)%strike !* DEG_TO_RAD
            dip    = plane(pl)%dip    !* DEG_TO_RAD

            CALL interpolate(plane(pl)%u, plane(pl)%v, plane(pl)%rake, u, v, rake)

            CALL perturbed_mechanism(strike, dip, rake, nrl)            !< get new strike, dip, rake (all values in radians)

            strike = strike * RAD_TO_DEG
            dip    = dip * RAD_TO_DEG

            ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
            ! ----------------------------------------------- corner parameters ----------------------------------------------------

            DO icr = 1, 3
              rupture(icr) = MINVAL(nodes(iuc(icr), ivc(icr))%rupture)
              beta(icr)    = vinterp(input%velocity(vel)%depth, input%velocity(vel)%vs, input%velocity(vel)%vsgrad, z(icr))
              rho(icr)     = vinterp(input%velocity(vel)%depth, input%velocity(vel)%rho, input%velocity(vel)%rhograd, z(icr))
              rake(icr)    = rake(icr) * RAD_TO_DEG
            ENDDO

            ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
            ! ---------------------------------------- moment rate function of triangle --------------------------------------------

            nt = npts

            IF (input%source%is_point) THEN
              mrf = dbrune(time, input%source%freq)
            ELSE
              CALL mrf_rik(iuc, ivc, dt, rupture, mrf, nt)     !< moment rate function (averaged over corners)
            ENDIF

            CALL normalize(mrf, dt)

            DO icr = 1, 3
              CALL utm2geo(y(icr), x(icr), input%origin%lon, input%origin%lat, lon(icr), lat(icr))     !< from UTM to geo coordinates
            ENDDO

            READ(lu) lon_, lat_, z_, strike_, dip_, area_, rupture_, dt_, beta_, rho_, rake_, slip_, nt_, slip2_, nt2_, slip3_, nt3_

            IF (assert([lon_], mean(lon))) THEN
              CALL report_error('longitude mismatch: ' + num2char(lon_) + ', ' + num2char(mean(lon)))
              ok = 1
              RETURN
            ENDIF

            IF (assert([lat_], mean(lat))) THEN
              CALL report_error('latitude mismatch: ' + num2char(lat_) + ', ' + num2char(mean(lat)))
              ok = 1
              RETURN
            ENDIF

            IF (assert([z_], mean(z)*1.E-03_r32)) THEN
              CALL report_error('depth mismatch: ' + num2char(z_) + ', ' + num2char(mean(z)*1.E-03_r32))
              ok = 1
              RETURN
            ENDIF

            IF (assert([strike_], strike)) THEN
              CALL report_error('strike mismatch: ' + num2char(strike_) + ', ' + num2char(strike))
              ok = 1
              RETURN
            ENDIF

            IF (assert([dip_], dip)) THEN
              CALL report_error('dip mismatch: ' + num2char(dip_) + ', ' + num2char(dip))
              ok = 1
              RETURN
            ENDIF

            IF (assert([area_], area)) THEN
              CALL report_error('area mismatch: ' + num2char(area_) + ', ' + num2char(area))
              ok = 1
              RETURN
            ENDIF

            IF (assert([rupture_], mean(rupture))) THEN
              CALL report_error('rupture mismatch: ' + num2char(rupture_) + ', ' + num2char(mean(rupture)))
              ok = 1
              RETURN
            ENDIF

            IF (assert([dt_], dt)) THEN
              CALL report_error('dt mismatch: ' + num2char(dt_) + ', ' + num2char(dt))
              ok = 1
              RETURN
            ENDIF

            IF (assert([beta_], mean(beta)*1.E+02_r32)) THEN
              CALL report_error('vs mismatch: ' + num2char(beta_) + ', ' + num2char(mean(beta)*1.E+02_r32))
              ok = 1
              RETURN
            ENDIF

            IF (assert([rho_], mean(rho)*1.E-02_r32)) THEN
              CALL report_error('dens mismatch: ' + num2char(rho_) + ', ' + num2char(mean(rho)*1.E-02_r32))
              ok = 1
              RETURN
            ENDIF

            IF (assert([rake_], mean(rake))) THEN
              CALL report_error('rake mismatch: ' + num2char(rake_) + ', ' + num2char(mean(rake)))
              ok = 1
              RETURN
            ENDIF

            IF (assert([slip_], mean(slip)*1.E+02_r32)) THEN
              CALL report_error('slip mismatch: ' + num2char(slip_) + ', ' + num2char(mean(slip)*1.E+02_r32))
              ok = 1
              RETURN
            ENDIF

            IF (assert([nt_], nt)) THEN
              CALL report_error('mrf samples ("nt") mismatch: ' + num2char(nt_) + ', ' + num2char(nt))
              ok = 1
              RETURN
            ENDIF

            IF (assert([slip2_], 0._r32)) THEN
              CALL report_error('slip2 mismatch: ' + num2char(slip2_) + ', 0')
              ok = 1
              RETURN
            ENDIF

            IF (assert([nt2_], 0)) THEN
              CALL report_error('mrf samples ("nt2") mismatch: ' + num2char(nt2_) + ', 0')
              ok = 1
              RETURN
            ENDIF

            IF (assert([slip3_], 0._r32)) THEN
              CALL report_error('slip3 mismatch: ' + num2char(slip3_) + ', 0')
              ok = 1
              RETURN
            ENDIF

            IF (assert([nt3_], 0)) THEN
              CALL report_error('mrf samples ("nt3") mismatch: ' + num2char(nt3_) + ', 0')
              ok = 1
              RETURN
            ENDIF

            DO it = 1, nt_
              READ(lu) mrf_(it)
              mrf_(it) = mrf_(it) * 1.E+02_r32
              IF (assert([mrf_(it)], mrf(it))) THEN
                CALL report_error('mrf mismatch: ' + num2char(mrf_(it)) + ', ' + num2char(mrf(it)))
                ok = 1
                RETURN
              ENDIF
            ENDDO

          ENDDO
        ENDDO

        CALL dealloc_nodes()

        DEALLOCATE(time, mrf, mrf_)

      ENDDO        !< end loop over mesh refinements

      CLOSE(lu, iostat = ok)

      IF (ok .ne. 0) THEN
        CALL report_error('Error while closing file ' + fo)
        RETURN
      ENDIF

#endif

    END SUBROUTINE write_srf

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

      delta = MAXVAL(ABS(roughness)) * COS(plane(pl)%dip)! * DEG_TO_RAD)

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
                        num2char(rmax*1.E-03_r32, width=15, notation='f', precision=3, justify='r') + '|' + &
                        num2char(MINVAL(shooting)*1.E-03_r32, width=15, notation='f', precision=3, justify='r') + '|' + &
                        num2char(MAXVAL(shooting)*1.E-03_r32, width=15, notation='f', precision=3, justify='r') + '|' +  &
                        ! num2char(num2char(zmin, notation='f', width=6, precision=2) + ', ' +   &
                        !          num2char(zmax, notation='f', width=6, precision=2), width=15, justify='r') + '|',blankline=.false.)
                        num2char(mean([zmin, zmax])*1.E-03_r32, notation='f', width=15, precision=3, justify='r') + '|',    &
                        blankline=.false.)

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

      ! Purpose:
      !   to correct the amplitude of synthetic timeseries by removing the total impulse response of multiple bandpass filters used
      !   to compute synthetic seismograms with frequency-dependent scattering properties.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   08/03/21                  original version
      !

      USE, NON_INTRINSIC :: m_fft_real

      INTEGER(i32),                           INTENT(OUT) :: ok
      INTEGER(i32),                           INTENT(IN)  :: iter
      CHARACTER(:), ALLOCATABLE                           :: fo
      COMPLEX(r32), ALLOCATABLE, DIMENSION(:)             :: spectrum
      INTEGER(i32)                                        :: band, rec, i, lu, npts, ic
      REAL(r32),    ALLOCATABLE, DIMENSION(:)             :: freq, amp, phase, respz, blob

      !-----------------------------------------------------------------------------------------------------------------------------

      ASSOCIATE(model => input%attenuation(1), fmax => input%coda%fmax, dt => timeseries%sp%dt)

        ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --
        ! ----------------------------------------- overall filter response --------------------------------------------------------

        DO band = 1, SIZE(model%lcut)

          IF (model%lcut(band) .ge. fmax) EXIT

          IF (band .eq. 1) THEN
            CALL make_iir_plan(ok, 'butter', dt, [0.25_r32, MIN(model%hcut(band), fmax)], 'pass', 2, zphase = .true.)
          ELSE
            CALL make_iir_plan(ok, 'butter', dt, [model%lcut(band), MIN(model%hcut(band), fmax)], 'pass', 2, zphase = .true.)
          ENDIF

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
          CALL make_iir_plan(ok, 'butter', dt, [0.25_r32, MIN(model%hcut(band), fmax)], 'pass', 2, zphase = .true.)

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
            DO ic = 1, 3
              CALL fft(timeseries%sp%xyz(:, ic, rec), spectrum)
              spectrum(:) = spectrum(:) * respz(:)
              CALL ifft(timeseries%sp%xyz(:, ic, rec), spectrum)
            ENDDO
          ENDDO

          CALL destroy_fftw_plan([npts])

        ENDIF

      END ASSOCIATE


    END SUBROUTINE correct4impz

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE lookup_coda(ok, ref, rec, band, pl, vel, time)

      ! Purpose:
      !   to define a lookup table of coda envelopes for mesh refinement "ref", receiver "rec", frequency band "band", fault plane
      !   "pl" and velocity model "vel". Envelopes are computed at time values given by "time" and stored in module variable "coda".
      !   This subroutine computes also the amplitude of direct wave(s) for scaling purposes in "solve_isochron_integral".
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   08/03/21                  original version
      !   15/03/21                  compute P- and S-coda at once
      !

      INTEGER(i32),                INTENT(OUT) :: ok
      INTEGER(i32),                INTENT(IN)  :: ref, rec, band, vel, pl
      REAL(r32),    DIMENSION(:),  INTENT(IN)  :: time
      INTEGER(i32)                             :: totnutr, i, j, icr, src, wave, sheet, n, np, nt
      INTEGER(i32), DIMENSION(3)               :: iuc, ivc, shot
      REAL(r32)                                :: r, tp, ts, vs, wp, ws, gsp, vratio, pmin, pmax, tmin, tmax
      REAL(r32),    DIMENSION(2)               :: d0
      REAL(r32),    DIMENSION(3)               :: u, v, w, x, y, z, repi
      REAL(r32),    DIMENSION(2,2)             :: minbounds, maxbounds
      REAL(r32),    DIMENSION(3,2)             :: p, path, trvt, q

      !-----------------------------------------------------------------------------------------------------------------------------

      ok = 0

      totnutr = 2 * nutr(ref) - 1                     !< total triangles in a row

      minbounds(:,:) = BIG
      maxbounds(:,:) = -BIG

      vratio = 0._r32
      n      = 0

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! ---------------------------------- min/max path and traveltime over all triangles ------------------------------------------

      !$omp parallel do default(shared) private(i, j, iuc, ivc, u, v, w, icr, src, x, y, z, shot, repi, sheet, wave, p, path)   &
      !$omp private(trvt, q) reduction(max: maxbounds) reduction(min: minbounds) reduction(+: vratio, n)
      DO j = 1, nvtr(ref)
        DO i = 1, totnutr

          CALL cornr(j, i, iuc, ivc)                 !< corner indices for current triangle
          CALL cornr2uv(iuc, ivc, ref, u, v)         !< on-fault coordinates

          DO icr = 1, 3
            w(icr) = roughness(iuc(icr), ivc(icr))
          ENDDO

          CALL uvw2xyz(pl, u, v, w, x, y, z)            !< get cartesian coordinates

          ! this makes triangle degenerate (all "z" = MIN_DEPTH) if all "z" are above MIN_DEPTH
          ! DO icr = 1, 3
          !   z(icr) = MAX(MIN_DEPTH, z(icr))             !< make sure we never breach the free-surface (clip roughness)
          ! ENDDO

          IF (ANY(z .lt. MIN_DEPTH)) THEN
            w(:) = 0._r32                          !< simply set roughness to zero and recompute x, y, z
            CALL uvw2xyz(pl, u, v, w, x, y, z)
          ENDIF

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

              CALL raycorner(sheet, shooting, repi, z, shot, wave, p(:, wave), path(:, wave), trvt(:, wave), q(:, wave))

              IF (ALL(trvt(:, wave) .ne. 0._r32)) THEN
                minbounds(1, wave) = MIN(minbounds(1, wave), mean(path(:, wave)))
                maxbounds(1, wave) = MAX(maxbounds(1, wave), mean(path(:, wave)))
                minbounds(2, wave) = MIN(minbounds(2, wave), mean(trvt(:, wave)))
                maxbounds(2, wave) = MAX(maxbounds(2, wave), mean(trvt(:, wave)))
              ENDIF

            ENDDO

            IF (ALL(trvt .ne. 0._r32)) THEN
              vratio = vratio + mean(path(:, 1)) / mean(trvt(:, 1)) * mean(trvt(:, 2)) / mean(path(:, 2))   !< alpha/beta
              n      = n + 1
            ENDIF

          ENDDO

        ENDDO
      ENDDO  !< end loop over triangles
     !$omp end parallel do

      vratio = vratio / n         !< average Vp/Vs ratio

      coda%lopath = minbounds(1, :)
      coda%uppath = maxbounds(1, :)

      coda%lotrvt = minbounds(2, :)
      coda%uptrvt = maxbounds(2, :)

      IF (input%advanced%verbose .eq. 2) THEN

        ! CALL update_log(num2char('<coda P>', justify='c', width=30) +    &
        !                 num2char('Min trvt', width=15, justify='r')  + '|' + num2char('Max trvt', width=15, justify='r') + '|' +  &
        !                 num2char('Min path', width=15, justify='r')  + '|' + num2char('Max path', width=15, justify='r') + '|' +  &
        !                 num2char('Envelopes', width=15, justify='r') + '|' + num2char('<vp/vs>', width=15, justify='r')  + '|')
        !
        ! np = NINT( (coda%uppath(1) - coda%lopath(1)) / DCP ) + 1
        ! nt = NINT( (coda%uptrvt(1) - coda%lotrvt(1)) / DCT ) + 1
        !
        ! CALL update_log(num2char('', width=30) +  &
        !                 num2char(coda%lotrvt(1), width=15, justify='r', notation='f', precision=2) + '|' +    &
        !                 num2char(coda%uptrvt(1), width=15, justify='r', notation='f', precision=2) + '|' +    &
        !                 num2char(coda%lopath(1), width=15, justify='r', notation='f', precision=2) + '|' +    &
        !                 num2char(coda%uppath(1), width=15, justify='r', notation='f', precision=2) + '|' +    &
        !                 num2char(np*nt, width=15, justify='r') + '|' + &
        !                 num2char(vratio, width=15, justify='r', notation='f', precision=2) + '|', blankline=.false.)

        CALL update_log(num2char('<coda S>', justify='c', width=30) +    &
                        num2char('Min trvt', width=15, justify='r') + '|' + num2char('Max trvt', width=15, justify='r') + '|' + &
                        num2char('Min path', width=15, justify='r') + '|' + num2char('Max path', width=15, justify='r') + '|' + &
                        num2char('Envelopes', width=15, justify='r') + '|' + num2char('<vp/vs>', width=15, justify='r') + '|')

        np = NINT( (coda%uppath(2) - coda%lopath(2)) / DCP ) + 1
        nt = NINT( (coda%uptrvt(2) - coda%lotrvt(2)) / DCT ) + 1

        CALL update_log(num2char('', width=30) +  &
                        num2char(coda%lotrvt(2), width=15, justify='r', notation='f', precision=2) + '|' +    &
                        num2char(coda%uptrvt(2), width=15, justify='r', notation='f', precision=2) + '|' +    &
                        num2char(coda%lopath(2), width=15, justify='r', notation='f', precision=2) + '|' +    &
                        num2char(coda%uppath(2), width=15, justify='r', notation='f', precision=2) + '|' +    &
                        num2char(np*nt, width=15, justify='r')   + '|' + &
                        num2char(vratio, width=15, justify='r', notation='f', precision=2) + '|', blankline=.false.)

      ENDIF

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! -------------------------------------------- compute direct and coda waves -------------------------------------------------


      ! IF (ALLOCATED(coda%penvelope)) DEALLOCATE(coda%penvelope, coda%pdirect)
      !
      ! np = NINT( (coda%uppath(1) - coda%lopath(1)) / DCP ) + 1
      ! nt = NINT( (coda%uptrvt(1) - coda%lotrvt(1)) / DCT ) + 1
      !
      ! ALLOCATE(coda%penvelope(SIZE(time), nt, np), coda%pdirect(nt, np))

      ASSOCIATE(model => input%attenuation(input%receiver(rec)%attenuation))

        gsp = model%gps(band) / 6._r32

        ! wp = 1._r32
        ! ws = 0._r32
        !
        ! !$omp parallel do default(shared) private(i, j, tp, ts, r, vs, d0) reduction(max: ok)
        ! DO j = 1, np
        !   DO i = 1, nt
        !
        !     tp = coda%lotrvt(1) + (i - 1) * DCT
        !     r  = coda%lopath(1) + (j - 1) * DCP
        !
        !     ts = tp * vratio
        !
        !     vs = r / ts
        !
        !     IF (model%gss(band) .ge. GSS_MIN) THEN
        !       CALL rtt(time, tp, ts, model%gpp(band), model%gps(band), gsp, model%gss(band), model%b(band), vs, wp, ws,   &
        !                0.25_r32, d0, coda%penvelope(:, i, j), ok)
        !
        !       coda%pdirect(i, j) = d0(1)
        !     ELSE
        !       coda%penvelope(:, i, j) = 0._r32
        !       coda%pdirect(i, j)      = 1._r32
        !     ENDIF
        !
        !     ! print*, 'P ', tp, r, ts, vs, ok, SIZE(time)
        !
        !   ENDDO
        ! ENDDO
        ! !$omp end parallel do

        np = NINT( (coda%uppath(2) - coda%lopath(2)) / DCP ) + 1
        nt = NINT( (coda%uptrvt(2) - coda%lotrvt(2)) / DCT ) + 1

        IF (ALLOCATED(coda%senvelope)) DEALLOCATE(coda%senvelope, coda%sdirect)

        ALLOCATE(coda%senvelope(SIZE(time), nt, np), coda%sdirect(nt, np))

        ! wp = 0._r32
        ! ws = 1._r32
        wp = 1._r32               !< generic Ep/Es for double-couple in Poisson solids
        ws = 23.4_r32

        ! parallelization here creates problems to "rtt", not sure why (variables "alpha" and "beta" seems overwritten inside "rtt")
        ! As a solution, I commented parallel loop here and enabled parallelization inside "rtt" -> performance penalty
        !$omp parallel do default(shared) private(i, j, tp, ts, r, vs, d0) reduction(max: ok) collapse(2)
        DO j = 1, np
          DO i = 1, nt

            ts = coda%lotrvt(2) + (i - 1) * DCT
            r  = coda%lopath(2) + (j - 1) * DCP

            tp = ts / vratio

            vs = r / ts

            IF (model%gss(band) .ge. GSS_MIN) THEN
              CALL rtt(time, tp, ts, model%gpp(band), model%gps(band), gsp, model%gss(band), model%b(band), vs, wp, ws,   &
                       0.25_r32, d0, coda%senvelope(:, i, j), ok)

              coda%sdirect(i, j) = d0(2)
            ELSE
              coda%senvelope(:, i, j) = 0._r32
              coda%sdirect(i, j)      = 1._r32
            ENDIF

            ! print*, 'S ', ts, r, tp, vs, ok

          ENDDO
        ENDDO
        !$omp end parallel do

      END ASSOCIATE

print*, 'ok coda ', ok

    END SUBROUTINE lookup_coda

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE integration_step(ok, ref, rec, pl, vel, iter, step)

      ! Purpose:
      !   to define the integration time-step "step" such that the desired number of (average) cuts per triangle are obtained. Value
      !   is searched for mesh refinement "ref", receiver "rec", fault plane "pl", velocity model "vel" and iteration "iter".
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   08/03/21                  original version
      !

      INTEGER(i32),              INTENT(OUT) :: ok
      INTEGER(i32),              INTENT(IN)  :: ref, rec, pl, vel, iter
      REAL(r32),                 INTENT(OUT) :: step
      INTEGER(i32)                           :: seed, totnutr, i, j, icr, src, wtp, sheet
      INTEGER(i32), DIMENSION(2)             :: n
      INTEGER(i32), DIMENSION(3)             :: iuc, ivc, shot
      REAL(r32)                              :: strike, dip, urec, vrec, wrec, taumin, taumax
      REAL(r32),    DIMENSION(2)             :: dt
      REAL(r32),    DIMENSION(3)             :: slip, u, v, w, x, y, z, nrl, rake, repi, rupture
      REAL(r32),    DIMENSION(3)             :: p, path, trvt, q, tau

      !-----------------------------------------------------------------------------------------------------------------------------

      ok = 0

      dt(:) = 0._r32
      n(:)  = 0

      ! set "seed" such that random numbers depend on fault plane number and iteration
      seed = input%source%seed + (iter - 1) * SIZE(plane) + pl

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --
      ! ------------------------------------------ define RIK model on mesh nodes ------------------------------------------------

      ALLOCATE(nodes(nugr(ref), nvgr(ref)))

      CALL nodefun(ref, pl, vel, seed)                !< define slip, rupture time and rise time on mesh nodes

      totnutr = 2 * nutr(ref) - 1                     !< total triangles in a row

      ! generate fault plane roughness
      CALL fault_roughness(ok, ref, pl, iter)

      ! shoot rays between receivers and a set of sources spanning (vertically) current mesh
      CALL rayshooting(ok, ref, pl, vel)

      !$omp parallel do default(shared) private(i, j, iuc, ivc, icr, slip, u, v, w, x, y, z, nrl, strike, dip, rake, rupture)   &
      !$omp private(urec, vrec, wrec, repi, wtp, sheet, src, shot, p, path, trvt, q, tau, taumin, taumax)   &
      !$omp reduction(+: dt, n)
      DO j = 1, nvtr(ref)
        DO i = 1, totnutr

          ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
          ! ---------------------------------------- check corners have non-zero slip --------------------------------------------

          CALL cornr(j, i, iuc, ivc)        !< corner indices for current triangle

          DO icr = 1, 3
            slip(icr) = SUM(nodes(iuc(icr), ivc(icr))%slip)
          ENDDO

          IF (ALL(slip .eq. 0._r32)) CYCLE          !< jump to next triangle if current has zero slip

          ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
          ! --------------------------------------------- get perturbed mechanism ------------------------------------------------

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

          strike = plane(pl)%strike !* DEG_TO_RAD
          dip    = plane(pl)%dip    !* DEG_TO_RAD

          CALL interpolate(plane(pl)%u, plane(pl)%v, plane(pl)%rake, u, v, rake)

          CALL perturbed_mechanism(strike, dip, rake, nrl)            !< get new strike, dip, rake (all values in radians)

          ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
          ! ----------------------------------------------- corner parameters ----------------------------------------------------

          DO icr = 1, 3
            rupture(icr) = MINVAL(nodes(iuc(icr), ivc(icr))%rupture)
          ENDDO

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

          ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
          ! -------------------------------------- triangle contribution to ground motion ----------------------------------------

          urec = input%receiver(rec)%x
          vrec = input%receiver(rec)%y
          wrec = 0._r32

          CALL rotate(urec, vrec, wrec, 0._r32, 0._r32, -strike)       !< move to "u-t" coordinates for receiver

          DO icr = 1, 3
            repi(icr) = HYPOT(urec - u(icr), vrec - v(icr)) / 1000._r32        !< epicentral distance (corner-receiver)
          ENDDO

          ! loop over wave types
          DO wtp = 1, 2

            DO sheet = 1, maxsheets(wtp)

              CALL raycorner(sheet, shooting, repi, z, shot, wtp, p, path, trvt, q)

              IF (ANY(trvt .eq. 0._r32)) CYCLE       !< all corners must belong to same sheet

              tau(:) = trvt(:) + rupture(:)          !< sum travel-time and rupture time

              taumax = MAXVAL(tau)
              taumin = MINVAL(tau)

              dt(wtp) = dt(wtp) + (taumax - taumin) / (input%advanced%avecuts + 1)
              n(wtp)  = n(wtp) + 1

            ENDDO     !< end loop over sheets

          ENDDO    !< end loop over wave types

        ENDDO
      ENDDO
      !$omp end parallel do

      CALL dealloc_nodes()

      step = MAXVAL(dt / n)

    END SUBROUTINE integration_step

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

!     FUNCTION assert_i32(par, v1, v2, strict) RESULT(assert)
!
!       ! Purpose:
!       !   To determine if variable "par" belongs to range [v1, v2]. Square brackets are replaced by round brackets if "strict=.true."
!       !
!       ! Revisions:
!       !     Date                    Description of change
!       !     ====                    =====================
!       !   18/12/20                  original version
!       !
!
!       INTEGER(i32), DIMENSION(:),           INTENT(IN) :: par
!       INTEGER(i32),               OPTIONAL, INTENT(IN) :: v1, v2
!       LOGICAL,                    OPTIONAL, INTENT(IN) :: strict
!       LOGICAL                                          :: flag, assert
!
!       !-----------------------------------------------------------------------------------------------------------------------------
!
! #include "assert_incl.f90"
!
!     END FUNCTION assert_i32
!
!     ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
!     !===============================================================================================================================
!     ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
!
!     FUNCTION assert_r32(par, v1, v2, strict) RESULT(assert)
!
!       ! Purpose:
!       !   To determine if variable "par" belongs to range [v1, v2]. Square brackets are replaced by round brackets if "strict=.true."
!       !
!       ! Revisions:
!       !     Date                    Description of change
!       !     ====                    =====================
!       !   18/12/20                  original version
!       !
!
!       REAL(r32), DIMENSION(:),           INTENT(IN) :: par
!       REAL(r32),               OPTIONAL, INTENT(IN) :: v1, v2
!       LOGICAL,                 OPTIONAL, INTENT(IN) :: strict
!       LOGICAL                                       :: flag, assert
!
!       !-----------------------------------------------------------------------------------------------------------------------------
!
! #include "assert_incl.f90"
!
!     END FUNCTION assert_r32

END MODULE m_isochron
