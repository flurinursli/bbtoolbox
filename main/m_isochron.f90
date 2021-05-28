MODULE m_isochron

  USE                :: omp_lib
  USE, NON_INTRINSIC :: m_precisions
  USE, NON_INTRINSIC :: m_interpolation_r32
  USE, NON_INTRINSIC :: m_compgeo
  USE, NON_INTRINSIC :: m_stat
  USE, NON_INTRINSIC :: m_source
  USE, NON_INTRINSIC :: m_rik
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

  PUBLIC :: solve_isochron_integral, node2disk

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

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  COMPLEX(r32), ALLOCATABLE, DIMENSION(:) :: fxp, fzp                 !< free-surface amplification coefficients
  INTEGER(i32),              DIMENSION(2) :: minsheets, maxsheets
  REAL(r32),    ALLOCATABLE, DIMENSION(:) :: shooting                 !< depth of shooting points

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  PROCEDURE(rik_at_nodes), POINTER :: nodefun

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  CONTAINS

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE solve_isochron_integral(ok, pl, vel, iter)

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

      INTEGER(i32),                             INTENT(OUT) :: ok
      INTEGER(i32),                             INTENT(IN)  :: pl, vel, iter
      COMPLEX(r32)                                          :: fxpfzs, fzpfxs, g31, g21, g23, ga, gb
      COMPLEX(r32),              DIMENSION(3,3)             :: fp, fs, g
      COMPLEX(r32), ALLOCATABLE, DIMENSION(:)               :: fxp, fzp
      INTEGER(i32)                                          :: n, i, j, ref, icr, totnutr, seed, src, rec, sheet, wtp, it1, it2, ic
      INTEGER(i32)                                          :: it, icut, iat, band
      INTEGER(i32),              DIMENSION(3)               :: iuc, ivc, abv, blw, shot
      INTEGER(i32), ALLOCATABLE, DIMENSION(:,:)             :: ncuts
      REAL(r32)                                             :: m0, moment, strike, dip, urec, vrec, wrec, taumin, taumax, sd, cd
      REAL(r32)                                             :: srfvel, velloc, stho, r, dx, dy, dr, cpsi, spsi, sloloc, signp, sthf
      REAL(r32)                                             :: cthf, rn, rv, ru, bn, cn, bu, cu, bv, cv, rnbu, rnbv, rncu, rncv, guk
      REAL(r32)                                             :: ctho, tau31, tau21, tau23, c, du, dv, p13, p12, p32, dl, scale
      REAL(r32)                                             :: avepath, eta
      REAL(r32),                 DIMENSION(2)               :: po, ro, to, qo, zshots
      REAL(r32),                 DIMENSION(3)               :: u, v, w, x, y, z, slip, rise, rupture, rake, nrl, rho, alpha, beta
      REAL(r32),                 DIMENSION(3)               :: mu, p, q, path, repi, tau, sr, cr, rvec, thvec, phvec, ir, itheta
      REAL(r32),                 DIMENSION(3)               :: iphi, dtau
      REAL(r32),    ALLOCATABLE, DIMENSION(:)               :: fsp

      !-----------------------------------------------------------------------------------------------------------------------------

      CALL setup_interpolation('linear', 'zero', ok)

      IF (ok .ne. 0) RETURN

      CALL make_fftw_plan([n])
      CALL make_iir_plan(ok, 'butter', dt, fbands(:, k), 'pass', 2, zphase = .true.)

      CALL destroy_fftw_plan([n])
      CALL destroy_iir_plan()

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! ---------------------------------------- lookup table free-surface coefficients --------------------------------------------

      alpha = BIG
      beta  = BIG

      ! find minimum velocities amongst all receivers
      DO rec = 1, SIZE(input%receiver)
        alpha(1) = MIN(alpha(1), input%velocity(input%receiver(rec)%velocity)%vp(1))
        beta(1)  = MIN(beta(1), input%velocity(input%receiver(rec)%velocity)%vs(1))
      ENDDO

      n = NINT(1._r32 / beta(1) / DP) + 1      !< coefficients are computed at "n" slowness points

      ALLOCATE(fxp(n), fzp(n), fsp(n))

      ! determine coefficients valid for all receivers
      CALL setf(ok, alpha(1), beta(1), 1, fsp, fxp, fzp)

      IF (ok .ne. 0) RETURN

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! --------------------------------------- solve integral for each mesh refinement --------------------------------------------

      moment = 0._r32      !< keep track of actual moment

      nodefun => rik_at_nodes

      IF (input%source%is_point) nodefun => ptrsrc_at_nodes

      ! set "seed" such that random numbers depend on fault plane number and iteration
      seed = input%source%seed + (iter - 1) * SIZE(plane) + pl

      ALLOCATE(ncuts(2, SIZE(input%receiver)))

      ncuts(:,:) = 0        !< keep track of mean isochron cuts for each wave type and receiver

      DO ref = 1, SIZE(nvtr)         !< loop over mesh refinements

        ALLOCATE(nodes(nugr(ref), nvgr(ref)))

        CALL nodefun(ref, pl, vel, seed)                !< define slip, rupture time and rise time on mesh nodes

        totnutr = 2 * nutr(ref) - 1                     !< total triangles in a row

        m0 = 0._r32

        ! generate fault plane roughness
        CALL fault_roughness(ok, ref, pl, iter)

        ! shoot rays between receivers and a set of sources spanning (vertically) current mesh
        CALL rayshooting(ok, ref, pl, vel)

        !$omp parallel do default(shared) private(i, j, iuc, ivc, icr, slip, rise, rupture, u, v, w, x, y, z, nrl, strike, dip)  &
        !$omp private(rake, alpha, beta, rho, mu, wtp, sheet, rec, urec, vrec, wrec, src, shot, iab, ibl, zshots, po, xo, ro, to) &
        !$omp private(qo, p, q, r, path, trvt, repi, it1, it2, taumin, taumax, tau, cd, sd, sr, cr) reduction(+: m0)
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
              w(icr) = roughness(iuc(icr), ivc(icr))      !< add roughness
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

            CALL perturbed_mechanism(strike, dip, rake, nrl)            !< new strike, dip, rake (all values in radians)

            sd = SIN(dip)
            cd = COS(dip)

            DO icr = 1, 3
              sr(icr) = SIN(rake(icr))
              cr(icr) = COS(rake(icr))
            ENDDO

            ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
            ! ----------------------------------------------- corner parameters ----------------------------------------------------

            DO icr = 1, 3
              ! rise(icr)    = mean(nodes(iuc(icr), ivc(icr))%rise)
              rupture(icr) = MINVAL(nodes(iuc(icr), ivc(icr))%rupture)
              alpha(icr)   = vinterp(input%velocity(vel)%depth, input%velocity(vel)%vp, input%velocity(vel)%vpgrad, z(icr))
              beta(icr)    = vinterp(input%velocity(vel)%depth, input%velocity(vel)%vs, input%velocity(vel)%vsgrad, z(icr))
              rho(icr)     = vinterp(input%velocity(vel)%depth, input%velocity(vel)%rho, input%velocity(vel)%rhograd, z(icr))
              mu(icr)      = rho(icr) * beta(icr)**2         !< rigidity
            ENDDO

            ! move from "x-y" to "u-t" coordinates for corners
            u = x
            v = y
            w = 0._r32
            CALL rotate(u, v, w, 0._r32, 0._r32, -strike)

            ! find shooting points above/below
            DO icr = 1, 3
              DO src = 1, SIZE(shooting) - 1
                IF ( (z(icr) .ge. shooting(src)) .and. (z(icr) .lt. shooting(src + 1)) ) shot(icr) = src
              ENDDO
            ENDDO

            ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
            ! ---------------------------------------- moment rate function of triangle --------------------------------------------

            CALL mrf_rik(iuc, ivc, mrf)     !< return normalized mrf (area ~ 1)

            ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
            ! -------------------------------------- triangle contribution to ground motion ----------------------------------------

            ! loop over receivers
            DO rec = 1, SIZE(input%receiver)

              IF (input%receiver(rec)%velocity .ne. vel) CYCLE         !< consider only receivers for current velocity model "vel"

              urec = input%receiver(rec)%x
              vrec = input%receiver(rec)%y
              wrec = 0._r32

              CALL rotate(urec, vrec, wrec, 0._r32, 0._r32, -strike)       !< move to "u-t" coordinates for receiver

              DO icr = 1, 3
                repi(icr) = HYPOT(urec - u(icr), vrec - v(icr))        !< epicentral distance (corner-receiver)
              ENDDO

              iat = input%receiver(rec)%attenuation      !< attenuation model

              ! loop over wave types (i.e. 1=P, 2=S)
              DO wtp = 1, 2

                DO icr = 1, 3
                  abv(icr) = 1      !< reset sheet-related index
                  blw(icr) = 1
                ENDDO

                ! velocity at receiver and source
                IF (wtp .eq. 1) THEN
                  srfvel = input%velocity(input%receiver(rec)r%velocity)%vp(1)
                  velloc = alpha(icr)
                ELSE
                  srfvel = input%velocity(input%receiver(rec)%velocity)%vs(1)
                  velloc = beta(icr)
                ENDIF

                ! loop over sheets
                DO sheet = 1, maxsheets(wtp)    !< loop over sheets for each wave type

                  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
                  ! ----------------------------------------- ray parameters at corners --------------------------------------------

                  DO icr = 1, 3

                    CALL interpray(abv(icr), sr(icr), shot(icr), wtp, po(1), ro(1), to(1), qo(1))        !< shot-point above
                    CALL interpray(blw(icr), sr(icr), shot(icr) + 1, wtp, po(1), ro(2), to(2), qo(2))    !< shot-point below

                    tau(icr) = -1._r32

                    ! interpolate values at depth "z(icr)" only if corner is fully inside sheet
                    IF (ALL(to .ne. 0._r32)) THEN

                      zshots = [shooting(shot(icr)), shooting(shot(icr) + 1)]      !< shooting points depth vector

                      CALL interpolate(zshots, po, z(icr), p(icr))
                      CALL interpolate(zshots, ro, z(icr), path(icr))
                      CALL interpolate(zshots, to, z(icr), trvt)
                      CALL interpolate(zshots, qo, z(icr), q(icr))

                      ! sum travel-time and rupture time
                      tau(icr) = trvt + rupture(icr)

                    ENDIF

                  ENDDO

                  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
                  ! -------------------------------------------- time integration limits -------------------------------------------

                  IF (ANY(tau .lt. 0._r32)) CYCLE       !< all corners must belong to same sheet

                  taumax = MAXVAL(tau)
                  taumin = MINVAL(tau)

                  it1 = NINT(taumin / timeseries%sp%dt) + 1

                  it2 = it1 + (taumax - taumin) / timeseries%sp%dt - 1      !< when (taumax-taumin) < dt, it2 < it1
                  it2 = MIN(SIZE(timeseries%sp%time), it2)                  !< limit to max number of time points

                  ! jump to next sheet if no isochron is spanned
                  IF (it2 .lt. it1) CYCLE

                  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
                  ! ---------------------------------------------- integral kernel -------------------------------------------------

                  DO icr = 1, 3

                    dx = urec - u(icr)
                    dy = vrec - v(icr)

                    dr = HYPOT(dx, dy)

                    ! cos and sin of angle psi
                    cpsi = dx / dr
                    spsi = dy / dr

                    sloloc = 1._r32 / velloc        !< slowness at corner

                    ! "p" is "-p" for downgoing rays
                    IF (p(icr) .gt. sloloc) p(icr) = -(2._r32 * sloloc - p(icr))

                    ! get sign of ray (+ if upgoing, - if downgoing)
                    signp = SIGN(1._r32, p(icr))

                    p(icr) = ABS(p(icr))         !< from here onwards we need only its absolute value

                    ! sine of angle eta at source
                    sthf = velloc * p(icr)

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
                      ir(icr) = 2._r32 * mu(icr) * q(icr) * rn * (ru * cr(icr)) - rv * sr(icr))) * slip(icr)
                    ELSE
                      ! dot products
                      bn = (cthf * spsi * sd - sthf * cd)
                      cn = (cpsi * sd)
                      bu = -cthf * cpsi
                      cu = spsi
                      bv = -cthf * spsi * cd - sthf * sd
                      cv = -cpsi * cd

                      rnbu = (rn * bu + bn*ru)
                      rnbv = (rn * bv + bn*rv)
                      rncu = (rn * cu + cn*ru)
                      rncv = (rn * cv + cn*rv)

                      guk = mu(icr) * q(icr)

                      itheta(icr) = slip(icr) * guk * (rnbu * cr(icr)) - rnbv * sr(icr)))
                      iphi(icr)   = slip(icr) * guk * (rncu * cr(icr)) - rncv * sr(icr)))
                    ENDIF

                    stho = p(icr) * sfrvel

                    #ifdef ERROR_TRAP
                    IF (stho .gt. 1.01_r32) CALL report_error('solve_isochron_integral - ERROR: sin(eta) sensibly larger than 1')
                    #endif

                    stho = MIN(stho, 1._r32)     !< handle roundoff errors
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
                    CALL interpolate(fsp, fxp, p(icr), fxpfzs)
                    CALL interpolate(fsp, fzp, p(icr), fzpfxs)

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

                  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
                  ! --------------------------------------------- time integration -------------------------------------------------

                  ! differences in tau between corners
                  tau31 = tau(3) - tau(1)
                  tau21 = tau(2) - tau(1)
                  tau23 = tau(2) - tau(3)

                  ! apparent velocity of the tau curves
                  c = HYPOT(tau31 / dutr(ref), (tau(2) - 0.5_r32 * (tau(1) + tau(3))) / dvtr(ref))
                  c = 1._r32 / c

                  ! path averaged over corners
                  avepath = mean(path)

                  DO ic = 1, 3

                    g31 = g(ic, 3) - g(ic, 1)
                    g21 = g(ic, 2) - g(ic, 1)
                    g23 = g(ic, 2) - g(ic, 3)

                    DO it = it1, it2

                      ! t = (it - 1) * timeseries%sp%dt
                      ! t = timeseries%sp%time(it)
                      dtau(1) = timeseries%sp%time(it) - tau(1)
                      dtau(2) = timeseries%sp%time(it) - tau(2)
                      dtau(3) = timeseries%sp%time(it) - tau(3)

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

                      ! keep track of total cuts
                      IF ( (band .eq. 1) .and. (ic .eq. 1) ) THEN
                        IF (icut .ne. 0) ncuts(wtp, rec) = ncuts(wtp, rec) + (it(2) - it(1) + 1)
                      ENDIF

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

                      triseis(it, ic) = c * (ga + gb) * dl * 0.5_r32

                    ENDDO      !< end loop over time (cuts)

                  ENDDO     !< end loop over components of motion

                  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
                  ! --------------------------------------------- attenuation & coda -----------------------------------------------

                  ! loop over frequency bands
                  DO band = 1, SIZE(input%attenuation(iat)%gss)

                    eta = input%attenuation(iat)%b(band)

                    IF (wtp .eq. 1) THEN
                      eta = eta + input%attenuation(iat)%gpp(band) + input%attenuation(iat)%gps(band)       !< P
                    ELSE
                      eta = eta + input%attenuation(iat)%gps(band) / 6._r32 + input%attenuation(iat)%gss(band)    !< S
                    ENDIF

                    attenuation = EXP(-mean(path) * eta(band, wtp) / 2._r32)        !< attenuation factor direct wave (amplitude, not energy)

                    DO ic = 1, 3
                      DO it = it1, it2
                        direct(it, ic, band, rec) = direct(it, ic, band, rec) + triseis(it, ic) * attenuation
                      ENDDO
                    ENDDO

                    ! stack envelopes here (p-only, s-only, depending on wtp), weighted by "Q" and "dl"

                  ENDDO      !< end loop over frequency bands

                ENDDO     !< end loop over sheets

              ENDDO    !< end loop over wave types

              ! add hilbert transform here
              CALL fft(AIMAG(triseis(:, ic, band)), spectrum)

              spectrum(:) = spectrum(:) * hlb(:)

              CALL ifft(analytic, spectrum)

              seis = CMPLX(REAL(rseis, r32), -AIMAG(analytic), r32)




              ! move back from "u/t" (fault-parallel/fault-normal) to "x/y" coordinates
              CALL rotate(spx, spy, spz, 0._r32, 0._r32, strike)

              timeseries%sp%x(it, rec) = timeseries%sp%x(it, rec) + spx(it)
              timeseries%sp%y(it, rec) = timeseries%sp%y(it, rec) + spy(it)
              timeseries%sp%z(it, rec) = timeseries%sp%z(it, rec) + spz(it)

              ! add coda here

              ! convolve with mrf here

              ! add bandpass filter here

            ENDDO   !< end loop over receivers

            m0 = m0 + mean(mu * slip)        !< average moment contribution (area is added later)

          ENDDO
        ENDDO
        !$omp end parallel do

        CALL dealloc_nodes()

        moment = moment + m0 * dutr(ref) * dvtr(ref) * 0.5_r32      !< sum is over mesh refinements

      ENDDO        !< end loop over mesh refinements

      scale = plane(pl)%targetm0 / moment

      ! rescale synthetics to target moment
      DO rec = 1, SIZE(input%receiver)
        DO it = 1, SIZE(timeseries%sp%time)
          timeseries%sp%x(it, rec) = timeseries%sp%x(it, rec) * scale
          timeseries%sp%y(it, rec) = timeseries%sp%y(it, rec) * scale
          timeseries%sp%z(it, rec) = timeseries%sp%z(it, rec) * scale
        ENDDO
      ENDDO



    END SUBROUTINE solve_isochron_integral

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE setfs(ok, vp, vs, irort, fsp, fxp, fzp)

      INTEGER(i32),               INTENT(OUT) :: ok
      REAL(r32),                  INTENT(IN)  :: vp, vs
      INTEGER(i32),               INTENT(IN)  :: irort
      COMPLEX(r32), DIMENSION(:), INTENT(OUT) :: fsp, fxp, fzp
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
        IF ( (AIMAG(fxp(i)) .ne. 0._r32) .or. (AIMAG(fzp(i)) .ne. 0._r32) ) THEN
          CALL report_error('setfs - ERROR: coefficient fxp and/or fzp not purely real for p < 1/vp')
          RETURN
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

      nrmn = dpsm / dp + 2          !< points in the running mean window

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
              w(icr) = roughness(iuc(icr), ivc(icr))
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
          CALL cornr2uv(iuc, ivc, ref, u, v)        !< on-fault coordinates

          DO icr = 1, 3
            w(icr) = roughness(iuc(icr), ivc(icr))
          ENDDO

          CALL uvw2xyz(pl, u, v, w, x, y, z)            !< get cartesian coordinates

          ! compute max source-receiver distance
          DO icr = 1, 3
            DO rec = 1, SIZE(input%receiver)
              rmax = MAX(rmax, HYPOT(x(icr) - input%receiver(rec)%x, y(icr) - input%receiver(rec)%y) / 1000._r32)
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
                        num2char('Separation', width=15, justify='r') + '|')

        CALL update_log(num2char('', width=30, justify='c')  +  &
                        num2char(SIZE(shooting), width=15, justify='r') + '|' + &
                        num2char(rmax, width=15, notation='f', precision=1, justify='r') + '|' + &
                        num2char(MINVAL(shooting), width=15, notation='f', precision=2, justify='r') + '|' + &
                        num2char(MAXVAL(shooting), width=15, notation='f', precision=2, justify='r') + '|' +  &
                        num2char(num2char(zmin, notation='f', width=6, precision=2) + ', ' +   &
                                 num2char(zmax, notation='f', width=6, precision=2), width=15, justify='r') + '|',blankline=.false.)

      ENDIF

      rmax = 1.2 * rmax            !< slightly increase maximum distance

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

        CALL update_log(num2char('<sheets>', justify='c', width=30) + num2char('P-wave', width=15, justify='r') + '|' +  &
                        num2char('S-wave', width=15, justify='r') + '|')

        CALL update_log(num2char('', width=30, justify='c')  +  &
                        num2char(num2char(minsheets(1), width=2) + ', ' +   &
                                 num2char(maxsheets(1), width=2), width=15, justify='r') + '|' + &
                        num2char(num2char(minsheets(2), width=2) + ', ' +   &
                                 num2char(maxsheets(2), width=2), width=15, justify='r') + '|',  blankline=.false.)

      ENDIF

#ifdef PERF
      IF (input%advanced%verbose .eq. 2) THEN
        CALL update_log(num2char('<<elapsed time>>', justify='c', width=30) + num2char(tictoc(1), width=15, notation='f',   &
                        precision=3, justify='r') + '|', blankline=.false.)
      ENDIF
#endif

    END SUBROUTINE rayshooting

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *


    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *


    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *


END MODULE m_isochron
