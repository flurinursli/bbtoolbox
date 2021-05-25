MODULE m_wkbj

  USE, NON_INTRINSIC :: m_precisions
  USE, NON_INTRINSIC :: m_logfile
  USE, NON_INTRINSIC :: m_toolbox, ONLY: input

  IMPLICIT none

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  PRIVATE

  PUBLIC :: spred, interpray
  PUBLIC :: wkbj

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  INTEGER(i32), PARAMETER :: MAX_P_RAYS = 100, MAX_S_RAYS = 50
  REAL(r64),    PARAMETER :: PI = 3.14159265358979323846_r64
  REAL(r64),    PARAMETER :: DMAX = 0.2_r64                                !< step damping

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  TYPE :: ray
    REAL(r64), ALLOCATABLE, DIMENSION(:) :: pn, x, t, q, r
  END TYPE ray

  TYPE(ray), ALLOCATABLE, DIMENSION(:,:) :: wkbj

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  CONTAINS

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE spred(ok, phys, xmax, zsrc, src, wtp, lu, sheets)

      INTEGER(i32),                                          INTENT(OUT) :: ok
      REAL(r32),                 DIMENSION(:,:),             INTENT(IN)  :: phys
      REAL(r32),                                             INTENT(IN)  :: xmax, zsrc
      INTEGER(i32),                                          INTENT(IN)  :: src, wtp, lu
      INTEGER(i32),                                          INTENT(OUT) :: sheets

      COMPLEX(r64)                                                       :: frp, fthp, frs, fths, fr, fth, supd
      INTEGER(i32)                                                       :: imth, nn, nmod, j, ilay, i, ii, is, ir, num, lays
      INTEGER(i32)                                                       :: icns, icnr, isheet, itypef, itypeo, nveldim
      INTEGER(i32)                                                       :: npdif, nsdif, itry, ierr, jray, nrays
      INTEGER(i32),              DIMENSION(SIZE(phys,2)+2)               :: icn
      INTEGER(i32),              DIMENSION(20,2)                         :: nps
      INTEGER(i32),              DIMENSION(2)                            :: jmin, jmax
      INTEGER(i32), ALLOCATABLE, DIMENSION(:)                            :: itmp
      LOGICAL                                                            :: upgoin, spconv
      REAL(r64)                                                          :: rsur, dum, sygn, z, za, zb, dz, x, xx, xi, zr, zs, r
      REAL(r64)                                                          :: dpinit, ulast, dudp, dumin, sdudp, du, pmax, pmin, domda
      REAL(r64)                                                          :: dadzs, dbdzs, dadzr, dbdzr, drhdzs, drhdzr, sdomda
      REAL(r64)                                                          :: rsor, rrec, dp, pp, umin, umax, pplast, u, ustep, psourc
      REAL(r64)                                                          :: alff, alfo, betf, beto, rhof, rhoo, pta, ptb, cmin, cmax
      REAL(r64)                                                          :: velo, velf, tho, thf, qp, travtm, cvel, pn, vp, vs, tau
      REAL(r64),                                             PARAMETER   :: one = 1.000000000001_r64
      REAL(r64),                 DIMENSION(10000)                        :: vpn, vxx, vtt, vqp, vr, vsheet, vtypef, vtypeo
      REAL(r64),                 DIMENSION(SIZE(phys,2)+2)               :: rho
      REAL(r64),                 DIMENSION(6,SIZE(phys,2)+2)             :: vel
      REAL(r64),    ALLOCATABLE, DIMENSION(:)                            :: tmp

      !-----------------------------------------------------------------------------------------------------------------------------

      ok = 0

      sheets = 0

      nps(:,:) = 0
      jmin(:)  = 0
      jmax(:)  = 0

      rsur = 0._r64
      imth = 2
      nn   = 0

      frp  = (1._r64, 0._r64)
      fthp = (1._r64, 0._r64)
      frs  = (1._r64, 0._r64)
      fths = (1._r64, 0._r64)

      ! current version supports only receivers at free-surface
      zr = 0._r64
      zs = zsrc

      ! assign 1d velocity model values to local variables
      vel(:,:) = 0._r64
      rho(:) = 0._r64

      nveldim = SIZE(phys, 2)

      DO i = 1, nveldim
        vel(1:3, i) = phys(1:3, i)
        rho(i)      = phys(4, i)
      ENDDO

      nmod = nveldim

      ! determine indicies of model points which are tops of layers. Put these in icn. Ihis is needed IF we have layers with
      ! constant velocity values
      ilay = 0
      z = vel(1,1)
      DO i = 1,nmod
        IF (vel(1,i) .eq. z) THEN
          ilay = ilay + 1
          icn(ilay) = i
        ELSE
          z = vel(1,i)
        ENDIF
      ENDDO

      ! declare bottom of model an interface
      ilay = ilay + 1
      icn(ilay) = nmod

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! -------------------------------------- find layer indices for source and receiver ------------------------------------------

      ! "is, ir" is the model index at or above the (source,receiver).
      is = 0
      ir = 0
      dadzs = 0._r64; dbdzs = 0._r64; dadzr = 0._r64; dbdzr = 0._r64; drhdzs = 0._r64; drhdzr = 0._r64

      DO i = 1, nmod-1

        za = vel(1,i)
        zb = vel(1,i+1)
        dz = zb - za

        IF (dz .eq. 0.) CYCLE

        IF ( .not. ( (zs .ge. za .and. zs .lt. zb) .or. (zs .le. za .and. zs .gt. zb) )) THEN
          IF ( .not. ((zr .ge. za .and. zr .lt. zb) .or. (zr .le. za .and. zr .gt. zb)) ) CYCLE

          ! receiver is between za and zb
          ir = i
          dadzr = (vel(2,i+1) - vel(2,i)) / dz
          dbdzr = (vel(3,i+1) - vel(3,i)) / dz
          drhdzr = (rho(i+1) - rho(i)) / dz

          CYCLE
        ENDIF

        ! source is between za and zb
        is = i
        drhdzs = (rho(i+1) - rho(i)) / dz
        dadzs = (vel(2,i+1) - vel(2,i)) / dz
        dbdzs = (vel(3,i+1) - vel(3,i)) / dz

        IF ( .not. ((zr .ge. za .and. zr .lt. zb) .or. (zr .le. za .and. zr .gt. zb)) ) CYCLE

        ! receiver is between za and zb
        ir = i
        dadzr = (vel(2,i+1) - vel(2,i)) / dz
        dbdzr = (vel(3,i+1) - vel(3,i)) / dz
        drhdzr = (rho(i+1) - rho(i)) / dz

      ENDDO

      IF (is .eq. 0) THEN
        CALL report_error('spred - ERROR: is=0=model index above source! Check model and source depth and ray description')
        ok = 1
        RETURN
      ENDIF

      IF (ir .eq. 0) THEN
        CALL report_error('spred - ERROR: ir=0=model index above source! Check model and source depth and ray description')
        ok = 1
        RETURN
      ENDIF

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! ------------------------------ add fictitious interface at source and receiver depth ---------------------------------------

      ! icnr and icns are the model indices corresponding to the source and receiver
      icns = 0
      icnr = 0

      ! check for preexisting interfaces at zs or zr
      DO i = 1, ilay
         IF (zs .eq. vel(1,icn(i))) icns = icn(i)
         IF (zr .eq. vel(1,icn(i))) icnr = icn(i)
      ENDDO

      ! zs > zr
      IF (zs .gt. zr) THEN

         IF ( icns .eq. 0 ) THEN
            CALL inser2d (nmod, zs, is, vel, rho, dadzs, dbdzs, drhdzs)
            icns = is + 2
         ENDIF

         IF (icnr .eq. 0) THEN
            CALL inser2d (nmod, zr, ir, vel, rho, dadzr, dbdzr, drhdzr)
            icns = icns + 2
            icnr = ir + 2
         ENDIF

      ! zr > zs
      ELSEIF (zr .gt. zs) THEN

         IF (icnr .eq. 0) THEN
            CALL inser2d (nmod, zr, ir, vel, rho, dadzr, dbdzr, drhdzr)
            icnr = ir + 2
         ENDIF

         IF (icns .eq. 0) THEN
            CALL inser2d (nmod, zs, is, vel, rho, dadzs, dbdzs, drhdzs)
            icnr = icnr + 2
            icns = is + 2
         ENDIF

      ! zr = zs
      ELSE

         IF (icnr .eq. 0) THEN
            CALL inser2d (nmod, zr, ir, vel, rho, dadzs, dbdzs, drhdzs)
            icnr = ir + 2
            icns = icnr
         ENDIF

      ENDIF

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! ---------------------------------------- define velocities at source and receiver ------------------------------------------

      alff = vel(2,icns); alfo = vel(2,icnr)
      betf = vel(3,icns); beto = vel(3,icnr)
      rhof = rho(icns);   rhoo = rho(icnr)
      pta = 1._r64/ alff; ptb = 1._r64/ betf

#ifdef DEBUG
      WRITE(lu,*)
      WRITE(lu,103) nmod, ((vel(j, i), j=1,3), rho(i), i=1,nmod)
      103 FORMAT (T8, 'Z', T23, 'VP', T38, 'VS', T51, 'DENS', T60, 'NMOD=', T70, I5, / (4F15.8))
#endif

      ilay = 0

      za = vel(1,1)

      ! change velocity model to slowness model and fill "vel(4-6,1-nmod)" with interpolation coefficients.
      DO i = 1, nmod

        z = za
        za = vel(1, i)
        vp = vel(2, i)
        vs = vel(3, i)

        vel(2, i) = 1._r64 / vel(2,i)
        IF (vel(3, i) .ne. 0.) vel(3,i) = 1._r64 / vel(3,i)

        IF (z .eq. za) THEN
          ilay = ilay + 1
          icn(ilay) = i
          zb = vel(1, i)
          CYCLE
        ENDIF

        ! interpolation coefficients
        ii = i - 1
        z = zb
        zb = vel(1,i)

        DO j = 2, 3
          x = vel(j, ii)
          xi = vel(j, i)

          IF (x .ne. xi) THEN
            vel(j + 3, ii) = (xi - x) / (z - zb) / x / xi
          ELSE
            vel(j + 3, ii) = 0.5_r64 * (x - xi) * (x + xi) / (z - zb)
          ENDIF
        ENDDO

      ENDDO

#ifdef DEBUG
      WRITE(lu,104) (vel(i,icns), i=1,3), rho(icns)
      104 FORMAT ('*********' / 'Z, PSLO, SSLO, DENS=', T22, 4F15.8 / T7, 'PN', T22, 'X', T37, 'T', T52, 'Q', T60, 'ISHEET')
#endif

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! ------------------------------------------ loop over four rays: two up, two down -------------------------------------------

      nrays = 0

      ray: DO jray = 1, 4

        spconv = .false.
        isheet = 1
        sdomda = 1._r64

        nps(:, :) = 0
        nn = 0

        ! parameters for p-waves and s-waves
        IF (wtp .eq. 1) THEN

          num = MAX_P_RAYS
          itypef = 1
          itypeo = 1
          cmin = 0._r64
          nps(1, 1) = 1

          IF (jray .eq. 1) THEN
            cmax = 10000._r64
            nps(2, 1) = 0
            nn = 1
          ELSEIF (jray .eq. 2) THEN
            cmax = 0._r64
            nps(2, 1) = 2
            nn = 2
          ENDIF

        ELSEIF (wtp .eq. 2) THEN

          num = MAX_S_RAYS
          itypef = 2
          itypeo = 2
          cmin = -1._r64
          nps(1, 2) = 1

          IF (jray .eq. 1) THEN
            cmax = 10000._r64
            nn = 1
          ELSEIF (jray .eq. 2) THEN
            cmax = -1._r64
            nps(2, 2) = 2
            nn = 2
          ENDIF

        ENDIF

        upgoin = .true.

        IF (nn .eq. 0) CYCLE

        IF (cmax .eq. -9998._r64) cmax = MIN(1._r64 / vel(2, icns), 1._r64 / vel(3, nmod))
        IF (cmax .eq.  9998._r64) cmax = 1._r64 / vel(2, icns)
        IF (cmin .eq. -9999._r64) cmin = MAX ( (1._r64 / vel(2, 1) + 0.01_r64), 1._r64 / vel(3, icns) )
        IF (cmin .eq. 9999._r64)  cmin = MIN(1._r64 / vel(2, icns), 1._r64 / vel(3, nmod))

        IF (cmin .lt. 0._r64) THEN
          ! s source
          pmax = vel(3, icns)
          psourc = pmax
        ELSEIF (cmin .eq. 0._r64) THEN
          ! p source
          pmax = vel(2, icns)
          psourc = pmax
        ELSE
          pmax = 1._r64 / cmin
        ENDIF

        cmin = 1._r64 / pmax

        IF (cmax .lt. 0.) THEN
          ! rays bottom in s mode
          pmin = vel(3, nmod) * one
        ELSEIF (cmax .eq. 0.) THEN
          ! rays bottom in p mode
          pmin = vel(2, nmod) * one
        ELSEIF ( (cmax .gt. 0._r64) .and. (cmax .le. 9999._r64) ) THEN
          pmin = 1._r64 / cmax
        ELSE
          pmin = 0._r64
        ENDIF
        IF (pmin .ne. 0._r64) cmax = 1._r64 / pmin

        IF (jmin(1) .ne. 0) pmin = vel(jmin(1) + 1, jmin(2)) * one
        IF (jmax(1) .ne. 0) pmax = vel(jmax(1) + 1, jmax(2)) / one

        IF (cmax .lt. cmin) CYCLE

        ! this is a fudge to reverse the order of the table for rays leaving  the source travelling downward. First determine what
        ! layer the source is in
        lays = 1
        DO i = 1,icns
          IF ( (vel(1, i + 1) - vel(1, i)) .eq. 0._r64) lays = lays + 1
        ENDDO

        !..the ray type at the source is determined by which ray type has one more ray segment
        !  above the source than below (or vice versa). the direction of ray propagation (up or
        !  downgoing) is determined by which of the vice or versa is true above.
        !  note.... this algorithm assumes that it is impossible for a ray to change direction
        !  (reflect) or mode convert at the source interface (i.e. the source interface is
        !  completely transparent to all rays). thus, you are not allowed to put the source at a
        !  velocity discontinuity.

        npdif = nps(lays - 1, 1) - nps(lays, 1)
        nsdif = nps(lays - 1, 2) - nps(lays, 2)

        IF (npdif .eq. 1) THEN
          itypef = 1
          upgoin = .true.
        ELSEIF (npdif .eq. -1) THEN
          itypef = 1
          upgoin = .false.
        ELSEIF (nsdif .eq. 1) THEN
          itypef = 2
          upgoin = .true.
        ELSEIF (nsdif .eq. -1) THEN
          itypef = 2
          upgoin = .false.
        ELSE
          CALL report_error('spred - ERROR: bad npdif or nsdif. Remove discontinuity from v(z)')
          ok = 1
          RETURN
        ENDIF

        ! this assumes the conversion is at the free surface
        IF ( (itypef .eq. 2) .and. (itypeo .eq. 1) ) spconv = .true.

        !  downgoing rays
        IF (.not.upgoin) THEN
          dum = pmax
          pmax = pmin
          pmin = dum
        ENDIF

        ! determine an average step in ray parameter to use
        dp = (pmax - pmin) / REAL(num - 1, r64)
        rsor = vel(1, icns)
        rrec = vel(1, icnr)

        pp = pmin - dp
        itry = 0
        dpinit = dp
        ulast = -1._r64
        dudp = -1._r64
        umin = 0._r64
        umax = xmax

        ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --
        ! ---------------------------------------------- loop over ray parameter ---------------------------------------------------

        ! here we explore all rays from "pmin" to "pmax"
        raypar: DO i = 1, 10000

          pplast = pp
          pp = pp + dp

          ! IF upgoing, pmax>pmin, IF downgoing, pmax < pmin
          IF (upgoin .and. (pp .gt. pmax)) pp = pmax
          IF (.not.upgoin .and. (pp .lt. pmax)) pp = pmax

          IF (.not.upgoin .and. (pp .gt. pmin)) pp = pmin
          dp = pp - pplast

          ! calculate x(p) and tau(p)
          CALL pxint2d(nmod, pp, tau, xx, r, domda, vel, icn, nn, nps, itypef, itypeo, icns, icnr, ok)

          IF (ok .ne. 0) RETURN

          u = xx

          IF (i .eq. 1) GO TO 21

          ustep = u - ulast
          dudp = ustep / dp

          ! IF we have been inside (umin,umax) for 2 steps, or IF we have just stepped out of
          ! the interval, turn on step length damping and variable step length
          IF (bothin(u, ulast, umin, umax)) GO TO 26
          IF (stpout(u, ulast, umin, umax)) GO TO 26

          ! IF we are outside the interval don't DO anything
          IF (.not.stpin(u, ulast, umin, umax)) THEN
            dp = dpinit
            GO TO 21
          ENDIF

          ! we have just stepped into the interval. IF our first step into the interval has
          ! left us more than dumin from the boundary, shorten the step size and repeat the step
          du = MIN (ABS(u - umin), ABS(u - umax))
          IF (du .lt. dumin) GO TO 21
          itry = 0
          GO TO 27

          ! IF ustep too large, set dp smaller and recompute step. try this tactic a max of 9
          ! times to get ustep to a reasonable size
          26 IF (ABS(ustep) .lt. dumin) GO TO 23
          27 itry = itry + 1

          ! too many iterations, give up
          IF (itry .ge. 100) THEN
            IF (upgoin .and. ((pplast + dp) .gt. pmax))  CYCLE ray
            dp = dpinit
            GO TO 21
          ENDIF

          ! try another iteration
          pp = pplast
          dp = dp / 2._r64

          CYCLE raypar

          !  adjust ray parameter step for next step depending on du/dp
          23 itry = 0
          IF (dudp .ne. 0.) THEN
            IF (upgoin)      dp =  ABS(dumin / dudp)
            IF (.not.upgoin) dp = -ABS(dumin / dudp)
          ELSE
            dp = dpinit
          ENDIF

          21 CONTINUE

          ! save this point
          IF (itypeo .eq. 1) velo = alfo
          IF (itypeo .eq. 2) velo = beto
          IF (itypef .eq. 1) velf = alff
          IF (itypef .eq. 2) velf = betf
          tho = ASIN(pp * velo)
          thf = ASIN(pp * velf)
          IF (.not.upgoin) thf = PI - thf

          ! correct sign of domda for upgoing or downgoing
          sdudp = dudp
          IF (.not.upgoin) THEN
            domda = -domda
            sdudp = -sdudp
          ENDIF

          ! determine which sheet we're on
          IF (upgoin) isheet = 1
          IF (.not.upgoin) THEN
            IF (sdomda*sdudp .lt. 0._r64) THEN
              sdomda = -sdomda
              isheet = isheet + 1
            ENDIF
          ENDIF

          ! q factor
          sygn = 1._r64
          IF (domda .lt. 0._r64) sygn = -1._r64

          ! qp contains spreading factor and ratio of impendences at ends of the ray
          qp = sygn * SQRT( ABS(domda) / (rhoo*rhof*velf**3 * velo)) /  (4._r64 * PI * velf)

          IF (spconv) THEN
            fr = fr * supd
            fth = fth * supd
          ENDIF

          travtm = tau + pp * xx
          cvel = 1._r64 / MAX(pp, 1.d-5)
          pn = pp
          IF (.not. upgoin) pn = 2. * psourc - pp

          ulast = u
          itry = 0

          ! IF we are close to source keep max x step small
          IF (xx .lt. 5._r64) dumin = dmax / 4._r64
          IF ( (xx .lt. 10._r64) .and. (xx .ge. 5._r64)) dumin = dmax / 2._r64
          IF (xx .ge. 10._r64) dumin = dmax

          ! jump out of loop IF at max ray parameter
          IF (pp .eq. pmax) CYCLE ray
          IF ((jray .gt. 1) .and. (i .eq. 1)) CYCLE raypar

          nrays = nrays + 1           !< add values for current ray

          vpn(nrays) = pn
          vxx(nrays) = xx
          vtt(nrays) = travtm
          vqp(nrays) = qp
          vr(nrays)  = r

#ifdef DEBUG
          WRITE(lu, 105) pn, xx, r, travtm, qp, isheet, itypef, itypeo
          105      FORMAT (5G15.8, I5, 1X, I1, 1X, I1, 1X, 4E15.8)
#endif

          sheets = MAX(sheets, isheet)       !< store highest number of sheets amongst wavetype

          IF (itypef .ne. itypeo) THEN
            CALL report_error('spred - ERROR: wave conversion detected')
            ok = 1
          ENDIF

        ENDDO raypar
        ! ----------------end of ray parameter loop------------------------

      ENDDO ray
      ! ----------------------end of loop over rays------------------------

      IF (nrays .le. 1) THEN
        CALL report_error('spred - ERROR: no rays were found. Check velocity model, source and receivers location')
        ok = 1
        RETURN
      ENDIF

      ALLOCATE(wkbj(src, wtp)%pn(nrays), wkbj(src, wtp)%x(nrays), wkbj(src, wtp)%t(nrays), wkbj(src, wtp)%q(nrays))
      ALLOCATE(wkbj(src, wtp)%r(nrays))

      DO i = 1, nrays
        wkbj(src, wtp)%pn(i) = vpn(i)
        wkbj(src, wtp)%x(i)  = vxx(i)
        wkbj(src, wtp)%t(i)  = vtt(i)
        wkbj(src, wtp)%q(i)  = vqp(i)
        wkbj(src, wtp)%r(i)  = vr(i)
      ENDDO


    END SUBROUTINE spred

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE inser2d(nmod, z, ir, vel, rho, dadz, dbdz, drdz)

      ! Purpose:
      !   to insert an interface at depth "z" in the velocity model "vel". Velocities are continuous across the interface and the
      !   velocities are determined by linear interpolation. This subroutine was originally written by P.Spudich and modified by
      !   W.Imperatori.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   08/03/21                  original version
      !

      INTEGER(i32),                      INTENT(INOUT) :: nmod
      REAL(r64),                         INTENT(IN)    :: z
      INTEGER(i32),                      INTENT(IN)    :: ir
      REAL(r64),    DIMENSION(6,nmod+2), INTENT(INOUT) :: vel
      REAL(r64),    DIMENSION(nmod+2),   INTENT(INOUT) :: rho
      REAL(r64),                         INTENT(IN)    :: dadz, dbdz, drdz
      INTEGER(i32)                                     :: i, j, k

      !-----------------------------------------------------------------------------------------------------------------------------

      ! open up two slots in the vel array for insertion of interfaces
      DO i = ir + 1, nmod
        j = nmod - i + ir + 1
        rho(j + 2) = rho(j)
        DO k = 1, 3
          vel(k, j + 2) = vel(k, j)
        ENDDO
      ENDDO

      ! insert interface parameters
      vel(1, ir + 1) = z
      vel(1, ir + 2) = z
      vel(2, ir + 1) = vel(2, ir) + (z - vel(1, ir)) * dadz
      vel(2, ir + 2) = vel(2, ir + 1)
      vel(3, ir + 1) = vel(3, ir) + (z - vel(1, ir)) * dbdz
      vel(3, ir + 2) = vel(3, ir + 1)
      rho(ir + 1)    = rho(ir) + (z - vel(1, ir)) * drdz
      rho(ir + 2)    = rho(ir + 1)

      ! boost count of model points
      nmod = nmod + 2

    END SUBROUTINE inser2d

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE pxint2d(nmod, p, tau, x, r, domda, vel, icn, lln, nps, itypes, ityper, icns, icnr, ierr)

      ! Purpose:
      !   to compute delay time "tau", distance "x" and spreading factor "domda" (D(omega)/D(area)) for a geometric ray. This
      !   subroutine was originally  part of the program WKBJ3 written by S.K. Dey-Sarkar and C.H. Chapman as described in the paper
      !   "A simple method for the computation of body-wave seismograms", BSSA, Vol. 68, pp. 1577-1593, 1978. This subroutine was
      !   modified by P. Spudich and W. Imperatori.
      !
      !   p    -  ray parameter
      !  vel   -  array containing velocity structure
      !            vel(1,j) = depth of jth model point
      !             " (2,j) = 1./vp "    "    "     "
      !             " (3,j) = 1./vs "    "    "     "
      !             " (4,j) = density of "    "     "
      !             " (5,j) = (vp(j+1) - vp(j)) / (z(j+1) - z(j)) - p vel grad
      !             " (6,j) = (vs(j+1) - vs(j)) / (z(j+1) - z(j)) - s vel grad
      ! nmod   -  number of points in vel, i.e. vel(i,j) for
      !            1 .le. j .le. nmod defined
      !  icn   -  integer array containing indices of model points at top of
      !            layers. e.g. vel(1,icn(2)) is the depth to the top of
      !            layer 2
      !  lln   -  largest layer # in ray description
      !  nps   -  array describing number of p and s ray segments in layer.
      !            nps(i,j), i=layer #, j=1 is p segs, j=2 is s segs.
      ! itypes -  =1 if wave is in p mode at the source
      !           =2 if wave is in s mode at the source
      ! ityper -  =1 if wave is in p mode at the receiver
      !           =2 if wave is in s mode at the receiver
      !  icns  -  index of model point of top of the layer bounded above by
      !            the source (e.g. depth of source =vel(1,icns-1)=vel(1,icns)
      !  icnr  -  index of model point of top of the layer bounded above by
      !            the receiver.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   08/03/21                  original version
      !

      INTEGER(i32),                    INTENT(IN)  :: nmod
      REAL(r64),                       INTENT(IN)  :: p
      REAL(r64),                       INTENT(OUT) :: tau, x, domda, r
      REAL(r64),    DIMENSION(6,nmod), INTENT(IN)  :: vel
      INTEGER(i32),                    INTENT(IN)  :: lln, itypes, ityper, icnr, icns
      INTEGER(i32),                    INTENT(OUT) :: ierr
      INTEGER(i32), DIMENSION(20,2),   INTENT(IN)  :: nps
      INTEGER(i32), DIMENSION(lln),    INTENT(IN)  :: icn
      INTEGER(i32)                                 :: is, ig, ips, nsegs, iz, layer
      LOGICAL                                      :: tp
      REAL(r64)                                    :: sum, dx, dtau, dsum, factor, cs, dz, fact, dadom, dl, dr, rad
      REAL(r64)                                    :: sths, cths, cthr, cthl, cthg, slowg, sthl

      !-----------------------------------------------------------------------------------------------------------------------------

      sum = 0._r64; tau = 0._r64; x = 0._r64; domda = 0._r64
      dx = 0._r64; dtau = 0._r64; dsum = 0._r64
      ierr = 0
      factor = 0._r64
      dl = 0._r64

      r = 0._r64
      dr = 0._r64

      !...  cosines and sines of incidence angles at source and receiver
      cs = 1._r64 /  vel(itypes + 1, icns)

      IF (p*cs .gt. 1._r64) THEN
        ierr = 4
        CALL report_error('pxint2d - ERROR: p greater than slowness at source')
        RETURN
      ENDIF

      sths = p * cs
      cths = SQRT(1._r64 - sths**2)
      cthr = SQRT(1._r64 - (p / vel(ityper + 1, icnr))**2 )

      !...  note, we use abs below because we don't want sign of domda to change from downgoing
      !     to upgoing rays
      IF (sths .ne. 0.) factor = cths * cthr / sths**2

      !...  loop over wtps (p and s)
      wtp: DO ips = 1, 2

        tp = .false.
        is = ips + 1
        ig = ips + 4
        layer = 1

        IF (p .gt. vel(is, 1)) tp = .true.

        !  loop over model points from top to bottom doing x, tau, and domda integrals
        mpt: DO iz = 1, nmod-1

          ! thickness of current layer
          dz = vel(1, iz + 1) - vel(1, iz)

          IF (dz .eq. 0.) THEN
            layer = layer + 1
            IF (layer .gt. lln) CYCLE wtp
            CYCLE mpt
          ENDIF

          nsegs = nps(layer,ips)
          fact = factor

          !...  DO the following for waves IF they exist between model points iz and iz+1
          IF ( (nsegs .ne. 0) .and. .not.tp ) THEN
            sthl = p / vel(is, iz)
            !...     check for turning point above this region
            IF (sthl .gt. 1._r64) THEN
              tp = .true.
              CYCLE wtp
            ENDIF

            cthl = SQRT(1._r64 - (p/vel(is, iz))**2)

            !...     check for zero velocity gradient
            IF (vel(ig, iz) .eq. 0._r64) THEN
              !...        zero gradient. is the ray horizontal in this zone?
              IF (p .eq. vel(is, iz)) THEN
                !...           ray is horizontal
                dx = 10000._r64
                dtau = 0._r64
                ierr = 3
                dsum = 0._r64
              ELSE
                !...           ray passes through zone
                dx = dz * sthl / cthl
                dtau = dz * SQRT(vel(is, iz)**2 - p**2)
                cthg = cthl
                dsum = fact * dx / cthl / cthg
              ENDIF

              dr = SQRT(dz**2 + dx**2)

            ELSE
              !...        non-zero gradient. check for turning point between model points iz and iz+1
              !           calculate dx and dtau

              IF (p .ge. vel(is, iz + 1)) THEN
                !...           turning point exists
                tp = .true.
                IF (p .eq. 0.) THEN
                  ierr = 6
                  RETURN
                ENDIF

                dx = cthl / p / vel(ig,iz)
                dtau = (-LOG(p) - (cthl - LOG(vel(is, iz) + vel(is, iz) * cthl))) / vel(ig, iz)

                ! dr = (1._r64 / p - 1._r64 / vel(is,iz)) / vel(ig,iz)
                ! dr = nsegs * SQRT(dr**2 + dx**2)

                rad = 1._r64 / (p * vel(ig, iz))
                dr = ACOS(1._r64 - 2._r64 * dx**2 / rad**2) * rad / 2._r64

                !...           calculate spreading factor stuff.
                IF (cths .ne. 0.) THEN
                  !...              turning point not at source
                  dsum = -fact * dx / cthl**2
                ELSE
                  !...              turning point at source
                  fact = cthr / sths**2
                  dsum = fact * dx / cthl
                ENDIF

              ELSE

                !...           no turning point
                cthg = SQRT(1._r64 - (p / vel(is, iz + 1))**2)
                dx = 0._r64
                IF (p .ne. 0._r64) dx = -cthg / p / vel(ig, iz)
                slowg = vel(is, iz + 1)
                dtau = (cthg - LOG(slowg + slowg*cthg)) / vel(ig, iz)

                !...           add lower limit of integral
                IF (p .eq. 0._r64) THEN
                  dx = 0._r64
                  dr = dz
                ELSE
                  dx = dx + cthl / p / vel(ig,iz)

                  rad = 1._r64 / (p * vel(ig, iz))
                  dr  = SQRT(dz**2 + dx**2)
                  dr  = ACOS(1._r64 - 0.5_r64 * dr**2 / rad**2) * rad

                ENDIF
                dtau = dtau - (cthl - LOG(vel(is, iz) + vel(is, iz)*cthl)) / vel(ig, iz)
                dsum = fact * dx / cthl / cthg

                ! dr = SQRT(dz**2 + dx**2)

              ENDIF

            ENDIF

            dl = dl + nsegs * dz
            x = x + nsegs * dx
            tau = tau + nsegs * dtau
            sum = sum + nsegs * dsum

            ! r = r + dr
            r = r + dr * nsegs

          ENDIF

        ENDDO mpt
      ENDDO wtp

      dadom = x * sum

      !...  dadom for a vertical ray.  may not be correct for a complicated ray with mode
      !     conversions

      IF (sths .eq. 0._r64) dadom = dl**2
      IF (dadom .eq. 0._r64) THEN
        ierr = 5
      ELSE
        domda = 1._r64 / dadom
      ENDIF


      SELECT CASE (ierr)
        CASE(2)
          CALL report_error('pxint2d - ERROR: p larger than local slowness')
        CASE(3)
          CALL report_error('pxint2d - ERROR: spreading requested for horizontal ray in uniform medium')
        CASE(4)
          CALL report_error('pxint2d - ERROR: p greater than slowness at source')
        CASE(5)
          CALL report_error('pxint2d - ERROR: dadom=0 for some reason')
        CASE(6)
          CALL report_error('pxint2d - ERROR: vel(is, iz+1) .le. 0')
      END SELECT

    END SUBROUTINE pxint2d

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    LOGICAL FUNCTION bothin(u, ulast, a, b)

      REAL(r64), INTENT(IN) :: u, ulast, a, b

      !-----------------------------------------------------------------------------------------------------------------------------

      bothin = inside(u,a,b) .and. inside(ulast,a,b)

    END FUNCTION bothin

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    LOGICAL FUNCTION stpin(u, ulast, a, b)

      REAL(r64), INTENT(IN) :: u, ulast, a, b

      !-----------------------------------------------------------------------------------------------------------------------------

      stpin = (.not.inside(ulast,a,b)) .and. inside(u,a,b)

    END FUNCTION stpin

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    LOGICAL FUNCTION stpout(u, ulast, a, b)

      REAL(r64), INTENT(IN) :: u, ulast, a, b

      !-----------------------------------------------------------------------------------------------------------------------------

      stpout = inside(ulast,a,b) .and. .not.inside(u,a,b)

    END FUNCTION stpout

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    LOGICAL FUNCTION inside(u, a, b)

      REAL(r64), INTENT(IN) :: u, a, b

      !-----------------------------------------------------------------------------------------------------------------------------

      inside  = (u .ge. a) .and. (u .le. b)

    END FUNCTION inside

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE interpray(is, sr, src, wtp, po, ro, to, qo)

      INTEGER(i32),  INTENT(INOUT) :: is
      REAL(r32),     INTENT(IN)    :: sr
      INTEGER(i32),  INTENT(IN)    :: src, wtp
      REAL(r32),     INTENT(OUT)   :: po, ro, to, qo
      INTEGER(i32)                 :: i
      REAL(r64)                    :: delta

      !-----------------------------------------------------------------------------------------------------------------------------

      ASSOCIATE(p => wkbj(src, wtp)%pn, q => wkbj(src, wtp)%q, t => wkbj(src, wtp)%t, x => wkbj(src, wtp)%x, r => wkbj(src, wtp)%r)

        po = 0._r32; ro = 0._r32; to = 0._r32; qo = 0._r32

        DO i = is, SIZE(p) - 1

          IF ( ((sr .ge. x(i)) .and. (sr .lt. x(i + 1))) .or. ((sr .le. x(i)) .and. (sr .gt. x(i + 1))) ) THEN

            delta = (sr - x(i)) / (x(i + 1) - x(i))

            po = p(i) + (p(i + 1) - p(i)) * delta
            ! xo = x(i) + (x(i + 1) - x(i)) * delta
            ro = r(i) + (r(i + 1) - r(i)) * delta
            to = t(i) + (t(i + 1) - t(i)) * delta
            qo = q(i) + (q(i + 1) - q(i)) * delta

            is = i + 1       !< update index for next call in order to skip current branch

            EXIT

          ENDIF

        ENDDO

        is = MIN(is, SIZE(p))

      END ASSOCIATE

    END SUBROUTINE interpray

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *


    ! SUBROUTINE interpray(src, wtp, sheet, sr, pl, ql, tl, rl)
    ! !
    ! ! DESCRIPTION:
    ! !
    ! !   INTERPOLATE RAY-RELATED QUANTITIES (P,Q,T,R) AT SPECIFIC SOURCE-RECEIVER HORIZONTAL
    ! !   DISTANCE SR FOR GIVEN SOURCE (SRC), WAVE TYPE (WTP) AND SHEET (SHEET). VALUES ARE
    ! !   RETURNED IN PL, QL, TL, RL
    ! !
    ! ! AUTHOR: W.IMPERATORI (SED@ETH)
    ! !
    ! ! HISTORY: CREATED ON MAY 2017 (V1.0)
    !
    !   ! INDEX FOR SOURCE, WAVE TYPE AND SHEET
    !   INTEGER(i32),INTENT(IN) :: src, wtp, sheet
    !   ! HORIZONTAL SOURCE-RECEIVER DISTANCE
    !   REAL(r32),INTENT(IN)    :: sr
    !   REAL(r32),INTENT(OUT)   :: pl, ql, tl, rl
    !   INTEGER(i32)            :: nf, k
    !   REAL(r32)               :: delta, dp
    !
    !   !-----------------------------------------------------------------------------------------------------------------------------
    !
    !   ASSOCIATE(p => wkbj(src, wtp)%pn, q => wkbj(src, wtp)%q, t => wkbj(src, wtp)%t, x => wkbj(src, wtp)%x, r => wkbj(src, wtp)%r)
    !
    !     nf = 0
    !
    !     pl = 0._r32
    !     ql = 0._r32
    !     tl = 0._r32
    !     rl = 0._r32
    !
    !     IF (sr .gt. maxval(x)) THEN
    !       print*,'sr > x ',sr,maxval(x), src, wtp, sheet
    !       stop
    !     ENDIF
    !
    !     ! cycle over rays
    !     DO k = 1, SIZE(p) - 1
    !
    !       IF (x(k) .eq. x(k+1)) CYCLE
    !
    !       IF ( ((sr .ge. x(k)) .and. (sr .lt. x(k+1))) .or. ((sr .ge. x(k+1)) .and. (sr .lt. x(k))) ) THEN
    !
    !         ! update sheet number
    !         nf = nf + 1
    !
    !         ! are we in the right sheet?
    !         IF (nf .ne. sheet) CYCLE
    !
    !         delta = (sr - x(k)) / (x(k + 1) - x(k))
    !
    !         ! interpolate values at current sr distance
    !         pl = p(k) + (p(k + 1) - p(k)) * delta
    !         ql = q(k) + (q(k + 1) - q(k)) * delta
    !         rl = r(k) + (r(k + 1) - r(k)) * delta
    !         tl = t(k) + (t(k + 1) - t(k)) * delta
    !
    !         ! check if there is a problem with ray param. to be remove once
    !         ! triplicatins are tested
    !         dp = p(k + 1) - p(k)
    !
    !         IF (dp .eq. 0.) THEN
    !           print*,'inside process ray: ',p(k+1),p(k),x(k+1),x(k)
    !           !call exe_stop('iso.f90','ghost_fault',4)
    !         ENDIF
    !
    !       ENDIF
    !
    !     ENDDO ! end loop over rays
    !
    !     ! handle case where rays do not span distance range (i.e. r .gt. max(x)). this should not
    !     ! happen once the code is ready
    !     IF (nf .eq. 0) THEN
    !
    !       nf = 1
    !       k = SIZE(p) - 1
    !
    !       delta = (sr - x(k)) / (x(k + 1) - x(k))
    !
    !       pl = p(k) + (p(k + 1) - p(k)) * delta
    !       ql = q(k) + (q(k + 1) - q(k)) * delta
    !       rl = r(k) + (r(k + 1) - r(k)) * delta
    !       tl = t(k) + (t(k + 1) - t(k)) * delta
    !
    !       dp = p(k + 1) - p(k)
    !
    !       ! check if there is a problem with ray param. to be remove once
    !       ! triplicatins are tested
    !       dp = p(k + 1) - p(k)
    !
    !       IF (dp .eq. 0.) THEN
    !         print*,'inside process ray: ',p(k+1),p(k),x(k+1),x(k)
    !         !call exe_stop('iso.f90','ghost_fault',4)
    !       ENDIF
    !
    !     ENDIF
    !
    !   END ASSOCIATE
    !
    ! END SUBROUTINE interpray

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *



END MODULE m_wkbj
