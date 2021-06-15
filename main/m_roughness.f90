MODULE m_roughness

  USE, NON_INTRINSIC :: m_precisions
  USE, NON_INTRINSIC :: m_logfile
  USE, NON_INTRINSIC :: m_scarflib
  USE, NON_INTRINSIC :: m_source
  USE, NON_INTRINSIC :: m_stat
  USE, NON_INTRINSIC :: m_strings
  USE, NON_INTRINSIC :: m_toolbox, ONLY: input

  IMPLICIT none

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  PRIVATE

  PUBLIC :: fault_roughness
  PUBLIC :: roughness

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  REAL(r32), ALLOCATABLE, DIMENSION(:,:), TARGET :: roughness        !< "target" attribute used only inside this module

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  CONTAINS

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE fault_roughness(ok, ref, pl, iter)

      INTEGER(i32),                              INTENT(OUT) :: ok
      INTEGER(i32),                              INTENT(IN)  :: ref, pl, iter
      CHARACTER(:), ALLOCATABLE                              :: fo
      INTEGER(i32)                                           :: i, j, seed, lu, icr, totnutr
      INTEGER(i32),              DIMENSION(3)                :: iuc, ivc
      INTEGER(i32), ALLOCATABLE, DIMENSION(:)                :: npts
      REAL(r32)                                              :: dh, du, rms, sigma, mu
      REAL(r32),                 DIMENSION(2)                :: nc, fc
      REAL(r32),                 DIMENSION(3)                :: uc, vc
      REAL(r32),                 DIMENSION(8)                :: stats
      REAL(r32),                 DIMENSION(2),   PARAMETER   :: cl = [1._r32, 1._r32]
      REAL(r32),                 DIMENSION(:),   POINTER     :: u1 => NULL(), v1 => NULL(), r1 => NULL()
      REAL(r32),    ALLOCATABLE, DIMENSION(:)                :: avg, var
      REAL(r32),    ALLOCATABLE, DIMENSION(:,:), TARGET      :: u, v

      !-----------------------------------------------------------------------------------------------------------------------------

      ok = 0

      IF (ALLOCATED(roughness)) DEALLOCATE(roughness)

      IF (input%source%add_roughness .and. (input%advanced%verbose .eq. 2) .and. (ref .eq. 1) )  THEN
        ALLOCATE(avg(SIZE(nugr)), var(SIZE(nugr)), npts(SIZE(nugr)))
      ENDIF

      ! set "seed" such that roughness depends on iteration and fault plane
      seed = input%source%seed + (iter - 1) * SIZE(plane) + pl

      ! random field must be always generated on a grid whose extension is given by fault plane and whose resolution equals the
      ! minimum grid-step amongst all refinements
      ! nc = [MIN(umingr, vmingr(1)), MIN(umingr, vmingr(1))]
      nc = [umingr, vmingr(1)]
      ! fc = [MAX(umaxgr, vmaxgr(SIZE(vmaxgr))), MAX(umaxgr, vmaxgr(SIZE(vmaxgr)))]
      fc = [umaxgr, vmaxgr(SIZE(vmaxgr))]
      dh = MIN(MINVAL(dutr), MINVAL(dvtr))

      ALLOCATE(roughness(nugr(ref), nvgr(ref)), u(nugr(ref), nvgr(ref)), v(nugr(ref), nvgr(ref)))

      roughness(:,:) = 0._r32

      IF (input%source%add_roughness) THEN

        ! define points where the random field will be computed
        DO j = 1, nvgr(ref)
          du = (1 - MOD(j, 2)) * dutr(ref) / 2._r32                !< "du" is 0 for odd "j", dutr/2 for even "j"
          DO i = 1, nugr(ref) + MOD(j,2) - 1                       !< "i" tops at "nugr" for odd "j", at "nugr-1" for even "j"
            u(i, j) = umingr + (i - 1) * dutr(ref) + du
            v(i, j) = vmingr(ref) + (j - 1) * dvtr(ref)
          ENDDO
          IF (i .eq. nugr(ref)) THEN
            u(nugr(ref), j) = u(i - 1, j)      !< last u-point for even "j" is simply a copy
            v(nugr(ref), j) = v(i - 1, j)
          ENDIF
        ENDDO

        u1(1:nugr(ref)*nvgr(ref)) => u
        v1(1:nugr(ref)*nvgr(ref)) => v

        r1(1:nugr(ref)*nvgr(ref)) => roughness

        sigma = SUM(plane(:)%length) * 10**input%source%roughness         !< assume zero-mean, such that rms = std.dev.

        ! use user-defined PSD, unstructured grid, set "ds" always larger than "dh" to avoid aliasing
        CALL scarf_initialize(dh, 2, cl, sigma, u1, v1, nc = nc, fc = fc, ds = 2*MAX(dutr(ref), dvtr(ref)), rescale=1)

        CALL scarf_execute(seed, r1, stats)

        CALL scarf_finalize()

        NULLIFY(u1, v1, r1)

        ! rescale roughness such that alpha = 10**input%source%roughness = RMS / LENGTH
        ! rms = SQRT(mean(roughness)**2 + variance(roughness))
        ! rms = SUM(plane(:)%length) * (10**input%source%roughness) / rms
        !
        ! roughness = roughness * rms

        IF (input%advanced%verbose .eq. 2) THEN

          var(ref) = variance(roughness)
          avg(ref) = mean(roughness)
          npts(ref) = nugr(ref) * nvgr(ref)

          IF (ref .eq. SIZE(nugr)) THEN

            !rms = SQRT(mean(roughness)**2 + variance(roughness))
            CALL parallel_variance(var, avg, npts, rms, mu)

            rms = SQRT(mu**2 + rms)

            CALL update_log(num2char('<fault roughness>', justify='c', width=30) +    &
            num2char('Actual RMS', width=15, justify='r') + '|' +  &
            num2char('Expected RMS', width=15, justify='r') + '|')

            CALL update_log(num2char('', width=30, justify='c')  +  &
            num2char(rms, width=15, notation='f', precision=2, justify='r') + '|' + &
            num2char(sigma, width=15, notation='f', precision=2, justify='r') + '|', blankline=.false.)

            DEALLOCATE(var, avg, npts)

          ENDIF
        ENDIF

#ifdef DEBUG

        fo = 'roughness_' + num2char(ref) + '_' + num2char(pl) + '_' + num2char(iter) + '.bin'

        OPEN(newunit = lu, file = fo, status = 'replace', form = 'unformatted', access = 'stream', action = 'write', IOSTAT = ok)

        WRITE(lu, POS=1) nutr(ref)
        WRITE(lu)        nvtr(ref)

        totnutr = 2*nutr(ref) - 1

        DO j = 1, nvtr(ref)
          DO i = 1, totnutr

            CALL cornr(j, i, iuc, ivc)                  !< corner indices for current triangle

            CALL cornr2uv(iuc, ivc, ref, uc, vc)        !< on-fault coordinates

            DO icr = 1, 3
              WRITE(lu) uc(icr), vc(icr), roughness(iuc(icr), ivc(icr))
            ENDDO

          ENDDO
        ENDDO

        WRITE(lu) sigma

        CLOSE(lu, iostat = ok)

        IF (ok .ne. 0) THEN
          CALL report_error('Error while closing file ' + fo)
          RETURN
        ENDIF

#endif

      ENDIF

    END SUBROUTINE fault_roughness

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

END MODULE m_roughness
