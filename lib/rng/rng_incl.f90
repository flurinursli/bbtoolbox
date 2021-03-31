INTEGER(i32),                OPTIONAL, INTENT(IN) :: nt
INTEGER(i32)                                      :: threads
INTEGER(c_i32)                                    :: c_npts, i0, i1
REAL(r__),     DIMENSION(:), POINTER              :: x_flattened

!-----------------------------------------------------------------------------------------------------------------------------------

! set the number of threads to desired value
!$ IF (PRESENT(nt)) THEN
!$    threads = omp_get_max_threads()
!$    CALL omp_set_num_threads(nt)
!$ ENDIF

c_npts = SIZE(x)

i0 = 1
i1 = c_npts

x_flattened(1:c_npts) => x

! parallel region below is meant to parallelize the calculation of random numbers when "rng" is called from a serial region. In this
! case, "i0", "i1", "c_npts" and "c_skip" are modified and "setup_trng" called to update number of random deviates to be skipped.
! Note that in such configuration leapfrogging is disabled.
!$omp parallel default(shared) private(i0, i1) firstprivate(c_npts) copyin(c_skip)
!$ c_skip = c_skip + omp_get_thread_num() * (c_npts / omp_get_num_threads())
!$ IF (omp_get_num_threads() .gt. 1) CALL setup_trng(c_seed, c_skip, 1, 0)
!$ i0 = omp_get_thread_num() * c_npts / omp_get_num_threads() + 1
!$ i1 = (omp_get_thread_num() + 1) * c_npts / omp_get_num_threads()
!$ c_npts = i1 - i0 + 1
CALL fun(c_par1, c_par2, c_npts, x_flattened(i0:i1))
!$omp end parallel

! increase amount of random numbers to be skipped
c_skip = c_skip + SIZE(x)

NULLIFY(x_flattened)

! restore initial number of threads
!$ IF (PRESENT(nt)) CALL omp_set_num_threads(threads)
