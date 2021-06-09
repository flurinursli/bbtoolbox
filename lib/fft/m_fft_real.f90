MODULE m_fft_real

  ! Purpose:
  !   To provide routines to compute forward and inverse Fourier transform of real vectors and matrices by taking advantage of
  !   Hermitian symmetry. In other words, for N (NxM) real input data, one obtains N/2 + 1 (N/2 + 1, M) complex values in the
  !   transformed domain.
  !   By default, all transforms are in-place and data are copied to memory-aligned arrays for better performance and to avoid side
  !   effects (i.e. data overwritten). However, copying large arrays (especially in double-precision) may hide any performance gain
  !   associated to SIMD instructions exploitation.
  !
  !   The calling sequence is as follows:
  !     a) subroutine "make_fftw_plan" (called only once)
  !     b) subroutine "fft"/"ifft" (called as many times as necessary, as long as the array size doesn't change)
  !     c) subroutine "destroy_fftw_plan" (called only once)
  !
  !   Parallel execution can be achieved in two ways:
  !      i) for multithreaded 2D transforms, compile with -DFFTW_OMP and indicate the number of threads when creating the plan.
  !     ii) for concurrent 1D transforms of an array columns/rows, define an OPENMP parallel loop where each thread operates on
  !         different columns/rows.
  !
  !   When the library is compiled in double-precision (-DDOUBLE_PREC), both single- and double-precision input data are accepted
  !   (but transforms are always in double-precision). On the other hand, if the library is compiled in single-precision, only
  !   single-precision input data can be provided.
  !
  !   Routines rely on the FFTW library (version 3.x)
  !
  ! Revisions:
  !     Date                    Description of change
  !     ====                    =====================
  !   02/09/20                  original version
  !

  !$ USE             :: omp_lib
  USE, INTRINSIC     :: iso_c_binding
  USE, NON_INTRINSIC :: m_precisions

  IMPLICIT none

  INCLUDE 'fftw3.f03'

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  PRIVATE

  PUBLIC :: make_fftw_plan, destroy_fftw_plan, fft, ifft

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  INTERFACE fft
    MODULE PROCEDURE fft_1d_r32, fft_2d_r32
#ifdef DOUBLE_PREC
    MODULE PROCEDURE fft_1d_r64, fft_2d_r64
#endif
  END INTERFACE fft

  INTERFACE ifft
    MODULE PROCEDURE ifft_1d_r32, ifft_2d_r32
#ifdef DOUBLE_PREC
    MODULE PROCEDURE ifft_1d_r64, ifft_2d_r64
#endif
  END INTERFACE ifft

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  ! define FFTW plan
  INTEGER(i32)                            :: plan_type = FFTW_ESTIMATE

#ifdef DOUBLE_PREC
  INTEGER, PARAMETER :: r__   = r64
  INTEGER, PARAMETER :: c_r__ = c_r64
  INTEGER, PARAMETER :: c_c__ = c_c64
#else
  INTEGER, PARAMETER :: r__   = r32
  INTEGER, PARAMETER :: c_r__ = c_r32
  INTEGER, PARAMETER :: c_c__ = c_c32
#endif

  ! FFTW variables
  COMPLEX(c_c__), DIMENSION(:),   POINTER :: tu => NULL()
  COMPLEX(c_c__), DIMENSION(:,:), POINTER :: tv => NULL()
  REAL(c_r__),    DIMENSION(:),   POINTER :: u  => NULL()
  REAL(c_r__),    DIMENSION(:,:), POINTER :: v  => NULL()

  TYPE(c_ptr)                             :: pu, fplan, iplan

  !$omp threadprivate (u, tu, pu, fplan, iplan)

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  CONTAINS

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE make_fftw_plan(n, nt)

      INTEGER(i32), DIMENSION(:),           INTENT(IN) :: n               !< points to be transformed along each dimension
      INTEGER(i32),               OPTIONAL, INTENT(IN) :: nt              !< multithreaded transform (2D arrays only)
      INTEGER(i32)                                     :: npts, status

      !-----------------------------------------------------------------------------------------------------------------------------

      ! allow implicit (FFTW internal) parallelization only when transforming matrices
#ifdef FFTW_OMP
      IF (PRESENT(nt) .and. (SIZE(n) .eq. 2)) THEN
#ifdef DOUBLE_PREC
        status = fftw_init_threads()
        CALL fftw_plan_with_nthreads(nt)
#else
        status = fftwf_init_threads()
        CALL fftwf_plan_with_nthreads(nt)
#endif
      ENDIF
#endif

      ! points in transformed domain
      npts = n(1) / 2 + 1

      ! allow user OPENMP parallelization to transform vectors
      IF (SIZE(n) .eq. 1) THEN

        !$omp parallel default(shared)

        ! allocate contiguous memory
#ifdef DOUBLE_PREC
        pu = fftw_alloc_complex(INT(npts, c_size_t))
#else
        pu = fftwf_alloc_complex(INT(npts, c_size_t))
#endif

        ! map C-to-Fortran pointers
        CALL c_f_pointer(pu, u, [2*npts])
        CALL c_f_pointer(pu, tu, [npts])

        ! ALLOCATE(u(2*npts), tu(npts))

        !$omp critical (make_plan)
#ifdef DOUBLE_PREC
        fplan = fftw_plan_dft_r2c_1d(n(1), u, tu, plan_type)
        iplan = fftw_plan_dft_c2r_1d(n(1), tu, u, plan_type)
#else
        fplan = fftwf_plan_dft_r2c_1d(n(1), u, tu, plan_type)
        iplan = fftwf_plan_dft_c2r_1d(n(1), tu, u, plan_type)
#endif
        !$omp end critical (make_plan)

        ! DEALLOCATE(u, tu)

        !$omp end parallel

      ELSE

        ! allocate contiguous memory
#ifdef DOUBLE_PREC
        pu = fftw_alloc_complex(INT(npts * n(2), c_size_t))
#else
        pu = fftwf_alloc_complex(INT(npts * n(2), c_size_t))
#endif

        ! map C-to-Fortran pointers
        CALL c_f_pointer(pu, v, [2*npts, n(2)])
        CALL c_f_pointer(pu, tv, [npts, n(2)])

        ! ALLOCATE(v(2*npts, n(2)), tv(npts, n(2)))

#ifdef DOUBLE_PREC
        fplan = fftw_plan_dft_r2c_2d(n(2), n(1), v, tv, plan_type)
        iplan = fftw_plan_dft_c2r_2d(n(2), n(1), tv, v, plan_type)
#else
        fplan = fftwf_plan_dft_r2c_2d(n(2), n(1), v, tv, plan_type)
        iplan = fftwf_plan_dft_c2r_2d(n(2), n(1), tv, v, plan_type)
#endif
        ! DEALLOCATE(v, tv)

      ENDIF

    END SUBROUTINE make_fftw_plan

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE destroy_fftw_plan(n, nt)

      INTEGER(i32), DIMENSION(:),           INTENT(IN) :: n               !< points to be transformed along each dimension
      INTEGER(i32),               OPTIONAL, INTENT(IN) :: nt              !< multithreaded transform (2D arrays only)

      !-----------------------------------------------------------------------------------------------------------------------------

      IF (SIZE(n) .eq. 1) THEN

        NULLIFY(u, tu)

        !$omp parallel default(shared)

        !$omp critical (destroy_plan)
#ifdef DOUBLE_PREC
        CALL fftw_destroy_plan(fplan)
        CALL fftw_destroy_plan(iplan)

        CALL fftw_free(pu)
#else
        CALL fftwf_destroy_plan(fplan)
        CALL fftwf_destroy_plan(iplan)

        CALL fftwf_free(pu)
#endif
        !$omp end critical (destroy_plan)

        !$omp end parallel

      ELSE

        NULLIFY(v, tv)

#ifdef DOUBLE_PREC
        CALL fftw_destroy_plan(fplan)
        CALL fftw_destroy_plan(iplan)

        CALL fftw_free(pu)

#ifdef FFTW_OMP
        IF (PRESENT(nt)) CALL fftw_cleanup_threads()
#endif
#else
        CALL fftwf_destroy_plan(fplan)
        CALL fftwf_destroy_plan(iplan)

        CALL fftwf_free(pu)

#ifdef FFTW_OMP
        IF (PRESENT(nt)) CALL fftwf_cleanup_threads()
#endif
#endif
      ENDIF

    END SUBROUTINE destroy_fftw_plan

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE fft_1d_r32(x, z)

      REAL(r32),    DIMENSION(:),             INTENT(IN)  :: x
      COMPLEX(r32), DIMENSION(SIZE(x)/2 + 1), INTENT(OUT) :: z
      INTEGER(i32)                                        :: i

      !-----------------------------------------------------------------------------------------------------------------------------

      u(:) = 0._r__

      DO i = 1, SIZE(x)
        u(i) = x(i)
      ENDDO

#ifdef DOUBLE_PREC
      CALL fftw_execute_dft_r2c(fplan, u, tu)
#else
      CALL fftwf_execute_dft_r2c(fplan, u, tu)
#endif

      DO i = 1, SIZE(z)
        z(i) = tu(i)
      ENDDO

    END SUBROUTINE fft_1d_r32

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE ifft_1d_r32(x, z)

      REAL(r32),    DIMENSION(:),             INTENT(OUT) :: x
      COMPLEX(r32), DIMENSION(SIZE(x)/2 + 1), INTENT(IN)  :: z
      INTEGER(i32)                                        :: n, i

      !-----------------------------------------------------------------------------------------------------------------------------

      DO i = 1, SIZE(z)
        tu(i) = z(i)
      ENDDO

#ifdef DOUBLE_PREC
      CALL fftw_execute_dft_c2r(iplan, tu, u)
#else
      CALL fftwf_execute_dft_c2r(iplan, tu, u)
#endif

      n = SIZE(x)

      DO i = 1, n
        x(i) = u(i) / n
      ENDDO

    END SUBROUTINE ifft_1d_r32

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE fft_2d_r32(x, z)

      REAL(r32),    DIMENSION(:,:),                        INTENT(IN)  :: x
      COMPLEX(r32), DIMENSION(SIZE(x,1)/2 + 1, SIZE(x,2)), INTENT(OUT) :: z
      INTEGER(i32)                                                     :: i, j

      !-----------------------------------------------------------------------------------------------------------------------------

      v(:,:) = 0._r__

      DO j = 1, SIZE(x, 2)
        DO i = 1, SIZE(x, 1)
          v(i, j) = x(i, j)
        ENDDO
      ENDDO

#ifdef DOUBLE_PREC
      CALL fftw_execute_dft_r2c(fplan, v, tv)
#else
      CALL fftwf_execute_dft_r2c(fplan, v, tv)
#endif

      DO j = 1, SIZE(z, 2)
        DO i = 1, SIZE(z, 1)
          z(i, j) = tv(i, j)
        ENDDO
      ENDDO

    END SUBROUTINE fft_2d_r32

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE ifft_2d_r32(x, z)

      REAL(r32),    DIMENSION(:,:),                        INTENT(OUT) :: x
      COMPLEX(r32), DIMENSION(SIZE(x,1)/2 + 1, SIZE(x,2)), INTENT(IN)  :: z
      INTEGER(i32)                                                     :: n, i, j

      !-----------------------------------------------------------------------------------------------------------------------------

      DO j = 1, SIZE(z, 2)
        DO i = 1, SIZE(z, 1)
          tv(i, j) = z(i, j)
        ENDDO
      ENDDO

#ifdef DOUBLE_PREC
      CALL fftw_execute_dft_c2r(iplan, tv, v)
#else
      CALL fftwf_execute_dft_c2r(iplan, tv, v)
#endif

      n = SIZE(x)

      DO j = 1, SIZE(x, 2)
        DO i = 1, SIZE(x, 1)
          x(i, j) = v(i, j) / n
        ENDDO
      ENDDO

    END SUBROUTINE ifft_2d_r32

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

#ifdef DOUBLE_PREC

    SUBROUTINE fft_1d_r64(x, z)

      REAL(r64),    DIMENSION(:),             INTENT(IN)  :: x
      COMPLEX(r64), DIMENSION(SIZE(x)/2 + 1), INTENT(OUT) :: z
      INTEGER(i32)                                        :: i

      !-----------------------------------------------------------------------------------------------------------------------------

      u(:) = 0._r64

      DO i = 1, SIZE(x)
        u(i) = x(i)
      ENDDO

      CALL fftw_execute_dft_r2c(fplan, u, tu)

      DO i = 1, SIZE(z)
        z(i) = tu(i)
      ENDDO

    END SUBROUTINE fft_1d_r64

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE ifft_1d_r64(x, z)

      REAL(r64),    DIMENSION(:),             INTENT(OUT) :: x
      ! COMPLEX(r64), DIMENSION(SIZE(x)/2 + 1), INTENT(INOUT)  :: z
      COMPLEX(r64), DIMENSION(SIZE(x)/2 + 1), INTENT(IN)  :: z
      INTEGER(i32)                                        :: n, i

      !-----------------------------------------------------------------------------------------------------------------------------

      DO i = 1, SIZE(z)
        tu(i) = z(i)
      ENDDO

      CALL fftw_execute_dft_c2r(iplan, tu, u)

      n = SIZE(x)

      DO i = 1, n
        x(i) = u(i) / n
      ENDDO

      ! CALL fftw_execute_dft_c2r(iplan, z, x)
      !
      ! n = SIZE(x)
      !
      ! DO i = 1, n
      !   x(i) = REAL(x(i), r64) / n
      ! ENDDO

    END SUBROUTINE ifft_1d_r64

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE fft_2d_r64(x, z)

      REAL(r64),    DIMENSION(:,:),                        INTENT(IN)  :: x
      ! REAL(r64),    DIMENSION(:,:),                        INTENT(INOUT)  :: x        !< this is not overwritten
      COMPLEX(r64), DIMENSION(SIZE(x,1)/2 + 1, SIZE(x,2)), INTENT(OUT) :: z
      INTEGER(i32)                                                     :: i, j

      !-----------------------------------------------------------------------------------------------------------------------------

      v(:,:) = 0._r64

      DO j = 1, SIZE(x, 2)
        DO i = 1, SIZE(x, 1)
          v(i, j) = x(i, j)
        ENDDO
      ENDDO

      CALL fftw_execute_dft_r2c(fplan, v, tv)

      DO j = 1, SIZE(z, 2)
        DO i = 1, SIZE(z, 1)
          z(i, j) = tv(i, j)
        ENDDO
      ENDDO

      ! CALL fftw_execute_dft_r2c(fplan, x, z)

    END SUBROUTINE fft_2d_r64

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE ifft_2d_r64(x, z)

      REAL(r64),    DIMENSION(:,:),                        INTENT(OUT) :: x
      COMPLEX(r64), DIMENSION(SIZE(x,1)/2 + 1, SIZE(x,2)), INTENT(IN)  :: z
      ! COMPLEX(r64), DIMENSION(SIZE(x,1)/2 + 1, SIZE(x,2)), INTENT(INOUT)  :: z            !< this is overwritten
      INTEGER(i32)                                                     :: n, i, j

      !-----------------------------------------------------------------------------------------------------------------------------

      DO j = 1, SIZE(z, 2)
        DO i = 1, SIZE(z, 1)
          tv(i, j) = z(i, j)
        ENDDO
      ENDDO

      CALL fftw_execute_dft_c2r(iplan, tv, v)

      n = SIZE(x)

      DO j = 1, SIZE(x, 2)
        DO i = 1, SIZE(x, 1)
          x(i, j) = v(i, j) / n
        ENDDO
      ENDDO

      ! CALL fftw_execute_dft_c2r(iplan, z, x)
      !
      ! n = SIZE(x)
      !
      ! DO j = 1, SIZE(x, 2)
      !   DO i = 1, SIZE(x, 1)
      !     x(i, j) = x(i, j) / n
      !   ENDDO
      ! ENDDO

    END SUBROUTINE ifft_2d_r64

#endif
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

END MODULE m_fft_real
