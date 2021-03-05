#ifdef DOUBLE_PREC
MODULE m_llsq_r64
#else
MODULE m_llsq_r32
#endif

  ! Purpose:
  !   to compute the minimum norm solution to a linear least squares problem. For instance, given two sets of M observations {X, Y}
  !   and {X, Z} and a linear regression model Y(X) = alpha*X + beta, the problem can be cast in the following matrix form:
  !
  !     A = [x1 1]   B = [y1 z1]   X = [alpha(Y) alpha(Z)]
  !         [x2 1]       [y2 z2]       [beta(Y)  beta(Z)]
  !         [....]       [.....]
  !         [xM 1]       [yM zM]
  !
  !   where A is [m x n], B [m x nrhs] and X [n x nrhs]. In this case n = 2 (linear regression model) and nrhs = 2 (two sets of
  !   observations). The 2-norm |B - A*X| is solved for matrix X by using Singular Value Decomposition (SVD).
  !
  !   For all routines, arguments "A" and "B" are modified (intent is "inout"). In particular, the solution for each r.h.s. is
  !   contained in the first "n" rows of "B". Error flag "OK" returns the absolute highest error (0 = no error) and it must be
  !   initialised to 0 first time the routine is called. This approach is required when routine is inside OPENMP parallelised loops:
  !   in such cases, "OK" must appear in a "!$omp reduction(max:ok)" statement. Argument "rsq" (coefficient of determination) is
  !   optional, since it increases significantly the computation time.
  !
  !   Routines rely on the LAPACK library ("?gelss"). Precision is determined at compile time: use "-DDOUBLE_PREC" to
  !   enable double-precision. OPENMP parallelization is available for nrhs > 1. MKL LAPACK version is enabled by specifying the
  !   "-DMKL" preprocessor flag. The OPENMP parallelization is not recommended when the MKL library is used.
  !
  !   To compile:
  !       gfortran -c -O2 -march=native precisions.f90 llsq.f90 -cpp (-DDOUBLE_PREC, -DMKL) (-fopenmp)
  !       ifort -c -O3 -xHost precisions.f90 llsq.f90 -cpp (-DOUBLE_PREC, -DMKL) (-qopenmp)
  !
  ! Revisions:
  !     Date                    Description of change
  !     ====                    =====================
  !   02/09/20                  original version
  !

  USE                :: omp_lib
  USE, NON_INTRINSIC :: m_precisions
  USE, NON_INTRINSIC :: m_strings

  IMPLICIT none

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  PRIVATE

  PUBLIC :: llsq_solver, llsq_error

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  INTERFACE llsq_solver
    MODULE PROCEDURE llsq_matrix, llsq_vector
  END INTERFACE llsq_solver

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  ! set precisions of input/output arguments at compile-time
#ifdef DOUBLE_PREC
  INTEGER, PARAMETER :: r__ = r64
#else
  INTEGER, PARAMETER :: r__ = r32
#endif

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  CONTAINS

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE llsq_vector(a, b, ok, rsq)

      ! Purpose:
      !   to compute the minimum norm solution to a linear least squares problem (see module header) when nrhs = 1.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   02/09/20                  original version
      !

      REAL(r__),    DIMENSION(:,:),                   INTENT(INOUT) :: a                    !< m-by-n matrix
      REAL(r__),    DIMENSION(:),   TARGET,           INTENT(INOUT) :: b                    !< m-by-nrhs matrix
      INTEGER(i32),                                   INTENT(INOUT) :: ok                   !< error flag
      REAL(r__),                            OPTIONAL, INTENT(OUT)   :: rsq                  !< r-squared value
      REAL(r__),    DIMENSION(1)                                    :: r
      REAL(r__),    DIMENSION(:,:), POINTER                         :: m

      !-----------------------------------------------------------------------------------------------------------------------------

      m(1:SIZE(b), 1:1) => b(:)

      IF (PRESENT(rsq)) THEN
        CALL llsq_matrix(a, m, ok, r)
        rsq = r(1)
      ELSE
        CALL llsq_matrix(a, m, ok)
      ENDIF

    END SUBROUTINE llsq_vector

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE llsq_matrix(a, b, ok, rsq)

      ! Purpose:
      !   to compute the minimum norm solution to a linear least squares problem (see module header) when nrhs > 1.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   02/09/20                  original version
      !

      REAL(r__),                 DIMENSION(:,:),                 INTENT(INOUT) :: a                    !< m-by-n matrix
      REAL(r__),    TARGET,      DIMENSION(:,:),                 INTENT(INOUT) :: b                    !< m-by-nrhs matrix
      INTEGER(i32),                                              INTENT(INOUT) :: ok                   !< error flag
      REAL(r__),                 DIMENSION(SIZE(b,2)), OPTIONAL, INTENT(OUT)   :: rsq                  !< r-squared values
      INTEGER(i32)                                                             :: i, j, m, n
      INTEGER(i32)                                                             :: i0, i1, nrhs
      INTEGER(i32)                                                             :: ierr1, ierr2
      INTEGER(i32)                                                             :: l, lwork, rank, mthreads
      REAL(r__)                                                                :: rcond, sres
      REAL(r__),    ALLOCATABLE, DIMENSION(:)                                  :: s, work
      REAL(r__),    POINTER,     DIMENSION(:,:)                                :: b_ptr

      !-----------------------------------------------------------------------------------------------------------------------------

      IF (PRESENT(rsq)) rsq(:) = 0._r__

      m    = SIZE(a, 1)
      n    = SIZE(a, 2)
      nrhs = SIZE(b, 2)

      IF (m .ne. SIZE(b, 1)) THEN
        ok = 1
        RETURN
      ENDIF

      i0 = 1
      i1 = nrhs

      !$ mthreads = omp_get_max_threads()

      ! there cannot be more threads than r.h.s.
      !$ IF (nrhs .lt. mthreads) CALL omp_set_num_threads(nrhs)

      !$omp parallel default(shared) private(i, j, i0, i1, b_ptr, work, l, rank, rcond, lwork, s, sres) firstprivate(a, nrhs)  &
      !$omp reduction(max:ierr1, ierr2)

      ! split r.h.s. amongst OPENMP threads
      !$ i0 = 1 + NINT( (REAL(nrhs, r__) / omp_get_num_threads()) * omp_get_thread_num() )
      !$ i1 =     NINT( (REAL(nrhs, r__) / omp_get_num_threads()) * (omp_get_thread_num() + 1) )

      nrhs  = i1 - i0 + 1

      b_ptr => b(:, i0:i1)

      ! compute total sum of squares before "b" gets altered
      IF ( PRESENT(rsq) .and. (m .ge. n) ) THEN
        DO i = i0, i1
          rsq(i) = stot(b(:, i))
        ENDDO
      ENDIF

      ! workspace query to calculate optimal size for 'work'
      ALLOCATE(work(1))

      lwork = -1

      l = MIN(m, n)

      ALLOCATE(s(l))

      l = MAX(m, n)

      ! set flag for machine precision
      rcond = -1._r__

      ! first call meant to get correct size for "work"
#ifdef DOUBLE_PREC
#ifdef MKL
      CALL dgels('N', m, n, nrhs, a, m, b_ptr, l, work, lwork, ierr1)
#else
      CALL dgelss(m, n, nrhs, a, m, b_ptr, l, s, rcond, rank, work, lwork, ierr1)
#endif
#else
#ifdef MKL
      CALL sgels('N', m, n, nrhs, a, m, b_ptr, l, work, lwork, ierr1)
#else
      CALL sgelss(m, n, nrhs, a, m, b_ptr, l, s, rcond, rank, work, lwork, ierr1)
#endif
#endif

      lwork = NINT(work(1))

      DEALLOCATE(work)

      ALLOCATE(work(lwork))

#ifdef DOUBLE_PREC
#ifdef MKL
      CALL dgels('N', m, n, nrhs, a, m, b_ptr, l, work, lwork, ierr2)
#else
      CALL dgelss(m, n, nrhs, a, m, b_ptr, l, s, rcond, rank, work, lwork, ierr2)
#endif
#else
#ifdef MKL
      CALL sgels('N', m, n, nrhs, a, m, b_ptr, l, work, lwork, ierr2)
#else
      CALL sgelss(m, n, nrhs, a, m, b_ptr, l, s, rcond, rank, work, lwork, ierr2)
#endif
#endif

      DEALLOCATE(work, s)

#ifdef MKL
      rank = n
#endif

      ! return r-squared for each r.h.s.
      IF ( PRESENT(rsq) .and. (rank .eq. n) .and. (m .ge. n) ) THEN

        DO j = i0, i1

          sres  = 0._r__

          DO i = n + 1, m
            sres = sres + b(i, j)**2
          ENDDO

          rsq(j) = 1._r__ - sres / rsq(j)

        ENDDO

      ELSEIF ( PRESENT(rsq) .and. ((rank .ne. n) .or. (m .ne. n)) ) THEN

        ierr1 = 2

      ENDIF

      !$omp end parallel

      ! revert to initial number of threads
      !$ CALL omp_set_num_threads(mthreads)

      ok = MAX(ok, MAX(ABS(ierr1), ABS(ierr2)))

    END SUBROUTINE llsq_matrix

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    FUNCTION llsq_error(ierr) RESULT(msg)

      ! Purpose:
      !   to translate an error code into a text for all routines in current module.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   02/09/20                  original version
      !

      INTEGER(i32),             INTENT(IN) :: ierr
      CHARACTER(:), ALLOCATABLE            :: msg

      !-----------------------------------------------------------------------------------------------------------------------------

      SELECT CASE(ierr)
        CASE(1)
          msg = 'size mismatch in llsq_solver: size(a,1) must be equal to size(b,1)'

        CASE(2)
          msg = 'r-squared in llq_solver cannot be computed because (rank .ne. n) or (m .ne. n)'

        CASE DEFAULT
          msg = 'dgelss/sgells in llsq_solver failed with exit code: ' // num2char(ierr)

      END SELECT

    END FUNCTION llsq_error

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    REAL(r__) PURE FUNCTION stot(r)

      ! Purpose:
      !   To compute the total sum of squares according to the compensated-summation version of the two-pass algorithm. Calculations
      !   in double-precision to limit roundoff errors.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   02/09/20                  original version
      !

      REAL(r__),  DIMENSION(:), INTENT(IN) :: r
      INTEGER(i32)                         :: i, n
      REAL(r64)                            :: mu, s1, s2, x

      !-----------------------------------------------------------------------------------------------------------------------------

      n = SIZE(r)

      mu = mean(r)

      s1 = 0._r64
      s2 = 0._r64

      DO i = 1, n
        x = r(i) - mu
        s1 = s1 + x
        s2 = s2 + x**2
      ENDDO

      s1 = (s1**2) / n

      stot = (s2 - s1)

    END FUNCTION stot

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    REAL(r__) PURE FUNCTION mean(r)

      ! Purpose:
      !   To compute average. Calculations in double-precision to reduce risk of cancellation.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   02/09/20                  original version
      !

      REAL(r__),   DIMENSION(:), INTENT(IN) :: r
      INTEGER(i32)                          :: i
      INTEGER(i32)                          :: n
      REAL(r64)                             :: v, c

      !-----------------------------------------------------------------------------------------------------------------------------

      n = SIZE(r)

      v = 0._r64
      c = 1._r64

      DO i = 1, n
        v = v + (r(i) - v) / c
        c = c + 1._r64
      ENDDO

      mean = v

    END FUNCTION mean

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

#ifdef DOUBLE_PREC
END MODULE m_llsq_r64
#else
END MODULE m_llsq_r32
#endif
