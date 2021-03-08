MODULE m_strings

  USE, NON_INTRINSIC :: m_precisions

  IMPLICIT none

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  PRIVATE

  PUBLIC :: num2char, operator (+), uppercase, lowercase, char2num, split

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  INTERFACE OPERATOR (+)
    MODULE PROCEDURE add
  END INTERFACE

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  CONTAINS

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    FUNCTION add(a, b)

      ! Purpose:
      !   to concatenate two strings of arbitrary length into a new one
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   02/09/20                  original version
      !

      CHARACTER(*),              INTENT(IN) :: a, b
      CHARACTER(LEN(A) + LEN(B))            :: add

      !----------------------------------------------------------------------------------------------------------------------------

      add = a // b

    END FUNCTION add

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    RECURSIVE FUNCTION num2char(x, notation, width, digits, precision, sign, separator, justify, fill) RESULT(char)

      ! Purpose:
      !   to convert a number (integer or real) into a string or to adjust the length of a string. Available options are:
      !
      !                   notation    width   digits    precision   sign    separator   justify   fill
      !   integers                      x       x                    x                     x        x
      !   floats            x           x                   x        x         x           x        x
      !   characters                    x                                                  x        x
      !   logical                       x                                                  x        x
      !
      !   Valid arguments for optional parameters:
      !     notation  = 'd' (default), 'f' (fixed), 's' (scientific)
      !     width     = any positive integer (default = 0)
      !     digits    = any positive integer (can be larger than "width", default = 0)
      !     precision = any positive integer such that (width - precision) > 6 (default = 0)
      !     sign      = .true. or .false. (default is .false.)
      !     separator = ',' or '.' (default is '.')
      !     justify   = 'l' or 'r' (default is 'r')
      !     fill      = any character, it will be repeated to fill "width-len(char)" spaces
      !
      !   If no options are specified, default format descriptors are "A" (characters), "I0.0" (integers), "G0" (floats)
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   02/09/20                  original version
      !

      CLASS(*),                           INTENT(IN) :: x
      CHARACTER(1),             OPTIONAL, INTENT(IN) :: notation, justify, separator, fill
      INTEGER(i32),             OPTIONAL, INTENT(IN) :: width, digits, precision
      LOGICAL,                  OPTIONAL, INTENT(IN) :: sign
      CHARACTER(:), ALLOCATABLE                      :: char, pat, temp
      CHARACTER(2)                                   :: n, pls, sep
      INTEGER(i32)                                   :: w, m, p, i, r

      !-----------------------------------------------------------------------------------------------------------------------------

      w = 0
      m = 0
      p = 0

      ! by default, use point as decimal separator and don't add a plus in front of positive numbers
      sep = 'DP'
      pls = 'SS'

      IF (PRESENT(separator)) THEN
        IF (separator .eq. ',') sep = 'DC'
      ENDIF

      IF (PRESENT(sign)) THEN
        IF (sign .eqv. .true.) pls = 'SP'
      ENDIF

      SELECT TYPE(x)

        TYPE IS (CHARACTER(*))

          IF (PRESENT(width)) w = width

          SELECT CASE (w)
            CASE(0)
              char = x
            CASE DEFAULT
              ALLOCATE(character(w) :: char)
              WRITE(char, '(A)') x
          END SELECT

        TYPE IS (LOGICAL)

          IF (PRESENT(width)) w = width

          SELECT CASE (w)
            CASE(0)
              ALLOCATE(character(1) :: char)
            CASE DEFAULT
              ALLOCATE(character(w) :: char)
          END SELECT

          WRITE(char, '(L)') x .eqv. .true.

        ! handle integers. Note: "i0.0" or "i0.2" is accepted as format
        TYPE IS (INTEGER(i32))

          n = 'I'

          IF (.not.PRESENT(width) .and. .not.PRESENT(digits)) THEN

            pat = '(' + pls + TRIM(n) + '0' + ')'

            w = RANGE(x) + 2

          ELSE

            IF (PRESENT(width))  w = width
            IF (PRESENT(digits)) m = digits

            IF (w .eq. 0) w = RANGE(x) + 2

            pat = '(' + pls + TRIM(n) + num2char(w) + '.' + num2char(m) + ')'

          ENDIF

          ALLOCATE(character(w) :: char)

          WRITE(char, (pat)) x

        TYPE IS (INTEGER(i64))

          n = 'I'

          IF (PRESENT(width)) w = width

          IF (PRESENT(digits)) m = digits

          IF (w .eq. 0) w = RANGE(x) + 2

          pat = '(' + pls + TRIM(n) + num2char(w) + '.' + num2char(m) + ')'

          ALLOCATE(character(w) :: char)

          WRITE(char, (pat)) x

        ! handle real numbers. Note: width=0 is not accepted for scientific notation
        TYPE IS (REAL(r32))

          n = 'G'

          IF (PRESENT(notation)) THEN
            SELECT CASE (notation)
              CASE('s')
                n = 'ES'
              CASE('f')
                n = 'F'
              CASE DEFAULT
                n = 'G'
            END SELECT
          ENDIF

          IF (PRESENT(width)) w = width

          IF (TRIM(n) .eq. 'G') THEN
            IF (PRESENT(digits)) p = digits
          ELSE
            IF (PRESENT(precision)) p = precision
          ENDIF

          IF (w .eq. 0) w = RANGE(x) + 2

          ! format is "G0" if precision is not specified
          IF (p .eq. 0) THEN
            pat = '(' + pls + ',' + sep + ',' + 'G0' + ')'
          ELSE
            pat = '(' + pls + ',' + sep + ',' + TRIM(n) + num2char(w) + '.' + num2char(p) + ')'
          ENDIF

          ALLOCATE(character(w) :: char)

          WRITE(char, (pat)) x

        TYPE IS (REAL(r64))

          n = 'G'

          IF (PRESENT(notation)) THEN
            SELECT CASE (notation)
              CASE('s')
                n = 'ES'
              CASE('f')
                n = 'F'
              CASE DEFAULT
                n = 'G'
            END SELECT
          ENDIF

          IF (PRESENT(width)) w = width

          IF (TRIM(n) .eq. 'G') THEN
            IF (PRESENT(digits)) p = digits
          ELSE
            IF (PRESENT(precision)) p = precision
          ENDIF

          IF (w .eq. 0) w = RANGE(x) + 2

          ! format is "G0" if precision is not specified
          IF (p .eq. 0) THEN
            pat = '(' + pls + ',' + sep + ',' + 'G0' + ')'
          ELSE
            pat = '(' + pls + ',' + sep + ',' + TRIM(n) + num2char(w) + '.' + num2char(p) + ')'
          ENDIF

          ALLOCATE(character(w) :: char)

          WRITE(char, (pat)) x

      END SELECT

      IF (.not.PRESENT(width)) char = TRIM(ADJUSTL(char))

      ! justify string
      IF (PRESENT(justify)) THEN
        SELECT CASE (justify)
          CASE('l')
            char = ADJUSTL(char)
          CASE('r')
            char = ADJUSTR(char)
        END SELECT
      ENDIF

      ! add fill element (make sense only if "width" is given as argument)
      IF (PRESENT(fill) .and. PRESENT(width)) THEN

        r = LEN(char) - LEN_TRIM(char)

        IF (PRESENT(justify)) THEN
          SELECT CASE (justify)
            CASE('l')
              pat = '(A,' + num2char(r) + 'A)'
              WRITE(char, (pat)) TRIM(char), [(fill, i = 1,r)]
            CASE('r')
              r = LEN(char) - LEN_TRIM(ADJUSTL(char))
              pat = '(' + num2char(r) + 'A,A)'
              temp = TRIM(ADJUSTL(char))                       !< adjustl required because trim remove trailing blanks
              WRITE(char, (pat)) [(fill, i = 1,r)], temp
          END SELECT
        ELSE
          pat = '(A,' + num2char(r) + 'A)'
          WRITE(char, (pat)) TRIM(char), [(fill, i = 1,r)]
        ENDIF

      ENDIF

    END FUNCTION num2char

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    FUNCTION uppercase(x) RESULT(y)

      ! Purpose:
      !   to convert any lower-case character in a string into upper-case
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   02/09/20                  original version
      !

      CHARACTER(*),     INTENT(IN) :: x
      CHARACTER(LEN(x))            :: y
      INTEGER(i32)                 :: i

      !-----------------------------------------------------------------------------------------------------------------------------

      DO i = 1, LEN(x)
        IF ( (x(i:i) .ge. 'a') .and. (x(i:i) .le. 'z') ) THEN
          y(i:i) = ACHAR(IACHAR(x(i:i)) - 32)
        ELSE
          y(i:i) = x(i:i)
        ENDIF
      ENDDO

    END FUNCTION uppercase

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    FUNCTION lowercase(x) RESULT(y)

      ! Purpose:
      !   to convert any upper-case character in a string into lower-case
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   02/09/20                  original version
      !

      CHARACTER(*),     INTENT(IN) :: x
      CHARACTER(LEN(x))            :: y
      INTEGER(i32)                 :: i

      !-----------------------------------------------------------------------------------------------------------------------------

      DO i = 1, LEN(x)
        IF ( (x(i:i) .ge. 'A') .and. (x(i:i) .le. 'Z') ) THEN
          y(i:i) = ACHAR(IACHAR(x(i:i)) + 32)
        ELSE
          y(i:i) = x(i:i)
        ENDIF
      ENDDO

    END FUNCTION lowercase

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE char2num(char, num)

      ! Purpose:
      !   to convert a string into a number.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   02/09/20                  original version
      !

      CHARACTER(*), INTENT(IN)  :: char
      CLASS(*),     INTENT(OUT) :: num

      !-----------------------------------------------------------------------------------------------------------------------------

      SELECT TYPE(num)

        TYPE IS (CHARACTER(*))
          num = char

        TYPE IS (INTEGER(i32))
          READ(char, *) num

        TYPE IS (INTEGER(i64))
          READ(char, *) num

        TYPE IS (REAL(r32))
          READ(char, *) num

        TYPE IS (REAL(r64))
          READ(char, *) num

      END SELECT

    END SUBROUTINE char2num

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    RECURSIVE FUNCTION split(str, char, n) RESULT(chunk)

      ! Purpose:
      !   to split a string according to a specific character (including blank), returning the n-th chunk. If there are less chunks
      !   than queried, an empty string is returned
      !
      !   Example: given string = "this is a sample string", calling "split(string, ' ', n)" returns:
      !     "this",   if n = 1
      !     "string", if n = 5
      !     "",       if n > 5
      !
      !   Note that the output array results always allocated, even if empty.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   02/09/20                  original version
      !

      CHARACTER(*),             INTENT(IN) :: str, char
      INTEGER(i32),             INTENT(IN) :: n
      CHARACTER(:), ALLOCATABLE            :: chunk
      INTEGER(i32)                         :: pos, i

      !-----------------------------------------------------------------------------------------------------------------------------

      ! find first occurrance of "char" into "str"
      pos = INDEX(str, char)

      IF (pos .ne. 0) THEN
        chunk = str(1:pos - 1)
        i = n - 1
      ELSE
        IF (n .eq. 1) THEN
          chunk = str
        ELSE
          chunk = ""
        ENDIF
        i = 0
      ENDIF

      IF (i .ne. 0) chunk = split(str(pos + 1:LEN(str)), char, i)

    END FUNCTION split

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *


    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

END MODULE m_strings
