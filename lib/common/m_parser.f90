MODULE m_parser

  ! Purpose:
  !   To provide routines to parse generic sequential ASCII files, including input and miniseed files. A line of a generic input
  !   file can contain a keyword followed by one or more fields, each introducing a one or multiple strings, integers or real numbers,
  !   and must end with a comma, e.g.:
  !
  !     KEY field_a = 48.2, field_b = {4,2,5}, field_c = {1.,3.; 5.,-2.3; 7.3, 4.5}, field_d = 'any possible string',
  !
  !   To parse each field in such a line, one would open the file and then use:
  !
  !     parse(scalar_real, file_unit, 'field_a', ['=', ','], key = 'KEY')
  !     parse(vector_integer, file_unit, 'field_b', ['{', '}'], key = 'KEY', elm = ',')
  !     parse(matrix_real, file_unit, 'field_c', ['{', '}'], key = 'KEY', elm = ',', grp = ';')
  !     parse(string, file_unit, ["'", "'"], key = 'KEY')
  !
  !   Lines may be commented out, e.g.:
  !
  !     # KEY field_a = 48.2, [...]
  !
  !   In this case the command
  !
  !     parse(scalar_real, file_unit, 'field_a', ['=', ','], com = '#')
  !
  !   allows to skip commented lines. In general, a file may contain the same key more than once:
  !
  !     KEY field_a = 48.2, field_b = {4,2,5}, field_c = {1.,3.; 5.,-2.3; 7.3, 4.5}, field_d = 'any possible string',
  !     KEY field_a = 48.2, field_b = {4,2,5}, field_c = {1.,3.; 5.,-2.3; 7.3, 4.5}, field_d = 'any possible string',
  !
  !   In this case, one can read fields introduced by the a specific key-number, e.g.:
  !
  !     parse(scalar_real, file_unit, 'field_a', ['=', ','], key = 'KEY', nkey = 2)
  !
  !   A key may introduce more than a single field of the same kind, but these must occur on different (continued) lines, e.g.:
  !
  !     KEY field_d = 'string 1', field_a = 48.2, +
  !         field_d = 'string 2', field_a = 31.8, +
  !         field_d = 'string 3', field_a = 7.1,
  !
  !   Under these circumstances, one must indicate also the field number. For instance:
  !
  !     parse(scalar_real, file_unit, 'field_a', ['=', ','], key = KEY, nfield = 3)
  !
  !   would lead to scalar_real = 7.1.
  !   Note that the line continuation symbol is handled automatically but, since it is hard-coded, must be changed (if necessary)
  !   before compiling the module (see parameter "cont" below). Also, a single field cannot be expanded into two or more lines.
  !
  !   Parsing routines can handle single- and double-precision real numbers.
  !
  ! Revisions:
  !     Date                    Description of change
  !     ====                    =====================
  !   02/09/20                  original version
  !

  USE, NON_INTRINSIC :: m_precisions
  USE, NON_INTRINSIC :: m_strings

  IMPLICIT none

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  PRIVATE

  PUBLIC :: parse, parser_error, recurrences, is_empty, read_miniseed, parse_split, parse_split_index

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * -

  INTERFACE parse
    MODULE PROCEDURE parse_char, parse_char_vec, parse_char_mat, parse_int, parse_int_vec, parse_int_mat, parse_real,  &
                     parse_real_vec, parse_real_mat
  END INTERFACE parse

  INTERFACE parse_split
    MODULE PROCEDURE parse_split_r32
  END INTERFACE parse_split

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * -

  CHARACTER(1), PARAMETER :: comma = ',', equal = '=', void = ' '   !< constants for comma, equal and void
  CHARACTER(1), PARAMETER :: cont = '+'                             !< constant for line continuation symbol
  CHARACTER(3), PARAMETER :: char_empty = '???'                     !< constant for "no character found"
  CHARACTER(6), PARAMETER :: fmt = '(A200)'                         !< format to scan generic sequential text file
  INTEGER(i32), PARAMETER :: ch = 200                               !< max expected string length of a numeric field
  INTEGER(i32), PARAMETER :: tab = 9                                !< ASCII code for horizontal tab
  INTEGER(i32), PARAMETER :: int_empty = HUGE(0)                    !< constant for "no integer number found"
  REAL(r32),    PARAMETER :: real_empty = HUGE(0._r32)              !< constant for "no real number found"

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---

  CONTAINS

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    FUNCTION position(buffer, marker, del, rev)

      ! Purpose:
      !   if delimiter "del" is not given, it returns the index of first and last letter of word "marker" inside "buffer". Word means
      !   that "marker" is recognized as a string if preceded only by void, comma or tab and followed by void, equal and tab. If "del"
      !   is specified, it returns the index of the first and last letter of the string introduced by "marker" and delimited by the
      !   characters in "del". If "rev" is set to true, it is expected that "marker" follows the string of interest: also, in this
      !   case, "marker" is considered as single word if followed only by void, comma or tab and preceeded by void, equal and tab.
      !   Indices are set to zero if nothing is found.
      !
      !   Examples:
      !     buffer = "f1 = 20.3, f2 = -5.2, f3 = 3.97", marker = "f2", del = ["=", ","], returns [16, 20]
      !     buffer = "20.3 f1, -5.2 f2, 3.97 f3", marker = "f2", del = [" ", ","], rev = .true., returns [9, 13]
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   02/09/20                  original version
      !

      CHARACTER(*),                         INTENT(IN) :: buffer, marker
      CHARACTER(1), DIMENSION(2), OPTIONAL, INTENT(IN) :: del
      LOGICAL,                    OPTIONAL, INTENT(IN) :: rev
      CHARACTER(1)                                     :: tmp
      INTEGER(i32)                                     :: loc
      INTEGER(i32), DIMENSION(2)                       :: position
      LOGICAL                                          :: backward

      !-----------------------------------------------------------------------------------------------------------------------------

      loc = 1                           !< a string index cannot start from zero

      backward = .false.

      ! set result to zero, implying that "marker" was not found in "buffer"
      position(:) = 0

      ! set-up flag for rev search: IF enabled, it is expected that "marker" follows the numeric value it refers to
      IF (PRESENT(rev)) backward = rev

      ! scan "buffer" letter by letter to find first and last occurrance of "marker"
      DO

        ! quit loop immediately if "marker" is not inside "buffer"
        IF (INDEX(buffer(loc:), marker) .eq. 0) EXIT

        IF (backward) THEN
          loc = INDEX(buffer(loc:), marker)
        ELSE
          loc = loc + INDEX(buffer(loc:), marker) - 1
        ENDIF

        ! if something was found, we verify that "marker" is a single word, i.e. it is not part of a larger string. "marker" is
        ! considered as a single word if preceded only by void, comma or tab and followed by void, equal and tab. For reversed
        ! search, these conditions are swapped.
        IF (loc .ne. 0) THEN

          IF (backward) THEN

            tmp = buffer(loc + LEN(marker):loc + LEN(marker))

            ! IF char after "marker" is not void, comma or tab, move one position forward and repeat search
            IF ( (tmp .ne. void) .and. (tmp .ne. comma) .and. (IACHAR(tmp) .ne. 9) ) THEN
              loc = loc + 1
              CYCLE
            ENDIF

            IF (loc .gt. 1) THEN

              tmp = buffer(loc-1:loc-1)

              ! IF char before is not void, equal or tab, move one position forward and repeat search
              IF ( (tmp .ne. void) .and. (tmp .ne. equal) .and. (IACHAR(tmp) .ne. 9) ) THEN
                loc = loc + 1
                CYCLE
              ENDIF

            ENDIF

          ELSE

            IF (loc .gt. 1) THEN

              tmp = buffer(loc-1:loc-1)

              ! IF char before is not void, comma or tab, move one position forward and repeat search
              IF ( (tmp .ne. void) .and. (tmp .ne. comma) .and. (IACHAR(tmp) .ne. 9) ) THEN
                loc = loc + 1
                CYCLE
              ENDIF

            ENDIF

            tmp = buffer(loc + LEN(marker):loc + LEN(marker))

            ! IF char after "marker" is not void, equal or tab, move one position forward and repeat search
            IF ( (tmp .ne. void) .and. (tmp .ne. equal) .and. (IACHAR(tmp) .ne. 9) ) THEN
              loc = loc + 1
              CYCLE
            ENDIF

          ENDIF

          ! if execution arrives at this point, it means "marker" was recognized as a single word in "buffer".

          ! assign index of first and last letter of "marker" in "buffer"
          position(:) = [loc, loc + LEN(marker)]

          ! move on to return first/last index of a value introduced by "marker" and delimited by "del"
          IF (PRESENT(del)) THEN

            IF (backward) THEN

              !loc = loc - INDEX(reverse_string(buffer(1:loc-1)), del(1))
              loc = INDEX(buffer(1:loc - 1), del(1), back = .true.) - 1

              ! index of last letter
              position(2) = loc

              ! take whole line IF second delimiter is empty
              IF (del(2) .eq. void) THEN
                loc = LEN(buffer)
              ELSE
                !loc = loc - INDEX(reverse_string(buffer(1:loc-1)),del2) + 1
                loc = INDEX(buffer(1:loc - 1), del(2), back = .true.) + 1
              ENDIF

              ! index of first letter
              position(1) = loc

            ELSE

              ! fix case when one delimiter is missing and uncorrectly assigns next one (e.g. avoid returning "y" when "field1 x, field2= y,")
              IF ((del(2) .ne. void) .and. (INDEX(buffer(loc:), del(1)) .gt. INDEX(buffer(loc:), del(2)))) THEN

                position(:) = 0

              ELSE

                loc = loc + INDEX(buffer(loc:), del(1))

                ! index of first letter
                position(1) = loc

                ! take whole line if second delimiter is empty
                IF (del(2) .eq. void) THEN
                  loc = LEN(buffer)
                ELSE
                  loc = loc + INDEX(buffer(loc:), del(2)) - 2
                ENDIF

                ! index of first letter
                position(2) = loc

              ENDIF

            ENDIF

          ENDIF

        ENDIF

        EXIT

      ENDDO

    END FUNCTION position

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    LOGICAL FUNCTION is_commented(buffer, com)

      ! Purpose:
      !   to determine whether string "buffer" is introduced by comment marker "com", i.e. if a line is commented or not.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   02/09/20                  original version
      !

      CHARACTER(*),           INTENT(IN) :: buffer
      CHARACTER(1), OPTIONAL, INTENT(IN) :: com

      !-----------------------------------------------------------------------------------------------------------------------------

      is_commented = .false.

      IF (PRESENT(com)) is_commented = INDEX(ADJUSTL(TRIM(buffer)), com) .eq. 1

    END FUNCTION is_commented

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    FUNCTION recurrences(ok, fo, field, key, nkey, com)

      ! Purpose:
      !   to count how many times "field" occurs in file unit "fo". If requested, only "field" strings introduced by the "nkey"-th
      !   keyword "key" are considered.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   02/09/20                  original version
      !

      INTEGER(i32),                         INTENT(OUT) :: ok
      INTEGER(i32),                         INTENT(IN)  :: fo
      CHARACTER(*),                         INTENT(IN)  :: field
      CHARACTER(*),               OPTIONAL, INTENT(IN)  :: key
      INTEGER(i32),               OPTIONAL, INTENT(IN)  :: nkey
      CHARACTER(1),               OPTIONAL, INTENT(IN)  :: com
      CHARACTER(200)                                    :: buffer
      INTEGER(i32)                                      :: recurrences
      INTEGER(i32)                                      :: ios, keycount
      INTEGER(i32),  DIMENSION(2)                       :: pos
      LOGICAL                                           :: skipkey

      !-----------------------------------------------------------------------------------------------------------------------------

      ok          = 0
      recurrences = 0
      keycount    = 0
      skipkey     = .false.

      DO

        READ(fo,(fmt), iostat = ios) buffer

        IF (ios .eq. -1) EXIT

        IF (ios .ne. 0) THEN
          ok = 6
          RETURN
        ENDIF

        ! jump to next line if current one is empty
        IF (LEN_TRIM(buffer) .eq. 0) CYCLE

        ! jump to next line if current one is commented
        IF (is_commented(buffer, com)) CYCLE

        ! jump to next line if we stumbled over wrong key or key number
        IF ( PRESENT(key) .and. (skipkey .eqv. .false.) ) THEN
          IF (verify_key(buffer, key, keycount, nkey) .eqv. .false.) CYCLE
        ENDIF

        ! find if and where "field" is a word inside "buffer"
        pos = position(buffer, field)

        ! if so, increase counter
        IF (ALL(pos .ne. 0)) recurrences = recurrences + 1

        ! no need to check key for next line as long as continuation symbol is found at end of current line
        skipkey = is_continued(buffer)

      ENDDO

      REWIND(fo, iostat = ios)

      IF (ios .ne. 0) ok = 6

    END FUNCTION recurrences

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE parse_char(ok, v, fo, field, del, key, nkey, nfield, com, rev)

      ! Purpose:
      !   to parse file attached to unit "fo" to return the string preceeded by "nfield"-th "field", delimited by "del" and
      !   associated to "nkey"-th "key". A record is neglected if introduced by "com". The string is expected to preceed "field" if
      !   condition "rev" is true. All other parse routines are built on this one.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   02/09/20                  original version
      !

      INTEGER(i32),                                       INTENT(OUT) :: ok
      CHARACTER(:),  ALLOCATABLE,                         INTENT(OUT) :: v
      INTEGER(i32),                                       INTENT(IN)  :: fo
      CHARACTER(*),                                       INTENT(IN)  :: field
      CHARACTER(1),               DIMENSION(2),           INTENT(IN)  :: del
      CHARACTER(*),                             OPTIONAL, INTENT(IN)  :: key
      INTEGER(i32),                             OPTIONAL, INTENT(IN)  :: nkey, nfield
      CHARACTER(1),                             OPTIONAL, INTENT(IN)  :: com
      LOGICAL,                                  OPTIONAL, INTENT(IN)  :: rev
      CHARACTER(200)                                                  :: buffer
      INTEGER(i32)                                                    :: keycount, fieldcount, ios
      INTEGER(i32),               DIMENSION(2)                        :: pos
      LOGICAL                                                         :: skipkey

      !-----------------------------------------------------------------------------------------------------------------------------

      ok         = 0
      keycount   = 0
      fieldcount = 0
      skipkey    = .false.

      ! set default result to empty
      v = char_empty

      DO

        READ(fo,(fmt), iostat = ios) buffer

        IF (ios .eq. -1) EXIT                       !< EOF reached: quit loop

        IF (ios .ne. 0) THEN
          ok = 7
          RETURN
        ENDIF

        ! jump to next line if current one is empty
        IF (LEN_TRIM(buffer) .eq. 0) CYCLE

        ! jump to next line if current one is commented
        IF (is_commented(buffer, com)) CYCLE

        ! jump to next line if we stumbled over wrong key or key number
        IF ( PRESENT(key) .and. (skipkey .eqv. .false.) ) THEN
          IF (verify_key(buffer, key, keycount, nkey) .eqv. .false.) CYCLE
        ENDIF

        ! find if "field" is a word inside "buffer" and, if so, return position of first and last value associated to "field"
        pos = position(buffer, field, del, rev)

        IF (ALL(pos .ne. 0)) THEN

          fieldcount = fieldcount + 1                  !< update field counter

          ! jump to next line if we stumbled over wrong "field" number
          IF (PRESENT(nfield)) THEN
            IF (fieldcount .ne. nfield) THEN
              skipkey = is_continued(buffer)          !< disable "key" checking if line ends with continuation symbol
              CYCLE
            ENDIF
          ENDIF

          v = buffer(pos(1):pos(2))

          EXIT                                        !< stop reading file

        ENDIF

        ! disable "key" checking if line ends with continuation symbol
        skipkey = is_continued(buffer)

      ENDDO

      REWIND(fo, iostat = ios)

      IF (ios .ne. 0) ok = 7

    END SUBROUTINE parse_char

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE parse_int(ok, v, fo, field, del, key, nkey, nfield, com, rev)

      ! Purpose:
      !   to parse file attached to unit "fo" to return the numeric value preceeded by "nfield"-th "field", delimited by "del" and
      !   associated to "nkey"-th "key". A record is neglected if introduced by "com". The numeric value is expected to preceed
      !   "field" if condition "rev" is true.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   02/09/20                  original version
      !

      INTEGER(i32),                                      INTENT(OUT) :: ok, v
      INTEGER(i32),                                      INTENT(IN)  :: fo
      CHARACTER(*),                                      INTENT(IN)  :: field
      CHARACTER(1),              DIMENSION(2),           INTENT(IN)  :: del
      CHARACTER(*),                            OPTIONAL, INTENT(IN)  :: key
      INTEGER(i32),                            OPTIONAL, INTENT(IN)  :: nkey, nfield
      CHARACTER(1),                            OPTIONAL, INTENT(IN)  :: com
      LOGICAL,                                 OPTIONAL, INTENT(IN)  :: rev
      INTEGER(i32), ALLOCATABLE, DIMENSION(:)                        :: x

      !-----------------------------------------------------------------------------------------------------------------------------

      CALL parse(ok, x, fo, field, del, '?', key, nkey, nfield, com, rev)

      IF (ok .ne. 0) RETURN

      v = x(1)

    END SUBROUTINE parse_int

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE parse_real(ok, v, fo, field, del, key, nkey, nfield, com, rev)

      ! Purpose:
      !   to parse file attached to unit "fo" to return the numeric value preceeded by "nfield"-th "field", delimited by "del" and
      !   associated to "nkey"-th "key". A record is neglected if introduced by "com". The numeric value is expected to preceed
      !   "field" if condition "rev" is true.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   02/09/20                  original version
      !

      INTEGER(i32),                                      INTENT(OUT) :: ok
      REAL(r32),                                         INTENT(OUT) :: v
      INTEGER(i32),                                      INTENT(IN)  :: fo
      CHARACTER(*),                                      INTENT(IN)  :: field
      CHARACTER(1),              DIMENSION(2),           INTENT(IN)  :: del
      CHARACTER(*),                            OPTIONAL, INTENT(IN)  :: key
      INTEGER(i32),                            OPTIONAL, INTENT(IN)  :: nkey, nfield
      CHARACTER(1),                            OPTIONAL, INTENT(IN)  :: com
      LOGICAL,                                 OPTIONAL, INTENT(IN)  :: rev
      REAL(r32),    ALLOCATABLE, DIMENSION(:)                        :: x

      !-----------------------------------------------------------------------------------------------------------------------------

      CALL parse(ok, x, fo, field, del, '?', key, nkey, nfield, com, rev)

      IF (ok .ne. 0) RETURN

      v = x(1)

    END SUBROUTINE parse_real

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE parse_char_vec(ok, v, fo, field, del, elm, key, nkey, nfield, com, rev)

      ! Purpose:
      !   to parse file attached to unit "fo" to return the set of strings preceeded by "nfield"-th "field", delimited by "del" and
      !   associated to "nkey"-th "key". Each string is separated according to "elm". A record is neglected if introduced by "com".
      !   Strings are expected to preceed "field" if condition "rev" is true.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   02/09/20                  original version
      !

      INTEGER(i32),                                        INTENT(OUT) :: ok
      CHARACTER(*),  ALLOCATABLE, DIMENSION(:),            INTENT(OUT) :: v
      INTEGER(i32),                                        INTENT(IN)  :: fo
      CHARACTER(*),                                        INTENT(IN)  :: field
      CHARACTER(1),               DIMENSION(2),            INTENT(IN)  :: del
      CHARACTER(1),                                        INTENT(IN)  :: elm
      CHARACTER(*),                              OPTIONAL, INTENT(IN)  :: key
      INTEGER(i32),                              OPTIONAL, INTENT(IN)  :: nkey, nfield
      CHARACTER(1),                              OPTIONAL, INTENT(IN)  :: com
      LOGICAL,                                   OPTIONAL, INTENT(IN)  :: rev
      CHARACTER(ch), ALLOCATABLE, DIMENSION(:,:)                       :: x

      !-----------------------------------------------------------------------------------------------------------------------------

      CALL parse(ok, x, fo, field, del, elm, '?', key, nkey, nfield, com, rev)

      IF (ok .ne. 0) RETURN

      v = x(:, 1)

    END SUBROUTINE parse_char_vec

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE parse_int_vec(ok, v, fo, field, del, elm, key, nkey, nfield, com, rev)

      ! Purpose:
      !   to parse file attached to unit "fo" to return the numeric values preceeded by "nfield"-th "field", delimited by "del" and
      !   associated to "nkey"-th "key". Vector elements are separated according to "elm". A record is neglected if introduced by
      !   "com". The numeric value is expected to preceed "field" if condition "rev" is true.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   02/09/20                  original version
      !

      INTEGER(i32),                                       INTENT(OUT) :: ok
      INTEGER(i32), ALLOCATABLE, DIMENSION(:),            INTENT(OUT) :: v
      INTEGER(i32),                                       INTENT(IN)  :: fo
      CHARACTER(*),                                       INTENT(IN)  :: field
      CHARACTER(1),              DIMENSION(2),            INTENT(IN)  :: del
      CHARACTER(1),                                       INTENT(IN)  :: elm
      CHARACTER(*),                             OPTIONAL, INTENT(IN)  :: key
      INTEGER(i32),                             OPTIONAL, INTENT(IN)  :: nkey, nfield
      CHARACTER(1),                             OPTIONAL, INTENT(IN)  :: com
      LOGICAL,                                  OPTIONAL, INTENT(IN)  :: rev
      INTEGER(i32), ALLOCATABLE, DIMENSION(:,:)                       :: x

      !-----------------------------------------------------------------------------------------------------------------------------

      CALL parse(ok, x, fo, field, del, elm, '?', key, nkey, nfield, com, rev)

      IF (ok .ne. 0) RETURN

      v = x(:, 1)

    END SUBROUTINE parse_int_vec

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE parse_real_vec(ok, v, fo, field, del, elm, key, nkey, nfield, com, rev)

      ! Purpose:
      !   to parse file attached to unit "fo" to return the numeric values preceeded by "nfield"-th "field", delimited by "del" and
      !   associated to "nkey"-th "key". Vector elements are separated according to "elm". A record is neglected if introduced by
      !   "com". The numeric value is expected to preceed "field" if condition "rev" is true.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   02/09/20                  original version
      !

      INTEGER(i32),                                       INTENT(OUT) :: ok
      REAL(r32),    ALLOCATABLE, DIMENSION(:),            INTENT(OUT) :: v
      INTEGER(i32),                                       INTENT(IN)  :: fo
      CHARACTER(*),                                       INTENT(IN)  :: field
      CHARACTER(1),              DIMENSION(2),            INTENT(IN)  :: del
      CHARACTER(1),                                       INTENT(IN)  :: elm
      CHARACTER(*),                             OPTIONAL, INTENT(IN)  :: key
      INTEGER(i32),                             OPTIONAL, INTENT(IN)  :: nkey, nfield
      CHARACTER(1),                             OPTIONAL, INTENT(IN)  :: com
      LOGICAL,                                  OPTIONAL, INTENT(IN)  :: rev
      REAL(r32),    ALLOCATABLE, DIMENSION(:,:)                       :: x

      !-----------------------------------------------------------------------------------------------------------------------------

      CALL parse(ok, x, fo, field, del, elm, '?', key, nkey, nfield, com, rev)

      IF (ok .ne. 0) RETURN

      v = x(:, 1)

    END SUBROUTINE parse_real_vec

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE parse_char_mat(ok, v, fo, field, del, elm, grp, key, nkey, nfield, com, rev)

      ! Purpose:
      !   to parse file attached to unit "fo" to return the set of strings preceeded by "nfield"-th "field", delimited by "del" and
      !   associated to "nkey"-th "key". Column and row elements are separated according to "elm" and "grp", respectively. A record
      !   is neglected if introduced by "com". Strings are expected to preceed "field" if condition "rev" is true.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   02/09/20                  original version
      !

      INTEGER(i32),                                         INTENT(OUT) :: ok
      CHARACTER(*),  ALLOCATABLE, DIMENSION(:,:),           INTENT(OUT) :: v
      INTEGER(i32),                                         INTENT(IN)  :: fo
      CHARACTER(*),                                         INTENT(IN)  :: field
      CHARACTER(1),               DIMENSION(2),             INTENT(IN)  :: del
      CHARACTER(1),                                         INTENT(IN)  :: elm, grp
      CHARACTER(*),                               OPTIONAL, INTENT(IN)  :: key
      INTEGER(i32),                               OPTIONAL, INTENT(IN)  :: nkey, nfield
      CHARACTER(1),                               OPTIONAL, INTENT(IN)  :: com
      LOGICAL,                                    OPTIONAL, INTENT(IN)  :: rev
      CHARACTER(ch)                                                     :: u, z
      CHARACTER(:),  ALLOCATABLE                                        :: x
      CHARACTER(ch), ALLOCATABLE, DIMENSION(:)                          :: w
      INTEGER(i32)                                                      :: i, j

      !-----------------------------------------------------------------------------------------------------------------------------

      CALL parse(ok, x, fo, field, del, key, nkey, nfield, com, rev)

      IF (ok .ne. 0) RETURN

      IF (is_empty(x)) THEN

        ALLOCATE(v(1, 1))

        v = char_empty

      ELSE

        j = 1

        DO

          z = split(x, grp, j)

          IF (z .eq. '') EXIT

          i = 1

          DO

            u = split(z, elm, i)

            IF (u .eq. '') EXIT

            IF (ALLOCATED(w)) THEN
              w = [w, u]
            ELSE
              w = [u]
            ENDIF

            i = i + 1

          ENDDO

          j = j + 1

        ENDDO

        v = RESHAPE(w, [i-1, j-1])

      ENDIF

    END SUBROUTINE parse_char_mat

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE parse_int_mat(ok, v, fo, field, del, elm, grp, key, nkey, nfield, com, rev)

      ! Purpose:
      !   to parse file attached to unit "fo" to return the numeric values preceeded by "nfield"-th "field", delimited by "del" and
      !   associated to "nkey"-th "key". Column and row elements are separated according to "elm" and "grp", respectively. A record
      !   is neglected if introduced by "com". The numeric value is expected to preceed "field" if condition "rev" is true.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   02/09/20                  original version
      !

      INTEGER(i32),                                         INTENT(OUT) :: ok
      INTEGER(i32),  ALLOCATABLE, DIMENSION(:,:),           INTENT(OUT) :: v
      INTEGER(i32),                                         INTENT(IN)  :: fo
      CHARACTER(*),                                         INTENT(IN)  :: field
      CHARACTER(1),               DIMENSION(2),             INTENT(IN)  :: del
      CHARACTER(1),                                         INTENT(IN)  :: elm, grp
      CHARACTER(*),                               OPTIONAL, INTENT(IN)  :: key
      INTEGER(i32),                               OPTIONAL, INTENT(IN)  :: nkey, nfield
      CHARACTER(1),                               OPTIONAL, INTENT(IN)  :: com
      LOGICAL,                                    OPTIONAL, INTENT(IN)  :: rev
      CHARACTER(ch), ALLOCATABLE, DIMENSION(:,:)                        :: x
      INTEGER(i32)                                                      :: i, j, ierr

      !-----------------------------------------------------------------------------------------------------------------------------

      ierr = 0

      CALL parse(ok, x, fo, field, del, elm, grp, key, nkey, nfield, com, rev)

      IF (ok .ne. 0) RETURN

      ALLOCATE(v(SIZE(x,1),SIZE(x,2)))

      DO j = 1, SIZE(x, 2)
        DO i = 1, SIZE(x, 1)
          IF (TRIM(x(i, j)) .eq. char_empty) THEN
            v(i, j) = int_empty
          ELSE
            READ(x(i, j), *, iostat = ierr) v(i, j)
            IF (ierr .ne. 0) EXIT
          ENDIF
        ENDDO
      ENDDO

      IF (ierr .ne. 0) ok = 8

    END SUBROUTINE parse_int_mat

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE parse_real_mat(ok, v, fo, field, del, elm, grp, key, nkey, nfield, com, rev)

      ! Purpose:
      !   to parse file attached to unit "fo" to return the numeric values preceeded by "nfield"-th "field", delimited by "del" and
      !   associated to "nkey"-th "key". Column and row elements are separated according to "elm" and "grp", respectively. A record
      !   is neglected if introduced by "com". The numeric value is expected to preceed "field" if condition "rev" is true.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   02/09/20                  original version
      !

      INTEGER(i32),                                         INTENT(OUT) :: ok
      REAL(r32),     ALLOCATABLE, DIMENSION(:,:),           INTENT(OUT) :: v
      INTEGER(i32),                                         INTENT(IN)  :: fo
      CHARACTER(*),                                         INTENT(IN)  :: field
      CHARACTER(1),               DIMENSION(2),             INTENT(IN)  :: del
      CHARACTER(1),                                         INTENT(IN)  :: elm, grp
      CHARACTER(*),                               OPTIONAL, INTENT(IN)  :: key
      INTEGER(i32),                               OPTIONAL, INTENT(IN)  :: nkey, nfield
      CHARACTER(1),                               OPTIONAL, INTENT(IN)  :: com
      LOGICAL,                                    OPTIONAL, INTENT(IN)  :: rev
      CHARACTER(ch), ALLOCATABLE, DIMENSION(:,:)                        :: x
      INTEGER(i32)                                                      :: i, j, ierr

      !-----------------------------------------------------------------------------------------------------------------------------

      CALL parse(ok, x, fo, field, del, elm, grp, key, nkey, nfield, com, rev)

      IF (ok .ne. 0) RETURN

      ALLOCATE(v(SIZE(x,1), SIZE(x,2)))

      DO j = 1, SIZE(x, 2)
        DO i = 1, SIZE(x, 1)
          IF (TRIM(x(i, j)) .eq. char_empty) THEN
            v(i, j) = real_empty
          ELSE
            READ(x(i, j), *, iostat = ierr) v(i, j)
            IF (ierr .ne. 0) EXIT
          ENDIF
        ENDDO
      ENDDO

      IF (ierr .ne. 0) ok = 9

    END SUBROUTINE parse_real_mat

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    LOGICAL FUNCTION verify_key(buffer, key, keycount, nkey, com)

      ! Purpose:
      !   to determine if string "key" is present inside string "buffer" as the first single (isolated) word and, if so, if it
      !   corresponds to the "nkey"-th entry.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   02/09/20                  original version
      !   08/03/21                  added ADJUSTL since "key" must be first word
      !

      CHARACTER(*),                        INTENT(IN)    :: buffer, key
      INTEGER(i32),                        INTENT(INOUT) :: keycount
      INTEGER(i32),              OPTIONAL, INTENT(IN)    :: nkey
      CHARACTER(1),              OPTIONAL, INTENT(IN)    :: com
      INTEGER(i32), DIMENSION(2)                         :: pos

      !-----------------------------------------------------------------------------------------------------------------------------

      verify_key = .true.

      ! find "key" inside buffer (remove leading spaces)
      pos = position(ADJUSTL(buffer), key)

      IF (ALL(pos .eq. 0) .or. (pos(1) .ne. 1)) THEN             !< "key" not present or not first word
        verify_key = .false.
      ELSE
        keycount = keycount + 1                                  !< update keyword counter
        IF (PRESENT(nkey)) THEN
          IF (keycount .ne. nkey) verify_key = .false.           !< verify key number matching
        ENDIF
      ENDIF

    END FUNCTION verify_key

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    LOGICAL FUNCTION is_continued(buffer)

      ! Purpose:
      !   to determined if string "buffer" terminates with the predefined continuation symbol.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   02/09/20                  original version
      !

      CHARACTER(*), INTENT(IN) :: buffer
      INTEGER(i32)             :: i

      !-----------------------------------------------------------------------------------------------------------------------------

      i = LEN_TRIM(buffer)

      ! return "true" if last letter is continue-to-next-line symbol "cont"
      is_continued = buffer(i:i) .eq. cont

    END FUNCTION is_continued

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    ! ELEMENTAL SUBROUTINE set_to_empty(x)
    !
    !   CLASS(*), INTENT(INOUT) :: x
    !
    !   !-----------------------------------------------------------------------------------------------------------------------------
    !
    !   SELECT TYPE(x)
    !     TYPE IS (CHARACTER(*))
    !       x = char_empty
    !     TYPE IS (INTEGER(i32))
    !       x = int_empty
    !     TYPE IS (REAL(r32))
    !       x = real_empty
    !     TYPE IS (REAL(r64))
    !       x = real_empty
    !   END SELECT
    !
    ! END SUBROUTINE set_to_empty

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    LOGICAL ELEMENTAL FUNCTION is_empty(x)

      ! Purpose:
      !   to determine if input value "x" is empty according to the predefined is-empty constants.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   02/09/20                  original version
      !

      CLASS(*), INTENT(IN) :: x

      !-----------------------------------------------------------------------------------------------------------------------------

      SELECT TYPE(x)
        TYPE IS (CHARACTER(*))
          is_empty = x .eq. char_empty
        TYPE IS (INTEGER(i32))
          is_empty = x .eq. int_empty
        TYPE IS (REAL(r32))
          is_empty = x .eq. real_empty
        TYPE IS (REAL(r64))
          is_empty = x .eq. real_empty
      END SELECT

    END FUNCTION is_empty

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE read_miniseed(ok, fo, dt, timeseries, tmax)

      ! Purpose:
      !   to read a miniseed ASCII file "fo" and return its content in terms of time-step "dt" and values "timeseries" up to "tmax"
      !   (take whole record if the latter is not specified). Although each miniseed file contains one component of motion, multiple-
      !   component can be specified by using square brackets, e.g. "file_name_[NEZ]_xxx".
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   02/09/20                  original version
      !   21/01/21                  added option to read up "tmax"
      !

      INTEGER(i32),                                        INTENT(OUT) :: ok
      CHARACTER(*),                                        INTENT(IN)  :: fo
      REAL(r32),                                           INTENT(OUT) :: dt
      REAL(r32),    ALLOCATABLE, DIMENSION(:,:),           INTENT(OUT) :: timeseries
      REAL(r32),                                 OPTIONAL, INTENT(IN)  :: tmax
      CHARACTER(:), ALLOCATABLE                                        :: buffer
      INTEGER(i32)                                                     :: i, j, k, ierr
      INTEGER(i32)                                                     :: npts, ncomp, lu
      INTEGER(i32)                                                     :: sbrack, ebrack

      !-----------------------------------------------------------------------------------------------------------------------------

      ok = 0

      ncomp = 1

      ! determine if "fo" indicates one or more components of motion (e.g. ...[nez]...)

      sbrack = INDEX(fo, '[')
      ebrack = INDEX(fo, ']')

      IF (sbrack .ne. 0) THEN
        ncomp = ebrack - sbrack - 1
      ENDIF

      DO k = 1, ncomp

        IF (ncomp .eq. 1) THEN
          buffer = fo
        ELSE
          buffer = fo(1:sbrack - 1) + fo(sbrack + k:sbrack + k) + fo(ebrack + 1:)
        ENDIF

        OPEN(newunit = lu, file = TRIM(buffer), status = 'old', form = 'formatted', access = 'sequential', action = 'read', &
             iostat = ierr)

        IF (ierr .ne. 0) THEN
          ok = 1
          RETURN
        ENDIF

        ! get number of points
        CALL parse(ierr, npts, lu, 'samples', [' ', ','], rev = .true.)

        IF (is_empty(npts) .or. (ierr .ne. 0)) THEN
          ok = 2
          RETURN
        ENDIF

        ! get time sampling
        CALL parse(ierr, dt, lu, 'sps', [' ', ','], rev = .true.)

        IF (is_empty(dt) .or. (ierr .ne. 0)) THEN
          ok = 3
          RETURN
        ENDIF

        dt = 1._r32 / dt

        ! determine number of points to be read
        IF (PRESENT(tmax)) THEN
          IF (tmax .gt. 0._r32) npts = MIN(npts, NINT(tmax / dt))
        ENDIF

        ! read only multiple of 6 lines
        npts = 6 * (npts / 6)

        IF (.not.ALLOCATED(timeseries)) ALLOCATE(timeseries(npts, ncomp))

        ! skip one-line header
        READ(lu, (fmt), iostat = ierr)

        IF (ierr .ne. 0) THEN
          ok = 4
          RETURN
        ENDIF

        ! now read values
        DO i = 1, npts / 6
          READ(lu, *, iostat = ierr) (timeseries((i - 1)*6 + j, k), j = 1, 6)
          IF (ierr .ne. 0) EXIT
        ENDDO

        IF (ierr .ne. 0) THEN
          ok = 4
          RETURN
        ENDIF

        CLOSE(lu, iostat = ierr)

        IF (ierr .ne. 0) THEN
          ok = 5
          RETURN
        ENDIF

      ENDDO

    END SUBROUTINE read_miniseed

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE parse_split_r32(line, numeric)

      ! Purpose:
      !   to read a line containing an arbitrary number of numeric values and returning each entry separately.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   08/03/21                  original version
      !

      CHARACTER(*),                                        INTENT(IN)  :: line
      REAL(r32),    ALLOCATABLE, DIMENSION(:),             INTENT(OUT) :: numeric
      INTEGER(i32)                                                     :: n
      REAL(r32),                 DIMENSION(LEN_TRIM(line))             :: x

      !-----------------------------------------------------------------------------------------------------------------------------

      n = 1

      DO
        READ(line, *, ERR = 10, END = 10) x(1:n)
        n = n + 1
      ENDDO

      10 CONTINUE

      n = n - 1

      numeric = x(1:n)

    END SUBROUTINE parse_split_r32

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    FUNCTION parse_split_index(line, word, separator) RESULT(col)

      ! Purpose:
      !   to split "line" according to "separator" and return at which index "word" is found (0 otherwise)
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   08/03/21                  original version
      !

      CHARACTER(*),                                        INTENT(IN)  :: line, word, separator
      INTEGER(i32)                                                     :: i, w, column, col
      CHARACTER(:), ALLOCATABLE                                        :: buffer

      !-----------------------------------------------------------------------------------------------------------------------------

      w = LEN(word)

      column = 0
      col    = 0

      buffer = separator + TRIM(line)

      DO i = 2, LEN_TRIM(buffer) - w + 1

        IF (buffer(i:i) .ne. separator) THEN          !< we have some text

          IF (buffer(i-1:i-1) .eq. separator) THEN    !< previous character was a separator

            column = column + 1

            IF (buffer(i:i + w - 1) .eq. word) THEN   !< word found: stop search
              col = column
              EXIT
            ENDIF

          ENDIF

        ENDIF

      ENDDO

    END FUNCTION parse_split_index

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    FUNCTION parser_error(ierr) RESULT(msg)

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
          msg = 'IO error in read_miniseed: could not open file'

        CASE(2)
          msg = 'IO error in read_miniseed: number of samples not found'

        CASE(3)
          msg = 'IO error in read_miniseed: sampling interval not found'

        CASE(4)
          msg = 'IO error in read_miniseed: error while reading file'

        CASE(5)
          msg = 'IO error in read_miniseed: could not close file'

        CASE(6)
          msg = 'IO error in recurrences: error while reading file'

        CASE(7)
          msg = 'IO error in parse (parse_char): error while reading file'

        CASE(8)
          msg = 'error in parse (parse_int_mat): string-to-integer conversion error'

        CASE(9)
          msg = 'error in parse (parse_real_mat): string-to-real conversion error'

      END SELECT


    END FUNCTION parser_error

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *



END MODULE m_parser
