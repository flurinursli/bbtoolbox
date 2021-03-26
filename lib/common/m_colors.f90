MODULE m_colors

  ! Purpose:
  !   to provide facilities to alter font color, background color and style of a string by adding ANSI/VT100 instructions. These
  !   require compatiable terminals and/or emulators (e.g. "less -R")
  !
  ! Revisions:
  !     Date                    Description of change
  !     ====                    =====================
  !   25/03/21                  original version
  !

  USE, NON_INTRINSIC :: m_precisions

  IMPLICIT none

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  PRIVATE

  PUBLIC :: colorize

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  CHARACTER(2), PARAMETER :: cs = ACHAR(27) // '['
  CHARACTER(1), PARAMETER :: ce = 'm', join = ';'
  CHARACTER(4), PARAMETER :: clearall = cs // '0' // ce

  ! foreground color codes
  CHARACTER(7), PARAMETER, DIMENSION(2,9) :: colors_foreground = RESHAPE(['black  ','30     ',      &
                                                                          'red    ','31     ',      &
                                                                          'green  ','32     ',      &
                                                                          'yellow ','33     ',      &
                                                                          'blue   ','34     ',      &
                                                                          'magenta','35     ',      &
                                                                          'cyan   ','36     ',      &
                                                                          'white  ','37     ',      &
                                                                          'default','39     '], [2,9])

  ! background color codes
  CHARACTER(7), PARAMETER, DIMENSION(2,9) :: colors_background = RESHAPE(['black  ','40     ',      &
                                                                          'red    ','41     ',      &
                                                                          'green  ','42     ',      &
                                                                          'yellow ','43     ',      &
                                                                          'blue   ','44     ',      &
                                                                          'magenta','45     ',      &
                                                                          'cyan   ','46     ',      &
                                                                          'white  ','47     ',      &
                                                                          'default','49     '], [2,9])

  ! formatting styles
  CHARACTER(13), PARAMETER, DIMENSION(2,6) :: styles = RESHAPE(['bold         ','1            ',      &
                                                                'italics      ','2            ',      &
                                                                'underlined   ','4            ',      &
                                                                'blink        ','5            ',      &
                                                                'reverse      ','7            ',      &
                                                                'strikethrough','9            '], [2,6])

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  CONTAINS

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    FUNCTION colorize(string, foreground, background, style)

      ! Purpose:
      !   to modify a string by adding ANSI/VT100 instructions for (foreground/background) color and formatting. If desired color/
      !   style is not available, the input string is returned unmodified.
      !
      !   Example: string = colorize(string, foreground = 'red', style = 'bold')
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   25/03/21                  original version
      !

      CHARACTER(*),                       INTENT(IN) :: string
      CHARACTER(*),             OPTIONAL, INTENT(IN) :: foreground, background, style
      CHARACTER(:), ALLOCATABLE                      :: colorize, buffer, set

      !-----------------------------------------------------------------------------------------------------------------------------

      buffer = cs

      set = match(colors_foreground, foreground)

      IF (set .ne. '') buffer = buffer // set

      set = match(colors_background, background)

      IF (set .ne. '') buffer = buffer // join // set

      set = match(styles, style)

      IF (set .ne. '') buffer = buffer // join // set

      buffer = buffer // ce

      IF (buffer .ne. (cs // ce)) THEN
        colorize = buffer // string // ' ' // clearall
      ELSE
        colorize = string
      ENDIF

      !=============================================================================================================================

      CONTAINS

        FUNCTION match(dict, key)

          CHARACTER(*),             DIMENSION(:,:),           INTENT(IN) :: dict
          CHARACTER(*),                             OPTIONAL, INTENT(IN) :: key
          CHARACTER(:), ALLOCATABLE                                      :: match
          INTEGER(i32)                                                   :: i

          !-------------------------------------------------------------------------------------------------------------------------

          match = ''

          IF (PRESENT(key)) THEN
            DO i = 1, SIZE(dict, 2)
              IF (INDEX(dict(1, i), key) .ne. 0) match = TRIM(dict(2, i))
            ENDDO
          ENDIF

        END FUNCTION match

    END FUNCTION colorize

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *


END MODULE m_colors
