MODULE m_toolbox

  USE, NON_INTRINSIC :: m_precisions
  USE, NON_INTRINSIC :: m_logfile
  USE, NON_INTRINSIC :: m_parser
  USE, NON_INTRINSIC :: m_strings

  IMPLICIT none

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  PRIVATE

  PUBLIC :: input, read_input_file, echo_input, broadcast, watch_start, watch_stop

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  INTERFACE assert
    MODULE PROCEDURE assert_i32, assert_r32
  END INTERFACE assert

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  TYPE :: src
    CHARACTER(8)   :: type
    CHARACTER(256) :: file
    INTEGER(i32)   :: seed, samples
    LOGICAL        :: is_point = .false., add_roughness = .true., add_rik = .true.
    REAL(r32)      :: x, y, z, m0, strike, dip, rake, freq, azimuth
    REAL(r32)      :: roughness, correlation, l0, aparam, vrfact
  END TYPE src

  TYPE :: rcv
    CHARACTER(32) :: file
    INTEGER(i32)  :: attenuation, velocity
    REAL(r32)     :: x, y, z
  END TYPE rcv

  TYPE :: att
    REAL(r32), ALLOCATABLE, DIMENSION(:) ::  gpp, gps, gss, b, lcut, hcut
  END TYPE att

  TYPE :: mdl
    REAL(r32), ALLOCATABLE, DIMENSION(:) :: vp, vs, rho, depth
  END TYPE mdl

  TYPE :: io
    CHARACTER(8)   :: format
    CHARACTER(16)  :: variable
    CHARACTER(256) :: folder
  END TYPE io

  TYPE :: hf
    CHARACTER(4) :: model
    INTEGER(i32) :: seed, samples
    LOGICAL      :: add_coherency = .true.
    REAL(r32)    :: fmax, matching, bandwidth, alpha, threshold
  END TYPE hf

  TYPE :: adv
    INTEGER(i32) :: pmw = 4, avecuts = 2, sheets = 0, waves = 2, verbose = 1
    REAL(r32)    :: vrfact = 0.85_r32
  END TYPE adv

  TYPE dict
    TYPE(adv)                            :: advanced
    TYPE(hf)                             :: coda
    TYPE(io)                             :: input, output
    TYPE(src)                            :: source
    TYPE(rcv), ALLOCATABLE, DIMENSION(:) :: receiver
    TYPE(att), ALLOCATABLE, DIMENSION(:) :: attenuation
    TYPE(mdl), ALLOCATABLE, DIMENSION(:) :: velocity
  END TYPE dict

  ! sort of dictionary contaning all input parameters
  TYPE(dict) :: input

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  CONTAINS

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE read_input_file(ok)

      ! Purpose:
      !   to read and check the bbtolbox input file. Name of input file is retrieved from command line. On exit, "ok" is not zero if
      !   an error occurred.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   08/03/21                  original version
      !

      INTEGER(i32),                             INTENT(OUT) :: ok
      CHARACTER(64)                                         :: fo
      CHARACTER(:), ALLOCATABLE                             :: str
      INTEGER(i32)                                          :: lu, n, i, p
      REAL(r32)                                             :: z, vp, vs, rho
      REAL(r32),    ALLOCATABLE, DIMENSION(:,:)             :: fbands

      !-----------------------------------------------------------------------------------------------------------------------------

      ok = 0

      CALL get_command_argument(1, fo, status = ok)                        !< get input file name from command line

      IF (ok .ne. 0) THEN
        CALL report_error('Could not read from command line')
        RETURN
      ENDIF

      OPEN(newunit = lu, file = TRIM(fo), status = 'old', form = 'formatted', access = 'sequential', action = 'read', iostat = ok)

      IF (ok .ne. 0) THEN
        CALL report_error('Error while opening file' + TRIM(fo))
        RETURN
      ENDIF

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! --------------------------------------------------- attenuation ------------------------------------------------------------

      n = recurrences(ok, lu, 'attenuation', key = 'attenuation', com = '#')                        !< number of attenuation models
      CALL missing_arg(ok, is_empty(n) .or. (n .le. 0), 'No attenuation model(s) found in input file')

      IF (ok .ne. 0) RETURN

      ALLOCATE(input%attenuation(n))


      DO i = 1, n

        CALL parse(ok, input%attenuation(i)%gpp, lu, 'gpp', ['[', ']'], ',', 'attenuation', nkey = i, com = '#')     !< gpp
        CALL missing_arg(ok, ALL(is_empty(input%attenuation(i)%gpp)), 'Argument "gpp" for attenuation model #' + num2char(i) +  &
        'not found')

        IF (ok .ne. 0) RETURN

        CALL parse(ok, input%attenuation(i)%gps, lu, 'gps', ['[', ']'], ',', 'attenuation', nkey = i, com = '#')     !< gps
        CALL missing_arg(ok, ALL(is_empty(input%attenuation(i)%gps)), 'Argument "gps" for attenuation model #' + num2char(i) +  &
        'not found')

        IF (ok .ne. 0) RETURN

        CALL parse(ok, input%attenuation(i)%gss, lu, 'gss', ['[', ']'], ',', 'attenuation', nkey = i, com = '#')     !< gss
        CALL missing_arg(ok, ALL(is_empty(input%attenuation(i)%gss)), 'Argument "gss" for attenuation model #' + num2char(i) +  &
        'not found')

        IF (ok .ne. 0) RETURN

        CALL parse(ok, input%attenuation(i)%b, lu, 'b', ['[', ']'], ',', 'attenuation', nkey = i, com = '#')     !< b
        CALL missing_arg(ok, ALL(is_empty(input%attenuation(i)%b)), 'Argument "b" for attenuation model #' + num2char(i) +  &
        'not found')

        IF (ok .ne. 0) RETURN

        CALL parse(ok, fbands, lu, 'frequency', ['[', ']'], ',', ';', 'attenuation', nkey = i, com = '#')
        CALL missing_arg(ok, ALL(is_empty(fbands)), 'Argument "frequency" for attenuation model #' + num2char(i) + 'not found')

        IF (ok .ne. 0) RETURN

        input%attenuation(i)%lcut = fbands(1, :)
        input%attenuation(i)%hcut = fbands(2, :)

      ENDDO

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! ----------------------------------------------------- velocity -------------------------------------------------------------

      n = recurrences(ok, lu, 'layer', key = 'layer', com = '#')                        !< number of layers
      CALL missing_arg(ok, is_empty(n) .or. (n .le. 0), 'No velocity model(s) found in input file')

      IF (ok .ne. 0) RETURN

      p = 0

      ! find number of velocity models by counting how many times a "layer" command without "depth" field is present
      DO i = 1, n
        CALL parse(ok, z, lu, 'depth', ['=', ','], 'layer', nkey = i, com = '#')
        IF (is_empty(z)) p = p + 1
      ENDDO

      ALLOCATE(input%velocity(p))

      p = 0

      DO i = 1, n

        CALL parse(ok, z, lu, 'depth', ['=', ','], 'layer', nkey = i, com = '#')

        IF (is_empty(z)) THEN
          p = p + 1
          z = 0._r32
        ENDIF

        CALL parse(ok, vp, lu, 'vp', ['=', ','], 'layer', nkey = i, com = '#')     !< vp
        CALL missing_arg(ok, is_empty(vp), 'Argument "vp" for layer keyword #' + num2char(i) + 'not found')

        IF (ok .ne. 0) RETURN

        CALL parse(ok, vs, lu, 'vs', ['=', ','], 'layer', nkey = i, com = '#')     !< vs
        CALL missing_arg(ok, is_empty(vs), 'Argument "vs" for layer keyword #' + num2char(i) + 'not found')

        IF (ok .ne. 0) RETURN

        CALL parse(ok, rho, lu, 'rho', ['=', ','], 'layer', nkey = i, com = '#')    !< rho
        CALL missing_arg(ok, is_empty(rho), 'Argument "rho" for layer keyword #' + num2char(i) + 'not found')

        IF (ok .ne. 0) RETURN

        IF (.not.ALLOCATED(input%velocity(p)%depth)) THEN
          input%velocity(p)%depth = [z]
          input%velocity(p)%vp    = [vp]
          input%velocity(p)%vs    = [vs]
          input%velocity(p)%rho   = [rho]
        ELSE
          input%velocity(p)%depth = [input%velocity(p)%depth, z]
          input%velocity(p)%vp    = [input%velocity(p)%vp, vp]
          input%velocity(p)%vs    = [input%velocity(p)%vs, vs]
          input%velocity(p)%rho   = [input%velocity(p)%rho, rho]
        ENDIF

      ENDDO

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! ----------------------------------------------------- receivers ------------------------------------------------------------

      n = recurrences(ok, lu, 'rec', key = 'rec', com = '#')                        !< number of receivers

      CALL missing_arg(ok, is_empty(n) .or. (n .le. 0), 'No stations found in input file')

      IF (ok .ne. 0) RETURN

      ALLOCATE(input%receiver(n))

      ! loop over stations to store parameters
      DO i = 1, n

        CALL parse(ok, str, lu, 'file', ["'", "'"], 'rec', nkey = i, com = '#')
        CALL missing_arg(ok, is_empty(str), 'Argument "file" for receiver #' + num2char(i) + ' not found')

        IF (ok .ne. 0) RETURN

        input%receiver(i)%file = TRIM(str)      !< file

        CALL parse(ok, input%receiver(i)%velocity, lu, 'velocity', ['=', ','], 'rec', nkey = i, com = '#')   !< velocity
        CALL missing_arg(ok, is_empty(input%receiver(i)%velocity), 'Argument "velocity" for receiver #' + num2char(i) +  &
        ' not found')

        IF (ok .ne. 0) RETURN

        CALL parse(ok, input%receiver(i)%attenuation, lu, 'attenuation', ['=', ','], 'rec', nkey = i, com = '#')  !< attenuation
        CALL missing_arg(ok, is_empty(input%receiver(i)%attenuation), 'Argument "attenuation" for receiver #' + num2char(i) +  &
        ' not found')

        IF (ok .ne. 0) RETURN

        CALL parse(ok, input%receiver(i)%x, lu, 'x', ['=', ','], 'rec', nkey = i, com = '#')       !< x
        CALL missing_arg(ok, is_empty(input%receiver(i)%x), 'Argument "x" for receiver #' + num2char(i) + ' not found')

        IF (ok .ne. 0) RETURN

        CALL parse(ok, input%receiver(i)%y, lu, 'y', ['=', ','], 'rec', nkey = i, com = '#')       !< y
        CALL missing_arg(ok, is_empty(input%receiver(i)%y), 'Argument "y" for receiver #' + num2char(i) + ' not found')

        IF (ok .ne. 0) RETURN

        CALL parse(ok, input%receiver(i)%z, lu, 'z', ['=', ','], 'rec', nkey = i, com = '#')    !< z
        CALL missing_arg(ok, is_empty(input%receiver(i)%z), 'Argument "z" for receiver #' + num2char(i) + ' not found')

        IF (ok .ne. 0) RETURN

      ENDDO

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! ------------------------------------------------ coda parameters -----------------------------------------------------------

      CALL parse(ok, input%coda%fmax, lu, 'fmax', ['=', ','], 'coda', com = '#')               !< fmax
      CALL missing_arg(ok, is_empty(input%coda%fmax), 'Argument "fmax" for keyword "coda" not found')

      IF (ok .ne. 0) RETURN

      CALL parse(ok, input%coda%matching, lu, 'matching frequency', ['=', ','], 'coda', com = '#')    !< matching
      CALL missing_arg(ok, is_empty(input%coda%matching), 'Argument "matching frequency" for keyword "coda" not found')

      IF (ok .ne. 0) RETURN

      CALL parse(ok, input%coda%bandwidth, lu, 'bandwidth', ['=', ','], 'coda', com = '#')     !< bandwidth
      CALL missing_arg(ok, is_empty(input%coda%bandwidth), 'Argument "bandwidth" for keyword "coda" not found')

      IF (ok .ne. 0) RETURN

      CALL parse(ok, input%coda%seed, lu, 'seed', ['=', ','], 'coda', com = '#')     !< seed
      CALL missing_arg(ok, is_empty(input%coda%seed), 'Argument "seed" for keyword "coda" not found')

      IF (ok .ne. 0) RETURN

      CALL parse(ok, input%coda%samples, lu, 'samples', ['=', ','], 'coda', com = '#')    !< samples
      CALL missing_arg(ok, is_empty(input%coda%samples), 'Argument "samples" for keyword "coda" not found')

      IF (ok .ne. 0) RETURN

      CALL parse(ok, str, lu, 'model', ["'", "'"], 'coda', com = '#')
      CALL missing_arg(ok, .false., 'Argument "model" for keyword "coda" not found')

      IF (ok .ne. 0) RETURN

      input%coda%model = TRIM(str)       !< model

      IF (is_empty(input%coda%model)) input%coda%add_coherency = .false.

      IF (input%coda%add_coherency) THEN
        CALL parse(ok, input%coda%alpha, lu, 'alpha', ['=', ','], 'coda', com = '#')     !< alpha
        CALL missing_arg(ok, is_empty(input%coda%alpha), 'Argument "alpha" for keyword "coda" not found')

        IF (ok .ne. 0) RETURN

        CALL parse(ok, input%coda%threshold, lu, 'threshold', ['=', ','], 'coda', com = '#')     !< threshold
        CALL missing_arg(ok, is_empty(input%coda%threshold), 'Argument "threshold" for keyword "coda" not found')

        IF (ok .ne. 0) RETURN
      ENDIF

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! --------------------------------------------------- input files ------------------------------------------------------------

      CALL parse(ok, str, lu, 'folder', ["'", "'"], 'input', com = '#')
      CALL missing_arg(ok, is_empty(str), 'Argument "folder" for keyword "input" not found')

      IF (ok .ne. 0) RETURN

      input%input%folder = TRIM(str)     !< folder

      CALL parse(ok, str, lu, 'format', ["'", "'"], 'input', com = '#')
      CALL missing_arg(ok, is_empty(str), 'Argument "format" for keyword "input" not found')

      IF (ok .ne. 0) RETURN

      input%input%format = TRIM(str)     !< format

      CALL parse(ok, str, lu, 'variable', ["'", "'"], 'input', com = '#')
      CALL missing_arg(ok, is_empty(str), 'Argument "variable" for keyword "input" not found')

      IF (ok .ne. 0) RETURN

      input%input%variable = TRIM(str)     !< variable

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! -------------------------------------------------- output files ------------------------------------------------------------

      CALL parse(ok, str, lu, 'folder', ["'", "'"], 'output', com = '#')
      CALL missing_arg(ok, is_empty(str), 'Argument "folder" for keyword "output" not found')

      IF (ok .ne. 0) RETURN

      input%output%folder = TRIM(str)     !< folder

      CALL parse(ok, str, lu, 'variable', ["'", "'"], 'output', com = '#')
      CALL missing_arg(ok, is_empty(str), 'Argument "variable" for keyword "output" not found')

      IF (ok .ne. 0) RETURN

      input%output%variable = TRIM(str)     !< variable

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! ------------------------------------------------------ source --------------------------------------------------------------

      n = recurrences(ok, lu, 'source', key = 'source', com = '#')                        !< find number of sources

      IF (ok .ne. 0) THEN
        CALL report_error(parser_error(ok))
        RETURN
      ENDIF

      IF (n .ge. 2) THEN
        CALL report_error('Only one point-source can be specified')
        ok = 1
        RETURN
      ENDIF

      IF (n .eq. 1) input%source%is_point = .true.

      IF (input%source%is_point) THEN

        CALL parse(ok, input%source%x, lu, 'x', ['=', ','], 'source', com = '#')     !< x
        CALL missing_arg(ok, is_empty(input%source%x), 'Argument "x" for keyword "source" not found')

        IF (ok .ne. 0) RETURN

        CALL parse(ok, input%source%y, lu, 'y', ['=', ','], 'source', com = '#')     !< y
        CALL missing_arg(ok, is_empty(input%source%y), 'Argument "y" for keyword "source" not found')

        IF (ok .ne. 0) RETURN

        CALL parse(ok, input%source%z, lu, 'z', ['=', ','], 'source', com = '#')     !< z
        CALL missing_arg(ok, is_empty(input%source%z), 'Argument "z" for keyword "source" not found')

        IF (ok .ne. 0) RETURN

        CALL parse(ok, input%source%m0, lu, 'm0', ['=', ','], 'source', com = '#')     !< m0
        CALL missing_arg(ok, is_empty(input%source%m0), 'Argument "m0" for keyword "source" not found')

        IF (ok .ne. 0) RETURN

        CALL parse(ok, input%source%m0, lu, 'm0', ['=', ','], 'source', com = '#')     !< m0
        CALL missing_arg(ok, is_empty(input%source%m0), 'Argument "m0" for keyword "source" not found')

        IF (ok .ne. 0) RETURN

        CALL parse(ok, input%source%strike, lu, 'strike', ['=', ','], 'source', com = '#')     !< strike
        CALL missing_arg(ok, is_empty(input%source%strike), 'Argument "strike" for keyword "source" not found')

        IF (ok .ne. 0) RETURN

        CALL parse(ok, input%source%dip, lu, 'dip', ['=', ','], 'source', com = '#')     !< dip
        CALL missing_arg(ok, is_empty(input%source%dip), 'Argument "dip" for keyword "source" not found')

        IF (ok .ne. 0) RETURN

        CALL parse(ok, input%source%rake, lu, 'rake', ['=', ','], 'source', com = '#')     !< rake
        CALL missing_arg(ok, is_empty(input%source%rake), 'Argument "rake" for keyword "source" not found')

        IF (ok .ne. 0) RETURN

        CALL parse(ok, input%source%freq, lu, 'freq', ['=', ','], 'source', com = '#')     !< freq
        CALL missing_arg(ok, is_empty(input%source%freq), 'Argument "freq" for keyword "source" not found')

        IF (ok .ne. 0) RETURN

        CALL parse(ok, str, lu, 'type', ["'", "'"], 'source', com = '#')     !< type
        CALL missing_arg(ok, is_empty(str), 'Argument "type" for keyword "source" not found')

        IF (ok .ne. 0) RETURN

        input%source%type = TRIM(str)

        input%source%add_rik = .false.
        input%source%add_roughness = .false.

      ENDIF

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! ----------------------------------------------------- rupture --------------------------------------------------------------

      n = recurrences(ok, lu, 'rupture', key = 'rupture', com = '#')

      IF ((n .gt. 0) .and. input%source%is_point) THEN
        CALL report_error('Either point-source or extended-source can be specified')
        ok = 1
        RETURN
      ENDIF

      IF (.not.input%source%is_point) THEN

        CALL missing_arg(ok, is_empty(n) .or. (n .le. 0), 'No source found in input file')

        IF (ok .ne. 0) RETURN

        IF (n .ge. 2) THEN
          CALL report_error('Only one rupture model can be specified')
          ok = 1
          RETURN
        ENDIF

        CALL parse(ok, str, lu, 'file', ["'", "'"], 'rupture', com = '#')     !< file
        CALL missing_arg(ok, is_empty(str), 'Argument "file" for keyword "rupture" not found')

        IF (ok .ne. 0) RETURN

        input%source%file = TRIM(str)

        CALL parse(ok, str, lu, 'type', ["'", "'"], 'rupture', com = '#')     !< type
        CALL missing_arg(ok, is_empty(str), 'Argument "type" for keyword "rupture" not found')

        IF (ok .ne. 0) RETURN

        input%source%type = TRIM(str)

        CALL parse(ok, input%source%x, lu, 'x', ['=', ','], 'rupture', com = '#')     !< x
        CALL missing_arg(ok, is_empty(input%source%x), 'Argument "x" for keyword "rupture" not found')

        IF (ok .ne. 0) RETURN

        CALL parse(ok, input%source%y, lu, 'y', ['=', ','], 'rupture', com = '#')     !< y
        CALL missing_arg(ok, is_empty(input%source%y), 'Argument "y" for keyword "rupture" not found')

        IF (ok .ne. 0) RETURN

        CALL parse(ok, input%source%azimuth, lu, 'azimuth', ['=', ','], 'rupture', com = '#')     !< azimuth
        CALL missing_arg(ok, is_empty(input%source%azimuth), 'Argument "azimuth" for keyword "rupture" not found')

        IF (ok .ne. 0) RETURN

        CALL parse(ok, input%source%roughness, lu, 'roughness', ['=', ','], 'rupture', com = '#')     !< roughness
        CALL missing_arg(ok, .false., '')

        IF (ok .ne. 0) RETURN

        IF (is_empty(input%source%roughness)) input%source%add_roughness = .false.

        CALL parse(ok, input%source%correlation, lu, 'corr', ['=', ','], 'rupture', com = '#')     !< correlation
        CALL missing_arg(ok, .false., '')

        IF (ok .ne. 0) RETURN

        IF (is_empty(input%source%correlation)) input%source%add_rik = .false.

        IF (input%source%add_rik) THEN

          CALL parse(ok, input%source%l0, lu, 'l0', ['=', ','], 'rupture', com = '#')     !< l0
          CALL missing_arg(ok, is_empty(input%source%l0), 'Argument "l0" for keyword "rupture" not found')

          IF (ok .ne. 0) RETURN

          CALL parse(ok, input%source%aparam, lu, 'aparam', ['=', ','], 'rupture', com = '#')     !< aparam
          CALL missing_arg(ok, is_empty(input%source%aparam), 'Argument "aparam" for keyword "rupture" not found')

          IF (ok .ne. 0) RETURN

          CALL parse(ok, input%source%vrfact, lu, 'vrfact', ['=', ','], 'rupture', com = '#')     !< vrfact
          CALL missing_arg(ok, is_empty(input%source%vrfact), 'Argument "vrfact" for keyword "rupture" not found')

          IF (ok .ne. 0) RETURN

          CALL parse(ok, input%source%seed, lu, 'seed', ['=', ','], 'rupture', com = '#')     !< seed
          CALL missing_arg(ok, is_empty(input%source%seed), 'Argument "seed" for keyword "rupture" not found')

          IF (ok .ne. 0) RETURN

          CALL parse(ok, input%source%samples, lu, 'samples', ['=', ','], 'rupture', com = '#')     !< samples
          CALL missing_arg(ok, is_empty(input%source%samples), 'Argument "samples" for keyword "rupture" not found')

          IF (ok .ne. 0) RETURN

        ENDIF

      ENDIF

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! ----------------------------------------------------- advanced -------------------------------------------------------------

      CALL parse(ok, z, lu, 'pmw', ['=', ','], 'advanced', com = '#')     !< pmw
      CALL missing_arg(ok, .false., '')

      IF (ok .ne. 0) RETURN

      IF (.not.is_empty(z)) input%advanced%pmw = z

      CALL parse(ok, z, lu, 'vrfact', ['=', ','], 'advanced', com = '#')     !< vrfact
      CALL missing_arg(ok, .false., '')

      IF (ok .ne. 0) RETURN

      IF (.not.is_empty(z)) input%advanced%vrfact = z

      CALL parse(ok, p, lu, 'avecuts', ['=', ','], 'advanced', com = '#')     !< avecuts
      CALL missing_arg(ok, .false., '')

      IF (ok .ne. 0) RETURN

      IF (.not.is_empty(p)) input%advanced%avecuts = p

      CALL parse(ok, p, lu, 'sheets', ['=', ','], 'advanced', com = '#')     !< sheets
      CALL missing_arg(ok, .false., '')

      IF (ok .ne. 0) RETURN

      IF (.not.is_empty(p)) input%advanced%sheets = p

      CALL parse(ok, p, lu, 'waves', ['=', ','], 'advanced', com = '#')     !< waves
      CALL missing_arg(ok, .false., '')

      IF (ok .ne. 0) RETURN

      IF (.not.is_empty(p)) input%advanced%waves = p

      CALL parse(ok, p, lu, 'verbose', ['=', ','], 'advanced', com = '#')     !< verbose
      CALL missing_arg(ok, .false., '')

      IF (ok .ne. 0) RETURN

      IF (.not.is_empty(p)) input%advanced%verbose = p

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! -------------------------------------- check whether all input parameters make sense ---------------------------------------

      DO i = 1, SIZE(input%attenuation)
        IF (assert(input%attenuation(i)%gpp, 0._r32) .or. assert(input%attenuation(i)%gps, 0._r32) .or.   &
            assert(input%attenuation(i)%gss, 0._r32) .or. assert(input%attenuation(i)%b, 0._r32)) THEN
          CALL report_error('All attenuation parameters must be larger than or equal to zero')
          ok = 1
          RETURN
        ENDIF

        IF (assert(input%attenuation(i)%lcut, 0._r32, strict = .true.) .or. assert(input%attenuation(i)%hcut, 0._r32,  &
            strict = .true.)) THEN
          CALL report_error('All attenuation frequency values must be larger than zero')
          ok = 1
          RETURN
        ENDIF

        IF (assert(input%attenuation(i)%hcut - input%attenuation(i)%lcut, 0._r32, strict = .true.)) THEN
          CALL report_error('Upper limit of attenuation frequency bands cannot be lower than lower limit')
          ok = 1
          RETURN
        ENDIF
      ENDDO

      DO i = 1, SIZE(input%velocity)
        IF (assert(input%velocity(i)%vp, 0._r32, strict = .true.) .or. assert(input%velocity(i)%vs, 0._r32, strict = .true.) .or. &
            assert(input%velocity(i)%rho, 0._r32, strict = .true.)) THEN
          CALL report_error('Physical parameters for velocity model(s) must be larger than zero')
          ok = 1
          RETURN
        ENDIF

        IF (assert(input%velocity(i)%depth, 0._r32)) THEN
          CALL report_error('Depth of layer(s) must be equal to or larger than zero')
          ok = 1
          RETURN
        ENDIF
      ENDDO

      IF (assert(input%receiver(:)%z, 0._r32)) THEN
        CALL report_error('Receiver(s) depth must be equal to or larger than zero')
        ok = 1
        RETURN
      ENDIF

      IF (assert(input%receiver(:)%attenuation, 1, SIZE(input%attenuation))) THEN
        CALL report_error('Wrong attenuation model specified for receiver(s)')
        ok = 1
        RETURN
      ENDIF

      IF (assert(input%receiver(:)%velocity, 1, SIZE(input%velocity))) THEN
        CALL report_error('Wrong velocity model specified for receiver(s)')
        ok = 1
        RETURN
      ENDIF

      IF (assert([input%coda%fmax, input%coda%matching, input%coda%bandwidth], 0._r32, strict = .true.)) THEN
        CALL report_error('Frequency parameters for coda must be larger than zero')
        ok = 1
        RETURN
      ENDIF

      IF (input%source%is_point) THEN
        IF (assert([input%source%z], 0._r32)) THEN
          CALL report_error('Source depth must be equal to or larger than zero')
          ok = 1
          RETURN
        ENDIF

        IF (assert([input%source%m0], 0._r32, strict = .true.)) THEN
          CALL report_error('Source moment must be larger than zero')
          ok = 1
          RETURN
        ENDIF

        IF (assert([input%source%freq], 0._r32, strict = .true.)) THEN
          CALL report_error('Source frequency parameter must be larger than zero')
          ok = 1
          RETURN
        ENDIF
      ENDIF

      ! missing checks for rik, roughness, coherency and advanced

    END SUBROUTINE read_input_file

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !==============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE missing_arg(ok, flag, text)

      INTEGER(i32), INTENT(INOUT) :: ok
      LOGICAL,      INTENT(IN)    :: flag
      CHARACTER(*), INTENT(IN)    :: text

      !-----------------------------------------------------------------------------------------------------------------------------

      IF (ok .ne. 0) THEN
        CALL report_error(parser_error(ok))
      ELSE
        IF (flag) THEN
          CALL report_error(text)
          ok = 1
        ENDIF
      ENDIF

    END SUBROUTINE missing_arg

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !==============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    FUNCTION assert_i32(par, v1, v2, strict) RESULT(assert)

      ! Purpose:
      !   To determine if variable "par" belongs to range [v1, v2]. Square brackets are replaced by round brackets if "strict=.true."
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   18/12/20                  original version
      !

      INTEGER(i32), DIMENSION(:),           INTENT(IN) :: par
      INTEGER(i32),               OPTIONAL, INTENT(IN) :: v1, v2
      LOGICAL,                    OPTIONAL, INTENT(IN) :: strict
      LOGICAL                                          :: flag, assert

      !-----------------------------------------------------------------------------------------------------------------------------

#include "assert_incl.f90"

    END FUNCTION assert_i32

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    FUNCTION assert_r32(par, v1, v2, strict) RESULT(assert)

      ! Purpose:
      !   To determine if variable "par" belongs to range [v1, v2]. Square brackets are replaced by round brackets if "strict=.true."
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   18/12/20                  original version
      !

      REAL(r32), DIMENSION(:),           INTENT(IN) :: par
      REAL(r32),               OPTIONAL, INTENT(IN) :: v1, v2
      LOGICAL,                 OPTIONAL, INTENT(IN) :: strict
      LOGICAL                                       :: flag, assert

      !-----------------------------------------------------------------------------------------------------------------------------

#include "assert_incl.f90"

    END FUNCTION assert_r32

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE echo_input()

      ! Purpose:
      !   to report input parameters to standard output.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   08/03/21                  original version
      !

      USE, INTRINSIC     :: iso_fortran_env, only: compiler_version
      USE, NON_INTRINSIC :: m_logfile

      INTEGER(i32) :: i, j

      !---------------------------------------------------------------------------------------------------------------------------------

      CALL update_log('***********************************************************************************************************')
      CALL update_log('BBTOOLBOX v1.0')
      CALL update_log('Compiled with ' + COMPILER_VERSION())
      CALL update_log('-----------------------------------------------------------------------------------------------------------')
      CALL update_log('Summary of input parameters', blankline = .false.)

      CALL update_log(num2char('Input folder LF seismograms', width=30, fill='.') +  &
                      num2char(TRIM(input%input%folder), width=79, justify = 'c') + '|')
      CALL update_log(num2char('Format LF seismograms', width=30, fill='.') +  &
                      num2char(TRIM(input%input%format), width=79, justify='c') + '|', blankline = .false.)
      CALL update_log(num2char('Variable LF seismograms', width=30, fill='.') +  &
                      num2char(TRIM(input%input%variable), width=79, justify='c') + '|', blankline = .false.)

      CALL update_log(num2char('Output folder HF seismograms', width=30, fill='.') +  &
                      num2char(TRIM(input%output%folder), width=79, justify = 'c') + '|')
      CALL update_log(num2char('Variable HF seismograms', width=30, fill='.') +  &
                      num2char(TRIM(input%output%variable), width=79, justify='c') + '|', blankline = .false.)

      CALL update_log(num2char('Coda parameters', width=30, fill='.') + num2char('Fmax', width=15, justify='r') + '|' +  &
                      num2char('Matching fr.', width=15, justify='r') + '|' + num2char('Bandwidth', width=15, justify='r') + '|' + &
                      num2char('Seed', width=15, justify='r') + '|' + num2char('Samples', width=15, justify='r') + '|')
      CALL update_log(num2char('', width=30) + &
                      num2char(input%coda%fmax, notation='f', precision=3, width=15, justify='r')      + '|' +  &
                      num2char(input%coda%matching, notation='f', precision=3, width=15, justify='r')  + '|' +  &
                      num2char(input%coda%bandwidth, notation='f', precision=3, width=15, justify='r') + '|' +  &
                      num2char(input%coda%seed, width=15, justify='r')                                 + '|' +  &
                      num2char(input%coda%samples, width=15, justify='r')                              + '|', blankline = .false.)

      IF (input%coda%add_coherency) THEN
        CALL update_log(num2char('Coherency parameters', width=30, fill='.') +     &
                        num2char('Model', width=15, justify='r')       + '|' +     &
                        num2char('Alpha', width=15, justify='r')       + '|' +     &
                        num2char('Threshold', width=15, justify='r')   + '|')
        CALL update_log(num2char('', width=30) +  &
                        num2char(input%coda%model, width=15, justify='r') + '|' + &
                        num2char(input%coda%alpha, notation='s', precision=3, width=15, justify='r') + '|' + &
                        num2char(input%coda%threshold, notation='f', precision=3, width=15, justify='r') + '|', blankline = .false.)
      ENDIF

      IF (input%source%is_point) THEN
        CALL update_log(num2char('Source parameters', width=30, fill='.') +     &
                        num2char('x', width=15, justify='r') + '|' +  &
                        num2char('y', width=15, justify='r') + '|' +  &
                        num2char('z', width=15, justify='r') + '|' +  &
                        num2char('Strike', width=15, justify='r') + '|' +  &
                        num2char('Dip', width=15, justify='r')    + '|')
        CALL update_log(num2char('', width=30) + &
                        num2char(input%source%x, notation='f', width=15, precision=3, justify='r') + '|' + &
                        num2char(input%source%y, notation='f', width=15, precision=3, justify='r') + '|' + &
                        num2char(input%source%z, notation='f', width=15, precision=3, justify='r') + '|' + &
                        num2char(input%source%strike, notation='f', width=15, precision=1, justify='r') + '|' + &
                        num2char(input%source%dip, notation='f', width=15, precision=1, justify='r')    + '|', blankline = .false.)
        CALL update_log(num2char('(continued)', width=30, fill='.', justify='c') +  &
                        num2char('Rake', width=15, justify='r')   + '|' + &
                        num2char('Moment', width=15, justify='r')     + '|' +  &
                        num2char('Type', width=15, justify='r') + '|' +  &
                        num2char('Frequency', width=15, justify='r') + '|', blankline = .false.)
        CALL update_log(num2char('', width=30) +  &
                        num2char(input%source%rake, notation='f', width=15, precision=1, justify='r')   + '|' + &
                        num2char(input%source%m0, notation='s', width=15, precision=3, justify='r')     + '|' + &
                        num2char(input%source%type, width=15, justify='r') + '|' + &
                        num2char(input%source%freq, notation='f', width=15, precision=3, justify='r') + '|', blankline = .false.)
      ELSE
        CALL update_log(num2char('Rupture file', width=30, fill='.')   +   &
                        num2char(TRIM(input%source%file), width=79, justify='c') + '|')
        CALL update_log(num2char('Rupture parameters', width=30, fill='.')  + &
                        num2char('x', width=15, justify='r') + '|' + &
                        num2char('y', width=15, justify='r') + '|' + &
                        num2char('Azimuth', width=15, justify='r') + '|' + &
                        num2char('Type', width=15, justify='r') + '|' + &
                        num2char('Roughness', width=15, justify='r') + '|', blankline = .false.)
        IF (input%source%add_roughness) THEN
          CALL update_log(num2char('', width=30) +  &
                          num2char(input%source%x, width=15, notation='f', precision=3, justify='r') + '|' + &
                          num2char(input%source%y, width=15, notation='f', precision=3, justify='r') + '|' + &
                          num2char(input%source%azimuth, width=15, notation='f', precision=3, justify='r') + '|' + &
                          num2char(input%source%type, width=15, justify='r') + '|' + &
                          num2char(input%source%roughness, width=15, notation='f', precision=1, justify='r') + '|',  &
                          blankline = .false.)
        ELSE
          CALL update_log(num2char('', width=30) +  &
                          num2char(input%source%x, width=15, notation='f', precision=3, justify='r') + '|' + &
                          num2char(input%source%y, width=15, notation='f', precision=3, justify='r') + '|' + &
                          num2char(input%source%azimuth, width=15, notation='f', precision=3, justify='r') + '|' + &
                          num2char(input%source%type, width=15, justify='r') + '|' + &
                          num2char('None', width=15, justify='r') + '|', blankline = .false.)
        ENDIF
        IF (input%source%add_rik) THEN
          CALL update_log(num2char('RIK parameters', width=30, fill='.') +  &
                          num2char('Correlation', width=15, justify='r') + '|' + &
                          num2char('L0', width=15, justify='r')      + '|' + &
                          num2char('Aparam', width=15, justify='r')  + '|' + &
                          num2char('Vrfact', width=15, justify='r')  + '|' + &
                          num2char('Seed', width=15, justify='r')    + '|' + &
                          num2char('Samples', width=15, justify='r') + '|')
          CALL update_log(num2char('', width=30) +  &
                          num2char(input%source%correlation, width=15, notation='f', precision=2, justify='r') + '|' + &
                          num2char(input%source%l0, width=15, notation='f', precision=1, justify='r')          + '|' + &
                          num2char(input%source%aparam, width=15, notation='f', precision=1, justify='r')      + '|' + &
                          num2char(input%source%vrfact, width=15, notation='f', precision=2, justify='r')      + '|' + &
                          num2char(input%source%seed, width=15, justify='r')                                   + '|' + &
                          num2char(input%source%samples, width=15, justify='r') + '|', blankline = .false.)
        ENDIF

      ENDIF

      DO i = 1, SIZE(input%velocity)
        CALL update_log(num2char('Velocity Model ' + num2char(i), width=30, fill='.') +   &
                        num2char('Vp', width=15, justify='r') + '|' + &
                        num2char('Vs', width=15, justify='r') + '|' + &
                        num2char('Rho', width=15, justify='r') + '|' + &
                        num2char('Depth', width=15, justify='r') + '|')
        DO j = 1, SIZE(input%velocity(i)%vp)
          CALL update_log(num2char('', width=30)  +  &
                          num2char(input%velocity(i)%vp(j), width=15, notation='f', precision=3, justify='r') + '|' + &
                          num2char(input%velocity(i)%vs(j), width=15, notation='f', precision=3, justify='r') + '|' + &
                          num2char(input%velocity(i)%rho(j), width=15, notation='f', precision=3, justify='r') + '|' + &
                          num2char(input%velocity(i)%depth(j), width=15, notation='f', precision=3, justify='r') + '|',  &
                          blankline = .false.)
        ENDDO
      ENDDO

      DO i = 1, SIZE(input%attenuation)
        CALL update_log(num2char('Attenuation Model ' + num2char(i), width=30, fill='.') +   &
                        num2char('Low Freq', width=15, justify='r') + '|' + &
                        num2char('High Freq', width=15, justify='r') + '|' + &
                        num2char('Gpp', width=15, justify='r') + '|' + &
                        num2char('Gps', width=15, justify='r') + '|' + &
                        num2char('Gss', width=15, justify='r') + '|' + &
                        num2char('b', width=15, justify='r') + '|')
        DO j = 1, SIZE(input%attenuation(i)%gpp)
          CALL update_log(num2char('', width=30)  +  &
                          num2char(input%attenuation(i)%lcut(j), width=15, notation='f', precision=1, justify='r') + '|' + &
                          num2char(input%attenuation(i)%hcut(j), width=15, notation='f', precision=1, justify='r') + '|' + &
                          num2char(input%attenuation(i)%gpp(j), width=15, notation='f', precision=3, justify='r') + '|' + &
                          num2char(input%attenuation(i)%gps(j), width=15, notation='f', precision=3, justify='r') + '|' + &
                          num2char(input%attenuation(i)%gss(j), width=15, notation='f', precision=3, justify='r') + '|' + &
                          num2char(input%attenuation(i)%b(j), width=15, notation='f', precision=3, justify='r') + '|',  &
                          blankline = .false.)
        ENDDO
      ENDDO

      CALL update_log(num2char('Receivers List', width=30, fill='.') + &
                      num2char('x', width=15, justify='r') + '|' +  &
                      num2char('y', width=15, justify='r') + '|' +  &
                      num2char('z', width=15, justify='r') + '|' +  &
                      num2char('File', width=15, justify='r') + '|' + &
                      num2char('Vel Model', width=15, justify='r') + '|' + &
                      num2char('Att Model', width=15, justify='r') + '|')

      DO i = 1, SIZE(input%receiver)
        CALL update_log(num2char('', width=30)  +   &
                        num2char(input%receiver(i)%x, width=15, notation='f', precision=3, justify='r') + '|' + &
                        num2char(input%receiver(i)%y, width=15, notation='f', precision=3, justify='r') + '|' + &
                        num2char(input%receiver(i)%z, width=15, notation='f', precision=3, justify='r') + '|' + &
                        num2char(TRIM(input%receiver(i)%file), width=15, justify='r') + '|' + &
                        num2char(input%receiver(i)%velocity, width=15, justify='r')   + '|' + &
                        num2char(input%receiver(i)%attenuation, width=15, justify='r') + '|', blankline = .false.)
      ENDDO

      CALL update_log(num2char('Advanced settings', width=30, fill='.') +   &
                      num2char('Pmw', width=15, justify = 'r') + '|' +   &
                      num2char('Avg cuts', width=15, justify = 'r')           + '|' +   &
                      num2char('Vr factor', width=15, justify = 'r')          + '|' +   &
                      num2char('Sheets', width=15, justify = 'r')             + '|' +   &
                      num2char('Waves', width=15, justify = 'r')              + '|' +   &
                      num2char('Verbose', width=15, justify = 'r')            + '|')

      CALL update_log(num2char('', width=30) + &
                      num2char(input%advanced%pmw, width=15, justify='r')     + '|' + &
                      num2char(input%advanced%avecuts, width=15, justify='r') + '|' + &
                      num2char(input%advanced%vrfact, width=15, notation='f', precision=2, justify='r') + '|' +  &
                      num2char(input%advanced%sheets, width=15, justify='r')  + '|' + &
                      num2char(input%advanced%waves, width=15, justify='r') + '|' + &
                      num2char(input%advanced%verbose, width=15, justify='r') + '|', blankline=.false.)


      CALL update_log('-----------------------------------------------------------------------------------------------------------')

    END SUBROUTINE echo_input

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE broadcast()

      ! Purpose:
      !   to broadcast all input parameters to all processes.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   08/03/21                  original version
      !

#ifdef MPI
      USE                :: mpi
      USE, NON_INTRINSIC :: m_precisions
      USE, NON_INTRINSIC :: m_logfile

      INTEGER(i32) :: ierr

      !-----------------------------------------------------------------------------------------------------------------------------

      ! coda
      CALL mpi_bcast(input%coda%fmax, 1, mpi_real, 0, mpi_comm_world, ierr)
      CALL mpi_bcast(input%coda%matching, 1, mpi_real, 0, mpi_comm_world, ierr)
      CALL mpi_bcast(input%coda%bandwidth, 1, mpi_real, 0, mpi_comm_world, ierr)
      CALL mpi_bcast(input%coda%seed, 1, mpi_int, 0, mpi_comm_world, ierr)
      CALL mpi_bcast(input%coda%samples, 1, mpi_int, 0, mpi_comm_world, ierr)
      CALL mpi_bcast(input%coda%model, 4, mpi_character, 0, mpi_comm_world, ierr)
      CALL mpi_bcast(input%coda%alpha, 1, mpi_real, 0, mpi_comm_world, ierr)
      CALL mpi_bcast(input%coda%threshold, 1, mpi_real, 0, mpi_comm_world, ierr)

      ! input
      CALL mpi_bcast(input%input%folder, 256, mpi_character, 0, mpi_comm_world, ierr)
      CALL mpi_bcast(input%input%format, 8, mpi_character, 0, mpi_comm_world, ierr)
      CALL mpi_bcast(input%input%variable, 16, mpi_character, 0, mpi_comm_world, ierr)

      ! output
      CALL mpi_bcast(input%output%folder, 256, mpi_character, 0, mpi_comm_world, ierr)
      CALL mpi_bcast(input%input%variable, 16, mpi_character, 0, mpi_comm_world, ierr)

      ! source
      CALL mpi_bcast(input%source%file, 256, mpi_character, 0, mpi_comm_world, ierr)
      CALL mpi_bcast(input%source%type, 8, mpi_character, 0, mpi_comm_world, ierr)

#endif

    END SUBROUTINE broadcast

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE watch_start(tictoc, comm)

      ! Purpose:
      !   To start the stopwatch. Timing is in double-precision. For mpi, if specific communicator handle not given, mpi_comm_world
      !   is used.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   08/03/21                  original version
      !

#ifdef MPI
      USE                :: mpi
#endif
      USE, NON_INTRINSIC :: m_precisions

      REAL(r64),              INTENT(OUT) :: tictoc                            !< initial time
      INTEGER(i32), OPTIONAL, INTENT(IN)  :: comm                              !< communicator handle (mpi-only)
      INTEGER(i32)                        :: ierr
      INTEGER(i64)                        :: t0, rate

      !-----------------------------------------------------------------------------------------------------------------------------

#ifdef MPI
      IF (PRESENT(comm)) THEN
        CALL mpi_barrier(comm, ierr)
      ELSE
        CALL mpi_barrier(mpi_comm_world, ierr)
      ENDIF
      tictoc = mpi_wtime()
#else
      CALL SYSTEM_CLOCK(t0, rate)
      tictoc = REAL(t0, r64) / rate
#endif

    END SUBROUTINE watch_start

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE watch_stop(tictoc, comm)

      ! Purpose:
      !   To stop the stopwatch and return elapsed time. Timing is in double-precision. For mpi, if specific communicator handle not
      !   given, mpi_comm_world is used.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   08/03/21                  original version
      !

#ifdef MPI
      USE                :: mpi
#endif
      USE, NON_INTRINSIC :: m_precisions

      REAL(r64),              INTENT(INOUT) :: tictoc                          !< elapsed time
      INTEGER(i32), OPTIONAL, INTENT(IN)    :: comm                            !< communicator handle
      INTEGER(i32)                          :: ierr
      INTEGER(i64)                          :: t0, rate

      !-----------------------------------------------------------------------------------------------------------------------------

#ifdef MPI
      IF (PRESENT(comm)) THEN
        CALL mpi_barrier(comm, ierr)
      ELSE
        CALL mpi_barrier(mpi_comm_world, ierr)
      ENDIF
      tictoc = mpi_wtime() - tictoc
#else
      CALL SYSTEM_CLOCK(t0, rate)
      tictoc = REAL(t0, r64) / rate - tictoc
#endif

    END SUBROUTINE watch_stop

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

END MODULE m_toolbox
