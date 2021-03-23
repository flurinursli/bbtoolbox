MODULE m_toolbox

  USE, INTRINSIC     :: iso_c_binding
  USE, NON_INTRINSIC :: m_precisions
  USE, NON_INTRINSIC :: m_logfile
  USE, NON_INTRINSIC :: m_parser
  USE, NON_INTRINSIC :: m_strings
#ifdef MPI
  USE                :: mpi
#endif

  IMPLICIT none

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  PRIVATE

  PUBLIC :: input
  PUBLIC :: read_input_file, missing_arg, watch_start, watch_stop, geo2utm, split_task


  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  INTERFACE assert
    MODULE PROCEDURE assert_i32, assert_r32
  END INTERFACE assert

  INTERFACE
    SUBROUTINE geo2utm(long, lat, long_0, lat_0, easting, northing) bind(C, name="fun_c")
      USE, INTRINSIC     :: iso_c_binding
      USE, NON_INTRINSIC :: m_precisions
      IMPLICIT none
      REAL(c_r32), VALUE       :: long, lat, long_0, lat_0
      REAL(c_r32), INTENT(OUT) :: easting, northing
    END SUBROUTINE
  END INTERFACE

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  REAL(r32), PARAMETER :: DEFAULT_VPGRAD = 0.01_r32, DEFAULT_VSGRAD = 0.01_r32

  TYPE :: src
    CHARACTER(8)   :: type = "Brune"
    CHARACTER(256) :: file
    INTEGER(i32)   :: seed, samples
    LOGICAL        :: is_point = .false., add_roughness = .true., add_rik = .true.
    REAL(r32)      :: x, y, lon, lat, z, m0, strike, dip, rake
    REAL(r32)      :: freq = 6.28_r32
    REAL(r32)      :: roughness, correlation, l0, aparam, vrfact
  END TYPE src

  TYPE :: rcv
    CHARACTER(32) :: file
    INTEGER(i32)  :: attenuation, velocity
    REAL(r32)     :: x, y, z, lon, lat
  END TYPE rcv

  TYPE :: att
    REAL(r32), ALLOCATABLE, DIMENSION(:) :: gpp, gps, gss, b, lcut, hcut
  END TYPE att

  TYPE :: mdl
    REAL(r32), ALLOCATABLE, DIMENSION(:) :: vp, vs, rho, depth, vpgrad, vsgrad
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
    REAL(r32)    :: fmax, matching, bandwidth
    REAL(r32)    :: alpha                   !< Luco-Wong coherence
    REAL(r32)    :: a, ak, f0, b            !< Harichandran-Vanmacke coherence
  END TYPE hf

  TYPE :: adv
    INTEGER(i32) :: pmw = 4, avecuts = 2, sheets = 0, waves = 2, verbose = 1
    REAL(r32)    :: vrfact = 0.85_r32
  END TYPE adv

  TYPE :: org
    LOGICAL   :: is_geo = .true.
    REAL(r32) :: lon, lat
  END TYPE org

  TYPE dict
    TYPE(org)                            :: origin
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

    SUBROUTINE read_input_file(ok, rank, ntasks)

      ! Purpose:
      !   to read and check the bbtolbox input file. Name of input file is retrieved from command line. On exit, "ok" is not zero if
      !   an error occurred.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   08/03/21                  original version
      !

      INTEGER(i32), INTENT(OUT) :: ok
      INTEGER(i32), INTENT(IN)  :: rank, ntasks
      INTEGER(i32)              :: ierr

      !-----------------------------------------------------------------------------------------------------------------------------

      IF (rank .eq. 0) CALL parse_input_file(ok)

#ifdef MPI
      CALL mpi_bcast(ok, 1, mpi_int, 0, mpi_comm_world, ierr)
#endif

      IF (ok .ne. 0) RETURN

#ifdef MPI
      CALL broadcast()
#endif

      IF (rank .eq. ntasks - 1) CALL echo()

    END SUBROUTINE read_input_file

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE parse_input_file(ok)

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
      REAL(r32)                                             :: z, lon, lat, vp, vs, rho, vpgrad, vsgrad, freq
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
      ! ----------------------------------------------------- origin ---------------------------------------------------------------

      CALL parse(ok, input%origin%lon, lu, 'lon', ['=',' '], 'origin', com='#')
      CALL missing_arg(ok, .false., '')

      IF (ok .ne. 0) RETURN

      CALL parse(ok, input%origin%lat, lu, 'lat', ['=',' '], 'origin', com='#')
      CALL missing_arg(ok, .false., '')

      IF (ok .ne. 0) RETURN

      IF (is_empty(input%origin%lon) .and. is_empty(input%origin%lat)) input%origin%is_geo = .false.

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
        CALL parse(ok, z, lu, 'depth', ['=', ' '], 'layer', nkey = i, com = '#')
        IF (is_empty(z)) p = p + 1
      ENDDO

      ALLOCATE(input%velocity(p))

      p = 0

      DO i = 1, n

        CALL parse(ok, z, lu, 'depth', ['=', ' '], 'layer', nkey = i, com = '#')

        IF (is_empty(z)) THEN
          z = 0._r32
          p = p + 1
        ENDIF

        CALL parse(ok, vpgrad, lu, 'vpgrad', ['=', ' '], 'layer', nkey = i, com = '#')
        CALL missing_arg(ok, .false., '')

        IF (ok .ne. 0) RETURN

        IF (is_empty(vpgrad)) vpgrad = DEFAULT_VPGRAD

        CALL parse(ok, vsgrad, lu, 'vsgrad', ['=', ' '], 'layer', nkey = i, com = '#')
        CALL missing_arg(ok, .false., '')

        IF (ok .ne. 0) RETURN

        IF (is_empty(vsgrad)) vsgrad = DEFAULT_VSGRAD

        CALL parse(ok, vp, lu, 'vp', ['=', ' '], 'layer', nkey = i, com = '#')     !< vp
        CALL missing_arg(ok, is_empty(vp), 'Argument "vp" for layer keyword #' + num2char(i) + 'not found')

        IF (ok .ne. 0) RETURN

        CALL parse(ok, vs, lu, 'vs', ['=', ' '], 'layer', nkey = i, com = '#')     !< vs
        CALL missing_arg(ok, is_empty(vs), 'Argument "vs" for layer keyword #' + num2char(i) + 'not found')

        IF (ok .ne. 0) RETURN

        CALL parse(ok, rho, lu, 'rho', ['=', ' '], 'layer', nkey = i, com = '#')    !< rho
        CALL missing_arg(ok, is_empty(rho), 'Argument "rho" for layer keyword #' + num2char(i) + 'not found')

        IF (ok .ne. 0) RETURN

        IF (.not.ALLOCATED(input%velocity(p)%depth)) THEN
          input%velocity(p)%vp  = [vp]
          input%velocity(p)%vs  = [vs]
          input%velocity(p)%rho = [rho]
          input%velocity(p)%vpgrad = [vpgrad]
          input%velocity(p)%vsgrad = [vsgrad]
          IF (.not.is_empty(z)) input%velocity(p)%depth  = [z]

        ELSE
          input%velocity(p)%vp  = [input%velocity(p)%vp, vp]
          input%velocity(p)%vs  = [input%velocity(p)%vs, vs]
          input%velocity(p)%rho = [input%velocity(p)%rho, rho]
          input%velocity(p)%vpgrad = [input%velocity(p)%vpgrad, vpgrad]
          input%velocity(p)%vsgrad = [input%velocity(p)%vsgrad, vsgrad]
          IF (.not.is_empty(z)) input%velocity(p)%depth = [input%velocity(p)%depth, z]
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

        CALL parse(ok, input%receiver(i)%velocity, lu, 'velocity', ['=', ' '], 'rec', nkey = i, com = '#')   !< velocity
        CALL missing_arg(ok, .false., '')

        IF (ok .ne. 0) RETURN

        IF (is_empty(input%receiver(i)%velocity)) input%receiver(i)%velocity = 1

        CALL parse(ok, input%receiver(i)%attenuation, lu, 'attenuation', ['=', ' '], 'rec', nkey = i, com = '#')  !< attenuation
        CALL missing_arg(ok, .false., '')

        IF (ok .ne. 0) RETURN

        IF (is_empty(input%receiver(i)%attenuation)) input%receiver(i)%attenuation = 1

        IF (input%origin%is_geo) THEN
          CALL parse(ok, input%receiver(i)%lon, lu, 'lon', ['=', ' '], 'rec', nkey = i, com = '#')       !< lon
          CALL missing_arg(ok, is_empty(input%receiver(i)%lon), 'Argument "lon" for receiver #' + num2char(i) + ' not found')

          IF (ok .ne. 0) RETURN

          CALL parse(ok, input%receiver(i)%lat, lu, 'lat', ['=', ' '], 'rec', nkey = i, com = '#')       !< lat
          CALL missing_arg(ok, is_empty(input%receiver(i)%lat), 'Argument "lat" for receiver #' + num2char(i) + ' not found')

          IF (ok .ne. 0) RETURN

          ! move from geographical to utm coordinates (y is "east-west" as lon, x is "north-south" as lat)
          CALL geo2utm(input%receiver(i)%lon, input%receiver(i)%lat, input%origin%lon, input%origin%lat, input%receiver(i)%y,  &
                       input%receiver(i)%x)

        ELSE
          CALL parse(ok, input%receiver(i)%x, lu, 'x', ['=', ' '], 'rec', nkey = i, com = '#')       !< x
          CALL missing_arg(ok, is_empty(input%receiver(i)%x), 'Argument "x" for receiver #' + num2char(i) + ' not found')

          IF (ok .ne. 0) RETURN

          CALL parse(ok, input%receiver(i)%y, lu, 'y', ['=', ' '], 'rec', nkey = i, com = '#')       !< y
          CALL missing_arg(ok, is_empty(input%receiver(i)%y), 'Argument "y" for receiver #' + num2char(i) + ' not found')

          IF (ok .ne. 0) RETURN
        ENDIF

        CALL parse(ok, input%receiver(i)%z, lu, 'z', ['=', ' '], 'rec', nkey = i, com = '#')    !< z
        CALL missing_arg(ok, is_empty(input%receiver(i)%z), 'Argument "z" for receiver #' + num2char(i) + ' not found')

        IF (ok .ne. 0) RETURN

      ENDDO

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! ------------------------------------------------ coda parameters -----------------------------------------------------------

      CALL parse(ok, input%coda%fmax, lu, 'fmax', ['=', ' '], 'coda', com = '#')               !< fmax
      CALL missing_arg(ok, is_empty(input%coda%fmax), 'Argument "fmax" for keyword "coda" not found')

      IF (ok .ne. 0) RETURN

      CALL parse(ok, input%coda%matching, lu, 'matching frequency', ['=', ' '], 'coda', com = '#')    !< matching
      CALL missing_arg(ok, is_empty(input%coda%matching), 'Argument "matching frequency" for keyword "coda" not found')

      IF (ok .ne. 0) RETURN

      CALL parse(ok, input%coda%bandwidth, lu, 'bandwidth', ['=', ' '], 'coda', com = '#')     !< bandwidth
      CALL missing_arg(ok, is_empty(input%coda%bandwidth), 'Argument "bandwidth" for keyword "coda" not found')

      IF (ok .ne. 0) RETURN

      CALL parse(ok, input%coda%seed, lu, 'seed', ['=', ' '], 'coda', com = '#')     !< seed
      CALL missing_arg(ok, is_empty(input%coda%seed), 'Argument "seed" for keyword "coda" not found')

      IF (ok .ne. 0) RETURN

      CALL parse(ok, input%coda%samples, lu, 'samples', ['=', ' '], 'coda', com = '#')    !< samples
      CALL missing_arg(ok, is_empty(input%coda%samples), 'Argument "samples" for keyword "coda" not found')

      IF (ok .ne. 0) RETURN

      CALL parse(ok, str, lu, 'model', ["'", "'"], 'coda', com = '#')
      CALL missing_arg(ok, .false., 'Argument "model" for keyword "coda" not found')

      IF (ok .ne. 0) RETURN

      input%coda%model = TRIM(str)       !< model

      IF (is_empty(input%coda%model)) input%coda%add_coherency = .false.

      IF (input%coda%add_coherency) THEN

        IF (input%coda%model .eq. 'lw') THEN

          CALL parse(ok, input%coda%alpha, lu, 'alpha', ['=', ' '], 'coda', com = '#')     !< alpha
          CALL missing_arg(ok, is_empty(input%coda%alpha), 'Argument "alpha" for keyword "coda" not found')

          IF (ok .ne. 0) RETURN

        ELSEIF (input%coda%model .eq. 'hv') THEN

          CALL parse(ok, input%coda%a, lu, 'a', ['=', ' '], 'coda', com = '#')     !< a
          CALL missing_arg(ok, is_empty(input%coda%a), 'Argument "a" for keyword "coda" not found')

          IF (ok .ne. 0) RETURN

          CALL parse(ok, input%coda%ak, lu, 'ak', ['=', ' '], 'coda', com = '#')     !< ak
          CALL missing_arg(ok, is_empty(input%coda%ak), 'Argument "ak" for keyword "coda" not found')

          IF (ok .ne. 0) RETURN

          CALL parse(ok, input%coda%b, lu, 'b', ['=', ' '], 'coda', com = '#')     !< b
          CALL missing_arg(ok, is_empty(input%coda%b), 'Argument "b" for keyword "coda" not found')

          IF (ok .ne. 0) RETURN

          CALL parse(ok, input%coda%f0, lu, 'f0', ['=', ' '], 'coda', com = '#')     !< f0
          CALL missing_arg(ok, is_empty(input%coda%f0), 'Argument "f0" for keyword "coda" not found')

          IF (ok .ne. 0) RETURN

        ENDIF

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

        IF (input%origin%is_geo) THEN

          CALL parse(ok, input%source%lon, lu, 'lon', ['=', ' '], 'source', com = '#')     !< lon
          CALL missing_arg(ok, is_empty(input%source%lon), 'Argument "lon" for keyword "source" not found')

          IF (ok .ne. 0) RETURN

          CALL parse(ok, input%source%lat, lu, 'lat', ['=', ' '], 'source', com = '#')     !< lat
          CALL missing_arg(ok, is_empty(input%source%lat), 'Argument "lat" for keyword "source" not found')

          IF (ok .ne. 0) RETURN

          ! move from geographical to utm coordinates (y is "east-west" as lon, x is "north-south" as lat)
          CALL geo2utm(input%source%lon, input%source%lat, input%origin%lon, input%origin%lat, input%source%y, input%source%x)

        ELSE

          CALL parse(ok, input%source%x, lu, 'x', ['=', ' '], 'source', com = '#')     !< x
          CALL missing_arg(ok, is_empty(input%source%x), 'Argument "x" for keyword "source" not found')

          IF (ok .ne. 0) RETURN

          CALL parse(ok, input%source%y, lu, 'y', ['=', ' '], 'source', com = '#')     !< y
          CALL missing_arg(ok, is_empty(input%source%y), 'Argument "y" for keyword "source" not found')

          IF (ok .ne. 0) RETURN

        ENDIF

        CALL parse(ok, input%source%z, lu, 'z', ['=', ' '], 'source', com = '#')     !< z
        CALL missing_arg(ok, is_empty(input%source%z), 'Argument "z" for keyword "source" not found')

        IF (ok .ne. 0) RETURN

        CALL parse(ok, input%source%m0, lu, 'm0', ['=', ' '], 'source', com = '#')     !< m0
        CALL missing_arg(ok, is_empty(input%source%m0), 'Argument "m0" for keyword "source" not found')

        IF (ok .ne. 0) RETURN

        CALL parse(ok, input%source%strike, lu, 'strike', ['=', ' '], 'source', com = '#')     !< strike
        CALL missing_arg(ok, is_empty(input%source%strike), 'Argument "strike" for keyword "source" not found')

        IF (ok .ne. 0) RETURN

        CALL parse(ok, input%source%dip, lu, 'dip', ['=', ' '], 'source', com = '#')     !< dip
        CALL missing_arg(ok, is_empty(input%source%dip), 'Argument "dip" for keyword "source" not found')

        IF (ok .ne. 0) RETURN

        CALL parse(ok, input%source%rake, lu, 'rake', ['=', ' '], 'source', com = '#')     !< rake
        CALL missing_arg(ok, is_empty(input%source%rake), 'Argument "rake" for keyword "source" not found')

        IF (ok .ne. 0) RETURN

        CALL parse(ok, freq, lu, 'freq', ['=', ' '], 'source', com = '#')     !< freq
        CALL missing_arg(ok, .false., '')

        IF (ok .ne. 0) RETURN

        IF (.not.is_empty(freq)) input%source%freq = freq

        CALL parse(ok, str, lu, 'type', ["'", "'"], 'source', com = '#')     !< type
        CALL missing_arg(ok, .false., '')

        IF (ok .ne. 0) RETURN

        IF (.not.is_empty(str)) input%source%type = TRIM(str)

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

        CALL parse(ok, freq, lu, 'freq', ['=', ' '], 'rupture', com = '#')     !< freq
        CALL missing_arg(ok, .false., '')

        IF (ok .ne. 0) RETURN

        IF (.not.is_empty(freq)) input%source%freq = freq

        CALL parse(ok, str, lu, 'type', ["'", "'"], 'rupture', com = '#')     !< type
        CALL missing_arg(ok, .false., '')

        IF (ok .ne. 0) RETURN

        IF (.not.is_empty(str)) input%source%type = TRIM(str)
        CALL parse(ok, input%source%roughness, lu, 'roughness', ['=', ' '], 'rupture', com = '#')     !< roughness
        CALL missing_arg(ok, .false., '')

        IF (ok .ne. 0) RETURN

        IF (is_empty(input%source%roughness)) input%source%add_roughness = .false.

        CALL parse(ok, input%source%correlation, lu, 'corr', ['=', ' '], 'rupture', com = '#')     !< correlation
        CALL missing_arg(ok, .false., '')

        IF (ok .ne. 0) RETURN

        IF (is_empty(input%source%correlation)) input%source%add_rik = .false.

        IF (input%source%add_rik) THEN

          CALL parse(ok, input%source%l0, lu, 'l0', ['=', ' '], 'rupture', com = '#')     !< l0
          CALL missing_arg(ok, is_empty(input%source%l0), 'Argument "l0" for keyword "rupture" not found')

          IF (ok .ne. 0) RETURN

          CALL parse(ok, input%source%aparam, lu, 'aparam', ['=', ' '], 'rupture', com = '#')     !< aparam
          CALL missing_arg(ok, is_empty(input%source%aparam), 'Argument "aparam" for keyword "rupture" not found')

          IF (ok .ne. 0) RETURN

          CALL parse(ok, input%source%vrfact, lu, 'vrfact', ['=', ' '], 'rupture', com = '#')     !< vrfact
          CALL missing_arg(ok, is_empty(input%source%vrfact), 'Argument "vrfact" for keyword "rupture" not found')

          IF (ok .ne. 0) RETURN

          CALL parse(ok, input%source%seed, lu, 'seed', ['=', ' '], 'rupture', com = '#')     !< seed
          CALL missing_arg(ok, is_empty(input%source%seed), 'Argument "seed" for keyword "rupture" not found')

          IF (ok .ne. 0) RETURN

          CALL parse(ok, input%source%samples, lu, 'samples', ['=', ' '], 'rupture', com = '#')     !< samples
          CALL missing_arg(ok, is_empty(input%source%samples), 'Argument "samples" for keyword "rupture" not found')

          IF (ok .ne. 0) RETURN

        ENDIF

      ENDIF

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
      ! ----------------------------------------------------- advanced -------------------------------------------------------------

      CALL parse(ok, z, lu, 'pmw', ['=', ' '], 'advanced', com = '#')     !< pmw
      CALL missing_arg(ok, .false., '')

      IF (ok .ne. 0) RETURN

      IF (.not.is_empty(z)) input%advanced%pmw = z

      CALL parse(ok, z, lu, 'vrfact', ['=', ' '], 'advanced', com = '#')     !< vrfact
      CALL missing_arg(ok, .false., '')

      IF (ok .ne. 0) RETURN

      IF (.not.is_empty(z)) input%advanced%vrfact = z

      CALL parse(ok, p, lu, 'avecuts', ['=', ' '], 'advanced', com = '#')     !< avecuts
      CALL missing_arg(ok, .false., '')

      IF (ok .ne. 0) RETURN

      IF (.not.is_empty(p)) input%advanced%avecuts = p

      CALL parse(ok, p, lu, 'sheets', ['=', ' '], 'advanced', com = '#')     !< sheets
      CALL missing_arg(ok, .false., '')

      IF (ok .ne. 0) RETURN

      IF (.not.is_empty(p)) input%advanced%sheets = p

      CALL parse(ok, p, lu, 'waves', ['=', ' '], 'advanced', com = '#')     !< waves
      CALL missing_arg(ok, .false., '')

      IF (ok .ne. 0) RETURN

      IF (.not.is_empty(p)) input%advanced%waves = p

      CALL parse(ok, p, lu, 'verbose', ['=', ' '], 'advanced', com = '#')     !< verbose
      CALL missing_arg(ok, .false., '')

      IF (ok .ne. 0) RETURN

      IF (.not.is_empty(p)) input%advanced%verbose = p

      CLOSE(lu, iostat = ok)

      IF (ok .ne. 0) THEN
        CALL report_error('Error while closing file ' + fo)
        RETURN
      ENDIF

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

    END SUBROUTINE parse_input_file

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

    SUBROUTINE echo()

      ! Purpose:
      !   to report input parameters to standard output.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   08/03/21                  original version
      !

      USE, INTRINSIC     :: iso_fortran_env, only: compiler_version

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
        IF (input%coda%model .eq. 'lw') THEN
          CALL update_log(num2char('Coherency parameters', width=30, fill='.') +     &
                          num2char('Model', width=15, justify='r')       + '|' +     &
                          num2char('Alpha', width=15, justify='r')       + '|' +     &
                          num2char('Threshold', width=15, justify='r')   + '|')
          CALL update_log(num2char('', width=30) +  &
                          num2char(input%coda%model, width=15, justify='r') + '|' + &
                          num2char(input%coda%alpha, notation='s', precision=3, width=15, justify='r') + '|', blankline=.false.)
        ELSEIF (input%coda%model .eq. 'hv') THEN
          CALL update_log(num2char('Coherency parameters', width=30, fill='.') +     &
                          num2char('Model', width=15, justify='r')       + '|' +     &
                          num2char('a', width=15, justify='r')       + '|' +     &
                          num2char('b', width=15, justify='r')   + '|' +  &
                          num2char('ak', width=15, justify='r') + '|' +  &
                          num2char('f0', width=15, justify='r') + '|')
          CALL update_log(num2char('', width=30) +  &
                          num2char(input%coda%model, width=15, justify='r') + '|' + &
                          num2char(input%coda%a, notation='f', precision=3, width=15, justify='r') + '|' + &
                          num2char(input%coda%b, notation='f', precision=3, width=15, justify='r') + '|' + &
                          num2char(input%coda%ak, notation='f', precision=3, width=15, justify='r') + '|' + &
                          num2char(input%coda%f0, notation='f', precision=3, width=15, justify='r') + '|', blankline=.false.)
        ENDIF
      ENDIF

      IF (input%source%is_point) THEN
        IF (input%origin%is_geo) THEN
          CALL update_log(num2char('Source parameters', width=30, fill='.') +     &
                          num2char('lon', width=15, justify='r') + '|' +  &
                          num2char('lat', width=15, justify='r') + '|' +  &
                          num2char('z', width=15, justify='r') + '|' +  &
                          num2char('Strike', width=15, justify='r') + '|' +  &
                          num2char('Dip', width=15, justify='r')    + '|')
          CALL update_log(num2char('', width=30) + &
                          num2char(input%source%lon, notation='f', width=15, precision=3, justify='r') + '|' + &
                          num2char(input%source%lat, notation='f', width=15, precision=3, justify='r') + '|' + &
                          num2char(input%source%z, notation='f', width=15, precision=3, justify='r') + '|' + &
                          num2char(input%source%strike, notation='f', width=15, precision=1, justify='r') + '|' + &
                          num2char(input%source%dip, notation='f', width=15, precision=1, justify='r')    + '|', blankline=.false.)
        ELSE
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
                          num2char(input%source%dip, notation='f', width=15, precision=1, justify='r')    + '|', blankline=.false.)
        ENDIF
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
                        num2char('Type', width=15, justify='r') + '|' + &
                        num2char('Frequency', width=15, justify='r') + '|' + &
                        num2char('Roughness', width=15, justify='r') + '|', blankline = .false.)
        IF (input%source%add_roughness) THEN
          CALL update_log(num2char('', width=30) +  &
                          num2char(input%source%type, width=15, justify='r') + '|' + &
                          num2char(input%source%freq, width=15, notation='f', precision=3, justify='r') + '|' + &
                          num2char(input%source%roughness, width=15, notation='f', precision=1, justify='r') + '|',  &
                          blankline = .false.)
        ELSE
          CALL update_log(num2char('', width=30) +  &
                          num2char(input%source%type, width=15, justify='r') + '|' + &
                          num2char(input%source%freq, width=15, notation='f', precision=3, justify='r') + '|' + &
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

      IF (input%origin%is_geo) THEN
        CALL update_log(num2char('Receivers List', width=30, fill='.') + &
                        num2char('lon', width=15, justify='r') + '|' +  &
                        num2char('lat', width=15, justify='r') + '|' +  &
                        num2char('z', width=15, justify='r') + '|' +  &
                        num2char('File', width=15, justify='r') + '|' + &
                        num2char('Vel Model', width=15, justify='r') + '|' + &
                        num2char('Att Model', width=15, justify='r') + '|')
      ELSE
        CALL update_log(num2char('Receivers List', width=30, fill='.') + &
                        num2char('x', width=15, justify='r') + '|' +  &
                        num2char('y', width=15, justify='r') + '|' +  &
                        num2char('z', width=15, justify='r') + '|' +  &
                        num2char('File', width=15, justify='r') + '|' + &
                        num2char('Vel Model', width=15, justify='r') + '|' + &
                        num2char('Att Model', width=15, justify='r') + '|')
      ENDIF

      IF (input%origin%is_geo) THEN
        DO i = 1, SIZE(input%receiver)
          CALL update_log(num2char('', width=30)  +   &
                          num2char(input%receiver(i)%lon, width=15, notation='f', precision=3, justify='r') + '|' + &
                          num2char(input%receiver(i)%lat, width=15, notation='f', precision=3, justify='r') + '|' + &
                          num2char(input%receiver(i)%z, width=15, notation='f', precision=3, justify='r') + '|' + &
                          num2char(TRIM(input%receiver(i)%file), width=15, justify='r') + '|' + &
                          num2char(input%receiver(i)%velocity, width=15, justify='r')   + '|' + &
                          num2char(input%receiver(i)%attenuation, width=15, justify='r') + '|', blankline = .false.)
        ENDDO
      ELSE
        DO i = 1, SIZE(input%receiver)
          CALL update_log(num2char('', width=30)  +   &
                          num2char(input%receiver(i)%x, width=15, notation='f', precision=3, justify='r') + '|' + &
                          num2char(input%receiver(i)%y, width=15, notation='f', precision=3, justify='r') + '|' + &
                          num2char(input%receiver(i)%z, width=15, notation='f', precision=3, justify='r') + '|' + &
                          num2char(TRIM(input%receiver(i)%file), width=15, justify='r') + '|' + &
                          num2char(input%receiver(i)%velocity, width=15, justify='r')   + '|' + &
                          num2char(input%receiver(i)%attenuation, width=15, justify='r') + '|', blankline = .false.)
        ENDDO
      ENDIF

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

    END SUBROUTINE echo

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

    SUBROUTINE broadcast()

      ! Purpose:
      !   to broadcast all input parameters to all processes. Make use of some manual packing/unpacking to minimize mpi calls.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   08/03/21                  original version
      !

#ifdef MPI

      INTEGER(i32)                            :: ierr, i, n
      CHARACTER(:), ALLOCATABLE               :: string
      INTEGER(i32), ALLOCATABLE, DIMENSION(:) :: intg
      LOGICAL,      ALLOCATABLE, DIMENSION(:) :: lg
      REAL(r32),    ALLOCATABLE, DIMENSION(:) :: float

      !-----------------------------------------------------------------------------------------------------------------------------

      float = [input%coda%fmax, input%coda%matching, input%coda%bandwidth, input%coda%alpha, input%source%x, &
               input%source%y, input%source%z, input%source%strike, input%source%dip, input%source%rake, input%source%m0,  &
               input%source%freq, input%source%roughness, input%source%correlation, input%source%l0, input%source%aparam, &
               input%source%vrfact, input%advanced%vrfact, input%source%lon, input%source%lat, input%origin%lon, input%origin%lat, &
               input%coda%a, input%coda%b, input%coda%ak, input%coda%f0]

      CALL mpi_bcast(float, SIZE(float), mpi_real, 0, mpi_comm_world, ierr)

      input%coda%fmax = float(1); input%coda%matching = float(2); input%coda%bandwidth = float(3); input%coda%alpha = float(4)
      input%source%x = float(5); input%source%y = float(6); input%source%z = float(7)
      input%source%strike = float(8); input%source%dip = float(9); input%source%rake = float(10); input%source%m0 = float(11)
      input%source%freq = float(12); input%source%roughness = float(13); input%source%correlation = float(14)
      input%source%l0 = float(15); input%source%aparam = float(16); input%source%vrfact = float(17)
      input%advanced%vrfact = float(18); input%source%lon = float(19); input%source%lat = float(20); input%origin%lon = float(21)
      input%origin%lat = float(22); input%coda%a = float(23); input%coda%b = float(24); input%coda%ak = float(25)
      input%coda%f0 = float(26)

      ALLOCATE(intg(12))

      IF (ALLOCATED(input%receiver)) THEN
        intg = [SIZE(input%velocity), SIZE(input%attenuation), SIZE(input%receiver), input%coda%seed, input%coda%samples,   &
                input%source%seed, input%source%samples, input%advanced%pmw, input%advanced%avecuts, input%advanced%sheets,  &
                input%advanced%waves, input%advanced%verbose]
      ENDIF

      CALL mpi_bcast(intg, 12, mpi_int, 0, mpi_comm_world, ierr)

      IF (.not.ALLOCATED(input%receiver)) THEN
        ALLOCATE(input%velocity(intg(1)), input%attenuation(intg(2)), input%receiver(intg(3)))
        input%coda%seed = intg(4); input%coda%samples = intg(5); input%source%seed = intg(6); input%source%samples = intg(7)
        input%advanced%pmw = intg(8); input%advanced%avecuts = intg(9); input%advanced%sheets = intg(10)
        input%advanced%waves = intg(11); input%advanced%verbose = intg(12)
      ENDIF

      DO i = 1, SIZE(input%velocity)
        IF (ALLOCATED(input%velocity(i)%vp)) n = SIZE(input%velocity(i)%vp)
        CALL mpi_bcast(n, 1, mpi_int, 0, mpi_comm_world, ierr)
        IF (.not.ALLOCATED(input%velocity(i)%vp)) THEN
          ALLOCATE(input%velocity(i)%vp(n), input%velocity(i)%vs(n), input%velocity(i)%rho(n), input%velocity(i)%depth(n))
          ALLOCATE(input%velocity(i)%vpgrad(n), input%velocity(i)%vsgrad(n))
        ENDIF
        float = [input%velocity(i)%vp(:), input%velocity(i)%vs(:), input%velocity(i)%rho(:), input%velocity(i)%depth(:),  &
                 input%velocity(i)%vpgrad(:), input%velocity(i)%vsgrad(:)]
        CALL mpi_bcast(float, SIZE(float), mpi_real, 0, mpi_comm_world, ierr)
        input%velocity(i)%vp = float(1:n)
        input%velocity(i)%vs = float(n+1:2*n)
        input%velocity(i)%rho = float(2*n+1:3*n)
        input%velocity(i)%depth = float(3*n+1:4*n)
        input%velocity(i)%vpgrad = float(4*n+1:5*n)
        input%velocity(i)%vsgrad = float(5*n+1:6*n)
      ENDDO

      DO i = 1, SIZE(input%attenuation)
        IF (ALLOCATED(input%attenuation(i)%gpp)) n = SIZE(input%attenuation(i)%gpp)
        CALL mpi_bcast(n, 1, mpi_int, 0, mpi_comm_world, ierr)
        IF (.not.ALLOCATED(input%attenuation(i)%gpp)) THEN
          ALLOCATE(input%attenuation(i)%gpp(n), input%attenuation(i)%gps(n), input%attenuation(i)%gss(n), input%attenuation(i)%b(n))
          ALLOCATE(input%attenuation(i)%lcut(n), input%attenuation(i)%hcut(n))
        ENDIF
        float = [input%attenuation(i)%gpp(:), input%attenuation(i)%gps(:), input%attenuation(i)%gss(:), input%attenuation(i)%b(:), &
                 input%attenuation(i)%lcut(:), input%attenuation(i)%hcut(:)]
        CALL mpi_bcast(float, SIZE(float), mpi_real, 0, mpi_comm_world, ierr)
        input%attenuation(i)%gpp  = float(1:n)
        input%attenuation(i)%gps  = float(n+1:2*n)
        input%attenuation(i)%gss  = float(2*n+1:3*n)
        input%attenuation(i)%b    = float(3*n+1:4*n)
        input%attenuation(i)%lcut = float(4*n+1:5*n)
        input%attenuation(i)%hcut = float(5*n+1:6*n)
      ENDDO

      n = SIZE(input%receiver)

      float = [input%receiver(:)%x, input%receiver(:)%y, input%receiver(:)%z, input%receiver(:)%lon, input%receiver(:)%lat]

      CALL mpi_bcast(float, SIZE(float), mpi_real, 0, mpi_comm_world, ierr)

      input%receiver(:)%x = float(1:n)
      input%receiver(:)%y = float(n+1:2*n)
      input%receiver(:)%z = float(2*n+1:3*n)
      input%receiver(:)%lon = float(3*n+1:4*n)
      input%receiver(:)%lat = float(4*n+1:5*n)

      intg = [input%receiver(:)%velocity, input%receiver(:)%attenuation]

      CALL mpi_bcast(intg, SIZE(intg), mpi_int, 0, mpi_comm_world, ierr)

      input%receiver(:)%velocity = intg(1:n)
      input%receiver(:)%attenuation = intg(n+1:2*n)

      CALL mpi_bcast(input%receiver(:)%file, n*32, mpi_character, 0, mpi_comm_world, ierr)

      string = input%source%file // input%source%type // input%input%folder // input%input%variable // input%input%format //  &
               input%coda%model // input%output%folder // input%output%variable

      CALL mpi_bcast(string, LEN(string), mpi_character, 0, mpi_comm_world, ierr)

      input%source%file = string(1:256); n = 257
      input%source%type = string(n:n+8-1); n = n + 8
      input%input%folder = string(n:n+256-1); n = n + 256
      input%input%variable = string(n:n+16-1); n = n + 16
      input%input%format = string(n:n+8-1); n = n + 8
      input%coda%model = string(n:n+4-1); n = n + 4
      input%output%folder = string(n:n+256-1); n = n + 256
      input%output%variable = string(n:n+16-1)

      lg = [input%source%is_point, input%source%add_roughness, input%source%add_rik, input%coda%add_coherency]

      CALL mpi_bcast(lg, SIZE(lg), mpi_logical, 0, mpi_comm_world, ierr)

      input%source%is_point = lg(1); input%source%add_roughness = lg(2); input%source%add_rik = lg(3)
      input%coda%add_coherency = lg(4)

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

    SUBROUTINE split_task(npts, ntasks, rank, i0, i1)

      ! Purpose:
      !   To evenly distribute a "npts" processes amongst "ntasks" tasks, returning first "i0" and last "i1" index for each "rank"
      !   process.
      !
      ! Revisions:
      !     Date                    Description of change
      !     ====                    =====================
      !   08/03/21                  original version
      !

      INTEGER(i32), INTENT(IN)  :: npts, rank, ntasks
      INTEGER(i32), INTENT(OUT) :: i0, i1

      !-----------------------------------------------------------------------------------------------------------------------------

      i0 = 1 + INT( REAL(npts, r32) / REAL(ntasks, r32) * REAL(rank, r32) )
      i1 = INT( REAL(npts, r32) / REAL(ntasks, r32) * REAL(rank + 1, r32) )

    END SUBROUTINE split_task

    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *
    !===============================================================================================================================
    ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- *

END MODULE m_toolbox
