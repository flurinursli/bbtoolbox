MODULE m_toolbox

  USE, NON_INTRINSIC :: m_precisions
  USE, NON_INTRINSIC :: m_logfile
  USE, NON_INTRINSIC :: m_parser
  USE, NON_INTRINSIC :: m_strings

  IMPLICIT none

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  INTERFACE wrong_par
    MODULE PROCEDURE wrong_par_i32, wrong_par_r32
  END INTERFACE wrong_par

  ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --

  TYPE :: src
    CHARACTER(16)  :: mrf
    CHARACTER(256) :: file
    INTEGER(i32)   :: seed, samples
    LOGICAL        :: is_point, add_roughness
    REAL(r32)      :: x, y, z, m0, strike, dip, rake, freq
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
    CHARACTER(256) :: file
  END TYPE io

  TYPE :: hf
    CHARACTER(4) :: model
    INTEGER(i32) :: seed, samples
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

      INTEGER(i32),             INTENT(OUT) :: ok
      CHARACTER(64)                         :: fo
      CHARACTER(:), ALLOCATABLE             :: str
      INTEGER(i32)                          :: lu, n, i, p

      !---------------------------------------------------------------------------------------------------------------------------------

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

      n = recurrences(ok, lu, 'rec', com = '#')                        !< find number of receivers

      IF (ok .ne. 0) THEN
        CALL report_error(parser_error(ok))
        RETURN
      ENDIF

      IF (is_empty(n) .or. (n .le. 0)) THEN
        CALL report_error('No stations found in input file')
        ok = 1
        RETURN
      ENDIF

      ALLOCATE(input%receiver(n))

      ! loop over stations to store parameters
      DO i = 1, n

        CALL parse(ok, str, lu, 'file', ["'", "'"], 'rec', nkey = i, com = '#')

        IF (ok .ne. 0) THEN
          CALL report_error(parser_error(ok))
          RETURN
        ENDIF

        IF (is_empty(str)) THEN
          CALL report_error('Argument "file" for receiver #' + num2char(i) + ' not found')
          ok = 1
          RETURN
        ENDIF

        input%receiver(i)%file = TRIM(str)                                         !< file

        CALL parse(ok, input%receiver(i)%velocity, lu, 'velocity', ['=', ','], 'rec', nkey = i, com = '#')   !< velocity

        IF (ok .ne. 0) THEN
          CALL report_error(parser_error(ok))
          RETURN
        ENDIF

        IF (is_empty(input%receiver(i)%velocity) .eqv. .false.) input%receiver(i)%velocity = 1

        CALL parse(ok, input%receiver(i)%attenuation, lu, 'attenuation', ['=', ','], 'rec', nkey = i, com = '#')  !< attenuation

        IF (ok .ne. 0) THEN
          CALL report_error(parser_error(ok))
          RETURN
        ENDIF

        IF (is_empty(input%receiver(i)%attenuation)) input%receiver(i)%attenuation = 1

        CALL parse(ok, input%receiver(i)%x, lu, 'x', ['=', ','], 'rec', nkey = i, com = '#')       !< x

        IF (ok .ne. 0) THEN
          CALL report_error(parser_error(ok))
          RETURN
        ENDIF

        IF (is_empty(input%receiver(i)%x)) THEN
          CALL report_error('Argument "x" for receiver #' + num2char(i) + ' not found')
          ok = 1
          RETURN
        ENDIF

        CALL parse(ok, input%receiver(i)%y, lu, 'y', ['=', ','], 'rec', nkey = i, com = '#')       !< y

        IF (ok .ne. 0) THEN
          CALL report_error(parser_error(ok))
          RETURN
        ENDIF

        IF (is_empty(input%receiver(i)%y)) THEN
          CALL report_error('Argument "y" for receiver #' + num2char(i) + ' not found')
          ok = 1
          RETURN
        ENDIF

        CALL parse(ok, input%receiver(i)%z, lu, 'z', ['=', ','], 'rec', nkey = i, com = '#')    !< z

        IF (ok .ne. 0) THEN
          CALL report_error(parser_error(ok))
          RETURN
        ENDIF

        IF (is_empty(input%receiver(i)%z)) THEN
          CALL report_error('Argument "z" for receiver #' + num2char(i) + ' not found')
          ok = 1
          RETURN
        ENDIF

      ENDDO

      CALL parse(ok, input%coda%fmax, lu, 'fmax', ['=', ','], 'coda', com = '#')               !< fmax
      IF (ok .ne. 0) THEN
        CALL report_error(parser_error(ok))
        RETURN
      ENDIF

      IF (is_empty(input%coda%fmax)) THEN
        CALL report_error('Argument "fmax" for keyword "coda" not found')
        ok = 1
        RETURN
      ENDIF

      CALL parse(ok, input%coda%matching, lu, 'matching frequency', ['=', ','], 'coda', com = '#')    !< matching
      IF (ok .ne. 0) THEN
        CALL report_error(parser_error(ok))
        RETURN
      ENDIF

      IF (is_empty(input%coda%matching)) THEN
        CALL report_error('Argument "matching frequency" for keyword "coda" not found')
        ok = 1
        RETURN
      ENDIF

      CALL parse(ok, input%coda%bandwidth, lu, 'bandwidth', ['=', ','], 'coda', com = '#')     !< bandwidth
      IF (ok .ne. 0) THEN
        CALL report_error(parser_error(ok))
        RETURN
      ENDIF

      IF (is_empty(input%coda%bandwidth)) THEN
        CALL report_error('Argument "bandwidth" for keyword "coda" not found')
        ok = 1
        RETURN
      ENDIF

      CALL parse(ok, input%coda%seed, lu, 'seed', ['=', ','], 'coda', com = '#')     !< seed
      IF (ok .ne. 0) THEN
        CALL report_error(parser_error(ok))
        RETURN
      ENDIF

      IF (is_empty(input%coda%seed)) THEN
        CALL report_error('Argument "seed" for keyword "coda" not found')
        ok = 1
        RETURN
      ENDIF

      CALL parse(ok, input%coda%samples, lu, 'samples', ['=', ','], 'coda', com = '#')    !< samples
      IF (ok .ne. 0) THEN
        CALL report_error(parser_error(ok))
        RETURN
      ENDIF

      IF (is_empty(input%coda%samples)) THEN
        CALL report_error('Argument "samples" for keyword "coda" not found')
        ok = 1
        RETURN
      ENDIF

      CALL parse(ok, str, lu, 'model', ["'", "'"], 'coda', com = '#')
      IF (ok .ne. 0) THEN
        CALL report_error(parser_error(ok))
        RETURN
      ENDIF

      input%coda%model = TRIM(str)       !< model

      CALL parse(ok, input%coda%alpha, lu, 'alpha', ['=', ','], 'coda', com = '#')     !< alpha
      IF (ok .ne. 0) THEN
        CALL report_error(parser_error(ok))
        RETURN
      ENDIF

      CALL parse(ok, input%coda%threshold, lu, 'threshold', ['=', ','], 'coda', com = '#')     !< threshold
      IF (ok .ne. 0) THEN
        CALL report_error(parser_error(ok))
        RETURN
      ENDIF

      CALL parse(ok, str, lu, 'folder', ["'", "'"], 'input', com = '#')
      IF (ok .ne. 0) THEN
        CALL report_error(parser_error(ok))
        RETURN
      ENDIF

      IF (is_empty(str)) THEN
        CALL report_error('Argument "folder" for keyword "input" not found')
        ok = 1
        RETURN
      ENDIF

      input%input%folder = TRIM(str)     !< folder

      CALL parse(ok, str, lu, 'format', ["'", "'"], 'input', com = '#')
      IF (ok .ne. 0) THEN
        CALL report_error(parser_error(ok))
        RETURN
      ENDIF

      IF (is_empty(str)) THEN
        CALL report_error('Argument "format" for keyword "input" not found')
        ok = 1
        RETURN
      ENDIF

      input%input%format = TRIM(str)     !< format

      CALL parse(ok, str, lu, 'variable', ["'", "'"], 'input', com = '#')
      IF (ok .ne. 0) THEN
        CALL report_error(parser_error(ok))
        RETURN
      ENDIF

      IF (is_empty(str)) THEN
        CALL report_error('Argument "variable" for keyword "input" not found')
        ok = 1
        RETURN
      ENDIF

      input%input%variable = TRIM(str)     !< variable

      CALL parse(ok, str, lu, 'folder', ["'", "'"], 'output', com = '#')
      IF (ok .ne. 0) THEN
        CALL report_error(parser_error(ok))
        RETURN
      ENDIF

      IF (is_empty(str)) THEN
        CALL report_error('Argument "folder" for keyword "output" not found')
        ok = 1
        RETURN
      ENDIF

      input%output%folder = TRIM(str)     !< folder

      CALL parse(ok, str, lu, 'variable', ["'", "'"], 'output', com = '#')
      IF (ok .ne. 0) THEN
        CALL report_error(parser_error(ok))
        RETURN
      ENDIF

      IF (is_empty(str)) THEN
        CALL report_error('Argument "variable" for keyword "output" not found')
        ok = 1
        RETURN
      ENDIF

      input%output%variable = TRIM(str)     !< variable


! down here!!!!!






      ! scan for frequency bands
      CALL parse(ok, fbands, lu, 'Bands', ['{', '}'], ',', ';', com = '#')                       !< fbands`
      IF (ok .ne. 0) CALL report_error(parser_error(ok))

      IF (ALL(is_empty(fbands))) THEN
        CALL report_error('Argument "Bands" not found')
        ok = 1
        RETURN
      ENDIF

      ! scan for inversion parameters
      CALL parse(ok, etass, lu, 'EtaSS', ['{', '}'], ',', com = '#')                             !< etass
      IF (ok .ne. 0) CALL report_error(parser_error(ok))

      IF (ALL(is_empty(etass))) THEN
        CALL report_error('Argument "EtaS" not found')
        ok = 1
        RETURN
      ENDIF

      CALL parse(ok, nu, lu, 'Nu', ['{', '}'], ',', com = '#')                                   !< nu
      IF (ok .ne. 0) CALL report_error(parser_error(ok))

      CALL parse(ok, mode, lu, 'Mode', ['=', ','], com = '#')                                    !< mode
      IF (ok .ne. 0) CALL report_error(parser_error(ok))

      IF (is_empty(mode)) THEN
        CALL report_error('Argument "Mode" not found')
        ok = 1
        RETURN
      ENDIF

      IF (mode .eq. 2) THEN

        CALL parse(ok, threshold, lu, 'Threshold', ['=', ','], com = '#')                       !< threshold
        IF (ok .ne. 0) CALL report_error(parser_error(ok))

        IF (is_empty(threshold)) THEN
          CALL report_error('Argument "Threshold" not found')
          ok = 1
          RETURN
        ENDIF

      ENDIF

      CALL parse(ok, sdwindow, lu, 'Swin', ['=', ','], 'DIRECT', com = '#')               !< sdwindow (for direct)
      IF (ok .ne. 0) CALL report_error(parser_error(ok))

      IF (is_empty(sdwindow)) THEN
        CALL report_error('Argument "Swin" for direct waves not found')
        ok = 1
        RETURN
      ENDIF

      CALL parse(ok, fwin, lu, 'Factor', ['=', '%'], 'DIRECT', com = '#')                !< fwin
      IF (ok .ne. 0) CALL report_error(parser_error(ok))

      ! by default, window starts at Tp and Ts
      IF (is_empty(fwin)) fwin = 100._r32

      fwin = fwin / 100._r32

      CALL parse(ok, scwindow, lu, 'Swin', ['=', ','], 'CODA', com = '#')                 !< scwindow (for coda)
      IF (ok .ne. 0) CALL report_error(parser_error(ok))

      IF (is_empty(scwindow)) THEN
        CALL report_error('Argument "Swin" for coda waves not found')
        ok = 1
        RETURN
      ENDIF

      CALL parse(ok, tlim, lu, 'Tlim', ['=', ','], 'CODA', com = '#')                 !< tlim (for coda)
      IF (ok .ne. 0) CALL report_error(parser_error(ok))

      ! by default, take whole time-series
      IF (is_empty(tlim)) tlim = 0._r32

      CALL parse(ok, nsi, lu, 'InitialModels', ['=', ','], com = '#')                           !< nsi
      IF (ok .ne. 0) CALL report_error(parser_error(ok))

      IF (is_empty(nsi)) THEN
        CALL report_error('Argument "InitialModels" not found')
        ok = 1
        RETURN
      ENDIF

      CALL parse(ok, ns, lu, 'Models', ['=', ','], com = '#')                                   !< ns
      IF (ok .ne. 0) CALL report_error(parser_error(ok))

      IF (is_empty(ns)) THEN
        CALL report_error('Argument "Models" not found')
        ok = 1
        RETURN
      ENDIF

      CALL parse(ok, nr, lu, 'Resampled', ['=', ','], com = '#')                                !< nr
      IF (ok .ne. 0) CALL report_error(parser_error(ok))

      IF (is_empty(nr)) THEN
        CALL report_error('Argument "Resampled" not found')
        ok = 1
        RETURN
      ENDIF

      CALL parse(ok, itermax, lu, 'Iterations', ['=', ','], com = '#')                          !< itermax
      IF (ok .ne. 0) CALL report_error(parser_error(ok))

      IF (is_empty(itermax)) THEN
        CALL report_error('Argument "Iterations" not found')
        ok = 1
        RETURN
      ENDIF

      CALL parse(ok, seed, lu, 'Seed', ['=', ','], com = '#')                                  !< seed
      IF (ok .ne. 0) CALL report_error(parser_error(ok))

      IF (is_empty(seed)) THEN
        CALL report_error('Argument "Seed" not found')
        ok = 1
        RETURN
      ENDIF

      CALL parse(ok, beta, lu, 'Beta', ['=', ','], com = '#')                                  !< beta
      IF (ok .ne. 0) CALL report_error(parser_error(ok))

      IF (is_empty(beta)) THEN
        CALL report_error('Argument "Beta" not found')
        ok = 1
        RETURN
      ENDIF

      CALL parse(ok, str, lu, 'Weight', ["'", "'"], com = '#')                                 !< noweight
      IF (ok .ne. 0) CALL report_error(parser_error(ok))

      IF (is_empty(str)) THEN
        CALL report_error('Argument "Weight" not found')
        ok = 1
        RETURN
      ENDIF

      noweight = .false.
      noweight = lowercase(str) .eq. 'n'

      ! if "nu" is not present, we are dealing with elastic RTT
      elastic = ALL(is_empty(nu))

      IF (elastic) THEN

        CALL parse(ok, etass2pp, lu, 'EtaSS/PP', ['{', '}'], ',', com = '#')                !< etass2pp
        IF (ok .ne. 0) CALL report_error(parser_error(ok))

        IF (ALL(is_empty(etass2pp))) THEN
          CALL report_error('Argument "EtaSS/PP" not found')
          ok = 1
          RETURN
        ENDIF

        CALL parse(ok, etaps2pp, lu, 'EtaPS/PP', ['{', '}'], ',', com = '#')                !< etaps2pp
        IF (ok .ne. 0) CALL report_error(parser_error(ok))

        IF (ALL(is_empty(etaps2pp))) THEN
          CALL report_error('Argument "EtaPS" not found')
          ok = 1
          RETURN
        ENDIF

        CALL parse(ok, pdwindow, lu, 'Pwin', ['=', ','], 'DIRECT', com = '#')               !< pdwindow (for direct)
        IF (ok .ne. 0) CALL report_error(parser_error(ok))

        IF (is_empty(pdwindow)) THEN
          CALL report_error('Argument "Pwin" for direct waves not found')
          ok = 1
          RETURN
        ENDIF

        CALL parse(ok, pcwindow, lu, 'Pwin', ['=', ','], 'CODA', com = '#')                 !< pcwindow (for coda)
        IF (ok .ne. 0) CALL report_error(parser_error(ok))

        IF (is_empty(pcwindow)) THEN
          CALL report_error('Argument "Pwin" for coda waves not found')
          ok = 1
          RETURN
        ENDIF

      ELSE

        ! assign some values to "acf" and "hurst" for acoustic isotropic RTT. Note: in such case, these are irrelevant
        IF (ALL(nu .eq. 0._r32)) THEN

          acf   = 'vk'
          hurst = [0.5_r32, 0.5_r32]

        ELSE

          CALL parse(ok, acf, lu, 'acf', ["'", "'"], com = '#')                                 !< acf
          IF (ok .ne. 0) CALL report_error(parser_error(ok))

          IF (is_empty(acf)) THEN
            CALL report_error('Argument "acf" not found')
            ok = 1
            RETURN
          ENDIF

          IF (acf .eq. 'vk') THEN
            CALL parse(ok, hurst, lu, 'Hurst', ['{', '}'], ',', com = '#')                     !< hurst
            IF (ok .ne. 0) CALL report_error(parser_error(ok))

            IF (ALL(is_empty(hurst))) THEN
              CALL report_error('Argument "hurst" not found')
              ok = 1
              RETURN
            ENDIF
          ELSE
            ! for Gaussian media, set Hurst exponent to any value to avoid triggering errors
            hurst = [0.5_r32, 0.5_r32]
          ENDIF

        ENDIF

      ENDIF

      ! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --
      ! -------------------------------------- check whether all input parameters make sense -------------------------------------------

      DO i = 1, SIZE(fbands, 2)
        IF (wrong_par(fbands(:, i), 0._r32, strict = .true.)) THEN
          CALL report_error('All frequency bands must be positive')
          ok = 1
          RETURN
        ENDIF
        IF (wrong_par(fbands(:, i), strict = .true.)) THEN
          CALL report_error('Lower bound for frequency bands must be smaller than upper bound')
          ok = 1
          RETURN
        ENDIF
      ENDDO

      IF (wrong_par([beta], 0._r32, strict = .true.)) THEN
        CALL report_error('Shear wave velocity ("beta") must be positive')
        ok = 1
        RETURN
      ENDIF

      IF (wrong_par([nsi, ns, nr, itermax], 0, strict = .true.)) THEN
        CALL report_error('NA parameters must be positive')
        ok = 1
        RETURN
      ENDIF

      IF (wrong_par([mode], 0, 2)) THEN
        CALL report_error('Mode parameter must lie in the range [0, 2]')
        ok = 1
        RETURN
      ENDIF

      IF (mode .eq. 2) THEN
        IF (wrong_par([threshold], 1)) THEN
          CALL report_error('Parameter "Threshold" must be larger than or equal to 1')
          ok = 1
          RETURN
        ENDIF
      ENDIF

      IF (wrong_par([sdwindow], 0._r32)) THEN
        CALL report_error('Parameter "Swin" for direct wave must be positive')
        ok = 1
        RETURN
      ENDIF
      IF (wrong_par([scwindow], 0._r32)) THEN
        CALL report_error('Parameter "Swin" for coda wave must be positive')
        ok = 1
        RETURN
      ENDIF
      IF (wrong_par([tlim], 0._r32)) THEN
        CALL report_error('Parameter "Tlim" for coda wave must be positive')
        ok = 1
        RETURN
      ENDIF

      IF (wrong_par([fwin], 0._r32)) THEN
        CALL report_error('Parameter "Factor" for direct wave must be positive')
        ok = 1
        RETURN
      ENDIF

      IF (wrong_par(etass, 0._r32, strict = .true.)) THEN
        CALL report_error('Parameter "EtaSS" must be positive')
        ok = 1
        RETURN
      ENDIF
      IF (wrong_par(etass)) THEN
        CALL report_error('Lower bound for parameter "EtaSS" must be smaller than or equal to upper bound')
        ok = 1
        RETURN
      ENDIF

      IF (elastic) THEN
        IF (wrong_par(etaps2pp, 0._r32, strict = .true.)) THEN
          CALL report_error('Parameter "EtaPS/PP" must be positive')
          ok = 1
          RETURN
        ENDIF
        IF (wrong_par(etaps2pp)) THEN
          CALL report_error('Lower bound for parameter "EtaPS/PP" must be smaller than or equal to upper bound')
          ok = 1
          RETURN
        ENDIF

        IF (wrong_par(etass2pp, 0._r32, strict = .true.)) THEN
          CALL report_error('Parameter "EtaSS/PP" must be positive')
          ok = 1
          RETURN
        ENDIF
        IF (wrong_par(etass2pp)) THEN
          CALL report_error('Lower bound for parameter "EtaSS/PP" must be smaller than or equal to upper bound')
          ok = 1
          RETURN
        ENDIF

        IF (wrong_par([pdwindow], 0._r32)) THEN
          CALL report_error('Parameter "Pwin" for direct wave must be positive')
          ok = 1
          RETURN
        ENDIF
        IF (wrong_par([pcwindow], 0._r32)) THEN
          CALL report_error('Parameter "Pwin" for coda wave must be positive')
          ok = 1
          RETURN
        ENDIF

      ELSE
        IF (wrong_par(nu, 0._r32)) THEN
          CALL report_error('"Nu" must be positive')
          ok = 1
          RETURN
        ENDIF
        IF (wrong_par(nu)) THEN
          CALL report_error('Lower bound for "Nu" must be smaller than or equal to upper bound')
          ok = 1
          RETURN
        ENDIF
        IF (wrong_par(hurst, 0._r32, 1._r32, .true.)) THEN
          CALL report_error('Hurst exponent must be in the interval (0, 1)')
          ok = 1
          RETURN
        ENDIF
        IF ( (acf .ne. 'vk') .and. (acf .ne. 'gs') ) THEN
          CALL report_error('Unknown value for "acf": ' + acf)
          ok = 1
          RETURN
        ENDIF
      ENDIF

    END SUBROUTINE read_input_file

! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
!===================================================================================================================================
! --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * --- * ---
