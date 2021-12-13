
#define DIRNIST '../NIST/'

! nist_isotope.txt: fetched from
! Atomic Weights and Isotopic Compositions with Relative Atomic Masses
! using "get_nist_atomic.sh"
!
! *.states and *.ionize: fetched from
! NIST Atomic Spectra Database Levels Form
! using "get_nist_level.sh"

MODULE nist_module
  
  USE ISO_FORTRAN_ENV
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: get_nist_mass_comp, get_nist_eneion, get_nist_states, get_nist_isofrac
  
CONTAINS

  SUBROUTINE get_nist_mass_comp(an_max, amass, comp_iso)
    INTEGER(INT32), INTENT(IN) :: an_max ! atomic number maximum
    REAL(REAL64), INTENT(INOUT) :: amass(:)
    REAL(REAL64), INTENT(INOUT) :: comp_iso(:,:)

    REAL(REAL64) :: mass, comp
    INTEGER(INT32) :: mn, an, is, iostat
    CHARACTER :: dummy*256, mass_c*20, comp_c*20, smass_c*20
    CHARACTER*20, PARAMETER :: a20 = '                    '

    comp_iso = 0d0
    
    OPEN(10, FILE=DIRNIST//'nist_isotope.txt', STATUS='OLD')
    DO
       READ(10,'(A16,I4)',IOSTAT=iostat) dummy, an
       IF(an > an_max) EXIT
       READ(10,'(A16,A2)') 
       READ(10,'(A14,I3)') dummy, mn
       READ(10,'(A23,A20)') dummy, mass_c
       IF(mass_c /= a20) THEN
          is = INDEX(mass_c,'(')
          IF(is == 0) THEN
             READ(mass_c,*) mass
          ELSE
             READ(mass_c(1:is-1),*) mass
          END IF
       END IF
       READ(10,'(A23,A20)') dummy, comp_c
       IF(comp_c /= a20) THEN
          is = INDEX(comp_c,'(')
          IF(is == 0) THEN
             READ(comp_c,*) comp
          ELSE
             READ(comp_c(1:is-1),*) comp
          END IF
       ELSE
          comp = 0d0
       END IF
       comp_iso(an,mn) = comp
       READ(10,'(A25,A20)') dummy, smass_c
       IF(smass_c /= a20) THEN
          IF(INDEX(smass_c,'[') == 0) THEN
             READ(smass_c(1:INDEX(smass_c,'(')-1),*) amass(an)
          ELSE
             IF(INDEX(smass_c,',') == 0) THEN
                READ(smass_c(INDEX(smass_c,'[')+1:INDEX(smass_c,']')-1),*) amass(an)
             ELSE
                amass(an) = amass(an) + mass * comp
             END IF
          END IF
       END IF
       READ(10,*)
       READ(10,*)
    END DO
    CLOSE(10)
    RETURN
  END SUBROUTINE get_nist_mass_comp
  

  SUBROUTINE get_nist_eneion(natm, nion, eneion)
    CHARACTER*3, INTENT(IN) :: natm
    CHARACTER*2, INTENT(IN) :: nion
    REAL(REAL64), INTENT(OUT) :: eneion

    INTEGER(INT32) :: iostat
      
    OPEN(1, FILE=DIRNIST//'levels/nist_'//natm//'.'//nion//'.ionize', STATUS='OLD')
    READ(1, *, IOSTAT=iostat) eneion
    IF(iostat /= 0) THEN
       PRINT *, '*** ERROR: could not get ionization energy for '//natm//'.'//nion
       STOP
    END IF
    CLOSE(1)
    
    RETURN
  END SUBROUTINE get_nist_eneion

  
  SUBROUTINE get_nist_states(natm, nion, gtot, ene, term, conf)
    CHARACTER*3, INTENT(IN) :: natm
    CHARACTER*2, INTENT(IN) :: nion
    INTEGER(INT32), ALLOCATABLE, INTENT(OUT) :: gtot(:)
    REAL(REAL64), ALLOCATABLE, INTENT(OUT) :: ene(:)
    CHARACTER*8, ALLOCATABLE, INTENT(INOUT) :: term(:)
    CHARACTER*32, ALLOCATABLE, INTENT(INOUT) :: conf(:)
    
    INTEGER(INT32) :: nstates, l, iostat
    INTEGER(INT32), PARAMETER :: max_buf = 256
    CHARACTER(max_buf):: linebuf
    
    OPEN(1, FILE=DIRNIST//'levels/nist_'//natm//'.'//nion//'.states', STATUS='OLD')
    nstates = 0
    DO
       READ(1,*,IOSTAT=iostat)
       IF(iostat /= 0) EXIT
       nstates = nstates + 1
    END DO
    REWIND(1)
    ALLOCATE(gtot(nstates), ene(nstates), term(nstates), conf(nstates))
    DO l = 1, nstates
       READ(1,'(A)',IOSTAT=iostat) linebuf
       READ(linebuf,*,IOSTAT=iostat) gtot(l), ene(l), term(l), conf(l)
       IF(iostat /= 0) THEN
          BACKSPACE(1)
          READ(1,*,IOSTAT=iostat) gtot(l), ene(l), term(l)
          IF(iostat /= 0) THEN
             term(l) = 'N/A'
             conf(l) = 'N/A'
          ELSE
             conf(l) = 'N/A'
          END IF
       END IF
    END DO
    CLOSE(1)
    
    RETURN
  END SUBROUTINE get_nist_states
  
  
  SUBROUTINE get_nist_isofrac(frac)
    REAL(REAL64), INTENT(OUT) :: frac
    
    REAL(REAL64) :: comp
    INTEGER(INT32) mass1, number, mass0, is, iostat
    CHARACTER*2 :: sym, sym0
    CHARACTER :: dummy*256, mass_c*20, comp_c*20
    CHARACTER*20, PARAMETER :: a20 = '                    '

    print '(a)', 'read composition:'
    
    OPEN(50, FILE='./composition.txt', STATUS='OLD', IOSTAT=iostat)
    IF(iostat /= 0) THEN
       print '(a)', '*** error: no composition file: ./composition.txt'
       STOP
    END IF
    OPEN(51, FILE=DIRNIST//'nist_isotope.txt', STATUS='OLD')
    frac = 1d0
    DO
       READ(50,*,IOSTAT=iostat) mass0, sym0, number
       IF(iostat /= 0) EXIT
       DO
          READ(51,*,IOSTAT=iostat)
          IF(iostat /= 0) THEN
             print '(a)', '*** error *** could not find symbol: ', sym0
             STOP
          END IF
          READ(51,'(A16,A2)') dummy, sym
          READ(51,'(A14,I3)') dummy, mass1
          READ(51,'(A23,A20)') dummy, mass_c
          READ(51,'(A23,A20)') dummy, comp_c
          IF(comp_c /= a20) THEN
             is = INDEX(comp_c,'(')
             IF(is == 0) THEN
                READ(comp_c,*) comp
             ELSE
                READ(comp_c(1:is-1),*) comp
             END IF
          END IF
          READ(51,*)
          READ(51,*)
          READ(51,*)
          IF(TRIM(sym) == TRIM(sym0) .AND. mass1 == mass0) EXIT
       END DO
       print '(a4,f9.7,i2)', ' '//sym//' ', comp, number
       frac = frac * comp ** number
       REWIND(51)
    END DO
    CLOSE(51)
    CLOSE(50)

    print '(a,f9.7)', ' isotopic fraction = ', frac
    print *, ''
    
    RETURN
  END SUBROUTINE get_nist_isofrac

  
END MODULE nist_module
