
! MAKE CONFIGURATION FILE, T-P GRID AND TABLE OF (T, LNK(T)) FOR FASTCHEM

#define ATOMIC_IONS
#define FASTCHEM_ORIGINAL
#define OPTAB_DATABASE_DIR "/S/data00/G5106/y0582/database/"
#define FASTCHEM_INPUT_DIR "/S/home01/G5106/y0582/FastChem/input/"

PROGRAM prep_FastChem

  USE ISO_FORTRAN_ENV
  USE const_module, ONLY : k_bol, clight, m_ele, pi, h, ev_to_erg
  IMPLICIT NONE

  REAL(REAL64), PARAMETER :: bar = 1d6
  REAL(REAL64) :: p, pmin, pmax, dp
  REAL(REAL64) :: a1, a2, a3, a4, a5
  REAL(REAL128) :: tmp_min, tmp_max, dtmp
  REAL(REAL128), ALLOCATABLE, DIMENSION(:) :: tmp, q_i, q_n, lnk_save
  REAL(REAL128) :: lnk, kt
  REAL(REAL64) :: eneion, aff(92)
  REAL(REAL64), PARAMETER :: tmp_threshold = 1d4
  INTEGER :: lmax, l, na, ni, ab, lnk_max, iostat, j, jmax
  INTEGER, PARAMETER :: max_line_len = 1024
  CHARACTER(MAX_LINE_LEN) :: linebuf
  CHARACTER :: species*16
  CHARACTER*64 :: fname_config, fname_tpgrid, fname_output, fname_output_monitor, fname_species, dname_lnK, id
  CHARACTER :: cna*3, cni*2, fname*32, cnia*3, cnib*3, csign*1, fname_pf*32
  CHARACTER(LEN=2), PARAMETER, DIMENSION(92) :: &
       celem = ['H ', 'He', 'Li', 'Be', 'B ', 'C ', 'N ', 'O ', 'F ', 'Ne', &
       'Na', 'Mg', 'Al', 'Si', 'P ', 'S ', 'Cl', 'Ar', 'K ', 'Ca', &
       'Sc', 'Ti', 'V ', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', &
       'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y ', 'Zr', &
       'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', &
       'Sb', 'Te', 'I ', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', &
       'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', &
       'Lu', 'Hf', 'Ta', 'W ', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', &
       'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', &
       'Pa', 'U ']
  CHARACTER(LEN=2), PARAMETER, DIMENSION(79) :: &
       abund = ['H ', 'He', 'O ', 'C ', 'Ne', 'N ', 'Mg', 'Si', 'Fe', 'S ', &
       'Al', 'Ar', 'Ca', 'Na', 'Ni', 'Cr', 'Cl', 'Mn', 'P ', 'K ', &
       'Co', 'Ti', 'Zn', 'F ', 'Cu', 'V ', 'Ge', 'Se', 'Sc', 'Ga', &
       'Sr', 'B ', 'Zr', 'Rb', 'As', 'Y ', 'Te', 'Ba', 'Sn', 'Mo', &
       'Ru', 'Pb', 'Cd', 'Pt', 'Ce', 'Pd', 'Nb', 'Nd', 'Os', 'Ir', &
       'Be', 'Hg', 'La', 'Dy', 'Cs', 'Gd', 'Li', 'Sb', 'Sm', 'Ag', &
       'Er', 'Au', 'Rh', 'Tl', 'W ', 'Hf', 'Yb', 'In', 'Pr', 'Bi', &
       'Eu', 'Ho', 'Tb', 'Re', 'Tm', 'Lu', 'Th', 'Ta', 'U ']


  OPEN(1,FILE='prep_FastChem.dat',STATUS='OLD', RECL=max_line_len)
  READ(1,*) id
  READ(1,*) tmp_min, tmp_max, lmax
  READ(1,*) pmin, pmax, jmax
  CLOSE(1)
!  id = 'tab40'

  ! T-P GRID
  ! log10(P/Ba)
!  jmax = 40 + 1
!  pmin = -4d0
!  pmax = 9d0

  ! log10(T/K)
!  lmax = 40 + 1
!  tmp_min = 2.5d0
!  tmp_max = 6d0
  
  CALL SYSTEM('mkdir '//FASTCHEM_INPUT_DIR//'lnK_'//id)
  CALL SYSTEM('mkdir '//FASTCHEM_INPUT_DIR//'../output')

  fname_config = 'config.input_'//TRIM(id)
  fname_tpgrid = 'tpgrid_'//TRIM(id)//'.dat'
  fname_output = TRIM(id)//'.dat'
  fname_output_monitor = TRIM(id)//'_monitor.dat'
  fname_species = 'lnK_'//TRIM(id)//'.dat'
  dname_lnk = 'lnK_'//TRIM(id)

  PRINT *, 'configuration file:', TRIM(fname_config)
  ! CONFIGURATION FILE
  OPEN(1, FILE=FASTCHEM_INPUT_DIR//TRIM(fname_config), RECL=max_line_len)
  WRITE(1,*) '#Atmospheric profile input file'
  WRITE(1,*) 'input/'//TRIM(fname_tpgrid)
  WRITE(1,*) ''
  WRITE(1,*) '#Chemistry output file'
  WRITE(1,*) 'output/'//TRIM(fname_output)
  WRITE(1,*) ''
  WRITE(1,*) '#Monitor output file'
  WRITE(1,*) 'output/'//TRIM(fname_output_monitor)
  WRITE(1,*) ''
  WRITE(1,*) '#FastChem console verbose level (1 - 4); 1 = almost silent, 4 = detailed console output'
  WRITE(1,*) '4'
  WRITE(1,*) ''
  WRITE(1,*) '#Output mixing ratios (MR) or particle number densities (ND, default)'
  WRITE(1,*) 'ND'
  WRITE(1,*) ''
  WRITE(1,*) '#Element abundance file'
  WRITE(1,*) 'input/element_abundances_solar.dat'
  WRITE(1,*) ''
  WRITE(1,*) '#Species data file'
  WRITE(1,*) 'input/'//TRIM(fname_species)
  WRITE(1,*) ''
  WRITE(1,*) '#Accuracy of chemistry iteration'
  WRITE(1,*) '1.0e-4'
  WRITE(1,*) ''
  WRITE(1,*) '#Max number of chemistry iterations'
  WRITE(1,*) '80000'
  WRITE(1,*) ''
  WRITE(1,*) '#Max number internal solver iterations'
  WRITE(1,*) '20000'
  WRITE(1,*) ''
  CLOSE(1)

  PRINT *, 'profile input file:', TRIM(fname_tpgrid)
  ALLOCATE(tmp(lmax), q_i(lmax), q_n(lmax), lnk_save(lmax))

  OPEN(1,FILE=FASTCHEM_INPUT_DIR//TRIM(fname_tpgrid), STATUS='REPLACE', RECL=max_line_len)
  WRITE(1,'(A)') '# T[K]-P[bar] grid'
  WRITE(1,'(A,F6.2,A,F6.2,A,I3,A)') '# logT = [',REAL(tmp_min),',',REAL(tmp_max),',',lmax,']'
  WRITE(1,'(A,F6.2,A,F6.2,A,I3,A)') '# logP = [',REAL(pmin),',',REAL(pmax),',',jmax,']'

  dtmp = (tmp_max - tmp_min) / (lmax - 1)
  dp = (pmax - pmin) / (jmax - 1)
  
  DO l = 1, lmax
     tmp(l) = tmp_min + (l - 1) * dtmp
     tmp(l) = 10d0**tmp(l)
     DO j = 1, jmax
        p = pmin + (j - 1) * dp
        !FASTCHEM ACCEPTS PRESSURE IN BAR
        WRITE(1,*) 10d0**p / bar, tmp(l)
     END DO
  END DO
  CLOSE(1) 
  

  ! read electron_affinity data
  OPEN(1,FILE='electron_affinity.tsv',STATUS='OLD', RECL=max_line_len)
  READ(1,*)
  READ(1,*)
  DO na = 1, 92
     READ(1,*) aff(na)
  END DO
  CLOSE(1)
  aff = aff * ev_to_erg / (h * clight) ! eV -> cm^-1
  

  PRINT *, 'species data file:', TRIM(fname_species)
  PRINT *, 'lnK directory:', TRIM(dname_lnK)

  OPEN(2,FILE=FASTCHEM_INPUT_DIR//TRIM(fname_species), RECL=max_line_len)
  WRITE(2,*) '#'
  WRITE(2,*) '#'
  WRITE(2,*) '#'

  
#ifdef FASTCHEM_ORIGINAL
  PRINT *, 'molecule lnK'
  OPEN(1, FILE=FASTCHEM_INPUT_DIR//'logK_ext.dat', STATUS='OLD', RECL=max_line_len)
  READ(1,*)!HEADER
  READ(1,*)!HEADER
  READ(1,*)!HEADER
  DO
     READ(1,'(A)',IOSTAT=iostat) linebuf
     IF(iostat /= 0) EXIT
     species = linebuf(1:INDEX(linebuf,' ')-1)
     READ(1,*) a1, a2, a3, a4, a5
!     PRINT *, TRIM(species)
     READ(1,*)

     fname = TRIM(dname_lnK)//'/'//TRIM(species)//'.dat'
#ifdef ATOMIC_IONS
     IF(exclude(species)) CYCLE
#endif     
     WRITE(2,*) TRIM(linebuf)
     WRITE(2,*) 'f input/'//TRIM(fname)
     WRITE(2,*) ' '
     
     OPEN(3, FILE=FASTCHEM_INPUT_DIR//TRIM(fname), RECL=max_line_len)
     WRITE(3,*) DBLE(tmp_min), DBLE(dtmp), ' log'
     DO l = 1, lmax
        lnK = a1/tmp(l) + a2*LOG(tmp(l)) + a3 + a4*tmp(l) + a5*tmp(l)**2
        IF(tmp(l) >= tmp_threshold) THEN
!           WRITE(3,*) DBLE(LOG(TINY(lnK)))
           WRITE(3,*) LOG(TINY(DBLE(lnK)))
        ELSE
           WRITE(3,*) DBLE(lnK)
        END IF
     END DO
     CLOSE(3)

!     PRINT *, 'lnK/'//TRIM(species)//'.dat'
  END DO
  CLOSE(1)
#endif  
#ifdef ATOMIC_IONS
  PRINT *, 'atom lnK'
  DO ab = 1, 79 ! abundance order
     DO na = 1, 92
        IF(celem(na) == abund(ab)) EXIT
     END DO
     WRITE(*,'(A3)',ADVANCE='NO') celem(na)//' '
     lnk_max = 0d0
     DO ni = -1, na-1 ! ionization level (-1: anion, 0: neutral, 1: 1+, 2:, 2+, ....)

        IF(ni+1 .GT. na) CYCLE ! veto
        IF(ni .EQ. -1 .AND. na .NE. 1) CYCLE ! consider anion only for H
        
        IF(ni .LE. 0) lnk_save(:) = 0d0

        WRITE(cna,'(I3.3)') na
        WRITE(cni,'(I2.2)') ni+1
        IF(ni   .EQ. -1) THEN ! ANION
           CALL atomic_pf(q_n, na+1, 0, tmp)
           eneion = aff(na)
           csign = '-'
           cnib = '1'
           cnia = ' 1'
        ELSE
           CALL atomic_pf(q_n, na, ni, tmp, eneion)
           csign = '+'
           IF(ni+1 .LT. 10) THEN
              WRITE(cnib,'(I1)') ni+1
              WRITE(cnia,'(I2)') -(ni+1)
           ELSE
              WRITE(cnib,'(I2)') ni+1
              WRITE(cnia,'(I3)') -(ni+1)
           END IF
        END IF

        IF(ni+1 .EQ. na) THEN
           q_i(:) = 1d0
        ELSE
           CALL atomic_pf(q_i, na, ni+1, tmp)
        END IF

        fname = TRIM(dname_lnk)//'/k'//cna//cni//'.dat'
        fname_pf = TRIM(dname_lnk)//'/pf'//cna//cni//'.dat'
        WRITE(2,*) TRIM(celem(na))//'1'//csign//TRIM(cnib)//' '//TRIM(celem(na))//'_Ion : '//celem(na)//' 1 e- '//cnia//' #'
        WRITE(2,*) 'f input/'//fname
        WRITE(2,*) ' '
        OPEN(1, FILE=FASTCHEM_INPUT_DIR//fname, RECL=max_line_len)
        OPEN(3, FILE=FASTCHEM_INPUT_DIR//fname_pf, RECL=max_line_len)
        WRITE(1,*) DBLE(tmp_min), DBLE(dtmp), ' log'
        DO l = 1, lmax
           kt = k_bol * tmp(l)
           lnk = 2d0 * (2d0 * pi * m_ele * kt / h**2)**1.5d0 * &
                EXP(-eneion * h * clight / kt) * q_i(l) / q_n(l) * (kt / 1d6)
           lnk = LOG(MAX(lnk,TINY(lnk)))
           IF(ni .EQ. -1) lnk = -lnk ! ANION
           WRITE(3,*) q_n(l), q_i(l), lnk
           lnk_save(l) = lnk + lnk_save(l)
           WRITE(1,*) DBLE(lnk_save(l))
        END DO
        lnk_max = MAX(INT(MAXVAL(DBLE(lnk_save))), lnk_max)
        CLOSE(1)
        CLOSE(3)
     END DO
!     PRINT *, lnk_max
  END DO
#endif
  CLOSE(2)
  PRINT *, ''
  PRINT *, 'Go to FastChem directory and run FastChem:'
  PRINT *, './fastchem input/'//TRIM(fname_config)

  STOP

CONTAINS

  SUBROUTINE atomic_pf(q,na,ni,tmp,eneion)
    USE HDF5
    USE H5LT

    INTEGER(HID_T) :: file_spec, grp_prop, grp_ion
    INTEGER(HSIZE_T), DIMENSION(1) :: dims

    REAL(REAL128), INTENT(OUT) :: q(:)
    INTEGER(INT32), INTENT(IN) :: na ! atomic number
    INTEGER(INT32), INTENT(IN) :: ni ! ionization level (0 = neutral)
    REAL(REAL128), INTENT(IN) :: tmp(:) ! temperature
    REAL(REAL64), INTENT(OUT), OPTIONAL :: eneion
    INTEGER(INT32), ALLOCATABLE :: gtot(:)
    REAL(REAL64), ALLOCATABLE :: ene(:)
    REAL(REAL64) :: mass
    INTEGER(INT32) :: error, j, nstates
    LOGICAL :: exists
    CHARACTER :: natm*3, nion*2
    
    
    ! 読み込みファイルを開く
    CALL h5open_f(error)
    
    !-----------------------------------------------------------------------------
    ! 同位体置換体の質量・存在比・分配関数
    CALL h5Fopen_f(OPTAB_DATABASE_DIR//'h5/NIST.h5', H5F_ACC_RDONLY_F, file_spec, error)
    CALL h5Gopen_f(file_spec, 'prop', grp_prop, error) 
    WRITE(natm,'(I3.3)') na
    WRITE(nion,'(I2.2)') ni
    CALL h5Gopen_f(grp_prop, natm//'.'//nion, grp_ion, error)
!    PRINT *, natm//'.'//nion
    dims = [1]
    CALL h5LTread_dataset_f(grp_ion, 'mass'   , H5T_IEEE_F64LE, mass, dims, error)
    CALL h5LTread_dataset_f(grp_ion, 'nstates', H5T_STD_I32LE , nstates, dims, error)
    IF(PRESENT(eneion)) THEN
       CALL h5LTread_dataset_f(grp_ion, 'eneion', H5T_IEEE_F64LE, eneion, dims, error)
    END IF
    ALLOCATE(ene(nstates), gtot(nstates))
    dims = [nstates]
    CALL h5LTread_dataset_f(grp_ion, 'ene' , H5T_IEEE_F64LE, ene , dims, error)
    CALL h5LTread_dataset_f(grp_ion, 'gtot', H5T_STD_I32LE, gtot, dims, error)
    gtot(:) = MAX(1, gtot(:))
    DO j = 1, UBOUND(tmp,1)
       CALL compute_q(tmp(j), ene, gtot, q(j))
    END DO
    DEALLOCATE(ene, gtot)
    CALL h5Gclose_f(grp_ion, error)
    CALL h5Gclose_f(grp_prop, error)
    CALL h5Fclose_f(file_spec, error)

#ifndef UNUSE_GFGAM
    CALL h5Fopen_f(OPTAB_DATABASE_DIR//'h5/gfgam.h5', H5F_ACC_RDONLY_F, file_spec, error)
    WRITE(natm,'(I3.3)') na
    WRITE(nion,'(I2.2)') ni
    ! ENERGY AND GTOT
    CALL h5Lexists_f(file_spec, natm//nion, exists, error)
    IF(exists) THEN
       CALL h5Gopen_f(file_spec, natm//nion, grp_ion, error)
       dims = [1]
       CALL h5LTread_dataset_f(grp_ion, 'nstates', H5T_STD_I32LE , nstates, dims, error)
       ALLOCATE(ene(nstates), gtot(nstates))
       dims = [nstates]
       CALL h5LTread_dataset_f(grp_ion, 'ene' , H5T_IEEE_F64LE, ene , dims, error)
       CALL h5LTread_dataset_f(grp_ion, 'gtot', H5T_STD_I32LE, gtot, dims, error)
       DO j = 1, UBOUND(tmp,1)
          CALL compute_q(tmp(j), ene, gtot, q(j))
       END DO
       DEALLOCATE(ene, gtot)
       CALL h5Gclose_f(grp_ion, error)
    END IF
    CALL h5Fclose_f(file_spec, error)
#endif

    RETURN
  END SUBROUTINE atomic_pf


  SUBROUTINE compute_q(temp, ene, gtot, q)

    USE const_module, ONLY : c2
    IMPLICIT NONE
    REAL(REAL128), INTENT(IN) :: temp
    REAL(REAL64), INTENT(IN) :: ene(:)
    INTEGER(INT32), INTENT(IN) :: gtot(:)
    REAL(REAL128), INTENT(OUT) :: q
    
    INTEGER(INT64) :: l

    q = 0d0
    DO l = 1, UBOUND(ene,1)
       q = q + gtot(l) * EXP(-c2 * ene(l) / temp)
    END DO

    RETURN
  END SUBROUTINE compute_q
  
         
  LOGICAL FUNCTION exclude(species)
    CHARACTER :: species*16
    exclude = .FALSE.
    IF(&
         TRIM(species) == 'Al1+' .OR. &
         TRIM(species) == 'Ar1+' .OR. &
         TRIM(species) == 'C1+' .OR. &
         TRIM(species) == 'Ca1+' .OR. &
         TRIM(species) == 'Cl1+' .OR. &
         TRIM(species) == 'Co1+' .OR. &
         TRIM(species) == 'Cr1+' .OR. &
         TRIM(species) == 'Cu1+' .OR. &
         TRIM(species) == 'F1+' .OR. &
         TRIM(species) == 'Fe1+' .OR. &
         TRIM(species) == 'Ge1+' .OR. &
         TRIM(species) == 'H1+' .OR. &
         TRIM(species) == 'H1-' .OR. &
         TRIM(species) == 'He1+' .OR. &
         TRIM(species) == 'K1+' .OR. &
         TRIM(species) == 'Mg1+' .OR. &
         TRIM(species) == 'Mn1+' .OR. &
         TRIM(species) == 'N1+' .OR. &
         TRIM(species) == 'Na1+' .OR. &
         TRIM(species) == 'Ne1+' .OR. &
         TRIM(species) == 'Ni1+' .OR. &
         TRIM(species) == 'O1+' .OR. &
         TRIM(species) == 'P1+' .OR. &
         TRIM(species) == 'S1+' .OR. &
         TRIM(species) == 'Si1+' .OR. &
         TRIM(species) == 'Ti1+' .OR. &
         TRIM(species) == 'V1+' .OR. &
         TRIM(species) == 'Zn1+' .OR. &
         TRIM(species) == 'Li1+' .OR. &
         TRIM(species) == 'Be1+' .OR. &
         TRIM(species) == 'B1+' .OR. &
         TRIM(species) == 'Sc1+' .OR. &
         TRIM(species) == 'Ga1+' .OR. &
         TRIM(species) == 'As1+' .OR. &
         TRIM(species) == 'Se1+' .OR. &
         TRIM(species) == 'Rb1+' .OR. &
         TRIM(species) == 'Sr1+' .OR. &
         TRIM(species) == 'Y1+' .OR. &
         TRIM(species) == 'Zr1+' .OR. &
         TRIM(species) == 'Nb1+' .OR. &
         TRIM(species) == 'Mo1+' .OR. &
         TRIM(species) == 'Ru1+' .OR. &
         TRIM(species) == 'Rh1+' .OR. &
         TRIM(species) == 'Pd1+' .OR. &
         TRIM(species) == 'Ag1+' .OR. &
         TRIM(species) == 'Cd1+' .OR. &
         TRIM(species) == 'In1+' .OR. &
         TRIM(species) == 'Sn1+' .OR. &
         TRIM(species) == 'Sb1+' .OR. &
         TRIM(species) == 'Te1+' .OR. &
         TRIM(species) == 'Cs1+' .OR. &
         TRIM(species) == 'Ba1+' .OR. &
         TRIM(species) == 'La1+' .OR. &
         TRIM(species) == 'Ce1+' .OR. &
         TRIM(species) == 'Pr1+' .OR. &
         TRIM(species) == 'Nd1+' .OR. &
         TRIM(species) == 'Sm1+' .OR. &
         TRIM(species) == 'Eu1+' .OR. &
         TRIM(species) == 'Gd1+' .OR. &
         TRIM(species) == 'Tb1+' .OR. &
         TRIM(species) == 'Dy1+' .OR. &
         TRIM(species) == 'Ho1+' .OR. &
         TRIM(species) == 'Er1+' .OR. &
         TRIM(species) == 'Tm1+' .OR. &
         TRIM(species) == 'Yb1+' .OR. &
         TRIM(species) == 'Lu1+' .OR. &
         TRIM(species) == 'Tm1+' .OR. &
         TRIM(species) == 'Hf1+' .OR. &
         TRIM(species) == 'Ta1+' .OR. &
         TRIM(species) == 'W1+' .OR. &
         TRIM(species) == 'Re1+' .OR. &
         TRIM(species) == 'Os1+' .OR. &
         TRIM(species) == 'Ir1+' .OR. &
         TRIM(species) == 'Pt1+' .OR. &
         TRIM(species) == 'Au1+' .OR. &
         TRIM(species) == 'Hg1+' .OR. &
         TRIM(species) == 'Tl1+' .OR. &
         TRIM(species) == 'Pb1+' .OR. &
         TRIM(species) == 'Bi1+' .OR. &
         TRIM(species) == 'Th1+' .OR. &
         TRIM(species) == 'U1+' .OR. &
         TRIM(species) == 'Al1++' .OR. &
         TRIM(species) == 'Ar1++' .OR. &
         TRIM(species) == 'C1++' .OR. &
         TRIM(species) == 'Ca1++' .OR. &
         TRIM(species) == 'Cl1++' .OR. &
         TRIM(species) == 'Co1++' .OR. &
         TRIM(species) == 'Cr1++' .OR. &
         TRIM(species) == 'Cu1++' .OR. &
         TRIM(species) == 'F1++' .OR. &
         TRIM(species) == 'Fe1++' .OR. &
         TRIM(species) == 'Ge1++' .OR. &
         TRIM(species) == 'He1++' .OR. &
         TRIM(species) == 'K1++' .OR. &
         TRIM(species) == 'Mg1++' .OR. &
         TRIM(species) == 'Mn1++' .OR. &
         TRIM(species) == 'N1++' .OR. &
         TRIM(species) == 'Na1++' .OR. &
         TRIM(species) == 'Ne1++' .OR. &
         TRIM(species) == 'Ni1++' .OR. &
         TRIM(species) == 'O1++' .OR. &
         TRIM(species) == 'P1++' .OR. &
         TRIM(species) == 'S1++' .OR. &
         TRIM(species) == 'Si1++' .OR. &
         TRIM(species) == 'Ti1++' .OR. &
         TRIM(species) == 'V1++' .OR. &
         TRIM(species) == 'Zn1++' .OR. &
         TRIM(species) == 'Li1++' .OR. &
         TRIM(species) == 'Be1++' .OR. &
         TRIM(species) == 'B1++' .OR. &
         TRIM(species) == 'Sc1++' .OR. &
         TRIM(species) == 'Ga1++' .OR. &
         TRIM(species) == 'As1++' .OR. &
         TRIM(species) == 'Se1++' .OR. &
         TRIM(species) == 'Rb1++' .OR. &
         TRIM(species) == 'Sr1++' .OR. &
         TRIM(species) == 'Y1++' .OR. &
         TRIM(species) == 'Zr1++' .OR. &
         TRIM(species) == 'Nb1++' .OR. &
         TRIM(species) == 'Mo1++' .OR. &
         TRIM(species) == 'Ru1++' .OR. &
         TRIM(species) == 'Rh1++' .OR. &
         TRIM(species) == 'Pd1++' .OR. &
         TRIM(species) == 'Ag1++' .OR. &
         TRIM(species) == 'Cd1++' .OR. &
         TRIM(species) == 'In1++' .OR. &
         TRIM(species) == 'Sn1++' .OR. &
         TRIM(species) == 'Sb1++' .OR. &
         TRIM(species) == 'Te1++' .OR. &
         TRIM(species) == 'Cs1++' .OR. &
         TRIM(species) == 'Ba1++' .OR. &
         TRIM(species) == 'La1++' .OR. &
         TRIM(species) == 'Ce1++' .OR. &
         TRIM(species) == 'Pr1++' .OR. &
         TRIM(species) == 'Nd1++' .OR. &
         TRIM(species) == 'Sm1++' .OR. &
         TRIM(species) == 'Eu1++' .OR. &
         TRIM(species) == 'Gd1++' .OR. &
         TRIM(species) == 'Tb1++' .OR. &
         TRIM(species) == 'Dy1++' .OR. &
         TRIM(species) == 'Ho1++' .OR. &
         TRIM(species) == 'Er1++' .OR. &
         TRIM(species) == 'Tm1++' .OR. &
         TRIM(species) == 'Yb1++' .OR. &
         TRIM(species) == 'Lu1++' .OR. &
         TRIM(species) == 'Tm1++' .OR. &
         TRIM(species) == 'Hf1++' .OR. &
         TRIM(species) == 'Ta1++' .OR. &
         TRIM(species) == 'W1++' .OR. &
         TRIM(species) == 'Re1++' .OR. &
         TRIM(species) == 'Os1++' .OR. &
         TRIM(species) == 'Ir1++' .OR. &
         TRIM(species) == 'Pt1++' .OR. &
         TRIM(species) == 'Au1++' .OR. &
         TRIM(species) == 'Hg1++' .OR. &
         TRIM(species) == 'Tl1++' .OR. &
         TRIM(species) == 'Pb1++' .OR. &
         TRIM(species) == 'Bi1++' .OR. &
         TRIM(species) == 'Th1++' .OR. &
         TRIM(species) == 'U1++') exclude = .TRUE.
    RETURN
  END FUNCTION exclude


END PROGRAM prep_FastChem
