
MODULE input_module

  USE ISO_FORTRAN_ENV
  IMPLICIT NONE

CONTAINS  

  SUBROUTINE read_mol_source_list(n_species)
    INTEGER(INT32), INTENT(OUT) :: n_species
    INTEGER(INT32) :: n, i, snum, check_total
    INTEGER(INT32), ALLOCATABLE :: check(:)
    CHARACTER*256, ALLOCATABLE :: source(:)
    CHARACTER*256 :: species, isotopol, linebuf, source0, source1
    INTEGER(INT32) :: iostat

    OPEN(2, FILE='input/mol_source.dat')
    OPEN(1, FILE='input/species_id.dat', STATUS='OLD')
    READ(1,*)
    n_species = 0
    main_loop: DO
       n = 1
       check_total = 0
       DO
          ALLOCATE(check(n), source(n))
          check(:) = 0
          READ(1,FMT='(A)') linebuf
          READ(linebuf,*,IOSTAT=iostat) snum, species, isotopol, (check(i), source(i), i = 1, n)
          IF(iostat /= 0) THEN
             SELECT CASE (SUM(check))
             CASE (0)
!                PRINT *, '+++ WARNING: no source is selected for ', TRIM(isotopol)
                EXIT
             CASE(1)
                source1 = TRIM(isotopol)//'__'//TRIM(source0)
                WRITE(2,*) TRIM(source1)
                n_species = n_species + 1
                EXIT
             CASE(2:)
                PRINT *, '*** ERROR: multiple sources are selected for ', TRIM(species)
                PRINT *, TRIM(linebuf)
                STOP
             END SELECT
          END IF
          IF(snum == 999) EXIT main_loop
          IF(check(n) == 1) THEN
             source0 = source(n)
          END IF
          BACKSPACE(1)
          DEALLOCATE(check, source)
          n = n + 1
       END DO
       DEALLOCATE(check, source)
    END DO main_loop
    DEALLOCATE(check, source)
    CLOSE(1)
    CLOSE(2)
    
    RETURN

  END SUBROUTINE read_mol_source_list


!#define GRID_CONST
#define GRID_LOG_CONST
#ifdef PHOENIX_GRID
#undef GRID_LOG_CONST
#endif
  SUBROUTINE read_grid(grd, dgrd)
    USE ISO_FORTRAN_ENV
    USE mpi_module, ONLY : para_range, kprc, myrk_k
    IMPLICIT NONE
    REAL(REAL64), ALLOCATABLE, INTENT(OUT) :: grd(:)
    REAL(REAL64), ALLOCATABLE, INTENT(OUT) :: dgrd(:)
    
    INTEGER(INT32) :: k, ks, ke, k_total

#ifdef GRID_CONST
    REAL(REAL64) :: grd_min, grd_max, dgrd0
#endif    
#ifdef GRID_LOG_CONST
    REAL(REAL64) :: grd_min, grd_max, dgrd0
    NAMELIST /grid_log_const/ k_total, grd_min, grd_max
#endif    
#ifdef PHOENIX_GRID
    REAL(REAL64), ALLOCATABLE :: grd_tot(:), dgrd_tot(:)
    REAL(REAL64) :: lam
#endif

#ifdef GRID_CONST
    grd_min = 1d2
    grd_max = 1d7
    dgrd0 = 10d0
    k_total = NINT((grd_max - grd_min) / dgrd0) + 1
    CALL para_range(1, k_total, kprc, myrk_k, ks, ke)
    ALLOCATE(grd(ks:ke), dgrd(ks:ke))

    DO k = ks, ke
       grd(k) = grd_min + dgrd0 * (k - 1)
    END DO
    DO k = ks, ke
       dgrd(k) = dgrd0
    END DO
#endif
#ifdef GRID_LOG_CONST
    OPEN(5, FILE='input/fort.5', STATUS='OLD')
    READ(5, NML=grid_log_const)
    CLOSE(5)
    k_total = k_total + 1
!    k_total = 10000000 + 1
!    grd_min = 1d0
!    grd_max = 7d0
    dgrd0 = (grd_max - grd_min) / (k_total - 1)
    CALL para_range(1, k_total, kprc, myrk_k, ks, ke)
    ALLOCATE(grd(k_total), dgrd(k_total))

    DO k = ks, ke
       grd(k) = 10d0**(grd_min + dgrd0 * (k - 1) + 0.5d0 * dgrd0)
    END DO
    DO k = ks, ke
       dgrd(k) = 10d0**(grd_min + dgrd0 * (k - 1) + dgrd0) - 10d0**(grd_min + dgrd0 * (k - 1))
    END DO
#endif    
#ifdef PHOENIX_GRID
    OPEN(2, FILE='input/grid_phoenix.dat', STATUS='OLD')
    READ(2,*) k_total
    IF(k_total <= 1) THEN
       PRINT *, '*** ERROR: phoenix grid points number <= 1'
       STOP
    END IF
    CALL para_range(1, k_total, kprc, myrk_k, ks, ke)
    ALLOCATE(grd(ks:ke), dgrd(ks:ke), grd_tot(k_total), dgrd_tot(k_total))
    DO k = 1, k_total
       READ(2,*) lam
       grd_tot(k) = 1d0 / (lam * 1d-8)
    END DO
    CLOSE(2)
    
    DO k = 1, k_total-1
       dgrd_tot(k) = grd_tot(k+1) - grd_tot(k)
    END DO
    dgrd_tot(k_total) = dgrd_tot(k_total-1)
    DO k = 1, k_total
       IF(k >= ks .AND. k <= ke) THEN
          grd(k) = grd_tot(k)
          dgrd(k) = dgrd_tot(k)
       END IF
    END DO

    DEALLOCATE(grd_tot, dgrd_tot)
#endif
    RETURN
    
  END SUBROUTINE read_grid
  

#define NUMBER_SPECIES_ARRAY 11000
  SUBROUTINE read_eos(np, temp, rho)
    USE ISO_FORTRAN_ENV
    USE HDF5
    USE H5LT
    USE code_module
    USE const_module, ONLY : k_bol
    IMPLICIT NONE
    REAL(REAL64), ALLOCATABLE, INTENT(OUT) :: np(:,:) ! partial density [cm^{-3}]
    REAL(REAL64), ALLOCATABLE, INTENT(OUT) :: temp(:) ! temperature
    REAL(REAL64), ALLOCATABLE, INTENT(OUT) :: rho(:) ! mass density [g cm^{-2}]
    INTEGER(INT32) :: jmax
#ifndef CROSS_CHECK    
    INTEGER(INT32) :: imax, error
    INTEGER(HID_T) :: file_id
    INTEGER(HSIZE_T), DIMENSION(1) :: dim1
    INTEGER(HSIZE_T), DIMENSION(2) :: dim2
#else
    INTEGER(INT32) :: code0, j
    INTEGER(INT32), PARAMETER :: jmax0 = 128
    REAL(REAL64) :: temp0(jmax0) = 0d0
    NAMELIST /cross_check/ code0, temp0
#endif
    INTEGER(INT32) :: codes, i, ios
    NAMELIST /select_species/ codes
    
#ifdef CROSS_CHECK
    OPEN(5, FILE='input/fort.5', STATUS='OLD')
    READ(5, NML=cross_check)
    CLOSE(5)

    jmax = 0
    DO j = 1, jmax0
       IF(temp0(j) == 0d0) EXIT
       jmax = jmax + 1
    END DO
    
    ALLOCATE(temp(jmax), rho(jmax), np(NUMBER_SPECIES_ARRAY,jmax))
    temp(1:jmax) = temp0(1:jmax)
    rho(1:jmax) = 1d0
    np(:,1:jmax) = 0d0

    np(code0,1:jmax) = 1d0
#else
    CALL h5open_f(error)
    CALL h5Fopen_f('input/eos.h5', H5F_ACC_RDONLY_F, file_id, error)
    dim1 = [1]
    CALL h5LTread_dataset_f(file_id, 'n_layer', H5T_STD_I32LE, jmax, dim1, error)
    CALL h5LTread_dataset_f(file_id, 'n_species', H5T_STD_I32LE, imax, dim1, error)

    ALLOCATE(np(imax,jmax), temp(jmax), rho(jmax))
    dim1 = [jmax]
    CALL h5LTread_dataset_f(file_id, 'temp', H5T_IEEE_F64LE, temp, dim1, error)

    CALL h5Eset_auto_f(0,error)
    CALL h5LTread_dataset_f(file_id, 'rho', H5T_IEEE_F64LE, rho, dim1, error)
    IF(error /= 0) THEN
       print *, 'WARNING: no rho data in EOS => zero is substituted'
       rho(:) = 0d0
    END IF
    CALL h5Eset_auto_f(1,error)
    
    dim2 = [imax,jmax]
    CALL h5LTread_dataset_f(file_id, 'ndens', H5T_IEEE_F64LE, np, dim2, error)
    CALL h5Fclose_f(file_id, error)
    CALL h5close_f(error)

    OPEN(5, FILE='input/fort.5', STATUS='OLD')
    READ(5, NML=select_species, IOSTAT=ios)
    CLOSE(5)
    IF(ios == 0) THEN
       DO i = 1, imax
          IF(i /= codes) np(codes,1:jmax) = 0d0
       END DO
       PRINT *, '*** CAUTION: ONLY A SINGLE SPECIES IS SELECTED'
    END IF
#endif
    RETURN
    
  END SUBROUTINE read_eos

END MODULE input_module
