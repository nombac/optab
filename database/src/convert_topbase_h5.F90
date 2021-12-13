
#define DIRH5 '../h5/'
#define DIRTOPBASE './'

! elevel files MUST BE downloaded with
! - General energy order
! - ONLY Energy (Ryd) wrt ground state and Statistical weight BEING CHECKED

! xsectn files MUST BE downloaded with
! - General energy order

PROGRAM convert_topbase_h5

  USE const_module, ONLY : megabarn_to_cm2, ryd_to_wnum
  USE ISO_FORTRAN_ENV
  USE ISO_C_BINDING
  USE HDF5
  USE H5LT
  
  IMPLICIT NONE
  
  INTEGER(HID_T) :: grp_lvl, grp_spec, file_h5
  INTEGER(HSIZE_T), DIMENSION(1) :: dims
  INTEGER(INT32) :: error
  INTEGER(INT32) :: nz, ne, nl, ios, np, nl_max, np_max
  INTEGER(INT32) :: i0, nz0, ne0, islp0, ilv0
  REAL(REAL64) :: ene0, gtot0
  REAL(REAL64), ALLOCATABLE :: ene(:), gtot(:), wnum(:), cross(:)
  INTEGER(INT32), ALLOCATABLE :: islp(:), ilv(:)
  CHARACTER(LEN=64), PARAMETER :: efmt = "(1X,I6,1X,I2,1X,I2,1X,I4,1X,I3,1X,E12.5,1X,F4.1)"
  CHARACTER(LEN=72) :: filename, command
  
  !***********************************************************
  CALL h5open_f(error)
  
  ! NIST データ
  filename = DIRH5//'TOPbase.h5'
  CALL h5Fcreate_f(TRIM(filename), H5F_ACC_TRUNC_F, file_h5, error)
  
  DO nz = 1, 26
     IF(nz==15 .OR. nz==17 .OR. nz==19 .OR. (nz >=21 .AND. nz <=25)) CYCLE ! NO DATA IN TOPBASE
     
     DO ne = 1, nz

        ! energy level
        OPEN(1,FILE=DIRTOPBASE//'/elevel.'//c02(nz)//'.'//c02(ne)//'.dat', STATUS='OLD',IOSTAT=ios)
        IF(ios /= 0) CYCLE
        
        CALL h5Gcreate_f(file_h5, c02(nz)//'.'//c02(ne), grp_spec, error)
        
        READ(1,*)
        READ(1,*)
        READ(1,*)
        nl_max = 0
        DO
           READ(1,efmt,IOSTAT=ios) i0, nz0, ne0, islp0, ilv0, ene0, gtot0
           IF(ios /= 0) EXIT
           nl_max = nl_max + 1
        END DO
        
        dims = [1]
        CALL h5LTmake_dataset_f(grp_spec, 'nlevel', 1, dims, H5T_STD_I32LE, nl_max, error)
        ALLOCATE(ene(nl_max), gtot(nl_max), islp(nl_max), ilv(nl_max))
        PRINT *, 'nz = ', nz, ' ne = ', ne, ' nlevel = ', nl_max
        
        REWIND(1)
        
        READ(1,*)
        READ(1,*)
        READ(1,*)
        DO nl = 1, nl_max
           READ(1,efmt) i0, nz0, ne0, islp(nl), ilv(nl), ene(nl), gtot(nl)
        END DO
        CLOSE(1)
        
        !CONVERSION FROM RYD TO WAVENUMBER
        ene = ene * ryd_to_wnum
        
        dims = [nl_max]
        CALL h5LTmake_dataset_f(grp_spec, 'ene'  , 1, dims, H5T_IEEE_F64LE, ene      , error)
        CALL h5LTmake_dataset_f(grp_spec, 'gtot' , 1, dims, H5T_STD_I32LE , INT(gtot), error)
        
        ! cross sections
        OPEN(1,FILE=DIRTOPBASE//'/xsectn.'//c02(nz)//'.'//c02(ne)//'.dat', STATUS='OLD',IOSTAT=ios)
        IF(ios /= 0) GOTO 999
        
        READ(1,*)
        READ(1,*)
        READ(1,*)
        DO nl = 1, nl_max
           CALL h5Gcreate_f(grp_spec, c04(nl), grp_lvl, error)
           READ(1,*) i0, nz0, ne0, islp0, ilv0, ene0, np_max
           IF((islp0 /= islp(nl)) .OR. (ilv0 /= ilv(nl))) THEN
              PRINT *, 'error: ', islp0, islp(nl), ilv0, ilv(nl), nz, ne
              STOP
           END IF
           dims = [1]
           CALL h5LTmake_dataset_f(grp_lvl, 'np', 1, dims, H5T_STD_I32LE, np_max, error)
           IF(np_max == 0) CYCLE
           ALLOCATE(wnum(np_max), cross(np_max))
           DO np = 1, np_max
              READ(1,*) wnum(np), cross(np)
           END DO
           !CONVERSION FROM RYD TO WAVENUMBER
           wnum = wnum * ryd_to_wnum
           !CONVERSION FROM MBARN TO CM^2
           cross = cross * megabarn_to_cm2
           dims = [np_max]
           CALL h5LTmake_dataset_f(grp_lvl, 'wnum'  , 1, dims, H5T_IEEE_F64LE, wnum, error)
           CALL h5LTmake_dataset_f(grp_lvl, 'cross' , 1, dims, H5T_IEEE_F64LE, cross, error)
           CALL h5Gclose_f(grp_lvl, error)
           DEALLOCATE(wnum, cross)
        END DO
        
        CLOSE(1)

999     CONTINUE
        DEALLOCATE(ene, gtot, islp, ilv)
        CALL h5Gclose_f(grp_spec, error)
        
     END DO
  END DO
  CALL h5Fclose_f(file_h5, error)

  print *, 'output file: '
  WRITE(command,"('ls -l ',A)") TRIM(filename)
  CALL SYSTEM(TRIM(command))
  
  CALL h5close_f(error)
  
  STOP
  
CONTAINS
  
  FUNCTION c02(n)
    CHARACTER(LEN=2) :: c02
    INTEGER(INT32), INTENT(IN) :: n
    
    WRITE(c02,'(I2.2)') n
    
    RETURN
  END FUNCTION c02
  
  FUNCTION c04(n)
    CHARACTER(LEN=4) :: c04
    INTEGER(INT32), INTENT(IN) :: n
    
    WRITE(c04,'(I4.4)') n
    
    RETURN
  END FUNCTION c04
  
END PROGRAM convert_topbase_h5
