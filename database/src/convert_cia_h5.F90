
#define DIRH5 '../h5/'

PROGRAM convert_cia_h5

  USE ISO_FORTRAN_ENV
  USE HDF5
  USE H5LT

  IMPLICIT NONE
  
  INTEGER(HID_T) :: file_id, group_id
  INTEGER(HSIZE_T) :: dims(1), dims2(2)
  INTEGER(INT32) :: error, num_te, num_wn0, num_wn, i, j, iostat
  REAL(REAL64), ALLOCATABLE :: wn(:), te(:), coeff(:,:)
  REAL(REAL64) :: temp_old, temp, wn_min, wn_max
  CHARACTER*32 :: species, filename

  CALL h5open_f(error)

  CALL h5Fcreate_f(DIRH5//'CIA_HITRAN.h5', H5F_ACC_TRUNC_F, file_id, error)

  OPEN(1, FILE='list_convert.txt', STATUS='OLD')
  DO
     READ(1,FMT='(A)',IOSTAT=iostat) filename
!     READ(1,*,IOSTAT=iostat) filename
     IF(iostat /= 0) THEN
        PRINT *, 'iostat for list_convert.txt = ', iostat
        EXIT
     END IF
     
     OPEN(2, FILE=TRIM(filename), STATUS='OLD', IOSTAT=iostat)
     IF(iostat /= 0) CYCLE
     PRINT *, TRIM(filename)

     READ(2,*) species, wn_min, wn_max, num_wn, temp
     DO j = 1, num_wn
        READ(2,*)
     END DO
     num_te = 1
     num_wn0 = num_wn
     temp_old = temp
     DO
        READ(2,*,IOSTAT=iostat) species, wn_min, wn_max, num_wn, temp
        IF(iostat /= 0) EXIT
        IF(num_wn /= num_wn0) EXIT
        IF(temp < temp_old) THEN
           PRINT *, '*** ERROR: temp not sorted ***', temp, temp_old
           STOP
        END IF
        num_te = num_te + 1
        DO j = 1, num_wn
           READ(2,*)
        END DO
     END DO
     REWIND(2)
     ALLOCATE(coeff(num_te, num_wn0), te(num_te), wn(num_wn0))
     DO i = 1, num_te
        READ(2,*) species, wn_min, wn_max, num_wn, temp
        te(i) = temp
        DO j = 1, num_wn
           READ(2,*) wn(j), coeff(i,j)
        END DO
     END DO
     CLOSE(2)     

     CALL h5Gcreate_f(file_id, TRIM(species), group_id, error)
     dims = [1]
     CALL h5LTmake_dataset_f(group_id, 'num_wn', 1, dims, H5T_STD_I32LE, num_wn0, error)
     CALL h5LTmake_dataset_f(group_id, 'num_temp', 1, dims, H5T_STD_I32LE, num_te, error)
     dims = [num_te]
     CALL h5LTmake_dataset_f(group_id, 'temp', 1, dims, H5T_IEEE_F64LE, te, error)
     dims = [num_wn]
     CALL h5LTmake_dataset_f(group_id, 'wn', 1, dims, H5T_IEEE_F64LE, wn, error)
     dims2 = [num_te, num_wn]
     CALL h5LTmake_dataset_f(group_id, 'coeff', 2, dims2, H5T_IEEE_F64LE, coeff, error)
     CALL h5Gclose_f(group_id, error)

     DEALLOCATE(coeff, te, wn)
  END DO

  CALL h5Fclose_f(file_id, error)

  CALL h5close_f(error)
  
  STOP

END PROGRAM convert_cia_h5
