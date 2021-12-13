
#define DIRH5 '../h5/'

PROGRAM convert_gfgam_h5

  USE const_module, ONLY : megabarn_to_cm2, ryd_to_wnum
  USE ISO_FORTRAN_ENV
  USE ISO_C_BINDING
  USE HDF5
  USE H5LT
  
  IMPLICIT NONE
  
  INTEGER(HID_T) :: grp_spec, file_h5
  INTEGER(HSIZE_T), DIMENSION(1) :: dims
  INTEGER(INT32) :: error, flag, count
  INTEGER(INT32) :: nz, ne, l, ios, l_max, nc
  REAL(REAL64) :: ene0, gtot0
  REAL(REAL64), ALLOCATABLE :: ene(:), gtot(:)
  CHARACTER(LEN=1024) :: linebuf
  CHARACTER(LEN=19), PARAMETER :: fmt = "(13X,F11.3,1X,F4.1)"
  CHARACTER(LEN=72) :: command, output
  
  !***********************************************************
  CALL h5open_f(error)
  
  ! NIST データ
  output = DIRH5//'gfgam.h5'
  CALL h5Fcreate_f(TRIM(output), H5F_ACC_TRUNC_F, file_h5, error)
  
  NZ_LOOP: DO nz = 1, 56
!  DO nz = 24, 24
     
     NE_LOOP: DO ne = nz, 1, -1
!     DO ne = 17, 17

        nc = nz - ne

        ! energy level
        OPEN(1,FILE='atoms/gf'//c02(nz)//c02(nc)//'.gam', STATUS='OLD',IOSTAT=ios)
        PRINT *, 'atoms/gf'//c02(nz)//c02(nc)//'.gam', ios
        IF(ios /= 0) CYCLE

        count = 0
        DO
           count = count + 1
           READ(1,'(a)',IOSTAT=ios) linebuf
           IF(ios /= 0) THEN 
              PRINT *, 'ERROR: ios=', ios
              CYCLE NE_LOOP
           END IF
           IF(INDEX(linebuf,'ELEM') /= 0) EXIT
        END DO

        l_max = 0
        flag = 0
        DO
           READ(1,fmt,IOSTAT=ios) ene0, gtot0
           IF(ios /= 0) THEN
              PRINT *, 'ERROR:', ios, c02(nz)//c02(nc), ene0, gtot0
              CYCLE NE_LOOP
           END IF
           
           IF(ene0 == 0d0 .AND. flag == 0) THEN
              flag = 1
           ELSE IF(ene0 == 0d0 .AND. flag == 1) THEN
              EXIT
           END IF
           l_max = l_max + 1
        END DO
        PRINT *, c02(nz)//c02(nc), ': l_max=', l_max, ene0, gtot0

        CALL h5Gcreate_f(file_h5, c03(nz)//c02(nc), grp_spec, error)
        dims = [1]
        CALL h5LTmake_dataset_f(grp_spec, 'nstates', 1, dims, H5T_STD_I32LE, l_max, error)
        
        ALLOCATE(ene(l_max), gtot(l_max))

        REWIND(1)

        DO
           READ(1,'(a)',IOSTAT=ios) linebuf
           IF(INDEX(linebuf,'ELEM') /= 0) EXIT
        END DO
        DO l = 1, l_max
           READ(1,fmt,IOSTAT=ios) ene0, gtot0
           ene(l) = ABS(ene0)
           gtot(l) = 2d0 * gtot0 + 1d0
!!$           IF(nz == 3 .AND. ne == 3) THEN
!!$              PRINT *, nz, ne, l, ene0, gtot0, gtot(l)
!!$           END IF
        END DO
        dims = [l_max]
        CALL h5LTmake_dataset_f(grp_spec, 'ene'  , 1, dims, H5T_IEEE_F64LE, ene      , error)
        CALL h5LTmake_dataset_f(grp_spec, 'gtot' , 1, dims, H5T_STD_I32LE , INT(gtot), error)

        CLOSE(1)

        DEALLOCATE(ene, gtot)

        CALL h5Gclose_f(grp_spec, error)

     END DO NE_LOOP
  END DO NZ_LOOP
  CALL h5Fclose_f(file_h5, error)
  print *, 'output file: '
  WRITE(command,"('ls -l ',A)") TRIM(output)
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
  
  FUNCTION c03(n)
    CHARACTER(LEN=3) :: c03
    INTEGER(INT32), INTENT(IN) :: n
    
    WRITE(c03,'(I3.3)') n
    
    RETURN
  END FUNCTION c03
  
  FUNCTION c04(n)
    CHARACTER(LEN=4) :: c04
    INTEGER(INT32), INTENT(IN) :: n
    
    WRITE(c04,'(I4.4)') n
    
    RETURN
  END FUNCTION c04
  
END PROGRAM convert_gfgam_h5
