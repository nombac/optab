
PROGRAM convert_Phoenix

  USE ISO_FORTRAN_ENV
  USE HDF5
  USE H5LT
  USE const_module, ONLY : k_bol
  
  IMPLICIT NONE

  CHARACTER*32768 :: linebuf
  INTEGER, PARAMETER :: imax = 2048
  CHARACTER*16 :: b(imax)
  CHARACTER*16, ALLOCATABLE :: species(:)
  INTEGER(INT32) :: i, i0, iostat, j, jmax, k, kmax, error
  INTEGER(INT32), ALLOCATABLE :: id(:)
  REAL(REAL64) :: x(imax), temp0, pres0, temp0_old, pres0_old
  REAL(REAL64), ALLOCATABLE :: n(:,:), pres(:), temp(:), rho(:)
  INTEGER(HID_T) :: file_id
  INTEGER(HSIZE_T), DIMENSION(1) :: dim1
  INTEGER(HSIZE_T), DIMENSION(2) :: dim2
  INTEGER(HID_T) :: filetype
  INTEGER(HID_T) :: data_id, space_id
  INTEGER(SIZE_T), PARAMETER :: sdim = 16
  CHARACTER(LEN=sdim), ALLOCATABLE, TARGET :: term(:)
  TYPE(C_PTR) :: f_ptr
!!$  INTEGER(INT32) :: count, it, ip, itsize, ipsize, modeos
!!$  REAL(REAL32) :: t2
!!$  REAL(REAL64), ALLOCATABLE :: ttab(:), ptab(:), rhotab(:), rtab(:), &
!!$       rosst(:,:), cap_planck(:,:), cap_planck2(:,:), drho(:,:)
  INTEGER(INT32) :: length, status
  CHARACTER(:), ALLOCATABLE :: filename
  INTRINSIC :: COMMAND_ARGUMENT_COUNT, GET_COMMAND_ARGUMENT
  
  CALL GET_COMMAND_ARGUMENT(1, LENGTH=length, STATUS=status)
  IF(status == 0) THEN
     ALLOCATE(CHARACTER(length) :: filename)
     CALL GET_COMMAND_ARGUMENT(1, filename, STATUS=status)
     IF(status == 0) THEN
        PRINT *, 'input file = ', filename
     END IF
  ELSE
     PRINT *, 'argument error'
     STOP
  END IF
  
  OPEN(1, FILE=filename, STATUS='OLD')
  READ(1,'(A)') linebuf

  ! count # of layers
  jmax = 0
  temp0_old = -1d0
  pres0_old = -1d0
  DO
     READ(1,*,IOSTAT=iostat) temp0, pres0
     IF(temp0 == temp0_old .AND. pres0 == pres0_old) EXIT
     temp0_old = temp0
     pres0_old = pres0
     IF(iostat /= 0) EXIT
     jmax = jmax + 1
  END DO

  PRINT *, 'Number of layers = ', jmax

  ! count # of columns
  i0 = 1
  DO
     READ(linebuf,*,IOSTAT=iostat) (b(i), i = 1, i0)
     IF(iostat /= 0) EXIT
     i0 = i0 + 1
  END DO
  i0 = i0 - 1
!!$  OPEN(10,FILE='species_aces_list.dat')
!!$  OPEN(11,FILE='species_aces_list2.dat')
!!$  write(10,*) '# of species = ', i0
!!$  DO i = 1, i0
!!$     write (10,fmt='(a)', advance='no') TRIM(b(i))//', '
!!$     write (11,fmt='(a)', advance='yes') TRIM(b(i))//', '
!!$  END DO
!!$  CLOSE(10)
!!$  CLOSE(11)

  ! read OPAC id list
  OPEN(2, FILE='list_id_Phoenix.txt', STATUS='OLD')
  kmax = 0
  DO
     READ(2,*,IOSTAT=iostat)
     IF(iostat /= 0) EXIT
     kmax = kmax + 1
  END DO
  REWIND(2)
  ALLOCATE(species(kmax), id(kmax))
  DO k = 1, kmax
     READ(2,*) species(k), id(k)
  END DO
  CLOSE(2)


  ! fill output arrays
#define NUMBER_SPECIES_ARRAY 11000  
  ALLOCATE(pres(jmax), temp(jmax), rho(jmax), n(NUMBER_SPECIES_ARRAY,jmax), term(NUMBER_SPECIES_ARRAY))
  term(:) = ""
  term(1) = "total"
  DO k = 1, kmax
     term(id(k)) = species(k)
  END DO
  n = 0d0
  REWIND(1)
  READ(1,*)
  DO j = 1, jmax
     READ(1,*) (x(i), i = 1, i0)
     pres(j) = x(2)
     temp(j) = x(1)
     rho(j) = x(3)
     print *, 'layer: ', j, ' T=', temp(j), ' P=', pres(j), ' rho=', rho(j)
     n(1,j) = pres(j) / (k_bol * temp(j))
     DO k = 1, kmax
        DO i = 1, i0
           IF(TRIM(b(i)) == TRIM(species(k))) THEN
              n(id(k),j) = x(i)
           END IF
        END DO
     END DO
  END DO

  CLOSE(1)

!!$  OPEN(19, FILE='fort.19', FORM='UNFORMATTED', STATUS='OLD', CONVERT='BIG_ENDIAN')
!!$  READ(19) itsize, ipsize, modeos
!!$  READ(19) t2
!!$  ALLOCATE(ttab(itsize))
!!$  ALLOCATE(ptab(ipsize))
!!$  ALLOCATE(rhotab(ipsize))
!!$  ALLOCATE(rtab(ipsize))
!!$  ALLOCATE(rosst(ipsize,itsize))
!!$  ALLOCATE(cap_planck(ipsize,itsize))
!!$  ALLOCATE(cap_planck2(ipsize,itsize))
!!$  ALLOCATE(drho(ipsize,itsize))
!!$
!!$  READ(19) ttab
!!$  READ(19) ptab
!!$  READ(19) rhotab
!!$  READ(19) rtab
!!$  READ(19) rosst
!!$  READ(19) cap_planck
!!$  READ(19) cap_planck2
!!$  READ(19) drho
!!$
!!$  count = 1
!!$  DO it = 1, itsize
!!$     DO ip = 1, ipsize
!!$        PRINT *, 'b: ', count, ttab(it), ptab(ip), rhotab(ip), rtab(ip), drho(ip,it)
!!$        count = count + 1
!!$     END DO
!!$  ENDDO
!!$
!!$  CLOSE(19)

  ! hdf5 output
  CALL h5open_f(error)
  CALL h5Fcreate_f('eos_'//filename(1:INDEX(filename,'_')-1)//'.h5', H5F_ACC_TRUNC_F, file_id, error)
  dim1 = [1]
  CALL h5LTmake_dataset_f(file_id, 'n_layer', 1, dim1, H5T_STD_I32LE, jmax, error)
  CALL h5LTmake_dataset_f(file_id, 'n_species', 1, dim1, H5T_STD_I32LE, NUMBER_SPECIES_ARRAY, error)
  dim1 = [jmax]
  CALL h5LTmake_dataset_f(file_id, 'temp', 1, dim1, H5T_IEEE_F64LE, temp, error)
  CALL h5LTmake_dataset_f(file_id, 'pres', 1, dim1, H5T_IEEE_F64LE, pres, error)
  CALL h5LTmake_dataset_f(file_id, 'rho', 1, dim1, H5T_IEEE_F64LE, rho, error)
  dim2 = [NUMBER_SPECIES_ARRAY,jmax]
  CALL h5LTmake_dataset_f(file_id, 'ndens', 2, dim2, H5T_IEEE_F64LE, n, error)

  dim1 = [NUMBER_SPECIES_ARRAY]
  CALL H5Tcopy_f(H5T_FORTRAN_S1, filetype, error)
  CALL H5Tset_size_f(filetype, sdim, error)
  CALL h5screate_simple_f(1, dim1, space_id, error)
  CALL h5dcreate_f(file_id, 'species', filetype, space_id, data_id, error)
  f_ptr = C_LOC(term(1)(1:1))
  CALL H5Dwrite_f(data_id, filetype, f_ptr, error);
  CALL h5Dclose_f(data_id, error)
  CALL h5Sclose_f(space_id, error)
  CALL H5Tclose_f(filetype, error)
  
  CALL h5Fclose_f(file_id, error)
  CALL h5close_f(error)
    
  DEALLOCATE(pres, temp, rho, n, species, id)

  STOP

END PROGRAM convert_Phoenix
