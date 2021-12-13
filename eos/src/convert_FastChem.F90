
PROGRAM convert_FastChem

  USE ISO_FORTRAN_ENV
  USE HDF5
  USE H5LT
  USE const_module, ONLY : k_bol, amu
  
  IMPLICIT NONE

  CHARACTER*32768 :: linebuf, linebuf0
  INTEGER, PARAMETER :: imax = 2048
  CHARACTER*16 :: b(imax)
  CHARACTER*16, ALLOCATABLE :: species(:)
  INTEGER(INT32) :: i, i0, iostat, j, jmax, k, kmax, error
  INTEGER(INT32), ALLOCATABLE :: id(:)
  REAL(REAL64) :: x(imax), dummy
  REAL(REAL64), ALLOCATABLE :: n(:,:), pres(:), temp(:), rho(:)
  INTEGER(HID_T) :: file_id
  INTEGER(HSIZE_T), DIMENSION(1) :: dim1
  INTEGER(HSIZE_T), DIMENSION(2) :: dim2
  INTEGER(HID_T) :: filetype
  INTEGER(HID_T) :: data_id, space_id
  INTEGER(SIZE_T), PARAMETER :: sdim = 16
  CHARACTER(LEN=sdim), ALLOCATABLE, TARGET :: term(:)
  TYPE(C_PTR) :: f_ptr
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
  READ(1,'(A)') linebuf0

  ! count # of layers
  jmax = 0
  DO
     READ(1,*,IOSTAT=iostat)
     IF(iostat /= 0) EXIT
     jmax = jmax + 1
  END DO

  linebuf = linebuf0(INDEX(linebuf0,"e-"):)

  ! count # of columns
  i0 = 1
  DO
     READ(linebuf,*,IOSTAT=iostat) (b(i), i = 1, i0)
     IF(iostat /= 0) EXIT
     i0 = i0 + 1
  END DO
  i0 = i0 - 1

  ! read OPAC id list
  OPEN(2, FILE='list_id_FastChem.txt', STATUS='OLD')
  kmax = 0
  DO
     READ(2,*,IOSTAT=iostat)
     IF(iostat /= 0) EXIT
     kmax = kmax + 1
  END DO
  kmax = kmax - 1
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
     print *, 'layer: ', j
     READ(1,*) pres(j), temp(j), dummy, n(1,j), rho(j), (x(i), i = 1, i0)
#define BAR (1e6)     
     pres(j) = pres(j) * BAR
     rho(j) = (amu * rho(j)) * n(1,j)
     DO k = 1, kmax
        DO i = 1, i0
           IF(TRIM(b(i)) == TRIM(species(k))) THEN
#ifdef MIXING_RATIO_INSTEAD_OF_NUMBER_DENSITY
              n(id(k),j) = x(i) * x(4)
#else              
              n(id(k),j) = x(i)
#endif              
           END IF
        END DO
     END DO
  END DO

  CLOSE(1)

  ! hdf5 output
  CALL h5open_f(error)
  CALL h5Fcreate_f(filename(INDEX(filename,'/',.TRUE.)+1:INDEX(filename,'.',.TRUE.))//'h5', H5F_ACC_TRUNC_F, file_id, error)
  dim1 = [1]
  CALL h5LTmake_dataset_f(file_id, 'n_layer', 1, dim1, H5T_STD_I32LE, jmax, error)
  CALL h5LTmake_dataset_f(file_id, 'n_species', 1, dim1, H5T_STD_I32LE, NUMBER_SPECIES_ARRAY, error)
  dim1 = [jmax]
  CALL h5LTmake_dataset_f(file_id, 'temp', 1, dim1, H5T_IEEE_F64LE, temp, error)
  CALL h5LTmake_dataset_f(file_id, 'pres', 1, dim1, H5T_IEEE_F64LE, pres, error)
  CALL h5LTmake_dataset_f(file_id, 'rho' , 1, dim1, H5T_IEEE_F64LE, rho , error)
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

END PROGRAM convert_FastChem
