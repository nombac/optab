
#define DIRH5 '../h5/'

PROGRAM convert_nist_h5

  USE ISO_FORTRAN_ENV
  USE ISO_C_BINDING
  USE HDF5
  USE H5LT
  USE nist_module, ONLY : get_nist_eneion, get_nist_states, get_nist_mass_comp

  IMPLICIT NONE
  
  INTEGER(HID_T) :: grp_ion, grp_prop, nist_h5, filetype
  INTEGER(HID_T) :: data_id, space_id
  INTEGER(SIZE_T), PARAMETER :: sdim = 8
  INTEGER(SIZE_T), PARAMETER :: sdim2 = 32
  INTEGER(HSIZE_T), DIMENSION(1) :: dims
  ! .states
  INTEGER(INT32), ALLOCATABLE :: gtot(:)
  REAL(REAL64), ALLOCATABLE :: ene(:)
  CHARACTER(LEN=sdim), ALLOCATABLE, TARGET :: term(:)
  CHARACTER(LEN=sdim2), ALLOCATABLE, TARGET :: conf(:)
  TYPE(C_PTR) :: f_ptr
  INTEGER(INT32) :: nstates
  INTEGER(INT32) :: error
  INTEGER(INT32) :: ni, na, indx
  CHARACTER :: natm*3, nion*2, filename*72, command*72
  INTEGER(INT32), PARAMETER :: na_uran = 92, mn_uran = 243
  REAL(REAL64),ALLOCATABLE :: comp_iso(:,:), mass(:), eneion(:)
  
  ALLOCATE(comp_iso(na_uran,mn_uran), mass(na_uran), eneion((na_uran + 1) * 100))
  comp_iso(:,:) = -1d0
  mass(:) = 0d0

  ! 平均質量を計算する
  CALL get_nist_mass_comp(na_uran, mass, comp_iso)
  
  !***********************************************************
  CALL h5open_f(error)
  
  ! NIST データ
  filename = DIRH5//'NIST.h5'
  CALL h5Fcreate_f(TRIM(filename), H5F_ACC_TRUNC_F, nist_h5, error)
  
  CALL h5Gcreate_f(nist_h5, 'prop', grp_prop, error)
  
  ! イオンごとの質量とレベルデータ
  DO na = 1, na_uran ! loop atom
     WRITE(natm,'(I3.3)') na
     DO ni = 0, na - 1 ! loop ion
        WRITE(nion,'(I2.2)') ni
        
        indx = na * 100 + ni
        CALL get_nist_eneion(natm, nion, eneion(indx))
        CALL get_nist_states(natm, nion, gtot, ene, term, conf)
        nstates = UBOUND(gtot,1)
        
        CALL h5Gcreate_f(grp_prop, natm//'.'//nion, grp_ion, error)
        dims = [1]
        CALL h5LTmake_dataset_f(grp_ion, 'mass'   , 1, dims, H5T_IEEE_F64LE, mass(na)    , error)
        CALL h5LTmake_dataset_f(grp_ion, 'nstates', 1, dims, H5T_STD_I32LE , nstates     , error)
        CALL h5LTmake_dataset_f(grp_ion, 'eneion' , 1, dims, H5T_IEEE_F64LE, eneion(indx), error)
        dims = [nstates]
        CALL h5LTmake_dataset_f(grp_ion, 'ene'    , 1, dims, H5T_IEEE_F64LE, ene     , error)
        CALL h5LTmake_dataset_f(grp_ion, 'gtot'   , 1, dims, H5T_STD_I32LE , gtot    , error)
        !---- write an array of strings ---
        ! https://support.hdfgroup.org/ftp/HDF5/examples/examples-by-api/hdf5-examples/1_8/FORTRAN/H5T/h5ex_t_stringC_F03.f90
        CALL H5Tcopy_f(H5T_FORTRAN_S1, filetype, error)
        CALL H5Tset_size_f(filetype, sdim, error)
        CALL h5screate_simple_f(1, dims, space_id, error)
        CALL h5dcreate_f(grp_ion, 'term', filetype, space_id, data_id, error)
        f_ptr = C_LOC(term(1)(1:1))
        CALL H5Dwrite_f(data_id, filetype, f_ptr, error);
        CALL h5Dclose_f(data_id, error)
        CALL h5Sclose_f(space_id, error)
        CALL H5Tclose_f(filetype, error)
        !----
        CALL H5Tcopy_f(H5T_FORTRAN_S1, filetype, error)
        CALL H5Tset_size_f(filetype, sdim2, error)
        CALL h5screate_simple_f(1, dims, space_id, error)
        CALL h5dcreate_f(grp_ion, 'conf', filetype, space_id, data_id, error)
        f_ptr = C_LOC(conf(1)(1:1))
        CALL H5Dwrite_f(data_id, filetype, f_ptr, error);
        CALL h5Dclose_f(data_id, error)
        CALL h5Sclose_f(space_id, error)
        CALL H5Tclose_f(filetype, error)
        !----
        CALL h5Gclose_f(grp_ion, error)
        
        DEALLOCATE(gtot, ene, term, conf)
     END DO
  END DO
  CALL h5Gclose_f(grp_prop, error)
  
  CALL h5Fclose_f(nist_h5, error)
  print *, 'output file: '
  WRITE(command,"('ls -l ',A)") TRIM(filename)
  CALL SYSTEM(TRIM(command))
  
  CALL h5close_f(error)
  
  DEALLOCATE(comp_iso, mass, eneion)
  STOP

END PROGRAM convert_nist_h5
