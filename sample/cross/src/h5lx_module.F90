#define DATA_XFER_MODE H5FD_MPIO_INDEPENDENT_F
!#define DATA_XFER_MODE H5FD_MPIO_COLLECTIVE_F
MODULE h5LX_module

  USE ISO_FORTRAN_ENV
  USE HDF5
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: h5LXset_extendible_dataset, h5LXset_dataset, &
       h5LXwrite_extendible_dataset_F64, h5LXwrite_extendible_dataset_F32, h5LXwrite_extendible_dataset_I32, &
       h5LXread_hyperslab_dataset_F64, h5LXread_hyperslab_dataset_F32, h5LXread_hyperslab_dataset_I32, &
       h5LXwrite_hyperslab_dataset_F64, h5LXwrite_hyperslab_dataset_F32, h5LXwrite_hyperslab_dataset_I32, &
       h5LXread_hyperslab_dataset_I16, h5LXopen_file

CONTAINS

  SUBROUTINE h5LXopen_file(filename, access_flag, file_spec, error, comm, info)
    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER(INT32), INTENT(IN) :: access_flag     ! File access flags  
    INTEGER(HID_T), INTENT(OUT) :: file_spec
    INTEGER(INT32), INTENT(INOUT) :: error
    INTEGER(INT32), INTENT(IN), OPTIONAL :: comm
    INTEGER(INT32), INTENT(IN), OPTIONAL :: info

    INTEGER(HID_T) :: plist_id

    CALL h5pcreate_F(H5P_FILE_ACCESS_F, plist_id, error)
    CALL h5pset_fapl_mpio_f(plist_id, comm, info, error)
    CALL h5Fopen_f(filename, access_flag, file_spec, error, plist_id)
    CALL h5Pclose_f(plist_id, error)
    
    RETURN    
  END SUBROUTINE h5LXopen_file


  SUBROUTINE h5LXset_extendible_dataset(dataset, location, name, datatype, chunk, gzip)
   
    INTEGER(HID_T), INTENT(OUT) :: dataset
    INTEGER(HID_T), INTENT(IN) :: location
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER(HID_T), INTENT(IN) :: datatype
    INTEGER(HSIZE_T), INTENT(IN) :: chunk(:)
    INTEGER(INT32), INTENT(IN) :: gzip
    
    INTEGER(HID_T) :: plist_id, space_id
    INTEGER(INT32) :: rank, error
    INTEGER(HSIZE_T), ALLOCATABLE :: maxdims(:)

    rank = UBOUND(chunk,1)

    ALLOCATE(maxdims(rank))
    maxdims = H5S_UNLIMITED_F
    
    CALL h5Pcreate_f(H5P_DATASET_CREATE_F, plist_id, error)
    CALL h5Pset_deflate_f(plist_id, gzip, error)
    CALL h5Pset_chunk_f(plist_id, rank, chunk, error)
    CALL h5Screate_simple_f(rank, chunk, space_id, error, maxdims)
    CALL h5Dcreate_f(location, TRIM(name), datatype, space_id, dataset, error, plist_id)
    CALL h5Sclose_f(space_id, error)
    CALL h5Pclose_f(plist_id, error)
    
    RETURN
  END SUBROUTINE h5LXset_extendible_dataset


  SUBROUTINE h5LXwrite_extendible_dataset_F64(dataset, var, rank, size, start, count)
    
    INTEGER(HID_T), INTENT(IN) :: dataset
    REAL(REAL64), DIMENSION(*), INTENT(IN) :: var
    INTEGER(INT32), INTENT(IN) :: rank
    INTEGER(HSIZE_T), DIMENSION(*), INTENT(IN) :: size
    INTEGER(HSIZE_T), DIMENSION(*), INTENT(IN) :: start
    INTEGER(HSIZE_T), DIMENSION(*), INTENT(IN) :: count
    
    INTEGER(HID_T) :: dataspace, space_mem
#ifdef MPIO    
    INTEGER(HID_T) :: xferlist
#endif    
    INTEGER(INT32) :: error

#ifdef MPIO    
    CALL h5Pcreate_f(H5P_DATASET_XFER_F, xferlist, error)
    CALL h5Pset_dxpl_mpio_f(xferlist, DATA_XFER_MODE, error) 
#endif
    CALL h5Dset_extent_f(dataset, size, error)
    CALL h5Screate_simple_f(rank, count, space_mem, error)
    CALL h5Dget_space_f(dataset, dataspace, error)
    CALL h5Sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
#ifdef MPIO    
    CALL h5Dwrite_f(dataset, H5T_IEEE_F64LE, var, count, error, space_mem, dataspace, xferlist)
    CALL h5Pclose_f(xferlist, error)
#else    
    CALL h5Dwrite_f(dataset, H5T_IEEE_F64LE, var, count, error, space_mem, dataspace)
#endif
    CALL h5Sclose_f(dataspace, error)
    CALL h5Sclose_f(space_mem, error)
    
    RETURN
  END SUBROUTINE h5LXwrite_extendible_dataset_F64



  SUBROUTINE h5LXwrite_extendible_dataset_F32(dataset, var, rank, size, start, count)
    
    INTEGER(HID_T), INTENT(IN) :: dataset
    REAL(REAL32), DIMENSION(*), INTENT(IN) :: var
    INTEGER(INT32), INTENT(IN) :: rank
    INTEGER(HSIZE_T), DIMENSION(*), INTENT(IN) :: size
    INTEGER(HSIZE_T), DIMENSION(*), INTENT(IN) :: start
    INTEGER(HSIZE_T), DIMENSION(*), INTENT(IN) :: count
    
    INTEGER(HID_T) :: dataspace, space_mem
#ifdef MPIO    
    INTEGER(HID_T) :: xferlist
#endif    
    INTEGER(INT32) :: error

#ifdef MPIO    
    CALL h5Pcreate_f(H5P_DATASET_XFER_F, xferlist, error)
    CALL h5Pset_dxpl_mpio_f(xferlist, DATA_XFER_MODE, error) 
#endif
    CALL h5Dset_extent_f(dataset, size, error)
    CALL h5Screate_simple_f(rank, count, space_mem, error)
    CALL h5Dget_space_f(dataset, dataspace, error)
    CALL h5Sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
#ifdef MPIO    
    CALL h5Dwrite_f(dataset, H5T_IEEE_F32LE, var, count, error, space_mem, dataspace, xferlist)
    CALL h5Pclose_f(xferlist, error)
#else   
    CALL h5Dwrite_f(dataset, H5T_IEEE_F32LE, var, count, error, space_mem, dataspace)
#endif
    CALL h5Sclose_f(dataspace, error)
    CALL h5Sclose_f(space_mem, error)
    
    RETURN
  END SUBROUTINE h5LXwrite_extendible_dataset_F32



  SUBROUTINE h5LXwrite_extendible_dataset_I32(dataset, var, rank, size, start, count)
    
    INTEGER(HID_T), INTENT(IN) :: dataset
    INTEGER(INT32), DIMENSION(*), INTENT(IN) :: var
    INTEGER(INT32), INTENT(IN) :: rank
    INTEGER(HSIZE_T), DIMENSION(*), INTENT(IN) :: size
    INTEGER(HSIZE_T), DIMENSION(*), INTENT(IN) :: start
    INTEGER(HSIZE_T), DIMENSION(*), INTENT(IN) :: count
    
    INTEGER(HID_T) :: space_file, space_mem
#ifdef MPIO    
    INTEGER(HID_T) :: xferlist
#endif    
    INTEGER(INT32) :: error

#ifdef MPIO    
    CALL h5Pcreate_f(H5P_DATASET_XFER_F, xferlist, error)
    CALL h5Pset_dxpl_mpio_f(xferlist, DATA_XFER_MODE, error) 
#endif
    CALL h5Dset_extent_f(dataset, size, error)
    CALL h5Dget_space_f(dataset, space_file, error)
    CALL h5Sselect_hyperslab_f(space_file, H5S_SELECT_SET_F, start, count, error)
    CALL h5Screate_simple_f(rank, count, space_mem, error)
#ifdef MPIO    
    CALL h5Dwrite_f(dataset, H5T_STD_I32LE, var, count, error, space_mem, space_file, xferlist)
    CALL h5Pclose_f(xferlist, error)
#else    
    CALL h5Dwrite_f(dataset, H5T_STD_I32LE, var, count, error, space_mem, space_file)
#endif
    CALL h5Sclose_f(space_file, error)
    CALL h5Sclose_f(space_mem, error)
    
    RETURN
  END SUBROUTINE h5LXwrite_extendible_dataset_I32



  SUBROUTINE h5LXset_dataset(location, name, datatype, rank, dims, dataset)
   
    INTEGER(HID_T), INTENT(IN) :: location
    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER(HID_T), INTENT(IN) :: datatype
    INTEGER(INT32), INTENT(IN) :: rank
    INTEGER(HSIZE_T), DIMENSION(*), INTENT(IN) :: dims
    INTEGER(HID_T), INTENT(OUT) :: dataset
    
    INTEGER(INT32) :: error
    INTEGER(HID_T) :: space_id
    
    CALL h5Screate_simple_f(rank, dims, space_id, error)
    CALL h5Dcreate_f(location, name, datatype, space_id, dataset, error)
    CALL h5Sclose_f(space_id, error)

    RETURN
  END SUBROUTINE h5LXset_dataset

  
  
  SUBROUTINE h5LXwrite_hyperslab_dataset_F64(dataset, var, rank, dims, start, count)
    
    INTEGER(HID_T), INTENT(IN) :: dataset
    REAL(REAL64), DIMENSION(*), INTENT(IN) :: var
    INTEGER(INT32), INTENT(IN) :: rank
    INTEGER(HSIZE_T), DIMENSION(*), INTENT(IN) :: dims
    INTEGER(HSIZE_T), DIMENSION(*), INTENT(IN) :: start
    INTEGER(HSIZE_T), DIMENSION(*), INTENT(IN) :: count
    
    INTEGER(HID_T) :: space_file, space_mem
#ifdef MPIO    
    INTEGER(HID_T) :: xferlist
#endif    
    INTEGER(INT32) :: error

#ifdef MPIO    
    CALL h5Pcreate_f(H5P_DATASET_XFER_F, xferlist, error)
    CALL h5Pset_dxpl_mpio_f(xferlist, DATA_XFER_MODE, error) 
#endif
    ! dataset用にfile dataspaceを確保する
    CALL h5Dget_space_f(dataset, space_file, error)
    ! 確保したfile dataspaceのうち、一部（subset）だけを選ぶ
    CALL h5Sselect_hyperslab_f(space_file, H5S_SELECT_SET_F, start, count, error)
    ! 値varを格納するmemory dataspaceを確保する
    CALL h5Screate_simple_f(rank, dims, space_mem, error)
#ifdef MPIO
    CALL h5Dwrite_f(dataset, H5T_IEEE_F64LE, var, dims, error, space_mem, space_file, xferlist)
    CALL h5Pclose_f(xferlist, error)
#else    
    CALL h5Dwrite_f(dataset, H5T_IEEE_F64LE, var, dims, error, space_mem, space_file)
#endif
    CALL h5Sclose_f(space_file, error)
    CALL h5Sclose_f(space_mem, error)
    
    RETURN
  END SUBROUTINE h5LXwrite_hyperslab_dataset_F64
  
  
  SUBROUTINE h5LXwrite_hyperslab_dataset_F32(dataset, var, rank, dims, start, count)
    
    INTEGER(HID_T), INTENT(IN) :: dataset
    REAL(REAL32), DIMENSION(*), INTENT(IN) :: var
    INTEGER(INT32), INTENT(IN) :: rank
    INTEGER(HSIZE_T), DIMENSION(*), INTENT(IN) :: dims
    INTEGER(HSIZE_T), DIMENSION(*), INTENT(IN) :: start
    INTEGER(HSIZE_T), DIMENSION(*), INTENT(IN) :: count
    
    INTEGER(HID_T) :: dataspace, space_mem
#ifdef MPIO    
    INTEGER(HID_T) :: xferlist
#endif    
    INTEGER(INT32) :: error

#ifdef MPIO    
    CALL h5Pcreate_f(H5P_DATASET_XFER_F, xferlist, error)
    CALL h5Pset_dxpl_mpio_f(xferlist, DATA_XFER_MODE, error) 
#endif
    CALL h5Screate_simple_f(rank, dims, space_mem, error)
    CALL h5Dget_space_f(dataset, dataspace, error)
    CALL h5Sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
#ifdef MPIO
    CALL h5Dwrite_f(dataset, H5T_IEEE_F32LE, var, dims, error, space_mem, dataspace, xferlist)
    CALL h5Pclose_f(xferlist, error)
#else    
    CALL h5Dwrite_f(dataset, H5T_IEEE_F32LE, var, dims, error, space_mem, dataspace)
#endif
    CALL h5Sclose_f(dataspace, error)
    CALL h5Sclose_f(space_mem, error)
    
    RETURN
  END SUBROUTINE h5LXwrite_hyperslab_dataset_F32
  
  
  SUBROUTINE h5LXwrite_hyperslab_dataset_I32(dataset, var, rank, dims, start, count)
    
    INTEGER(HID_T), INTENT(IN) :: dataset
    INTEGER(INT32), DIMENSION(*), INTENT(IN) :: var
    INTEGER(INT32), INTENT(IN) :: rank
    INTEGER(HSIZE_T), DIMENSION(*), INTENT(IN) :: dims
    INTEGER(HSIZE_T), DIMENSION(*), INTENT(IN) :: start
    INTEGER(HSIZE_T), DIMENSION(*), INTENT(IN) :: count
    
    INTEGER(HID_T) :: dataspace, space_mem
#ifdef MPIO    
    INTEGER(HID_T) :: xferlist
#endif    
    INTEGER(INT32) :: error

#ifdef MPIO    
    CALL h5Pcreate_f(H5P_DATASET_XFER_F, xferlist, error)
    CALL h5Pset_dxpl_mpio_f(xferlist, DATA_XFER_MODE, error) 
#endif
    CALL h5Screate_simple_f(rank, dims, space_mem, error)
    CALL h5Dget_space_f(dataset, dataspace, error)
    CALL h5Sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
#ifdef MPIO
    CALL h5Dwrite_f(dataset, H5T_STD_I32LE, var, dims, error, space_mem, dataspace, xferlist)
    CALL h5Pclose_f(xferlist, error)
#else    
    CALL h5Dwrite_f(dataset, H5T_STD_I32LE, var, dims, error, space_mem, dataspace)
#endif
    CALL h5Sclose_f(dataspace, error)
    CALL h5Sclose_f(space_mem, error)
    
    RETURN
  END SUBROUTINE h5LXwrite_hyperslab_dataset_I32
  
  
  SUBROUTINE h5LXread_hyperslab_dataset_F64(dataset, var, rank, dims, start, count)
    
    INTEGER(HID_T), INTENT(IN) :: dataset
    REAL(REAL64), DIMENSION(*), INTENT(OUT) :: var
    INTEGER(INT32) :: rank
    INTEGER(HSIZE_T), DIMENSION(*), INTENT(IN) :: dims
    INTEGER(HSIZE_T), DIMENSION(*), INTENT(IN) :: start
    INTEGER(HSIZE_T), DIMENSION(*), INTENT(IN) :: count

    INTEGER(HID_T) :: dataspace, space_mem
#ifdef MPIO    
    INTEGER(HID_T) :: xferlist
#endif    
    INTEGER(INT32) :: error

#ifdef MPIO    
    CALL h5pcreate_f(H5P_DATASET_XFER_F, xferlist, error)
    CALL h5Pset_dxpl_mpio_f(xferlist, DATA_XFER_MODE, error)
#endif    
    CALL h5Screate_simple_f(rank, dims, space_mem, error)
    CALL h5Dget_space_f(dataset, dataspace, error)
    CALL h5Sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
#ifdef MPIO
    CALL h5Dread_f(dataset, H5T_IEEE_F64LE, var, dims, error, space_mem, dataspace, xferlist)
    CALL h5Pclose_f(xferlist, error)
#else    
    CALL h5Dread_f(dataset, H5T_IEEE_F64LE, var, dims, error, space_mem, dataspace)
#endif    
    CALL h5Sclose_f(dataspace, error)
    CALL h5Sclose_f(space_mem, error)
    
  END SUBROUTINE h5LXread_hyperslab_dataset_F64
  

  
  SUBROUTINE h5LXread_hyperslab_dataset_F32(dataset, var, rank, dims, start, count)
    
    INTEGER(HID_T), INTENT(IN) :: dataset
    REAL(REAL32), DIMENSION(*), INTENT(OUT) :: var
    INTEGER(INT32) :: rank
    INTEGER(HSIZE_T), DIMENSION(*), INTENT(IN) :: dims
    INTEGER(HSIZE_T), DIMENSION(*), INTENT(IN) :: start
    INTEGER(HSIZE_T), DIMENSION(*), INTENT(IN) :: count
    
    INTEGER(HID_T) :: dataspace, space_mem
#ifdef MPIO
    INTEGER(HID_T) :: xferlist
#endif    
    INTEGER(INT32) :: error

#ifdef MPIO    
    CALL h5pcreate_f(H5P_DATASET_XFER_F, xferlist, error)
    CALL h5Pset_dxpl_mpio_f(xferlist, DATA_XFER_MODE, error)
#endif    
    CALL h5Screate_simple_f(rank, dims, space_mem, error)
    CALL h5Dget_space_f(dataset, dataspace, error)
    CALL h5Sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
#ifdef MPIO
    CALL h5Dread_f(dataset, H5T_IEEE_F32LE, var, dims, error, space_mem, dataspace, xferlist)
    CALL h5Pclose_f(xferlist, error)
#else    
    CALL h5Dread_f(dataset, H5T_IEEE_F32LE, var, dims, error, space_mem, dataspace)
#endif    
    CALL h5Sclose_f(dataspace, error)
    CALL h5Sclose_f(space_mem, error)
    
  END SUBROUTINE h5LXread_hyperslab_dataset_F32
  

  
  SUBROUTINE h5LXread_hyperslab_dataset_I32(dataset, var, rank, dims, start, count)
    
    INTEGER(HID_T), INTENT(IN) :: dataset
    INTEGER(INT32), DIMENSION(*), INTENT(OUT) :: var
    INTEGER(INT32) :: rank
    INTEGER(HSIZE_T), DIMENSION(*), INTENT(IN) :: dims
    INTEGER(HSIZE_T), DIMENSION(*), INTENT(IN) :: start
    INTEGER(HSIZE_T), DIMENSION(*), INTENT(IN) :: count

    INTEGER(HID_T) :: dataspace, space_mem
#ifdef MPIO    
    INTEGER(HID_T) :: xferlist
#endif    
    INTEGER(INT32) :: error

#ifdef MPIO    
    CALL h5pcreate_f(H5P_DATASET_XFER_F, xferlist, error)
    CALL h5Pset_dxpl_mpio_f(xferlist, DATA_XFER_MODE, error)
#endif    
    CALL h5Screate_simple_f(rank, dims, space_mem, error)
    CALL h5Dget_space_f(dataset, dataspace, error)
    CALL h5Sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
#ifdef MPIO
    CALL h5Dread_f(dataset, H5T_STD_I32LE, var, dims, error, space_mem, dataspace, xferlist)
    CALL h5Pclose_f(xferlist, error)
#else    
    CALL h5Dread_f(dataset, H5T_STD_I32LE, var, dims, error, space_mem, dataspace)
#endif    
    CALL h5Sclose_f(dataspace, error)
    CALL h5Sclose_f(space_mem, error)
    
  END SUBROUTINE h5LXread_hyperslab_dataset_I32
  

  
  SUBROUTINE h5LXread_hyperslab_dataset_I16(dataset, var, rank, dims, start, count)
    
    INTEGER(HID_T), INTENT(IN) :: dataset
    INTEGER(INT16), DIMENSION(*), INTENT(OUT) :: var
    INTEGER(INT32) :: rank
    INTEGER(HSIZE_T), DIMENSION(*), INTENT(IN) :: dims
    INTEGER(HSIZE_T), DIMENSION(*), INTENT(IN) :: start
    INTEGER(HSIZE_T), DIMENSION(*), INTENT(IN) :: count

    INTEGER(HID_T) :: dataspace, space_mem
#ifdef MPIO    
    INTEGER(HID_T) :: xferlist
#endif    
    INTEGER(INT32) :: error

#ifdef MPIO    
    CALL h5pcreate_f(H5P_DATASET_XFER_F, xferlist, error)
    CALL h5Pset_dxpl_mpio_f(xferlist, DATA_XFER_MODE, error)
#endif    
    CALL h5Screate_simple_f(rank, dims, space_mem, error)
    CALL h5Dget_space_f(dataset, dataspace, error)
    CALL h5Sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, start, count, error)
#ifdef MPIO
    CALL h5Dread_f(dataset, H5T_STD_I16LE, var, dims, error, space_mem, dataspace, xferlist)
    CALL h5Pclose_f(xferlist, error)
#else    
    CALL h5Dread_f(dataset, H5T_STD_I16LE, var, dims, error, space_mem, dataspace)
#endif    
    CALL h5Sclose_f(dataspace, error)
    CALL h5Sclose_f(space_mem, error)
    
  END SUBROUTINE h5LXread_hyperslab_dataset_I16
  

  
END MODULE h5LX_module
