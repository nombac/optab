
PROGRAM h5_create_Exomol2
  
  USE ISO_FORTRAN_ENV
  USE HDF5
  USE H5LT

  USE sort_module, ONLY : sorti

  IMPLICIT NONE

  REAL(8), PARAMETER :: m_ele = 9.10938d-28, clight = 2.99792d+10
  REAL(8), PARAMETER :: pi = 3.1415926535898d+00
  REAL(8), PARAMETER :: hbar = 1.0545726397789448d-27
  REAL(8), PARAMETER :: e2 = (1d0/137.035999084d0) * hbar * clight
  REAL(8), PARAMETER :: h = 6.6260755D-27

  INTEGER(INT64), PARAMETER :: chunk0 = 1024**2/8
  INTEGER(HID_T) :: space_unlimited, space_extend, space_written, space_slab, space_mem
  INTEGER(HID_T) :: prop_chunk
  INTEGER(HID_T) :: data_freq, data_lowe, data_gfva
  INTEGER(HID_T) :: dset_freq, dset_lowe, dset_gfva
  INTEGER(HID_T) :: file_id, descmp, st
  INTEGER(HID_T) :: grp_states, grp_trans, grp_species
  INTEGER(HID_T), ALLOCATABLE, DIMENSION(:) :: desc
  INTEGER(INT64) :: size_written, nprog, chunk_size

  CHARACTER*256 :: fn_states, fn_trans, dir, species
!  CHARACTER*256, PARAMETER :: filename_h5 = 'test.data.h5'
  CHARACTER*256 :: filename_h5, tmp_h5, test
  
  INTEGER(INT32), PARAMETER :: file_states = 10, descrans = 11, file_input = 12
  REAL(REAL64), ALLOCATABLE :: wn(:), min_finder(:)
  INTEGER(INT32), ALLOCATABLE :: gtot(:), m(:), read_flag(:)
  INTEGER(INT64) :: id, nc, n_rest, nc_max, nf_min(1), nn
  INTEGER(INT32) :: iostat
  INTEGER(INT32) :: id_u, id_l, gtot_tmp, ntimes
  REAL(REAL64) :: acoeff, wn_tmp
  REAL(REAL64), ALLOCATABLE, DIMENSION(:) :: freq, lowe, gfva, freq_tmp, lowe_tmp, gfva_tmp, freq1, lowe1, gfva1
  REAL(REAL64), ALLOCATABLE, DIMENSION(:,:) :: freq0, lowe0, gfva0
  INTEGER(INT64), ALLOCATABLE, DIMENSION(:) :: indx, nlines, chunk, offset, nlines_read, chunk1
  
  INTEGER(INT32) :: error
  INTEGER(INT64) :: l, nstates, nlines_all, nf, nlines_written, size
  CHARACTER*5 :: str_wns, str_wne
  INTEGER(INT32):: wnc, wn0, wn1, dwn, wrange

  dir = 'data/Exomol/'
!  dir = '/Users/shirose/Downloads/'
  
  CALL h5open_f(error)

  ! species loop
  OPEN(file_input, FILE='input', STATUS='OLD')
  DO
     READ(file_input, *, IOSTAT=iostat) species, wn0, wn1, dwn, wrange
     IF(iostat /= 0) EXIT
     
     dir = TRIM(dir)//TRIM(species)//'/'
     PRINT *, TRIM(dir)
     
     ! create HDF5 file
     filename_h5 = TRIM(dir)//TRIM(species)//'.h5'
     CALL h5Fcreate_f(TRIM(filename_h5), H5F_ACC_TRUNC_F, file_id, error)

     CALL h5Gcreate_f(file_id, TRIM(species)           , grp_species, error)
     CALL h5Gcreate_f(file_id, TRIM(species)//'/states', grp_states , error)
     CALL h5Gcreate_f(file_id, TRIM(species)//'/trans' , grp_trans  , error)

     !%%% states %%%
     !read states data
     fn_states = TRIM(dir)//TRIM(species)//'.states'
     OPEN(file_states, FILE=fn_states, STATUS='OLD')
     nstates = 0
     DO 
        READ(file_states, *, IOSTAT=iostat)
        IF(iostat /= 0) EXIT
        nstates = nstates + 1
     END DO
     CLOSE(file_states)
!     PRINT *, TRIM(fn_states), nstates

     ALLOCATE(gtot(nstates), wn(nstates))
     OPEN(file_states, file=TRIM(fn_states), status='OLD')
     DO l = 1, nstates
        READ(file_states, *) id, wn_tmp, gtot_tmp
        wn(id) = wn_tmp
        gtot(id) = gtot_tmp
     END DO
     CLOSE(file_states)

     !write states data
     CALL h5LTmake_dataset_f(grp_states, 'nstates', 1, (/1_HID_T/), H5T_STD_I64LE , nstates , error)
     CALL h5LTmake_dataset_f(grp_states, 'freq'   , 1, (/nstates/), H5T_IEEE_F64LE, wn      , error)
     CALL h5LTmake_dataset_f(grp_states, 'gtot'   , 1, (/nstates/), H5T_STD_I32LE , gtot    , error)

     !SORT
     nf = 0
     DO wnc = wn0, wn1, dwn
        ! READ TRANS DATA
        WRITE(str_wns,'(I5.5)') wnc
        WRITE(str_wne,'(I5.5)') wnc + dwn
        IF(wnc == wnc + dwn) THEN
           fn_trans = TRIM(species)//'.trans'
        ELSE
           fn_trans = TRIM(species)//'__'//TRIM(str_wns//'-'//str_wne)//'.trans'
        ENDIF
        fn_trans = TRIM(dir)//TRIM(fn_trans)
        
        PRINT *, TRIM(fn_trans)

        OPEN(descrans, FILE=TRIM(fn_trans), STATUS='OLD', IOSTAT=iostat)
        nlines = 0
        DO
           READ(descrans, *, IOSTAT=iostat)
           IF(iostat /= 0) EXIT
           nlines = nlines + 1
        END DO
        CLOSE(descrans)
        PRINT *, 'nlines =', nlines

        ALLOCATE(freq(nlines), lowe(nlines), gfva(nlines), indx(nlines), &
             freq_tmp(nlines), lowe_tmp(nlines), gfva_tmp(nlines))

        ! COMPUTE GF VALUE, ...
        OPEN(descrans, file=TRIM(fn_trans), status='OLD')
        PRINT *, 'begin reading ...'
        DO l = 1, nlines
           READ(descrans, *) id_u, id_l, acoeff
           freq_tmp(l) = wn(id_u) - wn(id_l)
           lowe_tmp(l) = wn(id_l)
           gfva_tmp(l) = gtot(id_u) * acoeff / (wn(id_u) - wn(id_l))**2
           indx(l) = l
        END DO
        PRINT *, 'done'
        CLOSE(descrans)
        PRINT *, 'begin sorting ...'
        CALL sorti(nlines, freq_tmp, indx)
        PRINT *, 'done'
        DO l = 1, nlines
           freq(l) = freq_tmp(indx(l))
           lowe(l) = lowe_tmp(indx(l))
           gfva(l) = gfva_tmp(indx(l))
        END DO
        lowe = lowe * h * clight
        gfva = gfva * (m_ele * clight / (8d0 * pi**2 * e2))

        tmp_h5 = TRIM(fn_trans)//'.h5'
        CALL h5Fcreate_f(TRIM(tmp_h5), H5F_ACC_TRUNC_F, descmp, error)
        CALL h5Screate_simple_f(1, (/nlines/), space_written, error)
        chunk_size = MIN(chunk0, nlines)
        CALL h5Pcreate_f(H5P_DATASET_CREATE_F, prop_chunk, error)
        CALL h5Pset_chunk_f(prop_chunk, 1, [chunk_size], error)
        
        CALL h5Dcreate_f(descmp, 'freq', H5T_IEEE_F64LE, space_written, data_freq, error, prop_chunk)
        CALL h5Dcreate_f(descmp, 'lowe', H5T_IEEE_F64LE, space_written, data_lowe, error, prop_chunk)
        CALL h5Dcreate_f(descmp, 'gfva', H5T_IEEE_F64LE, space_written, data_gfva, error, prop_chunk)
        CALL h5Dwrite_f(data_freq, H5T_IEEE_F64LE, freq, (/nlines/), error)
        CALL h5Dwrite_f(data_lowe, H5T_IEEE_F64LE, lowe, (/nlines/), error)
        CALL h5Dwrite_f(data_gfva, H5T_IEEE_F64LE, gfva, (/nlines/), error)

        CALL h5LTmake_dataset_f(descmp, 'chunk0', 1, (/1_HID_T/), H5T_STD_I64LE, chunk_size, error)
        CALL h5LTmake_dataset_f(descmp, 'nlines' , 1, (/1_HID_T/), H5T_STD_I64LE, nlines , error)
        
        CALL h5Dclose_f(data_freq, error)
        CALL h5Dclose_f(data_lowe, error)
        CALL h5Dclose_f(data_gfva, error)
        CALL h5Pclose_f(prop_chunk, error)
        CALL h5Sclose_f(space_written, error)
        CALL h5Fclose_f(descmp, error)
        DEALLOCATE(freq, lowe, gfva, indx, freq_tmp, lowe_tmp, gfva_tmp)
        nf = nf + 1
        
     END DO
     DEALLOCATE(gtot, wn)
     
     PRINT *, 'nf=', nf
     
     ALLOCATE(desc(nf), nlines(nf), chunk(nf), offset(nf), chunk1(nf))
     ALLOCATE(freq0(chunk0,nf), lowe0(chunk0,nf), gfva0(chunk0,nf), min_finder(nf))
     ALLOCATE(m(nf), read_flag(nf), nlines_read(nf))
     ALLOCATE(freq1(chunk0), lowe1(chunk0), gfva1(chunk0))

     ! MERGED FILE
     CALL h5Pcreate_f(H5P_DATASET_CREATE_F, prop_chunk, error)
     CALL h5Pset_chunk_f(prop_chunk, 1, (/chunk0/), error)
     CALL h5Screate_simple_f(1, (/chunk0/), space_unlimited, error, (/H5S_UNLIMITED_F/))
     CALL h5Dcreate_f(grp_trans, 'freq', H5T_IEEE_F64LE, space_unlimited, dset_freq, error, prop_chunk)
     CALL h5Dcreate_f(grp_trans, 'lowe', H5T_IEEE_F64LE, space_unlimited, dset_lowe, error, prop_chunk)
     CALL h5Dcreate_f(grp_trans, 'gfva', H5T_IEEE_F64LE, space_unlimited, dset_gfva, error, prop_chunk)
     CALL h5Sclose_f(space_unlimited, error)
     CALL h5Pclose_f(prop_chunk, error)
     

     ! MERGE

     read_flag(1:nf) = 0
     nlines_read(1:nf) = 0
     offset(1:nf) = 0
     freq0(1:chunk0,1:nf) = HUGE(0d0)
     ntimes = 0

     ! OPEN FILES AND READ INITIAL CHUNK
     nf = 0
     DO wnc = wn0, wn1, dwn
        nf = nf + 1
       
        WRITE(tmp_h5,'(A,"__",I5.5,"-",I5.5,".trans.h5")') TRIM(dir)//TRIM(species), wnc, wnc+dwn
        CALL h5Fopen_f(TRIM(tmp_h5), H5F_ACC_RDONLY_F, desc(nf), error)
        CALL h5ltread_dataset_f(desc(nf), 'nlines', H5T_STD_I64LE, nlines(nf), [1_HID_T], error)
        CALL h5ltread_dataset_f(desc(nf), 'n_chunk', H5T_STD_I64LE, chunk(nf), [1_HID_T], error)

        chunk1(nf) = MIN(chunk(nf), nlines(nf) - offset(nf))
        CALL h5Dopen_f(desc(nf), 'freq', data_freq, error)
        CALL h5Dopen_f(desc(nf), 'lowe', data_lowe, error)
        CALL h5Dopen_f(desc(nf), 'gfva', data_gfva, error)
        CALL h5Dget_space_f(data_freq, space_slab, error)
        CALL h5Dget_space_f(data_gfva, space_slab, error)
        CALL h5Dget_space_f(data_lowe, space_slab, error)
        CALL h5Screate_simple_f(1, [chunk1(nf)], space_mem, error)
        CALL h5Sselect_hyperslab_f(space_slab, H5S_SELECT_SET_F, [offset(nf)], [1_HID_T], error, [1_HID_T], [chunk1(nf)])
        CALL H5DREAD_F(data_freq, H5T_IEEE_F64LE, freq0(1:chunk1(nf),nf), [chunk1(nf)], error, space_mem, space_slab)
        CALL H5DREAD_F(data_lowe, H5T_IEEE_F64LE, lowe0(1:chunk1(nf),nf), [chunk1(nf)], error, space_mem, space_slab)
        CALL H5DREAD_F(data_gfva, H5T_IEEE_F64LE, gfva0(1:chunk1(nf),nf), [chunk1(nf)], error, space_mem, space_slab)
        CALL h5Sclose_f(space_mem, error)
        CALL h5Dclose_f(data_freq, error)
        CALL h5Dclose_f(data_lowe, error)
        CALL h5Dclose_f(data_gfva, error)
        offset(nf) = offset(nf) + chunk1(nf)
        
        nlines_read(nf) = nlines_read(nf) + chunk1(nf)
        IF(nlines_read(nf) == nlines(nf)) THEN
           read_flag(nf) = 1
        END IF

        m(nf) = 1
        min_finder(nf) = freq0(m(nf),nf)
     END DO

     nlines_all = SUM(nlines)
     PRINT *, 'nlines_all = ', nlines_all

     m(1:nf) = 1

     nprog = 1

     nlines_written = 0
     size = 0

     st = 0
     DO l = 1, nlines_all
        ! FIND FILE NUMBER HAVING MINIMUM VALUE
        nf_min = MINLOC(min_finder)
        nf = nf_min(1)

        ! STORE THE MINIMUM VALUE IN THE RECORD ARRAY
        st = st + 1
        freq1(st) = freq0(m(nf),nf)
        lowe1(st) = lowe0(m(nf),nf)
        gfva1(st) = gfva0(m(nf),nf)

        ! IF THE RECORD ARRAY IS FULL, WRITE IT OUT
        IF(st == chunk0 .OR. l == nlines_all) THEN
           CALL h5screate_simple_f(1, [st], space_written, error)
           size = size + st
           CALL h5Dset_extent_f(dset_freq, [size], error)
           CALL h5Dset_extent_f(dset_lowe, [size], error)
           CALL h5Dset_extent_f(dset_gfva, [size], error)
           CALL h5Dget_space_f(dset_freq, space_extend, error)
           CALL h5Dget_space_f(dset_lowe, space_extend, error)
           CALL h5Dget_space_f(dset_gfva, space_extend, error)
           CALL h5Sselect_hyperslab_f(space_extend, H5S_SELECT_SET_F, [nlines_written], (/1_HID_T/), error, BLOCK=[st])
           ! write the line data
           CALL h5Dwrite_f(dset_freq, H5T_IEEE_F64LE, freq1(1:st), [st], error, space_written, space_extend)
           CALL h5Dwrite_f(dset_lowe, H5T_IEEE_F64LE, lowe1(1:st), [st], error, space_written, space_extend)
           CALL h5Dwrite_f(dset_gfva, H5T_IEEE_F64LE, gfva1(1:st), [st], error, space_written, space_extend)
           ! close extended space
           CALL h5Sclose_f(space_extend, error)
           ! close dump space
           CALL h5Sclose_f(space_written, error)
           nlines_written = nlines_written + st
           st = 0
        END IF

        ! PRINT PROGRESS
        IF(REAL(l)/REAL(nlines_all)*100. > REAL(nprog)) THEN
           PRINT *, freq1(st), l, nf, INT(REAL(l)/REAL(nlines_all)*100.), '% done'
           nprog = nprog + 1
        END IF

        ! PROCEED TO THE NEXT ELEMENT
        m(nf) = m(nf) + 1

        IF(m(nf) <= chunk1(nf)) THEN
           min_finder(nf) = freq0(m(nf),nf)
        ELSE
           IF(read_flag(nf) == 1) THEN
              min_finder(nf) = HUGE(0d0)
              PRINT *, 'search completed: nf = ', nf
           ELSE
              chunk1(nf) = MIN(chunk(nf), nlines(nf) - offset(nf))
              CALL h5Dopen_f(desc(nf), 'freq', data_freq, error)
              CALL h5Dopen_f(desc(nf), 'lowe', data_lowe, error)
              CALL h5Dopen_f(desc(nf), 'gfva', data_gfva, error)
              CALL h5Dget_space_f(data_freq, space_slab, error)
              CALL h5Dget_space_f(data_gfva, space_slab, error)
              CALL h5Dget_space_f(data_lowe, space_slab, error)
              CALL h5Screate_simple_f(1, [chunk1(nf)], space_mem, error)
              CALL h5Sselect_hyperslab_f(space_slab, H5S_SELECT_SET_F, [offset(nf)], [1_HID_T], error, [1_HID_T], [chunk1(nf)])
              CALL H5DREAD_F(data_freq, H5T_IEEE_F64LE, freq0(1:chunk1(nf),nf), [chunk1(nf)], error, space_mem, space_slab)
              CALL H5DREAD_F(data_lowe, H5T_IEEE_F64LE, lowe0(1:chunk1(nf),nf), [chunk1(nf)], error, space_mem, space_slab)
              CALL H5DREAD_F(data_gfva, H5T_IEEE_F64LE, gfva0(1:chunk1(nf),nf), [chunk1(nf)], error, space_mem, space_slab)
              CALL h5Sclose_f(space_mem, error)
              CALL h5Dclose_f(data_freq, error)
              CALL h5Dclose_f(data_lowe, error)
              CALL h5Dclose_f(data_gfva, error)
              offset(nf) = offset(nf) + chunk1(nf)
              
              nlines_read(nf) = nlines_read(nf) + chunk1(nf)
              IF(nlines_read(nf) == nlines(nf)) THEN
                 read_flag(nf) = 1
              END IF

              m(nf) = 1
              min_finder(nf) = freq0(m(nf),nf)
           END IF
        END IF

     END DO


     PRINT *, 'total number of lines read', SUM(nlines_read), '/', nlines_all
     PRINT *, 'total number of lines written', nlines_written, '/', nlines_all

     nf = 0
     DO wnc = wn0, wn1, dwn
        nf = nf + 1
        PRINT *, 'min_finder(nf)=', nf, min_finder(nf)
        CALL h5Fclose_f(desc(nf), error)
     END DO
     
     CALL h5LTmake_dataset_f(grp_trans, 'nlines' , 1, (/1_HID_T/), H5T_STD_I64LE, nlines_all, error)
     CALL h5LTmake_dataset_f(grp_trans, 'n_chunk', 1, (/1_HID_T/), H5T_STD_I64LE, chunk0    , error)

     CALL h5Gclose_f(grp_trans  , error)
     CALL h5Gclose_f(grp_states , error)
     CALL h5Gclose_f(grp_species, error)
     CALL h5Fclose_f(file_id, error)

     DEALLOCATE(desc, nlines, chunk)
  END DO
  
  CALL h5close_f(error)

  STOP

END PROGRAM h5_create_Exomol2
