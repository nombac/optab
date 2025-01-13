
#define DIRH5 '../h5/'
#define DIRNIST '../NIST/'

PROGRAM convert_lines_h5
  
  USE ISO_FORTRAN_ENV
  IMPLICIT NONE
  
  CHARACTER*256, ALlOCATABLE :: fname(:)
  CHARACTER*256 :: linebuf
  INTEGER(INT32) :: nfiles, n, iostat
  
  ! コンバートするファイル名を読み込み
  OPEN(1, FILE='list_convert.txt', STATUS='OLD')
  nfiles = 0
  DO
     READ(1,*,IOSTAT=iostat)
     IF(iostat /= 0) EXIT
     nfiles = nfiles + 1
  END DO
  REWIND(1)
  ALLOCATE(fname(nfiles))
  DO n = 1, nfiles
     READ(1,FMT='(A)') linebuf
     fname(n) = TRIM(linebuf)
  END DO
  CLOSE(1)

  SELECT CASE (TRIM(fname(1)(INDEX(fname(1),'.',.TRUE.):)))
  CASE ('.trans')
     CALL convert_exomol_h5(fname, gzip=1)
  CASE ('.par')
     CALL convert_hitran_h5(fname, gzip=1)
  CASE ('.dat')
     CALL convert_kurucz_h5(fname(1), gzip=1)
  CASE DEFAULT
     PRINT *, '*** ERROR *** : could not figure out the format'
  END SELECT

  DEALLOCATE(fname)

  STOP

CONTAINS

  
  SUBROUTINE convert_hitran_h5(filename, gzip)
    USE HDF5
    USE H5LT
    USE H5LX_module, ONLY : h5LXset_extendible_dataset, h5LXwrite_extendible_dataset_F64, h5LXwrite_extendible_dataset_F32
    USE code_module, ONLY : AIR
    USE const_module, ONLY : m_ele, pi, e2, clight, atm_to_bar

    CHARACTER(*), INTENT(IN), DIMENSION(:) :: filename
    INTEGER(INT32), INTENT(IN) :: gzip

    INTEGER(HID_T) :: data_wnum, data_lowe, data_gfva
    INTEGER(HID_T) :: data_gm_1, data_nx_1, data_gm_2, data_nx_2
    INTEGER(HID_T) :: grp_trans, grp_prop, file_h5
    INTEGER(HSIZE_T), DIMENSION(1) :: size, offset, block_size, dims
    INTEGER(HSIZE_T), PARAMETER :: chunk0 = 32768
    
    REAL(REAL64), ALLOCATABLE, DIMENSION(:) :: wnum, lowe
    REAL(REAL32), ALLOCATABLE, DIMENSION(:) :: gm_1, gm_2, nx_1, nx_2, gfva

    REAL(REAL64) :: mass, frac
    INTEGER(INT64) :: l, nlines, n_rest, nc, nc_max
    INTEGER(INT32) :: n, ntemp
    INTEGER(INT32) :: id
    INTEGER(INT32) :: error
    CHARACTER(LEN=256) :: prefix
    CHARACTER(LEN=72) :: output, command
    
    ! .prop
    REAL(REAL64), ALLOCATABLE :: temp(:), pf(:)

    ALLOCATE(wnum(chunk0), lowe(chunk0), gm_1(chunk0), gm_2(chunk0), nx_1(chunk0), nx_2(chunk0), gfva(chunk0))

    CALL h5open_f(error)

    DO n = 1, UBOUND(filename,1)

       WRITE(*,FMT='(I3.3,A1,I3.3,X,A36)', ADVANCE='NO') n, '/', UBOUND(filename,1), trim(filename(n))

       ! 変換先のHDF5ファイルを開く
       prefix = filename(n)(INDEX(filename(n),'/')+1:INDEX(filename(n),'.',.TRUE.)-1)
       output = DIRH5//TRIM(prefix)//'.h5'
       CALL h5Fcreate_f(TRIM(output), H5F_ACC_TRUNC_F, file_h5, error)
     
       !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       ! 同位体置換体ごとの質量・存在比・分配関数
       CALL get_hitran_prop(filename(n), id, mass, frac, temp, pf, ntemp)
       CALL h5Gcreate_f(file_h5, 'prop', grp_prop, error)
       dims = [1]
       CALL h5LTmake_dataset_f(grp_prop, 'id'  , 1, dims, H5T_STD_I32LE , id   , error)
       CALL h5LTmake_dataset_f(grp_prop, 'mass', 1, dims, H5T_IEEE_F64LE, mass , error)
       CALL h5LTmake_dataset_f(grp_prop, 'frac', 1, dims, H5T_IEEE_F64LE, frac , error)
       CALL h5LTmake_dataset_f(grp_prop, 'ntemp',1, dims, H5T_STD_I32LE , ntemp, error)
       dims = [ntemp]
       CALL h5LTmake_dataset_f(grp_prop, 'temp', 1, dims, H5T_IEEE_F64LE, temp , error)
       CALL h5LTmake_dataset_f(grp_prop, 'pf'  , 1, dims, H5T_IEEE_F64LE, pf   , error)
       CALL h5Gclose_f(grp_prop, error)
       DEALLOCATE(temp,pf)

       ! ラインデータ書き込みデータセット
       CALL h5Gcreate_f(file_h5, 'trans', grp_trans, error)
       dims = [1]
       CALL h5LTmake_dataset_f(grp_trans, 'id_1', 1, dims, H5T_STD_I32LE, AIR, error)
       CALL h5LTmake_dataset_f(grp_trans, 'id_2', 1, dims, H5T_STD_I32LE, id , error)
       
       CALL h5LXset_extendible_dataset(data_wnum, grp_trans, 'wnum', H5T_IEEE_F64LE, [chunk0], gzip)
       CALL h5LXset_extendible_dataset(data_lowe, grp_trans, 'lowe', H5T_IEEE_F64LE, [chunk0], gzip)
       CALL h5LXset_extendible_dataset(data_gfva, grp_trans, 'gfva', H5T_IEEE_F32LE, [chunk0], gzip)
       CALL h5LXset_extendible_dataset(data_gm_1, grp_trans, 'gm_1', H5T_IEEE_F32LE, [chunk0], gzip)
       CALL h5LXset_extendible_dataset(data_nx_1, grp_trans, 'nx_1', H5T_IEEE_F32LE, [chunk0], gzip)
       CALL h5LXset_extendible_dataset(data_gm_2, grp_trans, 'gm_2', H5T_IEEE_F32LE, [chunk0], gzip)
       CALL h5LXset_extendible_dataset(data_nx_2, grp_trans, 'nx_2', H5T_IEEE_F32LE, [chunk0], gzip)

       CALL count_lines(filename(n), nlines)
       
       ! チャンクサイズで書き込み
       nc_max = nlines / chunk0
       n_rest = nlines - nc_max * chunk0

       block_size = [chunk0]
       offset = [0]
       size = [0]
       
       OPEN(2, FILE=filename(n), STATUS='OLD')
       
       DO nc = 1, nc_max + 1
          IF(nc == nc_max + 1) THEN
             block_size = [n_rest]
          END IF
          
          !***********************************************************
          ! ラインデータ読み込み
          DO l = 1, block_size(1)
             READ(2,*) wnum(l), lowe(l), gfva(l), gm_1(l), gm_2(l), nx_1(l)
             nx_2(l) = nx_1(l)
             gm_1(l) = gm_1(l) / REAL(atm_to_bar)
             gm_2(l) = gm_2(l) / REAL(atm_to_bar)
          END DO
          !***********************************************************
          
          size = size + block_size

          CALL h5LXwrite_extendible_dataset_F64(data_wnum, wnum, 1, size, offset, block_size) 
          CALL h5LXwrite_extendible_dataset_F64(data_lowe, lowe, 1, size, offset, block_size)
          CALL h5LXwrite_extendible_dataset_F32(data_gfva, gfva, 1, size, offset, block_size)
          CALL h5LXwrite_extendible_dataset_F32(data_gm_1, gm_1, 1, size, offset, block_size)
          CALL h5LXwrite_extendible_dataset_F32(data_nx_1, nx_1, 1, size, offset, block_size)
          CALL h5LXwrite_extendible_dataset_F32(data_gm_2, gm_2, 1, size, offset, block_size)
          CALL h5LXwrite_extendible_dataset_F32(data_nx_2, nx_2, 1, size, offset, block_size)
         
          offset = offset + block_size
       END DO
       
       CLOSE(2)
       
       ! ラインの本数とチャンクサイズを記録
       dims = [1]
       CALL h5LTmake_dataset_f(grp_trans, 'chunk0', 1, dims, H5T_STD_I32LE, chunk0, error)
       CALL h5LTmake_dataset_f(grp_trans, 'nlines', 1, dims, H5T_STD_I64LE, size(1)   , error)

       CALL h5Dclose_f(data_wnum, error)
       CALL h5Dclose_f(data_lowe, error)
       CALL h5Dclose_f(data_gfva, error)
       CALL h5Dclose_f(data_gm_1, error)
       CALL h5Dclose_f(data_nx_1, error)
       CALL h5Dclose_f(data_gm_2, error)
       CALL h5Dclose_f(data_nx_2, error)
       CALL h5Gclose_f(grp_trans, error)
       CALL h5Fclose_f(file_h5, error)
       
       WRITE(*,FMT='(A2,I13,A8)', ADVANCE='NO') ': ', nlines, ' lines: '
       WRITE(*,FMT='(A5)',ADVANCE='NO') '  => '
       WRITE(command,"('ls -sh ',A)") TRIM(output)
       CALL SYSTEM(TRIM(command))
    END DO
    
    !***********************************************************
    
    CALL h5close_f(error)
    
    DEALLOCATE(wnum, lowe, gm_1, gm_2, nx_1, nx_2, gfva)

    RETURN    
  END SUBROUTINE convert_hitran_h5



  SUBROUTINE convert_exomol_h5(filename, gzip)
    USE HDF5
    USE H5LT
    USE H5LX_module, ONLY : h5LXset_extendible_dataset, h5LXwrite_extendible_dataset_F64, h5LXwrite_extendible_dataset_F32
    USE const_module, ONLY : m_ele, pi, e2, clight
    USE nist_module, ONLY : get_nist_isofrac

    CHARACTER(*), INTENT(IN), DIMENSION(:) :: filename
    INTEGER(INT32), INTENT(IN) :: gzip

    INTEGER(HID_T) :: data_wnum, data_lowe, data_gfva
    INTEGER(HID_T) :: data_gm_1, data_nx_1, data_gm_2, data_nx_2
    INTEGER(HID_T) :: grp_trans, grp_prop, file_h5
    INTEGER(HSIZE_T), DIMENSION(1) :: size, offset, block_size, dims
    INTEGER(HSIZE_T), PARAMETER :: chunk0 = 32768

    REAL(REAL64), ALLOCATABLE, DIMENSION(:) :: wnum, lowe
    REAL(REAL32), ALLOCATABLE, DIMENSION(:) :: gm_1, gm_2, nx_1, nx_2, gfva

    INTEGER(INT64) :: l, nlines, nc_max, nc
    INTEGER(INT64) :: nlines_total, nlines_def
    INTEGER(INT64) :: id_u, id_l
    INTEGER(INT32) :: n, ntemp
    INTEGER(INT32) :: id, jt_1, jt_2
    INTEGER(INT32) :: error

    REAL(REAL64) :: mass, frac
    REAL(REAL32) :: acoeff
    
    CHARACTER :: prefix*256, command*72, output*72
    ! .states
    INTEGER(INT32), ALLOCATABLE, DIMENSION(:) :: gtot
    REAL(REAL32), ALLOCATABLE, DIMENSION(:) :: jtot
    REAL(REAL64), ALLOCATABLE, DIMENSION(:) :: ene, temp, pf
    ! .broad
    INTEGER(INT32) :: id_10, id_20
    REAL(REAL32), ALLOCATABLE, DIMENSION(:) :: jt_10, jt_20
    REAL(REAL32), ALLOCATABLE, DIMENSION(:) :: gm_10, gm_20, nx_10, nx_20

    ALLOCATE(wnum(chunk0), lowe(chunk0), gm_1(chunk0), gm_2(chunk0), nx_1(chunk0), nx_2(chunk0), gfva(chunk0))
    
    ! 同位体置換体比
!    CALL get_nist_isofrac(frac)
    !***********************************************************
    CALL get_exomol_prefix(filename(1), prefix, id)
    CALL get_exomol_isofrac(prefix, frac)
    ! レベルデータ読み込み
    CALL read_exomol_states(prefix, ene, gtot, jtot)
    ! 定義ファイル読み込み、質量とラインの本数
    CALL read_exomol_def(prefix, mass, nlines_def)
    ! 分配関数読み込み
    CALL read_exomol_pf(prefix, temp, pf)
    ntemp = UBOUND(temp,1)
    ! 圧力広がり読み込み
    CALL read_exomol_broad(prefix, id_10, gm_10, nx_10, jt_10, id_20, gm_20, nx_20, jt_20)
    !***********************************************************
    
    CALL h5open_f(error)

    ! 変換先のHDF5ファイルを開く
    output = DIRH5//TRIM(prefix(INDEX(prefix,'/')+1:))//'.h5'
    CALL h5Fcreate_f(TRIM(output), H5F_ACC_TRUNC_F, file_h5, error)

    ! 同位体置換体ごとの質量・存在比・分配関数
    CALL h5Gcreate_f(file_h5, 'prop', grp_prop, error)
    dims = [1]
    CALL h5LTmake_dataset_f(grp_prop, 'id'  , 1, dims, H5T_STD_I32LE , id   , error)
    CALL h5LTmake_dataset_f(grp_prop, 'mass', 1, dims, H5T_IEEE_F64LE, mass , error)
    CALL h5LTmake_dataset_f(grp_prop, 'frac', 1, dims, H5T_IEEE_F64LE, frac , error)
    CALL h5LTmake_dataset_f(grp_prop, 'ntemp',1, dims, H5T_STD_I32LE , ntemp, error)
    dims = [ntemp]
    CALL h5LTmake_dataset_f(grp_prop, 'temp', 1, dims, H5T_IEEE_F64LE, temp , error)
    CALL h5LTmake_dataset_f(grp_prop, 'pf'  , 1, dims, H5T_IEEE_F64LE, pf   , error)
    CALL h5Gclose_f(grp_prop, error)
       
    ! ラインデータ書き込みデータセット
    CALL h5Gcreate_f(file_h5, 'trans', grp_trans, error)
    dims = [1]
    CALL h5LTmake_dataset_f(grp_trans, 'id_1', 1, dims, H5T_STD_I32LE, id_10, error)
    CALL h5LTmake_dataset_f(grp_trans, 'id_2', 1, dims, H5T_STD_I32LE, id_20, error)

    CALL h5LXset_extendible_dataset(data_wnum, grp_trans, 'wnum', H5T_IEEE_F64LE, [chunk0], gzip)
    CALL h5LXset_extendible_dataset(data_lowe, grp_trans, 'lowe', H5T_IEEE_F64LE, [chunk0], gzip)
    CALL h5LXset_extendible_dataset(data_gfva, grp_trans, 'gfva', H5T_IEEE_F32LE, [chunk0], gzip)
    CALL h5LXset_extendible_dataset(data_gm_1, grp_trans, 'gm_1', H5T_IEEE_F32LE, [chunk0], gzip)
    CALL h5LXset_extendible_dataset(data_nx_1, grp_trans, 'nx_1', H5T_IEEE_F32LE, [chunk0], gzip)
    CALL h5LXset_extendible_dataset(data_gm_2, grp_trans, 'gm_2', H5T_IEEE_F32LE, [chunk0], gzip)
    CALL h5LXset_extendible_dataset(data_nx_2, grp_trans, 'nx_2', H5T_IEEE_F32LE, [chunk0], gzip)

    offset = [0]
    size = [0]
    
    nlines_total = 0

    ! 波数ごとに複数のファイルに分かれている場合
    DO n = 1, UBOUND(filename,1)
       write(*,fmt='(a,i3.3,a,i3.3,a)', advance='no') ' ', n, '/', UBOUND(filename,1), ': convert '//trim(filename(n))

       CALL count_lines(filename(n), nlines)
       
       nlines_total = nlines_total + nlines
       ! チャンクサイズで書き込み
       nc_max = nlines / chunk0

       OPEN(2, FILE=filename(n), STATUS='OLD')
       DO nc = 1, nc_max + 1
!          PRINT *, 'nc = ', nc, ' / ', nc_max+1
          IF(nc == nc_max + 1) THEN
             block_size = [nlines - nc_max * chunk0]
          ELSE
             block_size = [chunk0]
          END IF

          !***********************************************************
          ! ラインデータ読み込み
          DO l = 1, block_size(1)
             READ(2, *) id_u, id_l, acoeff
             wnum(l) = ene(id_u) - ene(id_l)
             lowe(l) = ene(id_l)
             ! Einstein A => gf
             gfva(l) = REAL(gtot(id_u) * acoeff / (ene(id_u) - ene(id_l))**2 &
                  * (m_ele * clight / (8d0 * pi**2 * e2)))
             
             jt_1 = FINDLOC(jt_10, MAX(MIN(jtot(id_l), MAXVAL(jt_10)), MINVAL(jt_10)), 1)
             gm_1(l) = gm_10(jt_1)
             nx_1(l) = nx_10(jt_1)

             jt_2 = FINDLOC(jt_20, MAX(MIN(jtot(id_l), MAXVAL(jt_20)), MINVAL(jt_20)), 1)
             gm_2(l) = gm_20(jt_2)
             nx_2(l) = nx_20(jt_2)

             IF(wnum(l) <= 0d0) gfva(l) = 0d0
          END DO
          !***********************************************************
          
          size = size + block_size

          CALL h5LXwrite_extendible_dataset_F64(data_wnum, wnum, 1, size, offset, block_size) 
          CALL h5LXwrite_extendible_dataset_F64(data_lowe, lowe, 1, size, offset, block_size)
          CALL h5LXwrite_extendible_dataset_F32(data_gfva, gfva, 1, size, offset, block_size)
          CALL h5LXwrite_extendible_dataset_F32(data_gm_1, gm_1, 1, size, offset, block_size)
          CALL h5LXwrite_extendible_dataset_F32(data_nx_1, nx_1, 1, size, offset, block_size)
          CALL h5LXwrite_extendible_dataset_F32(data_gm_2, gm_2, 1, size, offset, block_size)
          CALL h5LXwrite_extendible_dataset_F32(data_nx_2, nx_2, 1, size, offset, block_size)
         
          offset = offset + block_size
       END DO
       CLOSE(2)
       
       print '(a,i12,a)', ': ', nlines, ' lines done'
       IF(size(1) /= nlines_total) THEN
          print '(a, 2i10)', '*** error *** size /= nlines_total: ', size(1), nlines_total
          STOP
       END IF
    END DO
   
    DEALLOCATE(gtot, jtot, ene, temp, pf, jt_10, jt_20, gm_10, gm_20, nx_10, nx_20)
       
    ! ラインの本数とチャンクサイズを記録
    dims = [1]
    CALL h5LTmake_dataset_f(grp_trans, 'chunk0', 1, dims, H5T_STD_I32LE, chunk0 , error)
    CALL h5LTmake_dataset_f(grp_trans, 'nlines', 1, dims, H5T_STD_I64LE, size(1), error)

    !***********************************************************
    IF(nlines_total /= nlines_def) THEN
       print '(a,i12,a,i12)', '+++ warning: nlines_total = ', nlines_total, ' does not match nlines in .def: ', nlines_def
    ELSE
       print '(a,i12,a,i12)', ' nlines_total = ', nlines_total, ' matches nlines in .def: ', nlines_def
    END IF
    !***********************************************************
    DEALLOCATE(wnum, lowe, gm_1, gm_2, nx_1, nx_2, gfva)
    
    CALL h5Dclose_f(data_wnum, error)
    CALL h5Dclose_f(data_lowe, error)
    CALL h5Dclose_f(data_gfva, error)
    CALL h5Dclose_f(data_gm_1, error)
    CALL h5Dclose_f(data_nx_1, error)
    CALL h5Dclose_f(data_gm_2, error)
    CALL h5Dclose_f(data_nx_2, error)
    CALL h5Gclose_f(grp_trans, error)
    CALL h5Fclose_f(file_h5, error)
    CALL h5close_f(error)
    
    WRITE(*,FMT='(A13)',ADVANCE='NO') 'output file: '
    WRITE(command,"('ls -sh ',A)") TRIM(output)
    CALL SYSTEM(TRIM(command))

    RETURN

  END SUBROUTINE convert_exomol_h5


  SUBROUTINE convert_kurucz_h5(filename, gzip)
    USE HDF5
    USE H5LT
    USE H5LX_module, ONLY : h5LXset_extendible_dataset, &
         h5LXwrite_extendible_dataset_F64, h5LXwrite_extendible_dataset_F32, h5LXwrite_extendible_dataset_I32
    USE code_module, ONLY : HI, ELECTRON
    USE const_module, ONLY : ene_hyd, h, c
    USE nist_module, ONLY : get_nist_mass_comp, get_nist_eneion
    IMPLICIT NONE
    
    CHARACTER(*), INTENT(IN) :: filename
    INTEGER(INT32), INTENT(IN) :: gzip

    INTEGER(HID_T) :: data_wnum, data_lowe, data_gfva, data_code
    INTEGER(HID_T) :: data_gm_0, data_gm_1, data_gm_2
    INTEGER(HID_T) :: grp_trans, file_h5
    INTEGER(HSIZE_T), DIMENSION(1) :: size, offset, block_size, dims

    INTEGER(HSIZE_T), PARAMETER :: chunk0 = 32768
    REAL(REAL64), ALLOCATABLE, DIMENSION(:) :: wnum, lowe
    REAL(REAL32), ALLOCATABLE, DIMENSION(:) :: gm_0, gm_1, gm_2, gfva
    INTEGER(INT32), ALLOCATABLE, DIMENSION(:) :: code

    INTEGER(INT64) :: l, nlines, lcount, nc, nc_max, nlines_total
    INTEGER(INT32) :: ni, na, indx
    
    INTEGER(INT32) :: error
    CHARACTER :: lowlabel*10, highlabel*10, cdum*256
    INTEGER(INT32) :: isot, EJE_HF_shift, EJO_HF_shift, idum, isot2
    REAL(REAL32) :: lambda, loggf, code0, lowj, gamr, gam4, gam6, highj, lg_HF_frc, lg_frc, frc_nist
    REAL(REAL64) :: highe0, lowe0

    INTEGER(INT32), PARAMETER :: na_uran = 92, mn_uran = 243
    CHARACTER :: natm*3, nion*2, output*72, command*72
    REAL(REAL64), ALLOCATABLE :: comp_iso(:,:), mass(:), eneion(:)
    CHARACTER :: prefix*256

    ALLOCATE(wnum(chunk0), lowe(chunk0), gm_0(chunk0), gm_1(chunk0), gm_2(chunk0), gfva(chunk0), code(chunk0))

    ALLOCATE(comp_iso(na_uran,mn_uran), mass(na_uran), eneion((na_uran + 1) * 100))
    comp_iso(:,:) = -1d0
    mass(:) = 0d0

    ! 平均質量を計算する
    CALL get_nist_mass_comp(na_uran, mass, comp_iso)

    ! イオン化エネルギー
    DO na = 1, na_uran ! loop atom
       WRITE(natm,'(I3.3)') na
       DO ni = 0, na - 1 ! loop ion
          WRITE(nion,'(I2.2)') ni
          indx = na * 100 + ni
          CALL get_nist_eneion(natm, nion, eneion(indx))
       END DO
    END DO
    
    CALL h5open_f(error)
    !***********************************************************
    ! Kurucz
    !***********************************************************
    ! ファイル名を決める
    prefix = filename(1:INDEX(filename,'.')-1)
    prefix = prefix(INDEX(prefix,'/')+1:)
    output = DIRH5//TRIM(prefix)//'.h5'
    CALL h5Fcreate_f(TRIM(output), H5F_ACC_TRUNC_F, file_h5, error)

    ! ラインデータ書き込みテンプレート
    CALL h5Gcreate_f(file_h5, 'trans', grp_trans, error)
    dims = [1]
    CALL h5LTmake_dataset_f(grp_trans, 'id_1', 1, dims, H5T_STD_I32LE, ELECTRON, error)
    CALL h5LTmake_dataset_f(grp_trans, 'id_2', 1, dims, H5T_STD_I32LE, HI      , error)
    CALL count_lines(filename, nlines)
    dims = [nlines]
    CALL h5LXset_extendible_dataset(data_wnum, grp_trans, 'wnum', H5T_IEEE_F64LE, [chunk0], gzip)
    CALL h5LXset_extendible_dataset(data_lowe, grp_trans, 'lowe', H5T_IEEE_F64LE, [chunk0], gzip)
    CALL h5LXset_extendible_dataset(data_gfva, grp_trans, 'gfva', H5T_IEEE_F32LE, [chunk0], gzip)
    CALL h5LXset_extendible_dataset(data_code, grp_trans, 'code', H5T_STD_I32LE , [chunk0], gzip)
    CALL h5LXset_extendible_dataset(data_gm_0, grp_trans, 'gm_0', H5T_IEEE_F32LE, [chunk0], gzip)
    CALL h5LXset_extendible_dataset(data_gm_1, grp_trans, 'gm_1', H5T_IEEE_F32LE, [chunk0], gzip)
    CALL h5LXset_extendible_dataset(data_gm_2, grp_trans, 'gm_2', H5T_IEEE_F32LE, [chunk0], gzip)

    nlines_total = nlines_total + nlines
    ! チャンクサイズで書き込み
    nc_max = INT(nlines / chunk0,INT32)
    
    offset = [0]
    size = [0]
    
    OPEN(2, FILE=TRIM(filename), STATUS='OLD')
    lcount = 0
    DO nc = 1, nc_max+1
       !          print *, 'nc = ', nc, ' / ', nc_max+1
       IF(nc == nc_max + 1) THEN
          block_size = [nlines - nc_max * chunk0]
       ELSE
          block_size = [chunk0]
       END IF
       
       !***********************************************************
       ! ラインデータ読み込み
       DO l = 1, block_size(1)
          READ(2,'(F11.4,F7.3,F6.2,F12.3,F5.2,1X,A10,F12.3,F5.2,1X,A10, 3F6.2,A4,2I2,I3,F6.3,I3,F6.3,2I5)') &
               lambda, loggf, code0, lowe0, lowj, lowlabel, highe0, highj, highlabel, &
               gamr, gam4, gam6, cdum, idum, idum, isot2, lg_HF_frc, isot, lg_frc, EJE_HF_shift, EJO_HF_shift
          lcount = lcount + 1
          
          ! CORRECTION OF MASS NUMBER
          IF(isot == 864) isot = 138

          ! ISOTOPIC CORRECTION OF loggf
          IF(lg_frc /= 0) THEN
             frc_nist = REAL(comp_iso(INT(code0),isot))
             IF(frc_nist /= 0d0) THEN ! employ NIST value
                loggf = loggf + LOG10(frc_nist)
             ELSE
                loggf = loggf + lg_frc
             END IF
#ifdef WARNING             
             IF(lg_frc /= 0d0 .AND. ABS((lg_frc - LOG10(frc_nist))/lg_frc) > 0.10d0) THEN
                PRINT *, 'Strange isotopic fraction in Kurucz:', lg_frc, ' compared with NIST:', LOG10(frc_nist), &
                     ' atomic number:', INT(code0*100), ' mass number:', isot
             ENDIF
#endif             
          END IF

          ! Fine Stricture correction of loggf
          loggf = loggf + lg_HF_frc
          
          ! SHIFT CORRECTION OF ENERGIES
          lowe0 = ABS(lowe0) + DBLE(EJE_HF_shift) / 1d3
          highe0= ABS(highe0)+ DBLE(EJO_HF_shift) / 1d3
          
          wnum(l) = ABS(highe0 - lowe0)
          lowe(l) = MIN(lowe0, highe0)
          gfva(l) = 10e0 ** loggf
          code(l) = INT(code0*100)

          ! natural broadening
          gm_0(l) = 10e0 ** gamr
          ! pressure broadening by electron (quadratic Stark)
          gm_1(l) = 10e0 ** gam4
          IF(gam4 == 0.) gm_1(l) = 0.
          ! pressure broadening by HI (van der Waals)
          gm_2(l) = 10e0 ** gam6
          IF(gam6 == 0.) gm_2(l) = 0.
          
          ! IGNORE 'CONTINUUM'
          IF(TRIM(highlabel) == 'CONTINUUM') gfva(l) = 0d0
          ! IGNORE THE LINE IF THE TWO MASS NUMBERS DO NOT MATCH
          
          IF(wnum(l) <= 0d0) gfva(l) = 0d0
       END DO
          !***********************************************************
          
       size = size + block_size
          
       CALL h5LXwrite_extendible_dataset_F64(data_wnum, wnum, 1, size, offset, block_size) 
       CALL h5LXwrite_extendible_dataset_F64(data_lowe, lowe, 1, size, offset, block_size)
       CALL h5LXwrite_extendible_dataset_F32(data_gfva, gfva, 1, size, offset, block_size)
       CALL h5LXwrite_extendible_dataset_I32(data_code, code, 1, size, offset, block_size)
       CALL h5LXwrite_extendible_dataset_F32(data_gm_0, gm_0, 1, size, offset, block_size)
       CALL h5LXwrite_extendible_dataset_F32(data_gm_1, gm_1, 1, size, offset, block_size)
       CALL h5LXwrite_extendible_dataset_F32(data_gm_2, gm_2, 1, size, offset, block_size)

       offset = offset + block_size
       
    END DO
    CLOSE(2)
    
    IF(lcount /= nlines) THEN
       PRINT '(A, 2I12)', '*** ERROR *** lcount /= nlines: ', lcount, nlines
       STOP
    END IF
        
    dims = [1]
    CALL h5LTmake_dataset_f(grp_trans, 'chunk0', 1, dims, H5T_STD_I32LE, chunk0 , error)
    CALL h5LTmake_dataset_f(grp_trans, 'nlines', 1, dims, H5T_STD_I64LE, size(1), error)

    CALL h5Dclose_f(data_wnum, error)
    CALL h5Dclose_f(data_lowe, error)
    CALL h5Dclose_f(data_gfva, error)
    CALL h5Dclose_f(data_code, error)
    CALL h5Dclose_f(data_gm_0, error)
    CALL h5Dclose_f(data_gm_1, error)
    CALL h5Dclose_f(data_gm_2, error)
    CALL h5Gclose_f(grp_trans, error)
    CALL h5Fclose_f(file_h5, error)
    CALL h5close_f(error)

    WRITE(*,FMT='(A36,A1,I13,A6)',ADVANCE='NO') TRIM(filename), ':', nlines, ' lines'
    WRITE(*,FMT='(A5)',ADVANCE='NO') '  => '
    WRITE(command,"('ls -sh ',A)") TRIM(output)
    CALL SYSTEM(TRIM(command))

    DEALLOCATE(comp_iso, mass, eneion)
    DEALLOCATE(wnum, lowe, gm_0, gm_1, gm_2, gfva, code)

    RETURN    
  END SUBROUTINE convert_kurucz_h5

  
  SUBROUTINE count_lines(filename, nlines)
    CHARACTER(*), INTENT(IN) :: filename
    INTEGER(INT64), INTENT(OUT) :: nlines
    
    INTEGER(INT32) :: iostat
    
    OPEN(1, FILE=TRIM(filename), STATUS='OLD')
    nlines = 0
    DO
       READ(1,*,IOSTAT=iostat)
       IF(iostat /= 0) EXIT
       nlines = nlines + 1
    END DO
    CLOSE(1)
    
    RETURN
  END SUBROUTINE count_lines

  
  SUBROUTINE get_hitran_prop(filename, code, mass, frac, temp, pf, ntemp)
    USE code_module, ONLY : HITRAN_TO_INTERNAL
    CHARACTER(*), INTENT(IN) :: filename
    INTEGER(INT32), INTENT(OUT) :: code
    REAL(REAL64), INTENT(OUT) :: mass
    REAL(REAL64), INTENT(OUT) :: frac
    REAL(REAL64), ALLOCATABLE, INTENT(OUT) :: temp(:)
    REAL(REAL64), ALLOCATABLE, INTENT(OUT) :: pf(:)
    INTEGER(INT32), INTENT(OUT) :: ntemp
    
    INTEGER(INT32) :: nt
    
    OPEN(1, FILE=TRIM(filename(1:INDEX(filename,'.',.TRUE.)-1))//'.prop', STATUS='OLD')
    READ(1,*) code
    READ(1,*) frac
    READ(1,*) mass
    READ(1,*) ntemp
    ALLOCATE(temp(ntemp), pf(ntemp))
    DO nt = 1, ntemp
       READ(1,*) temp(nt), pf(nt)
    END DO
    CLOSE(1) 

    code = code + HITRAN_TO_INTERNAL
   
    RETURN
  END SUBROUTINE get_hitran_prop
  

  SUBROUTINE get_exomol_prefix(filename, prefix, code)
    USE code_module, ONLY : HITRAN_TO_INTERNAL
    CHARACTER(*), INTENT(IN) :: filename
    CHARACTER*256, INTENT(OUT) :: prefix
    INTEGER(INT32), INTENT(OUT) :: code
    
    CHARACTER*256 :: prefix0, prefix1, species, isotopologue
    
    prefix0 = filename(INDEX(filename,'/',.TRUE.)+1:)

    IF(INDEX(TRIM(prefix0),'_',.TRUE.) - INDEX(TRIM(prefix0),'_') == 1) THEN
       prefix = filename(1:INDEX(filename, '.', .TRUE.)-1)
    ELSE
       IF(INDEX(TRIM(prefix0),'_p') /= 0) THEN
          prefix = filename(1:INDEX(filename, '.', .TRUE.)-1)
       ELSE
          prefix = filename(1:INDEX(filename, '_', .TRUE.)-2)
       END IF
    END IF

    prefix1 = prefix(1:INDEX(prefix,'__')-1)
    
    IF(INDEX(prefix1,'cis-') /= 0) THEN
       prefix1 = prefix1(INDEX(prefix1,'cis-')+4:)
    END IF
    IF(INDEX(prefix1,'trans-') /= 0) THEN
       prefix1 = prefix1(INDEX(prefix1,'trans-')+6:)
    END IF

    OPEN(1,FILE='molecular_id.tsv',STATUS='OLD')
    READ(1,*)
    DO
       READ(1,*) code, species, isotopologue
       IF(code == 999) THEN
          print *, '*** ERROR: could not find species id for ', TRIM(prefix1)
          STOP
       END IF
       IF(TRIM(isotopologue) == TRIM(prefix1)) EXIT
    END DO
    CLOSE(1)

    code = code + HITRAN_TO_INTERNAL

    print '(a,i5)', 'isotopologue:'//TRIM(prefix1)//' species: '//TRIM(species)//' id: ', code
    print *, ''

    RETURN
  END SUBROUTINE get_exomol_prefix


  SUBROUTINE read_exomol_pf(prefix, temp, pf)
    CHARACTER(*), INTENT(IN) :: prefix
    REAL(REAL64), ALLOCATABLE, INTENT(OUT) :: temp(:)
    REAL(REAL64), ALLOCATABLE, INTENT(OUT) :: pf(:)

    INTEGER(INT32) :: iostat, ntemp, l
    
    print '(a)', 'read partition function: '//trim(prefix)//'.pf:'

    OPEN(1, FILE=TRIM(prefix)//'.pf', STATUS='OLD')
    ntemp = 0
    DO
       READ(1, *, IOSTAT=iostat)
       IF(iostat /= 0) EXIT
       ntemp = ntemp + 1
    END DO
    ALLOCATE(temp(ntemp), pf(ntemp))
    REWIND(1)
    DO l = 1, ntemp
       READ(1,*) temp(l), pf(l)
    END DO
    CLOSE(1)

    print '(a,i10)', ' ntemp =', ntemp
    print *, ''
    
    RETURN
  END SUBROUTINE read_exomol_pf
  

  SUBROUTINE read_exomol_broad(prefix, id_10, gm_10, nx_10, jt_10, id_20, gm_20, nx_20, jt_20)
    USE code_module, ONLY : H2, HeI
    CHARACTER(*), INTENT(IN) :: prefix
    INTEGER(INT32), INTENT(OUT) :: id_10, id_20
    REAL(REAL32), ALLOCATABLE, INTENT(OUT) :: gm_10(:), gm_20(:)
    REAL(REAL32), ALLOCATABLE, INTENT(OUT) :: nx_10(:), nx_20(:)
    REAL(REAL32), ALLOCATABLE, INTENT(OUT) :: jt_10(:), jt_20(:) ! J quantum number integer/half integer

    INTEGER(INT32) :: iostat, jmax_1, jmax_2, l, l_skip
    CHARACTER*256 :: linebuf, prefix0, bcode

    prefix0 = prefix(1:INDEX(prefix,'__')-1)

    ! broadener = H2
    id_10 = H2
    OPEN(1, FILE=TRIM(prefix0)//'__H2.broad', STATUS='OLD', IOSTAT=iostat)
    IF(iostat /= 0) THEN
       print '(a)', '+++ warning: '//trim(prefix0)//'__H2.broad does not exist'
       jmax_1 = 1
       ALLOCATE(gm_10(jmax_1), nx_10(jmax_1), jt_10(jmax_1))
       gm_10(1) = 0d0
       nx_10(1) = 0d0
       jt_10(1) = 0
    ELSE
       PRINT '(a)', 'read broadening data: '//TRIM(prefix0)//'__H2.broad'
       jmax_1 = 0
       l_skip = 0
       DO
          READ(1,FMT='(A)',IOSTAT=iostat) linebuf
          IF(iostat /= 0) EXIT
          IF(linebuf(1:2) /= 'a0') THEN
             l_skip = l_skip + 1
          ELSE
             jmax_1 = jmax_1 + 1
          END IF
       END DO
       ALLOCATE(gm_10(jmax_1), nx_10(jmax_1), jt_10(jmax_1))
       REWIND(1)
       DO l = 1, l_skip
          READ(1,*)
       END DO
       DO l = 1, jmax_1
          READ(1,*) bcode, gm_10(l), nx_10(l), jt_10(l)
       END DO
       CLOSE(1)
    END IF
    
    ! broadener = He
    id_20 = HeI
    OPEN(1, FILE=TRIM(prefix0)//'__He.broad', STATUS='OLD', IOSTAT=iostat)
    IF(iostat /= 0) THEN
       print '(a)', '+++ warning: '//trim(prefix0)//'__He.broad does not exist'
       jmax_2 = 1
       ALLOCATE(gm_20(jmax_2), nx_20(jmax_2), jt_20(jmax_2))
       gm_20(1) = 0d0
       nx_20(1) = 0d0
    ELSE
       PRINT '(a)', 'read broadening data: '//TRIM(prefix0)//'__He.broad'
       jmax_2 = 0
       l_skip = 0
       DO
          READ(1,FMT='(A)',IOSTAT=iostat) linebuf
          IF(iostat /= 0) EXIT
          IF(linebuf(1:2) /= 'a0') THEN
             l_skip = l_skip + 1
          ELSE
             jmax_2 = jmax_2 + 1
          END IF
       END DO
       ALLOCATE(gm_20(jmax_2), nx_20(jmax_2), jt_20(jmax_2))
       REWIND(1)
       DO l = 1, l_skip
          READ(1,*)
       END DO
       DO l = 1, jmax_2
          READ(1,*) bcode, gm_20(l), nx_20(l), jt_20(l)
       END DO
       CLOSE(1)
    END IF
    print *, ''
    
    RETURN
  END SUBROUTINE read_exomol_broad
  

  SUBROUTINE read_exomol_def(prefix, mass, nlines_def)
    USE const_module, ONLY : amu
   
    CHARACTER(*), INTENT(IN) :: prefix
    REAL(REAL64), INTENT(OUT) :: mass
    INTEGER(INT64), INTENT(OUT) :: nlines_def
    
    CHARACTER*256 :: linebuf
    INTEGER(INT32) :: iostat

    print '(a)', 'read definition: '//trim(prefix)//'.def: '

    OPEN(1, FILE=TRIM(prefix)//'.def', STATUS='OLD')

    ! read mass
    DO
       READ(1,FMT='(A)',IOSTAT=iostat) linebuf
       IF(iostat /= 0) THEN
          print *, '*** read error: '//trim(prefix)//'.def'
          STOP
       END IF
       IF(INDEX(linebuf,'Isotopologue mass') /= 0) EXIT
    END DO
    BACKSPACE(1)
    READ(1,*) mass
    print '(a,f6.2)', ' mass (Da) = ', mass
    mass = mass * amu

    REWIND(1)

    ! read # of transitions
    DO
       READ(1,FMT='(A)',IOSTAT=iostat) linebuf
       IF(iostat /= 0) THEN
          print *, '*** read error: '//trim(prefix)//'.def'
          STOP
       END IF
       IF(INDEX(linebuf,'Total number of transitions') /= 0) EXIT
    END DO
    BACKSPACE(1)
    READ(1,*) nlines_def
    
    CLOSE(1)
    
    print '(a1, a, i12)', ' ', 'nlines_def = ', nlines_def
    print *, ''

    RETURN
  END SUBROUTINE read_exomol_def
  
      
  SUBROUTINE read_exomol_states(prefix, ene, gtot, jtot)
    CHARACTER(*), INTENT(IN) :: prefix
    REAL(REAL64), ALLOCATABLE, INTENT(OUT) :: ene(:)
    INTEGER(INT32), ALLOCATABLE, INTENT(OUT) :: gtot(:)
    REAL(REAL32), ALLOCATABLE, INTENT(OUT) :: jtot(:)
    
    INTEGER(INT64) :: l, id, nstates, iostat
    
    print '(a)', 'read level data: '//trim(prefix)//'.states:'

    OPEN(1, FILE=TRIM(prefix)//'.states', STATUS='OLD')
    nstates = 0
    DO
       READ(1, *, IOSTAT=iostat)
       IF(iostat /= 0) EXIT
       nstates = nstates + 1
    END DO
    ALLOCATE(ene(nstates), gtot(nstates), jtot(nstates))
    REWIND(1)
    DO l = 1, nstates
       READ(1,*) id, ene(id), gtot(id), jtot(id)
    END DO
    CLOSE(1)

    print '(a,i10)', ' nstates =', nstates
    print *, ''
    
    RETURN
  END SUBROUTINE read_exomol_states
  
  
  SUBROUTINE get_exomol_isofrac(prefix, frac)
    
    CHARACTER(LEN=*), INTENT(IN) :: prefix
    REAL(REAL64), INTENT(OUT) :: frac

    CHARACTER(len=32), DIMENSION(10) :: form

    INTEGER(INT32), PARAMETER :: almax = 26, elmax = 10
    CHARACTER(LEN=2), DIMENSION(elmax) :: sym
    INTEGER(INT32), DIMENSION(elmax) :: num, iso
    
    CHARACTER(len=1), DIMENSION(almax), PARAMETER :: alpha = [&
         'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', &
         'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z']
    
    INTEGER :: i1, i2, i3, el0, el, al

    REAL(REAL64) :: comp
    INTEGER(INT32) mass1, number, mass0, is, iostat
    CHARACTER*2 :: sym1, sym0
    CHARACTER :: dummy*256, mass_c*20, comp_c*20
    CHARACTER*20, PARAMETER :: a20 = '                    '

    el0 = 1
    
    OPEN(1, FILE=TRIM(prefix)//'.def', STATUS='OLD')
    READ(1,*)
    READ(1,*) form(1)
    CLOSE(1)

    IF(INDEX(form(1),'cis-') /= 0) THEN
       form(1) = form(1)(INDEX(form(1),'cis-')+4:)
    END IF
    IF(INDEX(form(1),'trans-') /= 0) THEN
       form(1) = form(1)(INDEX(form(1),'trans-')+6:)
    END IF

    IF(INDEX(form(1),'+') /= 0) THEN
       form(1) = form(1)(1:INDEX(form(1),'+')-1)
    END IF

    DO el = 1, elmax - 1
       i1 = INDEX(form(el),'(',BACK=.TRUE.)
       IF(i1 /= 1) THEN
          form(el+1) = TRIM(form(el)(:i1-1))
          form(el) = TRIM(form(el)(i1:))
          el0 = el + 1
       ELSE
          EXIT
       END IF
    END DO

    OPEN(51, FILE=DIRNIST//'nist_isotope.txt', STATUS='OLD')
    frac = 1d0
    DO el = 1, el0
       DO al = 1, almax
          i2 = INDEX(form(el),alpha(al))
          IF(i2 /= 0) EXIT
       END DO
       i3 = INDEX(form(el),')')
       
       sym(el) = TRIM(form(el)(i2:i3-1))
       READ(form(el)(2:i2-1),*) iso(el)

       IF(form(el)(LEN(TRIM(form(el))):) == ')') THEN
          num(el) = 1
       ELSE
          read(form(el)(INDEX(form(el),')')+1:),*) num(el)
       ENDIF

       mass0 = iso(el)
       sym0 = TRIM(sym(el))
       number = num(el)

       DO
          READ(51,*,IOSTAT=iostat)
          IF(iostat /= 0) THEN
             print '(a)', '*** error *** could not find symbol: ', sym0
             STOP
          END IF
          READ(51,'(A16,A2)') dummy, sym1
          READ(51,'(A14,I3)') dummy, mass1
          READ(51,'(A23,A20)') dummy, mass_c
          READ(51,'(A23,A20)') dummy, comp_c
          IF(comp_c /= a20) THEN
             is = INDEX(comp_c,'(')
             IF(is == 0) THEN
                READ(comp_c,*) comp
             ELSE
                READ(comp_c(1:is-1),*) comp
             END IF
          END IF
          READ(51,*)
          READ(51,*)
          READ(51,*)
          IF(TRIM(sym1) == TRIM(sym0) .AND. mass1 == mass0) EXIT
       END DO
       print '(a4,f9.7,i2)', ' '//sym1//' ', comp, number
       frac = frac * comp ** number
       REWIND(51)
    END DO

    print *, 'frac=', frac

    CLOSE(51)
    
    RETURN
  END SUBROUTINE get_exomol_isofrac
  

END PROGRAM convert_lines_h5
