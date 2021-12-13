
PROGRAM preprop_hitran

  USE ISO_FORTRAN_ENV
  USE const_module, ONLY : amu, clight, m_ele, e2, pi
  IMPLICIT NONE

  INTEGER(INT32) :: n, nfiles, iostat, num, iostat0
  CHARACTER*256, ALLOCATABLE, DIMENSION(:) :: filename
  INTEGER(INT32) :: code, code0, code00, count, ntemp, nt, co
  INTEGER(INT32) :: global, local, afgl, gi, cnt, temp
  REAL(REAL64) :: frac1, mass1, q, q1
  CHARACTER*256 :: pf1, form1, dec_prop, dec_par
  REAL(REAL64) :: wnum, int, acoeff, aw, sw, lowe, td, ps, gu, gl, gfva
  INTEGER(INT64) :: ref, nlines, nlines_decomp
  INTEGER(INT32) :: mol, iso, err
  CHARACTER :: uq*15, lq*15, flg*1, ulq*15, llq*15, colon*1
  INTEGER(INT32), PARAMETER :: max_buf = 256
  CHARACTER(max_buf):: linebuf

  INTEGER(INT32) :: length, status
  CHARACTER(:), ALLOCATABLE :: filename0
  INTRINSIC :: COMMAND_ARGUMENT_COUNT, GET_COMMAND_ARGUMENT

  CALL GET_COMMAND_ARGUMENT(1, LENGTH=length, STATUS=status)
  IF(status == 0) THEN
     ALLOCATE(CHARACTER(length) :: filename0)
     CALL GET_COMMAND_ARGUMENT(1, filename0, STATUS=status)
     IF(status == 0) THEN
        PRINT *, 'input file = ', filename0
     END IF
  ELSE
     PRINT *, 'argument error'
     STOP
  END IF

  ! コンバートするファイル名を読み込み
!!$  OPEN(1, FILE='list_decomp.txt', STATUS='OLD')
!!$  nfiles = 0
!!$  DO
!!$     READ(1,*,IOSTAT=iostat)
!!$     IF(iostat /= 0) EXIT
!!$     nfiles = nfiles + 1
!!$  END DO
!!$  REWIND(1)
!!$  PRINT *, '# of original par files: ', nfiles
!!$  ALLOCATE(filename(nfiles))
!!$  DO n = 1, nfiles
!!$     READ(1,FMT='(a)') filename(n)
!!$  END DO
!!$  CLOSE(1)

  nfiles = 1
  ALLOCATE(filename(nfiles))
  filename(1) = filename0

  ! 分子種コード
  num = INDEX(TRIM(filename(1)),'_')
  READ(filename(1)(num-2:num-1),*) code0
  DO n = 2, nfiles
     READ(filename(n)(num-2:num-1),*) code00
     IF(code00 /= code0) THEN
        PRINT *, '*** ERROR: # OF SPECIES MUST BE ONE: ', code00, code0
        STOP
     END IF
  END DO

  ! 当該分子種の同位体置換体の分だけファイルを用意する
  OPEN(12, FILE='hitran_meta.txt', STATUS='OLD')
  DO 
     READ(12,'(A)',IOSTAT=iostat0) linebuf
     IF(iostat0 /= 0) THEN
        PRINT *, 'END OF FILE'
        STOP
     END IF
     READ(linebuf,*,IOSTAT=iostat) code, colon, form1
     IF(iostat /= 0) CONTINUE
     IF(code == code0 .AND. colon == ':') EXIT     
  END DO
  PRINT '(A11,I3,A2,A)', ' SPECIES: (', code, ') ',  TRIM(form1)
  READ(12,*)
  READ(12,*)
  count = 0
  DO
     READ(12,*,IOSTAT=iostat) global, local
     IF(iostat /= 0) EXIT
     count = count + 1
  END DO
  PRINT *, '# of isotopologues: ', count

  REWIND(12)
  OPEN(4,FILE='list_convert.txt',STATUS='REPLACE')
  DO 
     READ(12,'(A)',IOSTAT=iostat0) linebuf
     IF(iostat0 /= 0) THEN
        PRINT *, 'END OF FILE'
        STOP
     END IF
     READ(linebuf,*,IOSTAT=iostat) code, colon
     IF(iostat /= 0) CONTINUE
     IF(code == code0 .AND. colon == ':') EXIT     
  END DO
  READ(12,*)
  READ(12,*)

  DO cnt = 1, count
     READ(12,*) global, local, form1, afgl, frac1, mass1, q1, pf1, gi
     !     IF(code0 == 10 .AND. global == 130) CYCLE ! for NO2, q130.txt not exist
     IF(local == 0) local = 10
     PRINT '(A10,I4,A2,A)', 'local ID: ', local, ': ', TRIM(form1)
     dec_prop = TRIM(form1)//'__'//filename(1)(INDEX(filename(1),'_',.TRUE.)+1:INDEX(filename(1),'.',.TRUE.)-1)//'.prop'
     dec_par  = TRIM(form1)//'__'//filename(1)(INDEX(filename(1),'_',.TRUE.)+1:INDEX(filename(1),'.',.TRUE.)-1)//'.par'
     !        PRINT *, 'output: ', 'decomp/'//TRIM(dec_prop)
!     IF(local == 10) local = 0 ! for HITRAN CO2
     OPEN(100+local, FILE='decomp/'//TRIM(dec_prop), STATUS='REPLACE')
     WRITE(100+local,*) code, ' HITRAN species code'
     WRITE(100+local,*) frac1, ' isotopic fraction'
     WRITE(100+local,*) mass1 * amu, ' mass [g]'
     OPEN(3, FILE='Q/'//TRIM(pf1), STATUS='OLD')
     ntemp = 0
     DO
        READ(3,*,IOSTAT=iostat)
        IF(iostat /= 0) EXIT
        ntemp = ntemp + 1
     END DO
     REWIND(3)
     WRITE(100+local,*) ntemp, '# of entries in partition function table'
     DO nt = 1, ntemp
        READ(3,'(I4,F22.8)') temp, q
        WRITE(100+local,*) temp, q
     END DO
     CLOSE(3)
     CLOSE(100+local)
     !        PRINT *, 'output: ', 'decomp/'//TRIM(dec_par)
     
     OPEN(100+local, FILE='decomp/'//TRIM(dec_par), STATUS='REPLACE')
     WRITE(4,FMT='(A)') 'decomp/'//TRIM(dec_par)
  END DO
  CLOSE(12)
  CLOSE(4)
  
  ! オリジナルファイル読み込み
  nlines = 0
  DO n = 1, nfiles
     OPEN(99, FILE=TRIM(filename(n)), STATUS='OLD')
     DO
        READ(99,FMT='(I2,I1,F12.6,E10.3,E10.3,F5.4,F5.3,F10.4,F4.2,F8.6,A15,A15,A15,A15,I6,I12,A1,F7.1,F7.1)', IOSTAT=iostat) &
             mol, iso, wnum, int, acoeff, aw, sw, lowe, td, ps, uq, lq, ulq, llq, err, ref, flg, gu, gl
!       IF(code0 == 10 .AND. iso == 2) CYCLE ! for NO2, q130.txt not exist
        IF(iostat /= 0) EXIT
        gfva = gu * acoeff / wnum**2 * (clight * m_ele / (8d0 * pi**2 * e2))
        WRITE(100+iso,*) wnum, lowe, gfva, aw, sw, td
        nlines = nlines + 1
     END DO
     print '(i2.2,a,i2.2,a)', n, '/', nfiles, ' decomposed '//trim(filename(n))
     CLOSE(99)
  END DO

  nlines_decomp = 0
  DO co = 1, count
     REWIND(100+co)
     DO
        READ(100+co,*,IOSTAT=iostat)
        IF(iostat /= 0) EXIT
        nlines_decomp = nlines_decomp + 1
     END DO
     CLOSE(100+co)
  END DO

  IF(nlines_decomp /= nlines) THEN
     PRINT *, '*** ERROR: processed numbers of lines do not match: ', nlines, nlines_decomp
  ELSE
     PRINT *, 'processed numbers of lines match: ', nlines, nlines_decomp
  END IF

  STOP
CONTAINS  
  
END PROGRAM preprop_hitran
