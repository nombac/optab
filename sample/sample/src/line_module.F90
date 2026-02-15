
MODULE line_module

  USE ISO_FORTRAN_ENV
  USE HDF5
  USE H5LT
  USE MPI
  USE mpi_module, ONLY : myrk
  
  PRIVATE
  PUBLIC :: line_molec, line_kurucz, line_kurucz_p

CONTAINS

  SUBROUTINE line_molec(source, temp0, np, grd, dgrd, alp, cnt, pmean, temp2, pmean2, slines)
    USE code_module, ONLY : TOTAL
    USE const_module, ONLY : c2, k_bol, clight, e2, m_ele, pi, h
    USE h5lx_module, ONLY : h5LXread_hyperslab_dataset_F64, h5LXread_hyperslab_dataset_F32, &
         h5LXopen_file
    USE mpi_module, ONLY : myrk_line, mprc, myrk_m, lprc, myrk_l, mpi_line_world, para_range, &
         kprc, mpi_jconst_world, myrk_jconst, nprc_jconst
    IMPLICIT NONE

    CHARACTER(*), INTENT(IN) :: source ! line source
    REAL(REAL64), INTENT(IN) :: temp0(:) ! temperature [K]
    REAL(REAL64), INTENT(IN) :: np(:,:) ! partial density
    REAL(REAL64), INTENT(IN) :: grd(:), dgrd(:) ! wavenumber grid [cm^{-1}]
    REAL(REAL64), INTENT(INOUT) :: alp(:,:) ! opacity [cm^{-1}]
    REAL(REAL64), INTENT(IN) :: cnt(:,:) ! continuum opacity [cm^{-1}]
    REAL(REAL64), INTENT(INOUT) :: pmean(:)
    REAL(REAL64), INTENT(IN) :: temp2 ! radiation temperature
    REAL(REAL64), INTENT(INOUT) :: pmean2(:) ! two-temp planck_mean opacity
    INTEGER(INT64), INTENT(INOUT), OPTIONAL :: slines(:,:)
    
    ! HDF5 variables
    INTEGER(HID_T) :: file_spec, grp_trans, grp_prop
    INTEGER(HID_T) :: data_wnum, data_lowe, data_gfva, data_gm_1, data_nx_1, data_gm_2, data_nx_2
    INTEGER(HSIZE_T), DIMENSION(1) :: offset, block_size, dims
    
    ! grid index
    INTEGER(INT32) :: ks, ke, kmin, kmax
    ! lines
    INTEGER(INT64) :: l, nlines, n_rest
    INTEGER(INT32) :: nc, nc_max, chunk0, ls, le
    ! main variables
    REAL(REAL64), ALLOCATABLE, DIMENSION(:) :: wnum, lowe
    REAL(REAL32), ALLOCATABLE, DIMENSION(:) :: gfva, gm_1, nx_1, gm_2, nx_2
    ! Read from .prop
    REAL(REAL64) :: mass, frac
    REAL(REAL64), ALLOCATABLE :: temp(:), pf(:), q(:)
    INTEGER(INT32) :: id, ntemp
    ! Lorentz profile
    REAL(REAL64), PARAMETER :: temp_ref = 296d0 ! [K]
    REAL(REAL64), PARAMETER :: pres_ref = 1d6   ! [1 bar = 10^6 dyne]
    INTEGER(INT32) :: id_1, id_2
    REAL(REAL64) gn, gamma
#ifndef CROSS_CHECK    
    REAL(REAL64) g6, pres_1, pres_2
    INTEGER(INT32) :: k_wnum
#endif    
    ! Gaussian profile
    REAL(REAL64), ALLOCATABLE :: prof(:)
    REAL(REAL64) :: sigma
    REAL(REAL64) :: strength, strength0
    ! parameters read from files
    INTEGER(INT32) :: icutoff
    CHARACTER(LEN=2) :: iprof = 'Vo'
    CHARACTER :: fname*128
    REAL(REAL64) cutoff0_voigt, cutoff0_gauss, cutoff_p, cutoff_m, &
         strength_cutoff0, wnum_crit, delta, delta_crit, delta_voigt
    ! others
    INTEGER(INT32) :: error, j, js, je, prog
    REAL(REAL64), ALLOCATABLE :: dummy(:)
#ifdef DEBUG
    INTEGER(INT32), PARAMETER :: countmax = 13
    REAL(REAL64) :: t_count(0:countmax), t0, t1, t2
    INTEGER(INT64) :: n_count(0:countmax)
#endif
#ifdef PLANCK_MEAN    
    REAL(REAL64) a, ap
#endif
    NAMELIST /cutoffs/ icutoff, cutoff0_voigt, cutoff0_gauss, strength_cutoff0, wnum_crit, &
         delta_crit, delta_voigt
#ifdef DEBUG
    t1 = mpi_wtime()
    t0 = mpi_wtime()
    t_count(:) = 0d0
    n_count(:) = 0
    t2 = 0d0
#endif

    alp = 0d0

    js = 1
    je = UBOUND(temp0,1)

    ks = 1
    ke = UBOUND(grd,1)

    ALLOCATE(prof(ks:ke))

    !-----------------------------------------------------------------------------
    ! read cutoff parameters
    OPEN(5, FILE='input/fort.5', STATUS='OLD')
    READ(5, NML=cutoffs)
    CLOSE(5)

    !-----------------------------------------------------------------------------
    ! start HDF5
    CALL h5open_f(error)
    CALL h5Eset_auto_f(0,error)
    ! open line list file
    fname = 'input/h5/'//TRIM(source)//'.h5'
#ifndef MPIO    
    CALL h5Fopen_f(TRIM(fname), H5F_ACC_RDONLY_F, file_spec, error)
#else
    CALL h5LXopen_file(TRIM(fname), H5F_ACC_RDONLY_F, file_spec, error, &
         MPI_COMM_WORLD, MPI_INFO_NULL)
#endif
    IF(error /= 0) THEN
       PRINT *, '*** '//TRIM(fname)//' DOES NOT EXISTS ***'
       CALL MPI_FINALIZE(ERROR); STOP
    END IF
    !-----------------------------------------------------------------------------
    ! read id, frac, mass, and compute q
    CALL h5Gopen_f(file_spec, 'prop', grp_prop, error)
    dims = [1]
    CALL h5LTread_dataset_f(grp_prop, 'id'   , H5T_STD_I32LE , id   , dims, error)

    CALL h5LTread_dataset_f(grp_prop, 'frac' , H5T_IEEE_F64LE, frac , dims, error)
    CALL h5LTread_dataset_f(grp_prop, 'mass' , H5T_IEEE_F64LE, mass , dims, error)
    CALL h5LTread_dataset_f(grp_prop, 'ntemp', H5T_STD_I32LE , ntemp, dims, error)
    ALLOCATE(temp(ntemp), pf(ntemp))
    dims = [ntemp]
    CALL h5LTread_dataset_f(grp_prop, 'temp' , H5T_IEEE_F64LE, temp , dims, error)
    CALL h5LTread_dataset_f(grp_prop, 'pf'   , H5T_IEEE_F64LE, pf   , dims, error)
    ALLOCATE(q(UBOUND(temp0,1)))
    CALL interpolate_q(temp0(js:je), temp, pf, q(js:je))
    DEALLOCATE(temp, pf)
    CALL h5Gclose_f(grp_prop, error)

    !-----------------------------------------------------------------------------
    ! prepare for reading line list
    CALL h5Gopen_f(file_spec, 'trans', grp_trans, error)
    dims = [1]
    CALL h5LTread_dataset_f(grp_trans, 'nlines', H5T_STD_I64LE, nlines, dims, error)
    IF(nlines == 0) GOTO 999
    CALL h5LTread_dataset_f(grp_trans, 'chunk0', H5T_STD_I32LE, chunk0, dims, error)
    CALL h5LTread_dataset_f(grp_trans, 'id_1'  , H5T_STD_I32LE, id_1  , dims, error)
    CALL h5LTread_dataset_f(grp_trans, 'id_2'  , H5T_STD_I32LE, id_2  , dims, error)
    CALL h5Dopen_f(grp_trans, 'wnum', data_wnum, error)
    CALL h5Dopen_f(grp_trans, 'lowe', data_lowe, error)
    CALL h5Dopen_f(grp_trans, 'gfva', data_gfva, error)
    CALL h5Dopen_f(grp_trans, 'gm_1', data_gm_1, error)
    CALL h5Dopen_f(grp_trans, 'nx_1', data_nx_1, error)
    CALL h5Dopen_f(grp_trans, 'gm_2', data_gm_2, error)
    CALL h5Dopen_f(grp_trans, 'nx_2', data_nx_2, error)

!    IF(myrk == 0) write(*,FMT='(a25,i11)') TRIM(source), nlines
     ! prepare arrays
    
    !-----------------------------------------------------------------------------
    ! main part
    block_size = [MIN(chunk0, nlines)]
    nc_max = INT(nlines / block_size(1))
    n_rest = nlines - nc_max * block_size(1)
    ALLOCATE(&
         gfva(block_size(1)), wnum(block_size(1)), lowe(block_size(1)), &
         gm_1(block_size(1)), nx_1(block_size(1)), gm_2(block_size(1)), nx_2(block_size(1)))
    CALL para_range(1, INT(block_size(1)), lprc, myrk_l, ls, le)
    offset = myrk_m * block_size
    
    prog = 0
#ifdef DEBUG    
    t_count(0) = t_count(0) + (mpi_wtime() - t0)
    n_count(0) = n_count(0) + 1
#endif
   
    BLOCK: DO nc = 1 + myrk_m, nc_max + 1, mprc
!       if(myrk == 0) call progress(nc, nc_max, prog, dprog=1)
       ! ----- READ DATA BLOCK -----
#ifdef DEBUG       
       t0 = mpi_wtime()
#endif       
       IF(nc == nc_max + 1) THEN
          IF(n_rest == 0) EXIT
          DEALLOCATE(gfva, wnum, lowe, gm_1, nx_1, gm_2, nx_2)
          block_size = [n_rest]
          ALLOCATE(&
               gfva(block_size(1)), wnum(block_size(1)), lowe(block_size(1)), &
               gm_1(block_size(1)), nx_1(block_size(1)), gm_2(block_size(1)), nx_2(block_size(1)))
          CALL para_range(1, INT(block_size(1)), lprc, myrk_l, ls, le)
       END IF
       
       CALL h5LXread_hyperslab_dataset_F64(data_wnum, wnum, 1, block_size, offset, block_size)
       CALL h5LXread_hyperslab_dataset_F64(data_lowe, lowe, 1, block_size, offset, block_size)
       CALL h5LXread_hyperslab_dataset_F32(data_gfva, gfva, 1, block_size, offset, block_size)
       CALL h5LXread_hyperslab_dataset_F32(data_gm_1, gm_1, 1, block_size, offset, block_size)
       CALL h5LXread_hyperslab_dataset_F32(data_nx_1, nx_1, 1, block_size, offset, block_size)
       CALL h5LXread_hyperslab_dataset_F32(data_gm_2, gm_2, 1, block_size, offset, block_size)
       CALL h5LXread_hyperslab_dataset_F32(data_nx_2, nx_2, 1, block_size, offset, block_size)
       offset = offset + block_size * mprc
#ifdef LINE_STATISTICS
       DO j = js, je
          slines(j,1) = slines(j,1) + block_size(1)
       END DO
#endif    
#ifdef DEBUG       
       t_count(1) = t_count(1) + (mpi_wtime() - t0)
       n_count(1) = n_count(1) + 1
#endif       
       ! ---------------------------
       
       LAYER: DO j = js, je
          ! *** skip layer if Q is note defined OR abundance is zero
#ifdef DEBUG             
             t0 = mpi_wtime()
#endif             
          IF(q(j) * np(id,j) == 0d0) THEN
#ifdef LINE_STATISTICS
             slines(j,2) = slines(j,2) + (le - ls + 1)
#endif
#ifdef DEBUG
             t2 = t2 + (mpi_wtime() - t0)
#endif             
             CYCLE
          ELSE
#ifdef DEBUG
             t2 = t2 + (mpi_wtime() - t0)
#endif             
          END IF
          ! --------------------------------------------------------------
          LINE: DO l = ls, le
             ! *** skip lines of zero-wavenumber
#ifdef DEBUG             
             t0 = mpi_wtime()
#endif             
             IF(wnum(l) == 0d0) THEN
#ifdef LINE_STATISTICS                
                slines(j,2) = slines(j,2) + 1
#endif                
#ifdef DEBUG
                t2 = t2 + (mpi_wtime() - t0)
#endif             
                CYCLE
             ELSE
#ifdef DEBUG
                t2 = t2 + (mpi_wtime() - t0)
#endif             
             END IF
             
             ! --- COMPUTE LINE STRENGTH [cm] -----
#ifdef DEBUG             
             t0 = mpi_wtime()
#endif             
             strength0 = (pi * e2 / (m_ele * clight**2)) * gfva(l) * EXP(-c2 * lowe(l) / temp0(j)) / q(j)
             strength = strength0 * (1d0 - EXP(-c2 * wnum(l) / temp0(j))) 
#ifdef DEBUG             
             t_count(2) = t_count(2) + (mpi_wtime() - t0)
             n_count(2) = n_count(2) + 1
#endif             
             ! ------------------------------------
#ifdef STRENGTH_CUTOFF             
#ifndef CROSS_CHECK
             ! *** IGNORE LINES WHOSE STRENGTH IS BELOW CUTOFF
             IF(strength < strength_cutoff(wnum(l), temp0(j), wnum_crit, strength_cutoff0)) THEN
#ifdef LINE_STATISTICS                
                slines(j,2) = slines(j,2) + 1
#endif                
                CYCLE
             END IF
#endif             
#endif    
             ! --- COMPUTE  DOPPLER WIDTH [cm^-1] -----
#ifdef DEBUG             
             t0 = mpi_wtime()
#endif             
             sigma = SQRT(2d0 * k_bol * temp0(j) / mass) / clight * wnum(l)

             ! --- COMPUTE PLANCK-MEAN OPACITY -----
#ifdef PLANCK_MEAN
             a = wnum(l) - c2 * sigma**2 / temp0(j)
             ap = a / sigma
             pmean(j) = pmean(j) + (np(id,j) * frac) * strength0 * &
                  EXP(-c2 / temp0(j) * (wnum(l) + a) / 2d0) / SQRT(pi) / sigma * &
                  0.25d0 * sigma**4 * (SQRT(pi) * ap * (2d0 * ap**2 + 3d0) * (ERF(ap) + 1d0) + &
                  2d0 * (ap**2 + 1d0) * EXP(-ap**2))
#else             
             pmean(j) = pmean(j) + (np(id,j) * frac) * strength0 * wnum(l)**3 * EXP(-c2 * wnum(l) / temp0(j))
#endif
             pmean2(j) = pmean2(j) + (np(id,j) * frac) * strength * wnum(l)**3 / (EXP(c2 * wnum(l) / temp2) - 1d0)
#ifdef DEBUG             
             t_count(3) = t_count(3) + (mpi_wtime() - t0)
             n_count(3) = n_count(3) + 1
#endif             
             ! ---------------------------------------
             
             ! --- LINE SORT BASED ON RATIO OF LINE CENTER TO CONTINUUM -----
#ifdef CROSS_CHECK
             iprof = 'bG'
#else             
#ifdef DEBUG             
             t0 = mpi_wtime()
#endif             
             k_wnum = MAX(find_grid(grd, wnum(l)), ks)
             delta = ((np(id,j) * frac) * strength / (sigma * SQRT(pi))) / (cnt(k_wnum,j) + TINY(0d0))
             IF(delta > delta_voigt) THEN ! VOIGT LINES
                iprof = 'Vo'
             ELSE ! GAUSS LINES
                iprof = 'Ga'
             END IF
#ifdef DEBUG             
             t_count(4) = t_count(4) + (mpi_wtime() - t0)
             n_count(4) = n_count(4) + 1
#endif
             ! --------------------------------------------------------------
#ifdef DEBUG             
             t0 = mpi_wtime()
#endif             
             IF(delta < delta_crit) THEN
#ifdef LINE_STATISTICS          
                slines(j,2) = slines(j,2) + 1
#endif
#ifdef DEBUG
                t2 = t2 + (mpi_wtime() - t0)
#endif             
                CYCLE
             ELSE
#ifdef DEBUG
                t2 = t2 + (mpi_wtime() - t0)
#endif             
             END IF
#endif             
             ! --------------------------------------------------------------
             
             ! --- SEARCH WINDOW IN WAVELENGTH -----
#ifdef DEBUG             
             t0 = mpi_wtime()
#endif
             IF(iprof == 'Ga') THEN
                cutoff_p = wnum(l) + cutoff0_gauss * sigma
                cutoff_m = wnum(l) - cutoff0_gauss * sigma
             ELSE
                cutoff_p = wnum(l) + cutoff0_voigt
                cutoff_m = wnum(l) - cutoff0_voigt
             END IF
#ifdef DEBUG             
             t_count(5) = t_count(5) + (mpi_wtime() - t0)
             n_count(5) = n_count(5) + 1
#endif             
             !--------------------------------------
             
             ! *** グリッド範囲にないものは無視する
#ifdef DEBUG             
             t0 = mpi_wtime()
#endif             
             IF(cutoff_m > grd(ke) .OR. cutoff_p < grd(ks)) THEN
#ifdef LINE_STATISTICS                
                slines(j,2) = slines(j,2) + 1
#endif                
#ifdef DEBUG
                t2 = t2 + (mpi_wtime() - t0)
#endif             
                CYCLE
             ELSE
#ifdef DEBUG
                t2 = t2 + (mpi_wtime() - t0)
#endif             
             END IF
             
             ! --- SEARCH WINDOW in GRID POINTS -----
#ifdef DEBUG             
             t0 = mpi_wtime()
#endif             
             kmin = find_grid(grd, cutoff_m) + 1
             kmax = find_grid(grd, cutoff_p)
#ifdef DEBUG             
             t_count(6) = t_count(6) + (mpi_wtime() - t0)
             n_count(6) = n_count(6) + 1
#endif             
             ! --------------------------------------

             ! *** ライン幅がグリッド幅に比べて狭すぎるものは無視する
#ifdef DEBUG             
             t0 = mpi_wtime()
#endif             
             IF(kmax < kmin) THEN
#ifdef LINE_STATISTICS                
                slines(j,2) = slines(j,2) + 1
#endif                
#ifdef DEBUG
                t2 = t2 + (mpi_wtime() - t0)
#endif             
                CYCLE
             ELSE
#ifdef DEBUG
                t2 = t2 + (mpi_wtime() - t0)
#endif             
             END IF
             
             ! ----------------------------------------   
#ifdef DEBUG             
             t0 = mpi_wtime()
#endif             
             ! Lorentz width [cm^-1]
             IF(iprof == 'Vo' .OR. iprof == 'Lo') THEN
                ! classical broadening
                gn = g_nat(wnum(l)) / (4d0 * pi * clight) ![cm^-1]
                gamma = gn
#ifndef CROSS_CHECK
                ! van der Waals broadening by species 1 and species 2
                ! HITRAN: 1 = air, 2 = self
                ! Exomol: 1 = H2 , 2 = He
                IF(id_1 == 0) THEN !HITRAN (AIR)
                   pres_1 = (np(TOTAL,j) - np(id_2,j)) * k_bol * temp0(j)
                ELSE ! EXOMOL(H2)
                   pres_1 = np(id_1,j) * k_bol * temp0(j)
                END IF
                pres_2 = np(id_2,j) * k_bol * temp0(j) ! HITRAN(SELF), EXOMOL(He)
                g6 = gm_1(l) * (temp_ref / temp0(j)) ** nx_1(l) * (pres_1 / pres_ref) + &
                     gm_2(l) * (temp_ref / temp0(j)) ** nx_2(l) * (pres_2 / pres_ref)
                gamma = gamma + g6 ![cm^-1]
#endif                
#ifdef LINE_STATISTICS
                slines(j,3) = slines(j,3) + 1
#endif                
             ELSE
                gamma = 0d0
#ifdef LINE_STATISTICS
                slines(j,4) = slines(j,4) + 1
#endif                
             END IF
#ifdef DEBUG             
             t_count(7) = t_count(7) + (mpi_wtime() - t0)
             n_count(7) = n_count(7) + 1
#endif             
             ! ------------------------------------------------------
             ! line profile [cm]
#ifdef DEBUG             
             t0 = mpi_wtime()
#endif             
             CALL line_prof(grd(kmin:kmax), dgrd(kmin:kmax), wnum(l), gamma, sigma, iprof, prof(kmin:kmax))
#ifdef DEBUG
             IF(iprof == 'Vo') THEN
                t_count( 8) = t_count( 8) + (mpi_wtime() - t0)
                n_count( 8) = n_count( 8) + (kmax - kmin + 1)
                n_count(12) = n_count(12) + 1
             ELSE
                t_count(11) = t_count(11) + (mpi_wtime() - t0)
                n_count(11) = n_count(11) + (kmax - kmin + 1)
                n_count(13) = n_count(13) + 1
             END IF
#endif             
             ! -----------------------------
             ! absorption coefficient [cm^-1] = number density [cm^-3] * cross section [cm^2]
#ifdef DEBUG             
             t0 = mpi_wtime()
#endif
             alp(kmin:kmax,j) = alp(kmin:kmax,j) + (np(id,j) * frac) * (strength * prof(kmin:kmax))

#ifdef DEBUG             
             t_count(9) = t_count(9) + (mpi_wtime() - t0)
             n_count(9) = n_count(9) + 1
#endif
          END DO LINE
       END DO LAYER
    END DO BLOCK
#ifdef DEBUG             
    t0 = mpi_wtime()
#endif
    DEALLOCATE(gfva, wnum, lowe, gm_1, nx_1, gm_2, nx_2)
    IF(nprc_jconst /= 1) THEN
       IF(myrk_jconst /= 0) THEN
          ALLOCATE(dummy(ke*je))
          CALL MPI_REDUCE(alp, dummy, ke*je, MPI_DOUBLE_PRECISION, MPI_SUM, 0, mpi_jconst_world, error)
          DEALLOCATE(dummy)
       ELSE
          CALL MPI_REDUCE(MPI_IN_PLACE, alp, ke*je, MPI_DOUBLE_PRECISION, MPI_SUM, 0, mpi_jconst_world, error)
       END IF
    END IF
    IF(kprc /= 1) THEN
       IF(myrk_line /= 0) THEN
          ALLOCATE(dummy(ke*je))
          CALL MPI_REDUCE(alp, dummy, ke*je, MPI_DOUBLE_PRECISION, MPI_SUM, 0, mpi_line_world, error)
          DEALLOCATE(dummy)
       ELSE
          CALL MPI_REDUCE(MPI_IN_PLACE, alp, ke*je, MPI_DOUBLE_PRECISION, MPI_SUM, 0, mpi_line_world, error)
       END IF
    END IF
    
    ! CLOSE location identifiers
    CALL h5Dclose_f(data_wnum, error)
    CALL h5Dclose_f(data_lowe, error)
    CALL h5Dclose_f(data_gfva, error)
    CALL h5Dclose_f(data_gm_1, error)
    CALL h5Dclose_f(data_nx_1, error)
    CALL h5Dclose_f(data_gm_2, error)
    CALL h5Dclose_f(data_nx_2, error)
999 CONTINUE
    CALL h5Gclose_f(grp_trans, error)
    CALL h5Fclose_f(file_spec, error)
    CALL h5close_f(error)
    DEALLOCATE(q, prof)

!!$#ifdef LINE_STATISTICS    
!!$    DO j = js, je
!!$       IF(slines(j,1) /= SUM(slines(j,2:4))) THEN
!!$          PRINT *, '***** LINE STATISTICS ERROR: ', TRIM(source), j, slines(j,1), slines(j,2), slines(j,3), slines(j,4)
!!$          STOP
!!$       END IF
!!$    END DO
!!$#endif    
#ifdef DEBUG    
    t_count(10) = t_count(10) + (mpi_wtime() - t0)
    n_count(10) = n_count(10) + 1
    t1 = (mpi_wtime() - t1)
    CALL record_count(TRIM(source), t1, t2, t_count, n_count, je)
#endif

    RETURN
  END SUBROUTINE line_molec
  

  SUBROUTINE line_kurucz(source, temp0, np, grd, dgrd, alp, cnt, pmean, temp2, pmean2, slines)
    USE code_module, ONLY : ELECTRON, HI, HeI, H2
    USE const_module, ONLY : c2, k_bol, clight, e2, m_ele, pi, ene_hyd, a_0, hbar, alpha, m_hyd, amu, h, cm_to_ang
    USE h5lx_module, ONLY : h5LXread_hyperslab_dataset_F64, h5LXread_hyperslab_dataset_F32, h5LXread_hyperslab_dataset_I32, &
         h5LXopen_file
    USE mpi_module, ONLY : myrk_line, mprc, myrk_m, lprc, myrk_l, mpi_line_world, para_range, &
         kprc, mpi_jconst_world, myrk_jconst, nprc_jconst
    IMPLICIT NONE

    !Peter Schwerdtfeger & Jeffrey K. Nagle (2018)
    REAL(REAL64), PARAMETER :: pol_HI = 4.507d0 * a_0**3 ![cm^3] polarizability of HI
    REAL(REAL64), PARAMETER :: pol_HeI = 1.38375d0 * a_0**3 ![cm^3] polarizability of HeI
    !Wilkins & Taylor 1968
    REAL(REAL64), PARAMETER :: pol_H2 = 0.8d-24 ![cm^3] polarizability of H2
    
    CHARACTER(*), INTENT(IN) :: source ! line source
    REAL(REAL64), INTENT(IN) :: temp0(:) ! temperature [K]
    REAL(REAL64), INTENT(IN) :: np(:,:) ! partial density
    REAL(REAL64), INTENT(IN) :: grd(:), dgrd(:) ! wavenumber grid [cm^{-1}]
    REAL(REAL64), INTENT(OUT) :: alp(:,:) ! opacity [cm^{-1}]
    REAL(REAL64), INTENT(IN) :: cnt(:,:) ! opacity [cm^{-1}]
    REAL(REAL64), INTENT(INOUT) :: pmean(:) ! planck_mean opacity
    REAL(REAL64), INTENT(IN) :: temp2 ! radiation temperature
    REAL(REAL64), INTENT(INOUT) :: pmean2(:) ! two-temp planck_mean opacity
    INTEGER(INT64), INTENT(INOUT), OPTIONAL :: slines(:,:)

    ! HDF5 variables
    INTEGER(HID_T) :: file_spec, grp_trans, grp_prop, grp_ion
    INTEGER(HID_T) :: data_wnum, data_lowe, data_gfva, data_code, data_gm_0, data_gm_1, data_gm_2
    INTEGER(HSIZE_T), DIMENSION(1) :: offset, block_size, dims
    
    ! grid index
    INTEGER(INT32) :: ks, ke, kmin, kmax
    ! lines
    INTEGER(INT64) :: l, nlines, n_rest
    INTEGER(INT32) :: nc, nc_max, ls, le
!    INTEGER(INT32), PARAMETER :: chunk0 = 128**2 *8
    ! main variables
    REAL(REAL64), ALLOCATABLE, DIMENSION(:) :: wnum, lowe
    REAL(REAL32), ALLOCATABLE, DIMENSION(:) :: gfva, gm_0, gm_1, gm_2
    INTEGER(INT32), ALLOCATABLE, DIMENSION(:) :: code
    ! Read from .prop
    INTEGER(INT32), ALLOCATABLE :: gtot(:)
    INTEGER(INT32) :: na, ni, nstates, indx, chunk0
    REAL(REAL64), ALLOCATABLE :: ene(:), q(:,:)
    INTEGER(INT32), PARAMETER :: na_uran = 92, mn_uran = 238
    REAL(REAL64), ALLOCATABLE :: mass(:), wn_ion(:)
    CHARACTER :: natm*3, nion*2, fname*128
    ! Lorentz profile
    ! constants
    REAL(REAL64) g6, g4, gn, gamma
    REAL(REAL64), PARAMETER :: temp_ref = 1d4
    ! Gaussian profile
    REAL(REAL64), ALLOCATABLE :: prof(:)
    REAL(REAL64) :: sigma
    REAL(REAL64) :: strength, strength0
    ! parameters read from files
    INTEGER(INT32) :: icutoff, k_wnum
    CHARACTER(LEN=2) :: iprof = 'Vo'
    REAL(REAL64) cutoff0_voigt, cutoff0_gauss, cutoff_p, cutoff_m, strength_cutoff0, wnum_crit, &
         delta_crit, delta_voigt, delta
    ! others
    INTEGER(INT32) :: error, j, js, je, prog, z, k
    REAL(REAL64) :: r2u, r2l
    REAL(REAL64), ALLOCATABLE :: dummy(:)
#ifndef UNUSE_GFGAM
    LOGICAL :: exists
#endif    
#ifdef PLANCK_MEAN    
    REAL(REAL64) :: a, ap
#endif
#ifdef DEBUG
    INTEGER(INT32), PARAMETER :: countmax = 13
    REAL(REAL64) :: t_count(0:countmax), t0, t1, t2
    INTEGER(INT64) :: n_count(0:countmax)
#endif
    NAMELIST /cutoffs/ icutoff, cutoff0_voigt, cutoff0_gauss, strength_cutoff0, wnum_crit, &
         delta_crit, delta_voigt
    
    ALLOCATE(mass((na_uran+1)*100))
    ALLOCATE(wn_ion((na_uran+1)*100))
#ifdef DEBUG
    t1 = mpi_wtime()
    t0 = mpi_wtime()
    t_count(:) = 0d0
    n_count(:) = 0
    t2 = 0d0
#endif
    alp = 0d0

    js = 1
    je = UBOUND(temp0,1)

    ks = 1
    ke = UBOUND(grd,1)
    
    ALLOCATE(prof(ks:ke))

    ! パラメータ読み込み
    OPEN(5, FILE='input/fort.5', STATUS='OLD')
    READ(5, NML=cutoffs)
    CLOSE(5)

    ! 読み込みファイルを開く
    CALL h5open_f(error)
    CALL h5Eset_auto_f(0,error)

    !-----------------------------------------------------------------------------
    ! 同位体置換体の質量・存在比・分配関数
    ALLOCATE(q((na_uran+1)*100,UBOUND(temp0,1)))
    fname = 'input/h5/NIST.h5'
    CALL h5Fopen_f(TRIM(fname), H5F_ACC_RDONLY_F, file_spec, error)
    IF(error /= 0) THEN
       PRINT *, '*** '//TRIM(fname)//' DOES NOT EXISTS ***'
       CALL MPI_FINALIZE(ERROR); STOP
    END IF
    CALL h5Gopen_f(file_spec, 'prop', grp_prop, error) 
    DO na = 1, na_uran
       WRITE(natm,'(I3.3)') na
       DO ni = 0, na - 1 ! loop ion
          indx = na * 100 + ni
          WRITE(nion,'(I2.2)') ni
          CALL h5Gopen_f(grp_prop, natm//'.'//nion, grp_ion, error)
          dims = [1]
          CALL h5LTread_dataset_f(grp_ion, 'mass'   , H5T_IEEE_F64LE, mass(indx) , dims, error)
          CALL h5LTread_dataset_f(grp_ion, 'nstates', H5T_STD_I32LE , nstates, dims, error)
          CALL h5LTread_dataset_f(grp_ion, 'eneion', H5T_IEEE_F64LE, wn_ion(indx), dims, error)
          ALLOCATE(ene(nstates), gtot(nstates))
          dims = [nstates]
          CALL h5LTread_dataset_f(grp_ion, 'ene' , H5T_IEEE_F64LE, ene , dims, error)
          CALL h5LTread_dataset_f(grp_ion, 'gtot', H5T_STD_I32LE, gtot, dims, error)
          DO j = js, je
             CALL compute_q(temp0(j), ene, gtot, q(indx,j))
          END DO
          DEALLOCATE(ene, gtot)
          CALL h5Gclose_f(grp_ion, error)
       END DO
    END DO
    CALL h5Gclose_f(grp_prop, error)
    CALL h5Fclose_f(file_spec, error)
#ifndef UNUSE_GFGAM
    fname = 'input/h5/gfgam.h5'
    CALL h5Fopen_f(TRIM(fname), H5F_ACC_RDONLY_F, file_spec, error)
    IF(error /= 0) THEN
       PRINT *, '*** '//TRIM(fname)//' DOES NOT EXISTS ***'
       CALL MPI_FINALIZE(ERROR); STOP
    END IF
    DO na = 1, na_uran
       WRITE(natm,'(I3.3)') na
       DO ni = 0, na - 1 ! loop ion
          indx = na * 100 + ni
          WRITE(nion,'(I2.2)') ni
          ! ENERGY AND GTOT
          CALL h5Lexists_f(file_spec, natm//nion, exists, error)
          IF(exists) THEN
             CALL h5Gopen_f(file_spec, natm//nion, grp_ion, error)
             dims = [1]
             CALL h5LTread_dataset_f(grp_ion, 'nstates', H5T_STD_I32LE , nstates, dims, error)
             ALLOCATE(ene(nstates), gtot(nstates))
             dims = [nstates]
             CALL h5LTread_dataset_f(grp_ion, 'ene' , H5T_IEEE_F64LE, ene , dims, error)
             CALL h5LTread_dataset_f(grp_ion, 'gtot', H5T_STD_I32LE, gtot, dims, error)
             DO j = js, je
                CALL compute_q(temp0(j), ene, gtot, q(indx,j))
             END DO
             DEALLOCATE(ene, gtot)
             CALL h5Gclose_f(grp_ion, error)
          END IF
       END DO
    END DO
    CALL h5Fclose_f(file_spec, error)
#endif    
    !-----------------------------------------------------------------------------

    ! HDF5ファイルを開く
    fname = 'input/h5/'//TRIM(source)//'.h5'
#ifndef MPIO    
    CALL h5Fopen_f(TRIM(fname), H5F_ACC_RDONLY_F, file_spec, error)
#else    
    CALL h5LXopen_file(TRIM(fname), H5F_ACC_RDONLY_F, file_spec, error, &
         MPI_COMM_WORLD, MPI_INFO_NULL)
#endif    
    IF(error /= 0) THEN
       PRINT *, '*** '//TRIM(fname)//' DOES NOT EXISTS ***'
       CALL MPI_FINALIZE(ERROR); STOP
    END IF
    !-----------------------------------------------------------------------------
    CALL h5Gopen_f(file_spec, 'trans', grp_trans, error)
    dims = [1]
    CALL h5LTread_dataset_f(grp_trans, 'nlines', H5T_STD_I64LE, nlines, dims, error)
    CALL h5LTread_dataset_f(grp_trans, 'chunk0', H5T_STD_I32LE, chunk0, dims, error)
    CALL h5Dopen_f(grp_trans, 'wnum', data_wnum, error)
    CALL h5Dopen_f(grp_trans, 'lowe', data_lowe, error)
    CALL h5Dopen_f(grp_trans, 'gfva', data_gfva, error)
    CALL h5Dopen_f(grp_trans, 'code', data_code, error)
    CALL h5Dopen_f(grp_trans, 'gm_0', data_gm_0, error)
    CALL h5Dopen_f(grp_trans, 'gm_1', data_gm_1, error)
    CALL h5Dopen_f(grp_trans, 'gm_2', data_gm_2, error)

!    if(myrk == 0) print *, 'nlines = ', nlines

    !-----------------------------------------------------------------------------

    block_size = [MIN(chunk0, nlines)]
    nc_max = INT(nlines / block_size(1))
    n_rest = nlines - nc_max * block_size(1)
    ALLOCATE(gfva(block_size(1)), wnum(block_size(1)), lowe(block_size(1)), code(block_size(1)), &
             gm_0(block_size(1)), gm_1(block_size(1)), gm_2(block_size(1)))
    CALL para_range(1, INT(block_size(1)), lprc, myrk_l, ls, le)
    offset = myrk_m * block_size

    prog = 0
#ifdef DEBUG
    t_count(0) = t_count(0) + (mpi_wtime() - t0)
    n_count(0) = n_count(0) + 1
#endif
    BLOCK: DO nc = 1 + myrk_m, nc_max + 1, mprc
       ! データブロックを読み込む -------
#ifdef DEBUG       
       t0 = mpi_wtime()
#endif       
       IF(nc == nc_max + 1) THEN
          IF(n_rest == 0) EXIT
          DEALLOCATE(gfva, wnum, lowe, code, gm_0, gm_1, gm_2)
          block_size = [n_rest]
          ALLOCATE(gfva(block_size(1)), wnum(block_size(1)), lowe(block_size(1)), code(block_size(1)), &
                   gm_0(block_size(1)), gm_1(block_size(1)), gm_2(block_size(1)))
          CALL para_range(1, INT(block_size(1)), lprc, myrk_l, ls, le)
       END IF

       CALL h5LXread_hyperslab_dataset_F64(data_wnum, wnum, 1, block_size, offset, block_size)
       CALL h5LXread_hyperslab_dataset_F64(data_lowe, lowe, 1, block_size, offset, block_size)
       CALL h5LXread_hyperslab_dataset_F32(data_gfva, gfva, 1, block_size, offset, block_size)
       CALL h5LXread_hyperslab_dataset_I32(data_code, code, 1, block_size, offset, block_size)
       CALL h5LXread_hyperslab_dataset_F32(data_gm_0, gm_0, 1, block_size, offset, block_size)
       CALL h5LXread_hyperslab_dataset_F32(data_gm_1, gm_1, 1, block_size, offset, block_size)
       CALL h5LXread_hyperslab_dataset_F32(data_gm_2, gm_2, 1, block_size, offset, block_size)
       offset = offset + block_size * mprc
#ifdef LINE_STATISTICS
       DO j = js, je
          slines(j,1) = slines(j,1) + block_size(1)
       END DO
#endif    
#ifdef DEBUG
       t_count(1) = t_count(1) + (mpi_wtime() - t0)
       n_count(1) = n_count(1) + 1
#endif       
       ! -------------------------------
       LAYER: DO j = js, je
          LINE: DO l = ls, le
             ! --------------------------------------------------------------
#ifdef DEBUG             
             t0 = mpi_wtime()
#endif             
             ! *** 分圧0の化学種のラインは無視する
             IF(np(code(l),j) == 0d0 .OR. wnum(l) == 0d0) THEN
#ifdef LINE_STATISTICS          
                slines(j,2) = slines(j,2) + 1
#endif
#ifdef DEBUG
                t2 = t2 + (mpi_wtime() - t0)
#endif             
                CYCLE
             ELSE
#ifdef DEBUG             
                t2 = t2 + (mpi_wtime() - t0)
#endif             
             END IF
             ! --------------------------------------------------------------

!#define SELECT_SINGLE_SPECIES 2600
#ifdef SELECT_SINGLE_SPECIES         
             IF(code(l) /= SELECT_SINGLE_SPECIES) THEN
#ifdef LINE_STATISTICS          
                slines(j,2) = slines(j,2) + 1
#endif
                CYCLE
             END IF
#endif
             ! --- COMPUTE LINE STRENGTH [cm] -----
#ifdef DEBUG
             t0 = mpi_wtime()
#endif             
             strength0 = (pi * e2 / (m_ele * clight**2)) * gfva(l) * EXP(-c2 * lowe(l) / temp0(j)) / q(code(l),j)
             strength = strength0 * (1d0 - EXP(-c2 * wnum(l) / temp0(j)))
#ifdef DEBUG
             t_count(2) = t_count(2) + (mpi_wtime() - t0)
             n_count(2) = n_count(2) + 1
#endif
             ! ------------------------------------
#ifdef STRENGTH_CUTOFF             
#ifndef CROSS_CHECK
             ! *** 強さが弱いラインは無視する
             IF(strength < strength_cutoff(wnum(l), temp0(j), wnum_crit, strength_cutoff0)) THEN
#ifdef LINE_STATISTICS          
                slines(j,2) = slines(j,2) + 1
#endif
                CYCLE
             END IF
#endif
#endif             
             ! --- COMPUTE DOPPLER WIDTH [cm^-1] -----
#ifdef DEBUG             
             t0 = mpi_wtime()
#endif             
             sigma = SQRT(2d0 * k_bol * temp0(j) / (mass(code(l)) * amu)) / clight * wnum(l)

             ! --- COMPUTE PLANCK-MEAN OPACITY -----
#ifdef PLANCK_MEAN
             a = wnum(l) - c2 * sigma**2 / temp0(j)
             ap = a / sigma
             pmean(j) = pmean(j) + np(code(l),j) * strength * &
                  EXP(-c2 / temp0(j) * (wnum(l) + a) / 2d0) / SQRT(pi) / sigma * &
                  0.25d0 * sigma**4 * (SQRT(pi) * ap * (2d0 * ap**2 + 3d0) * (ERF(ap) + 1d0) + &
                  2d0 * (ap**2 + 1d0) * EXP(-ap**2))
#else             
             pmean(j) = pmean(j) + np(code(l),j) * strength0 * wnum(l)**3 * EXP(-c2 * wnum(l) / temp0(j))
#endif
             pmean2(j) = pmean2(j) + np(code(l),j) * strength * wnum(l)**3 / (EXP(c2 * wnum(l) / temp2) - 1d0)
#ifdef DEBUG             
             t_count(3) = t_count(3) + (mpi_wtime() - t0)
             n_count(3) = n_count(3) + 1
#endif             
             ! --------------------------------------

             ! --- LINE SORT BASED ON RATIO OF LINE CENTER TO CONTINUUM -----
#ifdef CROSS_CHECK
             iprof = 'bG'
#else        
#ifdef DEBUG             
             t0 = mpi_wtime()
#endif
             k_wnum = MAX(find_grid(grd, wnum(l)), ks)
             delta = (np(code(l),j) * strength / (sigma * SQRT(pi))) / (cnt(k_wnum,j) + TINY(0d0))
             IF(delta > delta_voigt) THEN ! VOIGT LINES
                iprof = 'Vo'
             ELSE ! GAUSS LINES
                iprof = 'Ga'
             END IF
#ifdef DEBUG             
             t_count(4) = t_count(4) + (mpi_wtime() - t0)
             n_count(4) = n_count(4) + 1
#endif
             ! --------------------------------------------------------------
#ifdef DEBUG             
             t0 = mpi_wtime()
#endif             
             IF(delta < delta_crit) THEN
#ifdef LINE_STATISTICS          
                slines(j,2) = slines(j,2) + 1
#endif
#ifdef DEBUG             
                t2 = t2 + (mpi_wtime() - t0)
#endif             
                CYCLE
             ELSE
#ifdef DEBUG             
                t2 = t2 + (mpi_wtime() - t0)
#endif             
             END IF
#endif       
             ! --------------------------------------------------------------

             ! --- SEARCH WINDOW IN WAVELENGTH -----
#ifdef DEBUG             
             t0 = mpi_wtime()
#endif             
#ifndef WAVELENGTH_BASE
             IF(iprof == 'Ga') THEN
                cutoff_p = wnum(l) + cutoff0_gauss * sigma
                cutoff_m = wnum(l) - cutoff0_gauss * sigma
             ELSE
                cutoff_p = wnum(l) + cutoff0_voigt
                cutoff_m = wnum(l) - cutoff0_voigt
             END IF
#else             
             IF(iprof == 'Ga') THEN
                cutoff_m = 1d0 / (1d0 / wnum(l) + (cutoff0_gauss * sigma / wnum(l)**2))
                cutoff_p = 1d0 / (1d0 / wnum(l) - (cutoff0_gauss * sigma / wnum(l)**2))
                IF(cutoff_p <= 0d0) cutoff_p = grd(ke)
             ELSE
                cutoff_m = 1d0 / (1d0 / wnum(l) + (cutoff0_voigt / cm_to_ang))
                cutoff_p = 1d0 / (1d0 / wnum(l) - (cutoff0_voigt / cm_to_ang))
                IF(cutoff_p <= 0d0) cutoff_p = grd(ke)
             END IF
#endif             
#ifdef DEBUG             
             t_count(5) = t_count(5) + (mpi_wtime() - t0)
             n_count(5) = n_count(5) + 1
#endif             

             !--------------------------------------

             ! *** グリッド範囲にないものは無視する
#ifdef DEBUG             
             t0 = mpi_wtime()
#endif             
             IF(cutoff_m > grd(ke) .OR. cutoff_p < grd(ks)) THEN
#ifdef LINE_STATISTICS          
                slines(j,2) = slines(j,2) + 1
#endif
#ifdef DEBUG             
                t2 = t2 + (mpi_wtime() - t0)
#endif             
                CYCLE
             ELSE
#ifdef DEBUG             
                t2 = t2 + (mpi_wtime() - t0)
#endif             
             END IF

             ! --- SEARCH WINDOW in GRID POINTS -----
#ifdef DEBUG             
             t0 = mpi_wtime()
#endif             
             kmin = find_grid(grd, cutoff_m) + 1
             kmax = find_grid(grd, cutoff_p)
#ifdef DEBUG             
             t_count(6) = t_count(6) + (mpi_wtime() - t0)
             n_count(6) = n_count(6) + 1
#endif             
             ! --------------------------------------
             ! *** ライン幅がグリッド幅に比べて狭すぎるものは無視する
#ifdef DEBUG             
             t0 = mpi_wtime()
#endif             
             IF(kmax < kmin) THEN
#ifdef LINE_STATISTICS          
                slines(j,2) = slines(j,2) + 1
#endif
#ifdef DEBUG             
                t2 = t2 + (mpi_wtime() - t0)
#endif             
                CYCLE
             ELSE
#ifdef DEBUG             
                t2 = t2 + (mpi_wtime() - t0)
#endif             
             END IF

             ! ---- Lorentz width [s^-1] <= WATCH OUT FOR UNIT! -----
#ifdef DEBUG             
             t0 = mpi_wtime()
#endif
             IF(iprof == 'Vo' .OR. iprof == 'Lo') THEN
                ! natural broadening
                gn = gm_0(l)
                IF(gn == 0d0) THEN ! classical broadening
                   gn = g_nat(wnum(l))
                END IF
                gamma =gn / (4d0 * pi * clight)
#ifndef CROSS_CHECK
                ! radiator radius (divided by a_0) squared 
                z = (code(l) - (code(l) / 100) * 100) + 1
                r2u = rsqd(z, wn_ion(code(l)), lowe(l)+wnum(l))
                r2l = rsqd(z, wn_ion(code(l)), lowe(l)        )
                ! quadratic Stark broadening by electron
                g4 = gm_1(l) * np(ELECTRON,j) ![s^-1]
                g4 = g4 * (temp0(j) / temp_ref)**(-0.5d0)
                ! employ Griem if not provided by Kurucz and <R^2> > 0 (exclude autoionization lines)
                IF(g4 == 0d0 .or. code(l) == HI) THEN
#ifndef PHOENIX_BROADENING
                   g4 = g_stark(r2u, r2l, temp0(j), np(ELECTRON,j))
#else                   
                   g4 = g_stark_phoenix(lowe(l)+wnum(l), code(l), wn_ion(code(l)))
#endif                   
                END IF
                ! van der Waals broadening by hydrogen atom
                g6 = gm_2(l)
                g6 = g6 * (temp0(j) / temp_ref)**0.3d0
                ! employ C6 if not provided by Kurucz and <R^2> > 0 (exclude autoionization lines),
                IF(g6 == 0d0) THEN
#ifndef PHOENIX_BROADENING
                   g6 = g_vdw(r2u, r2l, temp0(j), pol_HI, mass(HI)*amu, mass(code(l))*amu)
#else                   
                   g6 = g_vdw_phoenix(lowe(l), wnum(l), code(l), wn_ion(code(l)), temp0(j), mass(HI), mass(code(l)))
#endif                   
                END IF
                ! perturbers: HI, HeI, and H2
                ! Note that atomic size dependence in C6 is ignored.
#ifndef PHOENIX_BROADENING
                g6 = g6 * (np(HI,j) + &
                     (pol_HeI / pol_HI)**0.4d0 * (mass(HeI)     / mass(HI))**0.3d0 * np(HeI,j) + &
                     (pol_H2  / pol_HI)**0.4d0 * (mass(HI)*2d0  / mass(HI))**0.3d0 * np(H2 ,j))
#else                
                g6 = g6 * (np(HI,j) + 0.42d0*np(HeI,j) + 0.85d0*np(H2,j))
#endif                
                ! total Lorentz width [cm^-1]
                gamma = gamma + (g4 + g6) / (4d0 * pi * clight)
#endif                
#ifdef LINE_STATISTICS
                slines(j,3) = slines(j,3) + 1
#endif                
             ELSE
                gamma = 0d0
#ifdef LINE_STATISTICS
                slines(j,4) = slines(j,4) + 1
#endif                
             END IF
#ifdef DEBUG             
             t_count(7) = t_count(7) + (mpi_wtime() - t0)
             n_count(7) = n_count(7) + 1
#endif
             ! ------------------------------------------------------

             ! ----- LINE PROFILE [cm] -----
#ifdef DEBUG             
             t0 = mpi_wtime()
#endif
             DO k = kmin, kmax
                CALL line_prof(grd(k:k), dgrd(k:k), wnum(l), gamma, sigma, iprof, prof=prof(k:k))
             END DO
!             CALL line_prof(grd(kmin:kmax), dgrd(kmin:kmax), wnum(l), gamma, sigma, iprof, prof=prof(kmin:kmax))
#ifdef DEBUG
             IF(iprof == 'Vo') THEN
                t_count( 8) = t_count( 8) + (mpi_wtime() - t0)
                n_count( 8) = n_count( 8) + (kmax - kmin + 1)
                n_count(12) = n_count(12) + 1
             ELSE
                t_count(11) = t_count(11) + (mpi_wtime() - t0)
                n_count(11) = n_count(11) + (kmax - kmin + 1)
                n_count(13) = n_count(13) + 1
             END IF
#endif             
             ! -----------------------------
             ! ----- ABSORPTION COEFFICIENT [CM^-1] = number density [cm^-3] * cross section [cm^2] -----
#ifdef DEBUG             
             t0 = mpi_wtime()
#endif
             alp(kmin:kmax,j) = alp(kmin:kmax,j) + np(code(l),j) * (strength * prof(kmin:kmax))
#ifdef DEBUG             
             t_count(9) = t_count(9) + (mpi_wtime() - t0)
             n_count(9) = n_count(9) + 1
#endif             
          END DO LINE
       END DO LAYER
    END DO BLOCK
#ifdef DEBUG
    t0 = mpi_wtime()
#endif
    DEALLOCATE(gfva, wnum, lowe, code, gm_0, gm_1, gm_2)
    IF(nprc_jconst /= 1) THEN
       IF(myrk_jconst /= 0) THEN
          ALLOCATE(dummy(ke*je))
          CALL MPI_REDUCE(alp, dummy, ke*je, MPI_DOUBLE_PRECISION, MPI_SUM, 0, mpi_jconst_world, error)
          DEALLOCATE(dummy)
       ELSE
          CALL MPI_REDUCE(MPI_IN_PLACE, alp, ke*je, MPI_DOUBLE_PRECISION, MPI_SUM, 0, mpi_jconst_world, error)
       END IF
    END IF
    IF(kprc /= 1) THEN
       IF(myrk_line /= 0) THEN
          ALLOCATE(dummy(ke*je))
          CALL MPI_REDUCE(alp, dummy, ke*je, MPI_DOUBLE_PRECISION, MPI_SUM, 0, mpi_line_world, error)
          DEALLOCATE(dummy)
       ELSE
          CALL MPI_REDUCE(MPI_IN_PLACE, alp, ke*je, MPI_DOUBLE_PRECISION, MPI_SUM, 0, mpi_line_world, error)
       END IF
    END IF
   
    ! CLOSE location identifiers
    CALL h5Dclose_f(data_wnum, error)
    CALL h5Dclose_f(data_lowe, error)
    CALL h5Dclose_f(data_gfva, error)
    CALL h5Dclose_f(data_code, error)
    CALL h5Dclose_f(data_gm_0, error)
    CALL h5Dclose_f(data_gm_1, error)
    CALL h5Dclose_f(data_gm_2, error)
    
    CALL h5Gclose_f(grp_trans, error)
    CALL h5Fclose_f(file_spec, error)
    CALL h5close_f(error)
    
    DEALLOCATE(q, mass, prof)
#ifdef DEBUG
    t_count(10) = t_count(10) + (mpi_wtime() - t0)
    n_count(10) = n_count(10) + 1
    t1 = (mpi_wtime() - t1)
    CALL record_count(TRIM(source), t1, t2, t_count, n_count, je)
#endif

    RETURN
  END SUBROUTINE line_kurucz

#ifdef DEBUG
  SUBROUTINE record_count(source, t_tot, t_cyc, t, n, jmax)
    USE mpi_module, ONLY : unit
    CHARACTER(LEN=*), INTENT(IN) :: source
    REAL(REAL64), INTENT(IN) :: t_tot
    REAL(REAL64), INTENT(IN) :: t_cyc
    REAL(REAL64), INTENT(IN) :: t(0:)
    INTEGER(INT64), INTENT(IN) :: n(0:)
    INTEGER(INT32), INTENT(IN) :: jmax
    INTEGER(INT32) :: u    
    CHARACTER(LEN=4) :: id

    u = unit + 3000

    WRITE(id,'(I4)') unit
    OPEN(u,FILE='output/debug_'//id//'.data',STATUS='UNKNOWN',POSITION='APPEND')
    WRITE(u,'(A10,X,A,F8.2,X,I3)') source, ': TOTAL', t_tot, jmax
    WRITE(u,'(A10,X,A,F8.2)') source, ': CYCLE', t_cyc 
    WRITE(u,'(A10,X,A,F8.2,X,I13)') source, ': init ', t( 0), n( 0)
    WRITE(u,'(A10,X,A,F8.2,X,I13)') source, ': readd', t( 1), n( 1)
    WRITE(u,'(A10,X,A,F8.2,X,I13)') source, ': stren', t( 2), n( 2)
    WRITE(u,'(A10,X,A,F8.2,X,I13)') source, ': sigma', t( 3), n( 3)
    WRITE(u,'(A10,X,A,F8.2,X,I13)') source, ': delta', t( 4), n( 4)
    WRITE(u,'(A10,X,A,F8.2,X,I13)') source, ': cutof', t( 5), n( 5)
    WRITE(u,'(A10,X,A,F8.2,X,I13)') source, ': ksear', t( 6), n( 6)
    WRITE(u,'(A10,X,A,F8.2,X,I13)') source, ': gamma', t( 7), n( 7)
    WRITE(u,'(A10,X,A,F8.2,X,I13,X,E12.4,X,F8.2,X,I13)') &
         source, ': profv', t( 8), n( 8), t( 8)/DBLE(n( 8)), DBLE(n( 8))/DBLE(n(12)), n(12)
    WRITE(u,'(A10,X,A,F8.2,X,I13,X,E12.4,X,F8.2,X,I13)') &
         source, ': profg', t(11), n(11), t(11)/DBLE(n(11)), DBLE(n(11))/DBLE(n(13)), n(13)
    WRITE(u,'(A10,X,A,F8.2,X,I13)') source, ': alpha', t( 9), n( 9)
    WRITE(u,'(A10,X,A,F8.2,X,I13)') source, ': final', t(10), n(10)
    WRITE(u,'(A10,X,A,F8.2,X,I13)') source, ': total', SUM(t(0:10)), SUM(n(0:10))
    CLOSE(u)

    RETURN    
  END SUBROUTINE record_count
#endif

  REAL(REAL64) FUNCTION rsqd(z, ion, ene)
    USE const_module, ONLY : ene_hyd, h, c
    INTEGER(INT32) :: z
    REAL(REAL64) :: ion, ene
    REAL(REAL64) :: nsqd
    nsqd = MAX(DBLE(z)**2 * ene_hyd/(h*c) / (ion - ene), 0d0)
    rsqd = 2.5d0 * nsqd**2 / z**2
    RETURN
  END FUNCTION rsqd


  SUBROUTINE line_kurucz_p(source, temp0, np, grd, dgrd, alp, cnt, pmean, temp2, pmean2, slines)
    ! compute line absorption coefficient from Kurucz database, following Phoenix convension
    USE code_module, ONLY : ELECTRON, HI, HeI, H2
    USE const_module, ONLY : c2, k_bol, clight, e2, m_ele, pi, ene_hyd, a_0, hbar, alpha, m_hyd, amu, h, cm_to_ang
    USE h5lx_module, ONLY : h5LXread_hyperslab_dataset_F64, h5LXread_hyperslab_dataset_F32, &
         h5LXread_hyperslab_dataset_I16, h5LXopen_file
    USE mpi_module, ONLY : myrk_line, mprc, myrk_m, lprc, myrk_l, mpi_line_world, para_range, kprc, myrk
    IMPLICIT NONE

    CHARACTER(*), INTENT(IN) :: source ! line source
    REAL(REAL64), INTENT(IN) :: temp0(:) ! temperature [K]
    REAL(REAL64), INTENT(IN) :: np(:,:) ! partial density
    REAL(REAL64), INTENT(IN) :: grd(:), dgrd(:) ! wavenumber grid [cm^{-1}]
    REAL(REAL64), INTENT(OUT) :: alp(:,:) ! opacity [cm^{-1}]
    REAL(REAL64), INTENT(IN) :: cnt(:,:) ! opacity [cm^{-1}]
    REAL(REAL64), INTENT(INOUT) :: pmean(:) ! planck_mean opacity
    REAL(REAL64), INTENT(IN) :: temp2 ! radiation temperature
    REAL(REAL64), INTENT(INOUT) :: pmean2(:) ! two-temp planck_mean opacity
    INTEGER(INT64), INTENT(INOUT), OPTIONAL :: slines(:,:)

    !Peter Schwerdtfeger & Jeffrey K. Nagle (2018)
    REAL(REAL64), PARAMETER :: pol_HI = 4.507d0 * a_0**3 ![cm^3] polarizability of HI
    REAL(REAL64), PARAMETER :: pol_HeI = 1.38375d0 * a_0**3 ![cm^3] polarizability of HeI
    !Wilkins & Taylor 1968
    REAL(REAL64), PARAMETER :: pol_H2 = 0.8d-24 ![cm^3] polarizability of H2
    
    ! HDF5 variables
    INTEGER(HID_T) :: file_spec, grp_trans, grp_prop, grp_ion
    INTEGER(HID_T) :: data_wnum, data_lowe, data_gfva, data_code, data_gm_0, data_gm_1, data_gm_2
    INTEGER(HSIZE_T), DIMENSION(1) :: offset, block_size, dims
    
    ! grid index
    INTEGER(INT32) :: ks, ke, kmin, kmax
    ! lines
    INTEGER(INT64) :: l, nlines, n_rest
    INTEGER(INT32) :: nc, nc_max, ls, le
    INTEGER(INT32) :: chunk0! = 128**2 *8
    ! main variables
    REAL(REAL64), ALLOCATABLE, DIMENSION(:) :: wnum, lowe
    REAL(REAL32), ALLOCATABLE, DIMENSION(:) :: gfva, gm_0, gm_1, gm_2
    INTEGER(INT16), ALLOCATABLE, DIMENSION(:) :: code
    ! Read from .prop
    INTEGER(INT32), ALLOCATABLE :: gtot(:)
    INTEGER(INT32) :: na, ni, nstates, indx
    REAL(REAL64), ALLOCATABLE :: ene(:), q(:,:)
    INTEGER(INT32), PARAMETER :: na_uran = 92, mn_uran = 238
    REAL(REAL64), ALLOCATABLE :: mass(:), wn_ion(:)
    CHARACTER :: natm*3, nion*2, fname*128
    ! Lorentz profile
    ! constants
    REAL(REAL64), PARAMETER :: temp_ref = 1d4
    REAL(REAL64) g6, g4, gn, gamma
    ! Gaussian profile
    REAL(REAL64), ALLOCATABLE :: prof(:)
    REAL(REAL64) :: sigma, sigma_h
    REAL(REAL64) :: strength, strength0
    ! parameters read from files
    INTEGER(INT32) :: icutoff, k_wnum
    CHARACTER(LEN=2) :: iprof = 'Vo'
    REAL(REAL64) cutoff0_voigt, cutoff0_gauss, cutoff_p, cutoff_m, strength_cutoff0, wnum_crit, &
         delta_crit, delta_voigt, delta
    ! others
    INTEGER(INT32) :: error, j, js, je, prog, code0
    REAL(REAL64), ALLOCATABLE :: dummy(:)
!    INTEGER(INT32) :: k
    NAMELIST /cutoffs/ icutoff, cutoff0_voigt, cutoff0_gauss, strength_cutoff0, wnum_crit, &
         delta_crit, delta_voigt

    ALLOCATE(mass((na_uran+1)*100), wn_ion((na_uran+1)*100))

    ! パラメータ読み込み
    OPEN(5, FILE='input/fort.5', STATUS='OLD')
    READ(5, NML=cutoffs)
    CLOSE(5)

    js = 1
    je = UBOUND(temp0,1)

    ks = 1
    ke = UBOUND(grd,1)
    
    ALLOCATE(prof(ks:ke))

    ! 読み込みファイルを開く
    CALL h5open_f(error)
    CALL h5Eset_auto_f(0,error)

    !-----------------------------------------------------------------------------
    ! 同位体置換体の質量・存在比・分配関数
    ALLOCATE(q((na_uran+1)*100,UBOUND(temp0,1)))
    fname = 'input/h5/NIST.h5'
    CALL h5Fopen_f(TRIM(fname), H5F_ACC_RDONLY_F, file_spec, error)
    CALL h5Gopen_f(file_spec, 'prop', grp_prop, error)
    IF(error /= 0) THEN
       PRINT *, '*** '//TRIM(fname)//' DOES NOT EXISTS ***'
       CALL MPI_FINALIZE(ERROR); STOP
    END IF
    DO na = 1, na_uran
       WRITE(natm,'(I3.3)') na
       DO ni = 0, na - 1 ! loop ion
          indx = na * 100 + ni
          WRITE(nion,'(I2.2)') ni
          CALL h5Gopen_f(grp_prop, natm//'.'//nion, grp_ion, error)
          dims = [1]
          CALL h5LTread_dataset_f(grp_ion, 'mass'   , H5T_IEEE_F64LE, mass(indx) , dims, error)
          CALL h5LTread_dataset_f(grp_ion, 'nstates', H5T_STD_I32LE , nstates, dims, error)
          CALL h5LTread_dataset_f(grp_ion, 'eneion', H5T_IEEE_F64LE, wn_ion(indx), dims, error)
          ALLOCATE(ene(nstates), gtot(nstates))
          dims = [nstates]
          CALL h5LTread_dataset_f(grp_ion, 'ene' , H5T_IEEE_F64LE, ene , dims, error)
          CALL h5LTread_dataset_f(grp_ion, 'gtot', H5T_STD_I32LE, gtot, dims, error)
          DO j = js, je
             CALL compute_q(temp0(j), ene, gtot, q(indx,j))
          END DO
          DEALLOCATE(ene, gtot)
          CALL h5Gclose_f(grp_ion, error)
       END DO
    END DO
    CALL h5Gclose_f(grp_prop, error)
    CALL h5Fclose_f(file_spec, error)
    !-----------------------------------------------------------------------------

    ! HDF5ファイルを開く
    CALL h5LXopen_file('input/h5/'//TRIM(source)//'.h5', H5F_ACC_RDONLY_F, file_spec, error, &
         MPI_COMM_WORLD, MPI_INFO_NULL)
    CALL h5Gopen_f(file_spec, 'ATOMIC_LINES', grp_trans, error)
    dims = [1]
    CALL h5LTread_dataset_f(grp_trans, 'nlines', H5T_STD_I64LE, nlines, dims, error)
    !if(myrk == 0) print *, 'nlines = ', nlines

    CALL h5LTread_dataset_f(grp_trans, 'chunky', H5T_STD_I32LE, chunk0, dims, error)
    CALL h5Dopen_f(grp_trans, 'wl', data_wnum, error)
    CALL h5Dopen_f(grp_trans, 'chi', data_lowe, error)
    CALL h5Dopen_f(grp_trans, 'log_gf', data_gfva, error)
    CALL h5Dopen_f(grp_trans, 'phx_code', data_code, error)
    CALL h5Dopen_f(grp_trans, 'gam_rad', data_gm_0, error)
    CALL h5Dopen_f(grp_trans, 'gam_stark', data_gm_1, error)
    CALL h5Dopen_f(grp_trans, 'gam_vdW', data_gm_2, error)
    !-----------------------------------------------------------------------------

    block_size = [MIN(chunk0, nlines)]
    nc_max = INT(nlines / block_size(1), INT32)
    n_rest = nlines - nc_max * block_size(1)
    ALLOCATE(gfva(block_size(1)), wnum(block_size(1)), lowe(block_size(1)), code(block_size(1)), &
             gm_0(block_size(1)), gm_1(block_size(1)), gm_2(block_size(1)))
    CALL para_range(1, INT(block_size(1)), lprc, myrk_l, ls, le)
    offset = myrk_m * block_size

    prog = 0
    BLOCK: DO nc = 1 + myrk_m, nc_max + 1, mprc
       ! データブロックを読み込む -------
       IF(nc == nc_max + 1) THEN
          IF(n_rest == 0) EXIT
          DEALLOCATE(gfva, wnum, lowe, code, gm_0, gm_1, gm_2)
          block_size = [n_rest]
          ALLOCATE(gfva(block_size(1)), wnum(block_size(1)), lowe(block_size(1)), code(block_size(1)), &
               gm_0(block_size(1)), gm_1(block_size(1)), gm_2(block_size(1)))
          CALL para_range(1, INT(block_size(1)), lprc, myrk_l, ls, le)
       END IF

       CALL h5LXread_hyperslab_dataset_F64(data_wnum, wnum, 1, block_size, offset, block_size)
       wnum = 1d8 / wnum ! convert to lambda in A to wavenumber in cm^-1
       CALL h5LXread_hyperslab_dataset_F64(data_lowe, lowe, 1, block_size, offset, block_size)
       CALL h5LXread_hyperslab_dataset_F32(data_gfva, gfva, 1, block_size, offset, block_size)
       gfva = 10e0**gfva
       CALL h5LXread_hyperslab_dataset_I16(data_code, code, 1, block_size, offset, block_size)
       CALL h5LXread_hyperslab_dataset_F32(data_gm_0, gm_0, 1, block_size, offset, block_size)
       CALL h5LXread_hyperslab_dataset_F32(data_gm_1, gm_1, 1, block_size, offset, block_size)
       CALL h5LXread_hyperslab_dataset_F32(data_gm_2, gm_2, 1, block_size, offset, block_size)
       offset = offset + block_size * mprc
       ! -------------------------------

#ifdef LINE_STATISTICS
       DO j = js, je
          slines(j,1) = slines(j,1) + block_size(1)
       END DO
#endif    
       LAYER: DO j = js, je
          LINE: DO l = ls, le
             code0 = map(code(l))

             ! *** 分圧0の化学種のラインは無視する
             IF(np(code0,j) == 0d0 .OR. wnum(l) == 0d0) THEN
#ifdef LINE_STATISTICS          
                slines(j,2) = slines(j,2) + 1
#endif
                CYCLE
             END IF

             ! --- COMPUTE LINE STRENGTH [cm] -----
             strength0 = (pi * e2 / (m_ele * clight**2)) * gfva(l) * EXP(-c2 * lowe(l) / temp0(j)) / q(code0, j)
             strength = strength0 * (1d0 - EXP(-c2 * wnum(l) / temp0(j)))
             ! ------------------------------------
             ! --- COMPUTE DOPPLER WIDTH [cm^-1] -----
             sigma = SQRT(2d0 * k_bol * temp0(j) / (mass(code0) * amu)) / clight * wnum(l)
             sigma_h = SQRT(2d0 * k_bol * temp0(j) / (mass(HI) * amu)) / clight * wnum(l)
             ! ---------------------------------------

#ifdef PLANCK_MEAN
             a = wnum(l) - c2 * sigma**2 / temp0(j)
             ap = a / sigma
             pmean(j) = pmean(j) + np(code0,j) * strength * &
                  EXP(-c2 / temp0(j) * (wnum(l) + a) / 2d0) / SQRT(pi) / sigma * &
                  0.25d0 * sigma**4 * (SQRT(pi) * ap * (2d0 * ap**2 + 3d0) * (ERF(ap) + 1d0) + &
                  2d0 * (ap**2 + 1d0) * EXP(-ap**2))
#else             
             pmean(j) = pmean(j) + np(code0,j) * strength0 * wnum(l)**3 * EXP(-c2 * wnum(l) / temp0(j))
#endif
             pmean2(j) = pmean2(j) + np(code0,j) * strength * wnum(l)**3 / (EXP(c2 * wnum(l) / temp2) - 1d0)

             ! --- LINE SORT BASED ON RATIO OF LINE CENTER TO CONTINUUM -----
#ifdef CROSS_CHECK
             iprof = 'bG'
#else             
             k_wnum = MAX(find_grid(grd, wnum(l)), ks)
             delta = (np(code0,j) * strength / (sigma_h * SQRT(pi))) / MAX(cnt(k_wnum,j), TINY(0d0)) &
                  / (1d0 - EXP(-c2 * wnum(l) / temp0(j))) !REMOVE STIMULATED EMISSION FACTOR
             IF(delta < delta_crit) THEN
#ifdef LINE_STATISTICS          
                slines(j,2) = slines(j,2) + 1
#endif
                CYCLE
             ELSE IF(delta > delta_voigt) THEN ! VOIGT LINES
                iprof = 'Vo'
             ELSE
                iprof = 'Ga'
             END IF
#endif
             ! --------------------------------------------------------------
             ! --- SEARCH WINDOW IN WAVELENGTH -----
             IF(iprof == 'Ga') THEN
                cutoff_m = 1d0 / (1d0 / wnum(l) + (cutoff0_gauss * sigma / wnum(l)**2))
                cutoff_p = 1d0 / (1d0 / wnum(l) - (cutoff0_gauss * sigma / wnum(l)**2))
                IF(cutoff_p <= 0d0) cutoff_p = grd(ke)
             ELSE
                cutoff_m = 1d0 / (1d0 / wnum(l) + (cutoff0_voigt / cm_to_ang))
                cutoff_p = 1d0 / (1d0 / wnum(l) - (cutoff0_voigt / cm_to_ang))
                IF(cutoff_p <= 0d0) cutoff_p = grd(ke)
             END IF

             !--------------------------------------

             ! *** グリッド範囲にないものは無視する
             IF(cutoff_m > grd(ke) .OR. cutoff_p < grd(ks)) THEN
#ifdef LINE_STATISTICS          
                slines(j,2) = slines(j,2) + 1
#endif
                CYCLE
             END IF

             ! --- SEARCH WINDOW in GRID POINTS -----
             kmin = find_grid(grd, cutoff_m) + 1
             kmax = find_grid(grd, cutoff_p)
             ! --------------------------------------

             ! *** ライン幅がグリッド幅に比べて狭すぎるものは無視する
             IF(kmax < kmin) THEN
#ifdef LINE_STATISTICS          
                slines(j,2) = slines(j,2) + 1
#endif
                CYCLE
             END IF

             ! ---- Lorentz width [s^-1] <= WATCH OUT FOR UNIT! -----
             IF(iprof == 'Vo' .OR. iprof == 'Lo') THEN
             ! natural broadening
                gn = gm_0(l)
                IF(gn <= 0d0) THEN ! classical broadening
                   gn = g_nat(wnum(l))
                END IF
                gamma = (gn) / (2d0 * pi * clight) / 2d0
#ifndef CROSS_CHECK             
                ! quadratic Stark broadening by electron
                ! employ Griem if not provided by Kurucz and <R^2> > 0 (exclude autoionization lines)
                g4 = gm_1(l)
                IF(g4 <= 0d0 .OR. code0 == HI) THEN
                   g4 = g_stark_phoenix(lowe(l)+wnum(l), code0, wn_ion(code0))
                END IF
                g4 = g4 * np(ELECTRON,j)
                
                ! van der Waals broadening by hydrogen atom
                g6 = gm_2(l)
                g6 = g6 * (temp0(j) / temp_ref)**0.3d0 ! SCALE TO temp_ref
                ! employ C6 if not provided by Kurucz and <R^2> > 0 (exclude autoionization lines)
                IF(g6 <= 0d0) THEN
                   g6 = g_vdw_phoenix(lowe(l), wnum(l), code0, wn_ion(code0), temp0(j), mass(HI), mass(code0))
                END IF
                !                g6 =g6 * np(HI,j)
                g6 = g6 * (np(HI,j) + &
                     (pol_HeI / pol_HI)**0.4d0 * (mass(HeI)     / mass(HI))**0.3d0 * np(HeI,j) + &
                     (pol_H2  / pol_HI)**0.4d0 * (mass(HI)*2d0  / mass(HI))**0.3d0 * np(H2 ,j))
          ! total Lorentz width [cm^-1]
                ! / (2 pi c): angular frequncy -> wavenumber
                ! / 2: FWHM -> HWHM
                gamma = gamma + (g4 + g6) / (2d0 * pi * clight) / 2d0
#endif
#ifdef LINE_STATISTICS
                slines(j,3) = slines(j,3) + 1
#endif                
             ELSE
                gamma = 0d0
#ifdef LINE_STATISTICS
                slines(j,4) = slines(j,4) + 1
#endif                
             END IF
             ! ------------------------------------------------------

             ! ----- LINE PROFILE [cm] -----
             CALL line_prof(grd(kmin:kmax), dgrd(kmin:kmax), wnum(l), gamma, sigma, iprof, prof=prof(kmin:kmax))
             ! -----------------------------
             ! ----- ABSORPTION COEFFICIENT [CM^-1] = number density [cm^-3] * cross section [cm^2] -----
             alp(kmin:kmax,j) = alp(kmin:kmax,j) + np(code0,j) * (strength * prof(kmin:kmax))
             ! ------------------------------------------------------------------------------------------
          END DO LINE
       END DO LAYER
!       if(myrk == 0) call progress(nc, nc_max, prog, dprog=1)
    END DO BLOCK

    DEALLOCATE(gfva, wnum, lowe, code, gm_0, gm_1, gm_2, q, prof)

    IF(kprc /= 1) THEN
       IF(myrk_line /= 0) THEN
          ALLOCATE(dummy(ke*je))
          CALL MPI_REDUCE(alp, dummy, ke*je, MPI_DOUBLE_PRECISION, MPI_SUM, 0, mpi_line_world, error)
          DEALLOCATE(dummy)
       ELSE
          CALL MPI_REDUCE(MPI_IN_PLACE, alp, ke*je, MPI_DOUBLE_PRECISION, MPI_SUM, 0, mpi_line_world, error)
       END IF
    END IF
   
    ! CLOSE location identifiers
    CALL h5Dclose_f(data_wnum, error)
    CALL h5Dclose_f(data_lowe, error)
    CALL h5Dclose_f(data_gfva, error)
    CALL h5Dclose_f(data_code, error)
    CALL h5Dclose_f(data_gm_0, error)
    CALL h5Dclose_f(data_gm_1, error)
    CALL h5Dclose_f(data_gm_2, error)
    CALL h5Gclose_f(grp_trans, error)
    CALL h5Fclose_f(file_spec, error)
    CALL h5close_f(error)

    DEALLOCATE(mass, wn_ion)

    RETURN
  END SUBROUTINE line_kurucz_p


  REAL(REAL64) FUNCTION g_nat(wnum)
    USE const_module, ONLY : pi, e2, m_ele, clight
    REAL(REAL64), INTENT(IN) :: wnum
    REAL(REAL64) :: omega
    omega = 2d0 * pi * clight * wnum
    g_nat = (2d0 / 3d0) * (e2 * omega**2) / (m_ele * clight**3)
    RETURN
  END FUNCTION g_nat

  
  REAL(REAL64) FUNCTION g_stark(r2u, r2l, temp, n_ele)
    USE const_module, ONLY : k_bol, hbar, pi, a_0, m_ele, ene_hyd
    REAL(REAL64), INTENT(IN) :: r2u, r2l
    REAL(REAL64), INTENT(IN) :: temp, n_ele
    REAL(REAL64), PARAMETER :: gaunt = 0.2d0
    
    g_stark = 2d0 * 8d0 * (pi / 3d0)**1.5d0 * (hbar * a_0 / m_ele) * &
         n_ele * (ene_hyd / (k_bol * temp))**0.5d0 * gaunt * &
         (r2u + r2l)
    
    RETURN
  END FUNCTION g_stark


  REAL(REAL64) FUNCTION g_vdw(r2u, r2l, temp, pol_pert, mass_pert, mass_rad)
    USE const_module, ONLY : k_bol, e2, a_0, hbar
    REAL(REAL64), INTENT(IN) :: r2u, r2l
    REAL(REAL64), INTENT(IN) :: temp, pol_pert, mass_pert, mass_rad
    REAL(REAL64) :: mass_pert2, v_pert, c6
    REAL(REAL64), PARAMETER :: g6_corr = 2.5d0 ! Valenti & Piskunov (1996)
    
    c6 = (e2 * a_0**2) / hbar * pol_pert * MAX((r2u - r2l), 0d0)
    mass_pert2 = mass_pert * mass_rad / (mass_pert + mass_rad)
    v_pert = SQRT(3d0 * k_bol * temp / mass_pert2)
    g_vdw = g6_corr * 8.08d0 * c6**0.4d0 * v_pert**0.6d0
    
    RETURN
  END FUNCTION g_vdw

  
  REAL(REAL64) FUNCTION g_vdw_phoenix(lowe, wnum, code, wn_ion, temp, mass_HI, mass)
    USE const_module, ONLY : ene_hyd, h, clight, amu, k_bol
    REAL(REAL64) :: lowe
    REAL(REAL64) :: wnum
    INTEGER(INT32) :: code, z
    REAL(REAL64) :: wn_ion
    REAL(REAL64) :: temp
    REAL(REAL64) :: mass_HI, mass
    REAL(REAL64), PARAMETER :: g6_corr = 2.5d0
    REAL(REAL64) :: c6cons, c6, rmass, vv1, pert6, gam6

    z = (code - (code / 100) * 100) + 1
    c6cons = g6_corr * 1.01d-32 * z**2 * (ene_hyd / (h * clight))**2
    c6 = c6cons * (1d0 / (wn_ion - (lowe + wnum))**2 - 1d0 / (wn_ion - lowe)**2) ! Unsoeld (1955)
    rmass = mass_HI * (1d0 / (1d0 + mass_HI / mass))
    vv1 = 0.958d0 * (2d0 / (rmass * amu))**0.30d0
    IF(c6 > 0d0) THEN
       gam6 = 8.08d0 * vv1 * c6**0.4d0
    ELSE
       gam6 = 0d0
    END IF
    pert6 = (k_bol * temp)**0.3d0
    g_vdw_phoenix = gam6 * pert6
    
    RETURN
  END FUNCTION g_vdw_phoenix
  
  
  REAL(REAL64) FUNCTION g_stark_phoenix(highe, code, wn_ion)
    USE const_module, ONLY : ene_hyd, h, clight
    REAL(REAL64) :: highe
    INTEGER(INT32) :: code, z
    REAL(REAL64) :: wn_ion
    REAL(REAL64) :: nsqd
    
    z = (code - (code / 100) * 100) + 1
    ! effective atomic radius squared for evaluating pressure broadening
    nsqd = (z**2) * (ene_hyd / (h * clight)) / (wn_ion - highe)
    IF(nsqd < 0d0) nsqd = 25d0
    g_stark_phoenix = 1d-8 * nsqd**1.5d0
    
    RETURN
  END FUNCTION g_stark_phoenix
                

  FUNCTION map(code)
    INTEGER(INT32) :: map
    INTEGER(INT16) :: code
    INTEGER(INT32), PARAMETER :: imax = 840
    INTEGER(INT32), PARAMETER, DIMENSION(imax) :: elem = [&
         100, 101, &
         200, 201, 202, &
         300, 301, 302, 303, &
         400, 401, 402, 403, 404, &
         500, 501, 502, 503, 504, 505, &
         600, 601, 602, 603, 604, 605, 606, &
         700, 701, 702, 703, 704, 705, 706, 707, &
         800, 801, 802, 803, 804, 805, 806, 807, 808, &
         900, 901, 902, 903, 904, 905, 906, 907, 908, 909,&
         1000,1001,1002,1003,1004,1005,1006,1007,1008,1009,1010, &
         1100,1101,1102,1103,1104,1105,1106,1107,1108,1109,1110,1111, &
         1200,1201,1202,1203,1204,1205,1206,1207,1208,1209,1210,1211,1212, &
         1300,1301,1302,1303,1304,1305,1306,1307,1308,1309,1310,1311,1312,1313, &
         1400,1401,1402,1403,1404,1405,1406,1407,1408,1409,1410,1411,1412,1413,1414, &
         1500,1501,1502,1503,1504,1505,1506,1507,1508,1509,1510,1511,1512,1513,1514,1515, &
         1600,1601,1602,1603,1604,1605,1606,1607,1608,1609,1610,1611,1612,1613,1614,1615,1616, &
         1700,1701,1702,1703,1704,1705,1706,1707,1708,1709,1710,1711,1712,1713,1714,1715,1716,1717, &
         1800,1801,1802,1803,1804,1805,1806,1807,1808,1809,1810,1811,1812,1813,1814,1815,1816,1817,1818,&
         1900,1901,1902,1903,1904,1905,1906,1907,1908,1909,1910,1911,1912,1913,1914,1915,1916,1917,1918,1919,&
         2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018,2019,2020,&
         2100,2101,2102,2103,2104,2105,2106,2107,2108,2109,2110,2111,2112,2113,2114,2115,2116,2117,2118,2119,2120,&
         2121, &
         2200,2201,2202,2203,2204,2205,2206,2207,2208,2209,2210,2211,2212,2213,2214,2215,2216,2217,2218,2219,2220,&
         2221,2222, &
         2300,2301,2302,2303,2304,2305,2306,2307,2308,2309,2310,2311,2312,2313,2314,2315,2316,2317,2318,2319,2320,&
         2321,2322,2323,&
         2400,2401,2402,2403,2404,2405,2406,2407,2408,2409,2410,2411,2412,2413,2414,2415,2416,2417,2418,2419,2420,&
         2421,2422,2423,2424,&
         2500,2501,2502,2503,2504,2505,2506,2507,2508,2509,2510,2511,2512,2513,2514,2515,2516,2517,2518,2519,2520,&
         2521,2522,2523,2524,2525,&
         2600,2601,2602,2603,2604,2605,2606,2607,2608,2609,2610,2611,2612,2613,2614,2615,2616,2617,2618,2619,2620,&
         2621,2622,2623,2624,2625,2626,&
         2700,2701,2702,2703,2704,2705,2706,2707,2708,2709,2710,2711,2712,2713,2714,2715,2716,2717,2718,2719,2720,&
         2721,2722,2723,2724,2725,2726,2727,&
         2800,2801,2802,2803,2804,2805,2806,2807,2808,2809,2810,2811,2812,2813,2814,2815,2816,2817,2818,2819,2820,&
         2821,2822,2823,2824,2825,2826,2827,2828,&
         2900,2901,2902,2903,2904,2905,2906,2907,2908,2909,2910,2911,2912,2913,2914,2915,2916,2917,2918,2919,2920,&
         2921,2922,2923,2924,2925,2926,2927,2928,2929,&
         3000,3001,3002,3003,3004,3005,3006,3007,3008,3009,3010,3011,3012,3013,3014,3015,3016,3017,3018,3019,3020,&
         3021,3022,3023,3024,3025,3026,3027,3028,3029,3030,&
         3100,3101,3102,3103,3104,&
         3200,3201,3202,3203,3204,&
         3300,3301,3302,3303,3304,&
         3400,3401,3402,3403,3404,&
         3500,3501,3502,3503,3504,&
         3600,3601,3602,3603,3604,&
         3700,3701,3702,3703,3704,&
         3800,3801,3802,3803,3804,&
         3900,3901,3902,3903,3904,&
         4000,4001,4002,4003,4004,&
         4100,4101,4102,4103,4104,&
         4200,4201,4202,4203,4204,&
         4300,4301,4302,4303,4304,&
         4400,4401,4402,4403,4404,&
         4500,4501,4502,4503,4504,&
         4600,4601,4602,4603,4604,&
         4700,4701,4702,4703,4704,&
         4800,4801,4802,4803,4804,&
         4900,4901,4902,4903,4904,&
         5000,5001,5002,5003,5004,&
         5100,5101,5102,5103,5104,&
         5200,5201,5202,5203,5204,&
         5300,5301,5302,5303,5304,&
         5400,5401,5402,5403,5404,&
         5500,5501,5502,5503,5504,&
         5600,5601,5602,5603,5604,&
         5700,5701,5702,5703,5704,&
         5800,5801,5802,5803,5804,&
         5900,5901,5902,5903,5904,&
         6000,6001,6002,6003,6004,&
         6100,6101,6102,6103,6104,&
         6200,6201,6202,6203,6204,&
         6300,6301,6302,6303,6304,&
         6400,6401,6402,6403,6404,&
         6500,6501,6502,6503,6504,&
         6600,6601,6602,6603,6604,&
         6700,6701,6702,6703,6704,&
         6800,6801,6802,6803,6804,&
         6900,6901,6902,6903,6904,&
         7000,7001,7002,7003,7004,&
         7100,7101,7102,7103,7104,&
         7200,7201,7202,7203,7204,&
         7300,7301,7302,7303,7304,&
         7400,7401,7402,7403,7404,&
         7500,7501,7502,7503,7504,&
         7600,7601,7602,7603,7604,&
         7700,7701,7702,7703,7704,&
         7800,7801,7802,7803,7804,&
         7900,7901,7902,7903,7904,&
         8000,8001,8002,8003,8004,&
         8100,8101,8102,8103,8104,&
         8200,8201,8202,8203,8204,&
         8300,8301,8302,8303,8304,&
         8400,8401,8402,8403,8404,&
         8500,8501,8502,8503,8504,&
         8600,8601,8602,8603,8604,&
         8700,8701,8702,8703,8704,&
         8800,8801,8802,8803,8804,&
         8900,8901,8902,8903,8904,&
         9000,9001,9002,9003,9004,&
         9100,9101,9102,9103,9104,&
         9200,9201,9202,9203,9204,&
         9300,9301,9302,9303,9304,&
         9400,9401,9402,9403,9404,&
         9500,9501,9502,9503,9504,&
         9600,9601,9602,9603,9604,&
         9700,9701,9702,9703,9704,&
         9800,9801,9802,9803,9804,&
         9900,9901,9902,9903,9904]

    map = elem(code - 84)
    RETURN
  END FUNCTION map



  FUNCTION find_grid(grd, wnum)

    IMPLICIT NONE
    INTEGER(INT32) :: find_grid
    REAL(REAL64), INTENT(IN) :: grd(:)
    REAL(REAL64), INTENT(IN) :: wnum
    INTEGER(INT32) :: m, m1, m2

    m1 = 1
    m2 = UBOUND(grd,1)

    IF(wnum > grd(m2)) THEN
       find_grid = m2
       RETURN
    END IF
    IF(wnum < grd(m1)) THEN
       find_grid = m1 - 1
       RETURN
    END IF
    
    DO WHILE((m2 - m1) > 1)
       m = (m1 + m2) / 2
       IF(wnum > grd(m)) THEN
          m1 = m
       ELSE
          m2 = m
       END IF
    END DO
    
    find_grid = m1
    
    RETURN
  END FUNCTION find_grid
  
  
  SUBROUTINE compute_q(temp, ene, gtot, q)

    USE const_module, ONLY : c2
    IMPLICIT NONE
    REAL(REAL64), INTENT(IN) :: temp
    REAL(REAL64), INTENT(IN) :: ene(:)
    INTEGER(INT32), INTENT(IN) :: gtot(:)
    REAL(REAL64), INTENT(OUT) :: q
    
    INTEGER(INT64) :: l

    q = 0d0
    DO l = 1, UBOUND(ene,1)
       q = q + gtot(l) * EXP(-c2 * ene(l) / temp)
    END DO

    RETURN
  END SUBROUTINE compute_q
  
         
  SUBROUTINE interpolate_q(temp0, temp, pf, q)

    IMPLICIT NONE
    REAL(REAL64), INTENT(IN) :: temp0(:)
    REAL(REAL64), INTENT(IN) :: temp(:)
    REAL(REAL64), INTENT(IN) :: pf(:)
    REAL(REAL64), INTENT(OUT) :: q(:)
    INTEGER(INT32) :: ntemp, n, j
    
    ntemp = UBOUND(temp,1)

    DO j = 1, UBOUND(temp0,1)
       IF(temp0(j) < temp(1)) THEN
          q(j) = 0d0
       ELSE IF(temp0(j) > temp(ntemp)) THEN
          q(j) = 0d0
       ELSE
          DO n = 1, ntemp - 1
             IF(temp0(j) <= temp(n+1)) EXIT
          END DO
          q(j) = ((temp0(j) - temp(n)) * pf(n+1) + (temp(n+1) - temp0(j)) * pf(n)) / (temp(n+1) - temp(n))
       END IF
    END DO
    
    RETURN
  END SUBROUTINE interpolate_q

  
  SUBROUTINE line_prof(grd, dgrd, wnum, gamma, sigma, iprof, prof)

    USE const_module, ONLY : pi
    USE voigt_module, ONLY : hum1zpf16, hum1wei24, hum1wei32, voigt_bs
    IMPLICIT NONE
    REAL(REAL64), INTENT(IN) :: grd(:)
    REAL(REAL64), INTENT(IN) :: dgrd(:)
    REAL(REAL64), INTENT(IN) :: wnum
    REAL(REAL64), INTENT(IN) :: gamma
    REAL(REAL64), INTENT(IN) :: sigma
    CHARACTER(LEN=2), INTENT(IN) :: iprof
    REAL(REAL64), INTENT(OUT) :: prof(:)
    COMPLEX(REAL64), ALLOCATABLE :: vgtKL(:)

    SELECT CASE (iprof)
       CASE('Ga') ! Gaussian
          prof = EXP(-((grd - wnum) / sigma)**2) / (SQRT(pi) * sigma)
       CASE('Lo') ! Lorentizan
          prof = gamma / ((grd - wnum)**2 + gamma**2) / pi
       CASE('bG') ! binned Gaussian
          ! Yurchenko et al. 2018, A&A 614, A131
          prof = (ERF(((grd + dgrd / 2d0) - wnum) / sigma) - ERF(((grd - dgrd / 2d0) - wnum) / sigma)) / (2d0 * dgrd)
       CASE('Vo') ! Voigt
#ifndef WAVELENGTH_BASE
          CALL voigt_bs(ABS(grd - wnum)/sigma, gamma/sigma, prof)
#else          
          CALL voigt_bs(ABS(1d0/grd - 1d0/wnum)*wnum**2/sigma, gamma/sigma, prof)
#endif          
          prof = prof / (SQRT(pi) * sigma)
       CASE('Vs')
          IF(gamma / sigma > 1d-8) THEN !Voigt
          ! Schreier 2018, MNRAS, 479, 3068
             ALLOCATE(vgtKL(UBOUND(grd,1)))
             CALL hum1wei24(UBOUND(grd,1) - LBOUND(grd,1), (grd - wnum) / sigma, gamma / sigma, vgtKL)
             !    CALL hum1wei32(UBOUND(grd,1) - LBOUND(grd,1), (grd - wnum) / sigma, gamma / sigma, vgtKL)
             !    CALL hum1zpf16(UBOUND(grd,1) - LBOUND(grd,1), (grd - wnum) / sigma, gamma / sigma, vgtKL)
             prof = DBLE(vgtKL) / (sigma * SQRT(pi))
             DEALLOCATE(vgtKL)
          ELSE ! Gaussian
             prof = EXP(-((grd - wnum) / sigma)**2) / (SQRT(pi) * sigma)
          ENDIF
       CASE DEFAULT
          PRINT *, '*** ERROR in line prof: ', iprof
          STOP
    END SELECT
    
    RETURN
  END SUBROUTINE line_prof
  

  REAL(REAL64) FUNCTION strength_cutoff(wnum, temp, wnum_crit, strength_cutoff0)

    USE const_module, ONLY : c2
    IMPLICIT NONE
    REAL(REAL64) :: wnum
    REAL(REAL64) :: temp
    REAL(REAL64) :: wnum_crit
    REAL(REAL64) :: strength_cutoff0
    IF(wnum >  wnum_crit) THEN
       strength_cutoff = strength_cutoff0
    ELSE
       strength_cutoff = strength_cutoff0 * (wnum / wnum_crit) * TANH(c2 * wnum / (2d0 * temp))
    END IF
    RETURN
  END FUNCTION strength_cutoff
  
#ifdef PRINT_PROGRESS  
  SUBROUTINE progress(nc, nc_max, prog, dprog)

    IMPLICIT NONE
    INTEGER(INT32), INTENT(IN) :: nc
    INTEGER(INT32), INTENT(IN) :: nc_max
    INTEGER(INT32), INTENT(INOUT) :: prog
    INTEGER(INT32), INTENT(IN) :: dprog
    
    if(DBLE(nc * 100) / DBLE(nc_max + 1) >= DBLE(prog + dprog)) then
       print '(I3,A)', INT(DBLE(nc * 100) / DBLE(nc_max + 1)), ' % DONE'
       prog = prog + dprog
    end if
    
    RETURN
  END SUBROUTINE progress
#endif  
END MODULE line_module
