
MODULE output_module

  USE ISO_FORTRAN_ENV
  IMPLICIT NONE

CONTAINS  

  SUBROUTINE output_table(np, temp, rho, mean_p, mean_r, temp2, mean_p2)
    USE HDF5
    USE H5LT
    USE MPI
    USE h5LX_module, ONLY : h5LXset_dataset, h5LXwrite_hyperslab_dataset_F64
    USE mpi_module, ONLY : myrk
    REAL(REAL64), INTENT(IN) :: np(:,:)
    REAL(REAL64), INTENT(IN) :: temp(:)
    REAL(REAL64), INTENT(IN) :: rho(:)
    REAL(REAL64), INTENT(IN) :: mean_p(:)
    REAL(REAL64), INTENT(IN) :: mean_r(:)
    REAL(REAL64), INTENT(IN) :: temp2
    REAL(REAL64), INTENT(IN) :: mean_p2(:)
    
    INTEGER(HID_T) :: file!, data_temp, data_rho, data_mean_p, data_mean_r, data_mean_p2, data_ndens
    INTEGER(HSIZE_T) :: dims(1), dims2(2)!, dims_global(1), dims2_global(2), start(1), count(1), start2(2), count2(2)
    INTEGER(INT32) :: error, nmax, jmax

    jmax = UBOUND(temp,1)
    
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, mean_p, jmax, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, error)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, mean_r, jmax, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, error)
    CALL MPI_ALLREDUCE(MPI_IN_PLACE, mean_p2, jmax, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, error)
    
    IF(myrk == 0) THEN
       CALL h5open_f(error)
       CALL h5Fcreate_f('output/opac_table.h5', H5F_ACC_TRUNC_F, file, error)
       dims = [1]
       CALL h5LTmake_dataset_f(file, 'jmax', 1, dims, H5T_STD_I32LE, jmax, error)
       CALL h5LTmake_dataset_f(file, 'temp2', 1, dims, H5T_IEEE_F64LE, temp2, error)
       
       dims = [jmax]
       CALL h5LTmake_dataset_f(file, 'temp', 1, dims, H5T_IEEE_F64LE, temp, error)
       CALL h5LTmake_dataset_f(file, 'rho', 1, dims, H5T_IEEE_F64LE, rho, error)
       CALL h5LTmake_dataset_f(file, 'mean_p', 1, dims, H5T_IEEE_F64LE, mean_p, error)
       CALL h5LTmake_dataset_f(file, 'mean_r', 1, dims, H5T_IEEE_F64LE, mean_r, error)
       CALL h5LTmake_dataset_f(file, 'mean_p2', 1, dims, H5T_IEEE_F64LE, mean_p2, error)
       nmax = UBOUND(np,1)
       dims2 = [nmax, jmax]
       CALL h5LTmake_dataset_f(file, 'ndens', 2, dims2, H5T_IEEE_F64LE, np(1:nmax,1:jmax), error)
       
       CALL h5Fclose_f(file, error)
       CALL h5close_f(error)
    END IF
    
    RETURN
  END SUBROUTINE output_table

  
  SUBROUTINE output_mono(j, temp, nden, grd, alp, sca, cnt, temp2, rho, pmean, pmean2)
    USE HDF5
    USE H5LT
    USE const_module, ONLY : c2, h, c, pi

    INTEGER(INT32), INTENT(IN) :: j
    REAL(REAL64), INTENT(IN) :: temp
    REAL(REAL64), INTENT(IN) :: nden(:)
    REAL(REAL64), INTENT(IN) :: grd(:)
    REAL(REAL64), INTENT(IN) :: alp(:)
    REAL(REAL64), INTENT(IN) :: sca(:)
    REAL(REAL64), INTENT(IN) :: cnt(:)
    REAL(REAL64), INTENT(IN) :: temp2
    REAL(REAL64), INTENT(IN) :: rho
    REAL(REAL64), INTENT(IN) :: pmean
    REAL(REAL64), INTENT(IN) :: pmean2
    INTEGER(INT32) :: ns
    INTEGER :: error
    INTEGER(HID_T) :: file_id!, plist_id, data_grid, data_abs, data_sca
    INTEGER(HSIZE_T), DIMENSION(1) :: dims1
    REAL(REAL64), ALLOCATABLE :: b(:), dbdt(:) ! Planck function
    REAL(REAL64) :: nume, deno, weight, pla, pla2, ros, plac, plal, numec, plac2, plal2
    INTEGER(INT64) :: k, ks, ke, k_total

    ks = LBOUND(alp,1)
    ke = UBOUND(alp,1)

    ! ------------------------------------------
    ! MEAN OPACITIES

    ALLOCATE(b(ks:ke), dbdt(ks:ke))
    
    ! PLANCK FUNCTION
    DO k = ks, ke
       b(k) = grd(k)**3 / (EXP(c2 * grd(k) / temp) - 1d0)
       dbdt(k) = grd(k)**4 / temp**2 * EXP(-c2 * grd(k) / temp) / (EXP(-c2 * grd(k) / temp) - 1d0)**2
    END DO
    b = b * (2d0 * h * c**2)
    dbdt = dbdt * (2d0 * h * c**2 * c2)

    ! PLANCK
    nume = 0d0
    numec= 0d0
    deno = 0d0
    DO k = ks, ke
       IF(k == ks) THEN
          weight = (grd(k+1) - grd(k  )) * 0.5d0
       ELSE IF(k == ke) THEN
          weight = (grd(k  ) - grd(k-1)) * 0.5d0
       ELSE
          weight = (grd(k+1) - grd(k-1)) * 0.5d0
       END IF
       nume = nume + b(k) * DBLE(alp(k)) * weight
       numec= numec+ b(k) * DBLE(cnt(k) - sca(k)) * weight
       deno = deno + b(k)                * weight
    END DO
    pla = REAL(nume / deno)
    plac= REAL(numec/ deno)
    plal= REAL(pmean / ((pi**4 / 15d0) * (temp / c2)**4))

    ! ROSSELAND
    nume = 0d0
    deno = 0d0
    DO k = ks, ke
       IF(alp(k) + sca(k) == 0d0) CYCLE
       IF(k == ks) THEN
          weight = (grd(k+1) - grd(k  )) * 0.5d0
       ELSE IF(k == ke) THEN
          weight = (grd(k  ) - grd(k-1)) * 0.5d0
       ELSE
          weight = (grd(k+1) - grd(k-1)) * 0.5d0
       END IF
       nume = nume + dbdt(k)                           * weight
       deno = deno + dbdt(k) / DBLE((alp(k) + sca(k))) * weight
    END DO
    ros = REAL(nume / deno)

    ! PLANCK FUNCTION with temp2
    DO k = ks, ke
       b(k) = grd(k)**3 / (EXP(c2 * grd(k) / temp2) - 1d0)
    END DO
    b = b * (2d0 * h * c**2)

    ! PLANCK
    nume = 0d0
    numec = 0d0
    deno = 0d0
    DO k = ks, ke
       IF(k == ks) THEN
          weight = (grd(k+1) - grd(k  )) * 0.5d0
       ELSE IF(k == ke) THEN
          weight = (grd(k  ) - grd(k-1)) * 0.5d0
       ELSE
          weight = (grd(k+1) - grd(k-1)) * 0.5d0
       END IF
       nume = nume + DBLE(alp(k)) * b(k) * weight
       numec= numec+ b(k) * DBLE(cnt(k) - sca(k)) * weight
       deno = deno + b(k)                * weight
    END DO
    pla2 = REAL(nume / deno)
    plac2= REAL(numec/ deno)
    plal2= REAL(pmean2 / ((pi**4 / 15d0) * (temp2 / c2)**4))

    DEALLOCATE(b, dbdt)
    ! ------------------------------------------

    k_total = ke - ks + 1
    ns = UBOUND(nden,1) - LBOUND(nden,1) + 1

    CALL h5open_f(error)
    CALL h5Fcreate_f('output/mono_'//i2c(j)//'.h5', H5F_ACC_TRUNC_F, file_id, error)
    
    dims1 = [1]
    CALL h5LTmake_dataset_f(file_id, 'temp', 1, dims1, H5T_IEEE_F64LE, temp, error)
    CALL h5LTmake_dataset_f(file_id, 'temp2', 1, dims1, H5T_IEEE_F64LE, temp2, error)
    CALL h5LTmake_dataset_f(file_id, 'rho', 1, dims1, H5T_IEEE_F64LE, rho, error)
    ! mean opacities
    ! grid-based Planck means
    CALL h5LTmake_dataset_f(file_id, 'pla', 1, dims1, H5T_IEEE_F64LE, pla, error)
    CALL h5LTmake_dataset_f(file_id, 'pla2', 1, dims1, H5T_IEEE_F64LE, pla2, error)
    ! line-based Planck means and continuum Planck means
    CALL h5LTmake_dataset_f(file_id, 'plac', 1, dims1, H5T_IEEE_F64LE, plac, error)
    CALL h5LTmake_dataset_f(file_id, 'plal', 1, dims1, H5T_IEEE_F64LE, plal, error)
    CALL h5LTmake_dataset_f(file_id, 'plac2', 1, dims1, H5T_IEEE_F64LE, plac2, error)
    CALL h5LTmake_dataset_f(file_id, 'plal2', 1, dims1, H5T_IEEE_F64LE, plal2, error)
    ! Rosseland means
    CALL h5LTmake_dataset_f(file_id, 'ros', 1, dims1, H5T_IEEE_F64LE, ros, error)
    dims1 = [ns]
    CALL h5LTmake_dataset_f(file_id, 'nden', 1, dims1, H5T_IEEE_F64LE, nden, error)
    ! monochromatic opacities
    dims1 = [k_total]
    CALL h5LTmake_dataset_f(file_id, 'grd', 1, dims1, H5T_IEEE_F64LE, grd, error)
    CALL h5LTmake_dataset_f(file_id, 'abs', 1, dims1, H5T_IEEE_F64LE, alp, error)
    CALL h5LTmake_dataset_f(file_id, 'sca', 1, dims1, H5T_IEEE_F64LE, sca, error)
    CALL h5LTmake_dataset_f(file_id, 'cnt', 1, dims1, H5T_IEEE_F64LE, cnt, error)

    CALL h5Fclose_f(file_id, error)
    CALL h5close_f(error)

    RETURN
  END SUBROUTINE output_mono


  FUNCTION i2c(j)
    IMPLICIT NONE
    CHARACTER*5 :: i2c
    INTEGER(INT32) :: j
    WRITE(i2c,FMT='(I5.5)') j
    RETURN
  END FUNCTION i2c

END MODULE output_module
