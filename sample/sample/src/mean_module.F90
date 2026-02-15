MODULE mean_module
  
  USE ISO_FORTRAN_ENV
  IMPLICIT NONE

CONTAINS  

  SUBROUTINE mean_opac(pla, ros, alp, sca, temp, grd, temp2, pla2)
    USE MPI
    USE const_module, ONLY : c2, h, clight
    USE mpi_module, ONLY : mpi_grid_world, kprc
    
    REAL(REAL64), INTENT(OUT) :: pla   ! Planck-mean opacity
    REAL(REAL64), INTENT(OUT) :: ros   ! Rosseland-mean opacity
    REAL(REAL64), INTENT(IN) :: alp(:) ! monochromatic opacity
    REAL(REAL64), INTENT(IN) :: sca(:) ! monochromatic opacity
    REAL(REAL64), INTENT(IN) :: temp   ! temperature
    REAL(REAL64), INTENT(IN) :: grd(:) !wavenumber grid
    REAL(REAL64), INTENT(IN) :: temp2  ! radiation temperature
    REAL(REAL64), INTENT(OUT) :: pla2  ! Planck-mean opacity with temp2

    REAL(REAL64), ALLOCATABLE :: b(:), dbdt(:) ! Planck function
    REAL(REAL64) :: nume, deno, weight
    INTEGER(INT64) :: k, ks, ke
    INTEGER(INT32) :: error

    ks = LBOUND(alp,1)
    ke = UBOUND(alp,1)

    ALLOCATE(b(ks:ke), dbdt(ks:ke))
    
    ! PLANCK FUNCTION
    DO k = ks, ke
       b(k) = grd(k)**3 / (EXP(c2 * grd(k) / temp) - 1d0)
       dbdt(k) = grd(k)**4 / temp**2 * EXP(-c2 * grd(k) / temp) / (EXP(-c2 * grd(k) / temp) - 1d0)**2
    END DO
    b = b * (2d0 * h * clight**2)
    dbdt = dbdt * (2d0 * h * clight**2 * c2)

    ! PLANCK
    nume = 0d0
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
       deno = deno + b(k)                * weight
    END DO
    IF(kprc /= 1) THEN
       CALL MPI_ALLREDUCE(MPI_IN_PLACE, nume, 1, MPI_DOUBLE_PRECISION, MPI_SUM, mpi_grid_world, error)
       CALL MPI_ALLREDUCE(MPI_IN_PLACE, deno, 1, MPI_DOUBLE_PRECISION, MPI_SUM, mpi_grid_world, error)
    END IF
    pla = REAL(nume / deno)

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
    IF(kprc /= 1) THEN
       CALL MPI_ALLREDUCE(MPI_IN_PLACE, nume, 1, MPI_DOUBLE_PRECISION, MPI_SUM, mpi_grid_world, error)
       CALL MPI_ALLREDUCE(MPI_IN_PLACE, deno, 1, MPI_DOUBLE_PRECISION, MPI_SUM, mpi_grid_world, error)
    END IF
    ros = REAL(nume / deno)

    ! PLANCK FUNCTION with temp2
    DO k = ks, ke
       b(k) = grd(k)**3 / (EXP(c2 * grd(k) / temp2) - 1d0)
    END DO
    b = b * (2d0 * h * clight**2)

    ! PLANCK
    nume = 0d0
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
       deno = deno + b(k)                * weight
    END DO
    IF(kprc /= 1) THEN
       CALL MPI_ALLREDUCE(MPI_IN_PLACE, nume, 1, MPI_DOUBLE_PRECISION, MPI_SUM, mpi_grid_world, error)
       CALL MPI_ALLREDUCE(MPI_IN_PLACE, deno, 1, MPI_DOUBLE_PRECISION, MPI_SUM, mpi_grid_world, error)
    END IF
    pla2 = REAL(nume / deno)

    DEALLOCATE(b, dbdt)
    
    RETURN
  END SUBROUTINE mean_opac


  SUBROUTINE mean_opac_lambda(pla, ros, alp, sca, temp, grd, temp2, pla2)
    USE ISO_FORTRAN_ENV
    USE MPI
    USE const_module, ONLY : c2, h, clight, k_bol
    USE mpi_module, ONLY : mpi_grid_world, kprc
    IMPLICIT NONE
    
    REAL(REAL64), INTENT(OUT) :: pla ! Planck-mean opacity
    REAL(REAL64), INTENT(OUT) :: ros ! Rosseland-mean opacity
    REAL(REAL64), INTENT(IN) :: alp(:) ! monochromatic opacity
    REAL(REAL64), INTENT(IN) :: sca(:) ! monochromatic opacity
    REAL(REAL64), INTENT(IN) :: temp ! temperature
    REAL(REAL64), INTENT(IN) :: grd(:) !wavenumber grid
    REAL(REAL64), INTENT(IN) :: temp2 ! radiation temperature
    REAL(REAL64), INTENT(OUT) :: pla2 ! Planck-mean opacity with temp2

    REAL(REAL64), ALLOCATABLE :: b(:), dbdt(:) ! Planck function
    REAL(REAL64) :: nume, deno, lambda, weight
    INTEGER(INT64) :: k, ks, ke
    INTEGER(INT32) :: error
    REAL(REAL64) :: cons, a, c1_locln, dc1_locln, am1i
    REAL(REAL64), parameter :: c1_loc = 2d0 * h * clight**2, &
         c2_loc = h * clight / k_bol

    ks = LBOUND(alp,1)
    ke = UBOUND(alp,1)

    ALLOCATE(b(ks:ke), dbdt(ks:ke))

    ! PLANCK FUNCTION
    DO k = ks, ke
       lambda = 1d0 / grd(k)
       cons = c2_loc / lambda
       a = cons / temp
       c1_locln = log(c1_loc) - 5d0 * log(lambda)
       dc1_locln = log(c1_loc * c2_loc) - 6d0 * log(lambda)
       IF(a > 40d0) THEN
          b(k) = exp(c1_locln - a)
          dbdt(k) = exp(dc1_locln - (a + 2d0 * log(temp)))
       ELSE
          a = exp(a)
          am1i = 1d0 / (a - 1d0)
          b(k) = c1_loc / lambda**5 * am1i
          dbdt(k) = (c1_loc * c2_loc) / lambda**6 * a * (am1i / temp)**2
       END IF
    END DO

    ! PLANCK
    nume = 0d0
    deno = 0d0
    DO k = ks, ke
       IF(k == ks) THEN
          weight = 0.5d0 * (1d0 / grd(k  ) - 1d0 / grd(k+1))
       ELSE IF(k == ke) THEN
          weight = 0.5d0 * (1d0 / grd(k-1) - 1d0 / grd(k  ))
       ELSE
          weight = 0.5d0 * (1d0 / grd(k-1) - 1d0 / grd(k+1))
       END IF
       nume = nume + b(k) * DBLE(alp(k)) * weight
       deno = deno + b(k)                * weight
    END DO
    IF(kprc /= 1) THEN
       CALL MPI_ALLREDUCE(MPI_IN_PLACE, nume, 1, MPI_DOUBLE_PRECISION, MPI_SUM, mpi_grid_world, error)
       CALL MPI_ALLREDUCE(MPI_IN_PLACE, deno, 1, MPI_DOUBLE_PRECISION, MPI_SUM, mpi_grid_world, error)
    END IF
    pla = REAL(nume / deno)

    ! ROSSELAND
    nume = 0d0
    deno = 0d0
    DO k = ks, ke
       IF(alp(k) + sca(k) == 0d0) CYCLE
       IF(k == ks) THEN
          weight = 0.5d0 * (1d0 / grd(k  ) - 1d0 / grd(k+1))
       ELSE IF(k == ke) THEN
          weight = 0.5d0 * (1d0 / grd(k-1) - 1d0 / grd(k  ))
       ELSE
          weight = 0.5d0 * (1d0 / grd(k-1) - 1d0 / grd(k+1))
       END IF
       nume = nume + dbdt(k)                           * weight
       deno = deno + dbdt(k) / DBLE((alp(k) + sca(k))) * weight
    END DO

    IF(kprc /= 1) THEN
       CALL MPI_ALLREDUCE(MPI_IN_PLACE, nume, 1, MPI_DOUBLE_PRECISION, MPI_SUM, mpi_grid_world, error)
       CALL MPI_ALLREDUCE(MPI_IN_PLACE, deno, 1, MPI_DOUBLE_PRECISION, MPI_SUM, mpi_grid_world, error)
    END IF
    ros = REAL(nume / deno)

    ! PLANCK FUNCTION with temp2
    DO k = ks, ke
       lambda = 1d0 / grd(k)
       cons = c2_loc / lambda
       c1_locln = log(c1_loc) - 5d0 * log(lambda)
       a = cons / temp2
       if(a > 40d0) then
          b(k) = exp(c1_locln - a)
       else
          a = exp(a)
          am1i = 1d0 / (a - 1d0)
          b(k) = c1_loc / lambda**5 * am1i
       endif
    END DO

    ! PLANCK
    nume = 0d0
    deno = 0d0
    DO k = ks, ke
       IF(k == ks) THEN
          weight = 0.5d0 * (1d0 / grd(k  ) - 1d0 / grd(k+1))
       ELSE IF(k == ke) THEN
          weight = 0.5d0 * (1d0 / grd(k-1) - 1d0 / grd(k  ))
       ELSE
          weight = 0.5d0 * (1d0 / grd(k-1) - 1d0 / grd(k+1))
       END IF
       nume = nume + DBLE(alp(k)) * b(k) * weight
       deno = deno + b(k)                * weight
    END DO
    IF(kprc /= 1) THEN
       CALL MPI_ALLREDUCE(MPI_IN_PLACE, nume, 1, MPI_DOUBLE_PRECISION, MPI_SUM, mpi_grid_world, error)
       CALL MPI_ALLREDUCE(MPI_IN_PLACE, deno, 1, MPI_DOUBLE_PRECISION, MPI_SUM, mpi_grid_world, error)
    END IF
    pla2 = REAL(nume / deno)

    DEALLOCATE(b, dbdt)
    
    RETURN
  END SUBROUTINE mean_opac_lambda

END MODULE mean_module
