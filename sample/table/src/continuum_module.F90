
! Use Karzas & Latter's table in place of Mathisen's formula for H I and He II
#define KARZAS_IN_PLACE_OF_MATHISEN

MODULE continuum_module

  USE ISO_FORTRAN_ENV
  USE HDF5
  USE H5LT
  USE MPI
  USE mpi_module, ONLY : myrk

CONTAINS

  SUBROUTINE brems_atom_vanhoof_2014(temp, grd, gaunt)
    ! Gaunt factors for bremsstrahlung of positive atomic ions
    ! van Hoof et al. 2014, MNRAS, 444, 420
    USE const_module, ONLY : alpha, clight, k_bol, m_ele, h, rydberg
    USE code_module, ONLY : ZnI
    IMPLICIT NONE

    REAL(REAL64), INTENT(IN) :: temp(:)
    REAL(REAL64), INTENT(IN) :: grd(:)
    REAL(REAL64), INTENT(INOUT) :: gaunt(:,:,:)

    INTEGER(INT32) :: n, i, j, imax, jmax, ls, le, l, z
    REAL(REAL64), ALLOCATABLE :: g(:,:), gam2(:), u(:), u0(:)
    REAL(REAL64) :: step, gam2_min, u_min, gam20

    gaunt = 0d0

    OPEN(1, FILE='input/1016620_Supplementary_Data/gauntff.dat', STATUS='OLD')
    DO n = 1, 29
       READ(1,*)
    END DO
    READ(1,*) imax, jmax
    ALLOCATE(g(imax,jmax), gam2(imax), u(jmax))
    ALLOCATE(u0(UBOUND(grd,1)))

    READ(1,*) gam2_min
    READ(1,*) u_min
    READ(1,*) step
    DO n = 1, 9
       READ(1,*)
    END DO
    READ(1,*) ((g(i,j), i = 1, imax), j = 1, jmax)
    CLOSE(1)
    
    DO i = 1, imax
       gam2(i) = gam2_min + (i - 1) * step
    END DO
    DO j = 1, jmax
       u(j) = u_min + (j - 1) * step
    END DO

    ls = LBOUND(temp,1)
    le = UBOUND(temp,1)

    DO l = ls, le
       u0(:) = LOG10((h * clight * grd(:)) / (k_bol * temp(l)))
       DO z = 1, (ZnI/100)
          gam20 = LOG10((z**2 * rydberg) / (k_bol * temp(l)))
          CALL interp2d_yarr(gam20, u0(:), gam2(:), u(:), g(:,:), gaunt(:,l,z))
       END DO
    END DO

    DEALLOCATE(g, gam2, u)
    DEALLOCATE(u0)
    RETURN
    
  END SUBROUTINE brems_atom_vanhoof_2014


  SUBROUTINE brems_atom_ferland(temp, grd, gaunt)
    ! Gaunt factors for bremsstrahlung of positive atomic ions
    USE const_module, ONLY : alpha, clight, k_bol, m_ele, h, rydberg
    USE code_module, ONLY : ZnI
    USE gffree_module, ONLY : gffree
    IMPLICIT NONE

    REAL(REAL64), INTENT(IN) :: temp(:)
    REAL(REAL64), INTENT(IN) :: grd(:)
    REAL(REAL64), INTENT(INOUT) :: gaunt(:,:,:)

    INTEGER(INT32) :: ls, le, l, z
    REAL(REAL64), ALLOCATABLE :: u0(:)
    REAL(REAL64) :: gam20
    INTEGER(INT32) :: k, ks, ke

    gaunt = 0d0

    ks = LBOUND(grd,1)
    ke = UBOUND(grd,1)

    ALLOCATE(u0(UBOUND(grd,1)))

    ls = LBOUND(temp,1)
    le = UBOUND(temp,1)

    DO l = ls, le
       u0(ks:ke) = 1.43883d8 * (grd(ks:ke)/1d8) / temp(l)
       DO z = 1, (ZnI/100)
          gam20 = z**2 * 1.5780787d5 / temp(l)
          DO k = ks, ke
             gaunt(k,l,z) = gffree(gam20, u0(k))
          END DO
       END DO
    END DO

    DEALLOCATE(u0)
    RETURN
    
  END SUBROUTINE brems_atom_ferland


  SUBROUTINE brems_hm_john_1988(temp, grd, coeff)
    ! coefficients of bremsstrahlung of H-
    ! John 1988, A&A, 193, 189
    USE const_module, ONLY : alpha, hbar, pi, clight, k_bol, m_ele, ev_to_erg, h_electron_aff_ev, h, c2, megabarn_to_cm2
    IMPLICIT NONE
    REAL(REAL64), INTENT(IN) :: temp(:)
    REAL(REAL64), INTENT(IN) :: grd(:)
    REAL(REAL64), INTENT(OUT) :: coeff(:,:)

    INTEGER(INT32) :: ks, ke, k, n, js, je, j
    REAL(REAL64), PARAMETER, DIMENSION(6,6) :: &
         a1 = RESHAPE([&
         0d0, 0d0, 0d0, 0d0, 0d0, 0d0,&
         2483.3460d0, 285.8270d0, -2054.2910d0, 2827.7760d0, -1341.5370d0, 208.9520d0,&
         -3449.8890d0, -1158.3820d0, 8746.5230d0, -11485.6320d0, 5303.6090d0, -812.9390d0, &
         2200.0400d0, 2427.7190d0, -13651.1050d0, 16755.5240d0, -7510.4940d0, 1132.7380d0, &
         -696.2710d0, -1841.4000d0, 8624.9700d0, -10051.5300d0, 4400.0670d0, -655.0200d0, &
         88.2830d0, 444.5170d0, -1863.8640d0, 2095.2880d0, -901.7880d0, 132.9850d0], SHAPE(a1)), &
         a2 = RESHAPE([&
         518.1021d0, -734.8666d0, 1021.1775d0, -479.0721d0, 93.1373d0, -6.4285d0, &
         473.2636d0, 1443.4137d0, -1977.3395d0, 922.3575d0, -178.9275d0, 12.3600d0, &
         -482.2089d0, -737.1616d0, 1096.8827d0, -521.1341d0, 101.7963d0, - 7.0571d0, &
         115.5291d0, 169.6374d0, -245.6490d0, 114.2430d0, -21.9972d0, 1.5097d0, &
         0d0, 0d0, 0d0, 0d0, 0d0, 0d0,&
         0d0, 0d0, 0d0, 0d0, 0d0, 0d0], SHAPE(a2))
    REAL(REAL64), PARAMETER :: cm_to_micron = 1d4
    REAL(REAL64) :: theta, lambda, f, a(6,6)

    coeff(:,:) = 0d0

    ks = 1
    ke = UBOUND(grd,1)
    js = 1
    je = UBOUND(temp,1)

    DO j = js, je
       theta = 5040d0 / temp(j)
       ! Following Phoenix, extend the temperature range.
       !      IF(theta < 0.5d0 .OR. theta > 3.6d0) CYCLE

       DO k = ks, ke
          lambda = 1d0 / grd(k) * cm_to_micron

          ! Following Phoenix, extrapolate to shorter wavelengths.
          lambda = MAX(lambda, 0.1823d0)
          IF(lambda < 0.3645d0) THEN
             a = a2
          ELSE
             a = a1
          END IF

          f = 0
          DO n = LBOUND(a,2), UBOUND(a,2)
             f = f + theta**((DBLE(n) + 1d0) / 2d0) * &
                  (lambda**2 * a(1,n) + a(2,n) + a(3,n) / lambda + a(4,n) / lambda**2 + &
                  a(5,n) / lambda**3 + a(6,n) / lambda**4)
          END DO
          ! f can be negative when temp(j) is as low as 300K.
          coeff(k,j) = MAX(1d-29 * f, 0d0)
       END DO
    END DO

    RETURN
    
  END SUBROUTINE brems_hm_john_1988


  SUBROUTINE brems_h2m_john_1975(temp, grd, coeff)
    ! coefficients of bremsstrahlung of H2-
    ! John 1975, MNRAS, 172, 305, eq(3)
    USE const_module, ONLY : alpha, hbar, pi, clight, k_bol, m_ele
    IMPLICIT NONE
    REAL(REAL64), INTENT(IN) :: temp(:)
    REAL(REAL64), INTENT(IN) :: grd(:)
    REAL(REAL64), INTENT(OUT) :: coeff(:,:)

    INTEGER(INT32) :: ks, ke, k, n, js, je, j
    REAL(REAL64), PARAMETER, DIMENSION(4,3) :: &
         a = RESHAPE([&
         0.4612d0, 0.2403d0, 0.03423d0, -0.06184d0, &
         0.006916d0, -1.0119d0, 0.9102d0, -0.2181d0, &
         -0.03981d0, 0.3508d0, -0.3048d0, 0.07262d0], SHAPE(a))
    REAL(REAL64), PARAMETER :: q = 1d0         
    REAL(REAL64), PARAMETER :: cm_to_micron = 1d4
    REAL(REAL64) :: theta, lambda, f

    coeff(:,:) = 0d0

    ks = 1
    ke = UBOUND(grd,1)
    js = 1
    je = UBOUND(temp,1)

    DO j = js, je
       IF(temp(j) < 100d0 .OR. temp(j) > 15000d0) CYCLE
       theta = 5040.2d0 / temp(j)
       DO k = ks, ke
          lambda = 1d0 / grd(k) * cm_to_micron
          IF(lambda < 0.5d0 .OR. lambda > 10d0) CYCLE
          
          f = 0
          DO n = LBOUND(a,2), UBOUND(a,2)
             f = f + theta**((DBLE(n) + q) / 2d0) * &
                  (lambda**2 * a(1,n) + a(2,n) + a(3,n) / lambda + a(4,n) / lambda**2)
          END DO
          coeff(k,j) = MAX(1d-26 * f, 0d0)
       END DO
    END DO
    
    RETURN
    
  END SUBROUTINE brems_h2m_john_1975



  SUBROUTINE brems_h2m_john_1978(temp, grd, coeff)
    ! coefficients of bremsstrahlung of H2-
    ! John 1978, A&A, 67, 395, eq.7
    USE const_module, ONLY : alpha, hbar, pi, clight, k_bol, m_ele
    IMPLICIT NONE
    REAL(REAL64), INTENT(IN) :: temp(:)
    REAL(REAL64), INTENT(IN) :: grd(:)
    REAL(REAL64), INTENT(OUT) :: coeff(:,:)

    INTEGER(INT32) :: ks, ke, k, n, js, je, j
    REAL(REAL64), PARAMETER, DIMENSION(0:2) :: a = [-9.451d-4, 4.750d-3,-5.468d-4]
    REAL(REAL64), PARAMETER, DIMENSION(0:2) :: b = [ 6.556d-2,-1.104d-1, 6.784d-2]
    REAL(REAL64), PARAMETER, DIMENSION(0:2) :: c = [ 1.712d0 ,-5.033d0 , 3.691d0 ]
    REAL(REAL64), PARAMETER, DIMENSION(0:2) :: d = [-6.416d0 , 1.769d1 ,-1.223d1 ]
    REAL(REAL64), PARAMETER :: cm_to_micron = 1d4
    REAL(REAL64) :: theta, lambda, f, dk2

    coeff(:,:) = 0d0

    ks = 1
    ke = UBOUND(grd,1)
    js = 1
    je = UBOUND(temp,1)

    DO j = js, je
       theta = 5040d0 / temp(j)
       DO k = ks, ke
          lambda = 1d0 / grd(k) * cm_to_micron
          dk2 = 0.09113d0 / lambda
          IF(dk2 > 0.15d0) CYCLE
          
          f = 0
          DO n = 0, 2
             f = f + theta**(DBLE(n) / 2d0) * &
                  (dk2**(-2) * a(n) + b(n) + c(n) * dk2 + c(n) * dk2**2)
          END DO
          coeff(k,j) = MAX(1d-26 * f, 0d0)
       END DO
    END DO
    
    RETURN
    
  END SUBROUTINE brems_h2m_john_1978



  SUBROUTINE coeff_cia(temp, np, grd, alp)
    USE code_module, ONLY : HI, HeI, H2, N2, H2O, CO2, CH4, O2
    IMPLICIT NONE
    
    REAL(REAL64), INTENT(IN) :: temp(:)
    REAL(REAL64), INTENT(IN) :: np(:,:)
    REAL(REAL64), INTENT(IN) :: grd(:)
    REAL(REAL64), INTENT(OUT) :: alp(:,:)

    REAL(REAL64), ALLOCATABLE :: out(:,:)

    INTEGER(INT32), PARAMETER :: &
         cia_h2_h2 = 1, &
         cia_h2_he = 1, &
         cia_n2_n2 = 1, &
         cia_ch4_he = 1, &
         cia_co2_ch4 = 1, &
         cia_co2_co2 = 1, &
         cia_co2_h2 = 1, &
         cia_co2_he = 1, &
         cia_h2_ch4 = 1, &
         cia_h2_h = 1, &
         cia_he_h = 1, &
         cia_n2_h2o = 1, &
         cia_n2_h2 = 1, &
         cia_n2_he = 1, &
         cia_o2_co2 = 1, &
         cia_o2_o2 = 1, &
         cia_o2_n2 = 1
    INTEGER(INT32) :: j, js, je

    alp(:,:) = 0d0

    js = LBOUND(temp,1)
    je = UBOUND(temp,1)
    ALLOCATE(out(UBOUND(grd,1),js:je))
    
    IF(CIA_H2_H2 == 1) THEN
       CALL coeff_cia_hitran('H2-H2', temp(:), grd(:), out(:,:))
       DO j = js, je
          alp(:,j) = alp(:,j) + (np(H2,j) * np(H2,j)) * out(:,j)
       END DO
    END IF
    IF(CIA_H2_He == 1) THEN
       CALL coeff_cia_hitran('H2-He', temp(:), grd(:), out(:,:))
       DO j = js, je
          alp(:,j) = alp(:,j) + (np(H2,j) * np(HeI,j)) * out(:,j)
       END DO
    END IF
    IF(CIA_N2_N2 == 1) THEN
       CALL coeff_cia_hitran('N2-N2', temp(:), grd(:), out(:,:))
       DO j = js, je
          alp(:,j) = alp(:,j) + (np(N2,j) * np(N2,j)) * out(:,j)
       END DO
    END IF
    IF(CIA_CH4_He == 1) THEN
       CALL coeff_cia_hitran('CH4-He', temp(:), grd(:), out(:,:))
       DO j = js, je
          alp(:,j) = alp(:,j) + (np(CH4,j) * np(HeI,j)) * out(:,j)
       END DO
    END IF
    IF(CIA_CO2_CH4 == 1) THEN
       CALL coeff_cia_hitran('CO2-CH4', temp(:), grd(:), out(:,:))
       DO j = js, je
          alp(:,j) = alp(:,j) + (np(CO2,j) * np(CH4,j)) * out(:,j)
       END DO
    END IF
    IF(CIA_CO2_CO2 == 1) THEN
       CALL coeff_cia_hitran('CO2-CO2', temp(:), grd(:), out(:,:))
       DO j = js, je
          alp(:,j) = alp(:,j) + (np(CO2,j) * np(CO2,j)) * out(:,j)
       END DO
    END IF
    IF(CIA_CO2_H2 == 1) THEN
       CALL coeff_cia_hitran('CO2-H2', temp(:), grd(:), out(:,:))
       DO j = js, je
          alp(:,j) = alp(:,j) + (np(CO2,j) * np(H2,j)) * out(:,j)
       END DO
    END IF
    IF(CIA_CO2_He == 1) THEN
       CALL coeff_cia_hitran('CO2-He', temp(:), grd(:), out(:,:))
       DO j = js, je
          alp(:,j) = alp(:,j) + (np(CO2,j) * np(HeI,j)) * out(:,j)
       END DO
    END IF
    IF(CIA_H2_CH4 == 1) THEN
       CALL coeff_cia_hitran('H2-CH4', temp(:), grd(:), out(:,:))
       DO j = js, je
          alp(:,j) = alp(:,j) + (np(H2,j) * np(CH4,j)) * out(:,j)
       END DO
    END IF
    IF(CIA_H2_H == 1) THEN
       CALL coeff_cia_hitran('H2-H', temp(:), grd(:), out(:,:))
       DO j = js, je
          alp(:,j) = alp(:,j) + (np(H2,j) * np(HI,j)) * out(:,j)
       END DO
    END IF
    IF(CIA_He_H == 1) THEN
       CALL coeff_cia_hitran('He-H', temp(:), grd(:), out(:,:))
       DO j = js, je
          alp(:,j) = alp(:,j) + (np(HeI,j) * np(HI,j)) * out(:,j)
       END DO
    END IF
    IF(CIA_N2_H2O == 1) THEN
       CALL coeff_cia_hitran('N2-H2O', temp(:), grd(:), out(:,:))
       DO j = js, je
          alp(:,j) = alp(:,j) + (np(N2,j) * np(H2O,j)) * out(:,j)
       END DO
    END IF
    IF(CIA_N2_H2 == 1) THEN
       CALL coeff_cia_hitran('N2-H2', temp(:), grd(:), out(:,:))
       DO j = js, je
          alp(:,j) = alp(:,j) + (np(N2,j) * np(H2,j)) * out(:,j)
       END DO
    END IF
    IF(CIA_N2_He == 1) THEN
       CALL coeff_cia_hitran('N2-He', temp(:), grd(:), out(:,:))
       DO j = js, je
          alp(:,j) = alp(:,j) + (np(N2,j) * np(HeI,j)) * out(:,j)
       END DO
    END IF
    IF(CIA_O2_CO2 == 1) THEN
       CALL coeff_cia_hitran('O2-CO2', temp(:), grd(:), out(:,:))
       DO j = js, je
          alp(:,j) = alp(:,j) + (np(O2,j) * np(CO2,j)) * out(:,j)
       END DO
    END IF
    IF(CIA_O2_O2 == 1) THEN
       CALL coeff_cia_hitran('O2-O2', temp(:), grd(:), out(:,:))
       DO j = js, je
          alp(:,j) = alp(:,j) + (np(O2,j) * np(O2,j)) * out(:,j)
       END DO
    END IF
    IF(CIA_O2_N2 == 1) THEN
       CALL coeff_cia_hitran('O2-N2', temp(:), grd(:), out(:,:))
       DO j = js, je
          alp(:,j) = alp(:,j) + (np(O2,j) * np(N2,j)) * out(:,j)
       END DO
    END IF

    DEALLOCATE(out)
    
    RETURN
    
  END SUBROUTINE coeff_cia

  
  
  SUBROUTINE coeff_cia_hitran(species, temp, grd, coeff)
    ! coefficients for CIA
    ! https://hitran.org/cia/
    USE HDF5
    USE H5LT
    IMPLICIT NONE
    CHARACTER(*), INTENT(IN) :: species
    REAL(REAL64), INTENT(IN) :: temp(:)
    REAL(REAL64), INTENT(IN) :: grd(:)
    REAL(REAL64), INTENT(OUT) :: coeff(:,:)

    INTEGER(HID_T) :: file_id, group_id
    INTEGER(HSIZE_T) :: dims(1), dims2(2)
    INTEGER(INT32) :: error, num_wn, num_te, j
    REAL(REAL64), ALLOCATABLE :: coeff0(:,:), wn(:), te(:)
    LOGICAL :: exists

    coeff(:,:) = 0d0
    
    CALL h5open_f(error)
    CALL h5Fopen_f('input/h5/CIA_HITRAN.h5', H5F_ACC_RDONLY_F, file_id, error)
    
    CALL h5Lexists_f(file_id, TRIM(species), exists, error)
    IF(exists) THEN
!       if(myrk == 0) PRINT *, 'collision-induced absorption: '//TRIM(species)
       CALL h5Gopen_f(file_id, TRIM(species), group_id, error)
       dims = [1]
       CALL h5LTread_dataset_f(group_id, 'num_wn', H5T_STD_I32LE, num_wn, dims, error)
       CALL h5LTread_dataset_f(group_id, 'num_temp', H5T_STD_I32LE, num_te, dims, error)
       ALLOCATE(coeff0(num_te, num_wn), wn(num_wn), te(num_te))
       dims = [num_te]
       CALL h5LTread_dataset_f(group_id, 'temp', H5T_IEEE_F64LE, te, dims, error)
       dims = [num_wn]
       CALL h5LTread_dataset_f(group_id, 'wn', H5T_IEEE_F64LE, wn, dims, error)
       dims2 = [num_te, num_wn]
       CALL h5LTread_dataset_f(group_id, 'coeff', H5T_IEEE_F64LE, coeff0, dims2, error)
       DO j = LBOUND(temp,1), UBOUND(temp,1)
          CALL interp2d_yarr(temp(j), grd, te, wn, coeff0, coeff(:,j))
       END DO
       CALL h5Gclose_f(group_id, error)
       DEALLOCATE(coeff0, wn, te)
    ELSE
       if(myrk == 0) PRINT *, '+++ WARNING: NO DATA FOR '//TRIM(species)//' +++'
    END IF

    CALL h5Fclose_f(file_id, error)
    CALL h5close_f(error)

    RETURN
  END SUBROUTINE coeff_cia_hitran


  SUBROUTINE photo_topbase(temp, grd, np, alp, photoion_verner)

    USE const_module, ONLY : c2
    IMPLICIT NONE
    
    REAL(REAL64), INTENT(IN) :: temp(:) ! TEMPERATUER [K]
    REAL(REAL64), INTENT(IN) :: grd(:)  ! WAVENUMBER GRID GRID [cm^-1]
    REAL(REAL64), INTENT(IN) :: np(:,:) ! PARTIAL NUMBER DENSITY [cm^-3]
    REAL(REAL64), INTENT(OUT) :: alp(:,:) ! ABSORPTION COEFFICIENT [cm^-1]
    INTEGER(INT32), INTENT(IN) :: photoion_verner

    ! HDF5 variables
    INTEGER(HID_T) :: file_spec, grp_spec, grp_lvl
    INTEGER(HSIZE_T), DIMENSION(1) :: dims
    
    REAL(REAL64), ALLOCATABLE :: crs(:,:), crs0(:)
    REAL(REAL64), ALLOCATABLE :: wnum(:), cross(:), ene(:)
    INTEGER(INT32), ALLOCATABLE :: gtot(:)

    INTEGER(INT32) :: n_min, n, nlevel, error, npo, npo_max, nz, ne
    INTEGER(INT32) :: j, js, je, k, ks, ke
    REAL(REAL64) :: zeta

    IF(photoion_verner == 0) THEN
#ifdef CROSS_CHECK
       n_min = 2
#else    
       n_min = 1
#endif
    ELSE
       n_min = 2
    END IF
    
    js = LBOUND(temp,1)
    je = UBOUND(temp,1)
    ks = LBOUND(grd,1)
    ke = UBOUND(grd,1)

    alp(:,:) = 0d0

    ALLOCATE(crs(LBOUND(grd,1):UBOUND(grd,1),js:je), crs0(LBOUND(grd,1):UBOUND(grd,1)))

    CALL h5open_f(error)
    CALL h5Fopen_f('input/h5/TOPbase.h5', H5F_ACC_RDONLY_F, file_spec, error)
    CALL h5Eset_auto_f(0, error)

    DO nz = 1, 26
       IF(nz==15 .OR. nz==17 .OR. nz==19 .OR. (nz >=21 .AND. nz <=25)) CYCLE ! NO DATA IN TOPBASE
       DO ne = 1, nz
          crs(:,:) = 0d0
          CALL h5Gopen_f(file_spec, c02(nz)//'.'//c02(ne), grp_spec, error)
          IF(error /= 0) CYCLE
          dims = [1]
          CALL h5LTread_dataset_f(grp_spec, 'nlevel', H5T_STD_I32LE, nlevel, dims, error)
          ALLOCATE(ene(nlevel), gtot(nlevel))
          dims = [nlevel]
          CALL h5LTread_dataset_f(grp_spec, 'ene', H5T_IEEE_F64LE, ene, dims, error)
          CALL h5LTread_dataset_f(grp_spec, 'gtot', H5T_STD_I32LE, gtot, dims, error)
          DO n = n_min, nlevel
             CALL h5Gopen_f(grp_spec, c04(n), grp_lvl, error)
             IF(error /= 0) THEN
                PRINT *, error, c04(n)
                STOP
             END IF
             dims = [1]
             CALL h5LTread_dataset_f(grp_lvl, 'np', H5T_STD_I32LE, npo_max, dims, error)
             IF(npo_max == 0) THEN
                CALL h5Gclose_f(grp_lvl, error)
                CYCLE
             END IF
             ALLOCATE(wnum(npo_max), cross(npo_max))
             dims = [npo_max]
             CALL h5LTread_dataset_f(grp_lvl, 'wnum', H5T_IEEE_F64LE, wnum, dims, error)
             CALL h5LTread_dataset_f(grp_lvl, 'cross', H5T_IEEE_F64LE, cross, dims, error)
             CALL h5Gclose_f(grp_lvl, error)
             crs0(:) = 0d0
             DO k = ks, ke
                IF(grd(k) < wnum(1) .OR. grd(k) >= wnum(npo_max)) CYCLE
                DO npo = 1, npo_max
                   IF(wnum(npo) > grd(k)) EXIT
                END DO
                crs0(k) = ((wnum(npo) - grd(k)) * cross(npo-1) + (grd(k) - wnum(npo-1)) * cross(npo)) / &
                     (wnum(npo) - wnum(npo-1))
             END DO
             DO j = js, je
                crs(:,j) = crs(:,j) + crs0(:) * gtot(n) * EXP(-c2 * ene(n) / temp(j))
             END DO
             DEALLOCATE(wnum, cross)
          END DO
          DO j = js, je
             zeta = 0d0
             DO n = 1, nlevel
                zeta = zeta + gtot(n) * EXP(-c2 * ene(n) / temp(j))
             END DO
             alp(:,j) = alp(:,j) + (1d0 - EXP(-c2 * grd(:) / temp(j))) * np(nz*100+(nz-ne),j) / zeta * crs(:,j)
          END DO
          CALL h5Gclose_f(grp_spec, error)
          DEALLOCATE(ene, gtot)
       END DO
    END DO
    
    CALL h5Fclose_f(file_spec, error)
    
    DEALLOCATE(crs, crs0)

    CALL h5close_f(error)

    RETURN
  END SUBROUTINE photo_topbase

  FUNCTION c02(n)
    CHARACTER(LEN=2) :: c02
    INTEGER(INT32), INTENT(IN) :: n
    
    WRITE(c02,'(I2.2)') n
    
    RETURN
  END FUNCTION c02
  
  FUNCTION c04(n)
    CHARACTER(LEN=4) :: c04
    INTEGER(INT32), INTENT(IN) :: n
    
    WRITE(c04,'(I4.4)') n
    
    RETURN
  END FUNCTION c04

  
  SUBROUTINE photo_mathisen(temp, grd, np, zeta, alp, photoion_verner)

    USE code_module, ONLY : HI, HeI, HeII, CI, NI, OI, NaI, MgI, MgII, AlI, SiI, SiII, SI, CaI, CaII, FeI
    USE const_module, ONLY : c2
    IMPLICIT NONE
    
    REAL(REAL64), INTENT(IN) :: temp(:)
    REAL(REAL64), INTENT(IN) :: grd(:)
    REAL(REAL64), INTENT(IN) :: np(:,:)
    REAL(REAL64), INTENT(IN) :: zeta(:,:)
    REAL(REAL64), INTENT(OUT) :: alp(:,:)
    INTEGER(INT32), INTENT(IN) :: photoion_verner

    REAL(REAL64), ALLOCATABLE :: crs(:,:)

    INTEGER(INT32) :: n_min
    INTEGER(INT32) :: j, js, je
    INTEGER(INT32), PARAMETER :: &
         photo_HI_mathisen = 1, &
         photo_HeI_mathisen = 1, &
         photo_HeII_mathisen = 1, &
         photo_CI_mathisen = 1, &
         photo_NI_mathisen = 1, &
         photo_OI_mathisen = 1, &
         photo_NaI_mathisen = 1, &
         photo_MgI_mathisen = 1, &
         photo_MgII_mathisen = 1, &
         photo_AlI_mathisen = 1, &
         photo_SiI_mathisen = 1, &
         photo_SiII_mathisen = 1, &
         photo_SI_mathisen = 1, &
         photo_CaI_mathisen = 1, &
         photo_CaII_mathisen = 1, &
         photo_FeI_mathisen = 1

    IF(photoion_verner == 0) THEN
#ifdef CROSS_CHECK
       n_min = 2
#else    
       n_min = 1
#endif
    ELSE
       n_min = 2
    END IF
    
    js = LBOUND(temp,1)
    je = UBOUND(temp,1)

    alp(:,:) = 0d0

    ALLOCATE(crs(LBOUND(grd,1):UBOUND(grd,1),js:je))

    IF(photo_HI_mathisen == 1) THEN
       CALL cross_hydrogenic(temp(:), grd(:), crs(:,:), 1, n_min)
       DO j = js, je
          alp(:,j) = alp(:,j) + (1d0 - EXP(-c2 * grd(:) / temp(j))) * np(HI,j) / zeta(j,HI) * crs(:,j)
       END DO
    END IF
    IF(photo_HeI_mathisen == 1) THEN
       CALL cross_HeI_mathisen(temp(:), grd(:), crs(:,:), n_min)
       DO j = js, je
          alp(:,j) = alp(:,j) + (1d0 - EXP(-c2 * grd(:) / temp(j))) * np(HeI,j) / zeta(j,HeI) * crs(:,j)
       END DO
    END IF
    IF(photo_HeII_mathisen == 1) THEN
       CALL cross_hydrogenic(temp(:), grd(:), crs(:,:), 2, n_min)
       DO j = js, je
          alp(:,j) = alp(:,j) + (1d0 - EXP(-c2 * grd(:) / temp(j))) * np(HeII,j) / zeta(j,HeII) * crs(:,j)
       END DO
    END IF
    IF(photo_CI_mathisen == 1) THEN
       CALL cross_CI_mathisen(temp(:), grd(:), crs(:,:), n_min)
       DO j = js, je
          alp(:,j) = alp(:,j) + (1d0 - EXP(-c2 * grd(:) / temp(j))) * np(CI,j) / zeta(j,CI) * crs(:,j)
       END DO
    END IF
    IF(photo_NI_mathisen == 1) THEN
       CALL cross_NI_mathisen(temp(:), grd(:), crs(:,:), n_min)
       DO j = js, je
          alp(:,j) = alp(:,j) + (1d0 - EXP(-c2 * grd(:) / temp(j))) * np(NI,j) / zeta(j,NI) * crs(:,j)
       END DO
    END IF
    IF(photo_OI_mathisen == 1) THEN
       CALL cross_OI_mathisen(temp(:), grd(:), crs(:,:), n_min)
       DO j = js, je
          alp(:,j) = alp(:,j) + (1d0 - EXP(-c2 * grd(:) / temp(j))) * np(OI,j) / zeta(j,OI) * crs(:,j)
       END DO
    END IF
    IF(photo_NaI_mathisen == 1) THEN
       CALL cross_NaI_mathisen(temp(:), grd(:), crs(:,:), n_min)
       DO j = js, je
          alp(:,j) = alp(:,j) + (1d0 - EXP(-c2 * grd(:) / temp(j))) * np(NaI,j) / zeta(j,NaI) * crs(:,j)
       END DO
    END IF
    IF(photo_MgI_mathisen == 1) THEN
       CALL cross_MgI_mathisen(temp(:), grd(:), crs(:,:), n_min)
       DO j = js, je
          alp(:,j) = alp(:,j) + (1d0 - EXP(-c2 * grd(:) / temp(j))) * np(MgI,j) / zeta(j,MgI) * crs(:,j)
       END DO
    END IF
    IF(photo_MgII_mathisen == 1) THEN
       CALL cross_MgII_mathisen(temp(:), grd(:), crs(:,:), n_min)
       DO j = js, je
          alp(:,j) = alp(:,j) + (1d0 - EXP(-c2 * grd(:) / temp(j))) * np(MgII,j) / zeta(j,MgII) * crs(:,j)
       END DO
    END IF
    IF(photo_AlI_mathisen == 1) THEN
       CALL cross_AlI_mathisen(temp(:), grd(:), crs(:,:), n_min)
       DO j = js, je
          alp(:,j) = alp(:,j) + (1d0 - EXP(-c2 * grd(:) / temp(j))) * np(AlI,j) / zeta(j,AlI) * crs(:,j)
       END DO
    END IF
    IF(photo_SiI_mathisen == 1) THEN
       CALL cross_SiI_mathisen(temp(:), grd(:), crs(:,:), n_min)
       DO j = js, je
          alp(:,j) = alp(:,j) + (1d0 - EXP(-c2 * grd(:) / temp(j))) * np(SiI,j) / zeta(j,SiI) * crs(:,j)
       END DO
    END IF
    IF(photo_SiII_mathisen == 1) THEN
       CALL cross_SiII_mathisen(temp(:), grd(:), crs(:,:), n_min)
       DO j = js, je
          alp(:,j) = alp(:,j) + (1d0 - EXP(-c2 * grd(:) / temp(j))) * np(SiII,j) / zeta(j,SiII) * crs(:,j)
       END DO
    END IF
    IF(photo_SI_mathisen == 1) THEN
       CALL cross_SI_mathisen(temp(:), grd(:), crs(:,:), n_min)
       DO j = js, je
          alp(:,j) = alp(:,j) + (1d0 - EXP(-c2 * grd(:) / temp(j))) * np(SI,j) / zeta(j,SI) * crs(:,j)
       END DO
    END IF
    IF(photo_CaI_mathisen == 1) THEN
       CALL cross_CaI_mathisen(temp(:), grd(:), crs(:,:), n_min)
       DO j = js, je
          alp(:,j) = alp(:,j) + (1d0 - EXP(-c2 * grd(:) / temp(j))) * np(CaI,j) / zeta(j,CaI) * crs(:,j)
       END DO
    END IF
    IF(photo_CaII_mathisen == 1) THEN
       CALL cross_CaII_mathisen(temp(:), grd(:), crs(:,:), n_min)
       DO j = js, je
          alp(:,j) = alp(:,j) + (1d0 - EXP(-c2 * grd(:) / temp(j))) * np(CaII,j) / zeta(j,CaII) * crs(:,j)
       END DO
    END IF
    IF(photo_FeI_mathisen == 1) THEN
       CALL cross_FeI_mathisen(temp(:), grd(:), crs(:,:), n_min)
       DO j = js, je
          alp(:,j) = alp(:,j) + (1d0 - EXP(-c2 * grd(:) / temp(j))) * np(FeI,j) / zeta(j,FeI) * crs(:,j)
       END DO
    END IF
    
    DEALLOCATE(crs)

    RETURN
  END SUBROUTINE photo_mathisen

  
  SUBROUTINE cross_HI_mathisen(temp, grd, crs, n_min)
    USE const_module, ONLY : cm_to_ang, ev_to_erg, k_bol, pi, m_ele, e2, clight, h
    IMPLICIT NONE

    REAL(REAL64), INTENT(IN) :: temp(:)
    REAL(REAL64), INTENT(IN) :: grd(:)
    REAL(REAL64), INTENT(OUT) :: crs(:,:)
    INTEGER(INT32), INTENT(IN) :: n_min

    INTEGER(INT32), PARAMETER :: n_max = 15
    REAL(REAL64), DIMENSION(n_max) :: g, chi, crs0, wa
    REAL(REAL64), PARAMETER :: z = 1
    REAL(REAL64), PARAMETER :: fac0 = 64d0 * pi**3 / (3d0 * SQRT(3d0)) * (m_ele * e2**5 * clight**2) / h**6
    REAL(REAL64) :: lambda
    INTEGER(INT32) :: k, j, n

    DO n = 1, n_max
       g(n) = 2d0 * DBLE(n)**2
       chi(n) = 13.595d0 * (1d0 - 1d0 / DBLE(n)**2) * ev_to_erg
       wa(n) = 911.8d0 * DBLE(n)**2
    END DO
    
    DO k = LBOUND(grd,1), UBOUND(grd,1)
       lambda = (1d0 / grd(k)) * cm_to_ang
       DO n = 1, n_max
#ifdef KARZAS_IN_PLACE_OF_MATHISEN          
          crs0(n) = fac0 * (z**4 / DBLE(n)**5) / grd(k)**3 * gaunt_karzas(grd(k)/z**2,n)
#else          
          crs0(n) = fac0 * (z**4 / DBLE(n)**5) * (lambda * 1d-8)**3 * gaunt_bf(lambda*z**2,n)
#endif          
       END DO

       DO j = LBOUND(temp,1), UBOUND(temp,1)
          crs(k,j) = 0d0
          DO n = n_min, n_max
             IF(lambda >= wa(n)) CYCLE
             crs(k,j) = crs(k,j) + crs0(n) * g(n) * EXP(-chi(n)/ (k_bol * temp(j)))
          END DO
       END DO
    END DO

    RETURN
  END SUBROUTINE cross_HI_mathisen


  FUNCTION gaunt_karzas(wn,n)
    USE const_module, ONLY : clight, h, rydberg, e2, m_ele
    IMPLICIT NONE
    REAL(REAL64) :: gaunt_karzas
    REAL(REAL64), INTENT(IN) :: wn
    INTEGER(INT32), INTENT(IN) :: n
    INTEGER(INT32), PARAMETER :: n_max = 15, imax1 = 51, imax2 = 37, imax3 = 29
    INTEGER(INT32), PARAMETER, DIMENSION(n_max) :: imax = [&
         imax1, imax1, imax1, imax1, imax1, imax1, imax1, imax1, imax2, imax1, imax3, imax3, imax3, imax3, imax3]
    REAL(REAL64), SAVE :: ele(imax1) = 0d0, pho(imax1,n_max) = 0d0, gaunt(imax1,n_max) = 0d0
    REAL(REAL64) :: x
    LOGICAL, SAVE :: first = .TRUE.
    INTEGER(INT32) :: i, nn

    IF(first) THEN
       OPEN(1,FILE='input/Karzas_Latter_1961.tsv',STATUS='OLD')
       READ(1,*) !SKIP
       READ(1,*) !SKIP
       READ(1,*) !SKIP
       READ(1,*) !SKIP
       READ(1,*) !SKIP
       DO i = 1, imax1 ! n= 1 - 8, 10
          READ(1,*) ele(i), (pho(i,nn), gaunt(i,nn), nn = 1, 8), (pho(i,nn), gaunt(i,nn), nn = 10, 10)
       END DO
       READ(1,*) !SKIP
       DO i = 1, imax2 ! n = 9
          READ(1,*) ele(i), (pho(i,nn), gaunt(i,nn), nn = 9, 9)
       END DO
       READ(1,*) !SKIP
       DO i = 1, imax3 ! n = 11 - 15
          READ(1,*) ele(i), (pho(i,nn), gaunt(i,nn), nn = 11, 15)
       END DO
       CLOSE(1)
       DO nn = 1, n_max
          DO i = 1, imax(nn)-1
             IF(pho(i,nn) == pho(imax(nn),nn)) THEN
                pho(i,nn) = pho(imax(nn),nn) + ele(i)
             END IF
          END DO
       END DO
       first = .FALSE.
    END IF

    x = h * clight * wn / rydberg
    gaunt_karzas = sig_interp(x, pho(1:imax(n),n), gaunt(1:imax(n),n))
    
    RETURN
  END FUNCTION gaunt_karzas


  SUBROUTINE cross_HeI_mathisen(temp, grd, crs, n_min)
    USE const_module, ONLY : cm_to_ang, ev_to_erg, k_bol, megabarn_to_cm2
    IMPLICIT NONE

    REAL(REAL64), INTENT(IN) :: temp(:)
    REAL(REAL64), INTENT(IN) :: grd(:)
    REAL(REAL64), INTENT(OUT) :: crs(:,:)
    INTEGER(INT32), INTENT(IN) :: n_min

    INTEGER(INT32), PARAMETER :: n_max = 10
    REAL(REAL64), PARAMETER, DIMENSION(n_max) :: &
         wa = [&
         504.3d0, 2600.1d0, 3121.2d0, 3421.2d0, 3679.0d0, &
         6631.4d0, 7434.0d0, 7842.1d0, 8187.3d0, 8259.6d0], &
         chi = [&
         0.000d0, 19.819d0, 20.615d0, 20.964d0, 21.218d0, &
         22.718d0, 22.920d0, 23.007d0, 23.073d0,23.087d0] * ev_to_erg
    INTEGER(INT32), PARAMETER, DIMENSION(n_max) :: &
         g = [1, 3, 1, 9, 3, 3, 1, 9, 20, 3]

    INTEGER(INT32), PARAMETER :: n1 = 48, n2 = 19, n3 = 20
    REAL(REAL64), PARAMETER, DIMENSION(n1) :: &
         lam1 = [&
         504.3d0, 501.5d0, 498.7d0, 493.3d0, 488.1d0, 480.3d0, &
         477.8d0, 454.0d0, 443.0d0, 395.0d0, 356.4d0, 348.2d0, &
         324.6d0, 302.0d0, 298.1d0, 275.6d0, 260.6d0, 256.2d0, &
         239.4d0, 224.6d0, 220.0d0, 215.0d0, 210.0d0, 205.0d0, &
         200.0d0, 195.0d0, 190.0d0, 185.0d0, 180.0d0, 175.0d0, &
         170.0d0, 165.0d0, 160.0d0, 155.0d0, 150.0d0, 145.0d0, &
         135.0d0, 130.0d0, 125.0d0, 120.0d0, 115.0d0, 110.0d0, &
         105.0d0, 100.0d0,  95.0d0,  90.0d0,  85.0d0,  80.0d0], &
         sig1 = [&
         7.376d0, 7.317d0, 7.259d0, 7.143d0, 7.030d0, 6.857d0, &
         6.800d0, 6.284d0, 6.041d0, 4.977d0, 4.138d0, 3.961d0, &
         3.474d0, 3.025d0, 2.945d0, 2.522d0, 2.259d0, 2.179d0, &
         1.901d0, 1.684d0, 1.610d0, 1.530d0, 1.450d0, 1.380d0, &
         1.300d0, 1.220d0, 1.140d0, 1.080d0, 1.020d0, 0.961d0, &
         0.903d0, 0.847d0, 0.792d0, 0.738d0, 0.687d0, 0.637d0, &
         0.542d0, 0.497d0, 0.454d0, 0.412d0, 0.373d0, 0.335d0, &
         0.299d0, 0.265d0, 0.233d0, 0.202d0, 0.174d0, 0.147d0]
    REAL(REAL64), PARAMETER, DIMENSION(n2) :: &
         lam2 = [&
         2600.6d0, 2427.0d0, 2275.5d0, 2141.8d0, 2023.0d0, &
         1820.9d0, 1655.5d0, 1517.7d0, 1401.0d0, 1214.3d0, &
         958.8d0,  792.2d0,  674.8d0,  587.8d0,  520.6d0, &
         467.3d0,  423.8d0,  387.7d0,  357.3d0], &
         sig2 = [&
         5.496d0, 5.144d0, 4.808d0, 4.489d0, 4.191d0, 3.655d0, &
         3.199d0, 2.812d0, 2.485d0, 1.968d0, 1.301d0, 0.913d0, &
         0.671d0, 0.512d0, 0.404d0, 0.328d0, 0.275d0, 0.240d0, 0.227d0]
    REAL(REAL64), PARAMETER, DIMENSION(n3) :: &
         lam3 = [&
         3122.6d0, 2875.0d0, 2664.9d0, 2483.3d0, 2325.0d0, &
         2062.0d0, 1852.4d0, 1681.5d0, 1539.5d0, 1317.0d0, &
         1021.7d0,  834.6d0,  705.4d0,  610.9d0,  538.6d0, &
         481.7d0,  435.7d0,  397.6d0,  365.7d0,  338.6d0], &
         sig3 = [&
         9.237d0, 8.172d0, 7.255d0, 6.474d0, 5.803d0, &
         4.724d0, 3.904d0, 3.270d0, 2.770d0, 2.047d0, &
         1.223d0, 0.793d0, 0.545d0, 0.389d0, 0.285d0, &
         0.211d0, 0.157d0, 0.115d0, 0.075d0, 0.026d0]

    REAL(REAL64), DIMENSION(n_max) :: crs0 = 0d0
    REAL(REAL64) :: lambda, loglam
    INTEGER(INT32) :: k, j, n

    DO k = LBOUND(grd,1), UBOUND(grd,1)
       lambda = (1d0 / grd(k)) * cm_to_ang
       loglam = LOG10(lambda)

       ! 1s2 1S: Stewart 1978a, Marr & West 1976
       crs0(1) = sig_interp(lambda, lam1, sig1)
       ! 2s 3S: Stewart 1979
       crs0(2) = sig_interp(lambda, lam2, sig2)
       ! 2s 1S: Stewart 1978c
       crs0(3) = sig_interp(lambda, lam3, sig3)
       ! 2p 1Po: Gingerich 1964
       crs0(4) = 10d0**(2.90d0 * loglam - 9.01d0) + 10d0**(3.3d0 * loglam - 11.91d0)
       ! 3s 3S: Hunger & van Blerkom 1967
       crs0(5) = 10d0**(3.50d0 * loglam -11.36d0) + 10d0**(3.6d0 * loglam - 13.04d0)
       ! 2p 3Po: Gingerich 1964
       crs0(6) = 10d0**(1.54d0 * loglam - 4.94d0)
       ! 3s 1S: Hunger & van Blerkom 1967
       crs0(7) = 10d0**(1.86d0 * loglam - 6.01d0)
       ! 3p 3Po: Hunger & van Blerkom 1967
       crs0(8) = 10d0**(2.60d0 * loglam - 8.63d0)
       ! 3d 3D,1D: Hunger & van Blerkom 1967
       crs0(9) = 10d0**(3.69d0 * loglam -13.18d0)
       ! 3p 1Po: Hunger & van Blerkom 1967
       crs0(10)= 10d0**(2.89d0 * loglam - 9.86d0)

       crs0 = crs0 * megabarn_to_cm2
       
       DO j = LBOUND(temp,1), UBOUND(temp,1)
          crs(k,j) = 0d0
          DO n = n_min, n_max
             IF(lambda > wa(n)) CYCLE
             crs(k,j) = crs(k,j) + crs0(n) * g(n) * EXP(-chi(n)/ (k_bol * temp(j)))
          END DO
       END DO
    END DO

    RETURN
  END SUBROUTINE cross_HeI_mathisen

  
  SUBROUTINE cross_HeII_mathisen(temp, grd, crs, n_min)
    USE const_module, ONLY : cm_to_ang, ev_to_erg, k_bol, pi, m_ele, e2, clight, h
    IMPLICIT NONE

    REAL(REAL64), INTENT(IN) :: temp(:)
    REAL(REAL64), INTENT(IN) :: grd(:)
    REAL(REAL64), INTENT(OUT) :: crs(:,:)
    INTEGER(INT32), INTENT(IN) :: n_min

    INTEGER(INT32), PARAMETER :: n_max = 15
    REAL(REAL64), DIMENSION(n_max) :: g, chi, crs0, wa
    REAL(REAL64), PARAMETER :: z = 2
    REAL(REAL64), PARAMETER :: fac0 = 64d0 * pi**3 / (3d0 * SQRT(3d0)) * (m_ele * e2**5 * clight**2) / h**6
    REAL(REAL64) :: lambda
    INTEGER(INT32) :: k, j, n

    DO n = 1, n_max
       g(n) = 2d0 * DBLE(n)**2
       chi(n) = 54.416d0 * (1d0 - 1d0 / DBLE(n)**2) * ev_to_erg
       wa(n) = 227.8d0 * DBLE(n)**2
    END DO
    
    DO k = LBOUND(grd,1), UBOUND(grd,1)
       lambda = (1d0 / grd(k)) * cm_to_ang
       DO n = 1, n_max
#ifdef KARZAS_IN_PLACE_OF_MATHISEN          
          crs0(n) = fac0 * (z**4 / DBLE(n)**5) / grd(k)**3 * gaunt_karzas(grd(k)/z**2,n)
#else          
          crs0(n) = fac0 * (z**4 / DBLE(n)**5) * (lambda * 1d-8)**3 * gaunt_bf(lambda*z**2,n)
#endif          
       END DO
       DO j = LBOUND(temp,1), UBOUND(temp,1)
          crs(k,j) = 0d0
          DO n = n_min, n_max
             IF(lambda > wa(n)) CYCLE
             crs(k,j) = crs(k,j) + crs0(n) * g(n) * EXP(-chi(n)/ (k_bol * temp(j)))
          END DO
       END DO
    END DO

    RETURN
  END SUBROUTINE cross_HeII_mathisen

  
  SUBROUTINE cross_hydrogenic(temp, grd, crs, z, n_min)
    USE const_module, ONLY : cm_to_ang, k_bol, rydberg, h, clight, pi, e2, m_ele
    IMPLICIT NONE

    REAL(REAL64), INTENT(IN) :: temp(:)
    REAL(REAL64), INTENT(IN) :: grd(:)
    REAL(REAL64), INTENT(OUT) :: crs(:,:)
    INTEGER(INT32), INTENT(IN) :: z
    INTEGER(INT32), INTENT(IN) :: n_min

    INTEGER(INT32), PARAMETER :: n_max = 15
    REAL(REAL64), PARAMETER :: fac0 = 64d0 * pi**4 / (3d0 * SQRT(3d0)) * (m_ele * e2**5 / (clight**4 * h**6))
    REAL(REAL64), DIMENSION(n_max) :: g, chi, crs0, wa
    REAL(REAL64) :: lambda
    INTEGER(INT32) :: k, j, n

    DO n = 1, n_max
       g(n) = 2d0 * DBLE(n)**2
       chi(n) = rydberg * DBLE(z)**2 * (1d0 - 1d0 / DBLE(n)**2)
       wa(n) = (h * clight / rydberg) * DBLE(n)**2 / DBLE(z)**2
    END DO

    DO k = LBOUND(grd,1), UBOUND(grd,1)
       lambda = (1d0 / grd(k))
       DO n = 1, n_max
#ifdef KARZAS_IN_PLACE_OF_MATHISEN
          crs0(n) = fac0 * (z**4 / DBLE(n)**5) / grd(k)**3 * gaunt_karzas(grd(k)/z**2,n)
#else          
          crs0(n) = fac0 * (z**4 / DBLE(n)**5) * lambda**3 * gaunt_bf((lambda*cm_to_ang)*z**2,n)
#endif          
       END DO
       
       DO j = LBOUND(temp,1), UBOUND(temp,1)
          crs(k,j) = 0d0
          DO n = n_min, n_max
             IF(lambda > wa(n)) CYCLE
             crs(k,j) = crs(k,j) + crs0(n) * g(n) * EXP(-chi(n)/ (k_bol * temp(j)))
          END DO
       END DO
    END DO

    RETURN
  END SUBROUTINE cross_hydrogenic

  
  FUNCTION gaunt_bf(lambda,n)
    IMPLICIT NONE
    REAL(REAL64) :: gaunt_bf
    REAL(REAL64), INTENT(IN) :: lambda
    INTEGER(INT32), INTENT(IN) :: n
    INTEGER(INT32), PARAMETER :: n_max = 15
    REAL(REAL64), PARAMETER :: safe_fac = 3d0
    REAL(REAL64), PARAMETER, DIMENSION(n_max) :: &
         a = [&
         -2.02848d+05, 1.105d0, 1.101d0, 1.101d0, 1.102d0, &
         1.0986d0, -2.06411d+07, -2.91537d+07, -5.25160d+07, -4.82529d+07, &
         -5.74871d+07, -8.81478d+07, -1.30735d+08, -1.51854d+08, -1.90977d+08], &
         b = [&
         -1.31854d+04, -7.922d-05, -3.290d-05, -1.923d-05, -1.304d-05, &
         -9.020d-06, 8.32740d+04, 1.02091d+05, 1.92771d+05, 1.39700d+05, &
         1.61021d+05, 1.95139d+05, 2.47876d+05, 2.74494d+05, 3.36368d+05], &
         c = [&
         -2.06536d-01, 4.536d-09, 1.152d-09, 5.110d-10, 2.638d-10, &
         1.367d-10, 7.79095d-01, 8.01148d-01, 7.46637d-01, 8.35914d-01, &
         8.47916d-01, 8.55662d-01, 8.56825d-01, 8.66258d-01, 8.65170d-01], &
         d = [9.36244d+04, 0.0d+00, 0.0d+00, 0.0d+00, 0.0d+00, &
         0.0d+00, -1.57600d+07, -2.30675d+07, -3.97851d+07, -4.02115d+07, &
         -5.10197d+07, -7.94763d+07, -1.18616d+08, -1.39892d+08, -1.76141d+08], &
         e = [2.21980d+03, 0.0d+00, 0.0d+00, 0.0d+00, 0.0d+00, &
         0.0d+00, 6.02788d+04, 7.63484d+04, 1.33609d+05, 1.09690d+05, &
         1.28543d+05, 1.58058d+05, 2.02013d+05, 2.26509d+05, 2.77724d+05]

    IF(n >= 2 .AND. n <= 6) THEN ! n = 2-6
       gaunt_bf = a(n) + lambda * (b(n) + lambda * c(n))
    ELSE ! n = 1, 7-15
       gaunt_bf = c(n) * (a(n) + lambda * (b(n) + lambda)) / (d(n) + lambda * (e(n) + lambda))
       ! exclude the range where the denominator becomes very small
       IF(lambda <= (-e(n) + SQRT(e(n)**2 - 4d0 * d(n))) / 2d0 * safe_fac) THEN
          gaunt_bf = 0d0
       END IF
    END IF

    RETURN    
  END FUNCTION gaunt_bf

  
  SUBROUTINE cross_CI_mathisen(temp, grd, crs, n_min)
    USE const_module, ONLY : cm_to_ang, ev_to_erg, k_bol, megabarn_to_cm2
    IMPLICIT NONE

    REAL(REAL64), INTENT(IN) :: temp(:)
    REAL(REAL64), INTENT(IN) :: grd(:)
    REAL(REAL64), INTENT(OUT) :: crs(:,:)
    INTEGER(INT32), INTENT(IN) :: n_min

    INTEGER(INT32), PARAMETER :: n_max = 6
    REAL(REAL64), PARAMETER, DIMENSION(n_max) :: &
         wa = [1101.1d0, 1240.3d0, 1444.6d0, 3284.0d0, 3468.0d0, 4894.0d0], &
         chi = [0.000d0, 1.264d0, 2.684d0, 7.485d0, 7.685d0, 8.727d0] * ev_to_erg
    INTEGER(INT32), PARAMETER, DIMENSION(n_max) :: &
         g = [9, 5, 1, 9, 3, 27]

    INTEGER(INT32), PARAMETER :: n1 = 19, n2 = 18, n3 = 13
    REAL(REAL64), PARAMETER, DIMENSION(n1) :: &
         lam1 = [&
         1101.d0, 1013.d0, 911.d0, 881.d0, 869.d0, 860.d0, 842.d0, &
         814.d0,  804.d0, 798.d0, 780.d0, 777.d0, 753.d0, 701.d0, &
         651.d0,  597.d0, 490.d0, 466.d0, 429.d0], &
         sig1 = [&
         17.2d0, 16.7d0, 16.4d0, 17.0d0, 18.2d0, 11.9d0, 14.0d0, &
         14.8d0, 17.7d0, 13.4d0, 16.6d0, 13.2d0, 14.9d0, 12.6d0, &
         10.9d0,  9.8d0,  7.1d0,  6.6d0,  5.3d0]
    REAL(REAL64), PARAMETER, DIMENSION(n2) :: &
         lam2 = [&
         1240.d0, 1171.d0, 1139.d0, 1070.d0, 1044.d0, 1000.d0, &
         930.d0,  911.d0,  895.d0,  893.d0,  828.d0,  759.d0, &
         701.d0,  651.d0,  570.d0,  506.d0,  414.d0,  351.d0], &
         sig2 = [&
         10.5d0,  9.7d0, 10.8d0, 25.0d0, 45.1d0, 35.0d0, 22.6d0, &
         51.8d0, 69.1d0, 33.0d0, 15.6d0, 10.9d0, 10.7d0,  9.7d0, &
         7.7d0,  6.5d0,  5.2d0,  3.9d0]
    REAL(REAL64), PARAMETER, DIMENSION(n3) :: &
         lam3 = [&
         1445.d0, 1302.d0, 1139.d0, 1040.d0, 1013.d0, 997.d0, &
         959.d0,  911.d0,  828.d0,  701.d0,  570.d0, 456.d0, 351.d0], &
         sig3 = [&
         14.3d0, 12.8d0,  8.5d0,  3.6d0, 30.2d0, 139.4d0, &
         39.6d0, 21.6d0, 14.4d0, 10.2d0,  7.2d0,   5.9d0, 3.6d0]

    REAL(REAL64), DIMENSION(n_max) :: crs0
    REAL(REAL64) :: lambda, loglam
    INTEGER(INT32) :: k, j, n

    DO k = LBOUND(grd,1), UBOUND(grd,1)
       lambda = (1d0 / grd(k)) * cm_to_ang
       loglam = LOG10(lambda)

       ! 2p2 3P: Hoffman et al. 1983
       crs0(1) = sig_interp(lambda, lam1, sig1)
       ! 2p2 1D: Hoffman et al. 1983
       crs0(2) = sig_interp(lambda, lam2, sig2)
       ! 2p2 1S: Hoffman et al. 1983
       crs0(3) = sig_interp(lambda, lam3, sig3)
       ! 3s 3Po: Vernazza et al. 1976
       crs0(4) = 0.20d0 * (lambda / 3284d0)**1.2d0
       ! 3s 1Po: Vernazza et al. 1976
       crs0(5) = 1.54d0 * (lambda / 3468d0)**1.2d0
       ! 3p 3D,3S,3P: Vernazza et al. 1976
       crs0(6) = 2.10d0 * (lambda / 4894d0)**1.5d0

       crs0 = crs0 * megabarn_to_cm2
       
       DO j = LBOUND(temp,1), UBOUND(temp,1)
          crs(k,j) = 0d0
          DO n = n_min, n_max
             IF(lambda > wa(n)) CYCLE
             crs(k,j) = crs(k,j) + crs0(n) * g(n) * EXP(-chi(n)/ (k_bol * temp(j)))
          END DO
       END DO
    END DO

    RETURN
  END SUBROUTINE cross_CI_mathisen

  
  SUBROUTINE cross_NI_mathisen(temp, grd, crs, n_min)
    USE const_module, ONLY : cm_to_ang, ev_to_erg, k_bol, megabarn_to_cm2
    IMPLICIT NONE

    REAL(REAL64), INTENT(IN) :: temp(:)
    REAL(REAL64), INTENT(IN) :: grd(:)
    REAL(REAL64), INTENT(OUT) :: crs(:,:)
    INTEGER(INT32), INTENT(IN) :: n_min

    INTEGER(INT32), PARAMETER :: n_max = 3
    REAL(REAL64), PARAMETER, DIMENSION(n_max) :: &
         wa = [853.1d0, 1020.4d0, 1131.4d0], &
         chi = [0.000d0, 2.384d0, 3.576d0] * ev_to_erg
    INTEGER(INT32), PARAMETER, DIMENSION(n_max) :: &
         g = [4, 10, 6]

    INTEGER(INT32), PARAMETER :: n1 = 27, n2 = 16, n3 = 16
    REAL(REAL64), PARAMETER, DIMENSION(n1) :: &
         lam1 = [&
         843.8d0, 734.9d0, 706.4d0, 695.6d0, 690.4d0, 687.7d0, &
         680.1d0, 665.2d0, 650.9d0, 647.4d0, 646.3d0, 637.2d0, &
         621.0d0, 617.8d0, 607.5d0, 562.5d0, 529.8d0, 511.9d0, &
         500.7d0, 492.6d0, 474.6d0, 440.2d0, 416.1d0, 387.8d0, &
         364.5d0, 337.5d0, 303.8d0], &
         sig1 = [&
         12.980d0, 14.331d0, 14.968d0, 16.397d0, 23.826d0, 5.784d0, &
         12.516d0, 13.692d0, 14.854d0, 19.342d0, 7.284d0, 13.566d0, &
         13.516d0, 13.499d0, 13.694d0, 13.044d0, 12.477d0, 12.113d0, &
         11.884d0, 11.713d0, 11.244d0, 10.342d0,  9.633d0,  8.774d0, &
         8.035d0,  7.149d0,  6.046d0]
    REAL(REAL64), PARAMETER, DIMENSION(n2) :: &
         lam2 = [&
         1012.5d0, 911.3d0, 828.4d0, 759.4d0, 701.0d0, 690.4d0, &
         650.9d0, 607.5d0, 569.5d0, 570.7d0, 479.6d0, 433.9d0, &
         405.0d0, 364.5d0, 331.4d0, 303.8d0], &
         sig2 = [&
         5.3d0,  6.8d0, 13.8d0, 14.3d0, 29.5d0, 33.2d0, &
         16.2d0, 12.8d0, 12.0d0, 11.0d0, 10.1d0, 8.9d0, &
         8.2d0,  7.1d0, 6.2d0,  5.5d0]
    REAL(REAL64), PARAMETER, DIMENSION(n3) :: &
         lam3 = [&
         1104.6d0, 1012.5d0, 959.2d0, 867.9d0, 810.0d0, &
         759.4d0,  743.9d0, 701.0d0, 650.9d0, 607.5d0, &
         552.3d0,  492.6d0, 444.5d0, 414.2d0, 364.5d0, 303.8d0], &
         sig3 = [&
         7.9d0,  6.0d0, 10.2d0,  9.3d0, 15.0d0, &
         20.4d0, 28.6d0, 18.2d0, 14.0d0, 13.5d0, &
         11.4d0, 10.4d0,  9.3d0,  8.5d0, 8.2d0,  5.5d0]

    REAL(REAL64), DIMENSION(n_max) :: crs0
    REAL(REAL64) :: lambda, loglam
    INTEGER(INT32) :: k, j, n

    DO k = LBOUND(grd,1), UBOUND(grd,1)
       lambda = (1d0 / grd(k)) * cm_to_ang
       loglam = LOG10(lambda)

       ! 2p3 4So: Le Dourneuf et al. 1979
       crs0(1) = sig_interp(lambda, lam1, sig1)
       ! 2p3 2Do: Zeippen et al. 1980
       crs0(2) = sig_interp(lambda, lam2, sig2)
       ! 2P3 2Po: Zeippen et al. 1980
       crs0(3) = sig_interp(lambda, lam3, sig3)

       crs0 = crs0 * megabarn_to_cm2
       
       DO j = LBOUND(temp,1), UBOUND(temp,1)
          crs(k,j) = 0d0
          DO n = n_min, n_max
             IF(lambda > wa(n)) CYCLE
             crs(k,j) = crs(k,j) + crs0(n) * g(n) * EXP(-chi(n)/ (k_bol * temp(j)))
          END DO
       END DO
    END DO

    RETURN
  END SUBROUTINE cross_NI_mathisen

  
  SUBROUTINE cross_OI_mathisen(temp, grd, crs, n_min)
    USE const_module, ONLY : cm_to_ang, ev_to_erg, k_bol, megabarn_to_cm2
    IMPLICIT NONE

    REAL(REAL64), INTENT(IN) :: temp(:)
    REAL(REAL64), INTENT(IN) :: grd(:)
    REAL(REAL64), INTENT(OUT) :: crs(:,:)
    INTEGER(INT32), INTENT(IN) :: n_min

    INTEGER(INT32), PARAMETER :: n_max = 3
    
    REAL(REAL64), PARAMETER, DIMENSION(n_max) :: chi = [0.0d0, 1.967d0, 4.190d0] * ev_to_erg
    INTEGER(INT32), PARAMETER, DIMENSION(n_max) :: g = [9, 5, 1]

    INTEGER(INT32), PARAMETER :: n1 = 23
    REAL(REAL64), PARAMETER, DIMENSION(n1) :: &
         lam1 = [&
         910.d0, 800.d0, 733.d0, 732.9d0, 732.1d0, 732.d0, 700.d0, 666.d0, &
         665.9d0, 665.1d0, 665.d0, 650.d0, 625.d0, 600.d0, 510.d0, 480.d0, &
         465.d0, 455.d0, 450.d0, 400.d0, 300.d0, 200.d0, 125.d0], &
         sig1 = [&
         3.9d0,  4.1d0,  4.8d0,  5.1d0,  7.4d0,  7.7d0, 8.3d0, 8.9d0, &
         9.2d0, 11.3d0, 11.6d0, 12.6d0, 13.3d0, 13.1d0, 12.0d0, 12.7d0, &
         11.2d0, 12.2d0, 11.3d0, 10.3d0, 7.4d0,  3.8d0,  1.5d0]

    REAL(REAL64), DIMENSION(n_max) :: crs0
    REAL(REAL64) :: lambda, ratio, a, s, sigma_l, lambda0
    INTEGER(INT32) :: k, j, n

    DO k = LBOUND(grd,1), UBOUND(grd,1)
       lambda = (1d0 / grd(k)) * cm_to_ang

       ! 2p4 3P: Samson et al. 1983
       lambda0 = 910.4d0
       IF(lambda <= lambda0) THEN
          crs0(1) = sig_interp(lambda, lam1, sig1)
       END IF

       crs0(2) = 0d0
       ! 2p4 1D -> 2Do: Henry 1970
       lambda0 = 828d0
       IF(lambda <= lambda0) THEN
          a = 6.829d0
          s = 1.5d0
          sigma_l = 4.64d0
          ratio = (lambda / lambda0)
          crs0(2) = crs0(2) + sigma_l * (a * ratio**s + (1d0 - a) * ratio**(s + 1d0))
       END IF
       ! 2p4 1D -> 2Po: Henry 1970
       lambda0 = 744d0
       IF(lambda <= lambda0) THEN
          a = 4.800d0
          s = 1.5d0
          sigma_l = 1.95d0
          ratio = (lambda / lambda0)
          crs0(2) = crs0(2) + sigma_l * (a * ratio**s + (1d0 - a) * ratio**(s + 1d0))
       END IF

       ! 2p4 1S: Henry 1970
       lambda0 = 858d0
       IF(lambda <= lambda0) THEN
          a = 5.124d0
          s = 1.5d0
          sigma_l = 7.65d0
          ratio = (lambda / lambda0)
          crs0(3) = sigma_l * (a * ratio**s + (1d0 - a) * ratio**(s + 1d0))
       END IF
       
       crs0 = MAX(0d0,crs0 * megabarn_to_cm2)
       
       DO j = LBOUND(temp,1), UBOUND(temp,1)
          crs(k,j) = 0d0
          DO n = n_min, n_max
             crs(k,j) = crs(k,j) + crs0(n) * g(n) * EXP(-chi(n)/ (k_bol * temp(j)))
          END DO
       END DO
    END DO

    RETURN
  END SUBROUTINE cross_OI_mathisen

  
  SUBROUTINE cross_NaI_mathisen(temp, grd, crs, n_min)
    USE const_module, ONLY : cm_to_ang, ev_to_erg, k_bol, megabarn_to_cm2
    IMPLICIT NONE

    REAL(REAL64), INTENT(IN) :: temp(:)
    REAL(REAL64), INTENT(IN) :: grd(:)
    REAL(REAL64), INTENT(OUT) :: crs(:,:)
    INTEGER(INT32), INTENT(IN) :: n_min

    INTEGER(INT32), PARAMETER :: n_max = 2
    
    REAL(REAL64), PARAMETER, DIMENSION(n_max) :: &
         wa = [2412.66d0, 4086.62d0], &
         chi = [0.00d0, 2.104d0] * ev_to_erg
    INTEGER(INT32), PARAMETER, DIMENSION(n_max) :: &
         g = [2, 6]

    INTEGER(INT32), PARAMETER :: n1 = 34
    REAL(REAL64), PARAMETER, DIMENSION(n1) :: &
         lam1 = [&
         2412.7d0, 2400.d0, 2350.d0, 2300.d0, 2250.d0, 2200.d0, 2150.d0, 2100.d0, 2050.d0, &
         2000.d0, 1950.d0, 1900.d0, 1850.d0, 1800.d0, 1750.d0, 1700.d0, 1650.d0, 1600.d0, &
         1550.d0, 1500.d0, 1450.d0, 1400.d0, 1350.d0, 1300.d0, 1250.d0, 1200.d0, 1150.d0, &
         1100.d0, 1050.d0, 1000.d0,  900.d0, 775.d0,  660.d0,  575.d0], &
         sig1 = [&
         0.092d0, 0.089d0, 0.078d0, 0.065d0, 0.050d0, 0.032d0, 0.016d0, 0.006d0, 0.001d0, &
         0.000d0, 0.000d0, 0.000d0, 0.001d0, 0.004d0, 0.013d0, 0.027d0, 0.040d0, 0.051d0, &
         0.060d0, 0.069d0, 0.078d0, 0.087d0, 0.095d0, 0.103d0, 0.110d0, 0.118d0, 0.124d0, &
         0.129d0, 0.134d0, 0.136d0, 0.136d0, 0.132d0, 0.122d0, 0.110d0]
    
    REAL(REAL64), DIMENSION(n_max) :: crs0
    REAL(REAL64) :: lambda
    INTEGER(INT32) :: k, j, n

    DO k = LBOUND(grd,1), UBOUND(grd,1)
       lambda = (1d0 / grd(k)) * cm_to_ang

       ! 3s 2S: Butler & Mendoza 1983
       crs0(1) = sig_interp(lambda, lam1, sig1)
       ! 3p 2Po: Laughlin 1978, Rothe 1969
       crs0(2) = 7.95d0 * (lambda / wa(2))**3
       
       crs0 = crs0 * megabarn_to_cm2
       
       DO j = LBOUND(temp,1), UBOUND(temp,1)
          crs(k,j) = 0d0
          DO n = n_min, n_max
             IF(lambda >= wa(n)) CYCLE
             crs(k,j) = crs(k,j) + crs0(n) * g(n) * EXP(-chi(n)/ (k_bol * temp(j)))
          END DO
       END DO
    END DO

    RETURN
  END SUBROUTINE cross_NaI_mathisen

  
  SUBROUTINE cross_MgI_mathisen(temp, grd, crs, n_min)
    USE const_module, ONLY : cm_to_ang, ev_to_erg, k_bol, megabarn_to_cm2
    IMPLICIT NONE

    REAL(REAL64), INTENT(IN) :: temp(:)
    REAL(REAL64), INTENT(IN) :: grd(:)
    REAL(REAL64), INTENT(OUT) :: crs(:,:)
    INTEGER(INT32), INTENT(IN) :: n_min

    INTEGER(INT32), PARAMETER :: n_max = 8
    
    REAL(REAL64), PARAMETER, DIMENSION(n_max) :: &
         wa = [1621.5d0, 2513.8d0, 3756.6d0, 4884.0d0, 5504.d0, 6550.0d0, 7234.0d0, 7292.0d0], &
         chi = [0.0d0, 2.714d0, 4.346d0, 5.108d0, 5.394d0, 5.753d0, 5.932d0, 5.946d0] * ev_to_erg
    INTEGER(INT32), PARAMETER, DIMENSION(n_max) :: &
         g = [1, 9, 3, 3, 1, 5, 9, 15]

    INTEGER(INT32), PARAMETER :: n1 = 8, n2 = 7, n3 = 7
    REAL(REAL64), PARAMETER, DIMENSION(n1) :: &
         lam1 = [1622.d0, 1588.d0, 1519.d0, 1504.d0, 1419.d0, 1317.d0, 1262.d0, 1180.d0], &
         sig1 = [2.55d0, 2.00d0, 1.20d0, 1.00d0, 0.45d0, 0.00d0, 0.95d0, 0.01d0]
    REAL(REAL64), PARAMETER, DIMENSION(n2) :: &
         lam2 = [3757.d0, 2953.d0, 2676.d0, 2466.d0, 2304.d0, 2172.d0, 2059.d0], &
         sig2 = [125.d0, 100.d0, 80.d0, 60.d0, 40.d0, 20.d0, 0.1d0]
    REAL(REAL64), PARAMETER, DIMENSION(n3) :: &
         lam3 = [3100.d0, 3050.d0, 3025.d0, 3009.d0, 2990.d0, 2950.d0, 2900.d0], &
         sig3 = [0.1d0, 71.d0, 353.d0, 800.d0, 200.d0, 33.d0, 0.1d0]

    REAL(REAL64), DIMENSION(n_max) :: crs0
    REAL(REAL64) :: lambda
    INTEGER(INT32) :: k, j, n

    DO k = LBOUND(grd,1), UBOUND(grd,1)
       lambda = (1d0 / grd(k)) * cm_to_ang

       ! 3s2 1S: Bates & Altick 1973
       crs0(1) = sig_interp(lambda, lam1, sig1)
       ! 3p 3Po: Vernazza et al. 1976
       crs0(2) = 20d0 * (lambda / wa(2))**2.7d0
       ! 3p 1Po: Thomson et al. 1974, Bradley et al. 1976
       crs0(3) = sig_interp(lambda, lam2, sig2) + sig_interp(lambda, lam3, sig3)
       ! Vernazza et al. 1981
       ! 4s 3S
       crs0(4) = 2.2d0 * (lambda / wa(4))**2.3d0
       ! 4s 1S
       crs0(5) = 0.5d0 * (lambda / wa(5))**9.1d0
       ! 3d 1D
       crs0(6) =40.1d0 * (lambda / wa(6))**1.7d0
       ! 4p 3Po
       crs0(7) =26.3d0 * (lambda / wa(7))**3.0d0
       ! 3d 3D
       crs0(8) =34.2d0 * (lambda / wa(8))**2.0d0
       
       crs0 = crs0 * megabarn_to_cm2
       
       DO j = LBOUND(temp,1), UBOUND(temp,1)
          crs(k,j) = 0d0
          DO n = n_min, n_max
             IF(lambda >= wa(n)) CYCLE
             crs(k,j) = crs(k,j) + crs0(n) * g(n) * EXP(-chi(n)/ (k_bol * temp(j)))
          END DO
       END DO
    END DO

    RETURN
  END SUBROUTINE cross_MgI_mathisen

  
  SUBROUTINE cross_MgII_mathisen(temp, grd, crs, n_min)
    USE const_module, ONLY : cm_to_ang, ev_to_erg, k_bol, megabarn_to_cm2
    IMPLICIT NONE

    REAL(REAL64), INTENT(IN) :: temp(:)
    REAL(REAL64), INTENT(IN) :: grd(:)
    REAL(REAL64), INTENT(OUT) :: crs(:,:)
    INTEGER(INT32), INTENT(IN) :: n_min

    INTEGER(INT32), PARAMETER :: n_max = 13
    
    REAL(REAL64), PARAMETER, DIMENSION(n_max) :: &
         wa = [&
         824.6d0, 1169.d0, 1943.d0, 2009.d0, 2461.d0, 3512.d0, &
         3577.d0 , 3641.d0, 4201.d0, 5534.d0, 5603.d0, 5689.d0, 5694.d0], &
         chi = [&
         0.00d0, 4.430d0, 8.655d0, 8.864d0, 9.998d0, 11.505d0, &
         11.569d0, 11.603d0, 12.084d0, 12.795d0, 12.823d0, 12.856d0, 12.858d0] * ev_to_erg
    INTEGER(INT32), PARAMETER, DIMENSION(n_max) :: &
         g = [2, 6, 2, 10, 6, 2, 10, 14, 6, 2, 10, 14, 18]

    INTEGER(INT32), PARAMETER :: n1 = 5
    REAL(REAL64), PARAMETER, DIMENSION(n1) :: &
         lam1 = [825.d0, 568.d0, 433.d0, 350.d0, 293.d0], &
         sig1 = [0.210d0, 0.243d0, 0.211d0, 0.175d0, 0.150d0]

    REAL(REAL64), PARAMETER, DIMENSION(2:13) :: &
         a = [&
         0.474d0, 0.161d0, 5.335d0, 1.660d0, 0.102d0, 9.555d0, 4.800d0, &
         3.483d0, 0.0797d0, 13.809d0, 10.30d0, 4.400d0], &
         b = [&
         0.477d0, 4.305d0, 0.790d0, 0.510d0, 6.546d0, 2.060d0, 1.000d0, &
         0.636d0, 5.878d0, 1.733d0, 1.000d0, 1.000d0], &
         s = [&
         2.0d0, 1.5d0, 3.0d0, 3.0d0, 1.5d0, 3.5d0, 3.0d0, &
         3.0d0, 1.0d0, 3.0d0, 3.0d0, 3.0d0]
    
    REAL(REAL64), DIMENSION(n_max) :: crs0
    REAL(REAL64) :: lambda, ratio
    INTEGER(INT32) :: k, j, n

    DO k = LBOUND(grd,1), UBOUND(grd,1)
       lambda = (1d0 / grd(k)) * cm_to_ang

       ! 3s 2S: Butler et al. 1984
       crs0(1) = sig_interp(lambda, lam1, sig1)
       ! Mihalas 1972
       ! 02: 3p 2Po
       ! 03: 4s 2S
       ! 04: 3d 2D
       ! 05: 4p 2Po
       ! 06: 5s 2S
       ! 07: 4d 2D
       ! 08: 4f 2F
       ! 09: 5p 2Po
       ! 10: 6s 2S
       ! 11: 5d 2D
       ! 12: 5f 2Fo
       ! 13: 5g 2G
       DO n = 2, 13
          ratio = lambda / wa(n)
          crs0(n) = a(n) * (b(n) * ratio**s(n) + (1d0 - b(n)) * ratio**(s(n) + 1d0))
       END DO
       
       crs0 = MAX(0d0, crs0 * megabarn_to_cm2)
       
       DO j = LBOUND(temp,1), UBOUND(temp,1)
          crs(k,j) = 0d0
          DO n = n_min, n_max
             IF(lambda > wa(n)) CYCLE
             crs(k,j) = crs(k,j) + crs0(n) * g(n) * EXP(-chi(n)/ (k_bol * temp(j)))
          END DO
       END DO
    END DO

    RETURN
  END SUBROUTINE cross_MgII_mathisen

  
  SUBROUTINE cross_AlI_mathisen(temp, grd, crs, n_min)
    USE const_module, ONLY : cm_to_ang, ev_to_erg, k_bol, megabarn_to_cm2
    IMPLICIT NONE

    REAL(REAL64), INTENT(IN) :: temp(:)
    REAL(REAL64), INTENT(IN) :: grd(:)
    REAL(REAL64), INTENT(OUT) :: crs(:,:)
    INTEGER(INT32), INTENT(IN) :: n_min

    INTEGER(INT32), PARAMETER :: n_max = 7
    
    REAL(REAL64), PARAMETER, DIMENSION(n_max) :: &
         wa = [2071.3d0, 4361.0d0, 6312.0d0, 6528.0d0, 9444.0d0, 10519.d0, 12496.d0], &
         chi = [0.000d0, 3.143d0, 4.022d0, 4.087d0, 4.673d0, 4.807d0, 4.994d0] * ev_to_erg
    INTEGER(INT32), PARAMETER, DIMENSION(n_max) :: &
         g = [6, 2, 10, 6, 2, 10, 6]

    REAL(REAL64), DIMENSION(n_max) :: crs0
    REAL(REAL64) :: lambda
    INTEGER(INT32) :: k, j, n

    DO k = LBOUND(grd,1), UBOUND(grd,1)
       lambda = (1d0 / grd(k)) * cm_to_ang

       ! 3p 2Po: Vernazza et al. 1976
       crs0(1) = 65d0 * (lambda / wa(1))**4.4d0
       ! 4s 2S: Dragon & Mutschlecner 1980
       crs0(2) = 1.22d0 * (lambda / wa(2))**3d0
       ! Vernezza et al. 1976
       ! 3d 2D
       crs0(3) = 47.0d0 * (lambda / wa(3))**1.83d0
       ! 4p 2Po
       crs0(4) = 14.5d0 * (lambda / wa(4))**1.00d0
       ! 5s 2S
       crs0(5) = 56.7d0 * (lambda / wa(5))**1.90d0
       ! 4d 2D
       crs0(6) = 50.0d0 * (lambda / wa(6))**3.00d0
       ! 5p 2Po
       crs0(7) = 50.0d0 * (lambda / wa(7))**3.00d0
       
       crs0 = crs0 * megabarn_to_cm2
       
       DO j = LBOUND(temp,1), UBOUND(temp,1)
          crs(k,j) = 0d0
          DO n = n_min, n_max
             IF(lambda > wa(n)) CYCLE
             crs(k,j) = crs(k,j) + crs0(n) * g(n) * EXP(-chi(n)/ (k_bol * temp(j)))
          END DO
       END DO
    END DO

    RETURN
  END SUBROUTINE cross_AlI_mathisen

  
  SUBROUTINE cross_SiI_mathisen(temp, grd, crs, n_min)
    USE const_module, ONLY : cm_to_ang, ev_to_erg, k_bol, megabarn_to_cm2
    IMPLICIT NONE

    REAL(REAL64), INTENT(IN) :: temp(:)
    REAL(REAL64), INTENT(IN) :: grd(:)
    REAL(REAL64), INTENT(OUT) :: crs(:,:)
    INTEGER(INT32), INTENT(IN) :: n_min

    INTEGER(INT32), PARAMETER :: n_max = 7
    
    REAL(REAL64), PARAMETER, DIMENSION(n_max) :: &
         wa = [1521.d0, 1682.d0, 1986.d0, 3863.d0, 4040.d0, 5416.d0, 5839.d0], &
         chi = [0.000d0,0.781d0, 1.909d0, 4.942d0, 5.082d0, 5.863d0, 6.028d0] * ev_to_erg
    INTEGER(INT32), PARAMETER, DIMENSION(n_max) :: &
         g = [9, 5, 1, 9, 3, 3, 27]

    REAL(REAL64), PARAMETER, DIMENSION(1:3) :: &
         sig_l = [39.16d0, 34.49d0, 33.56d0], &
         a = [4.420d0, 6.459d0, 10.013d0], &
         b = [8.934d0, 5.142d0, 5.500d0], &
         s = [5d0, 3d0, 3d0]
    
    REAL(REAL64), DIMENSION(n_max) :: crs0
    REAL(REAL64) :: lambda, ratio
    INTEGER(INT32) :: k, j, n

    DO k = LBOUND(grd,1), UBOUND(grd,1)
       lambda = (1d0 / grd(k)) * cm_to_ang
       
       DO n = 1, 3
          ratio = lambda / wa(n)
          crs0(n) = sig_l(n) * (a(n) * ratio**s(n) + (b(n) - 2d0 * a(n)) * ratio**(s(n) + 1d0) + &
               (1d0 + a(n) - b(n)) * ratio**(s(n) + 2d0))
       END DO
       crs0(4) = 1.25d0 * (lambda / wa(4))**2d0
       crs0(5) = 4.09d0 * (lambda / wa(5))**2d0
       crs0(6) =20.40d0 * (lambda / wa(6))**3d0
       crs0(7) =14.10d0 * (lambda / wa(7))**3d0
       
       crs0 = MAX(0d0, crs0 * megabarn_to_cm2)
       
       DO j = LBOUND(temp,1), UBOUND(temp,1)
          crs(k,j) = 0d0
          DO n = n_min, n_max
             IF(lambda > wa(n)) CYCLE
             crs(k,j) = crs(k,j) + crs0(n) * g(n) * EXP(-chi(n)/ (k_bol * temp(j)))
          END DO
       END DO
    END DO

    RETURN
  END SUBROUTINE cross_SiI_mathisen

  
  SUBROUTINE cross_SiII_mathisen(temp, grd, crs, n_min)
    USE const_module, ONLY : cm_to_ang, ev_to_erg, k_bol, megabarn_to_cm2
    IMPLICIT NONE

    REAL(REAL64), INTENT(IN) :: temp(:)
    REAL(REAL64), INTENT(IN) :: grd(:)
    REAL(REAL64), INTENT(OUT) :: crs(:,:)
    INTEGER(INT32), INTENT(IN) :: n_min

    INTEGER(INT32), PARAMETER :: n_max = 1
    
    REAL(REAL64), PARAMETER, DIMENSION(n_max) :: &
         wa = [758.5d0], &
         chi = [0.00d0] * ev_to_erg
    INTEGER(INT32), PARAMETER, DIMENSION(n_max) :: &
         g = [6]

    INTEGER(INT32), PARAMETER :: n1 = 15
    REAL(REAL64), PARAMETER, DIMENSION(n1) :: &
         lam1 = [&
         759.d0, 747.d0, 742.d0, 741.d0, 735.d0, 729.d0, 723.d0, 712.d0, &
         690.d0, 633.d0, 616.d0, 610.d0, 600.d0, 556.d0, 542.d0], &
         sig1 = [&
         2.8d0, 1.7d0, 0.8d0, 2.0d0, 10.5d0, 6.1d0, 3.9d0, 3.0d0, &
         2.2d0, 0.9d0, 0.5d0, 1.6d0, 1.1d0, 0.5d0, 0.4d0]
    
    REAL(REAL64), DIMENSION(n_max) :: crs0
    REAL(REAL64) :: lambda
    INTEGER(INT32) :: k, j, n

    DO k = LBOUND(grd,1), UBOUND(grd,1)
       lambda = (1d0 / grd(k)) * cm_to_ang

       crs0(1) = sig_interp(lambda, lam1, sig1)
       
       crs0 = crs0 * megabarn_to_cm2
       
       DO j = LBOUND(temp,1), UBOUND(temp,1)
          crs(k,j) = 0d0
          DO n = n_min, n_max
             IF(lambda > wa(n)) CYCLE
             crs(k,j) = crs(k,j) + crs0(n) * g(n) * EXP(-chi(n)/ (k_bol * temp(j)))
          END DO
       END DO
    END DO

    RETURN
  END SUBROUTINE cross_SiII_mathisen

  
  SUBROUTINE cross_SI_mathisen(temp, grd, crs, n_min)
    USE const_module, ONLY : cm_to_ang, ev_to_erg, k_bol, megabarn_to_cm2
    IMPLICIT NONE

    REAL(REAL64), INTENT(IN) :: temp(:)
    REAL(REAL64), INTENT(IN) :: grd(:)
    REAL(REAL64), INTENT(OUT) :: crs(:,:)
    INTEGER(INT32), INTENT(IN) :: n_min

    INTEGER(INT32), PARAMETER :: n_max = 3
  
    REAL(REAL64), PARAMETER, DIMENSION(n_max) :: chi = [0.000d0, 1.145d0, 2.750d0] * ev_to_erg
    INTEGER(INT32), PARAMETER, DIMENSION(n_max) :: g = [9, 5, 1]

    INTEGER(INT32), PARAMETER :: n1 = 23, n2 = 10, n3 = 49
    REAL(REAL64), PARAMETER, DIMENSION(n1) :: &
         lam1 = [&
         1195d0, 1190.d0, 1177.d0, 1176.5d0, 1176.d0, 1174.5d0, 1174.d0, 1173.5d0, &
         1173.d0, 1171.d0, 1170.d0, 1169.d0, 1168.5d0, 1168.d0, 1167.d0, 1165.d0, &
         1164.d0, 1160.d0, 1150.d0 , 1140.d0 , 1130.d0, 1124.d0, 1121.d0], &
         sig1 = [&
         140.d0, 110.d0, 61.d0, 80.d0, 130.d0, 65.d0, 55.d0, 65.d0, &
         140.d0, 106.d0, 180.d0, 60.d0, 50.d0, 60.d0, 175.d0, 90.d0, &
         70.d0, 54.d0, 38.d0, 30.d0, 18.d0, 8.d0, 15.d0]
    REAL(REAL64), PARAMETER, DIMENSION(n2) :: &
         lam2a = [1197.d0, 1123.d0, 1058.d0, 948.d0, 785.d0, 669.d0, 584.d0,  517.d0,  453.d0, 403.d0], &
         sig2a = [13.82d0, 14.61d0, 15.03d0, 15.04d0, 12.85d0, 9.47d0, 6.27d0, 3.84d0, 1.96d0, 0.97d0], &
         lam2b = [1016.d0, 962.d0, 914.d0, 831.d0, 703.d0, 609.d0, 537.d0, 480.d0, 424.d0, 380.d0], &
         sig2b = [27.43d0, 28.00d0, 27.46d0, 24.21d0, 15.42d0, 8.70d0, 4.67d0, 2.47d0, 1.14d0, 0.56d0], &
         lam2c = [925.d0, 880.d0, 840.d0, 769.d0, 658.d0, 575.d0, 510.d0, 459.d0, 408.d0, 367.d0], &
         sig2c = [20.87d0, 20.78d0, 19.73d0, 16.03d0, 8.43d0, 4.01d0, 1.92d0, 0.95d0, 0.44d0, 0.23d0]
    REAL(REAL64), PARAMETER, DIMENSION(n3) :: &
         lam3 = [&
         1629.d0, 1613.d0, 1600.d0, 1587.d0, 1578.d0, 1575.d0, 1570.d0, 1562.d0, 1538.d0, 1504.d0, &
         1389.d0, 1379.d0, 1375.d0, 1370.d0, 1364.d0, 1360.d0, 1357.d0, 1355.5d0, 1354.d0, 1351.d0, &
         1349.d0, 1316.d0, 1313.d0, 1312.1d0, 1312.d0, 1300.d0, 1292.d0, 1282.d0, 1277.d0, 1270.d0, &
         1266.d0, 1264.d0, 1262.d0, 1256.d0, 1250.d0, 1235.d0, 1227.d0, 1224.d0, 1221.d0, 1219.d0, &
         1212.d0, 1208.d0, 1205.d0, 1202.d0, 1198.d0, 1191.d0, 1188.d0, 1176.d0, 1163.d0], &
         sig3 = [&
         100.d0, 18.d0,   5.d0, 0.3d0, 0.2d0,  6.d0, 2.d0, 1.d+0, 0.3d0, 0.d0, &
         0.d0,  2.d0, 100.d0, 1.d0 , 0.d0 , 30.d0, 5.d0, 2.d0, 5.d0, 100.d0, &
         1.d0,  0.d0,   0.d0, 2.9d0, 3.d0,   4.d0, 10.d0, 70.d0, 800.d0, 100.d0, &
         2.d0, 40.d0, 500.d0, 100.d0, 80.d0, 250.d0, 50.d0, 500.d0, 14.d0, 12.d0, &
         40.d0, 800.d0, 20.d0, 13.d0, 25.d0, 800.d0, 40.d0, 6.d0, 5.d0]
    
    REAL(REAL64), PARAMETER, DIMENSION(1:3) :: &
         sig_l = [20.99d0, 7.87d0, 22.56d0], &
         a = [0.659d0, 0.543d0, -1.148d0], &
         b = [6.193d0, 7.433d0, 6.973d0], &
         s = [1.5d0, 2d0, 1.0d0]
    
    REAL(REAL64), DIMENSION(n_max) :: crs0
    REAL(REAL64) :: lambda, ratio, lambda0
    INTEGER(INT32) :: k, j, n, m

    ! C1a: 3p4 3P (g=9, chi=    0) -> ?     (wa=1196.8) : Tondello 72
    ! C1b: 3p4 3P (g=9, chi=    0) -> 4S    (wa=1197  ) : Dill+ 75
    ! C1c: 3p4 3P (g=9, chi=    0) -> 2D    (wa=1016  ) : Dill+ 75
    ! C1d: 3p4 3P (g=9, chi=    0) -> 2P    (wa= 925  ) : Dill+ 75
    ! C2a: 3p4 1D (g=5, chi=1.145) -> 2D    (wa=1121  ) : Chapman & Henry 71
    ! C2b: 3p4 1D (g=5, chi=1.145) -> 2P    (wa=1011  ) : Chapman & Henry 71
    ! C3a: 3p4 1S (g=1, chi=2.750) -> 4S,2D (wa=1629  ) : McGuire 79
    ! C3b: 3p4 1S (g=1, chi=2.750) -> 2P    (wa=1164  ) : Chapman & Henry 71

    DO k = LBOUND(grd,1), UBOUND(grd,1)
       lambda = (1d0 / grd(k)) * cm_to_ang
       crs0 = 0d0

       IF(lambda <= 1196.8d0) THEN ! C1a
          crs0(1) = crs0(1) + sig_interp(lambda, lam1, sig1)
       END IF
       IF(lambda <= 1197d0) THEN ! C1b
          crs0(1) = crs0(1) + sig_interp(lambda, lam2a, sig2a)
       ENDIF
       IF(lambda <= 1016d0) THEN ! C1c
          crs0(1) = crs0(1) + sig_interp(lambda, lam2b, sig2b)
       ENDIF
       IF(lambda <=  925d0) THEN ! C1d
          crs0(1) = crs0(1) + sig_interp(lambda, lam2c, sig2c)
       ENDIF

       lambda0 = 1121d0
       IF(lambda <= lambda0) THEN ! C2a
          m = 1
          ratio = lambda / lambda0
          crs0(2) = crs0(2) + sig_l(m) * (&
                            (a(m)) * ratio**(s(m)) + &
               (b(m) - 2d0 * a(m)) * ratio**(s(m) + 1d0) + &
               (1d0 + a(m) - b(m)) * ratio**(s(m) + 2d0))
       ENDIF
       lambda0 = 1011d0
       IF(lambda <= lambda0) THEN ! C2b
          m = 2
          ratio = lambda / lambda0
          crs0(2) = crs0(2) + sig_l(m) * (&
                            (a(m)) * ratio**(s(m)) + &
               (b(m) - 2d0 * a(m)) * ratio**(s(m) + 1d0) + &
               (1d0 + a(m) - b(m)) * ratio**(s(m) + 2d0))
       ENDIF

       IF(lambda <= 1629d0) THEN ! C3a
          crs0(3) = crs0(3) + sig_interp(lambda, lam3, sig3)
       ENDIF
       lambda0 = 1164d0
       IF(lambda <= lambda0) THEN ! C3b
          m = 3
          ratio = lambda / lambda0
          crs0(3) = sig_l(m) * (&
                            (a(m)) * ratio**(s(m)) + &
               (b(m) - 2d0 * a(m)) * ratio**(s(m) + 1d0) + &
               (1d0 + a(m) - b(m)) * ratio**(s(m) + 2d0))
       END IF

       crs0 = MAX(0d0,crs0) * megabarn_to_cm2
       
       DO j = LBOUND(temp,1), UBOUND(temp,1)
          crs(k,j) = 0d0
          DO n = n_min, n_max
             crs(k,j) = crs(k,j) + crs0(n) * g(n) * EXP(-chi(n)/ (k_bol * temp(j)))
          END DO
       END DO
    END DO

    RETURN
  END SUBROUTINE cross_SI_mathisen

  
  SUBROUTINE cross_CaI_mathisen(temp, grd, crs, n_min)
    USE const_module, ONLY : cm_to_ang, ev_to_erg, k_bol, megabarn_to_cm2
    IMPLICIT NONE

    REAL(REAL64), INTENT(IN) :: temp(:)
    REAL(REAL64), INTENT(IN) :: grd(:)
    REAL(REAL64), INTENT(OUT) :: crs(:,:)
    INTEGER(INT32), INTENT(IN) :: n_min

    INTEGER(INT32), PARAMETER :: n_max = 3
  
    REAL(REAL64), PARAMETER, DIMENSION(n_max) :: chi = [0.000d0, 1.892d0, 2.570d0] * ev_to_erg, &
         wa = [2028.2d0, 2938d0, 3500d0]
    INTEGER(INT32), PARAMETER, DIMENSION(n_max) :: g = [1, 9, 20]

    INTEGER(INT32), PARAMETER :: n1 = 30
    REAL(REAL64), PARAMETER, DIMENSION(n1) :: &
         lam1 = [&
         2029.d0, 2012.d0, 2001.d0, 1974.d0, 1943.d0, 1922.d0, 1912.d0, 1892.d0, 1873.d0, 1868.d0, &
         1863.d0, 1844.d0, 1835.d0, 1817.d0, 1799.d0, 1773.d0, 1764.d0, 1754.d0, 1752.1d0, 1751.3d0, &
         1750.4d0, 1749.d0, 1739.d0, 1722.d0, 1709.d0, 1690.d0, 1587.d0, 1540.d0, 1415.d0, 1342.d0], &
         sig1 = [&
         3.50d0, 2.05d0, 1.48d0, 0.50d0, 0.00d0, 0.29d0, 0.88d0, 5.45d0, 33.70d0, 38.40d0, &
         32.28d0, 8.04d0, 4.57d0, 1.76d0, 0.69d0, 0.02d0, 0.07d0, 10.69d0, 52.64d0, 130.21d0, &
         100.71d0, 19.95d0, 1.18d0, 0.14d0, 0.61d0, 1.06d0, 2.08d0, 2.22d0, 1.30d0, 0.90d0]
    
    REAL(REAL64), DIMENSION(n_max) :: crs0
    REAL(REAL64) :: lambda
    INTEGER(INT32) :: k, j, n

    ! C1: 4s2 1S  (g= 1, chi=    0) -> ?     (wa=2028.2) : Scott+ 83
    ! C2: 4p 3Po  (g= 9, chi=1.892) -> ?     (wa=2938  ) : Kelm & Schluter 1962
    ! C3: 3d 1D,3D(g=20, chi=2.570) -> ?     (wa=3500  ) : Kelm & Schluter 1962

    DO k = LBOUND(grd,1), UBOUND(grd,1)
       lambda = (1d0 / grd(k)) * cm_to_ang
       
       crs0(1) = sig_interp(lambda, lam1, sig1)
       crs0(2) = 6.1d0 * (lambda / wa(2))**20d0
       crs0(3) = 5.0d0 * (lambda / wa(3))**15d0

       crs0 = crs0 * megabarn_to_cm2
       
       DO j = LBOUND(temp,1), UBOUND(temp,1)
          crs(k,j) = 0d0
          DO n = n_min, n_max
             IF(lambda > wa(n)) CYCLE
             crs(k,j) = crs(k,j) + crs0(n) * g(n) * EXP(-chi(n)/ (k_bol * temp(j)))
          END DO
       END DO
    END DO

    RETURN
  END SUBROUTINE cross_CaI_mathisen

  
  SUBROUTINE cross_CaII_mathisen(temp, grd, crs, n_min)
    USE const_module, ONLY : cm_to_ang, ev_to_erg, k_bol, megabarn_to_cm2
    IMPLICIT NONE

    REAL(REAL64), INTENT(IN) :: temp(:)
    REAL(REAL64), INTENT(IN) :: grd(:)
    REAL(REAL64), INTENT(OUT) :: crs(:,:)
    INTEGER(INT32), INTENT(IN) :: n_min

    INTEGER(INT32), PARAMETER :: n_max = 12
  
    REAL(REAL64), PARAMETER, DIMENSION(n_max) :: &
         wa = [&
         1044.4d0, 1219.d0, 1420.d0, 2294.d0, 2571.d0, 2844.d0, &
         3611.0d0, 3988.d0, 4343.d0, 4708.d0, 5644.d0, 5686.d0], &
         chi = [&
         0.000d0, 1.697d0, 3.142d0, 6.468d0, 7.049d0, 7.512d0, &
         8.438d0, 8.763d0, 9.017d0, 9.238d0, 9.675d0, 9.691d0] * ev_to_erg
    INTEGER(INT32), PARAMETER, DIMENSION(n_max) :: &
         g = [2, 10, 6, 2, 10, 6, 14, 2, 10, 6, 14, 18]

    INTEGER(INT32), PARAMETER :: n1 = 11
    REAL(REAL64), PARAMETER, DIMENSION(n1) :: &
         lam1 = [1045.d0, 1033.d0, 999.d0, 947.d0, 883.d0, 812.d0, 740.d0,  669.d0, 603.d0, 542.d0, 487.d0], &
         sig1 = [.2094d0, .2119d0, .2179d0, .2252d0, .2299d0, .2295d0, .2224d0, .2106d0, .1944d0, .1753d0, .1554d0]
    REAL(REAL64), PARAMETER, DIMENSION(2:12) :: &
         a = [6.14d0, 2.40d0, 0.08d0, 21.50d0, 6.40d0, 7.07d0, 0.064d0, 23.78d0, 9.63d0, 14.91d0, 4.39d0], &
         b = [1.919d0, 0.302d0, 5.052d0, 2.904d0, 0.099d0, 2.163d0, 5.129d0, 3.291d0, 1.446d0, 2.427d0, 1.275d0], &
         s = [1.5d0, 3.0d0, 2.0d0, 3.75d0, 2.0d0, 4.5d0, 1.25d0, 3.5d0, 3.0d0, 4.0d0, 4.5d0]
    
    REAL(REAL64), DIMENSION(n_max) :: crs0
    REAL(REAL64) :: lambda, ratio
    INTEGER(INT32) :: k, j, n

    !  C1: 4s 2S  (g= 2, chi=    0) -> ?     (wa=1044.4) : Black+ 72
    !  C2: 3d 2D  (g=10, chi=1.697) -> ?     (wa=1219  ) : Mihalas 73
    !  C3: 4p 2Po (g= 6, chi=3.142) -> ?     (wa=1420  ) : Mihalas 73
    !  C4: 5s 2S  (g= 2, chi=6.468) -> ?     (wa=2294  ) : Mihalas 73
    !  C5: 4d 2D  (g=10, chi=7.049) -> ?     (wa=2571  ) : Mihalas 73
    !  C6: 5p 2Po (g= 6, chi=7.512) -> ?     (wa=2844  ) : Mihalas 73
    !  C7: 4f 2Fo (g=14, chi=8.438) -> ?     (wa=3611  ) : Mihalas 73
    !  C8: 6s 2S  (g= 2, chi=8.763) -> ?     (wa=3988  ) : Mihalas 73
    !  C9: 5d 2D  (g=10, chi=9.017) -> ?     (wa=4343  ) : Mihalas 73
    ! C10: 6p 2Po (g= 6, chi=9.238) -> ?     (wa=4708  ) : Mihalas 73
    ! C11: 5f 2Fo (g=14, chi=9.675) -> ?     (wa=5644  ) : Mihalas 73
    ! C12: 5g 2G  (g=18, chi=9.691) -> ?     (wa=5686  ) : Mihalas 73
    DO k = LBOUND(grd,1), UBOUND(grd,1)
       lambda = (1d0 / grd(k)) * cm_to_ang

       crs0(1) = sig_interp(lambda, lam1, sig1)
       DO n = 2, n_max
          ratio = lambda / wa(n)
          crs0(n) = a(n) * (b(n) * ratio**s(n) + (1d0 - b(n)) * ratio**(s(n) + 1d0))
       END DO

       crs0 = crs0 * megabarn_to_cm2
       
       DO j = LBOUND(temp,1), UBOUND(temp,1)
          crs(k,j) = 0d0
          DO n = n_min, n_max
             IF(lambda > wa(n)) CYCLE
             crs(k,j) = crs(k,j) + crs0(n) * g(n) * EXP(-chi(n)/ (k_bol * temp(j)))
          END DO
       END DO
    END DO

    RETURN
  END SUBROUTINE cross_CaII_mathisen

  
  SUBROUTINE cross_FeI_mathisen(temp, grd, crs, n_min)
    USE const_module, ONLY : cm_to_ang, ev_to_erg, k_bol, megabarn_to_cm2
    IMPLICIT NONE

    REAL(REAL64), INTENT(IN) :: temp(:)
    REAL(REAL64), INTENT(IN) :: grd(:)
    REAL(REAL64), INTENT(OUT) :: crs(:,:)
    INTEGER(INT32), INTENT(IN) :: n_min

    INTEGER(INT32), PARAMETER :: n_max = 26
  
    REAL(REAL64), PARAMETER, DIMENSION(n_max) :: &
         wa = [&
         1569.d0, 1785.d0, 1958.d0, 2184.d0, 2246.d0, 2278.d0, 2280.d0, 2344.d0, 2408.d0, &
         2464.d0, 2466.d0, 2536.d0, 2570.d0, 2679.d0, 2693.d0, 2735.d0, 2864.d0, 2936.d0, &
         3139.d0, 3310.d0, 3376.d0, 3430.d0, 3555.d0, 3779.d0, 3780.d0, 4004.d0], &
         chi = [&
         0.000d0, 0.925d0, 1.538d0, 2.193d0, 2.350d0, 2.427d0, 2.433d0, 2.581d0, 2.721d0, &
         2.839d0, 2.843d0, 2.982d0, 3.046d0, 3.242d0, 3.267d0, 3.375d0, 3.541d0, 3.647d0, &
         3.921d0, 4.124d0, 4.198d0, 4.256d0, 4.383d0, 4.590d0, 4.591d0, 4.774d0] * ev_to_erg
    INTEGER(INT32), PARAMETER, DIMENSION(n_max) :: &
         g = [&
         25, 35, 21, 15,  9, 33, 35, 21, 27,  9, 49, 48, 18, &
         25, 48, 35, 19, 52, 36, 21, 60, 14, 72, 42, 15, 15], &
         indx = [&
         1, 14, 15, 16, 2, 3, 5, 4, 17, 18, 6, 11, 19, 7, 20, 8, 21, 12, 9, 26, 22, 13, 23, 25, 10, 24]

    INTEGER(INT32), PARAMETER :: n1 = 20
    REAL(REAL64), PARAMETER, DIMENSION(n1) :: &
         lam1 = [&
         1569.d0, 1560.d0, 1552.d0, 1545.d0, 1540.d0, 1535.d0, 1530.d0, 1510.d0, 1430.d0, 1400.d0, &
         1373.d0, 1355.d0, 1350.d0, 1334.d0, 1325.d0, 1305.d0, 1280.d0, 1250.d0, 1200.d0, 1150.d0], &
         sig1 = [&
         1.0d0, 2.9d0, 7.7d0, 4.4d0, 5.0d0, 5.3d0, 3.8d0, 0.9d0, 1.2d0, 2.5d0, &
         11.3d0, 5.4d0, 9.2d0, 6.3d0, 10.8d0, 2.8d0, 6.7d0, 2.3d0,  0.6d0, 0.4d0]

    REAL(REAL64), PARAMETER, DIMENSION(2:26) :: &
         a = [&
         1.007d0, 1.007d0, 0.402d0, 1.007d0, 1.006d0, 1.006d0, 0.442d0, 0.462d0, 0.483d0, &
         0.424d0, 1.005d0, 1.004d0, 1.140d0, 1.149d0, 0.498d0, 1.177d0, 0.561d0, 1.188d0, &
         0.597d0, 1.207d0, 0.673d0, 1.255d0, 1.285d0, 0.709d0, 0.667d0], &
         cs = [&
         3.331d0, 3.331d0, 4.117d0, 3.331d0, 3.609d0, 3.918d0, 4.840d0, 5.505d0, 6.618d0, &
         4.453d0, 4.292d0, 5.009d0, 2.007d0, 2.192d0, 2.914d0, 2.665d0, 3.268d0, 2.833d0, &
         3.552d0, 3.134d0, 4.379d0, 3.825d0, 4.263d0, 4.855d0, 4.299d0]
    
    REAL(REAL64), DIMENSION(n_max) :: crs0 = 0d0
    REAL(REAL64) :: lambda, ratio
    INTEGER(INT32) :: k, j, n

    DO k = LBOUND(grd,1), UBOUND(grd,1)
       lambda = (1d0 / grd(k)) * cm_to_ang

       crs0(1) = sig_interp(lambda, lam1, sig1)
       DO n = 2, n_max
          ratio = lambda / wa(n)
          crs0(n) = cs(indx(n)) * (a(indx(n)) * ratio**3 + (1d0 - a(indx(n))) * ratio**4)
       END DO

       crs0 = crs0 * megabarn_to_cm2
       
       DO j = LBOUND(temp,1), UBOUND(temp,1)
          crs(k,j) = 0d0
          DO n = n_min, n_max
             IF(lambda > wa(n)) CYCLE
             crs(k,j) = crs(k,j) + crs0(n) * g(n) * EXP(-chi(n)/ (k_bol * temp(j)))
          END DO
       END DO
    END DO

    RETURN
  END SUBROUTINE cross_FeI_mathisen

  
  REAL(REAL64) FUNCTION sig_interp(lambda, lam, sig)
    IMPLICIT NONE
    REAL(REAL64), INTENT(IN) :: lambda
    REAL(REAL64), INTENT(IN) :: lam(:)
    REAL(REAL64), INTENT(IN) :: sig(:)
    INTEGER(INT32) :: nlam, n
    
    nlam = UBOUND(lam,1)

    IF(lam(1) < lam(nlam)) THEN
       PRINT *, 'invalid order in sig_interp'
       STOP
    END IF
    IF(lambda > lam(1) .OR. lambda <= lam(nlam)) THEN
       sig_interp = 0d0
    ELSE
       DO n = 2, nlam
          IF(lambda > lam(n)) EXIT
       END DO
       sig_interp = ((lambda - lam(n)) * sig(n-1) + (lam(n-1) - lambda) * sig(n)) / (lam(n-1) - lam(n))
    END IF
    
    RETURN
  END FUNCTION sig_interp

  
  SUBROUTINE cross_rayleigh_he_tarafdar_1969(grd, crs)
    ! Rayleigh scattering cross sections for He
    ! Tarafdar & Vardya (1969), MNRAS, 145, 171
    IMPLICIT NONE
    REAL(REAL64), INTENT(IN) :: grd(:)    ! wavenumber grid [cm^-1]
    REAL(REAL64), INTENT(OUT) :: crs(:) ! cross section [cm^-2]

    REAL(REAL64), PARAMETER, DIMENSION(3) :: a = [5.616d-46, 0.475d-10, 0.185d-20] ! for He
!   REAL(REAL64), PARAMETER, DIMENSION(3) :: a = [4.426d-44, 2.743d-10, 6.960d-20] ! for C
!   REAL(REAL64), PARAMETER, DIMENSION(3) :: a = [2.275d-44, 1.992d-10, 3.362d-20] ! for N
!   REAL(REAL64), PARAMETER, DIMENSION(3) :: a = [0.762d-44, 1.008d-10, 1.038d-20] ! for O
    REAL(REAL64), PARAMETER :: lambda_max_he = 584d0 ! 1st resonance line of He I
    REAL(REAL64) :: lambda
    INTEGER(INT32) :: k, ks, ke

    ks = 1
    ke = UBOUND(grd,1)

    crs = 0d0
    DO k = ks, ke
       lambda = 1d0 / grd(k)
       IF(lambda < lambda_max_he) CYCLE

       crs(k) = a(1) / lambda**4 * (1d0 + a(2) / lambda**2 + a(3) / lambda**4)
    END DO

    RETURN    
  END SUBROUTINE cross_rayleigh_he_tarafdar_1969
  
    
  SUBROUTINE cross_rayleigh_he_rohrmann_2018(grd, crs)
    ! Rayleigh scattering cross sections for He
    ! Rohrmann (2018), MNRAS, 473, 457
    USE const_module, ONLY : clight, pi, hartree, hbar, a_0
    IMPLICIT NONE
    REAL(REAL64), INTENT(IN) :: grd(:)    ! wavenumber grid [cm^-1]
    REAL(REAL64), INTENT(OUT) :: crs(:) ! cross section [cm^-2]

    REAL(REAL64), PARAMETER, DIMENSION(14) :: &
         omegaj = [0.779748d0, 0.848433d0, 0.872505d0, 0.883667d0, 0.889738d0, 0.893403d0, 0.895782d0, &
         0.897415d0, 0.898583d0, 0.899448d0, 0.900105d0, 0.900617d0, 0.901024d0, 0.901352d0], &
         omega0 = [0.821611d0, 0.862814d0, 0.879035d0, 0.887165d0, 0.891823d0, 0.894740d0, 0.896694d0, &
         0.898068d0, 0.899066d0, 0.899812d0, 0.900386d0, 0.900845d0, 0.901212d0, 0.901522d0], &
         alphad = [186.9873d0, 465.6428d0, 1003.9436d0, 1728.4265d0, 2702.5630d0, 4452.0236d0, 6436.1707d0, &
         9005.2470d0, 12444.8560d0, 17253.3920d0, 21009.2542d0, 30141.8518d0, 39888.7213d0, 71608.2446d0]
    REAL(REAL64), PARAMETER, DIMENSION(0:3) :: pj = [1.383817d0, -2.7631796d0, 1.2159289d0, -0.0188510d0]
    REAL(REAL64), PARAMETER, DIMENSION(1:3) :: qj = [-3.1117233d0, 2.8718963d0, -0.7404416d0]
    REAL(REAL64), PARAMETER :: omegau = hartree/ hbar
    REAL(REAL64) :: omega, b, z, alpha = 0d0
    INTEGER(INT32) :: k, j

    crs = 0d0
    DO k = LBOUND(grd,1), UBOUND(grd,1)
       omega = 2d0 * pi * clight * grd(k) / omegau
       alpha = 0d0
       IF(omega < 0.7339d0) THEN
          alpha = (pj(0) + pj(1)*omega**2 + pj(2)*omega**4 + pj(3)*omega**6) / &
                    (1d0 + qj(1)*omega**2 + qj(2)*omega**4 + qj(3)*omega**6)
       ELSE IF(omega < omegaj(1)) THEN
          z = omega - 0.7282d0
          alpha = EXP(1.6520099d0 + 22.475155d0 * z + 1.0197470d12 * z**9)
       ELSE
          DO j = 13, 1, -1
             IF(omega < omegaj(j+1)) THEN
                IF(omega < omega0(j)) THEN
                   b = pi / (2d0 * (omega0(j) - omegaj(j)))
                ELSE
                   b = pi / (2d0 * (omegaj(j+1) - omega0(j)))
                END IF
                alpha = alphad(j) / b * TAN(b * (omega - omega0(j)))
             END IF
          END DO
       END IF
       alpha = alpha * a_0**3
       crs(k) = 8d0 * pi / (3d0 * clight**4) * (omega * omegau)**4 * alpha**2
    END DO

    RETURN    
  END SUBROUTINE cross_rayleigh_he_rohrmann_2018
  
    
  SUBROUTINE cross_rayleigh_h2_tarafdar_1973(grd, crs)
    ! Rayleigh scattering cross sections for H2
    ! Tarafdar & Vardya (1973), MNRAS, 163, 261
    USE const_module, ONLY : clight, h, ev_to_erg
    IMPLICIT NONE
    REAL(REAL64), INTENT(IN) :: grd(:)    ! wavenumber grid [cm^-1]
    REAL(REAL64), INTENT(OUT) :: crs(:) ! cross section [cm^-2]

    REAL(REAL64), PARAMETER, DIMENSION(3) :: a = [8.779d-13, 1.323d-6, 2.245d0] ! for H2
!   REAL(REAL64), PARAMETER, DIMENSION(3) :: a = [3.947d-12, 4.229d-6, 5.414d0] ! for N2
    REAL(REAL64), PARAMETER :: cm_to_angstrom = 1d8
    REAL(REAL64), PARAMETER :: lambda_max_h2 = 1215d0 ! maximum lambda in Victor & Dalgarno, which Tarafda & Vardya is based on
    REAL(REAL64) :: lambda
    INTEGER(INT32) :: k, ks, ke

    ks = 1
    ke = UBOUND(grd,1)

    crs = 0d0
    DO k = ks, ke
       lambda = 1d0 / grd(k) * cm_to_angstrom
       IF(lambda < lambda_max_h2) CYCLE

       crs(k) = a(1) / lambda**4 + a(2) / lambda**6 + a(3) / lambda**8
    END DO

    RETURN    
  END SUBROUTINE cross_rayleigh_h2_tarafdar_1973
    

  SUBROUTINE cross_rayleigh_h_lee_2005(grd, crs)
    ! Rayleigh scattering cross sections for H
    ! Lee (2005), MNRAS, 358, 1472
    USE const_module, ONLY : alpha, m_ele, hbar, clight, erg_to_ev, megabarn_to_cm2, pi
    IMPLICIT NONE
    REAL(REAL64), INTENT(IN) :: grd(:)    ! wavenumber grid [cm^-1]
    REAL(REAL64), INTENT(OUT) :: crs(:) ! cross section [cm^-2]

    REAL(REAL64), PARAMETER, DIMENSION(0:5) :: &
         c = [1.26563d0, 3.73828125d0, 8.813930935d0, 19.15379502d0, 39.92303232d0, 81.10881152d0]
    REAL(REAL64), PARAMETER :: omega_lyman = (alpha * clight)**2 * m_ele / (2d0 * hbar)
    REAL(REAL64), PARAMETER :: crs_thomson = 8d0 * pi / 3d0 * (alpha * hbar / (m_ele * clight))**2
    REAL(REAL64) :: omega, f, x
    INTEGER(INT32) :: k, ks, ke, n
    
    ks = 1
    ke = UBOUND(grd,1)

    crs = 0d0
    DO k = ks, ke
       omega = 2d0 * pi * clight * grd(k)
       x = omega / omega_lyman
       IF(x >= 0.75d0) CYCLE ! discard blueward of Lyman alpha

       f = 0d0
       DO n = LBOUND(c,1), UBOUND(c,1)
          f = f + c(n) * x**(2d0 * DBLE(n))
       END DO
       crs(k) = x**4 * f
    END DO

    crs = crs * crs_thomson

    RETURN    
  END SUBROUTINE cross_rayleigh_h_lee_2005
  
    
  SUBROUTINE cross_photo_hm_john_1988(grd, crs)
    ! Photoionization cross sections for H-
    ! John 1988, A&A, 193, 189
    USE const_module, ONLY : h, clight, ev_to_erg, megabarn_to_cm2, h_electron_aff_ev
    IMPLICIT NONE
    REAL(REAL64), INTENT(IN) :: grd(:)    ! wavenumber grid [cm^-1]
    REAL(REAL64), INTENT(OUT) :: crs(:)   ! cross section [cm^-2]

    REAL(REAL64), PARAMETER, DIMENSION(6) :: c = [152.519d0, 49.534d0, -118.858d0, 92.536d0, -34.194d0, 4.982d0]
    REAL(REAL64), PARAMETER :: lambda_min = 0.125d0
    REAL(REAL64), PARAMETER :: cm_to_micron = 1d4
    REAL(REAL64), PARAMETER :: lambda_h_electron_aff = h * clight / (h_electron_aff_ev * ev_to_erg) * cm_to_micron
    REAL(REAL64) :: lambda, f, x
    INTEGER(INT32) :: k, ks, ke, n
    
    ks = 1
    ke = UBOUND(grd,1)

    crs(:) = 0d0

    DO k = ks, ke
       lambda = 1d0 / grd(k) * cm_to_micron
       IF(lambda < lambda_min .OR. lambda > lambda_h_electron_aff) CYCLE
       f = 0d0
       x = 1d0 / lambda - 1d0 / lambda_h_electron_aff
       DO n = 1, 6
          f = f + c(n) * x**((DBLE(n) - 1d0)/2d0)
       END DO
       crs(k) = lambda**3 * x**1.5d0 * f
    END DO

    crs = crs * megabarn_to_cm2 ! megabarn -> cm^2

    RETURN    
  END SUBROUTINE cross_photo_hm_john_1988
  
    
  SUBROUTINE cross_photo_hm_ohmura_1962(grd, crs)
    ! Photoionization cross sections for H-
    ! Ohmura & Ohmura 1962, Physical Review, 118, 154
    USE const_module, ONLY : h, clight, ev_to_erg, h_electron_aff_ev, pi, m_ele, hbar, a_0
    IMPLICIT NONE
    REAL(REAL64), INTENT(IN) :: grd(:)    ! wavenumber grid [cm^-1]
    REAL(REAL64), INTENT(OUT) :: crs(:)   ! cross section [cm^-2]

    REAL(REAL64), PARAMETER :: gam = 0.2355883d0, rho = 2.646d0
    REAL(REAL64), PARAMETER :: lambda_h_electron_aff = h * clight / (h_electron_aff_ev * ev_to_erg)
    REAL(REAL64) :: k_e, lambda
    INTEGER(INT32) :: k, ks, ke
    
    ks = 1
    ke = UBOUND(grd,1)

    crs(:) = 0d0

    DO k = ks, ke
       lambda = 1d0 / grd(k)
       IF(lambda > lambda_h_electron_aff) CYCLE
       
       k_e = SQRT(4d0 * pi * m_ele * clight / hbar * (grd(k) - 1d0 / lambda_h_electron_aff))
       k_e = k_e * a_0
       crs(k) = 6.8475d-18 * gam * k_e**3 / (1d0 - gam * rho) / (gam**2 + k_e**2)**3
    END DO

    RETURN    
  END SUBROUTINE cross_photo_hm_ohmura_1962
  
    
  SUBROUTINE cross_photo_h2_yan_2001(grd, crs)
    ! Photoionization cross sections for H2
    ! Yan et al. 2001, ApJ, 559, 1194
    USE const_module, ONLY : h, clight, erg_to_ev, barn_to_cm2, h2_ion_ev
    IMPLICIT NONE
    REAL(REAL64), INTENT(IN) :: grd(:)    ! wavenumber grid [cm^-1]
    REAL(REAL64), INTENT(OUT) :: crs(:) ! cross section [cm^-2]

    REAL(REAL64) :: x, ene
    INTEGER(INT32) :: k, ks, ke

    ks = 1
    ke = UBOUND(grd,1)

    crs = 0d0
    DO k = ks, ke
       ene = h * clight * grd(k) * erg_to_ev
       x = ene / h2_ion_ev
       IF(ene < h2_ion_ev) THEN
          CYCLE
       ELSE IF(ene < 18d0) THEN
          crs(k) = 1d7 * (1d0 - 197.448d0 * x**(-0.5d0) + 438.823d0 * x**(-1d0) - 260.481d0 * x**(-1.5d0) + 17.915d0 * x**(-2d0))
       ELSE IF(ene < 30d0) THEN
          crs(k) = (-145.528d0 + 351.394d0 * x**0.5 - 274.294d0 * x + 74.320d0 * x**1.5d0) / (ene / 1d3)**3.5d0
       ELSE IF(ene < 85d0) THEN
          crs(k) = (65.304d0 - 91.762d0 * x**0.5d0 + 51.778d0 * x - 9.364d0 * x**1.5d0) / (ene / 1d3)**3.5d0
       ELSE
          crs(k) = 45.57d0 * (1d0 - 2.003d0 * x**(-0.5d0) - 4.806d0 * x**(-1d0) + 50.577d0 * x**(-1.5d0) &
               - 171.044d0 * x**(-2d0) + 231.608d0 * x**(-2.5d0) - 81.885d0 * x**(-3d0)) / (ene / 1d3)**3.5d0
       END IF
       crs(k) = MAX(crs(k), 0d0)
    END DO

    crs = crs * barn_to_cm2 ! barn -> cm^2

    RETURN    
  END SUBROUTINE cross_photo_h2_yan_2001
  
    
  SUBROUTINE cross_photo_he_yan_2001(grd, crs)
    ! Photoionization cross sections for H2
    ! Yan et al. 2001, ApJ, 559, 1194
    USE const_module, ONLY : h, clight, erg_to_ev, barn_to_cm2, helium_ion_wavenum
    IMPLICIT NONE
    REAL(REAL64), INTENT(IN) :: grd(:)    ! wavenumber grid [cm^-1]
    REAL(REAL64), INTENT(OUT) :: crs(:) ! cross section [cm^-2]

    REAL(REAL64), PARAMETER, DIMENSION(6) :: a = [-4.7416d0, 14.8200d0, -30.8678d0, 37.3584d0, -23.4585d0, 5.9133d0]
    REAL(REAL64), PARAMETER :: helium_ion_ev = helium_ion_wavenum * h * clight * erg_to_ev
    REAL(REAL64) :: x, ene, f
    INTEGER(INT32) :: k, ks, ke, n

    ks = 1
    ke = UBOUND(grd,1)

    crs = 0d0
    DO k = ks, ke
       ene = h * clight * grd(k) * erg_to_ev
       x = ene / helium_ion_ev
       IF(x < 1d0) CYCLE
       f = 1d0
       DO n = 1, 6
          f = f + a(n) / x ** (DBLE(n) / 2d0)
       END DO
       crs(k) = 733d0 / (ene / 1d3)**3.5d0 * f
    END DO

    crs = crs * barn_to_cm2 ! barn -> cm^2

    RETURN    
  END SUBROUTINE cross_photo_he_yan_2001

  
  SUBROUTINE cross_photo_verner(grd, crs, z0, ne0)
#define VERNER_ALLSPECIES
    ! Photoionization cross sections of atoms and ions in the ground state
    ! Verner et al. (1995) A&AS, 109, 125, Table 1
    ! Verner et al. (1996) ApJ, 465, 487, Table 1
    ! http://www.pa.uky.edu/~verner/photo.html
    USE code_module, ONLY : HI, HeI, HeII, CI, NI, OI, NaI, MgI, MgII, AlI, SiI, SiII, SI, CaI, CaII, FeI
    USE const_module, ONLY : megabarn_to_cm2, wn_to_ev
    IMPLICIT NONE
    REAL(REAL64), INTENT(IN) :: grd(:)    ! wavenumber grid [cm^-1]
    REAL(REAL64), INTENT(OUT) :: crs(:) ! cross section [cm^-2]
    INTEGER(INT32), INTENT(IN) :: z0 ! atomic number
    INTEGER(INT32), INTENT(IN) :: ne0 ! electron number
    
    INTEGER(INT32) :: z, ne, n, l, iostat, ks, ke, k
#ifndef VERNER_ALLSPECIES    
    INTEGER(INT32) :: code
#endif    
    REAL(REAL64) :: eth, e0, sigma0, ya, p, yw, emax, y0, y1, e, f, q, x, y
    REAL(REAL64), ALLOCATABLE :: crs95(:)

    crs(:) = 0d0
 
    ks = 1
    ke = UBOUND(grd,1)

    ALLOCATE(crs95(ks:ke))

    ! Verner et al. (1995) A&AS, 109, 125, Table 1
    crs95(:) = 0d0
    OPEN(1, FILE='input/photo/table1.dat', STATUS='OLD')
    DO
       READ(1,*,IOSTAT=iostat) z, ne, n, l, eth, e0, sigma0, ya, p, yw
       IF(iostat /= 0) EXIT
       IF(z == z0 .AND. ne == ne0) THEN
          DO k = ks, ke
#ifndef VERNER_ALLSPECIES
             code = z * 100 + (z - ne)
             IF(&
                  code == HI .OR. &
                  code == HeI .OR. &
                  code == HeII .OR. &
                  code == CI .OR. &
                  code == NI .OR. &
                  code == OI .OR. &
                  code == NaI .OR. &
                  code == MgI .OR. &
                  code == MgII .OR. &
                  code == AlI .OR. &
                  code == SiI .OR. &
                  code == SiII .OR. &
                  code == SI .OR. &
                  code == CaI .OR. &
                  code == CaII .OR. &
                  code == FeI) THEN
                CYCLE
             END IF
#endif       
             e = grd(k) * wn_to_ev
             IF(e < eth) CYCLE
             y = e / e0
             q = 5.5d0 + l - 0.5d0 * p
             f = ((y - 1d0)**2 + yw**2) * y ** (-q) * (1d0 + SQRT(y / ya)) ** (-p)
             crs95(k) = crs95(k) + (sigma0 * f)
          END DO
       END IF
    END DO
    CLOSE(1)
#ifdef VERNER95_ONLY
    crs = crs95
#else
    ! Verner et al. (1996) ApJ, 465, 487, Table 1
    OPEN(2, FILE='input/photo/photo.dat', STATUS='OLD')
    DO
       READ(2,*,IOSTAT=iostat) z, ne, eth, emax, e0, sigma0, ya, p, yw, y0, y1
       IF(iostat /= 0) EXIT
       IF(z == z0 .AND. ne == ne0) THEN
          DO k = ks, ke
#ifndef VERNER_ALLSPECIES
             code = z * 100 + (z - ne)
             IF(&
                  code == HI .OR. &
                  code == HeI .OR. &
                  code == HeII .OR. &
                  code == CI .OR. &
                  code == NI .OR. &
                  code == OI .OR. &
                  code == NaI .OR. &
                  code == MgI .OR. &
                  code == MgII .OR. &
                  code == AlI .OR. &
                  code == SiI .OR. &
                  code == SiII .OR. &
                  code == SI .OR. &
                  code == CaI .OR. &
                  code == CaII .OR. &
                  code == FeI) THEN
                CYCLE
             END IF
#endif       
             e = grd(k) * wn_to_ev
             IF(e < eth) THEN
                CYCLE
#ifndef VERNER96_ONLY                
             ELSE IF(e > emax) THEN ! adopt Verner 1995 for inner shells
                crs(k) = crs95(k)
#endif                
             ELSE ! adopt Verner 1996 for outer shells
                x = e / e0 - y0
                y = SQRT(x**2 + y1**2)
                q = 5.5d0 - 0.5d0 * p
                f = ((x - 1d0)**2 + yw**2) * y ** (-q) * (1d0 + SQRT(y / ya)) ** (-p)
                crs(k) = (sigma0 * f)
             END IF
          END DO
       END IF
    END DO
    CLOSE(2)
#endif
    DEALLOCATE(crs95)

    crs(:) = crs(:) * megabarn_to_cm2 ! megabarn -> cm^2
    
    RETURN
    
  END SUBROUTINE cross_photo_verner

  
  SUBROUTINE compute_frac_g(temp, frac_g, zeta)
    ! compute the fraction of the ground state configuration
    USE ISO_C_BINDING
    USE HDF5
    USE H5LT
    USE const_module, ONLY : c2
    USE code_module, ONLY : ZnI
    IMPLICIT NONE
    REAL(REAL64), INTENT(IN) :: temp(:)
    REAL(REAL64), INTENT(INOUT) :: frac_g(:,:)
    REAL(REAL64), INTENT(INOUT) :: zeta(:,:)
    
    INTEGER(SIZE_T), PARAMETER :: sdim = 8, sdim2 = 32
    INTEGER(HID_T) :: file_id, grp_prop, grp_ion, data_term, data_conf, space_id, memtype
    INTEGER(HSIZE_T), DIMENSION(1:1) :: dims

    INTEGER(INT64) :: l
    INTEGER(INT32) :: nstates, error, z, ne, j, code
    INTEGER(INT32), ALLOCATABLE :: gtot(:)
    REAL(REAL64), ALLOCATABLE :: ene(:)
    REAL(REAL64) :: q, q_g, dq
    CHARACTER(LEN=sdim), ALLOCATABLE, TARGET :: term(:)
    CHARACTER(LEN=sdim) :: term_ground
    CHARACTER(LEN=sdim2), ALLOCATABLE, TARGET :: conf(:)
    CHARACTER(LEN=sdim2) :: conf_ground
    CHARACTER :: natm*3, nion*2
    TYPE(c_ptr) :: f_ptr

    CALL h5open_f(error)
    CALL h5Fopen_f('input/h5/NIST.h5', H5F_ACC_RDONLY_F, file_id, error)
    CALL h5Gopen_f(file_id, 'prop', grp_prop, error)

    DO z = 1, (ZnI/100)
       DO ne = 1, z
          code = 100 * z + (z - ne)
          WRITE(natm,'(I3.3)') z
          WRITE(nion,'(I2.2)') z - ne
          CALL h5Gopen_f(grp_prop, natm//'.'//nion, grp_ion, error)
          dims = [1]
          CALL h5LTread_dataset_f(grp_ion, 'nstates', H5T_STD_I32LE , nstates, dims, error)
          ALLOCATE(ene(nstates), gtot(nstates), term(nstates), conf(nstates))
          dims = [nstates]
          CALL h5LTread_dataset_f(grp_ion, 'ene' , H5T_IEEE_F64LE, ene , dims, error)
          CALL h5LTread_dataset_f(grp_ion, 'gtot', H5T_STD_I32LE, gtot, dims, error)
          !---- read an array of strings ---
          ! https://support.hdfgroup.org/ftp/HDF5/examples/examples-by-api/hdf5-examples/1_8/FORTRAN/H5T/h5ex_t_stringC_F03.f90
          CALL h5Dopen_f(grp_ion, 'term', data_term, error)
          CALL H5Dget_space_f(data_term, space_id, error)
          CALL H5Tcopy_f(H5T_FORTRAN_S1, memtype, error)
          CALL H5Tset_size_f(memtype, sdim, error)
          f_ptr = C_LOC(term(1)(1:1))
          CALL H5Dread_f(data_term, memtype, f_ptr, error, space_id)
          CALL H5Dclose_f(data_term, error)
          CALL H5Sclose_f(space_id, error)
          CALL H5Tclose_f(memtype, error)
          !----
          CALL h5Dopen_f(grp_ion, 'conf', data_conf, error)
          CALL H5Dget_space_f(data_conf, space_id, error)
          CALL H5Tcopy_f(H5T_FORTRAN_S1, memtype, error)
          CALL H5Tset_size_f(memtype, sdim2, error)
          f_ptr = C_LOC(conf(1)(1:1))
          CALL H5Dread_f(data_conf, memtype, f_ptr, error, space_id)
          CALL H5Dclose_f(data_conf, error)
          CALL H5Sclose_f(space_id, error)
          CALL H5Tclose_f(memtype, error)
          !----
          
          DO j = 1, UBOUND(temp,1)
             term_ground = term(1)
             conf_ground = conf(1)
             q = 0d0
             q_g = 0d0
             DO l = 1, UBOUND(ene,1)
                dq = gtot(l) * EXP(-c2 * ene(l) / temp(j))
                q = q + dq
                IF(term(l) == term_ground .AND. conf(l) == conf_ground) THEN
                   q_g = q_g + dq
                END IF
             END DO
!#define FRACG_INVERSE_Q             
#ifndef FRACG_INVERSE_Q             
             frac_g(j,code) = q_g / q
#else
             frac_g(j,code) = 1d0 / q
#endif             
             zeta(j,code) = q
          END DO
          
          DEALLOCATE(ene, gtot, term, conf)
          CALL h5Gclose_f(grp_ion, error)
       END DO
    END DO
    
    CALL h5Gclose_f(grp_prop, error)
    CALL h5Fclose_f(file_id, error)
    CALL h5close_f(error)

    RETURN
  END SUBROUTINE compute_frac_g


  SUBROUTINE get_eneion(code, wn_ion)
    ! retrieve ionization energy of atom/ion from NIST data
    USE HDF5
    USE H5LT
    USE const_module, ONLY : c2
    IMPLICIT NONE
    INTEGER(INT32), INTENT(IN) :: code
    REAL(REAL64), INTENT(OUT) :: wn_ion
    
    INTEGER(HID_T) :: file_id, grp_prop, grp_ion
    INTEGER(HSIZE_T), DIMENSION(1:1) :: dims

    INTEGER(INT32) :: error
    CHARACTER :: natm*3, nion*2

    CALL h5open_f(error)
    CALL h5Fopen_f('input/h5/NIST.h5', H5F_ACC_RDONLY_F, file_id, error)
    CALL h5Gopen_f(file_id, 'prop', grp_prop, error)
    WRITE(natm,'(I3.3)') code / 100
    WRITE(nion,'(I2.2)') code - (code / 100) * 100
    CALL h5Gopen_f(grp_prop, natm//'.'//nion, grp_ion, error)
    dims = [1]
    CALL h5LTread_dataset_f(grp_ion, 'eneion', H5T_IEEE_F64LE , wn_ion, dims, error)
    CALL h5Gclose_f(grp_ion, error)
    CALL h5Gclose_f(grp_prop, error)
    CALL h5Fclose_f(file_id, error)
    CALL h5close_f(error)
    
    RETURN
  END SUBROUTINE get_eneion


  SUBROUTINE interp2d_yarr(x0, y0, x, y, f, f0)
    ! bilinear interpolation
    IMPLICIT NONE
    REAL(REAL64), INTENT(IN) :: x0
    REAL(REAL64), INTENT(IN) :: y0(:)
    REAL(REAL64), INTENT(IN) :: x(:)
    REAL(REAL64), INTENT(IN) :: y(:)
    REAL(REAL64), INTENT(IN) :: f(:,:)
    REAL(REAL64), INTENT(OUT) :: f0(:)

    INTEGER(INT32) :: i, j, k
    REAL(REAL64) :: x_min, x_max, y_min, y_max, f_u, f_d

    f0 = 0d0
    x_min = MINVAL(x)
    x_max = MAXVAL(x)
    y_min = MINVAL(y)
    y_max = MAXVAL(y)

    IF(x0 >= x_max .OR. x0 < x_min) RETURN

    DO i = 1, UBOUND(x,1)
       IF(x(i) > x0) EXIT
    END DO

    DO k = 1, UBOUND(y0,1)
       IF(y0(k) >= y_max .OR. y0(k) < y_min) CYCLE

       DO j = 1, UBOUND(y,1)
          IF(y(j) > y0(k)) EXIT
       END DO
       f_u = ((x(i) - x0) * f(i-1,j  ) + (x0 - x(i-1)) * f(i,j  )) / (x(i) - x(i-1))
       f_d = ((x(i) - x0) * f(i-1,j-1) + (x0 - x(i-1)) * f(i,j-1)) / (x(i) - x(i-1))
       f0(k) = ((y(j) - y0(k)) * f_d + (y0(k) - y(j-1)) * f_u) / (y(j) - y(j-1))
    END DO

    RETURN
  END SUBROUTINE interp2d_yarr

  
END MODULE continuum_module
