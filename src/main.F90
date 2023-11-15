
#define NA -1

PROGRAM main
  
  USE ISO_FORTRAN_ENV
  USE MPI
  USE code_module, ONLY : HI, HeI, HeII, CI, NI, OI, NaI, MgI, MgII, AlI, SiI, SiII, SI, CaI, CaII, FeI, ELECTRON, H2, ZnI, TOTAL
  USE const_module, ONLY : k_bol, c2, alpha, hbar, pi, clight, m_ele, crs_thomson, h, ev_to_erg, h_electron_aff_ev
  USE continuum_module, ONLY : compute_frac_g, &
       brems_atom_vanHoof_2014, &
       brems_atom_ferland, &
       brems_hm_john_1988, &
       brems_h2m_john_1978, &
       cross_photo_hm_john_1988, &
       cross_photo_hm_ohmura_1962, &
       cross_photo_verner, &
       cross_photo_h2_yan_2001, &
       cross_rayleigh_h_lee_2005, &
       cross_rayleigh_h_rohrmann_2022, &
       cross_rayleigh_he_tarafdar_1969, &
       cross_rayleigh_he_rohrmann_2018, &
       cross_rayleigh_h2_tarafdar_1973, &
       coeff_cia, &
       photo_mathisen, photo_topbase
  USE input_module, ONLY : read_mol_source_list, read_grid, read_eos
  USE output_module, ONLY : output_table, output_mono
  USE line_module, ONLY : line_molec, line_kurucz, line_kurucz_p
  USE mpi_module, ONLY : mpi_first, myrk, para_range, myrk_j, jprc, myrk_jconst, unit
  USE mean_module, ONLY : mean_opac, mean_opac_lambda
  
  IMPLICIT NONE

  REAL(REAL64), ALLOCATABLE :: np(:,:) ! partial density (species,layer)
  REAL(REAL64), ALLOCATABLE :: temp(:) ! temperature (layer)
  REAL(REAL64), ALLOCATABLE :: rho(:)  ! mass density (layer)
  REAL(REAL64), ALLOCATABLE :: grd(:) ! grid (wavenumber)
  REAL(REAL64), ALLOCATABLE :: dgrd(:) ! grid (wavenumber)
  REAL(REAL64), ALLOCATABLE :: alp(:,:) ! continuum opacity (wavenumber,layer)
  REAL(REAL64), ALLOCATABLE :: line(:,:) ! line opacity (wavenumber,layer)
  REAL(REAL64), ALLOCATABLE :: sca(:,:) ! scattering opacity (wavenumber,layer)
  REAL(REAL64), ALLOCATABLE :: out(:,:,:) ! scattering opacity (wavenumber,layer)
  REAL(REAL64), ALLOCATABLE :: pmean(:) ! Planck-mean opacity
  REAL(REAL64), ALLOCATABLE :: pmean2(:) ! two-temp Planck-mean
  REAL(REAL64), ALLOCATABLE :: p_e(:)
  REAL(REAL64), ALLOCATABLE :: frac_g(:,:)
  REAL(REAL64), ALLOCATABLE :: zeta(:,:)
  INTEGER(INT64), ALLOCATABLE :: nlines(:,:)
  INTEGER(INT32), ALLOCATABLE :: map(:)

  REAL(REAL64) :: fac, wn_ion, fac_0, temp2
  INTEGER(INT32) :: j, jmax, iostat, code, z, ne, n_species, js, je, error, ns
  INTEGER(INT32) :: ks, ke
  CHARACTER(LEN=64), ALLOCATABLE :: source(:)
  INTEGER(INT32) :: jj, iblock, count
  CHARACTER(LEN=4) :: id
  REAL(REAL64), ALLOCATABLE :: temp1(:), p_e1(:), np1(:,:), rho1(:)
  INTEGER(INT32) :: &
       line_molecules = 0, &
       line_kurucz_phoenix = 0, &
       line_kurucz_gfall = 0, &
       line_kurucz_gfpred = 0, & 
       rayleigh_scattering_h2 = 0, &
       rayleigh_scattering_he = 0, &
       rayleigh_scattering_h = 0, &
       electron_scattering = 0, &
       cia = 0, &
       photoion_h2 = 0, &
       photoion_verner = 0, &
       photoion_h_minus = 0, &
       brems_h2_minus = 0, &
       brems_h_minus = 0, &
       brems_atomicions = 0, &
       photoion_mathisen = 0, &
       photoion_topbase = 0
  NAMELIST /switches/ &
       line_molecules, &
       line_kurucz_phoenix, &
       line_kurucz_gfall, &
       line_kurucz_gfpred, & 
       rayleigh_scattering_h2, &
       rayleigh_scattering_he, &
       rayleigh_scattering_h, &
       electron_scattering, &
       cia, &
       photoion_h2, &
       photoion_verner, &
       photoion_h_minus, &
       brems_h2_minus, &
       brems_h_minus, &
       brems_atomicions, &
       photoion_mathisen, &
       photoion_topbase
  NAMELIST /radtemp/ temp2
  NAMELIST /block_cyclic/ iblock
  
  CALL MPI_INIT(error)

  CALL mpi_first

  OPEN(5, FILE='input/fort.5', STATUS='OLD')
  READ(5, NML=switches)
  READ(5, NML=radtemp)
  READ(5, NML=block_cyclic)
  CLOSE(5)

  ! SET MOLECULAR LINE SOURCES
  IF(myrk == 0) THEN
     CALL read_mol_source_list(n_species)
  END IF
  CALL MPI_BCAST(n_species, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, error)
  ALLOCATE(source(n_species))
  OPEN(1, FILE='input/mol_source.dat', STATUS='OLD')
  DO ns = 1, n_species
     READ(1,*,IOSTAT=iostat) source(ns)
  END DO
  CLOSE(1)

  ! SET GRID
  CALL read_grid(grd, dgrd)
  ks = LBOUND(grd,1)
  ke = UBOUND(grd,1)

  ! SET PARTIAL PRESSURES
  CALL read_eos(np, temp, rho)
  jmax = UBOUND(np,2)
  ALLOCATE(p_e(jmax))
  ! electron pressure
#ifndef CROSS_CHECK
  p_e(:) = np(ELECTRON,:) * k_bol * temp(:)
#else
  p_e(:) = 1d0
#endif

  ! CYCLIC PARALLEL
  count = 0
  DO j = 1 + myrk_j, jmax, jprc
     count = count + 1
  END DO
  ALLOCATE(temp1(count), p_e1(count), np1(UBOUND(np,1),count), rho1(count), map(count))

  count = 0
  WRITE(id,'(I4)') unit
  OPEN(unit+3000,FILE='output/parallel_'//id//'.data',STATUS='UNKNOWN')
  DO j = 1 + myrk_j, jmax, jprc
     count = count + 1
     temp1(count) = temp(j)
     p_e1(count) = p_e(j)
     np1(:,count) = np(:,j)
     rho1(count) = rho(j)
     map(count) = j
     WRITE(unit+3000,'(2(I3,x),2(F8.2,x))') count, j, LOG10(temp(j)), LOG10(temp(j)*k_bol*np(TOTAL,j))
  END DO
  CLOSE(unit+3000)

  LAYER: DO jj = 0, iblock-1

     ! BLOCK PARALLEL
     CALL para_range(1, count, iblock, jj, js, je)

     ALLOCATE(sca(ks:ke,js:je), alp(ks:ke,js:je), line(ks:ke,js:je), frac_g(js:je,UBOUND(np,1)), &
          zeta(js:je,UBOUND(np,1)), pmean(js:je), pmean2(js:je))
     ALLOCATE(nlines(js:je,4))

     alp(:,:) = 0d0
     line(:,:) = 0d0
     sca(:,:) = 0d0
     pmean(:) = 0d0
     pmean2(:) = 0d0

     ! *********************
     ! CONTINUUM  ABSORPTION
     ! *********************
     
     call wtime(name='ground state frac:')
     CALL compute_frac_g(temp1(js:je), frac_g(js:je,:), zeta(js:je,:))
     call wtime()

     ! BREMSSTRAHLUNG
     
     IF(brems_h_minus == 1) THEN
        call wtime(name='brems(H-):')
        ALLOCATE(out(ks:ke,js:je,NA:NA))
        CALL brems_hm_john_1988(temp1(js:je), grd(:), out(:,js:je,NA))
        DO j = js, je
           alp(:,j) = alp(:,j) + (p_e1(j) * np1(HI,j)) * out(:,j,NA)
        END DO
        DEALLOCATE(out)
        call wtime()
     END IF
     
     IF(brems_h2_minus == 1) THEN
        call wtime(name='brems(H2-):')
        ALLOCATE(out(ks:ke,js:je,NA:NA))
        CALL brems_h2m_john_1978(temp1(js:je), grd(:), out(:,js:je,NA))
        DO j = js, je
           alp(:,j) = alp(:,j) + (p_e1(j) * np1(H2,j)) * out(:,j,NA)
        END DO
        DEALLOCATE(out)
        call wtime()
     END IF
     
     IF(brems_atomicions == 1) THEN
#ifndef PHOENIX_BREMS
        call wtime(name='brems(atomic ion):')
        ALLOCATE(out(ks:ke,js:je,NA:NA))
        fac_0 = SQRT(8d0 * alpha**6 * hbar**4 / (27d0 * pi * clight**2 * k_bol**3 * m_ele**3))
        DO z = 1, (ZnI/100)
           CALL brems_atom_vanHoof_2014(temp1(js:je), grd(:), out(:,js:je,NA), z - ne)
           DO ne = 0, z - 1
              code = z * 100 + (z - ne)
              DO j = js, je
                 alp(:,j) = alp(:,j) + p_e1(j) * np1(code,j) * (z - ne)**2 / temp1(j)**1.5d0 * fac_0 * &
                      (1d0 - EXP(-c2 * grd(:) / temp1(j))) * out(:,j,NA) / grd(:)**3
              END DO
           END DO
        END DO
        DEALLOCATE(out)
        call wtime()
#else        
        call wtime(name='brems(H only, ferland):')
        ALLOCATE(out(ks:ke,js:je,1:(ZnI/100)))
        CALL brems_atom_ferland(temp1(js:je), grd(:), out(:,js:je,:))
        DO z = 1, 1!(ZnI/100)
           DO ne = 0, z - 1
              code = z * 100 + (z - ne)
              DO j = js, je
                 alp(:,j) = alp(:,j) + np1(code,j) * np1(ELECTRON,j) * z**2 / &
                      SQRT(temp1(j)) * 1.37040445d-47 * &
                      (1d0 - EXP(-c2 * grd(:) / temp1(j))) * out(:,j,z) * (1d8 / grd(:))**3
              END DO
           END DO
        END DO
        DEALLOCATE(out)
        call wtime()
#endif        
     END IF

     ! PHOTOIONIZATION
     
     IF(photoion_verner == 1) THEN
        call wtime(name='photo(Verner):')
        ALLOCATE(out(ks:ke,NA:NA,NA:NA))
        DO z = 1, (ZnI/100)
           DO ne = 1, z
              code = z * 100 + (z - ne)
              CALL cross_photo_verner(grd(:), out(:,NA,NA), z, ne)
              DO j = js, je
                 alp(ks:ke,j) = alp(ks:ke,j) + (np1(code,j) * frac_g(j,code)) * out(ks:ke,NA,NA) &
                      * (1d0 - EXP(-c2 * grd(ks:ke) / temp1(j)))
              END DO
           END DO
        END DO
        DEALLOCATE(out)
        call wtime()
     END IF
     
     IF(photoion_mathisen == 1) THEN
        call wtime(name='photo(Mathisen):')
        ALLOCATE(out(ks:ke,js:je,NA:NA))
        CALL photo_mathisen(temp1(js:je), grd(:), np1(:,js:je), zeta(js:je,:), out(:,js:je,NA), photoion_verner)
        DO j = js, je
           alp(:,j) = alp(:,j) + out(:,j,NA)
        END DO
        DEALLOCATE(out)
        call wtime()
     END IF
     
     IF(photoion_topbase == 1) THEN
        call wtime(name='photo(TOPbase):')
        ALLOCATE(out(ks:ke,js:je,NA:NA))
        CALL photo_topbase(temp1(js:je), grd(:), np1(:,js:je), out(:,js:je,NA), photoion_verner)
        DO j = js, je
           alp(:,j) = alp(:,j) + out(:,j,NA)
        END DO
        DEALLOCATE(out)
        call wtime()
     END IF
     
     IF(photoion_h_minus == 1) THEN
        call wtime(name='photo(H-):')
        ALLOCATE(out(ks:ke,NA:NA,NA:NA))
!!$        CALL cross_photo_hm_john_1988(grd(:), out(:,NA,NA))
        CALL cross_photo_hm_ohmura_1962(grd(:), out(:,NA,NA))
        wn_ion = h_electron_aff_ev * ev_to_erg / (h * clight)
        DO j = js, je
           fac = p_e1(j) * np1(HI,j) * EXP(c2 * wn_ion / temp1(j)) * 0.750d0 * temp1(j)**(-2.5d0)
           fac = fac * (2d0 / zeta(j,HI)) ! correction for implicit assumption of Z_HI = 2
           alp(:,j) = alp(:,j) + fac * (1d0 - EXP(-c2 * grd(:) / temp1(j))) * out(:,NA,NA)
        END DO
        DEALLOCATE(out)
        call wtime()
     END IF
        
     IF(photoion_h2 == 1) THEN
        call wtime(name='photo(H2):')
        ALLOCATE(out(ks:ke,NA:NA,NA:NA))
        CALL cross_photo_h2_yan_2001(grd(:), out(:,NA,NA))
        DO j = js, je
           alp(:,j) = alp(:,j) + np1(H2,j) * out(:,NA,NA) * (1d0 - EXP(-c2 * grd(ks:ke) / temp1(j)))
        END DO
        DEALLOCATE(out)
        call wtime()
     END IF
     
     ! COLLISION-INDUCED ABSORPTION (CIA)
     IF(cia == 1) THEN
        call wtime(name='CIA:')
        ALLOCATE(out(ks:ke,js:je,NA:NA))
        CALL coeff_cia(temp1(js:je), np1(:,js:je), grd(:), out(:,js:je,NA))
        DO j = js, je
           alp(:,j) = alp(:,j) + out(:,j,NA)
        END DO
        DEALLOCATE(out)
        call wtime()
     END IF
     
     ! **********
     ! SCATTERING
     ! **********

     IF(electron_scattering == 1) THEN    
        call wtime(name='electron sca:')
        DO j = js, je
           sca(:,j) = sca(:,j) + np1(ELECTRON,j) * crs_thomson
        END DO
        call wtime()
     END IF
     
     IF(rayleigh_scattering_h == 1) THEN
        call wtime(name='Rayleigh sca(H):')
        ALLOCATE(out(ks:ke,NA:NA,NA:NA))
!!$        CALL cross_rayleigh_h_lee_2005(grd(:), out(:,NA,NA))
        CALL cross_rayleigh_h_rohrmann_2022(grd(:), out(:,NA,NA))
        DO j = js, je
           sca(:,j) = sca(:,j) + (np1(HI,j) * frac_g(j, HI)) * out(:,NA,NA)
        END DO
        DEALLOCATE(out)
        call wtime()
     END IF
     
     IF(rayleigh_scattering_he == 1) THEN
        call wtime(name='Rayleigh sca(He):')
        ALLOCATE(out(ks:ke,NA:NA,NA:NA))
!!$        CALL cross_rayleigh_he_tarafdar_1969(grd(:), out(:,NA,NA))
        CALL cross_rayleigh_he_rohrmann_2018(grd(:), out(:,NA,NA))
        DO j = js, je
           sca(:,j) = sca(:,j) + (np1(HeI,j) * frac_g(j,HeI)) * out(:,NA,NA)
        END DO
        DEALLOCATE(out)
        call wtime()
     END IF
     
     IF(rayleigh_scattering_h2 == 1) THEN
        call wtime(name='Rayleigh sca(H2):')
        ALLOCATE(out(ks:ke,NA:NA,NA:NA))
        CALL cross_rayleigh_h2_tarafdar_1973(grd(:), out(:,NA,NA))
        DO j = js, je
           sca(:,j) = sca(:,j) + np1(H2,j) * out(:,NA,NA)
        END DO
        DEALLOCATE(out)
        call wtime()
     END IF

     ! ***************
     ! LINE ABSORPTION
     ! ***************
     IF(line_molecules == 1) THEN
        call wtime(name='molecular lines:',nlines=nlines)
        DO ns = 1, n_species
           ALLOCATE(out(ks:ke,js:je,NA:NA))
           CALL line_molec(source(ns), temp1(js:je), np1(:,js:je), grd(:), dgrd(:), out(:,js:je,NA), &
                alp(:,js:je)+sca(:,js:je), pmean(js:je), temp2, pmean2(js:je), nlines(js:je,:))
           line(:,js:je) = line(:,js:je) + out(:,js:je,NA)
           DEALLOCATE(out)
        END DO
        call wtime(nlines=nlines)
     END IF

     IF(line_kurucz_phoenix == 1) THEN
        call wtime(name='Kurucz phoenix:',nlines=nlines)
        ALLOCATE(out(ks:ke,js:je,NA:NA))
        CALL line_kurucz_p('kurucz_phoenix', temp1(js:je), np1(:,js:je), grd(:), dgrd(:), out(:,js:je,NA), &
             alp(:,js:je)+sca(:,js:je), pmean(js:je), temp2, pmean2(js:je), nlines(js:je,:))
        line(:,js:je) = line(:,js:je) + out(:,js:je,NA)
        DEALLOCATE(out)
        call wtime(nlines=nlines)
     END IF
     
     IF(line_kurucz_gfall == 1) THEN
        call wtime(name='Kurucz gfall:',nlines=nlines)
        ALLOCATE(out(ks:ke,js:je,NA:NA))
        CALL line_kurucz('gfall08oct17', temp1(js:je), np1(:,js:je), grd(:), dgrd(:), out(:,js:je,NA), &
             alp(:,js:je)+sca(:,js:je), pmean(js:je), temp2, pmean2(js:je), nlines(js:je,:))
        line(:,js:je) = line(:,js:je) + out(:,js:je,NA)
        DEALLOCATE(out)
        call wtime(nlines=nlines)
     END IF
     
     IF(line_kurucz_gfpred == 1) THEN
        call wtime(name='Kurucz gfpred:',nlines=nlines)
        ALLOCATE(out(ks:ke,js:je,NA:NA))
        CALL line_kurucz('gfpred26apr18', temp1(js:je), np1(:,js:je), grd(:), dgrd(:), out(:,js:je,NA), &
             alp(:,js:je)+sca(:,js:je), pmean(js:je), temp2, pmean2(js:je), nlines(js:je,:))
        line(:,js:je) = line(:,js:je) + out(:,js:je,NA)
        DEALLOCATE(out)
        call wtime(nlines=nlines)
     END IF 
     ! ----------------
     ! output monochromatic opacities
     ! compute mean opacities
     DO j = js, je
        IF(myrk_jconst == 0) THEN
           CALL output_mono(map(j), temp1(j), np1(:,j), grd(:), alp(:,j), sca(:,j), line(:,j), temp2, rho1(j), &
                pmean(j), pmean2(j))
        END IF
     END DO
     DEALLOCATE(alp, sca, line, frac_g, zeta, nlines, pmean, pmean2)
     
  END DO layer

  !  DEALLOCATE(np, temp, grd, dgrd, mean_p, mean_r, mean_p2)
  DEALLOCATE(temp1, p_e1, np1, rho1, map)
  DEALLOCATE(np, temp, grd, dgrd)
  IF(line_molecules == 1 .AND. n_species /= 0) THEN
     DEALLOCATE(source)
  END IF
     
  CALL MPI_FINALIZE(error)

  STOP

CONTAINS

  SUBROUTINE wtime(name,nlines)
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: name
    INTEGER(INT64), OPTIONAL, INTENT(INOUT) :: nlines(:,:)
    REAL(REAL64), SAVE :: t0
    CHARACTER(LEN=4) :: id
    
    WRITE(id,'(I4)') unit
    IF(PRESENT(name)) THEN
       OPEN(unit,FILE='output/timing_'//id//'.data',STATUS='UNKNOWN',POSITION='APPEND')
       WRITE(unit,FMT='(A24)',ADVANCE='NO') name
       t0 = MPI_WTIME()
       IF(PRESENT(nlines)) THEN
          OPEN(unit+1000,FILE='output/nlines_'//id//'.data',STATUS='UNKNOWN',POSITION='APPEND')
          WRITE(unit+1000,FMT='(A24)',ADVANCE='NO') name
          nlines(:,:) = 0
       END IF
    ELSE
       WRITE(unit,FMT='(F8.2,A)') MPI_WTIME() - t0, ' [s]'
       CLOSE(unit)
       IF(PRESENT(nlines)) THEN
          WRITE(unit+1000,FMT='(4I13)') SUM(nlines(:,1)), SUM(nlines(:,2)), SUM(nlines(:,3)), SUM(nlines(:,4))
       END IF
       CLOSE(unit+1000)
    END IF

  END SUBROUTINE wtime
  
END PROGRAM main
