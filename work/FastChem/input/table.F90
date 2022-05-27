PROGRAM table
  IMPLICIT NONE
  REAL(8) :: t, p, tmin, tmax, dt, pmin, pmax, dp
  INTEGER :: i, j, imax, jmax

  ! pressure unit
  ! 1 bar = 10^6 Ba

  ! temperature in theta
  imax = 11
  tmin = 2d0  ! T = 5040/2 = 2520 K
  tmax = 50d0 ! T = 5040/50= 100.8K

  ! pressure in Ba
  jmax = 11
  pmin = -6d0 ! P = 10^-12 bar
  pmax = 9d0  ! P = 10^  3 bar

  ! temperature
  ! imax = 8
  ! tmin = 3d0
  ! tmax = 4.26d0

  ! pressure in Ba
  ! jmax = 8
  ! pmin = -8d0
  ! pmax = 11.6d0

  dt = (tmax - tmin) / (imax - 1)
  dp = (pmax - pmin) / (jmax - 1)
  
  OPEN(1,FILE='table.dat',STATUS='REPLACE')
  WRITE(1,*) '# table.dat'
  WRITE(1,*) '# pressure in bar, temperature in K'
  DO i = 1, imax
     t = tmin + (i - 1) * dt
     DO j = 1, jmax
        p = pmin + (j - 1) * dp
        WRITE(1,*) 10d0**p / 1d6, 5040d0/t
     END DO
  END DO
  CLOSE(1)  

END PROGRAM table
