PROGRAM table
  IMPLICIT NONE
  REAL(8) :: t, p, tmin, tmax, dt, pmin, pmax, dp
  INTEGER :: i, j, imax, jmax

  imax = 11
  tmin = 2d0
  tmax = 50d0

  jmax = 11
  pmin = -6d0
  pmax = 9d0

  ! imax = 8
  ! tmin = 3d0
  ! tmax = 4.26d0

  ! jmax = 8
  ! pmin = -8d0
  ! pmax = 11.6d0

  dt = (tmax - tmin) / (imax - 1)
  dp = (pmax - pmin) / (jmax - 1)
  
  OPEN(1,FILE='table.dat',STATUS='REPLACE')
  WRITE(1,*) '# table.dat'
  DO i = 1, imax
     t = tmin + (i - 1) * dt
     DO j = 1, jmax
        p = pmin + (j - 1) * dp
        PRINT *, i, j, t, p
!        WRITE(1,*) 10d0**t, 10d0**p / 1d6
        WRITE(1,*) 5040d0/t, 10d0**p / 1d6
     END DO
  END DO
  CLOSE(1)  

END PROGRAM table
